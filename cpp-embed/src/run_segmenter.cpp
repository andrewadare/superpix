// See http://julia.readthedocs.org/en/latest/manual/embedding/
// for more info and examples on embedding Julia in C/C++.

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include <julia.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

bool checkImageSizes(cv::Mat &img1, cv::Mat &img2);
void cvMat2Array(cv::Mat &bgr_img, cv::Mat &dep_img, unsigned int* pbuff);
int loadJuliaFunctions(string filename);

int main(int argc, char *argv[])
{
  jl_init(JULIA_INIT_DIR);

  // 1st arg (required): Name of input image
  if (argc < 3)
  {
    cout << "Usage: " << argv[0] << " colorimg.jpg depthimg.jpg" << endl;
    return 1;
  }
  string rgbfile = string(argv[1]);
  string depfile = string(argv[2]);

  // Read image and get dimensions
  cv::Mat rgb_img = cv::imread(rgbfile);
  cv::Mat dep_img = cv::imread(depfile);

  if (!checkImageSizes(rgb_img, dep_img))
    return -1;

  int nrows = rgb_img.rows, ncols = rgb_img.cols;
  const int sz = nrows*ncols;

  cout << "Read color cv::Mat from " << rgbfile
       << " (" << nrows << "x" << ncols << " rows x cols)" << endl;
  cout << "Read depth cv::Mat from " << depfile
       << " (" << nrows << "x" << ncols << " rows x cols)" << endl;

  // Copy image data into an array
  unsigned int* pbuff = new unsigned[sz];
  cvMat2Array(rgb_img, dep_img, pbuff);

  // TEMP: write pbuff to a text file
  ofstream fs("pbuff.txt");
  for (int i = 0; i < sz; i++)
    fs << pbuff[i] << endl;
  fs.close();

  // Wrap the array and its 2d dimensions for julia.
  jl_value_t *array_type = jl_apply_array_type(jl_uint32_type, 1);
  jl_array_t *pbuff_jl = jl_ptr_to_array_1d(array_type, pbuff, sz, false);

  jl_value_t *nrows_jl = jl_box_int64(nrows);
  jl_value_t *ncols_jl = jl_box_int64(ncols);

  // Tell Julia not to trash our stuff.
  jl_value_t **args;
  JL_GC_PUSHARGS(args, 3);
  args[0] = (jl_value_t*)pbuff_jl;
  args[1] = nrows_jl;
  args[2] = ncols_jl;

  // Read Julia function definitions into memory
  char jlfile[256];
  sprintf(jlfile, "%s/repos/superpix/cpp-embed/src/segmenter.jl", getenv("HOME"));
  int ret = loadJuliaFunctions(string(jlfile));
  if (ret != 0)
    return -1;

  jl_function_t *f = jl_get_function(jl_current_module, "main");
  jl_call(f, args, 3);

  JL_GC_POP();
  jl_atexit_hook();
  return 0;
}

bool checkImageSizes(cv::Mat &img1, cv::Mat &img2)
{
  if (img1.rows != img2.rows || img1.cols != img2.cols)
  {
    cout << "Error: Image sizes do not match:" << endl;
    cout << "  " << img1.rows << "x" << img1.cols << endl;
    cout << "  " << img2.rows << "x" << img2.cols << endl;
    return false;
  }
  return true;
}

void cvMat2Array(cv::Mat &bgr_img, cv::Mat &dep_img, unsigned int* pbuff)
{
  // Use 32 bit unsigned int to hold a pixel in DRGB format as follows:
  // from left to right,
  // the first 8 bits are for the depth channel
  // the next 8 bits are for the red channel
  // the next 8 bits are for the green channel
  // the last 8 bits are for the blue channel
  //
  // To use this function, allocate a C++ buffer like so:
  // unsigned int* pbuff = new unsigned[img.rows*img.cols];

  if (!checkImageSizes(bgr_img, dep_img))
    return;

  unsigned char *cdat = (unsigned char*)(bgr_img.data);
  unsigned char *ddat = (unsigned char*)(dep_img.data);
  int r=0, g=0, b=0, d=0, counter=0;

  int ncc = bgr_img.channels();
  int ndc = dep_img.channels();
  for (int i = 0; i < bgr_img.rows; i++)
  {
    for (int j = 0; j < bgr_img.cols; j++)
    {
      b = cdat[bgr_img.step*i + ncc*j    ];
      g = cdat[bgr_img.step*i + ncc*j + 1];
      r = cdat[bgr_img.step*i + ncc*j + 2];

      d = ddat[dep_img.step*i + ndc*j];

      pbuff[counter++] = (d << 24) | (r << 16) | (g << 8) | b;
    }
  }
}

// This function is intended to read in a julia source file containing
// definitions for functions, types, etc.
// Any function calls in the file will be triggered!
int loadJuliaFunctions(string filename)
{
  ifstream file(filename);
  if (file)
  {
    string str;
    string file_contents;
    while (getline(file, str))
    {
      file_contents += str;
      file_contents.push_back('\n');
    }

    jl_eval_string(file_contents.c_str());
  }
  else
  {
    printf("Problem finding or opening %s.\n", filename.c_str());
    return -1;
  }
  return 0;
}
