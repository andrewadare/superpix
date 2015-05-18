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

void cvMat2Array(cv::Mat &img, unsigned int* pbuff);
int loadJuliaFunctions(string filename);

int main(int argc, char *argv[])
{
  jl_init(JULIA_INIT_DIR);

  // 1st arg (required): Name of input image
  if (argc < 2)
  {
    cout << "Usage: " << argv[0] << " in.jpg" << endl;
    return 1;
  }
  string imfilename = string(argv[1]);

  // Read image and get dimensions
  cv::Mat img = cv::imread(imfilename);
  cout << "Read cv::Mat from " << imfilename
       << " (" << img.rows << "x" << img.cols << " rows x cols)" << endl;

  // Copy image data into an array  
  const int sz = img.rows*img.cols;
  unsigned int* pbuff = new unsigned[sz];
  cvMat2Array(img, pbuff);

  // TEMP: write pbuff to a text file
  ofstream fs("pbuff.txt");
  for (int i = 0; i < sz; i++)
    fs << pbuff[i] << endl;
  fs.close();

  // Wrap the array and its 2d dimensions for julia.
  jl_value_t *array_type = jl_apply_array_type(jl_uint32_type, 1);
  jl_array_t *pbuff_jl = jl_ptr_to_array_1d(array_type, pbuff, sz, false);

  jl_value_t *nrows_jl = jl_box_int64(img.rows);
  jl_value_t *ncols_jl = jl_box_int64(img.cols);

  // Tell Julia not to trash our stuff.
  jl_value_t **args;
  JL_GC_PUSHARGS(args, 3);
  args[0] = (jl_value_t*)pbuff_jl;
  args[1] = nrows_jl;
  args[2] = ncols_jl;

  // Read Julia function definitions into memory
  int ret = loadJuliaFunctions("/Users/adare/cv/mat2jl/src/segmenter.jl");
  if (ret != 0)
    return -1;

  jl_function_t *f = jl_get_function(jl_current_module, "main");
  jl_call(f, args, 3);

  JL_GC_POP();
  jl_atexit_hook();
  return 0;
}

void cvMat2Array(cv::Mat &img, unsigned int* pbuff)
{
  // Use 32 bit unsigned int to hold a pixel in ARGB format as follows:
  // from left to right,
  // the first 8 bits are for the alpha channel (currently unused)
  // the next 8 bits are for the red channel
  // the next 8 bits are for the green channel
  // the last 8 bits are for the blue channel
  //
  // To use this function, allocate a C++ buffer like so:
  // unsigned int* pbuff = new unsigned[img.rows*img.cols];

  unsigned char *dat = (unsigned char*)(img.data);
  int r=0, g=0, b=0, counter=0; // add int a=0 ?
  int c = img.channels();
  for (int i = 0; i < img.rows; i++)
  {
    for (int j = 0; j < img.cols; j++)
    {
      b = dat[img.step*i + c*j    ];
      g = dat[img.step*i + c*j + 1];
      r = dat[img.step*i + c*j + 2];
      
      // a = dat[img.step*i + c*j + 3];
      // Or, assign a from a second depth image?

      pbuff[counter++] = (r << 16) | (g << 8) | b;
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
