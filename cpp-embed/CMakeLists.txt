cmake_minimum_required(VERSION 2.8)

project(segmenter)

# Julia
#######

execute_process(COMMAND which julia OUTPUT_VARIABLE JL_SHARE)
string(REPLACE "bin" "share" JL_SHARE ${JL_SHARE})
string(REPLACE "\n"  ""      JL_SHARE ${JL_SHARE})
string(REPLACE "share" "lib" JL_LIB_DIR ${JL_SHARE})
execute_process(COMMAND ${JL_SHARE}/julia-config.jl --cflags OUTPUT_VARIABLE CMAKE_CXX_FLAGS)
execute_process(COMMAND ${JL_SHARE}/julia-config.jl --ldlibs OUTPUT_VARIABLE CMAKE_SHARED_LINKER_FLAGS)
string(REPLACE "\n"  "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
string(REPLACE "\n"  "" CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall")

link_directories(${JL_LIB_DIR})

message("JL_SHARE:   " ${JL_SHARE})
message("JL_LIB_DIR: " ${JL_LIB_DIR})
message("CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})
message("CMAKE_SHARED_LINKER_FLAGS:   " ${CMAKE_SHARED_LINKER_FLAGS})

# OpenCV
########

# For system-wide install location (default is /usr/local)
if(EXISTS "${CMAKE_INSTALL_PREFIX}/share/OpenCV/OpenCVConfig.cmake")
  set(OpenCV_DIR ${CMAKE_INSTALL_PREFIX}/share/OpenCV)

# For local/custom OpenCV install location (modify/add as needed)
elseif(EXISTS "$ENV{HOME}/opencv/opencv-3.0.0-rc1/build/OpenCVConfig.cmake")
  set(OpenCV_DIR "$ENV{HOME}/opencv/opencv-3.0.0-rc1/build")
endif()

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

message("\nOpenCV_DIR: " ${OpenCV_DIR})

# SLIC
######

# set(SLIC ${PROJECT_SOURCE_DIR}/slic)
# include_directories(${SLIC})
# add_library(slic SHARED 
#             ${SLIC}/SLIC.cpp)

# This project
##############

set(SRC ${PROJECT_SOURCE_DIR}/src)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

# Get list of .cpp files. CMake must be re-run whenever new files are added.
file(GLOB MY_CPPS RELATIVE ${SRC} ${SRC}/*.cpp)

# Make an executable from each .cpp file
message("\nExecutables to be built:")
foreach(CPPFILE ${MY_CPPS})
  # Extract foo from foo.cpp as exename
  string(REPLACE ".cpp" "" exename ${CPPFILE})
  # Generate foo executable
  add_executable(${exename} ${SRC}/${CPPFILE})
  target_link_libraries(${exename} 
                        ${OpenCV_LIBS} 
                        julia)
  # Print executable list
  message("  ${exename}")
endforeach()

message("\nExecutables will be installed in ${EXECUTABLE_OUTPUT_PATH}.\n")
