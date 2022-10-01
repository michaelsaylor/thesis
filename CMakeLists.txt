cmake_minimum_required(VERSION 3.23)
project (thesis)

set(CMAKE_CXX_STANDARD 20)

add_executable(thesis main.cpp crosssection_1_v3-W.cpp dNdW.cpp sigmavst_v2.cpp sigmavsw_v2.cpp sigmavsx_v2.cpp tbfk1.hpp TGreenFun.cpp TGreenFun.hpp)
