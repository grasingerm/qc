cmake_minimum_required(VERSION 3.1)

enable_language(C Fortran)
file(GLOB srcs *.f)
add_library(libbasic3d SHARED ${srcs})
