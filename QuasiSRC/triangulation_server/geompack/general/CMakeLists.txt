cmake_minimum_required(VERSION 3.1)

enable_language(C Fortran)
file(GLOB srcs *.f)
add_library(libgeneral SHARED ${srcs})
