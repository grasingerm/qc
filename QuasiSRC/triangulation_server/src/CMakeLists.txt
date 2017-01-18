cmake_minimum_required(VERSION 3.1)

include_directories(../include)
add_library(libtriang SHARED triangulation_client.c util.c)
add_executable(triang_server triangulation_server.c error.f)
target_link_libraries(triang_server libdeltr3d libbasic3d libgeneral libtriang)
