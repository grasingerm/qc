cmake_minimum_required(VERSION 3.1)

include_directories(../../include)
include_directories(../../triangulation_server/include)
include_directories(../../triangulation_server/src)

add_executable(quasi main.cc stamp.cc)
add_executable(test-quasi test.cc stamp.cc)

target_link_libraries(quasi libquasi libtriang CGAL CGAL_Core gmp boost_thread mpfr)
target_link_libraries(test-quasi libquasi libtriang CGAL CGAL_Core gmp boost_thread mpfr)
