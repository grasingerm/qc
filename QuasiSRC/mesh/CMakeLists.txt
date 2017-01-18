cmake_minimum_required(VERSION 3.1)
project(qcmesh)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

include_directories(../include)
add_executable(mesh mesh.cc)

if (UNIX)
    target_link_libraries(mesh m)
endif (UNIX)

target_link_libraries(mesh libquasi)
