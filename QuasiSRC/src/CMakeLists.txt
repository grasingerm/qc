cmake_minimum_required(VERSION 3.1)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

include_directories(../triangulation_server/include)
include_directories(../include)
include_directories(/usr/local/include)
include_directories(/usr/include)

file(GLOB libsources *.f *.c *.cc)

add_library(libquasi SHARED ${libsources})

if (UNIX)
    target_link_libraries(libquasi m)
endif (UNIX)

target_link_libraries(libquasi Threads::Threads)
