cmake_minimum_required(VERSION 3.01)
project(int2048)

set(CMAKE_CXX_STANDARD 14)

include_directories(${project_source_dir}/)

add_executable(code int2048.cpp)
