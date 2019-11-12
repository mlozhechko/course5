cmake_minimum_required(VERSION 3.11)
project(course)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g")

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(course main.cpp)
target_link_libraries(course ${VTK_LIBRARIES})