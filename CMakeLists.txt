cmake_minimum_required(VERSION 3.11)
project(course)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O2")

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories("./project/src" "./project/include")

add_executable(
        course
        project/src/main.cpp
        project/src/lite_tetrahedron.cpp
        project/src/line.cpp
        project/src/plane.cpp
        project/src/application_core.cpp)

target_link_libraries(course ${VTK_LIBRARIES})