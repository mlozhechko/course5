cmake_minimum_required(VERSION 3.11)
project(course)

set(CMAKE_CPP_COMPILER gcc)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -fopenmp -ftree-vectorize -mavx -march=native -mfpmath=sse")

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(Boost 1.61.0 COMPONENTS container REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})


include_directories("./project/src" "./project/include")

add_executable(
        course
        project/src/main.cpp
        project/src/tetra.cpp
        project/src/line.cpp
        project/src/plane.cpp
        project/src/app.cpp)

target_link_libraries(course ${VTK_LIBRARIES})
target_link_libraries(course ${Boost_LIBRARIES})