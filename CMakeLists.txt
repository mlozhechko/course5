cmake_minimum_required(VERSION 3.9)
project(course)

set(CMAKE_CXX_STANDARD 17)
if(UNIX AND NOT APPLE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -fopenmp")
endif()

if(APPLE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -Xpreprocessor -fopenmp")
        set(CMAKE_EXE_LINKER_FLAGS "-I /usr/local/include -L/usr/local/lib -lomp")
endif()

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(Boost 1.61.0 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include_directories("./project/src" "./project/include")

add_executable(
        course
        project/src/main.cpp
        project/src/tetra.cpp
        project/src/line.cpp
        project/src/plane.cpp
        project/src/object3d_base.cpp
        project/src/object3d_accretion_disk.cpp
        project/src/object3d_roche_lobe.cpp
        project/src/object3d_sphere.cpp
        project/src/object2d.cpp)

target_link_libraries(course ${VTK_LIBRARIES})
target_link_libraries(course ${Boost_LIBRARIES})