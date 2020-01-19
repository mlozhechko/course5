C++17 parallel ray-tracing application. 

Requirements:
1. Locally installed VTK library (version >= 7.0)
2. Locally installed boost::program_options library
3. CMake (version >= 3.11)
4. C++17 compatible compiler
5. OpenMP compatible compiler

Source data format is VTK unstructured grid of tetrahedrons

Build:
``` 
cmake build .
make -j8
```
 
Run:
```
./course --help
Allowed options:
  -h [ --help ]                         produce help message
  -f [ --file ] arg                     source file
  -d [ --destination ] arg              destination file
  -j [ --threads ] arg                  number of parallel threads
  -x [ --resolution_x ] arg (=1200)     set x axis resolution
  -y [ --resolution_y ] arg (=900)      set y axis resolution
  -X [ --angle_around_x ] arg (=0)      rotate view plane by angle around x 
                                        axis
  -Y [ --angle_around_y ] arg (=0)      rotate view plane by angle around y 
                                        axis
  -D [ --donor_angle ] arg (=0)         initial donor angle around y axis
  -I [ --initial system angle ] arg (=0)
                                        initial angle of system y axis

```
For example
```
./course -f filename.vtk -d result.vti -j16 -x 2400 -y 1800
```

Known problems:
1. For know application can not limit RAM usage. The most part of RAM usage fall on lines objects (pixels in future).
Therefore if you have std::bad_alloc or similar errors just lower resolution of plane (f.e. 600 x 450)
2. 3d object hit outside of domain causes segmentation fault (should be fixed as soon as possible)
3. For some reasons there are MAX_NUMBER_OF_THREADS = 32 constant in config.hpp;
If you want to use more than 32 threads just change constant, it won't cause any problem. 