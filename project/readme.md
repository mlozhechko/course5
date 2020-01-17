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
./course
```