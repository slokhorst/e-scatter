# e-scatter

Fast CUDA-enabled simulator for electron scattering processes in materials.

## Building

Dependencies of the core simulator:
* C++11 compiler (tested with GCC 4.8 - 6)
* CUDA 7.5 or newer
* CUB

Build the program:
````sh
cd e-scatter
mkdir build
cd build
cmake ..
make
````

To also build all of the utilities, append `-DWITH_UTILS=1` to the `cmake` command.
Dependencies required to build the utilities:
* Boost
* muParser
* ImageMagick
* FFTW
* SDL2
* OpenGL
* GLEW

## Usage

The main program is located at `bin/csdem`. It requires a geometry file, which can be generated with `bin/line-gen`, and an exposure file, which can be generated with `scripts/pri-gen.py`.
