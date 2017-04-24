# e-scatter

Fast CUDA-enabled simulator for electron scattering processes in materials.

[![Build Status](https://travis-ci.org/eScatter/e-scatter.svg?branch=master)](https://travis-ci.org/eScatter/e-scatter)

## Building

Dependencies of the core simulator:
* C++14 compiler (tested with GCC 5)
* CUDA 8.0 or newer
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

The main program is located at `bin/csdem`. It requires a geometry file, which
can be generated with `bin/line-gen`, and an exposure file, which can be
generated with `scripts/pri-gen.py`.

## Testing

To run the unit tests, the googletest library should be available in the 
3rdparty directory; clone it from github and run `cmake` with `-DWITH_TEST=ON`.

````sh
cd 3rdparty
git clone git@github.com:google/googletest.git
mkdir -p ../build
cd ../build
cmake .. -DWITH_TEST=ON
make
make test
````

## Creating materials

The material files contain most of the physics involved.

* [ELSEPA](http://adsabs.harvard.edu/abs/2005CoPhC.165..157S) can be downloaded from the
[Computer Physics communications Program library](http://www.cpc.cs.qub.ac.uk/). It has an
attribute-only license for non-commercial use. We use ELSEPA to compute Mott cross-sections.

* Livermore database [EPDL97](https://www-nds.iaea.org/epdl97/) or 
[ENDF/B-VII.1](http://www.nndc.bnl.gov/endf/b7.1/download.html)


