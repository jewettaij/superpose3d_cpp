[![Build Status](https://travis-ci.org/jewettaij/superpose3d_cpp.svg?branch=master)](https://travis-ci.org/jewettaij/superpose3d_cpp.svg?branch=master)
[![codecov](https://codecov.io/gh/jewettaij/superpose3d_cpp/branch/master/graph/badge.svg)](https://codecov.io/gh/jewettaij/superpose3d_cpp)
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()
[![C++11](https://img.shields.io/badge/C%2B%2B-11-blue.svg)](https://isocpp.org/std/the-standard)
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/jewettaij/superpose3d_cpp)]()


superpose3d_cpp
===========

Note: There is a ***python version*** of this repository
[***here***](https://github.com/jewettaij/superpose3d).

**superpose3d_cpp** is a header-only C++ library containing the definition
of a class whose single public member function, *Superpose()*,
takes two N×3 arrays representing coordinates of points
from a point cloud (denoted *X<sub>ni</sub>* and *x<sub>ni</sub>*) as arguments.

*Superpose3D::Superpose()* attempts to superimpose
them using **rotations**, **translations**, and (optionally) **scale**
transformations in order to minimize the root-mean-squared-distance (RMSD)
between corresponding points from either point cloud, where RMSD is defined as:

<img src="http://latex.codecogs.com/gif.latex?\large&space;RMSD=\sqrt\left\sum_{n=1}^N\,w_n\,\sum_{i=1}^3 \left|X_{ni}-\left(\sum_{j=1}^3 c R_{ij}x_{nj}+T_i\right)\right|^2\quad\middle/\quad\sum_{n=1}^N w_n}\right}"/>

If *w<sub>n</sub>* are omitted, then equal weights are used.  In that case:

<img src="http://latex.codecogs.com/gif.latex?\large&space;RMSD=\sqrt{\,\frac{1}{n}\,\sum_{n=1}^N\,\,\sum_{i=1}^3 \left|X_{ni}-\left(\sum_{j=1}^3 cR_{ij}x_{nj}+T_i\right)\right|^2}"/>

...where:
```
   T = a translation vector (a 1-D array containing x,y,z displacements),
   R = a rotation matrix    (a 3x3 array whose determinant = 1),
   c = a scale factor       (a number, optional, 1 by default)
```

After invoking Superpose3D::Superpose(), the optimal translation, rotation and
scale factor are stored in data members named *T*, *R*, and *c*, respectively.
(*T* is implementad as a C-style array and
 *R* is implemented as a pointer-to-pointer.)

The coordinate arrays (*X<sub>ni</sub>* and *x<sub>ni</sub>*)
can be implemented as T\*\* (pointer-to-pointer),
vector\<vector\<T\>\>&, fixed-size arrays,
or any other C or C++ object which supports \[\]\[\].
(Here **T** is any real numeric type.  Complex numbers are not supported.)
Likewise, the weights (*w*, if specified) can be implemented as arrays
or any other C++ container supporting \[\].



##  Example usage

```cpp
#include "superpose3d.hpp"
using namespace superpose3d;

// ...

double **X;   // 1st point cloud (note: use "double **X" not "double (*X)[3]")
double **x;   // 2nd point cloud (the mobile point cloud)
double *w;    // optional weights used in calculation of RMSD

// Allocate space for X and x, and load their coordinates (omitted)
// ...

// Create an instance of the "Superpose3D" class.

Superpose3D<double, double const* const*, double const*> superposer(N, w);

// This will allocate memory to store temporary arrays used later.
// Once created, it can be used multiple times on different point clouds of the
// same size (without incurring the cost of memory allocation on the heap).
// Notes:
// -"N" is the number of points in either point cloud.
// -"double const* const*" is the type of array for storing coordinates in this
//   example.  However you can use vectors or fixed sized arrays. For example:
//  Superpose3D<double, vector<vector<double>>, vector<double>> superposer(N,w);
// -"double const*" is the type of array for storing the weights in this
//   example.  If you are content to assign equal weight to each point, you
//   can omit the 3rd template argument and also omit "w" from the constructor.
//   (Most users will use it this way.)  In that case instantiate it this way:
//  Superpose3D<double, double const* const*> superposer(N);

// Calculate the optimal supperposition between the two point clouds (X and x)

double rmsd =
  superposer.Superpose(X, x);

// Note: The optimal rotation, translation, and scale factor will be stored in
//       superposer.R, superposer.T, and superposer.c, respectively.
```
*(A complete working example can be found [here](tests/test_superpose3d.cpp).)*

Each point in the point cloud will be given equal weights when calculating RMSD.
If you want to specify the weights (*w<sub>n</sub>* in the formula above),
then use:
```
superposer.Superpose(X, x, w);
```
This function implements a more general variant of the method from this paper:
R. Diamond, (1988)
"A Note on the Rotational Superposition Problem",
 Acta Cryst. A44, pp. 211-216.

This version has been augmented slightly to support scale transformations.
(I.E. multiplication by scalars.  This can be useful for the registration
of two different annotated volumetric 3-D images of the same object taken
at different magnifications.)

By default scale transformations are disabled.  (By default *c=1*.)
If you want to allow scale transformations, then use:
```
superposer.Superpose(X, x, w, true);
```

Note that if you enable scale transformations, you should be wary if the function returns a negative **c** value.  Negative **c** values correspond to inversions (reflections).  For this reason, if you are using this function to compare the conformations of molecules, you should probably set the fourth argument to *false*.  This will prevent matching a molecule with its stereoisomer.


## Downloading

This repository has a [dependency](https://github.com/jewettaij/jacobi_pd)
so you must use the **--recursive** argument when cloning it.  For example:

```
git clone --recursive https://github.com/jewettaij/superpose3d_cpp ~/superpose3d_cpp
```
## Installation

This is a header-only library.

Copy the *hpp* file(s) in the [include](include) subdirectory,
and the *hpp* files in the
[jacobi_pd/include](https://github.com/jewettaij/jacobi_pd/tree/master/include)
subdirectory to a location in your
[include path](https://www.rapidtables.com/code/linux/gcc/gcc-i.html).
*(Both repositories share the same license.)*

## Development Status: *Beta*

## License

*superpose3d_cpp* is available under the terms of the [MIT license](LICENSE.md).

