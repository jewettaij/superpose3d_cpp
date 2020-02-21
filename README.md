[![Build Status](https://travis-ci.org/jewettaij/superpose3d_cpp.svg?branch=master)](https://travis-ci.org/jewettaij/superpose3d_cpp.svg?branch=master)
[![codecov](https://codecov.io/gh/jewettaij/superpose3d_cpp/branch/master/graph/badge.svg)](https://codecov.io/gh/jewettaij/superpose3d_cpp)
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()
[![C++11](https://img.shields.io/badge/C%2B%2B-11-blue.svg)](https://isocpp.org/std/the-standard)
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/jewettaij/superpose3d_cpp)]()


superpose3d_cpp
===========

Note: There is a ***python version*** of this repository
[***here***](https://github.com/jewettaij/superpose3d).

**superpose3d_cpp** contains a [header file](include/superpose3d.hpp)
containing the definition of a class, *Superpose3D*.  It's single public member
function, *Superpose()*, takes two N×3 arrays representing coordinates of points
from a point cloud (denoted *X<sub>ni</sub>* and *x<sub>ni</sub>*) as arguments,
and attempts to superimpose them (treating them as rigid bodies).
More specifically, *Superpose3D::Superpose()* attempts to superimpose
them using **rotations**, **translations**, and (optionally) **scale**
transformations in order to minimize the root-mean-squared-distance (RMSD)
between corresponding points from either point cloud, where RMSD is defined as:

<img src="http://latex.codecogs.com/gif.latex?\large&space;RMSD=\sqrt{\,\frac{1}{N}\,\sum_{n=1}^N\,\,\sum_{i=1}^3 \left|X_{ni}-\left(\sum_{j=1}^3 cR_{ij}x_{nj}+T_i\right)\right|^2}"/>

...where:
```
   T = a translation vector (a 1-D array containing x,y,z displacements),
   R = a rotation matrix    (a 3x3 array whose determinant = 1),
   c = a scalar             (a number, optional, 1 by default)
```
After invoking Superpose3D::Superpose(), the optimal translation, rotation and
scale factor are stored in Superpose3D data members named
*T*, *R*, and *c*, respectively.
(*T* is implemented as a C-style array, and
 *R* is implemented as a C-style 3x3 array in pointer-to-pointer format.)

A weighted version of the RMSD minimization algorithm is also available
if the caller supplies an extra argument specifying the weight of every
point in the cloud (*w<sub>n</sub>*).  In that case, RMSD is defined as:
<img src="http://latex.codecogs.com/gif.latex?\large&space;RMSD=\sqrt\left\sum_{n=1}^N\,w_n\,\sum_{i=1}^3 \left|X_{ni}-\left(\sum_{j=1}^3 c R_{ij}x_{nj}+T_i\right)\right|^2\quad\middle/\quad\sum_{n=1}^N w_n}\right}"/>

The coordinate arrays (*X<sub>ni</sub>* and *x<sub>ni</sub>*)
can be implemented as T\*\* (pointer-to-pointer),
vector\<vector\<T\>\>&, fixed-size arrays,
or any other C or C++ object which supports \[\]\[\].
(Here **T** is any real numeric type.  Complex numbers are not supported.)
Similarly, the weights (*w*, if specified) can be implemented as arrays
or any other C++ container supporting \[\].

### Algorithm

This function implements a more general variant of the method from this paper:
R. Diamond, (1988)
"A Note on the Rotational Superposition Problem",
 Acta Cryst. A44, pp. 211-216.

### Scale transformations

This version has been augmented slightly to support scale transformations.
(I.E. multiplication by scalars.  This can be useful for the registration
of two different annotated volumetric 3-D images of the same object taken
at different magnifications.)

Note that if you enable scale transformations, you should be wary if the function returns a negative **c** value.  Negative **c** values correspond to inversions (reflections).  For this reason, if you are using this function to compare the conformations of molecules, you should probably set the fourth argument to *false*.  This will prevent matching a molecule with its stereoisomer.

### Rotation angle, axis, and quaternion
If the corresponding rotation angle and rotation axis are also needed, they
can be inferred from the ***q*** data member ("*superposer.q*" in the
example below). After invoking Superpose(), the *q* member will store the
[quaternion corresponding to rotation *R*](https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation).
The first element of *q* will store *cos(θ/2)* (where *θ* is the
rotation angle).  The remaining 3 elements of *q* will store the
axis of rotation (with length *sin(θ/2)*).



##  Example usage

```cpp
#include "superpose3d.hpp"
using namespace superpose3d;

// ...
int N = 147913; // the number of points in each cloud
double **X;     // <-- Nx3 array of coordinates for the "frozen" point cloud
double **x;     // <-- Nx3 array of coordinates for the "mobile" point cloud

// Allocate space for X and x, and load their coordinates (omitted) ...

// Create an instance of the "Superpose3D" class.

Superpose3D<double, double **> superposer(N);

// "double **" is the type of array for storing coordinates in this example.
// (If the arrays are read-only, then you can use "double const* const*".)
// You can also use vectors or other objects which support [][]. For example:
// Superpose3D<double, vector<vector<double>>&> superposer(N);
// This will allocate memory to store temporary arrays used later.
// Once created, it can be used multiple times on different point clouds of the
// same size (without incurring the cost of memory allocation on the heap).

// Calculate the optimal supperposition between the two point clouds (X and x)

double rmsd = superposer.Superpose(X, x);

// Note: The optimal rotation, translation, and scale factor will be stored in
//       superposer.R, superposer.T, and superposer.c, respectively.
//       (A quaternion describing the rotation is stored in superpose.q.)
```
*(A complete working example can be found [here](tests/test_superpose3d.cpp).)*

By default scale transformations are disabled.  (By default *c=1*.)
If you want to allow scale transformations, then use:
```
superposer.Superpose(X, x, true);
```

### Example using non-uniform weights

By default point in the point cloud will be given equal weights when
calculating RMSD.  If you want to specify different weights for each point
(ie. *w<sub>n</sub>* in the formula above), then see the following example:

```cpp
// ...

double *w;    // optional: weights used in calculation of RMSD
// Fill the "w" array (omitted)

Superpose3D<double, double **, double*> superposer(N, w);
double rmsd = superposer.Superpose(X, x);
// "double*" is the type of array for the weights in this example ("w").
// (For read-only arrays, you can use use "double const*".)
```


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

## Benchmarking

The performance of the algorithm is *O(N)*.
For large *N*, the computation time required (per point in the cloud)
is approximately 4.0e-08 seconds.

<sub>
(Details: This was measured on a single 1.7GHz i5-4210U CPU core.
For this test, the [tests/test_superpose3d.cpp](tests/test_superpose3d.cpp)
file was compiled using g++ with the "-Ofast" compiler flag, and then run with
and without with the line invoking Superpose3D::Superpose() commented out.)
</sub>

## Development Status: *Stable*

The code in [superpose3d.hpp](include/superpose3d.hpp) has been
[tested](.travis.yml) for accuracy, coverage, memory safety and speed
using thousands of large randomly generated point clouds.

## License

*superpose3d_cpp* is available under the terms of the [MIT license](LICENSE.md).

