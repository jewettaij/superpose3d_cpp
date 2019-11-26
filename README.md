[![Build Status](https://travis-ci.org/jewettaij/superpose3d_cpp.svg?branch=master)](https://travis-ci.org/jewettaij/superpose3d_cpp.svg?branch=master)
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()
[![GitHub All Releases](https://img.shields.io/github/downloads/jewettaij/superpose3d_cpp/total)]()
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/jewettaij/superpose3d_cpp)]()


superpose3d_cpp
===========

## WARNING: THIS CODE IS NOT WORKING YET (-andrew 2019-11-22)

**superpose3d_cpp** is a header-only C++ library containing a defninition
of a class whose single public member function, *Superpose()*,
takes two N x 3 multidimensional arrays
(*of the same length*, **N**) representing points
from a point cloud (**X_i** and **x_i**) as arguments.
Treating them as rigid objects,
*Superpose3D::Superpose()* attempts to superimpose
them using **rotations**, **translations**, and (optionally) **scale**
transformations in order to minimize the root-mean-squared-distance (RMSD)
between corresponding points from either point cloud, where RMSD is defined as:
```
   RMSD = sqrt((Σ_i  w_i * |X_i - Σ_j(c*R_ij*x_j + T_i))|^2) / (Σ_j w_j))
```
If *w_i* are omitted (ie. if *w_i = nullptr*),
then equal weights are used.  In that case:
```
   RMSD = sqrt(( Σ_i |X_i - Σ_j(c*R_ij*x_j + T_i) )|^2 ) / N)
```
...where:
```
   T = a translation vector (a 1-D numpy array containing x,y,z displacements),
   R = a rotation matrix    (a 3x3 numpy array whose determinant = 1),
   c = a scale factor       (a number)
n```
After invoking Superpose3D::Superpose(), the optimal translation, rotation and
scale factor are stored in data members named *T*, *R*, and *c*, respectively.

##  Example usage

```
#include "superpose3d.hpp"
using namespace superpose3d_lammps;

int main(int argc, char **argv) {

  double **X;   // (note: use "double **X" not "double (*X)[3]")
  double **x;
  double **w;

  // Allocate space for X and x, and load their coordinates (omitted)

  Superpose3D superposer(N);

  // Calculate the optimal supperposition between the two point clouds (X and x)

  double rmsd =
    superposer.Superpose(X, x);

  // Note: The optimal rotation, translation, and scale factor will be stored in
  //       superposer.R, superposer.T, and superposer.c, respectively.
}
```
If you want to specify the weights (*w_i* in the formula above), then use:
```
superposer.Superpose(X, x, w);
```
This function implements a more general variant of the method from this paper:
R. Diamond, (1988)
"A Note on the Rotational Superposition Problem",
 Acta Cryst. A44, pp. 211-216.

This version has been augmented slightly to support scale transformations.  (I.E. multiplication by scalars.  This can be useful for the registration of two different annotated volumetric 3-D images of the same object taken at different magnifications.)

By default scale transformations are disabled.  (By default $c=1$.)
If you want to allow scale transformations, then use:
```
superposer.Superpose(X, x, w, true);
```

Note that if you enable scale transformations (i.e. if the fourth argument is *true*), you should be wary if the function returns a negative **c** value.  Negative **c** values correspond to inversions (reflections).  For this reason, if you are using this function to compare the conformations of molecules, you should probably set the fourth argument to *false*.  This will prevent matching a molecule with its stereoisomer.


## Downloading

This repository has a [dependency](https://github.com/mrcdr/lambda-lanczos)
so you must use the **--recursive** argument when cloaning it.  For example:

```
git clone --recursive https://github.com/jewettaij/superpose3d_cpp ~/superpose3d_cpp
```

## Installation

This is a header-only library.

Copy the files *"include/superpose3d.hpp"*, and all of the *hpp* files in the
*"lambda-lanczos/include/lambda_lanczos/" directory to a location in your
[include path](https://www.rapidtables.com/code/linux/gcc/gcc-i.html).


## License

superpose3d is available under the terms of the [MIT license](LICENSE.md).
