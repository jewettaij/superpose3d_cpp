[![Build Status](https://travis-ci.org/jewettaij/superpose3d_cpp.svg?branch=master)](https://travis-ci.org/jewettaij/superpose3d_cpp.svg?branch=master)
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/jewettaij/superpose3d_cpp)]()


superpose3d_cpp
===========

Note: There is a python version of this repository
[here](https://github.com/jewettaij/superpose3d).

**superpose3d_cpp** is a header-only C++ library containing the definition
of a class whose single public member function, *Superpose()*,
takes two N×3 arrays representing coordinates of points
from a point cloud (denoted *X<sub>ni</sub>* and *x<sub>ni</sub>*) as arguments.
(Both *X<sub>ni</sub* and *x<sub>ni</sub>* should be implemented as C style
pointer-to-pointer arrays.)
Treating them as rigid objects,
*Superpose3D::Superpose()* attempts to superimpose
them using **rotations**, **translations**, and (optionally) **scale**
transformations in order to minimize the root-mean-squared-distance (RMSD)
between corresponding points from either point cloud, where RMSD is defined as:

<img src="http://latex.codecogs.com/gif.latex?\large&space;RMSD=\sqrt\left\sum_{n=1}^N\,w_n\,\sum_{i=1}^3 \left|X_{ni}-\left(\sum_{j=1}^3 c R_{ij}x_{nj}+T_i\right)\right|^2\quad\middle/\quad\sum_{n=1}^N w_n}\right}"/>

If *w<sub>n</sub>* are omitted (ie. if *w<sub>n</sub> = nullptr*),
then equal weights are used.  In that case:

<img src="http://latex.codecogs.com/gif.latex?\large&space;RMSD=\sqrt{\,\frac{1}{n}\,\sum_{n=1}^N\,\,\sum_{i=1}^3 \left|X_{ni}-\left(\sum_{j=1}^3 cR_{ij}x_{nj}+T_i\right)\right|^2}"/>

...where:
```
   T = a translation vector (a 1-D numpy array containing x,y,z displacements),
   R = a rotation matrix    (a 3x3 numpy array whose determinant = 1),
   c = a scale factor       (a number)
```

After invoking Superpose3D::Superpose(), the optimal translation, rotation and
scale factor are stored in data members named *T*, *R*, and *c*, respectively.


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

Superpose3D superposer(N);
// (N is the number of points in either point cloud.)

// Calculate the optimal supperposition between the two point clouds (X and x)

double rmsd =
  superposer.Superpose(X, x);

// Note: The optimal rotation, translation, and scale factor will be stored in
//       superposer.R, superposer.T, and superposer.c, respectively.
```
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

This repository has a [dependency](https://github.com/mrcdr/lambda-lanczos)
so you must use the **--recursive** argument when cloaning it.  For example:

```
git clone --recursive https://github.com/jewettaij/superpose3d_cpp ~/superpose3d_cpp
```
## Installation

This is a header-only library.

Copy the files [include/superpose3d.hpp](include),
and all of the *hpp* files in the
[lambda-lanczos/include/lambda_lanczos](lambda-lanczos/include/lambda_lanczos)
directory to a location in your
[include path](https://www.rapidtables.com/code/linux/gcc/gcc-i.html).

## *Additional Modifications Needed*

*Currently, if you copy these files into the same directory,
you will need to modify your include path
(ie. your "-I" compiler arguments)
and the *#include* statements in your header files
([here](./include/superpose3d.hpp) and
 [here](https://github.com/mrcdr/lambda-lanczos/blob/master/include/lambda_lanczos/lambda_lanczos.hpp))
to delete "lambda_lanczos/" from these paths where it appears.
I will try to persuade the author of "lambda_lanczos" to remove
explicit references to this subdirectory in his library code.
-Andrew 2019-12-04*


## Development Status: *Beta*

The source code in the ".hpp" header files are unlikely to change,
but *include paths* could change in the future.
(See [above](#Additional-Modifications-Needed).)

## License

*superpose3d_cpp* is available under the terms of the [MIT license](LICENSE.md).

