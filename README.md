<<<<<<< HEAD
[![Build Status](https://travis-ci.org/jewettaij/superpose3d_cpp.svg?branch=master)](https://travis-ci.org/jewettaij/superpose3d_cpp.svg?branch=master)
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()
[![GitHub All Releases](https://img.shields.io/github/downloads/jewettaij/superpose3d_cpp/total)]()
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/jewettaij/superpose3d_cpp)]()
=======
[![Build Status](https://travis-ci.org/jewettaij/superpose3d.svg?branch=master)](./.travis.yml)
[![GitHub](https://img.shields.io/github/license/jewettaij/superpose3d)](./LICENSE.md)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/superpose3d)](https://pypistats.org/packages/moltemplate)
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/jewettaij/superpose3d)]()

>>>>>>> 87d35ac1174a02406891be4a71e5dbe4c136902f


superpose3d_cpp
===========


<<<<<<< HEAD
**superpose3d** is a header-only C++ library containing a defninition
of a class whose single public member function, *Superpose3d()*,
takes two N x 3 multidimensional arrays
(*of the same length*, **N**) representing points
from a point cloud (**X_i** and **x_i**)
as arguments.
Treating them as rigid objects,
*Superpose3D::Superpose3D()* attempts to superimpose
=======
```python
def Superpose3D(X_i,    # <-- Nx3 array of coords for the "frozen" point cloud
                x_i,    # <-- Nx3 array of coords for the "mobile" point cloud
                w_i=None, #<- optional weights for the calculation of RMSD
                          #   (default w_i = 1 for all i)
                allow_rescale=False)  #<--attempt to rescale mobile point cloud?
```

Superpose3D() takes two ordered lists (or numpy arrays) of xyz coordinates
(*of the same length*, **N**) representing points in a point cloud (**X_i** and
**x_i**). Treating them as rigid objects, "Superpose3D()" attempts to superimpose
>>>>>>> 87d35ac1174a02406891be4a71e5dbe4c136902f
them using **rotations**, **translations**, and (optionally) **scale**
transformations in order to minimize the root-mean-squared-distance (RMSD)
between corresponding points from either point cloud, where RMSD is defined as:
```
   RMSD = sqrt((Σ_i  w_i * |X_i - Σ_j(c*R_ij*x_j + T_i))|^2) / (Σ_j w_j))
```
If *w_i=nullptr*, equal weights are used.  In that case:
```
   RMSD = sqrt(( Σ_i |X_i - Σ_j(c*R_ij*x_j + T_i) )|^2 ) / N)
```
...where:
```
   T_j  = a translation vector (a 1-D numpy array containing x,y,z displacements),
   R_ij = a rotation matrix    (a 3x3 numpy array whose determinant = 1),
    c   = a scalar             (a number)
```

##  Example usage

```
#include <superpose3d.hpp>
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

  // Note: The optimal rotation, translation, and scale factor are stored in
  //       superposer.R, superposer.T, and superposer.c, respectively.

  // NOTE: If you want to specify the weights, then invoke it this way
  // superposer.Superpose(X, x, w);
  // If you want to allow scale transformations, then use
  // superposer.Superpose(X, x, w, true);
}
```

This function implements a more general variant of the method from this paper:
R. Diamond, (1988)
"A Note on the Rotational Superposition Problem",
 Acta Cryst. A44, pp. 211-216.

This version has been augmented slightly to support scale transformations.  (I.E. multiplication by scalars.  This can be useful for the registration of two different annotated volumetric 3-D images of the same object taken at different magnifications.)

Note that if you enable scale transformations (i.e. if *allow_rescale=true*), you should be wary if the function returns a negative **c** value.  Negative **c** values correspond to inversions (reflections).  For this reason, if you are using this function to compare the conformations of molecules, you should probably set *allow_rescale=false*.  This will prevent matching a molecule with its stereoisomer.


## Installation

This is a header-only library.
Just copy the directory containing "superpose3d.hpp" to a location in
your *include path*.

## License

superpose3d is available under the terms of the [MIT license](LICENSE.md).
