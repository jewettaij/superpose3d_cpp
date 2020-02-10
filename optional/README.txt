This directory contains an alternative version the "peigencalc.hpp" file
distributed with "superpose3d.hpp".
It calculates the principle eigenvalue and eigenvector of a 4x4 matrix,
using 3rd-party code to do the hard work.  It uses the Lanczos algorithm from:
https://github.com/mrcdr/lambda-lanczos
This was the original eigenvalue/vector calculator used by superpose3d.hpp.
This implementation turned out to be slow for 4x4 matrices.  I currently
use a faster method for calculating eigenvalues, but I kept the original
code here in case anyone prefers this implementation.
There is no reason to keep this file.  I may delete it in future versions.
