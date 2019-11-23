///   @file  superpose3d.hpp
///   @brief Calculate the optimal rotation, translation and scale needed to
///          optimally fit two different point clouds containing n points.
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Andrew Jewett (Scripps Research)
------------------------------------------------------------------------- */


#ifndef _SUPERPOSE3D_HPP
#define _SUPERPOSE3D_HPP



namespace superpose3d_lammps {


/// @brief
/// Allocate a 2-dimensional table row-major order

template<typename Entry, typename Integer>
static inline void
Alloc2D(Integer const size[2], //!< size of the array in x,y directions
        Entry **paX,           //!< pointer to 1-D contiguous-memory array
        Entry ***paaX)         //!< pointer to 2-D multidimensional array
{
  assert(paX && paaX);

  *paX = new Entry [size[0] * size[1]];

  // Allocate a conventional 2-dimensional
  // pointer-to-a-pointer data structure.
  *paaX = new Entry* [size[1]];
  for(Integer iy=0; iy<size[1]; iy++)
    (*paaX)[iy] = &((*paX)[ iy*size[0] ]);
  // The caller can access the contents of *paX using (*paaX)[i][j] notation.
}


/// @brief
/// Slightly different version of Alloc2D()
/// In this version, the the size of the array specified by 2 integer arguments.
template<typename Entry>
static inline void
Alloc2D(size_t M,              //!< size of the array (outer)
        size_t N,              //!< size of the array (inner)
        Entry **paX,           //!< pointer to 1-D contiguous-memory array
        Entry ***paaX)         //!< pointer to 2-D multidimensional array
{
  size_t size[2];
  size[0] = M;
  size[1] = N;
  Alloc2D(size, paX, paaX);
}


/// @brief
/// This function is the corresponding way to dellocate arrays
/// that were created using Alloc2D()

template<typename Entry, typename Integer>
static inline void
Dealloc2D(Integer const size[2], //!< size of the array in x,y directions
          Entry **paX,          //!< pointer to 1-D contiguous-memory array
          Entry ***paaX)        //!< pointer to 2-D multidimensional array
{
  if (paaX && *paaX) {
    delete [] (*paaX);
    *paaX = nullptr;
  }
  if (paX && *paX) {
    delete [] *paX;
    *paX = nullptr;
  }
}

/// @brief
/// Slightly different version of Dealloc2D()
/// In this version, the the size of the array specified by 2 integer arguments.
template<typename Entry>
static inline void
Dealloc2D(size_t M,           //!< size of the array (outer)
          size_t N,           //!< size of the array (inner)
          Entry **paX,        //!< pointer to 1-D contiguous-memory array
          Entry ***paaX)      //!< pointer to 2-D multidimensional array
{
  size_t size[2];
  size[0] = M;
  size[1] = N;
  Dealloc2D(size, paX, paaX);
}


/// @brief
/// Superpose() is the stand-alone function invoked by
/// Superpose3D::Superpose().  You can use this function by itself without
/// creating an instance of the "Superpose3D" class, however you incur a
/// slow down if you do not preallocate the temporary arrys it needs.
  

template<typename Scalar>
static inline Scalar
Superpose(Scalar **aaRotate,  //!< store rotation here
          Scalar *aTranslate, //!< store translation here
          size_t N,           //!< number of points in both point clouds
          Scalar const *const *aaXf_orig, //!< coords for the "frozen" object
          Scalar const *const *aaXm_orig, //!< coords for the "mobile" object
          Scalar const *aWeights=nullptr, //!< optional weights
          Scalar *pC=nullptr, //!< rescale mobile object? if so store "c"here
          Scalar **aaXf, //optional preallocated aaXf array (Nx3 array)
          Scalar **aaXm  //optional preallocated aaXm array (Nx3 array)
          )
{
  assert(aaRotate && aTranslate);
  assert(aaXf_orig && aaXm_orig);

  bool alloc_aWeights = false;
  if (aWeights == nulptr) {
    aWeights = new Scalar[N];
    alloc_aWeights = true;
  }

  // Find the center of mass of each object:
  Scalar aCenter_f[3] = {0.0, 0.0, 0.0};
  Scalar aCenter_m[3] = {0.0, 0.0, 0.0};
  sum_weights = 0.0;
  for (size_t n=0; n < N; n++) {
    for (int d=0; d < 3; d++) {
      aCenter_f[d] += aaXf_orig[n][d]*aWeights[n];
      aCenter_m[d] += aaXm_orig[n][d]*aWeights[n];
      sum_weights += aWeights[n];
    }
  }
  for (int d=0; d < 3; d++) {
    aCenter_f[d] /= sum_weights;
    aCenter_m[d] /= sum_weights;
  }

  //Subtract the centers-of-mass from the original coordinates for each object

  // Create some temporary arrays we will need aaXf, aaXm
  //aaXf[i][d] = dth coord of ith particle in "fixed" object
  //aaXf[i][d] = dth coord of ith particle in "mobile" object

  // Ugly details:  we might need to allocate space for these arrays
  Scalar *_aXf = nullptr; //storage space for aaXf; aaXf[i][d] = _aXf[3*i+d]
  Scalar *_aXm = nullptr; //storage space for aaXm; aaXm[i][d] = _aXm[3*i+d]

  // decide whether to allocate the aaXf array
  if (! aaXf) {
    Alloc2D(N, 3, &_aXf, &aaXf);
    assert(_aXf);
  }
  assert(aaXf);
  // decide whether to allocate the aaXm array
  if (! aaXm) {
    Alloc2D(N, 3, &_aXm, &aaXm);
    assert(_aXm);
  }
  assert(aaXm);
  for (size_t n=0; n < N; n++) {
    for (int d=0; d < 3; d++) {
      aaXf[n][d] = aaXf_orig[n][d] - aCenter_f[d];
      aaXm[n][d] = aaXm_orig[n][d] - aCenter_m[d];
    }
  }

  if (allow_rescale) {
    // Optional: For numerical stability, we might as well rescale the
    // coordinates initially to make sure they have the same approximate
    // scale before we attempt to superimpose them.
    // This is only necessary if one object is much bigger than the other
    // (ie. by several orders of magnitude).
    // Note: This is NOT the optimal scale factor.
    //       (That must be determined later.)
    Rgf = 0.0;  // <-- the RMS size of the particles in the frozen object aaXf
    Rgm = 0.0;  // <-- the RMS size of the particles in the mobile object aaXm
    for (size_t n=0; n < N; n++) {
      for (int d=0; d < 3; d++) {
        Rgf += aWeights[n]*((aaXf[n][d])**2);
        Rgm += aWeights[n]*((aaXm[n][d])**2);
      }
    }
    Rgf = sqrt(Rgf / sum_weights);
    Rgm = sqrt(Rgm / sum_weights);

    for (size_t n=0; n < N; n++) {
      for (int d=0; d < 3; d++) {
        aaXf[n][d] /= Rgf;
        aaXm[n][d] /= Rgm;
      }
    }
  } //if (allow_rescale)

  // Calculate the "M" array from the Diamond paper (equation 16)
  Scalar M[3][3];
  
  for (size_t n=0; n < N; n++)
    for (int i=0; i < 3; i++)
      for (int j=0; j < 3; j++)
        M[i][j] += aWeights[n] * aaXm[n][i] * aaXf[n][j];

  // Calculate Q (equation 17)
  Scalar traceM = 0.0;
  for (int i=0; i < 3; i++)
    traceM += M[i][i];
  Scalar Q[3][3];
  for (int i=0; i < 3; i++) {
    for (int j=0; j < 3; j++) {
      Q[i][j] = M[i][j] + M[j][i];
      if (i==j)
        Q[i][j] -= 2.0 * traceM;
    }
  }

  // Calculate V (equation 18)
  Scalar V[3];
  V[0] = M[1][2] - M[2][1];
  V[1] = M[2][0] - M[0][2];
  V[2] = M[0][1] - M[1][0];

  // Calculate "P" (equation 22)
  Scalar *_P[4*4]; // contiguous 1D array for storing contents of P
  Scalar **P;      // 2D array (in a format compatible with matrix solvers)
  for (int i=0; i < 3; i++)
    P[i] = &(_P[4*i]);

  for (int i=0; i < 3; i++)
    for (int j=0; j < 3; j++)
      P[i][j] = Q[i][j];

  P[0][3] = V[0];
  P[3][0] = V[0];
  P[1][3] = V[1];
  P[3][1] = V[1];
  P[2][3] = V[2];
  P[3][2] = V[2];
  P[3][3] = 0.0;

  aEigenvals, aaEigenvects = LA.eigh(P);
  #http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.eigh.html

  eval_max = aEigenvals[0];
  i_eval_max = 0;
  for (int i=1; i < 4; i++) {
    if (aEigenvals[i] > eval_max) {
      eval_max = aEigenvals[i];
      i_eval_max = i;
    }
  }

  // The vector "p" contains the optimal rotation (in quaternion format)
  Scalar p[4];
  p[0] = aaEigenvects[0][ i_eval_max ];
  p[1] = aaEigenvects[1][ i_eval_max ];
  p[2] = aaEigenvects[2][ i_eval_max ];
  p[3] = aaEigenvects[3][ i_eval_max ];

  Scalar pnorm = 0.0;
  for (int i=0; i < 4; i++)
    pnorm += p[i]*p[i];
  pnorm = sqrt(pnorm);
  for (int i=0; i < 4; i++)
    p[i] /= pnorm;

  // Finally, calculate the rotation matrix corresponding to "p"
  // (convert a quaternion into a 3x3 rotation matrix)
  Scalar aaRotate[3][3];

  aaRotate[0][0] =  (p[0]*p[0])-(p[1]*p[1])-(p[2]*p[2])+(p[3]*p[3]);
  aaRotate[1][1] = -(p[0]*p[0])+(p[1]*p[1])-(p[2]*p[2])+(p[3]*p[3]);
  aaRotate[2][2] = -(p[0]*p[0])-(p[1]*p[1])+(p[2]*p[2])+(p[3]*p[3]);
  aaRotate[0][1] = 2*(p[0]*p[1] - p[2]*p[3]);
  aaRotate[1][0] = 2*(p[0]*p[1] + p[2]*p[3]);
  aaRotate[1][2] = 2*(p[1]*p[2] - p[0]*p[3]);
  aaRotate[2][1] = 2*(p[1]*p[2] + p[0]*p[3]);
  aaRotate[0][2] = 2*(p[0]*p[2] + p[1]*p[3]);
  aaRotate[2][0] = 2*(p[0]*p[2] - p[1]*p[3]);
    
  Scalar pPp = eval_max;

  // Optional: Decide the scale factor, c
  Scalar c = 1.0;   // by default, don't rescale the coordinates

  if (*pC) {
    // If the user supplies a non-NULL pC argument, then they want
    // to calculate the optimal scaling factor, c, (and store it in *pC).
    Scalar Waxaixai = 0.0;
    Scalar WaxaiXai = 0.0;
    for (size_t a=0; a < N; a++) {
      for (int i=0; i < 3; i++) {
        Waxaixai += aWeights[a] * aaXm[a][i] * aaXm[a][i];
        WaxaiXai += aWeights[a] * aaXm[a][i] * aaXf[a][i];
      }
    }
    Scalar c = (WaxaiXai + pPp) / Waxaixai;

    // Recall that we previously divided the two sets of coordinates by Rgm
    // and Rgf respectively.(I thought it might improve numerical stability)
    // Before returning "c" to the caller, we need to incorporate those
    // factors into "c" as well.
    c *= Rgf / Rgm;
        // And, lastly, undo this before calculating E0 below
    for (size_t n=0; n < N; n++) {
      for (int d=0; d < 3; d++) {
        aaXf[n][d] *= Rgf;
        aaXm[n][d] *= Rgm;
      }
    }
    *pC = c;
  } // if (*pC)


  // Finally compute the RMSD between the two coordinate sets:
  // First compute E0 from equation 24 of the paper
  Scalar E0 = 0.0;
  for (size_t n=0; n < N; n++)
    for (int d=0; d < 3; d++)
      // (remember to include the scale factor "c" that we inserted)
      E0 += aWeights[n] * ((aaXf[n][d] - c*aaXm[n][d])**2);
  Scalar sum_sqr_dist = E0 - 2.0*pPp;
  if (sum_sqr_dist < 0.0)
    sum_sqr_dist = 0.0;
  Scalar rmsd = sqrt(sum_sqr_dist/sum_weights);

  // Lastly, calculate the translational offset:
  // Recall that:
  //RMSD=sqrt((Sum_i  w_i * |X_i - Sum_j(c*R_ij*x_j + T_i))|^2) / (Sum_j w_j))
  //    =sqrt((Sum_i  w_i * |X_i - x_i')|^2) / (Sum_j w_j))
  //  where
  // x_i' = Sum_j(c*R_ij*x_j) + T_i
  //      = Xcm_i + c*R_ij*(x_j - xcm_j)
  //  and Xcm and xcm = center_of_mass for the frozen and mobile point clouds
  //                  = aCenter_f[]       and       aCenter_m[],  respectively
  // Hence:
  //  T_i = Xcm_i - Sum_j c*R_ij*xcm_j  =  aTranslate[i]

  Scalar aTranslate[3];
  for (int i=0; i < 3; i++) {
    aTranslate[i] = aCenter_f[i];
    for (int j=0; j < 3; j++) {
      aTranslate[i] -= c*aaRotate[i][j]*aCenter_m[j];
    }
  }

  // An alternate method to compute "aTranslate" using matrices:
  //Rmatrix = np.matrix(aaRotate)
  //TcolumnVec = np.matrix(np.empty((3,1))) //3x1 matrix<->[[0],[0],[0]]
  //for d in range(0,3):
  //    TcolumnVec[d][0] = -aCenter_m[d]
  //TcolumnVec = c * Rmatrix * TcolumnVec
  //for d in range(0,3):
  //   TcolumnVec[d][0] += aCenter_f[d]
  // #Turn the column vector back into an ordinary numpy array of size 3:
  //aTranslate = np.array(TcolumnVec.transpose())[0]

  // Deallocate the temporary arrays we created earlier (if necessary).
  if (_aXf)
    Dealloc2D(N, 3, &_aXf, &aaXf);
  if (_aXm)
    Dealloc2D(N, 3, &_aXm, &aaXm);
  if (alloc_aWeights)
    delete [] aWeights;

  return rmsd;
}



/// @brief  Superpose3d is a class with only one important member function
///         Superpose().  It is useful for repeatedly calculating the optimal
///         superposition (rotations, translations, and scale transformations)
///         between two point clouds of the same size.
template<typename Scalar>
class Superposer3D {
private:
  size_t N;              //number of points in the point clouds
  Scalar **aaXf_shifted; //preallocated space for fixed point cloud (Nx3 array)
  Scalar **aaXm_shifted; //preallocated space for mobile point cloud (Nx3 array)
  // preallocated space for 2D arrays:
  Scalar *aXf_shifted;   //contiguous memory allocated for aaXf_shifted
  Scalar *aXm_shifted;   //contiguous memory allocated for aaXm_shifted
  Scalar _R[9];          //contiguous memory allocated for R

  void Alloc(_N) {
    N = _N;
    Alloc2D(N, 3, &aXf_shifted, &aaXf_shifted);
    Alloc2D(N, 3, &aXm_shifted, &aaXm_shifted);
    assert(aXf && aXm);
  }

  void Dealloc() {
    Dealloc2D(aXf_shifted, aaXf_shifted);
    Dealloc2D(aXm_shifted, aaXm_shifted);
  }

public:
  // The next 3 data members store the rotation, translation and scale
  // after optimal superposition
  Scalar **R;  // store optimal rotation here
  Scalar T[3]; // store optimal translation here
  Scalar c;    // store optimal scale (typically 1,unless requested by the user)

  Superpose3D(size_t _N) //!< number of points in both point clouds
  {
    Alloc(N);
    R[0] = &(_R[0]);
    R[1] = &(_R[3]);
    R[2] = &(_R[6]);
  }  

  ~Superpose3D() {
    Dealloc();
  }

  size_t size() {
    return N;
  }

  // copy constructor
  Superpose3D(const Superpose3D<Scalar>& source) {
    Init();
    Resize(source.size()); // allocates and initializes afXf and afXm
    assert(N == source.size());
    std::copy(source.afXf,
              source.afXf + N*3,
              afXf);
    std::copy(source.afXm,
              source.afXm + N*3,
              afXm);
  }

  /// @brief
  /// Takes two lists of xyz coordinates, (of the same length)
  /// and attempts to superimpose them using rotations, translations, and 
  /// (optionally) rescale operations are applied to the coordinates in the
  /// "aaXm_orig" array in order to minimize the root-mean-squared-distance
  /// (RMSD) between them.
  /// This function implements a more general variant of the method from:
  /// R. Diamond, (1988)
  /// "A Note on the Rotational Superposition Problem", 
  /// Acta Cryst. A44, pp. 211-216
  /// This version has been augmented slightly.  The version in the original 
  /// paper only considers rotation and translation and does not allow the 
  /// coordinates of either object to be rescaled (multiplication by a scalar).
  ///
  /// @returns
  /// The RMSD between the 2 pointclouds after optimal rotation, translation
  /// (and scaling if requested) was applied to the "mobile" point cloud.
  /// After this function is called, the optimal rotation, translation,
  /// and scale (if requested) will be stored in the "R", "T", and "c"
  /// public data members.
  Scalar Superpose(
          Scalar const *const *aaXf_orig, //!< coords for the "frozen" object
          Scalar const *const *aaXm_orig, //!< coords for the "mobile" object
          Scalar const *aWeights=nullptr, //!< optional weights
          bool allow_rescale=false        //!< rescale mobile object?
                   )
  {
    return 
      Superpose(Scalar **aaRotate,  //!< store rotation here
          Scalar *aTranslate, //!< store translation here
          N,
          Scalar const *const *aaXf_orig, //!< coords for the "frozen" object
          Scalar const *const *aaXm_orig, //!< coords for the "mobile" object
          Scalar const *aWeights=nullptr, //!< optional weights
          Scalar *pC=nullptr, //!< rescale mobile object? if so store "c"here
                aaXf_shifted,
                aaXm_shifted);
  }

}; // class Superpose3D

} //namespace superpose3d_lammps



#endif //#ifndef _SUPERPOSE3D_HPP
