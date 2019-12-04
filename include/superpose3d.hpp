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

#include "peigencalc.hpp"
using namespace lambda_lanczos;


namespace superpose3d_lammps {


// ----------------------------------------------------------------
// ---------- some utilities I need to declare beforehand ---------
// ----------------------------------------------------------------


// Because I allocate 2-dimensional arrays frequently, I created a 
// few functions that make this more convenient.

/// @brief  Allocate a 2-dimensional table row-major order
template<typename Entry, typename Integer>
void Alloc2D(Integer const size[2], //!< size of the array in x,y directions
             Entry **paX,           //!< pointer to 1-D contiguous-memory array
             Entry ***paaX);        //!< pointer to 2-D multidimensional array

/// @brief
/// Slightly different version of Alloc2D()
/// In this version, the the size of the array specified by 2 integer arguments.
template<typename Entry>
void Alloc2D(size_t M,              //!< size of the array (outer)
             size_t N,              //!< size of the array (inner)
             Entry **paX,           //!< pointer to 1-D contiguous-memory array
             Entry ***paaX);        //!< pointer to 2-D multidimensional array

/// @brief
/// This function is the corresponding way to dellocate arrays
/// that were created using Alloc2D()
template<typename Entry>
void Dealloc2D(Entry **paX,          //!< pointer to 1-D contiguous-memory array
               Entry ***paaX);       //!< pointer to 2-D multidimensional array


template<typename Scalar>
static inline Scalar SQR(Scalar x) {return x*x;}


// This is a stand-alone function invoked by Superpose3D::Superpose()
// but it does all the work.  This function was not intended for public use.
// (See implementation below for a description of the variables.)
template<typename Scalar>
inline static Scalar
_Superpose3D(size_t N, Scalar **aaRotate, Scalar *aTranslate,
             Scalar const *const *aaXf_o, Scalar const *const *aaXm_o,
             Scalar const *aWeights=nullptr, Scalar *pC=nullptr,
             PEigenCalculator<Scalar> *pPE=nullptr,
             Scalar **aaXf_s=nullptr, Scalar **aaXm_s=nullptr);



// -----------------------------------------------------------
// ------------------------ INTERFACE ------------------------
// -----------------------------------------------------------


/// @brief  Superpose3d is a class with only one important member function
///         Superpose().  It is useful for repeatedly calculating the optimal
///         superposition (rotations, translations, and scale transformations)
///         between two point clouds of the same size.
template<typename Scalar>
class Superpose3D {
private:
  size_t N;              //number of points in the point clouds
  PEigenCalculator<Scalar> eigen_calculator; // calculates principal eigenvalues
  // (contiguous) preallocated space for 2D arrays:
  Scalar **aaXf_shifted; //preallocated space for fixed point cloud (Nx3 array)
  Scalar **aaXm_shifted; //preallocated space for mobile point cloud (Nx3 array)
  Scalar *aXf_shifted;   //contiguous memory allocated for aaXf_shifted
  Scalar *aXm_shifted;   //contiguous memory allocated for aaXm_shifted
  Scalar *_R;            //contiguous memory allocated for R

public:
  // The next 3 data members store the rotation, translation and scale
  // after optimal superposition
  Scalar **R;  //!< store optimal rotation here
  Scalar T[3]; //!< store optimal translation here
  Scalar c;  //!< store optimal scale (typically 1 unless requested by the user)

  Superpose3D(size_t N=0);  //!< N = number of points in both point clouds

  ~Superpose3D();

  /// @brief specify he number of points in both point clouds
  void SetNumPoints(size_t N) {
    Dealloc();
    Alloc(N);
    assert(this->N == N);
  }

  /// @brief return the number of points in both point clouds
  size_t GetNumPoints() {
    return N;
  }

  /// @brief
  /// Takes two lists of xyz coordinates (of the same length, specified earlier)
  /// and attempts to superimpose them using rotations, translations, and 
  /// (optionally) rescale operations are applied to the coordinates in the
  /// "aaXm_orig" array in order to minimize the root-mean-squared-distance
  /// (RMSD) between them, where RMSD is defined as:
  ///    sqrt((Sum_i  w_i * |X_i - (Sum_jc*R_ij*x_j + T_i))|^2) / (Sum_j w_j))
  /// The "X_i" and "x_i" are coordinates of the ith fixed and mobile point,
  /// (represented by "aaXf" and "aaXm" below), and "w_i" are weights
  /// (represented by "aWeights", which, if omitted, are assumed to be equal).
  /// This function implements a more general variant of the method from:
  /// R. Diamond, (1988)
  /// "A Note on the Rotational Superposition Problem", 
  /// Acta Cryst. A44, pp. 211-216
  /// This version has been augmented slightly.  The version in the original 
  /// paper only considers rotation and translation and does not allow the 
  /// coordinates of either object to be rescaled (multiplication by a scalar).
  /// You can enable the ability to rescale the coordinates by setting
  /// "allow_rescale" to true.  (By default, this feature is disabled.)
  ///
  /// @returns
  /// The RMSD between the 2 pointclouds after optimal rotation, translation
  /// (and scaling if requested) was applied to the "mobile" point cloud.
  /// After this function is called, the optimal rotation, translation,
  /// and scale (if requested) will be stored in the "R", "T", and "c"
  /// public data members.
  Scalar Superpose(
          Scalar const *const *aaXf,      //!< coords for the "frozen" object
          Scalar const *const *aaXm,      //!< coords for the "mobile" object
          Scalar const *aWeights=nullptr, //!< optional weights
          bool allow_rescale=false        //!< rescale mobile object? (c!=1?)
                   )
  {
    return
      _Superpose3D(N,
                   R,
                   T,
                   aaXf,
                   aaXm,
                   aWeights,
                   &c,
                   &eigen_calculator,
                   aaXf_shifted,
                   aaXm_shifted);
  }

private:

  void Alloc(size_t N);
  void Init();
  void Dealloc();

  // memory management: copy constructor, swap, and 
  Superpose3D(const Superpose3D<Scalar>& source);
  void swap(Superpose3D<Scalar> &other);
  Superpose3D<Scalar>& operator = (Superpose3D<Scalar> source);
}; // class Superpose3D





// -------------- IMPLEMENTATION --------------

template<typename Entry, typename Integer>
void Alloc2D(Integer const size[2], //!< size of the array in x,y directions
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

template<typename Entry>
void Alloc2D(size_t nrows,          //!< size of the array (outer)
             size_t ncolumns,       //!< size of the array (inner)
             Entry **paX,           //!< pointer to 1-D contiguous-memory array
             Entry ***paaX)         //!< pointer to 2-D multidimensional array
{
  size_t size[2];
  size[0] = ncolumns;
  size[1] = nrows;
  Alloc2D(size, paX, paaX);
}


template<typename Entry>
void Dealloc2D(Entry **paX,          //!< pointer to 1-D contiguous-memory array
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


template<typename Scalar>
Superpose3D<Scalar>::Superpose3D(size_t N)
  :eigen_calculator(4)
{
  Init();
  Alloc(N);
}

template<typename Scalar>
Superpose3D<Scalar>::~Superpose3D() {
  Dealloc();
}





template<typename Scalar>
static inline Scalar
_Superpose3D(size_t N,             // number of points in both point clouds
             Scalar **aaRotate,    // store rotation here
             Scalar *aTranslate,   // store translation here
             Scalar const *const *aaXf_o, // coords for the "frozen" object
             Scalar const *const *aaXm_o, // coords for the "mobile" object
             Scalar const *aWeights,         // optional weights
             Scalar *pC,      // rescale mobile object? if so store "c"here
             PEigenCalculator<Scalar> *pPE, //!< eigenvalue calculator
             Scalar **aaXf_s,   // optional preallocated (Nx3) temporary array
             Scalar **aaXm_s)   // optional preallocated (Nx3) temporary array
{
  assert(aaRotate && aTranslate);
  assert(aaXf_o && aaXm_o);

  bool alloc_pPE = false;
  if (! pPE) {
    alloc_pPE = true;
    pPE = new PEigenCalculator<Scalar>(4);
  }

  // Find the center of mass of each object:
  Scalar aCenter_f[3] = {0.0, 0.0, 0.0};
  Scalar aCenter_m[3] = {0.0, 0.0, 0.0};
  Scalar sum_weights = 0.0;
  for (size_t n=0; n < N; n++) {
    Scalar weight = 1.0;
    if (aWeights)
      weight = aWeights[n];
    for (int d=0; d < 3; d++) {
      aCenter_f[d] += aaXf_o[n][d]*weight;
      aCenter_m[d] += aaXm_o[n][d]*weight;
    }
    sum_weights += weight;
  }
  for (int d=0; d < 3; d++) {
    aCenter_f[d] /= sum_weights;
    aCenter_m[d] /= sum_weights;
  }

  //Subtract the centers-of-mass from the original coordinates for each object

  // Create some temporary arrays we will need aaXf_s, aaXm_s
  //aaXf_s[i][d] = dth coord of ith particle in "fixed" object
  //aaXm_s[i][d] = dth coord of ith particle in "mobile" object

  // Ugly details:  we might need to allocate space for these arrays
  Scalar *aXf_s=nullptr; //storage space for aaXf_s; aaXf_s[i][d] = aXf_s[3*i+d]
  Scalar *aXm_s=nullptr; //storage space for aaXm_s; aaXm_s[i][d] = aXm_s[3*i+d]

  // decide whether to allocate the aaXf_s array
  if (! aaXf_s) {
    Alloc2D(N, 3, &aXf_s, &aaXf_s);
    assert(aXf_s);
  }
  assert(aaXf_s);
  // decide whether to allocate the aaXm_s array
  if (! aaXm_s) {
    Alloc2D(N, 3, &aXm_s, &aaXm_s);
    assert(aXm_s);
  }
  assert(aaXm_s);
  for (size_t n=0; n < N; n++) {
    for (int d=0; d < 3; d++) {
      // shift the coordinates so that the new center of mass is at the origin
      aaXf_s[n][d] = aaXf_o[n][d] - aCenter_f[d];
      aaXm_s[n][d] = aaXm_o[n][d] - aCenter_m[d];
    }
  }

  bool allow_rescale = pC != nullptr;

  Scalar Rgf=0.0;// <--the RMS size of the particles in the frozen object aaXf_s
  Scalar Rgm=0.0;// <--the RMS size of the particles in the mobile object aaXm_s

  if (allow_rescale) {
    // Optional: For numerical stability, we might as well rescale the
    // coordinates initially to make sure they have the same approximate
    // scale before we attempt to superimpose them.
    // This is probably only useful if one object is much bigger than the other.
    // Note: This is NOT the optimal scale factor.
    //       (That must be determined later.)
    for (size_t n=0; n < N; n++) {
      Scalar weight = 1.0;
      if (aWeights)
        weight = aWeights[n];
      for (int d=0; d < 3; d++) {
        Rgf += weight * SQR(aaXf_s[n][d]);
        Rgm += weight * SQR(aaXm_s[n][d]);
      }
    }
    Rgf = sqrt(Rgf / sum_weights);
    Rgm = sqrt(Rgm / sum_weights);

    for (size_t n=0; n < N; n++) {
      for (int d=0; d < 3; d++) {
        aaXf_s[n][d] /= Rgf;
        aaXm_s[n][d] /= Rgm;
      }
    }
  } //if (allow_rescale)

  // Calculate the "M" array from the Diamond paper (equation 16)
  Scalar M[3][3];
  for (int i=0; i < 3; i++)
    for (int j=0; j < 3; j++)
      M[i][j] = 0.0;

  for (size_t n=0; n < N; n++) {
    Scalar weight = 1.0;
    if (aWeights)
      weight = aWeights[n];
    for (int i=0; i < 3; i++) {
      for (int j=0; j < 3; j++) {
        M[i][j] += weight * aaXm_s[n][i] * aaXf_s[n][j];
      }
    }
  }

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
  Scalar _P[4*4]; // contiguous 1D array for storing contents of P
  Scalar *P[4];      // 2D array (in a format compatible with matrix solvers)
  for (int i=0; i < 4; i++)
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

  // The vector "p" will contain the optimal rotation (in quaternion format)
  Scalar p[4];
  assert(pPE);
  Scalar eval_max = pPE->PrincipalEigen(P, p, true);

  // Now normalize p
  Scalar pnorm = 0.0;
  for (int i=0; i < 4; i++)
    pnorm += p[i]*p[i];
  pnorm = sqrt(pnorm);
  for (int i=0; i < 4; i++)
    p[i] /= pnorm;

  // Finally, calculate the rotation matrix corresponding to "p"
  // (convert a quaternion into a 3x3 rotation matrix)

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

  if (pC) {
    // If the user supplies a non-NULL pC argument, then they want
    // to calculate the optimal scaling factor, c, (and store it in *pC).
    Scalar Waxaixai = 0.0;
    Scalar WaxaiXai = 0.0;
    for (size_t a=0; a < N; a++) {
      Scalar weight = 1.0;
      if (aWeights)
        weight = aWeights[a];
      for (int i=0; i < 3; i++) {
        Waxaixai += weight * aaXm_s[a][i] * aaXm_s[a][i];
        WaxaiXai += weight * aaXm_s[a][i] * aaXf_s[a][i];
      }
    }
    c = (WaxaiXai + pPp) / Waxaixai;

    // Recall that we previously divided the two sets of coordinates by Rgm
    // and Rgf respectively.(I thought it might improve numerical stability)
    // Before returning "c" to the caller, we need to incorporate those
    // factors into "c" as well.
    c *= Rgf / Rgm;
        // And, lastly, undo this before calculating E0 below
    for (size_t n=0; n < N; n++) {
      for (int d=0; d < 3; d++) {
        aaXf_s[n][d] *= Rgf;
        aaXm_s[n][d] *= Rgm;
      }
    }
    *pC = c;
  } // if (pC)


  // Finally compute the RMSD between the two coordinate sets:
  // First compute E0 from equation 24 of the paper
  Scalar E0 = 0.0;
  for (size_t n=0; n < N; n++) {
    Scalar weight = 1.0;
    if (aWeights)
      weight = aWeights[n];
    for (int d=0; d < 3; d++)
      // (remember to include the scale factor "c" that we inserted)
      E0 += weight * (SQR(aaXf_s[n][d] - c*aaXm_s[n][d]));
  }
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
  if (aXf_s)
    Dealloc2D(&aXf_s, &aaXf_s);
  if (aXm_s)
    Dealloc2D(&aXm_s, &aaXm_s);
  if (alloc_pPE)
    delete pPE;
  return rmsd;
}


template<typename Scalar>
void Superpose3D<Scalar>::Init() {
  _R = nullptr;
  R = nullptr;
  aaXf_shifted = nullptr;
  aaXm_shifted = nullptr;
  aXf_shifted = nullptr;
  aXm_shifted = nullptr;
}

template<typename Scalar>
void Superpose3D<Scalar>::Alloc(size_t N) {
  this->N = N;
  Alloc2D(3, 3, &_R, &R);
  Alloc2D(N, 3, &aXf_shifted, &aaXf_shifted);
  Alloc2D(N, 3, &aXm_shifted, &aaXm_shifted);
  assert(aXf_shifted && aXm_shifted);
}

template<typename Scalar>
void Superpose3D<Scalar>::Dealloc() {
  if (R)
    Dealloc2D(&_R, &R);
  if (aaXf_shifted) {
    assert(aXf_shifted && aaXm_shifted && aXm_shifted);
    Dealloc2D(&aXf_shifted, &aaXf_shifted);
    Dealloc2D(&aXm_shifted, &aaXm_shifted);
  }
}

template<typename Scalar>
Superpose3D<Scalar>::Superpose3D(const Superpose3D<Scalar>& source)
  :eigen_calculator(4,true)
{
  Init();
  Alloc(source.N);
  assert(N == source.npoints());
  std::copy(source.aXf_shifted,
            source.aXf_shifted + N*3,
            aXf_shifted);
  std::copy(source.aXm_shifted,
            source.aXm_shifted + N*3,
            aXm_shifted);
}

template<typename Scalar>
void Superpose3D<Scalar>::swap(Superpose3D<Scalar> &other) {
  std::swap(aXf_shifted, other.aXf_shifted);
  std::swap(aXm_shifted, other.aXm_shifted);
  std::swap(aaXf_shifted, other.aaXf_shifted);
  std::swap(aaXm_shifted, other.aaXm_shifted);
  std::swap(N, other.N);
}

template<typename Scalar>
Superpose3D<Scalar>&
Superpose3D<Scalar>::operator = (Superpose3D<Scalar> source) {
  this->swap(source);
  return *this;
}


} //namespace superposed3d_lammps



#endif //#ifndef _SUPERPOSE3D_HPP
