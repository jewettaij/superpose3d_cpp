/// @file     superpose3d.hpp
/// @brief    Calculate the optimal rotation, translation and scale needed to
///           optimally fit two different point clouds containing n points.
/// @author   Andrew Jewett
/// @license  MIT

#ifndef _SUPERPOSE3D_HPP
#define _SUPERPOSE3D_HPP



#include "matrix_alloc.hpp"

// Note: The Superpose3D::Superpose() function need to calculate the eigenvalues
// and eigenvectors of a 4x4 matrix.  Two methods: Lanczos or Jacobi:
#ifdef SUPERPOSE3D_USES_LANCZOS
// If you select this version, you must download "lambda_lanczos.hpp" from
// https://github.com/mrcdr/lambda-lanczos.  The "peigencalc_lanczos.hpp"
// file should be located in a different directory with this repository.
#include "peigencalc_lanczos.hpp"
#else
// DEFAULT:
// The ordinary Jacobi diagonalization code turned out to be much faster
// than LambdaLanczos for 4x4 matrices, so we use this method by default.
// If you select this version, you must download "jacobi.hpp" and
// "matrix_alloc.hpp" from https://github.com/jewettaij/jacobi_pd
#include "peigencalc.hpp"
#endif


namespace superpose3d {

using namespace matrix_alloc;

// -----------------------------------------------------------
// ------- some utilities I need to declare beforehand -------
// -----------------------------------------------------------
template<typename Scalar>
static inline Scalar SQR(Scalar x) {return x*x;}

// This is a stand-alone function invoked by Superpose3D::Superpose()
// but it does all the work.  This function was not intended for public use.
// (See implementation below for a description of the variables.)
template<typename Scalar, typename ConstArrayOfCoords>
inline static Scalar
_Superpose3D(size_t N, Scalar **aaRotate, Scalar *aTranslate,
             ConstArrayOfCoords aaXf_o, ConstArrayOfCoords aaXm_o,
             Scalar const *aWeights=nullptr, Scalar *pC=nullptr,
             PEigenCalculator<Scalar, Scalar*, Scalar const* const*> *pPE=nullptr,
             Scalar **aaXf_s=nullptr, Scalar **aaXm_s=nullptr);

// -----------------------------------------------------------
// ------------------------ INTERFACE ------------------------
// -----------------------------------------------------------


/// @brief  Superpose3d is a class with only one important member function
///         Superpose().  It is useful for repeatedly calculating the optimal
///         superposition (rotations, translations, and scale transformations)
///         between two point clouds of the same size.
template<typename Scalar, typename ConstArrayOfCoords>
class Superpose3D {
private:
  size_t N;              //number of points in the point clouds
  PEigenCalculator<Scalar, Scalar*, Scalar const* const*>
       eigen_calculator; // calculates principal eigenvalues
  // (contiguous) preallocated space for 2D arrays:
  Scalar **aaXf_shifted; //preallocated space for fixed point cloud (Nx3 array)
  Scalar **aaXm_shifted; //preallocated space for mobile point cloud (Nx3 array)

public:
  // The next 3 data members store the rotation, translation and scale
  // after optimal superposition
  Scalar **R;  //!< store optimal rotation here
  Scalar T[3]; //!< store optimal translation here
  Scalar c;  //!< store optimal scale (typically 1 unless requested by the user)

  Superpose3D(size_t N = 0);  //!< N = number of points in both point clouds

  ~Superpose3D();

  /// @brief specify he number of points in both point clouds
  void SetNumPoints(size_t N) {
    Dealloc();
    Alloc(N);
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
          ConstArrayOfCoords aaXf,        //!< coords for the "frozen" object
          ConstArrayOfCoords aaXm,        //!< coords for the "mobile" object
          Scalar const *aWeights=nullptr, //!< optional weights
          bool allow_rescale=false        //!< rescale mobile object? (c!=1?)
                   )
  {
    return
      _Superpose3D(N, R, T, aaXf, aaXm, aWeights, &c,
                   &eigen_calculator, aaXf_shifted, aaXm_shifted);
  }

  // memory management: copy and move constructor, swap, and assignment operator
  Superpose3D(const Superpose3D<Scalar, ConstArrayOfCoords>& source);
  Superpose3D(Superpose3D<Scalar, ConstArrayOfCoords>&& other);
  void swap(Superpose3D<Scalar, ConstArrayOfCoords> &other);
  Superpose3D<Scalar, ConstArrayOfCoords>& operator = (Superpose3D<Scalar, ConstArrayOfCoords> source);

private:

  // memory management:
  void Alloc(size_t N);
  void Init();
  void Dealloc();

}; // class Superpose3D





// -------------- IMPLEMENTATION --------------

template<typename Scalar, typename ConstArrayOfCoords>
static inline Scalar
_Superpose3D(size_t N,             // number of points in both point clouds
             Scalar **aaRotate,    // store rotation here
             Scalar *aTranslate,   // store translation here
             ConstArrayOfCoords aaXf_o, // coords for the "frozen" object
             ConstArrayOfCoords aaXm_o, // coords for the "mobile" object
             Scalar const *aWeights,    // optional weights
             Scalar *pC,           // rescale mobile object? if so store "c"here
             PEigenCalculator<Scalar, Scalar*, Scalar const* const*> *pPE, //!< eigenvalue calculator
             Scalar **aaXf_s,   // optional preallocated (Nx3) temporary array
             Scalar **aaXm_s)   // optional preallocated (Nx3) temporary array
{
  assert(aaRotate && aTranslate);
  assert(aaXf_o && aaXm_o);

  bool alloc_pPE = false;
  if (! pPE) {
    alloc_pPE = true;
    pPE = new PEigenCalculator<Scalar, Scalar*, Scalar const* const*>(4);
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

  bool alloc_aaXf_s = false;
  bool alloc_aaXm_s = false;
  if (! aaXf_s) {   // need to allocate the aaXf_s array?
    Alloc2D(N, 3, &aaXf_s);
    alloc_aaXf_s = true;
  }
  if (! aaXm_s) {   // need to allocate the aaXm_s array?
    Alloc2D(N, 3, &aaXm_s);
    alloc_aaXm_s = true;
  }
  assert(aaXf_s && aaXm_s);

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

  // Deallocate the temporary arrays we created earlier (if necessary).
  if (alloc_aaXf_s)
    Dealloc2D(&aaXf_s);
  if (alloc_aaXm_s)
    Dealloc2D(&aaXm_s);
  if (alloc_pPE)
    delete pPE;
  return rmsd;
}


template<typename Scalar, typename ConstArrayOfCoords>
Superpose3D<Scalar, ConstArrayOfCoords>::Superpose3D(size_t N)
  :eigen_calculator(4)
{
  Init();
  Alloc(N);
}

template<typename Scalar, typename ConstArrayOfCoords>
Superpose3D<Scalar, ConstArrayOfCoords>::~Superpose3D() {
  Dealloc();
}

template<typename Scalar, typename ConstArrayOfCoords>
void Superpose3D<Scalar, ConstArrayOfCoords>::
Init() {
  R = nullptr;
  aaXf_shifted = nullptr;
  aaXm_shifted = nullptr;
}

template<typename Scalar, typename ConstArrayOfCoords>
void Superpose3D<Scalar, ConstArrayOfCoords>::
Alloc(size_t N) {
  this->N = N;
  Alloc2D(3, 3, &R);
  Alloc2D(N, 3, &aaXf_shifted);
  Alloc2D(N, 3, &aaXm_shifted);
}

template<typename Scalar, typename ConstArrayOfCoords>
void Superpose3D<Scalar, ConstArrayOfCoords>::
Dealloc() {
  if (R)
    Dealloc2D(&R);
  if (aaXf_shifted)
    Dealloc2D(&aaXf_shifted);
  if (aaXm_shifted)
    Dealloc2D(&aaXm_shifted);
}

template<typename Scalar, typename ConstArrayOfCoords>
Superpose3D<Scalar, ConstArrayOfCoords>::
Superpose3D(const Superpose3D<Scalar, ConstArrayOfCoords>& source)
  :eigen_calculator(4, true)
{
  Init();
  Alloc(source.N);
  assert(N == source.npoints());
  for (int i = 0; i < N; i++) {
    std::copy(source.aaXf_shifted[i],
              source.aaXf_shifted[i] + N,
              aaXf_shifted[i]);
    std::copy(source.aaXm_shifted[i],
              source.aaXm_shifted[i] + N,
              aaXm_shifted[i]);
  }
}

template<typename Scalar, typename ConstArrayOfCoords>
void Superpose3D<Scalar, ConstArrayOfCoords>::
swap(Superpose3D<Scalar, ConstArrayOfCoords> &other) {
  std::swap(N, other.N);
  std::swap(R, other.R);
  std::swap(aaXf_shifted, other.aaXf_shifted);
  std::swap(aaXm_shifted, other.aaXm_shifted);
}

// Move constructor (C++11)
template<typename Scalar, typename ConstArrayOfCoords>
Superpose3D<Scalar, ConstArrayOfCoords>::
Superpose3D(Superpose3D<Scalar, ConstArrayOfCoords>&& other) {
  Init();
  swap(*this, other);
}

// Using the "copy-swap" idiom for the assignment operator
template<typename Scalar, typename ConstArrayOfCoords>
Superpose3D<Scalar, ConstArrayOfCoords>&
Superpose3D<Scalar, ConstArrayOfCoords>::
operator = (Superpose3D<Scalar, ConstArrayOfCoords> source) {
  this->swap(source);
  return *this;
}


} //namespace superposed3d



#endif //#ifndef _SUPERPOSE3D_HPP
