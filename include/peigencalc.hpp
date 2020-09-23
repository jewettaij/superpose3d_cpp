/// @file     peigencalc.hpp
/// @brief    Provide a means to calculate the largest (or smallest) eigenvalue
///           (and corresponding eigenvector) of a small dense square matrix.
/// @author   Andrew Jewett
/// @license  MIT

#ifndef _PEIGENCALC_HPP
#define _PEIGENCALC_HPP

#include <vector>
#include "jacobi_pd.hpp"
using namespace jacobi_pd;

/// @brief PEigenCalculator is a class containing only one useful member
/// function PrincipalEigen().  This function calculates the principal (largest
/// or smallest) eigenvalue and corresponding eigenvector of an n x n matrix.
/// In some cases, this is faster than diagionalizing the entire matrix.
/// @note:
/// This particular version of the PEigenCalculator is simply a wrapper for the
/// "jacobi_pd::Jacobi" class which calculates ALL the eigenvalues.
/// I only intended to apply it to 4x4 matrices, and for this size, it turns out
/// the ordinary Jacobi eigenvalue algorithm is actually faster than some
/// approaches that generate only one eigenvalue, such as "LambdaLancos"
/// (https://github.com/mrcdr/lambda-lanczos)
/// Other code (like "superpose3d.hpp") which needs to calculate
/// eigenvalues and eigenvectors should not need to be concerned with how
/// the calculation is implemented.  (Later on, if I choose to implement it
/// a different way, I won't have to modify "superpose3d.hpp".)

template<typename Scalar, typename Vector, typename ConstMatrix>

class PEigenCalculator
{
  size_t n;                         // the size of the matrices to be analyzed
  std::vector<Scalar> evals;          // preallocated array for the eigenvalues
  std::vector<std::vector<Scalar> > evecs; //preallocated array for the evectors
  Jacobi<Scalar,
         std::vector<Scalar>&,
         std::vector<std::vector<Scalar> >&,
         ConstMatrix> ecalc;
public:

  /// @brief  Specify the size of the matrices you want to diagonalize later.
  /// @param  matrix_size  the number of rows (or columns) of the matrix
  void SetSize(int matrix_size);

  PEigenCalculator(int matrix_size=0) {
    SetSize(matrix_size);
  }

  /// @brief  Calculate the principal eigenvalue and eigenvector of a matrix.
  /// @return Return the principal eigenvalue of the matrix.
  ///         If you want the eigenvector, pass a non-null "evector" argument.
  Scalar
  PrincipalEigen(ConstMatrix matrix,   //!< the input patrix
                 Vector evector,       //!< the eigenvector is stored here
                 bool find_max=false); //!< want the max or min eigenvalue?

}; // class PEigenCalculator



// -------- IMPLEMENTATION --------

template<typename Scalar, typename Vector, typename ConstMatrix>
Scalar PEigenCalculator<Scalar, Vector, ConstMatrix>::
  PrincipalEigen(ConstMatrix matrix,
                 Vector eigenvector,
                 bool find_max)
{
  assert(n > 0);
  ecalc.Diagonalize(matrix,
                    evals,
                    evecs,
                    Jacobi<Scalar,          //<--specify the sorting criteria
                           std::vector<Scalar>&,
                           std::vector<std::vector<Scalar> >&,
                           ConstMatrix>::SORT_DECREASING_EVALS);

  int which_eigenvalue = n-1;
  if (find_max)
    which_eigenvalue = 0;

  if (eigenvector) {
    // If the caller requested the eigenvector as well, then
    // return it to the caller by copying the data into eigenvector[].
    for (int i = 0; i < n; i++)
      eigenvector[i] = evecs[which_eigenvalue][i];
  }

  return evals[which_eigenvalue];
}


template<typename Scalar, typename Vector, typename ConstMatrix>
void PEigenCalculator<Scalar, Vector, ConstMatrix>::
SetSize(int matrix_size) {
  n = matrix_size;
  ecalc.SetSize(n);
  evals.resize(n);
  evecs.resize(n);
  for (int i = 0; i < n; i++)
    evecs[i].resize(n);
}

#endif //#ifndef _PEIGENCALC_HPP
