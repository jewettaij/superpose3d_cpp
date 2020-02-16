/// @file     peigencalc_lanczos.hpp
/// @brief    Provide a means to calculate the largest (or smallest) eigenvalue
///           (and corresponding eigenvector) of a small dense square matrix.
/// @author   Andrew Jewett
/// @license  MIT
#ifndef _PEIGENCALC_HPP
#define _PEIGENCALC_HPP

#include <vector>
using std::vector;


#include "lambda_lanczos.hpp"

/// @brief PEigenCalculator is a class containing only one useful member
/// function PrincipalEigen().  This function calculates the principal (largest
/// or smallest) eigenvalue and corresponding eigenvector of an n x n matrix.
/// This can be faster than diagionalizing the entire matrix.
/// This code is a wrapper.  Internally, it uses the "LambdaLanczos" class.
/// (which is more general and can work with large sparse matrices.  You can
///  download that code here: https://github.com/mrcdr/lambda-lanczos)

template<typename Scalar, typename Vector, typename ConstMatrix>
class PEigenCalculator
{
  size_t n;            // the size of the matrix
  vector<Scalar> evec; // preallocated vector (lambda_lanzcos does not use ptrs)

public:
  void SetSize(int matrix_size) {
    n = matrix_size;
    evec.resize(n);
  }

  PEigenCalculator(int matrix_size=0):evec(matrix_size) {
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
  auto matmul = [&](const vector<Scalar>& in, vector<Scalar>& out) {
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        out[i] += matrix[i][j]*in[j];
      }
    } 
  };
  auto init_vec = [&](vector<Scalar>& vec) {
    for(int i = 0; i < n; i++)
      vec[i] = 0.0;
    vec[0] = 1.0;
  };

  Scalar eval;
  // (The next two lines do all the hard work.)
  lambda_lanczos::LambdaLanczos<Scalar> ll_engine(matmul, n, find_max);
  ll_engine.init_vector = init_vec;
  size_t itern = ll_engine.run(eval, evec);

  if (eigenvector) {
    // If the caller requested the eigenvector as well, then
    // return it to the caller by copying the data into eigenvector[].
    for (int i = 0; i < n; i++)
      eigenvector[i] = evec[i];
  }

  return eval;
}


#endif //#ifndef _PEIGENCALC_HPP
