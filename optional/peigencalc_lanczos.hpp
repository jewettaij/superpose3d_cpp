/// @file     peigencalc_lanczos.hpp
/// @brief    Provide a means to calculate the largest (or smallest) eigenvalue
///           (and corresponding eigenvector) of a small dense square matrix.
/// @author   Andrew Jewett
/// @license  MIT
#ifndef _PEIGENCALC_HPP
#define _PEIGENCALC_HPP

#include <vector>

#include "lambda_lanczos.hpp"

/// @brief
///   PEigenCalculator is a class containing only one useful member function:
///   PrincipalEigen().  This function calculates the principal (largest
///   or smallest) eigenvalue and corresponding eigenvector of a square
///   n x n matrix.  This can be faster than diagionalizing the entire matrix.
///   (For example by using the Lanczos algorithm or something similar.)
/// @note
///   This code is a wrapper. Internally, it uses the "LambdaLanczos" class.
/// @note
///   For matrices larger than 13x13, PEigenCalculator::PrincipleEigen()
///   is usually faster than Jacobi::Diagonalize().)

template<typename Scalar, typename Vector, typename ConstMatrix>
class PEigenCalculator
{
  size_t n;                 // the size of the matrix
  std::vector<Scalar> evec; // preallocated vector

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
  auto matmul = [&](const std::vector<Scalar>& in, std::vector<Scalar>& out) {
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        out[i] += matrix[i][j]*in[j];
      }
    } 
  };
  auto init_vec = [&](std::vector<Scalar>& vec) {
    for(int i = 0; i < n; i++)
      vec[i] = 0.0;
    vec[0] = 1.0;
  };

  // "ll_engine" calculates the eigenvalue and eigenvector.
  LambdaLanczos<Scalar> ll_engine(matmul, n, find_max);
  
  // The Lanczos algorithm selects the eigenvalue with the largest magnitude.
  // In order to insure that this is the one we want (maxima or minima), we can
  // add a constant to all of the eigenvalues by setting "eigenvalue_offset".
  Scalar eval_upper_bound = 0.0;
  for (int i = 0; i < n; i++) {
    Scalar sum_row = 0.0;
    for (int j = 0; j < n; i++)
      sum_row += std::abs(matrix[i][j]);
    if (eval_upper_bound < sum_row)
      eval_upper_bound = sum_row;
  }
  if (find_max)
    ll_engine.eigenvalue_offset = eval_upper_bound;
  else
    ll_engine.eigenvalue_offset = -eval_upper_bound;

  ll_engine.init_vector = init_vec;

  Scalar eval;

  // This line does all of the hard work:
  size_t itern = ll_engine.run(eval, evec);

  for (int i = 0; i < n; i++)
    eigenvector[i] = evec[i];

  return eval;
}


#endif //#ifndef _PEIGENCALC_HPP
