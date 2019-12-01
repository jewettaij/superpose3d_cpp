#ifndef _PEIGENCALC_HPP
#define _PEIGENCALC_HPP

#include <vector>
using std::vector;

#include "lambda_lanczos/lambda_lanczos.hpp"
//using lambda_lanczos::LambdaLanczos;


namespace lambda_lanczos {

/// @brief PEigenCalculator is a class containing only one useful member
/// function PrincipalEigen().  This function calculates the principal (largest
/// or smallest) eigenvalue and corresponding eigenvector of an n x n matrix.

template<typename Scalar>
class PEigenCalculator
{
  size_t n;            // the size of the matrix
  vector<Scalar> evec; // preallocated vector (lambda_lanzcos does not use ptrs)

public:

  PEigenCalculator(int matrix_size):evec(matrix_size) //!< must specify matrix size in advance
  {
    n = matrix_size;
    assert(n > 0);
  }

  /// @brief  Calculate the principal eigenvalue and eigenvector of a matrix.
  /// @return Return the principal eigenvalue of the matrix.
  ///         If you want the eigenvector, pass a non-null "evector" argument.
  Scalar
  PrincipalEigen(Scalar const* const *matrix, //!< the input patrix
                 Scalar *evector=nullptr, //!< optional: eigenvector stored here
                 bool find_max=false);    //!< want the max or min eigenvalue?

}; // class PEigenCalculator


template<typename Scalar>
Scalar PEigenCalculator<Scalar>::
  PrincipalEigen(Scalar const* const *matrix,
                 Scalar *eigenvector,
                 bool find_max)
{
  auto matmul = [&](const vector<double>& in, vector<double>& out) {
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
  LambdaLanczos<Scalar> ll_engine(matmul, n, find_max);
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


} // namespace lambda_lanczos {

#endif //#ifndef _PEIGENCALC_HPP
