#ifndef _PEIGENCALC_HPP
#define _PEIGENCALC_HPP

#include <vector>
using std::vector;

#include "lambda_lanczos/lambda_lanczos.hpp"
using lambda_lanczos::LambdaLanczos;


namespace superpose3d_lammps {


/// @brief PEigenCalculator caluclates the principal (largest)
/// eigenvalue and corresponding eigenvector of an n x n matrix.
/// Right now it is just a wrapper enclosing "lambda-lanczos".
/// (That might change of other developers want to swap it with 
///  the "Eigen" library or something else.)

template<typename Scalar>
class PEigenCalculator {
  size_t n;                  // the size of the matrix
  vector<vector<Scalar> > M; // the matrix
  vector<Scalar> evec;       // preallocated vector (needed by lambda_lanzcos)
  LambdaLanczos<Scalar> ll_engine;  // this is the object that does the work

  void Init();
  void SetSize(int n);
  void SetFindMax(bool findmax);

public:

  PEigenCalculator(int n=0, bool findmax = false):ll_engine(), evec()
  {
    Init();
    SetSize(n);
    SetFindMax(findmax);
  }
  
  /// @brief  Calculate the principal eigenvalue and eigenvector of a matrix.
  /// @return Return the principal eigenvalue of the matrix.
  ///         If you want the eigenvector, pass a non-null "evect" argument.
  Scalar
  PrincipalEigen(Scalar const* const *matrix,  //!< the input patrix
                 Scalar *evect=nullptr   //!< optional: store eigenvector
                 );

}; // class PEigenCalculator


// --- PEigenCalculator IMPLEMENTATION ---
template<typename Scalar>
void PEigenCalculator<Scalar>::Init() {
  n = 0;
  auto matmul = [&](const vector<double>& in, vector<double>& out) {
    for(int i = 0;i < n;i++) {
      for(int j = 0;j < n;j++) {
        out[i] += M[i][j]*in[j];
      }
    } 
  };
  ll_engine.SetMul(matmul);
  auto init_vec = [&](vector<double>& vec) {
    for(int i = 0;i < n;i++) {
      vec[i] = 0.0;
    }
    vec[0] = 1.0;
  };
  ll_engine.SetInitVec(init_vec);
}

template<typename Scalar>
void PEigenCalculator<Scalar>::SetSize(int n) {
  if (this->n != n) {
    this->n = n;
    M.resize(n);
    for (int i = 0; i < n; i++)
      M[i].resize(n);
    evec.resize(n);
  }
  ll_engine.SetSize(n);
}

template<typename Scalar>
void PEigenCalculator<Scalar>::SetFindMax(bool findmax) {
  ll_engine.SetFindMax(findmax);
}

template<typename Scalar>
Scalar PEigenCalculator<Scalar>::
  PrincipalEigen(Scalar const* const *matrix,  //!< the input patrix
                 Scalar *eigenvector)          //!< optional: store eigenvector
{
  assert(n > 0);

  Scalar eval;

  // We must copy the data from matrix into M.
  // (Because "matmul" refers to M.)
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      M[i][j] = matrix[i][j];

  size_t itern = ll_engine.run(eval, evec);

  if (eigenvector) {
    // If the caller requested the eigenvector as well, then
    // return it to the caller by copying the data into eigenvector[].
    for (int i = 0; i < n; i++)
      eigenvector[i] = evec[i];
  }
 
  return eval;
}

}  // namespace superpose3d_lammps {


#endif //#ifndef _PEIGENCALC_HPP
