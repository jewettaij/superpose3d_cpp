#include <iostream>
#include <cmath>
using namespace std;

#include "superpose3d.hpp"
using namespace superpose3d;


int main(int argc, char **argv) {
  double *_X; //the immobile point cloud (in continguous memory)
  double **X; //the immobile point cloud (as a pointer->pointer)
  double *_x; //the mobile point cloud (in continguous memory)
  double **x; //the mobile point cloud (as a pointer->pointer)

  const int N = 4;

  // Allocate the immobile point cloud array (X) and fill it with coordinates
  Alloc2D(N, 3, &_X, &X);
  X[0][0] = 0.0;
  X[0][1] = 0.0;
  X[0][2] = 0.0;

  X[1][0] = 0.0;
  X[1][1] = 0.0;
  X[1][2] = 1.0;

  X[2][0] = 0.0;
  X[2][1] = 1.0;
  X[2][2] = 1.0;

  X[3][0] = 1.0;
  X[3][1] = 1.0;
  X[3][2] = 1.0;

  // Allocate the mobile point cloud array (x) and fill it with coordinates
  Alloc2D(N, 3, &_x, &x);
  x[0][0] = 0.0;
  x[0][1] = 0.0;
  x[0][2] = 0.0;

  x[1][0] = 0.0;
  x[1][1] = 1.0;
  x[1][2] = 0.0;

  x[2][0] = 0.0;
  x[2][1] = 1.0;
  x[2][2] = -1.0;

  x[3][0] = 1.0;
  x[3][1] = 1.0;
  x[3][2] = -1.0;

  // ----------------- main code -----------------
  // Note: Superpose only works on 2D arrays which are implemented as C-style
  // pointer->pointer arrays.  (vector<vector>> or fixed size arrays like
  // "double _X[4][3]" will not work and cannot be sent to
  // Superpose3D::Superpose())

  // Now superimpose the two point clouds:
  double rmsd;
  Superpose3D<double> s(N);

  rmsd = s.Superpose(X, x);

  // Print the optimal rotation and translations

  cout << "Optimal superposition rmsd = " << rmsd << "\n";
  cout << "optimal translation, Ti =\n"<<s.T[0]<<" "<<s.T[1]<<" "<<s.T[2]<<"\n";
  cout << "Optimal Rotation, Rij = \n";
  for (int iy = 0; iy < 3; iy++) {
    for (int ix = 0; ix < 3; ix++)
      cout << s.R[iy][ix] << " ";
    cout << "\n";
  }
  // -------------------------------------------

  // Since one point cloud is just a rotated version of the other, we expect
  // that the RMSD between them should be 0.  Insure that this is so.
  assert(std::abs(rmsd) < 1.0e-6);

  // Now create some versions of "X" that have been modified in some way
  // and try again:
  double *_Xscshift;
  double **Xscshift;
  Alloc2D(N, 3, &_Xscshift, &Xscshift);
  for (int i = 0; i < N; i++) {
    for (int d = 0; d < 3; d++) {
      Xscshift[i][d] = 2.0 * x[i][d];
      if (d==1) Xscshift[i][d] += 200; //coords are scaled and shifted along y
    }
  }

  // Now try superposition again using these new coordinates
  rmsd = s.Superpose(X, Xscshift, nullptr, true);
  
  cout << "Optimal Scale factor, C = " << s.c << "\n";
  cout << "Optimal superposition rmsd = " << rmsd << "\n";
  cout << "optimal translation, Ti =\n"<<s.T[0]<<" "<<s.T[1]<<" "<<s.T[2]<<"\n";
  cout << "Optimal Rotation, Rij = \n";
  for (int iy = 0; iy < 3; iy++) {
    for (int ix = 0; ix < 3; ix++)
      cout << s.R[iy][ix] << " ";
    cout << "\n";
  }

  // Now apply this transformation to the new coordinates
  // and recalculate the RMSD manually (by comparing coordinate sets).
  // Check to see if the RMSD computed by two different methods agrees.
  double **Xprime;
  double *_aXprime;
  Alloc2D(N, 3, &_aXprime, &Xprime);

  // Now apply this transformation to the mobile point cloud. Save in "aaXprime"
  for (size_t n = 0; n < N; n++) {
    for (int i = 0; i < 3; i++) {
      Xprime[n][i] = 0.0;
      for (int j = 0; j < 3; j++)
        Xprime[n][i] += s.c * s.R[i][j] * Xscshift[n][j];
      Xprime[n][i] += s.T[i];
    }
  }

  // Manually calculate the sum-squared distance between points in the original
  // point cloud, and the mobile point cloud (after transformation applied).
  double RMSD = 0.0;

  for (size_t n = 0; n < N; n++)
    RMSD += (SQR(X[n][0] - Xprime[n][0]) +
             SQR(X[n][1] - Xprime[n][1]) +
             SQR(X[n][2] - Xprime[n][2]));

  RMSD = sqrt(RMSD / N);

  // Compare the manual method of computing RMSD with the
  // prediction from Superpose3d::Superpose() (currently stored in "rmsd")
  assert(abs(RMSD - rmsd) < 1.0e-6);

  Dealloc2D(&_Xscshift, &Xscshift);
  Dealloc2D(&_aXprime, &Xprime);
  Dealloc2D(&_X, &X);
  Dealloc2D(&_x, &x);

  return EXIT_SUCCESS;
} // main()
