#include "../src/superpose3d.hpp"

void main(int argc, char **argv) {
  const int N = 4;

  // Allocate the arrays for the coordinates for this quick and dirty example.
  // (To make it easier to initialize them, use fixed-size arrays.)
  double _X[N][3] = {{0,0,0},{0,0,1},{0,1,1},{1,1,1}};
  double _x[N][3] = {{0,0,0},{0,1,0},{0,1,-1},{1,1,-1}};
  // Must convert these 2D arrays into to pointer-to-pointer-to-double
  double *X[N] = {_X[0], _X[1], _X[2], _X[3]};
  double *x[N] = {_x[0], _x[1], _x[2], _x[3]};
  // (Note: x=X after rotation around the Z axis)

  // Allocate space for X and x, and load the coordinates (omitted)

  Superpose3D s(N);
  s.Superpose(X, x);

  // Print the optimal rotation and translations

  cout << "Optimal superposition rmsd = " << s.rmsd << "\n";
  cout << "optimal translation, Ti =\n"<<s.T[0]<<" "<<s.T[1]<<" "<<s.T[2]<<"\n";
  cout << "Optimal Rotation, Rij = \n";
  for (int iy = 0; iy < 3; iy++) {
    for (int ix = 0; ix < 3; ix++)
      cout << s.R[iy][ix] << " ";
    cout << "\n";
  }

  // Since one point cloud is just a rotate version of the other, we expect
  // that the RMSD between them should be 0.  Insure that this is so.
  assert(abs(s.rmsd) < 1.0e-6);

  // Now create some versions of "X" that have been modified in some way
  // and try again:
  //double _Xscaled[N][3]; // allocate space
  //double _Xshifted[N][3];
  double _Xscshift[N][3];
  //double *Xscaled[N]= {_Xcaled[0], _Xscaled[1], _Xscaled[2], _Xscaled[3]};
  //double *Xshifted[N]={_Xshifted[0], _Xshifted[1], _Xshifted[2],_Xshifted[3]};
  double *Xscshift[N]= {_Xscshift[0], _Xscshift[1], _Xscshift[2], _Xscshift[3]};
  for (int i = 0; i < N; i++) {
    for (int d = 0; d < 3; d++ {
      //Xscaled[i][d]  = 2.0 * X[i][d];  //coords are scaled
      // Xshifted[i][d] = X[i][d];
      //if (d==0) Xshifted += 100.0;     //coords are shifted in the x direction
      Xscshift[i][d] = 2.0 * X[i][d];
      if (d==1) Xscshift[i][d] += 200; //coords are scaled and shifted along y
    }
  }

  // Now try superposition again using these new coordinates
  result = s.Superpose(X, Xscshift, nullptr, true);
  
  cout << "Optimal Scale factor, C = " << s.C << "\n";
  cout << "Optimal superposition rmsd = " << s.rmsd << "\n";
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
  double **aaXprime;
  double *aXprime;
  Alloc2D(N, 3, &aXprime, &aaXprime);
  for (size_t n = 0; n < N; n++) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++)
        cout << s.R[i][j] << " ";
      cout << "\n";
    }
  }

  // Now apply this transformation to the mobile point cloud. Save in "aaXprime"
  for (size_t n = 0; n < N; n++) {
    for (int i = 0; i < 3; i++) {
      aaXprime[n][i] = 0.0;
      for (int j = 0; j < 3; j++)
        aaXprime[n][i] += s.c * s.R[i][j] * xscshift[j];
      aaXprime[n][i] += s.T[i];
    }
  }

  // Calculate the sum-squared distance between points in the original
  // point cloud, and the mobile point cloud (after transformation applied).
  double RMSD = 0.0;

  for (size_t n = 0; n < N; n++)
    RMSD += ((aaXf[i][0] - aaXprime[i][0])**2 +
             (aaXf[i][1] - aaXprime[i][1])**2 +
             (aaXf[i][2] - aaXprime[i][2])**2);

  RMSD = sqrt(RMSD / N);

  assert(abs(RMSD - result[0]) < 1.0e-6);

  Dealloc2D(&aXprime, &aaXprime);

} // main()
