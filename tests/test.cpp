#include <alloc2d.hpp>
#include <superpose3d.hpp>

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
  double _Xscaled[N][3]; // allocate space
  double _Xshifted[N][3];
  double _Xscshift[N][3];
  double *Xscaled[N]= {_Xcaled[0], _Xscaled[1], _Xscaled[2], _Xscaled[3]};
  double *Xshifted[N]= {_Xshifted[0], _Xshifted[1], _Xshifted[2], _Xshifted[3]};
  double *Xscshift[N]= {_Xscshift[0], _Xscshift[1], _Xscshift[2], _Xscshift[3]};
  for (int i = 0; i < N; i++) {
    for (int d = 0; d < 3; d++ {
      Xscaled[i][d] = 2.0 * X[i][d];  //coords are scaled
      Xshifted[i][d] = X[i][d];
      if (d==0) Xshifted += 100.0;    //coords are shifted in the x direction
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

  CONTINUEHERE: the following code needs to be converted from python to c++

  # Does the RMSD returned in result[0] match the RMSD calculated manually?
  R = np.matrix(result[1])              # rotation matrix
  T = np.matrix(result[2]).transpose()  # translation vector (3x1 matrix)
  c = result[3]                         # scalar
  _x = np.matrix(xscshift).transpose()
  _xprime = c*R*_x + T
  xprime = np.array(_xprime.transpose()) # convert to length 3 numpy array
  RMSD = 0.0

  for i in range(0, len(X)):
      RMSD += ((X[i][0] - xprime[i][0])**2 +
               (X[i][1] - xprime[i][1])**2 +
               (X[i][2] - xprime[i][2])**2)

  assert(abs(RMSD - result[0]) < 1.0e-6)

} // main()
