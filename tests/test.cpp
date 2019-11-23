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

  Superpose3D superposer(N);
  superposer.Superpose(X, x);

  // Print the optimal rotation and translations

  cout << "Optimal superposition rmsd = " << superposer.rmsd << "\n";
  cout << "optimal translation, Ti =\n"<<superposer.T[0]<<" "<<superposer.T[1]<<" "<<superposer.T[2]<<"\n";
  cout << "Optimal Rotation, Rij = \n";
  for (int iy = 0; iy < 3; iy++) {
    for (int ix = 0; ix < 3; ix++)
      cout << superpose.R[iy][ix] << " ";
    cout << "\n";
  }

  // Since one point cloud is just a rotate version of the other, we expect
  // that the RMSD between them should be 0.  Insure that this is so.
  assert(abs(superposer.rmsd) < 1.0e-6);

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
  result = Superpose3D(X, Xscshift, nullptr, true);
  
  cout << "Optimal Scale factor, C = " << superposer.C << "\n";
  cout << "Optimal superposition rmsd = " << superposer.rmsd << "\n";
  cout << "optimal translation, Ti =\n"<<superposer.T[0]<<" "<<superposer.T[1]<<" "<<superposer.T[2]<<"\n";
  cout << "Optimal Rotation, Rij = \n";
  for (int iy = 0; iy < 3; iy++) {
    for (int ix = 0; ix < 3; ix++)
      cout << superpose.R[iy][ix] << " ";
    cout << "\n";
  }

} // main()
