#include <iostream>
#include <cmath>
#include <chrono>
#include <random>
using namespace std;

#include "superpose3d.hpp"
using namespace superpose3d;


int main(int argc, char **argv) {
  double **X; //the immobile point cloud
  double **x; //the mobile point cloud
  double **xprime; //the mobile point cloud after transformations were applied
  double *w;  //weights to apply to each point when calculating RMSD

  int n_points = 4;
  if (argc > 1)
    n_points = atoi(argv[1]);
  int n_point_clouds = 1;
  if (argc > 2)
    n_point_clouds = atoi(argv[2]);

  // Superpose3D<double,double const*const*,double const*> s(n_points,weights);
  Superpose3D<double, double const* const*> superposer(n_points);
  // Now test the copy constructor and = operator
  Superpose3D<double, double const* const*> superposer_cpy(superposer);
  Superpose3D<double, double const* const*> su = superposer_cpy;
  // Test SetNumPoints()
  su.SetNumPoints(n_points);

  w = new double [n_points];
  for (int i=0; i < n_points; i++)
    w[i] = 1.0 / n_points;

  // Allocate the immobile point cloud array (X) and fill it with coordinates
  Alloc2D(n_points, 3, &X);
  // Allocate the mobile point cloud array (x) and fill it with coordinates
  Alloc2D(n_points, 3, &x);
  // Allocate space for the transformed mobile point cloud
  Alloc2D(n_points, 3, &xprime);

  int seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine rand_generator(seed);
  std::normal_distribution<double> gaussian_distribution(0, 1.0);

  for (int i_cl = 0; i_cl < n_point_clouds; i_cl++) {

    for (int i = 0; i < n_points; i++)
      for (int d = 0; d < 3; d++)
        X[i][d] = gaussian_distribution(rand_generator);
    for (int i = 0; i < n_points; i++)
      for (int d = 0; d < 3; d++)
        x[i][d] = gaussian_distribution(rand_generator);

    // Now generate some random weights
    double w_norm2 = 0.0;
    for (int i = 0; i < n_points; i++) {
      w[i] = gaussian_distribution(rand_generator);
      w_norm2 += w[i]*w[i];
    }
    for (int i = 0; i < n_points; i++)
      w[i] = w[i]*w[i] / w_norm2;
    su.SetWeights(w); // inform su that we want to use these weights
  
    // ----------------- main code -----------------

    // Now superimpose the two point clouds:
    double rmsd;

    bool allow_rescale = false;
    if (i_cl % 2 == 0) // (<-- test both with and without allowing rescale)
      allow_rescale = true;

    rmsd = su.Superpose(X, x, allow_rescale);

    if (n_point_clouds == 1) {
      cout << "X (frozen) = \n";
      for (size_t n = 0; n < n_points; n++) {
        for (int d = 0; d < 3; d++) {
          cout << X[n][d];
          if (d+1 < 3)
            cout << " ";
          else
            cout << "\n";
        }
      }
      cout << "x (mobile) = \n";
      for (size_t n = 0; n < n_points; n++) {
        for (int d = 0; d < 3; d++) {
          cout << x[n][d];
          if (d+1 < 3)
            cout << " ";
          else
            cout << "\n";
        }
      }
      cout << "Optimal Scale factor, C = " << su.c << "\n";
      cout << "Optimal superposition rmsd = " << rmsd << "\n";
      cout << "optimal translation, Ti =\n"
           << su.T[0] << " " << su.T[1] << " " << su.T[2] << "\n";
      cout << "Optimal Rotation, Rij = \n";
      for (int iy = 0; iy < 3; iy++) {
        for (int ix = 0; ix < 3; ix++)
          cout << su.R[iy][ix] << " ";
        cout << "\n";
      }
    }

    // Now apply this transformation to the mobile point cloud, save the
    // new coordinates in "xprime".  Then calculate the RMSD by summing the
    // (squared) distances between points in "X" and points in "xprime".
    // Then compare this RMSD with the rmsd calculated by Superpose().
    // They should agree.  First store the new coordinates in "xprime":
    for (size_t n = 0; n < n_points; n++) {
      for (int i = 0; i < 3; i++) {
        xprime[n][i] = 0.0;
        for (int j = 0; j < 3; j++)
          xprime[n][i] += su.c * su.R[i][j] * x[n][j];
        xprime[n][i] += su.T[i];
      }
    }
    // Now calculate the RMSD beween X and xprime
    double RMSD = 0.0;
    double sum_w = 0.0;
    for (size_t n = 0; n < n_points; n++)
      sum_w += w[n];
    for (size_t n = 0; n < n_points; n++)
      RMSD += w[n] * (SQR(X[n][0] - xprime[n][0]) +
                      SQR(X[n][1] - xprime[n][1]) +
                      SQR(X[n][2] - xprime[n][2]));
    RMSD = sqrt(RMSD / sum_w);

    // Compare this (direct) method for computing RMSD with the prediction
    // from Superpose3d::Superpose() (which is currently stored in "rmsd")
    assert(abs(RMSD - rmsd) < 1.0e-6);

  } // for (i_cl = 0; i_cl < n_point_clouds; i_cl++)

  Dealloc2D(&xprime);
  Dealloc2D(&X);
  Dealloc2D(&x);
  delete [] w;

  cout << "\n" << "test passed." << endl;
  return EXIT_SUCCESS;
} // main()
