#include <iostream>
#include <cmath>
#include <iomanip>
//#include <random> 
#include <vector>
#include "meep.hpp"
#include "newton_raphson.hpp"
#include <unistd.h>


/*
Algorithm for solving set of three coupled nonlinear quadratic forms for
deriving E from D field in cases where full Chi2 tensor is required:
Ex' = Bx + Cx^2 + Dy^2 + Ez^2 + Fyz + Gxz + Hxy,
Ey' = Jy + Kx^2 + Ly^2 + Mz^2 + Nyz + Oxz + Pxy,
Ez' = Rz + Sx^2 + Ty^2 + Uz^2 + Vyz + Wxz + Xxy,

where En' = Dn - Pn_dispersive (as implemented in update_eh.cpp, and meep's lorentzian material
model in which Pn_dispersive are the polarisations due to 1st order material dispersion (i.e. not
including nonlinear 'chi type' polarisations)
*/

using namespace std;


namespace meep {



const double TOLERANCE = 1e-10;
const int MAX_ITERATIONS = 500;

//// Struct to hold constant parameters for each equation
//struct Parameters { // moved to .hpp??!
//  double A, B, C, D, E, F, G, H;
//};

// Function definitions
vector<double> equations(realnum x, realnum y, realnum z, const Parameters &p1, const Parameters &p2,
                         const Parameters &p3);
vector<vector<double> > jacobian(realnum x, realnum y, realnum z, const Parameters &p1,
                                 const Parameters &p2, const Parameters &p3);
vector<double> solveLinearSystem(const vector<vector<double> > &J, const vector<double> &F);

// Compute the coefficient matrix (ignoring constant terms)
vector<vector<double> > computeCoefficientMatrix(const Parameters &p1, const Parameters &p2,
                                                 const Parameters &p3) {
  vector<vector<double> > M(3, vector<double>(3));
  M[0][0] = -p1.B;
  M[0][1] = -p1.H;
  M[0][2] = -p1.G;
  M[1][0] = -p2.B;
  M[1][1] = -p2.H;
  M[1][2] = -p2.G;
  M[2][0] = -p3.B;
  M[2][1] = -p3.H;
  M[2][2] = -p3.G;
  return M;
}

// Compute determinant of a 3x3 matrix
double determinant3x3(const vector<vector<double> > &M) {
  return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) -
         M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
         M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
}

// Compute rank of a 3x3 matrix (naive approach)
int rank3x3(const vector<vector<double> > &M) {
  double detM = determinant3x3(M);
  if (detM != 0) return 3;

  // Check if any 2x2 minor is nonzero (indicating rank 2)
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      vector<vector<double> > minorM;
      for (int m = 0; m < 3; m++) {
        if (m == i) continue;
        vector<double> row;
        for (int n = 0; n < 3; n++) {
          if (n == j) continue;
          row.push_back(M[m][n]);
        }
        minorM.push_back(row);
      }
      if (minorM[0][0] * minorM[1][1] - minorM[0][1] * minorM[1][0] != 0) { return 2; }
    }
  }
  return 1;
}

bool newtonRaphson(realnum x, realnum y, realnum z, const Parameters &p1,
                             const Parameters &p2, const Parameters &p3,  realnum* fw,
                             realnum* fw_2,
                             realnum* fw_3) {
  for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
    vector<double> F = equations(x, y, z, p1, p2, p3);
    vector<vector<double> > J = jacobian(x, y, z, p1, p2, p3);

    // CHECK. The next calc of detJ and check are only necessary for checking purposes, not the NR
    // alg.
    double detJ = determinant3x3(J);
    if (detJ == 0) {
      cout << "Jacobian determinant is 0: The system has dependent or singular solutions at (" << x
           << ", " << y << ", " << z << ")!" << endl;
    }
    // END CHECK

    vector<double> delta = solveLinearSystem(J, F);

    x -= delta[0];
    y -= delta[1];
    z -= delta[2];

    if (fabs(delta[0]) < TOLERANCE && fabs(delta[1]) < TOLERANCE && fabs(delta[2]) < TOLERANCE) {
   //   cout << "Converged after " << iter + 1 << " iterations:   "; //.\n";
   //   cout << "x = " << x << ", y = " << y << ", z = " << z << "\n";
      *fw = x; /// Update E field values!
   *fw_2 = y;
      *fw_3 = z;
      return true;
    }
  }
  
  cout << "Newton's method did not converge within " << MAX_ITERATIONS << " iterations.\n";
  return false;
}

vector<double> equations(double x, double y, double z, const Parameters &p1, const Parameters &p2,
                         const Parameters &p3) {
  return {p1.A - (p1.B * x + p1.C * x * x + p1.D * y * y + p1.E * z * z + p1.F * y * z +
                  p1.G * x * z + p1.H * x * y),
          p2.A - (p2.B * y + p2.C * x * x + p2.D * y * y + p2.E * z * z + p2.F * y * z +
                  p2.G * x * z + p2.H * x * y),
          p3.A - (p3.B * z + p3.C * x * x + p3.D * y * y + p3.E * z * z + p3.F * y * z +
                  p3.G * x * z + p3.H * x * y)};
}

vector<vector<double> > jacobian(double x, double y, double z, const Parameters &p1,
                                 const Parameters &p2, const Parameters &p3) {
  return {{-p1.B - 2 * p1.C * x - p1.G * z - p1.H * y, -2 * p1.D * y - p1.F * z - p1.H * x,
           -2 * p1.E * z - p1.F * y - p1.G * x},
          {-2 * p2.C * x - p2.G * z - p2.H * y, -p2.B - 2 * p2.D * y - p2.F * z - p2.H * x,
           -2 * p2.E * z - p2.F * y - p2.G * x},
          {-2 * p3.C * x - p3.G * z - p3.H * y, -2 * p3.D * y - p3.F * z - p3.H * x,
           -p3.B - 2 * p3.E * z - p3.F * y - p3.G * x}};
}

vector<double> solveLinearSystem(const vector<vector<double> > &J, const vector<double> &F) {
  vector<vector<double> > A = J;
  vector<double> b = F;
  int n = 3;

  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      double factor = A[j][i] / A[i][i];
      for (int k = i; k < n; k++) {
        A[j][k] -= factor * A[i][k];
      }
      b[j] -= factor * b[i];
    }
  }

  vector<double> x(n);
  for (int i = n - 1; i >= 0; i--) {
    x[i] = b[i];
    for (int j = i + 1; j < n; j++) {
      x[i] -= A[i][j] * x[j];
    }
    x[i] /= A[i][i];
  }
  return x;
}


void runNR(realnum seed1, realnum seed2, realnum seed3, realnum* fw, realnum* fw_2, realnum* fw_3, const Parameters &p1, const Parameters &p2,
           const Parameters &p3) { // TODO need to confirm that passing fw through as a ref like this actually works...
        
  //  cout << "Doing NR" << endl;
       
      // CHECK 2:
      vector<vector<double> > M = computeCoefficientMatrix(p1, p2, p3);
      int rankM = rank3x3(M);
      if (rankM < 3) {
        cout << "Coefficient matrix has rank < 3: The system is globally dependent!" << endl;
        cout << " s1: " << seed1 << " s2: " << seed2 << " s3: " << seed3 << " f1: " << *fw << " fw2: "
             << *fw_2 << "fw3: " << *fw_3 << endl;
        cout << " p1A: " << p1.A << " p1B: " << p1.B << " p1F: " << p1.F << endl;
        cout << " p2A: " << p2.A << " p2B: " << p2.B << " p2F: " << p2.G << endl;
        cout << " p3A: " << p3.A << " p3B: " << p3.B << " p3F: " << p3.H << endl;
        sleep(1);

      return;
      }
      // END CHECK 2

      //Parameters p1 = {a, 3.2, 0, 0, 0, 0.000002, 0, 0}; // p1.A is the D-P field p1.B is epsilon, then chi2
      //Parameters p2 = {b, 3.2, 0, 0, 0, 0, 0.000002, 0};
      //Parameters p3 = {c, 3.2, 0, 0, 0, 0, 0, 0.000002};
      bool counter = false;
      for (int i = 2, imax = 999; i < imax; i *= i) {
       if (newtonRaphson(seed1+3.06*i, seed2-2.43*i, seed3+1.277*i, p1, p2, p3, fw, fw_2, fw_3)) {
     //   if (newtonRaphson(seed1, seed2, seed3, p1, p2, p3, fw, fw_2, fw_3)) {
         if (counter) { cout << "True " << i << endl; 
         sleep(19);
         }
          counter = false;  
          break;
        }
        else {
          cout << "NR didn't converge: " << i << endl;
          cout << " s1: " << seed1 << " s2: " << seed2 << " s3: " << seed3 << " f1: " << *fw
               << " fw2: " << *fw_2 << "fw3: " << *fw_3 << endl;
          cout << " p1A: " << p1.A << " p1B: " << p1.B << " p1F: " << p1.F << endl;
          cout << " p2A: " << p2.A << " p2B: " << p2.B << " p2F: " << p2.G << endl;
          cout << " p3A: " << p3.A << " p3B: " << p3.B << " p3F: " << p3.H << endl;
          
          counter = true;
        }
      }
      if (counter) {
        cout << "FALSE "  << endl;
        sleep(19);
      }
  }

}


