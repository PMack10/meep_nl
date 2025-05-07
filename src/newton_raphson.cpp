#include <iostream>
#include <cmath>
#include <iomanip>
#include <random> 
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



const double TOLERANCE = 1e-8;
int MAX_ITERATIONS = 500;
const double FIELDCHECKPERCENT = 1e-4;

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
                             realnum* fw_2, realnum *fw_3, double tol1,
                   double tol2, double tol3) {
  for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
    vector<double> F = equations(x, y, z, p1, p2, p3);
    vector<vector<double> > J = jacobian(x, y, z, p1, p2, p3);

    // CHECK. The next calc of detJ and check are only necessary for checking purposes, not the NR
    // alg.
    // double detJ = determinant3x3(J);
    // if (detJ == 0) {
    //  cout << "Jacobian determinant is 0: The system has dependent or singular solutions at (" <<
    //  x
    //       << ", " << y << ", " << z << ")!" << endl;
    //}
    // END CHECK

    vector<double> delta = solveLinearSystem(J, F);

    x -= delta[0];
    y -= delta[1];
    z -= delta[2];

    // if (fabs(delta[0]) < TOLERANCE && fabs(delta[1]) < TOLERANCE && fabs(delta[2]) < TOLERANCE) {
    if (fabs(delta[0]) < tol1 && fabs(delta[1]) < tol2 && fabs(delta[2]) < tol3) {
      //   cout << "Converged after " << iter + 1 << " iterations:   "; //.\n";
      //   cout << "x = " << x << ", y = " << y << ", z = " << z << "\n";
    //  double fieldCheckPercent = fmax(fabs(x), fmax(fabs(y), fabs(z))) * FIELDCHECKPERCENT;
      double ax = std::abs(x), ay = std::abs(y), az = std::abs(z);
      double fieldCheckPercent = (ax > ay ? (ax > az ? ax : az) : (ay > az ? ay : az)) * FIELDCHECKPERCENT;
      vector<double> fCheck = equations(
          x, y, z, p1, p2, p3); // call eqns here with the final values to verify... (because
                                // current stopping condition is just NR algo step delta)
      if (fCheck[0] <= fieldCheckPercent && fCheck[1] <= fieldCheckPercent && fCheck[2] <= fieldCheckPercent) {
        *fw = x; /// Update E field values!
        *fw_2 = y;
        *fw_3 = z;
        return true;
      }
      else {
        cout << "Failed field tolerance check" << endl;
        cout << "x, y, z: " << x << " " << y << " " << z << endl;
       // sleep(10);
      }
    }
  }
 // cout << "Newton's method did not converge within " << MAX_ITERATIONS << " iterations.\n";
  return false;
}

vector<double> equations(double x, double y, double z, const Parameters &p1, const Parameters &p2,
                         const Parameters &p3) {
  return {p1.A - (p1.B * x + p1.F * y * z + p1.G * x * z + p1.H * x * y),
          p2.A - (p2.B * y + p2.F * y * z + p2.G * x * z + p2.H * x * y),
          p3.A - (p3.B * z + p3.F * y * z + p3.G * x * z + p3.H * x * y)}; //TODO, deleted the n^2 terms for efficiency (as not relevant for 43m)
  //return {p1.A - (p1.B * x + p1.C * x * x + p1.D * y * y + p1.E * z * z + p1.F * y * z +
  //                p1.G * x * z + p1.H * x * y),
  //        p2.A - (p2.B * y + p2.C * x * x + p2.D * y * y + p2.E * z * z + p2.F * y * z +
  //                p2.G * x * z + p2.H * x * y),
  //        p3.A - (p3.B * z + p3.C * x * x + p3.D * y * y + p3.E * z * z + p3.F * y * z +
  //                p3.G * x * z + p3.H * x * y)};
}

vector<vector<double> > jacobian(double x, double y, double z, const Parameters &p1,
                                 const Parameters &p2, const Parameters &p3) {
  return {{-p1.B - p1.G * z - p1.H * y,  - p1.F * z - p1.H * x, - p1.F * y - p1.G * x},
          {- p2.G * z - p2.H * y, -p2.B - p2.F * z - p2.H * x, - p2.F * y - p2.G * x},
          {- p3.G * z - p3.H * y, - p3.F * z - p3.H * x, -p3.B - p3.F * y - p3.G * x}};//TODO, deleted the n^2 terms for efficiency (as not relevant for 43m)
  //return {{-p1.B - 2 * p1.C * x - p1.G * z - p1.H * y, -2 * p1.D * y - p1.F * z - p1.H * x,
  //         -2 * p1.E * z - p1.F * y - p1.G * x},
  //        {-2 * p2.C * x - p2.G * z - p2.H * y, -p2.B - 2 * p2.D * y - p2.F * z - p2.H * x,
  //         -2 * p2.E * z - p2.F * y - p2.G * x},
  //        {-2 * p3.C * x - p3.G * z - p3.H * y, -2 * p3.D * y - p3.F * z - p3.H * x,
  //         -p3.B - 2 * p3.E * z - p3.F * y - p3.G * x}};
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

//static uniform_real_distribution<double> dist(-1, 1); // Define range
static uniform_real_distribution<double> distSign(-1.0, 1.0); // Define range
std::lognormal_distribution<double> dist(1.0, 90);
static random_device rd;                              // Obtain a random seed
static double seedMax = 1e33; //TODO reduced from e125

double getRandomNumber() {
  static mt19937 gen(rd()); // Seed the generator 
  return dist(gen)*distSign(gen);// *1e80;
 // return dist(gen) * 1e80;
}


void runNR(realnum seed1, realnum seed2, realnum seed3, realnum* fw, realnum* fw_2, realnum* fw_3, const Parameters &p1, const Parameters &p2,
           const Parameters &p3) { // TODO need to confirm that passing fw through as a ref like this actually works...
  std::cout << std::fixed << std::setprecision(12);
  //  cout << "Doing NR" << endl;
      
      double fwxInitial = *fw;
      MAX_ITERATIONS = 250;
      double tol1 = fmax(fabs(TOLERANCE * (*fw))*0.0001, TOLERANCE); // TODO now added equation percentage accuracy check so may be able to ditch the fmax here and just use a fixed tolerance
      double tol2 = fmax(fabs(TOLERANCE * (*fw_2))*0.0001, TOLERANCE);
      double tol3 = fmax(fabs(TOLERANCE * (*fw_3)) * 0.0001, TOLERANCE);
      double _seed1 = seed1;
      double _seed2 = seed2;
      double _seed3 = seed3;

      // CHECK 2:
      //vector<vector<double> > M = computeCoefficientMatrix(p1, p2, p3);
      //int rankM = rank3x3(M);
      //if (rankM < 3) {
      //  cout << "Coefficient matrix has rank < 3: The system is globally dependent!" << endl;
      //  cout << " s1: " << seed1 << " s2: " << seed2 << " s3: " << seed3 << " f1: " << *fw << " fw2: "
      //       << *fw_2 << " fw3: " << *fw_3 << endl;
      //  cout << " p1A: " << p1.A << " p1B: " << p1.B << " p1F: " << p1.F << endl;
      //  cout << " p2A: " << p2.A << " p2B: " << p2.B << " p2G: " << p2.G << endl;
      //  cout << " p3A: " << p3.A << " p3B: " << p3.B << " p3H: " << p3.H << endl;
      //  sleep(4);

      //return;
      //}
      // END CHECK 2

      //Parameters p1 = {a, 3.2, 0, 0, 0, 0.000002, 0, 0}; // p1.A is the D-P field p1.B is epsilon, then chi2
      //Parameters p2 = {b, 3.2, 0, 0, 0, 0, 0.000002, 0};
      //Parameters p3 = {c, 3.2, 0, 0, 0, 0, 0, 0.000002};
      bool counter = false;
      for (int i = 0, imax = 100; i < imax; ++i) {
        if (newtonRaphson(_seed1, _seed2, _seed3, p1, p2, p3, fw, fw_2, fw_3, tol1, tol2, tol3)) {
          //   if (newtonRaphson(seed1, seed2, seed3, p1, p2, p3, fw, fw_2, fw_3)) {
       cout << "Converged " << i << endl;
        /*     cout << "tols: " << tol1 << "  " << tol2 << "  " << tol3 << endl;
          cout << " s1: " << _seed1 << " s2: " << _seed2 << " s3: " << _seed3 << " f1: " << *fw
               << " fw2: " << *fw_2 << "fw3: " << *fw_3 << endl;
          cout << " p1A: " << p1.A << " p1B: " << p1.B << " p1F: " << p1.F << endl;
          cout << " p2A: " << p2.A << " p2B: " << p2.B << " p2F: " << p2.G << endl;
          cout << " p3A: " << p3.A << " p3B: " << p3.B << " p3F: " << p3.H << endl;*/
          counter = false;
          break;
        }
        else {
          cout << "NR didn't converge: " << i << endl;
          //cout << "tols: " << tol1 << "  " << tol2 << "  " << tol3 << endl;
          //cout << " s1: " << _seed1 << " s2: " << _seed2 << " s3: " << _seed3 << " f1: " << *fw
          //     << " fw2: " << *fw_2 << "fw3: " << *fw_3 << endl;
          //cout << " p1A: " << p1.A << " p1B: " << p1.B << " p1F: " << p1.F << endl;
          //cout << " p2A: " << p2.A << " p2B: " << p2.B << " p2F: " << p2.G << endl;
          //cout << " p3A: " << p3.A << " p3B: " << p3.B << " p3F: " << p3.H << endl;

          /// try adjusting seeds to see if it helps NR to converge:
          switch (i) {
            case 0:
              _seed1 = seed1 * seedMax;
              MAX_ITERATIONS = 600;
              break;
            case 1:
              _seed1 = seed1;
              _seed2 = seed2 * seedMax;
              break;
            case 2:
              _seed2 = seed2;
              _seed3 = seed3 * seedMax;
              break;
            case 3:
              _seed3 = seed3;
              _seed1 = -seed1 * seedMax;
              break;
            case 4:
              _seed1 = seed1;
              _seed2 = -seed2 * seedMax;
              break;
            case 5:
              _seed2 = seed2;
              _seed3 = -seed3 * seedMax;
              break;
            case 6:
              _seed1 = seed1 * seedMax;
              _seed2 = seed2 * seedMax;
              _seed3 = seed3;
              break;
            case 7:
              _seed1 = seed1 * seedMax;
              _seed2 = seed2;
              _seed3 = seed3 * seedMax;
              break;
            case 8:
              _seed1 = seed1;
              _seed2 = seed2 * seedMax;
              _seed3 = seed3 * seedMax;
              break;
            case 9:
              _seed1 = -seed1 * seedMax;
              _seed2 = -seed2 * seedMax;
              _seed3 = seed3;
              break;
            case 10:
              _seed1 = -seed1 * seedMax;
              _seed2 = seed2;
              _seed3 = -seed3 * seedMax;
              break;
            case 11:
              _seed1 = seed1;
              _seed2 = -seed2 * seedMax;
              _seed3 = -seed3 * seedMax;
              break;
            case 12:
              _seed1 = seed1 * seedMax;
              _seed2 = seed2 * seedMax;
              _seed3 = seed3 * seedMax;
              break;
            case 13:
              _seed1 = -seed1 * seedMax;
              _seed2 = -seed2 * seedMax;
              _seed3 = -seed3 * seedMax;
              break;
            default:
              cout << "NR Random Seeds !" << endl;
              _seed1 = getRandomNumber();
              _seed2 = getRandomNumber();
              _seed3 = getRandomNumber();
          //    cout << "Getting Random seeds " << i<< endl;
              break;
          }
          counter = true;
        }
      }

      if (counter) { ///TODO if it still doesn't converge, consider not updating the fields or something...
        cout << "FALSE "  << endl;
        cout << "FIz: " << fwxInitial << " FOz: " << *fw << endl;
        cout << "NR didn't converge: " << endl;
        cout << "tols: " << tol1 << "  " << tol2 << "  " << tol3 << endl;
        cout << " s1: " << _seed1 << " s2: " << _seed2 << " s3: " << _seed3 << " f1: " << *fw
             << " fw2: " << *fw_2 << " fw3: " << *fw_3 << endl;
        cout << " p1A: " << p1.A << " p1B: " << p1.B << " p1F: " << p1.F << endl;
        cout << " p2A: " << p2.A << " p2B: " << p2.B << " p2F: " << p2.G << endl;
        cout << " p3A: " << p3.A << " p3B: " << p3.B << " p3F: " << p3.H << endl;
        sleep(7);
      }

     // cout << "FIz: " << fwxInitial << " FOz: " << *fw<<endl;


  }

}




/// OLD 'for' loop:

// for (int i = 2, imax = 270000; i < imax; i*=3) {
//  if (newtonRaphson(seed1+3.06*i, seed2-2.43*i, seed3+1.277*i, p1, p2, p3, fw, fw_2, fw_3, tol1,
//  tol2, tol3)) {
////   if (newtonRaphson(seed1, seed2, seed3, p1, p2, p3, fw, fw_2, fw_3)) {
///*    if (counter) { cout << "True " << i << endl;
//    sleep(12);
//    }*/
//     counter = false;
//     break;
//   }
//   else {
//     cout << "NR didn't converge: " << i << endl;
//     cout << "tols: " << tol1 << "  " << tol2 << "  " << tol3 << endl;
//     cout << " s1: " << seed1 << " s2: " << seed2 << " s3: " << seed3 << " f1: " << *fw
//          << " fw2: " << *fw_2 << "fw3: " << *fw_3 << endl;
//     cout << " p1A: " << p1.A << " p1B: " << p1.B << " p1F: " << p1.F << endl;
//     cout << " p2A: " << p2.A << " p2B: " << p2.B << " p2F: " << p2.G << endl;
//     cout << " p3A: " << p3.A << " p3B: " << p3.B << " p3F: " << p3.H << endl;
//
//     counter = true;
//   }
// }




////OLD full NR code before stripping stuff for efficiency:
//namespace meep {
//
//const double TOLERANCE = 1e-8;
//int MAX_ITERATIONS = 500;
//
////// Struct to hold constant parameters for each equation
//// struct Parameters { // moved to .hpp??!
////   double A, B, C, D, E, F, G, H;
//// };
//
//// Function definitions
//vector<double> equations(realnum x, realnum y, realnum z, const Parameters &p1,
//                         const Parameters &p2, const Parameters &p3);
//vector<vector<double> > jacobian(realnum x, realnum y, realnum z, const Parameters &p1,
//                                 const Parameters &p2, const Parameters &p3);
//vector<double> solveLinearSystem(const vector<vector<double> > &J, const vector<double> &F);
//
//// Compute the coefficient matrix (ignoring constant terms)
//vector<vector<double> > computeCoefficientMatrix(const Parameters &p1, const Parameters &p2,
//                                                 const Parameters &p3) {
//  vector<vector<double> > M(3, vector<double>(3));
//  M[0][0] = -p1.B;
//  M[0][1] = -p1.H;
//  M[0][2] = -p1.G;
//  M[1][0] = -p2.B;
//  M[1][1] = -p2.H;
//  M[1][2] = -p2.G;
//  M[2][0] = -p3.B;
//  M[2][1] = -p3.H;
//  M[2][2] = -p3.G;
//  return M;
//}
//
//// Compute determinant of a 3x3 matrix
//double determinant3x3(const vector<vector<double> > &M) {
//  return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) -
//         M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
//         M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
//}
//
//// Compute rank of a 3x3 matrix (naive approach)
//int rank3x3(const vector<vector<double> > &M) {
//  double detM = determinant3x3(M);
//  if (detM != 0) return 3;
//
//  // Check if any 2x2 minor is nonzero (indicating rank 2)
//  for (int i = 0; i < 3; i++) {
//    for (int j = 0; j < 3; j++) {
//      vector<vector<double> > minorM;
//      for (int m = 0; m < 3; m++) {
//        if (m == i) continue;
//        vector<double> row;
//        for (int n = 0; n < 3; n++) {
//          if (n == j) continue;
//          row.push_back(M[m][n]);
//        }
//        minorM.push_back(row);
//      }
//      if (minorM[0][0] * minorM[1][1] - minorM[0][1] * minorM[1][0] != 0) { return 2; }
//    }
//  }
//  return 1;
//}
//
//bool newtonRaphson(realnum x, realnum y, realnum z, const Parameters &p1, const Parameters &p2,
//                   const Parameters &p3, realnum *fw, realnum *fw_2, realnum *fw_3, double tol1,
//                   double tol2, double tol3) {
//  for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
//    vector<double> F = equations(x, y, z, p1, p2, p3);
//    vector<vector<double> > J = jacobian(x, y, z, p1, p2, p3);
//
//    // CHECK. The next calc of detJ and check are only necessary for checking purposes, not the NR
//    // alg.
//     double detJ = determinant3x3(J);
//     if (detJ == 0) {
//      cout << "Jacobian determinant is 0: The system has dependent or singular solutions at (" <<
//      x
//           << ", " << y << ", " << z << ")!" << endl;
//    }
//    // END CHECK
//
//    vector<double> delta = solveLinearSystem(J, F);
//
//    x -= delta[0];
//    y -= delta[1];
//    z -= delta[2];
//
//    // if (fabs(delta[0]) < TOLERANCE && fabs(delta[1]) < TOLERANCE && fabs(delta[2]) < TOLERANCE) {
//    if (fabs(delta[0]) < tol1 && fabs(delta[1]) < tol2 && fabs(delta[2]) < tol3) {
//      //   cout << "Converged after " << iter + 1 << " iterations:   "; //.\n";
//      //   cout << "x = " << x << ", y = " << y << ", z = " << z << "\n";
//      *fw = x; /// Update E field values!
//      *fw_2 = y;
//      *fw_3 = z;
//
//      // TODO call eqns here with the final values to verify... (because current stopping condition
//      // is just NR algo step delta)
//      return true;
//    }
//  }
//
//  // cout << "Newton's method did not converge within " << MAX_ITERATIONS << " iterations.\n";
//  return false;
//}
//
//vector<double> equations(double x, double y, double z, const Parameters &p1, const Parameters &p2,
//                         const Parameters &p3) {
//   return {p1.A - (p1.B * x + p1.C * x * x + p1.D * y * y + p1.E * z * z + p1.F * y * z +
//                   p1.G * x * z + p1.H * x * y),
//           p2.A - (p2.B * y + p2.C * x * x + p2.D * y * y + p2.E * z * z + p2.F * y * z +
//                   p2.G * x * z + p2.H * x * y),
//           p3.A - (p3.B * z + p3.C * x * x + p3.D * y * y + p3.E * z * z + p3.F * y * z +
//                   p3.G * x * z + p3.H * x * y)};
//}
//
//vector<vector<double> > jacobian(double x, double y, double z, const Parameters &p1,
//                                 const Parameters &p2, const Parameters &p3) {
//  return {{-p1.B - 2 * p1.C * x - p1.G * z - p1.H * y, -2 * p1.D * y - p1.F * z - p1.H * x,
//           -2 * p1.E * z - p1.F * y - p1.G * x},
//          {-2 * p2.C * x - p2.G * z - p2.H * y, -p2.B - 2 * p2.D * y - p2.F * z - p2.H * x,
//           -2 * p2.E * z - p2.F * y - p2.G * x},
//          {-2 * p3.C * x - p3.G * z - p3.H * y, -2 * p3.D * y - p3.F * z - p3.H * x,
//           -p3.B - 2 * p3.E * z - p3.F * y - p3.G * x}};
//}
//
//vector<double> solveLinearSystem(const vector<vector<double> > &J, const vector<double> &F) {
//  vector<vector<double> > A = J;
//  vector<double> b = F;
//  int n = 3;
//
//  for (int i = 0; i < n; i++) {
//    for (int j = i + 1; j < n; j++) {
//      double factor = A[j][i] / A[i][i];
//      for (int k = i; k < n; k++) {
//        A[j][k] -= factor * A[i][k];
//      }
//      b[j] -= factor * b[i];
//    }
//  }
//
//  vector<double> x(n);
//  for (int i = n - 1; i >= 0; i--) {
//    x[i] = b[i];
//    for (int j = i + 1; j < n; j++) {
//      x[i] -= A[i][j] * x[j];
//    }
//    x[i] /= A[i][i];
//  }
//  return x;
//}
//
//// static uniform_real_distribution<double> dist(-1, 1); // Define range
//static uniform_real_distribution<double> distSign(-1.0, 1.0); // Define range
//std::lognormal_distribution<double> dist(1.0, 90);
//static random_device rd; // Obtain a random seed
//static double seedMax = 1e125;
//
//double getRandomNumber() {
//  static mt19937 gen(rd());         // Seed the generator
//  return dist(gen) * distSign(gen); // *1e80;
//                                    // return dist(gen) * 1e80;
//}
//
//void runNR(realnum seed1, realnum seed2, realnum seed3, realnum *fw, realnum *fw_2, realnum *fw_3,
//           const Parameters &p1, const Parameters &p2,
//           const Parameters &p3) { // TODO need to confirm that passing fw through as a ref like
//                                   // this actually works...
//  std::cout << std::fixed << std::setprecision(12);
//  //  cout << "Doing NR" << endl;
//
//  double fwxInitial = *fw;
//  MAX_ITERATIONS = 250;
//  double tol1 = fmax(fabs(TOLERANCE * (*fw)) * 0.0001, TOLERANCE);
//  double tol2 = fmax(fabs(TOLERANCE * (*fw_2)) * 0.0001, TOLERANCE);
//  double tol3 = fmax(fabs(TOLERANCE * (*fw_3)) * 0.0001, TOLERANCE);
//  double _seed1 = seed1;
//  double _seed2 = seed2;
//  double _seed3 = seed3;
//
//  // CHECK 2:
//   vector<vector<double> > M = computeCoefficientMatrix(p1, p2, p3);
//   int rankM = rank3x3(M);
//   if (rankM < 3) {
//    cout << "Coefficient matrix has rank < 3: The system is globally dependent!" << endl;
//    cout << " s1: " << seed1 << " s2: " << seed2 << " s3: " << seed3 << " f1: " << *fw << " fw2: "
//         << *fw_2 << " fw3: " << *fw_3 << endl;
//    cout << " p1A: " << p1.A << " p1B: " << p1.B << " p1F: " << p1.F << endl;
//    cout << " p2A: " << p2.A << " p2B: " << p2.B << " p2G: " << p2.G << endl;
//    cout << " p3A: " << p3.A << " p3B: " << p3.B << " p3H: " << p3.H << endl;
//    sleep(4);
//
//   return;
//   }
//  //  END CHECK 2
//
//  // Parameters p1 = {a, 3.2, 0, 0, 0, 0.000002, 0, 0}; // p1.A is the D-P field p1.B is epsilon,
//  // then chi2 Parameters p2 = {b, 3.2, 0, 0, 0, 0, 0.000002, 0}; Parameters p3 = {c, 3.2, 0, 0, 0,
//  // 0, 0, 0.000002};
//  bool counter = false;
//  for (int i = 0, imax = 100; i < imax; ++i) {
//    if (newtonRaphson(_seed1, _seed2, _seed3, p1, p2, p3, fw, fw_2, fw_3, tol1, tol2, tol3)) {
//      //   if (newtonRaphson(seed1, seed2, seed3, p1, p2, p3, fw, fw_2, fw_3)) {
//      /*     cout << "Converged " << i << endl;
//           cout << "tols: " << tol1 << "  " << tol2 << "  " << tol3 << endl;
//           cout << " s1: " << _seed1 << " s2: " << _seed2 << " s3: " << _seed3 << " f1: " << *fw
//                << " fw2: " << *fw_2 << "fw3: " << *fw_3 << endl;
//           cout << " p1A: " << p1.A << " p1B: " << p1.B << " p1F: " << p1.F << endl;
//           cout << " p2A: " << p2.A << " p2B: " << p2.B << " p2F: " << p2.G << endl;
//           cout << " p3A: " << p3.A << " p3B: " << p3.B << " p3F: " << p3.H << endl;*/
//      counter = false;
//      break;
//    }
//    else {
//      // cout << "NR didn't converge: " << i << endl;
//      // cout << "tols: " << tol1 << "  " << tol2 << "  " << tol3 << endl;
//      // cout << " s1: " << _seed1 << " s2: " << _seed2 << " s3: " << _seed3 << " f1: " << *fw
//      //      << " fw2: " << *fw_2 << "fw3: " << *fw_3 << endl;
//      // cout << " p1A: " << p1.A << " p1B: " << p1.B << " p1F: " << p1.F << endl;
//      // cout << " p2A: " << p2.A << " p2B: " << p2.B << " p2F: " << p2.G << endl;
//      // cout << " p3A: " << p3.A << " p3B: " << p3.B << " p3F: " << p3.H << endl;
//
//      /// try adjusting seeds to see if it helps NR to converge:
//      switch (i) {
//        case 0:
//          _seed1 = seed1 * seedMax;
//          MAX_ITERATIONS = 600;
//          break;
//        case 1:
//          _seed1 = seed1;
//          _seed2 = seed2 * seedMax;
//          break;
//        case 2:
//          _seed2 = seed2;
//          _seed3 = seed3 * seedMax;
//          break;
//        case 3:
//          _seed3 = seed3;
//          _seed1 = -seed1 * seedMax;
//          break;
//        case 4:
//          _seed1 = seed1;
//          _seed2 = -seed2 * seedMax;
//          break;
//        case 5:
//          _seed2 = seed2;
//          _seed3 = -seed3 * seedMax;
//          break;
//        case 6:
//          _seed1 = seed1 * seedMax;
//          _seed2 = seed2 * seedMax;
//          _seed3 = seed3;
//          break;
//        case 7:
//          _seed1 = seed1 * seedMax;
//          _seed2 = seed2;
//          _seed3 = seed3 * seedMax;
//          break;
//        case 8:
//          _seed1 = seed1;
//          _seed2 = seed2 * seedMax;
//          _seed3 = seed3 * seedMax;
//          break;
//        case 9:
//          _seed1 = -seed1 * seedMax;
//          _seed2 = -seed2 * seedMax;
//          _seed3 = seed3;
//          break;
//        case 10:
//          _seed1 = -seed1 * seedMax;
//          _seed2 = seed2;
//          _seed3 = -seed3 * seedMax;
//          break;
//        case 11:
//          _seed1 = seed1;
//          _seed2 = -seed2 * seedMax;
//          _seed3 = -seed3 * seedMax;
//          break;
//        case 12:
//          _seed1 = seed1 * seedMax;
//          _seed2 = seed2 * seedMax;
//          _seed3 = seed3 * seedMax;
//          break;
//        case 13:
//          _seed1 = -seed1 * seedMax;
//          _seed2 = -seed2 * seedMax;
//          _seed3 = -seed3 * seedMax;
//          break;
//        default:
//          cout << "NR Random Seeds !" << endl;
//          _seed1 = getRandomNumber();
//          _seed2 = getRandomNumber();
//          _seed3 = getRandomNumber();
//          //    cout << "Getting Random seeds " << i<< endl;
//          break;
//      }
//      counter = true;
//    }
//  }
//
//  if (counter) { /// TODO if it still doesn't converge, consider not updating the fields or
//                 /// something...
//    cout << "FALSE " << endl;
//    cout << "FIz: " << fwxInitial << " FOz: " << *fw << endl;
//    cout << "NR didn't converge: " << endl;
//    cout << "tols: " << tol1 << "  " << tol2 << "  " << tol3 << endl;
//    cout << " s1: " << _seed1 << " s2: " << _seed2 << " s3: " << _seed3 << " f1: " << *fw
//         << " fw2: " << *fw_2 << " fw3: " << *fw_3 << endl;
//    cout << " p1A: " << p1.A << " p1B: " << p1.B << " p1F: " << p1.F << endl;
//    cout << " p2A: " << p2.A << " p2B: " << p2.B << " p2F: " << p2.G << endl;
//    cout << " p3A: " << p3.A << " p3B: " << p3.B << " p3F: " << p3.H << endl;
//    sleep(7);
//  }
//
//  // cout << "FIz: " << fwxInitial << " FOz: " << *fw<<endl;
//}
//
//} // namespace meep