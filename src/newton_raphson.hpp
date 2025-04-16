#ifndef NONLINEAR_SOLVER_HPP
#define NONLINEAR_SOLVER_HPP
#include "meep.hpp"
#include <vector>



namespace meep {

// Struct to hold constant parameters for each equation
struct Parameters {
  double A, B, C, D, E, F, G, H;
};  


// Top-level function to run the Newton-Raphson solver
void runNR(realnum seed1, realnum seed2, realnum seed3, realnum*  fw, realnum*  fw_2, realnum*  fw_3,
           const Parameters &p1, const Parameters &p2, const Parameters &p3);


} // namespace meep



#endif // NONLINEAR_SOLVER_HPP
