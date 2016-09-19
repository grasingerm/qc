//
// File:      NonLinearSolver.cc
// Package:   solvers
//
// Various generic solvers.
//
#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif // HAVE_CONFIG_H

#include "NonLinearSolver.h"

//
//
//

namespace quasicontinuum {

//
// Constructor.
//
NonLinearSolver::NonLinearSolver(double tolerance, int maxNumberIterations,
                                 int debugLevel)
    : d_debugLevel(debugLevel), d_maxNumberIterations(maxNumberIterations),
      d_tolerance(tolerance) {

  //
  //
  //
  return;
}

//
// Destructor.
//
NonLinearSolver::~NonLinearSolver() {

  //
  //
  //
  return;
}

//
// Solve non-linear algebraic equation.
//
NonLinearSolver::ReturnValueType
NonLinearSolver::solve(SolverFunction &function) {

  //
  //
  //
  return this->solve(function, d_tolerance, d_maxNumberIterations,
                     d_debugLevel);
}
}
