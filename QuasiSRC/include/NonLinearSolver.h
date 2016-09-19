//
// File:      NonLinearSolver.h
// Package:   solvers
//
// Various generic solvers.
//
#if !defined(solver_non_linear_NonLinearSolver_h)
#define solver_non_linear_NonLinearSolver_h

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif // HAVE_CONFIG_H

//
//
//

namespace quasicontinuum {

//
// forward declarations
//
class SolverFunction;

/**
 * @brief Base class for non-linear algebraic solver.
 */
class NonLinearSolver {

  //
  // types
  //
public:
  enum ReturnValueType {
    SUCCESS = 0,
    FAILURE,
    LINESEARCH_FAILED,
    MAX_ITER_REACHED
  };

  //
  // methods
  //
public:
  /**
   * @brief Destructor.
   */
  virtual ~NonLinearSolver() = 0;

  /**
   * @brief Solver non-linear algebraic equation.
   *
   * @param function SolverFunction object (functor) to compute energy and
   *
   * @return Return value indicating success or failure.
   */
  ReturnValueType solve(SolverFunction &function);

  /**
   * @brief Solve non-linear algebraic equation.
   *
   * @param function SolverFunction object (functor) to compute energy and
   * forces.
   * @param tolerance Tolerance (relative) required for convergence.
   * @param maxNumberIterations Maximum number of iterations.
   * @param debugLevel Debug output level:
   *                   0 - no debug output
   *                   1 - limited debug output
   *                   2 - all debug output.
   *
   * @return Return value indicating success or failure.
   */
  virtual ReturnValueType solve(SolverFunction &function, double tolerance,
                                int maxNumberIterations, int debugLevel = 0,
                                int output_flag = 0, int loadNumber = 0,
                                int minMethod = 0,
                                int minMethodMaxIterations = 1,
                                double minMethodTolerance = 1.0e-6) = 0;

  //
  //  Free energy minimization wrt u and w simultaneously
  //
  virtual ReturnValueType
  CGMinimization(SolverFunction &function, double tolerance,
                 int maxNumberIterations, int debugLevel = 0,
                 int output_flag = 0, int loadNumber = 0) = 0;

  //
  //  Free energy minimization wrt u and w alternatively
  //
  virtual ReturnValueType CGAlternateMinimization(
      SolverFunction &function, double tolerance, int maxNumberIterations,
      int debugLevel = 0, int output_flag = 0, int loadNumber = 0,
      int minMethodMaxIterations = 1, double minMethodTolerance = 1.0e-6) = 0;

protected:
  /**
   * @brief Constructor.
   *
   * @param tolerance Tolerance (relative) required for convergence.
   * @param maxNumberIterations Maximum number of iterations.
   * @param debugLevel Debug output level:
   *                   0 - no debug output
   *                   1 - limited debug output
   *                   2 - all debug output.
   */
  NonLinearSolver(double tolerance, int maxNumberIterations, int debugLevel);

private:
  //
  // copy constructor/assignment operator
  //
  NonLinearSolver(const NonLinearSolver &);            // not implemented
  NonLinearSolver &operator=(const NonLinearSolver &); // not implemented

  //
  // data
  //
private:
  const int d_debugLevel;
  const int d_maxNumberIterations;
  const double d_tolerance;
};
}

#endif // solver_non_linear_NonLinearSolver_h
