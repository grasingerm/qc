//
// File:      CGNonLinearSolver.h
// Package:   solvers
//
// Various generic solvers.
//
#if !defined(CGNONLINEARSOLVER_H)
#define CGNONLINEARSOLVER_H

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif // HAVE_CONFIG_H

#if defined(HAVE_VECTOR)
#include <vector>
#else
#error vector header file not available
#endif // HAVE_VECTOR

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#else
#error No pthread.h available.
#endif /* HAVE_PTHREAD_H */

#include <utility>

#include "DataTypes.h"
#include "NonLinearSolver.h"

//
//
//

namespace quasicontinuum {

/**
 * @brief Concrete class implementing Conjugate Gradient non-linear
 * algebraic solver.
 */
class CGNonLinearSolver : public NonLinearSolver {

  //
  // types
  //
public:
  //
  // methods
  //
public:
  /**
   * @brief Constructor.
   *
   * @param tolerance Tolerance (relative) required for convergence.
   * @param maxNumberIterations Maximum number of iterations.
   * @param debugLevel Debug output level:
   *                   0 - no debug output
   *                   1 - limited debug output
   *                   2 - all debug output.
   * @param lineSearchTolerance Tolereance required for line search
   * convergence.
   * @param lineSearchMaxIterations Maximum number of iterations for the
   * line search.
   */
  CGNonLinearSolver(double tolerance, int maxNumberIterations, int debugLevel,
                    double lineSearchTolerance = 1.0e-6,
                    int lineSearchMaxIterations = 10,
                    double rebuild_cutoff = 0.0, int minMethod = 0,
                    int minMethodMaxIterations = 1,
                    double minMethodTolerance = 1.0e-6);

  /**
   * @brief Destructor.
   */
  virtual ~CGNonLinearSolver();

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
   * @param output_flag To output the data in between iteration
   *                    0 - no output
   *                    1 - output
   *
   * @return Return value indicating success or failure.
   */
  virtual NonLinearSolver::ReturnValueType
  solve(SolverFunction &function, double tolerance, int maxNumberIterations,
        int debugLevel = 0, int output_flag = 0, int loadNumber = 0,
        int minMethod = 0, int minMethodMaxIterations = 1,
        double minMethodTolerance = 1.0e-6);

  //
  //  Free energy minimization wrt u and w simultaneously
  //
  virtual NonLinearSolver::ReturnValueType
  CGMinimization(SolverFunction &function, double tolerance,
                 int maxNumberIterations, int debugLevel = 0,
                 int output_flag = 0, int loadNumber = 0);

  //
  //  Free energy minimization wrt u and w alternatively
  //
  virtual NonLinearSolver::ReturnValueType CGAlternateMinimization(
      SolverFunction &function, double tolerance, int maxNumberIterations,
      int debugLevel = 0, int output_flag = 0, int loadNumber = 0,
      int minMethodMaxIterations = 1, double minMethodTolerance = 1.0e-6);

private:
  //
  // copy constructor/assignment operator
  //
  CGNonLinearSolver(const CGNonLinearSolver &);            // not implemented
  CGNonLinearSolver &operator=(const CGNonLinearSolver &); // not implemented

  /**
   * @brief Initialize direction.
   */
  void initializeDirection();

  /**
   * @brief Perform line search.
   *
   * @param function SolverFunction object (functor) to compute energy and
   *                 forces.
   * @param tolerance Tolerance (relative) required for convergence.
   * @param maxNumberIterations Maximum number of iterations.
   * @param debugLevel Debug output level:
   *                   0 - no debug output
   *                   1 - limited debug output
   *                   2 - all debug output.
   *
   * @return Return value indicating success or failure.
   */
  NonLinearSolver::ReturnValueType lineSearch(SolverFunction &function,
                                              double tolerance,
                                              int maxNumberIterations,
                                              int debugLevel);

  /**
   * @brief Compute delta_d and eta_p.
   *
   & @return Pair containing delta_d and eta_p.
  */
  std::pair<double, double> computeDeltaD();

  /**
   * @brief Compute eta.
   *
   * @return Value of eta.
   */
  double computeEta();

  /**
   * @brief Compute delta new and delta mid.
   */
  void computeDeltas();

  /**
   * @brief Update direction.
   */
  void updateDirection();

  /**
   * @brief Compute residual L2-norm.
   *
   * @return Value of the residual L2-norm.
   */
  double computeResidualL2Norm() const;

  /**
   * @brief Compute the total number of unknowns in all
   * processors.
   *
   * @return Number of unknowns in all processors.
   */
  int computeTotalNumberUnknowns() const;

  /**
   * @brief Update solution x -> x + \alpha direction.
   *
   * @param alpha Scalar value.
   * @param direction Direction vector.
   * @param function Solver function object.
   */
  void updateSolution(double alpha, const std::vector<double> &direction,
                      SolverFunction &function);

  //
  // data
  //
private:
  std::vector<double> d_solution;
  std::vector<double> d_direction;
  mutable std::vector<double> d_gradient;
  std::vector<double> d_s;
  std::vector<double> d_sOld;
  double d_deltaNew;
  double d_deltaMid;
  double d_deltaOld;
  double d_beta;
  int d_numberUnknowns;
  int d_iter;
  mutable std::vector<int> d_unknownCountFlag;
  double d_lineSearchTolerance;
  int d_lineSearchMaxIterations;
  mutable pthread_mutex_t d_lock;
  enum quasicontinuum::mt_version_t d_version;
  const int REBUILD;
  const int NO_REBUILD;
  double d_accumulated_r;
  double d_rebuild_cutoff;
  int d_minMethod;
  int d_minMethodMaxIterations;
  double d_minMethodTolerance;
  int d_minimizingMethod;

  //
  // data
  //
private:
};
}

#endif // CGNONLINEARSOLVER_H
