//
// File:      SolverFunction.h
// Package:   solvers
//
// Various generic solvers.
//
#if !defined(solvers_non_linear_SolverFunction_h)
#define solvers_non_linear_SolverFunction_h

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif // HAVE_CONFIG_H

#include <vector>

namespace quasicontinuum {

/**
 * @brief Abstract class for solver functions.
 */
class SolverFunction {

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
   */
  SolverFunction();

  /**
   * @brief Destructor.
   */
  virtual ~SolverFunction() = 0;

  /**
   * @brief Obtain number of unknowns.
   *
   * @return Number of unknowns.
   */
  virtual int getNumberUnknowns() const = 0;

  /**
   * @brief Compute function value (aka energy).
   *
   *
   * @return Function value.
   */
  virtual double value() const = 0;

  /**
   * @brief Compute function gradient (aka forces).
   *
   * @param gradient STL vector for gradient values.
   */
  // virtual void gradient(std::vector<double> & gradient) const = 0;

  /**
   * @brief Compute function gradient (aka forces).
   *
   * @param gradient STL vector for gradient values.
   * @param flag Flag for function, necessary for rebuiling neighbor lists
   */
  virtual void gradient(std::vector<double> &gradient,
                        const int flag) const = 0;

  /**
   * @brief Apply preconditioner to function gradient.
   *
   * @param s STL vector for s=-M^{-1} gradient.
   * @param gradient STL vector for gradient values.
   */
  virtual void precondition(std::vector<double> &s,
                            const std::vector<double> &gradient) const = 0;

  /**
   * @brief Update solution.
   *
   * @param solution Updated solution.
   */
  virtual void update(const std::vector<double> &solution) = 0;

  /**
   * @brief Obtain current solution.
   *
   * @param solution Space for current solution.
   */
  virtual void solution(std::vector<double> &solution) const = 0;

  /**
   * @brief For each unknown indicate whether that unknown should
   * be accumulated. This functionality is needed in the case of
   * parallel execution when domain decomposition is
   * employed. Unknowns residing on processor boundary should only
   * be accumulated once when dot products of vertex fields are
   * computed (e.g. residual).
   *
   * @return A vector of int values for each unknown. Value of 1
   * indicates that the unknown should be counted and 0 otherwise.
   */
  virtual std::vector<int> getUnknownCountFlag() const = 0;

  /**
   * @brief changes position fixity mask
   *
   * fixity_flag = 0 : reset to original state
   #               1 : fix all
   *               2 : free all
   */
  virtual void setPositionFixity(const int fixity_flag) const = 0;

  /**
   * @brief changes frequency fixity mask
   *
   * fixity_flag = 0 : reset to original state
   #               1 : fix all
   *               2 : free all
   */
  virtual void setFrequencyFixity(const int fixity_flag) const = 0;

  /**
   * @brief recomputes the data of SolverFunction instance
   */
  virtual void recomputeData() const = 0;

private:
  //
  // copy constructor/assignment operator
  //
  SolverFunction(const SolverFunction &);            // not implemented
  SolverFunction &operator=(const SolverFunction &); // not implemented

  //
  // data
  //
private:
};
}

#endif // solvers_non_linear_SolverFunction_h
