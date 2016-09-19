//
// File:      QuadraticFunction.h
// Package:   tests/solvers
//
// Solvers tests.
//
#if !defined(QUASICONTINUAFORCEFUNCTION_H)
#define QUASICONTINUAFORCEFUNCTION_H

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif // HAVE_CONFIG_H

#include "SolverFunction.h"

//
//
//

namespace quasicontinuum {

/**
 * @brief Implement a function of
 */
class QuasicontinuaForceFunction : public SolverFunction {

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
   * @param x Initial guess.
   */
  QuasicontinuaForceFunction();

  /**
   * @brief Destructor.
   */
  virtual ~QuasicontinuaForceFunction();

  /**
   * @brief Obtain number of unknowns.
   *
   * @return Number of unknowns.
   */
  virtual int getNumberUnknowns() const;

  /**
   * @brief Compute function value (aka energy).
   *
   *
   * @return Function value.
   */
  virtual double value() const;

  /**
   * @brief Compute function gradient (aka forces).
   *
   * @param gradient STL vector for gradient values.
   */
  // virtual void gradient(std::vector<double> & gradient) const;

  /**
   * @brief Compute function gradient (aka forces).
   *
   * @param gradient STL vector for gradient values.
   * @param flag Flag for rebuilding neighbor lists.
   */
  virtual void gradient(std::vector<double> &gradient, const int flag) const;

  /**
   * @brief Apply preconditioner to function gradient.
   *
   * @param s STL vector for s=-M^{-1} gradient.
   * @param gradient STL vector for gradient values.
   */
  virtual void precondition(std::vector<double> &s,
                            const std::vector<double> &gradient) const;

  /**
   * @brief Update solution.
   *
   * @param solution Updated solution.
   */
  virtual void update(const std::vector<double> &solution);

  /**
   * @brief Obtain current solution.
   *
   * @param solution Space for current solution.
   */
  virtual void solution(std::vector<double> &solution) const;

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
  virtual std::vector<int> getUnknownCountFlag() const;

  /**
   * @brief Get current solution size.
   */
  virtual int getSolutionSize() const;

  /**
   * @brief changes position fixity mask
   *
   * fixity_flag = 0 : reset to original state
   #               1 : fix all
   *               2 : free all
   */
  virtual void setPositionFixity(const int fixity_flag) const;

  /**
   * @brief changes frequency fixity mask
   *
   * fixity_flag = 0 : reset to original state
   #               1 : fix all
   *               2 : free all
   */
  virtual void setFrequencyFixity(const int fixity_flag) const;

  /**
   * @brief recomputes the data of SolverFunction instance
   */
  virtual void recomputeData() const;

private:
  //
  // copy constructor/assignment operator
  //
  QuasicontinuaForceFunction(
      const QuasicontinuaForceFunction &); // not implemented
  QuasicontinuaForceFunction &
  operator=(const QuasicontinuaForceFunction &); // not implemented

  /**
   * @brief Function to get current Quasicontinua solutions.
   *
   * @return Vector of current Quasicontinua solutions.
   */
  virtual void GetQuasicontinuaSolution() const;
  //
  // data
  //
private:
  mutable std::vector<double> d_solution;
  mutable pthread_mutex_t d_lock;
  quasicontinuum::mt_version_t d_version;
  mutable std::vector<int> d_numQuasicontinuumUnknowns;
};
}

#endif // QUASICONTINUAFORCEFUNCTION_H
