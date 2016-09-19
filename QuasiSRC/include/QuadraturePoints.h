//
// QuadraturePoints.h
//

#if !defined(QUADRATUREPOINTS_H)
#define QUADRATUREPOINTS_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

// vector
#ifdef HAVE_VECTOR
#include <vector>
#else
#ifdef HAVE_VECTOR_H
#include <vector.h>
#else
#error No vector or vector.h available
#endif // HAVE_VECTOR_H
#endif // HAVE_VECTOR

//
//
//

namespace quasicontinuum {

/**
 * @brief Singleton container for Quasicontinua instance
 */
class QuadraturePoints {

  //
  // public methods
  //

public:
  /**
   * @brief getInstance.
   */
  static QuadraturePoints *getInstance();

  /**
   * @brief destroyInstance.
   */
  static void destroyInstance();

  /**
   * @brief Write the quadrature points on matrix quadVectors.
   *
   * @param pointer to matrix a of size n BY 3.
   * @param number of 3 dimensional points we need : n.
   */
  void getQuadraturePoints(
      std::vector<std::vector<std::vector<double>>> &quadVectors,
      std::vector<double> &quadWeights, int n);

  //
  //  get Quadrature vectors and weights
  //  dimension : it is the space in which individual vectors belong to.
  //  number : integration space is in R^(number * dimension)
  //
  std::vector<std::pair<std::vector<std::vector<double>>, double>>
  getQuadratureVectors(int number, int dimension);

  //
  // private methods
  //

private:
  /**
   * @brief Constructor.
   */
  QuadraturePoints();

  /**
   * @brief Copy constructor.
   */
  QuadraturePoints(QuadraturePoints const &);

  /**
   * @brief Assignment operator.
   */
  const QuadraturePoints &operator=(const QuadraturePoints &);

  /**
   * @brief Destructor.
   */
  ~QuadraturePoints();

  //
  // private data types
  //
private:
  static QuadraturePoints *_instance;
};
}

#endif // QUADRATUREPOINTS_H
