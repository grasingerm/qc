//
// Shape.h
//

#if !defined(SHAPE_H)
#define SHAPE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include "DataTypes.h"

//
//
//

namespace quasicontinuum {

/**
 * @brief Singleton container for Quasicontinua instance
 */
class Shape {

  //
  // public methods
  //

public:
  /**
   * @brief getInstance.
   */
  static Shape *getInstance();

  /**
   * @brief destroyInstance.
   */
  static void destroyInstance();

  //
  //  computeShapeFunction()
  //
  double computeShapeFunction(const double p1[3], const double p2[3],
                              const double p3[3], const double p4[3],
                              const double p[3]);

  //
  //  computeSiteShapeFunction()
  //
  double computeSiteShapeFunction(const struct node_t *P_node_0,
                                  const struct node_t *P_node_1,
                                  const struct node_t *P_node_2,
                                  const struct node_t *P_node_3, const int l[3],
                                  struct lattice_t *P_lattice,
                                  const int iQuasi);

  //
  //  computeShapeFunctionGradient()
  //
  void computeShapeFunctionGradient(const double p1[3], const double p2[3],
                                    const double p3[3], const double p4[3],
                                    double gradient[3]);

  //
  // private methods
  //

private:
  /**
   * @brief Constructor.
   */
  Shape();

  /**
   * @brief Copy constructor.
   */
  Shape(Shape const &);

  /**
   * @brief Assignment operator.
   */
  const Shape &operator=(const Shape &);

  /**
   * @brief Destructor.
   */
  ~Shape();

  //
  // private data types
  //
private:
  static Shape *_instance;
};
}

#endif // SHAPE_H
