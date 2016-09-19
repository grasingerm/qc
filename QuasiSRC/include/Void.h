//
// Void.h
//

#if !defined(VOID_H)
#define VOID_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#ifdef HAVE_VECTOR
#include <vector>
#else
#ifdef HAVE_VECTOR_H
#include <vector.h>
#else
#error No vector or vector.h available
#endif // HAVE_VECTOR_H
#endif // HAVE_VECTOR

#include "DataTypes.h"

//
//
//

namespace quasicontinuum {

/**
 * @brief Singleton container for Quasicontinua instance
 */
class Void {

  //
  // public methods
  //

public:
  /**
   * @brief getInstance.
   */
  static Void *getInstance();

  /**
   * @brief destroyInstance.
   */
  static void destroyInstance();

  //
  //  createVoid()
  //
  void createVoid(struct node_list_t *P_node_list,
                  const struct lattice_t *P_lattice, const int iQuasi);

  //
  //  removeVoidElements()
  //
  void removeVoidElements(struct element_list_t *P_element_list);

  //
  //  removeVoidNodes()
  //
  void removeVoidNodes(struct node_list_t *P_node_list);

  //
  //  insertVoidNodes()
  //
  void insertVoidNodes(struct node_list_t *P_node_list);

  //
  //  setVoid()
  //
  void setVoid(const int enable_void_flag);

  //
  //  isVoidEnable()
  //
  int isVoidEnable();

  //
  //  setVoidParameters()
  //
  void setVoidParameters(const double center[3], const enum void_t voidType,
                         const int number_parameters,
                         std::vector<double> param_list);

  //
  //  insertVoidSiteCache()
  //
  void insertVoidSiteCache(const std::vector<int> site, const int iQuasi);

  //
  //  clearVoidCache()
  //
  void clearVoidSiteCache();

  //
  //  findSiteInVoidCache()
  //  0 : site not a void
  //  1 : site is void
  //
  int findSiteInVoidCache(const std::vector<int> site, const int iQuasi);

  //
  //  printVoidCache()
  //
  void printVoidCache();

  //
  //  voidAtomsMissed()
  //
  std::vector<std::pair<int, std::vector<int>>> voidAtomsMissed();

  //
  //  fixVoidBoundary()
  //
  void fixVoidBoundary(struct node_list_t *P_node_list,
                       struct element_list_t *P_element_list);

  //
  // private methods
  //

private:
  /**
   * @brief Constructor.
   */
  Void();

  /**
   * @brief Copy constructor.
   */
  Void(Void const &);

  /**
   * @brief Assignment operator.
   */
  const Void &operator=(const Void &);

  /**
   * @brief Destructor.
   */
  ~Void();

  //
  //  constructor that will be called by getInstance()
  //
  Void(const int enable_void_default);

  //
  // private data types
  //
private:
  static Void *_instance;
  int d_enableVoid;
  std::vector<std::vector<std::vector<int>>> d_voidSiteCache;
  std::vector<double> d_voidCenter;
  double d_voidRadius;
  std::vector<double> d_voidParameters;
  enum void_t d_voidType;
};
}

#endif // VOID_H
