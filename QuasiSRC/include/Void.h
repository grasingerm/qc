#if !defined(VOID_H)
#define VOID_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <vector>
#include "DataTypes.h"

namespace quasicontinuum {

/**
 * @brief Singleton container for Quasicontinua instance
 */
class Void {

public:
  /**
   * @brief getInstance.
   */
  static Void *getInstance();

  /**
   * @brief destroyInstance.
   */
  static void destroyInstance();

  void createVoid(struct node_list_t *P_node_list,
                  const struct lattice_t *P_lattice, const int iQuasi);

  void removeVoidElements(struct element_list_t *P_element_list);

  void removeVoidNodes(struct node_list_t *P_node_list);

  void insertVoidNodes(struct node_list_t *P_node_list);

  void setVoid(const int enable_void_flag);

  int isVoidEnable();

  void setVoidParameters(const double center[3], const enum void_t voidType,
                         const int number_parameters,
                         std::vector<double> param_list);

  void insertVoidSiteCache(const std::vector<int> site, const int iQuasi);

  void clearVoidSiteCache();

  //
  //  findSiteInVoidCache()
  //  0 : site not a void
  //  1 : site is void
  //
  int findSiteInVoidCache(const std::vector<int> site, const int iQuasi);

  void printVoidCache();

  std::vector<std::pair<int, std::vector<int>>> voidAtomsMissed();

  void fixVoidBoundary(struct node_list_t *P_node_list,
                       struct element_list_t *P_element_list);

private:
  Void();
  Void(Void const &);
  const Void &operator=(const Void &);
  ~Void();
  Void(const int enable_void_default);

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
