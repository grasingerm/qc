//
// CreateMesh.h
//

#if !defined(CREATEMESH_H)
#define CREATEMESH_H

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
class CreateMesh {

  //
  // public methods
  //

public:
  //
  // struct data type
  struct edge_t {
    double length;
    struct node_t *node[2];
  };

  //
  //  score_list_t
  struct score_list_t {
    double energy_drop;
    struct element_t *P_element;
  };

  /**
   * @brief getInstance.
   */
  static CreateMesh *getInstance();

  /**
   * @brief destroyInstance.
   */
  static void destroyInstance();

  //
  //  setRemeshTolerance()
  //
  void setRemeshTolerance(double tolerance);

  //
  //  createElementsCentral()
  //
  //  computeForce : 0 - ON, 1 - OFF
  //
  void createElementsCentral(int computeForce);

  //
  //  createElements()
  //
  void createElements(int iQuasi, struct element_list_t *P_element_list,
                      struct node_list_t *P_node_list,
                      struct lattice_t *P_lattice, double r_cluster_cutoff,
                      int restartFlag);

  //
  //  remeshCentral()
  //
  int remeshCentral(int i_load);

  //
  //  remesh()
  //
  int remesh(int iQuasi, struct element_list_t *P_element_list,
             struct all_node_list_t *P_node_list, struct lattice_t *P_lattice,
             struct indentor_t *P_indentor, double r_cluster_cutoff,
             int i_load);
  //
  //  fixBoundary()
  //
  void fixBoundary(struct node_list_t *P_node_list,
                   struct element_list_t *P_element_list,
                   struct lattice_t *P_lattice, int restartFlag);

  //
  //  private definitions
  //
private:
  //
  // private methods
  //

private:
  /**
   * @brief Constructor.
   */
  CreateMesh();

  //
  //  Constructor
  //
  CreateMesh(const double remesh_tolerance);

  /**
   * @brief Copy constructor.
   */
  CreateMesh(CreateMesh const &);

  /**
   * @brief Assignment operator.
   */
  const CreateMesh &operator=(const CreateMesh &);

  /**
   * @brief Destructor.
   */
  ~CreateMesh();

  //
  //  remeshElementBisection()
  //
  void RemeshElementBisection(struct element_t *P_element,
                              struct all_node_list_t *P_node_list,
                              struct lattice_t *P_lattice, const int iQuasi);

  //
  //  FindBisectionSite()
  //
  int FindBisectionSite(int l[3], struct lattice_t *P_lattice,
                        const struct element_t *P_element, const int iQuasi);

  //
  //  FindCommonEdgeElements()
  //
  void FindCommonEdgeElements(std::vector<struct element_t *> *P_neigh_elements,
                              const struct node_t *P_node_1,
                              const struct node_t *P_node_2);

  //
  //  CompareEdges()
  //
  bool CompareEdges(const edge_t &lhs, const edge_t &rhs);

  //
  //  CheckEnergyRemesh()
  //
  void CheckEnergyRemesh(std::vector<score_list_t> *P_score_lists,
                         struct element_list_t *P_element_list,
                         struct all_node_list_t *P_node_list,
                         struct indentor_t *P_indentor, int i_load, int iQuasi);

  //
  //  CompareScoreList()
  //
  bool CompareScoreList(const score_list_t &lhs, const score_list_t &rhs);

  //
  //  DumpScoreList()
  //
  void DumpScoreList(std::vector<score_list_t> score_list);

  //
  // private data types
  //
private:
  static CreateMesh *_instance;
  double d_remeshTolerance;
};
}

#endif // CREATEMESH_H
