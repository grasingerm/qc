#if !defined(CREATEMESH_H)
#define CREATEMESH_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <vector>
#include "DataTypes.h"

namespace quasicontinuum {

/**
 * @brief Singleton container for Quasicontinua instance
 */
class CreateMesh {

public:
  struct edge_t {
    double length;
    struct node_t *node[2];
  };

  struct score_list_t {
    double energy_drop;
    struct element_t *P_element;
  };

  static CreateMesh *getInstance();

  static void destroyInstance();

  void setRemeshTolerance(double tolerance);

  //
  //  createElementsCentral()
  //
  //  computeForce : 0 - ON, 1 - OFF
  //
  void createElementsCentral(int computeForce);

  void createElements(int iQuasi, struct element_list_t *P_element_list,
                      struct node_list_t *P_node_list,
                      struct lattice_t *P_lattice, double r_cluster_cutoff,
                      int restartFlag);

  int remeshCentral(int i_load);

  int remesh(int iQuasi, struct element_list_t *P_element_list,
             struct all_node_list_t *P_node_list, struct lattice_t *P_lattice,
             struct indentor_t *P_indentor, double r_cluster_cutoff,
             int i_load);
  void fixBoundary(struct node_list_t *P_node_list,
                   struct element_list_t *P_element_list,
                   struct lattice_t *P_lattice, int restartFlag);

private:
  CreateMesh();
  CreateMesh(const double remesh_tolerance);
  CreateMesh(CreateMesh const &);
  const CreateMesh &operator=(const CreateMesh &);
  ~CreateMesh();
  void RemeshElementBisection(struct element_t *P_element,
                              struct all_node_list_t *P_node_list,
                              struct lattice_t *P_lattice, const int iQuasi);
  int FindBisectionSite(int l[3], struct lattice_t *P_lattice,
                        const struct element_t *P_element, const int iQuasi);
  void FindCommonEdgeElements(std::vector<struct element_t *> *P_neigh_elements,
                              const struct node_t *P_node_1,
                              const struct node_t *P_node_2);
  bool CompareEdges(const edge_t &lhs, const edge_t &rhs);
  void CheckEnergyRemesh(std::vector<score_list_t> *P_score_lists,
                         struct element_list_t *P_element_list,
                         struct all_node_list_t *P_node_list,
                         struct indentor_t *P_indentor, int i_load, int iQuasi);
  bool CompareScoreList(const score_list_t &lhs, const score_list_t &rhs);
  void DumpScoreList(std::vector<score_list_t> score_list);

private:
  static CreateMesh *_instance;
  double d_remeshTolerance;
};
}

#endif // CREATEMESH_H
