//
// Node.h
//

#if !defined(NODE_H)
#define NODE_H

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
class Node {

  //
  // public methods
  //

public:
  /**
   * @brief getInstance.
   */
  static Node *getInstance();

  /**
   * @brief destroyInstance.
   */
  static void destroyInstance();

  // ***************************************//
  // I = boundingbox.c

  //
  //  getBoundingBox()
  //
  void getBoundingBox(struct boundingbox_t *P_boundingbox,
                      struct node_list_t *P_node_list,
                      enum mt_version_t mt_version);

  // ***************************************//
  // I = clean_nodes.c

  //
  //  cleanNodes()
  //
  void cleanNodes(struct all_node_list_t *P_node_list);

  // ***************************************//
  // I = node.c

  //
  //	makeActiveNode()
  //
  void makeActiveNode(struct node_list_t *P_node_list, const int i_node);

  //
  //	makeInactiveNode()
  //
  void makeInactiveNode(struct node_list_t *P_node_list, const int i_node);

  // ***************************************//
  // I = node_clusters.c

  //
  //	isNodeAtSurface()
  //
  int isNodeAtSurface(const struct node_t *P_node);

  //
  //	createNodalClusters()
  //
  void createNodalClusters(int iQuasi, struct node_list_t *P_node_list,
                           struct lattice_t *P_lattice, int *cluster_size,
                           enum mt_version_t mt_version);

  // ***************************************//
  // I = node_clusters_weight.c

  //
  //	nodalClusterWeightTest()
  //
  int nodalClusterWeightTest(const int iQuasi,
                             const struct node_list_t *P_node_list,
                             struct lattice_t *P_lattice);

  //
  //	computeNodalClusterSitesAndWeights()
  //
  void computeNodalClusterSitesAndWeights(int iQuasi,
                                          struct lattice_t *P_lattice,
                                          struct node_list_t *P_node_list,
                                          double r_cluster,
                                          enum mt_version_t mt_version);

  // ***************************************//
  // I = function constructed from triangulation_client.c

  //
  //  checkForDuplicateNodes()
  //
  void checkForDuplicateNodes(struct node_list_t *P_node_list);

  // ***************************************//
  // I = function of remesh_elements.c

  //
  //  addNewNode()
  //
  void addNewNode(const int l[3], struct all_node_list_t *P_node_list,
                  struct lattice_t *P_lattice, const int iQuasi);

  //
  //  isNodeAtomistic()
  //
  bool isNodeAtomistic(const struct node_t *P_node,
                       const double atomistic_element_size);

  //
  // private methods
  //

private:
  /**
   * @brief Constructor.
   */
  Node();

  /**
   * @brief Copy constructor.
   */
  Node(Node const &);

  /**
   * @brief Assignment operator.
   */
  const Node &operator=(const Node &);

  /**
   * @brief Destructor.
   */
  ~Node();

  //
  //	BuildNodalClusters()
  //
  void BuildNodalClusters(int iQuasi, struct node_list_t *P_node_list,
                          struct lattice_t *P_element_list, int *cluster_size,
                          enum mt_version_t mt_version);

  //
  //	ResolveOverlaps()
  //
  void ResolveOverlaps(int iQuasi, struct node_list_t *P_node_list,
                       struct lattice_t *P_lattice, int *cluster_size,
                       enum mt_version_t mt_version);

  //
  //	CheckClusters()
  //
  void CheckClusters(const struct node_list_t *P_node_list,
                     const struct lattice_t *P_lattice);

  //
  //	ComputeLatticeExactSum()
  //
  double ComputeLatticeExactSum(const int iQuasi, struct lattice_t *P_lattice);

  //
  //	ComputeLatticeApproxSum()
  //
  double ComputeLatticeApproxSum(const int iQuasi, struct lattice_t *P_lattice,
                                 const struct node_list_t *P_node_list);

  //
  //	ComputeClusterNodalWeightsLamped()
  //
  void ComputeClusterNodalWeightsLamped(int iQuasi,
                                        struct node_list_t *P_node_list,
                                        struct lattice_t *P_lattice);

  //
  // private data types
  //
private:
  static Node *_instance;
};
}

#endif // NODE_H
