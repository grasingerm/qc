//
// Node.cc
//

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <alloca.h>
#include <search.h>
#include <math.h>
#include <unistd.h>
#include <assert.h>

#include <errno.h>

#include "C_Interface.h"
#include "DataTypes.h"
#include "Element.h"
#include "Error.h"
#include "Lattice.h"
#include "MiscFunctions.h"
#include "Node.h"
#include "PairPotentials.h"
#include "Quasicontinua.h"
#include "Shape.h"

#include "monitor.h"
#include "numa.h"
#include "threads.h"

// used in getBoundingBox() function.
static struct quasicontinuum::boundingbox_t *thread_bbox =
    (struct quasicontinuum::boundingbox_t *)NULL;

// globale variable used in makeActiveNode()
#define N_ALLOC_PAGES 2
static long int page_size;
static int node_bucket;

// global variables used in createNodalClusters()
#define SITE_ALLOC_CHUNK 10
#define CLUSTER_BUCKET 5
static int current_cluster = 0;
static pthread_mutex_t resolve_overlaps_lock = PTHREAD_MUTEX_INITIALIZER;

// global variable used in computeNodalClusterSitesAndWeights()
#define CLUSTER_SIZE_INCR 1.5
#define CLUSTER_INCR_MAX 100
#define CLUSTER_SIZE 3
#define CLUSTER_SURFACE_SIZE 3

// used in ComputeClusterNodalWeightsLampedWorker()
#define LAMPED_NODE_BUCKET 2
#define NODE_BUCKET 2
#define EXPLICIT_SHAPE_SUM_NUMBER_SITES 500
static int lamped_current_node = 0;
static pthread_mutex_t lamped_lock = PTHREAD_MUTEX_INITIALIZER;

//
//
//

namespace quasicontinuum {

//
//  thread block
//

// namesapce for getBoundingBox(), makeActiveNode(), makeInactiveNode()
namespace {

//
// data type for getBoundingBoxWorker()
//
struct getBoundingBoxData_t {
  struct boundingbox_t *thread_bbox;
  struct node_list_t *P_node_list;
};

//
//  getBoundingBoxProcessingLocal(_
//
void getBoundingBoxProcessingLocal(struct boundingbox_t *P_boundingbox,
                                   const struct node_list_t *P_node_list,
                                   const int i_start, const int i_end) {
  struct node_t *P_node;

  int i_node;

  // initialize
  P_node = P_node_list->nodes[i_start];

  memcpy(P_boundingbox->vertex_1, P_node->position, sizeof(double[3]));
  memcpy(P_boundingbox->vertex_2, P_node->position, sizeof(double[3]));

  // scan all nodes
  for (i_node = i_start; i_node <= i_end; i_node++) {
    P_node = P_node_list->nodes[i_node];

    P_boundingbox->vertex_1[0] =
        MIN(P_boundingbox->vertex_1[0], P_node->position[0]);
    P_boundingbox->vertex_1[1] =
        MIN(P_boundingbox->vertex_1[1], P_node->position[1]);
    P_boundingbox->vertex_1[2] =
        MIN(P_boundingbox->vertex_1[2], P_node->position[2]);

    P_boundingbox->vertex_2[0] =
        MAX(P_boundingbox->vertex_2[0], P_node->position[0]);
    P_boundingbox->vertex_2[1] =
        MAX(P_boundingbox->vertex_2[1], P_node->position[1]);
    P_boundingbox->vertex_2[2] =
        MAX(P_boundingbox->vertex_2[2], P_node->position[2]);
  }

  return;
} // end of getBoundingBoxProcessingLocal(),

//
//  getBoundingBoxWorker()
//
void *getBoundingBoxWorker(void *arg) {
  struct boundingbox_t *thread_bbox =
      ((struct getBoundingBoxData_t *)arg)->thread_bbox;
  struct node_list_t *P_node_list =
      ((struct getBoundingBoxData_t *)arg)->P_node_list;

  int number_threads = get_number_threads();
  int my_id = get_my_tid();
  int number_node_thread;
  int number_node_start;
  int number_node_end;

  /**
   * get my share of data
   */

  get_share(my_id, number_threads, P_node_list->number_nodes,
            &number_node_thread, &number_node_start, &number_node_end);

  /**
   * process
   */

  getBoundingBoxProcessingLocal(&(thread_bbox[my_id]), P_node_list,
                                number_node_start, number_node_end);

  return (NULL);
} // end of getBoundingBoxWorker()

} // end of namespace for getBoudingBox(),makeActiveNode()

// namespace for createNodalClusters()
namespace {
struct BuildNodalClustersData_t {
  int iQuasi;
  struct node_list_t *P_node_list;
  struct lattice_t *P_lattice;
  int *cluster_size;
  Lattice *latticeClass;
};

struct ResolveOverlapsData_t {
  int iQuasi;
  struct node_list_t *P_node_list;
  struct lattice_t *P_lattice;
  int *cluster_size;
};

struct ComputeClusterNodalWeightsLampedData_t {
  int iQuasi;
  struct node_list_t *P_node_list;
  struct lattice_t *P_lattice;
};

//
//	WeightTestFunctionLocal()
//
double WeightTestFunctionLocal(const double position[3]) {
  static double a;
  static double b;
  static double c;
  static double d;

  static int init_done = 0;

  if (init_done == 0) {

    a = 10 * drand48();
    b = 10 * drand48();
    c = drand48();
    d = drand48();
    a = 1.0;
    b = c = d = 0.0;
    init_done = 1;
  }

  return (a + b * position[0] + c * position[1] + d * position[2]);
} // end of WeightTestFunctionLocal()

//
//	clusterAddSiteLocal()
//
void clusterAddSiteLocal(struct node_t *P_node, int l[3]) {
  int i_site;

  // acquire lock
  pthread_mutex_lock(&P_node->lock);

  i_site = P_node->site_cluster.number_sites + 1;

  // check allocated space
  if (P_node->site_cluster.alloc == 0) {

    if ((P_node->site_cluster.sites =
             (int(*)[3])malloc(SITE_ALLOC_CHUNK * sizeof(int[3]))) == NULL) {
      PERROR_MALLOC();
      ERROR("malloc()");
      exit(EXIT_FAILURE);
    }

    P_node->site_cluster.alloc = SITE_ALLOC_CHUNK;
  } else if (i_site > P_node->site_cluster.alloc) {
    P_node->site_cluster.alloc += SITE_ALLOC_CHUNK;

    P_node->site_cluster.sites =
        (int(*)[3])realloc(P_node->site_cluster.sites,
                           sizeof(int[3]) * P_node->site_cluster.alloc);

    if (P_node->site_cluster.sites == NULL) {
      PERROR_MALLOC();
      ERROR("realloc()");
      exit(EXIT_FAILURE);
    }
  }

#if defined(__QC_SGI)
  change_memory_placement(P_node->site_cluster.sites,
                          sizeof(int[3]) * P_node->site_cluster.alloc, 0);
#endif /* sgi */

  // insert site
  P_node->site_cluster.number_sites++;

  i_site--;
  P_node->site_cluster.sites[i_site][0] = l[0];
  P_node->site_cluster.sites[i_site][1] = l[1];
  P_node->site_cluster.sites[i_site][2] = l[2];

  pthread_mutex_unlock(&P_node->lock);

  return;
} // end of clusterAddSiteLocal()

//
//	BuildNodalClustersProcessingLocal()
//
void BuildNodalClustersProcessingLocal(const int iQuasi,
                                       struct node_list_t *P_node_list,
                                       struct lattice_t *P_lattice,
                                       const int *cluster_size, int node_start,
                                       int node_end) {
  struct shell_t *shell;

  int l[3];
  int i_node;
  int shell_number;

  // loop over all nodes and add all sites that are within r_cluster
  for (i_node = node_start; i_node <= node_end; i_node++) {
    struct node_t *P_node = P_node_list->nodes[i_node];
    int local_cluster_size = cluster_size[i_node];

    // add node itself
    clusterAddSiteLocal(P_node, P_node->l);

    // process shells between 1 and local_cluster_size
    for (shell_number = 1; shell_number <= local_cluster_size; shell_number++) {
      int i_site;
      int i_elem;

      // get shell
      shell = Lattice::getInstance()->getShell(shell_number, P_lattice);

      // process all sites in a current shell
      for (i_site = 0; i_site < shell->number_sites; i_site++) {
        int in_cluster = 0;

        l[0] = P_node->l[0] + shell->site[i_site][0];
        l[1] = P_node->l[1] + shell->site[i_site][1];
        l[2] = P_node->l[2] + shell->site[i_site][2];

        // check if lattice site is inside the lattice
        if (Lattice::getInstance()->isSiteInsideLattice(P_lattice, l, iQuasi) ==
            RETURN_FAILURE)
          continue;

        // check if the site belongs to the element cluster at the node
        for (i_elem = 0; i_elem < P_node->element_list.number_elements;
             i_elem++) {
          struct node_t *P_node_0 =
              P_node->element_list.elements[i_elem]->node[0];
          struct node_t *P_node_1 =
              P_node->element_list.elements[i_elem]->node[1];
          struct node_t *P_node_2 =
              P_node->element_list.elements[i_elem]->node[2];
          struct node_t *P_node_3 =
              P_node->element_list.elements[i_elem]->node[3];

          if (MiscFunctions::getInstance()->checkSiteInTetra(
                  P_node_0->initial_position, P_node_1->initial_position,
                  P_node_2->initial_position, P_node_3->initial_position, l,
                  P_lattice, iQuasi) == INSIDE) {
            in_cluster = 1;
            break;
          }
        }

        // site does not belong to any of the elements in the cluster;
        // bail out
        if (in_cluster == 0)
          continue;

        // got this far, add site
        clusterAddSiteLocal(P_node, l);
      }
    }
  }

  return;
} // end of BuildNodalClustersProcessingLocal()

//
//	BuildNodalClustersWorker()
//
void *BuildNodalClustersWorker(void *arg) {
  const int iQuasi = ((struct BuildNodalClustersData_t *)arg)->iQuasi;
  struct node_list_t *P_node_list =
      ((struct BuildNodalClustersData_t *)arg)->P_node_list;
  struct lattice_t *P_lattice =
      ((struct BuildNodalClustersData_t *)arg)->P_lattice;
  const int *cluster_size =
      ((struct BuildNodalClustersData_t *)arg)->cluster_size;

  int number_threads = get_number_threads();
  int my_id = get_my_tid();
  int number_node_thread;
  int number_node_start;
  int number_node_end;

  /**
   * get my share of elements
   */

  get_share(my_id, number_threads, P_node_list->number_nodes,
            &number_node_thread, &number_node_start, &number_node_end);

  /**
   *  process elements
   */

  BuildNodalClustersProcessingLocal(iQuasi, P_node_list, P_lattice,
                                    cluster_size, number_node_start,
                                    number_node_end);

  return (NULL);
} // end of BuildNodalClustersWorker()

//
//	RemoveOverlapsLocal()
//
void RemoveOverlapsLocal(struct node_t *P_node_i, struct node_t *P_node_j) {

  pthread_mutex_t *P_lock1;
  pthread_mutex_t *P_lock2;

  int(*i_sites)[3];
  int(*j_sites)[3];

  int *mask_i;
  int *mask_j;

  int r_i[3];
  int r_j[3];

  int i_site;
  int j_site;

  int i_number_sites = 0;
  int j_number_sites = 0;

  // to avoid deadlock mutexes have to be aquired in order; here the
  // order will be set by the node number
  if (P_node_i->number < P_node_j->number) {
    P_lock1 = &(P_node_i->lock);
    P_lock2 = &(P_node_j->lock);
  } else {
    P_lock1 = &(P_node_j->lock);
    P_lock2 = &(P_node_i->lock);
  }

  // lock both locks in the right order
  pthread_mutex_lock(P_lock1);
  pthread_mutex_lock(P_lock2);

  // allocate space for i_sites and j_sites
  i_sites =
      (int(*)[3])alloca(P_node_i->site_cluster.number_sites * sizeof(int[3]));

  memcpy(&(i_sites[0]), &(P_node_i->site_cluster.sites[0]),
         P_node_i->site_cluster.number_sites * sizeof(int[3]));

  j_sites =
      (int(*)[3])alloca(P_node_j->site_cluster.number_sites * sizeof(int[3]));

  memcpy(&(j_sites[0]), &(P_node_j->site_cluster.sites[0]),
         P_node_j->site_cluster.number_sites * sizeof(int[3]));

  // allocate space for mask_i and mask_j
  mask_i = (int *)alloca(P_node_i->site_cluster.number_sites * sizeof(int));

  memset(mask_i, '\0', P_node_i->site_cluster.number_sites * sizeof(int));

  mask_j = (int *)alloca(P_node_j->site_cluster.number_sites * sizeof(int));

  memset(mask_j, '\0', P_node_j->site_cluster.number_sites * sizeof(int));

  // scan both nodes for same sites; simple linear search here
  for (i_site = 0; i_site < P_node_i->site_cluster.number_sites; i_site++)
    for (j_site = 0; j_site < P_node_j->site_cluster.number_sites; j_site++)
      if ((P_node_i->site_cluster.sites[i_site][0] ==
           P_node_j->site_cluster.sites[j_site][0]) &&
          (P_node_i->site_cluster.sites[i_site][1] ==
           P_node_j->site_cluster.sites[j_site][1]) &&
          (P_node_i->site_cluster.sites[i_site][2] ==
           P_node_j->site_cluster.sites[j_site][2])) { /* same site found */

        // compute distance between both nodes and the site; use
        // lattice cooridinates
        r_i[0] = P_node_i->l[0] - P_node_i->site_cluster.sites[i_site][0];
        r_i[1] = P_node_i->l[1] - P_node_i->site_cluster.sites[i_site][1];
        r_i[2] = P_node_i->l[2] - P_node_i->site_cluster.sites[i_site][2];

        r_j[0] = P_node_j->l[0] - P_node_i->site_cluster.sites[i_site][0];
        r_j[1] = P_node_j->l[1] - P_node_i->site_cluster.sites[i_site][1];
        r_j[2] = P_node_j->l[2] - P_node_i->site_cluster.sites[i_site][2];

        // flag the node that the site should be removed from
        if (r_i[0] * r_i[0] + r_i[1] * r_i[1] + r_i[2] * r_i[2] <
            r_j[0] * r_j[0] + r_j[1] * r_j[1] + r_j[2] * r_j[2]) {
          // P_node_i is closer
          mask_j[j_site] = 1;
        } else {
          // P_node_j is closer
          mask_i[i_site] = 1;
        }
      }

  // remove all sites that have been flagged
  for (i_site = 0; i_site < P_node_i->site_cluster.number_sites; i_site++)
    if (!mask_i[i_site]) {
      P_node_i->site_cluster.sites[i_number_sites][0] = i_sites[i_site][0];
      P_node_i->site_cluster.sites[i_number_sites][1] = i_sites[i_site][1];
      P_node_i->site_cluster.sites[i_number_sites][2] = i_sites[i_site][2];

      i_number_sites++;
    }

  P_node_i->site_cluster.number_sites = i_number_sites;

  for (j_site = 0; j_site < P_node_j->site_cluster.number_sites; j_site++)
    if (!mask_j[j_site]) {
      P_node_j->site_cluster.sites[j_number_sites][0] = j_sites[j_site][0];
      P_node_j->site_cluster.sites[j_number_sites][1] = j_sites[j_site][1];
      P_node_j->site_cluster.sites[j_number_sites][2] = j_sites[j_site][2];

      j_number_sites++;
    }

  P_node_j->site_cluster.number_sites = j_number_sites;

  // unlock locks
  pthread_mutex_unlock(P_lock2);
  pthread_mutex_unlock(P_lock1);

  // cleanup
  return;
} // end of RemoveOverlapsLocal()

//
//	ResolveOverlapsProcessingLocal()
//
void ResolveOverlapsProcessingLocal(const int iQuasi,
                                    struct node_list_t *P_node_list,
                                    struct lattice_t *P_lattice,
                                    const int *cluster_size,
                                    const int node_start, const int node_end) {
  struct shell_t *P_shell;

  double r[3];
  double r_cluster_i;
  double r_cluster_j;

  int i_node;
  int j_node;

  Lattice *latticeC = Lattice::getInstance();

  // loop over all nodes
  for (i_node = node_start; i_node <= node_end; i_node++) {
    struct node_t *P_node_i = P_node_list->nodes[i_node];

    P_shell = latticeC->getShell(cluster_size[i_node], P_lattice);

    r_cluster_i = latticeC->getShellRadius(P_lattice, P_shell->site[0], iQuasi);

    r_cluster_i = sqrt(r_cluster_i);

    for (j_node = 0; j_node < P_node_list->number_nodes; j_node++) {
      struct node_t *P_node_j = P_node_list->nodes[j_node];

      if (i_node == j_node)
        continue;

      P_shell = latticeC->getShell(cluster_size[j_node], P_lattice);
      r_cluster_j =
          latticeC->getShellRadius(P_lattice, P_shell->site[0], iQuasi);
      r_cluster_j = sqrt(r_cluster_j);

      // check the distance between nodes; if it is larger then
      // 2*r_cluster there is no chance for clusters to overlap
      r[0] = P_node_i->initial_position[0] - P_node_j->initial_position[0];
      r[1] = P_node_i->initial_position[1] - P_node_j->initial_position[1];
      r[2] = P_node_i->initial_position[2] - P_node_j->initial_position[2];

      if (r[0] * r[0] + r[1] * r[1] + r[2] * r[2] >
          (r_cluster_i + r_cluster_j) * (r_cluster_i + r_cluster_j))
        continue;

      // clusters overlap--resolve conflicts
      RemoveOverlapsLocal(P_node_i, P_node_j);
    }
  }

  return;
} // end of ResolveOverlapsProcessingLocal()

//
//	ResolveOverlapsWorker()
//
void *ResolveOverlapsWorker(void *arg) {
  const int iQuasi = ((struct ResolveOverlapsData_t *)arg)->iQuasi;
  struct node_list_t *P_node_list =
      ((struct ResolveOverlapsData_t *)arg)->P_node_list;
  struct lattice_t *P_lattice =
      ((struct ResolveOverlapsData_t *)arg)->P_lattice;
  const int *cluster_size = ((struct ResolveOverlapsData_t *)arg)->cluster_size;

  int i_node;
  int j_node;
  int node_start;
  int node_end;

  // loop until work

  /* CONSTCOND */
  while (1) {
    pthread_mutex_lock(&resolve_overlaps_lock);

    if (current_cluster == P_node_list->number_nodes) {
      pthread_mutex_unlock(&resolve_overlaps_lock);
      break;
    }

    node_start = current_cluster;

    // update current_cluster
    current_cluster += CLUSTER_BUCKET;

    if (current_cluster >= P_node_list->number_nodes)
      current_cluster = P_node_list->number_nodes;

    node_end = current_cluster - 1;

    pthread_mutex_unlock(&resolve_overlaps_lock);

    // do the work
    ResolveOverlapsProcessingLocal(iQuasi, P_node_list, P_lattice, cluster_size,
                                   node_start, node_end);
  }

  return (NULL);
} // end of ResolveOverlapsWorker()

//
//	ComputeElementClusterRadiusLocal()
//
//	compute radius of the cluster; take into account elements that
//  contain less than EXPLICIT_SHAPE_SUM_NUMBER_SITES sites inside
//  the element;
//
double ComputeElementClusterRadiusLocal(const struct node_t *P_node) {
  double r[3];
  double dr;
  double radius = 0.0;

  int i_elem;

  // loop over all elements attached to node
  for (i_elem = 0; i_elem < P_node->element_list.number_elements; i_elem++) {
    const struct element_t *P_element = P_node->element_list.elements[i_elem];
    int j_node;

    // check if the element is included in the explicit sum
    if (P_element->number_sites > EXPLICIT_SHAPE_SUM_NUMBER_SITES)
      continue;

    // process nodes attached to the element
    for (j_node = 0; j_node < 4; j_node++) {
      const struct node_t *P_node_j = P_element->node[j_node];

      r[0] = P_node->initial_position[0] - P_node_j->initial_position[0];
      r[1] = P_node->initial_position[1] - P_node_j->initial_position[1];
      r[2] = P_node->initial_position[2] - P_node_j->initial_position[2];

      dr = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];

      if (dr > radius)
        radius = dr;
    }
  }

  return (radius);
} // end of ComputeElementClusterRadiusLocal()

//
//	ComputeLargeElements_ShapeFunctionSumLocal()
//
//	process all elements attached to the node and include their
//  contribution to the shape function sum
//
double ComputeLargeElements_ShapeFunctionSumLocal(const struct node_t *P_node) {
  double sum = 0.0;

  int i_elem;

  for (i_elem = 0; i_elem < P_node->element_list.number_elements; i_elem++) {
    const struct element_t *P_element = P_node->element_list.elements[i_elem];

    if (P_element->number_sites > EXPLICIT_SHAPE_SUM_NUMBER_SITES)
      sum += P_element->number_sites / 4.0;
  }

  return (sum);
} // end of ComputeLargeElements_ShapeFunctionSumLocal()

//
//	ComputeShapeFunctionSumLocal()
//
//	for a given node compute the value of the sum of shape function
//  anchored at the node over the domain
//
double ComputeShapeFunctionSumLocal(const int iQuasi,
                                    const struct node_t *P_node,
                                    struct lattice_t *P_lattice) {
  // initialize sum to 1 to account for the site that is occupied by
  // the node
  double sum = 1.0;
  double element_cluster_radius = 0.0;

  int shell_number = 1;
  int break_flag = 1;
  int i_node;

  // compute element cluster radius
  element_cluster_radius = ComputeElementClusterRadiusLocal(P_node);

  // compute contribution from large elements
  sum += ComputeLargeElements_ShapeFunctionSumLocal(P_node);

  /**
   * the general idea is as follows:
   * start from a node and proceed on the shell-by-shell basis to
   * identify lattice sites; each of the lattice sites should belong
   * to an element that is attached to the node;
   * the calculation is stopped if the radius of the
   * current shell is larger than the max distance between the central
   * node and all other nodes in its element cluster
   */
  do {
    Lattice *latticeC = Lattice::getInstance();

    struct shell_t *P_shell = latticeC->getShell(shell_number, P_lattice);
    int l[3];
    int i_site;

    // compute radius of the current shell and check if it is smaller
    // than element_cluster_radius
    if (latticeC->getShellRadius(P_lattice, P_shell->site[0], iQuasi) >
        element_cluster_radius)
      break;

    // process all sites in the current shell
    for (i_site = 0; i_site < P_shell->number_sites; i_site++) {
      int i_elem;
      l[0] = P_node->l[0] + P_shell->site[i_site][0];
      l[1] = P_node->l[1] + P_shell->site[i_site][1];
      l[2] = P_node->l[2] + P_shell->site[i_site][2];

      // compute shape function for the site using all elements attached
      // to the node
      for (i_elem = 0; i_elem < P_node->element_list.number_elements;
           i_elem++) {
        struct element_t *P_element = P_node->element_list.elements[i_elem];
        struct node_t *P_nodes[3];
        int j_node = 0;

        // skip large elements for which contribution is computed
        // approximately
        if (P_element->number_sites > EXPLICIT_SHAPE_SUM_NUMBER_SITES)
          continue;

        // check if the site is located inside an element
        if (MiscFunctions::getInstance()->checkSiteInTetra(
                P_element->node[0]->initial_position,
                P_element->node[1]->initial_position,
                P_element->node[2]->initial_position,
                P_element->node[3]->initial_position, l, P_lattice,
                iQuasi) == OUTSIDE)
          continue;

        for (i_node = 0; i_node < 4; i_node++)
          if (P_element->node[i_node] != P_node)
            P_nodes[j_node++] = P_element->node[i_node];

        // compute shape function at current site given an element;
        // make sure that it is anchored at the current node (ie. P_node)
        sum += Shape::getInstance()->computeSiteShapeFunction(
            P_node, P_nodes[0], P_nodes[1], P_nodes[2], l, P_lattice, iQuasi);

        // if we got this far than the element for the site has been
        // located;

        break;
      }
    }

    // switch to the next shell
    shell_number++;
  } while (break_flag);

  return (sum);
} // end of ComputeShapeFunctionSumLocal()

//
//	ComputeClusterShapeSumSelfLocal()
//
double ComputeClusterShapeSumSelfLocal(const int iQuasi,
                                       const struct node_t *P_node,
                                       struct lattice_t *P_lattice) {
  double sum = 0.0;

  int i_site;

  // cluster sites all belong to the cluster of elements around the node
  // for each site in the cluster get the shape function
  for (i_site = 0; i_site < P_node->site_cluster.number_sites; i_site++) {
    struct node_t *P_nodes[3];
    struct element_t *P_element;

    int *l = P_node->site_cluster.sites[i_site];
    int i_elem;
    int i_node;
    int j_node = 0;

    // locate the element the site is in
    for (i_elem = 0; i_elem < P_node->element_list.number_elements; i_elem++) {
      P_element = P_node->element_list.elements[i_elem];
      if (MiscFunctions::getInstance()->checkSiteInTetra(
              P_element->node[0]->initial_position,
              P_element->node[1]->initial_position,
              P_element->node[2]->initial_position,
              P_element->node[3]->initial_position, l, P_lattice,
              iQuasi) == INSIDE)
        break; /* element found */
    }

// check if the site has been found in the cluster of elements that
// surrounds the node
#if defined(_QC_NODE_CLUSTER_WEIGHT_DEBUG_)

    assert(i_elem < P_node->element_list.number_elements);

#endif /* _QC_NODE_CLUSTER_WEIGHT_DEBUG_ */

    // find three other nodes on the element list
    for (i_node = 0; i_node < 4; i_node++)
      if (P_element->node[i_node] != P_node)
        P_nodes[j_node++] = P_element->node[i_node];

    // compute shape function and add it to the sum
    sum += Shape::getInstance()->computeSiteShapeFunction(
        P_node, P_nodes[0], P_nodes[1], P_nodes[2], l, P_lattice, iQuasi);
  }

  return (sum);
} // end of ComputeClusterShapeSumSelfLocal()

//
//	ComputeClusterShapeSumOtherLocal()
//
double ComputeClusterShapeSumOtherLocal(const int iQuasi,
                                        const struct node_t *P_cluster_node,
                                        const struct node_t *P_node,
                                        struct lattice_t *P_lattice) {
  double sum = 0.0;

  int i_site;

  // some error checks might be inserted:
  // 1) check if P_node_i and P_node_j belong to a common element
  // 2) check if P_node_i != P_node_j
  for (i_site = 0; i_site < P_cluster_node->site_cluster.number_sites;
       i_site++) {
    struct node_t *P_nodes[3];
    struct element_t *P_element;

    int *l = P_cluster_node->site_cluster.sites[i_site];
    int i_elem;
    int i_node;
    int j_node = 0;

    // locate the element the site is in
    for (i_elem = 0; i_elem < P_node->element_list.number_elements; i_elem++) {
      P_element = P_node->element_list.elements[i_elem];
      if (MiscFunctions::getInstance()->checkSiteInTetra(
              P_element->node[0]->initial_position,
              P_element->node[1]->initial_position,
              P_element->node[2]->initial_position,
              P_element->node[3]->initial_position, l, P_lattice,
              iQuasi) == INSIDE)
        break; /* element found */
    }

    // check if the element has been found
    if (i_elem == P_node->element_list.number_elements)
      continue;

    // find three other nodes on the element list
    for (i_node = 0; i_node < 4; i_node++)
      if (P_element->node[i_node] != P_node)
        P_nodes[j_node++] = P_element->node[i_node];

    // compute shape function and add it to the sum
    sum += Shape::getInstance()->computeSiteShapeFunction(
        P_node, P_nodes[0], P_nodes[1], P_nodes[2], l, P_lattice, iQuasi);
  }

  return (sum);
} // end of ComputeClusterShapeSumOtherLocal()

//
//	ComputeClusterNodalWeightsLampedWorker()
//
void *ComputeClusterNodalWeightsLampedWorker(void *arg) {
  const int iQuasi =
      ((struct ComputeClusterNodalWeightsLampedData_t *)arg)->iQuasi;
  struct node_list_t *P_node_list =
      ((struct ComputeClusterNodalWeightsLampedData_t *)arg)->P_node_list;
  struct lattice_t *P_lattice =
      ((struct ComputeClusterNodalWeightsLampedData_t *)arg)->P_lattice;

  int node_start;
  int node_end;
  int i_node;
  volatile int spin = 1;

  while (spin) {
    pthread_mutex_lock(&lamped_lock);

    // check the value of current_atom; if its number_atoms-1
    // we are done here
    if (lamped_current_node == P_node_list->number_nodes) {
      pthread_mutex_unlock(&lamped_lock);
      break;
    }

    node_start = lamped_current_node;

    // update current_node
    lamped_current_node += NODE_BUCKET;

    if (lamped_current_node >= P_node_list->number_nodes)
      lamped_current_node = P_node_list->number_nodes;

    node_end = lamped_current_node - 1;

    pthread_mutex_unlock(&lamped_lock);

    // process node
    for (i_node = node_start; i_node <= node_end; i_node++) {
      struct node_t *P_node_i = P_node_list->nodes[i_node];
      double cluster_sum;
      int j_node;

      // add self
      cluster_sum =
          ComputeClusterShapeSumSelfLocal(iQuasi, P_node_i, P_lattice);

      for (j_node = 0; j_node < P_node_i->node_list.number_nodes; j_node++) {
        struct node_t *P_node_j = P_node_i->node_list.nodes[j_node];

        cluster_sum += ComputeClusterShapeSumOtherLocal(iQuasi, P_node_j,
                                                        P_node_i, P_lattice);
      }

      // compute exact sum and weight
      P_node_i->weight =
          ComputeShapeFunctionSumLocal(iQuasi, P_node_i, P_lattice) /
          cluster_sum;
    }
  }

  return (NULL);
} // end of ComputeClusterNodalWeightsLampedWorker()

//
//	checkWeightsLocal()
//
int checkWeightsLocal(const struct node_list_t *P_node_list) {
  int i_node;

  for (i_node = 0; i_node < P_node_list->number_nodes; i_node++)
    if (P_node_list->nodes[i_node]->weight < 0.0)
      return (1);

  return (0);
} // end of checkWeightsLocal()

//
//	negativeWeightInfoLocal()
//
void negativeWeightInfoLocal(const struct node_list_t *P_node_list,
                             struct lattice_t *P_lattice, const int iQuasi) {
  int i_node;

  for (i_node = 0; i_node < P_node_list->number_nodes; i_node++) {
    struct node_t *P_node = P_node_list->nodes[i_node];

    if (P_node->weight < 1.0) {
      int i_site;
      int i_elem;

      d_print("Node: %d\n", P_node->number);

      // loop over all elements incident to node and check if
      // cluster sites are located in all of them
      for (i_elem = 0; i_elem < P_node->element_list.number_elements;
           i_elem++) {
        struct element_t *P_element = P_node->element_list.elements[i_elem];
        int i_hit = 0;

        // loop over all sites in the cluster; check their location
        // against the element
        for (i_site = 0; i_site < P_node->site_cluster.number_sites; i_site++)
          if (P_element ==
              Element::getInstance()->locateSiteElement(
                  NULL, P_node->site_cluster.sites[i_site], P_lattice, iQuasi))
            i_hit++;

        d_print("element %d contains %d sites (of %d)\n", i_elem, i_hit,
                (int)P_element->number_sites);
      }
    }
  }

  return;
} // end of negativeWeightInfoLocal()

//
//	zeroClustersLocal()
//
void zeroClustersLocal(struct node_list_t *P_node_list) {
  int i_node;

  for (i_node = 0; i_node < P_node_list->number_nodes; i_node++)
    P_node_list->nodes[i_node]->site_cluster.number_sites = 0;

  return;
} // end of zeroClustersLocal()
} // end of namespace for createNodalClusters()

//
//	namespace for checkForDuplicateNodes()
namespace {
struct checkForDuplicateNodesData_t {
  struct node_list_t *P_node_list;
};

//
//	getShareDuplicateNodesLocal()
//
void getShareDuplicateNodesLocal(const int my_id, const int number_threads,
                                 const int number_data, int *number_data_thread,
                                 int *number_data_start1, int *number_data_end1,
                                 int *number_data_start2,
                                 int *number_data_end2) {
  // check if number_data >= number_threads; if not signal en error
  // and abort
  assert(number_data >= number_threads);
  assert(my_id < number_threads);

  // node matrix is symmetric, only check lower half,
  // give each thread top portion and bottom portion for equal work

  *number_data_thread = number_data / (number_threads * 2);

  *number_data_start1 = my_id * (*number_data_thread);

  *number_data_end1 = *number_data_start1 + *number_data_thread;

  *number_data_end2 = number_data - my_id * (*number_data_thread);

  // last thread
  if (my_id == number_threads - 1) {
    *number_data_start2 = *number_data_end1;
  } else {
    *number_data_start2 = *number_data_end2 - *number_data_thread;
  }

  return;
} // end of getShareDuplicateNodesLocal()

//
//	checkForDuplicateNodesProcessingLocal()
//
void checkForDuplicateNodesProcessingLocal(
    const struct node_list_t *P_node_list, int start_node, int end_node) {
  int i_node;
  int j_node;

  // dont need to check last node
  if (end_node == P_node_list->number_nodes) {
    end_node--;
  }

  // check position of given nodes against nodes below themselves
  for (i_node = start_node; i_node < end_node; i_node++) {
    struct node_t *P_node_i = P_node_list->nodes[i_node];

    for (j_node = i_node + 1; j_node < P_node_list->number_nodes; j_node++) {
      struct node_t *P_node_j = P_node_list->nodes[j_node];

      // check position and output if same
      if ((P_node_i->l[0] == P_node_j->l[0]) &&
          (P_node_i->l[1] == P_node_j->l[1]) &&
          (P_node_i->l[2] == P_node_j->l[2])) {
        d_print("nodes: %d and %d are identical\n", i_node, j_node);
        d_print("node %d: % e % e % e\n", i_node, P_node_i->initial_position[0],
                P_node_i->initial_position[1], P_node_i->initial_position[2]);
        d_print("node %d: % e % e % e\n", j_node, P_node_j->initial_position[0],
                P_node_j->initial_position[1], P_node_j->initial_position[2]);
      }
    }
  }

  return;
} // end of checkForDuplicateNodesProcessingLocal()

//
//	checkForDuplicateNodesWorker()
//
void *checkForDuplicateNodesWorker(void *arg) {
  const struct node_list_t *P_node_list =
      ((struct checkForDuplicateNodesData_t *)arg)->P_node_list;

  int number_threads = get_number_threads();
  int my_id = get_my_tid();
  int number_data_thread;
  int number_start1;
  int number_end1;
  int number_start2;
  int number_end2;

  getShareDuplicateNodesLocal(my_id, number_threads, P_node_list->number_nodes,
                              &number_data_thread, &number_start1, &number_end1,
                              &number_start2, &number_end2);

  // check first block of nodes
  checkForDuplicateNodesProcessingLocal(P_node_list, number_start1,
                                        number_end1);

  // check second block of nodes
  checkForDuplicateNodesProcessingLocal(P_node_list, number_start2,
                                        number_end2);

  return (NULL);
} // end of checkForDuplicateNodesWorker()
} // end of namespace for checkForDuplicateNodes()

//
//  initializing static data
//
Node *Node::_instance = NULL;

//
// constructor
//

Node::Node() {

  //
  //
  //
  return;
}

//
// destructor
//

Node::~Node() {

  //
  //
  //
  return;
}

//
// getInstance method
//

Node *Node::getInstance() {

  //
  // if not created, create
  //
  if (_instance == NULL) {
    _instance = new Node();
  }

  //
  // return instance
  //
  return _instance;
}

//
// destroy instance method
//

void Node::destroyInstance() {

  //
  // delete instance
  //
  delete _instance;

  //
  //
  //
  return;
}

// ***************************************//
// I = boundingbox.c

//
//  getBoundingBox()
//
void Node::getBoundingBox(struct boundingbox_t *P_boundingbox,
                          struct node_list_t *P_node_list,
                          enum mt_version_t mt_version) {
  struct getBoundingBoxData_t getBoundingBoxData;

  // static data thread_bbox is define inside thread block

  int number_threads;
  int i_thread;

  switch (mt_version) {
  case SINGLE_THREADED:
    getBoundingBoxProcessingLocal(P_boundingbox, P_node_list, 0,
                                  P_node_list->number_nodes - 1);

    break;

  case MULTI_THREADED:
    number_threads = get_max_number_threads();

    // check for space in thread_bbox
    if (thread_bbox == (struct boundingbox_t *)NULL) {
      thread_bbox = (struct boundingbox_t *)malloc(
          number_threads * sizeof(struct boundingbox_t));

      if (thread_bbox == NULL) {
        PERROR_MALLOC();
        ERROR("malloc()");
        exit(EXIT_FAILURE);
      }
    }

    getBoundingBoxData.thread_bbox = thread_bbox;
    getBoundingBoxData.P_node_list = P_node_list;

    thread_monitor(getBoundingBoxWorker, (void *)&getBoundingBoxData,
                   number_threads);

    // scan thread specific data and get P_boundingbox
    memcpy(P_boundingbox, &(thread_bbox[0]), sizeof(struct boundingbox_t));

    for (i_thread = 0; i_thread < number_threads; i_thread++) {
      P_boundingbox->vertex_1[0] =
          MIN(P_boundingbox->vertex_1[0], thread_bbox[i_thread].vertex_1[0]);

      P_boundingbox->vertex_1[1] =
          MIN(P_boundingbox->vertex_1[1], thread_bbox[i_thread].vertex_1[1]);

      P_boundingbox->vertex_1[2] =
          MIN(P_boundingbox->vertex_1[2], thread_bbox[i_thread].vertex_1[2]);

      P_boundingbox->vertex_2[0] =
          MAX(P_boundingbox->vertex_2[0], thread_bbox[i_thread].vertex_2[0]);

      P_boundingbox->vertex_2[1] =
          MAX(P_boundingbox->vertex_2[1], thread_bbox[i_thread].vertex_2[1]);

      P_boundingbox->vertex_2[2] =
          MAX(P_boundingbox->vertex_2[2], thread_bbox[i_thread].vertex_2[2]);
    }
    break;
  }

  return;
} // end of getBoundingBox()

// ***************************************//
// I = clean_nodes.c

//
//  cleanNodes()
//  Description : Reallocate space to include new nodes, add new nodes at
//                the end of continuous block.
//                This function must be called after all elements has been
//                cleaned; otherwise references to removed elements will
//                be present.
//
void Node::cleanNodes(struct all_node_list_t *P_node_list) {
  int number_nodes = P_node_list->node_list.number_nodes +
                     P_node_list->new_node_list.number_nodes;

  int i_node = P_node_list->node_list.number_nodes;
  int j_node;

  // put inactive nodes back on the list
  for (j_node = 0; j_node < P_node_list->node_list.number_inactive_nodes;
       j_node++)
    P_node_list->node_list.nodes[i_node + j_node] =
        P_node_list->node_list.inactive_nodes[j_node];

  P_node_list->node_list.number_inactive_nodes = 0;

  // allocate space for more nodes and copy contents
  for (j_node = 0; j_node < P_node_list->new_node_list.number_nodes; j_node++) {
    struct node_t *P_node;
    struct node_t *P_new_node;

    int k_node = i_node + j_node;

    makeActiveNode(&P_node_list->node_list, k_node);

    /**
     * copy contents of an appropriate new node; data copied:
     * -node number
     * -fix mask
     * -fix_w_mask
     * -initial position
     * -initial_temperature
     * -initial_frequency
     * -temperature
     * -frequency
     * -position
     * -velocity
     * -lattice site
     * -acceleration
     * -force_frequency
     */

    P_node = P_node_list->node_list.nodes[k_node];
    P_new_node = P_node_list->new_node_list.nodes[j_node];

    P_node->number = P_new_node->number;
    P_node->fix_mask = P_new_node->fix_mask;
    P_node->fix_w_mask = P_new_node->fix_w_mask;

    memcpy(P_node->initial_position, P_new_node->initial_position,
           sizeof(P_node->initial_position));

    // memcpy(P_node->initial_temperature, P_new_node->initial_temperature,
    // sizeof(double));
    P_node->initial_temperature = P_new_node->initial_temperature;

    // memcpy(P_node->initial_frequency, P_new_node->initial_frequency,
    // sizeof(double));
    P_node->initial_frequency = P_new_node->initial_frequency;

    P_node->initial_tau = P_new_node->initial_tau;

    memcpy(P_node->position, P_new_node->position, sizeof(P_node->position));

    // memcpy(P_node->temperature, P_new_node->temperature,
    // sizeof(double));
    P_node->temperature = P_new_node->temperature;

    // memcpy(P_node->frequency, P_new_node->frequency,
    // sizeof(double));
    P_node->frequency = P_new_node->frequency;

    P_node->tau = P_new_node->tau;

    memcpy(P_node->velocity, P_new_node->velocity, sizeof(P_node->velocity));

    memcpy(P_node->l, P_new_node->l, sizeof(P_node->l));

    memcpy(P_node->acceleration, P_new_node->acceleration,
           sizeof(P_node->acceleration));

    // memcpy(P_node->force_frequency, P_new_node->force_frequency,
    // 	sizeof(double));
    P_node->force_frequency = P_new_node->force_frequency;
  }

  // update number of nodes
  P_node_list->node_list.number_nodes = number_nodes;

  // clean new_nodes
  P_node_list->new_node_list.number_nodes = 0;

  return;
} // end of cleanNodes()

// ***************************************//
// I = node.c

//
//	makeActiveNode()
//
void Node::makeActiveNode(struct node_list_t *P_node_list, const int i_node) {
  static int init = 0;
  int j_node;

  if (init == 0) {
    // extern int errno;
    int old_errno;

    old_errno = errno;

    // for NUMA we really want to allocate more then a page aligned on
    // page boundary
    page_size = sysconf(_SC_PAGESIZE);
    if (page_size < 0 && errno != old_errno) {
      ERROR("sysconf()");
      exit(EXIT_FAILURE);
    }

    // compute node bucket; we want it to be as close to page_size as
    // possible
    node_bucket = N_ALLOC_PAGES * page_size / sizeof(struct node_t);

    if (node_bucket == 0)
      node_bucket = 1;

    init = 1;
  }

  // check inactive queue
  if (P_node_list->number_inactive_nodes > 0) {
    // check the buffer
    if (P_node_list->alloc_space < i_node + 1) {
      if (P_node_list->alloc_space == 0)
        P_node_list->nodes = (struct node_t **)malloc(sizeof(struct node_t *));
      else
        P_node_list->nodes = (struct node_t **)realloc(
            P_node_list->nodes,
            (P_node_list->alloc_space + 1) * sizeof(struct node_t *));

#if defined(__QC_SGI)
      change_memory_placement(P_node_list->nodes, P_node_list->alloc_space + 1,
                              0);
#endif /* sgi */

      if (P_node_list->nodes == NULL) {
        PERROR_MALLOC();
        ERROR("malloc()");
        exit(EXIT_FAILURE);
      }

      P_node_list->alloc_space++;
    }

    // insert
    j_node = --P_node_list->number_inactive_nodes;
    P_node_list->nodes[i_node] = P_node_list->inactive_nodes[j_node];

    return;
  }

  // there are no inactive nodes
  if (P_node_list->alloc_space < i_node + 1) {
    // there is not enough space in the arena-add some more
    if (P_node_list->alloc_space == 0)
      P_node_list->nodes =
          (struct node_t **)malloc(node_bucket * sizeof(struct node_t *));
    else
      P_node_list->nodes = (struct node_t **)realloc(
          P_node_list->nodes,
          (P_node_list->alloc_space + node_bucket) * sizeof(struct node_t *));

#if defined(__QC_SGI)
    change_memory_placement(P_node_list->nodes,
                            P_node_list->alloc_space + node_bucket, 0);
#endif /* sgi */

    if (P_node_list->nodes == NULL) {
      PERROR_MALLOC();
      ERROR("malloc()");
      exit(EXIT_FAILURE);
    }

    P_node_list->alloc_space += node_bucket;

// allocate bucket of space
#if defined(__QC_SGI)
    P_node_list->nodes[i_node] =
        memalign(page_size, node_bucket * sizeof(struct node_t));

    change_memory_placement(P_node_list->nodes[i_node],
                            node_bucket * sizeof(struct node_t), 0);
#else
    P_node_list->nodes[i_node] =
        (struct node_t *)calloc(node_bucket, sizeof(struct node_t));
#endif /* sgi */

    if (P_node_list->nodes[i_node] == NULL) {
      PERROR_MALLOC();
      exit(EXIT_FAILURE);
    }

#if defined(__QC_SGI)
    memset(P_node_list->nodes[i_node], 0, node_bucket * sizeof(struct node_t));
#endif /* sgi */

    pthread_mutex_init(&P_node_list->nodes[i_node]->lock, NULL);
    pthread_mutex_init(&P_node_list->nodes[i_node]->site_cluster.lock, NULL);

    // index into the bucket
    for (j_node = 1; j_node < node_bucket; j_node++) {
      P_node_list->nodes[i_node + j_node] =
          P_node_list->nodes[i_node + j_node - 1] + 1;
      pthread_mutex_init(&P_node_list->nodes[i_node + j_node]->lock, NULL);
      pthread_mutex_init(
          &P_node_list->nodes[i_node + j_node]->site_cluster.lock, NULL);
    }
  }

  return;
} // end of makeActiveNode()

//
//	makeInactiveNode()
//
void Node::makeInactiveNode(struct node_list_t *P_node_list, const int i_node) {
  // check space on the inactive queue
  if (P_node_list->alloc_inactive_space < i_node + 1) {
    if (P_node_list->alloc_inactive_space == 0)
      P_node_list->inactive_nodes =
          (struct node_t **)malloc(node_bucket * sizeof(struct node_t *));
    else
      P_node_list->inactive_nodes = (struct node_t **)realloc(
          P_node_list->inactive_nodes,
          (P_node_list->alloc_inactive_space + node_bucket) *
              sizeof(struct node_t *));

    if (P_node_list->inactive_nodes == NULL) {
      PERROR_MALLOC();
      ERROR("malloc()");
      exit(EXIT_FAILURE);
    }

    P_node_list->alloc_inactive_space += node_bucket;
  }

  P_node_list->inactive_nodes[P_node_list->number_inactive_nodes] =
      P_node_list->nodes[i_node];
  P_node_list->number_inactive_nodes++;

  P_node_list->nodes[i_node] = NULL;

  return;
} // end of makeInactiveNode()

// ***************************************//
// I = node_clusters.c

//
//	isNodeAtSurface()
//
int Node::isNodeAtSurface(const struct node_t *P_node) {
  if (P_node->fix_mask & SURFACE_MASK)
    return (1);

  return (0);
} // end of isNodeAtSurface()

//
//	createNodalClusters()
//
void Node::createNodalClusters(int iQuasi, struct node_list_t *P_node_list,
                               struct lattice_t *P_lattice, int *cluster_size,
                               enum mt_version_t mt_version) {
  // create clusters
  BuildNodalClusters(iQuasi, P_node_list, P_lattice, cluster_size, mt_version);

  // resolve overlaps
  ResolveOverlaps(iQuasi, P_node_list, P_lattice, cluster_size, mt_version);

  // check clusters for overlaps
  CheckClusters(P_node_list, P_lattice);

  return;
} // end of createNodalClusters()

// ***************************************//
// I = node_clusters_weight.c

//
//	nodalClusterWeightTest()
//
int Node::nodalClusterWeightTest(const int iQuasi,
                                 const struct node_list_t *P_node_list,
                                 struct lattice_t *P_lattice) {
  double rel_error;
  double exact_sum = 0.0;
  double approx_sum = 0.0;

  // compute exact sum
  exact_sum = ComputeLatticeExactSum(iQuasi, P_lattice);

  /**
   * compute approx sum
   */

  approx_sum = ComputeLatticeApproxSum(iQuasi, P_lattice, P_node_list);

  /**
   *
   */

  rel_error = fabs((approx_sum - exact_sum) / exact_sum);

  d_print("***** Weight test ******\n");
  d_print("exact sum = % e approx sum = % e rel. error = % e\n", exact_sum,
          approx_sum, rel_error);

  if (rel_error >= 1.0e-6)
    abort();

  /**
  for(i_node=0; i_node < P_node_list->number_nodes; i_node++)
    printf("[%d] % e\n", i_node, P_node_list->nodes[i_node]->weight);
  */

  return (0);
} // end of nodalClusterWeightTest()

//
//	computeNodalClusterSitesAndWeights()
//
void Node::computeNodalClusterSitesAndWeights(int iQuasi,
                                              struct lattice_t *P_lattice,
                                              struct node_list_t *P_node_list,
                                              double r_cluster,
                                              enum mt_version_t mt_version) {
  int *cluster_size;
  int *cluster_incr;

  int i_incr;
  int i_node;

  // create cluster_sizes
  if ((cluster_size = (int *)malloc(P_node_list->number_nodes * sizeof(int))) ==
      NULL) {
    ERROR("malloc()");
    exit(EXIT_FAILURE);
  }

  if ((cluster_incr = (int *)malloc(P_node_list->number_nodes * sizeof(int))) ==
      NULL) {
    ERROR("malloc()");
    exit(EXIT_FAILURE);
  }

#if defined(__QC_SGI)
  change_memory_placement(cluster_size, P_node_list->number_nodes * sizeof(int),
                          0);
  change_memory_placement(cluster_incr, P_node_list->number_nodes * sizeof(int),
                          0);
#endif /* sgi */

  for (i_node = 0; i_node < P_node_list->number_nodes; i_node++) {
    if (isNodeAtSurface(P_node_list->nodes[i_node]))
      cluster_size[i_node] = CLUSTER_SURFACE_SIZE;
    else
      cluster_size[i_node] = CLUSTER_SIZE;
  }

  // check for negative weights; if any are present increase
  // cluster size
  for (i_incr = 0; i_incr < CLUSTER_INCR_MAX; i_incr++) {
    // create initial clusters
    setbuf(stdout, NULL);
    d_print("Creating clusters...");
    // call Node::createNodalClusters()
    createNodalClusters(iQuasi, P_node_list, P_lattice, cluster_size,
                        mt_version);
    d_print("Completed\n");

    // compute weights; call Node::ComputeClusterNodalWeightsLamped()
    ComputeClusterNodalWeightsLamped(iQuasi, P_node_list, P_lattice);
    break;

    // check for negative weights
    if (checkWeightsLocal(P_node_list) == 0)
      break;

    d_print("Some weights are negative-increasing cluster radius.\n");

    // call negativeWeightInfoLocal(P_node_list, P_lattice, iQuasi);
    // to know more about nodes with negative weights.
    // exit(EXIT_SUCCESS);

    // clean clusters
    zeroClustersLocal(P_node_list);

    // increase cluster radius
    memset(cluster_incr, 0, sizeof(int) * P_node_list->number_nodes);

    for (i_node = 0; i_node < P_node_list->number_nodes; i_node++)
      if (P_node_list->nodes[i_node]->weight < 0.0) {
        struct node_t *P_node = P_node_list->nodes[i_node];
        int j_node;

        if (cluster_incr[i_node] == 0) {
          cluster_size[i_node] += 1;
          cluster_incr[i_node] = 1;
        }

        for (j_node = 0; j_node < P_node->node_list.number_nodes; j_node++) {
          int node_j_number = P_node->node_list.nodes[j_node]->number;

          if (cluster_incr[node_j_number] == 0) {
            cluster_size[node_j_number] += 1;
            cluster_incr[node_j_number] = 1;
          }
        }

        d_print("size of the cluster at node %d increased to %d\n", i_node,
                cluster_size[i_node]);
        d_print("number of sites in a cluster %d\n",
                P_node_list->nodes[i_node]->site_cluster.number_sites);
        d_print("weight % 8.6e\n", P_node_list->nodes[i_node]->weight);
      }
  }

  if (i_incr == CLUSTER_INCR_MAX) {
    d_print("Clusters has been incresed by %d increments\n", i_incr);
    d_print("but some weights are still negative\n");
    D_ERROR("computeNodalClusterSitesAndWeights()");
    exit(EXIT_FAILURE);
  }

  // cleanup
  free(cluster_size);
  free(cluster_incr);

  return;
} // end of computeNodalClusterSitesAndWeights()

// ***************************************//
// I = function constructed from triangulation_client.c

//
//  checkForDuplicateNodes()
//
void Node::checkForDuplicateNodes(struct node_list_t *P_node_list) {
  enum mt_version_t version = MULTI_THREADED;

  /*determine whether multithreaded*/
  if (get_max_number_threads() == 1) {
    version = SINGLE_THREADED;
  }

  switch (version) {
  case MULTI_THREADED: {
    struct checkForDuplicateNodesData_t data;

    data.P_node_list = P_node_list;

    /*pass data to multithreading routine*/
    thread_monitor(checkForDuplicateNodesWorker, (void *)&data,
                   get_max_number_threads());
  } break;

  case SINGLE_THREADED:
    /*runs over whole list*/
    checkForDuplicateNodesProcessingLocal(P_node_list, 0,
                                          P_node_list->number_nodes);
    break;
  }

  return;
} // end of checkForDuplicateNodes()

// ***************************************//
// I = function of remesh_elements.c

//
//  addNewNode()
//
//	add new node to the list of new nodes;
//  this is run from multiple threads so locking is required
//
void Node::addNewNode(const int l[3], struct all_node_list_t *P_node_list,
                      struct lattice_t *P_lattice, const int iQuasi) {
  double S[4]; // initial state: S[3] is initial frequency
  double s[4]; // current state: s[3] is current frequency

  double X[3];
  double x[3];

  double W;
  double T;
  double Tau;
  double w;
  double t;
  double tau;

  double atomicMass =
      Quasicontinua::getInstance()->getQuasi(iQuasi).getAtomicMass();
  double boltzmanConstant =
      PairPotentials::getInstance()->getBoltzmanConstant();

  int i_node;
  int number_new_nodes;

  // compute initial state

  Lattice *latticeC = Lattice::getInstance();

  latticeC->getSiteInitialState(S, l, P_lattice, iQuasi);

  // compute current state
  latticeC->getSiteCurrentState(s, l, P_lattice, iQuasi);

  // copy it to relevant data

  for (int i = 0; i < 3; i++) {
    X[i] = S[i];
    x[i] = s[i];
  }

  W = S[3];
  T = Quasicontinua::getInstance()->getTemperature();
  Tau = sqrt(atomicMass * boltzmanConstant * T) / W;

  w = s[3];
  t = T;
  tau = sqrt(atomicMass * boltzmanConstant * T) / w;

  // lock mutex to write on node data
  pthread_mutex_lock(&P_node_list->new_node_list.lock);

  // check if site has not already been added; scan the list of
  // new nodes if compare_sites() returns RETURN_SUCCESS bail out
  // (releasing a lock)
  for (i_node = 0; i_node < P_node_list->new_node_list.number_nodes; i_node++)
    if (latticeC->compareSites(
            l, P_node_list->new_node_list.nodes[i_node]->l) == RETURN_SUCCESS) {
      // site has already been added
      pthread_mutex_unlock(&P_node_list->new_node_list.lock);
      return;
    }

  number_new_nodes = P_node_list->new_node_list.number_nodes;
  P_node_list->new_node_list.number_nodes++;

  makeActiveNode(&P_node_list->new_node_list, number_new_nodes);

  // fill in the data; this CANNOT be moved outside lock
  // since realloc() may move contents of memory to some
  // other location, so we then would fill some garbage
  memcpy(&(P_node_list->new_node_list.nodes[number_new_nodes]
               ->initial_position[0]),
         &(X[0]), sizeof(X));

  memcpy(&(P_node_list->new_node_list.nodes[number_new_nodes]->position[0]),
         &(x[0]), sizeof(x));

  memcpy(
      &(P_node_list->new_node_list.nodes[number_new_nodes]->initial_frequency),
      &(W), sizeof(W));

  memcpy(&(P_node_list->new_node_list.nodes[number_new_nodes]->frequency), &(w),
         sizeof(w));

  memcpy(&(P_node_list->new_node_list.nodes[number_new_nodes]
               ->initial_temperature),
         &(T), sizeof(T));

  memcpy(&(P_node_list->new_node_list.nodes[number_new_nodes]->temperature),
         &(t), sizeof(t));

  memcpy(&(P_node_list->new_node_list.nodes[number_new_nodes]->initial_tau),
         &(Tau), sizeof(Tau));

  memcpy(&(P_node_list->new_node_list.nodes[number_new_nodes]->tau), &(tau),
         sizeof(tau));

  P_node_list->new_node_list.nodes[number_new_nodes]->number =
      P_node_list->node_list.number_nodes + number_new_nodes;

  // copy get the mask for the node
  P_node_list->new_node_list.nodes[number_new_nodes]->fix_mask = FREE_MASK;

  // fix_w_mask
  P_node_list->new_node_list.nodes[number_new_nodes]->fix_w_mask = FREE_MASK;

  // save site
  P_node_list->new_node_list.nodes[number_new_nodes]->l[0] = l[0];
  P_node_list->new_node_list.nodes[number_new_nodes]->l[1] = l[1];
  P_node_list->new_node_list.nodes[number_new_nodes]->l[2] = l[2];

  // done writing to node, so unlock the mutex
  pthread_mutex_unlock(&P_node_list->new_node_list.lock);

  return;
} // end of addNewNode()

// *********** private functions  ************ //

//
//	BuildNodalClusters()
//
void Node::BuildNodalClusters(int iQuasi, struct node_list_t *P_node_list,
                              struct lattice_t *P_lattice, int *cluster_size,
                              enum mt_version_t mt_version) {
  switch (mt_version) {
  case MULTI_THREADED: {
    struct BuildNodalClustersData_t data;

    data.iQuasi = iQuasi;
    data.P_node_list = P_node_list;
    data.P_lattice = P_lattice;
    data.cluster_size = cluster_size;

    thread_monitor(BuildNodalClustersWorker, (void *)&data,
                   get_max_number_threads());
  } break;

  case SINGLE_THREADED:

    BuildNodalClustersProcessingLocal(iQuasi, P_node_list, P_lattice,
                                      cluster_size, 0,
                                      P_node_list->number_nodes - 1);
    break;
  }

  return;
} // end of BuildNodalClusters()

//
//	ResolveOverlaps()
//
void Node::ResolveOverlaps(int iQuasi, struct node_list_t *P_node_list,
                           struct lattice_t *P_lattice, int *cluster_size,
                           enum mt_version_t mt_version) {
  struct ResolveOverlapsData_t ResolveOverlapsData;

  switch (mt_version) {
  case SINGLE_THREADED:
    ResolveOverlapsProcessingLocal(iQuasi, P_node_list, P_lattice, cluster_size,
                                   0, P_node_list->number_nodes - 1);
    break;

  case MULTI_THREADED:
    current_cluster = 0;

    ResolveOverlapsData.iQuasi = iQuasi;
    ResolveOverlapsData.P_node_list = P_node_list;
    ResolveOverlapsData.P_lattice = P_lattice;
    ResolveOverlapsData.cluster_size = cluster_size;

    thread_monitor(ResolveOverlapsWorker, (void *)&ResolveOverlapsData,
                   get_max_number_threads());

    break;
  }

  return;
} // end of ResolveOverlaps()

//
//	CheckClusters()
//
void Node::CheckClusters(const struct node_list_t *P_node_list,
                         const struct lattice_t *P_lattice) {
  // at this point I have not implemented this function.

  // Later if it's required, see node_clusters.c and implement the
  // check_clusters() function.

  return;
} // end of CheckClusters()

//
//	ComputeLatticeExactSum()
//
double Node::ComputeLatticeExactSum(const int iQuasi,
                                    struct lattice_t *P_lattice) {
  double sum = 0.0;

  int l[3];

  // generate sites
  for (l[0] = P_lattice->l_start[0]; l[0] <= P_lattice->l_end[0]; l[0]++)
    for (l[1] = P_lattice->l_start[1]; l[1] <= P_lattice->l_end[1]; l[1]++)
      for (l[2] = P_lattice->l_start[2]; l[2] <= P_lattice->l_end[2]; l[2]++) {
        double X[3];

        if (Lattice::getInstance()->isLatticeSite(l, P_lattice, iQuasi) ==
            RETURN_FAILURE)
          continue;

        // find material coordinates
        Lattice::getInstance()->getSiteInitialPosition(X, l, P_lattice, iQuasi);

        // compute value of the test function at the current site

        // sum += WeightTestFunctionLocal(X);
        sum += 1.0; // WeightTestFunctionLocal() returns 1.0 always.
      }

  return (sum);
} // end of ComputeLatticeExactSum()

//
//	ComputeLatticeApproxSum()
//
double Node::ComputeLatticeApproxSum(const int iQuasi,
                                     struct lattice_t *P_lattice,
                                     const struct node_list_t *P_node_list) {
  double sum = 0.0;

  int i_node;
  int i_site;

  // loop over all nodes
  for (i_node = 0; i_node < P_node_list->number_nodes; i_node++) {
    struct node_t *P_node = P_node_list->nodes[i_node];

    for (i_site = 0; i_site < P_node->site_cluster.number_sites; i_site++) {
      int *l = P_node->site_cluster.sites[i_site];
      double X[3];

      // find material coordinates
      Lattice::getInstance()->getSiteInitialPosition(X, l, P_lattice, iQuasi);

      // compute value of the test function
      // sum += P_node->weight*WeightTestFunctionLocal(X);
      sum += P_node->weight * 1.0; // WeightTestFunctionLocal() returns 1.0;
    }
  }

  return (sum);
} // end of ComputeLatticeApproxSum()

//
//	ComputeClusterNodalWeightsLamped()
//
void Node::ComputeClusterNodalWeightsLamped(int iQuasi,
                                            struct node_list_t *P_node_list,
                                            struct lattice_t *P_lattice) {
  struct ComputeClusterNodalWeightsLampedData_t data;

  // initialize current node
  lamped_current_node = 0;

  // fill in data
  data.iQuasi = iQuasi;
  data.P_node_list = P_node_list;
  data.P_lattice = P_lattice;

  // call thread function
  thread_monitor(ComputeClusterNodalWeightsLampedWorker, (void *)&data,
                 get_max_number_threads());

  return;
} // end of ComputeClusterNodalWeightsLamped()

//
//  isNodeAtomistic()
//
bool Node::isNodeAtomistic(const struct node_t *P_node,
                           const double atomistic_element_size) {
  const double tolerance = 0.0000001;

  // check weight
  if (P_node->weight > 1.0 + tolerance || P_node->weight < 1.0 - tolerance)
    return false;

  // check element list of node
  int iElement;
  struct element_list_t element_list = P_node->element_list;

  for (iElement = 0; iElement < element_list.number_elements; iElement++) {
    const struct element_t *P_element = element_list.elements[iElement];

    if (MiscFunctions::getInstance()->tetraVolume(
            P_element->node[0]->position, P_element->node[1]->position,
            P_element->node[2]->position,
            P_element->node[3]->position) > atomistic_element_size)

      return false;
  }

  // if reached here, node is atomistic
  return true;
} // end of isNodeAtomistic()
}
