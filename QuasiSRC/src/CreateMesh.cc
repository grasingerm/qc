#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <pthread.h>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>

#include <algorithm> // for std::sort()
#include <assert.h>

// C++ source
#include "CGNonLinearSolver.h"
#include "C_Interface.h"
#include "CreateMesh.h"
#include "DataTypes.h"
#include "Electrostatics.h"
#include "Element.h"
#include "ForceEnergyCalculation.h"
#include "Indent.h"
#include "Lattice.h"
#include "MiscFunctions.h"
#include "Node.h"
#include "Output.h"
#include "PairPotentials.h"
#include "Quasicontinua.h"
#include "QuasicontinuaForceFunction.h"
#include "Void.h"

// C source
#include "monitor.h"
#include "numa.h"
#include "threads.h"
#include "triangulation_client.h"

// global variables
#define REL_ENERGY_ERR_MAX 0.01
#define CYLINDER_RADIUS 2.0
#define CYLINDER_HEIGHT 10.0
#define VICINITY_RADIUS 1.5

namespace quasicontinuum {

CreateMesh *CreateMesh::_instance = NULL;

namespace {

bool CompareEdgesLocal(const CreateMesh::edge_t &lhs,
                       const CreateMesh::edge_t &rhs) {
  return lhs.length > rhs.length;
}

bool CompareScoreListLocal(const CreateMesh::score_list_t &lhs,
                           const CreateMesh::score_list_t &rhs) {
  return lhs.energy_drop < rhs.energy_drop;
}

double CheckEnergyStrainLocal(struct element_t *P_element, const int i_load) {
  double center[3];
  double err_norm;
  double radius;

  err_norm = Element::getInstance()->computeElementStrain(P_element, i_load);

  radius = MiscFunctions::getInstance()->findTetraCenter(
      &(P_element->node[0]->initial_position[0]),
      &(P_element->node[1]->initial_position[0]),
      &(P_element->node[2]->initial_position[0]),
      &(P_element->node[3]->initial_position[0]), center);

  if (radius < 0.0)
    abort();

  err_norm = sqrt(fabs(err_norm));

  return (err_norm);
} // end of CheckEnergyStrainLocal()

/**
  * check_energy worker thread
  */
// mutex lock, to write to P_score_list
pthread_mutex_t score_list_lock;

struct CheckEnergyData_t {
  std::vector<CreateMesh::score_list_t> *P_score_list;
  struct element_list_t *P_element_list;
  struct all_node_list_t *P_node_list;
  struct indentor_t *P_indentor;
  struct boundingbox_t bbox;
  double r_contact_zone;
  int i_load;
};

void *CheckEnergyWorker(void *arg) {

  std::vector<CreateMesh::score_list_t> *P_score_list =
      ((struct CheckEnergyData_t *)arg)->P_score_list;
  struct element_list_t *P_element_list =
      ((struct CheckEnergyData_t *)arg)->P_element_list;
  const struct indentor_t *P_indentor =
      ((struct CheckEnergyData_t *)arg)->P_indentor;

  double r_contact_zone = ((struct CheckEnergyData_t *)arg)->r_contact_zone *
                          VICINITY_RADIUS * VICINITY_RADIUS;

  const int i_load = ((struct CheckEnergyData_t *)arg)->i_load;

  int number_threads = get_number_threads();
  int my_id = get_my_tid();

  int number_elem_thread;
  int number_elem_start;
  int number_elem_end;

  int i_elem;

  get_share(my_id, number_threads, P_element_list->number_elements,
            &number_elem_thread, &number_elem_start, &number_elem_end);

  /**
    * loop over all elements
    */

  for (i_elem = number_elem_start; i_elem <= number_elem_end; i_elem++) {
    struct element_t *P_element = P_element_list->elements[i_elem];

    CreateMesh::score_list_t score_list;

    score_list.energy_drop = CheckEnergyStrainLocal(P_element, i_load);

    score_list.P_element = P_element;

    pthread_mutex_lock(&score_list_lock);
    P_score_list->push_back(score_list);
    pthread_mutex_unlock(&score_list_lock);
  }

  return (NULL);
} // end of CheckEnergyWorker()

} // end of namespace for CheckEnergyWorker()

CreateMesh::CreateMesh() {
  return;
}

CreateMesh::CreateMesh(const double remesh_tolerance) {
  d_remeshTolerance = remesh_tolerance;

  return;
}

CreateMesh::~CreateMesh() {
  return;
}

CreateMesh *CreateMesh::getInstance() {
  if (_instance == NULL) {
    _instance = new CreateMesh(0.0);
  }

  return _instance;
}

void CreateMesh::destroyInstance() {
  delete _instance;

  return;
}

//
//  createElementsCentral()
//
void CreateMesh::createElementsCentral(int computeForce) {
  // get quasicontinua class instance
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  int numQuasi = quasicontinua->size();

  // loop over quasi and create mesh
  for (int iQuasi = 0; iQuasi < numQuasi; iQuasi++) {
    // data from iQuasi
    struct element_list_t &element_list =
        quasicontinua->getQuasi(iQuasi).getElementList();

    struct all_node_list_t &all_node_list =
        quasicontinua->getQuasi(iQuasi).getNodeList();

    struct lattice_t lattice = quasicontinua->getQuasi(iQuasi).getLattice();

    int restartFlag;
    // if 1 - restart is on
    restartFlag = quasicontinua->getQuasi(iQuasi).isRestartOn();

    // get cutoff radius for cluster site calculation
    const int potentialNumber =
        PairPotentials::getInstance()->doQuasicontinuumInteract(iQuasi, iQuasi);

    double r_cluster_cutoff;

    if (potentialNumber != -1) {
      r_cluster_cutoff =
          PairPotentials::getInstance()->getCutoffRadiusClusterList(
              potentialNumber);
    } else {
      r_cluster_cutoff = 0.00001;
    }

    // create elements for iQuasi
    d_print("CreateMesh : Creating elements for quasi = %i\n", iQuasi);
    createElements(iQuasi, &element_list, &all_node_list.node_list, &lattice,
                   r_cluster_cutoff, restartFlag);

    // initialize hash table in site-element map cache
    Element::getInstance()->checkHashTableAllocation(iQuasi);
    // print message that meshing of quasicontinuum is complete

    d_print("CreateMesh : Done\n");
  }

  //  initialize CrossNeighborList and fill neighbor list for force
  //  calculation
  int rebuild_neighbor_flag = 1;
  int rebuild_cluster_flag = 1;
  // d_print("CrossNeighborList...");
  CrossNeighborList::getInstance()->computeCrossNeighborList(
      rebuild_neighbor_flag, rebuild_cluster_flag);
  d_print("complete\n");

  // compute electrostatics geometry for force calculation
  Electrostatics *electroC = Electrostatics::getInstance();

  if (electroC->isElectrostaticEnable() == 1) {
    electroC->computeGeometries();
    d_print("complete\n");
  }

  if (computeForce == 0) {
    // compute forces
    d_print("Residual Forces...");
    ForceEnergyCalculation::getInstance()->removeResidualForces();
    d_print("complete\n");

    bool compute_in_reference = false;
    rebuild_neighbor_flag = 1;
    d_print("Forces...");
    ForceEnergyCalculation::getInstance()->computeForceEnergy(
        rebuild_neighbor_flag, compute_in_reference);
    d_print("complete\n");

    //
    // in very first call to createElementsCentral(), copy the energy to
    // initial_energy of
    //  all node list data.
    //
    ForceEnergyCalculation::getInstance()->copyEnergyToInitialEnergy();
  }

  return;
}

//
//  createElements()
//
void CreateMesh::createElements(int iQuasi,
                                struct element_list_t *P_element_list,
                                struct node_list_t *P_node_list,
                                struct lattice_t *P_lattice,
                                double r_cluster_cutoff, int restartFlag) {
  
  Node *nodeC = Node::getInstance();
  Element *elementC = Element::getInstance();

  //  if restartFlag is 1 - elements are already created. We only need
  //  to compute nodal weights, site-element-map etc
  //
  if (restartFlag == 0) {
    // mesh data for meshing
    struct mesh_data_t tetra_data_mesh = MESHDATA_INITIALIZER;

    // check for duplicate node in node list before passing it to
    // triangulation client.
    nodeC->checkForDuplicateNodes(P_node_list);

    // call create_elements_external() of triangulation_client.c
    create_elements_external(P_node_list, &tetra_data_mesh);

    // from new mesh_data create the elements. This also calls the
    // function computeElementNumberSites() and fills the
    // P_element->number_sites
    elementC->createElementsFromMeshData(P_element_list, P_node_list, P_lattice,
                                         &tetra_data_mesh);
  }

  /**
    * create void
    */
  Void *voidClass = Void::getInstance();

  if (voidClass->isVoidEnable() != 0) {
    voidClass->createVoid(P_node_list, P_lattice, iQuasi);
  }

  // call fixBoundary() member of createMesh class
  fixBoundary(P_node_list, P_element_list, P_lattice, restartFlag);

  /**
    * build BSP tree for site -> element mapping
    */
  elementC->initializeSiteElementMapData(iQuasi, P_element_list);

  // call binElements()
  elementC->binElements(P_element_list, P_node_list, MULTI_THREADED);

  // NOTE: MULTI_THREADED changes problem (1 thread vs. 2 threads
  // will not match)
  nodeC->computeNodalClusterSitesAndWeights(iQuasi, P_lattice, P_node_list,
                                            r_cluster_cutoff,
                                            SINGLE_THREADED); // See note above

  /**
    * correct boundary conditions (if needed)
    */
  fixBoundary(P_node_list, P_element_list, P_lattice, restartFlag);

  return;
} // end of createElements()

// fix_boundary.c related functions

//
//  fixBoundary()
//
void CreateMesh::fixBoundary(struct node_list_t *P_node_list,
                             struct element_list_t *P_element_list,
                             struct lattice_t *P_lattice, int restartFlag) {
  if (restartFlag == 0) {
    // fix nodes for the indentor problem
    if (Indent::getInstance()->isIndentEnable() != 0) {
      Indent::getInstance()->fixIndentBoundary(P_node_list, P_lattice);
    }

    if (Void::getInstance()->isVoidEnable() != 0) {
      Void::getInstance()->fixVoidBoundary(P_node_list, P_element_list);
    }
  } else {
    Void::getInstance()->removeVoidNodes(P_node_list);
  }

  return;
}

void CreateMesh::setRemeshTolerance(double tolerance) {
  d_remeshTolerance = tolerance;

  return;
}

int CreateMesh::remeshCentral(int i_load) {
  // get quasicontinua class instance
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  int numQuasi = quasicontinua->size();

  int total_new_nodes = 0;

  for (int iQuasi = 0; iQuasi < numQuasi; iQuasi++) {
    // data from iQuasi
    struct element_list_t &element_list =
        quasicontinua->getQuasi(iQuasi).getElementList();

    struct all_node_list_t &all_node_list =
        quasicontinua->getQuasi(iQuasi).getNodeList();

    struct lattice_t lattice = quasicontinua->getQuasi(iQuasi).getLattice();

    struct indentor_t &indentor = quasicontinua->getQuasi(iQuasi).getIndentor();

    // get cutoff radius for cluster site calculation
    const int potentialNumber =
        PairPotentials::getInstance()->doQuasicontinuumInteract(iQuasi, iQuasi);

    double r_cluster_cutoff;

    if (potentialNumber != -1) {
      r_cluster_cutoff =
          PairPotentials::getInstance()->getCutoffRadiusClusterList(
              potentialNumber);
    } else {
      r_cluster_cutoff = 0.00001;
    }

    // create elements for iQuasi
    d_print("CreateMesh : Remeshing quasi = %i\n", iQuasi);
    total_new_nodes += remesh(iQuasi, &element_list, &all_node_list, &lattice,
                              &indentor, r_cluster_cutoff, i_load);
    d_print("CreateMesh : Done\n");

    //
    //  initialize trace of K data in quasicontinuum
    //
    quasicontinua->getQuasi(iQuasi).initializeTraceOfKData();
  }

  //  initialize CrossNeighborList and fill neighbor list for force
  //  calculation
  int rebuild_neighbor_flag = 1;
  int rebuild_cluster_flag = 1;
  // d_print("CreateMesh : CrossNeighborList...");
  CrossNeighborList::getInstance()->computeCrossNeighborList(
      rebuild_neighbor_flag, rebuild_cluster_flag);
  d_print("complete\n");

  // compute electrostatics geometry for force calculation
  Electrostatics *electroC = Electrostatics::getInstance();

  if (electroC->isElectrostaticEnable() == 1) {
    electroC->clearGeometries();
    electroC->computeGeometries();
    d_print("complete\n");
  }

  d_print("Residual Forces...");
  ForceEnergyCalculation::getInstance()->removeResidualForces();
  d_print("complete\n");

  bool compute_in_reference = false;
  rebuild_neighbor_flag = 1;
  d_print("Force...");
  ForceEnergyCalculation::getInstance()->computeForceEnergy(
      rebuild_neighbor_flag, compute_in_reference);
  d_print("complete\n");

  return (total_new_nodes);
}

int CreateMesh::remesh(int iQuasi, struct element_list_t *P_element_list,
                       struct all_node_list_t *P_node_list,
                       struct lattice_t *P_lattice,
                       struct indentor_t *P_indentor, double r_cluster_cutoff,
                       int i_load) {
  std::vector<score_list_t> score_lists;

  int i_elem;
  int n_new_nodes;

  // get pointer to Element and Node class
  Element *elementC = Element::getInstance();

  Node *nodeC = Node::getInstance();

  CheckEnergyRemesh(&score_lists, P_element_list, P_node_list, P_indentor,
                    i_load, iQuasi);

  /**
    * sort score list in descending order
    */
  assert(score_lists.size() == P_element_list->number_elements);

  std::sort(score_lists.begin(), score_lists.end(), CompareScoreListLocal);

#if defined(_REMESH_FORCE_QUICK_DEBUG_)
  DumpScoreList(score_lists);
  exit(0);
#endif /* _REMESH_FORCE_QUICK_DEBUG_ */

  /**
    * keep poping elements from sorted score list
    */

  i_elem = 0;

  while (score_lists[i_elem].energy_drop > d_remeshTolerance &&
         i_elem < P_element_list->number_elements) {
    RemeshElementBisection(score_lists[i_elem].P_element, P_node_list,
                           P_lattice, iQuasi);

    i_elem++;
  }

  /**
    * prepare for triangulation
    */

  n_new_nodes = P_node_list->new_node_list.number_nodes;

  d_print("CreateMesh : max error=%e\n", score_lists[0].energy_drop);
  d_print("CreateMesh : added %d nodes\n", n_new_nodes);

  // we are using CrossNeighborList for neighbor list calculation.

  elementC->cleanElements(P_element_list, MULTI_THREADED);

  Void::getInstance()->insertVoidNodes(&P_node_list->node_list);

  nodeC->cleanNodes(P_node_list);

  elementC->hashTableClean(iQuasi);
  elementC->cleanSiteElementMapData(iQuasi);

  // triangulate and create new elements.
  int restartFlag = 0; // 0 - OFF , 1 - ON
  createElements(iQuasi, P_element_list, &(P_node_list->node_list), P_lattice,
                 r_cluster_cutoff, restartFlag);

  return (n_new_nodes);
}

void CreateMesh::RemeshElementBisection(struct element_t *P_element,
                                        struct all_node_list_t *P_node_list,
                                        struct lattice_t *P_lattice,
                                        const int iQuasi) {
  int l[3];

  if ((FindBisectionSite(l, P_lattice, P_element, iQuasi)) == RETURN_SUCCESS) {
    Node::getInstance()->addNewNode(l, P_node_list, P_lattice, iQuasi);
  }

  return;
}

int CreateMesh::FindBisectionSite(int l[3], struct lattice_t *P_lattice,
                                  const struct element_t *P_element,
                                  const int iQuasi) {
  std::vector<struct element_t *> neigh_elements;

  struct node_t *P_node_0 = P_element->node[0];
  struct node_t *P_node_1 = P_element->node[1];
  struct node_t *P_node_2 = P_element->node[2];
  struct node_t *P_node_3 = P_element->node[3];

  std::vector<edge_t> edges;
  edge_t temp_edge;

  double center[3];
  double X[3];
  double r;
  double r_min = 0.0;

  int l_center[3];
  int l_i[3];
  int n_neigh_elements;
  int i_node;
  int i_elem;
  int i_site;
  int stop_flag = 0;
  int shell_number = 1;
  int site_found = 0;

  /**
    * fill in edges
    */
  MiscFunctions *miscC = MiscFunctions::getInstance();
  Lattice *latticeC = Lattice::getInstance();

  temp_edge.length = miscC->getDistanceSqrBetweenPoints(
      P_node_0->initial_position, P_node_1->initial_position);
  temp_edge.node[0] = P_node_0;
  temp_edge.node[1] = P_node_1;

  // push_back this to edges
  edges.push_back(temp_edge);

  temp_edge.length = miscC->getDistanceSqrBetweenPoints(
      P_node_0->initial_position, P_node_2->initial_position);
  temp_edge.node[0] = P_node_0;
  temp_edge.node[1] = P_node_2;

  // push_back this to edges
  edges.push_back(temp_edge);

  temp_edge.length = miscC->getDistanceSqrBetweenPoints(
      P_node_0->initial_position, P_node_3->initial_position);
  temp_edge.node[0] = P_node_0;
  temp_edge.node[1] = P_node_3;

  // push_back this to edges
  edges.push_back(temp_edge);

  temp_edge.length = miscC->getDistanceSqrBetweenPoints(
      P_node_1->initial_position, P_node_2->initial_position);
  temp_edge.node[0] = P_node_1;
  temp_edge.node[1] = P_node_2;

  // push_back this to edges
  edges.push_back(temp_edge);

  temp_edge.length = miscC->getDistanceSqrBetweenPoints(
      P_node_1->initial_position, P_node_3->initial_position);
  temp_edge.node[0] = P_node_1;
  temp_edge.node[1] = P_node_3;

  // push_back this to edges
  edges.push_back(temp_edge);

  temp_edge.length = miscC->getDistanceSqrBetweenPoints(
      P_node_2->initial_position, P_node_3->initial_position);
  temp_edge.node[0] = P_node_2;
  temp_edge.node[1] = P_node_3;

  // push_back this to edges
  edges.push_back(temp_edge);

  /**
    * sort edges
    */

  // qsort(&(edges[0]), 6, sizeof(struct edge_t), CompareEdges);

  std::sort(edges.begin(), edges.end(), CompareEdgesLocal);

  /**
    * find center of the longest edge
    */

  center[0] = (edges[0].node[0]->initial_position[0] +
               edges[0].node[1]->initial_position[0]) /
              2.0;
  center[1] = (edges[0].node[0]->initial_position[1] +
               edges[0].node[1]->initial_position[1]) /
              2.0;
  center[2] = (edges[0].node[0]->initial_position[2] +
               edges[0].node[1]->initial_position[2]) /
              2.0;

  /**
    * locate lattice site closest to the center
    */

  latticeC->findClosestSite(l_center, center, P_lattice, iQuasi);

  /**
    * find all elements with common edge
    */

  FindCommonEdgeElements(&neigh_elements, edges[0].node[0], edges[0].node[1]);

  /**
    * if l_center is inside a cluster of elements with the common
    * edge and it is NOT a node we have a hit.
    */

  if (latticeC->isSiteInsideLattice(P_lattice, l_center, iQuasi) ==
      RETURN_SUCCESS)
    for (i_elem = 0; i_elem < neigh_elements.size(); i_elem++) {
      P_node_0 = neigh_elements[i_elem]->node[0];
      P_node_1 = neigh_elements[i_elem]->node[1];
      P_node_2 = neigh_elements[i_elem]->node[2];
      P_node_3 = neigh_elements[i_elem]->node[3];

      /**
        * check if center site is inside an element
        */

      if (miscC->checkSiteInTetra(
              P_node_0->initial_position, P_node_1->initial_position,
              P_node_2->initial_position, P_node_3->initial_position, l_center,
              P_lattice, iQuasi) == INSIDE) {
        if ((latticeC->compareSites(l_center, P_node_0->l) == RETURN_FAILURE) &&
            (latticeC->compareSites(l_center, P_node_1->l) == RETURN_FAILURE) &&
            (latticeC->compareSites(l_center, P_node_2->l) == RETURN_FAILURE) &&
            (latticeC->compareSites(l_center, P_node_3->l) == RETURN_FAILURE))
          stop_flag = 1;
      }

      /**
        * the l_center is the site; cleanup and return
        */
      if (stop_flag) {
        l[0] = l_center[0];
        l[1] = l_center[1];
        l[2] = l_center[2];

        return (RETURN_SUCCESS);
      }
    }

  /**
    * l_center is not the right site; scan all sites in the vicinity
    * of l_center
    */

  /**
    * initialize r_min
    */

  for (i_elem = 0, r_min = 0.0; i_elem < neigh_elements.size(); i_elem++)
    for (i_node = 0; i_node < 4; i_node++) {
      r = miscC->getDistanceSqrBetweenPoints(
          center, neigh_elements[i_elem]->node[i_node]->initial_position);
      r_min = MAX(r_min, r);
    }

  /**
    * scan all sites within common edge elements
    */

  stop_flag = 0;

  while (!stop_flag) {
    struct shell_t *i_shell = latticeC->getShell(shell_number, P_lattice);

    int i_hit = 0;

    /**
      * loop over all sites in the current shell
      */

    for (i_site = 0; i_site < i_shell->number_sites; i_site++) {

      l_i[0] = l_center[0] + i_shell->site[i_site][0];
      l_i[1] = l_center[1] + i_shell->site[i_site][1];
      l_i[2] = l_center[2] + i_shell->site[i_site][2];

      if (latticeC->isSiteInsideLattice(P_lattice, l_i, iQuasi) ==
          RETURN_FAILURE)
        continue;

      /**
        * check if the site is localted inside any of the elements
        */

      for (i_elem = 0; i_elem < neigh_elements.size(); i_elem++) {
        P_node_0 = neigh_elements[i_elem]->node[0];
        P_node_1 = neigh_elements[i_elem]->node[1];
        P_node_2 = neigh_elements[i_elem]->node[2];
        P_node_3 = neigh_elements[i_elem]->node[3];

        /**
          * check if site is inside an element
          */
        if (miscC->checkSiteInTetra(
                P_node_0->initial_position, P_node_1->initial_position,
                P_node_2->initial_position, P_node_3->initial_position, l_i,
                P_lattice, iQuasi) == INSIDE) {
          i_hit++;

          /**
            * check if site is not a node
            */
          if (latticeC->compareSites(l_i, P_node_0->l) == RETURN_SUCCESS)
            continue;
          if (latticeC->compareSites(l_i, P_node_1->l) == RETURN_SUCCESS)
            continue;
          if (latticeC->compareSites(l_i, P_node_2->l) == RETURN_SUCCESS)
            continue;
          if (latticeC->compareSites(l_i, P_node_3->l) == RETURN_SUCCESS)
            continue;

          /**
            * get material coordinates
            */
          latticeC->getSiteInitialPosition(X, l_i, P_lattice, iQuasi);

          /**
            * get distance between l and center
            */
          r = miscC->getDistanceSqrBetweenPoints(center, X);

          if (r < r_min) {
            r_min = r;

            l[0] = l_i[0];
            l[1] = l_i[1];
            l[2] = l_i[2];

            site_found = 1;
          }
        }
      }
    }

    /**
      * if i_hit is still 0 no sites where found in the current shell
      * that belong to any of the elements
      */

    if (i_hit == 0)
      stop_flag = 1;

    /**
      * go to the next shell
      */

    shell_number++;
  }

  if (site_found)
    return (RETURN_SUCCESS);

  return (RETURN_FAILURE);
} // end of FindBisectionSite()

//
//  FindCommonEdgeElements()
//  given two nodes finad all elements that have common edge
//
void CreateMesh::FindCommonEdgeElements(
    std::vector<struct element_t *> *P_neigh_elements,
    const struct node_t *P_node_1, const struct node_t *P_node_2) {
  P_neigh_elements->clear();

  int i_elem;

  /**
    * for P_node_1 scan the list of all elements
    * for each element check if it has P_node_2 on the list of nodes
    */

  for (i_elem = 0; i_elem < P_node_1->element_list.number_elements; i_elem++) {
    struct element_t *P_element = P_node_1->element_list.elements[i_elem];

    if (P_node_2 == P_element->node[0] || P_node_2 == P_element->node[1] ||
        P_node_2 == P_element->node[2] || P_node_2 == P_element->node[3]) {
      /**
        * common element found
        */

      P_neigh_elements->push_back(P_element);
    }
  }

  return;
}

//
//  CompareEdges()
//
bool CreateMesh::CompareEdges(const edge_t &lhs, const edge_t &rhs) {
  return lhs.length > rhs.length;
}

//
//  CompareScoreList()
//
bool CreateMesh::CompareScoreList(const score_list_t &lhs,
                                  const score_list_t &rhs) {
  return lhs.energy_drop < rhs.energy_drop;
}

//
//  CheckEnergyRemesh()
//
void CreateMesh::CheckEnergyRemesh(std::vector<score_list_t> *P_score_lists,
                                   struct element_list_t *P_element_list,
                                   struct all_node_list_t *P_node_list,
                                   struct indentor_t *P_indentor, int i_load,
                                   int iQuasi) {
  // clear the score_lit data before inserting score_list
  P_score_lists->clear();

  struct boundingbox_t bbox;

  struct CheckEnergyData_t checkEnergyData;

  Node::getInstance()->getBoundingBox(&bbox, &P_node_list->node_list,
                                      MULTI_THREADED);

  // get radius of the contact zone between indentor and sample
  Indent::getInstance()->findIndentorContactRadius(
      &checkEnergyData.r_contact_zone, NULL, &P_node_list->node_list,
      P_indentor, iQuasi, SINGLE_THREADED);

  // fill checkEnergyData
  checkEnergyData.P_score_list = P_score_lists;
  checkEnergyData.P_element_list = P_element_list;
  checkEnergyData.P_node_list = P_node_list;
  checkEnergyData.P_indentor = P_indentor;
  checkEnergyData.bbox = bbox;
  checkEnergyData.r_contact_zone =
      checkEnergyData.r_contact_zone * checkEnergyData.r_contact_zone;
  checkEnergyData.i_load = i_load;

  thread_monitor(CheckEnergyWorker, (void *)&checkEnergyData,
                 get_max_number_threads());

  return;
}

//
//  DumpScoreList()
//
void CreateMesh::DumpScoreList(std::vector<score_list_t> score_list) {
  for (int i_elem = 0; i_elem < score_list.size(); i_elem++) {
    d_print("[elam:%d] %e \n", score_list[i_elem].P_element->number,
            score_list[i_elem].energy_drop);
  }

  return;
}
}
