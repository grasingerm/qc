//
// ForceEnergyCalculation.cc
//
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#ifdef HAVE_LIMITS
#include <limits>
#else
#ifdef HAVE_LIMITS_H
#include <limits.h>
#else
#error No limits or limits.h available
#endif
#endif

#ifdef STDC_HEADERS
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#else
#error No standard C library headers found
#endif // STDC_HEADERS

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
// no  _REENTRANT
#else
#error No pthread.h available.
#endif // HAVE_PTHREAD_H

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

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#else
#error unistd.h not found.
#endif /* HAVE_UNISTD_H */

#include <iostream>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#include "C_Interface.h"
#include "CrossNeighborList.h"
#include "DataTypes.h"
#include "DataTypes.h"
#include "Electrostatics.h"
#include "Error.h"
#include "ForceEnergyCalculation.h"
#include "Indent.h"
#include "Input.h"
#include "Lattice.h"
#include "Node.h"
#include "Output.h"
#include "PairPotentials.h"
#include "QuadraturePoints.h"
#include "Quasicontinua.h"
#include "Quasicontinuum.h"
#include "RunData.h"
#include "Shape.h"

#include "monitor.h"
#include "threads.h"

// static double energyDebug=0.0;
static double forceDebug[3];
static pthread_mutex_t debugLock = PTHREAD_MUTEX_INITIALIZER;
static int residual_debug = 1;

static int dummy = 0;

// mean separation tolerance : if mean distance between two atoms
// is below this tolerance, do not compute force and energy
//
// this is to exclude atoms which are very close
static double meanSeparationTol = 0.3;

// kinetic energy and entropy energy needs to be in eV unit
static double energy_factor_kinetic = 1.036427059e-4;
static double energy_factor_entropy = 1.036427059e-4;

// static int l_debug[3];

//
//
//

namespace quasicontinuum {
//
// namespace for bucket lock init and cluster force addition
//
namespace {
static std::vector<std::vector<pthread_mutex_t>> bucketLock;

// not using newBucketLock currently
// static std::vector<std::vector<std::vector<pthread_mutex_t>>> newBucketLock;
static pthread_mutex_t lock_print = PTHREAD_MUTEX_INITIALIZER;

void bucketLockInitLocal(int numQuasi, int numBuckets) {
  // // std::cout<<"here bucket lock"<<numQuasi<<"    "<<numBuckets<<std::endl;
  // 	bucketLock.resize(numQuasi);

  // for(int i=0; i < bucketLock.size(); i++)
  // 	for(int j=0; j < numBuckets; j++)
  // 			bucketLock[i].push_back(PTHREAD_MUTEX_INITIALIZER);

  // 	// std::cout<<"here end of bucket lock"<<std::endl;

  // return;

  // std::cout<<"here bucket lock"<<numQuasi<<"    "<<numBuckets<<std::endl;
  if (bucketLock.size() == 0 || bucketLock.size() > numQuasi)
    bucketLock.resize(numQuasi);

  for (int i = 0; i < bucketLock.size(); i++) {
    if (bucketLock[i].size() == 0 || bucketLock[i].size() > numBuckets) {
      bucketLock[i].resize(numBuckets);

      for (int j = 0; j < numBuckets; j++)
        bucketLock[i][j] = PTHREAD_MUTEX_INITIALIZER;
    }
  }

  // std::cout<<"here end of bucket lock"<<std::endl;

  return;
}

//
//
//
struct AddClusterSiteForcesToNodesData_t {
  int iQuasi;
  std::vector<std::vector<cluster_site_data_t>> *P_clusterData;
  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      *P_clusterForceData;
};

//
//
//
struct ComputeTraceOfKAtNodesData_t {
  int iQuasi;
  std::vector<std::vector<cluster_site_data_t>> *P_clusterData;
  std::vector<std::vector<double>> *P_traceOfK;
  std::vector<double> *P_traceOfKNode;
};

//
//
//
void *AddClusterSiteForcesToNodesWorker(void *arg) {
  int iQuasi = ((struct AddClusterSiteForcesToNodesData_t *)arg)->iQuasi;
  std::vector<std::vector<cluster_site_data_t>> *P_clusterData =
      ((struct AddClusterSiteForcesToNodesData_t *)arg)->P_clusterData;

  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      *P_clusterForceData =
          ((struct AddClusterSiteForcesToNodesData_t *)arg)->P_clusterForceData;

  int number_threads = get_number_threads();
  int my_id = get_my_tid();
  int number_thread;
  int bucket_start;
  int bucket_end;

  get_share(my_id, number_threads, P_clusterData->size(), &number_thread,
            &bucket_start, &bucket_end);

  struct node_list_t *P_node_list =
      &(Quasicontinua::getInstance()->getQuasi(iQuasi).getNodeList().node_list);

  struct lattice_t *P_lattice =
      &(Quasicontinua::getInstance()->getQuasi(iQuasi).getLattice());

  double boltzmanC = PairPotentials::getInstance()->getBoltzmanConstant();

  double temperature = Quasicontinua::getInstance()->getTemperature();

  // get Shape Instance
  Shape *shapeC = Shape::getInstance();

  for (int iB = bucket_start; iB <= bucket_end; iB++) {
    for (int iC = 0; iC < (*P_clusterData)[iB].size(); iC++) {
      cluster_site_data_t iC_data = (*P_clusterData)[iB][iC];

      int node_err = 0;
      if (iC_data.second.size() == 5)
        node_err = -1;

      struct node_t *P_node_0 = P_node_list->nodes[iC_data.second[0]];
      struct node_t *P_node_1 = P_node_list->nodes[iC_data.second[1]];
      struct node_t *P_node_2 = P_node_list->nodes[iC_data.second[2]];
      struct node_t *P_node_3 = P_node_list->nodes[iC_data.second[3]];

      double nodeWeight;
      if (node_err == -1)
        nodeWeight = P_node_list->nodes[iC_data.second[4]]->weight;
      else
        nodeWeight = P_node_0->weight;

      //
      //  get value of shape function of node_0 at jC
      //
      int l_iC[3];
      for (int dof = 0; dof < 3; dof++)
        l_iC[dof] = iC_data.first.first[dof];

      // process node 0
      double shape = shapeC->computeSiteShapeFunction(
          P_node_0, P_node_1, P_node_2, P_node_3, l_iC, P_lattice, iQuasi);

      // update potential and kinetic energy
      if (node_err == -1) {
        struct node_t *P_node_4 = P_node_list->nodes[iC_data.second[4]];
        pthread_mutex_lock(&(P_node_4->lock));
        P_node_4->energy.potential +=
            nodeWeight * (*P_clusterForceData)[iB][iC].second;

        // P_node_4->energy.kinetic += energy_factor_kinetic * nodeWeight * 3.0
        // * boltzmanC * temperature/2.0;
        pthread_mutex_unlock(&(P_node_4->lock));
      }

      // add contribution to node_0
      pthread_mutex_lock(&(P_node_0->lock));
      P_node_0->acceleration[0] +=
          nodeWeight * (*P_clusterForceData)[iB][iC].first[0] * shape;
      P_node_0->acceleration[1] +=
          nodeWeight * (*P_clusterForceData)[iB][iC].first[1] * shape;
      P_node_0->acceleration[2] +=
          nodeWeight * (*P_clusterForceData)[iB][iC].first[2] * shape;

      P_node_0->force_frequency +=
          nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;

      if (node_err == 0) {
        P_node_0->energy.potential +=
            nodeWeight * (*P_clusterForceData)[iB][iC].second;
        // P_node_0->energy.kinetic += nodeWeight * 3.0 * boltzmanC *
        // temperature;
      }
      pthread_mutex_unlock(&(P_node_0->lock));

      // process node 1
      shape = shapeC->computeSiteShapeFunction(
          P_node_1, P_node_0, P_node_2, P_node_3, l_iC, P_lattice, iQuasi);

      // add contribution to node_1
      pthread_mutex_lock(&(P_node_1->lock));
      P_node_1->acceleration[0] +=
          nodeWeight * (*P_clusterForceData)[iB][iC].first[0] * shape;
      P_node_1->acceleration[1] +=
          nodeWeight * (*P_clusterForceData)[iB][iC].first[1] * shape;
      P_node_1->acceleration[2] +=
          nodeWeight * (*P_clusterForceData)[iB][iC].first[2] * shape;

      P_node_1->force_frequency +=
          nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;
      pthread_mutex_unlock(&(P_node_1->lock));

      // process node 2
      shape = shapeC->computeSiteShapeFunction(
          P_node_2, P_node_0, P_node_1, P_node_3, l_iC, P_lattice, iQuasi);

      // add contribution to node_2
      pthread_mutex_lock(&(P_node_2->lock));
      P_node_2->acceleration[0] +=
          nodeWeight * (*P_clusterForceData)[iB][iC].first[0] * shape;
      P_node_2->acceleration[1] +=
          nodeWeight * (*P_clusterForceData)[iB][iC].first[1] * shape;
      P_node_2->acceleration[2] +=
          nodeWeight * (*P_clusterForceData)[iB][iC].first[2] * shape;

      P_node_2->force_frequency +=
          nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;
      pthread_mutex_unlock(&(P_node_2->lock));

      // process node 3
      shape = shapeC->computeSiteShapeFunction(
          P_node_3, P_node_0, P_node_1, P_node_2, l_iC, P_lattice, iQuasi);

      // add contribution to node_3
      pthread_mutex_lock(&(P_node_3->lock));
      P_node_3->acceleration[0] +=
          nodeWeight * (*P_clusterForceData)[iB][iC].first[0] * shape;
      P_node_3->acceleration[1] +=
          nodeWeight * (*P_clusterForceData)[iB][iC].first[1] * shape;
      P_node_3->acceleration[2] +=
          nodeWeight * (*P_clusterForceData)[iB][iC].first[2] * shape;

      P_node_3->force_frequency +=
          nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;
      pthread_mutex_unlock(&(P_node_3->lock));
    }
  }

  return NULL;
}

//
//
//
void *ComputeTraceOfKAtNodesWorker(void *arg) {
  int iQuasi = ((struct ComputeTraceOfKAtNodesData_t *)arg)->iQuasi;
  std::vector<std::vector<cluster_site_data_t>> *P_clusterData =
      ((struct ComputeTraceOfKAtNodesData_t *)arg)->P_clusterData;

  std::vector<std::vector<double>> *P_traceOfK =
      ((struct ComputeTraceOfKAtNodesData_t *)arg)->P_traceOfK;
  std::vector<double> *P_traceOfKNode =
      ((struct ComputeTraceOfKAtNodesData_t *)arg)->P_traceOfKNode;

  int number_threads = get_number_threads();
  int my_id = get_my_tid();
  int number_thread;
  int bucket_start;
  int bucket_end;

  get_share(my_id, number_threads, P_clusterData->size(), &number_thread,
            &bucket_start, &bucket_end);

  struct node_list_t *P_node_list =
      &(Quasicontinua::getInstance()->getQuasi(iQuasi).getNodeList().node_list);

  struct lattice_t *P_lattice =
      &(Quasicontinua::getInstance()->getQuasi(iQuasi).getLattice());

  double boltzmanC = PairPotentials::getInstance()->getBoltzmanConstant();

  double temperature = Quasicontinua::getInstance()->getTemperature();

  // get Shape Instance
  Shape *shapeC = Shape::getInstance();

  for (int iB = bucket_start; iB <= bucket_end; iB++) {
    for (int iC = 0; iC < (*P_clusterData)[iB].size(); iC++) {
      cluster_site_data_t iC_data = (*P_clusterData)[iB][iC];

      int node_err = 0;
      if (iC_data.second.size() == 5)
        node_err = -1;

      struct node_t *P_node_0 = P_node_list->nodes[iC_data.second[0]];
      struct node_t *P_node_1 = P_node_list->nodes[iC_data.second[1]];
      struct node_t *P_node_2 = P_node_list->nodes[iC_data.second[2]];
      struct node_t *P_node_3 = P_node_list->nodes[iC_data.second[3]];

      double nodeWeight;
      if (node_err == -1)
        nodeWeight = P_node_list->nodes[iC_data.second[4]]->weight;
      else
        nodeWeight = P_node_0->weight;

      //
      //  get value of shape function of node_0 at jC
      //
      int l_iC[3];
      for (int dof = 0; dof < 3; dof++)
        l_iC[dof] = iC_data.first.first[dof];

      // process node 0
      double shape = shapeC->computeSiteShapeFunction(
          P_node_0, P_node_1, P_node_2, P_node_3, l_iC, P_lattice, iQuasi);

      // add contribution to node_0
      pthread_mutex_lock(&(P_node_0->lock));
      (*P_traceOfKNode)[iC_data.second[0]] +=
          nodeWeight * (*P_traceOfK)[iB][iC] * shape;
      pthread_mutex_unlock(&(P_node_0->lock));

      // process node 1
      shape = shapeC->computeSiteShapeFunction(
          P_node_1, P_node_0, P_node_2, P_node_3, l_iC, P_lattice, iQuasi);

      // add contribution to node_1
      pthread_mutex_lock(&(P_node_1->lock));
      (*P_traceOfKNode)[iC_data.second[1]] +=
          nodeWeight * (*P_traceOfK)[iB][iC] * shape;
      pthread_mutex_unlock(&(P_node_1->lock));

      // process node 2
      shape = shapeC->computeSiteShapeFunction(
          P_node_2, P_node_0, P_node_1, P_node_3, l_iC, P_lattice, iQuasi);

      // add contribution to node_2
      pthread_mutex_lock(&(P_node_2->lock));
      (*P_traceOfKNode)[iC_data.second[2]] +=
          nodeWeight * (*P_traceOfK)[iB][iC] * shape;
      pthread_mutex_unlock(&(P_node_2->lock));

      // process node 3
      shape = shapeC->computeSiteShapeFunction(
          P_node_3, P_node_0, P_node_1, P_node_2, l_iC, P_lattice, iQuasi);

      // add contribution to node_3
      pthread_mutex_lock(&(P_node_3->lock));
      (*P_traceOfKNode)[iC_data.second[3]] +=
          nodeWeight * (*P_traceOfK)[iB][iC] * shape;
      pthread_mutex_unlock(&(P_node_3->lock));
    }
  }

  return NULL;
}
}

//
//  namespace for EAM
//
namespace {
struct iQuasiEAMForceEnergyData_t {
  int iQuasi;
  int numQuasi;
  std::vector<std::vector<neigh_site_data_t>> *P_iQuasi_neighbor_data;
  std::vector<std::vector<std::vector<cluster_site_data_t>>> *P_cluster_data;
  std::vector<std::vector<std::vector<std::pair<std::vector<double>, double>>>>
      *P_clusterForceEnergy;
};

//
//  findSiteInVectorLocal()
//
void findSiteInVectorLocal(
    std::vector<std::pair<std::vector<int>, std::vector<double>>> *P_vectors,
    std::vector<int> site, int &loc) {
  //
  //
  //
  loc = -1;

  for (int i = 0; i < (*P_vectors).size(); i++) {
    if ((*P_vectors)[i].first[0] == site[0] &&
        (*P_vectors)[i].first[1] == site[1] &&
        (*P_vectors)[i].first[2] == site[2])
      loc = i;
  }

  return;
}

//
//  EAMGetNeighborsOfSiteLocal()
//
void EAMGetNeighborsOfSiteLocal(
    int iQuasi, std::pair<std::vector<int>, std::vector<double>> siteAndState,
    std::vector<std::vector<std::pair<std::vector<int>, std::vector<double>>>>
        &neighbors,
    int &loc_site_in_quadVector, std::vector<int> &iQuasi_quasi_map,
    std::vector<int> &iQuasi_jQuasi_interaction,
    std::vector<int> &jQuasi_jQuasi_interaction, int &number_neighbors) {
  //
  //  Task:
  //  (1) put only interacting quasis neighbor site in neighbors
  //  (2) store the mapping between counter and quasi in iQuasi_quasi_map
  //  (2) also get iQuasi - jQuasi interaction potential number
  //  (4) get jQuasi -jQuasi interaction potential number
  //  (5) get the location of site "siteAndState" in neighbors vector and
  //      correspondingly find the location of "siteAndState" in quadVectors
  //      For this, we assume the "siteAndState" will be the first element
  //      in neighbors[iQuasi] vector. Thus, it's location in quadVectors
  //      will simply be \sum_{k < iQuasi} neighbors[k].size().
  //  (6) compute the total number of neighboring site neighbors vector
  //

  //
  //  get Quasicontinua instance
  //
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  //
  //  get pairPotential instance
  //
  PairPotentials *pairC = PairPotentials::getInstance();

  //
  //  get number of quasis
  //
  int quasiSize = Quasicontinua::getInstance()->size();

  //
  //  initialization
  //
  number_neighbors = 0;
  loc_site_in_quadVector = -1;

  //
  //  loop over quasis and compute data
  //
  for (int jQuasi = 0; jQuasi < quasiSize; jQuasi++) {
    //
    //  fill the interaction data
    //
    jQuasi_jQuasi_interaction.push_back(
        pairC->doQuasicontinuumInteract(jQuasi, jQuasi));
    iQuasi_jQuasi_interaction.push_back(
        pairC->doQuasicontinuumInteract(iQuasi, jQuasi));

    //
    //  condition for loc_site_in_quadVector
    //
    if (jQuasi == iQuasi)
      loc_site_in_quadVector = number_neighbors;

    //
    //  loop over jQuasi only if interacts with iQuasi
    //
    if (iQuasi_jQuasi_interaction[jQuasi] != -1) {
      // iQuasi and jQuasi interact with each other

      // get cutoff radius for neighbor list
      double r_cut =
          pairC->getCutoffRadiusNeighborList(iQuasi_jQuasi_interaction[jQuasi]);

      neighbors.push_back(
          quasicontinua->getQuasi(jQuasi).getNeighborSitesAndStateAroundPoint(
              siteAndState, r_cut));

      //
      //  store the quasi number in map
      //
      iQuasi_quasi_map.push_back(jQuasi);

      // increment the number_neighbors
      number_neighbors += neighbors[neighbors.size() - 1].size();
    }
  }

  return;
} // end of EAMGetNeighborsOfSiteLocal()

//
//  EAMNeighborSiteProcessingLocal()
//
void EAMNeighborSiteProcessingLocal(
    int iQuasi, int numQuasi, neigh_site_data_t iNSite_data,
    std::vector<std::vector<std::vector<cluster_site_data_t>>> *P_cluster_data,
    std::vector<
        std::vector<std::vector<std::pair<std::vector<double>, double>>>>
        *P_clusterForceEnergy) {
  //
  //  get Quasicontinua handle
  //
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  //
  //  get temperature
  //
  double temperature = quasicontinua->getTemperature();

  //
  //  get vector of sigma for all quasi, where sigma = sqrt(k_B * temp)
  //  boltzman constant vector, for all quasis
  //
  //  Note : since we mainly use sqrt(2)*sigma, we get sqrt(2*k_B*temp)
  //  from Quasicontinua
  //
  std::vector<double> sigmaVector = quasicontinua->getSigmaVector();

  //
  //  get handle of PairPotentials
  //
  PairPotentials *pairC = PairPotentials::getInstance();

  //
  //  get Instance of Lattice
  //
  Lattice *latticeC = Lattice::getInstance();

  //
  //  get instance of Shape class
  //
  Shape *shapeC = Shape::getInstance();

  //
  //  handle to keep neighbors of iNSite
  //
  std::vector<std::vector<std::pair<std::vector<int>, std::vector<double>>>>
      neighbors;
  std::vector<int> iQuasi_jQuasi_interaction;
  std::vector<int> jQuasi_jQuasi_interaction;
  std::vector<int> iQuasi_quasi_map;

  int number_neighbors;
  int loc_iNSite_in_quadVector;

  //
  // call function to compute neighbors and related data
  //
  EAMGetNeighborsOfSiteLocal(iQuasi, iNSite_data.first, neighbors,
                             loc_iNSite_in_quadVector, iQuasi_quasi_map,
                             iQuasi_jQuasi_interaction,
                             jQuasi_jQuasi_interaction, number_neighbors);

  //
  //  check if we have iNSite somewhere in neighbors vector, and we have
  //  successfully found it's place in quadVectors
  //
  if (loc_iNSite_in_quadVector == -1) {
    d_print("Problem finding iNSite, in it's neighbor sites vector \n");
    D_ERROR("neighborStatesAroundSite()");
    exit(EXIT_FAILURE);
  }

  //
  //  create handle for quadrature points and get all the quadrature points
  //
  std::vector<std::pair<std::vector<std::vector<double>>, double>> quadVectors =
      QuadraturePoints::getInstance()->getQuadratureVectors(number_neighbors,
                                                            3);

  //
  //  int used in loop
  //
  int quad;

  //
  //  factor to be multiplied in phase average
  //
  double factor_quad = pow(1 / sqrt(M_PI), 1.0 * 3 * number_neighbors);
  double factor_iN = sigmaVector[iQuasi] / (iNSite_data.first.second[3]);

  // for debugging "cluster site not found in neighbors of iNSite"
  int init = 0;

  for (quad = 0; quad < quadVectors.size(); quad++) {
    //
    //  desnity and EAM related data for iNSite
    //
    double density_iN = 0.0;

    //
    //  first compute density at iN and embedding function
    //
    int vector_count = 0;
    for (int it = 0; it < iQuasi_quasi_map.size(); it++) {
      int jQuasi = iQuasi_quasi_map[it];

      for (int jN = 0; jN < neighbors[it].size(); jN++) {
        //
        //  increment the vector_count
        //
        vector_count = vector_count + 1;

        // skip when the current site jN is same as iN
        if (vector_count - 1 == loc_iNSite_in_quadVector) {
          continue;
        } else {
          //
          //  first compute r_hat(jN, jN)
          //
          std::vector<double> r_hat_iN_jN;
          r_hat_iN_jN.resize(3);

          //
          //  compute r between two atoms iN and jN
          //
          for (int dof = 0; dof < 3; dof++)
            r_hat_iN_jN[dof] =
                iNSite_data.first.second[dof] - neighbors[it][jN].second[dof];

          //
          //  get magnitude of r_hat_iN_jN
          //
          double r_hat_mag = sqrt(r_hat_iN_jN[0] * r_hat_iN_jN[0] +
                                  r_hat_iN_jN[1] * r_hat_iN_jN[1] +
                                  r_hat_iN_jN[2] * r_hat_iN_jN[2]);

          // proceed only if mean distance between iN and jN
          // is above meanSeparationTol
          if (r_hat_mag > meanSeparationTol) {
            double factor_jN =
                sigmaVector[jQuasi] / (neighbors[it][jN].second[3]);

            //
            //  compute r between two atoms iN and jN
            //
            for (int dof = 0; dof < 3; dof++)
              r_hat_iN_jN[dof] =
                  factor_iN *
                      quadVectors[quad].first[loc_iNSite_in_quadVector][dof] +
                  iNSite_data.first.second[dof] -
                  factor_jN * quadVectors[quad].first[vector_count - 1][dof] -
                  neighbors[it][jN].second[dof];

            //
            //  get magnitude of r_hat_iN_jN
            //
            double r_hat_mag = sqrt(r_hat_iN_jN[0] * r_hat_iN_jN[0] +
                                    r_hat_iN_jN[1] * r_hat_iN_jN[1] +
                                    r_hat_iN_jN[2] * r_hat_iN_jN[2]);

            //
            //  get density data from PairPotentials class
            //
            std::pair<double, double> g_dg_jN = pairC->getDensityAndDerivative(
                jQuasi_jQuasi_interaction[jQuasi], r_hat_mag);

            //
            //  add the contribution to density
            //
            density_iN += g_dg_jN.first;
          } // meanSeparationTol test
        }   // if site is not iNSite
      }     // loop over site jN
    }       // loop over counter of iQuasi_quasi_map

    //
    //  now density is computed, so we can compute
    //  embedding function and it's derivative
    //  Note : potential number is iQuasi-iQuasi interaction
    //
    std::pair<double, double> F_dF_iN =
        pairC->getEmbeddingFunctionAndDerivative(
            jQuasi_jQuasi_interaction[iQuasi], density_iN);

    //
    //  loop over all associated cluster sites of site iN
    //
    int jC;
    int jQ;
    int it;
    for (it = 0; it < iQuasi_quasi_map.size(); it++) {
      // for debugging "cluster site not found in neighbors of iNSite"
      // init = 0;
      //
      //  get quasi number
      //
      jQ = iQuasi_quasi_map[it];

      //
      //  get the node list of this quasi
      //
      struct all_node_list_t &jQuasi_all_node_list =
          quasicontinua->getQuasi(jQ).getNodeList();

      struct node_list_t &jQuasi_node_list = jQuasi_all_node_list.node_list;

      // get lattice of jQ
      struct lattice_t jQuasi_lattice =
          quasicontinua->getQuasi(jQ).getLattice();

      for (jC = 0; jC < iNSite_data.second[jQ].size(); jC++) {
        //
        //  handle to store force on jC
        //
        std::vector<double> f_jC;
        f_jC.push_back(0.0);
        f_jC.push_back(0.0);
        f_jC.push_back(0.0);
        f_jC.push_back(0.0);

        //
        //  get jC cluster site data from cluster list
        //
        //	jB is bucket number, jLoc is location in bucket
        //
        int jB = iNSite_data.second[jQ][jC][0];
        int jLoc = iNSite_data.second[jQ][jC][1];
        cluster_site_data_t jC_data = (*P_cluster_data)[jQ][jB][jLoc];

        //
        //  check if current site is cluster site itself
        //
        //  If it's cluster site, than process it
        //
        if (jQ == iQuasi &&
            jC_data.first.first[0] == iNSite_data.first.first[0] &&
            jC_data.first.first[1] == iNSite_data.first.first[1] &&
            jC_data.first.first[2] == iNSite_data.first.first[2]) {
          //
          //  we are processing cluster site jC
          //
          std::pair<double, double> F_dF_jC = F_dF_iN;

          // //
          // //  first add the energy to the associated node of jC
          // //
          // int jC_node_0 = jC_data.second[0];

          // if(jC_data.second.size() == 4)
          // {
          //   struct node_t * P_node_0 = jQuasi_node_list.nodes[jC_node_0];

          //   pthread_mutex_lock(&(P_node_0->lock));
          //   P_node_0->energy.potential +=
          //     factor_quad * F_dF_jC.first * P_node_0->weight *
          //     quadVectors[quad].second;
          //   pthread_mutex_unlock(&(P_node_0->lock));
          // }
          // else
          // {
          //   struct node_t * P_node_0 =
          //   jQuasi_node_list.nodes[jC_data.second[4]];

          //   pthread_mutex_lock(&(P_node_0->lock));
          //   P_node_0->energy.potential +=
          //     factor_quad * F_dF_jC.first * P_node_0->weight *
          //     quadVectors[quad].second;
          //   pthread_mutex_unlock(&(P_node_0->lock));
          // }

          // also add the contribution to energy into Cluster Force data
          pthread_mutex_lock(&bucketLock[jQ][jB]);
          (*P_clusterForceEnergy)[jQ][jB][jLoc].second +=
              factor_quad * F_dF_jC.first * quadVectors[quad].second;
          pthread_mutex_unlock(&bucketLock[jQ][jB]);

          //
          //  loop over neighbors of iN (or jC), both are same
          //
          int vec_counter = 0;
          int kt;
          for (kt = 0; kt < iQuasi_quasi_map.size(); kt++) {
            int kQ = iQuasi_quasi_map[kt];
            int kN;

            for (kN = 0; kN < neighbors[kt].size(); kN++) {
              //
              //  increment the counter
              //
              vec_counter = vec_counter + 1;

              if (vec_counter - 1 == loc_iNSite_in_quadVector) {
                // do nothing
                continue;
              } else {
                //
                //  first compute r_hat(jC, kN)
                //
                std::vector<double> r_hat_jC_kN;
                r_hat_jC_kN.resize(3);

                // compute distance between mean position of
                // of atom kN and jC
                for (int dof = 0; dof < 3; dof++)
                  r_hat_jC_kN[dof] =
                      jC_data.first.second[dof] - neighbors[kt][kN].second[dof];

                double r_hat_mag = sqrt(r_hat_jC_kN[0] * r_hat_jC_kN[0] +
                                        r_hat_jC_kN[1] * r_hat_jC_kN[1] +
                                        r_hat_jC_kN[2] * r_hat_jC_kN[2]);
                // continue only if
                // mean separation between two atom kN and jC
                // is above 0.3
                if (r_hat_mag > meanSeparationTol) {
                  // compute force and energy as
                  // mean separation is above tolerance

                  double factor_kN =
                      sigmaVector[kQ] / (neighbors[kt][kN].second[3]);

                  for (int dof = 0; dof < 3; dof++)
                    r_hat_jC_kN[dof] =
                        factor_iN *
                            quadVectors[quad].first[loc_iNSite_in_quadVector]
                                                   [dof] +
                        jC_data.first.second[dof] -
                        factor_kN *
                            quadVectors[quad].first[vec_counter - 1][dof] -
                        neighbors[kt][kN].second[dof];

                  //
                  //  get magnitude of r_hat
                  //
                  double r_hat_mag = sqrt(r_hat_jC_kN[0] * r_hat_jC_kN[0] +
                                          r_hat_jC_kN[1] * r_hat_jC_kN[1] +
                                          r_hat_jC_kN[2] * r_hat_jC_kN[2]);

                  //
                  //  get g and dg
                  //
                  //  Note : vector iQuasi_quasi_interaction contains potential
                  //  number
                  //  for (iQuasi, kQuasi), for kQuasi =0 to numQuasi
                  //
                  //  Note : get function "g" for atom kN, i.e. use kN's
                  //  potential number
                  //  kN belongs to kQ quasi, so use kQ potential number
                  //
                  std::pair<double, double> g_dg_kN =
                      pairC->getDensityAndDerivative(
                          jQuasi_jQuasi_interaction[kQ], r_hat_mag);

                  //
                  //  add the contribution of kN to f_jC
                  //
                  //	Note : force due to potential = - derivative of
                  //potential
                  //
                  f_jC[0] += factor_quad * (-g_dg_kN.second) * r_hat_jC_kN[0] *
                             F_dF_jC.second * quadVectors[quad].second /
                             r_hat_mag;
                  f_jC[1] += factor_quad * (-g_dg_kN.second) * r_hat_jC_kN[1] *
                             F_dF_jC.second * quadVectors[quad].second /
                             r_hat_mag;
                  f_jC[2] += factor_quad * (-g_dg_kN.second) * r_hat_jC_kN[2] *
                             F_dF_jC.second * quadVectors[quad].second /
                             r_hat_mag;

                  // frequecny force
                  f_jC[3] +=
                      factor_quad * (-g_dg_kN.second) *
                      (r_hat_jC_kN[0] *
                           quadVectors[quad].first[loc_iNSite_in_quadVector]
                                                  [0] +
                       r_hat_jC_kN[1] *
                           quadVectors[quad].first[loc_iNSite_in_quadVector]
                                                  [1] +
                       r_hat_jC_kN[2] *
                           quadVectors[quad].first[loc_iNSite_in_quadVector]
                                                  [2]) *
                      F_dF_jC.second * quadVectors[quad].second *
                      (-factor_iN / jC_data.first.second[3]) / r_hat_mag;
                } // meanSeparationTol test
              }   // if site is not iNSite
            }     // loop over kN
          }
        } // iNSite is cluster site jC
        else {
          //
          //  find the location of jC in quadVectors
          //  For this, we first need to find the location in
          //  neighbors vector
          //
          int loc_jC = 0;
          int loc_jC_quadVector = 0;
          findSiteInVectorLocal(&(neighbors[it]), jC_data.first.first, loc_jC);

          //
          //  check if we found the jC in neighbors
          //
          if (loc_jC == -1) {
            pthread_mutex_lock(&lock_print);
            d_print("cluster site = (%d %d, %d) of Quasi %d is not found in "
                    "neighbors of site = (%i, %i, %i) of Quasi = %i \n",
                    jC_data.first.first[0], jC_data.first.first[1],
                    jC_data.first.first[2], jQ, iNSite_data.first.first[0],
                    iNSite_data.first.first[1], iNSite_data.first.first[2],
                    iQuasi);
            d_print("*********************** dumping neighbors of iNSite "
                    "******************\n ");
            for (int i = 0; i < neighbors.size(); i++) {
              d_print("						Quasi = %d\n",
                      i);
              for (int j = 0; j < neighbors[i].size(); j++) {
                d_print("(%d, %d, %d)\n", neighbors[i][j].first[0],
                        neighbors[i][j].first[1], neighbors[i][j].first[2]);
              }
            }

            d_print("**************** dumping associated cluster site of "
                    "iNSite *************\n");
            for (int i = 0; i < iNSite_data.second.size(); i++) {
              d_print("				Quasi = %d\n", i);
              for (int j = 0; j < iNSite_data.second[i].size(); j++) {
                d_print("(%d, %d, %d)\n",
                        (*P_cluster_data)[i][iNSite_data.second[i][j][0]]
                                         [iNSite_data.second[i][j][1]]
                                             .first.first[0],
                        (*P_cluster_data)[i][iNSite_data.second[i][j][0]]
                                         [iNSite_data.second[i][j][1]]
                                             .first.first[1],
                        (*P_cluster_data)[i][iNSite_data.second[i][j][0]]
                                         [iNSite_data.second[i][j][1]]
                                             .first.first[2]);
              }
            }

            d_print("*********************** size of data "
                    "***********************\n");
            d_print("number of cluster sites associated to iNSite for iQuasi "
                    "%d is = %d \n",
                    jQ, iNSite_data.second[jQ].size());
            d_print("number of neighbors around iNSite = %d\n",
                    neighbors[it].size());
            pthread_mutex_unlock(&lock_print);
            // exit(EXIT_FAILURE);
            continue;
          }

          for (int temp_it = 0; temp_it < it; temp_it++)
            loc_jC_quadVector += neighbors[temp_it].size();

          loc_jC_quadVector += loc_jC;

          //
          //  compute r_hat_jC_iN
          //
          std::vector<double> r_hat_jC_iN;
          r_hat_jC_iN.resize(3);

          // compute distance between mean position of
          // of atom iN and jC
          for (int dof = 0; dof < 3; dof++)
            r_hat_jC_iN[dof] =
                jC_data.first.second[dof] - iNSite_data.first.second[dof];

          //
          //  compute r_hat_mag
          //
          double r_hat_mag = sqrt(r_hat_jC_iN[0] * r_hat_jC_iN[0] +
                                  r_hat_jC_iN[1] * r_hat_jC_iN[1] +
                                  r_hat_jC_iN[2] * r_hat_jC_iN[2]);

          // continue only if
          // mean separation between two atom kN and jC
          // is above 0.3
          if (r_hat_mag > meanSeparationTol) {

            double factor_jC = sigmaVector[jQ] / (jC_data.first.second[3]);

            for (int dof = 0; dof < 3; dof++)
              r_hat_jC_iN[dof] =
                  factor_jC * quadVectors[quad].first[loc_jC_quadVector][dof] +
                  jC_data.first.second[dof] -
                  factor_iN *
                      quadVectors[quad].first[loc_iNSite_in_quadVector][dof] -
                  iNSite_data.first.second[dof];

            //
            //  compute r_hat_mag
            //
            double r_hat_mag = sqrt(r_hat_jC_iN[0] * r_hat_jC_iN[0] +
                                    r_hat_jC_iN[1] * r_hat_jC_iN[1] +
                                    r_hat_jC_iN[2] * r_hat_jC_iN[2]);

            //
            //  compute function "g" for atom jC
            //
            std::pair<double, double> g_dg_jC = pairC->getDensityAndDerivative(
                jQuasi_jQuasi_interaction[jQ], r_hat_mag);

            //
            //  add the contribution to force
            //
            f_jC[0] += factor_quad * (-g_dg_jC.second) * r_hat_jC_iN[0] *
                       F_dF_iN.second * quadVectors[quad].second / r_hat_mag;
            f_jC[1] += factor_quad * (-g_dg_jC.second) * r_hat_jC_iN[1] *
                       F_dF_iN.second * quadVectors[quad].second / r_hat_mag;
            f_jC[2] += factor_quad * (-g_dg_jC.second) * r_hat_jC_iN[2] *
                       F_dF_iN.second * quadVectors[quad].second / r_hat_mag;

            // frequecny force
            f_jC[3] += factor_quad * (-g_dg_jC.second) *
                       (r_hat_jC_iN[0] *
                            quadVectors[quad].first[loc_jC_quadVector][0] +
                        r_hat_jC_iN[1] *
                            quadVectors[quad].first[loc_jC_quadVector][1] +
                        r_hat_jC_iN[2] *
                            quadVectors[quad].first[loc_jC_quadVector][2]) *
                       F_dF_iN.second * quadVectors[quad].second *
                       (-factor_jC / jC_data.first.second[3]) / r_hat_mag;
          } // meanSeparationTol test
        }   // iNSite is not jC cluster site

        // add force to Force Data
        pthread_mutex_lock(&bucketLock[jQ][jB]);
        (*P_clusterForceEnergy)[jQ][jB][jLoc].first[0] += f_jC[0];
        (*P_clusterForceEnergy)[jQ][jB][jLoc].first[1] += f_jC[1];
        (*P_clusterForceEnergy)[jQ][jB][jLoc].first[2] += f_jC[2];
        (*P_clusterForceEnergy)[jQ][jB][jLoc].first[3] += f_jC[3];
        pthread_mutex_unlock(&bucketLock[jQ][jB]);
      } // loop over jC
    }   // loop over jQ
  }     // loop over quad

  return;
} // end of EAMNeighborSiteProcessingLocal()

//
//  EAMNeighborBucketProcessingLocal()
//
void EAMNeighborBucketProcessingLocal(
    int iQuasi, int numQuasi,
    std::vector<std::vector<neigh_site_data_t>> *P_iQuasi_neighbor_data,
    std::vector<std::vector<std::vector<cluster_site_data_t>>> *P_cluster_data,
    int bucket_start, int bucket_end,
    std::vector<
        std::vector<std::vector<std::pair<std::vector<double>, double>>>>
        *P_clusterForceEnergy) {
  int iBucket;
  int iNSite;

  for (iBucket = bucket_start; iBucket <= bucket_end; iBucket++) {
    for (iNSite = 0; iNSite < (*P_iQuasi_neighbor_data)[iBucket].size();
         iNSite++) {
      //
      //  call function to process iNSite
      //
      EAMNeighborSiteProcessingLocal(iQuasi, numQuasi,
                                     (*P_iQuasi_neighbor_data)[iBucket][iNSite],
                                     P_cluster_data, P_clusterForceEnergy);
    } // loop over site
  }   // loop over bucket

  return;
} // end of EAMNeighborBucketProcessingLocal()

//
//  iQuasiEAMForceEnergyWorker()
//
void *iQuasiEAMForceEnergyWorker(void *arg) {
  int iQuasi = ((struct iQuasiEAMForceEnergyData_t *)arg)->iQuasi;
  int numQuasi = ((struct iQuasiEAMForceEnergyData_t *)arg)->numQuasi;
  std::vector<std::vector<neigh_site_data_t>> *P_iQuasi_neighbor_data =
      ((struct iQuasiEAMForceEnergyData_t *)arg)->P_iQuasi_neighbor_data;
  std::vector<std::vector<std::vector<cluster_site_data_t>>> *P_cluster_data =
      ((struct iQuasiEAMForceEnergyData_t *)arg)->P_cluster_data;

  std::vector<std::vector<std::vector<std::pair<std::vector<double>, double>>>>
      *P_clusterForceEnergy =
          ((struct iQuasiEAMForceEnergyData_t *)arg)->P_clusterForceEnergy;

  int number_threads = get_number_threads();
  int my_id = get_my_tid();
  int number_thread;
  int bucket_start;
  int bucket_end;

  get_share(my_id, number_threads, P_iQuasi_neighbor_data->size(),
            &number_thread, &bucket_start, &bucket_end);

  EAMNeighborBucketProcessingLocal(iQuasi, numQuasi, P_iQuasi_neighbor_data,
                                   P_cluster_data, bucket_start, bucket_end,
                                   P_clusterForceEnergy);

  return ((void *)NULL);
} // end of iQuasiEAMForceEnergyWorker()
} // end of namespace for EAM

//
//  namespace for Pairwise
//
namespace {
struct iQuasiPairwiseForceEnergyData_t {
  int iQuasi;
  int numQuasi;
  std::pair<int, int> coreFlag;
  std::vector<std::vector<cluster_site_data_t>> *P_iQuasi_cluster_data;
  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      *P_clusterForceEnergy;
};

//
//  PairwiseGetNeighborsOfSiteLocal()
//
void PairwiseGetNeighborsOfSiteLocal(
    int iQuasi, std::pair<std::vector<int>, std::vector<double>> siteAndState,
    std::vector<std::vector<std::vector<double>>> &neighbors,
    std::vector<int> &iQuasi_quasi_map,
    std::vector<int> &iQuasi_jQuasi_interaction,
    std::vector<int> &jQuasi_jQuasi_interaction) {
  //
  //  Task:
  //  (1) put only interacting quasis neighbor site in neighbors
  //  (2) store the mapping between counter and quasi in iQuasi_quasi_map
  //  (2) also get iQuasi - jQuasi interaction potential number
  //  (4) get jQuasi -jQuasi interaction potential number
  //

  //
  //  get Quasicontinua instance
  //
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  //
  //  get pairPotential instance
  //
  PairPotentials *pairC = PairPotentials::getInstance();

  //
  //  get number of quasis
  //
  int numQuasi = quasicontinua->size();

  //
  //  loop over quasis and compute data
  //
  for (int jQuasi = 0; jQuasi < numQuasi; jQuasi++) {
    //
    //  fill the interaction data
    //
    jQuasi_jQuasi_interaction.push_back(
        pairC->doQuasicontinuumInteract(jQuasi, jQuasi));
    iQuasi_jQuasi_interaction.push_back(
        pairC->doQuasicontinuumInteract(iQuasi, jQuasi));

    //
    //  loop over jQuasi only if interacts with iQuasi
    //
    if (iQuasi_jQuasi_interaction[jQuasi] != -1) {
      // iQuasi and jQuasi interact with each other

      //  get cutoff radius
      double r_cut =
          pairC->getCutoffRadiusNeighborList(iQuasi_jQuasi_interaction[jQuasi]);

      neighbors.push_back(
          quasicontinua->getQuasi(jQuasi).getNeighborStateAroundPoint(
              siteAndState, r_cut));

      //
      //  store the quasi number in map
      //
      iQuasi_quasi_map.push_back(jQuasi);
    }
  }

  return;
} // end of PairwiseGetNeighborsOfSiteLocal()

//
//  PairwiseClusterSiteProcessingLocal()
//
void PairwiseShellQuasiClusterSiteProcessingLocal(
    int iQuasi, int coreQuasi, cluster_site_data_t iCSite_data, int iC_bucket,
    int iC_location,
    std::vector<std::vector<std::pair<std::vector<double>, double>>>
        *P_clusterForceEnergy) {
  //
  //  force and energy at cluster site due to pairwise interaction
  //
  std::vector<double> f_iC;
  double energy_iC = 0.0;
  f_iC.push_back(0.0);
  f_iC.push_back(0.0);
  f_iC.push_back(0.0);
  f_iC.push_back(0.0);

  //
  //  get Quasicontinua handle
  //
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  //
  //  get PairPotentials instance
  //
  PairPotentials *pairC = PairPotentials::getInstance();

  //
  //  get Instance of Lattice, Shape
  //
  Lattice *latticeC = Lattice::getInstance();
  Shape *shapeC = Shape::getInstance();

  //
  //  get Sigmavectors
  //
  std::vector<double> sigmaVector = quasicontinua->getSigmaVector();

  //
  // handle for neighbor sites
  //
  std::vector<std::vector<std::vector<double>>> neighbors;

  //
  //  handle foe interaction data and map
  //
  std::vector<int> iQuasi_quasi_map;
  std::vector<int> iQuasi_jQuasi_interaction;
  std::vector<int> jQuasi_jQuasi_interaction;

  //
  //  call function to get neighbors
  //
  PairwiseGetNeighborsOfSiteLocal(iQuasi, iCSite_data.first, neighbors,
                                  iQuasi_quasi_map, iQuasi_jQuasi_interaction,
                                  jQuasi_jQuasi_interaction);

  //
  //  get quadrature points for pairwise interaction
  //  Note : For pairwise, we have phase average integration over
  //         R^6 space, i.e. 2 vectors each on R^3 space.
  //
  std::vector<std::pair<std::vector<std::vector<double>>, double>> quadVectors =
      QuadraturePoints::getInstance()->getQuadratureVectors(2, 3);

  //
  int quad;

  //
  //  factor_quad = integration factor
  //  factor_iC = sqrt(2) /sigma / frequency for ICSite
  //
  double factor_quad = pow(1 / sqrt(M_PI), 1.0 * 3 * 2);
  double factor_iC = sigmaVector[iQuasi] / (iCSite_data.first.second[3]);

  //
  //
  //
  int it_coreQuasi = -1;

  //
  //  loop over quadrature points
  //  Note : first vector in R^6 quad space is corresponding to iCSite,
  //         and second vector is corresponding to neighbor site
  //
  //	Note : We will deal with iQuasi-coreQuasi interaction separately
  //				 Phase average is for interaction of shell-shell
  //type
  //
  for (quad = 0; quad < quadVectors.size(); quad++) {
    //
    //  get two vectors and weight
    //
    std::vector<double> quadV_1 = quadVectors[quad].first[0];
    std::vector<double> quadV_2 = quadVectors[quad].first[1];
    double quadWeight = quadVectors[quad].second;

    //
    //  loop over all interacting quasis
    //
    int it;
    int jQ;
    int jN;
    for (it = 0; it < iQuasi_quasi_map.size(); it++) {
      //
      //  get Quasi number
      //
      jQ = iQuasi_quasi_map[it];

      //
      // proceed only if jQ is not equal to coreQuasi
      //
      // if core quasi of shell is -1, meaning there is no core of this shell
      // then no need to check
      // else check
      if (coreQuasi != -1) {
        if (jQ == coreQuasi) {
          it_coreQuasi = it;
        }
      }

      if (jQ != coreQuasi) {
        //
        //  loop over neighboring site of jQ
        //
        for (jN = 0; jN < neighbors[it].size(); jN++) {
          //
          //  get the state of jN
          //
          std::vector<double> jN_state = neighbors[it][jN];

          //
          // get factor corresponding to jNSite
          //
          double factor_jN = sigmaVector[jQ] / jN_state[3];

          //
          //  compute r_hat_iC_jN and magnitude
          //
          std::vector<double> r_hat_iC_jN;
          double r_hat_mag = 0.0;

          r_hat_iC_jN.resize(3);

          //
          //  first get the difference of mean position of two sites
          //  If mean position difference is less than half of lattice size
          //  the atoms are most probably same atoms.
          //  For such atoms we will not proceed further
          //
          for (int dof = 0; dof < 3; dof++)
            r_hat_iC_jN[dof] = iCSite_data.first.second[dof] - jN_state[dof];

          r_hat_mag = sqrt(r_hat_iC_jN[0] * r_hat_iC_jN[0] +
                           r_hat_iC_jN[1] * r_hat_iC_jN[1] +
                           r_hat_iC_jN[2] * r_hat_iC_jN[2]);

          // check is  r_hat_mag above meanSeparationTol
          if (r_hat_mag > meanSeparationTol) {
            for (int dof = 0; dof < 3; dof++)
              r_hat_iC_jN[dof] = r_hat_iC_jN[dof] + factor_iC * quadV_1[dof] -
                                 factor_jN * quadV_2[dof];

            r_hat_mag = sqrt(r_hat_iC_jN[0] * r_hat_iC_jN[0] +
                             r_hat_iC_jN[1] * r_hat_iC_jN[1] +
                             r_hat_iC_jN[2] * r_hat_iC_jN[2]);

            //
            //  compute phi and derivative of phi at r_hat_mag
            //
            std::pair<double, double> dphi_phi = pairC->getPairForceAndEnergy(
                iQuasi_jQuasi_interaction[jQ], r_hat_mag);

            //
            //  add contribution to energy
            //
            energy_iC += 0.5 * factor_quad * dphi_phi.second * quadWeight;

            //
            //  add contribution to force
            //
            f_iC[0] += factor_quad * (dphi_phi.first) * r_hat_iC_jN[0] *
                       quadWeight / r_hat_mag;

            f_iC[1] += factor_quad * (dphi_phi.first) * r_hat_iC_jN[1] *
                       quadWeight / r_hat_mag;

            f_iC[2] += factor_quad * (dphi_phi.first) * r_hat_iC_jN[2] *
                       quadWeight / r_hat_mag;

            f_iC[3] +=
                factor_quad * (dphi_phi.first) *
                (r_hat_iC_jN[0] * quadV_1[0] + r_hat_iC_jN[1] * quadV_1[1] +
                 r_hat_iC_jN[2] * quadV_1[2]) *
                (-factor_iC / iCSite_data.first.second[3]) * quadWeight /
                r_hat_mag;
          } // meanSeparationTol test
        }   // jQ is not the core of shell iQuasi
      }     // loop over jN of jQ
    }       // loop over counter of iQuasi_quasi_map
  }         // loop over quad vectors

  //
  //	now process iQuasi-coreQuasi interaction
  //	Note : Verify first if indeed coreQuasi is core of shell
  if (coreQuasi != -1) {
    if (it_coreQuasi == -1) {
      // it should not be -1 if coreQuasi is not -1
      pthread_mutex_lock(&lock_print);
      d_print("error in finding coreQuasi\n");
      pthread_mutex_unlock(&lock_print);

      exit(EXIT_FAILURE);
    } else {
      // get neighbors of cluster site of iQuasi in coreQuasi
      for (int jN = 0; jN < neighbors[it_coreQuasi].size(); jN++) {
        // get state of jN
        std::vector<double> jN_state = neighbors[it_coreQuasi][jN];

        // find r and r_mag
        std::vector<double> r_hat_iC_jN;
        double r_hat_mag = 0.0;

        for (int dof = 0; dof < 3; dof++)
          r_hat_iC_jN.push_back(iCSite_data.first.second[dof] - jN_state[dof]);

        r_hat_mag = sqrt(r_hat_iC_jN[0] * r_hat_iC_jN[0] +
                         r_hat_iC_jN[1] * r_hat_iC_jN[1] +
                         r_hat_iC_jN[2] * r_hat_iC_jN[2]);

        // find the force and energy from PairPotentials
        // pthread_mutex_lock(&lock_print);
        // d_print("shell Quasi = %d, core Quasi = %d, potential number = %d\n",
        // iQuasi, coreQuasi, iQuasi_jQuasi_interaction[coreQuasi]);
        // pthread_mutex_unlock(&lock_print);

        //
        //	For shell-core interaction : do not test
        //	for meanSeparationTol
        //
        std::pair<double, double> dphi_phi = pairC->getPairForceAndEnergy(
            iQuasi_jQuasi_interaction[coreQuasi], r_hat_mag);

        // energy
        energy_iC += 0.5 * dphi_phi.second;

        // force
        f_iC[0] += (dphi_phi.first) * r_hat_iC_jN[0] / r_hat_mag;
        f_iC[1] += (dphi_phi.first) * r_hat_iC_jN[1] / r_hat_mag;
        f_iC[2] += (dphi_phi.first) * r_hat_iC_jN[2] / r_hat_mag;
      } // loop over jN
    }
  }

  //
  //  Note : Force is -f_iC
  //
  for (int dof = 0; dof < 4; dof++)
    f_iC[dof] = -f_iC[dof];

  //
  //  update the cluster site force energy data
  //
  pthread_mutex_lock(&bucketLock[iQuasi][iC_bucket]);
  (*P_clusterForceEnergy)[iC_bucket][iC_location].first[0] += f_iC[0];
  (*P_clusterForceEnergy)[iC_bucket][iC_location].first[1] += f_iC[1];
  (*P_clusterForceEnergy)[iC_bucket][iC_location].first[2] += f_iC[2];
  (*P_clusterForceEnergy)[iC_bucket][iC_location].first[3] += f_iC[3];

  (*P_clusterForceEnergy)[iC_bucket][iC_location].second += energy_iC;
  pthread_mutex_unlock(&bucketLock[iQuasi][iC_bucket]);

  return;
} // end of PairwiseClusterSiteProcessingLocal()

//
//  PairwiseClusterSiteProcessingLocal()
//
void PairwiseCoreQuasiClusterSiteProcessingLocal(
    int iQuasi, int coreQuasi, cluster_site_data_t iCSite_data, int iC_bucket,
    int iC_location,
    std::vector<std::vector<std::pair<std::vector<double>, double>>>
        *P_clusterForceEnergy) {
  //
  //	define local force and energy variables
  //
  double energy_iC = 0.0;
  std::vector<double> f_iC;
  f_iC.push_back(0.0);
  f_iC.push_back(0.0);
  f_iC.push_back(0.0);
  f_iC.push_back(0.0);

  //
  //	get the neighbors of iC in quasi coreQuasi
  //
  std::vector<std::vector<double>> neighbors;

  PairPotentials *pairC = PairPotentials::getInstance();
  int potentialNumber = pairC->doQuasicontinuumInteract(iQuasi, coreQuasi);

  if (potentialNumber != -1) {
    // get cutoff radius
    double r_cut = pairC->getCutoffRadiusNeighborList(potentialNumber);
    neighbors = Quasicontinua::getInstance()
                    ->getQuasi(coreQuasi)
                    .getNeighborStateAroundPoint(iCSite_data.first, r_cut);

    // loop over neighbor sites
    for (int jN = 0; jN < neighbors.size(); jN++) {
      // find r and r_mag
      std::vector<double> r;
      double r_mag = 0.0;

      for (int dof = 0; dof < 3; dof++)
        r.push_back(iCSite_data.first.second[dof] - neighbors[jN][dof]);

      r_mag = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

      // find the force and energy from PairPotentials
      std::pair<double, double> dphi_phi =
          pairC->getPairForceAndEnergy(potentialNumber, r_mag);

      // energy
      energy_iC += 0.5 * dphi_phi.second;

      // force
      f_iC[0] += (dphi_phi.first) * r[0] / r_mag;
      f_iC[1] += (dphi_phi.first) * r[1] / r_mag;
      f_iC[2] += (dphi_phi.first) * r[2] / r_mag;
    }
  }
  //
  //  Note : Force is -f_iC
  //
  for (int dof = 0; dof < 4; dof++)
    f_iC[dof] = -f_iC[dof];

  //
  //  update the cluster site force energy data
  //
  pthread_mutex_lock(&bucketLock[iQuasi][iC_bucket]);
  (*P_clusterForceEnergy)[iC_bucket][iC_location].first[0] += f_iC[0];
  (*P_clusterForceEnergy)[iC_bucket][iC_location].first[1] += f_iC[1];
  (*P_clusterForceEnergy)[iC_bucket][iC_location].first[2] += f_iC[2];
  (*P_clusterForceEnergy)[iC_bucket][iC_location].first[3] += f_iC[3];

  (*P_clusterForceEnergy)[iC_bucket][iC_location].second += energy_iC;
  pthread_mutex_unlock(&bucketLock[iQuasi][iC_bucket]);
}

//
//  PairwiseShellQuasiClusterBucketProcessingLocal()
//
void PairwiseClusterBucketProcessingLocal(
    int iQuasi, int numQuasi, std::pair<int, int> coreFlag,
    std::vector<std::vector<cluster_site_data_t>> *P_iQuasi_cluster_data,
    int bucket_start, int bucket_end,
    std::vector<std::vector<std::pair<std::vector<double>, double>>>
        *P_clusterForceEnergy) {
  //
  //  loop over bucket
  //
  int iBucket;
  int iCSite;
  for (iBucket = bucket_start; iBucket <= bucket_end; iBucket++) {
    for (iCSite = 0; iCSite < (*P_iQuasi_cluster_data)[iBucket].size();
         iCSite++) {
      //
      //  process cluster site
      //
      cluster_site_data_t iCSite_data =
          (*P_iQuasi_cluster_data)[iBucket][iCSite];

      //
      //  call function to process current cluster site of iQuasi
      //
      if (coreFlag.first == -1)
        PairwiseShellQuasiClusterSiteProcessingLocal(
            iQuasi, coreFlag.second, iCSite_data, iBucket, iCSite,
            P_clusterForceEnergy);
      else
        PairwiseCoreQuasiClusterSiteProcessingLocal(
            iQuasi, coreFlag.second, iCSite_data, iBucket, iCSite,
            P_clusterForceEnergy);
    } // loop over cluster site
  }   // loop over bucket

  return;
} // end of PairwiseClusterBucketProcessingLocal()

//
//  iQuasiPairwiseForceEnergyWorker()
//
void *iQuasiPairwiseForceEnergyWorker(void *arg) {
  int iQuasi = ((struct iQuasiPairwiseForceEnergyData_t *)arg)->iQuasi;
  int numQuasi = ((struct iQuasiPairwiseForceEnergyData_t *)arg)->numQuasi;
  std::pair<int, int> coreFlag =
      ((struct iQuasiPairwiseForceEnergyData_t *)arg)->coreFlag;
  std::vector<std::vector<cluster_site_data_t>> *P_iQuasi_cluster_data =
      ((struct iQuasiPairwiseForceEnergyData_t *)arg)->P_iQuasi_cluster_data;
  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      *P_clusterForceEnergy =
          ((struct iQuasiPairwiseForceEnergyData_t *)arg)->P_clusterForceEnergy;

  int number_threads = get_number_threads();
  int my_id = get_my_tid();
  int number_thread;
  int bucket_start;
  int bucket_end;

  get_share(my_id, number_threads, P_iQuasi_cluster_data->size(),
            &number_thread, &bucket_start, &bucket_end);

  // pthread_mutex_lock(&lock_print);
  // d_print("iQuasiPairwiseForceEnergy()\n");
  // pthread_mutex_unlock(&lock_print);

  PairwiseClusterBucketProcessingLocal(iQuasi, numQuasi, coreFlag,
                                       P_iQuasi_cluster_data, bucket_start,
                                       bucket_end, P_clusterForceEnergy);

  return ((void *)NULL);
} // end of iQuasiPairwiseForceEnergyWorker()
} // end of namespace for Pairwise

//
//  namespace for Pairwise - Harmonin approximation
//
namespace {
struct iQuasiPairwiseForceEnergyHarmonicData_t {
  int iQuasi;
  int numQuasi;
  std::pair<int, int> coreFlag;
  std::vector<std::vector<cluster_site_data_t>> *P_iQuasi_cluster_data;
  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      *P_clusterForceEnergy;
  std::vector<std::vector<double>> *P_traceOfK;
};

//
//  PairwiseClusterSiteProcessingLocal()
//
void PairwiseHarmonicShellQuasiClusterSiteProcessingLocal(
    int iQuasi, int coreQuasi, cluster_site_data_t iCSite_data, int iC_bucket,
    int iC_location,
    std::vector<std::vector<std::pair<std::vector<double>, double>>>
        *P_clusterForceEnergy,
    std::vector<std::vector<double>> *P_traceOfK) {
  //
  //  force and energy at cluster site due to pairwise interaction
  //
  std::vector<double> f_iC;
  double energy_iC = 0.0;
  f_iC.push_back(0.0);
  f_iC.push_back(0.0);
  f_iC.push_back(0.0);
  f_iC.push_back(0.0);

  double traceOfK_iC = 0.0;

  //
  //  get Quasicontinua handle
  //
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  //
  //  get PairPotentials instance
  //
  PairPotentials *pairC = PairPotentials::getInstance();

  //
  //  get Instance of Lattice, Shape
  //
  Lattice *latticeC = Lattice::getInstance();
  Shape *shapeC = Shape::getInstance();

  //
  //  get Sigmavectors
  //
  std::vector<double> sigmaVector = quasicontinua->getSigmaVector();

  //
  // handle for neighbor sites
  //
  std::vector<std::vector<std::vector<double>>> neighbors;

  //
  //  handle foe interaction data and map
  //
  std::vector<int> iQuasi_quasi_map;
  std::vector<int> iQuasi_jQuasi_interaction;
  std::vector<int> jQuasi_jQuasi_interaction;

  //
  //  call function to get neighbors
  //
  PairwiseGetNeighborsOfSiteLocal(iQuasi, iCSite_data.first, neighbors,
                                  iQuasi_quasi_map, iQuasi_jQuasi_interaction,
                                  jQuasi_jQuasi_interaction);

  //
  //  loop over all interacting quasis
  //
  int it;
  int jQ;
  int jN;
  for (it = 0; it < iQuasi_quasi_map.size(); it++) {
    //
    //  get Quasi number
    //
    jQ = iQuasi_quasi_map[it];

    //
    //  loop over neighboring site of jQ
    //
    for (jN = 0; jN < neighbors[it].size(); jN++) {
      //
      //  get the state of jN
      //
      std::vector<double> jN_state = neighbors[it][jN];

      //
      //  compute r_hat_iC_jN and magnitude
      //
      std::vector<double> r_hat_iC_jN;
      double r_hat_mag = 0.0;

      r_hat_iC_jN.resize(3);
      for (int dof = 0; dof < 3; dof++)
        r_hat_iC_jN[dof] = iCSite_data.first.second[dof] - jN_state[dof];

      r_hat_mag = sqrt(r_hat_iC_jN[0] * r_hat_iC_jN[0] +
                       r_hat_iC_jN[1] * r_hat_iC_jN[1] +
                       r_hat_iC_jN[2] * r_hat_iC_jN[2]);

      // check if r_hat_mag is above meanSeparationTol
      if (r_hat_mag > meanSeparationTol) {
        //
        //  compute phi and derivative of phi at r_hat_mag
        //
        double phi = 0.0;
        double dphi = 0.0;
        double d_dphi = 0.0;

        pairC->getHarmonicApproxData(phi, dphi, d_dphi,
                                     iQuasi_jQuasi_interaction[jQ], r_hat_mag);

        //
        //  add contribution to energy
        //
        energy_iC += 0.5 * phi;

        //
        //	add contribution to traceOfK
        //
        traceOfK_iC += d_dphi + 2.0 * dphi / (r_hat_mag);

        // pthread_mutex_lock(&bucketLock[iQuasi][iC_bucket]);
        // if(dummy < 100)
        // {
        // 	d_print("d_dphi + 2.0*dphi/(r_hat_mag) = %f\n", d_dphi +
        // 2.0*dphi/(r_hat_mag));
        // 	dummy +=1;
        // }
        // pthread_mutex_unlock(&bucketLock[iQuasi][iC_bucket]);

        //
        //  add contribution to force
        //
        f_iC[0] += dphi * r_hat_iC_jN[0] / r_hat_mag;
        f_iC[1] += dphi * r_hat_iC_jN[1] / r_hat_mag;
        f_iC[2] += dphi * r_hat_iC_jN[2] / r_hat_mag;
      } // meanSeparationTol test
    }   // loop over jN of jQ
  }     // loop over counter of iQuasi_quasi_map

  //
  //  Note : Force is -f_iC
  //
  for (int dof = 0; dof < 4; dof++)
    f_iC[dof] = -f_iC[dof];

  //
  //  update the cluster site force energy data
  //
  pthread_mutex_lock(&bucketLock[iQuasi][iC_bucket]);
  (*P_clusterForceEnergy)[iC_bucket][iC_location].first[0] += f_iC[0];
  (*P_clusterForceEnergy)[iC_bucket][iC_location].first[1] += f_iC[1];
  (*P_clusterForceEnergy)[iC_bucket][iC_location].first[2] += f_iC[2];

  (*P_clusterForceEnergy)[iC_bucket][iC_location].second += energy_iC;

  (*P_traceOfK)[iC_bucket][iC_location] += traceOfK_iC;
  // d_print("%f\n", traceOfK_iC);
  pthread_mutex_unlock(&bucketLock[iQuasi][iC_bucket]);

  return;
} // end of PairwiseClusterSiteProcessingLocal()

//
//  PairwiseClusterSiteProcessingLocal()
//
void PairwiseHarmonicCoreQuasiClusterSiteProcessingLocal(
    int iQuasi, int coreQuasi, cluster_site_data_t iCSite_data, int iC_bucket,
    int iC_location,
    std::vector<std::vector<std::pair<std::vector<double>, double>>>
        *P_clusterForceEnergy,
    std::vector<std::vector<double>> *P_traceOfK) {
  //
  //	define local force and energy variables
  //
  double energy_iC = 0.0;
  std::vector<double> f_iC;
  f_iC.push_back(0.0);
  f_iC.push_back(0.0);
  f_iC.push_back(0.0);
  f_iC.push_back(0.0);

  //
  //	get the neighbors of iC in quasi coreQuasi
  //
  std::vector<std::vector<double>> neighbors;

  PairPotentials *pairC = PairPotentials::getInstance();
  int potentialNumber = pairC->doQuasicontinuumInteract(iQuasi, coreQuasi);

  if (potentialNumber != -1) {
    // get cutoff radius
    double r_cut = pairC->getCutoffRadiusNeighborList(potentialNumber);
    neighbors = Quasicontinua::getInstance()
                    ->getQuasi(coreQuasi)
                    .getNeighborStateAroundPoint(iCSite_data.first, r_cut);

    // loop over neighbor sites
    for (int jN = 0; jN < neighbors.size(); jN++) {
      // find r and r_mag
      std::vector<double> r;
      double r_mag = 0.0;

      for (int dof = 0; dof < 3; dof++)
        r.push_back(iCSite_data.first.second[dof] - neighbors[jN][dof]);

      r_mag = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

      // find the force and energy from PairPotentials
      std::pair<double, double> dphi_phi =
          pairC->getPairForceAndEnergy(potentialNumber, r_mag);

      // energy
      energy_iC += 0.5 * dphi_phi.second;

      // force
      f_iC[0] += dphi_phi.first * r[0] / r_mag;
      f_iC[1] += dphi_phi.first * r[1] / r_mag;
      f_iC[2] += dphi_phi.first * r[2] / r_mag;
    }
  }

  //
  //  Note : Force is -f_iC
  //
  for (int dof = 0; dof < 4; dof++)
    f_iC[dof] = -f_iC[dof];

  //
  //  update the cluster site force energy data
  //
  pthread_mutex_lock(&bucketLock[iQuasi][iC_bucket]);
  (*P_clusterForceEnergy)[iC_bucket][iC_location].first[0] += f_iC[0];
  (*P_clusterForceEnergy)[iC_bucket][iC_location].first[1] += f_iC[1];
  (*P_clusterForceEnergy)[iC_bucket][iC_location].first[2] += f_iC[2];
  (*P_clusterForceEnergy)[iC_bucket][iC_location].first[3] += f_iC[3];

  (*P_clusterForceEnergy)[iC_bucket][iC_location].second += energy_iC;
  pthread_mutex_unlock(&bucketLock[iQuasi][iC_bucket]);
}

//
//  PairwiseShellQuasiClusterBucketProcessingLocal()
//
void PairwiseHarmonicClusterBucketProcessingLocal(
    int iQuasi, int numQuasi, std::pair<int, int> coreFlag,
    std::vector<std::vector<cluster_site_data_t>> *P_iQuasi_cluster_data,
    int bucket_start, int bucket_end,
    std::vector<std::vector<std::pair<std::vector<double>, double>>>
        *P_clusterForceEnergy,
    std::vector<std::vector<double>> *P_traceOfK) {
  //
  //  get Quasicontinua, PairPotentials
  //
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  PairPotentials *pairC = PairPotentials::getInstance();

  //
  //  loop over bucket
  //
  int iBucket;
  int iCSite;
  for (iBucket = bucket_start; iBucket <= bucket_end; iBucket++) {
    for (iCSite = 0; iCSite < (*P_iQuasi_cluster_data)[iBucket].size();
         iCSite++) {
      //
      //  process cluster site
      //
      cluster_site_data_t iCSite_data =
          (*P_iQuasi_cluster_data)[iBucket][iCSite];

      //
      //  call function to process current cluster site of iQuasi
      //
      if (coreFlag.first == -1)
        PairwiseHarmonicShellQuasiClusterSiteProcessingLocal(
            iQuasi, coreFlag.second, iCSite_data, iBucket, iCSite,
            P_clusterForceEnergy, P_traceOfK);
      else
        PairwiseHarmonicCoreQuasiClusterSiteProcessingLocal(
            iQuasi, coreFlag.second, iCSite_data, iBucket, iCSite,
            P_clusterForceEnergy, P_traceOfK);
    } // loop over cluster site
  }   // loop over bucket

  return;
} // end of PairwiseQuasiHarmonicClusterBucketProcessingLocal()

//
//  iQuasiPairwiseForceEnergyHarmonicWorker()
//
void *iQuasiPairwiseForceEnergyHarmonicWorker(void *arg) {
  int iQuasi = ((struct iQuasiPairwiseForceEnergyHarmonicData_t *)arg)->iQuasi;
  int numQuasi =
      ((struct iQuasiPairwiseForceEnergyHarmonicData_t *)arg)->numQuasi;
  std::pair<int, int> coreFlag =
      ((struct iQuasiPairwiseForceEnergyHarmonicData_t *)arg)->coreFlag;
  std::vector<std::vector<cluster_site_data_t>> *P_iQuasi_cluster_data =
      ((struct iQuasiPairwiseForceEnergyHarmonicData_t *)arg)
          ->P_iQuasi_cluster_data;
  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      *P_clusterForceEnergy =
          ((struct iQuasiPairwiseForceEnergyHarmonicData_t *)arg)
              ->P_clusterForceEnergy;
  std::vector<std::vector<double>> *P_traceOfK =
      ((struct iQuasiPairwiseForceEnergyHarmonicData_t *)arg)->P_traceOfK;

  int number_threads = get_number_threads();
  int my_id = get_my_tid();
  int number_thread;
  int bucket_start;
  int bucket_end;

  get_share(my_id, number_threads, P_iQuasi_cluster_data->size(),
            &number_thread, &bucket_start, &bucket_end);

  PairwiseHarmonicClusterBucketProcessingLocal(
      iQuasi, numQuasi, coreFlag, P_iQuasi_cluster_data, bucket_start,
      bucket_end, P_clusterForceEnergy, P_traceOfK);

  return ((void *)NULL);
} // end of iQuasiPairwiseForceEnergyHarmonicWorker()
} // end of namespace for Pairwise Quasi Harmonic

//
//  namespace for Entropy
//
namespace {
struct iQuasiEntropyForceEnergyData_t {
  int iQuasi;
  int numQuasi;
  std::vector<std::vector<cluster_site_data_t>> *P_iQuasi_cluster_data;
  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      *P_clusterForceEnergy;
};

struct iQuasiEntropyForceEnergyHarmonicData_t {
  int iQuasi;
  int numQuasi;
  std::vector<std::vector<cluster_site_data_t>> *P_iQuasi_cluster_data;
  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      *P_clusterForceEnergy;
  std::vector<std::vector<double>> *P_traceOfK;
};

void *iQuasiEntropyForceEnergyWorker(void *arg) {
  int iQuasi = ((struct iQuasiEntropyForceEnergyData_t *)arg)->iQuasi;
  int numQuasi = ((struct iQuasiEntropyForceEnergyData_t *)arg)->numQuasi;
  std::vector<std::vector<cluster_site_data_t>> *P_iQuasi_cluster_data =
      ((struct iQuasiEntropyForceEnergyData_t *)arg)->P_iQuasi_cluster_data;
  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      *P_clusterForceEnergy =
          ((struct iQuasiEntropyForceEnergyData_t *)arg)->P_clusterForceEnergy;

  int number_threads = get_number_threads();
  int my_id = get_my_tid();
  int number_thread;
  int bucket_start;
  int bucket_end;

  get_share(my_id, number_threads, P_iQuasi_cluster_data->size(),
            &number_thread, &bucket_start, &bucket_end);

  // (N_A, k_B, h, h/2pi, epsilon, epsioln*4pi)
  std::vector<double> universalConstants =
      PairPotentials::getInstance()->getUniversalConstants();

  double boltzmanC = universalConstants[1];
  double maxPlanck_pi = universalConstants[3];

  double temperature = Quasicontinua::getInstance()->getTemperature();

  std::vector<double> sigmaVector =
      Quasicontinua::getInstance()->getSigmaVector();

  int iBucket;
  int iCluster;
  for (iBucket = bucket_start; iBucket <= bucket_end; iBucket++) {
    for (iCluster = 0; iCluster < (*P_iQuasi_cluster_data)[iBucket].size();
         iCluster++) {
      // compute entropy of a system
      double freq = (*P_iQuasi_cluster_data)[iBucket][iCluster].first.second[3];

      // if(iBucket == 0)
      // 	d_print("log(2.5) = %f\n", log(2.5));

      double entropy =
          3.0 * boltzmanC * log(sigmaVector[iQuasi] * sigmaVector[iQuasi] /
                                (maxPlanck_pi * freq)) +
          4 * boltzmanC;
      entropy = entropy * energy_factor_entropy;

      double d_entropy = -3.0 * boltzmanC / freq;
      d_entropy = d_entropy * energy_factor_entropy;

      double kinetic_energy =
          3.0 * boltzmanC * temperature * energy_factor_entropy / 2.0;

      //
      //  update the cluster site force energy data
      //
      //	Note that if Free energy F = V + W
      //	V : interatomic potential
      //  W : entropic energy = -T * S
      //
      //	We define force fx = -partial F/ partial x
      //  Similarly for y,z, frequency
      //
      //	Hence, f_freq = - partial W/partial freq
      //	= - ( - T * partial S/partial freq)
      //  = T * partial S/partial freq
      //
      pthread_mutex_lock(&bucketLock[iQuasi][iBucket]);
      (*P_clusterForceEnergy)[iBucket][iCluster].first[3] +=
          temperature * d_entropy;

      (*P_clusterForceEnergy)[iBucket][iCluster].second +=
          -entropy * temperature;
      pthread_mutex_unlock(&bucketLock[iQuasi][iBucket]);
    }
  }

  return ((void *)NULL);
}

void *iQuasiEntropyForceEnergyHarmonicWorker(void *arg) {
  int iQuasi = ((struct iQuasiEntropyForceEnergyHarmonicData_t *)arg)->iQuasi;
  int numQuasi =
      ((struct iQuasiEntropyForceEnergyHarmonicData_t *)arg)->numQuasi;
  std::vector<std::vector<cluster_site_data_t>> *P_iQuasi_cluster_data =
      ((struct iQuasiEntropyForceEnergyHarmonicData_t *)arg)
          ->P_iQuasi_cluster_data;
  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      *P_clusterForceEnergy =
          ((struct iQuasiEntropyForceEnergyHarmonicData_t *)arg)
              ->P_clusterForceEnergy;
  std::vector<std::vector<double>> *P_traceOfK =
      ((struct iQuasiEntropyForceEnergyHarmonicData_t *)arg)->P_traceOfK;

  int number_threads = get_number_threads();
  int my_id = get_my_tid();
  int number_thread;
  int bucket_start;
  int bucket_end;

  get_share(my_id, number_threads, P_iQuasi_cluster_data->size(),
            &number_thread, &bucket_start, &bucket_end);

  // (N_A, k_B, h, h/2pi, epsilon, epsioln*4pi)
  std::vector<double> universalConstants =
      PairPotentials::getInstance()->getUniversalConstants();

  double boltzmanC = universalConstants[1];
  double maxPlanck_pi = universalConstants[3];

  double temperature = Quasicontinua::getInstance()->getTemperature();

  std::vector<double> sigmaVector =
      Quasicontinua::getInstance()->getSigmaVector();

  double mass = Quasicontinua::getInstance()->getQuasi(iQuasi).getAtomicMass();

  int iBucket;
  int iCluster;
  for (iBucket = bucket_start; iBucket <= bucket_end; iBucket++) {
    for (iCluster = 0; iCluster < (*P_iQuasi_cluster_data)[iBucket].size();
         iCluster++) {
      // compute entropy of a system
      double freq = (*P_iQuasi_cluster_data)[iBucket][iCluster].first.second[3];

      //	entropy for quasi harmonic approximation
      double traceK = (*P_traceOfK)[iBucket][iCluster];

      double a =
          6.0 * boltzmanC * temperature * energy_factor_entropy / maxPlanck_pi;
      double b = traceK / freq + 3.0 * freq * energy_factor_kinetic / mass;

      double entropy = 3.0 * boltzmanC * (log(a) - log(b)) + 4.0 * boltzmanC;

      entropy = entropy * energy_factor_entropy;

      //	for quasi harmonic we also need to add 3 K_b * T to energy
      double add_energy = 3.0 * boltzmanC * temperature * energy_factor_entropy;

      double d_entropy =
          -3.0 * boltzmanC *
          (-traceK / (freq * freq) + 3.0 * energy_factor_kinetic / mass) / (b);

      // scaling to keep entropy in correct unit.
      d_entropy = d_entropy * energy_factor_entropy;

      //
      //  update the cluster site force energy data
      //
      //
      //	Note that if Free energy F = V + W
      //	V : interatomic potential
      //  W : entropic energy = -T * S
      //
      //	We define force fx = -partial F/ partial x
      //  Similarly for y,z, frequency
      //
      //	Hence, f_freq = - partial W/partial freq
      //	= - ( - T * partial S/partial freq)
      //  = T * partial S/partial freq
      //
      pthread_mutex_lock(&bucketLock[iQuasi][iBucket]);
      (*P_clusterForceEnergy)[iBucket][iCluster].first[3] +=
          temperature * d_entropy;

      (*P_clusterForceEnergy)[iBucket][iCluster].second +=
          -entropy * temperature + add_energy;
      pthread_mutex_unlock(&bucketLock[iQuasi][iBucket]);
    }
  }

  return ((void *)NULL);
}
} // end of namespace for iQuasiEntropyForceEnergyWorker

//
//  namespace for Electrostatics
//
namespace {
struct iQuasiElectrostaticsForceEnergyData_t {
  int iQuasi;
  int numQuasi;
  bool compute_in_reference;
  std::vector<std::vector<cluster_site_data_t>> *P_iQuasi_cluster_data;
  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      *P_clusterForceEnergy;
};

//
//  ElectrostaticsClusterBucketProcessingLocal()
//
void ElectrostaticsClusterBucketProcessingLocal(
    int iQuasi, int numQuasi, bool compute_in_reference,
    std::vector<std::vector<cluster_site_data_t>> *P_iQuasi_cluster_data,
    int bucket_start, int bucket_end,
    std::vector<std::vector<std::pair<std::vector<double>, double>>>
        *P_clusterForceEnergy) {
  //
  //  get Quasicontinua, PairPotentials
  //
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  PairPotentials *pairC = PairPotentials::getInstance();

  //
  //  loop over bucket
  //
  int iBucket;
  int iCSite;
  for (iBucket = bucket_start; iBucket <= bucket_end; iBucket++) {
    //
    //	check if we are trying to access the memory not present
    //	in P_clusterForceEnergy
    //
    if ((*P_clusterForceEnergy)[iBucket].size() <
        (*P_iQuasi_cluster_data)[iBucket].size()) {
      pthread_mutex_lock(&(bucketLock[iQuasi][iBucket]));
      d_print("Check InitializeClusterForceEnergyData(), size of cluster force "
              "and clustyer site data must match\n");
      pthread_mutex_unlock(&(bucketLock[iQuasi][iBucket]));
    }

    for (iCSite = 0; iCSite < (*P_iQuasi_cluster_data)[iBucket].size();
         iCSite++) {
      //
      //  process cluster site
      //
      cluster_site_data_t iCSite_data =
          (*P_iQuasi_cluster_data)[iBucket][iCSite];

      //
      //
      //
      std::vector<double> force;
      for (int dof = 0; dof < 4; dof++)
        force.push_back(0.0);

      double energy;
      energy = 0.0;

      Electrostatics::getInstance()->computeForceOnSite(
          iQuasi, iCSite_data, force, energy, compute_in_reference);

      // //
      // //  call Electrostatics function to process cluster site
      // //
      // Electrostatics::getInstance()->computeForceOnNodeDueToSite(iQuasi,
      // 	iCSite_data,
      // 	compute_in_reference);

      // pthread_mutex_lock(&lock_print);
      // d_print("ElectrostaticsClusterBucketProcessingLocal\n");
      // pthread_mutex_unlock(&lock_print);

      // //
      // //  update the cluster site force energy data
      // //
      // if(forceEnergy.first.size() == 0)
      // 	std::cout<<"error"<<std::endl;

      pthread_mutex_lock(&(bucketLock[iQuasi][iBucket]));
      (*P_clusterForceEnergy)[iBucket][iCSite].first[0] += force[0];
      (*P_clusterForceEnergy)[iBucket][iCSite].first[1] += force[1];
      (*P_clusterForceEnergy)[iBucket][iCSite].first[2] += force[2];
      (*P_clusterForceEnergy)[iBucket][iCSite].first[3] += force[3];

      (*P_clusterForceEnergy)[iBucket][iCSite].second += energy;
      pthread_mutex_unlock(&(bucketLock[iQuasi][iBucket]));

    } // loop over cluster site
  }   // loop over bucket

  return;
} // end of ElectrostaticsClusterBucketProcessingLocal()

//
//  iQuasiElectrostaticsForceEnergyWorker()
//
void *iQuasiElectrostaticsForceEnergyWorker(void *arg) {
  int iQuasi = ((struct iQuasiElectrostaticsForceEnergyData_t *)arg)->iQuasi;
  int numQuasi =
      ((struct iQuasiElectrostaticsForceEnergyData_t *)arg)->numQuasi;
  bool compute_in_reference =
      ((struct iQuasiElectrostaticsForceEnergyData_t *)arg)
          ->compute_in_reference;
  std::vector<std::vector<cluster_site_data_t>> *P_iQuasi_cluster_data =
      ((struct iQuasiElectrostaticsForceEnergyData_t *)arg)
          ->P_iQuasi_cluster_data;
  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      *P_clusterForceEnergy =
          ((struct iQuasiElectrostaticsForceEnergyData_t *)arg)
              ->P_clusterForceEnergy;

  int number_threads = get_number_threads();
  int my_id = get_my_tid();
  int number_thread;
  int bucket_start;
  int bucket_end;

  get_share(my_id, number_threads, P_iQuasi_cluster_data->size(),
            &number_thread, &bucket_start, &bucket_end);

  // pthread_mutex_lock(&lock_print);
  // d_print("iQuasiElectrostaticsForceEnergyWorker\n");
  // pthread_mutex_unlock(&lock_print);

  ElectrostaticsClusterBucketProcessingLocal(
      iQuasi, numQuasi, compute_in_reference, P_iQuasi_cluster_data,
      bucket_start, bucket_end, P_clusterForceEnergy);

  return ((void *)NULL);
} // end of iQuasiElectrostaticsForceEnergyWorker()

} // end of namespace for Electrostatics

//
//  namespace for additional works
//
namespace {
// pthread lock used in total energy worker
static pthread_mutex_t totalEnergyLock = PTHREAD_MUTEX_INITIALIZER;

//
//  iQuasiAddExternalForcesWorker()
//
void *iQuasiAddExternalForcesWorker(void *arg) {
  struct node_list_t *P_node_list = (struct node_list_t *)arg;

  int i_node;

  int my_id = get_my_tid();
  int number_threads = get_number_threads();
  int node_start;
  int node_end;
  int node_thread;

  /**
   * compute share of load
   */

  get_share(my_id, number_threads, P_node_list->number_nodes, &node_thread,
            &node_start, &node_end);

  /**
   * loop over nodes
   */

  for (i_node = node_start; i_node <= node_end; i_node++) {

    struct node_t *P_node = P_node_list->nodes[i_node];

    // if(residual_debug == 1)
    // {
    // 	if(P_node->l[0] >= 10 && P_node->l[0] <= 15)
    // 		if(P_node->l[1] >= 10 && P_node->l[1] <= 15)
    // 			if(P_node->l[2] >= 10 && P_node->l[2] <= 15)
    // 			{
    // 				pthread_mutex_lock(&lock_print);
    // 				d_print("force = (%f, %f, %f, %f)\n",
    // 					P_node->acceleration[0],
    // 					P_node->acceleration[1],
    // 					P_node->acceleration[2],
    // 					P_node->force_frequency);
    // 				pthread_mutex_unlock(&lock_print);
    // 			}
    // }

    P_node->acceleration[0] += P_node->external_force[0];
    P_node->acceleration[1] += P_node->external_force[1];
    P_node->acceleration[2] += P_node->external_force[2];

    P_node->force_frequency += P_node->external_freq_force;

    // if(residual_debug == 1)
    // {
    // 	if(P_node->l[0] >= 10 && P_node->l[0] <= 15)
    // 		if(P_node->l[1] >= 10 && P_node->l[1] <= 15)
    // 			if(P_node->l[2] >= 10 && P_node->l[2] <= 15)
    // 			{
    // 				pthread_mutex_lock(&lock_print);
    // 				d_print("ext force = (%f, %f, %f, %f)\n",
    // 					P_node->external_force[0],
    // 					P_node->external_force[1],
    // 					P_node->external_force[2],
    // 					P_node->external_freq_force);
    // 				d_print("force = (%f, %f, %f, %f)\n",
    // 					P_node->acceleration[0],
    // 					P_node->acceleration[1],
    // 					P_node->acceleration[2],
    // 					P_node->force_frequency);
    // 				pthread_mutex_unlock(&lock_print);
    // 			}
    // }
  }

  return ((void *)NULL);
} // end of iQuasiAddExternalForcesWorker()

//
//  iQuasiComputeTotalEnergyWorker()
//
void *iQuasiComputeTotalEnergyWorker(void *arg) {
  struct all_node_list_t *P_all_node_list = (struct all_node_list_t *)arg;

  int number_threads = get_number_threads();
  int my_id = get_my_tid();
  int node_per_thread;
  int node_start;
  int node_end;

  /**
   * get share of clusters
   */

  get_share(my_id, number_threads, P_all_node_list->node_list.number_nodes,
            &node_per_thread, &node_start, &node_end);

  // process data
  struct node_list_t *P_node_list = &P_all_node_list->node_list;
  struct energy_t local_energy = ENERGY_INITIALIZER;

  int i_node;

  /**
   * loop over all clusters
   */

  for (i_node = node_start; i_node <= node_end; i_node++) {
    struct node_t *P_node = P_node_list->nodes[i_node];

    local_energy.kinetic += P_node->energy.kinetic;
    local_energy.potential += P_node->energy.potential;
  }

  /**
   * add local contribution
   */

  pthread_mutex_lock(&totalEnergyLock);

  P_all_node_list->energy.kinetic += local_energy.kinetic;
  P_all_node_list->energy.potential += local_energy.potential;

  pthread_mutex_unlock(&totalEnergyLock);

  return ((void *)NULL);
} // end of iQuasiComputeTotalEnergyWorker()

//
//  iQuasiComputeExternalEnergyWorker()
//
void *iQuasiComputeExternalEnergyWorker(void *arg) {
  struct all_node_list_t *P_all_node_list = (struct all_node_list_t *)arg;

  struct node_list_t *P_node_list = &(P_all_node_list->node_list);

  int i_node;
  int my_id = get_my_tid();
  int number_threads = get_number_threads();
  int node_start;
  int node_end;
  int node_thread;

  /**
   * compute share of load
   */

  get_share(my_id, number_threads, P_node_list->number_nodes, &node_thread,
            &node_start, &node_end);

  /**
   * loop over nodes
   */

  for (i_node = node_start; i_node <= node_end; i_node++) {
    struct node_t *P_node = P_node_list->nodes[i_node];

    pthread_mutex_lock(&totalEnergyLock);

    P_all_node_list->energy.potential -=
        P_node->external_force[0] *
            (P_node->position[0] - P_node->initial_position[0]) +
        P_node->external_force[1] *
            (P_node->position[1] - P_node->initial_position[1]) +
        P_node->external_force[2] *
            (P_node->position[2] - P_node->initial_position[2]) +
        P_node->external_freq_force *
            (P_node->frequency - P_node->initial_frequency);

    pthread_mutex_unlock(&totalEnergyLock);
  }

  return ((void *)NULL);
} // end of iQuasiComputeExternalEnergyWorker()

struct testAddForceEnergyData_t {
  int iQuasi;
  int flag;
  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      *P_clusterForceData;
  std::vector<std::vector<cluster_site_data_t>> *P_clusterSiteData;
  std::vector<std::vector<std::pair<std::vector<double>, std::vector<double>>>>
      *P_testClusterForceData;
  std::vector<std::pair<std::vector<double>, std::vector<double>>>
      *P_testNodeForceData;
  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      *P_testTotalClusterForceData;
  std::vector<std::pair<std::vector<double>, double>> *P_testTotalNodeForceData;
};

void *testAddForceEnergyWorker(void *arg) {
  int iQuasi = ((struct testAddForceEnergyData_t *)arg)->iQuasi;
  const int flag = ((struct testAddForceEnergyData_t *)arg)->flag;
  std::vector<std::vector<cluster_site_data_t>> *P_clusterSiteData =
      ((struct testAddForceEnergyData_t *)arg)->P_clusterSiteData;

  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      *P_clusterForceData =
          ((struct testAddForceEnergyData_t *)arg)->P_clusterForceData;

  std::vector<std::vector<std::pair<std::vector<double>, std::vector<double>>>>
      *P_testClusterForceData =
          ((struct testAddForceEnergyData_t *)arg)->P_testClusterForceData;

  std::vector<std::pair<std::vector<double>, std::vector<double>>>
      *P_testNodeForceData =
          ((struct testAddForceEnergyData_t *)arg)->P_testNodeForceData;

  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      *P_testTotalClusterForceData =
          ((struct testAddForceEnergyData_t *)arg)->P_testTotalClusterForceData;

  std::vector<std::pair<std::vector<double>, double>>
      *P_testTotalNodeForceData =
          ((struct testAddForceEnergyData_t *)arg)->P_testTotalNodeForceData;

  int number_threads = get_number_threads();
  int my_id = get_my_tid();
  int number_thread;
  int bucket_start;
  int bucket_end;

  get_share(my_id, number_threads, P_clusterSiteData->size(), &number_thread,
            &bucket_start, &bucket_end);

  struct node_list_t *P_node_list =
      &(Quasicontinua::getInstance()->getQuasi(iQuasi).getNodeList().node_list);

  struct lattice_t *P_lattice =
      &(Quasicontinua::getInstance()->getQuasi(iQuasi).getLattice());

  double boltzmanC = PairPotentials::getInstance()->getBoltzmanConstant();

  double temperature = Quasicontinua::getInstance()->getTemperature();

  // get Shape Instance
  Shape *shapeC = Shape::getInstance();

  for (int iB = bucket_start; iB <= bucket_end; iB++) {
    for (int iC = 0; iC < (*P_clusterSiteData)[iB].size(); iC++) {
      cluster_site_data_t iC_data = (*P_clusterSiteData)[iB][iC];

      //
      //	add data of P_clusterForceData to P_testClusterForceData
      //
      int forceIndex[4];
      int freqIndex = -1;
      int energyIndex = -1;
      if (flag == 0) {
        for (int i = 0; i < 4; i++)
          forceIndex[i] = i;

        energyIndex = 0;
      }

      if (flag == 1) {
        for (int i = 0; i < 4; i++)
          forceIndex[i] = i + 4;

        energyIndex = 1;
      }

      if (flag == 2) {
        freqIndex = 8;
        energyIndex = 2;
      }

      if (flag == 3) {
        for (int i = 0; i < 4; i++)
          forceIndex[i] = i + 9;
        energyIndex = 3;
      }

      // add
      // separate the entropy job with others
      if (flag == 2) {
        (*P_testClusterForceData)[iB][iC].first[freqIndex] =
            (*P_clusterForceData)[iB][iC].first[3];

        (*P_testClusterForceData)[iB][iC].second[energyIndex] =
            (*P_clusterForceData)[iB][iC].second;

        (*P_testTotalClusterForceData)[iB][iC].first[3] +=
            (*P_clusterForceData)[iB][iC].first[3];

        (*P_testTotalClusterForceData)[iB][iC].second +=
            (*P_clusterForceData)[iB][iC].second;
      } else {
        (*P_testClusterForceData)[iB][iC].first[forceIndex[0]] =
            (*P_clusterForceData)[iB][iC].first[0];
        (*P_testClusterForceData)[iB][iC].first[forceIndex[1]] =
            (*P_clusterForceData)[iB][iC].first[1];
        (*P_testClusterForceData)[iB][iC].first[forceIndex[2]] =
            (*P_clusterForceData)[iB][iC].first[2];
        (*P_testClusterForceData)[iB][iC].first[forceIndex[3]] =
            (*P_clusterForceData)[iB][iC].first[3];

        (*P_testClusterForceData)[iB][iC].second[energyIndex] =
            (*P_clusterForceData)[iB][iC].second;

        // total cluster data
        (*P_testTotalClusterForceData)[iB][iC].first[0] +=
            (*P_clusterForceData)[iB][iC].first[0];
        (*P_testTotalClusterForceData)[iB][iC].first[1] +=
            (*P_clusterForceData)[iB][iC].first[1];
        (*P_testTotalClusterForceData)[iB][iC].first[2] +=
            (*P_clusterForceData)[iB][iC].first[2];
        (*P_testTotalClusterForceData)[iB][iC].first[3] +=
            (*P_clusterForceData)[iB][iC].first[3];

        (*P_testTotalClusterForceData)[iB][iC].second +=
            (*P_clusterForceData)[iB][iC].second;
      }

      int node_err = 0;
      if (iC_data.second.size() == 5)
        node_err = -1;

      int node_0 = iC_data.second[0];
      int node_1 = iC_data.second[1];
      int node_2 = iC_data.second[2];
      int node_3 = iC_data.second[3];

      struct node_t *P_node_0 = P_node_list->nodes[iC_data.second[0]];
      struct node_t *P_node_1 = P_node_list->nodes[iC_data.second[1]];
      struct node_t *P_node_2 = P_node_list->nodes[iC_data.second[2]];
      struct node_t *P_node_3 = P_node_list->nodes[iC_data.second[3]];

      double nodeWeight;
      if (node_err == -1)
        nodeWeight = P_node_list->nodes[iC_data.second[4]]->weight;
      else
        nodeWeight = P_node_0->weight;

      //
      //  get value of shape function of node_0 at jC
      //
      int l_iC[3];
      for (int dof = 0; dof < 3; dof++)
        l_iC[dof] = iC_data.first.first[dof];

      // process node 0
      double shape = shapeC->computeSiteShapeFunction(
          P_node_0, P_node_1, P_node_2, P_node_3, l_iC, P_lattice, iQuasi);

      // update potential and kinetic energy
      if (node_err == -1) {
        int node_5 = iC_data.second[4];

        pthread_mutex_lock(&lock_print);
        (*P_testNodeForceData)[node_5].second[energyIndex] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].second;
        pthread_mutex_unlock(&lock_print);

        pthread_mutex_lock(&lock_print);
        (*P_testTotalNodeForceData)[node_5].second +=
            nodeWeight * (*P_clusterForceData)[iB][iC].second;
        pthread_mutex_unlock(&lock_print);
      }

      // add contribution to node_0
      pthread_mutex_lock(&lock_print);
      if (flag == 2) {
        (*P_testNodeForceData)[node_0].first[freqIndex] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;

        (*P_testTotalNodeForceData)[node_0].first[3] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;
      } else {
        (*P_testNodeForceData)[node_0].first[forceIndex[0]] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[0] * shape;
        (*P_testNodeForceData)[node_0].first[forceIndex[1]] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[1] * shape;
        (*P_testNodeForceData)[node_0].first[forceIndex[2]] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[2] * shape;
        (*P_testNodeForceData)[node_0].first[forceIndex[3]] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;

        (*P_testTotalNodeForceData)[node_0].first[0] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[0] * shape;
        (*P_testTotalNodeForceData)[node_0].first[1] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[1] * shape;
        (*P_testTotalNodeForceData)[node_0].first[2] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[2] * shape;
        (*P_testTotalNodeForceData)[node_0].first[3] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;
      }

      if (node_err == 0) {
        (*P_testNodeForceData)[node_0].second[energyIndex] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].second;

        (*P_testTotalNodeForceData)[node_0].second +=
            nodeWeight * (*P_clusterForceData)[iB][iC].second;
      }
      pthread_mutex_unlock(&lock_print);

      // process node 1
      shape = shapeC->computeSiteShapeFunction(
          P_node_1, P_node_0, P_node_2, P_node_3, l_iC, P_lattice, iQuasi);

      // add contribution to node_1
      pthread_mutex_lock(&lock_print);
      if (flag == 2) {
        (*P_testNodeForceData)[node_1].first[freqIndex] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;

        (*P_testTotalNodeForceData)[node_1].first[3] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;
      } else {
        (*P_testNodeForceData)[node_1].first[forceIndex[0]] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[0] * shape;
        (*P_testNodeForceData)[node_1].first[forceIndex[1]] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[1] * shape;
        (*P_testNodeForceData)[node_1].first[forceIndex[2]] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[2] * shape;
        (*P_testNodeForceData)[node_1].first[forceIndex[3]] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;

        (*P_testTotalNodeForceData)[node_1].first[0] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[0] * shape;
        (*P_testTotalNodeForceData)[node_1].first[1] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[1] * shape;
        (*P_testTotalNodeForceData)[node_1].first[2] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[2] * shape;
        (*P_testTotalNodeForceData)[node_1].first[3] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;
      }
      pthread_mutex_unlock(&lock_print);

      // process node 2
      shape = shapeC->computeSiteShapeFunction(
          P_node_2, P_node_0, P_node_1, P_node_3, l_iC, P_lattice, iQuasi);

      // add contribution to node_2
      pthread_mutex_lock(&lock_print);
      if (flag == 2) {
        (*P_testNodeForceData)[node_2].first[freqIndex] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;

        (*P_testTotalNodeForceData)[node_2].first[3] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;
      } else {
        (*P_testNodeForceData)[node_2].first[forceIndex[0]] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[0] * shape;
        (*P_testNodeForceData)[node_2].first[forceIndex[1]] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[1] * shape;
        (*P_testNodeForceData)[node_2].first[forceIndex[2]] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[2] * shape;
        (*P_testNodeForceData)[node_2].first[forceIndex[3]] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;

        (*P_testTotalNodeForceData)[node_2].first[0] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[0] * shape;
        (*P_testTotalNodeForceData)[node_2].first[1] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[1] * shape;
        (*P_testTotalNodeForceData)[node_2].first[2] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[2] * shape;
        (*P_testTotalNodeForceData)[node_2].first[3] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;
      }
      pthread_mutex_unlock(&lock_print);

      // process node 3
      shape = shapeC->computeSiteShapeFunction(
          P_node_3, P_node_0, P_node_1, P_node_2, l_iC, P_lattice, iQuasi);

      // add contribution to node_3
      pthread_mutex_lock(&lock_print);
      if (flag == 2) {
        (*P_testNodeForceData)[node_3].first[freqIndex] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;

        (*P_testTotalNodeForceData)[node_3].first[3] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;
      } else {
        (*P_testNodeForceData)[node_3].first[forceIndex[0]] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[0] * shape;
        (*P_testNodeForceData)[node_3].first[forceIndex[1]] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[1] * shape;
        (*P_testNodeForceData)[node_3].first[forceIndex[2]] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[2] * shape;
        (*P_testNodeForceData)[node_3].first[forceIndex[3]] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;

        (*P_testTotalNodeForceData)[node_3].first[0] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[0] * shape;
        (*P_testTotalNodeForceData)[node_3].first[1] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[1] * shape;
        (*P_testTotalNodeForceData)[node_3].first[2] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[2] * shape;
        (*P_testTotalNodeForceData)[node_3].first[3] +=
            nodeWeight * (*P_clusterForceData)[iB][iC].first[3] * shape;
      }
      pthread_mutex_unlock(&lock_print);
    }
  }

  return NULL;
}
} // end of namespace for additional works

ForceEnergyCalculation *ForceEnergyCalculation::_instance = NULL;

//
// constructor
//

ForceEnergyCalculation::ForceEnergyCalculation() {
  // by default : there are no atomistic loads.
  // call noAtomisticLoads() to initialize the data
  noAtomisticLoads();

  //
  //
  //
  return;
}

//
// destructor
//

ForceEnergyCalculation::~ForceEnergyCalculation() {

  //
  //
  //
  return;
}

//
// getInstance method
//

ForceEnergyCalculation *ForceEnergyCalculation::getInstance() {

  //
  // if not created, create
  //
  if (_instance == NULL) {
    _instance = new ForceEnergyCalculation();
  }

  //
  // return instance
  //
  return _instance;
}

//
// destroy instance method
//

void ForceEnergyCalculation::destroyInstance() {

  //
  // delete instance
  //
  delete _instance;

  //
  //
  //
  return;
}

//
//  noAtomisticLoads()
//
// 	this function is called by constructor of
// 	ForceEnergyCalculation
//
void ForceEnergyCalculation::noAtomisticLoads(void) {
  //
  //	we have no atomistic loads for this system
  //

  // set d_atomisticFlag to 0
  d_atomisticFlag = 0;

  d_atomisticElementSize = 0.0;

  d_atomisticDeadLoads.clear();

  return;
}

//
//  setAtomisticDeadLoads()
//
void ForceEnergyCalculation::setAtomisticDeadLoads(
    std::vector<std::vector<double>> forces, const double elementSize) {
  //
  //	set the atomistic flag to 1
  //
  d_atomisticFlag = 1;

  //
  // resize the d_atomisticDeadLoads
  //
  d_atomisticDeadLoads.clear();
  d_atomisticDeadLoads.resize(forces.size());

  // set each dof
  int iQuasi;
  int dof;
  for (iQuasi = 0; iQuasi < forces.size(); ++iQuasi)
    for (dof = 0; dof < 3; ++dof)
      d_atomisticDeadLoads[iQuasi].push_back(forces[iQuasi][dof]);

  // set atomistic element size
  d_atomisticElementSize = elementSize;

  return;
} // end of setAtomisticDeadLoads()

//
//  setRemoveResidualFlags()
//
void ForceEnergyCalculation::setRemoveResidualFlags(
    const std::vector<int> flags) {
  //
  //  resize th d_residualForcesFlags()
  //
  if (flags.size() != 0) {
    d_residualForcesFlags.resize(flags.size());
    for (int i = 0; i < flags.size(); ++i) {
      if (flags[i] == 0)
        d_residualForcesFlags[i] = false;
      else if (flags[i] == 1)
        d_residualForcesFlags[i] = true;

      else {
        d_print("Forces flag not properly set.\n");
        D_ERROR("setRemoveResidualFlags");
        exit(0);
      }
    }
  } else {
    d_print("Forces flags supplied are empty");
    D_ERROR("setRemoveResidualFlags()");
    exit(0);
  }

  return;
} // end of setRemoveResidualFlags()

//
//  removeResidualForces()
//
void ForceEnergyCalculation::removeResidualForces(void) {
  residual_debug = 0;

  //
  //	time data
  //
  timeval t1, t2;
  double elapsedTime;
  gettimeofday(&t1, NULL);

  //
  //	get Quasicontinua and it's size
  //
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();
  int numQuasi = quasicontinua->size();

  //
  //  vector of int to store flags for each node
  //
  std::vector<std::vector<bool>> nodeRemoveFlags;
  nodeRemoveFlags.resize(numQuasi);

  //
  //  reset to node state to it's initial state
  //
  std::vector<std::vector<std::vector<double>>> currentStates;

  for (int iQuasi = 0; iQuasi < numQuasi; iQuasi++) {
    std::vector<std::vector<double>> iQuasiStates =
        quasicontinua->getQuasi(iQuasi).resetToInitialState();

    currentStates.push_back(iQuasiStates);

    // reset the external forces to zero
    quasicontinua->getQuasi(iQuasi).initializeExternalForces();
  }

  //
  //  compute force
  //  Note : We need to rebuild the neighbor list as we have changed the
  //  position of nodes to initial configuration.
  //
  int rebuild_neighbor_flag = 1;
  bool compute_in_reference = true;
  computeForceEnergy(rebuild_neighbor_flag, compute_in_reference);
  residual_debug = 1;

  //
  // nodal forces for surface nodes will (in general) be non-zero.
  // copy them as -external forces for nodes that have weight > W
  // (i.e. nodes that are NOT in the atomistic region) with W > 1
  // default is to not remove for all atoms
  //
  for (int iQuasi = 0; iQuasi < numQuasi; iQuasi++) {
    Quasicontinuum iQuasicontinuum = quasicontinua->getQuasi(iQuasi);

    struct all_node_list_t *P_all_node_list = &(iQuasicontinuum.getNodeList());

    struct lattice_t lattice = iQuasicontinuum.getLattice();

    //
    //  resize the flag vector
    //
    nodeRemoveFlags[iQuasi].resize(P_all_node_list->node_list.number_nodes);

    // get Node instance
    Node *nodeC = Node::getInstance();

    for (int iNode = 0; iNode < P_all_node_list->node_list.number_nodes;
         iNode++) {
      nodeRemoveFlags[iQuasi][iNode] = false;

      struct node_t *P_node = P_all_node_list->node_list.nodes[iNode];

      //
      // skip atomistic nodes ( by default d_atomisticFlag=0)
      // therefore, by default it will not skip
      // any node.
      // See constructor of ForceEnergyCalculation class.
      //
      if (d_atomisticFlag == 1)
        if (nodeC->isNodeAtomistic(P_node, d_atomisticElementSize))
          continue;

      //
      // check for residual flags
      //
      if (d_residualForcesFlags[0] == true) {
        nodeRemoveFlags[iQuasi][iNode] = true;
        continue;
      }

      //
      // check for xBeg
      //
      if (d_residualForcesFlags[1] == true &&
          P_node->l[0] == lattice.l_start[0]) {
        nodeRemoveFlags[iQuasi][iNode] = true;
        SetResidualFlagsForAdjacentNodes(P_node, nodeRemoveFlags[iQuasi]);
        continue;
      }

      //
      // check for xEnd
      //
      if (d_residualForcesFlags[2] == true &&
          P_node->l[0] == lattice.l_end[0]) {
        nodeRemoveFlags[iQuasi][iNode] = true;
        SetResidualFlagsForAdjacentNodes(P_node, nodeRemoveFlags[iQuasi]);
        continue;
      }

      //
      // check for yBeg
      //
      if (d_residualForcesFlags[3] == true &&
          P_node->l[1] == lattice.l_start[1]) {
        nodeRemoveFlags[iQuasi][iNode] = true;
        SetResidualFlagsForAdjacentNodes(P_node, nodeRemoveFlags[iQuasi]);
        continue;
      }

      //
      // check for yEnd
      //
      if (d_residualForcesFlags[4] == true &&
          P_node->l[2] == lattice.l_end[2]) {
        nodeRemoveFlags[iQuasi][iNode] = true;
        SetResidualFlagsForAdjacentNodes(P_node, nodeRemoveFlags[iQuasi]);
        continue;
      }

      //
      // check for zBeg
      //
      if (d_residualForcesFlags[5] == true &&
          P_node->l[3] == lattice.l_start[3]) {
        nodeRemoveFlags[iQuasi][iNode] = true;
        SetResidualFlagsForAdjacentNodes(P_node, nodeRemoveFlags[iQuasi]);
        continue;
      }

      //
      // check for zEnd
      //
      if (d_residualForcesFlags[6] == true &&
          P_node->l[3] == lattice.l_end[3]) {
        nodeRemoveFlags[iQuasi][iNode] = true;
        SetResidualFlagsForAdjacentNodes(P_node, nodeRemoveFlags[iQuasi]);
        continue;
      }
    } // loop over nodes

    //
    //  remove forces and reset the current configuration of quasi
    //
    iQuasicontinuum.resetToCurrentState(currentStates[iQuasi]);

    for (int iNode = 0; iNode < P_all_node_list->node_list.number_nodes;
         iNode++) {
      struct node_t *P_node = P_all_node_list->node_list.nodes[iNode];

      //
      //  add atomistic dead loads
      //	by default d_atomisticFlag = 0
      //
      //	unless setAtomisticDeadLoads() function is called
      //  the value of d_atomisticFlag will be 0
      //
      if (d_atomisticFlag == 1) {
        if (nodeC->isNodeAtomistic(P_node, d_atomisticElementSize)) {
          P_node->external_force[0] += d_atomisticDeadLoads[iQuasi][0];
          P_node->external_force[1] += d_atomisticDeadLoads[iQuasi][1];
          P_node->external_force[2] += d_atomisticDeadLoads[iQuasi][2];

          P_node->external_freq_force += 0.0;

          continue;
        }
      }

      //
      //  set external forces equal to computed force in computeForceEnergy()
      //
      if (nodeRemoveFlags[iQuasi][iNode] == false)
        continue;

      P_node->external_force[0] = -P_node->acceleration[0];
      P_node->external_force[1] = -P_node->acceleration[1];
      P_node->external_force[2] = -P_node->acceleration[2];

      // P_node->external_freq_force = -P_node->force_frequency;

      // if(P_node->l[0] >= 10 && P_node->l[0] <= 15)
      // 	if(P_node->l[1] >= 10 && P_node->l[1] <= 15)
      // 		if(P_node->l[2] >= 10 && P_node->l[2] <= 15)
      // 		{
      // 			d_print("force = (%f, %f, %f, %f)\n",
      // 				P_node->acceleration[0],
      // 				P_node->acceleration[1],
      // 				P_node->acceleration[2],
      // 				P_node->force_frequency);
      // 			d_print("ext force = (%f, %f, %f, %f)\n",
      // 				P_node->external_force[0],
      // 				P_node->external_force[1],
      // 				P_node->external_force[2],
      // 				P_node->external_freq_force);
      // 		}
    } // loop over nodes
  }   // loop over quasi

  // get time
  gettimeofday(&t2, NULL);
  elapsedTime = (t2.tv_sec - t1.tv_sec);
  elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000000.0;
  if (RunData::getInstance()->d_timeOut == true)
    RunData::getInstance()->addComputeResidualForceTime(0, elapsedTime);
  return;
} // end of removeResidualForces()

//
// private functions
//

//
//  SetResidualFlagsForAdjacentNodes()
//
void ForceEnergyCalculation::SetResidualFlagsForAdjacentNodes(
    const struct node_t *P_node, std::vector<bool> &node_residual_flags) {
  int jNode;
  int nodeAdd;
  for (jNode = 0; jNode < P_node->node_list.number_nodes; jNode++) {
    nodeAdd = P_node->node_list.nodes[jNode]->number;

    node_residual_flags[nodeAdd] = true;
  }

  return;
} // end of SetResidualFlagsForAdjacentNodes()

//
//  computeForceEnergy()
//
void ForceEnergyCalculation::computeForceEnergy(int rebuild_neighbor_flag,
                                                bool compute_in_reference) {
  int numQuasi = Quasicontinua::getInstance()->size();

  // time data
  timeval t1, t2;
  double elapsedTime;

  //
  //  first check if neighbor list is required to rebuild
  //
  int rebuild_cluster_flag = 0;

  CrossNeighborList *neighborC = CrossNeighborList::getInstance();
  neighborC->computeCrossNeighborList(rebuild_neighbor_flag,
                                      rebuild_cluster_flag);

  // initialize bucket
  int numBuckets = neighborC->getNumberBuckets();
  bucketLockInitLocal(numQuasi, numBuckets);

  //
  //	get the problem info
  //
  //	if minMethod = 2 and also quadMethod = 1
  //	then we are minimizing the frequency using quasi harmonic
  //	approximation. In this case force calculation from
  //	electrostatics is not needed.
  //
  int minMethod = PairPotentials::getInstance()->d_minMethod;
  int quadMethod = PairPotentials::getInstance()->d_quadMethod;

  int flag_elctric = 0;
  if (minMethod == 2 && quadMethod == 1)
    flag_elctric = -1;

  //
  //  compute electric field if Electrostatics is enabled
  //
  if (Electrostatics::getInstance()->isElectrostaticEnable() == 1) {
    if (flag_elctric != -1) {
      gettimeofday(&t1, NULL);

      Electrostatics::getInstance()->computeElectricField(compute_in_reference);

      gettimeofday(&t2, NULL);
      elapsedTime = (t2.tv_sec - t1.tv_sec);
      elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000000.0;
      if (RunData::getInstance()->d_timeOut == true)
        RunData::getInstance()->addComputeElectricFieldTime(0, elapsedTime);
    }
  }

  // initialize the force data
  InitializeClusterForceEnergyData();

  //  // New : trying to see how code scales
  // // with newBucketLock.
  // // To initialize, this we need to first initialize
  // // d_clusterForceEnergy
  // InitializeNewBucketLock();

  //
  //	initialize data for trace(hessian) of potential
  //
  if (quadMethod == 1)
    InitializeTraceOfKData();

  //
  //  loop over each quasi and compute force at cluster sites and
  //  nodes(electrostatics
  //	forces are directly added to nodes)
  //
  for (int iQuasi = 0; iQuasi < numQuasi; iQuasi++) {
    if (quadMethod == 1) {
      iQuasiForceEnergyHarmonic(iQuasi, compute_in_reference);
    } else {
      iQuasiForceEnergy(iQuasi, compute_in_reference);
    }
  }

  //
  // once forces are computed at cluster sites, add it to nodes for all quasis
  //
  for (int iQuasi = 0; iQuasi < numQuasi; iQuasi++) {
    iQuasiAddForceAndEnergyToNodes(iQuasi);
  }

  // //
  // //	clear the cluster force data
  // //
  // d_clusterForceEnergy.clear();

  // if(quadMethod == 1)
  // {
  // 	d_traceOfK.clear();
  // }

  // bucketLock.clear();

  return;
}

//
//  iQuasiForceEnergy()
//
void ForceEnergyCalculation::iQuasiForceEnergy(int iQuasi,
                                               bool compute_in_reference) {
  //
  //  get Quasicontinua handle
  //
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  //
  //	find if quasi is core or shell
  //
  int coreFlag = quasicontinua->getCoreShell(iQuasi);

  //
  //  first initialize the energy of all the nodes of iQuasi
  //
  quasicontinua->getQuasi(iQuasi).initializeNodeEnergy();

  //
  //  initialize the force of all the nodes of iQuasi
  //
  quasicontinua->getQuasi(iQuasi).initializeNodeForces();

  //
  //  get numQuasi
  //
  int numQuasi = quasicontinua->size();

  // time data
  timeval t1, t2;
  double elapsedTime;

  //
  //  process EAM potential
  //
  if (PairPotentials::getInstance()->d_EAMFlag == true && coreFlag == -1) {
    d_print("EAM...");
    gettimeofday(&t1, NULL);

    iQuasiEAMForceEnergy(iQuasi);

    gettimeofday(&t2, NULL);
    elapsedTime = (t2.tv_sec - t1.tv_sec);
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000000.0;
    if (RunData::getInstance()->d_timeOut == true)
      RunData::getInstance()->addComputeEAMTime(iQuasi, elapsedTime);
  }

  //
  //  process pairwise interaction
  //
  d_print("Pairwise...");
  gettimeofday(&t1, NULL);

  iQuasiPairwiseForceEnergy(iQuasi);

  gettimeofday(&t2, NULL);
  elapsedTime = (t2.tv_sec - t1.tv_sec);
  elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000000.0;
  if (RunData::getInstance()->d_timeOut == true)
    RunData::getInstance()->addComputePairwiseTime(iQuasi, elapsedTime);

  //
  // 	process Entropy's contribution to force and energy
  //
  d_print("Entropy...");
  if (coreFlag == -1) {
    gettimeofday(&t1, NULL);

    iQuasiEntropyForceEnergy(iQuasi);

    gettimeofday(&t2, NULL);
    elapsedTime = (t2.tv_sec - t1.tv_sec);
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000000.0;
    if (RunData::getInstance()->d_timeOut == true)
      RunData::getInstance()->addComputeEntropyTime(iQuasi, elapsedTime);
  }

  //
  // process Electrostatics only if it's enabled
  //
  if (Electrostatics::getInstance()->isElectrostaticEnable() == 1) {
    d_print("Electrostatics...");
    gettimeofday(&t1, NULL);

    iQuasiElectrostaticsForceEnergy(iQuasi, compute_in_reference);

    gettimeofday(&t2, NULL);
    elapsedTime = (t2.tv_sec - t1.tv_sec);
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000000.0;
    if (RunData::getInstance()->d_timeOut == true)
      RunData::getInstance()->addComputeElectrostaticsTime(iQuasi, elapsedTime);
  }

  //
  //  compute node-indentor interaction force energy
  //
  if (Indent::getInstance()->isIndentEnable() == 1) {
    d_print("Node-Indentor Interaction...");
    Indent::getInstance()->indentorNodesInteraction(iQuasi);
  }

  return;
} // end of iQuasiForceEnergy()

//
//  iQuasiForceEnergy()
//
void ForceEnergyCalculation::iQuasiForceEnergyHarmonic(
    int iQuasi, bool compute_in_reference) {
  //
  //  get Quasicontinua handle
  //
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  //
  //	find if quasi is core or shell
  //
  int coreFlag = quasicontinua->getCoreShell(iQuasi);

  //
  //  first initialize the energy of all the nodes of iQuasi
  //
  quasicontinua->getQuasi(iQuasi).initializeNodeEnergy();

  //
  //  initialize the force of all the nodes of iQuasi
  //
  quasicontinua->getQuasi(iQuasi).initializeNodeForces();

  //
  //	get the problem info
  //
  //	if ninMethod = 2 and also quadMethod = 1
  //	then we are minimizing the frequency using quasi harmonic
  //	approximation. In this case force calculation from
  // 	electrostatics is not needed.
  //
  int minMethod = PairPotentials::getInstance()->d_minMethod;
  int quadMethod = PairPotentials::getInstance()->d_quadMethod;

  //
  //  get numQuasi
  //
  int numQuasi = quasicontinua->size();

  //
  //  process EAM potential
  //
  //	Note : EAM not implemented for quasi harmonic approximation
  //
  // if(PairPotentials::getInstance()->d_EAMFlag == true && coreFlag == -1)
  // {
  // 	d_print("EAM...");
  // 	iQuasiEAMForceEnergy(iQuasi);
  // }

  //
  //  process pairwise interaction
  //
  d_print("Pairwise...");
  iQuasiPairwiseForceEnergyHarmonic(iQuasi);

  //
  // 	process Entropy's contribution to force and energy
  //
  d_print("Entropy...");
  if (coreFlag == -1)
    iQuasiEntropyForceEnergyHarmonic(iQuasi);

  //
  // process Electrostatics only if it's enabled
  //
  if (Electrostatics::getInstance()->isElectrostaticEnable() == 1) {
    // 	no need to call this fucntion if minMethod =1
    //	and quadMethod = 1
    int elect = 0;
    if (minMethod == 2 && quadMethod == 1)
      elect = -1;

    if (elect != -1) {
      d_print("Electrostatics...");
      iQuasiElectrostaticsForceEnergy(iQuasi, compute_in_reference);
    }
  }

  //
  //	we need to get traceOfK at node
  //
  ComputeTraceOfKAtNodes(iQuasi);

  //
  //  compute node-indentor interaction force energy
  //
  if (Indent::getInstance()->isIndentEnable() == 1) {
    d_print("Node-Indentor Interaction...");
    Indent::getInstance()->indentorNodesInteraction(iQuasi);
  }

  return;
} // end of iQuasiForceEnergy()

//
//  iQuasiForceEnergy()
//
void ForceEnergyCalculation::iQuasiAddForceAndEnergyToNodes(int iQuasi) {
  //	time data
  timeval t1, t2;
  double elapsedTime;

  //
  // add forces of cluster site to ndoes
  //
  gettimeofday(&t1, NULL);

  AddClusterSiteForcesToNodes(iQuasi);

  gettimeofday(&t2, NULL);
  elapsedTime = (t2.tv_sec - t1.tv_sec);
  elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000000.0;
  if (RunData::getInstance()->d_timeOut == true)
    RunData::getInstance()->addClusterForceToNodeForceTime(iQuasi, elapsedTime);

  //
  //  add external forces to node forces
  //
  iQuasiAddExternalForces(iQuasi);

  //
  //  compute total energy from all nodes and add it to all_node_list data
  //
  iQuasiComputeTotalEnergy(iQuasi);

  //
  //  compute energy due to external forces and add it to all_node_list data
  //
  iQuasiComputeExternalEnergy(iQuasi);

  // add the potential and kinetic energy to get total energy
  // std::cout<<"total energy..."<<std::endl;
  iQuasiAddEnergy(iQuasi);

  return;
} // end of iQuasiAddForceAndEnergyToNodes()

//
//  iQuasiEAMForceEnergy()
//
void ForceEnergyCalculation::iQuasiEAMForceEnergy(int iQuasi) {
  // std::cout<<"ForceEnergyCalculation : Computing EAM force and Energy for
  // Quasi = "<<iQuasi<<std::endl;

  //
  //  get CrossNeighborList handle
  //
  CrossNeighborList *crossC = CrossNeighborList::getInstance();

  //
  //  get the data from CrossNeighborList
  //
  std::vector<std::vector<neigh_site_data_t>> &iQuasi_neighbor_data =
      crossC->getQuasiNeighborData(iQuasi);

  //
  //  get Cluster data
  //
  std::vector<std::vector<std::vector<cluster_site_data_t>>> &cluster_data =
      crossC->getClusterData();

  //
  //  prepare data for thread work
  //
  struct iQuasiEAMForceEnergyData_t data;
  data.iQuasi = iQuasi;
  data.numQuasi = Quasicontinua::getInstance()->size();
  data.P_iQuasi_neighbor_data = &iQuasi_neighbor_data;
  data.P_cluster_data = &cluster_data;
  data.P_clusterForceEnergy = &d_clusterForceEnergy;

  thread_monitor(iQuasiEAMForceEnergyWorker, (void *)&data,
                 get_max_number_threads());

  // std::cout<<"ForceEnergyCalculation : Done"<<std::endl;

  return;
} // end of iQuasiEAMForceEnergy()

//
//  iQuasiPairwiseForceEnergy()
//
void ForceEnergyCalculation::iQuasiPairwiseForceEnergy(int iQuasi) {
  //
  //  get cluster data for iQuasi
  //
  std::vector<std::vector<cluster_site_data_t>> &iQuasi_cluster_data =
      CrossNeighborList::getInstance()->getQuasiClusterData(iQuasi);

  //
  // coreFlag data
  //
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();
  int numQuasi = quasicontinua->size();

  std::pair<int, int> coreFlag;
  coreFlag.second = -1;

  int i = Quasicontinua::getInstance()->getCoreShell(iQuasi);
  if (i == -1) {
    // iQuasi is shell
    coreFlag.first = -1;

    // find the core of it's shell
    for (int jQuasi = 0; jQuasi < numQuasi; jQuasi++) {
      int j = quasicontinua->getCoreShell(jQuasi);
      if (j == iQuasi) {
        // jQuasi is the core lattice of iQuasi shell lattice
        coreFlag.second = jQuasi;
      }
    }
  } else {
    // iQuasi is core
    coreFlag.first = 0;

    // put it's shell in second
    coreFlag.second = i;
  }

  //
  // prepare thread data to process cluster sites
  //
  struct iQuasiPairwiseForceEnergyData_t data;

  data.iQuasi = iQuasi;
  data.numQuasi = Quasicontinua::getInstance()->size();
  data.coreFlag = coreFlag;
  data.P_iQuasi_cluster_data = &iQuasi_cluster_data;
  data.P_clusterForceEnergy = &(d_clusterForceEnergy[iQuasi]);

  thread_monitor(iQuasiPairwiseForceEnergyWorker, (void *)&data,
                 get_max_number_threads());

  // //
  // //	Checking if data d_clusterForceEnergy and
  // //	iQuasi_cluster_data are taking very large space
  // //
  // static int testCounter = 51;
  // // printf("HEre .... \n");
  // if(testCounter < 50)
  // {
  // 	char * filename = NULL;
  // 	long int path_max;
  // 	path_max = pathconf(".", _PC_PATH_MAX);

  // 	filename = (char *) malloc(path_max*sizeof(char));
  // 	sprintf(filename,"%d_debug_check_cluster_and_force_data", testCounter);

  // 	FILE * file = fopen(filename, "w");

  // 	// cluster data
  // 	fprintf(file, "*** testing iQuasi_cluster_data iQuasi = %d***\n",
  // iQuasi);

  // 	fprintf(file, "total size = %d\n",
  // CrossNeighborList::getInstance()->getSizeOfClusterData());
  // 	fprintf(file, "size = %d\n", iQuasi_cluster_data.size());

  // 	int n = 0;
  // 	for(int i=0; i < iQuasi_cluster_data.size(); i++)
  // 	{
  // 		n = n + iQuasi_cluster_data[i].size();

  // 		//
  // 		// check inside the each data of bucket
  // 		//
  // 		for(int j = 0; j < iQuasi_cluster_data[i].size(); j++)
  // 		{
  // 			if(iQuasi_cluster_data[i][j].first.first.size() > 5)
  // 				fprintf(file, "error in site number %d of bucket %d\n", j,
  // i);

  // 			if(iQuasi_cluster_data[i][j].first.second.size() > 5)
  // 				fprintf(file, "error in site number %d of bucket %d\n", j,
  // i);

  // 			if(iQuasi_cluster_data[i][j].second.size() > 7)
  // 				fprintf(file, "error in site number %d of bucket %d\n", j,
  // i);
  // 		}
  // 	}
  // 	fprintf(file, "total number of data in all buckets = %d\n", n);

  // 	// cluster force data
  // 	fprintf(file, "*** testing d_clusterForceEnergy[iQausi] iQuasi =
  // %d***\n", iQuasi);
  // 	fprintf(file, "total size = %d\n", d_clusterForceEnergy.size());
  // 	fprintf(file, "size = %d\n", d_clusterForceEnergy[iQuasi].size());

  // 	n = 0;
  // 	for(int i = 0; i < d_clusterForceEnergy[iQuasi].size(); i++)
  // 	{
  // 		n = n+ d_clusterForceEnergy[iQuasi][i].size();

  // 		// check inside
  // 		for(int j =0; j < d_clusterForceEnergy[iQuasi][i].size(); j++)
  // 		{
  // 			if(d_clusterForceEnergy[iQuasi][i][j].first.size() > 5)
  // 				fprintf(file, "error in site number %d of bucket %d\n", j,
  // i);

  // 			// if(d_clusterForceEnergy[iQuasi][i][j].second.size() >
  // 2)
  // 			// 	fprintf(file, "error in site number %d of bucket %d\n", j,
  // i);
  // 		}
  // 	}
  // 	fprintf(file, "total number of data in all buckets = %d\n", n);

  // 	fclose(file);
  // 	free(filename);

  // 	testCounter++;
  // }

  return;
} // end of iQuasiPairwiseForceEnergy()

//
//  iQuasiPairwiseForceEnergy()
//
void ForceEnergyCalculation::iQuasiPairwiseForceEnergyHarmonic(int iQuasi) {
  //
  //  get cluster data for iQuasi
  //
  std::vector<std::vector<cluster_site_data_t>> &iQuasi_cluster_data =
      CrossNeighborList::getInstance()->getQuasiClusterData(iQuasi);

  //
  // coreFlag data
  //
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();
  int numQuasi = quasicontinua->size();

  std::pair<int, int> coreFlag;
  coreFlag.second = -1;

  int i = Quasicontinua::getInstance()->getCoreShell(iQuasi);
  if (i == -1) {
    // iQuasi is shell
    coreFlag.first = -1;

    // find the core of it's shell
    for (int jQuasi = 0; jQuasi < numQuasi; jQuasi++) {
      int j = quasicontinua->getCoreShell(jQuasi);
      if (j == iQuasi) {
        // jQuasi is the core lattice of iQuasi shell lattice
        coreFlag.second = jQuasi;
      }
    }
  } else {
    // iQuasi is core
    coreFlag.first = 0;

    // put it's shell in second
    coreFlag.second = i;
  }

  //
  // prepare thread data to process cluster sites
  //
  struct iQuasiPairwiseForceEnergyHarmonicData_t data;

  data.iQuasi = iQuasi;
  data.numQuasi = Quasicontinua::getInstance()->size();
  data.coreFlag = coreFlag;
  data.P_iQuasi_cluster_data = &iQuasi_cluster_data;
  data.P_clusterForceEnergy = &(d_clusterForceEnergy[iQuasi]);
  data.P_traceOfK = &(d_traceOfK[iQuasi]);

  thread_monitor(iQuasiPairwiseForceEnergyHarmonicWorker, (void *)&data,
                 get_max_number_threads());

  return;
} // end of iQuasiPairwiseForceEnergyHarmonic()

//
//	iQuasiEntropyForceEnergy()
//
void ForceEnergyCalculation::iQuasiEntropyForceEnergy(int iQuasi) {
  //
  //	get cluster data for iQuasi
  //
  std::vector<std::vector<cluster_site_data_t>> &iQuasi_cluster_data =
      CrossNeighborList::getInstance()->getQuasiClusterData(iQuasi);

  //
  //	prepare data for thread function
  //
  struct iQuasiEntropyForceEnergyData_t data;

  data.iQuasi = iQuasi;
  data.numQuasi = Quasicontinua::getInstance()->size();
  data.P_iQuasi_cluster_data = &iQuasi_cluster_data;
  data.P_clusterForceEnergy = &(d_clusterForceEnergy[iQuasi]);

  thread_monitor(iQuasiEntropyForceEnergyWorker, (void *)&data,
                 get_max_number_threads());

  return;
}

//
//	iQuasiEntropyForceEnergy()
//
void ForceEnergyCalculation::iQuasiEntropyForceEnergyHarmonic(int iQuasi) {
  //
  //	get cluster data for iQuasi
  //
  std::vector<std::vector<cluster_site_data_t>> &iQuasi_cluster_data =
      CrossNeighborList::getInstance()->getQuasiClusterData(iQuasi);

  //
  //	prepare data for thread function
  //
  struct iQuasiEntropyForceEnergyHarmonicData_t data;

  data.iQuasi = iQuasi;
  data.numQuasi = Quasicontinua::getInstance()->size();
  data.P_iQuasi_cluster_data = &iQuasi_cluster_data;
  data.P_clusterForceEnergy = &(d_clusterForceEnergy[iQuasi]);
  data.P_traceOfK = &(d_traceOfK[iQuasi]);

  thread_monitor(iQuasiEntropyForceEnergyHarmonicWorker, (void *)&data,
                 get_max_number_threads());

  return;
}

//
//  iQuasiElectrostaticsForceEnergy()
//
void ForceEnergyCalculation::iQuasiElectrostaticsForceEnergy(
    int iQuasi, bool compute_in_reference) {
  // std::cout<<"ForceEnergyCalculation : Computing Electrostatics force and
  // Energy for Quasi = "<<iQuasi<<std::endl;

  //
  //  get cluster data for iQuasi
  //
  std::vector<std::vector<cluster_site_data_t>> &iQuasi_cluster_data =
      CrossNeighborList::getInstance()->getQuasiClusterData(iQuasi);

  //
  //  process each cluster site in a thread
  //
  struct iQuasiElectrostaticsForceEnergyData_t data;

  data.iQuasi = iQuasi;
  data.numQuasi = Quasicontinua::getInstance()->size();
  data.compute_in_reference = compute_in_reference;
  data.P_iQuasi_cluster_data = &iQuasi_cluster_data;
  data.P_clusterForceEnergy = &(d_clusterForceEnergy[iQuasi]);

  thread_monitor(iQuasiElectrostaticsForceEnergyWorker, (void *)&data,
                 get_max_number_threads());

  // std::cout<<"ForceEnergyCalculation : Done"<<std::endl;

  return;
} // end of iQuasiElectrostaticsForceEnergy()

//
//  AddClusterSiteForcesToNodes()
//
void ForceEnergyCalculation::AddClusterSiteForcesToNodes(int iQuasi) {
  // get Quasicontinua
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();
  int numQuasi = quasicontinua->size();

  struct AddClusterSiteForcesToNodesData_t data;
  data.iQuasi = iQuasi;
  data.P_clusterData =
      &(CrossNeighborList::getInstance()->getQuasiClusterData(iQuasi));
  data.P_clusterForceData = &(d_clusterForceEnergy[iQuasi]);

  thread_monitor(AddClusterSiteForcesToNodesWorker, (void *)&data,
                 get_max_number_threads());

  return;
}

//
//
//
void ForceEnergyCalculation::ComputeTraceOfKAtNodes(int iQuasi) {
  // get Quasicontinua
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();
  int numQuasi = quasicontinua->size();

  struct ComputeTraceOfKAtNodesData_t data;
  data.iQuasi = iQuasi;
  data.P_clusterData =
      &(CrossNeighborList::getInstance()->getQuasiClusterData(iQuasi));
  data.P_traceOfK = &(d_traceOfK[iQuasi]);
  data.P_traceOfKNode = &(quasicontinua->getQuasi(iQuasi).getTraceOfKData());

  // std::vector< std::vector< cluster_site_data_t > > clusterData =
  // 	CrossNeighborList::getInstance()->getQuasiClusterData(iQuasi);

  // d_print("\ndebugging traceOfK data for iquasi = %d\n",iQuasi);
  // d_print("size of d_traceOfK[iQuasi] = %u\n", d_traceOfK[iQuasi].size());
  // if(clusterData.size() != d_traceOfK[iQuasi].size())
  // 	d_print("number of buckets in d_traceOfK (%u) and clusterData (%u) don't
  // match\n", d_traceOfK[iQuasi].size(), clusterData.size());

  // for(int i=0; i < clusterData.size(); i++)
  // 	if(clusterData[i].size() != d_traceOfK[iQuasi][i].size())
  // 		d_print("in bucket = %d number of cluster data in clusterData (%u)
  // and d_traceOfK (%u) don't match\n", i, clusterData[i].size(),
  // d_traceOfK[iQuasi][i].size());
  // d_print("test done\n");

  thread_monitor(ComputeTraceOfKAtNodesWorker, (void *)&data,
                 get_max_number_threads());

  return;
}

//
//  iQuasiAddExternalForces()
//
void ForceEnergyCalculation::iQuasiAddExternalForces(int iQuasi) {
  // std::cout<<"ForceEnergyCalculation : Adding External Forces for Quasi =
  // "<<iQuasi<<std::endl;

  // call threaded function to add the external forces to node forces
  struct node_list_t *P_node_list =
      &(Quasicontinua::getInstance()->getQuasi(iQuasi).getNodeList().node_list);

  thread_monitor(iQuasiAddExternalForcesWorker, (void *)P_node_list,
                 get_max_number_threads());

  // std::cout<<"ForceEnergyCalculation : Done"<<std::endl;

  return;
} // end of iQuasiAddExternalForces()

//
//  iQuasiComputeExternalEnergy()
//
void ForceEnergyCalculation::iQuasiComputeExternalEnergy(int iQuasi) {
  // std::cout<<"ForceEnergyCalculation : Compute External Energy for Quasi =
  // "<<iQuasi<<std::endl;

  // get the all_node_list data of iQuasi
  struct all_node_list_t &all_node_list =
      Quasicontinua::getInstance()->getQuasi(iQuasi).getNodeList();

  // call thread worker to do the job
  thread_monitor(iQuasiComputeExternalEnergyWorker, (void *)&all_node_list,
                 get_max_number_threads());

  // std::cout<<"ForceEnergyCalculation : Done"<<std::endl;

  return;
} // end of iQuasiComputeExternalEnergy()

//
//  iQuasiComputeTotalEnergy()
//
void ForceEnergyCalculation::iQuasiComputeTotalEnergy(int iQuasi) {
  // std::cout<<"ForceEnergyCalculation : Compute Total Energy for Quasi =
  // "<<iQuasi<<std::endl;

  // get all_node_list data of iQuasi
  struct all_node_list_t &all_node_list =
      Quasicontinua::getInstance()->getQuasi(iQuasi).getNodeList();

  // call thread function to do the job
  thread_monitor(iQuasiComputeTotalEnergyWorker, (void *)&all_node_list,
                 get_max_number_threads());

  // std::cout<<"ForceEnergyCalculation : Done"<<std::endl;

  return;
} // end of iQuasiComputeTotalEnergy()

//
//  iQuasiAddEnergy()
//
void ForceEnergyCalculation::iQuasiAddEnergy(int iQuasi) {
  // std::cout<<"ForceEnergyCalculation : Adding all Energy to all_node_list of
  // Quasi = "<<iQuasi<<std::endl;

  // get all_node_list of iQuasi
  struct all_node_list_t &all_node_list =
      Quasicontinua::getInstance()->getQuasi(iQuasi).getNodeList();

  // get indentor of iQuasi
  struct indentor_t indentor =
      Quasicontinua::getInstance()->getQuasi(iQuasi).getIndentor();

  // add the indentor energy to all_node_list
  all_node_list.energy.kinetic += indentor.energy.kinetic;
  all_node_list.energy.potential += indentor.energy.potential;

  // update total energy
  all_node_list.energy.total =
      all_node_list.energy.kinetic + all_node_list.energy.potential;

  // std::cout<<"ForceEnergyCalculation : Done"<<std::endl;

  return;
}

//
//  copyEnergyToInitialEnergy
//
void ForceEnergyCalculation::copyEnergyToInitialEnergy(void) {
  // loop over quasis
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  for (int iQuasi = 0; iQuasi < quasicontinua->size(); iQuasi++) {
    // get all node list
    struct all_node_list_t &all_node_list =
        quasicontinua->getQuasi(iQuasi).getNodeList();

    // copy the energy to initial energy
    memcpy(&(all_node_list.initial_energy), &(all_node_list.energy),
           sizeof(all_node_list.energy));
  }

  return;
}

//
//  InitializeClusterForceEnergyData()
//
void ForceEnergyCalculation::InitializeClusterForceEnergyData() {
  // get neighbor data and cluster data
  CrossNeighborList *crossC = CrossNeighborList::getInstance();

  const std::vector<std::vector<std::vector<cluster_site_data_t>>>
      &clusterData = crossC->getClusterData();

  int dummy = clusterData.size();
  int numBuckets = clusterData[0].size();

  int numQuasi = Quasicontinua::getInstance()->size();

  // initialize the data
  if (d_clusterForceEnergy.size() == 0 ||
      d_clusterForceEnergy.size() > numQuasi)
    d_clusterForceEnergy.resize(numQuasi);

  for (int iQuasi = 0; iQuasi < numQuasi; iQuasi++) {
    if (d_clusterForceEnergy[iQuasi].size() == 0 ||
        d_clusterForceEnergy[iQuasi].size() > numBuckets)
      d_clusterForceEnergy[iQuasi].resize(numBuckets);

    for (int bucket = 0; bucket < clusterData[iQuasi].size(); bucket++) {
      // resize the data inside bucket
      if (d_clusterForceEnergy[iQuasi][bucket].size() == 0 ||
          d_clusterForceEnergy[iQuasi][bucket].size() >
              clusterData[iQuasi][bucket].size())
        d_clusterForceEnergy[iQuasi][bucket].resize(
            clusterData[iQuasi][bucket].size());

      for (int i_cluster = 0; i_cluster < clusterData[iQuasi][bucket].size();
           i_cluster++) {
        if (d_clusterForceEnergy[iQuasi][bucket][i_cluster].first.size() == 0 ||
            d_clusterForceEnergy[iQuasi][bucket][i_cluster].first.size() > 4)
          d_clusterForceEnergy[iQuasi][bucket][i_cluster].first.resize(4);

        d_clusterForceEnergy[iQuasi][bucket][i_cluster].first[0] = 0.0;
        d_clusterForceEnergy[iQuasi][bucket][i_cluster].first[1] = 0.0;
        d_clusterForceEnergy[iQuasi][bucket][i_cluster].first[2] = 0.0;
        d_clusterForceEnergy[iQuasi][bucket][i_cluster].first[3] = 0.0;
        d_clusterForceEnergy[iQuasi][bucket][i_cluster].second = 0.0;
      }
    }
  }

  return;
}

//
//  getQuasiClusterSiteForceData()
//
std::vector<std::vector<std::pair<std::vector<double>, double>>> &
ForceEnergyCalculation::getQuasiClusterSiteForceData(int iQuasi) {
  return d_clusterForceEnergy[iQuasi];
}

//
//	testForceEnergyCalculation()
//
void ForceEnergyCalculation::testForceEnergyCalculation(
    int rebuild_neighbor_flag, int rebuild_cluster_flag,
    std::vector<std::vector<
        std::vector<std::pair<std::vector<double>, std::vector<double>>>>>
        &test_cluster_force_energy,
    std::vector<
        std::vector<std::pair<std::vector<double>, std::vector<double>>>>
        &test_node_force_energy,
    std::vector<
        std::vector<std::vector<std::pair<std::vector<double>, double>>>>
        &test_total_cluster_force_energy,
    std::vector<std::vector<std::pair<std::vector<double>, double>>>
        &test_total_node_force_energy) {
  d_print("Force Calculation...");

  int flag;
  //
  //	compute crossneighbor list if necessary
  //
  CrossNeighborList *neighborC = CrossNeighborList::getInstance();
  neighborC->computeCrossNeighborList(rebuild_neighbor_flag,
                                      rebuild_cluster_flag);
  //
  //	Initialize the data
  //
  InitializeTestForceEnergyData(
      test_cluster_force_energy, test_node_force_energy,
      test_total_cluster_force_energy, test_total_node_force_energy);

  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  int numQuasi = quasicontinua->size();

  // initialize bucket
  int numBuckets = neighborC->getNumberBuckets();
  bucketLockInitLocal(numQuasi, numBuckets);

  //
  //  compute electric field if Electrostatics is enabled
  //
  if (Electrostatics::getInstance()->isElectrostaticEnable() == 1) {
    d_print("compute electric field...");
    Electrostatics::getInstance()->computeElectricField(false);
  }

  //
  //	EAM interaction
  //
  if (PairPotentials::getInstance()->d_EAMFlag == true) {
    InitializeClusterForceEnergyData();
    for (int iQuasi = 0; iQuasi < numQuasi; iQuasi++) {
      d_print("EAM...");
      iQuasiEAMForceEnergy(iQuasi);
    }

    //
    //	put the eam force energy data in d_clusterForceEnergy to
    //	test_cluster_force_energy and test_node_force_energy
    //
    flag = 0; // 0 for EAM
    testAddForceEnergy(test_cluster_force_energy, test_node_force_energy,
                       test_total_cluster_force_energy,
                       test_total_node_force_energy, flag);
  }

  //
  //	paiwrise interaction
  //
  InitializeClusterForceEnergyData();
  for (int iQuasi = 0; iQuasi < numQuasi; iQuasi++) {
    d_print("Pairwise...");
    iQuasiPairwiseForceEnergy(iQuasi);
  }

  //
  //	put the pairwise force energy data in d_clusterForceEnergy to
  //	test_cluster_force_energy and test_node_force_energy
  //
  flag = 1; // 1 for pairwise
  testAddForceEnergy(test_cluster_force_energy, test_node_force_energy,
                     test_total_cluster_force_energy,
                     test_total_node_force_energy, flag);

  //
  // 	entropy force and energy
  //
  InitializeClusterForceEnergyData();
  for (int iQuasi = 0; iQuasi < numQuasi; iQuasi++) {
    d_print("Entropy...");
    iQuasiEntropyForceEnergy(iQuasi);
  }

  //
  //
  //
  flag = 2; // 2 for entropy
  testAddForceEnergy(test_cluster_force_energy, test_node_force_energy,
                     test_total_cluster_force_energy,
                     test_total_node_force_energy, flag);

  //
  //	electrostatics force energy
  //
  if (Electrostatics::getInstance()->isElectrostaticEnable() == 1) {
    InitializeClusterForceEnergyData();
    for (int iQuasi = 0; iQuasi < numQuasi; iQuasi++) {
      d_print("Electrostatics...");
      iQuasiElectrostaticsForceEnergy(iQuasi, false);
    }

    //
    //
    //
    flag = 3; // 3 for electrostatics
    testAddForceEnergy(test_cluster_force_energy, test_node_force_energy,
                       test_total_cluster_force_energy,
                       test_total_node_force_energy, flag);
  }

  //
  //	clear the cluster data
  //
  d_clusterForceEnergy.clear();

  d_print("done\n");
  return;
}

//
//	InitializeTestForceEnergyData()
//
void ForceEnergyCalculation::InitializeTestForceEnergyData(
    std::vector<std::vector<
        std::vector<std::pair<std::vector<double>, std::vector<double>>>>>
        &test_cluster_force_energy,
    std::vector<
        std::vector<std::pair<std::vector<double>, std::vector<double>>>>
        &test_node_force_energy,
    std::vector<
        std::vector<std::vector<std::pair<std::vector<double>, double>>>>
        &test_total_cluster_force_energy,
    std::vector<std::vector<std::pair<std::vector<double>, double>>>
        &test_total_node_force_energy) {
  // get neighbor data and cluster data
  CrossNeighborList *crossC = CrossNeighborList::getInstance();

  const std::vector<std::vector<std::vector<cluster_site_data_t>>>
      &clusterData = crossC->getClusterData();
  d_print("******\n");
  d_print("testing cluster data after updating only the neighbors\n");
  d_print("quasi size = %d, numbuckets = %d\n", clusterData.size(),
          clusterData[0].size());
  for (int k = 0; k < 10; k++) {
    for (int l = 0; l < clusterData[0][k].size(); l++)
      d_print("(x,y,z) = (%f, %f, %f)\n", clusterData[0][k][l].first.second[0],
              clusterData[0][k][l].first.second[1],
              clusterData[0][k][l].first.second[2]);
  }

  int dummy = clusterData.size();
  int numBuckets = clusterData[0].size();

  int numQuasi = Quasicontinua::getInstance()->size();

  test_cluster_force_energy.clear();
  test_cluster_force_energy.resize(numQuasi);

  test_node_force_energy.clear();
  test_node_force_energy.resize(numQuasi);

  test_total_cluster_force_energy.clear();
  test_total_cluster_force_energy.resize(numQuasi);

  test_total_node_force_energy.clear();
  test_total_node_force_energy.resize(numQuasi);

  for (int iQuasi = 0; iQuasi < numQuasi; iQuasi++) {
    // first initialize cluster data
    test_cluster_force_energy[iQuasi].resize(numBuckets);
    test_total_cluster_force_energy[iQuasi].resize(numBuckets);

    for (int bucket = 0; bucket < clusterData[iQuasi].size(); bucket++) {
      for (int i_cluster = 0; i_cluster < clusterData[iQuasi][bucket].size();
           i_cluster++) {
        std::pair<std::vector<double>, std::vector<double>> data;

        for (int i = 0; i < 13; i++) {
          data.first.push_back(0.0);
          if (i < 3)
            data.second.push_back(0.0);
        }

        test_cluster_force_energy[iQuasi][bucket].push_back(data);

        std::pair<std::vector<double>, double> data_1;

        for (int i = 0; i < 4; i++) {
          data_1.first.push_back(0.0);
        }
        data_1.second = 0.0;

        test_total_cluster_force_energy[iQuasi][bucket].push_back(data_1);
      }
    }

    // get node list for iQuasi
    struct node_list_t node_list =
        Quasicontinua::getInstance()->getQuasi(iQuasi).getNodeList().node_list;

    for (int i = 0; i < node_list.number_nodes; i++) {
      std::pair<std::vector<double>, std::vector<double>> data;
      for (int i = 0; i < 13; i++) {
        data.first.push_back(0.0);
        if (i < 3)
          data.second.push_back(0.0);
      }

      test_node_force_energy[iQuasi].push_back(data);

      std::pair<std::vector<double>, double> data_1;
      for (int i = 0; i < 4; i++) {
        data_1.first.push_back(0.0);
      }
      data_1.second = 0.0;

      test_total_node_force_energy[iQuasi].push_back(data_1);
    }
  }

  return;
}

//
//  testAddClusterForceEnergy()
//
//  flag : 0 - for eam
//         1 - for pairwise
//         2 - for entropy
//         3 - electrostatics
//
void ForceEnergyCalculation::testAddForceEnergy(
    std::vector<std::vector<
        std::vector<std::pair<std::vector<double>, std::vector<double>>>>>
        &test_cluster_force_energy,
    std::vector<
        std::vector<std::pair<std::vector<double>, std::vector<double>>>>
        &test_node_force_energy,
    std::vector<
        std::vector<std::vector<std::pair<std::vector<double>, double>>>>
        &test_total_cluster_force_energy,
    std::vector<std::vector<std::pair<std::vector<double>, double>>>
        &test_total_node_force_energy,
    int flag) {
  for (int iQuasi = 0; iQuasi < d_clusterForceEnergy.size(); iQuasi++) {
    struct testAddForceEnergyData_t data;
    data.iQuasi = iQuasi;
    data.flag = flag;
    data.P_clusterForceData = &(d_clusterForceEnergy[iQuasi]);
    data.P_clusterSiteData =
        &(CrossNeighborList::getInstance()->getQuasiClusterData(iQuasi));
    data.P_testClusterForceData = &(test_cluster_force_energy[iQuasi]);
    data.P_testNodeForceData = &(test_node_force_energy[iQuasi]);
    data.P_testTotalClusterForceData =
        &(test_total_cluster_force_energy[iQuasi]);
    data.P_testTotalNodeForceData = &(test_total_node_force_energy[iQuasi]);

    thread_monitor(testAddForceEnergyWorker, (void *)&data,
                   get_max_number_threads());
  }

  return;
}

//
//	InitializeTraceOfKData()
//
void ForceEnergyCalculation::InitializeTraceOfKData() {
  // get neighbor data and cluster data
  CrossNeighborList *crossC = CrossNeighborList::getInstance();

  const std::vector<std::vector<std::vector<cluster_site_data_t>>>
      &clusterData = crossC->getClusterData();

  int dummy = clusterData.size();
  int numBuckets = clusterData[0].size();

  int numQuasi = Quasicontinua::getInstance()->size();

  if (d_traceOfK.size() == 0 || d_traceOfK.size() > numQuasi)
    d_traceOfK.resize(numQuasi);

  for (int iQuasi = 0; iQuasi < numQuasi; iQuasi++) {
    if (d_traceOfK[iQuasi].size() == 0 ||
        d_traceOfK[iQuasi].size() > numBuckets)
      d_traceOfK[iQuasi].resize(numBuckets);

    for (int bucket = 0; bucket < clusterData[iQuasi].size(); bucket++) {
      if (d_traceOfK[iQuasi][bucket].size() == 0 ||
          d_traceOfK[iQuasi][bucket].size() >
              clusterData[iQuasi][bucket].size())
        d_traceOfK[iQuasi][bucket].resize(clusterData[iQuasi][bucket].size());

      for (int i_cluster = 0; i_cluster < clusterData[iQuasi][bucket].size();
           i_cluster++) {
        d_traceOfK[iQuasi][bucket][i_cluster] = 0.0;
      }
    }
  }

  return;
}

//
//  getQuasiClusterTraceOfK()
//
std::vector<std::vector<double>> &
ForceEnergyCalculation::getQuasiClusterTraceOfK(int iQuasi) {
  return d_traceOfK[iQuasi];
}

//
//	InitializeNewBucketLock()
//
void ForceEnergyCalculation::InitializeNewBucketLock() {
  // 	newBucketLock.clear();
  // newBucketLock.resize(d_clusterForceEnergy.size());

  // for(int iQuasi=0; iQuasi< d_clusterForceEnergy.size(); iQuasi++)
  // {
  // 	newBucketLock[iQuasi].resize(d_clusterForceEnergy[iQuasi].size());

  // 	for(int bucket=0; bucket < d_clusterForceEnergy[iQuasi].size();
  // bucket++)
  // 	{
  // 		for(int i_cluster = 0; i_cluster <
  // d_clusterForceEnergy[iQuasi][bucket].size(); i_cluster++)
  // 		{
  // 			newBucketLock[iQuasi][bucket].push_back(PTHREAD_MUTEX_INITIALIZER);
  // 		}
  // 	}
  // }

  // if(newBucketLock.size() != d_clusterForceEnergy.size())
  // 	d_print("Check initialization of newBucketLock. Quasi size don't
  // match\n");

  return;
}
}
