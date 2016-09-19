//
// Quasicontinua implementation file
//
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

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
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
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

#include <cassert>
#include <cmath>

// c++ source files
#include "C_Interface.h"
#include "CreateMesh.h"
#include "DataTypes.h"
#include "Element.h"
#include "Error.h"
#include "Indent.h"
#include "Input.h"
#include "Lattice.h"
#include "MiscFunctions.h"
#include "Node.h"
#include "PairPotentials.h"
#include "QuadraturePoints.h"
#include "Quasicontinua.h"
#include "Quasicontinuum.h"
#include "Void.h"

// c source files
#include "monitor.h"
#include "threads.h"

// static int forceF = 0;

//
//
//
namespace quasicontinuum {

//
// Constructor.
//
Quasicontinuum::Quasicontinuum(const int id, int ac, char *av[])
    : d_element_list(ELEMENT_LIST_INITIALIZER),
      d_node_list(ALL_NODE_LIST_INITIALIZER),
      d_lattice(), // FIXME: Needs to use LATTICE_INITIALIZER
      d_indentor(INDENTOR_INITIALIZER), d_qc_options(QC_OPTIONS_INITIALIZER),
      d_force_flag(0), d_energy_flag(0), d_output_flag(0),
      d_mt_version(MULTI_THREADED), d_mass(0.0) {
  //
  //  set id
  //
  d_id = id;

  //
  // process input
  //
  d_print("Quasicontinuum : Calling quasiInput()\n");
  d_qc_options = Input::getInstance()->quasiInput(
      d_id, &d_lattice, &d_node_list, &d_element_list, &d_indentor,
      d_materialName, d_mass, ac, av);

  //
  //  after setting up the node data of quasi set the fixity data
  //
  d_positionFixityOriginal.clear();
  d_frequencyFixityOriginal.clear();

  for (int iNode = 0; iNode < d_node_list.node_list.number_nodes; iNode++) {
    d_positionFixityOriginal.push_back(
        d_node_list.node_list.nodes[iNode]->fix_mask);
    d_frequencyFixityOriginal.push_back(
        d_node_list.node_list.nodes[iNode]->fix_w_mask);
  }

  //
  // change to single threaded if required
  //
  if (get_max_number_threads() == 1)
    d_mt_version = SINGLE_THREADED;

  //
  // initialize shells
  //
  d_print("Quasicontinuum : Initializing lattice shells\n");
  Lattice::getInstance()->initializeShells(d_lattice);

  //
  //  initialize d_traceOfKNode data
  //
  initializeTraceOfKData();

  //
  //
  //
  return;
}

//
// Destructor.
//
Quasicontinuum::~Quasicontinuum() {}

//
//  isRestartOn()
//
//  return 1 - if restart is on
//         0 - if restart is off
//
int Quasicontinuum::isRestartOn() const {
  //
  // check d_qc_options
  //
  if (d_qc_options.restart == ON)
    return 1;
  else
    return 0;
}

int Quasicontinuum::isRestartOn() {
  //
  // check d_qc_options
  //
  if (d_qc_options.restart == ON)
    return 1;
  else
    return 0;
}

//
// Get forces from data structures.
//
const std::vector<double> Quasicontinuum::getForces() {
  //
  // vector to hold forces
  //
  std::vector<double> myForces;

  //
  // counter for each node
  //
  int i_node;

  //
  // get statistics : 0 means zero temperature problem, 1 means constant temp
  //
  int statistics;
  statistics = PairPotentials::getInstance()->d_statistics;

  //
  // loop over each node and populate solution
  //
  for (i_node = 0; i_node < d_node_list.node_list.number_nodes; i_node++) {

    //
    // get pointer to each node
    //
    struct node_t *P_node = d_node_list.node_list.nodes[i_node];

    //
    // check each node against fixity mask
    //
    if ((P_node->fix_mask & FIX_X_MASK) == FREE_MASK) {
      // if(forceF < 50)
      //   printf("Force X here\n");
      myForces.push_back(P_node->acceleration[0]);
    }
    if ((P_node->fix_mask & FIX_Y_MASK) == FREE_MASK) {
      // if(forceF < 50)
      //   printf("Force Y here\n");
      myForces.push_back(P_node->acceleration[1]);
    }
    if ((P_node->fix_mask & FIX_Z_MASK) == FREE_MASK) {
      // if(forceF < 50)
      //   printf("Force Z here\n");
      myForces.push_back(P_node->acceleration[2]);
    }

    // force due to frequency
    if (statistics == 1)
      if ((P_node->fix_w_mask & FIX_W_MASK) == FREE_MASK) {
        // if(forceF < 50)
        //   printf("Force W here\n");
        myForces.push_back(P_node->force_frequency);
      }

    // forceF += 1;
  }

  //
  //
  //
  return myForces;
}

//
//  updateState()
//
void Quasicontinuum::updateState(const std::vector<double> &update) {
  //
  // get statistics : 0 means zero temperature problem, 1 means constant temp
  //
  int statistics;
  statistics = PairPotentials::getInstance()->d_statistics;

  //
  // counter for each node
  //
  int i_node;

  //
  // counter for update
  //
  unsigned int updateCount = 0;

  //
  // loop over each node and populate solution
  //
  for (i_node = 0; i_node < d_node_list.node_list.number_nodes; i_node++) {
    //
    // get pointer to each node
    //
    struct node_t *P_node = d_node_list.node_list.nodes[i_node];

    //
    // check each node against fixity mask
    //
    if ((P_node->fix_mask & FIX_X_MASK) == FREE_MASK) {
      //
      // update x position
      //
      P_node->position[0] = update[updateCount];

      //
      // update count
      //
      ++updateCount;
    }

    if ((P_node->fix_mask & FIX_Y_MASK) == FREE_MASK) {
      //
      // update y position
      //
      P_node->position[1] = update[updateCount];

      //
      // update count
      //
      ++updateCount;
    }

    if ((P_node->fix_mask & FIX_Z_MASK) == FREE_MASK) {
      //
      // update z position
      //
      P_node->position[2] = update[updateCount];

      //
      // update count
      //
      ++updateCount;
    }

    // update frequency after node position if statistics = 1
    if (statistics == 1) {
      if ((P_node->fix_w_mask & FIX_W_MASK) == FREE_MASK) {
        //
        // update frequency
        //
        P_node->frequency = update[updateCount];

        //
        // update count
        //
        ++updateCount;
      }
    }
  }

  //
  // compare sizes to make sure update was correct
  //
  assert(update.size() == updateCount);

  //
  //
  //
  return;
}

//
//  applyDeformation()
//
void Quasicontinuum::applyDeformation(
    const int flag, const std::vector<std::vector<double>> &def,
    const std::vector<double> &shiftVector) {

  //
  // loop over each node and populate solution
  //
  for (int i_node = 0; i_node < d_node_list.node_list.number_nodes; i_node++) {
    //
    // get pointer to each node
    //
    struct node_t *P_node = d_node_list.node_list.nodes[i_node];

    //
    // switch between whole domain and just boundaries
    //
    switch (flag) {
    //
    // whole domain
    //
    case 0:
      //
      // apply deformation
      //
      for (int iIndex = 0; iIndex < 3; ++iIndex)
        for (int jIndex = 0; jIndex < 3; ++jIndex)
          P_node->position[iIndex] +=
              def[iIndex][jIndex] *
              (P_node->initial_position[jIndex] + shiftVector[jIndex]);

      break;

    //
    // boundaries only
    //
    case 1:
      //
      // check if a dof is fixed
      //
      if (P_node->fix_mask != 0) {
        //
        // apply deformation
        //
        for (int iIndex = 0; iIndex < 3; ++iIndex)
          for (int jIndex = 0; jIndex < 3; ++jIndex)
            P_node->position[iIndex] +=
                def[iIndex][jIndex] *
                (P_node->initial_position[jIndex] + shiftVector[jIndex]);
      }

      break;
    }
  }
  //
  //
  //
  return;
}

//
// Return position.
//
const std::vector<double> Quasicontinuum::getPosition() {
  //
  // vector to return for solution
  //
  std::vector<double> myPosition;

  //
  // counter for each node
  //
  int i_node;

  //
  // loop over each node and populate solution
  //
  for (i_node = 0; i_node < d_node_list.node_list.number_nodes; i_node++) {
    //
    // get pointer to each node
    //
    struct node_t *P_node = d_node_list.node_list.nodes[i_node];

    //
    // check each node against fixity mask
    //
    if ((P_node->fix_mask & FIX_X_MASK) == FREE_MASK)
      myPosition.push_back(P_node->position[0]);
    if ((P_node->fix_mask & FIX_Y_MASK) == FREE_MASK)
      myPosition.push_back(P_node->position[1]);
    if ((P_node->fix_mask & FIX_Z_MASK) == FREE_MASK)
      myPosition.push_back(P_node->position[2]);
  }

  //
  //
  //
  return myPosition;
}

//
// Return frequency.
//
const std::vector<double> Quasicontinuum::getFrequency() {
  //
  // vector to return for solution
  //
  std::vector<double> myFrequency;

  //
  // counter for each node
  //
  int i_node;

  //
  // loop over each node and populate solution
  //
  for (i_node = 0; i_node < d_node_list.node_list.number_nodes; i_node++) {
    //
    // get pointer to each node
    //
    struct node_t *P_node = d_node_list.node_list.nodes[i_node];

    //
    // check each node against fixity mask
    //
    if ((P_node->fix_w_mask & FIX_W_MASK) == FREE_MASK)
      myFrequency.push_back(P_node->frequency);
  }

  //
  //
  //
  return myFrequency;
}

//
// Return state.
//
const std::vector<double> Quasicontinuum::getState() {
  //
  // get statistics : 0 means zero temperature problem, 1 means constant temp
  //
  int statistics;
  statistics = PairPotentials::getInstance()->d_statistics;

  //
  // vector to return for solution
  //
  std::vector<double> myState;

  //
  // counter for each node
  //
  int i_node;

  //
  // loop over each node and populate solution
  //
  for (i_node = 0; i_node < d_node_list.node_list.number_nodes; i_node++) {
    //
    // get pointer to each node
    //
    struct node_t *P_node = d_node_list.node_list.nodes[i_node];

    //
    // check each node against fixity mask
    //
    if ((P_node->fix_mask & FIX_X_MASK) == FREE_MASK)
      myState.push_back(P_node->position[0]);
    if ((P_node->fix_mask & FIX_Y_MASK) == FREE_MASK)
      myState.push_back(P_node->position[1]);
    if ((P_node->fix_mask & FIX_Z_MASK) == FREE_MASK)
      myState.push_back(P_node->position[2]);

    if (statistics == 1) {
      if ((P_node->fix_w_mask & FIX_W_MASK) == FREE_MASK) {
        // check if we enter this for core type quasi
        if (Quasicontinua::getInstance()->getCoreShell(d_id) != -1)
          d_print("for core also passing freq variable\n");

        myState.push_back(P_node->frequency);
      }
    }
  }
  //
  //
  //
  return myState;
}

//
// Return precondition values.
//
const std::vector<double> Quasicontinuum::getPreconditionValues() {
  //
  // get statistics : 0 means zero temperature problem, 1 means constant temp
  //
  int statistics;
  statistics = PairPotentials::getInstance()->d_statistics;

  //
  // vector to return for solution
  //
  std::vector<double> myPrecondition;

  //
  // counter for each node
  //
  int i_node;

  //
  // loop over each node and populate solution
  //
  for (i_node = 0; i_node < d_node_list.node_list.number_nodes; i_node++) {
    //
    // get pointer to each node
    //
    struct node_t *P_node = d_node_list.node_list.nodes[i_node];

    //
    // check each node against fixity mask
    //
    if ((P_node->fix_mask & FIX_X_MASK) == FREE_MASK)
      myPrecondition.push_back(P_node->weight);
    if ((P_node->fix_mask & FIX_Y_MASK) == FREE_MASK)
      myPrecondition.push_back(P_node->weight);
    if ((P_node->fix_mask & FIX_Z_MASK) == FREE_MASK)
      myPrecondition.push_back(P_node->weight);

    if (statistics == 1) {
      if ((P_node->fix_w_mask & FIX_W_MASK) == FREE_MASK)
        myPrecondition.push_back(0.0001 * P_node->weight);
    }
  }
  //
  //
  //
  return myPrecondition;
}

//
// get site location
//
std::pair< std::vector<int>, std::vector<double> >
Quasicontinuum::getClusterSiteLocationInNodeCluster(int iNode, int iSite) {
  //
  // site location variables
  //
  double siteLocation[3];
  int *l_i =
      &(d_node_list.node_list.nodes[iNode]->site_cluster.sites[iSite][0]);

  //
  // get spatial coordinates
  //
  Lattice::getInstance()->getSiteCurrentPosition(siteLocation, l_i, &d_lattice,
                                                 d_id);

  //
  // site location variable return
  //
  std::pair<std::vector<int>, std::vector<double>> returnSiteLocation;

  //
  // populate vector
  //
  returnSiteLocation.first.push_back(l_i[0]);
  returnSiteLocation.first.push_back(l_i[1]);
  returnSiteLocation.first.push_back(l_i[2]);
  returnSiteLocation.second.push_back(siteLocation[0]);
  returnSiteLocation.second.push_back(siteLocation[1]);
  returnSiteLocation.second.push_back(siteLocation[2]);

  //
  //
  //
  return returnSiteLocation;
}

//
// get site location
//
std::pair<std::vector<int>, std::vector<double>>
Quasicontinuum::getClusterSiteStateInNodeCluster(int iNode, int iSite) {
  //
  // site location variables
  //
  double siteState[4];
  int *l_i =
      &(d_node_list.node_list.nodes[iNode]->site_cluster.sites[iSite][0]);

  //
  // get spatial coordinates
  //
  Lattice::getInstance()->getSiteCurrentState(siteState, l_i, &d_lattice, d_id);

  //
  // site location variable return
  //
  std::pair<std::vector<int>, std::vector<double>> returnSiteState;

  //
  // populate vector
  //
  returnSiteState.first.push_back(l_i[0]);
  returnSiteState.first.push_back(l_i[1]);
  returnSiteState.first.push_back(l_i[2]);
  returnSiteState.second.push_back(siteState[0]);
  returnSiteState.second.push_back(siteState[1]);
  returnSiteState.second.push_back(siteState[2]);
  returnSiteState.second.push_back(siteState[3]);

  //
  //
  //
  return returnSiteState;
}

//
// get lattice location
//
std::vector<double> Quasicontinuum::getLatticeLocation(std::vector<int> site) {
  //
  // site location variables
  //
  double siteLocation[3];
  int l_i[3];
  l_i[0] = site[0];
  l_i[1] = site[1];
  l_i[2] = site[2];

  //
  // get spatial coordinates
  //
  Lattice::getInstance()->getSiteCurrentPosition(siteLocation, l_i, &d_lattice,
                                                 d_id);

  //
  // site location variable return
  //
  std::vector<double> returnSiteLocation;

  //
  // populate vector
  //
  returnSiteLocation.push_back(siteLocation[0]);
  returnSiteLocation.push_back(siteLocation[1]);
  returnSiteLocation.push_back(siteLocation[2]);

  //
  //
  //
  return returnSiteLocation;
}

//
// get lattice location
//
std::vector<double> Quasicontinuum::getLatticeState(std::vector<int> site) {
  //
  // site location variables
  //
  double siteState[4];
  int l_i[3];
  l_i[0] = site[0];
  l_i[1] = site[1];
  l_i[2] = site[2];

  //
  // get spatial coordinates
  //
  Lattice::getInstance()->getSiteCurrentState(siteState, l_i, &d_lattice, d_id);

  //
  // site location variable return
  //
  std::vector<double> returnSiteState;

  //
  // populate vector
  //
  returnSiteState.push_back(siteState[0]);
  returnSiteState.push_back(siteState[1]);
  returnSiteState.push_back(siteState[2]);
  returnSiteState.push_back(siteState[3]);

  //
  //
  //
  return returnSiteState;
}

//
// get node cluster weights
//
const std::vector<double> Quasicontinuum::getWeights() {
  //
  // site location variables
  //
  std::vector<double> nodeWeights;

  //
  // loop over all nodes
  //
  for (int iNode = 0; iNode < d_node_list.node_list.number_nodes; ++iNode) {
    //
    // populate vector
    //
    nodeWeights.push_back(d_node_list.node_list.nodes[iNode]->weight);
  }
  //
  //
  //
  return nodeWeights;
}

//
//  set current positions to initial positions and get actual
//  current positions
//
const std::vector<std::vector<double>>
Quasicontinuum::resetToInitialPosition(void) {

  //
  // vector of current positions
  //
  std::vector<std::vector<double>> currentPositions;

  //
  // resize currentPositions
  //
  currentPositions.resize(d_node_list.node_list.number_nodes);

  //
  // loop over all nodes
  //
  for (int iNode = 0; iNode < d_node_list.node_list.number_nodes; ++iNode) {
    //
    // pointer to current node
    //
    struct node_t *P_node = d_node_list.node_list.nodes[iNode];

    //
    // get node position
    //
    currentPositions[iNode].push_back(P_node->position[0]);
    currentPositions[iNode].push_back(P_node->position[1]);
    currentPositions[iNode].push_back(P_node->position[2]);

    //
    // set current position to initial position
    //
    P_node->position[0] = P_node->initial_position[0];
    P_node->position[1] = P_node->initial_position[1];
    P_node->position[2] = P_node->initial_position[2];
  }
  //
  //
  //
  return currentPositions;
}

//
// set current state to initial state and get actual current state
//
const std::vector<std::vector<double>>
Quasicontinuum::resetToInitialState(void) {
  //
  // vector of current state
  //
  std::vector<std::vector<double>> currentStates;

  //
  // resize currentState
  //
  currentStates.resize(d_node_list.node_list.number_nodes);

  //
  // loop over all nodes
  //
  for (int iNode = 0; iNode < d_node_list.node_list.number_nodes; ++iNode) {
    //
    // pointer to current node
    //
    struct node_t *P_node = d_node_list.node_list.nodes[iNode];

    //
    // get node position
    //
    currentStates[iNode].push_back(P_node->position[0]);
    currentStates[iNode].push_back(P_node->position[1]);
    currentStates[iNode].push_back(P_node->position[2]);
    currentStates[iNode].push_back(P_node->frequency);
    currentStates[iNode].push_back(P_node->tau);

    //
    // set current position to initial position
    //
    P_node->position[0] = P_node->initial_position[0];
    P_node->position[1] = P_node->initial_position[1];
    P_node->position[2] = P_node->initial_position[2];
    P_node->frequency = P_node->initial_frequency;
    P_node->tau = P_node->initial_tau;
  }
  //
  //
  //
  return currentStates;
}

//
// reset current positions
//
void Quasicontinuum::resetToCurrentPosition(
    const std::vector<std::vector<double>> currentPositions) {
  //
  // loop over all nodes
  //
  for (int iNode = 0; iNode < d_node_list.node_list.number_nodes; ++iNode) {
    //
    // pointer to current node
    //
    struct node_t *P_node = d_node_list.node_list.nodes[iNode];

    //
    // reset current position
    //
    P_node->position[0] = currentPositions[iNode][0];
    P_node->position[1] = currentPositions[iNode][1];
    P_node->position[2] = currentPositions[iNode][2];
  }
  //
  //
  //
  return;
}

//
// reset current positions
//
void Quasicontinuum::resetToCurrentState(
    const std::vector<std::vector<double>> currentStates) {
  //
  // loop over all nodes
  //
  for (int iNode = 0; iNode < d_node_list.node_list.number_nodes; ++iNode) {
    //
    // pointer to current node
    //
    struct node_t *P_node = d_node_list.node_list.nodes[iNode];

    //
    // reset current position
    //
    P_node->position[0] = currentStates[iNode][0];
    P_node->position[1] = currentStates[iNode][1];
    P_node->position[2] = currentStates[iNode][2];
    P_node->frequency = currentStates[iNode][3];
  }
  //
  //
  //
  return;
}

//
// find site in lattice near point by expanding shell
//
void Quasicontinuum::FindSiteInLatticeNearPoint(
    const std::vector<int> &site_index_guess,
    std::vector<int> &return_site) const {
  return_site.clear();
  return_site.resize(3);

  // check if site is in void
  // assumes lattices are all the same across multiple quasicontinuum
  // except for atoms removed due to voids
  Void *voidC = Void::getInstance();
  Lattice *latticeC = Lattice::getInstance();

  if (voidC->isVoidEnable() == 1 &&
      voidC->findSiteInVoidCache(site_index_guess, d_id) == RETURN_SUCCESS) {
    std::vector<int> test_site(3, 0);
    bool break_variable = false;

    // loop over cube of lattice sites and find a site
    for (int i = -1; i <= 1; ++i) {
      // break loop if site found
      if (break_variable == true)
        break;

      for (int j = -1; j <= 1; ++j) {
        // break loop if site found
        if (break_variable == true)
          break;

        for (int k = -1; k <= 1; ++k) {
          // new site to check
          test_site[0] = site_index_guess[0] + i;
          test_site[1] = site_index_guess[1] + j;
          test_site[2] = site_index_guess[2] + k;

          // check if test site is inside/outside sample
          if (latticeC->isSiteInsideLattice(&d_lattice, &test_site[0], d_id) ==
              RETURN_FAILURE)
            continue;

          // check if test site is a lattice site
          if (latticeC->isLatticeSite(&test_site[0], &d_lattice, d_id) ==
              RETURN_FAILURE)
            continue;

          // if here, site is good, set it
          for (int i_dof = 0; i_dof < 3; ++i_dof)
            return_site[i_dof] = test_site[i_dof];

          // set break_variable
          break_variable = true;
          break;
        } // k loop
      }   // j loop
    }     // i loop

    // if break variable == 0, no sites found
    if (break_variable == false) {
      d_print("no sites found\n");
      d_print("Site: %d   %d  %d\n", site_index_guess[0], site_index_guess[1],
              site_index_guess[2]);
      D_ERROR("FindSiteInLatticeNearPoint()");
      exit(EXIT_FAILURE);
    }
  } // end if initial site is a void

  else
    for (int i_dof = 0; i_dof < 3; ++i_dof)
      return_site[i_dof] = site_index_guess[i_dof];

  // finish
  return;
} // end of Quasicontinuum::FindSiteInLatticeNearPoint

//
//  getNeighborSitesAndStateAroundPoint()
//
std::vector<std::pair<std::vector<int>, std::vector<double>>>
Quasicontinuum::getNeighborSitesAndStateAroundPoint(
    const std::pair<std::vector<int>, std::vector<double>> siteLocation,
    const double cutoffRadius) {
  std::vector<std::pair<std::vector<int>, std::vector<double>>>
      returnSitesAndState;

  //
  // lattice coordinates of i site
  //
  int iLatticeCoordinates[3];
  iLatticeCoordinates[0] = siteLocation.first[0];
  iLatticeCoordinates[1] = siteLocation.first[1];
  iLatticeCoordinates[2] = siteLocation.first[2];

  //
  // spatial coordinates of i site
  //
  double iSpatialCoordinates[3];
  iSpatialCoordinates[0] = siteLocation.second[0];
  iSpatialCoordinates[1] = siteLocation.second[1];
  iSpatialCoordinates[2] = siteLocation.second[2];

  //  get shift
  std::vector<double> iShift = Quasicontinua::getInstance()->getShift(d_id);

  // get pointer to Lattice class
  Lattice *latticeC = Lattice::getInstance();

  //
  // lattice coordinates of j site
  //
  int jLatticeCoordinates[3];

  //
  // check if i site is lattice point in this Quasicontinuum
  //
  if (latticeC->isLatticeSite(iLatticeCoordinates, &d_lattice, d_id) ==
      RETURN_SUCCESS) {
    // set j site base coordinates
    jLatticeCoordinates[0] = iLatticeCoordinates[0];
    jLatticeCoordinates[1] = iLatticeCoordinates[1];
    jLatticeCoordinates[2] = iLatticeCoordinates[2];

  } else {
    // test site array
    int testSite[3];

    // break variable
    int breakVariable = 0;

    // loop over cube of lattice sites and find a site
    for (int i = -1; i <= 1; i++) {
      // break loop if site found
      if (breakVariable == 1)
        break;

      for (int j = -1; j <= 1; j++) {
        // break loop if site found
        if (breakVariable == 1)
          break;

        for (int k = -1; k <= 1; k++) {
          testSite[0] = iLatticeCoordinates[0] + i;
          testSite[1] = iLatticeCoordinates[1] + j;
          testSite[2] = iLatticeCoordinates[2] + k;

          // check if test site is inside/outside sample
          if (latticeC->isSiteInsideLattice(&d_lattice, testSite, d_id) ==
              RETURN_FAILURE)
            continue;

          // check if test site is a lattice site
          if (latticeC->isLatticeSite(testSite, &d_lattice, d_id) ==
              RETURN_FAILURE)
            continue;

          // set j site location
          jLatticeCoordinates[0] = testSite[0];
          jLatticeCoordinates[1] = testSite[1];
          jLatticeCoordinates[2] = testSite[2];

          // set breakVariable
          breakVariable = 1;

          break;
        }
      }
    }

    // if break variable == 0, no sites found
    if (breakVariable == 0) {
      d_print("no sites found\n");
      d_print("Site: %d   %d  %d\n", iLatticeCoordinates[0],
              iLatticeCoordinates[1], iLatticeCoordinates[2]);
      D_ERROR("getNeighborSitesAndStateAroundPoint()");
      exit(EXIT_FAILURE);
    }
  }

  //
  // spatial coordinates of site
  //
  double siteState[4];

  //
  // lattice coordinates of site
  //
  int siteLattice[3];

  //
  // shell counter
  //
  int shell_number = 1;

  //
  // counter for number of sites added
  //
  int sitesAdded = 1;

  //
  // variable to hold separation between sites
  //
  double r[3];
  double rMag;

  //
  // set site lattice coordinates
  //
  siteLattice[0] = jLatticeCoordinates[0];
  siteLattice[1] = jLatticeCoordinates[1];
  siteLattice[2] = jLatticeCoordinates[2];

  //
  // check if lattice site is inside the lattice
  //
  if (latticeC->isSiteInsideLattice(&d_lattice, siteLattice, d_id) ==
      RETURN_SUCCESS) {
    // get spatial coordinates of site
    latticeC->getSiteCurrentState(siteState, siteLattice, &d_lattice, d_id);

    //
    // add shift and get rij
    //
    for (int dof = 0; dof < 3; dof++) {
      siteState[dof] += iShift[dof];

      r[dof] = siteState[dof] - siteLocation.second[dof];
    }

    //
    // find magnitude of separation
    //
    rMag = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

    //
    // compare separation to cutoff radius
    //
    if (rMag <= cutoffRadius) {
      // move site lattice to vector
      std::pair<std::vector<int>, std::vector<double>> data;

      for (int iDof = 0; iDof < 3; ++iDof) {
        data.first.push_back(siteLattice[iDof]);
        data.second.push_back(siteState[iDof]);
      }

      data.second.push_back(siteState[3]);

      //
      // push it to return data
      //
      returnSitesAndState.push_back(data);
    }
  }

  //
  // shell variable
  //
  struct shell_t *P_shell;

  //
  // continue looping over shells until no sites are added
  //
  while (sitesAdded != 0) {
    //
    // reset sitesAdded to 0
    //
    sitesAdded = 0;

    //
    // get new shell
    //
    P_shell = latticeC->getShell(shell_number, &d_lattice);

    //
    // make sure shell exists
    //
    if (P_shell == NULL)
      exit(EXIT_FAILURE);

    //
    // loop over all sites in shell
    //
    for (int siteCounter = 0; siteCounter < P_shell->number_sites;
         ++siteCounter) {
      //
      // set site lattice coordinates
      //
      siteLattice[0] = jLatticeCoordinates[0] + P_shell->site[siteCounter][0];
      siteLattice[1] = jLatticeCoordinates[1] + P_shell->site[siteCounter][1];
      siteLattice[2] = jLatticeCoordinates[2] + P_shell->site[siteCounter][2];

      //
      // check if lattice site is inside the lattice
      //
      if (latticeC->isSiteInsideLattice(&d_lattice, siteLattice, d_id) ==
          RETURN_FAILURE)
        continue;

      //
      // get spatial coordinates of site
      //
      latticeC->getSiteCurrentState(siteState, siteLattice, &d_lattice, d_id);

      //
      // add global shift
      //
      for (int dof = 0; dof < 3; dof++) {
        siteState[dof] += iShift[dof];

        r[dof] = siteState[dof] - siteLocation.second[dof];
      }

      //
      // find magnitude of separation
      //
      rMag = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

      //
      // compare separation to cutoff radius
      //
      if (rMag <= cutoffRadius) {
        // move site lattice to vector
        std::pair<std::vector<int>, std::vector<double>> data;

        for (int iDof = 0; iDof < 3; ++iDof) {
          data.first.push_back(siteLattice[iDof]);
          data.second.push_back(siteState[iDof]);
        }

        data.second.push_back(siteState[3]);

        //
        // push it to return data
        //
        returnSitesAndState.push_back(data);

        //
        // increment sites added
        //
        sitesAdded++;
      }
    }

    //
    // increment shell number
    //
    shell_number++;
  }

  //
  //
  //
  return returnSitesAndState;

} // end of getNeighborSitesAndStateAroundPoint()

//
//  initializeExternalForces()
//
void Quasicontinuum::initializeExternalForces(void) {
  struct node_list_t node_list = d_node_list.node_list;

  //
  // loop over nodes and put zero in external forces
  //
  for (int iNode = 0; iNode < node_list.number_nodes; iNode++) {
    struct node_t *P_node = node_list.nodes[iNode];

    P_node->external_force[0] = 0.0;
    P_node->external_force[1] = 0.0;
    P_node->external_force[2] = 0.0;

    P_node->external_freq_force = 0.0;
  }

  return;
}

//
//  initializeNodeEnergy()
//
void Quasicontinuum::initializeNodeEnergy(void) {
  //
  //  loop over nodes and put zero
  //
  for (int iNode = 0; iNode < d_node_list.node_list.number_nodes; iNode++) {
    struct node_t *P_node = d_node_list.node_list.nodes[iNode];

    P_node->energy.kinetic = 0.0;
    P_node->energy.potential = 0.0;
  }

  return;
}

//
//  initializeNodeForces()
//
void Quasicontinuum::initializeNodeForces(void) {
  //
  //  loop over nodes and put zero
  //
  for (int iNode = 0; iNode < d_node_list.node_list.number_nodes; iNode++) {
    struct node_t *P_node = d_node_list.node_list.nodes[iNode];

    P_node->acceleration[0] = 0.0;
    P_node->acceleration[1] = 0.0;
    P_node->acceleration[2] = 0.0;

    P_node->force_frequency = 0.0;

    //
    // also initialize trace of K data
    //
    d_traceOfKNode[iNode] = 0.0;
  }

  return;
}

//
//  getNeighborStateAroundPoint()
//
//  Returns Glabal position and frequency of sites within cutoff radius.
//
std::vector<std::vector<double>> Quasicontinuum::getNeighborStateAroundPoint(
    const std::pair<std::vector<int>, std::vector<double>> siteLocation,
    const double cutoffRadius) {
  std::vector<std::vector<double>> returnStates;

  //
  // lattice coordinates of i site
  //
  int iLatticeCoordinates[3];
  iLatticeCoordinates[0] = siteLocation.first[0];
  iLatticeCoordinates[1] = siteLocation.first[1];
  iLatticeCoordinates[2] = siteLocation.first[2];

  //
  // spatial coordinates of i site
  //
  double iSpatialCoordinates[3];
  iSpatialCoordinates[0] = siteLocation.second[0];
  iSpatialCoordinates[1] = siteLocation.second[1];
  iSpatialCoordinates[2] = siteLocation.second[2];

  //  get shift
  std::vector<double> iShift = Quasicontinua::getInstance()->getShift(d_id);

  // get pointer to Lattice class
  Lattice *latticeC = Lattice::getInstance();

  //
  // lattice coordinates of j site
  //
  int jLatticeCoordinates[3];

  //
  // check if i site is lattice point in this Quasicontinuum
  //
  if (latticeC->isLatticeSite(iLatticeCoordinates, &d_lattice, d_id) ==
      RETURN_SUCCESS) {
    // set j site base coordinates
    jLatticeCoordinates[0] = iLatticeCoordinates[0];
    jLatticeCoordinates[1] = iLatticeCoordinates[1];
    jLatticeCoordinates[2] = iLatticeCoordinates[2];

  } else {
    // test site array
    int testSite[3];

    // break variable
    int breakVariable = 0;

    // loop over cube of lattice sites and find a site
    for (int i = -1; i <= 1; i++) {
      // break loop if site found
      if (breakVariable == 1)
        break;

      for (int j = -1; j <= 1; j++) {
        // break loop if site found
        if (breakVariable == 1)
          break;

        for (int k = -1; k <= 1; k++) {
          testSite[0] = iLatticeCoordinates[0] + i;
          testSite[1] = iLatticeCoordinates[1] + j;
          testSite[2] = iLatticeCoordinates[2] + k;

          // check if test site is inside/outside sample
          if (latticeC->isSiteInsideLattice(&d_lattice, testSite, d_id) ==
              RETURN_FAILURE)
            continue;

          // check if test site is a lattice site
          if (latticeC->isLatticeSite(testSite, &d_lattice, d_id) ==
              RETURN_FAILURE)
            continue;

          // set j site location
          jLatticeCoordinates[0] = testSite[0];
          jLatticeCoordinates[1] = testSite[1];
          jLatticeCoordinates[2] = testSite[2];

          // set breakVariable
          breakVariable = 1;

          break;
        }
      }
    }

    // if break variable == 0, no sites found
    if (breakVariable == 0) {
      d_print("no sites found\n");
      d_print("Site: %d   %d  %d \n", iLatticeCoordinates[0],
              iLatticeCoordinates[1], iLatticeCoordinates[2]);
      exit(EXIT_FAILURE);
    }
  }

  //
  // spatial coordinates of site
  //
  double siteState[4];

  //
  // lattice coordinates of site
  //
  int siteLattice[3];

  //
  // shell counter
  //
  int shell_number = 1;

  //
  // counter for number of sites added
  //
  int sitesAdded = 1;

  //
  // variable to hold separation between sites
  //
  double r[3];
  double rMag;

  //
  // set site lattice coordinates
  //
  siteLattice[0] = jLatticeCoordinates[0];
  siteLattice[1] = jLatticeCoordinates[1];
  siteLattice[2] = jLatticeCoordinates[2];

  //
  // check if lattice site is inside the lattice
  //
  if (latticeC->isSiteInsideLattice(&d_lattice, siteLattice, d_id) ==
      RETURN_SUCCESS) {
    // get spatial coordinates of site
    latticeC->getSiteCurrentState(siteState, siteLattice, &d_lattice, d_id);

    //
    // add shift and get rij
    //
    for (int dof = 0; dof < 3; dof++) {
      siteState[dof] += iShift[dof];

      r[dof] = siteState[dof] - siteLocation.second[dof];
    }

    //
    // find magnitude of separation
    //
    rMag = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

    //
    // compare separation to cutoff radius
    //
    if (rMag <= cutoffRadius) {
      // move site lattice to vector
      std::vector<double> data;

      for (int iDof = 0; iDof < 3; ++iDof)
        data.push_back(siteState[iDof]);

      data.push_back(siteState[3]);

      //
      // push it to return data
      //
      returnStates.push_back(data);
    }
  }

  //
  // shell variable
  //
  struct shell_t *P_shell;

  //
  // continue looping over shells until no sites are added
  //
  while (sitesAdded != 0) {
    //
    // reset sitesAdded to 0
    //
    sitesAdded = 0;

    //
    // get new shell
    //
    P_shell = latticeC->getShell(shell_number, &d_lattice);

    //
    // make sure shell exists
    //
    if (P_shell == NULL)
      exit(EXIT_FAILURE);

    //
    // loop over all sites in shell
    //
    for (int siteCounter = 0; siteCounter < P_shell->number_sites;
         ++siteCounter) {
      //
      // set site lattice coordinates
      //
      siteLattice[0] = jLatticeCoordinates[0] + P_shell->site[siteCounter][0];
      siteLattice[1] = jLatticeCoordinates[1] + P_shell->site[siteCounter][1];
      siteLattice[2] = jLatticeCoordinates[2] + P_shell->site[siteCounter][2];

      //
      // check if lattice site is inside the lattice
      //
      if (latticeC->isSiteInsideLattice(&d_lattice, siteLattice, d_id) ==
          RETURN_FAILURE)
        continue;

      //
      // get spatial coordinates of site
      //
      latticeC->getSiteCurrentState(siteState, siteLattice, &d_lattice, d_id);

      //
      // add global shift
      //
      for (int dof = 0; dof < 3; dof++) {
        siteState[dof] += iShift[dof];

        r[dof] = siteState[dof] - siteLocation.second[dof];
      }

      //
      // find magnitude of separation
      //
      rMag = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

      //
      // compare separation to cutoff radius
      //
      if (rMag <= cutoffRadius) {
        // move site lattice to vector
        std::vector<double> data;

        for (int iDof = 0; iDof < 3; ++iDof)
          data.push_back(siteState[iDof]);

        data.push_back(siteState[3]);

        //
        // push it to return data
        //
        returnStates.push_back(data);

        //
        // increment sites added
        //
        sitesAdded++;
      }
    }

    //
    // increment shell number
    //
    shell_number++;
  }

  //
  //
  //
  return returnStates;

} // end of getNeighborSitesAndStateAroundPoint()

//
// get sites in cluster
//
void Quasicontinuum::getSitesInCluster(
    const std::vector<int> &site_index_guess,
    const std::vector<double> &site_location,
    const double &cutoff_radius_squared,
    const std::vector<double> global_shift_vector,
    std::vector<std::vector<int>> &vector_return_sites,
    std::vector<std::vector<double>> &vector_return_locations) {
  // clear return variables
  vector_return_sites.clear();
  vector_return_locations.clear();

  // sites and locations
  std::vector<int> return_site(3, 0);
  std::vector<double> return_location(3, 0.0);

  // get good site and location
  FindSiteInLatticeNearPoint(site_index_guess, return_site);

  // get Lattice instance
  Lattice *latticeC = Lattice::getInstance();

  latticeC->getSiteCurrentPosition(&return_location[0], &return_site[0],
                                   &d_lattice, d_id);

  for (int i_dof = 0; i_dof < 3; ++i_dof)
    return_location[i_dof] += global_shift_vector[i_dof];

  // add site and location to vector of returns if inside cutoff radius
  if (MiscFunctions::getInstance()->IsPointInSphere(
          return_location, site_location, cutoff_radius_squared) == true) {
    vector_return_sites.push_back(return_site);
    vector_return_locations.push_back(return_location);
  }

  // shell counter
  int shell_number = 1;

  // counter for number of sites added
  bool sites_added = true;

  // set base site coordinates
  const std::vector<int> base_site{return_site[0], return_site[1],
                                   return_site[2]};

  // shell variable
  struct shell_t *P_shell;

  // continue looping over shells until no sites are added
  while (sites_added == true) {
    // reset sites_added to false
    sites_added = false;

    // get new shell
    P_shell = latticeC->getShell(shell_number, &d_lattice);

    // make sure shell exists
    if (P_shell == NULL) {
      D_ERROR("getSitesInCluster()");
      exit(EXIT_FAILURE);
    }

    // loop over all sites in shell
    for (int site_counter = 0; site_counter < P_shell->number_sites;
         ++site_counter) {
      // set site lattice coordinates
      for (int i_dof = 0; i_dof < 3; ++i_dof)
        return_site[i_dof] =
            base_site[i_dof] + P_shell->site[site_counter][i_dof];

      // check if site is inside the lattice
      if (latticeC->isSiteInsideLattice(&d_lattice, &return_site[0], d_id) ==
          RETURN_FAILURE)
        continue;

      // get spatial coordinates of site
      latticeC->getSiteCurrentPosition(&return_location[0], &return_site[0],
                                       &d_lattice, d_id);

      for (int i_dof = 0; i_dof < 3; ++i_dof)
        return_location[i_dof] += global_shift_vector[i_dof];

      // add site and location to vector of returns if inside cutoff radius
      if (MiscFunctions::getInstance()->IsPointInSphere(
              return_location, site_location, cutoff_radius_squared) == true) {
        vector_return_sites.push_back(return_site);
        vector_return_locations.push_back(return_location);
        sites_added = true;
      }
    } // end loop over all sites in shell

    // increment shell number
    shell_number++;
  } // end of while loop

  //
  return;
} // Quasicontinuum::GetSitesInCluster

//
//  modifyInitialConfiguration()
//
void Quasicontinuum::modifyInitialConfiguration() {
  d_print("here modifyInitialConfiguration()\n");
  // loop over nodes
  for (int iNode = 0; iNode < d_node_list.node_list.number_nodes; iNode++) {
    struct node_t *P_node = d_node_list.node_list.nodes[iNode];

    P_node->initial_position[0] = P_node->position[0];
    P_node->initial_position[1] = P_node->position[1];
    P_node->initial_position[2] = P_node->position[2];
    P_node->initial_frequency = P_node->frequency;
    P_node->initial_tau = P_node->tau;
  }

  d_node_list.initial_energy.kinetic = d_node_list.energy.kinetic;
  d_node_list.initial_energy.potential = d_node_list.energy.potential;
  d_node_list.initial_energy.total =
      d_node_list.initial_energy.kinetic + d_node_list.initial_energy.potential;

  //
  return;
} // end of modifyInitialConfiguration()

//
//  set or reset position fixity mask of all nodes
//
//  fixity_flag = 0 - reset the fixity mask
//                1 - set fixity mask so that all
//                    position dof are fixed
//                2 - set fixity mask so that all
//                    position dof are free
//
void Quasicontinuum::setPositionFixity(const int fixity_flag) {
  switch (fixity_flag) {
  case 0: {
    // check size of d_positionFixityOriginal
    if (d_positionFixityOriginal.size() != 
        static_cast<unsigned>(d_node_list.node_list.number_nodes)) {
      d_print(
          "Check d_positionFixityOriginal and function setPositionFixity()\n");
      exit(EXIT_FAILURE);
    }

    // replace the fix_mask of all nodes with d_positionFixityOriginal

    struct node_list_t *P_node_list = &(d_node_list.node_list);

    for (int iNode = 0; iNode < P_node_list->number_nodes; iNode++) {
      struct node_t *P_node = P_node_list->nodes[iNode];

      P_node->fix_mask = d_positionFixityOriginal[iNode];
    }
  }

  break;

  case 1: {
    // fix all position dof of nodes
    int fixity;
    fixity = FREE_MASK;
    fixity |= FIX_X_MASK | FIX_Y_MASK | FIX_Z_MASK;

    //
    // also put the recent fix_mask of nodes to
    // d_positionFixityOriginal so that later
    // it can be used to revert the fix_mask
    // to original value
    //
    d_positionFixityOriginal.clear();

    struct node_list_t *P_node_list = &(d_node_list.node_list);

    for (int iNode = 0; iNode < P_node_list->number_nodes; iNode++) {
      struct node_t *P_node = P_node_list->nodes[iNode];

      // put it to backup
      d_positionFixityOriginal.push_back(P_node->fix_mask);

      //  modify it now
      P_node->fix_mask = fixity;
    }
  }

  break;

  case 2: {
    // free all position dof of nodes
    int fixity;
    fixity = FREE_MASK;

    //
    // also put the recent fix_mask of nodes to
    // d_positionFixityOriginal so that later
    // it can be used to revert the fix_mask
    // to original value
    //
    d_positionFixityOriginal.clear();

    struct node_list_t *P_node_list = &(d_node_list.node_list);

    for (int iNode = 0; iNode < P_node_list->number_nodes; iNode++) {
      struct node_t *P_node = P_node_list->nodes[iNode];

      // put it to backup
      d_positionFixityOriginal.push_back(P_node->fix_mask);

      //  modify it now
      P_node->fix_mask = fixity;
    }
  }

  break;

  default: {
    // it should not reach here
    d_print("check fixity_flag = %d in setPositionFixity()\n");
    exit(EXIT_FAILURE);
  }

  break;
  }

  return;
}

//
//  set or reset frequency fixity mask of all nodes
//
//  fixity_flag = 0 - reset the fixity mask
//                1 - set fixity mask so that all
//                    frequency dof are fixed
//                2 - set fixity mask so that all
//                    frequency dof are free
//
void Quasicontinuum::setFrequencyFixity(const int fixity_flag) {
  switch (fixity_flag) {
  case 0: {
    // check size of d_frequencyFixityOriginal
    if (d_frequencyFixityOriginal.size() !=
        static_cast<unsigned>(d_node_list.node_list.number_nodes)) {
      d_print(
          "Check d_frequencyFixityOriginal and function setPositionFixity()\n");
      exit(EXIT_FAILURE);
    }

    // replace the fix_mask of all nodes with d_frequencyFixityOriginal

    struct node_list_t *P_node_list = &(d_node_list.node_list);

    for (int iNode = 0; iNode < P_node_list->number_nodes; iNode++) {
      struct node_t *P_node = P_node_list->nodes[iNode];

      P_node->fix_w_mask = d_frequencyFixityOriginal[iNode];
    }
  }

  break;

  case 1: {
    // fix all frequency dof of nodes
    int fixity;
    fixity = FREE_MASK;
    fixity |= FIX_W_MASK;

    //
    // also put the recent fix_mask of nodes to
    // d_frequencyFixityOriginal so that later
    // it can be used to revert the fix_w_mask
    // to original value
    //
    d_frequencyFixityOriginal.clear();

    struct node_list_t *P_node_list = &(d_node_list.node_list);

    for (int iNode = 0; iNode < P_node_list->number_nodes; iNode++) {
      struct node_t *P_node = P_node_list->nodes[iNode];

      // put it to backup
      d_frequencyFixityOriginal.push_back(P_node->fix_w_mask);

      //  modify it now
      P_node->fix_w_mask = fixity;
    }
  }

  break;

  case 2: {
    // free all frequency dof of nodes
    int fixity;
    fixity = FREE_MASK;

    //
    // also put the recent fix_w_mask of nodes to
    // d_frequencyFixityOriginal so that later
    // it can be used to revert the fix_w_mask
    // to original value
    //
    d_frequencyFixityOriginal.clear();

    struct node_list_t *P_node_list = &(d_node_list.node_list);

    for (int iNode = 0; iNode < P_node_list->number_nodes; iNode++) {
      struct node_t *P_node = P_node_list->nodes[iNode];

      // put it to backup
      d_frequencyFixityOriginal.push_back(P_node->fix_w_mask);

      //  modify it now
      P_node->fix_w_mask = fixity;
    }
  }

  break;

  default: {
    // it should not reach here
    d_print("check fixity_flag = %d in setFrequencyFixity()\n");
    exit(EXIT_FAILURE);
  }

  break;
  }

  return;
}

//
//  getTraceOfKData()
//
const std::vector<double> &Quasicontinuum::getTraceOfKData() const {
  return d_traceOfKNode;
}

std::vector<double> &Quasicontinuum::getTraceOfKData() {
  return d_traceOfKNode;
}

//
// InitializeTraceOfKData()
//
void Quasicontinuum::initializeTraceOfKData() {
  d_traceOfKNode.clear();
  for (int i = 0; i < d_node_list.node_list.number_nodes; i++)
    d_traceOfKNode.push_back(0.0);

  return;
}
}
