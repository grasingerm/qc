//
// Electrostatics.cc
//

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

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

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#else
#error No pthread.h available.
#endif /* HAVE_PTHREAD_H */

#include <algorithm>
#include <cassert>
#include <cmath>
//#include <ctime>
#include <cstdio>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <tuple>
#include <utility>

// CGAL library
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Object.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/intersections.h>

#include <CGAL/Cartesian.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/convex_hull_traits_2.h>

//
#include "C_Interface.h"
#include "CrossNeighborList.h"
#include "DataTypes.h"
#include "Electrostatics.h"
#include "Element.h"
#include "Error.h"
#include "Input.h"
#include "Lattice.h"
#include "MiscFunctions.h"
#include "PairPotentials.h"
#include "QuadraturePoints.h"
#include "Quasicontinua.h"
#include "Quasicontinuum.h"
#include "Shape.h"
#include "Void.h"

#include "monitor.h"
#include "threads.h"

//
//
//

namespace quasicontinuum {

//
//  namespace for data used in threading
//
namespace {
typedef std::vector<double> My_point;
typedef std::vector<My_point> My_triangle;

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Point_3 Point_3;
typedef K::Segment_3 Segment_3;
typedef K::Triangle_3 Triangle_3;

static pthread_once_t bucket_once = PTHREAD_ONCE_INIT;
const static int numBuckets = 65535;
static std::vector<pthread_mutex_t> bucketLocks;
static pthread_mutex_t dataLock = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t sumLock = PTHREAD_MUTEX_INITIALIZER;

static double shiftTolerance = 0.0001;

static int sortVariable = 0;

static int currentElem = 0;
static int currentBucket = 0;
static int currentNode = 0;

const static int elemBucket = 5;
const static int chargeBucket = 5;
const static int nodeBucket = 5;

volatile static int whileWorker = 1;

//
//  thread data structures
//
struct computeFaces_argument_t {
  std::vector<std::vector<std::vector<int>>> *faces;
  std::vector<std::vector<std::pair<int, int>>> *faceElementList;
  struct element_list_t *elementList;
};

struct computePolarizationSites_argument_t {
  std::vector<std::vector<int>> *sites;
  struct element_list_t *elementList;
  struct lattice_t *lattice;
  int *quasi;
};

struct computeAtomisticElements_argument_t {
  std::vector<int> *atomisticElements;
  struct element_list_t *elementList;
  double *atomisticElementSize;
};

struct computeAtomisticElementLocations_argument_t {
  std::vector<int> *atomisticElements;
  std::vector<std::vector<std::vector<double>>> *atomisticElementLocations;
  struct element_list_t *elementList;
};

struct computeAtomisticNodes_argument_t {
  std::vector<int> *atomisticElements;
  std::vector<std::pair<int, std::vector<int>>> *atomisticNodes;
  struct node_list_t *nodeList;
};

struct computeAtomisticCharges_argument_t {
  std::vector<int> *atomisticElements;
  std::vector<std::vector<std::vector<double>>> *atomisticElementLocations;
  std::vector<std::vector<std::vector<int>>> *atomisticCharges;
  struct lattice_t *lattice;
  struct element_list_t *elementList;
  int *quasi;
};

struct computePolarizationLocations_argument_t {
  std::vector<int> *atomisticElements;
  std::vector<std::vector<int>> *elementPolarizationSites;
  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      *elementPolarizationLocations;
  int quasiId;
  double fixedCharge;
  struct lattice_t *lattice;
};

struct computeElementPolarization_argument_t {
  std::vector<int> *atomisticElements;
  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      *elementPolarizationLocations;
  double initialUnitCellVolume;
  struct element_list_t *elementList;
  std::vector<std::vector<double>> *elementPolarization;
};

struct computeChargeDensities_argument_t {
  std::vector<std::vector<std::vector<int>>> *faces;
  std::vector<std::vector<std::pair<int, int>>> *faceElementList;
  struct node_list_t *nodeList;
  std::vector<std::vector<std::pair<double, std::vector<std::vector<double>>>>>
      *chargeDensity;
  std::vector<std::vector<double>> *elementPolarization;
  struct lattice_t *lattice;
  std::vector<std::pair<int, double>> *electroBoundaries;
};

struct computeAtomisticChargeState_argument_t {
  std::vector<std::vector<std::vector<int>>> *atomisticCharges;
  std::vector<std::vector<std::vector<double>>> *atomisticChargeState;
  int quasiId;
  struct lattice_t *lattice;
};

struct computeFieldFromAtomisticCharges_argument_t {
  struct node_list_t *nodeList;
  std::vector<std::vector<double>> *electricField;
  std::vector<double> *electricPotential;
  std::vector<double> *fixedCharge;
  double *electricConstant;
  std::vector<std::vector<std::vector<double>>> *atomisticChargeState;
  std::vector<double> *iShift;
  int *iQuasi;
};

struct computeClusterCharges_argument_t {
  struct node_list_t *i_node_list;
  std::vector<std::vector<double>> *electricField;
  std::vector<double> *electricPotential;
  double *electricConstant;
  std::vector<std::vector<std::vector<int>>> *atomisticCharges;
  struct lattice_t *j_lattice;
  std::vector<double> *i_shift;
  std::vector<double> *j_shift;
  int *num_shells;
  double *j_charge;
  int *j_quasi;
  int *i_quasi;
  std::pair<int, int> *coreData;
};

struct compute_field_from_external_charges_argument_t {
  struct node_list_t *node_list;
  std::vector<std::vector<double>> *electric_field;
  std::vector<double> *electric_potential;
  double *electric_constant;
  std::vector<std::pair<double, std::vector<double>>> *external_charges;
  std::vector<double> *iShift;
  int *iQuasi;
  int *iCoreFlag;
};

struct computeChargedFace_argument_t {
  struct node_list_t *i_node_list;
  std::vector<std::vector<double>> *electricField;
  std::vector<double> *electricPotential;
  double *electricConstant;
  std::vector<double> *i_shift;
  std::vector<double> *zero_shift;
  int *num_shells;
  int *face_method;
  double *integration_error;
  struct node_list_t *zero_node_list;
  std::vector<std::vector<std::vector<int>>> *faces;
  std::vector<std::vector<std::pair<double, std::vector<std::vector<double>>>>>
      *charge_density;
  struct lattice_t *zero_lattice;
  std::vector<std::vector<std::vector<int>>> *atomistic_charges;
  int *zero_quasi;
};

struct computeClusterChargeSurface_argument_t {
  struct node_list_t *i_node_list;
  std::vector<std::vector<double>> *electricField;
  std::vector<double> *electricPotential;
  double *electricConstant;
  std::vector<std::vector<std::vector<int>>> *atomisticCharges;
  std::vector<std::vector<double>> *element_polarization;
  struct lattice_t *zero_lattice;
  std::vector<double> *i_shift;
  std::vector<double> *zero_shift;
  int *num_shells;
  int *face_method;
  double *integration_error;
  struct node_list_t *zero_node_list;
  std::vector<std::vector<std::pair<std::vector<int>, int>>>
      *base_node_hash_list;
  int *zero_quasi;
};

} // end of namespace for data used in threading

//
//  namespace for local functions used in thread functions
//
namespace {
bool sortFunction(const std::vector<int> i, const std::vector<int> j) {
  if (i[0] < j[0])
    return true;
  else if (i[0] == j[0] && i[1] < j[1])
    return true;
  else if (i[0] == j[0] && i[1] == j[1] && i[2] < j[2])
    return true;
  else
    return false;
}

bool uniqueFunction(const std::vector<int> i, const std::vector<int> j) {
  if (i[0] == j[0] && i[1] == j[1] && i[2] == j[2])
    return true;
  else
    return false;
}

void bucket_init(void) {
  // loop over all locks and initialize
  for (int i = 0; i < numBuckets; ++i)
    bucketLocks.push_back(PTHREAD_MUTEX_INITIALIZER);

  //
  return;
}

//
//  CantorPair()
//
int CantorPair(int a, int b) {
  //
  // key
  //
  int key = (a + b) * (a + b + 1) / 2 + b;

  //
  //
  //
  return key;
}

//
//  getTripletIntegersKey()
//
unsigned int getTripletIntegersKey(int a, int b, int c) {
  // check for negative values
  if (a < 0)
    a = -a;
  if (b < 0)
    a = -b;
  if (c < 0)
    a = -c;

  // get key for a and b
  const int d = CantorPair(a, b);

  // get key for c and d
  int key = CantorPair(c, d) % numBuckets;

  // check for negative value if key overruns (unlikely)
  if (key < 0)
    key = -key;

  //
  return key;
}

//
//  tetrahedronVolume()
//
double tetrahedronVolume(double a[3], double b[3], double c[3], double d[3]) {
  //
  // call MiscFunctions to compute volume
  //
  double volume = MiscFunctions::getInstance()->tetraVolume(a, b, c, d);

  //
  // make sure volume is positive
  //
  if (volume < 0.0)
    volume = volume * -1.0;

  //
  //
  //
  return volume;
}

//
// find normal
//
std::vector<double>
CalculateNormal(std::vector<std::vector<double>> nodeLocations) {
  //
  // get difference vectors
  //
  std::vector<double> vec1(3, 0.0);
  std::vector<double> vec2(3, 0.0);

  //
  // get vecs
  //
  for (int iDof = 0; iDof < 3; ++iDof) {
    vec1[iDof] = nodeLocations[1][iDof] - nodeLocations[0][iDof];
    vec2[iDof] = nodeLocations[2][iDof] - nodeLocations[0][iDof];
  }

  //
  // get cross product
  //
  std::vector<double> normal =
      MiscFunctions::getInstance()->crossProduct3x3(vec1, vec2);
  //
  // get length of normal
  //
  double norm = MiscFunctions::getInstance()->L2Norm(normal);

  //
  // make normal unit length
  //
  for (int iDof = 0; iDof < 3; ++iDof)
    normal[iDof] = normal[iDof] / norm;

  //
  //
  //
  return normal;
}

//
//  insertChargeIntoBucket()
//
void insertChargeIntoBucket(std::vector<int> site,
                            std::vector<std::vector<int>> &bucket) {
  //
  // bucket iterator
  //
  std::vector<std::vector<int>>::iterator bucketIterator;

  //
  // search bucket
  //
  for (unsigned int it = 0; it != bucket.size() + 1; ++it) {
    //
    // check if iterator is at end
    //
    if (it == bucket.size()) {
      //
      // insert charge into bucket
      //
      bucketIterator = bucket.end();
      bucket.insert(bucketIterator, site);

      //
      // break out of for loop
      //
      break;
    }

    //
    // check if charge matches
    //
    if (bucket[it][0] == site[0] && bucket[it][1] == site[1] &&
        bucket[it][2] == site[2]) {
      //
      // break out of for loop
      //
      break;
    }

    //
    // check if have passed point of where charge would be
    //
    if (bucket[it][0] > site[0]) {
      //
      // insert charge into bucket
      //
      bucketIterator = bucket.begin() + it;
      bucket.insert(bucketIterator, site);

      //
      // break out of for loop
      //
      break;
    }
  }

  //
  //
  //
  return;
}

//
// check if site is in element
//
bool isInElement(const int site[3],
                 std::vector<std::vector<double>> nodeLocations,
                 struct lattice_t &lattice, const int &quasi) {

  //
  // get nodes in arrays
  //
  const double node1[3] = {nodeLocations[0][0], nodeLocations[0][1],
                           nodeLocations[0][2]};
  const double node2[3] = {nodeLocations[1][0], nodeLocations[1][1],
                           nodeLocations[1][2]};
  const double node3[3] = {nodeLocations[2][0], nodeLocations[2][1],
                           nodeLocations[2][2]};
  const double node4[3] = {nodeLocations[3][0], nodeLocations[3][1],
                           nodeLocations[3][2]};

  //
  // location variable
  //
  double location[3];

  //
  // get site location
  //
  Lattice::getInstance()->getSiteInitialPosition(location, site, &lattice,
                                                 quasi);

  //
  // check tetrahedron
  //
  if (MiscFunctions::getInstance()->checkPointInTetra(
          node1, node2, node3, node4, location) == INSIDE)
    return (true);

  //
  //
  //
  return false;
}

//
//  sitesInElement(()
//
//  find all sites inside element, possibility of returning one atomistic site
//  that is not in element
//
std::vector<std::vector<int>>
sitesInElement(const std::vector<int> elementSite, const int iElem,
               const std::vector<std::vector<std::vector<double>>>
                   &atomisticElementLocations,
               std::vector<int> &atomisticElements, struct lattice_t &lattice,
               int &quasi) {
  //
  // variable to hold sites
  //
  std::vector<std::vector<int>> sitesInElement;

  //
  // array for site
  //
  const int site[3] = {elementSite[0], elementSite[1], elementSite[2]};

  //
  // add site
  //
  sitesInElement.push_back(elementSite);

  //
  // counter for number of sites added
  //
  int sitesAdded = 1;

  //
  // shell variable
  //
  struct shell_t *P_shell;

  //
  // lattice coordinates of site
  //
  int siteLattice[3];

  //
  // shell counter
  //
  int shell_number = 1;

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
    P_shell = Lattice::getInstance()->getShell(shell_number, &lattice);

    //
    // check to see if shell exists
    //
    if (P_shell == NULL) {
      d_print("Shell not found compute electrostatics: siteInElement()\n");
      exit(EXIT_FAILURE);
    }

    //
    // loop over all sites in shell
    //
    for (int siteCounter = 0; siteCounter < P_shell->number_sites;
         ++siteCounter) {
      //
      // set site lattice coordinates
      //
      siteLattice[0] = site[0] + P_shell->site[siteCounter][0];
      siteLattice[1] = site[1] + P_shell->site[siteCounter][1];
      siteLattice[2] = site[2] + P_shell->site[siteCounter][2];

      //
      // check if lattice site is inside the lattice
      //
      if (Lattice::getInstance()->isSiteInsideLattice(
              &lattice, siteLattice, quasi) == RETURN_SUCCESS) {
        //
        // see if site is in element
        //
        if (isInElement(siteLattice, atomisticElementLocations[iElem], lattice,
                        quasi) == true) {
          //
          // increment sites added
          //
          sitesAdded++;

          //
          // add site
          //
          std::vector<int> tempSite;
          for (int iSite = 0; iSite < 3; ++iSite)
            tempSite.push_back(siteLattice[iSite]);

          sitesInElement.push_back(tempSite);
        }
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
  return sitesInElement;
}

//
// insert sites into list
//
void InsertSitesIntoList(std::vector<std::vector<int>> &site_list,
                         const std::vector<std::vector<int>> &sites_to_insert) {
  // loop over all sites to insert
  for (int i_site = 0; i_site < sites_to_insert.size(); ++i_site) {
    // list iterator
    std::vector<std::vector<int>>::iterator site_list_iterator;

    // search list
    for (unsigned int it = 0; it != site_list.size() + 1; ++it) {
      // check if iterator is at end
      if (it == site_list.size()) {
        // insert site into list
        site_list_iterator = site_list.end();
        site_list.insert(site_list_iterator, sites_to_insert[i_site]);

        // break out of for loop
        break;
      }

      // check if site matches
      if (site_list[it][0] == sites_to_insert[i_site][0] &&
          site_list[it][1] == sites_to_insert[i_site][1] &&
          site_list[it][2] == sites_to_insert[i_site][2]) {
        // break out of for loop
        break;
      }

      // check if have passed point of where site would be
      if (site_list[it][0] > sites_to_insert[i_site][0]) {
        // insert charge into site_list
        site_list_iterator = site_list.begin() + it;
        site_list.insert(site_list_iterator, sites_to_insert[i_site]);

        // break out of for loop
        break;
      }
    } // end for loop over all sites in list
  }   // end for loop over sites to insert

  //
  return;
} // end InsertSitesIntoList

//
// insert node site and number into hash list
//
void InsertNodeIntoHashList(
    const std::vector<int> &site, const int node_number,
    std::vector<std::pair<std::vector<int>, int>> &bucket) {
  // bucket iterator
  std::vector<std::pair<std::vector<int>, int>>::iterator bucket_iterator;

  // search bucket
  for (unsigned int it = 0; it != bucket.size() + 1; ++it) {
    // check if iterator is at end
    if (it == bucket.size()) {
      // insert site into bucket
      bucket_iterator = bucket.end();
      const std::pair<std::vector<int>, int> temp_pair =
          std::make_pair(site, node_number);
      bucket.insert(bucket_iterator, temp_pair);

      // break out of for loop
      break;
    }

    // check if site matches
    if (bucket[it].first[0] == site[0] && bucket[it].first[1] == site[1] &&
        bucket[it].first[2] == site[2]) {
      // break out of for loop
      break;
    }

    // check if have passed point of where site would be
    if (bucket[it].first[0] > site[0]) {
      // insert site into bucket
      bucket_iterator = bucket.begin() + it;

      const std::pair<std::vector<int>, int> temp_pair =
          std::make_pair(site, node_number);

      bucket.insert(bucket_iterator, temp_pair);

      // break out of for loop
      break;
    }
  } // end of for loop over sites in bucket

  //
  return;
} // end of InsertNodeIntoHashList

//
// find outward normal of plane given fourth interior point
//
std::vector<double>
outwardNormal(std::vector<std::vector<double>> nodeLocations,
              std::vector<double> fourthNode) {
  //
  // get difference vectors
  //
  std::vector<double> vec1(3, 0.0);
  std::vector<double> vec2(3, 0.0);

  //
  // get vecs
  //
  for (int iDof = 0; iDof < 3; ++iDof) {
    vec1[iDof] = nodeLocations[1][iDof] - nodeLocations[0][iDof];
    vec2[iDof] = nodeLocations[2][iDof] - nodeLocations[0][iDof];
  }

  //
  // get cross product
  //
  std::vector<double> normal =
      MiscFunctions::getInstance()->crossProduct3x3(vec1, vec2);

  //
  // get length of normal
  //
  double norm = MiscFunctions::getInstance()->L2Norm(normal);

  //
  // make normal unit length
  //
  for (int iDof = 0; iDof < 3; ++iDof)
    normal[iDof] = normal[iDof] / norm;

  //
  // check normal against 4th point
  //
  double check = 0.0;
  for (int iDof = 0; iDof < 3; ++iDof)
    check += (fourthNode[iDof] - nodeLocations[0][iDof]) * normal[iDof];

  //
  // check if greater than 0
  //
  if (check > 0)
    for (int iDof = 0; iDof < 3; ++iDof)
      normal[iDof] *= -1.0;

  //
  //
  //
  return normal;
}

//
// add contribution from atom to electric field
//
void AddAtomicToFieldAndPotential(const std::vector<double> &iFieldState,
                                  const int &iQuasi,
                                  const std::vector<double> &jSourceState,
                                  const int &jQuasi,
                                  const double &electric_constant,
                                  const double &source_charge,
                                  std::vector<double> &field, double &potential,
                                  std::pair<int, int> coreData) {
  //
  //  compute field at iFieldState due to charge at jSourceState
  //

  //
  //  Four cases :
  //  Case 1 : shell - shell interaction coreData = (-1,-1)
  //  Case 2 : shell - core  interaction coreData = (-1,0)
  //  Case 3 : core  - shell  interaction coreData = (0,-1)
  //  Case 4 : core  - core  interaction coreData = (0,0)
  //

  int flag = -1;

  if (coreData.first == -1 && coreData.second == -1)
    flag = 0;
  if (coreData.first == -1 && coreData.second != -1)
    flag = 1;
  if (coreData.first != -1 && coreData.second == -1)
    flag = 2;
  if (coreData.first != -1 && coreData.second != -1)
    flag = 3;

  //
  //  if we are doing quasi harmonic approximation, then
  //  calculation will be done with respect to mean position
  //  and there will be freq_force as Trace of K for electrostatics
  //  interaction is zero. See NOTES.
  //
  int minMethod = PairPotentials::getInstance()->d_minMethod;
  int quadMethod = PairPotentials::getInstance()->d_quadMethod;
  if (quadMethod == 1)
    flag = 3; // only mean position will play role.

  if (flag == -1) {
    pthread_mutex_lock(&dataLock);
    d_print("check core data in AddAtomicToFieldAndPotential()\n");
    exit(EXIT_FAILURE);
  }

  switch (flag) {
  case 0: {
    // shell - shell

    //
    //  get Quadrature points for pairwise interaction between two charges
    //
    std::vector<std::pair<std::vector<std::vector<double>>, double>>
        quadVectors =
            QuadraturePoints::getInstance()->getQuadratureVectors(2, 3);

    //
    //  get Sigmavectors
    //
    std::vector<double> sigmaVector =
        Quasicontinua::getInstance()->getSigmaVector();

    //
    //  factor_quad = integration factor
    //  factor_iSite = sqrt(2) /sigma / frequency for iSite
    //
    double factor_quad = pow(1 / sqrt(M_PI), 1.0 * 3 * 2);
    double factor_iSite = sigmaVector[iQuasi] / (iFieldState[3]);
    double factor_jSite = sigmaVector[jQuasi] / (jSourceState[3]);

    int quad;
    for (quad = 0; quad < quadVectors.size(); quad++) {
      //
      //  get two vectors and weight
      //
      std::vector<double> quadV_1 = quadVectors[quad].first[0];
      std::vector<double> quadV_2 = quadVectors[quad].first[1];
      double quadWeight = quadVectors[quad].second;

      //
      //  compute r_hat_ij
      //
      std::vector<double> r_hat_ij;

      for (int dof = 0; dof < 3; dof++)
        r_hat_ij.push_back(iFieldState[dof] - jSourceState[dof]);

      // magnitude of r_hat_ij
      double r_hat_mag =
          sqrt(r_hat_ij[0] * r_hat_ij[0] + r_hat_ij[1] * r_hat_ij[1] +
               r_hat_ij[2] * r_hat_ij[2]);

      //
      //  write it to electric field and potential
      //
      double r_hat_3 = r_hat_mag * r_hat_mag * r_hat_mag;

      // process only if two atoms are not identical
      if (r_hat_mag > 0.001) {
        for (int dof = 0; dof < 3; dof++)
          r_hat_ij[dof] +=
              factor_iSite * quadV_1[dof] - factor_jSite * quadV_2[dof];

        r_hat_mag = sqrt(r_hat_ij[0] * r_hat_ij[0] + r_hat_ij[1] * r_hat_ij[1] +
                         r_hat_ij[2] * r_hat_ij[2]);

        for (int dof = 0; dof < 3; dof++)
          field[dof] += electric_constant * factor_quad * source_charge *
                        r_hat_ij[dof] * quadWeight / r_hat_3;

        double r_hat_dot = r_hat_ij[0] * quadV_1[0] + r_hat_ij[1] * quadV_1[1] +
                           r_hat_ij[2] * quadV_1[2];

        field[3] += electric_constant * factor_quad * source_charge *
                    (-sigmaVector[iQuasi] / (iFieldState[3] * iFieldState[3])) *
                    r_hat_dot * quadWeight / r_hat_3;

        potential += electric_constant * factor_quad * source_charge *
                     quadWeight / r_hat_mag;
      }
    }

    return;

    break;
  }

  case 1: {
    // shell - core

    //
    //  get Quadrature points for pairwise interaction between two charges
    //
    std::vector<std::pair<std::vector<std::vector<double>>, double>>
        quadVectors =
            QuadraturePoints::getInstance()->getQuadratureVectors(1, 3);

    //
    //  get Sigmavectors
    //
    std::vector<double> sigmaVector =
        Quasicontinua::getInstance()->getSigmaVector();

    //
    //  factor_quad = integration factor
    //  factor_iSite = sqrt(2) /sigma / frequency for iSite
    //
    double factor_quad = pow(1 / sqrt(M_PI), 1.0 * 3);
    double factor_iSite = sigmaVector[iQuasi] / (iFieldState[3]);

    int quad;
    for (quad = 0; quad < quadVectors.size(); quad++) {
      //
      //  get two vectors and weight
      //
      std::vector<double> quadV_1 = quadVectors[quad].first[0];
      double quadWeight = quadVectors[quad].second;

      //
      //  compute r_hat_ij
      //
      std::vector<double> r_hat_ij;

      for (int dof = 0; dof < 3; dof++)
        r_hat_ij.push_back(iFieldState[dof] - jSourceState[dof]);

      // magnitude of r_hat_ij
      double r_hat_mag =
          sqrt(r_hat_ij[0] * r_hat_ij[0] + r_hat_ij[1] * r_hat_ij[1] +
               r_hat_ij[2] * r_hat_ij[2]);

      //
      //  write it to electric field and potential
      //
      double r_hat_3 = r_hat_mag * r_hat_mag * r_hat_mag;

      // process only if two atoms are not identical
      // Also avoiding shell-core interaction of same atom
      if (r_hat_mag > 0.001) {
        for (int dof = 0; dof < 3; dof++)
          r_hat_ij[dof] += factor_iSite * quadV_1[dof];

        r_hat_mag = sqrt(r_hat_ij[0] * r_hat_ij[0] + r_hat_ij[1] * r_hat_ij[1] +
                         r_hat_ij[2] * r_hat_ij[2]);

        for (int dof = 0; dof < 3; dof++)
          field[dof] += electric_constant * factor_quad * source_charge *
                        r_hat_ij[dof] * quadWeight / r_hat_3;

        double r_hat_dot = r_hat_ij[0] * quadV_1[0] + r_hat_ij[1] * quadV_1[1] +
                           r_hat_ij[2] * quadV_1[2];

        field[3] += electric_constant * factor_quad * source_charge *
                    (-sigmaVector[iQuasi] / (iFieldState[3] * iFieldState[3])) *
                    r_hat_dot * quadWeight / r_hat_3;

        potential += electric_constant * factor_quad * source_charge *
                     quadWeight / r_hat_mag;
      }
    }
    return;

    break;
  }

  case 2: {
    // core - shell

    //
    //  get Quadrature points for pairwise interaction between two charges
    //
    std::vector<std::pair<std::vector<std::vector<double>>, double>>
        quadVectors =
            QuadraturePoints::getInstance()->getQuadratureVectors(1, 3);

    //
    //  get Sigmavectors
    //
    std::vector<double> sigmaVector =
        Quasicontinua::getInstance()->getSigmaVector();

    //
    //  factor_quad = integration factor
    //  factor_iSite = sqrt(2) /sigma / frequency for iSite
    //
    double factor_quad = pow(1 / sqrt(M_PI), 1.0 * 3);
    double factor_jSite = sigmaVector[jQuasi] / (jSourceState[3]);

    int quad;
    for (quad = 0; quad < quadVectors.size(); quad++) {
      //
      //  get two vectors and weight
      //
      std::vector<double> quadV_1 = quadVectors[quad].first[0];
      double quadWeight = quadVectors[quad].second;

      //
      //  compute r_hat_ij
      //
      std::vector<double> r_hat_ij;

      for (int dof = 0; dof < 3; dof++)
        r_hat_ij.push_back(iFieldState[dof] - jSourceState[dof]);

      // magnitude of r_hat_ij
      double r_hat_mag =
          sqrt(r_hat_ij[0] * r_hat_ij[0] + r_hat_ij[1] * r_hat_ij[1] +
               r_hat_ij[2] * r_hat_ij[2]);

      //
      //  write it to electric field and potential
      //
      double r_hat_3 = r_hat_mag * r_hat_mag * r_hat_mag;

      // process only if two atoms are not identical
      // Also avoiding shell-core interaction of same atom
      if (r_hat_mag > 0.001) {
        for (int dof = 0; dof < 3; dof++)
          r_hat_ij[dof] -= factor_jSite * quadV_1[dof];

        r_hat_mag = sqrt(r_hat_ij[0] * r_hat_ij[0] + r_hat_ij[1] * r_hat_ij[1] +
                         r_hat_ij[2] * r_hat_ij[2]);

        for (int dof = 0; dof < 3; dof++)
          field[dof] += electric_constant * factor_quad * source_charge *
                        r_hat_ij[dof] * quadWeight / r_hat_3;

        // Note : there is no force corresponding to frequency for core
        // atoms

        potential += electric_constant * factor_quad * source_charge *
                     quadWeight / r_hat_mag;
      }
    }

    return;

    break;
  }

  case 3: {
    // core - core

    //
    // for core - core interaction it's only nteraction between mean
    // position of each core
    //

    //
    //  compute r_hat_ij
    //
    std::vector<double> r_hat_ij;

    for (int dof = 0; dof < 3; dof++)
      r_hat_ij.push_back(iFieldState[dof] - jSourceState[dof]);

    // magnitude of r_hat_ij
    double r_hat_mag =
        sqrt(r_hat_ij[0] * r_hat_ij[0] + r_hat_ij[1] * r_hat_ij[1] +
             r_hat_ij[2] * r_hat_ij[2]);

    //
    //  write it to electric field and potential
    //
    double r_hat_3 = r_hat_mag * r_hat_mag * r_hat_mag;

    // process only if two atoms are not identical
    // Also avoiding shell-core interaction of same atom
    if (r_hat_mag > 0.1) {
      for (int dof = 0; dof < 3; dof++)
        field[dof] +=
            electric_constant * source_charge * r_hat_ij[dof] / r_hat_3;

      // Note : there is no force corresponding to frequency for core
      // atoms

      potential += electric_constant * source_charge / r_hat_mag;
    }

    return;

    break;
  }

  default: {
    //
    // error
    //
    d_print("Exit: AddAtomicToFieldAndPotential() failure\n");
    exit(EXIT_FAILURE);
  }
  }

  return;

  return;
} // end of AddAtomicToFieldAndPotential

//
// add contribution from atom to electric field
//
void AddAtomicToFieldAndPotentialExternal(
    const std::vector<double> &iFieldState, const int &iQuasi,
    const std::vector<double> &jExternalCharge, const double &electric_constant,
    const double &source_charge, std::vector<double> &field, double &potential,
    int iCoreFlag) {
  //
  // depending on iCoreFlag either do phase average of charge-charge
  //  interaction or just compute interaction of external charge
  //  with charge at mean position of core
  //
  if (iCoreFlag == -1) {
    // iQuasi is of shell type

    //
    //   get Quadrature points for 3 dimensional integration
    //
    std::vector<std::pair<std::vector<std::vector<double>>, double>>
        quadVectors =
            QuadraturePoints::getInstance()->getQuadratureVectors(1, 3);

    std::vector<double> sigmaVector =
        Quasicontinua::getInstance()->getSigmaVector();

    int quad;
    double factor_quad = pow(1 / sqrt(M_PI), 1.0 * 3);
    double factor_iSite = sigmaVector[iQuasi] / (iFieldState[3]);

    for (quad = 0; quad < quadVectors.size(); quad++) {
      //
      //  get two vectors and weight
      //
      std::vector<double> quadVector = quadVectors[quad].first[0];
      double quadWeight = quadVectors[quad].second;

      //
      //  compute r_hat_ij
      //
      std::vector<double> r_hat_ij;

      for (int dof = 0; dof < 3; dof++)
        r_hat_ij.push_back(iFieldState[dof] + factor_iSite * quadVector[dof] -
                           jExternalCharge[dof]);

      double r_hat_mag =
          sqrt(r_hat_ij[0] * r_hat_ij[0] + r_hat_ij[1] * r_hat_ij[1] +
               r_hat_ij[2] * r_hat_ij[2]);

      double r_hat_3 = r_hat_mag * r_hat_mag * r_hat_mag;

      if (r_hat_mag > 0.1) {
        // add contribution to field and potential
        for (int dof = 0; dof < 3; dof++)
          field[dof] += electric_constant * factor_quad * source_charge *
                        r_hat_ij[dof] * quadWeight / r_hat_3;

        double r_hat_dot = r_hat_ij[0] * quadVector[0] +
                           r_hat_ij[1] * quadVector[1] +
                           r_hat_ij[2] * quadVector[2];

        field[3] += electric_constant * factor_quad * source_charge *
                    (-sigmaVector[iQuasi] / (iFieldState[3] * iFieldState[3])) *
                    r_hat_dot * quadWeight / r_hat_3;

        potential += electric_constant * factor_quad * source_charge *
                     quadWeight / r_hat_mag;
      }
    }
  } else {
    //
    //  compute r_hat_ij
    //
    std::vector<double> r_hat_ij;

    for (int dof = 0; dof < 3; dof++)
      r_hat_ij.push_back(iFieldState[dof] - jExternalCharge[dof]);

    // magnitude of r_hat_ij
    double r_hat_mag =
        sqrt(r_hat_ij[0] * r_hat_ij[0] + r_hat_ij[1] * r_hat_ij[1] +
             r_hat_ij[2] * r_hat_ij[2]);

    //
    //  write it to electric field and potential
    //
    double r_hat_3 = r_hat_mag * r_hat_mag * r_hat_mag;

    // process only if two atoms are not identical
    // Also avoiding shell-core interaction of same atom
    if (r_hat_mag > 0.1) {
      for (int dof = 0; dof < 3; dof++)
        field[dof] +=
            electric_constant * source_charge * r_hat_ij[dof] / r_hat_3;

      // Note : there is no force corresponding to frequency for core
      // atoms

      potential += electric_constant * source_charge / r_hat_mag;
    }
  }

  return;
} // end of AddAtomicToFieldAndPotentialExternal

//
// check if site is in atomistic list
//
bool isSiteAtomistic(
    const std::vector<int> site,
    const std::vector<std::vector<std::vector<int>>> &atomisticCharges) {
  //
  // get bucket
  //
  unsigned int bucket = getTripletIntegersKey(site[0], site[1], site[2]);

  //
  // search bucket
  //
  for (unsigned int it = 0; it != atomisticCharges[bucket].size() + 1; ++it) {
    //
    // check if iterator is at end
    //
    if (it == atomisticCharges[bucket].size()) {
      //
      // return false
      //
      return false;
    }

    //
    // check if charge matches
    //
    if (atomisticCharges[bucket][it][0] == site[0] &&
        atomisticCharges[bucket][it][1] == site[1] &&
        atomisticCharges[bucket][it][2] == site[2]) {
      //
      // return true
      //
      return true;
    }

    //
    // check if have passed point of where charge would be
    //
    if (atomisticCharges[bucket][it][0] > site[0]) {
      //
      // return false
      //
      return false;
    }
  }

  //
  // not reached
  //
  return false;
}

//
// check if site is in node list
//
const bool
IsSiteANode(const std::vector<int> &site,
            const std::vector<std::pair<std::vector<int>, int>> &bucket,
            int &node_number) {
  // search list
  for (unsigned int i_node = 0; i_node < bucket.size(); ++i_node) {
    // check if site matches
    if (bucket[i_node].first[0] == site[0] &&
        bucket[i_node].first[1] == site[1] &&
        bucket[i_node].first[2] == site[2]) {
      //
      node_number = bucket[i_node].second;
      return true;
    }

    // check if have passed point of where element would be
    if (bucket[i_node].first[0] > site[0]) {
      //
      return false;
    }
  } // end of for loop over bucket of sites

  //
  return false;
} // end of IsSiteANode

//
// midpoint of line
//
std::vector<double> midPointLine(const std::vector<double> &location1,
                                 const std::vector<double> &location2) {
  //
  // midpoint
  //
  std::vector<double> midPoint(3, 0.0);
  for (int iDof = 0; iDof < 3; ++iDof)
    midPoint[iDof] = (location1[iDof] + location2[iDof]) / 2.0;

  //
  // return
  //
  return midPoint;
}

//
// divide triangle into 4 triangles on plane
//
std::vector<std::vector<std::vector<double>>>
divideTrianglePlane(const std::vector<double> location1,
                    const std::vector<double> location2,
                    const std::vector<double> location3) {
  //
  // m12
  //
  std::vector<double> midLocation12 = midPointLine(location1, location2);

  //
  // m23
  //
  std::vector<double> midLocation23 = midPointLine(location2, location3);

  //
  // m13
  //
  std::vector<double> midLocation13 = midPointLine(location1, location3);

  std::vector<std::vector<std::vector<double>>> triangles;
  std::vector<std::vector<double>> triangle;

  triangle.push_back(location1);
  triangle.push_back(midLocation12);
  triangle.push_back(midLocation13);
  triangles.push_back(triangle);

  triangle.clear();
  triangle.push_back(location2);
  triangle.push_back(midLocation12);
  triangle.push_back(midLocation23);
  triangles.push_back(triangle);

  triangle.clear();
  triangle.push_back(location3);
  triangle.push_back(midLocation13);
  triangle.push_back(midLocation23);
  triangles.push_back(triangle);

  triangle.clear();
  triangle.push_back(midLocation13);
  triangle.push_back(midLocation12);
  triangle.push_back(midLocation23);
  triangles.push_back(triangle);

  return triangles;
}

//
// calculate field from plane charge
//
std::pair<std::vector<double>, double> fieldFromTriangularFace(
    double chargeDensity, std::vector<std::vector<double>> planeLocations,
    std::vector<double> fieldLocation, double electricConstant, int method) {
#if 0
            std::vector<double> fieldDummy;
            double potentialDummy;
            std::ofstream outfile;
            outfile.open("quadratureConvergence.txt",std::ios_base::app);
            //outfile<<"Next Charge: "<<std::endl;
            for(method =0; method<6; ++method){
#endif

  //
  // field and potential
  //
  std::vector<double> field(3, 0.0);
  double potential = 0.0;
  std::vector<std::vector<double>> localQuadPoints;
  std::vector<double> localWeights;

  // get MiscFunctions instance
  MiscFunctions *miscC = MiscFunctions::getInstance();

  //
  // switch between different calculation methods
  //
  switch (method) {
  //
  // method 0, aggregate charge and treat as point charge
  //
  case 0: {
    //
    // get area of triangle
    //
    double area = miscC->TriangleArea(planeLocations[0], planeLocations[1],
                                      planeLocations[2]);

    //
    // get center of triangle
    //
    std::vector<double> center = miscC->TriangleCenter(
        planeLocations[0], planeLocations[1], planeLocations[2]);

    //
    // get r
    //
    std::vector<double> r;
    for (int iDof = 0; iDof < 3; ++iDof)
      r.push_back(fieldLocation[iDof] - center[iDof]);

    //
    // get r2
    //
    double r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];

    //
    // get r3
    //
    double r1 = std::sqrt(r2);
    double r3 = r2 * r1;

    //
    // add contribution to field
    //
    for (int iDof = 0; iDof < 3; ++iDof)
      field[iDof] = electricConstant * chargeDensity * area * r[iDof] / r3;

    //
    // add to potential
    //
    potential += electricConstant * chargeDensity * area / r1;
#if 0
            if(potential > 0.0001 || potential < -0.0001){
              outfile<<potential;
            }
#endif
    return std::make_pair(field, potential);
  }

  //
  // method 1, single quadrature point
  //
  case 1: {
    //
    // insert single quad point
    //
    std::vector<double> localLocation;
    localLocation.push_back(0.333333333333333);
    localLocation.push_back(0.333333333333333);
    localLocation.push_back(0.333333333333333);
    localQuadPoints.push_back(localLocation);

    //
    // insert local quad point weight
    //
    localWeights.push_back(1.0);

    break;
  }

  //
  // method 2, three quadrature points
  //
  case 2: {
    //
    // insert quad points
    //
    std::vector<double> localLocation;
    localLocation.push_back(0.666666666666667);
    localLocation.push_back(0.166666666666667);
    localLocation.push_back(0.166666666666667);
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.166666666666667;
    localLocation[1] = 0.666666666666667;
    localLocation[2] = 0.166666666666667;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.166666666666667;
    localLocation[1] = 0.166666666666667;
    localLocation[2] = 0.666666666666667;
    localQuadPoints.push_back(localLocation);

    //
    // insert local quad point weight
    //
    localWeights.push_back(0.333333333333333);
    localWeights.push_back(0.333333333333333);
    localWeights.push_back(0.333333333333333);

    break;
  }

  //
  // method 3, six quadrature points
  //
  case 3: {
    //
    // insert quad points
    //
    std::vector<double> localLocation;
    localLocation.push_back(0.816847572980459);
    localLocation.push_back(0.091576213509771);
    localLocation.push_back(0.091576213509771);
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.091576213509771;
    localLocation[1] = 0.816847572980459;
    localLocation[2] = 0.091576213509771;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.091576213509771;
    localLocation[1] = 0.091576213509771;
    localLocation[2] = 0.816847572980459;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.108103018168070;
    localLocation[1] = 0.445948490915965;
    localLocation[2] = 0.445948490915965;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.445948490915965;
    localLocation[1] = 0.108103018168070;
    localLocation[2] = 0.445948490915965;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.445948490915965;
    localLocation[1] = 0.445948490915965;
    localLocation[2] = 0.108103018168070;
    localQuadPoints.push_back(localLocation);

    //
    // insert local quad point weight
    //
    localWeights.push_back(0.109951743655322);
    localWeights.push_back(0.109951743655322);
    localWeights.push_back(0.109951743655322);
    localWeights.push_back(0.223381589678011);
    localWeights.push_back(0.223381589678011);
    localWeights.push_back(0.223381589678011);

    break;
  }

  //
  // method 4, seven quadrature points
  //
  case 4: {
    //
    // insert quad points
    //
    std::vector<double> localLocation;
    localLocation.push_back(0.333333333333333);
    localLocation.push_back(0.333333333333333);
    localLocation.push_back(0.333333333333333);
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.797426985353087;
    localLocation[1] = 0.101286507323456;
    localLocation[2] = 0.101286507323456;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.101286507323456;
    localLocation[1] = 0.797426985353087;
    localLocation[2] = 0.101286507323456;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.101286507323456;
    localLocation[1] = 0.101286507323456;
    localLocation[2] = 0.797426985353087;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.470142064105115;
    localLocation[1] = 0.470142064105115;
    localLocation[2] = 0.059715871789770;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.470142064105115;
    localLocation[1] = 0.059715871789770;
    localLocation[2] = 0.470142064105115;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.059715871789770;
    localLocation[1] = 0.470142064105115;
    localLocation[2] = 0.470142064105115;
    localQuadPoints.push_back(localLocation);

    //
    // insert local quad point weight
    //
    localWeights.push_back(0.225000000000000);
    localWeights.push_back(0.125939180544827);
    localWeights.push_back(0.125939180544827);
    localWeights.push_back(0.125939180544827);
    localWeights.push_back(0.132394152788506);
    localWeights.push_back(0.132394152788506);
    localWeights.push_back(0.132394152788506);

    break;
  }

  //
  // method 5, twelve quadrature points
  //
  case 5: {
    //
    // insert quad points
    //
    std::vector<double> localLocation;
    localLocation.push_back(0.873821971016996);
    localLocation.push_back(0.063089014491502);
    localLocation.push_back(0.063089014491502);
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.063089014491502;
    localLocation[1] = 0.873821971016996;
    localLocation[2] = 0.063089014491502;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.063089014491502;
    localLocation[1] = 0.063089014491502;
    localLocation[2] = 0.873821971016996;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.501426509658179;
    localLocation[1] = 0.249286745170910;
    localLocation[2] = 0.249286745170910;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.249286745170910;
    localLocation[1] = 0.501426509658179;
    localLocation[2] = 0.249286745170910;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.249286745170910;
    localLocation[1] = 0.249286745170910;
    localLocation[2] = 0.501426509658179;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.636502499121399;
    localLocation[1] = 0.310352451033785;
    localLocation[2] = 0.053145049844816;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.636502499121399;
    localLocation[1] = 0.053145049844816;
    localLocation[2] = 0.310352451033785;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.053145049844816;
    localLocation[1] = 0.636502499121399;
    localLocation[2] = 0.310352451033785;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.053145049844816;
    localLocation[1] = 0.310352451033785;
    localLocation[2] = 0.636502499121399;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.310352451033785;
    localLocation[1] = 0.053145049844816;
    localLocation[2] = 0.636502499121399;
    localQuadPoints.push_back(localLocation);

    localLocation[0] = 0.310352451033785;
    localLocation[1] = 0.636502499121399;
    localLocation[2] = 0.053145049844816;
    localQuadPoints.push_back(localLocation);

    //
    // insert local quad point weight
    //
    localWeights.push_back(0.050844906370207);
    localWeights.push_back(0.050844906370207);
    localWeights.push_back(0.050844906370207);
    localWeights.push_back(0.116786275726379);
    localWeights.push_back(0.116786275726379);
    localWeights.push_back(0.116786275726379);
    localWeights.push_back(0.082851075618374);
    localWeights.push_back(0.082851075618374);
    localWeights.push_back(0.082851075618374);
    localWeights.push_back(0.082851075618374);
    localWeights.push_back(0.082851075618374);
    localWeights.push_back(0.082851075618374);

    break;
  }

  default: { assert(method == -1); }
  }

  //
  // get area of triangle
  //
  double J = miscC->TriangleArea(planeLocations[0], planeLocations[1],
                                 planeLocations[2]);

  //
  // loop over all quad points
  //
  for (unsigned int iQuad = 0; iQuad < localQuadPoints.size(); ++iQuad) {
    //
    // global coordinates
    //
    std::vector<double> globalLocation(3, 0.0);

    //
    // get quad global coordinates
    //
    for (int iDof = 0; iDof < 3; ++iDof)
      for (int iShape = 0; iShape < 3; ++iShape)
        globalLocation[iDof] +=
            planeLocations[iShape][iDof] * localQuadPoints[iQuad][iShape];

    //
    // get r
    //
    std::vector<double> r;
    for (int iDof = 0; iDof < 3; ++iDof)
      r.push_back(fieldLocation[iDof] - globalLocation[iDof]);

    //
    // get r2
    //
    double r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];

    //
    // get r3
    //
    double r1 = std::sqrt(r2);
    double r3 = r2 * r1;

    //
    // add contribution to field
    //
    for (int iDof = 0; iDof < 3; ++iDof)
      field[iDof] += electricConstant * chargeDensity * J *
                     localWeights[iQuad] * r[iDof] / r3;

    //
    // add to potential
    //
    potential +=
        electricConstant * chargeDensity * J * localWeights[iQuad] / r1;
  }

  return std::make_pair(field, potential);
} // end of fieldFromTriangularFace()

//
//  static int shell_count = 1;
//  static int integration_level = 1;
//
//  calculate field from density face adaptively
//
std::pair<std::vector<double>, double> fieldFromDensityFaceAdaptive(
    const double &density,
    const std::vector<std::vector<double>> &planeLocations,
    const std::vector<double> &fieldLocation, const double &electricConstant,
    const int &faceMethod, const double &integrationError) {
  //
  // calculate field from nodes
  //
  std::vector<double> fieldOld(3, std::numeric_limits<double>::max());

  //
  // all levels of triangles
  //
  std::vector<std::vector<std::vector<std::vector<double>>>> trianglesAllLevels;

  //
  // triangles to loop over
  //
  std::vector<std::vector<std::vector<double>>> triangles;
  triangles.push_back(planeLocations);

  //
  // initial triangles
  //
  trianglesAllLevels.push_back(triangles);

  //
  // loop over subdivision levels
  //
  int nLevel = 0;
  std::vector<double> fieldFromPlanes(3, 0.0);
  double potentialFromPlanes = 0.0;

  while (std::abs(fieldFromPlanes[0] - fieldOld[0]) > integrationError ||
         std::abs(fieldFromPlanes[1] - fieldOld[1]) > integrationError ||
         std::abs(fieldFromPlanes[2] - fieldOld[2]) > integrationError) {

    if (nLevel == 5)
      break;

    // while(nLevel < integration_level){

    trianglesAllLevels.resize(nLevel + 1);

    //
    // subdivide if level greater than zero
    //
    if (nLevel > 0) {
      //
      // set old field values
      //
      for (int iDof = 0; iDof < 3; ++iDof)
        fieldOld[iDof] = fieldFromPlanes[iDof];

      //
      // subdivide triangles
      //
      for (int iTriangle = 0;
           iTriangle < int(trianglesAllLevels[nLevel - 1].size());
           ++iTriangle) {
        //
        // get new triangles
        //
        std::vector<std::vector<std::vector<double>>> tempTriangles =
            divideTrianglePlane(trianglesAllLevels[nLevel - 1][iTriangle][0],
                                trianglesAllLevels[nLevel - 1][iTriangle][1],
                                trianglesAllLevels[nLevel - 1][iTriangle][2]);

        for (int iTri = 0; iTri < int(tempTriangles.size()); ++iTri)
          trianglesAllLevels[nLevel].push_back(tempTriangles[iTri]);
      }
    }

    //
    // initialize variables to zero
    //
    fieldFromPlanes[0] = 0.0;
    fieldFromPlanes[1] = 0.0;
    fieldFromPlanes[2] = 0.0;
    potentialFromPlanes = 0.0;

    //
    // calculate values from subdivided triangles
    //
    for (int iTriangle = 0; iTriangle < int(trianglesAllLevels[nLevel].size());
         ++iTriangle) {
      //
      // field from triangular face
      //
      const std::pair<std::vector<double>, double> tempVar =
          fieldFromTriangularFace(density,
                                  trianglesAllLevels[nLevel][iTriangle],
                                  fieldLocation, electricConstant, faceMethod);

      for (int iDof = 0; iDof < 3; ++iDof)
        fieldFromPlanes[iDof] += tempVar.first[iDof];

      potentialFromPlanes += tempVar.second;
    }

    //
    // increment level
    //
    nLevel++;
  }

  //
  // return field
  //
  return std::make_pair(fieldFromPlanes, potentialFromPlanes);
}

//
//  AddTriangleCluster()
//
void AddTriangleCluster(
    const int quasi, const std::vector<double> &node_location,
    const std::vector<int> site_1, const std::vector<int> site_2,
    const std::vector<int> site_3, const double &electric_constant,
    struct lattice_t &lattice,
    const std::vector<std::vector<double>> &element_polarization,
    const std::vector<std::vector<std::vector<int>>> &atomistic_charges,
    const double &integration_error, const int &face_method,
    const double &test_shift, const std::vector<double> &face_shift,
    const int &const_index, const struct node_list_t &zero_node_list,
    const std::vector<std::vector<std::pair<std::vector<int>, int>>>
        &node_hash_list,
    std::vector<double> &electric_field, double &electric_potential) {
  // check if all sites are atomistic, return, contribution included in
  // atomistics
  if (isSiteAtomistic(site_1, atomistic_charges) == true &&
      isSiteAtomistic(site_2, atomistic_charges) == true &&
      isSiteAtomistic(site_3, atomistic_charges) == true)
    return;

  // get Lattice instance
  Lattice *latticeC = Lattice::getInstance();
  MiscFunctions *miscC = MiscFunctions::getInstance();

  // element to get polarization from
  int element_for_polarization = -1;

  // get site locations without shift
  std::vector<double> location_1(3, 0.0);
  std::vector<double> location_2(3, 0.0);
  std::vector<double> location_3(3, 0.0);
  latticeC->getSiteCurrentPosition(&location_1[0], &site_1[0], &lattice, quasi);
  latticeC->getSiteCurrentPosition(&location_2[0], &site_2[0], &lattice, quasi);
  latticeC->getSiteCurrentPosition(&location_3[0], &site_3[0], &lattice, quasi);

  // check if site is a node, if so find element that point is in
  const unsigned int bucket =
      getTripletIntegersKey(site_1[0], site_1[1], site_1[2]);

  int node_number;
  const bool site_is_a_node =
      IsSiteANode(site_1, node_hash_list[bucket], node_number);

  // get site point to check
  std::vector<double> point_to_check =
      miscC->vectorScale(MiscFunctions::getInstance()->vectorAddition(
                             location_1, location_2, location_3),
                         1.0 / 3.0);

  point_to_check[const_index] += test_shift;

  // find element if site is a node
  if (site_is_a_node == true) {
    // get element list
    const struct element_list_t &element_list =
        zero_node_list.nodes[node_number]->element_list;

    // loop over elements
    for (int i_element = 0; i_element < element_list.number_elements;
         ++i_element) {
      // get element
      const struct element_t *P_element = element_list.elements[i_element];

      // check if point is in element
      if (miscC->checkPointInTetra(
              P_element->node[0]->position, P_element->node[1]->position,
              P_element->node[2]->position, P_element->node[3]->position,
              &point_to_check[0]) == INSIDE) {
        element_for_polarization = P_element->number;
        break;
      }
    } // end of for loop over elements

    // if element is still -1, there was a problem
    if (element_for_polarization == -1) {
      if (node_location[0] > 58.12 && node_location[0] < 58.13 &&
          node_location[1] > 33.55 && node_location[1] < 33.56 &&
          node_location[2] > 61.90 && node_location[2] < 62.00) {
        d_print("Element not found: %d  %d  %d \n", site_1[0], site_1[1],
                site_1[2]);
      }
      return;
    }
  } // end of if site is a node block
  else {
    // get element site is in
    int l[3];
    for (int dof = 0; dof < 3; dof++)
      l[dof] = site_1[dof];

    double shape[4];
    const struct element_t *P_element =
        Element::getInstance()->locateSiteElement(&(shape[0]), &l[0], &lattice,
                                                  quasi);

    // set element for polarization
    element_for_polarization = P_element->number;
  }

  // insert locations into triangle
  std::vector<std::vector<double>> triangle_locations;
  triangle_locations.push_back(miscC->vectorAddition(location_1, face_shift));
  triangle_locations.push_back(miscC->vectorAddition(location_2, face_shift));
  triangle_locations.push_back(miscC->vectorAddition(location_3, face_shift));

  // calculate charge density
  double charge_density = 0.0;
  const std::vector<double> normal =
      outwardNormal(triangle_locations, point_to_check);
  for (int i_dof = 0; i_dof < 3; ++i_dof)
    charge_density +=
        element_polarization[element_for_polarization][i_dof] * normal[i_dof];

  // add to field
  const std::pair<std::vector<double>, double> return_field =
      fieldFromDensityFaceAdaptive(charge_density, triangle_locations,
                                   node_location, electric_constant,
                                   face_method, integration_error);

  // add to field
  for (int i_dof = 0; i_dof < 3; ++i_dof)
    electric_field[i_dof] += return_field.first[i_dof];

  // add to potential
  electric_potential += return_field.second;

  //
  return;
} // end of AddTriangleCluster

//
// compute contribution from single face of cluster
//
void ContributionFromSingleClusterSurface(
    const int &quasi, const std::vector<double> &node_location,
    const int &face_method, const double &integration_error,
    const double &electric_constant, const std::vector<int> &face_site,
    const int &const_index, const int &iterate_index_1,
    const int &iterate_index_2, const double &num_shells,
    const std::vector<std::vector<double>> &element_polarization,
    const std::vector<std::vector<std::vector<int>>> &atomistic_charges,
    const std::vector<double> &face_shift, const double &test_shift,
    struct lattice_t &lattice, const struct node_list_t &zero_node_list,
    const std::vector<std::vector<std::pair<std::vector<int>, int>>>
        &node_hash_list,
    std::vector<double> &electric_field, double &electric_potential) {
  // current site to get contribution from
  std::vector<int> current_site(3, 0);

  // additional sites for integration
  std::vector<int> site_2(3, 0);
  std::vector<int> site_3(3, 0);
  std::vector<int> site_4(3, 0);

  // set const index value
  current_site[const_index] = face_site[const_index];
  site_2[const_index] = face_site[const_index];
  site_3[const_index] = face_site[const_index];
  site_4[const_index] = face_site[const_index];

  // get Lattice Instance
  Lattice *latticeC = Lattice::getInstance();

  // iterate over 2 coordinates
  for (current_site[iterate_index_1] = face_site[iterate_index_1] - num_shells;
       current_site[iterate_index_1] <= face_site[iterate_index_1] + num_shells;
       ++current_site[iterate_index_1]) {
    for (current_site[iterate_index_2] =
             face_site[iterate_index_2] - num_shells;
         current_site[iterate_index_2] <=
         face_site[iterate_index_2] + num_shells;
         ++current_site[iterate_index_2]) {
      // flag for both triangles
      bool include_triangle_1 = true;
      bool include_triangle_2 = true;

      // set other sites
      site_2[iterate_index_1] = current_site[iterate_index_1] + 1;
      site_2[iterate_index_2] = current_site[iterate_index_2];
      site_3[iterate_index_1] = current_site[iterate_index_1];
      site_3[iterate_index_2] = current_site[iterate_index_2] + 1;
      site_4[iterate_index_1] = current_site[iterate_index_1] + 1;
      site_4[iterate_index_2] = current_site[iterate_index_2] + 1;

      // check if sites are in lattice
      if (latticeC->isSiteInsideLattice(&lattice, &current_site[0], quasi) ==
          RETURN_FAILURE)
        include_triangle_1 = false;
      if (latticeC->isSiteInsideLattice(&lattice, &site_2[0], quasi) ==
          RETURN_FAILURE) {
        include_triangle_1 = false;
        include_triangle_2 = false;
      }
      if (latticeC->isSiteInsideLattice(&lattice, &site_3[0], quasi) ==
          RETURN_FAILURE) {
        include_triangle_1 = false;
        include_triangle_2 = false;
      }
      if (latticeC->isSiteInsideLattice(&lattice, &site_4[0], quasi) ==
          RETURN_FAILURE)
        include_triangle_2 = false;

      // add contribution from triangle, potentially
      if (include_triangle_1 == true)
        AddTriangleCluster(quasi, node_location, current_site, site_2, site_3,
                           electric_constant, lattice, element_polarization,
                           atomistic_charges, integration_error, face_method,
                           test_shift, face_shift, const_index, zero_node_list,
                           node_hash_list, electric_field, electric_potential);
      if (include_triangle_2 == true)
        AddTriangleCluster(quasi, node_location, site_4, site_2, site_3,
                           electric_constant, lattice, element_polarization,
                           atomistic_charges, integration_error, face_method,
                           test_shift, face_shift, const_index, zero_node_list,
                           node_hash_list, electric_field, electric_potential);
    } // end of for loop over second iterated index
  }   // end of for loop over first iterated index

  //
  return;
} // end of ContributionFromSingleClusterSurface

//
// get shell of triangles
//
void AddShellTriangles(
    const std::vector<std::vector<int>> &coords, struct lattice_t &zero_lattice,
    const std::vector<double> &zero_shift, const int &const_index,
    const int &iterate_index_1, const int &iterate_index_2,
    const bool &negative_side,
    const std::vector<std::vector<std::vector<int>>> &atomistic_charges,
    std::vector<Triangle_3> &shell_triangles, const int &zero_quasi) {
  // current site to get contribution from
  std::vector<int> current_site(3, 0);

  // Get Lattice instance
  Lattice *latticeC = Lattice::getInstance();

  // additional sites for integration
  std::vector<int> site_2(3, 0);
  std::vector<int> site_3(3, 0);
  std::vector<int> site_4(3, 0);

  // set const index value
  if (negative_side == true) {
    current_site[const_index] = coords[const_index][0];
    site_2[const_index] = coords[const_index][0];
    site_3[const_index] = coords[const_index][0];
    site_4[const_index] = coords[const_index][0];
  } else {
    current_site[const_index] = coords[const_index][1];
    site_2[const_index] = coords[const_index][1];
    site_3[const_index] = coords[const_index][1];
    site_4[const_index] = coords[const_index][1];
  }

  // iterate over 2 coordinates
  for (current_site[iterate_index_1] = coords[iterate_index_1][0];
       current_site[iterate_index_1] <= coords[iterate_index_1][1];
       ++current_site[iterate_index_1]) {
    for (current_site[iterate_index_2] = coords[iterate_index_2][0];
         current_site[iterate_index_2] <= coords[iterate_index_2][1];
         ++current_site[iterate_index_2]) {
      // flag for both triangles
      bool include_triangle_1 = true;
      bool include_triangle_2 = true;

      // set other sites
      site_2[iterate_index_1] = current_site[iterate_index_1] + 1;
      site_2[iterate_index_2] = current_site[iterate_index_2];
      site_3[iterate_index_1] = current_site[iterate_index_1];
      site_3[iterate_index_2] = current_site[iterate_index_2] + 1;
      site_4[iterate_index_1] = current_site[iterate_index_1] + 1;
      site_4[iterate_index_2] = current_site[iterate_index_2] + 1;

      // check if sites are in lattice
      if (latticeC->isSiteInsideLattice(&zero_lattice, &current_site[0],
                                        zero_quasi) == RETURN_FAILURE)
        include_triangle_1 = false;
      if (latticeC->isSiteInsideLattice(&zero_lattice, &site_2[0],
                                        zero_quasi) == RETURN_FAILURE) {
        include_triangle_1 = false;
        include_triangle_2 = false;
      }
      if (latticeC->isSiteInsideLattice(&zero_lattice, &site_3[0],
                                        zero_quasi) == RETURN_FAILURE) {
        include_triangle_1 = false;
        include_triangle_2 = false;
      }
      if (latticeC->isSiteInsideLattice(&zero_lattice, &site_4[0],
                                        zero_quasi) == RETURN_FAILURE)
        include_triangle_2 = false;

      // add triangle 1
      if (include_triangle_1 == true) {
        // see if sites any sites are not atomistic
        if (isSiteAtomistic(current_site, atomistic_charges) == false ||
            isSiteAtomistic(site_2, atomistic_charges) == false ||
            isSiteAtomistic(site_3, atomistic_charges) == false) {
          // get site locations without shift
          std::vector<double> location_1(3, 0.0);
          std::vector<double> location_2(3, 0.0);
          std::vector<double> location_3(3, 0.0);
          latticeC->getSiteCurrentPosition(&location_1[0], &current_site[0],
                                           &zero_lattice, zero_quasi);

          latticeC->getSiteCurrentPosition(&location_2[0], &site_2[0],
                                           &zero_lattice, zero_quasi);

          latticeC->getSiteCurrentPosition(&location_3[0], &site_3[0],
                                           &zero_lattice, zero_quasi);

          // insert locations into triangle
          std::vector<Point_3> points;

          for (int i_dof = 0; i_dof < 3; ++i_dof) {
            location_1[i_dof] += zero_shift[i_dof];
            location_2[i_dof] += zero_shift[i_dof];
            location_3[i_dof] += zero_shift[i_dof];
          }

          points.push_back(
              Point_3(location_1[0], location_1[1], location_1[2]));
          points.push_back(
              Point_3(location_2[0], location_2[1], location_2[2]));
          points.push_back(
              Point_3(location_3[0], location_3[1], location_3[2]));
          shell_triangles.push_back(
              Triangle_3(points[0], points[1], points[2]));
        }
      }

      // add triangle 2
      if (include_triangle_2 == true) {
        // see if sites any sites are not atomistic
        if (isSiteAtomistic(site_4, atomistic_charges) == false ||
            isSiteAtomistic(site_2, atomistic_charges) == false ||
            isSiteAtomistic(site_3, atomistic_charges) == false) {
          // get site locations without shift
          std::vector<double> location_1(3, 0.0);
          std::vector<double> location_2(3, 0.0);
          std::vector<double> location_3(3, 0.0);

          latticeC->getSiteCurrentPosition(&location_1[0], &site_4[0],
                                           &zero_lattice, zero_quasi);

          latticeC->getSiteCurrentPosition(&location_2[0], &site_2[0],
                                           &zero_lattice, zero_quasi);

          latticeC->getSiteCurrentPosition(&location_3[0], &site_3[0],
                                           &zero_lattice, zero_quasi);

          // insert locations into triangle
          std::vector<Point_3> points;

          for (int i_dof = 0; i_dof < 3; ++i_dof) {
            location_1[i_dof] += zero_shift[i_dof];
            location_2[i_dof] += zero_shift[i_dof];
            location_3[i_dof] += zero_shift[i_dof];
          }

          points.push_back(
              Point_3(location_1[0], location_1[1], location_1[2]));
          points.push_back(
              Point_3(location_2[0], location_2[1], location_2[2]));
          points.push_back(
              Point_3(location_3[0], location_3[1], location_3[2]));
          shell_triangles.push_back(
              Triangle_3(points[0], points[1], points[2]));
        }
      }
    } // end of for loop over second iterated index
  }   // end of for loop over first iterated index

  //
  return;
} // end of AddShellTriangles

//
// get shell of triangles
//
void GetShellOfTriangles(
    const std::vector<std::vector<int>> &coords, struct lattice_t &zero_lattice,
    const std::vector<double> &zero_shift,
    const std::vector<std::vector<std::vector<int>>> &atomistic_charges,
    std::vector<Triangle_3> &shell_triangles, const int &zero_quasi) {
  // variables for adding faces
  int const_index;
  int iterate_index_1;
  int iterate_index_2;
  bool negative_side;

  // add triangles on -x face
  const_index = 0;
  iterate_index_1 = 1;
  iterate_index_2 = 2;
  negative_side = true;

  AddShellTriangles(coords, zero_lattice, zero_shift, const_index,
                    iterate_index_1, iterate_index_2, negative_side,
                    atomistic_charges, shell_triangles, zero_quasi);

  // add triangles on +x face
  const_index = 0;
  iterate_index_1 = 1;
  iterate_index_2 = 2;
  negative_side = false;

  AddShellTriangles(coords, zero_lattice, zero_shift, const_index,
                    iterate_index_1, iterate_index_2, negative_side,
                    atomistic_charges, shell_triangles, zero_quasi);

  // add triangles on -y face
  const_index = 1;
  iterate_index_1 = 0;
  iterate_index_2 = 2;
  negative_side = true;

  AddShellTriangles(coords, zero_lattice, zero_shift, const_index,
                    iterate_index_1, iterate_index_2, negative_side,
                    atomistic_charges, shell_triangles, zero_quasi);

  // add triangles on +y face
  const_index = 1;
  iterate_index_1 = 0;
  iterate_index_2 = 2;
  negative_side = false;

  AddShellTriangles(coords, zero_lattice, zero_shift, const_index,
                    iterate_index_1, iterate_index_2, negative_side,
                    atomistic_charges, shell_triangles, zero_quasi);

  // add triangles on -z face
  const_index = 2;
  iterate_index_1 = 0;
  iterate_index_2 = 1;
  negative_side = true;

  AddShellTriangles(coords, zero_lattice, zero_shift, const_index,
                    iterate_index_1, iterate_index_2, negative_side,
                    atomistic_charges, shell_triangles, zero_quasi);

  // add triangles on +z face
  const_index = 2;
  iterate_index_1 = 0;
  iterate_index_2 = 1;
  negative_side = false;

  AddShellTriangles(coords, zero_lattice, zero_shift, const_index,
                    iterate_index_1, iterate_index_2, negative_side,
                    atomistic_charges, shell_triangles, zero_quasi);

  //
  return;
} // end of GetShellOfTriangles

//
// determine if site is in cluster
//
bool NodeInCluster(const std::vector<int> &node_site, const int &num_shells,
                   const std::vector<int> &face_node_site) {
  // check coordinates
  for (int i_dof = 0; i_dof < 3; ++i_dof)
    if (face_node_site[i_dof] < node_site[i_dof] - num_shells &&
        face_node_site[i_dof] > node_site[i_dof] + num_shells + 1)
      return false;

  // if reached node is in cluster
  return true;
}

//
// break up intersection point polygon into triangles for integration
//
void GetTrianglesFromIntersectionPoints(
    std::vector<Point_3> &unique_intersection_points,
    std::vector<My_triangle> &integration_triangles) {
  // get normal
  std::vector<std::vector<double>> intersection_points_double;
  for (int i_point = 0; i_point < 3; ++i_point) {
    const My_point point = {
        CGAL::to_double(unique_intersection_points[i_point].x()),
        CGAL::to_double(unique_intersection_points[i_point].y()),
        CGAL::to_double(unique_intersection_points[i_point].z())};

    intersection_points_double.push_back(point);
  }

  // get test triangle
  My_triangle test_triangle;

  test_triangle.push_back(intersection_points_double[0]);
  test_triangle.push_back(intersection_points_double[1]);
  test_triangle.push_back(intersection_points_double[2]);

  // get normal
  std::vector<double> normal = CalculateNormal(test_triangle);

  for (int i_dof = 0; i_dof < 3; ++i_dof)
    if (normal[i_dof] < 0.0)
      normal[i_dof] *= -1.0;

  // ordered points container
  std::vector<My_point> ordered_points;

  // get convex hull on yz plane
  if (normal[0] >= normal[1] && normal[0] >= normal[2]) {
    // data
    std::vector<Point_2> ch_polygon;
    std::vector<Point_2> unique_intersection_points_2d;

    // get projection of points on yz plane
    for (int i_point = 0; i_point < unique_intersection_points.size();
         ++i_point)
      unique_intersection_points_2d.push_back(
          Point_2(unique_intersection_points[i_point].y(),
                  unique_intersection_points[i_point].z()));

    // make convex hull
    convex_hull_2(unique_intersection_points_2d.begin(),
                  unique_intersection_points_2d.end(),
                  std::back_inserter(ch_polygon));

    // loop over convex hull and insert points
    for (int i_point = 0; i_point < unique_intersection_points_2d.size();
         ++i_point) {
      // loop over all points
      for (int j_point = 0; j_point < intersection_points_double.size();
           ++j_point) {
        // check if point is within tolerance
        const double y_diff =
            CGAL::to_double(unique_intersection_points_2d[i_point][0]) -
            intersection_points_double[j_point][1];
        const double z_diff =
            CGAL::to_double(unique_intersection_points_2d[i_point][1]) -
            intersection_points_double[j_point][2];
        if (y_diff < 0.001 && y_diff > -0.001 && z_diff < 0.001 &&
            z_diff > -0.001) {
          // insert point into ordered list
          ordered_points.push_back(intersection_points_double[j_point]);

          // break to next convex hull point
          break;
        } // end of tolerance
        else if (j_point == unique_intersection_points_2d.size() - 1) {
          d_print("Error in convex hull algorithm\n");
          D_ERROR("GetTrianglesFromIntersectionPoints()");
          std::exit(EXIT_FAILURE);
        }
      } // end loop over j_points in intersection_points_double
    }   // end loop over i_points in unique_intersection_points_2d
  }     // end yz plane hull generation

  // get convex hull on xz plane
  else if (normal[1] >= normal[0] && normal[1] >= normal[2]) {
    // data
    std::vector<Point_2> ch_polygon;
    std::vector<Point_2> unique_intersection_points_2d;

    // get projection of points on yz plane
    for (int i_point = 0; i_point < unique_intersection_points.size();
         ++i_point)
      unique_intersection_points_2d.push_back(
          Point_2(unique_intersection_points[i_point].x(),
                  unique_intersection_points[i_point].z()));

    // make convex hull
    convex_hull_2(unique_intersection_points_2d.begin(),
                  unique_intersection_points_2d.end(),
                  std::back_inserter(ch_polygon));

    // loop over convex hull and insert points
    for (int i_point = 0; i_point < unique_intersection_points_2d.size();
         ++i_point) {
      // loop over all points
      for (int j_point = 0; j_point < intersection_points_double.size();
           ++j_point) {
        // check if point is within tolerance
        const double x_diff =
            CGAL::to_double(unique_intersection_points_2d[i_point][0]) -
            intersection_points_double[j_point][0];
        const double z_diff =
            CGAL::to_double(unique_intersection_points_2d[i_point][1]) -
            intersection_points_double[j_point][2];
        if (x_diff < 0.001 && x_diff > -0.001 && z_diff < 0.001 &&
            z_diff > -0.001) {
          // insert point into ordered list
          ordered_points.push_back(intersection_points_double[j_point]);

          // break to next convex hull point
          break;
        } // end of tolerance
        else if (j_point == unique_intersection_points_2d.size() - 1) {
          d_print("Error in convex hull algorithm\n");
          D_ERROR("GetTrianglesFromIntersectionPoints()");
          std::exit(EXIT_FAILURE);
        }
      } // end loop over j_points in intersection_points_double
    }   // end loop over i_points in unique_intersection_points_2d
  }     // end of xz convex hull

  // get convex hull on xy plane
  else {
    // data
    std::vector<Point_2> ch_polygon;
    std::vector<Point_2> unique_intersection_points_2d;

    // get projection of points on yz plane
    for (int i_point = 0; i_point < unique_intersection_points.size();
         ++i_point)
      unique_intersection_points_2d.push_back(
          Point_2(unique_intersection_points[i_point].x(),
                  unique_intersection_points[i_point].y()));

    // make convex hull
    convex_hull_2(unique_intersection_points_2d.begin(),
                  unique_intersection_points_2d.end(),
                  std::back_inserter(ch_polygon));

    // loop over convex hull and insert points
    for (int i_point = 0; i_point < unique_intersection_points_2d.size();
         ++i_point) {
      // loop over all points
      for (int j_point = 0; j_point < intersection_points_double.size();
           ++j_point) {
        // check if point is within tolerance
        const double x_diff =
            CGAL::to_double(unique_intersection_points_2d[i_point][0]) -
            intersection_points_double[j_point][0];
        const double y_diff =
            CGAL::to_double(unique_intersection_points_2d[i_point][1]) -
            intersection_points_double[j_point][1];
        if (x_diff < 0.001 && x_diff > -0.001 && y_diff < 0.001 &&
            y_diff > -0.001) {
          // insert point into ordered list
          ordered_points.push_back(intersection_points_double[j_point]);

          // break to next convex hull point
          break;
        } // end of tolerance
        else if (j_point == unique_intersection_points_2d.size() - 1) {
          d_print("Error in convex hull algorithm\n");
          D_ERROR("GetTrianglesFromIntersectionPoints()");
          std::exit(EXIT_FAILURE);
        }
      } // end loop over j_points in intersection_points_double
    }   // end loop over i_points in unique_intersection_points_2d
  }     // end of xy convex hull

  // make sure things went ok
  assert(ordered_points.size() == unique_intersection_points.size());

  // get average location
  My_point average_loc(3, 0.0);
  for (int i_point = 0; i_point < ordered_points.size(); ++i_point)
    for (int i_dof = 0; i_dof < 3; ++i_dof)
      average_loc[i_dof] += ordered_points[i_point][i_dof];

  for (int i_dof = 0; i_dof < 3; ++i_dof)
    average_loc[i_dof] /= static_cast<double>(ordered_points.size());

  // insert triangles
  for (int i_point = 0; i_point < ordered_points.size() - 1; ++i_point) {
    // get triangle
    const My_triangle temp_tri = {ordered_points[i_point],
                                  ordered_points[i_point + 1], average_loc};

    // insert triangle into return triangles
    integration_triangles.push_back(temp_tri);
  }

  // insert last triangle
  const My_triangle temp_tri = {ordered_points[0],
                                ordered_points[ordered_points.size() - 1],
                                average_loc};

  // insert triangle into return triangles
  integration_triangles.push_back(temp_tri);

  //
  return;
} // end of get_triangles_from_intersection_points

//
//  InterpolateValuesAtPoint()
//
void InterpolateValuesAtPoint(
    cluster_site_data_t iCSite_data, int iQuasi,
    const std::vector<std::vector<double>> electric_field,
    const std::vector<double> electric_potential, double charge,
    std::vector<double> &f_iC, double &energy_iC) {
  // matrix variable
  std::vector<std::vector<double>> matrix;
  matrix.resize(3);
  for (int i_dof = 0; i_dof < 3; ++i_dof)
    matrix[i_dof].resize(3, 0.0);

  f_iC.clear();
  f_iC.push_back(0.0);
  f_iC.push_back(0.0);
  f_iC.push_back(0.0);
  f_iC.push_back(0.0);
  energy_iC = 0.0;

  // set node numbers
  const int node_0 = iCSite_data.second[0];
  const int node_1 = iCSite_data.second[1];
  const int node_2 = iCSite_data.second[2];
  const int node_3 = iCSite_data.second[3];

  // get site state
  std::vector<double> site_state = iCSite_data.first.second;

  //
  //  get Node list for iQuasi
  //
  struct node_list_t node_list =
      Quasicontinua::getInstance()->getQuasi(iQuasi).getNodeList().node_list;

  // get node locations
  std::vector<double> node_0_location;
  std::vector<double> node_1_location;
  std::vector<double> node_2_location;
  std::vector<double> node_3_location;

  for (int i_dof = 0; i_dof < 3; ++i_dof) {
    node_0_location.push_back(node_list.nodes[node_0]->position[i_dof]);
    node_1_location.push_back(node_list.nodes[node_1]->position[i_dof]);
    node_2_location.push_back(node_list.nodes[node_2]->position[i_dof]);
    node_3_location.push_back(node_list.nodes[node_3]->position[i_dof]);
  }

  // get handle of MiscFunctions class
  MiscFunctions *miscC = MiscFunctions::getInstance();

  // get node_difference_1, node_difference_2, node_difference_3
  const std::vector<double> node_difference_1 =
      miscC->vectorDifference(node_1_location, node_0_location);
  const std::vector<double> node_difference_2 =
      miscC->vectorDifference(node_2_location, node_0_location);
  const std::vector<double> node_difference_3 =
      miscC->vectorDifference(node_3_location, node_0_location);

  // get J
  const double J = miscC->dotProduct(
      miscC->crossProduct3x3(node_difference_1, node_difference_2),
      node_difference_3);

  // variables for interpolation
  double a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4;

  // get a1
  matrix[0][0] = node_1_location[0];
  matrix[0][1] = node_1_location[1];
  matrix[0][2] = node_1_location[2];
  matrix[1][0] = node_2_location[0];
  matrix[1][1] = node_2_location[1];
  matrix[1][2] = node_2_location[2];
  matrix[2][0] = node_3_location[0];
  matrix[2][1] = node_3_location[1];
  matrix[2][2] = node_3_location[2];
  a1 = miscC->determinant3x3(matrix) / J;

  // get a2
  matrix[0][0] = node_0_location[0];
  matrix[0][1] = node_0_location[1];
  matrix[0][2] = node_0_location[2];
  matrix[1][0] = node_2_location[0];
  matrix[1][1] = node_2_location[1];
  matrix[1][2] = node_2_location[2];
  matrix[2][0] = node_3_location[0];
  matrix[2][1] = node_3_location[1];
  matrix[2][2] = node_3_location[2];
  a2 = -1.0 * miscC->determinant3x3(matrix) / J;

  // get a3
  matrix[0][0] = node_0_location[0];
  matrix[0][1] = node_0_location[1];
  matrix[0][2] = node_0_location[2];
  matrix[1][0] = node_1_location[0];
  matrix[1][1] = node_1_location[1];
  matrix[1][2] = node_1_location[2];
  matrix[2][0] = node_3_location[0];
  matrix[2][1] = node_3_location[1];
  matrix[2][2] = node_3_location[2];
  a3 = miscC->determinant3x3(matrix) / J;

  // get a4
  matrix[0][0] = node_0_location[0];
  matrix[0][1] = node_0_location[1];
  matrix[0][2] = node_0_location[2];
  matrix[1][0] = node_1_location[0];
  matrix[1][1] = node_1_location[1];
  matrix[1][2] = node_1_location[2];
  matrix[2][0] = node_2_location[0];
  matrix[2][1] = node_2_location[1];
  matrix[2][2] = node_2_location[2];
  a4 = -1.0 * miscC->determinant3x3(matrix) / J;

  // get b1
  matrix[0][0] = 1.0;
  matrix[0][1] = node_1_location[1];
  matrix[0][2] = node_1_location[2];
  matrix[1][0] = 1.0;
  matrix[1][1] = node_2_location[1];
  matrix[1][2] = node_2_location[2];
  matrix[2][0] = 1.0;
  matrix[2][1] = node_3_location[1];
  matrix[2][2] = node_3_location[2];
  b1 = -1.0 * miscC->determinant3x3(matrix) / J;

  // get b2
  matrix[0][0] = 1.0;
  matrix[0][1] = node_0_location[1];
  matrix[0][2] = node_0_location[2];
  matrix[1][0] = 1.0;
  matrix[1][1] = node_2_location[1];
  matrix[1][2] = node_2_location[2];
  matrix[2][0] = 1.0;
  matrix[2][1] = node_3_location[1];
  matrix[2][2] = node_3_location[2];
  b2 = miscC->determinant3x3(matrix) / J;

  // get b3
  matrix[0][0] = 1.0;
  matrix[0][1] = node_0_location[1];
  matrix[0][2] = node_0_location[2];
  matrix[1][0] = 1.0;
  matrix[1][1] = node_1_location[1];
  matrix[1][2] = node_1_location[2];
  matrix[2][0] = 1.0;
  matrix[2][1] = node_3_location[1];
  matrix[2][2] = node_3_location[2];
  b3 = -1.0 * miscC->determinant3x3(matrix) / J;

  // get b4
  matrix[0][0] = 1.0;
  matrix[0][1] = node_0_location[1];
  matrix[0][2] = node_0_location[2];
  matrix[1][0] = 1.0;
  matrix[1][1] = node_1_location[1];
  matrix[1][2] = node_1_location[2];
  matrix[2][0] = 1.0;
  matrix[2][1] = node_2_location[1];
  matrix[2][2] = node_2_location[2];
  b4 = miscC->determinant3x3(matrix) / J;

  // get c1
  matrix[0][0] = 1.0;
  matrix[0][1] = node_1_location[0];
  matrix[0][2] = node_1_location[2];
  matrix[1][0] = 1.0;
  matrix[1][1] = node_2_location[0];
  matrix[1][2] = node_2_location[2];
  matrix[2][0] = 1.0;
  matrix[2][1] = node_3_location[0];
  matrix[2][2] = node_3_location[2];
  c1 = miscC->determinant3x3(matrix) / J;

  // get c2
  matrix[0][0] = 1.0;
  matrix[0][1] = node_0_location[0];
  matrix[0][2] = node_0_location[2];
  matrix[1][0] = 1.0;
  matrix[1][1] = node_2_location[0];
  matrix[1][2] = node_2_location[2];
  matrix[2][0] = 1.0;
  matrix[2][1] = node_3_location[0];
  matrix[2][2] = node_3_location[2];
  c2 = -1.0 * miscC->determinant3x3(matrix) / J;

  // get c3
  matrix[0][0] = 1.0;
  matrix[0][1] = node_0_location[0];
  matrix[0][2] = node_0_location[2];
  matrix[1][0] = 1.0;
  matrix[1][1] = node_1_location[0];
  matrix[1][2] = node_1_location[2];
  matrix[2][0] = 1.0;
  matrix[2][1] = node_3_location[0];
  matrix[2][2] = node_3_location[2];
  c3 = miscC->determinant3x3(matrix) / J;

  // get c4
  matrix[0][0] = 1.0;
  matrix[0][1] = node_0_location[0];
  matrix[0][2] = node_0_location[2];
  matrix[1][0] = 1.0;
  matrix[1][1] = node_1_location[0];
  matrix[1][2] = node_1_location[2];
  matrix[2][0] = 1.0;
  matrix[2][1] = node_2_location[0];
  matrix[2][2] = node_2_location[2];
  c4 = -1.0 * miscC->determinant3x3(matrix) / J;

  // get d1
  matrix[0][0] = 1.0;
  matrix[0][1] = node_1_location[0];
  matrix[0][2] = node_1_location[1];
  matrix[1][0] = 1.0;
  matrix[1][1] = node_2_location[0];
  matrix[1][2] = node_2_location[1];
  matrix[2][0] = 1.0;
  matrix[2][1] = node_3_location[0];
  matrix[2][2] = node_3_location[1];
  d1 = -1.0 * miscC->determinant3x3(matrix) / J;

  // get d2
  matrix[0][0] = 1.0;
  matrix[0][1] = node_0_location[0];
  matrix[0][2] = node_0_location[1];
  matrix[1][0] = 1.0;
  matrix[1][1] = node_2_location[0];
  matrix[1][2] = node_2_location[1];
  matrix[2][0] = 1.0;
  matrix[2][1] = node_3_location[0];
  matrix[2][2] = node_3_location[1];
  d2 = miscC->determinant3x3(matrix) / J;

  // get d3
  matrix[0][0] = 1.0;
  matrix[0][1] = node_0_location[0];
  matrix[0][2] = node_0_location[1];
  matrix[1][0] = 1.0;
  matrix[1][1] = node_1_location[0];
  matrix[1][2] = node_1_location[1];
  matrix[2][0] = 1.0;
  matrix[2][1] = node_3_location[0];
  matrix[2][2] = node_3_location[1];
  d3 = -1.0 * miscC->determinant3x3(matrix) / J;

  // get d4
  matrix[0][0] = 1.0;
  matrix[0][1] = node_0_location[0];
  matrix[0][2] = node_0_location[1];
  matrix[1][0] = 1.0;
  matrix[1][1] = node_1_location[0];
  matrix[1][2] = node_1_location[1];
  matrix[2][0] = 1.0;
  matrix[2][1] = node_2_location[0];
  matrix[2][2] = node_2_location[1];
  d4 = miscC->determinant3x3(matrix) / J;

  // get force
  for (int i_dof = 0; i_dof < 4; ++i_dof)
    f_iC[i_dof] =
        (a1 + b1 * site_state[0] + c1 * site_state[1] + d1 * site_state[2]) *
            electric_field[node_0][i_dof] +
        (a2 + b2 * site_state[0] + c2 * site_state[1] + d2 * site_state[2]) *
            electric_field[node_1][i_dof] +
        (a3 + b3 * site_state[0] + c3 * site_state[1] + d3 * site_state[2]) *
            electric_field[node_2][i_dof] +
        (a4 + b4 * site_state[0] + c4 * site_state[1] + d4 * site_state[2]) *
            electric_field[node_3][i_dof];

  for (int i_dof = 0; i_dof < 4; ++i_dof)
    f_iC[i_dof] *= charge;

  // get energy
  energy_iC =
      0.5 * charge *
      ((a1 + b1 * site_state[0] + c1 * site_state[1] + d1 * site_state[2]) *
           electric_potential[node_0] +
       (a2 + b2 * site_state[0] + c2 * site_state[1] + d2 * site_state[2]) *
           electric_potential[node_1] +
       (a3 + b3 * site_state[0] + c3 * site_state[1] + d3 * site_state[2]) *
           electric_potential[node_2] +
       (a4 + b4 * site_state[0] + c4 * site_state[1] + d4 * site_state[2]) *
           electric_potential[node_3]);

  //
  return;
} // end of InterpolateValuesAtPoint()
} // end of namespace for local functions

//
//  namespace for computeGeometries related functions
//
namespace {
//
//  computeFacesWorker
//
void *computeFacesWorker(void *arg) {
  //
  // cast back to get computeFaces_argument_t
  //
  struct computeFaces_argument_t *computeFaces_argument =
      static_cast<struct computeFaces_argument_t *>(arg);

  //
  // get handles
  //
  std::vector<std::vector<std::vector<int>>> &faces =
      *(computeFaces_argument->faces);
  std::vector<std::vector<std::pair<int, int>>> &faceElementList =
      *(computeFaces_argument->faceElementList);
  struct element_list_t &elementList = *(computeFaces_argument->elementList);

  //
  // get share of data
  //
  const int number_threads = quasicontinuum::get_number_threads();
  const int my_id = quasicontinuum::get_my_tid();
  int number_data_thread;
  int number_start;
  int number_end;

  //
  // get share of work
  //
  quasicontinuum::get_share(my_id, number_threads, elementList.number_elements,
                            &number_data_thread, &number_start, &number_end);

  //
  // loop over elements
  //
  for (int iElem = number_start; iElem <= number_end; ++iElem) {
    //
    // get node numbers
    //
    std::vector<int> nodes;
    nodes.push_back(elementList.elements[iElem]->node[0]->number);
    nodes.push_back(elementList.elements[iElem]->node[1]->number);
    nodes.push_back(elementList.elements[iElem]->node[2]->number);
    nodes.push_back(elementList.elements[iElem]->node[3]->number);
    std::sort(nodes.begin(), nodes.begin() + 4);

    //
    // set faces
    //
    std::vector<std::vector<int>> elementFaces;
    std::vector<int> face(4, 0);

    //
    // insert node numbers into face
    //
    face[0] = nodes[0];
    face[1] = nodes[1];
    face[2] = nodes[2];
    face[3] = nodes[3];

    //
    // insert 1st face
    //
    elementFaces.push_back(face);

    //
    // insert node numbers into face
    //
    face[0] = nodes[0];
    face[1] = nodes[1];
    face[2] = nodes[3];
    face[3] = nodes[2];

    //
    // insert 2nd face
    //
    elementFaces.push_back(face);

    //
    // insert node numbers into face
    //
    face[0] = nodes[1];
    face[1] = nodes[2];
    face[2] = nodes[3];
    face[3] = nodes[0];

    //
    // insert 3rd face
    //
    elementFaces.push_back(face);

    //
    // insert node numbers into face
    //
    face[0] = nodes[0];
    face[1] = nodes[2];
    face[2] = nodes[3];
    face[3] = nodes[1];

    //
    // insert 4th face
    //
    elementFaces.push_back(face);

    //
    // loop over all 4 faces
    //
    for (int iFace = 0; iFace < 4; ++iFace) {
//
// for debugging
//
#if 0

            faceNum++;
            const int node1 = elementFaces[iFace][0];
            const int node2 = elementFaces[iFace][1];
            const int node3 = elementFaces[iFace][2];

            std::cout<<faceNum<<" Node 1 "
               <<quasiNodeList.nodes[node1]->l[0]<<" "
               <<quasiNodeList.nodes[node1]->l[1]<<" "
               <<quasiNodeList.nodes[node1]->l[2]<<" Node 2 "
               <<quasiNodeList.nodes[node2]->l[0]<<" "
               <<quasiNodeList.nodes[node2]->l[1]<<" "
               <<quasiNodeList.nodes[node2]->l[2]<<" Node 3 "
               <<quasiNodeList.nodes[node3]->l[0]<<" "
               <<quasiNodeList.nodes[node3]->l[1]<<" "
               <<quasiNodeList.nodes[node3]->l[2]<<std::endl;
#endif

      //
      // get bucket
      //
      unsigned int bucket =
          getTripletIntegersKey(elementFaces[iFace][0], elementFaces[iFace][1],
                                elementFaces[iFace][2]);

      //
      // lock bucket
      //
      pthread_mutex_lock(&bucketLocks[bucket]);

      //
      // iterators
      //
      std::vector<std::vector<int>>::iterator faceIt;
      std::vector<std::pair<int, int>>::iterator listIt;

      //
      // search bucket
      //
      for (unsigned int it = 0; it != faces[bucket].size() + 1; ++it) {
        //
        // check if iterator is at end
        //
        if (it == faces[bucket].size()) {
          //
          // insert face into global faces structure
          //
          faceIt = faces[bucket].end();
          faces[bucket].insert(faceIt, elementFaces[iFace]);

          //
          // insert iElem into face element list
          //
          listIt = faceElementList[bucket].end();
          faceElementList[bucket].insert(listIt, std::make_pair(iElem, -1));

//
//  for debugging purpose
//
#if 0
                faceNum++;
                const int node1 = elementFaces[iFace][0];
                const int node2 = elementFaces[iFace][1];
                const int node3 = elementFaces[iFace][2];
                std::cout<<faceNum<<" Node 1 "
                   <<quasiNodeList.nodes[node1]->l[0]<<" "
                   <<quasiNodeList.nodes[node1]->l[1]<<" "
                   <<quasiNodeList.nodes[node1]->l[2]<<" Node 2 "
                   <<quasiNodeList.nodes[node2]->l[0]<<" "
                   <<quasiNodeList.nodes[node2]->l[1]<<" "
                   <<quasiNodeList.nodes[node2]->l[2]<<" Node 3 "
                   <<quasiNodeList.nodes[node3]->l[0]<<" "
                   <<quasiNodeList.nodes[node3]->l[1]<<" "
                   <<quasiNodeList.nodes[node3]->l[2]<<std::endl;
#endif

          //
          // break out of for loop
          //
          break;
        }

        //
        // check if face matches and add element to face element list
        //
        if (faces[bucket][it][0] == elementFaces[iFace][0] &&
            faces[bucket][it][1] == elementFaces[iFace][1] &&
            faces[bucket][it][2] == elementFaces[iFace][2]) {
          //
          // add additional 4th node to list
          //
          faces[bucket][it].push_back(elementFaces[iFace][3]);

          //
          // add second element to list
          //
          faceElementList[bucket][it].second = iElem;

          //
          // break out of for loop
          //
          break;
        }

        //
        // check if have passed point of where face would be
        //
        if (faces[bucket][it][0] > elementFaces[iFace][0]) {
          //
          // insert face into global faces structure
          //
          faceIt = faces[bucket].begin() + it;
          faces[bucket].insert(faceIt, elementFaces[iFace]);

          //
          // insert iElem into face element list
          //
          listIt = faceElementList[bucket].begin() + it;
          faceElementList[bucket].insert(listIt, std::make_pair(iElem, -1));

          //
          // break out of for loop
          //
          break;
        }
      }

      //
      // unlock bucket
      //
      pthread_mutex_unlock(&bucketLocks[bucket]);
    }
  }

  //
  //
  //
  return NULL;
} // end of computeFacesWorker()

//
//  computePolarizationSitesWorker()
//
void *computePolarizationSitesWorker(void *arg) {
  //
  // cast back to get computePolarizationSites_argument_t
  //
  struct computePolarizationSites_argument_t
      *computePolarizationSites_argument =
          static_cast<struct computePolarizationSites_argument_t *>(arg);

  //
  // get handles
  //
  std::vector<std::vector<int>> &sites =
      *(computePolarizationSites_argument->sites);
  struct element_list_t &elementList =
      *(computePolarizationSites_argument->elementList);
  struct lattice_t &lattice = *(computePolarizationSites_argument->lattice);
  int &quasi = *(computePolarizationSites_argument->quasi);

  //
  // get share of data
  //
  const int number_threads = quasicontinuum::get_number_threads();
  const int my_id = quasicontinuum::get_my_tid();
  int number_data_thread;
  int number_start;
  int number_end;

  //
  // get share of work
  //
  quasicontinuum::get_share(my_id, number_threads, elementList.number_elements,
                            &number_data_thread, &number_start, &number_end);

  //
  // loop over elements
  //
  for (int iElem = number_start; iElem <= number_end; ++iElem) {
    //
    // array to hold center of tetrahedron
    //
    double center[3];
    center[0] = 0.0;
    center[1] = 0.0;
    center[2] = 0.0;

    //
    // loop over all directions and add initial positions
    //
    for (int iDof = 0; iDof < 3; ++iDof) {
      //
      // loop over all nodes
      //
      for (int iNode = 0; iNode < 4; ++iNode) {
        //
        // add initial positions
        //
        center[iDof] +=
            elementList.elements[iElem]->node[iNode]->initial_position[iDof];
      }

      //
      // divide center by 4 to get position
      //
      center[iDof] = center[iDof] / 4.0;
    }

    //
    // array to hold site
    //
    int site[3];

    //
    // find closest site
    //
    Lattice::getInstance()->findClosestSite(site, center, &lattice, quasi);

    //
    // resize to numQuasi
    //
    sites[iElem].clear();

    //
    // insert site into element list
    //
    sites[iElem].push_back(site[0]);
    sites[iElem].push_back(site[1]);
    sites[iElem].push_back(site[2]);
  }

  //
  //
  //
  return NULL;
} // end of computePolarizationSitesWorker()

//
//  computeAtomisticElementsWorker()
//
void *computeAtomisticElementsWorker(void *arg) {
  //
  // cast back to get computeAtomisticElements_argument_t
  //
  struct computeAtomisticElements_argument_t
      *computeAtomisticElements_argument =
          static_cast<struct computeAtomisticElements_argument_t *>(arg);

  //
  // get handles
  //
  std::vector<int> &atomisticElements =
      *(computeAtomisticElements_argument->atomisticElements);
  struct element_list_t &elementList =
      *(computeAtomisticElements_argument->elementList);
  double &atomisticElementSize =
      *(computeAtomisticElements_argument->atomisticElementSize);

  //
  // get share of data
  //
  const int number_threads = quasicontinuum::get_number_threads();
  const int my_id = quasicontinuum::get_my_tid();
  int number_data_thread;
  int number_start;
  int number_end;

  //
  // get share of work
  //
  quasicontinuum::get_share(my_id, number_threads, elementList.number_elements,
                            &number_data_thread, &number_start, &number_end);

  //
  // loop over elements
  //
  for (int iElem = number_start; iElem <= number_end; ++iElem) {
    //
    // get volume of element
    //
    double elementVolume = tetrahedronVolume(
        elementList.elements[iElem]->node[0]->initial_position,
        elementList.elements[iElem]->node[1]->initial_position,
        elementList.elements[iElem]->node[2]->initial_position,
        elementList.elements[iElem]->node[3]->initial_position);

    //
    // check if size of element is less then atomistic element size
    //
    if (elementVolume <= atomisticElementSize)
      atomisticElements[iElem] = 0;
    else
      atomisticElements[iElem] = 1;
  }

  //
  //
  //
  return NULL;
} // end of computeAtomisticElementsWorker()

//
//  computeAtomisticElementLocationsWorker()
//
void *computeAtomisticElementLocationsWorker(void *arg) {
  //
  // cast back to get computeAtomisticElementLocations_argument_t
  //
  struct computeAtomisticElementLocations_argument_t
      *computeAtomisticElementLocations_argument =
          static_cast<struct computeAtomisticElementLocations_argument_t *>(
              arg);

  //
  // get handles
  //
  std::vector<int> &atomisticElements =
      *(computeAtomisticElementLocations_argument->atomisticElements);
  std::vector<std::vector<std::vector<double>>> &atomisticElementLocations =
      *(computeAtomisticElementLocations_argument->atomisticElementLocations);
  struct element_list_t &elementList =
      *(computeAtomisticElementLocations_argument->elementList);

  //
  // variables for elements
  //
  int number_start;
  int number_end;

  //
  // keep getting more elements
  //
  while (whileWorker) {
    //
    // get element lock
    //
    pthread_mutex_lock(&dataLock);

    //
    // set variables
    //
    if (currentElem == int(atomisticElements.size())) {
      pthread_mutex_unlock(&dataLock);
      break;
    } else if (currentElem + elemBucket >= int(atomisticElements.size())) {
      //
      // set variables
      //
      number_start = currentElem;
      number_end = atomisticElements.size();
      currentElem = number_end;
    } else {
      //
      // set variables
      //
      number_start = currentElem;
      number_end = number_start + elemBucket;
      currentElem = number_end;
    }

    //
    // get element lock
    //
    pthread_mutex_unlock(&dataLock);

    //
    // loop over elements
    //
    for (int iElem = number_start; iElem < number_end; ++iElem) {
      //
      // check for atomistic elements
      //
      if (atomisticElements[iElem] == 0) {
        //
        // vector of position
        //
        std::vector<double> nodePosition(3, 0.0);

        //
        // clear node positions
        //
        atomisticElementLocations[iElem].clear();

        //
        // loop over all nodes
        //
        for (int iNode = 0; iNode < 4; ++iNode) {
          //
          // populate node
          //
          nodePosition[0] =
              elementList.elements[iElem]->node[iNode]->initial_position[0];
          nodePosition[1] =
              elementList.elements[iElem]->node[iNode]->initial_position[1];
          nodePosition[2] =
              elementList.elements[iElem]->node[iNode]->initial_position[2];

          //
          // insert position into data structure
          //
          atomisticElementLocations[iElem].push_back(nodePosition);
        }
      }
    }
  }

  //
  //
  //
  return NULL;
} // end of computeAtomisticElementLocationsWorker()

//
// computeAtomisticNodesWorker()
//
void *computeAtomisticNodesWorker(void *arg) {
  //
  // cast back to get computeAtomisticNodes_argument_t
  //
  struct computeAtomisticNodes_argument_t *computeAtomisticNodes_argument =
      static_cast<struct computeAtomisticNodes_argument_t *>(arg);

  //
  // get handles
  //
  std::vector<int> &atomisticElements =
      *(computeAtomisticNodes_argument->atomisticElements);
  struct node_list_t &nodeList = *(computeAtomisticNodes_argument->nodeList);
  std::vector<std::pair<int, std::vector<int>>> &atomisticNodes =
      *(computeAtomisticNodes_argument->atomisticNodes);

  //
  // get share of data
  //
  const int number_threads = quasicontinuum::get_number_threads();
  const int my_id = quasicontinuum::get_my_tid();
  int number_data_thread;
  int number_start;
  int number_end;

  //
  // get share of work
  //
  quasicontinuum::get_share(my_id, number_threads, nodeList.number_nodes,
                            &number_data_thread, &number_start, &number_end);

  //
  // local atomistic nodes
  //
  std::vector<std::pair<int, std::vector<int>>> localAtomisticNodes;

  //
  // loop over nodes
  //
  for (int iNode = number_start; iNode <= number_end; ++iNode) {
    //
    // get number of elements node is in
    //
    int numElem = nodeList.nodes[iNode]->element_list.number_elements;

    //
    // vector of element numbers
    //
    std::vector<int> elemNumbers;

    //
    // loop over number of elements and get each element number
    //
    for (int iElem = 0; iElem < numElem; ++iElem)
      elemNumbers.push_back(
          nodeList.nodes[iNode]->element_list.elements[iElem]->number);

    //
    // variable to hold node sum
    //
    int nodeSum = 0;

    //
    // determine if node should be included in atomistic
    //
    for (int iElem = 0; iElem < numElem; ++iElem)
      nodeSum += atomisticElements[elemNumbers[iElem]];

    if (nodeSum == 0) {
      //
      // vector to hold site location
      //
      std::vector<int> site(3, 0);

      //
      // populate site
      //
      site[0] = nodeList.nodes[iNode]->l[0];
      site[1] = nodeList.nodes[iNode]->l[1];
      site[2] = nodeList.nodes[iNode]->l[2];

      //
      // make pair
      //
      localAtomisticNodes.push_back(std::make_pair(iNode, site));
    }
  }

  //
  // keep looping until thread turn to populate global list
  //
  while (whileWorker) {
    //
    // lock for sort variable
    //
    pthread_mutex_lock(&dataLock);

    //
    // check for turn
    //
    if (sortVariable == my_id) {
      //
      // insert nodes on end
      //
      atomisticNodes.insert(atomisticNodes.end(), localAtomisticNodes.begin(),
                            localAtomisticNodes.end());

      //
      // update sort variable
      //
      sortVariable++;

      //
      // unlock
      //
      pthread_mutex_unlock(&dataLock);

      //
      // break while loop
      //
      break;
    }

    //
    // unlock
    //
    pthread_mutex_unlock(&dataLock);
  }

  //
  //
  //
  return NULL;
} // end of computeAtomisticNodesWorker()

//
//  computeAtomisticChargesWorker()
//
void *computeAtomisticChargesWorker(void *arg) {
  //
  // cast back to get computeAtomisticNodes_argument_t
  //
  struct computeAtomisticCharges_argument_t *computeAtomisticCharges_argument =
      static_cast<struct computeAtomisticCharges_argument_t *>(arg);

  //
  // get handles
  //
  std::vector<int> &atomisticElements =
      *(computeAtomisticCharges_argument->atomisticElements);
  std::vector<std::vector<std::vector<double>>> &atomisticElementLocations =
      *(computeAtomisticCharges_argument->atomisticElementLocations);
  std::vector<std::vector<std::vector<int>>> &atomisticCharges =
      *(computeAtomisticCharges_argument->atomisticCharges);
  struct lattice_t &lattice = *(computeAtomisticCharges_argument->lattice);
  struct element_list_t &elementList =
      *(computeAtomisticCharges_argument->elementList);
  int quasi = *(computeAtomisticCharges_argument->quasi);

  //
  // variables for sites
  //
  int number_start;
  int number_end;

  //
  // local atomistic charges
  //
  std::vector<std::vector<int>> localSites;

  //
  // keep getting more elements
  //
  while (whileWorker) {
    //
    // get element lock
    //
    pthread_mutex_lock(&dataLock);

    //
    // set variables
    //
    if (currentElem == int(atomisticElements.size())) {
      pthread_mutex_unlock(&dataLock);
      break;
    } else if (currentElem + elemBucket >= int(atomisticElements.size())) {
      //
      // set variables
      //
      number_start = currentElem;
      number_end = atomisticElements.size();
      currentElem = number_end;
    } else {
      //
      // set variables
      //
      number_start = currentElem;
      number_end = number_start + elemBucket;
      currentElem = number_end;
    }

    //
    // get element lock
    //
    pthread_mutex_unlock(&dataLock);

    //
    // loop over atomistic elements
    //
    for (int iElem = number_start; iElem < number_end; ++iElem) {
      //
      // check for atomistic elements
      //
      if (atomisticElements[iElem] == 0) {
        //
        // loop over all nodes of element
        //
        for (int iNode = 0; iNode < 4; ++iNode) {
          //
          // get site
          //
          std::vector<int> elementSite;
          elementSite.push_back(elementList.elements[iElem]->node[iNode]->l[0]);
          elementSite.push_back(elementList.elements[iElem]->node[iNode]->l[1]);
          elementSite.push_back(elementList.elements[iElem]->node[iNode]->l[2]);

          //
          // get sites inside element
          //
          std::vector<std::vector<int>> addSites =
              sitesInElement(elementSite, iElem, atomisticElementLocations,
                             atomisticElements, lattice, quasi);

          //
          // loop over all sites
          //
          for (int iSite = 0; iSite < int(addSites.size()); ++iSite) {

            //
            // get bucket
            //
            int bucket = getTripletIntegersKey(
                addSites[iSite][0], addSites[iSite][1], addSites[iSite][2]);

            //
            // get bucket lock
            //
            pthread_mutex_lock(&bucketLocks[bucket]);

            //
            // insert into bucket
            //
            insertChargeIntoBucket(addSites[iSite], atomisticCharges[bucket]);

            //
            // unlock bucket lock
            //
            pthread_mutex_unlock(&bucketLocks[bucket]);
          }
        }
      }
    }
  }

  //
  //
  //
  return NULL;
} // end of computeAtomisticChargesWorker()
} // end of namespace for computeGeomeries related functions

//
//  namespace for computeElectricField()
//
namespace {
//
// compute polarization locations
//
void *computePolarizationLocationsWorker(void *arg) {
  //
  // cast back to get computeAtomisticNodes_argument_t
  //
  struct computePolarizationLocations_argument_t
      *computePolarizationLocations_argument =
          static_cast<struct computePolarizationLocations_argument_t *>(arg);

  //
  // get handles
  //
  std::vector<int> &atomisticElements =
      *(computePolarizationLocations_argument->atomisticElements);
  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      &elementPolarizationLocations = *(
          computePolarizationLocations_argument->elementPolarizationLocations);
  std::vector<std::vector<int>> &elementPolarizationSites =
      *(computePolarizationLocations_argument->elementPolarizationSites);
  int iQuasi = computePolarizationLocations_argument->quasiId;
  double charge = computePolarizationLocations_argument->fixedCharge;
  struct lattice_t &lattice = *(computePolarizationLocations_argument->lattice);

  //
  // variables for sites
  //
  int number_start;
  int number_end;

  //
  // get shift vector
  //
  std::vector<double> shift = Quasicontinua::getInstance()->getShift(iQuasi);

  //
  // keep getting more elements
  //
  while (whileWorker) {
    //
    // get element lock
    //
    pthread_mutex_lock(&dataLock);

    //
    // set variables
    //
    if (currentElem == int(atomisticElements.size())) {
      pthread_mutex_unlock(&dataLock);
      break;
    } else if (currentElem + elemBucket >= int(atomisticElements.size())) {
      //
      // set variables
      //
      number_start = currentElem;
      number_end = atomisticElements.size();
      currentElem = number_end;
    } else {
      //
      // set variables
      //
      number_start = currentElem;
      number_end = number_start + elemBucket;
      currentElem = number_end;
    }

    //
    // get element lock
    //
    pthread_mutex_unlock(&dataLock);

    //
    // temp vector to hold location
    //
    std::vector<double> tempLocation(3, 0.0);

    //
    // loop over elements
    //
    for (int iElem = number_start; iElem < number_end; ++iElem) {
      //
      // clear list if first quasi
      //
      if (iQuasi == 0)
        elementPolarizationLocations[iElem].clear();

      //
      // check for non atomistic elements
      //
      if (atomisticElements[iElem] == 1) {
        //
        // temporary location
        //
        double location[3];

        //
        // temporary site
        //
        int site[3];
        site[0] = elementPolarizationSites[iElem][0];
        site[1] = elementPolarizationSites[iElem][1];
        site[2] = elementPolarizationSites[iElem][2];

        if (Lattice::getInstance()->isSiteInsideLattice(
                &lattice, &site[0], iQuasi) == RETURN_FAILURE) {
          d_print("site = (%d, %d, %d) is not inside lattice\n", site[0],
                  site[1], site[2]);
        }

        //
        // get site location
        //

        Lattice::getInstance()->getSiteCurrentPosition(location, site, &lattice,
                                                       iQuasi);

        //
        // add shift to temp vector
        //
        tempLocation[0] = location[0] + shift[0];
        tempLocation[1] = location[1] + shift[1];
        tempLocation[2] = location[2] + shift[2];

        //
        // add actual location and charge to data
        //
        elementPolarizationLocations[iElem].push_back(
            std::make_pair(tempLocation, charge));
      }
    }
  }

  //
  //
  //
  return NULL;
}

//
// compute element polarizations
//
void *computeElementPolarizationWorker(void *arg) {
  //
  // cast back to get computeAtomisticNodes_argument_t
  //
  struct computeElementPolarization_argument_t
      *computeElementPolarization_argument =
          static_cast<struct computeElementPolarization_argument_t *>(arg);

  //
  // get handles
  //
  std::vector<int> &atomisticElements =
      *(computeElementPolarization_argument->atomisticElements);
  std::vector<std::vector<std::pair<std::vector<double>, double>>>
      &elementPolarizationLocations =
          *(computeElementPolarization_argument->elementPolarizationLocations);
  double initialUnitCellVolume =
      computeElementPolarization_argument->initialUnitCellVolume;
  struct element_list_t &elementList =
      *(computeElementPolarization_argument->elementList);
  std::vector<std::vector<double>> &elementPolarization =
      *(computeElementPolarization_argument->elementPolarization);

  //
  // variables for sites
  //
  int number_start;
  int number_end;

  //
  // keep getting more elements
  //
  while (whileWorker) {
    //
    // get element lock
    //
    pthread_mutex_lock(&dataLock);

    //
    // set variables
    //
    if (currentElem == int(atomisticElements.size())) {
      pthread_mutex_unlock(&dataLock);
      break;
    } else if (currentElem + elemBucket >= int(atomisticElements.size())) {
      //
      // set variables
      //
      number_start = currentElem;
      number_end = atomisticElements.size();
      currentElem = number_end;
    } else {
      //
      // set variables
      //
      number_start = currentElem;
      number_end = number_start + elemBucket;
      currentElem = number_end;
    }

    //
    // get element lock
    //
    pthread_mutex_unlock(&dataLock);

    //
    // loop over elements
    //
    for (int iElem = number_start; iElem < number_end; ++iElem) {
      //
      // check to see if element is not atomistic
      //
      if (atomisticElements[iElem] == 1) {
        //
        // calculate initial element volume
        //
        double initialElementVolume = tetrahedronVolume(
            elementList.elements[iElem]->node[0]->initial_position,
            elementList.elements[iElem]->node[1]->initial_position,
            elementList.elements[iElem]->node[2]->initial_position,
            elementList.elements[iElem]->node[3]->initial_position);

        //
        // calculate current element volume
        //
        double currentElementVolume =
            tetrahedronVolume(elementList.elements[iElem]->node[0]->position,
                              elementList.elements[iElem]->node[1]->position,
                              elementList.elements[iElem]->node[2]->position,
                              elementList.elements[iElem]->node[3]->position);

        //
        // caculate change in volume
        //
        double changeInVolume = currentElementVolume / initialElementVolume;

        //
        // calculate current unit cell volume
        //
        double unitCellVolume = changeInVolume * initialUnitCellVolume;

        //
        // calculate polarization
        //
        std::vector<double> polarization(3, 0.0);

        for (int iAtom = 0;
             iAtom < int(elementPolarizationLocations[iElem].size()); ++iAtom)
          for (int iDof = 0; iDof < 3; ++iDof)
            polarization[iDof] +=
                (elementPolarizationLocations[iElem][iAtom].first[iDof] -
                 elementPolarizationLocations[iElem][0].first[iDof]) *
                elementPolarizationLocations[iElem][iAtom].second;

        //
        // calculate charge in unit cell
        //
        double chargeCheck = 0.0;
        for (int iAtom = 0;
             iAtom < int(elementPolarizationLocations[iElem].size()); ++iAtom)
          chargeCheck += elementPolarizationLocations[iElem][iAtom].second;

        //
        // check to make sure charge is 0
        //
        assert(chargeCheck > -shiftTolerance && chargeCheck < shiftTolerance);

        //
        // calculate polarization density
        //
        for (int iDof = 0; iDof < 3; ++iDof)
          polarization[iDof] = polarization[iDof] / unitCellVolume;

        //
        // set element polarization
        //
        elementPolarization[iElem].clear();
        for (int iDof = 0; iDof < 3; ++iDof)
          elementPolarization[iElem].push_back(polarization[iDof]);
      } else if (atomisticElements[iElem] == 0) {
        //
        // if atomistic element set polarization to 0
        //
        elementPolarization[iElem].clear();
        for (int iDof = 0; iDof < 3; ++iDof)
          elementPolarization[iElem].push_back(0.0);
      }
    }
  }

  //
  //
  //
  return NULL;
}

//
// compute charge densities
//
void *computeChargeDensitiesWorker(void *arg) {
  // cast back to get computeAtomisticElementLocations_argument_t
  struct computeChargeDensities_argument_t *computeChargeDensities_argument =
      static_cast<struct computeChargeDensities_argument_t *>(arg);

  // get handles
  std::vector<std::vector<std::vector<int>>> &faces =
      *(computeChargeDensities_argument->faces);
  std::vector<std::vector<std::pair<int, int>>> &faceElementList =
      *(computeChargeDensities_argument->faceElementList);
  struct node_list_t &nodeList = *(computeChargeDensities_argument->nodeList);
  std::vector<std::vector<double>> &elementPolarization =
      *(computeChargeDensities_argument->elementPolarization);
  std::vector<std::vector<std::pair<double, std::vector<std::vector<double>>>>>
      &globalChargeDensity = *(computeChargeDensities_argument->chargeDensity);
  struct lattice_t &lattice = *(computeChargeDensities_argument->lattice);
  std::vector<std::pair<int, double>> &electroBoundaries =
      *(computeChargeDensities_argument->electroBoundaries);

  // variables for buckets
  int number_start;
  int number_end;

  // keep getting more buckets
  while (whileWorker) {
    // get bucket count lock
    pthread_mutex_lock(&dataLock);

    // set variables
    if (currentBucket == numBuckets) {
      pthread_mutex_unlock(&dataLock);
      break;
    } else if (currentBucket + elemBucket >= numBuckets) {
      // set variables
      number_start = currentBucket;
      number_end = numBuckets;
      currentBucket = number_end;
    } else {
      // set variables
      number_start = currentBucket;
      number_end = number_start + elemBucket;
      currentBucket = number_end;
    }

    // unlock bucket lock
    pthread_mutex_unlock(&dataLock);

    // loop over elements
    for (int iBucket = number_start; iBucket < number_end; ++iBucket) {
      // clear charge density
      globalChargeDensity[iBucket].clear();

      // loop over faces in bucket
      for (int iFace = 0; iFace < int(faces[iBucket].size()); ++iFace) {
        // check sizes
        assert(faces[iBucket].size() == faceElementList[iBucket].size());

        // density to calculate
        double chargeDensity = 0.0;

        // get first element
        const int element1 = faceElementList[iBucket][iFace].first;

        // get second element
        const int element2 = faceElementList[iBucket][iFace].second;

        // get nodes
        std::vector<int> faceNodes;
        for (int iNode = 0; iNode < 3; ++iNode)
          faceNodes.push_back(faces[iBucket][iFace][iNode]);

        // 4th nodes of first element
        int fourthNode = faces[iBucket][iFace][3];

        // get node locations
        std::vector<std::vector<double>> nodeLocations;
        std::vector<double> location(3, 0.0);
        for (int iNode = 0; iNode < 3; ++iNode) {
          for (int iDof = 0; iDof < 3; ++iDof) {
            // get location
            location[iDof] = nodeList.nodes[faceNodes[iNode]]->position[iDof];
          }

          // populate node locations
          nodeLocations.push_back(location);
        }

        // fourth node location
        std::vector<double> fourthNodeLocation(3, 0.0);
        for (int iDof = 0; iDof < 3; ++iDof) {
          // get location
          fourthNodeLocation[iDof] = nodeList.nodes[fourthNode]->position[iDof];
        }

        // calculate outward normal between elements from first element
        std::vector<double> normal =
            outwardNormal(nodeLocations, fourthNodeLocation);

        // add fourth node locations to nodeLocations
        nodeLocations.push_back(fourthNodeLocation);

        // get polarization on side 1
        std::vector<double> polar1(3, 0.0);
        polar1[0] = elementPolarization[element1][0];
        polar1[1] = elementPolarization[element1][1];
        polar1[2] = elementPolarization[element1][2];

        // get polarization on side 2
        std::vector<double> polar2(3, 0.0);

        // if there is no second element, then face has a void or edge on the
        // other side
        // assume polar2 = 0, can also assume polar2 = polar1 for no
        // contribution, see commented lines below
        if (element2 != -1) {
          // get polarization in element 2
          polar2[0] = elementPolarization[element2][0];
          polar2[1] = elementPolarization[element2][1];
          polar2[2] = elementPolarization[element2][2];

          // add second fourth node location because element2 must exist
          nodeLocations.resize(5);
          nodeLocations[4].resize(3);

          // get 2nd fourth node
          int fourthNode2 = faces[iBucket][iFace][4];

          // add location of 2nd fourth node to data structure
          for (int iDof = 0; iDof < 3; ++iDof)
            nodeLocations[4][iDof] =
                nodeList.nodes[fourthNode2]->position[iDof];

          // calculate charge density on face
          for (int iDof = 0; iDof < 3; ++iDof)
            chargeDensity += (polar1[iDof] - polar2[iDof]) * normal[iDof];

          // add to local charge density structure
          globalChargeDensity[iBucket].push_back(
              std::make_pair(chargeDensity, nodeLocations));
        } // end if element 2 exists
        else {
          // need to check for atomistic nodes
          // need to check for electro boundary conditions
          polar2[0] = 0.0;
          polar2[1] = 0.0;
          polar2[2] = 0.0;

          // calculate charge density on face
          for (int iDof = 0; iDof < 3; ++iDof)
            chargeDensity += (polar1[iDof] - polar2[iDof]) * normal[iDof];

          // check xBegin surface
          if (nodeList.nodes[faceNodes[0]]->l[0] ==
                  nodeList.nodes[faceNodes[1]]->l[0] &&
              nodeList.nodes[faceNodes[0]]->l[0] ==
                  nodeList.nodes[faceNodes[2]]->l[0] &&
              nodeList.nodes[faceNodes[0]]->l[0] == lattice.l_start[0]) {
            if (electroBoundaries[0].first == 0) {
              chargeDensity = 0.0;
            } else if (electroBoundaries[0].first == 1) {
              chargeDensity = 0.0;
              for (int iDof = 0; iDof < 3; ++iDof)
                chargeDensity += polar1[iDof] * normal[iDof];
            } else if (electroBoundaries[0].first == 2) {
              chargeDensity = electroBoundaries[0].second;
            } else {
              d_print("Electrostatic BC not implemented\n");
              exit(0);
            }
          }

          // check yEnd surface
          if (nodeList.nodes[faceNodes[0]]->l[0] ==
                  nodeList.nodes[faceNodes[1]]->l[0] &&
              nodeList.nodes[faceNodes[0]]->l[0] ==
                  nodeList.nodes[faceNodes[2]]->l[0] &&
              nodeList.nodes[faceNodes[0]]->l[0] == lattice.l_end[0]) {
            if (electroBoundaries[1].first == 0) {
              chargeDensity = 0.0;
            } else if (electroBoundaries[1].first == 1) {
              chargeDensity = 0.0;
              for (int iDof = 0; iDof < 3; ++iDof)
                chargeDensity += polar1[iDof] * normal[iDof];
            } else if (electroBoundaries[1].first == 2) {
              chargeDensity = electroBoundaries[1].second;
            } else {
              d_print("Electrostatic BC not implemented\n");
              exit(0);
            }
          }

          // check yBegin surface
          if (nodeList.nodes[faceNodes[0]]->l[1] ==
                  nodeList.nodes[faceNodes[1]]->l[1] &&
              nodeList.nodes[faceNodes[0]]->l[1] ==
                  nodeList.nodes[faceNodes[2]]->l[1] &&
              nodeList.nodes[faceNodes[0]]->l[1] == lattice.l_start[1]) {
            if (electroBoundaries[2].first == 0) {
              chargeDensity = 0.0;
            } else if (electroBoundaries[2].first == 1) {
              chargeDensity = 0.0;
              for (int iDof = 0; iDof < 3; ++iDof)
                chargeDensity += polar1[iDof] * normal[iDof];
            } else if (electroBoundaries[2].first == 2) {
              chargeDensity = electroBoundaries[2].second;
            } else {
              d_print("Electrostatic BC not implemented\n");
              exit(0);
            }
          }

          // check yEnd surface
          if (nodeList.nodes[faceNodes[0]]->l[1] ==
                  nodeList.nodes[faceNodes[1]]->l[1] &&
              nodeList.nodes[faceNodes[0]]->l[1] ==
                  nodeList.nodes[faceNodes[2]]->l[1] &&
              nodeList.nodes[faceNodes[0]]->l[1] == lattice.l_end[1]) {
            if (electroBoundaries[3].first == 0) {
              chargeDensity = 0.0;
            } else if (electroBoundaries[3].first == 1) {
              chargeDensity = 0.0;
              for (int iDof = 0; iDof < 3; ++iDof)
                chargeDensity += polar1[iDof] * normal[iDof];
            } else if (electroBoundaries[3].first == 2) {
              chargeDensity = electroBoundaries[3].second;
            } else {
              d_print("Electrostatic BC not implemented\n");
              exit(0);
            }
          }

          // check zBegin surface
          if (nodeList.nodes[faceNodes[0]]->l[2] ==
                  nodeList.nodes[faceNodes[1]]->l[2] &&
              nodeList.nodes[faceNodes[0]]->l[2] ==
                  nodeList.nodes[faceNodes[2]]->l[2] &&
              nodeList.nodes[faceNodes[0]]->l[2] == lattice.l_start[2]) {
            if (electroBoundaries[4].first == 0) {
              chargeDensity = 0.0;
            } else if (electroBoundaries[4].first == 1) {
              chargeDensity = 0.0;
              for (int iDof = 0; iDof < 3; ++iDof)
                chargeDensity += polar1[iDof] * normal[iDof];
            } else if (electroBoundaries[4].first == 2) {
              chargeDensity = electroBoundaries[4].second;
            } else {
              d_print("Electrostatic BC not implemented\n");
              exit(0);
            }
          }

          // check zEnd surface
          if (nodeList.nodes[faceNodes[0]]->l[2] ==
                  nodeList.nodes[faceNodes[1]]->l[2] &&
              nodeList.nodes[faceNodes[0]]->l[2] ==
                  nodeList.nodes[faceNodes[2]]->l[2] &&
              nodeList.nodes[faceNodes[0]]->l[2] == lattice.l_end[2]) {
            if (electroBoundaries[5].first == 0) {
              chargeDensity = 0.0;
            } else if (electroBoundaries[5].first == 1) {
              chargeDensity = 0.0;
              for (int iDof = 0; iDof < 3; ++iDof)
                chargeDensity += polar1[iDof] * normal[iDof];
            } else if (electroBoundaries[5].first == 2) {
              chargeDensity = electroBoundaries[5].second;
            } else {
              d_print("Electrostatic BC not implemented\n");
              exit(0);
            }
          }

          // add to local charge density structure
          globalChargeDensity[iBucket].push_back(
              std::make_pair(chargeDensity, nodeLocations));
        } // end element 2 does not exist
      }   // end for loop over faces in bucket
    }     // end for loop over buckets
  }       // end while loop over start/stop buckets

  //
  return NULL;
} // end computeChargeDensitiesWorker

//
// compute atomistic charge locations
//
void *computeAtomisticChargeStateWorker(void *arg) {
  //
  // cast back to get computeAtomisticChargeLocations_argument_t
  //
  struct computeAtomisticChargeState_argument_t
      *computeAtomisticChargeState_argument =
          static_cast<struct computeAtomisticChargeState_argument_t *>(arg);

  //
  // get handles
  //
  std::vector<std::vector<std::vector<int>>> &atomisticCharges =
      *(computeAtomisticChargeState_argument->atomisticCharges);
  std::vector<std::vector<std::vector<double>>> &atomisticChargeState =
      *(computeAtomisticChargeState_argument->atomisticChargeState);
  int iQuasi = computeAtomisticChargeState_argument->quasiId;
  struct lattice_t &lattice = *(computeAtomisticChargeState_argument->lattice);

  //
  // variables for sites
  //
  int number_start;
  int number_end;

  //
  // local atomistic charge locations variable
  //
  std::vector<std::vector<double>> localChargeState;

  //
  // get shift vector
  //
  std::vector<double> shift = Quasicontinua::getInstance()->getShift(iQuasi);

  //
  // keep getting more elements
  //
  while (whileWorker) {
    //
    // get data lock
    //
    pthread_mutex_lock(&dataLock);

    //
    // set variables
    //
    if (currentBucket == numBuckets) {
      pthread_mutex_unlock(&dataLock);
      break;
    } else if (currentBucket + chargeBucket >= numBuckets) {
      //
      // set variables
      //
      number_start = currentBucket;
      number_end = numBuckets;
      currentBucket = number_end;
    } else {
      //
      // set variables
      //
      number_start = currentBucket;
      number_end = number_start + chargeBucket;
      currentBucket = number_end;
    }

    //
    // unlock data lock
    //
    pthread_mutex_unlock(&dataLock);

    //
    // temp vector to hold location
    //
    std::vector<double> tempState(4, 0.0);

    //
    // loop over buckets
    //
    for (int iBucket = number_start; iBucket < number_end; ++iBucket) {
      //
      // loop over charges
      //
      for (int iCharge = 0; iCharge < int(atomisticCharges[iBucket].size());
           ++iCharge) {
        //
        // temporary location
        //
        double state[4];

        //
        // temporary site
        //
        int site[3];
        site[0] = atomisticCharges[iBucket][iCharge][0];
        site[1] = atomisticCharges[iBucket][iCharge][1];
        site[2] = atomisticCharges[iBucket][iCharge][2];

        //
        // get site location
        //
        // if(quasicontinuum::is_site_inside_lattice(&lattice,site) ==
        // RETURN_FAILURE)
        if (Void::getInstance()->findSiteInVoidCache(
                atomisticCharges[iBucket][iCharge], iQuasi) == true)
          continue;

        Lattice::getInstance()->getSiteCurrentState(state, site, &lattice,
                                                    iQuasi);

        //
        // add shift to temp vector
        //
        tempState[0] = state[0] + shift[0];
        tempState[1] = state[1] + shift[1];
        tempState[2] = state[2] + shift[2];
        tempState[3] = state[3];

        //
        // add actual location to local data structure
        //
        localChargeState.push_back(tempState);
      }
    }
  }

  //
  // get sumLock
  //
  pthread_mutex_lock(&sumLock);

  //
  // add local charge locations to global data structure
  //
  atomisticChargeState[iQuasi].insert(atomisticChargeState[iQuasi].end(),
                                      localChargeState.begin(),
                                      localChargeState.end());

  //
  // unlock sumLock
  //
  pthread_mutex_unlock(&sumLock);

  //
  //
  //
  return NULL;
}

//
// compute potential and field from atomistic charges
//
void *computeFieldFromAtomisticChargesWorker(void *arg) {
  // cast back to get computeFieldFromAtomisticCharges_argument_t
  struct computeFieldFromAtomisticCharges_argument_t
      *computeFieldFromAtomisticCharges_argument =
          static_cast<struct computeFieldFromAtomisticCharges_argument_t *>(
              arg);

  // get handles
  const struct node_list_t &nodeList =
      *(computeFieldFromAtomisticCharges_argument->nodeList);
  std::vector<std::vector<double>> &electricField =
      *(computeFieldFromAtomisticCharges_argument->electricField);
  std::vector<double> &electricPotential =
      *(computeFieldFromAtomisticCharges_argument->electricPotential);
  const std::vector<double> &fixedCharge =
      *(computeFieldFromAtomisticCharges_argument->fixedCharge);
  const double &electricConstant =
      *(computeFieldFromAtomisticCharges_argument->electricConstant);
  const std::vector<std::vector<std::vector<double>>> &atomisticChargeState =
      *(computeFieldFromAtomisticCharges_argument->atomisticChargeState);
  const std::vector<double> &iShift =
      *(computeFieldFromAtomisticCharges_argument->iShift);
  const int &iQuasi = *(computeFieldFromAtomisticCharges_argument->iQuasi);

  // variables for nodes
  int number_start;
  int number_end;

  // keep getting more nodes
  while (whileWorker) {
    // get node lock
    pthread_mutex_lock(&dataLock);

    // set variables
    if (currentNode == nodeList.number_nodes) {
      pthread_mutex_unlock(&dataLock);
      break;
    } else if (currentNode + nodeBucket >= nodeList.number_nodes) {
      // set variables
      number_start = currentNode;
      number_end = nodeList.number_nodes;
      currentNode = number_end;
    } else {
      // set variables
      number_start = currentNode;
      number_end = number_start + nodeBucket;
      currentNode = number_end;
    }

    // unlock data lock
    pthread_mutex_unlock(&dataLock);

    //
    //  core flag
    //
    int iCoreFlag = Quasicontinua::getInstance()->getCoreShell(iQuasi);

    // node location
    std::vector<double> node_state(4, 0.0);

    // loop over nodes
    for (int i_node = number_start; i_node < number_end; ++i_node) {
      // get global node location
      for (int i_dof = 0; i_dof < 3; ++i_dof)
        node_state[i_dof] =
            nodeList.nodes[i_node]->position[i_dof] + iShift[i_dof];

      node_state[3] = nodeList.nodes[i_node]->frequency;

      // get field from atomistic charges
      for (int j_quasi = 0; j_quasi < int(atomisticChargeState.size());
           ++j_quasi) {
        //
        //  get core flag
        //
        int jCoreFlag = Quasicontinua::getInstance()->getCoreShell(j_quasi);

        std::pair<int, int> coreData;
        coreData.first = iCoreFlag;
        coreData.second = jCoreFlag;

        // set charge
        const double charge = fixedCharge[j_quasi];

        // loop over all charges in i_quasi
        for (int j_charge = 0; j_charge < atomisticChargeState[j_quasi].size();
             ++j_charge) {
          AddAtomicToFieldAndPotential(
              node_state, iQuasi, atomisticChargeState[j_quasi][j_charge],
              j_quasi, electricConstant, charge, electricField[i_node],
              electricPotential[i_node], coreData);
        }
      } // end of loop over quasicontinuum
    }   // end for loop over nodes in list
  }     // end while loop over nodes

  //
  return NULL;
} // computeFieldFromAtomisticChargesWorker

//
// compute contribution from cluster charges
//
void *computeClusterChargesWorker(void *arg) {
  // cast back to get computeClusterCharges_argument_t
  struct computeClusterCharges_argument_t *computeClusterCharges_argument =
      static_cast<struct computeClusterCharges_argument_t *>(arg);

  // get handles
  const struct node_list_t &i_node_list =
      *(computeClusterCharges_argument->i_node_list);
  std::vector<std::vector<double>> &electric_field =
      *(computeClusterCharges_argument->electricField);
  std::vector<double> &electric_potential =
      *(computeClusterCharges_argument->electricPotential);
  const double &electric_constant =
      *(computeClusterCharges_argument->electricConstant);
  const std::vector<std::vector<std::vector<int>>> &atomistic_charges =
      *(computeClusterCharges_argument->atomisticCharges);
  struct lattice_t &j_lattice = *(computeClusterCharges_argument->j_lattice);
  const std::vector<double> &i_shift =
      *(computeClusterCharges_argument->i_shift);
  const std::vector<double> &j_shift =
      *(computeClusterCharges_argument->j_shift);
  const int &num_shells = *(computeClusterCharges_argument->num_shells);
  const double &j_charge = *(computeClusterCharges_argument->j_charge);
  const int &i_quasi = *(computeClusterCharges_argument->i_quasi);

  std::pair<int, int> &coreData = *(computeClusterCharges_argument->coreData);

  int &j_quasi = *(computeClusterCharges_argument->j_quasi);

  // variables for nodes
  int number_start;
  int number_end;

  // keep getting more nodes
  while (whileWorker) {
    // get node lock
    pthread_mutex_lock(&dataLock);

    // set variables
    if (currentNode == i_node_list.number_nodes) {
      pthread_mutex_unlock(&dataLock);
      break;
    } else if (currentNode + nodeBucket >= i_node_list.number_nodes) {
      // set variables
      number_start = currentNode;
      number_end = i_node_list.number_nodes;
      currentNode = number_end;
    } else {
      // set variables
      number_start = currentNode;
      number_end = number_start + nodeBucket;
      currentNode = number_end;
    }

    // unlock data lock
    pthread_mutex_unlock(&dataLock);

    // node location
    std::vector<double> node_state(4, 0.0);

    // loop over nodes
    for (int i_node = number_start; i_node < number_end; ++i_node) {
      // get global node location
      for (int i_dof = 0; i_dof < 3; ++i_dof)
        node_state[i_dof] =
            i_node_list.nodes[i_node]->position[i_dof] + i_shift[i_dof];

      node_state[3] = i_node_list.nodes[i_node]->frequency;

      // set base site index
      const std::vector<int> base_site_index = {
          i_node_list.nodes[i_node]->l[0], i_node_list.nodes[i_node]->l[1],
          i_node_list.nodes[i_node]->l[2]};

      // loop over shells of unit cells
      for (int x_site = base_site_index[0] - num_shells;
           x_site <= base_site_index[0] + num_shells; ++x_site) {
        for (int y_site = base_site_index[1] - num_shells;
             y_site <= base_site_index[1] + num_shells; ++y_site) {
          for (int z_site = base_site_index[2] - num_shells;
               z_site <= base_site_index[2] + num_shells; ++z_site) {
            // set site
            const std::vector<int> site = {x_site, y_site, z_site};

            // check if site is in lattice
            if (Lattice::getInstance()->isSiteInsideLattice(
                    &j_lattice, &site[0], j_quasi) == RETURN_FAILURE)
              continue;

            // check if site is atomistic
            if (isSiteAtomistic(site, atomistic_charges) == true)
              continue;

            // get spatial coordinates
            std::vector<double> jSiteState(4, 0.0);

            Lattice::getInstance()->getSiteCurrentState(
                &jSiteState[0], &site[0], &j_lattice, j_quasi);

            for (int i_dof = 0; i_dof < 3; ++i_dof)
              jSiteState[i_dof] += j_shift[i_dof];

            // add contribution to field
            AddAtomicToFieldAndPotential(node_state, i_quasi, jSiteState,
                                         j_quasi, electric_constant, j_charge,
                                         electric_field[i_node],
                                         electric_potential[i_node], coreData);
          } // end for loop over z_sites
        }   // end for loop over y_sites
      }     // end for loop over x_sites
    }       // end for loop over local nodes
  }         // end while loop over global nodes

  //
  return NULL;
} // end computeClusterChargesWorker

//
// compute potential and field from atomistic charges
//
void *ComputeFieldFromExternalChargesWorker(void *arg) {
  // cast back to get compute_field_from_external_charges_argument_t
  struct compute_field_from_external_charges_argument_t
      *compute_field_from_external_charges_argument =
          static_cast<struct compute_field_from_external_charges_argument_t *>(
              arg);

  // get handles
  const struct node_list_t &node_list =
      *(compute_field_from_external_charges_argument->node_list);
  std::vector<std::vector<double>> &electric_field =
      *(compute_field_from_external_charges_argument->electric_field);
  std::vector<double> &electric_potential =
      *(compute_field_from_external_charges_argument->electric_potential);
  const double &electric_constant =
      *(compute_field_from_external_charges_argument->electric_constant);
  const std::vector<std::pair<double, std::vector<double>>> &external_charges =
      *(compute_field_from_external_charges_argument->external_charges);
  const std::vector<double> &iShift =
      *(compute_field_from_external_charges_argument->iShift);
  const int &iQuasi = *(compute_field_from_external_charges_argument->iQuasi);
  int &iCoreFlag = *(compute_field_from_external_charges_argument->iCoreFlag);

  // variables for nodes
  int number_start;
  int number_end;

  // keep getting more nodes
  while (whileWorker) {
    // get node lock
    pthread_mutex_lock(&dataLock);

    // set variables
    if (currentNode == node_list.number_nodes) {
      pthread_mutex_unlock(&dataLock);
      break;
    } else if (currentNode + nodeBucket >= node_list.number_nodes) {
      // set variables
      number_start = currentNode;
      number_end = node_list.number_nodes;
      currentNode = number_end;
    } else {
      // set variables
      number_start = currentNode;
      number_end = number_start + nodeBucket;
      currentNode = number_end;
    }

    // unlock data lock
    pthread_mutex_unlock(&dataLock);

    // node location
    std::vector<double> node_state(4, 0.0);

    // loop over nodes
    for (int i_node = number_start; i_node < number_end; ++i_node) {
      // get global node location
      for (int i_dof = 0; i_dof < 3; ++i_dof)
        node_state[i_dof] =
            node_list.nodes[i_node]->position[i_dof] + iShift[i_dof];

      node_state[3] = node_list.nodes[i_node]->frequency;

      // get field from external charges
      for (int j_charge = 0; j_charge < external_charges.size(); ++j_charge) {
        // add to field and potential
        AddAtomicToFieldAndPotentialExternal(
            node_state, iQuasi, external_charges[j_charge].second,
            electric_constant, external_charges[j_charge].first,
            electric_field[i_node], electric_potential[i_node], iCoreFlag);
      } // end of loop over charges
    }   // end for loop over nodes in list
  }     // end while loop over nodes

  //
  return NULL;
} // ComputeFieldFromExternalChargesWorker

//
// compute contribution from cluster charge surfaces
//
void *computeClusterChargeSurfaceWorker(void *arg) {
  // cast back to get computeClusterChargeSurface_argument_t
  struct computeClusterChargeSurface_argument_t
      *computeClusterChargeSurface_argument =
          static_cast<struct computeClusterChargeSurface_argument_t *>(arg);

  // get handles
  const struct node_list_t &i_node_list =
      *(computeClusterChargeSurface_argument->i_node_list);
  std::vector<std::vector<double>> &electric_field =
      *(computeClusterChargeSurface_argument->electricField);
  std::vector<double> &electric_potential =
      *(computeClusterChargeSurface_argument->electricPotential);
  const double &electric_constant =
      *(computeClusterChargeSurface_argument->electricConstant);
  const std::vector<std::vector<std::vector<int>>> &atomistic_charges =
      *(computeClusterChargeSurface_argument->atomisticCharges);
  const std::vector<std::vector<double>> &element_polarization =
      *(computeClusterChargeSurface_argument->element_polarization);
  struct lattice_t &zero_lattice =
      *(computeClusterChargeSurface_argument->zero_lattice);
  const std::vector<double> &i_shift =
      *(computeClusterChargeSurface_argument->i_shift);
  const std::vector<double> &zero_shift =
      *(computeClusterChargeSurface_argument->zero_shift);
  const int &num_shells = *(computeClusterChargeSurface_argument->num_shells);
  const int &face_method = *(computeClusterChargeSurface_argument->face_method);
  const double &integration_error =
      *(computeClusterChargeSurface_argument->integration_error);
  const struct node_list_t &zero_node_list =
      *(computeClusterChargeSurface_argument->zero_node_list);
  const std::vector<std::vector<std::pair<std::vector<int>, int>>>
      &node_hash_list =
          *(computeClusterChargeSurface_argument->base_node_hash_list);
  int zero_quasi = *(computeClusterChargeSurface_argument->zero_quasi);

  // variables for nodes
  int number_start;
  int number_end;

  // keep getting more nodes
  while (whileWorker) {
    // get node lock
    pthread_mutex_lock(&dataLock);

    // set variables
    if (currentNode == i_node_list.number_nodes) {
      pthread_mutex_unlock(&dataLock);
      break;
    } else if (currentNode + nodeBucket >= i_node_list.number_nodes) {
      // set variables
      number_start = currentNode;
      number_end = i_node_list.number_nodes;
      currentNode = number_end;
    } else {
      // set variables
      number_start = currentNode;
      number_end = number_start + nodeBucket;
      currentNode = number_end;
    }

    // unlock data lock
    pthread_mutex_unlock(&dataLock);

    // node location
    std::vector<double> node_location(3, 0.0);

    // loop over nodes
    for (int i_node = number_start; i_node < number_end; ++i_node) {
      // get global node location
      for (int i_dof = 0; i_dof < 3; ++i_dof)
        node_location[i_dof] =
            i_node_list.nodes[i_node]->position[i_dof] + i_shift[i_dof];

      // set base site index
      const std::vector<int> base_site_index = {
          i_node_list.nodes[i_node]->l[0], i_node_list.nodes[i_node]->l[1],
          i_node_list.nodes[i_node]->l[2]};

      // center of face to iterate over
      std::vector<int> face_site(3, 0);
      int const_index;
      int iterate_index_1;
      int iterate_index_2;
      double test_shift;

      // get contribution from -x face
      face_site[0] = base_site_index[0] - num_shells;
      face_site[1] = base_site_index[1];
      face_site[2] = base_site_index[2];
      const_index = 0;
      iterate_index_1 = 1;
      iterate_index_2 = 2;
      test_shift = -0.1;
      if (face_site[0] != zero_lattice.l_start[0])
        ContributionFromSingleClusterSurface(
            zero_quasi, node_location, face_method, integration_error,
            electric_constant, face_site, const_index, iterate_index_1,
            iterate_index_2, num_shells, element_polarization,
            atomistic_charges, zero_shift, test_shift, zero_lattice,
            zero_node_list, node_hash_list, electric_field[i_node],
            electric_potential[i_node]);

      // get contribution from +x face
      face_site[0] = base_site_index[0] + num_shells + 1;
      face_site[1] = base_site_index[1];
      face_site[2] = base_site_index[2];
      const_index = 0;
      iterate_index_1 = 1;
      iterate_index_2 = 2;
      test_shift = 0.1;
      if (face_site[0] != zero_lattice.l_end[0])
        ContributionFromSingleClusterSurface(
            zero_quasi, node_location, face_method, integration_error,
            electric_constant, face_site, const_index, iterate_index_1,
            iterate_index_2, num_shells, element_polarization,
            atomistic_charges, zero_shift, test_shift, zero_lattice,
            zero_node_list, node_hash_list, electric_field[i_node],
            electric_potential[i_node]);

      // get contribution from -y face
      face_site[0] = base_site_index[0];
      face_site[1] = base_site_index[1] - num_shells;
      face_site[2] = base_site_index[2];
      const_index = 1;
      iterate_index_1 = 0;
      iterate_index_2 = 2;
      test_shift = -0.1;
      if (face_site[1] != zero_lattice.l_start[1])
        ContributionFromSingleClusterSurface(
            zero_quasi, node_location, face_method, integration_error,
            electric_constant, face_site, const_index, iterate_index_1,
            iterate_index_2, num_shells, element_polarization,
            atomistic_charges, zero_shift, test_shift, zero_lattice,
            zero_node_list, node_hash_list, electric_field[i_node],
            electric_potential[i_node]);

      // get contribution from +y face
      face_site[0] = base_site_index[0];
      face_site[1] = base_site_index[1] + num_shells + 1;
      face_site[2] = base_site_index[2];
      const_index = 1;
      iterate_index_1 = 0;
      iterate_index_2 = 2;
      test_shift = 0.1;
      if (face_site[1] != zero_lattice.l_end[1])
        ContributionFromSingleClusterSurface(
            zero_quasi, node_location, face_method, integration_error,
            electric_constant, face_site, const_index, iterate_index_1,
            iterate_index_2, num_shells, element_polarization,
            atomistic_charges, zero_shift, test_shift, zero_lattice,
            zero_node_list, node_hash_list, electric_field[i_node],
            electric_potential[i_node]);

      // get contribution from -z face
      face_site[0] = base_site_index[0];
      face_site[1] = base_site_index[1];
      face_site[2] = base_site_index[2] - num_shells;
      const_index = 2;
      iterate_index_1 = 0;
      iterate_index_2 = 1;
      test_shift = -0.1;
      if (face_site[2] != zero_lattice.l_start[2])
        ContributionFromSingleClusterSurface(
            zero_quasi, node_location, face_method, integration_error,
            electric_constant, face_site, const_index, iterate_index_1,
            iterate_index_2, num_shells, element_polarization,
            atomistic_charges, zero_shift, test_shift, zero_lattice,
            zero_node_list, node_hash_list, electric_field[i_node],
            electric_potential[i_node]);

      // get contribution from +z face
      face_site[0] = base_site_index[0];
      face_site[1] = base_site_index[1];
      face_site[2] = base_site_index[2] + num_shells + 1;
      const_index = 2;
      iterate_index_1 = 0;
      iterate_index_2 = 1;
      test_shift = 0.1;
      if (face_site[2] != zero_lattice.l_end[2])
        ContributionFromSingleClusterSurface(
            zero_quasi, node_location, face_method, integration_error,
            electric_constant, face_site, const_index, iterate_index_1,
            iterate_index_2, num_shells, element_polarization,
            atomistic_charges, zero_shift, test_shift, zero_lattice,
            zero_node_list, node_hash_list, electric_field[i_node],
            electric_potential[i_node]);
    } // end for loop over local nodes
  }   // end while loop over global nodes

  //
  return NULL;
} // end computeClusterChargeSurfaceWorker

//
// compute contribution from charged faces
//
void *computeChargedFaceWorker(void *arg) {
  // cast back to get computeChargedFace_argument_t
  struct computeChargedFace_argument_t *computeChargedFace_argument =
      static_cast<struct computeChargedFace_argument_t *>(arg);

  // get handles
  const struct node_list_t &i_node_list =
      *(computeChargedFace_argument->i_node_list);
  std::vector<std::vector<double>> &electric_field =
      *(computeChargedFace_argument->electricField);
  std::vector<double> &electric_potential =
      *(computeChargedFace_argument->electricPotential);
  const double &electric_constant =
      *(computeChargedFace_argument->electricConstant);
  const std::vector<double> &i_shift = *(computeChargedFace_argument->i_shift);
  const std::vector<double> &zero_shift =
      *(computeChargedFace_argument->zero_shift);
  const int &num_shells = *(computeChargedFace_argument->num_shells);
  const int &face_method = *(computeChargedFace_argument->face_method);
  const double &integration_error =
      *(computeChargedFace_argument->integration_error);
  const struct node_list_t &zero_node_list =
      *(computeChargedFace_argument->zero_node_list);
  const std::vector<std::vector<std::vector<int>>> &faces =
      *(computeChargedFace_argument->faces);
  const std::vector<
      std::vector<std::pair<double, std::vector<std::vector<double>>>>>
      &charge_density = *(computeChargedFace_argument->charge_density);
  struct lattice_t &zero_lattice = *(computeChargedFace_argument->zero_lattice);
  const std::vector<std::vector<std::vector<int>>> &atomistic_charges =
      *(computeChargedFace_argument->atomistic_charges);
  const int zero_quasi = *(computeChargedFace_argument->zero_quasi);

  // variables for nodes
  int number_start;
  int number_end;

  // keep getting more nodes
  while (whileWorker) {
    // get node lock
    pthread_mutex_lock(&dataLock);

    // set variables
    if (currentNode == i_node_list.number_nodes) {
      pthread_mutex_unlock(&dataLock);
      break;
    } else if (currentNode + nodeBucket >= i_node_list.number_nodes) {
      // set variables
      number_start = currentNode;
      number_end = i_node_list.number_nodes;
      currentNode = number_end;
    } else {
      // set variables
      number_start = currentNode;
      number_end = number_start + nodeBucket;
      currentNode = number_end;
    }

    // unlock data lock
    pthread_mutex_unlock(&dataLock);

    // node site
    std::vector<int> node_site(3, 0);

    // loop over nodes
    for (int i_node = number_start; i_node < number_end; ++i_node) {
      // get global node site
      for (int i_dof = 0; i_dof < 3; ++i_dof)
        node_site[i_dof] = i_node_list.nodes[i_node]->l[i_dof];

      // get global node location
      std::vector<double> node_location(3, 0.0);
      for (int i_dof = 0; i_dof < 3; ++i_dof)
        node_location[i_dof] =
            i_node_list.nodes[i_node]->position[i_dof] + i_shift[i_dof];

      // get coordinates of shell box
      const std::vector<int> x_coords = {
          std::max(node_site[0] - num_shells, zero_lattice.l_start[0]),
          std::min(node_site[0] + num_shells + 1, zero_lattice.l_end[0])};
      const std::vector<int> y_coords = {
          std::max(node_site[1] - num_shells, zero_lattice.l_start[1]),
          std::min(node_site[1] + num_shells + 1, zero_lattice.l_end[1])};
      const std::vector<int> z_coords = {
          std::max(node_site[2] - num_shells, zero_lattice.l_start[2]),
          std::min(node_site[2] + num_shells + 1, zero_lattice.l_end[2])};
      const std::vector<std::vector<int>> coords = {x_coords, y_coords,
                                                    z_coords};

#if 0
          // get convex hull of shell surface
          Convex_hull_3 ch_shell_surface(3);
          for(int x_it = 0; x_it < 2; ++x_it)
            for(int y_it = 0; y_it < 2; ++y_it)
              for(int z_it = 0; z_it < 2; ++z_it)
          ch_shell_surface.insert(Point_3(x_coords[x_it],
                  y_coords[y_it],
                  z_coords[z_it]));
          assert(ch_shell_surface.is_valid());
          std::cout<<__LINE__<<std::endl;
          // define polyhedron to hold convex hull
          Polyhedron_3 ph_shell_surface;
          CGAL::convex_hull_d_to_polyhedron_3(ch_shell_surface,
                      ph_shell_surface);
          std::cout<<__LINE__<<std::endl;
          // get nef polyhedron
          Nef_polyhedron nef_shell_surface(ph_shell_surface);

          // get net polyhedron without atomistic section
          nef_shell_surface -= nef_atomistic;
          nef_shell_surface.convert_to_polyhedron(ph_shell_surface);
          std::cout<<__LINE__<<std::endl;
          // check if polyhedron is all triangles
          if(ph_shell_surface.is_pure_triangle() == false)
            for(Polyhedron_3::Halfedge_iterator edge_it = ph_shell_surface.halfedges_begin();
          edge_it != ph_shell_surface.halfedges_end();
          edge_it++){
              if(ph_shell_surface.is_triangle(edge_it) == false)
          ph_shell_surface.create_center_vertex(edge_it);

            }
          assert(ph_shell_surface.is_pure_triangle() == true);
          std::cout<<__LINE__<<std::endl;
          // get triangles from polyhedron
          std::vector<Triangle_3> shell_triangles;
          for(Polyhedron_3::Facet_iterator facet_it = ph_shell_surface.facets_begin();
              facet_it != ph_shell_surface.facets_end();
              ++facet_it){

            // get triangle locations
            std::vector<Point_3> triangle;
            HF_circulator edge = facet_it->facet_begin();
            CGAL_assertion(CGAL::circulator_size(edge) == 3);
            do{
              Point_3 point = edge->vertex()->point();
              triangle.push_back(point);
            }while(++edge != facet_it->facet_begin());

            // cast triangle to ints, should be correct
            // need to check
            std::vector<Point_3> triangle_vertexes;
            for(int i_point = 0; i_point < 3; ++i_point){
              std::vector<double> temp_location(3,0.0);
              std::vector<int> temp_site = {int(CGAL::to_double(triangle[i_point][0])),
                    int(CGAL::to_double(triangle[i_point][1])),
                    int(CGAL::to_double(triangle[i_point][2]))};

              // get site coordinates
              get_site_spatial_coordinates(&temp_location[0],
                   &temp_site[0],
                   &zero_lattice);
              for(int i_dof = 0; i_dof < 3; ++i_dof)
          temp_location[i_dof] += zero_shift[i_dof];

              // insert points into triangle
              triangle_vertexes.push_back(Point_3(temp_location[0],
                    temp_location[1],
                    temp_location[2]));
            }

            // insert triangle into list
            shell_triangles.push_back(Triangle_3(triangle_vertexes[0],
                   triangle_vertexes[1],
                   triangle_vertexes[2]));

          }
#endif

      // get triangles for intersections
      std::vector<Triangle_3> shell_triangles;

      GetShellOfTriangles(coords, zero_lattice, zero_shift, atomistic_charges,
                          shell_triangles, zero_quasi);

      // loop over faces
      for (int i_bucket = 0; i_bucket < numBuckets; ++i_bucket) {
        for (int i_face = 0; i_face < charge_density[i_bucket].size();
             ++i_face) {
          // get face node sites
          std::vector<int> face_node_numbers;
          for (int i_node = 0; i_node < 3; ++i_node)
            face_node_numbers.push_back(faces[i_bucket][i_face][i_node]);

          // get node sites
          std::vector<std::vector<int>> face_node_sites;
          for (int i_node = 0; i_node < 3; ++i_node) {
            std::vector<int> temp_site;
            for (int i_dof = 0; i_dof < 3; ++i_dof)
              temp_site.push_back(
                  zero_node_list.nodes[face_node_numbers[i_node]]->l[i_dof]);
            face_node_sites.push_back(temp_site);
          }

          // check if face nodes are inside cluster
          std::vector<bool> face_node_in_cluster;
          for (int i_node = 0; i_node < 3; ++i_node)
            face_node_in_cluster.push_back(
                NodeInCluster(node_site, num_shells, face_node_sites[i_node]));

          // if face is completely in cluster, continue onto next node
          if (face_node_in_cluster[0] == true &&
              face_node_in_cluster[1] == true &&
              face_node_in_cluster[2] == true)
            continue;

          // get face node global locations for CGAL
          std::vector<Point_3> face_node_locations_point;
          for (int i_node = 0; i_node < 3; ++i_node)
            face_node_locations_point.push_back(

                Point_3(zero_node_list.nodes[face_node_numbers[i_node]]
                                ->position[0] +
                            zero_shift[0],
                        zero_node_list.nodes[face_node_numbers[i_node]]
                                ->position[1] +
                            zero_shift[1],
                        zero_node_list.nodes[face_node_numbers[i_node]]
                                ->position[2] +
                            zero_shift[2]));

          Triangle_3 face_node_triangle = {face_node_locations_point[0],
                                           face_node_locations_point[1],
                                           face_node_locations_point[2]};

          // get face node global locations for calculation
          std::vector<std::vector<double>> face_node_locations_double;
          for (int i_node = 0; i_node < 3; ++i_node) {
            std::vector<double> temp_location;

            for (int i_dof = 0; i_dof < 3; ++i_dof)
              temp_location.push_back(
                  zero_node_list.nodes[face_node_numbers[i_node]]
                      ->position[i_dof] +
                  zero_shift[i_dof]);

            face_node_locations_double.push_back(temp_location);
          }

          // add full plane to field subtract parts inside below
          const std::pair<std::vector<double>, double> pos_return_field =
              fieldFromDensityFaceAdaptive(
                  charge_density[i_bucket][i_face].first,
                  face_node_locations_double, node_location, electric_constant,
                  face_method, integration_error);

          // add to field
          for (int i_dof = 0; i_dof < 3; ++i_dof)
            electric_field[i_node][i_dof] += pos_return_field.first[i_dof];

          // add to potential
          electric_potential[i_node] += pos_return_field.second;

          // get intersectons
          std::vector<Point_3> intersection_points;

          for (int i_triangle = 0; i_triangle < shell_triangles.size();
               ++i_triangle) {
            // check intersection
            CGAL::Object result =
                intersection(shell_triangles[i_triangle], face_node_triangle);

            // potential return values
            Point_3 return_point;
            Segment_3 return_segment;
            Triangle_3 return_triangle;
            std::vector<Point_3> return_vector_points;

            // see what is returned
            if (assign(return_point, result)) {
              intersection_points.push_back(return_point);
            } else if (assign(return_segment, result)) {
              intersection_points.push_back(return_segment[0]);
              intersection_points.push_back(return_segment[1]);
            } else if (assign(return_triangle, result)) {
              intersection_points.push_back(return_triangle[0]);
              intersection_points.push_back(return_triangle[1]);
              intersection_points.push_back(return_triangle[2]);
            } else if (assign(return_vector_points, result)) {
              for (int i_point = 0; i_point < return_vector_points.size();
                   ++i_point)
                intersection_points.push_back(return_vector_points[i_point]);
            }
          } // end for loop over triangles

          // if 0 or 1 points intersect go to next face
          if (intersection_points.size() < 2)
            continue;

          // unique points
          std::vector<Point_3> unique_intersection_points;

          // insert unique points
          for (int i_point = 0; i_point < intersection_points.size();
               ++i_point) {
            bool insert_point = true;
            for (int j_point = 0; j_point < unique_intersection_points.size();
                 ++j_point)
              if (CGAL::squared_distance(intersection_points[i_point],
                                         unique_intersection_points[j_point]) <
                  0.1)
                insert_point = false;

            if (insert_point == true)
              unique_intersection_points.push_back(
                  intersection_points[i_point]);
          }

          // check number of intersected points
          if (unique_intersection_points.size() < 2)
            continue;

          // insert nodes inside cluster if not already included
          for (int i_node = 0; i_node < 3; ++i_node) {
            if (face_node_in_cluster[i_node] == true) {
              bool insert_node = true;

              for (int j_point = 0; j_point < unique_intersection_points.size();
                   ++j_point)
                if (CGAL::squared_distance(
                        face_node_locations_point[i_node],
                        unique_intersection_points[j_point]) < 0.1)
                  insert_node = false;

              if (insert_node == true)
                unique_intersection_points.push_back(
                    face_node_locations_point[i_node]);
            }
          }

          // triangles for integration
          std::vector<My_triangle> integration_triangles;
          GetTrianglesFromIntersectionPoints(unique_intersection_points,
                                             integration_triangles);

          // loop over all triangles
          for (int i_triangle = 0; i_triangle < integration_triangles.size();
               ++i_triangle) {
            // subtract planes inside
            const double temp_charge = -charge_density[i_bucket][i_face].first;

            const std::pair<std::vector<double>, double> neg_return_field =
                fieldFromDensityFaceAdaptive(temp_charge,
                                             integration_triangles[i_triangle],
                                             node_location, electric_constant,
                                             face_method, integration_error);

            // add to field
            for (int i_dof = 0; i_dof < 3; ++i_dof)
              electric_field[i_node][i_dof] += neg_return_field.first[i_dof];

            // add to potential
            electric_potential[i_node] += neg_return_field.second;
          }
        } // end for loop over faces
      }   // end for loop over buckets of faces for charge densities
    }     // end for loop over local nodes
  }       // end while loop over global nodes

  //
  return NULL;
} // end computeChargedFaceWorker
} // end of namespace for computeElectricField()

//
//
//

Electrostatics *Electrostatics::_instance = NULL;

//
// constructor
//

Electrostatics::Electrostatics() {

  //
  //
  //
  return;
}

//
//  overloading constructor
//
Electrostatics::Electrostatics(int default_electrostatics_disabled) {
  // set d_electrostaticsEnabled to zero as default
  d_electrostaticsEnabled = default_electrostatics_disabled;

  // also set the base quasi to zero
  d_baseQuasi = 0;

  return;
}

//
// destructor
//

Electrostatics::~Electrostatics() {

  //
  //
  //
  return;
}

//
// getInstance method
//

Electrostatics *Electrostatics::getInstance() {

  //
  // if not created, create
  //
  if (_instance == NULL) {
    _instance = new Electrostatics(0);
  }

  //
  // return instance
  //
  return _instance;
}

//
// destroy instance method
//

void Electrostatics::destroyInstance() {

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
//  insertInputData()
//
void Electrostatics::insertInputData(
    const std::vector<double> charges, const double atomisticRadius,
    const double elementSize, const double electricConstant,
    const int faceMethod, const std::vector<std::pair<int, double>> boundaries,
    const double integration_error, const int num_shells) {
  //
  //  this funciton is being called means electrostatics is enabled in system
  //
  d_electrostaticsEnabled = 1;

  // set charges
  d_fixedCharge.clear();
  for (int i_charge = 0; i_charge < charges.size(); ++i_charge)
    d_fixedCharge.push_back(charges[i_charge]);

  // set atomistic radius
  d_atomisticRadius = atomisticRadius;
  d_atomisticRadiusSquared = atomisticRadius * atomisticRadius;

  // set element size to switch to atomistic
  d_atomisticElementSize = elementSize;

  // set electric constant (14.39964 eVxA/e^2)
  d_electricConstant = electricConstant;

  // set plane method
  d_faceMethod = faceMethod;

  // set boundary conditions
  d_electroBoundaries.clear();
  for (int i_bound = 0; i_bound < boundaries.size(); ++i_bound)
    d_electroBoundaries.push_back(boundaries[i_bound]);

  // set integration error
  d_integration_error = integration_error;

  // set number of shells
  d_num_shells = num_shells;

  //
  return;
} // end of insertInputData()

//
//  isElectrostaticEnable()
//
int Electrostatics::isElectrostaticEnable(void) {
  // return d_electrostaticsEnabled
  return d_electrostaticsEnabled;
}

//
//  clearExternalPointCharges()
//
void Electrostatics::clearExternalPointCharges(void) {
  d_externalPointCharges.clear();

  return;
}

//
//  insertExternalPointCharges()
//
void Electrostatics::insertExternalPointCharges(
    std::vector<std::pair<double, std::vector<double>>> externalCharges) {
  d_externalPointCharges = externalCharges;

  return;
}

//
//  computeGeometries()
//
void Electrostatics::computeGeometries(void) {
  d_print("Computing geometries...");

  //
  //  get data of zeroth quasicontinuum
  //
  Quasicontinuum zeroQuasicontinuum = Quasicontinua::getInstance()->getQuasi(0);

  struct node_list_t *P_node_list =
      &(zeroQuasicontinuum.getNodeList().node_list);

  struct lattice_t *P_lattice = &(zeroQuasicontinuum.getLattice());

  struct element_list_t *P_element_list =
      &(zeroQuasicontinuum.getElementList());

  // initialize face bucket locks
  pthread_once(&bucket_once, bucket_init);

  // set number of buckets in d_faces and d_faceElementList
  d_faces.resize(numBuckets);
  d_faceElementList.resize(numBuckets);

  //
  //  find faces and set face element list
  //
  //  stores all the faces of element, with both side of face and corresponding
  //  element.
  //
  computeFaces(P_element_list);

  //
  //  set polarization sites
  //
  //  assign site close to center of tetra for each element
  //
  computePolarizationSites(P_element_list, P_lattice);

  //
  // find atomistic elements
  //
  //  loop over each element and see if it's volume is less or greater
  //  than atomistic element limit.
  //
  //  put 0 if it's less , meaning it is atomistic element
  //
  computeAtomisticElements(P_element_list);

  //
  // get atomistic element node locations
  //
  //  store initial position of nodes of all atomistic elements in data
  //
  computeAtomisticElementLocations(P_element_list);

  //
  // find atomistic nodes
  //
  //  if all elements associated to node are atomistic then declare node as
  //  atomistics
  //
  computeAtomisticNodes(P_node_list);

  //
  // find all atomistic charges
  //
  //  puts all sites inside the atomistic elements, and including
  //  atomistic nodes.
  //
  computeAtomisticCharges(P_lattice, P_element_list);

  //
  // create node hash list
  //
  ComputeBaseNodeHashList(P_node_list);

  //
  return;
} // end of computeGeometries()

//
//  computeFaces()
//
void Electrostatics::computeFaces(struct element_list_t *P_element_list) {
  //
  // instantiate computeFaces_argument_t
  //
  computeFaces_argument_t computeFaces_argument;

  //
  // find faces and set face element list
  //
  computeFaces_argument.faces = &(d_faces);
  computeFaces_argument.faceElementList = &(d_faceElementList);
  computeFaces_argument.elementList = P_element_list;

  //
  // distribute work to threads
  //
  thread_monitor(computeFacesWorker,
                 static_cast<void *>(&computeFaces_argument),
                 get_max_number_threads());

  //   //
  //   // for debugging purpose
  //   //
  //   int faceNum = 0;
  //   for(int iBucket = 0; iBucket < numBuckets; ++iBucket){
  //     for(int iFace = 0; iFace < d_faces[iBucket].size(); ++iFace){
  // std::cout<<faceNum<<" "<<iBucket<<" "<<iFace<<"
  // "<<d_faces[iBucket][iFace][0]<<" "<<d_faces[iBucket][iFace][1]<<"
  // "<<d_faces[iBucket][iFace][2]<<std::endl;
  // faceNum++;
  //     }
  //   }

  //   std::exit(0);
  //   // end of debugging block

  //
  //
  //
  return;
} // end of computeFaces()

//
//  computePolarizationSites()
//
void Electrostatics::computePolarizationSites(
    struct element_list_t *P_element_list, struct lattice_t *P_lattice) {
  //
  // set d_elementPolarizationSites to correct size
  //
  if (P_element_list->number_elements != int(d_elementPolarizationSites.size()))
    d_elementPolarizationSites.resize(P_element_list->number_elements);

  //
  // instantiate computePolarizationSites_argument_t
  //
  computePolarizationSites_argument_t computePolarizationSites_argument;

  //
  // find polarization sites list
  //
  computePolarizationSites_argument.sites = &(d_elementPolarizationSites);
  computePolarizationSites_argument.elementList = P_element_list;
  computePolarizationSites_argument.lattice = P_lattice;
  computePolarizationSites_argument.quasi = &d_baseQuasi;

  //
  // distribute work to threads
  //
  thread_monitor(computePolarizationSitesWorker,
                 static_cast<void *>(&computePolarizationSites_argument),
                 get_max_number_threads());

//
//  for debug purpose
//
#if 0
    for(int iElement = 0; iElement < P_element_list->number_elements; ++iElement){
    std::cout<<d_elementPolarizationSites[iElement][0]<<"  "<<d_elementPolarizationSites[iElement][1]<<"  "<<d_elementPolarizationSites[iElement][2]<<std::endl;
    }

    std::cout<<P_element_list->number_elements<<std::endl;
    std::exit(0);
#endif

  //
  //
  //
  return;
} // end of computePolarizationSites()

//
// compute atomistic elements
//
void Electrostatics::computeAtomisticElements(
    struct element_list_t *P_element_list) {

  //
  // set d_atomisticElements to correct size
  //
  if (P_element_list->number_elements != int(d_atomisticElements.size()))
    d_atomisticElements.resize(P_element_list->number_elements);

  //
  // instantiate computeAtomisticElements_argument_t
  //
  computeAtomisticElements_argument_t computeAtomisticElements_argument;

  //
  // find atomistic element list
  //
  computeAtomisticElements_argument.atomisticElements = &(d_atomisticElements);
  computeAtomisticElements_argument.elementList = P_element_list;
  computeAtomisticElements_argument.atomisticElementSize =
      &(d_atomisticElementSize);

  //
  // distribute work to threads
  //
  thread_monitor(computeAtomisticElementsWorker,
                 static_cast<void *>(&computeAtomisticElements_argument),
                 get_max_number_threads());

//
// for debuggin
//
#if 0
    for(int iElement = 0; iElement < P_element_list->number_elements; ++iElement){
      std::cout<<d_atomisticElements[iElement]<<std::endl;
    }
    std::exit(0);
#endif

  //
  //
  //
  return;
}

//
// compute atomistic element locations
//
void Electrostatics::computeAtomisticElementLocations(
    struct element_list_t *P_element_list) {
  //
  // set current element to 0
  //
  currentElem = 0;

  //
  // set d_atomisticElements to correct size
  //
  if (P_element_list->number_elements !=
      int(d_atomisticElementLocations.size()))
    d_atomisticElementLocations.resize(P_element_list->number_elements);

  //
  // instantiate computeAtomisticElements_argument_t
  //
  computeAtomisticElementLocations_argument_t
      computeAtomisticElementLocations_argument;

  //
  // find atomistic element list
  //
  computeAtomisticElementLocations_argument.atomisticElements =
      &(d_atomisticElements);
  computeAtomisticElementLocations_argument.atomisticElementLocations =
      &(d_atomisticElementLocations);
  computeAtomisticElementLocations_argument.elementList = P_element_list;

  //
  // distribute work to threads
  //
  thread_monitor(
      computeAtomisticElementLocationsWorker,
      static_cast<void *>(&computeAtomisticElementLocations_argument),
      get_max_number_threads());

//
//  for debugging
//
#if 0
    for(int iElement = 0; iElement < P_element_list->number_elements; ++iElement){
      std::cout<<d_atomisticElements[iElement]<<std::endl;
    }
    std::exit(0);
#endif

  //
  //
  //
  return;
}

//
// compute atomistic nodes
//
void Electrostatics::computeAtomisticNodes(struct node_list_t *P_node_list) {
  //
  // set sort variable to 0
  //
  sortVariable = 0;

  //
  // instantiate computeAtomisticNodes_argument_t
  //
  computeAtomisticNodes_argument_t computeAtomisticNodes_argument;

  //
  // find atomistic nodes list
  //
  computeAtomisticNodes_argument.atomisticElements = &(d_atomisticElements);
  computeAtomisticNodes_argument.atomisticNodes = &(d_atomisticNodes);
  computeAtomisticNodes_argument.nodeList = P_node_list;

  //
  // distribute work to threads
  //
  thread_monitor(computeAtomisticNodesWorker,
                 static_cast<void *>(&computeAtomisticNodes_argument),
                 get_max_number_threads());

//
//  fir debugging
//
#if 0
    std::cout<<d_atomisticNodes.size()<<std::endl;
    for(int iNode = 0; iNode < d_atomisticNodes.size(); ++iNode){
      std::cout<<d_atomisticNodes[iNode].first<<"  "<<d_atomisticNodes[iNode].second[0]<<"  "<<d_atomisticNodes[iNode].second[1]<<"  "<<d_atomisticNodes[iNode].second[2]<<std::endl;
    }
    std::exit(0);
#endif

  //
  //
  //
  return;
} // end of computeAtomisticNodes()

//
// compute atomistic charges
//
void Electrostatics::computeAtomisticCharges(
    struct lattice_t *P_lattice, struct element_list_t *P_element_list) {
  //
  // resize atomistic charges for each quasi
  //
  if (int(d_atomisticCharges.size()) != numBuckets)
    d_atomisticCharges.resize(numBuckets);

  //
  // add atomistic nodes to charge list
  //
  for (int iCharge = 0; iCharge < int(d_atomisticNodes.size()); ++iCharge) {

    //
    // get bucket
    //
    int bucket = getTripletIntegersKey(d_atomisticNodes[iCharge].second[0],
                                       d_atomisticNodes[iCharge].second[1],
                                       d_atomisticNodes[iCharge].second[2]);

    //
    // insert into bucket
    //
    insertChargeIntoBucket(d_atomisticNodes[iCharge].second,
                           d_atomisticCharges[bucket]);
  }

  //
  // instantiate computeAtomisticCharges_argument_t
  //
  computeAtomisticCharges_argument_t computeAtomisticCharges_argument;

  //
  // set current site
  //
  currentElem = 0;

  //
  // find atomistic charges
  //
  computeAtomisticCharges_argument.atomisticElements = &(d_atomisticElements);
  computeAtomisticCharges_argument.atomisticElementLocations =
      &(d_atomisticElementLocations);
  computeAtomisticCharges_argument.atomisticCharges = &(d_atomisticCharges);
  computeAtomisticCharges_argument.lattice = P_lattice;
  computeAtomisticCharges_argument.elementList = P_element_list;
  computeAtomisticCharges_argument.quasi = &d_baseQuasi;

  //
  // distribute work to threads
  //
  thread_monitor(computeAtomisticChargesWorker,
                 static_cast<void *>(&computeAtomisticCharges_argument),
                 get_max_number_threads());

//
//  for debugging
//
#if 0
        int totSize = 0;
        int maxSize = 0;
        for(int iBucket = 0; iBucket < numBuckets; ++iBucket){
          totSize += d_atomisticCharges[iBucket].size();
          if(maxSize < d_atomisticCharges[iBucket].size())
      maxSize = d_atomisticCharges[iBucket].size();
          for(int iCharge = 0; iCharge < d_atomisticCharges[iBucket].size(); ++iCharge){
      //std::cout<<d_atomisticCharges[iBucket][iCharge][0]<<"  "<<d_atomisticCharges[iBucket][iCharge][1]<<"  "<<d_atomisticCharges[iBucket][iCharge][2]<<std::endl;
          }
        }

        std::cout<<totSize<<std::endl;
        std::cout<<maxSize<<std::endl;
        std::exit(0);
#endif

  //
  //
  //
  return;
} //  end of computeAtomisticCharges()

//
//  ComputeBaseNodeHashList()
//
void Electrostatics::ComputeBaseNodeHashList(struct node_list_t *P_node_list) {
  // resize number of buckets
  d_base_node_hash_list.resize(numBuckets);

  // loop over all nodes in list
  for (int i_node = 0; i_node < P_node_list->number_nodes; ++i_node) {
    // get site
    const std::vector<int> site = {P_node_list->nodes[i_node]->l[0],
                                   P_node_list->nodes[i_node]->l[1],
                                   P_node_list->nodes[i_node]->l[2]};

    // get bucket
    const unsigned int i_bucket =
        getTripletIntegersKey(site[0], site[1], site[2]);

    // add node to list
    InsertNodeIntoHashList(site, i_node, d_base_node_hash_list[i_bucket]);
  } // end of for loop over nodes

  return;
} // end of ComputeBaseNodeHashList

//
// clear geometries
//
void Electrostatics::clearGeometries(void) {
  // output geometries information
  // std::cout<<"Clearing electrostatic geometries"<<std::endl;

  // clear faces
  d_faces.clear();

  // clear face element list
  d_faceElementList.clear();

  // clear polarization sites
  d_elementPolarizationSites.clear();

  // clear element polarization locations
  d_elementPolarizationLocations.clear();

  // clear atomistic charges
  d_atomisticCharges.clear();

  // clear atomistic charges
  d_atomisticElements.clear();

  // clear atomistic nodes
  d_atomisticNodes.clear();

  // clear atomistic element locations
  d_atomisticElementLocations.clear();

  // clear base element list for electric field
  d_baseElementList.clear();

  // clear base node hash list
  d_base_node_hash_list.clear();

  //
  return;
}

//
// compute electric field in domain
//
void Electrostatics::computeElectricField(const bool compute_in_reference) {
  d_print("electric field calculation...");
  std::cout.setf(std::ios::unitbuf);
  // std::cout<<"Computing electric field...";

  // get quasicontinua and other data
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();
  const int num_quasi = quasicontinua->size();

  struct node_list_t &zero_node_list =
      quasicontinua->getQuasi(0).getNodeList().node_list;
  struct lattice_t lattice = quasicontinua->getQuasi(0).getLattice();
  struct element_list_t element_list =
      quasicontinua->getQuasi(0).getElementList();

  // compute necessary quantities
  computePolarization(&element_list);
  computeChargeDensities(&zero_node_list, &lattice);
  computeAtomisticChargeState();

  // check size of base element list
  if (d_baseElementList.size() != element_list.number_elements) {
    // clear and resize base element list
    d_baseElementList.clear();
    d_baseElementList.resize(element_list.number_elements);

    // temp node list
    std::vector<int> element_nodes(4, 0);

    // loop over elements and make node to element list
    for (int i_elem = 0; i_elem < element_list.number_elements; ++i_elem) {
      // get nodes
      for (int i_node = 0; i_node < 4; ++i_node)
        element_nodes[i_node] =
            element_list.elements[i_elem]->node[i_node]->number;

      // insert nodes into element list
      d_baseElementList[i_elem].insert(d_baseElementList[i_elem].begin(),
                                       element_nodes.begin(),
                                       element_nodes.end());
    } // end loop over elements
  }   // end if statement for d_baseElementList size

  // resize field and potential for num_quasi
  if (d_electricField.size() != num_quasi)
    d_electricField.resize(num_quasi);

  if (d_electricPotential.size() != num_quasi)
    d_electricPotential.resize(num_quasi);

  if (d_electricFieldReference.size() != num_quasi)
    d_electricFieldReference.resize(num_quasi);

  if (d_electricPotentialReference.size() != num_quasi)
    d_electricPotentialReference.resize(num_quasi);

  // loop over all quasicontinuum
  for (int i_quasi = 0; i_quasi < num_quasi; ++i_quasi) {
    // get number of nodes and node list
    struct node_list_t &node_list =
        quasicontinua->getQuasi(i_quasi).getNodeList().node_list;

    const int num_nodes = node_list.number_nodes;

    // resize for number of nodes
    if (d_electricField[i_quasi].size() != num_nodes)
      d_electricField[i_quasi].resize(num_nodes);

    if (d_electricPotential[i_quasi].size() != num_nodes)
      d_electricPotential[i_quasi].resize(num_nodes);

    if (d_electricFieldReference[i_quasi].size() != num_nodes)
      d_electricFieldReference[i_quasi].resize(num_nodes);

    if (d_electricPotentialReference[i_quasi].size() != num_nodes)
      d_electricPotentialReference[i_quasi].resize(num_nodes);

    // set values to 0
    for (int i_node = 0; i_node < num_nodes; ++i_node) {
      if (compute_in_reference == true) {
        d_electricFieldReference[i_quasi][i_node].clear();

        for (int i_dof = 0; i_dof < 4; ++i_dof)
          d_electricFieldReference[i_quasi][i_node].push_back(0.0);

        d_electricPotentialReference[i_quasi][i_node] = 0.0;
      } else {
        d_electricField[i_quasi][i_node].clear();

        for (int i_dof = 0; i_dof < 4; ++i_dof)
          d_electricField[i_quasi][i_node].push_back(0.0);

        d_electricPotential[i_quasi][i_node] = 0.0;
      }
    } // end of for loop over nodes
  }   // end of for loop over i_quasi

  // compute contribution from atomistic atoms
  ComputeContributionFromAtomistic(compute_in_reference);

  // compute contribution from atoms in cluster
  ComputeContributionFromClusters(compute_in_reference);

  // compute contribution from external charges
  ComputeContributionFromExternal(compute_in_reference);

  //
  ComputeContributionFromClusterSurface(compute_in_reference);

  //
  ComputeContributionFromChargedFaces(compute_in_reference);

  // std::cout<<"Completed."<<std::endl;
  return;
} // end of computeElectricField

//
// calculate polarization in domain
//
void Electrostatics::computePolarization(
    struct element_list_t *P_element_list) {
  //
  // check size of d_elementPolarization
  //
  if (d_atomisticElements.size() != d_elementPolarization.size())
    d_elementPolarization.resize(d_atomisticElements.size());

  //
  // check size of element polarization locations
  //
  if (d_atomisticElements.size() != d_elementPolarizationLocations.size())
    d_elementPolarizationLocations.resize(d_atomisticElements.size());

  //
  // compute element polarization locations
  //
  Electrostatics::computeElementPolarizationLocations();

  //
  // compute element polarizations
  //
  Electrostatics::computeElementPolarization(P_element_list);

  //
  //
  //
  return;
}

//
// compute charge densities
//
void Electrostatics::computeChargeDensities(struct node_list_t *P_node_list,
                                            struct lattice_t *P_lattice) {
  //
  // resize charge density
  //
  if (int(d_chargeDensity.size()) != numBuckets)
    d_chargeDensity.resize(numBuckets);

  //
  // set current bucket to 0
  //
  currentBucket = 0;

  //
  // instantiate computeChargeDensities_argument_t
  //
  computeChargeDensities_argument_t computeChargeDensities_argument;

  //
  // find faces and set face element list
  //
  computeChargeDensities_argument.faces = &(d_faces);
  computeChargeDensities_argument.faceElementList = &(d_faceElementList);
  computeChargeDensities_argument.nodeList = P_node_list;
  computeChargeDensities_argument.chargeDensity = &(d_chargeDensity);
  computeChargeDensities_argument.elementPolarization =
      &(d_elementPolarization);
  computeChargeDensities_argument.lattice = P_lattice;
  computeChargeDensities_argument.electroBoundaries = &(d_electroBoundaries);

  //
  // distribute work to threads
  //
  thread_monitor(computeChargeDensitiesWorker,
                 static_cast<void *>(&computeChargeDensities_argument),
                 get_max_number_threads());

  //
  //
  //
  return;
}

//
// compute atomistic charge locations
//
void Electrostatics::computeAtomisticChargeState() {
  //
  // get Quasicontinua
  //
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  //
  // get number of Quasicontinuums
  //
  const int numQuasi = quasicontinua->size();

  //
  // resize atomsisticChargeLocations
  //
  if (int(d_atomisticChargeState.size()) != numQuasi)
    d_atomisticChargeState.resize(numQuasi);

  //
  // loop over all quasicontinuum and add
  //
  for (int iQuasi = 0; iQuasi < numQuasi; ++iQuasi) {
    //
    // set current quasi
    //
    struct lattice_t lattice = quasicontinua->getQuasi(iQuasi).getLattice();

    //
    // clear charge locations
    //
    d_atomisticChargeState[iQuasi].clear();

    //
    // set current element
    //
    currentBucket = 0;

    //
    // instantiate computeAtomisticChargeLocations_argument_t
    //
    computeAtomisticChargeState_argument_t computeAtomisticChargeState_argument;

    //
    // compute polarization locations argument
    //
    computeAtomisticChargeState_argument.atomisticCharges =
        &(d_atomisticCharges);
    computeAtomisticChargeState_argument.atomisticChargeState =
        &(d_atomisticChargeState);
    computeAtomisticChargeState_argument.quasiId = iQuasi;
    computeAtomisticChargeState_argument.lattice = &lattice;

    //
    // distribute work to threads
    //
    thread_monitor(computeAtomisticChargeStateWorker,
                   static_cast<void *>(&computeAtomisticChargeState_argument),
                   get_max_number_threads());
  }

  //
  // add any atoms missed from void
  //
  const std::vector<std::pair<int, std::vector<int>>> addAtoms;

  //
  // loop over all atoms to add
  //
  for (int iAtom = 0; iAtom < addAtoms.size(); ++iAtom) {
    //
    // get quasi
    //
    const int iQuasi = addAtoms[iAtom].first;

    //
    // get lattice
    //
    struct lattice_t lattice = quasicontinua->getQuasi(iQuasi).getLattice();

    //
    // temporary location
    //
    double state[4];

    //
    // temporary site
    //
    int site[3];
    site[0] = addAtoms[iQuasi].second[0];
    site[1] = addAtoms[iQuasi].second[1];
    site[2] = addAtoms[iQuasi].second[2];

    //
    // get site location
    //
    Lattice::getInstance()->getSiteCurrentState(state, site, &lattice, iQuasi);

    //
    // add location to charge list
    //
    std::vector<double> vecState = {state[0], state[1], state[2], state[3]};
    d_atomisticChargeState[iQuasi].push_back(vecState);
  }

  //
  //
  //
  return;
}

//
// compute contribution to field and potential from atomistic atoms
//
void Electrostatics::ComputeContributionFromAtomistic(
    const bool &compute_in_reference) {
  // get number of quasicontinuums
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();
  const int num_quasi = quasicontinua->size();

  // loop over all quasicontinuum
  for (int i_quasi = 0; i_quasi < num_quasi; ++i_quasi) {
    // get number of nodes and node list
    struct node_list_t &node_list =
        quasicontinua->getQuasi(i_quasi).getNodeList().node_list;
    const int num_nodes = node_list.number_nodes;
    std::vector<double> shift = Quasicontinua::getInstance()->getShift(i_quasi);

    // set current node to 0
    currentNode = 0;

    // instantiate computeFieldFromAtomisticCharges_argument_t
    computeFieldFromAtomisticCharges_argument_t
        computeFieldFromAtomisticCharges_argument;

    // set data structures
    computeFieldFromAtomisticCharges_argument.nodeList = &node_list;

    if (compute_in_reference == true) {
      computeFieldFromAtomisticCharges_argument.electricField =
          &(d_electricFieldReference[i_quasi]);
      computeFieldFromAtomisticCharges_argument.electricPotential =
          &(d_electricPotentialReference[i_quasi]);
    } else {
      computeFieldFromAtomisticCharges_argument.electricField =
          &(d_electricField[i_quasi]);
      computeFieldFromAtomisticCharges_argument.electricPotential =
          &(d_electricPotential[i_quasi]);
    }
    computeFieldFromAtomisticCharges_argument.fixedCharge = &(d_fixedCharge);
    computeFieldFromAtomisticCharges_argument.electricConstant =
        &(d_electricConstant);
    computeFieldFromAtomisticCharges_argument.atomisticChargeState =
        &(d_atomisticChargeState);
    computeFieldFromAtomisticCharges_argument.iShift = &shift;
    computeFieldFromAtomisticCharges_argument.iQuasi = &i_quasi;

    // distribute work to threads
    thread_monitor(
        computeFieldFromAtomisticChargesWorker,
        static_cast<void *>(&computeFieldFromAtomisticCharges_argument),
        get_max_number_threads());
  } // end for loop over i_quasi

  return;
} // end of ComputeContributionFromAtomistic

//
// compute contribution to field and potential from clusters
//
void Electrostatics::ComputeContributionFromClusters(
    const bool &compute_in_reference) {
  // get number of quasicontinuums
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();
  const int num_quasi = quasicontinua->size();

  // loop over all quasicontinuum
  for (int i_quasi = 0; i_quasi < num_quasi; ++i_quasi) {
    // get number of nodes, node list, and shift
    struct node_list_t &node_list =
        quasicontinua->getQuasi(i_quasi).getNodeList().node_list;
    std::vector<double> i_shift =
        Quasicontinua::getInstance()->getShift(i_quasi);

    int iCoreFlag = quasicontinua->getCoreShell(i_quasi);

    // loop over all quasicontinuum
    for (int j_quasi = 0; j_quasi < num_quasi; ++j_quasi) {
      // get j_quasi lattice and shift
      struct lattice_t j_lattice =
          quasicontinua->getQuasi(j_quasi).getLattice();

      std::vector<double> j_shift =
          Quasicontinua::getInstance()->getShift(j_quasi);

      int jCoreFlag = quasicontinua->getCoreShell(j_quasi);

      std::pair<int, int> coreData;
      coreData.first = iCoreFlag;
      coreData.second = jCoreFlag;

      double j_charge = d_fixedCharge[j_quasi];

      // set current node to 0
      currentNode = 0;

      // instantiate computeClusterCharges_argument_t
      computeClusterCharges_argument_t computeClusterCharges_argument;

      // set data structures
      computeClusterCharges_argument.i_node_list = &node_list;

      if (compute_in_reference == true) {
        computeClusterCharges_argument.electricField =
            &(d_electricFieldReference[i_quasi]);
        computeClusterCharges_argument.electricPotential =
            &(d_electricPotentialReference[i_quasi]);
      } else {
        computeClusterCharges_argument.electricField =
            &(d_electricField[i_quasi]);
        computeClusterCharges_argument.electricPotential =
            &(d_electricPotential[i_quasi]);
      }

      computeClusterCharges_argument.electricConstant = &(d_electricConstant);
      computeClusterCharges_argument.atomisticCharges = &(d_atomisticCharges);
      computeClusterCharges_argument.j_lattice = &j_lattice;
      computeClusterCharges_argument.i_shift = &i_shift;
      computeClusterCharges_argument.j_shift = &j_shift;
      computeClusterCharges_argument.num_shells = &d_num_shells;
      computeClusterCharges_argument.j_charge = &j_charge;
      computeClusterCharges_argument.j_quasi = &j_quasi;
      computeClusterCharges_argument.i_quasi = &i_quasi;
      computeClusterCharges_argument.coreData = &coreData;

      // distribute work to threads
      thread_monitor(computeClusterChargesWorker,
                     static_cast<void *>(&computeClusterCharges_argument),
                     get_max_number_threads());
    } // end of for loop over j_quasi
  }   // end of for loop over i_quasi

  return;
} // end of ComputeContributionFromClusters

//
// compute contribution to field and potential from external loads
//
void Electrostatics::ComputeContributionFromExternal(
    const bool &compute_in_reference) {
  // if in reference do not compute, can be changed if wanted
  if (compute_in_reference == true)
    return;

  // get number of quasicontinuums
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();
  const int num_quasi = quasicontinua->size();

  // loop over all quasicontinuum
  for (int i_quasi = 0; i_quasi < num_quasi; ++i_quasi) {
    // get number of nodes and node list
    struct node_list_t &node_list =
        quasicontinua->getQuasi(i_quasi).getNodeList().node_list;
    const int num_nodes = node_list.number_nodes;
    std::vector<double> shift = Quasicontinua::getInstance()->getShift(i_quasi);

    int iCoreFlag = quasicontinua->getCoreShell(i_quasi);

    // set current node to 0
    currentNode = 0;

    // instantiate compute_field_from_external_charges_argument_t
    compute_field_from_external_charges_argument_t
        compute_field_from_external_charges_argument;

    // set data structures
    compute_field_from_external_charges_argument.node_list = &node_list;
    compute_field_from_external_charges_argument.electric_field =
        &(d_electricField[i_quasi]);
    compute_field_from_external_charges_argument.electric_potential =
        &(d_electricPotential[i_quasi]);
    compute_field_from_external_charges_argument.electric_constant =
        &(d_electricConstant);
    compute_field_from_external_charges_argument.external_charges =
        &(d_externalPointCharges);
    compute_field_from_external_charges_argument.iShift = &shift;
    compute_field_from_external_charges_argument.iQuasi = &i_quasi;
    compute_field_from_external_charges_argument.iCoreFlag = &iCoreFlag;

    // distribute work to threads
    thread_monitor(
        ComputeFieldFromExternalChargesWorker,
        static_cast<void *>(&compute_field_from_external_charges_argument),
        get_max_number_threads());
  } // end for loop over i_quasi

  return;
} // end of ComputeContributionFromExternal

//
// compute contribution to field and potential from cluster surface
//
void Electrostatics::ComputeContributionFromClusterSurface(
    const bool &compute_in_reference) {
  // get number of quasicontinuums
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();
  const int num_quasi = quasicontinua->size();

  // loop over all quasicontinuum
  for (int i_quasi = 0; i_quasi < num_quasi; ++i_quasi) {
    // get number of nodes, node list, and shift
    struct node_list_t &node_list =
        quasicontinua->getQuasi(i_quasi).getNodeList().node_list;
    std::vector<double> i_shift =
        Quasicontinua::getInstance()->getShift(i_quasi);

    // set j_quasi to 0
    const int j_quasi = 0;

    // get j_quasi lattice and shift
    struct lattice_t j_lattice = quasicontinua->getQuasi(j_quasi).getLattice();
    std::vector<double> j_shift =
        Quasicontinua::getInstance()->getShift(j_quasi);
    struct node_list_t &zero_node_list =
        quasicontinua->getQuasi(j_quasi).getNodeList().node_list;

    // set current node to 0
    currentNode = 0;

    // instantiate computeClusterChargeSurface_argument_t
    computeClusterChargeSurface_argument_t computeClusterChargeSurface_argument;

    // set data structures
    computeClusterChargeSurface_argument.i_node_list = &node_list;

    if (compute_in_reference == true) {
      computeClusterChargeSurface_argument.electricField =
          &(d_electricFieldReference[i_quasi]);
      computeClusterChargeSurface_argument.electricPotential =
          &(d_electricPotentialReference[i_quasi]);
    } else {
      computeClusterChargeSurface_argument.electricField =
          &(d_electricField[i_quasi]);
      computeClusterChargeSurface_argument.electricPotential =
          &(d_electricPotential[i_quasi]);
    }

    computeClusterChargeSurface_argument.electricConstant =
        &(d_electricConstant);
    computeClusterChargeSurface_argument.atomisticCharges =
        &(d_atomisticCharges);
    computeClusterChargeSurface_argument.element_polarization =
        &(d_elementPolarization);
    computeClusterChargeSurface_argument.zero_lattice = &j_lattice;
    computeClusterChargeSurface_argument.i_shift = &i_shift;
    computeClusterChargeSurface_argument.zero_shift = &j_shift;
    computeClusterChargeSurface_argument.num_shells = &d_num_shells;
    computeClusterChargeSurface_argument.face_method = &d_faceMethod;
    computeClusterChargeSurface_argument.integration_error =
        &d_integration_error;
    computeClusterChargeSurface_argument.zero_node_list = &zero_node_list;
    computeClusterChargeSurface_argument.base_node_hash_list =
        &d_base_node_hash_list;
    computeClusterChargeSurface_argument.zero_quasi = &d_baseQuasi;

    // distribute work to threads
    thread_monitor(computeClusterChargeSurfaceWorker,
                   static_cast<void *>(&computeClusterChargeSurface_argument),
                   get_max_number_threads());
  } // end of for loop over i_quasi

  return;
} // end of ComputeContributionFromClusterSurface

//
// compute contribution to field and potential from cluster surface
//
void Electrostatics::ComputeContributionFromChargedFaces(
    const bool &compute_in_reference) {
  // get number of quasicontinuums
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();
  const int num_quasi = quasicontinua->size();

  // loop over all quasicontinuum
  for (int i_quasi = 0; i_quasi < num_quasi; ++i_quasi) {
    // get number of nodes, node list, and shift
    struct node_list_t &node_list =
        quasicontinua->getQuasi(i_quasi).getNodeList().node_list;
    std::vector<double> i_shift =
        Quasicontinua::getInstance()->getShift(i_quasi);

    // set j_quasi to 0
    const int j_quasi = 0;

    // get j_quasi lattice and shift
    struct lattice_t j_lattice = quasicontinua->getQuasi(j_quasi).getLattice();
    std::vector<double> j_shift =
        Quasicontinua::getInstance()->getShift(j_quasi);
    struct node_list_t &zero_node_list =
        quasicontinua->getQuasi(j_quasi).getNodeList().node_list;

    // set current node to 0
    currentNode = 0;

    // instantiate computeChargedFace_argument_t
    computeChargedFace_argument_t computeChargedFace_argument;

    // set data structures
    computeChargedFace_argument.i_node_list = &node_list;

    if (compute_in_reference == true) {
      computeChargedFace_argument.electricField =
          &(d_electricFieldReference[i_quasi]);
      computeChargedFace_argument.electricPotential =
          &(d_electricPotentialReference[i_quasi]);
    } else {
      computeChargedFace_argument.electricField = &(d_electricField[i_quasi]);
      computeChargedFace_argument.electricPotential =
          &(d_electricPotential[i_quasi]);
    }

    computeChargedFace_argument.electricConstant = &(d_electricConstant);
    computeChargedFace_argument.i_shift = &i_shift;
    computeChargedFace_argument.zero_shift = &j_shift;
    computeChargedFace_argument.num_shells = &d_num_shells;
    computeChargedFace_argument.face_method = &d_faceMethod;
    computeChargedFace_argument.integration_error = &d_integration_error;
    computeChargedFace_argument.zero_node_list = &zero_node_list;
    computeChargedFace_argument.faces = &d_faces;
    computeChargedFace_argument.charge_density = &d_chargeDensity;
    computeChargedFace_argument.zero_lattice = &j_lattice;
    computeChargedFace_argument.atomistic_charges = &d_atomisticCharges;
    computeChargedFace_argument.zero_quasi = &d_baseQuasi;

    // distribute work to threads
    thread_monitor(computeChargedFaceWorker,
                   static_cast<void *>(&computeChargedFace_argument),
                   get_max_number_threads());
  } // end of for loop over i_quasi

  return;
} // end of ComputeContributionFromChargedFaces

//
// compute polarization locations
//
void Electrostatics::computeElementPolarizationLocations() {
  //
  // get Quasicontinua
  //
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  //
  // get number of Quasicontinuums
  //
  const int numQuasi = quasicontinua->size();

  //
  // loop over all quasicontinuum and add
  //
  for (int iQuasi = 0; iQuasi < numQuasi; ++iQuasi) {
    //
    // set current quasi
    //
    struct lattice_t lattice = quasicontinua->getQuasi(iQuasi).getLattice();

    //
    // set current element
    //
    currentElem = 0;

    //
    // instantiate computePolarizationLocations_argument_t
    //
    computePolarizationLocations_argument_t
        computePolarizationLocations_argument;

    //
    // compute polarization locations argument
    //
    computePolarizationLocations_argument.atomisticElements =
        &(d_atomisticElements);
    computePolarizationLocations_argument.elementPolarizationLocations =
        &(d_elementPolarizationLocations);
    computePolarizationLocations_argument.elementPolarizationSites =
        &(d_elementPolarizationSites);
    computePolarizationLocations_argument.quasiId = iQuasi;
    computePolarizationLocations_argument.fixedCharge = d_fixedCharge[iQuasi];
    computePolarizationLocations_argument.lattice = &lattice;

    //
    // distribute work to threads
    //
    thread_monitor(computePolarizationLocationsWorker,
                   static_cast<void *>(&computePolarizationLocations_argument),
                   get_max_number_threads());
  }

  return;
} // end of computeElementPolarizationLocations()

//
// compute polarization in elements
//
void Electrostatics::computeElementPolarization(
    struct element_list_t *P_element_list) {
  //
  // get Quasicontinua
  //
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  //
  // set current quasi
  //
  struct lattice_t lattice = quasicontinua->getQuasi(0).getLattice();

  //
  // initial unit cell volume
  //
  double unitCellVolume =
      (lattice.a1[1] * lattice.a2[2] - lattice.a1[2] * lattice.a2[1]) *
          lattice.a3[0] +
      (lattice.a1[2] * lattice.a2[0] - lattice.a1[0] * lattice.a2[2]) *
          lattice.a3[1] +
      (lattice.a1[0] * lattice.a2[1] - lattice.a1[1] * lattice.a2[0]) *
          lattice.a3[2];

  //
  // make sure volume is positive
  //
  if (unitCellVolume < 0.0)
    unitCellVolume = unitCellVolume * -1.0;

  //
  // set current element
  //
  currentElem = 0;

  //
  // instantiate computePolarizationLocations_argument_t
  //
  computeElementPolarization_argument_t computeElementPolarization_argument;

  //
  // compute polarization locations argument
  //
  computeElementPolarization_argument.atomisticElements =
      &(d_atomisticElements);
  computeElementPolarization_argument.elementPolarizationLocations =
      &(d_elementPolarizationLocations);
  computeElementPolarization_argument.elementPolarization =
      &(d_elementPolarization);
  computeElementPolarization_argument.elementList = P_element_list;
  computeElementPolarization_argument.initialUnitCellVolume = unitCellVolume;

  //
  // distribute work to threads
  //
  thread_monitor(computeElementPolarizationWorker,
                 static_cast<void *>(&computeElementPolarization_argument),
                 get_max_number_threads());

  //
  //
  //
  return;
}

//
// compute electric field in Reference domain
//
void Electrostatics::computeElectricFieldReference() {
  // get quasi holder
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  // get number of quasi
  const int numQuasi = quasicontinua->size();

  // node positions to hold current locations
  std::vector<std::vector<std::vector<double>>> nodePositions;

  // resize for numQuasi
  nodePositions.resize(numQuasi);

  // loop over all quasicontinuum except 0
  for (int iQuasi = 1; iQuasi < numQuasi; ++iQuasi) {
    // get nodeList
    struct node_list_t &quasiNodeList =
        quasicontinua->getQuasi(iQuasi).getNodeList().node_list;

    // get number of nodes in iQuasi
    int numberNodes = quasiNodeList.number_nodes;

    // resize nodePositions for number of nodes
    nodePositions[iQuasi].resize(numberNodes);

    // loop over nodes in each quasicontinuum
    for (int iNode = 0; iNode < numberNodes; ++iNode) {
      // save current position
      nodePositions[iQuasi][iNode].clear();
      for (int iDof = 0; iDof < 3; ++iDof)
        nodePositions[iQuasi][iNode].push_back(
            quasiNodeList.nodes[iNode]->position[iDof]);

      // set current position to initial position
      for (int iDof = 0; iDof < 3; ++iDof)
        quasiNodeList.nodes[iNode]->position[iDof] =
            quasiNodeList.nodes[iNode]->initial_position[iDof];
    } // end loop over nodes
  }   // end loop over quasi

  // save external forces and zero out current
  const std::vector<std::pair<double, std::vector<double>>>
      saved_external_point_charges{d_externalPointCharges};
  d_externalPointCharges.clear();

  // reset current quasicontinuum to 0
  Quasicontinuum quasicontinuum = quasicontinua->getQuasi(0);

  // compute electric field
  computeElectricField(true);

  // reset external charges
  for (int i_charge = 0; i_charge < saved_external_point_charges.size();
       ++i_charge)
    d_externalPointCharges.push_back(saved_external_point_charges[i_charge]);

  // reset saved current positions for all quasicontinuum except 0
  for (int iQuasi = 1; iQuasi < numQuasi; ++iQuasi) {
    // get nodeList
    struct node_list_t &quasiNodeList =
        quasicontinua->getQuasi(iQuasi).getNodeList().node_list;

    // get number of nodes in iQuasi
    int numberNodes = quasiNodeList.number_nodes;

    // loop over nodes in each quasicontinuum
    for (int iNode = 0; iNode < numberNodes; ++iNode) {
      // save current position
      for (int iDof = 0; iDof < 3; ++iDof)
        quasiNodeList.nodes[iNode]->position[iDof] =
            nodePositions[iQuasi][iNode][iDof];
    } // end loop over nodes
  }   // end loop over quasi

  //
  return;
} // end of Electrostatics::computeElectricFieldReference

//
// get electric field
//
const std::vector<std::vector<double>> &
Electrostatics::getElectricField(const int i_quasi) const {
  return d_electricField[i_quasi];
}

//
//  computeForceOnNodeDueToSite()
//
void Electrostatics::computeForceOnNodeDueToSite(
    int iQuasi, cluster_site_data_t iCSite_data, bool compute_in_reference) {
  //
  //  interpolate the field and potential at iCSite and get force and energy
  //
  std::vector<double> f_iC;
  for (int dof = 0; dof < 4; dof++)
    f_iC.push_back(0.0);

  double energy_iC = 0.0;

  // get Quasicontinua instance
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  // get Shape instance
  Shape *shapeC = Shape::getInstance();

  // interpolate appropriate values to site
  if (compute_in_reference == true) {
    InterpolateValuesAtPoint(iCSite_data, iQuasi,
                             d_electricFieldReference[iQuasi],
                             d_electricPotentialReference[iQuasi],
                             d_fixedCharge[iQuasi], f_iC, energy_iC);
  } else {
    InterpolateValuesAtPoint(iCSite_data, iQuasi, d_electricField[iQuasi],
                             d_electricPotential[iQuasi], d_fixedCharge[iQuasi],
                             f_iC, energy_iC);
  }

  //
  //  distribute the computed force on site to relevant nodes
  //
  struct all_node_list_t &iQuasi_all_node_list =
      quasicontinua->getQuasi(iQuasi).getNodeList();

  struct lattice_t iQuasi_lattice =
      quasicontinua->getQuasi(iQuasi).getLattice();

  struct node_list_t iQuasi_node_list = iQuasi_all_node_list.node_list;

  struct node_t *P_node_0 = iQuasi_node_list.nodes[iCSite_data.second[0]];
  struct node_t *P_node_1 = iQuasi_node_list.nodes[iCSite_data.second[1]];
  struct node_t *P_node_2 = iQuasi_node_list.nodes[iCSite_data.second[2]];
  struct node_t *P_node_3 = iQuasi_node_list.nodes[iCSite_data.second[3]];

  //
  //  cluster weight
  //
  double nodeWeight = P_node_0->weight;

  int l_iC[3];
  for (int dof = 0; dof < 3; dof++)
    l_iC[dof] = iCSite_data.first.first[dof];

  //
  //  process node_0
  //
  double shape = shapeC->computeSiteShapeFunction(
      P_node_0, P_node_1, P_node_2, P_node_3, l_iC, &iQuasi_lattice, iQuasi);

  //  add contribution to node 0
  //  Note : for this node, we will also add contribution to energy
  pthread_mutex_lock(&(P_node_0->lock));
  P_node_0->acceleration[0] += nodeWeight * f_iC[0] * shape;
  P_node_0->acceleration[1] += nodeWeight * f_iC[1] * shape;
  P_node_0->acceleration[2] += nodeWeight * f_iC[2] * shape;

  P_node_0->force_frequency += nodeWeight * f_iC[3] * shape;

  // update the energy
  P_node_0->energy.potential += nodeWeight * energy_iC;
  pthread_mutex_unlock(&(P_node_0->lock));

  //
  //  process node_1
  //
  shape = shapeC->computeSiteShapeFunction(
      P_node_1, P_node_0, P_node_2, P_node_3, l_iC, &iQuasi_lattice, iQuasi);

  //  add contribution to node 1
  pthread_mutex_lock(&(P_node_1->lock));
  P_node_1->acceleration[0] += nodeWeight * f_iC[0] * shape;
  P_node_1->acceleration[1] += nodeWeight * f_iC[1] * shape;
  P_node_1->acceleration[2] += nodeWeight * f_iC[2] * shape;

  P_node_1->force_frequency += nodeWeight * f_iC[3] * shape;
  pthread_mutex_unlock(&(P_node_1->lock));

  //
  //  process node_2
  //
  shape = shapeC->computeSiteShapeFunction(
      P_node_2, P_node_0, P_node_1, P_node_3, l_iC, &iQuasi_lattice, iQuasi);

  //  add contribution to node 2
  pthread_mutex_lock(&(P_node_2->lock));
  P_node_2->acceleration[0] += nodeWeight * f_iC[0] * shape;
  P_node_2->acceleration[1] += nodeWeight * f_iC[1] * shape;
  P_node_2->acceleration[2] += nodeWeight * f_iC[2] * shape;

  P_node_2->force_frequency += nodeWeight * f_iC[3] * shape;
  pthread_mutex_unlock(&(P_node_2->lock));

  //
  //  process node_3
  //
  shape = shapeC->computeSiteShapeFunction(
      P_node_3, P_node_0, P_node_1, P_node_2, l_iC, &iQuasi_lattice, iQuasi);

  //  add contribution to node 3
  pthread_mutex_lock(&(P_node_3->lock));
  P_node_3->acceleration[0] += nodeWeight * f_iC[0] * shape;
  P_node_3->acceleration[1] += nodeWeight * f_iC[1] * shape;
  P_node_3->acceleration[2] += nodeWeight * f_iC[2] * shape;

  P_node_3->force_frequency += nodeWeight * f_iC[3] * shape;
  pthread_mutex_unlock(&(P_node_3->lock));

  return;
}

//
//  computeForceOnSite()
//
void Electrostatics::computeForceOnSite(int iQuasi,
                                        cluster_site_data_t iCSite_data,
                                        std::vector<double> &force,
                                        double &energy,
                                        bool compute_in_reference) {
  //
  //  interpolate the field and potential at iCSite and get force and energy
  //

  // interpolate appropriate values to site
  if (compute_in_reference == true) {
    InterpolateValuesAtPoint(iCSite_data, iQuasi,
                             d_electricFieldReference[iQuasi],
                             d_electricPotentialReference[iQuasi],
                             d_fixedCharge[iQuasi], force, energy);
  } else {
    InterpolateValuesAtPoint(iCSite_data, iQuasi, d_electricField[iQuasi],
                             d_electricPotential[iQuasi], d_fixedCharge[iQuasi],
                             force, energy);
  }

  return;
}
}