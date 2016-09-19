//
// PairPotentials.cc
//

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

#ifdef STDC_HEADERS
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#else
#error No standard C library headers found
#endif /* STDC_HEADERS */

#ifdef HAVE_ASSERT_H
#include <assert.h>
#else
#error assert.h not found.
#endif /* HAVE_ASSERT_H */

#include <cmath>
#include <utility>

#include "DataTypes.h"
#include "Error.h"
#include "PairPotentials.h"

#include "spline.h" // for fortran subroutine

// gloabl variable
#define CUTOFF_NEIGHBOR_SCALE 1.1
#define CUTOFF_CLUSTER_SCALE 1.0

// universal physics constants
long double const c_avagradoNumber = 6.022140857e23;
long double const c_boltzmanConstant = 0.831445986;
long double const c_maxPlanckConstant = 39.903127109;
long double const c_maxPlanckPi = 6.350779924;
long double const c_dielectricConstantVacuum = 0.0055263488697;
long double const c_dielectricConstantVacuumWithPi = 0.069446148;

//
//
//

namespace quasicontinuum {

//
//
//

PairPotentials *PairPotentials::_instance = NULL;

//
// constructor
//

PairPotentials::PairPotentials() {
  //
  // put default values to universal constants
  //
  c_universalConst.clear();

  c_universalConst.push_back(6.022140857e23);  // avagrado number
  c_universalConst.push_back(0.831445986);     // boltzman const
  c_universalConst.push_back(39.903127109);    // max plack
  c_universalConst.push_back(6.350779924);     // max planck over 2pi
  c_universalConst.push_back(0.0055263488697); // dielectri const
  c_universalConst.push_back(0.069446148);     // dielectric * 4pi

  //
  //
  //
  return;
}

//
// destructor
//

PairPotentials::~PairPotentials() {

  //
  //
  //
  return;
}

//
// getInstance method
//

PairPotentials *PairPotentials::getInstance() {

  //
  // if not created, create
  //
  if (_instance == NULL) {
    _instance = new PairPotentials();
  }

  //
  // return instance
  //
  return _instance;
}

//
// destroy instance method
//

void PairPotentials::destroyInstance() {

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
// insert Potential interaction
//
void PairPotentials::insertNonEAMPotential(
    std::pair<int, int> quasicontinuumPair, int potentialType,
    std::vector<double> potentialParameters, const bool EAMFlag) {

  //
  // insert Potential interaction
  //
  d_PotentialInteractions.push_back(quasicontinuumPair);

  //
  // insert Potential type
  //
  d_PotentialType.push_back(potentialType);

  //
  // insert vector of Potential parameters
  //
  d_PotentialParameters.push_back(potentialParameters);

  //
  // set EAMFlag
  //
  d_EAMFlag = EAMFlag;

  //
  //
  //
  return;
}

//
// insert Potential interaction
//

void PairPotentials::insertEAMPotential(std::pair<int, int> quasicontinuumPair,
                                        int potentialType,
                                        std::vector<double> potentialParameters,
                                        const char *dataFile1,
                                        const int NGridPoints,
                                        const bool EAMFlag) {

  //
  // insert Potential interaction
  //
  d_PotentialInteractions.push_back(quasicontinuumPair);

  //
  // insert Potential type
  //
  d_PotentialType.push_back(potentialType);

  //
  // insert vector of Potential parameters (first parameter must be R_c)
  //
  d_PotentialParameters.push_back(potentialParameters);

  //
  // set EAMFlag
  //
  d_EAMFlag = EAMFlag;

  //
  // get current potential number
  //
  const unsigned int potNumber = d_PotentialType.size();

  //
  // resize data
  //
  d_RhoX.resize(potNumber);
  d_RhoY.resize(potNumber);
  d_RhoB.resize(potNumber);
  d_RhoC.resize(potNumber);
  d_RhoD.resize(potNumber);
  d_RhoX[potNumber - 1].resize(NGridPoints);
  d_RhoY[potNumber - 1].resize(NGridPoints);
  d_RhoB[potNumber - 1].resize(NGridPoints);
  d_RhoC[potNumber - 1].resize(NGridPoints);
  d_RhoD[potNumber - 1].resize(NGridPoints);

  d_EmbedX.resize(potNumber);
  d_EmbedY.resize(potNumber);
  d_EmbedB.resize(potNumber);
  d_EmbedC.resize(potNumber);
  d_EmbedD.resize(potNumber);
  d_EmbedX[potNumber - 1].resize(NGridPoints);
  d_EmbedY[potNumber - 1].resize(NGridPoints);
  d_EmbedB[potNumber - 1].resize(NGridPoints);
  d_EmbedC[potNumber - 1].resize(NGridPoints);
  d_EmbedD[potNumber - 1].resize(NGridPoints);

  d_PairX.resize(potNumber);
  d_PairY.resize(potNumber);
  d_PairB.resize(potNumber);
  d_PairC.resize(potNumber);
  d_PairD.resize(potNumber);
  d_PairX[potNumber - 1].resize(NGridPoints);
  d_PairY[potNumber - 1].resize(NGridPoints);
  d_PairB[potNumber - 1].resize(NGridPoints);
  d_PairC[potNumber - 1].resize(NGridPoints);
  d_PairD[potNumber - 1].resize(NGridPoints);

  //
  // read in data file
  //
  FILE *fp;
  fp = fopen(dataFile1, "r");
  assert(fp != NULL);

  //
  // read in rho
  //
  for (int iPoint = 0; iPoint < NGridPoints; ++iPoint)
    fscanf(fp, "%le %le %le %le %le", &(d_RhoX[potNumber - 1][iPoint]),
           &(d_RhoY[potNumber - 1][iPoint]), &(d_RhoB[potNumber - 1][iPoint]),
           &(d_RhoC[potNumber - 1][iPoint]), &(d_RhoD[potNumber - 1][iPoint]));

  for (int iPoint = 0; iPoint < NGridPoints; ++iPoint)
    fscanf(
        fp, "%le %le %le %le %le", &(d_EmbedX[potNumber - 1][iPoint]),
        &(d_EmbedY[potNumber - 1][iPoint]), &(d_EmbedB[potNumber - 1][iPoint]),
        &(d_EmbedC[potNumber - 1][iPoint]), &(d_EmbedD[potNumber - 1][iPoint]));

  for (int iPoint = 0; iPoint < NGridPoints; ++iPoint)
    fscanf(fp, "%le %le %le %le %le", &(d_PairX[potNumber - 1][iPoint]),
           &(d_PairY[potNumber - 1][iPoint]), &(d_PairB[potNumber - 1][iPoint]),
           &(d_PairC[potNumber - 1][iPoint]),
           &(d_PairD[potNumber - 1][iPoint]));

  //
  // close file
  //
  fclose(fp);

  //
  //
  //
  return;
}

//
// get forces and energies
//
std::pair<double, double>
PairPotentials::getPairForceAndEnergy(int potentialNumber, double separation) {
  // // d_print("potential number = %d\n", potentialNumber);
  // //
  // // variable to hold force and energy
  // //
  // std::pair<double,double> forceEnergy;
  // forceEnergy.first = 0.0;
  // forceEnergy.second = 0.0;

  // return forceEnergy;

  //
  // switch between cases
  //
  switch (d_PotentialType[potentialNumber]) {

  //
  // LJ case
  //
  case 0:

    //
    // return force and energy values
    //
    return LJ(potentialNumber, separation);

    /* NOTREACHED */
    break;

  //
  // EAM Johnson Case
  //
  case 1:

    //
    // return force and energy values
    //
    return EAMJohnson(potentialNumber, separation);

    /* NOTREACHED */
    break;

  //
  // EAM Data Case
  //
  case 2:

    //
    // return force and energy values
    //
    return EAMData(potentialNumber, separation);

    /* NOTREACHED */
    break;

  //
  // LJ NiMn case
  //
  case 3:

    //
    // return force and energy values
    //
    return LJNiMn(potentialNumber, separation);

    /* NOTREACHED */
    break;

  //
  // Buckingham
  //
  case 4:

    //
    // return force and energy values
    //
    return Buckingham(potentialNumber, separation);

    /* NOTREACHED */
    break;

  //
  // Harmonic
  //
  case 5:

    //
    // return force and energy values
    //
    return Harmonic(potentialNumber, separation);

    /* NOTREACHED */
    break;

  //
  // Anharmonic
  //
  case 6:

    //
    // return force and energy values
    //
    return Anharmonic(potentialNumber, separation);

    /* NOTREACHED */
    break;

  default:

    //
    // error
    //
    if (potentialNumber < 0)
      d_print("for negative potential number, getPairForce() is called\n");
    else if (potentialNumber >= d_PotentialType.size()) {
      d_print("potential number = %d is above d_PotentialType.size() = %d\n",
              potentialNumber, d_PotentialType.size());
    } else
      d_print("potentialNumber = %d, potentialType = %d\n", potentialNumber,
              d_PotentialType[potentialNumber]);

    d_print("potential number = %d, potentialType\n", potentialNumber);
    D_ERROR("Potential Type not implemented");
    exit(EXIT_FAILURE);
  }

  //
  // not reached
  //
  std::pair<double, double> returndummy;
  return returndummy;
}

//
// determine if pairs interact
//
int PairPotentials::doQuasicontinuumInteract(int quasicontinuumIdOne,
                                             int quasicontinuumIdTwo) {

  //
  // variable to determine if pairs interact
  //
  int pairInteract = -1;

  //
  // loop over all interaction pairs and see if they match input ids
  //
  for (unsigned int i = 0; i < d_PotentialInteractions.size(); ++i) {
    //
    // check if the pairs match
    //
    if (d_PotentialInteractions[i].first == quasicontinuumIdOne &&
        d_PotentialInteractions[i].second == quasicontinuumIdTwo) {

      //
      // set pair interact variable
      //
      pairInteract = i;

      //
      // return interaction
      //
      return pairInteract;
    }
  }

  //
  // check if the pairs match opposite numbers
  //
  for (unsigned int i = 0; i < d_PotentialInteractions.size(); ++i) {
    if (d_PotentialInteractions[i].first == quasicontinuumIdTwo &&
        d_PotentialInteractions[i].second == quasicontinuumIdOne) {
      //
      // set pair interact variable
      //
      pairInteract = i;

      //
      // return interaction potential number
      //
      return pairInteract;
    }
  }

  //
  //
  //
  return pairInteract;
}

//
// get cutoff radius
//
double PairPotentials::getCutoffRadius(int potentialNumber) {

  //
  //
  //
  return d_PotentialParameters[potentialNumber][0];
}

//
//  getCutoffRadiusNeighborList()
//
double PairPotentials::getCutoffRadiusNeighborList(const int potentialNumber) {
  // return cutoff radius to build cluster list
  return CUTOFF_NEIGHBOR_SCALE * d_PotentialParameters[potentialNumber][0];
}

//
//  getCutoffRadiusClusterList()
//
double PairPotentials::getCutoffRadiusClusterList(const int potentialNumber) {
  // return cutoff radius to build cluster list
  return CUTOFF_CLUSTER_SCALE * d_PotentialParameters[potentialNumber][0];
}

//
// get density
//
std::pair<double, double>
PairPotentials::getDensityAndDerivative(int potentialNumber, double r) {

  //
  // variable to hold force and energy
  //
  std::pair<double, double> density;

  //
  // switch between cases
  //
  switch (d_PotentialType[potentialNumber]) {
  //
  // not used
  //
  case 2: {
    //
    // get variables
    //
    const double rCut = d_PotentialParameters[potentialNumber][0];

    //
    // get number of grid points
    //
    int nGrid = d_RhoX[potentialNumber].size();

    //
    // check if less than rCut
    //
    if (r <= rCut) {
      //
      // find force and energy
      //
      SPLINE_F77(&nGrid, &r, &(d_RhoX[potentialNumber][0]),
                 &(d_RhoY[potentialNumber][0]), &(d_RhoB[potentialNumber][0]),
                 &(d_RhoC[potentialNumber][0]), &(d_RhoD[potentialNumber][0]),
                 &density.first, &density.second);
    }

    return density;
    /* NOTREACHED */
    break;
  }

  default:
    //
    // error
    //
    ERROR("Potential Type not implemented");
    exit(EXIT_FAILURE);
  }

  //
  //
  //
  return density;
}

//
// get density
//
std::pair<double, double>
PairPotentials::getEmbeddingFunctionAndDerivative(int potentialNumber,
                                                  double density) {
  //
  // variable to hold force and energy
  //
  std::pair<double, double> embeddingFunction;

  //
  // switch between cases
  //
  switch (d_PotentialType[potentialNumber]) {
  //
  // not used
  //
  case 2: {
    //
    // get number of grid points
    //
    int nGrid = d_EmbedX[potentialNumber].size();

    //
    // find force and energy
    //
    SPLINE_F77(&nGrid, &density, &(d_EmbedX[potentialNumber][0]),
               &(d_EmbedY[potentialNumber][0]), &(d_EmbedB[potentialNumber][0]),
               &(d_EmbedC[potentialNumber][0]), &(d_EmbedD[potentialNumber][0]),
               &embeddingFunction.first, &embeddingFunction.second);

    //
    //
    //
    return embeddingFunction;

    /* NOTREACHED */
    break;
  }

  default:
    //
    // error
    //
    ERROR("Potential Type not implemented");
    exit(EXIT_FAILURE);
  }
  //
  //
  //
  return embeddingFunction;
}

//
// LJ calculations
//
std::pair<double, double> PairPotentials::LJ(int potentialNumber, double r) {

  //
  // Parameter list:
  // 0 = cutoff radius (MUST BE CUTOFF Radius)
  // 1 = epsilon
  // 2 = epsilon x 12.0
  // 3 = sigma^12
  // 4 = sigma^6 x 2.0
  // 5 = sigma^6

  //
  // variable to hold force and energy
  //
  std::pair<double, double> forceEnergy;

  //
  // calculate force and energy
  //
  if (r > d_PotentialParameters[potentialNumber][0]) {

    //
    // no force or energy
    //
    forceEnergy.first = 0.0;
    forceEnergy.second = 0.0;

  }

  else {

    //
    // variables
    //
    const double r6 = r * r * r * r * r * r;
    const double r12 = r6 * r6;
    const double eps_4 = d_PotentialParameters[potentialNumber][1];
    const double eps_48 = d_PotentialParameters[potentialNumber][2];
    const double sigma_12 = d_PotentialParameters[potentialNumber][3];
    const double sigma_6_2 = d_PotentialParameters[potentialNumber][4];
    const double sigma_6 = d_PotentialParameters[potentialNumber][5];

    //
    // calculate force and energy
    //
    forceEnergy.first = eps_48 * (-sigma_12 / (r12 * r) + sigma_6 / (r6 * r));
    forceEnergy.second = eps_4 * (sigma_12 / r12 - sigma_6_2 / r6);
  }

  //
  //
  //
  return forceEnergy;
}

//
// EAMJohnson pair potential
//
std::pair<double, double> PairPotentials::EAMJohnson(int potentialNumber,
                                                     double r) {

  //
  // Parameter list:
  // 0 = r_cutoff (MUST BE CUTOFF Radius)
  // 1 = r_c
  // 2 = c
  // 3 = d
  // 4 = rgama
  // 5 = gama
  // 6 = r_e
  // 7 = gamma_phi
  // 8 = phie
  // 9 = a
  // 10 = b

  //
  // get variables
  //
  const double r_cutoff = d_PotentialParameters[potentialNumber][0];
  const double r_c = d_PotentialParameters[potentialNumber][1];
  const double c = d_PotentialParameters[potentialNumber][2];
  const double d = d_PotentialParameters[potentialNumber][3];
  const double rgama = d_PotentialParameters[potentialNumber][4];
  const double gama = d_PotentialParameters[potentialNumber][5];
  const double r_e = d_PotentialParameters[potentialNumber][6];
  const double gamma_phi = d_PotentialParameters[potentialNumber][7];
  const double phie = d_PotentialParameters[potentialNumber][8];
  const double a = d_PotentialParameters[potentialNumber][9];
  const double b = d_PotentialParameters[potentialNumber][10];

  //
  // variable to hold force and energy
  //
  std::pair<double, double> forceEnergy;

  //
  // check if less than cutoff radius
  //
  if (r <= r_cutoff) {

    //
    // check if less than r_c
    //
    if (r <= r_c) {

      //
      // local variable
      //
      double r_12 = r * r * r * r * r * r * r * r * r * r * r * r;

      //
      // set force and energy values
      //
      forceEnergy.first = -12.0 * c / (r_12 * r);
      forceEnergy.second = c / r_12 + d;

      //
      //
      //
      return forceEnergy;
    }

    //
    // check if less than rgama
    //
    if (r <= rgama) {

      //
      // local variable
      //
      double exponent = exp(-gama * (r / r_e - 1.0));

      //
      // set force and energy values
      //
      forceEnergy.first = gamma_phi * exponent;
      forceEnergy.second = phie * exponent;

      //
      //
      //
      return forceEnergy;

    }

    //
    // if not less than r_c or rgama
    //
    else {

      double dr = r_cutoff - r;

      //
      // set force and energy values
      //
      forceEnergy.first = (-3.0 * a * dr - 2.0 * b) * dr;
      forceEnergy.second = (a * dr + b) * dr * dr;

      //
      //
      //
      return forceEnergy;
    }

  } else {

    //
    // no force or energy
    //
    forceEnergy.first = 0.0;
    forceEnergy.second = 0.0;
  }

  //
  //
  //
  return forceEnergy;
}

//
// EAMData pair potential
//
std::pair<double, double> PairPotentials::EAMData(int potentialNumber,
                                                  double r) {

  //
  // Parameter list:
  // 0 = rCut (MUST BE CUTOFF Radius)
  //

  //
  // variable to hold force and energy
  //
  std::pair<double, double> forceEnergy;

  //
  // get variables
  //
  const double rCut = d_PotentialParameters[potentialNumber][0];

  //
  // get number of grid points
  //
  int nGrid = d_PairX[potentialNumber].size();

  //
  // check if less than rCut
  //
  if (r <= rCut) {

    //
    // find force and energy
    //
    SPLINE_F77(&nGrid, &r, &(d_PairX[potentialNumber][0]),
               &(d_PairY[potentialNumber][0]), &(d_PairB[potentialNumber][0]),
               &(d_PairC[potentialNumber][0]), &(d_PairD[potentialNumber][0]),
               &forceEnergy.second, &forceEnergy.first);

    //
    //
    //
    return forceEnergy;

  }

  //
  // if not less than rmax
  //
  else {

    //
    // no force or energy
    //
    forceEnergy.first = 0.0;
    forceEnergy.second = 0.0;
  }

  //
  //
  //
  return forceEnergy;
}

//
// LJNiMn calculations
//
std::pair<double, double> PairPotentials::LJNiMn(int potentialNumber,
                                                 double r) {

  //
  // Parameter list:
  // 0 = cutoff radius (MUST BE CUTOFF Radius)
  // 1 = sigma
  // 2 = sigma^6
  // 3 = sigma^12
  // 4 = epsilon

  //
  // variable to hold force and energy
  //
  std::pair<double, double> forceEnergy;

  //
  // calculate force and energy
  //
  if (r > d_PotentialParameters[potentialNumber][0]) {

    //
    // no force or energy
    //
    forceEnergy.first = 0.0;
    forceEnergy.second = 0.0;

  }

  else {

    //
    // variables
    //
    const double r2 = r * r;
    const double r6 = r2 * r2 * r2;
    const double r12 = r6 * r6;
    const double sigma = d_PotentialParameters[potentialNumber][1];
    const double sigma6 = d_PotentialParameters[potentialNumber][2];
    const double sigma12 = d_PotentialParameters[potentialNumber][3];
    const double epsilon = d_PotentialParameters[potentialNumber][4];

    //
    // calculate force and energy
    //
    forceEnergy.first =
        4.0 * epsilon * (-12.0 * sigma12 / (r12 * r) + 6.0 * sigma6 / (r6 * r));
    forceEnergy.second = 4.0 * epsilon * (sigma12 / r12 - sigma6 / r6);
  }

  //
  //
  //
  return forceEnergy;
}

//
// Buckingham calculations
//
std::pair<double, double> PairPotentials::Buckingham(int potentialNumber,
                                                     double r) {

  //
  // Parameter list:
  // 0 = cutoff radius (MUST BE CUTOFF Radius)
  // 1 = A
  // 2 = rho
  // 3 = C
  // 4 = shifted energy constant
  // 5 = shifted force constant

  //
  // variable to hold force and energy
  //
  std::pair<double, double> forceEnergy;

  const double cutoffRadius = d_PotentialParameters[potentialNumber][0];

  //
  // calculate force and energy
  //
  if (r > cutoffRadius) {

    //
    // no force or energy
    //
    forceEnergy.first = 0.0;
    forceEnergy.second = 0.0;

  }

  else {

    //
    // variables
    //
    const double r2 = r * r;
    const double r6 = r2 * r2 * r2;
    const double A = d_PotentialParameters[potentialNumber][1];
    const double rho = d_PotentialParameters[potentialNumber][2];
    const double C = d_PotentialParameters[potentialNumber][3];
    const double vC = d_PotentialParameters[potentialNumber][4];
    const double vF = d_PotentialParameters[potentialNumber][5];

    //
    // calculate force
    //
    forceEnergy.first = -A / rho * std::exp(-r / rho) + 6.0 * C / (r6 * r) - vF;

    //
    // calculate energy
    //
    forceEnergy.second =
        A * std::exp(-r / rho) - C / r6 - vC - vF * (r - cutoffRadius);
  }

  //
  //
  //
  return forceEnergy;
}

//
// Harmonic calculations
//
std::pair<double, double> PairPotentials::Harmonic(int potentialNumber,
                                                   double r) {

  //
  // Parameter list:
  // 0 = cutoff radius (MUST BE CUTOFF Radius)
  // 1 = k
  // 2 = C

  //
  // variable to hold force and energy
  //
  std::pair<double, double> forceEnergy;

  //
  // calculate force and energy
  //
  if (r > d_PotentialParameters[potentialNumber][0]) {

    //
    // no force or energy
    //
    forceEnergy.first = 0.0;
    forceEnergy.second = 0.0;

  }

  else {

    //
    // variables
    //
    const double k = d_PotentialParameters[potentialNumber][1];
    const double C = d_PotentialParameters[potentialNumber][2];

    //
    // calculate force
    //
    forceEnergy.first = k * (r - C);

    //
    // calculate energy
    //
    forceEnergy.second = 0.5 * k * (r - C) * (r - C);
  }

  //
  //
  //
  return forceEnergy;
}

//
// Anharmonic calculations
//
std::pair<double, double> PairPotentials::Anharmonic(int potentialNumber,
                                                     double r) {
  //
  // Parameter list:
  // 0 = r_cut
  // 1 = k2  // coefficient associated to r*r
  // 2 = k4  // coefficient associated to r*r*r*r
  //
  // e(r) = k2*r*r/2 + k4*r*r*r*r/24
  //

  //
  // variable to hold force and energy
  //
  std::pair<double, double> forceEnergy;

  double r_cut = d_PotentialParameters[potentialNumber][0];

  if (r < r_cut) {
    //
    // calculate force and energy
    //
    const double k2 = d_PotentialParameters[potentialNumber][1];
    const double k4 = d_PotentialParameters[potentialNumber][2];

    //
    // calculate force
    //
    forceEnergy.first = k2 * r + k4 * r * r * r / 6.0;

    //
    // calculate energy
    //
    forceEnergy.second = 0.5 * k2 * r * r + k4 * r * r * r * r / 24.0;
  } else {
    forceEnergy.first = 0.0;
    forceEnergy.second = 0.0;
  }
  //
  //
  //
  return forceEnergy;
}

//
//  setUniversalConstants()
//
void PairPotentials::setUniversalConstants(std::vector<double> constants) {
  //  c_universalConst[0] : avagrado number
  //  c_universalConst[0] : boltzman constant
  //  c_universalConst[0] : maxplanck constant
  //  c_universalConst[0] : maxplanck over 2pi constant
  //  c_universalConst[0] : electric constant
  //  c_universalConst[0] : electric constant * 4pi

  // boltzman constant
  c_universalConst[1] = constants[0];

  // max planck constant
  c_universalConst[2] = constants[1];

  // max planck / 2pi
  c_universalConst[3] = c_universalConst[2] / (2 * M_PI);

  //  electric constant
  c_universalConst[4] = constants[2];

  // electric constant * 4 Pi
  c_universalConst[5] = constants[2] * 4 * M_PI;

  return;
}

//
//  getBoltzmanConstant()
//
double PairPotentials::getBoltzmanConstant(void) { return c_universalConst[1]; }

//
//  getElectricConstant()
//
double PairPotentials::getElectricConstant() { return c_universalConst[4]; }

//
//  getUniversalConstants()
//
std::vector<double> PairPotentials::getUniversalConstants() {
  return c_universalConst;
}

//
//  getPotentialType()
//
int PairPotentials::getPotentialType(int potentialNumber) {
  if (potentialNumber >= 0 && potentialNumber < d_PotentialType.size())
    return d_PotentialType[potentialNumber];
  else
    return -1;
}

//
// get forces and energies
//
void PairPotentials::getHarmonicApproxData(double &phi, double &dphi,
                                           double &d_dphi, int potentialNumber,
                                           double separation) {
  //
  // switch between cases
  //
  switch (d_PotentialType[potentialNumber]) {

  //
  // LJ case
  //
  case 0:

    //
    // return force and energy values
    //
    LJ_Harmonic(phi, dphi, d_dphi, potentialNumber, separation);

    return;

    /* NOTREACHED */
    break;

  //
  // EAM Johnson Case
  //
  case 1:

    //
    // return force and energy values
    //
    d_print("harmonic approx for eam not implemented\n");
    exit(EXIT_SUCCESS);

    /* NOTREACHED */
    break;

  //
  // EAM Data Case
  //
  case 2:

    //
    // return force and energy values
    //
    // return EAMData(potentialNumber, separation);
    d_print("harmonic approx for eam not implemented\n");
    exit(EXIT_SUCCESS);

    /* NOTREACHED */
    break;

  //
  // LJ NiMn case
  //
  case 3:

    //
    // return force and energy values
    //
    LJNiMn_Harmonic(phi, dphi, d_dphi, potentialNumber, separation);

    return;

    /* NOTREACHED */
    break;

  //
  // Buckingham
  //
  case 4:

    //
    // return force and energy values
    //
    Buckingham_Harmonic(phi, dphi, d_dphi, potentialNumber, separation);

    return;

    /* NOTREACHED */
    break;

  //
  // Harmonic
  //
  case 5:

    //
    // return force and energy values
    //
    return Harmonic_Harmonic(phi, dphi, d_dphi, potentialNumber, separation);

    return;

    /* NOTREACHED */
    break;

  //
  // Anharmonic
  //
  case 6:

    //
    // return force and energy values
    //
    return Anharmonic_Harmonic(phi, dphi, d_dphi, potentialNumber, separation);

    return;

    /* NOTREACHED */
    break;

  default:

    //
    // error
    //
    if (potentialNumber < 0)
      d_print("for negative potential number, getPairForce() is called\n");
    else if (potentialNumber >= d_PotentialType.size()) {
      d_print("potential number = %d is above d_PotentialType.size() = %d\n",
              potentialNumber, d_PotentialType.size());
    } else
      d_print("potentialNumber = %d, potentialType = %d\n", potentialNumber,
              d_PotentialType[potentialNumber]);

    d_print("potential number = %d, potentialType\n", potentialNumber);
    D_ERROR("Potential Type not implemented");
    exit(EXIT_FAILURE);
  }

  //
  // not reached
  //
  return;
}

//
// LJ calculations
//
void PairPotentials::LJ_Harmonic(double &phi, double &dphi, double &d_dphi,
                                 int potentialNumber, double r) {

  //
  // Parameter list:
  // 0 = cutoff radius (MUST BE CUTOFF Radius)
  // 1 = epsilon
  // 2 = epsilon x 12.0
  // 3 = sigma^12
  // 4 = sigma^6 x 2.0
  // 5 = sigma^6

  //
  // calculate force and energy
  //
  if (r > d_PotentialParameters[potentialNumber][0]) {
    //
    // no force or energy
    //
    phi = 0.0;
    dphi = 0.0;
    d_dphi = 0.0;
  } else {
    //
    // variables
    //
    const double r6 = r * r * r * r * r * r;
    const double r12 = r6 * r6;
    const double eps_4 = d_PotentialParameters[potentialNumber][1];
    const double eps_48 = d_PotentialParameters[potentialNumber][2];
    const double sigma_12 = d_PotentialParameters[potentialNumber][3];
    const double sigma_6_2 = d_PotentialParameters[potentialNumber][4];
    const double sigma_6 = d_PotentialParameters[potentialNumber][5];

    //
    // calculate force and energy
    //
    phi = eps_4 * (sigma_12 / r12 - sigma_6_2 / r6);
    dphi = eps_48 * (-sigma_12 / (r12 * r) + sigma_6 / (r6 * r));
    d_dphi = eps_4 * 12.0 *
             (13.0 * sigma_12 / (r12 * r * r) - 7.0 * sigma_6 / (r6 * r * r));
  }

  //
  //
  //
  return;
}

//
// LJNiMn calculations
//
void PairPotentials::LJNiMn_Harmonic(double &phi, double &dphi, double &d_dphi,
                                     int potentialNumber, double r) {

  //
  // Parameter list:
  // 0 = cutoff radius (MUST BE CUTOFF Radius)
  // 1 = sigma
  // 2 = sigma^6
  // 3 = sigma^12
  // 4 = epsilon

  //
  // variable to hold force and energy
  //
  std::pair<double, double> forceEnergy;

  //
  // calculate force and energy
  //
  if (r > d_PotentialParameters[potentialNumber][0]) {
    //
    // no force or energy
    //
    phi = 0.0;
    dphi = 0.0;
    d_dphi = 0.0;
  } else {
    //
    // variables
    //
    const double r2 = r * r;
    const double r6 = r2 * r2 * r2;
    const double r12 = r6 * r6;
    const double sigma = d_PotentialParameters[potentialNumber][1];
    const double sigma6 = d_PotentialParameters[potentialNumber][2];
    const double sigma12 = d_PotentialParameters[potentialNumber][3];
    const double epsilon = d_PotentialParameters[potentialNumber][4];

    //
    // calculate force and energy
    //
    phi = 4.0 * epsilon * (sigma12 / r12 - sigma6 / r6);
    dphi =
        4.0 * epsilon * (-12.0 * sigma12 / (r12 * r) + 6.0 * sigma6 / (r6 * r));
    d_dphi = 4.0 * epsilon * (12.0 * 13.0 * sigma12 / (r12 * r * r) -
                              6.0 * 7.0 * sigma6 / (r6 * r * r));
  }

  //
  //
  //
  return;
}

//
// Buckingham calculations
//
void PairPotentials::Buckingham_Harmonic(double &phi, double &dphi,
                                         double &d_dphi, int potentialNumber,
                                         double r) {
  //
  // Parameter list:
  // 0 = cutoff radius (MUST BE CUTOFF Radius)
  // 1 = A
  // 2 = rho
  // 3 = C
  // 4 = shifted energy constant
  // 5 = shifted force constant

  //
  // variable to hold force and energy
  //
  std::pair<double, double> forceEnergy;

  const double cutoffRadius = d_PotentialParameters[potentialNumber][0];

  //
  // calculate force and energy
  //
  if (r > cutoffRadius) {
    //
    // no force or energy
    //
    phi = 0.0;
    dphi = 0.0;
    d_dphi = 0.0;
  } else {
    //
    // variables
    //
    const double r2 = r * r;
    const double r6 = r2 * r2 * r2;
    const double A = d_PotentialParameters[potentialNumber][1];
    const double rho = d_PotentialParameters[potentialNumber][2];
    const double C = d_PotentialParameters[potentialNumber][3];
    const double vC = d_PotentialParameters[potentialNumber][4];
    const double vF = d_PotentialParameters[potentialNumber][5];

    //
    //
    //
    phi = A * std::exp(-r / rho) - C / r6 - vC - vF * (r - cutoffRadius);
    dphi = -(A * std::exp(-r / rho)) / rho + 6.0 * C / (r6 * r) - vF;
    d_dphi =
        (A * std::exp(-r / rho)) / (rho * rho) - 6.0 * 7.0 * C / (r6 * r * r);
  }

  //
  //
  //
  return;
}

//
// Harmonic calculations
//
void PairPotentials::Harmonic_Harmonic(double &phi, double &dphi,
                                       double &d_dphi, int potentialNumber,
                                       double r) {
  //
  // Parameter list:
  // 0 = cutoff radius (MUST BE CUTOFF Radius)
  // 1 = k
  // 2 = C

  //
  // variable to hold force and energy
  //
  std::pair<double, double> forceEnergy;

  //
  // calculate force and energy
  //
  if (r > d_PotentialParameters[potentialNumber][0]) {
    //
    // no force or energy
    //
    phi = 0.0;
    dphi = 0.0;
    d_dphi = 0.0;
  } else {
    //
    // variables
    //
    const double k = d_PotentialParameters[potentialNumber][1];
    const double C = d_PotentialParameters[potentialNumber][2];

    phi = 0.5 * k * (r - C) * (r - C);
    dphi = k * (r - C);
    d_dphi = k;
  }

  //
  //
  //
  return;
}

//
// Anharmonic calculations
//
void PairPotentials::Anharmonic_Harmonic(double &phi, double &dphi,
                                         double &d_dphi, int potentialNumber,
                                         double r) {
  //
  // Parameter list:
  // 0 = r_cut
  // 1 = k2  // coefficient associated to r*r
  // 2 = k4  // coefficient associated to r*r*r*r
  //
  // e(r) = k2*r*r/2 + k4*r*r*r*r/24
  //

  //
  // variable to hold force and energy
  //
  std::pair<double, double> forceEnergy;

  double r_cut = d_PotentialParameters[potentialNumber][0];

  if (r < r_cut) {
    //
    // calculate force and energy
    //
    const double k2 = d_PotentialParameters[potentialNumber][1];
    const double k4 = d_PotentialParameters[potentialNumber][2];

    phi = 0.5 * k2 * r * r + k4 * r * r * r * r / 24.0;
    dphi = k2 * r + k4 * r * r * r / 6.0;
    d_dphi = k2 + 0.5 * k4 * r * r;
  } else {
    //
    // no force or energy
    //
    phi = 0.0;
    dphi = 0.0;
    d_dphi = 0.0;
  }
  //
  //
  //
  return;
}
}
