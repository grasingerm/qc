//
// Quasicontinua.cc
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

#include <cassert>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <utility>

#include "Error.h"
#include "PairPotentials.h"
#include "Quasicontinua.h"
#include "Quasicontinuum.h"

//
//
//

namespace quasicontinuum {

//
//
//

Quasicontinua *Quasicontinua::_instance = NULL;

//
// constructor
//

Quasicontinua::Quasicontinua() {

  //
  //
  //
  return;
}

//
// destructor
//

Quasicontinua::~Quasicontinua() {

  //
  //
  //
  return;
}

//
// getInstance method
//

Quasicontinua *Quasicontinua::getInstance() {

  //
  // if not created, create
  //
  if (_instance == NULL) {
    _instance = new Quasicontinua();
  }

  //
  // return instance
  //
  return _instance;
}

//
// destroy instance method
//

void Quasicontinua::destroyInstance() {

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
// insert instance of Quasicontinuum
//
void Quasicontinua::insert(int num_quasi, int ac, char *av[]) {

  //
  // insert quasi into quasicontinua
  //
  for (int i = 0; i < num_quasi; i++) {

    //
    // output current Quasicontinuum
    //
    d_print("Quasicontinua : Creating Quasicontinuum = %i\n", i);

    //
    // set current Quasicontinuum id
    //
    d_quasicontinuum_id = d_quasicontinua.size();

    //
    // instantiate a Quasicontinuum object
    //
    const Quasicontinuum quasicontinuum_instance(i, ac, av);

    //
    // insert Quasicontinuum object
    //
    Quasicontinua::d_quasicontinua.push_back(quasicontinuum_instance);
  }

  //
  // loop over each quasicontinuum and initialize energy and surface forces
  //
  // for (int iQuasi = 0; iQuasi < num_quasi; ++iQuasi){

  //   //
  //   // set current quasicontinuum id
  //   //
  //   d_quasicontinuum_id = iQuasi;

  //   //
  //   // initialize energy
  //   //
  //   d_quasicontinua[iQuasi].secondaryInitialization();

  // }

  //
  //
  //
  return;
}

//
// get Quasicontinua instance
//

const Quasicontinuum &Quasicontinua::get(id_type id) const {

  //
  // set current id
  //
  d_quasicontinuum_id = id; // FIXME:remove this dependence for calculations

  //
  // check validity of id
  //
  assert(static_cast<int>(id) >= 0 && id < d_quasicontinua.size());

  //
  // return Quasicontinuum instance
  //
  return (d_quasicontinua[id]);
}

//
// get Quasicontinuum instance
//

Quasicontinuum &Quasicontinua::get(id_type id) {

  //
  // set current id
  //
  d_quasicontinuum_id = id; // FIXME:remove this dependence for calculations

  //
  // check validity of id
  //
  assert(static_cast<int>(id) >= 0 && id < d_quasicontinua.size());

  //
  //
  //
  return (d_quasicontinua[id]);
}

//
// get Quasicontinua instance
//
// current id safe version of get()
//

const Quasicontinuum &Quasicontinua::getQuasi(id_type id) const {
  //
  // check validity of id
  //
  assert(static_cast<int>(id) >= 0 && id < d_quasicontinua.size());

  //
  // return Quasicontinuum instance
  //
  return (d_quasicontinua[id]);
}

//
// get Quasicontinuum instance
//

Quasicontinuum &Quasicontinua::getQuasi(id_type id) {
  //
  // check validity of id
  //
  assert(static_cast<int>(id) >= 0 && id < d_quasicontinua.size());

  //
  //
  //
  return (d_quasicontinua[id]);
}

//
// get size of quasicontinua
//

Quasicontinua::size_type Quasicontinua::size() const {

  //
  // return size of quasicontinua
  //
  return d_quasicontinua.size();
}

//
// get Id of current Quasicontinuum
//

Quasicontinua::id_type Quasicontinua::getCurrentId() const {

  //
  // return size of quasicontinua
  //
  return d_quasicontinuum_id;
}

//
// insert shift vector
//
void Quasicontinua::insertShift(int quasicontinuumId,
                                std::vector<double> shiftVector,
                                int coreShell) {

  //
  // pair to hold data
  //
  std::pair< int, std::vector<double> > data(quasicontinuumId, shiftVector);

  //
  // insert data into d_shifts
  //
  d_shifts.push_back(data);

  //
  // insert core shell type
  //
  d_coreShell.push_back(coreShell);

  //
  //
  //
  return;
}

//
// get relative shift between Quasicontinuum
//
std::vector<double> Quasicontinua::getRelativeShift(unsigned int iQuasi,
                                                    unsigned int jQuasi) {

  //
  // vector to hold relative shift
  //
  std::vector<double> shiftVector(3, 0.0);

  //
  // check Quasicontinuum ID's against size of d_shifts
  //
  assert(iQuasi < d_shifts.size() && static_cast<int>(iQuasi) >= 0);
  assert(jQuasi < d_shifts.size() && static_cast<int>(jQuasi) >= 0);

  //
  // vectors to hold shifts
  //
  std::vector<double> shiftOne;
  std::vector<double> shiftTwo;

  //
  // loop over all shifts and find
  //
  for (unsigned int iShift = 0; iShift < d_shifts.size(); ++iShift) {

    if (static_cast<int>(iQuasi) == d_shifts[iShift].first)
      shiftOne = d_shifts[iShift].second;

    if (static_cast<int>(jQuasi) == d_shifts[iShift].first)
      shiftTwo = d_shifts[iShift].second;
  
  }

  //
  // get shifts between Quasicontinuum
  //
  shiftVector[0] = shiftTwo[0] - shiftOne[0];
  shiftVector[1] = shiftTwo[1] - shiftOne[1];
  shiftVector[2] = shiftTwo[2] - shiftOne[2];

  //
  //
  //
  return shiftVector;
}

//
// get global shift vector of Quasicontinuum
//
std::vector<double> Quasicontinua::getShift(unsigned int QuasicontinuumID) {

  //
  // check QuasicontinuumID against size of d_shifts
  //
  assert(QuasicontinuumID < d_shifts.size() && 
         static_cast<int>(QuasicontinuumID) >= 0);

  //
  // loop over all shifts and find
  //
  for (unsigned int iShift = 0; iShift < d_shifts.size(); ++iShift) {
    //
    // check if QuasicontinuumID matches
    //
    if (int(QuasicontinuumID) == d_shifts[iShift].first) {
      //
      // return shift
      //
      return d_shifts[iShift].second;
    }
  }
}

//
// get whether core or shell type
//
int Quasicontinua::getCoreShell(unsigned int QuasicontinuumID) {

  //
  //
  //
  return d_coreShell[QuasicontinuumID];
}

//
//  setTemperature()
//
void Quasicontinua::setTemperature(const double temperature) {
  d_temperature = temperature;

  return;
}

//
//  getTemperature()
//
double Quasicontinua::getTemperature() { return d_temperature; }

//
//  getSigmaVector()
//
std::vector<double> Quasicontinua::getSigmaVector(void) {
  //
  //  vector to hold sigma for each quasi
  //
  std::vector<double> sigmaVector;

  //
  //  get handle for PairPotentials
  //
  PairPotentials *pairC = PairPotentials::getInstance();

  double k_B = pairC->getBoltzmanConstant();

  for (unsigned iQuasi = 0; iQuasi < d_quasicontinua.size(); ++iQuasi) {
    // get atomic mass for atoms of iQuasi
    long double mass = d_quasicontinua[iQuasi].getAtomicMass();

    // compute sigma
    double sigma = sqrt(2 * mass * k_B * d_temperature);

    // insert it into sigmaVector
    sigmaVector.push_back(sigma);
  }

  return sigmaVector;
}

}
