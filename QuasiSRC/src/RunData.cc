//
// RunData.cc
//

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <vector>

#include "C_Interface.h"
#include "Quasicontinua.h"
#include "RunData.h"
#include "threads.h"

namespace quasicontinuum {

RunData *RunData::_instance = NULL;

RunData::RunData() {
  d_timeOut = false;

  return;
}

RunData::~RunData() {
  return;
}

RunData *RunData::getInstance() {

  // if not created, create
  if (_instance == NULL) {
    _instance = new RunData();
  }

  return _instance;
}

void RunData::destroyInstance() {
  delete _instance;
  return;
}

//  initializeData()
void RunData::initializeData() {
  // get Quasicontinua class instance
  Quasicontinua *quasicontinua = Quasicontinua::getInstance();

  int numQuasi = quasicontinua->size();

  // clear and resize all the data
  // d_timeOut = false;
  d_numQuasi = numQuasi;
  d_numberThreads = get_max_number_threads();

  d_computeClusterAndNeighborListTime.clear();
  d_buildNeighborListTime.clear();
  d_updateClusterAndNeighborListTime.clear();
  d_computePairwiseTime.clear();
  d_computeEAMTime.clear();
  d_computeElectricFieldTime.clear();
  d_computeElectrostaticsTime.clear();
  d_computeEntropyTime.clear();
  d_clusterForceToNodeForceTime.clear();
  d_initialToFirstForceCalculationTime.clear();
  d_computeResidualForceTime.clear();
  d_numberOfClusterSites.clear();

  for (int i = 0; i < numQuasi; i++) {
    d_computeClusterAndNeighborListTime.push_back(0.0);
    d_buildNeighborListTime.push_back(0.0);
    d_updateClusterAndNeighborListTime.push_back(0.0);
    d_computePairwiseTime.push_back(0.0);
    d_computeEAMTime.push_back(0.0);
    d_computeElectricFieldTime.push_back(0.0);
    d_computeElectrostaticsTime.push_back(0.0);
    d_computeEntropyTime.push_back(0.0);
    d_clusterForceToNodeForceTime.push_back(0.0);
    d_initialToFirstForceCalculationTime.push_back(0.0);
    d_computeResidualForceTime.push_back(0.0);
    d_numberOfClusterSites.push_back(0);
  }

  return;
}

//
//  addComputeClusterAndNeighborListTime()
//
void RunData::addComputeClusterAndNeighborListTime(int iQuasi, double time) {
  if (iQuasi >= d_numQuasi) {
    d_print("check the value of iQuasi supplied to RunData\n");
    return;
  }

  d_computeClusterAndNeighborListTime[iQuasi] = time;

  return;
}

//
//  addBuildNeighborListTime()
//
void RunData::addBuildNeighborListTime(int iQuasi, double time) {
  if (iQuasi >= d_numQuasi) {
    d_print("check the value of iQuasi supplied to RunData\n");
    return;
  }

  d_buildNeighborListTime[iQuasi] = time;

  return;
}

//
//  addUpdateClusterAndNeighborListTime()
//
void RunData::addUpdateClusterAndNeighborListTime(int iQuasi, double time) {
  if (iQuasi >= d_numQuasi) {
    d_print("check the value of iQuasi supplied to RunData\n");
    return;
  }

  d_updateClusterAndNeighborListTime[iQuasi] = time;

  return;
}

//
//  addComputePairwiseTime()
//
void RunData::addComputePairwiseTime(int iQuasi, double time) {
  if (iQuasi >= d_numQuasi) {
    d_print("check the value of iQuasi supplied to RunData\n");
    return;
  }

  d_computePairwiseTime[iQuasi] = time;

  return;
}

//
//  addComputeEAMTime()
//
void RunData::addComputeEAMTime(int iQuasi, double time) {
  if (iQuasi >= d_numQuasi) {
    d_print("check the value of iQuasi supplied to RunData\n");
    return;
  }

  d_computeEAMTime[iQuasi] = time;

  return;
}

//
//  addComputeElectricFieldTime()
//
void RunData::addComputeElectricFieldTime(int iQuasi, double time) {
  if (iQuasi >= d_numQuasi) {
    d_print("check the value of iQuasi supplied to RunData\n");
    return;
  }

  d_computeElectricFieldTime[iQuasi] = time;

  return;
}

//
//  addComputeElectrostaticsTime()
//
void RunData::addComputeElectrostaticsTime(int iQuasi, double time) {
  if (iQuasi >= d_numQuasi) {
    d_print("check the value of iQuasi supplied to RunData\n");
    return;
  }

  d_computeElectrostaticsTime[iQuasi] = time;

  return;
}

//
//  addComputeEntropyTime()
//
void RunData::addComputeEntropyTime(int iQuasi, double time) {
  if (iQuasi >= d_numQuasi) {
    d_print("check the value of iQuasi supplied to RunData\n");
    return;
  }

  d_computeEntropyTime[iQuasi] = time;

  return;
}

//
//  addClusterForceToNodeForceTime()
//
void RunData::addClusterForceToNodeForceTime(int iQuasi, double time) {
  if (iQuasi >= d_numQuasi) {
    d_print("check the value of iQuasi supplied to RunData\n");
    return;
  }

  d_clusterForceToNodeForceTime[iQuasi] = time;

  return;
}

//
//  addInitialToFirstForceCalculationTime()
//
void RunData::addInitialToFirstForceCalculationTime(int iQuasi, double time) {
  if (iQuasi >= d_numQuasi) {
    d_print("check the value of iQuasi supplied to RunData\n");
    return;
  }

  d_initialToFirstForceCalculationTime[iQuasi] = time;

  return;
}

//
//  addComputeResidualForceTime()
//
void RunData::addComputeResidualForceTime(int iQuasi, double time) {
  if (iQuasi >= d_numQuasi) {
    d_print("check the value of iQuasi supplied to RunData\n");
    return;
  }

  d_computeResidualForceTime[iQuasi] = time;

  return;
}

//
//  addNumberOfClusterSites()
//
void RunData::addNumberOfClusterSites(int iQuasi, int n) {
  if (iQuasi >= d_numQuasi) {
    d_print("check the value of iQuasi supplied to RunData\n");
    return;
  }

  d_numberOfClusterSites[iQuasi] = n;

  return;
}

//
//  writeData()
//
std::vector<int> RunData::writeData(std::vector<std::vector<double>> &data) {
  data.clear();

  data.push_back(d_computeClusterAndNeighborListTime);
  data.push_back(d_buildNeighborListTime);
  data.push_back(d_updateClusterAndNeighborListTime);
  data.push_back(d_computePairwiseTime);
  data.push_back(d_computeEAMTime);
  data.push_back(d_computeElectricFieldTime);
  data.push_back(d_computeElectrostaticsTime);
  data.push_back(d_computeEntropyTime);
  data.push_back(d_clusterForceToNodeForceTime);
  data.push_back(d_initialToFirstForceCalculationTime);
  data.push_back(d_computeResidualForceTime);

  return d_numberOfClusterSites;
}
}
