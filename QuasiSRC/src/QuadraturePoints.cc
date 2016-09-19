//
// QuadraturePoints.cc
//

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include "C_Interface.h"
#include "QuadraturePoints.h"

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
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <utility>

#include "PairPotentials.h"

//
//
//

namespace quasicontinuum {
//
//
//

QuadraturePoints *QuadraturePoints::_instance = NULL;

//
// constructor
//

QuadraturePoints::QuadraturePoints() {

  //
  //
  //
  return;
}

//
// destructor
//

QuadraturePoints::~QuadraturePoints() {

  //
  //
  //
  return;
}

//
// getInstance method
//

QuadraturePoints *QuadraturePoints::getInstance() {

  //
  // if not created, create
  //
  if (_instance == NULL) {
    _instance = new QuadraturePoints();
  }

  //
  // return instance
  //
  return _instance;
}

//
// destroy instance method
//

void QuadraturePoints::destroyInstance() {

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
// get the quadrature point
//
void QuadraturePoints::getQuadraturePoints(
    std::vector<std::vector<std::vector<double>>> &quadVectors,
    std::vector<double> &quadWeights, int n) {

  //
  // Loop over method
  //
  int quadMethod = PairPotentials::getInstance()->d_quadMethod;

  switch (quadMethod) {

  //
  // 3 point quadrature method
  //
  case 3: {

    //
    // parameters for 3 point method
    //

    int dimension;
    long double V, r;

    dimension = 3 * n;

    V = std::pow(3.14159265358979323846, (double)dimension / 2);
    r = std::sqrt((double)n / 2);

    //
    // resize quadVectors and quadWeights
    //

    quadVectors.resize(3 * n * 2);
    quadWeights.resize(3 * n * 2);

    for (int i = 0; i < quadVectors.size(); i++) {
      quadVectors[i].resize(n);

      for (int j = 0; j < n; j++) {
        quadVectors[i][j].resize(3);
      }
    }

    //
    // put vectors and weight into quadVectors and quadWeights
    //

    for (int i = 0; i < (3 * n); i++) {
      for (int j = 0; j < (n - 1); j++) {
        for (int k = 0; k < 3; k++) {
          if (3 * j + k == i) {
            quadVectors[i][j][k] = r;
            quadVectors[i + 3 * n][j][k] = -r;
          }

          else {
            quadVectors[i][j][k] = 0;
            quadVectors[i + 3 * n][j][k] = 0;
          }
        }
      }

      quadWeights[i] = V / (3 * n);
      quadWeights[i + 3 * n] = -quadWeights[i];
    }

    return;
    break;
  }

  //
  // 5 point quadrature method
  //
  case 5: {

    //
    // NOT IMPLEMENTED YET
    //

    return;
    break;
  }

  default:

    //
    // error
    //
    std::cout << "Exit: Quadrature Point function" << std::endl;
    exit(EXIT_FAILURE);
  }

  return;
}

//
//  getQuadratureVectors
//
std::vector<std::pair<std::vector<std::vector<double>>, double>>
QuadraturePoints::getQuadratureVectors(int number, int dimension) {
  std::vector<std::pair<std::vector<std::vector<double>>, double>> returnData;

  int quadMethod = PairPotentials::getInstance()->d_quadMethod;

  //
  // Loop over method
  //
  switch (quadMethod) {

  //
  // 3 point quadrature method
  //
  case 3: {
    //
    // parameters for 3 point method
    //
    int n;
    long double V, r, W;
    int number_of_quad_vecs;

    n = dimension * number;
    number_of_quad_vecs = 2 * n; // for method 3

    V = std::pow(M_PI, (double)n / 2);
    W = V / (number_of_quad_vecs);
    r = std::sqrt((double)n / 2);

    //
    // resize quadVectors and quadWeights
    //
    returnData.resize(number_of_quad_vecs);

    for (int i = 0; i < returnData.size(); i++) {
      returnData[i].first.resize(number);

      for (int j = 0; j < number; j++) {
        returnData[i].first[j].resize(dimension);
      }
    }

    //
    // put vectors and weight into quadVectors and quadWeights
    //
    //  since we will process two vectors (0,...0,+r,0,..,0) and
    //  (0,...0,-r,0,..,0)
    //  we loop over only half of number of quad vectors.
    //
    for (int quad = 0; quad < number_of_quad_vecs / 2; quad++) {
      for (int vec = 0; vec < number; vec++) {
        for (int dof = 0; dof < dimension; dof++) {
          int loc = dimension * vec + dof;

          // put non-zero value in only one place in vector
          if (loc == quad) {
            returnData[quad].first[vec][dof] = r;
            returnData[quad + n].first[vec][dof] = -r;
          } else {
            returnData[quad].first[vec][dof] = 0.0;
            returnData[quad + n].first[vec][dof] = 0.0;
          }
        }
      }

      returnData[quad].second = W;
      returnData[quad + n].second = W;
    }
  } // end of method = 3
  break;

  //
  //  case 5
  //
  case 5: {
    // not implemented ye
    d_print("Quad method = 5 is not implemented yet.\n");
    exit(EXIT_FAILURE);
  } break;

  //
  //  default
  //
  default: {
    d_print("Not valid quad method. Check quasi input and put valid quad "
            "method.\n");
    exit(EXIT_FAILURE);
  }
  } // end of switch

  // if reached here then we have quad vectors and weight
  return returnData;
} // end of getQuadratureVectors()
}
