//
// Lattice.cc
//

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#ifdef STDC_HEADERS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#else
#error No standard C library headers found
#endif /* STDC_HEADERS */

#ifdef HAVE_MATH_H
#include <math.h>
#else
#error math.h not found.
#endif /* HAVE_MATH_H */

#include <iostream>
#include <vector>

#include "C_Interface.h"
#include "DataTypes.h"
#include "Element.h"
#include "Error.h"
#include "Lattice.h"
#include "MiscFunctions.h"
#include "Quasicontinua.h"
#include "Void.h"

// global variables
#define MAX_NUMBER_SHELLS 80
#define RETURN_OK 0
#define RETURN_ERR -1

//
//
//

namespace quasicontinuum {
//
//
//

int Lattice::d_initialize_bcc = 0;
int Lattice::d_initialize_fcc = 0;

int Lattice::d_bcc_number_shells = 0;
int Lattice::d_fcc_number_shells = 0;

int Lattice::d_bcc_p = 1;
int Lattice::d_fcc_p = 1;

int Lattice::d_initialConfigFlag = 0;

Lattice *Lattice::_instance = NULL;

//
// constructor
//

Lattice::Lattice() {

  //
  //
  //
  return;
}

//
// destructor
//

Lattice::~Lattice() {

  //
  //
  //
  return;
}

//
// getInstance method
//

Lattice *Lattice::getInstance() {

  //
  // if not created, create
  //
  if (_instance == NULL) {
    _instance = new Lattice();
  }

  //
  // return instance
  //
  return _instance;
}

//
// destroy instance method
//

void Lattice::destroyInstance() {

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
//  setInitialConfigFlag()
//
void Lattice::setInitialConfigFlag(int flag) {
  d_initialConfigFlag = flag;

  return;
}

// ***************************************//
// I = lattice.c

//
//  findClosestSite()
//
void Lattice::findClosestSite(int l[3], double X[3],
                              const struct lattice_t *P_lattice,
                              const int iQuasi) {
  double r[3];
  double d[3];
  double r_max = 12; // r_max is really squared distance
  double ld;

  int lp[3];
  int lr[3];
  int i;
  int j;
  int k;

  // find real indices of point
  GetLatticeCoordinates(r, X, P_lattice);

  //
  lp[0] = rint(r[0]);
  lp[1] = rint(r[1]);
  lp[2] = rint(r[2]);

  // check if lp is lattice site
  if (isLatticeSite(lp, P_lattice, iQuasi) == RETURN_SUCCESS) {
    l[0] = lp[0];
    l[1] = lp[1];
    l[2] = lp[2];

    return;
  }

  // lp is not lattice site. Create a cube of lattices sites centered at
  // lp and check all of them
  for (i = -1; i <= 1; i++)
    for (j = -1; j <= 1; j++)
      for (k = -1; k <= 1; k++) {
        lr[0] = lp[0] + i;
        lr[1] = lp[1] + j;
        lr[2] = lp[2] + k;

        /**
          * check if l is inside/outside sample
          */

        if ((lr[0] < P_lattice->l_start[0]) || (lr[0] > P_lattice->l_end[0]) ||
            (lr[1] < P_lattice->l_start[1]) || (lr[1] > P_lattice->l_end[1]) ||
            (lr[2] < P_lattice->l_start[2]) || (lr[2] > P_lattice->l_end[2]))
          continue;

        /**
          * check if l is a lattice site
          */

        if (isLatticeSite(lr, P_lattice, iQuasi) == RETURN_FAILURE)
          continue;

        /**
          * compute distance
          */

        d[0] = r[0] - lr[0];
        d[1] = r[1] - lr[1];
        d[2] = r[2] - lr[2];

        ld = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];

        if (ld < r_max) {

          r_max = ld;
          l[0] = lr[0];
          l[1] = lr[1];
          l[2] = lr[2];
        }
      }

  return;
}

//
//  isSiteInsideLattice()
//
int Lattice::isSiteInsideLattice(const struct lattice_t *P_lattice,
                                 const int l[3], const int iQuasi) {
  if (Void::getInstance()->isVoidEnable() != 0) {
    std::vector<int> site;
    for (int dof = 0; dof < 3; dof++)
      site.push_back(l[dof]);

    if (Void::getInstance()->findSiteInVoidCache(site, iQuasi) ==
        RETURN_SUCCESS)
      return RETURN_FAILURE;
  }

  if (l[0] >= P_lattice->l_start[0] && l[1] >= P_lattice->l_start[1] &&
      l[2] >= P_lattice->l_start[2] && l[0] <= P_lattice->l_end[0] &&
      l[1] <= P_lattice->l_end[1] && l[2] <= P_lattice->l_end[2])
    return (RETURN_SUCCESS);

  return (RETURN_FAILURE);
}

//
//  isLatticeSite()
//
int Lattice::isLatticeSite(const int l[3], const struct lattice_t *P_lattice,
                           const int iQuasi) {

  if (Void::getInstance()->isVoidEnable() != 0) {
    std::vector<int> site;
    for (int dof = 0; dof < 3; dof++)
      site.push_back(l[dof]);

    if (Void::getInstance()->findSiteInVoidCache(site, iQuasi) ==
        RETURN_SUCCESS)
      return RETURN_FAILURE;
  }

  switch (P_lattice->type) {
  case FCC:
    if ((l[0] + l[1] + l[2]) % 2 == 0)
      return RETURN_SUCCESS;
    break;

  case BCC:
    return RETURN_SUCCESS;
    break;
  }

  return (RETURN_FAILURE);
}

//
//  compareSites()
//
int Lattice::compareSites(const int l1[3], const int l2[3]) {
  if (memcmp(l1, l2, sizeof(int[3])) == 0)
    return (RETURN_SUCCESS);

  return (RETURN_FAILURE);
}

//
//  getSiteInitialPosition()
//
void Lattice::getSiteInitialPosition(double r[3], const int l[3],
                                     struct lattice_t *P_lattice,
                                     const int iQuasi) {
  if (d_initialConfigFlag == 0) {
    r[0] = l[0] * P_lattice->a1[0] + l[1] * P_lattice->a2[0] +
           l[2] * P_lattice->a3[0];

    r[1] = l[0] * P_lattice->a1[1] + l[1] * P_lattice->a2[1] +
           l[2] * P_lattice->a3[1];

    r[2] = l[0] * P_lattice->a1[2] + l[1] * P_lattice->a2[2] +
           l[2] * P_lattice->a3[2];

    return;
  } else {
    struct element_t *P_element;

    struct node_t *P_node_0;
    struct node_t *P_node_1;
    struct node_t *P_node_2;
    struct node_t *P_node_3;

    double shape[4];

    // get element
    P_element = Element::getInstance()->locateSiteElement(&(shape[0]), l,
                                                          P_lattice, iQuasi);

    P_node_0 = P_element->node[0];
    P_node_1 = P_element->node[1];
    P_node_2 = P_element->node[2];
    P_node_3 = P_element->node[3];

    // compute position
    r[0] = P_node_0->initial_position[0] * shape[0] +
           P_node_1->initial_position[0] * shape[1] +
           P_node_2->initial_position[0] * shape[2] +
           P_node_3->initial_position[0] * shape[3];

    r[1] = P_node_0->initial_position[1] * shape[0] +
           P_node_1->initial_position[1] * shape[1] +
           P_node_2->initial_position[1] * shape[2] +
           P_node_3->initial_position[1] * shape[3];

    r[2] = P_node_0->initial_position[2] * shape[0] +
           P_node_1->initial_position[2] * shape[1] +
           P_node_2->initial_position[2] * shape[2] +
           P_node_3->initial_position[2] * shape[3];

    return;
  }

  // if reached here then problem in the code
  d_print("Check getSiteInitialState() and initial flag of Lattice class \n");
  exit(EXIT_FAILURE);
}

//
//  getSiteInitialFrequency()
//
void Lattice::getSiteInitialFrequency(double &f, const int l[3],
                                      struct lattice_t *P_lattice,
                                      const int iQuasi) {
  struct element_t *P_element;

  struct node_t *P_node_0;
  struct node_t *P_node_1;
  struct node_t *P_node_2;
  struct node_t *P_node_3;

  double shape[4];

  // get element
  P_element = Element::getInstance()->locateSiteElement(&(shape[0]), l,
                                                        P_lattice, iQuasi);

  P_node_0 = P_element->node[0];
  P_node_1 = P_element->node[1];
  P_node_2 = P_element->node[2];
  P_node_3 = P_element->node[3];

  f = P_node_0->initial_frequency * shape[0] +
      P_node_1->initial_frequency * shape[1] +
      P_node_2->initial_frequency * shape[2] +
      P_node_3->initial_frequency * shape[3];

  return;
}

//
//  getSiteInitialTemperature()
//
void Lattice::getSiteInitialTemperature(double &T, const int l[3],
                                        struct lattice_t *P_lattice,
                                        const int iQuasi) {
  // for constant temperature quasi, call ghetTemperature() of
  // Quasicontinua to get temperature
  T = Quasicontinua::getInstance()->getTemperature();

  return;
}

//
//  getSiteInitialState()
//
void Lattice::getSiteInitialState(double S[4], const int l[3],
                                  struct lattice_t *P_lattice,
                                  const int iQuasi) {
  struct element_t *P_element;

  struct node_t *P_node_0;
  struct node_t *P_node_1;
  struct node_t *P_node_2;
  struct node_t *P_node_3;

  double shape[4];

  // get element
  P_element = Element::getInstance()->locateSiteElement(&(shape[0]), l,
                                                        P_lattice, iQuasi);

  P_node_0 = P_element->node[0];
  P_node_1 = P_element->node[1];
  P_node_2 = P_element->node[2];
  P_node_3 = P_element->node[3];

  // compute position
  if (d_initialConfigFlag == 0) {
    S[0] = l[0] * P_lattice->a1[0] + l[1] * P_lattice->a2[0] +
           l[2] * P_lattice->a3[0];

    S[1] = l[0] * P_lattice->a1[1] + l[1] * P_lattice->a2[1] +
           l[2] * P_lattice->a3[1];

    S[2] = l[0] * P_lattice->a1[2] + l[1] * P_lattice->a2[2] +
           l[2] * P_lattice->a3[2];
  } else {
    S[0] = P_node_0->initial_position[0] * shape[0] +
           P_node_1->initial_position[0] * shape[1] +
           P_node_2->initial_position[0] * shape[2] +
           P_node_3->initial_position[0] * shape[3];

    S[1] = P_node_0->initial_position[1] * shape[0] +
           P_node_1->initial_position[1] * shape[1] +
           P_node_2->initial_position[1] * shape[2] +
           P_node_3->initial_position[1] * shape[3];

    S[2] = P_node_0->initial_position[2] * shape[0] +
           P_node_1->initial_position[2] * shape[1] +
           P_node_2->initial_position[2] * shape[2] +
           P_node_3->initial_position[2] * shape[3];
  }

  S[3] = P_node_0->initial_frequency * shape[0] +
         P_node_1->initial_frequency * shape[1] +
         P_node_2->initial_frequency * shape[2] +
         P_node_3->initial_frequency * shape[3];

  return;
}

//
//  getSiteInitialStateInterpolate()
//
void Lattice::getSiteInitialStateInterpolate(double S[5], const int l[3],
                                             struct lattice_t *P_lattice,
                                             const int iQuasi) {
  struct element_t *P_element;

  struct node_t *P_node_0;
  struct node_t *P_node_1;
  struct node_t *P_node_2;
  struct node_t *P_node_3;

  double shape[4];

  // get element
  P_element = Element::getInstance()->locateSiteElement(&(shape[0]), l,
                                                        P_lattice, iQuasi);

  P_node_0 = P_element->node[0];
  P_node_1 = P_element->node[1];
  P_node_2 = P_element->node[2];
  P_node_3 = P_element->node[3];

  // compute position
  S[0] = P_node_0->initial_position[0] * shape[0] +
         P_node_1->initial_position[0] * shape[1] +
         P_node_2->initial_position[0] * shape[2] +
         P_node_3->initial_position[0] * shape[3];

  S[1] = P_node_0->initial_position[1] * shape[0] +
         P_node_1->initial_position[1] * shape[1] +
         P_node_2->initial_position[1] * shape[2] +
         P_node_3->initial_position[1] * shape[3];

  S[2] = P_node_0->initial_position[2] * shape[0] +
         P_node_1->initial_position[2] * shape[1] +
         P_node_2->initial_position[2] * shape[2] +
         P_node_3->initial_position[2] * shape[3];

  S[3] = P_node_0->initial_frequency * shape[0] +
         P_node_1->initial_frequency * shape[1] +
         P_node_2->initial_frequency * shape[2] +
         P_node_3->initial_frequency * shape[3];

  return;
}

//
//  getSiteCurrentPosition()
//
void Lattice::getSiteCurrentPosition(double r[3], const int l[3],
                                     struct lattice_t *P_lattice,
                                     const int iQuasi) {
  struct element_t *P_element;

  struct node_t *P_node_0;
  struct node_t *P_node_1;
  struct node_t *P_node_2;
  struct node_t *P_node_3;

  double shape[4];

  // get element
  P_element = Element::getInstance()->locateSiteElement(&(shape[0]), l,
                                                        P_lattice, iQuasi);

  P_node_0 = P_element->node[0];
  P_node_1 = P_element->node[1];
  P_node_2 = P_element->node[2];
  P_node_3 = P_element->node[3];

  // compute position
  r[0] = P_node_0->position[0] * shape[0] + P_node_1->position[0] * shape[1] +
         P_node_2->position[0] * shape[2] + P_node_3->position[0] * shape[3];

  r[1] = P_node_0->position[1] * shape[0] + P_node_1->position[1] * shape[1] +
         P_node_2->position[1] * shape[2] + P_node_3->position[1] * shape[3];

  r[2] = P_node_0->position[2] * shape[0] + P_node_1->position[2] * shape[1] +
         P_node_2->position[2] * shape[2] + P_node_3->position[2] * shape[3];

  return;
}

//
//  getSiteCurrentFrequency()
//
void Lattice::getSiteCurrentFrequency(double &f, const int l[3],
                                      struct lattice_t *P_lattice,
                                      const int iQuasi) {
  struct element_t *P_element;

  struct node_t *P_node_0;
  struct node_t *P_node_1;
  struct node_t *P_node_2;
  struct node_t *P_node_3;

  double shape[4];

  // get element
  P_element = Element::getInstance()->locateSiteElement(&(shape[0]), l,
                                                        P_lattice, iQuasi);

  P_node_0 = P_element->node[0];
  P_node_1 = P_element->node[1];
  P_node_2 = P_element->node[2];
  P_node_3 = P_element->node[3];

  // compute position
  f = P_node_0->frequency * shape[0] + P_node_1->frequency * shape[1] +
      P_node_2->frequency * shape[2] + P_node_3->frequency * shape[3];

  return;
}

//
//  getSiteCurrentTemperature()
//
void Lattice::getSiteCurrentTemperature(double &T, const int l[3],
                                        struct lattice_t *P_lattice,
                                        const int iQuasi) {
  // we are dealing with constant temperature multiscale method
  T = Quasicontinua::getInstance()->getTemperature();

  return;
}

//
//  getSiteCurrentPosition()
//
void Lattice::getSiteCurrentState(double s[4], const int l[3],
                                  struct lattice_t *P_lattice,
                                  const int iQuasi) {
  struct element_t *P_element;

  struct node_t *P_node_0;
  struct node_t *P_node_1;
  struct node_t *P_node_2;
  struct node_t *P_node_3;

  double shape[4];

  // get element
  P_element = Element::getInstance()->locateSiteElement(&(shape[0]), l,
                                                        P_lattice, iQuasi);

  P_node_0 = P_element->node[0];
  P_node_1 = P_element->node[1];
  P_node_2 = P_element->node[2];
  P_node_3 = P_element->node[3];

  // compute position
  s[0] = P_node_0->position[0] * shape[0] + P_node_1->position[0] * shape[1] +
         P_node_2->position[0] * shape[2] + P_node_3->position[0] * shape[3];

  s[1] = P_node_0->position[1] * shape[0] + P_node_1->position[1] * shape[1] +
         P_node_2->position[1] * shape[2] + P_node_3->position[1] * shape[3];

  s[2] = P_node_0->position[2] * shape[0] + P_node_1->position[2] * shape[1] +
         P_node_2->position[2] * shape[2] + P_node_3->position[2] * shape[3];

  s[3] = P_node_0->frequency * shape[0] + P_node_1->frequency * shape[1] +
         P_node_2->frequency * shape[2] + P_node_3->frequency * shape[3];

  return;
}

//
//  getSiteCurrentStateAndNodeInfo()
//
//  Note : this takes std::vector type data as input
//  This function writes position+freq dof on input data. it also
//  writes NodeInfo on input data
//
void Lattice::getSiteCurrentStateAndNodeInfo(
    std::pair<std::pair<std::vector<int>, std::vector<double>>,
              std::vector<int>> &iC_data,
    struct lattice_t *P_lattice, const int iQuasi) {
  // clear the data
  iC_data.second.clear();
  iC_data.second.resize(4);

  iC_data.first.second.clear();
  iC_data.first.second.resize(4);

  struct element_t *P_element;

  struct node_t *P_node_0;
  struct node_t *P_node_1;
  struct node_t *P_node_2;
  struct node_t *P_node_3;

  double shape[4];

  // get int[3] form of lattice
  int l[3];
  for (int dof = 0; dof < 3; dof++)
    l[dof] = iC_data.first.first[dof];

  // element
  P_element = Element::getInstance()->locateSiteElement(&(shape[0]), l,
                                                        P_lattice, iQuasi);

  P_node_0 = P_element->node[0];
  P_node_1 = P_element->node[1];
  P_node_2 = P_element->node[2];
  P_node_3 = P_element->node[3];

  // compute position+freq
  iC_data.first.second[0] =
      P_node_0->position[0] * shape[0] + P_node_1->position[0] * shape[1] +
      P_node_2->position[0] * shape[2] + P_node_3->position[0] * shape[3];

  iC_data.first.second[1] =
      P_node_0->position[1] * shape[0] + P_node_1->position[1] * shape[1] +
      P_node_2->position[1] * shape[2] + P_node_3->position[1] * shape[3];

  iC_data.first.second[2] =
      P_node_0->position[2] * shape[0] + P_node_1->position[2] * shape[1] +
      P_node_2->position[2] * shape[2] + P_node_3->position[2] * shape[3];

  iC_data.first.second[3] =
      P_node_0->frequency * shape[0] + P_node_1->frequency * shape[1] +
      P_node_2->frequency * shape[2] + P_node_3->frequency * shape[3];

  // get node numbers and add it to the iCNodeInfo
  iC_data.second[0] = P_node_0->number;
  iC_data.second[1] = P_node_1->number;
  iC_data.second[2] = P_node_2->number;
  iC_data.second[3] = P_node_3->number;
  return;
} // end of getSiteCurrentStateAndNodeInfo()

//
//  getSiteInitialState() : return (x,y,z,w,t)
//
void Lattice::getSiteCurrentStateNew(std::vector<double> &state,
                                     std::vector<int> site, const int iQuasi) {
  state.clear();
  state.resize(4);

  struct element_t *P_element;

  struct node_t *P_node_0;
  struct node_t *P_node_1;
  struct node_t *P_node_2;
  struct node_t *P_node_3;

  double shape[4];

  // get P_lattice
  struct lattice_t iLattice =
      Quasicontinua::getInstance()->getQuasi(iQuasi).getLattice();

  int l[3];

  for (int dof = 0; dof < 3; dof++)
    l[dof] = site[dof];

  P_element = Element::getInstance()->locateSiteElement(&(shape[0]), l,
                                                        &iLattice, iQuasi);

  P_node_0 = P_element->node[0];
  P_node_1 = P_element->node[1];
  P_node_2 = P_element->node[2];
  P_node_3 = P_element->node[3];

  // compute position
  state[0] =
      P_node_0->position[0] * shape[0] + P_node_1->position[0] * shape[1] +
      P_node_2->position[0] * shape[2] + P_node_3->position[0] * shape[3];

  state[1] =
      P_node_0->position[1] * shape[0] + P_node_1->position[1] * shape[1] +
      P_node_2->position[1] * shape[2] + P_node_3->position[1] * shape[3];

  state[2] =
      P_node_0->position[2] * shape[0] + P_node_1->position[2] * shape[1] +
      P_node_2->position[2] * shape[2] + P_node_3->position[2] * shape[3];

  state[3] = P_node_0->frequency * shape[0] + P_node_1->frequency * shape[1] +
             P_node_2->frequency * shape[2] + P_node_3->frequency * shape[3];

  return;
} // end of getSiteInitialStateNew()

//
// getShell()
//
struct shell_t *Lattice::getShell(int shell_number,
                                  const struct lattice_t *P_lattice) {
  struct shell_t *shell_return;
  if (shell_number > MAX_NUMBER_SHELLS) {
    d_print("Increase value of max_number of shells\n");
    exit(EXIT_FAILURE);
  }

  switch (P_lattice->type) {
  case FCC:
    d_print("FCC not implemented\n");
    exit(EXIT_FAILURE);
    break;

  case BCC:
    shell_return = &(d_bcc_shells[shell_number - 1]);
    return (shell_return);
    break;

  default:
    d_print("Not implemented\n");
    exit(EXIT_FAILURE);
  }

  d_print("Check lattice type\n");

  exit(EXIT_FAILURE);
}

//
//  initializeShells()
//
int Lattice::initializeShells(struct lattice_t lattice) {

  switch (lattice.type) {
  case FCC:
    // if fcc is not initialized already, initilize it
    if (d_initialize_fcc == 0) {
      //  counter for shells
      int i_shell;

      // generate shells up to maximum number
      for (i_shell = 1; i_shell <= MAX_NUMBER_SHELLS; i_shell++)
        if (GenerateNewShell(i_shell, lattice) == RETURN_ERR) {
          return (RETURN_ERR);
        }

      // set static variable that fcc has been initialized
      d_initialize_fcc = 1;

      return (RETURN_OK);
    }

    else {
      return (RETURN_OK);
    }
    break;

  case BCC:
    // if bcc is not initialized, initialize it
    if (d_initialize_bcc == 0) {
      int i_shell;

      // generate shells upto maximum number
      for (i_shell = 1; i_shell <= MAX_NUMBER_SHELLS; i_shell++)
        if (GenerateNewShell(i_shell, lattice) == RETURN_ERR) {
          return (RETURN_ERR);
        }

      // set static variable that bcc has been initialized
      d_initialize_bcc = 1;

      return (RETURN_OK);
    }

    else {
      return (RETURN_OK);
    }
    break;

  default:

    // initialize other shells if not done

    ERROR("lattice type not supported");
    exit(EXIT_FAILURE);
  } // end of switch()

  return (RETURN_ERR);
} // end of initializeShells()

//
//  getShellRadius()
//
double Lattice::getShellRadius(struct lattice_t *P_lattice, const int l[3],
                               const int iQuasi) {
  double X[3];

  // get initial coordinates
  getSiteInitialPosition(X, l, P_lattice, iQuasi);

  return (X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);
}

// private functions
//
//  GenerateNewShell() : used in initializeShells()
//
int Lattice::GenerateNewShell(int shell_number, struct lattice_t lattice) {

  // create shell_t to push into d_fcc_shells
  struct shell_t data_shell(SHELL_INITIALIZER);

  switch (lattice.type) {

  case FCC:
    // check if shell_number is already in d_fcc_shells
    if (shell_number <= d_fcc_shells.size())
      return (RETURN_OK);

    // fill sites into data_shell
    if (FillShell(&data_shell, shell_number, lattice) == RETURN_ERR)
      return (RETURN_ERR);

    // push back the data_shell to d_fcc_shells
    d_fcc_shells.push_back(data_shell);

    return (RETURN_OK);
    break;

  case BCC:
    // check if shell_number is already in d_bcc_shells
    if (shell_number <= d_bcc_shells.size())
      return (RETURN_OK);

    // fill sites into data_shell
    if (FillShell(&data_shell, shell_number, lattice) == RETURN_ERR)
      return (RETURN_ERR);

    // push back the data_shell to d_fcc_shells
    d_bcc_shells.push_back(data_shell);

    return (RETURN_OK);
    break;

  default:
    ERROR("lattice type not supported");
    exit(EXIT_FAILURE);
  }

  return (RETURN_ERR);
} // end of GenerateNewShell()

//
//  FillShell() : used in GenerateNewShells()
//
int Lattice::FillShell(struct shell_t *data_shell, int shell_number,
                       struct lattice_t lattice) {
  switch (lattice.type) {
  case FCC:
    do {
      int l_max;
      int l1;
      int l2;
      int l3;

      l_max = floor(sqrt(d_fcc_p));

      /**
        * loop over all integers inside [0,l_max]x[0,l_max]x[0,l_max]
        */

      for (l1 = 0; l1 <= l_max; l1++)
        for (l2 = 0; l2 <= l_max; l2++)
          for (l3 = 0; l3 <= l_max; l3++)
            if ((l1 * l1 + l2 * l2 + l3 * l3 == d_fcc_p) &&
                ((l1 + l2 + l3) % 2 == 0)) {
              // add (l1,l2,l3)
              AddNewSite(data_shell, l1, l2, l3);

              // l1
              if (l2 != 0)
                AddNewSite(data_shell, l1, -l2, l3);
              if (l3 != 0)
                AddNewSite(data_shell, l1, l2, -l3);
              if ((l2 != 0) && (l3 != 0))
                AddNewSite(data_shell, l1, -l2, -l3);

              // - l1
              if (l1 != 0) {
                AddNewSite(data_shell, -l1, l2, l3);

                if (l2 != 0)
                  AddNewSite(data_shell, -l1, -l2, l3);
                if (l3 != 0)
                  AddNewSite(data_shell, -l1, l2, -l3);
                if ((l2 != 0) && (l3 != 0))
                  AddNewSite(data_shell, -l1, -l2, -l3);
              }
            }

      d_fcc_p++;
    } while (data_shell->number_sites == 0);

    return (RETURN_OK);
    break;

  case BCC: {
    int x_count;
    int y_count;
    int z_count;

    // add +/- X faces
    for (y_count = -d_bcc_p; y_count <= d_bcc_p; ++y_count) {
      for (z_count = -d_bcc_p; z_count <= d_bcc_p; ++z_count) {
        // add faces
        if (AddNewSite(data_shell, -d_bcc_p, y_count, z_count) == RETURN_ERR)
          return (RETURN_ERR);

        if (AddNewSite(data_shell, d_bcc_p, y_count, z_count) == RETURN_ERR)
          return (RETURN_ERR);
      }
    }

    // add +/- Y faces
    for (x_count = -d_bcc_p + 1; x_count <= d_bcc_p - 1; ++x_count) {
      for (z_count = -d_bcc_p; z_count <= d_bcc_p; ++z_count) {
        // add faces
        if (AddNewSite(data_shell, x_count, -d_bcc_p, z_count) == RETURN_ERR)
          return (RETURN_ERR);

        if (AddNewSite(data_shell, x_count, d_bcc_p, z_count) == RETURN_ERR)
          return (RETURN_ERR);
      }
    }

    // add +/- Z faces
    for (x_count = -d_bcc_p + 1; x_count <= d_bcc_p - 1; ++x_count) {
      for (y_count = -d_bcc_p + 1; y_count <= d_bcc_p - 1; ++y_count) {
        // add faces
        if (AddNewSite(data_shell, x_count, y_count, -d_bcc_p) == RETURN_ERR)
          return (RETURN_ERR);

        if (AddNewSite(data_shell, x_count, y_count, d_bcc_p) == RETURN_ERR)
          return (RETURN_ERR);
      }
    }

    d_bcc_p++;

    return (RETURN_OK);
    break;
  }

  default:
    ERROR("lattice type not supported");
    exit(EXIT_FAILURE);
  }

  return (RETURN_ERR);
} // end of FillShell()

//
//  AddNewSite()
//
int Lattice::AddNewSite(struct shell_t *data_shell, int l1, int l2, int l3) {

  // check if this is the first site to be added to data_shell
  if (data_shell->number_sites == 0) {
    // allocate memory for first site to be added
    data_shell->site = (int(*)[3])malloc(sizeof(int[3]));
  } else {
    // this is not the first site, so increase space
    data_shell->site = (int(*)[3])realloc(
        data_shell->site, (data_shell->number_sites + 1) * sizeof(int[3]));
  }

  if (data_shell->site == NULL)
    return (RETURN_ERR);

  // add site to new allocated space
  data_shell->site[data_shell->number_sites][0] = l1;
  data_shell->site[data_shell->number_sites][1] = l2;
  data_shell->site[data_shell->number_sites][2] = l3;

  // increase the number_sites counter
  data_shell->number_sites++;

  // if(d_bcc_number_shells <= 4){
  //   std::cout<<"++++++"<<std::endl;
  //   std::cout<<l1<<" : "<<data_shell->site[data_shell->number_sites -
  //   1][0]<<std::endl;
  //   std::cout<<l2<<" : "<<data_shell->site[data_shell->number_sites
  //   -1][1]<<std::endl;
  //   std::cout<<l3<<" : "<<data_shell->site[data_shell->number_sites
  //   -1][2]<<std::endl;
  // }

  return (RETURN_OK);
}

//
//    GetLatticeCoordinates()
//
void Lattice::GetLatticeCoordinates(double r[3], double X[3],
                                    const struct lattice_t *P_lattice) {
  MiscFunctions *miscF = MiscFunctions::getInstance();
  double det =
      miscF->computeDeterminant(P_lattice->a1, P_lattice->a2, P_lattice->a3);

  r[0] = miscF->computeDeterminant(X, P_lattice->a2, P_lattice->a3) / det;

  r[1] = miscF->computeDeterminant(P_lattice->a1, X, P_lattice->a3) / det;

  r[2] = miscF->computeDeterminant(P_lattice->a1, P_lattice->a2, X) / det;

  return;
}
}