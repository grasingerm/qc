//
// C_Interface.cc
//

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "C_Interface.h"
#include "DataTypes.h"
#include "Error.h"
#include "Lattice.h"
#include "MiscFunctions.h"
#include "Output.h"
#include "Quasicontinua.h"
#include "Quasicontinuum.h"
#include "Shape.h"

//
//  f_print()
//
void f_print(FILE *file, const char *data, ...) {
  va_list f_vargs;

  va_start(f_vargs, data);

  // print only if file is not empty
  if (file != NULL) {
    vfprintf(file, data, f_vargs);
  }

  va_end(f_vargs);

  va_list vargs;
  va_start(vargs, data);
  vprintf(data, vargs);
  va_end(vargs);

  return;
}

//
//  d_print()
//
void d_print(const char *data, ...) {
  //
  //  get FILE from Output
  //
  FILE *file;
  file = quasicontinuum::Output::getInstance()->d_debugFile;

  va_list f_vargs;

  va_start(f_vargs, data);

  // print only if file is not empty
  if (file != NULL) {
    vfprintf(file, data, f_vargs);
  }

  va_end(f_vargs);

  va_list vargs;
  va_start(vargs, data);
  vprintf(data, vargs);
  va_end(vargs);

  return;
}

//
//  free_e()
//
void *free_e(void *p) {
  free(p);

  return ((void *)NULL);
}

//
//	isLatticeSite()
//
int isLatticeSite(const int l[3], struct lattice_t *P_lattice,
                  const int iQuasi) {
  struct quasicontinuum::lattice_t *pLattice =
      reinterpret_cast<quasicontinuum::lattice_t *>(P_lattice);
  // link it to relevant Lattice.cc function
  return (quasicontinuum::Lattice::getInstance()->isLatticeSite(l, pLattice,
                                                                iQuasi));
}

//
//	getSiteInitialPosition()
//
void getSiteInitialPosition(double r[3], const int l[3],
                            struct lattice_t *P_lattice, const int iQuasi) {
  struct quasicontinuum::lattice_t *pLattice =
      reinterpret_cast<quasicontinuum::lattice_t *>(P_lattice);
  // link it to relevant Lattice.cc function
  quasicontinuum::Lattice::getInstance()->getSiteInitialPosition(r, l, pLattice,
                                                                 iQuasi);

  return;
}

//
//	checkPointInTetra()
//
//  Return = 0 : Outside
//           1 : Inside
//
int checkPointInTetra(const double v1[3], const double v2[3],
                      const double v3[3], const double v4[3],
                      const double pos[3]) {
  if (quasicontinuum::MiscFunctions::getInstance()->checkPointInTetra(
          v1, v2, v3, v4, pos) == quasicontinuum::INSIDE)
    return 1;
  else {
    if (quasicontinuum::MiscFunctions::getInstance()->checkPointInTetra(
            v1, v2, v3, v4, pos) == quasicontinuum::OUTSIDE)
      return 0;
    else {
      ERROR("checkPointInTetra()");
      exit(EXIT_FAILURE);
    }
  }

  // return(quasicontinuum::MiscFunctions::getInstance()->checkPointInTetra(v1,
  // v2, v3, v4, pos) );
}

//
//	computeShapeFunction()
//
double computeShapeFunction(const double p1[3], const double p2[3],
                            const double p3[3], const double p4[3],
                            const double p[3]) {
  // link it to Shape.cc function
  return (quasicontinuum::Shape::getInstance()->computeShapeFunction(p1, p2, p3,
                                                                     p4, p));
}
