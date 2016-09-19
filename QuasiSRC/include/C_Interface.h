//
//  C_Interface.h
//
#if !defined(C_INTERFACE_H)
#define C_INTERFACE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#include <stdarg.h>
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "DataTypes.h"

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

void f_print(FILE *file, const char *data, ...);

void d_print(const char *data, ...);

void *free_e(void *p);

int isLatticeSite(const int l[3], struct lattice_t *p_lattice,
                  const int iQuasi);

void getSiteInitialPosition(double r[3], const int l[3],
                            struct lattice_t *p_lattice, const int iQuasi);

int checkPointInTetra(const double v1[3], const double v2[3],
                      const double v3[3], const double v4[3],
                      const double pos[3]);

double computeShapeFunction(const double p1[3], const double p2[3],
                            const double p3[3], const double p4[3],
                            const double p[3]);

#if defined(__cplusplus)
}
#endif /* __cplusplus */

#endif /* C_INTERFACE_H */