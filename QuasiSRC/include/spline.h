#ifndef _SPLINE_H_
#define _SPLINE_H_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#if defined(__cplusplus)
//namespace quasicontinuum {
  extern "C" {
#endif /* __cplusplus */

/**
 * Fortran prototypes
 */

#define SPLINE_F77 F77_FUNC(spline, SPLINE)
extern void SPLINE_F77(const int *, const  double *, const  double *, const  double *, const  double *,
		       const double *, const double *, double *, double *);

#if defined(__cplusplus)
  }
//}
#endif /* __cplusplus */

#endif /* _DENSITY_H_ */