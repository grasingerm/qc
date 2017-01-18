#ifndef _SPLINE_H_
#define _SPLINE_H_

#if defined(__cplusplus)
  extern "C" void SPLINE_F77(const int *, const double *, const double *,
                             const double *, const double *, const double *,
                             const double *, double *, double *);
#else
  fortran void SPLINE_F77(const int *, const double *, const double *,
                          const double *, const double *, const double *,
                          const double *, double *, double *);
#endif /* __cplusplus */

#endif
