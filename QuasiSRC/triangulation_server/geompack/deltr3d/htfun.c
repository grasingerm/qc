/**
 * $Log: htfun.c,v $
 * Revision 1.2  2002/03/07 23:38:52  knap
 * Updated build procedure to use automake/autoconf.
 *
 * Revision 1.1  2000/07/19 18:50:04  knap
 * Initial rev.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */

#if !defined(lint)
static const char rcsid[] = "$Id: htfun.c,v 1.2 2002/03/07 23:38:52 knap Exp $";
#endif /* !lint */

#ifdef __cplusplus
  extern "C" int HTFUN_F77(const int *, const int *, const int *,
                           const int *, const int *);
#else
  int HTFUN_F77(const int *, const int *, const int *,
                const int *, const int *);
#endif

int
HTFUN_F77(const int *A, 
       const int *B, 
       const int *C, 
       const int *N, 
       const int *P)
{

  long k;
  long a;
  long b;
  long c;
  long n;
  long p;

  a = *A;
  b = *B;
  c = *C;
  n = *N;
  p = *P;
  
  k = ((a)*(n) + (b))%(p);
  return((int) ((k*n + c)%p));

}
