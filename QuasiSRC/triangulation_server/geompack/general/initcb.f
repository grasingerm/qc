      subroutine initcb(tolin)
      implicit logical (a-z)
      double precision tolin
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: initialize global variables in common blocks
c        gerror, gconst, and gprint. the latter is used for
c        printing debugging information.
c
c     input parameters:
c	 tolin - relative tolerance used to determine tol
c
c     output parameters in common blocks:
c        ierr - error code, initialized to 0
c        pi - acos(-1.0d0)
c        tol - relative tolerance max(tolin,100.0d0*eps) where
c              eps is approximation to machine epsilon
c        iprt - standard output unit 6
c        msglvl - message level, initialized to 0
c
      integer ierr,iprt,msglvl
      double precision pi,tol
      common /gerror/ ierr
      common /gconst/ pi,tol
      common /gprint/ iprt,msglvl
      save /gerror/,/gconst/,/gprint/
c
      double precision eps,epsp1
c
      ierr = 0
      pi = acos(-1.0d0)
      eps = 1.0d0
   10 continue
	 eps = eps/2.0d0
	 epsp1 = 1.0d0 + eps
      if (epsp1 .gt. 1.0d0) go to 10
      tol = max(tolin,100.0d0*eps)
      iprt = 6
      msglvl = 0
      end
