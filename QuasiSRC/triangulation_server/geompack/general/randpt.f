      subroutine randpt(k,n,seed,axis,nptav,scale,trans,lda,a)
      implicit logical (a-z)
      integer axis,k,lda,n,nptav,seed
      double precision scale(k),trans(k),a(lda,n)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: generate n random k-dimensional points from the uniform
c        distribution.
c
c     input parameters:
c        k - dimension of points
c        n - number of random points
c        seed - seed for pseudo random number generator
c        axis,nptav - if axis < 1 or > k, then uniform random points are
c              generated; if 1 <= axis <= k then an average of nptav
c              uniform random points are generated with the same axis
c              coordinate on about n/nptav random parallel hyperplanes
c        scale(1:k),trans(1:k) - scale and translation factors for
c              coordinates 1 to k; ith coordinate of random point is
c              r*scale(i) + trans(i) where 0 < r < 1.
c        lda - leading dimension of array a in calling routine; should
c              be >= k
c
c     updated parameters:
c        seed - gets updated on each call to urand
c
c     output parameters:
c        a(1:k,1:n) - array of n uniform random k-d points
c
c     routines called:
c        urand
c
      integer i,j,m
      double precision r
      real urand
c
      if (axis .lt. 1 .or. axis .gt. k) then
         do 20 j = 1,n
	    do 10 i = 1,k
	       a(i,j) = urand(seed)*scale(i) + trans(i)
   10       continue
   20    continue
      else
	 m = int(urand(seed)*2.0*nptav + 0.5)
	 r = urand(seed)*scale(axis) + trans(axis)
         do 40 j = 1,n
	    do 30 i = 1,k
	       if (i .eq. axis) then
		  a(i,j) = r
	       else
	          a(i,j) = urand(seed)*scale(i) + trans(i)
	       endif
   30       continue
	    m = m - 1
	    if (m .le. 0) then
	       m = int(urand(seed)*2.0*nptav + 0.5)
	       r = urand(seed)*scale(axis) + trans(axis)
	    endif
   40    continue
      endif
      end
