      subroutine dsftdw(l,u,k,lda,a,map)
      implicit logical (a-z)
      integer k,l,lda,map(*),u
      double precision a(lda,*)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: sift a(*,map(l)) down a heap of size u.
c
c     input parameters:
c        l,u - lower and upper index of part of heap
c        k - dimension of points
c        lda - leading dimension of array a in calling routine
c        a(1:k,1:*),map(1:*) - see routine dhpsrt
c
c     updated parameters:
c        map(1:*) - see routine dhpsrt
c
c     routines called:
c        dless
c
      integer i,j,t
      logical dless
c
      i = l
      j = 2*i
      t = map(i)
   10 continue
      if (j .gt. u) go to 20
	 if (j .lt. u) then
	    if (dless(k,a(1,map(j)),a(1,map(j+1)))) j = j + 1
	 endif
	 if (dless(k,a(1,map(j)),a(1,t))) go to 20
	 map(i) = map(j)
	 i = j
	 j = 2*i
      go to 10
   20 continue
      map(i) = t
      end
