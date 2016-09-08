      subroutine isftdw(l,u,k,lda,a,map)
      implicit logical (a-z)
      integer k,l,lda,u
      integer a(lda,*),map(*)
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
c        a(1:k,1:*),map(1:*) - see routine ihpsrt
c
c     updated parameters:
c        map(1:*) - see routine ihpsrt
c
c     routines called:
c        iless
c
      integer i,j,t
      logical iless
c
      i = l
      j = 2*i
      t = map(i)
   10 continue
      if (j .gt. u) go to 20
	 if (j .lt. u) then
	    if (iless(k,a(1,map(j)),a(1,map(j+1)))) j = j + 1
	 endif
	 if (iless(k,a(1,map(j)),a(1,t))) go to 20
	 map(i) = map(j)
	 i = j
	 j = 2*i
      go to 10
   20 continue
      map(i) = t
      end
