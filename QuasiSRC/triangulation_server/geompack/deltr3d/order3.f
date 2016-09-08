      subroutine order3(i,j,k)
      implicit logical (a-z)
      integer i,j,k
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: order i, j, k so that i <= j <= k.
c
c     input parameters:
c        i,j,k - 3 integers
c
c     output parameters:
c        i,j,k - sorted in nondecreasing order
c
      integer t
c
      if (j .lt. i) then
	 if (k .lt. j) then
	    t = i
	    i = k
	    k = t
	 else if (k .lt. i) then
	    t = i
	    i = j
	    j = k
	    k = t
	 else
	    t = i
	    i = j
	    j = t
	 endif
      else
	 if (k .lt. i) then
	    t = i
	    i = k
	    k = j
	    j = t
	 else if (k .lt. j) then
	    t = j
	    j = k
	    k = t
	 endif
      endif
      end
