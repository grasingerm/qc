      subroutine rotiar(n,arr,shift)
      implicit logical (a-z)
      integer n,shift
      integer arr(0:n-1)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: rotate elements of integer array.
c
c     input parameters:
c	 n - number of elements of array
c	 arr(0:n-1) - integer array
c        shift - amount of (left) shift or rotation; arr(shift) on input
c              becomes arr(0) on output
c
c     updated parameters:
c	 arr(0:n-1) - rotated integer array
c
      integer a,b,i,j,k,l,m,r,sh,t
c
      sh = mod(shift,n)
      if (sh .lt. 0) sh = sh + n
      if (sh .eq. 0) return
      a = n
      b = sh
   20 continue
	 r = mod(a,b)
	 a = b
	 b = r
      if (r .gt. 0) go to 20
      m = n/a - 1
      do 40 i = 0,a-1
	 t = arr(i)
	 k = i
	 do 30 j = 1,m
	    l = k + sh
	    if (l .ge. n) l = l - n
	    arr(k) = arr(l)
	    k = l
   30    continue
	 arr(k) = t
   40 continue
      end
