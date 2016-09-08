      subroutine bnsrt3(binexp,n,a,map,bin,iwk)
      implicit logical (a-z)
      integer n
      integer bin(n),iwk(n),map(n)
      double precision a(3,*),binexp
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: use bin sort to obtain the permutation of n 3-dimensional
c        double precision points so that points are in increasing bin
c        order, where the n points are assigned to about n**binexp bins.
c
c     input parameters:
c        binexp - exponent for number of bins
c        n - number of points
c        a(1:3,1:*) - array of >= n 3-d double precision points
c        map(1:n) - the points of a with indices map(1), map(2), ...,
c              map(n) are to be sorted
c
c     updated parameters:
c        map(1:n) - elements are permuted so that bin of map(1) <=
c              bin of map(2) <= ... <= bin of map(n)
c
c     working parameters:
c        bin(1:n) - used for bin numbers and permutation of 1 to n
c        iwk(1:n) - used for copy of map array
c
c     routines called:
c        ihpsrt
c
      integer i,j,k,l,m,nside,nsidsq
      double precision dx,dy,dz,xmax,xmin,ymax,ymin,zmax,zmin
c
      nside = int(real(n)**(binexp/3.0) + 0.5)
      if (nside .le. 1) return
      nsidsq = nside**2
      xmin = a(1,map(1))
      ymin = a(2,map(1))
      zmin = a(3,map(1))
      xmax = xmin
      ymax = ymin
      zmax = zmin
      do 10 i = 2,n
         j = map(i)
	 xmin = min(xmin,a(1,j))
	 xmax = max(xmax,a(1,j))
	 ymin = min(ymin,a(2,j))
	 ymax = max(ymax,a(2,j))
	 zmin = min(zmin,a(3,j))
	 zmax = max(zmax,a(3,j))
   10 continue
      dx = 1.0001d0*(xmax - xmin)/dble(nside)
      dy = 1.0001d0*(ymax - ymin)/dble(nside)
      dz = 1.0001d0*(zmax - zmin)/dble(nside)
      if (dx .eq. 0.0d0) dx = 1.0d0
      if (dy .eq. 0.0d0) dy = 1.0d0
      if (dz .eq. 0.0d0) dz = 1.0d0
      do 20 i = 1,n
	 j = map(i)
	 iwk(i) = j
	 map(i) = i
	 k = int((a(1,j) - xmin)/dx)
	 l = int((a(2,j) - ymin)/dy)
	 m = int((a(3,j) - zmin)/dz)
	 if (mod(l,2) .eq. 0) then
	    bin(i) = l*nside + m
	 else
	    bin(i) = (l+1)*nside - m - 1
	 endif
	 if (mod(k,2) .eq. 0) then
	    bin(i) = k*nsidsq + bin(i)
	 else
	    bin(i) = (k+1)*nsidsq - bin(i) - 1
	 endif
   20 continue
      call ihpsrt(1,n,1,bin,map)
      do 30 i = 1,n
	 bin(i) = map(i)
   30 continue
      do 40 i = 1,n
	 map(i) = iwk(bin(i))
   40 continue
      end
