      integer function opside(a,b,c,d,e)
      implicit logical (a-z)
      double precision a(3),b(3),c(3),d(3),e(3)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: test if points d, e are on opposite sides of triangular
c              face with vertices a, b, c.
c
c     input parameters:
c        a(1:3),b(1:3),c(1:3),d(1:3),e(1:3) - five 3-d points
c
c     returned function value:
c        opside - +1 if d, e on opposite sides; -1 if on same side;
c                 2 if d is coplanar with face abc (abcd is degenerate
c                 tetra); 0 if e is coplanar with face abc
c
      double precision pi,tol
      common /gconst/ pi,tol
      save /gconst/
c
      integer i
      double precision ab(3),ac(3),ddp,dmax,edp,emax,nrml1,nrml2,nrml3
c
      do 10 i = 1,3
	 ab(i) = b(i) - a(i)
	 ac(i) = c(i) - a(i)
   10 continue
      emax = max(abs(a(1)),abs(a(2)),abs(a(3)),abs(b(1)),abs(b(2)),
     $   abs(b(3)),abs(c(1)),abs(c(2)),abs(c(3)))
      dmax = max(emax, abs(d(1)),abs(d(2)),abs(d(3)))
      nrml1 = ab(2)*ac(3) - ab(3)*ac(2)
      nrml2 = ab(3)*ac(1) - ab(1)*ac(3)
      nrml3 = ab(1)*ac(2) - ab(2)*ac(1)
      ddp = (d(1) - a(1))*nrml1 + (d(2) - a(2))*nrml2 +
     $   (d(3) - a(3))*nrml3
      if (abs(ddp) .le. tol*dmax) then
	 opside = 2
	 return
      endif
      emax = max(emax, abs(e(1)),abs(e(2)),abs(e(3)))
      edp = (e(1) - a(1))*nrml1 + (e(2) - a(2))*nrml2 +
     $   (e(3) - a(3))*nrml3
      if (abs(edp) .le. tol*emax) then
	 opside = 0
      else if (ddp*edp .lt. 0.0d0) then
	 opside = 1
      else
	 opside = -1
      endif
      end
