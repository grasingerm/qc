      subroutine ccsph(intest,a,b,c,d,e,centre,radsq,in)
      implicit logical (a-z)
      logical intest
      integer in
      double precision a(3),b(3),c(3),centre(3),d(3),e(3),radsq
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: find centre and square of radius of circumsphere through
c        four vertices of a tetrahedron, and possibly determine whether
c        a fifth 3-d point is inside sphere.
c
c     input parameters:
c        intest - .true. iff test for fifth point in sphere to be made
c        a(1:3),b(1:3),c(1:3),d(1:3) - 4 vertices of tetrahedron
c        e(1:3) - fifth point; referenced iff intest is .true.
c
c     output parameters:
c        centre(1:3) - centre of sphere; undefined if a,b,c,d coplanar
c        radsq - square of radius of sphere; -1 if a,b,c,d coplanar
c        in - contains following value if intest is .true.:
c              2 if a,b,c,d coplanar; 1 if e inside sphere;
c              0 if e on sphere; -1 if e outside sphere
c
      double precision pi,tol
      common /gconst/ pi,tol
      save /gconst/
c
      integer i
      double precision cmax,cp1,cp2,cp3,det,dsq
      double precision da(3),db(3),dc(3),rhs(3)
c
      do 10 i = 1,3
	 da(i) = a(i) - d(i)
	 db(i) = b(i) - d(i)
	 dc(i) = c(i) - d(i)
   10 continue
      rhs(1) = 0.5d0*(da(1)**2 + da(2)**2 + da(3)**2)
      rhs(2) = 0.5d0*(db(1)**2 + db(2)**2 + db(3)**2)
      rhs(3) = 0.5d0*(dc(1)**2 + dc(2)**2 + dc(3)**2)
      cmax = max(abs(a(1)),abs(a(2)),abs(a(3)),abs(b(1)),abs(b(2)),
     $   abs(b(3)),abs(c(1)),abs(c(2)),abs(c(3)),abs(d(1)),abs(d(2)),
     $   abs(d(3)))
      cp1 = db(2)*dc(3) - dc(2)*db(3)
      cp2 = dc(2)*da(3) - da(2)*dc(3)
      cp3 = da(2)*db(3) - db(2)*da(3)
      det = da(1)*cp1 + db(1)*cp2 + dc(1)*cp3
c     if (abs(det) .le. 0.01d0*tol*cmax) then
      if (abs(det) .le. tol*cmax) then
	 radsq = -1.0d0
	 in = 2
	 return
      endif
      centre(1) = (rhs(1)*cp1 + rhs(2)*cp2 + rhs(3)*cp3)/det
      cp1 = db(1)*rhs(3) - dc(1)*rhs(2)
      cp2 = dc(1)*rhs(1) - da(1)*rhs(3)
      cp3 = da(1)*rhs(2) - db(1)*rhs(1)
      centre(2) = (da(3)*cp1 + db(3)*cp2 + dc(3)*cp3)/det
      centre(3) = -(da(2)*cp1 + db(2)*cp2 + dc(2)*cp3)/det
      radsq = centre(1)**2 + centre(2)**2 + centre(3)**2
      do 20 i = 1,3
	 centre(i) = centre(i) + d(i)
   20 continue
      if (intest) then
         dsq = (e(1) - centre(1))**2 + (e(2) - centre(2))**2 +
     $      (e(3) - centre(3))**2
         if (dsq .gt. (1.0d0 + tol)*radsq) then
	    in = -1
         else if (dsq .lt. (1.0d0 - tol)*radsq) then
	    in = 1
         else
	    in = 0
         endif
      endif
      end
