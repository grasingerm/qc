      subroutine baryth(a,b,c,d,e,alpha,degen)
      implicit logical (a-z)
      double precision a(3),alpha(4),b(3),c(3),d(3),e(3)
      logical degen
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: compute barycentric coordinates of 3-d point with respect
c        to four vertices of a tetrahedron.
c
c     input parameters:
c        a(1:3),b(1:3),c(1:3),d(1:3) - 4 vertices of tetrahedron
c        e(1:3) - fifth point for which barycentric coordinates found
c
c     output parameters:
c        alpha(1:4) - scaled barycentric coord. (if degen = .false.)
c              such that e = (alpha(1)*a + alpha(2)*b + alpha(3)*c +
c              alpha(4)*d)/det where det = 6 * (volume of tetra abcd);
c              an alpha(i) may be set to 0 after tolerance test to
c              indicate that e is coplanar with a face, so sum of
c              alpha(i)/det may not be 1; if the actual barycentric
c              coordinates rather than just their signs are needed,
c              modify this routine to divide alpha(i) by det
c        degen - .true. iff a,b,c,d coplanar
c
      double precision pi,tol
      common /gconst/ pi,tol
      save /gconst/
c
      integer i
      double precision amax,bmax,cmax,cp1,cp2,cp3,det,dmax,emax
      double precision da(3),db(3),dc(3),de(3),ea(3),eb(3),ec(3)
c
      degen = .false.
      do 10 i = 1,3
	 da(i) = a(i) - d(i)
	 db(i) = b(i) - d(i)
	 dc(i) = c(i) - d(i)
   10 continue
      amax = max(abs(a(1)),abs(a(2)),abs(a(3)))
      bmax = max(abs(b(1)),abs(b(2)),abs(b(3)))
      cmax = max(abs(c(1)),abs(c(2)),abs(c(3)))
      dmax = max(abs(d(1)),abs(d(2)),abs(d(3)))
      cp1 = db(2)*dc(3) - db(3)*dc(2)
      cp2 = db(3)*dc(1) - db(1)*dc(3)
      cp3 = db(1)*dc(2) - db(2)*dc(1)
      det = da(1)*cp1 + da(2)*cp2 + da(3)*cp3
c     if (abs(det) .le. 0.01d0*tol*max(amax,bmax,cmax,dmax)) then
      if (abs(det) .le. tol*max(amax,bmax,cmax,dmax)) then
         degen = .true.
	 return
      endif
      do 20 i = 1,3
	 de(i) = e(i) - d(i)
	 ea(i) = a(i) - e(i)
	 eb(i) = b(i) - e(i)
	 ec(i) = c(i) - e(i)
   20 continue
      alpha(1) = de(1)*cp1 + de(2)*cp2 + de(3)*cp3
      cp1 = da(2)*de(3) - da(3)*de(2)
      cp2 = da(3)*de(1) - da(1)*de(3)
      cp3 = da(1)*de(2) - da(2)*de(1)
      alpha(2) = dc(1)*cp1 + dc(2)*cp2 + dc(3)*cp3
      alpha(3) = db(1)*cp1 + db(2)*cp2 + db(3)*cp3
      alpha(4) = ea(1)*(eb(2)*ec(3) - eb(3)*ec(2)) + ea(2)*(eb(3)*ec(1)
     $   - eb(1)*ec(3)) + ea(3)*(eb(1)*ec(2) - eb(2)*ec(1))
      if (det .lt. 0.0d0) then
	 alpha(1) = -alpha(1)
	 alpha(2) = -alpha(2)
	 alpha(4) = -alpha(4)
      else
	 alpha(3) = -alpha(3)
      endif
      emax = max(abs(e(1)),abs(e(2)),abs(e(3)))
      if (abs(alpha(1)) .le. tol*max(bmax,cmax,dmax,emax))
     $   alpha(1) = 0.0d0
      if (abs(alpha(2)) .le. tol*max(amax,cmax,dmax,emax))
     $   alpha(2) = 0.0d0
      if (abs(alpha(3)) .le. tol*max(amax,bmax,dmax,emax))
     $   alpha(3) = 0.0d0
      if (abs(alpha(4)) .le. tol*max(amax,bmax,cmax,emax))
     $   alpha(4) = 0.0d0
      end
