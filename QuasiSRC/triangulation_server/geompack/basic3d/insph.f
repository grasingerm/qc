      subroutine insph(a,b,c,d,centre,rad)
      implicit logical (a-z)
      double precision a(3),b(3),c(3),centre(3),d(3),rad
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: find centre and radius of insphere of tetrahedron.
c
c     input parameters:
c        a(1:3),b(1:3),c(1:3),d(1:3) - 4 vertices of tetrahedron
c
c     output parameters:
c        centre(1:3) - centre of insphere; undefined if a,b,c,d coplanar
c        rad - radius of insphere; 0 if a,b,c,d coplanar
c
      double precision pi,tol
      common /gconst/ pi,tol
      save /gconst/
c
      integer i,j,r
      double precision ab(3),ac(3),ad(3),bc(3),bd(3),mat(4,4),rtol,t
      logical singlr
c
c     compute unit outward (or inward) normals and equations of 4 faces.
c
      do 10 i = 1,3
	 ab(i) = b(i) - a(i)
	 ac(i) = c(i) - a(i)
	 ad(i) = d(i) - a(i)
	 bc(i) = c(i) - b(i)
	 bd(i) = d(i) - b(i)
   10 continue
      rtol = tol*max(abs(a(1)),abs(a(2)),abs(a(3)),abs(b(1)),abs(b(2)),
     $   abs(b(3)),abs(c(1)),abs(c(2)),abs(c(3)),abs(d(1)),abs(d(2)),
     $   abs(d(3)))
      mat(1,1) = ac(2)*ab(3) - ac(3)*ab(2)
      mat(1,2) = ac(3)*ab(1) - ac(1)*ab(3)
      mat(1,3) = ac(1)*ab(2) - ac(2)*ab(1)
      mat(2,1) = ab(2)*ad(3) - ab(3)*ad(2)
      mat(2,2) = ab(3)*ad(1) - ab(1)*ad(3)
      mat(2,3) = ab(1)*ad(2) - ab(2)*ad(1)
      mat(3,1) = ad(2)*ac(3) - ad(3)*ac(2)
      mat(3,2) = ad(3)*ac(1) - ad(1)*ac(3)
      mat(3,3) = ad(1)*ac(2) - ad(2)*ac(1)
      mat(4,1) = bc(2)*bd(3) - bc(3)*bd(2)
      mat(4,2) = bc(3)*bd(1) - bc(1)*bd(3)
      mat(4,3) = bc(1)*bd(2) - bc(2)*bd(1)
      singlr = .true.
      do 30 i = 4,1,-1
	 t = sqrt(mat(i,1)**2 + mat(i,2)**2 + mat(i,3)**2)
	 if (t .le. rtol) go to 100
	 mat(i,1) = mat(i,1)/t
	 mat(i,2) = mat(i,2)/t
	 mat(i,3) = mat(i,3)/t
	 if (i .eq. 4) then
	    mat(i,4) = mat(i,1)*b(1) + mat(i,2)*b(2) + mat(i,3)*b(3)
	 else
	    mat(i,4) = mat(i,1)*a(1) + mat(i,2)*a(2) + mat(i,3)*a(3)
	    do 20 j = 1,4
	       mat(i,j) = mat(i,j) - mat(4,j)
   20       continue
	 endif
   30 continue
c
c     use gaussian elimination with partial pivoting to solve 3 by 3
c     system of linear equations for centre of insphere.
c
      r = 1
      do 40 i = 2,3
	 if (abs(mat(i,1)) .gt. abs(mat(r,1))) r = i
   40 continue
      if (abs(mat(r,1)) .le. tol) go to 100
      if (r .ne. 1) then
	 do 50 j = 1,4
	    t = mat(1,j)
	    mat(1,j) = mat(r,j)
	    mat(r,j) = t
   50    continue
      endif
      do 70 i = 2,3
	 t = mat(i,1)/mat(1,1)
	 do 60 j = 2,4
	    mat(i,j) = mat(i,j) - t*mat(1,j)
   60    continue
   70 continue
      if (abs(mat(3,2)) .gt. abs(mat(2,2))) then
	 do 80 j = 2,4
	    t = mat(2,j)
	    mat(2,j) = mat(3,j)
	    mat(3,j) = t
   80    continue
      endif
      if (abs(mat(2,2)) .le. tol) go to 100
      t = mat(3,2)/mat(2,2)
      do 90 j = 3,4
	 mat(3,j) = mat(3,j) - t*mat(2,j)
   90 continue
      if (abs(mat(3,3)) .gt. tol) singlr = .false.
c
  100 continue
      if (singlr) then
	 rad = 0.0d0
      else
	 centre(3) = mat(3,4)/mat(3,3)
	 centre(2) = (mat(2,4) - mat(2,3)*centre(3))/mat(2,2)
	 centre(1) = (mat(1,4) - mat(1,3)*centre(3) -
     $      mat(1,2)*centre(2))/mat(1,1)
	 rad = abs(mat(4,1)*centre(1) + mat(4,2)*centre(2) +
     $      mat(4,3)*centre(3) - mat(4,4))
      endif
      end
