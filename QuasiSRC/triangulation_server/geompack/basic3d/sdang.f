      subroutine sdang(a,b,c,d,sang,dang)
      implicit logical (a-z)
      double precision a(3),b(3),c(3),d(3)
      double precision sang(4),dang(6)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: compute 4 solid angles (+ pi) and 6 dihedral angles of
c        tetrahedron in radians. solid angle at vertex a is
c           (dihedral angle at edge ab) + (dihedral angle at edge ac)
c           + (dihedral angle at edge ad) - pi
c
c     input parameters:
c        a(1:3),b(1:3),c(1:3),d(1:3) - 4 vertices of tetrahedron
c
c     output parameters:
c        sang(1:4) - 4 solid angles (+ pi) at a,b,c,d, respectively
c        dang(1:6) - 6 dihedral angles at ab,ac,ad,bc,bd,cd, resp.
c
c     routines called:
c        angle3
c
      double precision pi,tol
      common /gconst/ pi,tol
      save /gconst/
c
      integer i
      double precision angle3,ab(3),ac(3),ad(3),bc(3),bd(3)
      double precision nrmabc(3),nrmabd(3),nrmacd(3),nrmbcd(3),rtolsq
c
c     compute outward (or inward) normals at 4 faces to get angles.
c
      do 10 i = 1,3
	 ab(i) = b(i) - a(i)
	 ac(i) = c(i) - a(i)
	 ad(i) = d(i) - a(i)
	 bc(i) = c(i) - b(i)
	 bd(i) = d(i) - b(i)
   10 continue
      rtolsq = tol*max(abs(a(1)),abs(a(2)),abs(a(3)),abs(b(1)),
     $   abs(b(2)),abs(b(3)),abs(c(1)),abs(c(2)),abs(c(3)),abs(d(1)),
     $   abs(d(2)),abs(d(3)))
      rtolsq = rtolsq**2
      nrmabc(1) = ac(2)*ab(3) - ac(3)*ab(2)
      nrmabc(2) = ac(3)*ab(1) - ac(1)*ab(3)
      nrmabc(3) = ac(1)*ab(2) - ac(2)*ab(1)
      nrmabd(1) = ab(2)*ad(3) - ab(3)*ad(2)
      nrmabd(2) = ab(3)*ad(1) - ab(1)*ad(3)
      nrmabd(3) = ab(1)*ad(2) - ab(2)*ad(1)
      nrmacd(1) = ad(2)*ac(3) - ad(3)*ac(2)
      nrmacd(2) = ad(3)*ac(1) - ad(1)*ac(3)
      nrmacd(3) = ad(1)*ac(2) - ad(2)*ac(1)
      nrmbcd(1) = bc(2)*bd(3) - bc(3)*bd(2)
      nrmbcd(2) = bc(3)*bd(1) - bc(1)*bd(3)
      nrmbcd(3) = bc(1)*bd(2) - bc(2)*bd(1)
      dang(1) = pi - angle3(nrmabc,nrmabd,rtolsq)
      dang(2) = pi - angle3(nrmabc,nrmacd,rtolsq)
      dang(3) = pi - angle3(nrmabd,nrmacd,rtolsq)
      dang(4) = pi - angle3(nrmabc,nrmbcd,rtolsq)
      dang(5) = pi - angle3(nrmabd,nrmbcd,rtolsq)
      dang(6) = pi - angle3(nrmacd,nrmbcd,rtolsq)
      sang(1) = dang(1) + dang(2) + dang(3)
      sang(2) = dang(1) + dang(4) + dang(5)
      sang(3) = dang(2) + dang(4) + dang(6)
      sang(4) = dang(3) + dang(5) + dang(6)
      end
