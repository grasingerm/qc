      double precision function sangmn(a,b,c,d,sang)
      implicit logical (a-z)
      double precision a(3),b(3),c(3),d(3),sang(4)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: compute 4 solid angles of tetrahedron and their minimum.
c        actually, sin(solid angle / 2) is computed.
c
c     input parameters:
c        a(1:3),b(1:3),c(1:3),d(1:3) - 4 vertices of tetrahedron
c
c     output parameters:
c        sang(1:4) - 4 solid angles at a,b,c,d, resp (actually sin(a/2))
c
c     returned function value:
c        sangmn - sin(min solid angle / 2) = min(sang(1:4))
c              since sum of 4 solid angles <= 2*pi
c
      integer i
      double precision ab(3),ac(3),ad(3),bc(3),bd(3),cd(3)
      double precision denom,l1,l2,l3,lab,lac,lad,lbc,lbd,lcd,vol
c
      do 10 i = 1,3
	 ab(i) = b(i) - a(i)
	 ac(i) = c(i) - a(i)
	 ad(i) = d(i) - a(i)
	 bc(i) = c(i) - b(i)
	 bd(i) = d(i) - b(i)
	 cd(i) = d(i) - c(i)
   10 continue
      lab = sqrt(ab(1)**2 + ab(2)**2 + ab(3)**2)
      lac = sqrt(ac(1)**2 + ac(2)**2 + ac(3)**2)
      lad = sqrt(ad(1)**2 + ad(2)**2 + ad(3)**2)
      lbc = sqrt(bc(1)**2 + bc(2)**2 + bc(3)**2)
      lbd = sqrt(bd(1)**2 + bd(2)**2 + bd(3)**2)
      lcd = sqrt(cd(1)**2 + cd(2)**2 + cd(3)**2)
      vol = abs(ab(1)*(ac(2)*ad(3) - ac(3)*ad(2)) + ab(2)*(ac(3)*ad(1)
     $   - ac(1)*ad(3)) + ab(3)*(ac(1)*ad(2) - ac(2)*ad(1)))*2.0d0
      l1 = lab + lac
      l2 = lab + lad
      l3 = lac + lad
      denom = (l1+lbc)*(l1-lbc)*(l2+lbd)*(l2-lbd)*(l3+lcd)*(l3-lcd)
      if (denom .le. 0.0d0) then
	 sang(1) = 0.0d0
      else
         sang(1) = vol/sqrt(denom)
      endif
      l1 = lab + lbc
      l2 = lab + lbd
      l3 = lbc + lbd
      denom = (l1+lac)*(l1-lac)*(l2+lad)*(l2-lad)*(l3+lcd)*(l3-lcd)
      if (denom .le. 0.0d0) then
	 sang(2) = 0.0d0
      else
         sang(2) = vol/sqrt(denom)
      endif
      l1 = lac + lbc
      l2 = lac + lcd
      l3 = lbc + lcd
      denom = (l1+lab)*(l1-lab)*(l2+lad)*(l2-lad)*(l3+lbd)*(l3-lbd)
      if (denom .le. 0.0d0) then
	 sang(3) = 0.0d0
      else
         sang(3) = vol/sqrt(denom)
      endif
      l1 = lad + lbd
      l2 = lad + lcd
      l3 = lbd + lcd
      denom = (l1+lab)*(l1-lab)*(l2+lac)*(l2-lac)*(l3+lbc)*(l3-lbc)
      if (denom .le. 0.0d0) then
	 sang(4) = 0.0d0
      else
         sang(4) = vol/sqrt(denom)
      endif
      sangmn = min(sang(1),sang(2),sang(3),sang(4))
      end
