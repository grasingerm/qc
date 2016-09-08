      double precision function radrth(a,b,c,d)
      implicit logical (a-z)
      double precision a(3),b(3),c(3),d(3)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: compute radius (aspect) ratio of tetrahedron
c        = 3 * inradius / circumradius
c
c     input parameters:
c        a(1:3),b(1:3),c(1:3),d(1:3) - 4 vertices of tetrahedron
c
c     returned function value:
c        radrth - radius ratio of tetrahedron
c
      integer i
      double precision ab(3),ac(3),ad(3),bc(3),bd(3),cd(3)
      double precision cp1,cp2,cp3,denom,fa,fb,fc,fd
      double precision lab,lac,lad,lbc,lbd,lcd,pb,pc,pd,t1,t2,vol
c
      do 10 i = 1,3
	 ab(i) = b(i) - a(i)
	 ac(i) = c(i) - a(i)
	 ad(i) = d(i) - a(i)
	 bc(i) = c(i) - b(i)
	 bd(i) = d(i) - b(i)
	 cd(i) = d(i) - c(i)
   10 continue
      lab = ab(1)**2 + ab(2)**2 + ab(3)**2
      lac = ac(1)**2 + ac(2)**2 + ac(3)**2
      lad = ad(1)**2 + ad(2)**2 + ad(3)**2
      lbc = bc(1)**2 + bc(2)**2 + bc(3)**2
      lbd = bd(1)**2 + bd(2)**2 + bd(3)**2
      lcd = cd(1)**2 + cd(2)**2 + cd(3)**2
      pb = sqrt(lab*lcd)
      pc = sqrt(lac*lbd)
      pd = sqrt(lad*lbc)
      cp1 = ab(2)*ac(3) - ab(3)*ac(2)
      cp2 = ab(3)*ac(1) - ab(1)*ac(3)
      cp3 = ab(1)*ac(2) - ab(2)*ac(1)
      fd = sqrt(cp1**2 + cp2**2 + cp3**2)
      cp1 = ab(2)*ad(3) - ab(3)*ad(2)
      cp2 = ab(3)*ad(1) - ab(1)*ad(3)
      cp3 = ab(1)*ad(2) - ab(2)*ad(1)
      fc = sqrt(cp1**2 + cp2**2 + cp3**2)
      cp1 = bc(2)*bd(3) - bc(3)*bd(2)
      cp2 = bc(3)*bd(1) - bc(1)*bd(3)
      cp3 = bc(1)*bd(2) - bc(2)*bd(1)
      fa = sqrt(cp1**2 + cp2**2 + cp3**2)
      cp1 = ac(2)*ad(3) - ac(3)*ad(2)
      cp2 = ac(3)*ad(1) - ac(1)*ad(3)
      cp3 = ac(1)*ad(2) - ac(2)*ad(1)
      fb = sqrt(cp1**2 + cp2**2 + cp3**2)
      t1 = pb + pc
      t2 = pb - pc
      denom = (fa+fb+fc+fd)*sqrt(abs((t1+pd)*(t1-pd)*(pd+t2)*(pd-t2)))
      if (denom .eq. 0.0d0) then
	 radrth = 0.0d0
      else
         vol = ab(1)*cp1 + ab(2)*cp2 + ab(3)*cp3
         radrth = 12.0d0*vol**2/denom
      endif
      end
