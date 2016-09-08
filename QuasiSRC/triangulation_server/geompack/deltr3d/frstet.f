      subroutine frstet(shift,nv,vcl,map,i3,i4)
      implicit logical (a-z)
      logical shift
      integer i3,i4,nv
      integer map(nv)
      double precision vcl(3,*)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: shift or swap vertices if necessary so first 3 vertices
c        (according to map) are not collinear and first 4 vertices are
c        not coplanar (so that first tetrahedron is valid).
c
c     input parameters:
c        shift - if .true., map(3), map(4) may be updated due to shift,
c              else they may be updated due to swaps; in former case,
c              it is assumed map gives vertices in lexicographic order
c        nv - number of vertices
c        vcl(1:3,1:*) - vertex coordinate list
c        map(1:nv) - contains vertex indices of vcl
c
c     updated parameters:
c        map - shifted or 2 swaps applied if necessary so that vertices
c              indexed by map(1), map(2), map(3), map(4) not coplanar
c
c     output parameters:
c        i3,i4 - the indices such that map_in(i3) = map_out(3) and
c              map_in(i4) = map_out(4)
c
c     abnormal return:
c        ierr is set to 302, 303, or 304
c
      integer ierr
      double precision pi,tol
      common /gerror/ ierr
      common /gconst/ pi,tol
      save /gerror/,/gconst/
c
      integer i,k,l,m,m1,m2
      double precision dv2(3),dvk(3),dvl(3)
      double precision cmax,cp1,cp2,cp3,dmax,dotp
c
c     first check that consecutive vertices are not identical.
c
      if (shift) then
	 l = nv - 1
      else
	 l = 1
      endif
      m1 = map(1)
      do 20 i = 1,l
	 m = m1
	 m1 = map(i+1)
	 do 10 k = 1,3
	    cmax = max( abs(vcl(k,m)),abs(vcl(k,m1)) )
	    if (abs(vcl(k,m) - vcl(k,m1)) .gt. tol*cmax .and. cmax
     $         .gt. tol) go to 20
   10    continue
	 ierr = 302
	 return
   20 continue
c
c     find index k = i3 and l = i4.
c
      m1 = map(1)
      m2 = map(2)
      dv2(1) = vcl(1,m2) - vcl(1,m1)
      dv2(2) = vcl(2,m2) - vcl(2,m1)
      dv2(3) = vcl(3,m2) - vcl(3,m1)
      cmax = max( abs(vcl(1,m1)),abs(vcl(2,m1)),abs(vcl(3,m1)),
     $   abs(vcl(1,m2)),abs(vcl(2,m2)),abs(vcl(3,m2)) )
      k = 2
   30 continue
	 k = k + 1
	 if (k .gt. nv) then
	    ierr = 303
	    return
	 endif
	 m = map(k)
	 dvk(1) = vcl(1,m) - vcl(1,m1)
	 dvk(2) = vcl(2,m) - vcl(2,m1)
	 dvk(3) = vcl(3,m) - vcl(3,m1)
	 dmax = max(cmax, abs(vcl(1,m)),abs(vcl(2,m)),abs(vcl(3,m)) )
         cp1 = dv2(2)*dvk(3) - dv2(3)*dvk(2)
         cp2 = dv2(3)*dvk(1) - dv2(1)*dvk(3)
         cp3 = dv2(1)*dvk(2) - dv2(2)*dvk(1)
      if (max(abs(cp1),abs(cp2),abs(cp3)) .le. tol*dmax) go to 30
c
      cmax = dmax
      l = k
   40 continue
	 l = l + 1
	 if (l .gt. nv) then
	    ierr = 304
	    return
	 endif
	 m = map(l)
	 dvl(1) = vcl(1,m) - vcl(1,m1)
	 dvl(2) = vcl(2,m) - vcl(2,m1)
	 dvl(3) = vcl(3,m) - vcl(3,m1)
	 dmax = max(cmax, abs(vcl(1,m)),abs(vcl(2,m)),abs(vcl(3,m)) )
	 dotp = dvl(1)*cp1 + dvl(2)*cp2 + dvl(3)*cp3
      if (abs(dotp) .le. tol*dmax) go to 40
c
c     shift or swap elements of map if necessary.
c
      if (shift) then
         if (k .gt. 3) m1 = map(k)
         if (l .gt. 4) then
	    m2 = map(l)
	    do 50 i = l,k+2,-1
	       map(i) = map(i-1)
   50       continue
	    do 60 i = k+1,5,-1
	       map(i) = map(i-2)
   60       continue
	    map(4) = m2
         endif
         if (k .gt. 3) map(3) = m1
      else
         if (k .gt. 3) then
	    m = map(3)
	    map(3) = map(k)
	    map(k) = m
	 endif
         if (l .gt. 4) then
	    m = map(4)
	    map(4) = map(l)
	    map(l) = m
	 endif
      endif
      i3 = k
      i4 = l
      end
