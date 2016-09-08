      subroutine walkt3(pt,n,p,ntetra,vcl,vm,fc,ht,ifac,ivrt)
      implicit logical (a-z)
      integer ifac,ivrt,n,ntetra,p
      integer fc(7,*),ht(0:p-1),vm(n)
      double precision pt(3),vcl(3,*)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: walk through neighbouring tetrahedra of 3-d (delaunay)
c        triangulation until a tetrahedron is found containing point pt
c        or pt is found to be outside the convex hull. search is
c        guaranteed to terminate for a delaunay triangulation, else a
c        cycle may occur.
c
c     input parameters:
c        pt(1:3) - 3-d point
c        n - upper bound on vertex indices or size of vm
c        p - size of hash table
c        ntetra - number of tetrahedra in triangulation; used to detect
c              cycle
c        vcl(1:3,1:*) - vertex coordinate list
c        vm(1:n) - vertex mapping list (maps from local indices used in
c              fc to indices of vcl)
c        fc(1:7,1:*) - array of face records; see routine dtris3
c        ht(0:p-1) - hash table using direct chaining
c        ifac - index of face of fc to begin search at
c        ivrt - 4 or 5 to indicate tetrahedron to begin search at
c
c     updated parameters:
c        ifac - index of fc indicating tetrahedron or face containing pt
c        ivrt - 4 or 5 to indicate that fc(ivrt,ifac) is 4th vertex of
c              tetrahedron containing pt in its interior; 6 if pt lies
c              in interior of face fc(*,ifac); 0 if pt lies outside
c              convex hull on other side of face fc(*,ifac); 1, 2, or 3
c              if pt lies on interior of edge of face from vertices
c              fc(ivrt,ifac) to fc(ivrt mod 3 + 1,ifac); -1, -2, or -3
c              if pt is (nearly) vertex fc(-ivrt,ifac)
c
c     abnormal return:
c        ierr is set to 300, 301, or 307
c
c     routines called:
c        baryth,htsrc
c
      integer ierr
      common /gerror/ ierr
      save /gerror/
c
      integer a,aa,b,bb,c,cnt,d,htsrc,j,k,va,vb,vc,vd
      double precision t(4)
      logical degen,zero(4)
c
      cnt = 0
   10 continue
	 if (fc(ivrt,ifac) .le. 0) then
	    ivrt = 0
	    return
	 endif
	 cnt = cnt + 1
	 if (cnt .gt. ntetra) then
	    ierr = 307
	    return
	 endif
	 a = fc(1,ifac)
	 b = fc(2,ifac)
	 c = fc(3,ifac)
	 d = fc(ivrt,ifac)
	 va = vm(a)
	 vb = vm(b)
	 vc = vm(c)
	 vd = vm(d)
	 call baryth(vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,vd),pt,t,degen)
	 if (degen) then
	    ierr = 301
	    return
	 endif
	 if (t(1) .gt. 0.0d0 .and. t(2) .gt. 0.0d0 .and. t(3) .gt. 0.0d0
     $      .and. t(4) .gt. 0.0d0) then
	    return
	 else if (t(4) .lt. 0.0d0) then
	    ivrt = 9 - ivrt
	 else if (t(1) .lt. 0.0d0) then
	    ifac = htsrc(b,c,d,n,p,fc,ht)
	    if (ifac .le. 0) go to 30
	    if (fc(4,ifac) .eq. a) then
	       ivrt = 5
	    else
	       ivrt = 4
	    endif
	 else if (t(2) .lt. 0.0d0) then
	    ifac = htsrc(a,c,d,n,p,fc,ht)
	    if (ifac .le. 0) go to 30
	    if (fc(4,ifac) .eq. b) then
	       ivrt = 5
	    else
	       ivrt = 4
	    endif
	 else if (t(3) .lt. 0.0d0) then
	    ifac = htsrc(a,b,d,n,p,fc,ht)
	    if (ifac .le. 0) go to 30
	    if (fc(4,ifac) .eq. c) then
	       ivrt = 5
	    else
	       ivrt = 4
	    endif
	 else
c           all t(j) .ge. 0.0d0 and at least one t(j) .eq. 0.0d0
	    k = 0
	    do 20 j = 1,4
	       zero(j) = (t(j) .eq. 0.0d0)
	       if (zero(j)) k = k + 1
   20       continue
	    if (k .eq. 1) then
	       ivrt = 6
	       if (zero(1)) then
	          ifac = htsrc(b,c,d,n,p,fc,ht)
	          if (ifac .le. 0) go to 30
	       else if (zero(2)) then
	          ifac = htsrc(a,c,d,n,p,fc,ht)
	          if (ifac .le. 0) go to 30
	       else if (zero(3)) then
	          ifac = htsrc(a,b,d,n,p,fc,ht)
	          if (ifac .le. 0) go to 30
	       endif
	    else if (k .eq. 2) then
	       if (zero(4)) then
		  if (zero(3)) then
		     ivrt = 1
		  else if (zero(1)) then
		     ivrt = 2
		  else
		     ivrt = 3
		  endif
	       else
		  if (zero(3)) then
		     ifac = htsrc(a,b,d,n,p,fc,ht)
		     if (ifac .le. 0) go to 30
		     if (zero(2)) then
		        aa = a
		     else
		        aa = b
		     endif
	          else
		     ifac = htsrc(a,c,d,n,p,fc,ht)
		     if (ifac .le. 0) go to 30
		     aa = c
		  endif
		  bb = d
		  if (aa .gt. bb) then
		     j = aa
		     aa = bb
		     bb = j
		  endif
	          if (fc(1,ifac) .eq. aa) then
	             if (fc(2,ifac) .eq. bb) then
		        ivrt = 1
		     else
		        ivrt = 3
		     endif
	          else
		     ivrt = 2
	          endif
	       endif
	    else
c              k .eq. 3
	       if (.not. zero(1)) then
		  ivrt = -1
	       else if (.not. zero(2)) then
		  ivrt = -2
	       else if (.not. zero(3)) then
		  ivrt = -3
	       else
		  ifac = htsrc(a,b,d,n,p,fc,ht)
		  if (ifac .le. 0) go to 30
		  if (fc(1,ifac) .eq. d) then
		     ivrt = -1
		  else if (fc(2,ifac) .eq. d) then
		     ivrt = -2
		  else
		     ivrt = -3
		  endif
	       endif
	    endif
	    return
	 endif
      go to 10
   30 continue
      ierr = 300
      end
