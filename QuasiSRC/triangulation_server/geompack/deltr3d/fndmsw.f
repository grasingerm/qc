      subroutine fndmsw(crit,npt,sizht,vcl,vm,fc,ht,a,b,d,e,f,minbef,
     $   top,top2,impr)
      implicit logical (a-z)
      integer a,b,crit,d,e,f,npt,sizht,top,top2
      integer fc(7,*),ht(0:sizht-1),vm(npt)
      double precision minbef,vcl(3,*)
      logical impr
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: find a sequence of >= 3 local transformations to improve
c        3-d triangulation.
c
c     input parameters:
c        crit - criterion code; 1 for (local max-min) solid angle
c              criterion, 2 for radius ratio criterion, 3 for mean ratio
c              criterion, 0 (or anything else) for no swaps
c	 npt - number of 3-d vertices (points)
c	 sizht - size of hash table ht
c        vcl(1:3,1:*) - vertex coordinate list
c        vm(1:npt) - indices of vertices of vcl being triangulated
c        fc(1:7,1:*) - array of face records; see routine dtris3
c        ht(0:sizht-1) - hash table using direct chaining
c        a,b,d,e,f - vertices (local indices) in configuration with t23
c              face afb|de swapped out to produce only 1 tetra afde with
c              measure <= minbef; try to apply a t32 swap to remove afde
c              where af is the desired edge to be removed
c        minbef - (min tetra measure of an existing tetra of swap) + tol
c        top - pointer to list of 2 or 3 faces to be possibly swapped
c        top2 - pointer to stack of other faces of t32 or t44 swaps
c
c     updated parameters:
c        a,b,d,e,f - may be updated
c        fc(7,*) - some faces may be added to list
c        top2 - gets set to 0 and stack emptied
c
c     output parameters:
c        top - pointer to list of faces to be swapped if impr is
c              .true., else top = 0 and all faces removed from list
c        impr - .true. iff improvement is possible
c
c     abnormal return:
c        ierr is set to 300, 301, or 309
c
c     routines called:
c        baryth,htsrc,ifacty,tetmu
c
      integer ierr
      common /gerror/ ierr
      save /gerror/
c
      integer aa,bb,cc,c2,g,h,htsrc,i,ind,ind1,ind2,indx,indy,j,kf,kneg
      integer kzero,maxtf,ptr,va,vd,ve,vf,vg,vh,vi
      double precision nu1,nu2,nu3,nu4,nu5,tetmu,alpha(4)
      logical degen
      character type*3,typ2*3
      parameter (maxtf = 13)
c
      kf = 2
      impr = .false.
      ptr = top
   10 continue
      if (kf .ge. maxtf) go to 40
      indx = htsrc(a,f,d,npt,sizht,fc,ht)
      if (indx .le. 0) go to 60
      if (fc(4,indx) .eq. b) then
	 g = fc(5,indx)
      else
	 g = fc(4,indx)
      endif
      ind = htsrc(a,f,e,npt,sizht,fc,ht)
      if (ind .le. 0) go to 60
      if (fc(4,ind) .eq. b) then
	 h = fc(5,ind)
      else
	 h = fc(4,ind)
      endif
      if (g .le. 0 .or. h .le. 0) go to 40
      if (g .ne. h) then
	 ind1 = htsrc(a,f,h,npt,sizht,fc,ht)
         if (ind1 .le. 0) go to 60
	 if (fc(4,ind1) .ne. g .and. fc(5,ind1) .ne. g) go to 40
	 call ifacty(ind1,npt,sizht,vcl,vm,fc,ht,type,aa,bb,cc)
	 if (type .ne. 't23') then
	    ind2 = ind1
	    c2 = cc
	    typ2 = type
	    ind1 = htsrc(a,f,g,npt,sizht,fc,ht)
            if (ind1 .le. 0) go to 60
	    call ifacty(ind1,npt,sizht,vcl,vm,fc,ht,type,aa,bb,cc)
	    if (type .eq. 't23' .or. type .eq. 'n32') then
	       j = d
	       d = e
	       e = j
	       j = g
	       g = h
	       h = j
	    else if (typ2 .eq. 'n32') then
	       type = typ2
	       ind1 = ind2
	       cc = c2
	    else
	       go to 40
	    endif
	 endif
      endif
      va = vm(a)
      vd = vm(d)
      ve = vm(e)
      vf = vm(f)
      vg = vm(g)
      call baryth(vcl(1,va),vcl(1,vf),vcl(1,vd),vcl(1,ve),vcl(1,vg),
     $   alpha,degen)
      if (degen) then
	 ierr = 301
	 return
      else if (alpha(4) .ge. 0.0d0) then
         ierr = 309
         write(*,*) 'error 309 in fndmsw'
         write(*,*) vcl(1,va), vcl(2, va), vcl(3, va)
         write(*,*) vcl(1,vf), vcl(2, vf), vcl(3, vf)
         write(*,*) vcl(1,vd), vcl(2, vd), vcl(3, vd)
         write(*,*) vcl(1,ve), vcl(2, ve), vcl(3, ve)
         write(*,*) vcl(1,vg), vcl(2, vg), vcl(3, vg)
	 return
      endif
      kneg = 1
      kzero = 0
      do 20 j = 1,3
	 if (alpha(j) .lt. 0.0d0) then
	    kneg = kneg + 1
	 else if (alpha(j) .eq. 0.0d0) then
	    kzero = kzero + 1
	 endif
   20 continue
      if (kneg .ne. 2 .or. kzero .ne. 0 .or. alpha(3) .ge. 0.0d0)
     $   go to 40
      if (fc(7,ind) .ne. -1 .or. fc(7,indx) .ne. -1) go to 40
      indy = htsrc(a,f,g,npt,sizht,fc,ht)
      if (indy .le. 0) go to 60
      if (fc(7,indy) .ne. -1) go to 40
      nu1 = tetmu(crit,vcl(1,va),vcl(1,vd),vcl(1,ve),vcl(1,vg),alpha)
      nu2 = tetmu(crit,vcl(1,vf),vcl(1,vd),vcl(1,ve),vcl(1,vg),alpha)
      if (min(nu1,nu2) .le. minbef) go to 40
      fc(7,ind) = fc(7,ptr)
      fc(7,ptr) = ind
      kf = kf + 1
      if (g .eq. h) then
c        last face added to middle of list has type t32.
	 impr = .true.
         go to 50
      endif
      fc(7,indx) = indy
      fc(7,indy) = top2
      top2 = indx
      if (type .eq. 'n32') go to 30
      if (fc(7,ind1) .ne. -1) go to 40
      vh = vm(h)
      nu3 = tetmu(crit,vcl(1,va),vcl(1,vh),vcl(1,ve),vcl(1,vg),alpha)
      nu4 = tetmu(crit,vcl(1,vf),vcl(1,vh),vcl(1,ve),vcl(1,vg),alpha)
      if (max(nu3,nu4) .le. minbef) go to 40
      fc(7,ind1) = fc(7,ptr)
      fc(7,ptr) = ind1
      ptr = ind1
      kf = kf + 1
      if (min(nu3,nu4) .gt. minbef) then
c        last face added to middle of list has type t23.
	 impr = .true.
	 go to 50
      endif
      if (nu4 .lt. nu3) then
	 j = a
	 a = f
	 f = j
      endif
      b = f
      d = g
      f = h
      go to 10
c
   30 continue
      if (cc .eq. a) then
	 j = a
	 a = f
	 f = j
	 j = va
	 va = vf
	 vf = j
      endif
      indx = htsrc(a,h,e,npt,sizht,fc,ht)
      if (indx .le. 0) go to 60
      if (fc(7,indx) .ne. -1) go to 40
      if (fc(4,indx) .eq. f) then
	 i = fc(5,indx)
      else
	 i = fc(4,indx)
      endif
      ind2 = htsrc(a,h,i,npt,sizht,fc,ht)
      if (ind2 .le. 0) go to 60
      if (fc(4,ind2) .ne. g .and. fc(5,ind2) .ne. g) go to 40
      call ifacty(ind2,npt,sizht,vcl,vm,fc,ht,type,aa,bb,cc)
      if (type .ne. 't23') go to 40
      vh = vm(h)
      vi = vm(i)
      nu3 = tetmu(crit,vcl(1,ve),vcl(1,vf),vcl(1,vg),vcl(1,vh),alpha)
      if (nu3 .le. minbef) go to 40
      nu4 = tetmu(crit,vcl(1,va),vcl(1,vi),vcl(1,ve),vcl(1,vg),alpha)
      nu5 = tetmu(crit,vcl(1,vh),vcl(1,vi),vcl(1,ve),vcl(1,vg),alpha)
      if (max(nu4,nu5) .le. minbef) go to 40
      if (fc(7,ind1) .ne. -1 .or. fc(7,ind2) .ne. -1) go to 40
      indy = htsrc(a,h,g,npt,sizht,fc,ht)
      if (indy .le. 0) go to 60
      if (fc(7,indy) .ne. -1) go to 40
      fc(7,ind1) = fc(7,ptr)
      fc(7,ptr) = ind2
      fc(7,ind2) = ind1
      ptr = ind2
      kf = kf + 2
      if (min(nu4,nu5) .gt. minbef) then
c        last 2 faces added to middle of list have type t23, t32.
	 impr = .true.
	 go to 50
      endif
      fc(7,indx) = indy
      fc(7,indy) = top2
      top2 = indx
      if (nu5 .lt. nu4) then
	 j = a
	 a = h
	 h = j
      endif
      b = h
      d = g
      f = i
      go to 10
c
   40 continue
         ptr = top
         top = fc(7,ptr)
         fc(7,ptr) = -1
      if (top .ne. 0) go to 40
   50 continue
         ptr = top2
         top2 = fc(7,ptr)
         fc(7,ptr) = -1
      if (top2 .ne. 0) go to 50
      return
c
   60 continue
      ierr = 300
      end
