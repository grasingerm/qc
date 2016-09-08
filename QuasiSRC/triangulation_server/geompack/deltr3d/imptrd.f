      subroutine imptrd(bndcon,npt,sizht,maxfc,vcl,vm,nfc,ntetra,hdavfc,
     $   bf,fc,ht)
      implicit logical (a-z)
      logical bndcon
      integer hdavfc,maxfc,nfc,npt,ntetra,sizht
      integer bf(3,*),fc(7,maxfc),ht(0:sizht-1),vm(npt)
      double precision vcl(3,*)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: further improve given 3-d triangulation towards delaunay
c        one by using combination local transformations (not yet
c        guaranteed to produce delaunay triangulation).
c
c     input parameters:
c        bndcon - .true. iff boundary faces are constrained (i.e. not
c              swapped by local transformations)
c	 npt - number of 3-d vertices (points)
c	 sizht - size of hash table ht; a good choice is a prime number
c              which is about 1/8 * nface (or 3/2 * npt for random
c              points from the uniform distribution)
c        maxfc - maximum size available for fc array
c        vcl(1:3,1:*) - vertex coordinate list
c        vm(1:npt) - indices of vertices of vcl being triangulated
c        nfc - number of positions used in fc array
c        ntetra - number of tetrahedra in triangulation
c        hdavfc - head pointer to available fc records
c        bf(1:3,1:*) -  array of boundary face records; see dtris3
c        fc(1:7,1:maxfc) - array of face records; see routine dtris3
c        ht(0:sizht-1) - hash table using direct chaining
c        [note: bf, fc, ht should be as output by dtris3 or dtriw3,
c         except it is assumed fc(7,1:2) don't contain hdavbf, hdavfc.]
c
c     updated parameters:
c        nfc,ntetra,bf,fc,ht
c
c     abnormal return:
c        ierr is set to 11, 300, 301, 305, 308, or 309
c
c     routines called:
c        ccradi,ccsph,htsrc,ifacty,swapes,swaptf
c
      integer ierr,iprt,msglvl
      double precision pi,tol
      common /gerror/ ierr
      common /gconst/ pi,tol
      common /gprint/ iprt,msglvl
      save /gerror/,/gconst/,/gprint/
c
      integer a,aa,back,b,bb,c,cc,c2,d,e,f,front,g,h,htsrc,ibeg,in
      integer ind,ind0,ind1,ind2,ind3,j,top,va,vb,vc,vd,ve,vf,vg,vh
      double precision ccradi,maxaft,maxbef,mu1,mu2,mu3,mu4,mu5,mu6,mu7
      double precision mu8,nu1,nu2,nu3,nu4,nu5,nu6,nu7,tolp1,centre(3)
      logical first,t1,t2
      character type*3,typ2*3
c
      if (msglvl .eq. 4) write (iprt,600)
      front = 0
      back = 0
      ibeg = 1
      ind = 1
      tolp1 = 1.0d0 + tol
   10 continue
         if (fc(1,ind) .le. 0 .or. fc(5,ind) .le. 0) go to 70
	 call ifacty(ind,npt,sizht,vcl,vm,fc,ht,type,a,b,c)
         if (ierr .ne. 0) return
	 if (type .ne. 'n32' .and. type .ne. 'n44') go to 70
	 d = fc(4,ind)
	 e = fc(5,ind)
	 va = vm(a)
	 vb = vm(b)
	 vc = vm(c)
	 vd = vm(d)
	 ve = vm(e)
	 call ccsph(.true.,vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,vd),
     $      vcl(1,ve),centre,mu1,in)
	 if (in .eq. 2) then
	    ierr = 301
	    return
	 endif
	 if (in .le. 0) go to 70
         if (msglvl .eq. 4) write (iprt,610) type,ind,a,b,c,d,e
	 ind1 = htsrc(a,b,d,npt,sizht,fc,ht)
	 if (ind1 .le. 0) go to 80
	 if (fc(4,ind1) .eq. c) then
	    f = fc(5,ind1)
	 else
	    f = fc(4,ind1)
	 endif
	 ind2 = htsrc(a,b,f,npt,sizht,fc,ht)
	 if (ind2 .le. 0) go to 80
	 mu1 = 0.0625d0/mu1
	 if (type .eq. 'n44') go to 20
c
c        type .eq. 'n32'
 	 if (fc(4,ind2) .ne. e .and. fc(5,ind2) .ne. e) go to 70
	 call ifacty(ind2,npt,sizht,vcl,vm,fc,ht,type,aa,bb,cc)
         if (msglvl .eq. 4) write (iprt,620) aa,bb,cc,d,e,type
	 if (type .ne. 't23' .and. type .ne. 't32' .and. type .ne. 't44'
     $      .and. type .ne. 'n32') go to 70
	 vf = vm(f)
	 mu2 = ccradi(vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,ve))
	 mu3 = ccradi(vcl(1,va),vcl(1,vb),vcl(1,vf),vcl(1,vd))
	 mu4 = ccradi(vcl(1,va),vcl(1,vb),vcl(1,vf),vcl(1,ve))
	 nu1 = ccradi(vcl(1,va),vcl(1,vc),vcl(1,vd),vcl(1,ve))
	 nu2 = ccradi(vcl(1,vb),vcl(1,vc),vcl(1,vd),vcl(1,ve))
	 maxbef = max(mu1,mu2,mu3,mu4)
	 if (type .eq. 't23') then
	    nu3 = ccradi(vcl(1,va),vcl(1,vf),vcl(1,vd),vcl(1,ve))
	    nu4 = ccradi(vcl(1,vb),vcl(1,vf),vcl(1,vd),vcl(1,ve))
	    if (max(nu1,nu2,nu3,nu4) .le. maxbef*tolp1) go to 70
	    top = ind2
	    fc(7,ind2) = ind
	    fc(7,ind) = 0
	    go to 60
	 else if (type .eq. 't32') then
	    if (cc .eq. a) then
	       j = va
	       va = vb
	       vb = j
	    endif
	    mu5 = ccradi(vcl(1,vf),vcl(1,va),vcl(1,vd),vcl(1,ve))
	    nu3 = ccradi(vcl(1,vf),vcl(1,vb),vcl(1,vd),vcl(1,ve))
	    if (max(nu1,nu2,nu3) .le. max(maxbef,mu5)*tolp1) go to 70
	    top = ind2
	    fc(7,ind2) = ind
	    fc(7,ind) = 0
	    go to 60
	 else if (type .eq. 't44') then
	    if (cc .eq. a) then
	       j = a
	       a = b
	       b = j
	       j = va
	       va = vb
	       vb = j
	    endif
	    ind1 = htsrc(a,d,f,npt,sizht,fc,ht)
	    if (ind1 .le. 0) go to 80
	    if (fc(4,ind1) .eq. b) then
	       g = fc(5,ind1)
	    else
	       g = fc(4,ind1)
	    endif
	    vg = vm(g)
	    mu5 = ccradi(vcl(1,va),vcl(1,vf),vcl(1,vg),vcl(1,vd))
	    mu6 = ccradi(vcl(1,va),vcl(1,vf),vcl(1,vg),vcl(1,ve))
	    nu3 = ccradi(vcl(1,vf),vcl(1,vb),vcl(1,vd),vcl(1,ve))
	    nu4 = ccradi(vcl(1,va),vcl(1,vg),vcl(1,vd),vcl(1,ve))
	    nu5 = ccradi(vcl(1,vf),vcl(1,vg),vcl(1,vd),vcl(1,ve))
	    if (max(nu1,nu2,nu3,nu4,nu5) .le. max(maxbef,mu5,mu6)*tolp1)
     $         go to 70
	    top = ind2
	    fc(7,ind2) = ind
	    fc(7,ind) = 0
	    go to 60
	 else
	    if (cc .eq. a) then
	       j = a
	       a = b
	       b = j
	       j = va
	       va = vb
	       vb = j
	    endif
	    ind1 = htsrc(a,f,d,npt,sizht,fc,ht)
	    if (ind1 .le. 0) go to 80
	    if (fc(4,ind1) .eq. b) then
	       g = fc(5,ind1)
	    else
	       g = fc(4,ind1)
	    endif
	    ind1 = htsrc(a,f,g,npt,sizht,fc,ht)
	    if (ind1 .le. 0) go to 80
	    if (fc(4,ind1) .ne. e .and. fc(5,ind1) .ne. e) go to 70
	    call ifacty(ind1,npt,sizht,vcl,vm,fc,ht,type,aa,bb,cc)
	    if (type .ne. 't23' .and. type .ne. 't32' .and. type .ne.
     $         'n32') go to 70
	    vg = vm(g)
	    mu5 = ccradi(vcl(1,va),vcl(1,vf),vcl(1,vg),vcl(1,vd))
	    mu6 = ccradi(vcl(1,va),vcl(1,vf),vcl(1,vg),vcl(1,ve))
	    nu3 = ccradi(vcl(1,vb),vcl(1,vd),vcl(1,ve),vcl(1,vf))
	    if (type .eq. 't23') then
	       maxbef = max(maxbef,mu5,mu6)
	       nu4 = ccradi(vcl(1,va),vcl(1,vg),vcl(1,vd),vcl(1,ve))
	       nu5 = ccradi(vcl(1,vf),vcl(1,vg),vcl(1,vd),vcl(1,ve))
	       if (max(nu1,nu2,nu3,nu4,nu5) .le. maxbef*tolp1) go to 70
	       top = ind1
	    else if (type .eq. 't32') then
	       if (cc .eq. a) then
	          j = va
	          va = vf
	          vf = j
	       endif
	       mu7 = ccradi(vcl(1,va),vcl(1,vg),vcl(1,vd),vcl(1,ve))
	       maxbef = max(maxbef,mu5,mu6,mu7)
	       nu4 = ccradi(vcl(1,vf),vcl(1,vg),vcl(1,vd),vcl(1,ve))
	       if (max(nu1,nu2,nu3,nu4) .le. maxbef*tolp1) go to 70
	       top = ind1
	    else
	       if (g .eq. c) go to 70
	       if (cc .eq. a) then
	          j = a
	          a = f
	          f = j
	          j = va
	          va = vf
	          vf = j
	       endif
	       ind0 = htsrc(a,g,d,npt,sizht,fc,ht)
	       if (ind0 .le. 0) go to 80
	       if (fc(4,ind0) .eq. f) then
	          h = fc(5,ind0)
	       else
	          h = fc(4,ind0)
	       endif
	       ind0 = htsrc(a,g,h,npt,sizht,fc,ht)
	       if (ind0 .le. 0) go to 80
	       if (fc(4,ind0) .ne. e .and. fc(5,ind0) .ne. e) go to 70
	       call ifacty(ind0,npt,sizht,vcl,vm,fc,ht,type,aa,bb,cc)
	       if (type .ne. 't23') go to 70
	       vh = vm(h)
	       mu7 = ccradi(vcl(1,va),vcl(1,vg),vcl(1,vh),vcl(1,vd))
	       mu8 = ccradi(vcl(1,va),vcl(1,vg),vcl(1,vh),vcl(1,ve))
	       maxbef = max(maxbef,mu5,mu6,mu7,mu8)
	       nu4 = ccradi(vcl(1,vf),vcl(1,vg),vcl(1,vd),vcl(1,ve))
	       nu5 = ccradi(vcl(1,va),vcl(1,vh),vcl(1,vd),vcl(1,ve))
	       nu6 = ccradi(vcl(1,vg),vcl(1,vh),vcl(1,vd),vcl(1,ve))
	       if (max(nu1,nu2,nu3,nu4,nu5,nu6) .le. maxbef*tolp1)
     $            go to 70
	       top = ind0
	       fc(7,ind0) = ind1
	    endif
	    fc(7,ind1) = ind2
	    fc(7,ind2) = ind
	    fc(7,ind) = 0
	    go to 60
	 endif
c
c        type .eq. 'n44'
   20    continue
	 ind3 = ind2
	 if (fc(4,ind3) .eq. d) then
	    g = fc(5,ind3)
	 else
	    g = fc(4,ind3)
	 endif
	 ind2 = htsrc(a,b,g,npt,sizht,fc,ht)
	 if (ind2 .le. 0) go to 80
	 if (fc(4,ind2) .ne. e .and. fc(5,ind2) .ne. e) go to 70
	 call ifacty(ind2,npt,sizht,vcl,vm,fc,ht,type,aa,bb,cc)
         if (msglvl .eq. 4) write (iprt,620) aa,bb,cc,e,f,type
	 t1 = (type .eq. 't23' .or. type .eq. 't32' .or. type .eq.'t44')
	 call ifacty(ind3,npt,sizht,vcl,vm,fc,ht,typ2,aa,bb,c2)
         if (msglvl .eq. 4) write (iprt,620) aa,bb,c2,d,g,typ2
	 t2 = (typ2 .eq. 't23' .or. typ2 .eq. 't32' .or. typ2 .eq.'t44')
	 if (.not. t1 .and. .not. t2) go to 70
	 first = .true.
   30    continue
	 if (t1) go to 40
	    ind2 = ind3
	    type = typ2
	    cc = c2
	    j = d
	    d = e
	    e = j
	    j = vd
	    vd = ve
	    ve = j
	    j = f
	    f = g
	    g = j
   40    continue
	 vf = vm(f)
	 vg = vm(g)
         if (first) then
	    mu2 = ccradi(vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,ve))
	    mu3 = ccradi(vcl(1,va),vcl(1,vb),vcl(1,vf),vcl(1,vd))
	    mu4 = ccradi(vcl(1,va),vcl(1,vb),vcl(1,vg),vcl(1,ve))
	    mu5 = ccradi(vcl(1,va),vcl(1,vb),vcl(1,vg),vcl(1,vf))
	    nu1 = ccradi(vcl(1,va),vcl(1,vc),vcl(1,vd),vcl(1,ve))
	    nu2 = ccradi(vcl(1,vb),vcl(1,vc),vcl(1,vd),vcl(1,ve))
	    nu3 = ccradi(vcl(1,va),vcl(1,vf),vcl(1,vd),vcl(1,ve))
	    nu4 = ccradi(vcl(1,vb),vcl(1,vf),vcl(1,vd),vcl(1,ve))
         else
	    nu3 = ccradi(vcl(1,va),vcl(1,vf),vcl(1,vd),vcl(1,ve))
	    nu4 = ccradi(vcl(1,vb),vcl(1,vf),vcl(1,vd),vcl(1,ve))
	 endif
	 if (first) then
	    maxbef = max(mu1,mu2,mu3,mu4,mu5)
	    maxaft = max(nu1,nu2,nu3,nu4)
	 else
	    maxaft = max(nu1,nu2,nu3,nu4)
	 endif
	 if (type .eq. 't23') then
	    nu5 = ccradi(vcl(1,va),vcl(1,vg),vcl(1,ve),vcl(1,vf))
	    nu6 = ccradi(vcl(1,vb),vcl(1,vg),vcl(1,ve),vcl(1,vf))
	    if (max(maxaft,nu5,nu6) .le. maxbef*tolp1) go to 50
	    top = ind2
	    fc(7,ind2) = ind
	    fc(7,ind) = 0
	    go to 60
	 else if (type .eq. 't32') then
	    if (cc .eq. a) then
	       j = a
	       a = b
	       b = j
	       j = va
	       va = vb
	       vb = j
	    endif
	    mu6 = ccradi(vcl(1,va),vcl(1,vg),vcl(1,ve),vcl(1,vf))
	    nu5 = ccradi(vcl(1,vb),vcl(1,vg),vcl(1,ve),vcl(1,vf))
	    if (max(maxaft,nu5) .le. max(maxbef,mu6)*tolp1) go to 50
	    top = ind2
	    fc(7,ind2) = ind
	    fc(7,ind) = 0
	    go to 60
	 else
	    if (cc .eq. a) then
	       j = a
	       a = b
	       b = j
	       j = va
	       va = vb
	       vb = j
	    endif
	    ind1 = htsrc(a,e,g,npt,sizht,fc,ht)
	    if (ind1 .le. 0) go to 80
	    if (fc(4,ind1) .eq. b) then
	       h = fc(5,ind1)
	    else
	       h = fc(4,ind1)
	    endif
	    vh = vm(h)
	    mu6 = ccradi(vcl(1,va),vcl(1,vg),vcl(1,vh),vcl(1,ve))
	    mu7 = ccradi(vcl(1,va),vcl(1,vg),vcl(1,vh),vcl(1,vf))
	    nu5 = ccradi(vcl(1,vg),vcl(1,vb),vcl(1,ve),vcl(1,vf))
	    nu6 = ccradi(vcl(1,va),vcl(1,vh),vcl(1,ve),vcl(1,vf))
	    nu7 = ccradi(vcl(1,vg),vcl(1,vh),vcl(1,ve),vcl(1,vf))
	    if (max(maxaft,nu5,nu6,nu7) .le. max(maxbef,mu6,mu7)*tolp1)
     $         go to 50
	    top = ind2
	    fc(7,ind2) = ind
	    fc(7,ind) = 0
	    go to 60
	 endif
   50    continue
	 if (t1 .and. t2) then
	    t1 = .false.
	    first = .false.
	    go to 30
	 else
	    go to 70
	 endif
c
   60    continue
	 if (msglvl .eq. 4) write (iprt,630) 'combination swaps made'
	 call swaptf(top,npt,sizht,nfc,maxfc,vcl,vm,fc,ht,ntetra,hdavfc,
     $      front,back)
         if (ierr .ne. 0) return
	 call swapes(bndcon,0,npt,sizht,nfc,maxfc,vcl,vm,bf,fc,ht,
     $      ntetra,hdavfc,front,back,j)
         if (ierr .ne. 0) return
         ind = ind + 1
         if (ind .gt. nfc) ind = 1
	 ibeg = ind
	 go to 10
   70    continue
         ind = ind + 1
         if (ind .gt. nfc) ind = 1
      if (ind .ne. ibeg) go to 10
      return
c
   80 continue
      ierr = 300
c
  600 format (/1x,'imptrd')
  610 format (1x,'type ',a3,i7,' : ',5i7)
  620 format (4x,'face',3i7,' | ',2i7,' has type ',a3)
  630 format (4x,a)
      end
