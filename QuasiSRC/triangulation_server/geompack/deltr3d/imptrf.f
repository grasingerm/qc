      subroutine imptrf(bndcon,crit,npt,sizht,maxfc,vcl,vm,nfc,ntetra,
     $   hdavfc,bf,fc,ht)
      implicit logical (a-z)
      logical bndcon
      integer crit,hdavfc,maxfc,nfc,npt,ntetra,sizht
      integer bf(3,*),fc(7,maxfc),ht(0:sizht-1),vm(npt)
      double precision vcl(3,*)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: further improve a given 3-d triangulation by applying
c        local transformations based on a local criterion. combination
c        swaps are used to remove poorly shaped tetrahedra.
c
c     input parameters:
c        bndcon - .true. iff boundary faces are constrained (i.e. not
c              swapped by local transformations)
c        crit - criterion code; 1 for (local max-min) solid angle
c              criterion, 2 for radius ratio criterion, 3 for mean ratio
c              criterion, 0 (or anything else) for no swaps
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
c        ierr is set to 11, 300, 301, 308, or 309
c
c     routines called:
c        fndmsw,htsrc,ifacty,swapmu,swaptf,tetmu
c
      integer ierr,iprt,msglvl
      double precision pi,tol
      common /gerror/ ierr
      common /gconst/ pi,tol
      common /gprint/ iprt,msglvl
      save /gerror/,/gconst/,/gprint/
c
      integer a,aa,back,b,bb,c,cc,c2,d,dd,e,ee,f,ff,front,g,h,htsrc
      integer ibeg,ind,ind1,ind2,ind3,indx,indy,indz,j,top,top2,va,vb
      integer vc,vd,ve,vf,vg,vh
      double precision minaft,minbef,mu1,mu2,mu3,mu4,mu5,mu6,mu7
      double precision nu1,nu2,nu3,nu4,nu5,nu6,nu7,tetmu,s(4)
      logical first,impr,t1,t2
      character type*3,typ2*3
c
      if (msglvl .eq. 4) write (iprt,600) crit
      front = 0
      back = 0
      ibeg = 1
      ind = 1
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
         if (msglvl .eq. 4) write (iprt,610) type,ind,a,b,c,d,e
	 indx = htsrc(a,b,d,npt,sizht,fc,ht)
	 if (indx .le. 0) go to 80
	 if (fc(4,indx) .eq. c) then
	    f = fc(5,indx)
	 else
	    f = fc(4,indx)
	 endif
	 ind2 = htsrc(a,b,f,npt,sizht,fc,ht)
	 if (ind2 .le. 0) go to 80
	 if (type .eq. 'n44') go to 20
c
c        type .eq. 'n32'
	 if (fc(4,ind2) .ne. e .and. fc(5,ind2) .ne. e) go to 70
	 call ifacty(ind2,npt,sizht,vcl,vm,fc,ht,type,aa,bb,cc)
         if (msglvl .eq. 4) write (iprt,620) aa,bb,cc,d,e,type
	 if (type .ne. 't23' .and. type .ne. 't32' .and. type .ne.'t44'
     $      .and. type .ne. 'n32') go to 70
	 vf = vm(f)
	 mu1 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,vd),s)
	 mu2 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,ve),s)
	 mu3 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vf),vcl(1,vd),s)
	 mu4 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vf),vcl(1,ve),s)
	 nu1 = tetmu(crit,vcl(1,va),vcl(1,vc),vcl(1,vd),vcl(1,ve),s)
	 nu2 = tetmu(crit,vcl(1,vb),vcl(1,vc),vcl(1,vd),vcl(1,ve),s)
	 minbef = min(mu1,mu2,mu3,mu4)
	 if (type .eq. 't23') then
	    minbef = minbef + tol
	    if (min(nu1,nu2) .le. minbef) go to 70
	    nu3 = tetmu(crit,vcl(1,va),vcl(1,vf),vcl(1,vd),vcl(1,ve),s)
	    nu4 = tetmu(crit,vcl(1,vb),vcl(1,vf),vcl(1,vd),vcl(1,ve),s)
	    if (max(nu3,nu4) .le. minbef) go to 70
	    top = ind2
	    fc(7,ind2) = ind
	    fc(7,ind) = 0
	    if (min(nu3,nu4) .le. minbef) then
	       indy = htsrc(a,b,e,npt,sizht,fc,ht)
	       if (indy .le. 0) go to 80
	       top2 = indx
	       fc(7,indx) = indy
	       fc(7,indy) = 0
	       if (nu4 .lt. nu3) then
		  j = a
		  a = b
		  b = j
	       endif
 	       call fndmsw(crit,npt,sizht,vcl,vm,fc,ht,a,b,d,e,f,minbef,
     $            top,top2,impr)
	       if (ierr .ne. 0) return
	       if (.not. impr) go to 70
	    endif
	    go to 60
	 else if (type .eq. 't32') then
	    if (cc .eq. a) then
	       j = va
	       va = vb
	       vb = j
	    endif
	    mu5 = tetmu(crit,vcl(1,vf),vcl(1,va),vcl(1,vd),vcl(1,ve),s)
	    nu3 = tetmu(crit,vcl(1,vf),vcl(1,vb),vcl(1,vd),vcl(1,ve),s)
	    if (min(nu1,nu2,nu3) .le. min(minbef,mu5) + tol) go to 70
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
	    mu5 = tetmu(crit,vcl(1,va),vcl(1,vf),vcl(1,vg),vcl(1,vd),s)
	    mu6 = tetmu(crit,vcl(1,va),vcl(1,vf),vcl(1,vg),vcl(1,ve),s)
	    nu3 = tetmu(crit,vcl(1,vf),vcl(1,vb),vcl(1,vd),vcl(1,ve),s)
	    nu4 = tetmu(crit,vcl(1,va),vcl(1,vg),vcl(1,vd),vcl(1,ve),s)
	    nu5 = tetmu(crit,vcl(1,vf),vcl(1,vg),vcl(1,vd),vcl(1,ve),s)
	    if (min(nu1,nu2,nu3,nu4,nu5) .le. min(minbef,mu5,mu6) + tol)
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
	    ind3 = htsrc(a,f,d,npt,sizht,fc,ht)
	    if (ind3 .le. 0) go to 80
	    if (fc(4,ind3) .eq. b) then
	       g = fc(5,ind3)
	    else
	       g = fc(4,ind3)
	    endif
	    ind1 = htsrc(a,f,g,npt,sizht,fc,ht)
	    if (ind1 .le. 0) go to 80
	    if (fc(4,ind1) .ne. e .and. fc(5,ind1) .ne. e) go to 70
	    call ifacty(ind1,npt,sizht,vcl,vm,fc,ht,type,aa,bb,cc)
	    if (type .ne. 't23') go to 70
	    vg = vm(g)
	    mu5 = tetmu(crit,vcl(1,va),vcl(1,vf),vcl(1,vg),vcl(1,vd),s)
	    mu6 = tetmu(crit,vcl(1,va),vcl(1,vf),vcl(1,vg),vcl(1,ve),s)
	    minbef = min(minbef,mu5,mu6) + tol
	    nu3 = tetmu(crit,vcl(1,vb),vcl(1,vd),vcl(1,ve),vcl(1,vf),s)
	    if (min(nu1,nu2,nu3) .le. minbef) go to 70
	    nu4 = tetmu(crit,vcl(1,va),vcl(1,vg),vcl(1,vd),vcl(1,ve),s)
	    nu5 = tetmu(crit,vcl(1,vf),vcl(1,vg),vcl(1,vd),vcl(1,ve),s)
	    if (max(nu4,nu5) .le. minbef) go to 70
	    top = ind1
	    fc(7,ind1) = ind2
	    fc(7,ind2) = ind
	    fc(7,ind) = 0
	    if (min(nu4,nu5) .le. minbef) then
	       indy = htsrc(a,b,e,npt,sizht,fc,ht)
	       indz = htsrc(a,f,e,npt,sizht,fc,ht)
	       if (indy .le. 0 .or. indz .le. 0) go to 80
	       top2 = indx
	       fc(7,indx) = indy
	       fc(7,indy) = indz
	       fc(7,indz) = ind3
	       fc(7,ind3) = 0
	       if (nu5 .lt. nu4) then
		  j = a
		  a = f
		  f = j
	       endif
 	       call fndmsw(crit,npt,sizht,vcl,vm,fc,ht,a,f,d,e,g,minbef,
     $            top,top2,impr)
	       if (ierr .ne. 0) return
	       if (.not. impr) go to 70
	    endif
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
	    j = ind2
	    ind2 = ind3
	    ind3 = j
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
	    mu1 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,vd),s)
	    mu2 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,ve),s)
	    mu3 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vf),vcl(1,vd),s)
	    mu4 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vg),vcl(1,ve),s)
	    mu5 = tetmu(crit,vcl(1,va),vcl(1,vb),vcl(1,vg),vcl(1,vf),s)
	    nu1 = tetmu(crit,vcl(1,va),vcl(1,vc),vcl(1,vd),vcl(1,ve),s)
	    nu2 = tetmu(crit,vcl(1,vb),vcl(1,vc),vcl(1,vd),vcl(1,ve),s)
	    nu3 = tetmu(crit,vcl(1,va),vcl(1,vf),vcl(1,vd),vcl(1,ve),s)
	    nu4 = tetmu(crit,vcl(1,vb),vcl(1,vf),vcl(1,vd),vcl(1,ve),s)
         else
	    nu3 = tetmu(crit,vcl(1,va),vcl(1,vf),vcl(1,vd),vcl(1,ve),s)
	    nu4 = tetmu(crit,vcl(1,vb),vcl(1,vf),vcl(1,vd),vcl(1,ve),s)
	 endif
	 if (first) then
	    minbef = min(mu1,mu2,mu3,mu4,mu5)
	    minaft = min(nu1,nu2,nu3,nu4)
	 else
	    minaft = min(nu1,nu2,nu3,nu4)
	 endif
	 if (type .eq. 't23') then
	    if (minaft .le. minbef + tol) go to 50
	    nu5 = tetmu(crit,vcl(1,va),vcl(1,vg),vcl(1,ve),vcl(1,vf),s)
	    nu6 = tetmu(crit,vcl(1,vb),vcl(1,vg),vcl(1,ve),vcl(1,vf),s)
	    if (max(nu5,nu6) .le. minbef + tol) go to 50
	    top = ind2
	    fc(7,ind2) = ind
	    fc(7,ind) = 0
	    if (min(nu5,nu6) .le. minbef + tol) then
	       indy = htsrc(a,b,e,npt,sizht,fc,ht)
	       if (indy .le. 0) go to 80
	       top2 = ind3
	       fc(7,ind3) = indy
	       fc(7,indy) = 0
	       if (nu5 .le. nu6) then
		  aa = a
		  bb = b
	       else
		  aa = b
		  bb = a
	       endif
	       dd = e
	       ee = f
	       ff = g
 	       call fndmsw(crit,npt,sizht,vcl,vm,fc,ht,aa,bb,dd,ee,ff,
     $            minbef+tol,top,top2,impr)
	       if (ierr .ne. 0) return
	       if (.not. impr) go to 50
	    endif
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
	    mu6 = tetmu(crit,vcl(1,va),vcl(1,vg),vcl(1,ve),vcl(1,vf),s)
	    nu5 = tetmu(crit,vcl(1,vb),vcl(1,vg),vcl(1,ve),vcl(1,vf),s)
	    if (min(minaft,nu5) .le. min(minbef,mu6) + tol) go to 50
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
	    mu6 = tetmu(crit,vcl(1,va),vcl(1,vg),vcl(1,vh),vcl(1,ve),s)
	    mu7 = tetmu(crit,vcl(1,va),vcl(1,vg),vcl(1,vh),vcl(1,vf),s)
	    nu5 = tetmu(crit,vcl(1,vg),vcl(1,vb),vcl(1,ve),vcl(1,vf),s)
	    nu6 = tetmu(crit,vcl(1,va),vcl(1,vh),vcl(1,ve),vcl(1,vf),s)
	    nu7 = tetmu(crit,vcl(1,vg),vcl(1,vh),vcl(1,ve),vcl(1,vf),s)
	    if (min(minaft,nu5,nu6,nu7) .le. min(minbef,mu6,mu7) + tol)
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
	 call swapmu(bndcon,crit,npt,sizht,nfc,maxfc,vcl,vm,bf,fc,ht,
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
  600 format (/1x,'imptrf: criterion =',i3)
  610 format (1x,'type ',a3,i7,' : ',5i7)
  620 format (4x,'face',3i7,' | ',2i7,' has type ',a3)
  630 format (4x,a)
      end
