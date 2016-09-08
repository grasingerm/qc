      subroutine swaptf(top,npt,sizht,nfc,maxfc,vcl,vm,fc,ht,ntetra,
     $   hdavfc,front,back)
      implicit logical (a-z)
      integer back,front,hdavfc,maxfc,nfc,npt,ntetra,sizht,top
      integer fc(7,maxfc),ht(0:sizht-1),vm(npt)
      double precision vcl(3,*)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: swap transformable faces of type t23, t32, or t44 in 3-d
c        triangulation in the order given on stack.
c
c     input parameters:
c        top - pointer to stack of faces to be transformed; each face
c              should be transformable upon reaching top of stack
c        npt - number of 3-d points to be triangulated
c        sizht - size of hash table ht
c        nfc - number of positions used in fc array
c        maxfc - maximum size available for fc array
c        vcl(1:3,1:*) - vertex coordinate list
c        vm(1:npt) - indices of vertices of vcl being triangulated
c        fc(1:7,1:maxfc) - array of face records; see routine dtris3
c        ht(0:sizht-1) - hash table using direct chaining
c        ntetra - number of tetrahedra in triangulation
c        hdavfc - head pointer to available fc records
c        front,back - indices of front and back of queue of interior
c              faces for which local optimality test may be applied;
c              this routine adds faces to end of queue
c
c     updated parameters:
c        top,nfc,fc,ht,ntetra,hdavfc,front,back
c
c     abnormal return:
c        ierr is set to 11, 300, 301, 308, or 309
c
c     routines called:
c        availf,baryth,htdel,htins,htsrc,updatf
c
      integer ierr,iprt,msglvl
      common /gerror/ ierr
      common /gprint/ iprt,msglvl
      save /gerror/,/gprint/
c
      integer a,b,c,d,e,f,g,htsrc,ind,ind1,ind2,j,kneg,kzero
      integer va,vb,vc,vd,ve
      double precision alpha(4)
      logical degen
c
   10 continue
      if (top .eq. 0) return
	 ind = top
	 top = fc(7,ind)
	 if (fc(2,ind) .eq. 0) then
	    ierr = 308
	    return
	 endif
	 fc(7,ind) = -1
	 a = fc(1,ind)
	 b = fc(2,ind)
	 c = fc(3,ind)
	 d = fc(4,ind)
	 e = fc(5,ind)
	 va = vm(a)
	 vb = vm(b)
	 vc = vm(c)
	 vd = vm(d)
	 ve = vm(e)
	 if (msglvl .eq. 4) write (iprt,600) ind,a,b,c,d,e
	 call baryth(vcl(1,va),vcl(1,vb),vcl(1,vc),vcl(1,vd),vcl(1,ve),
     $      alpha,degen)
	 if (degen) then
	    ierr = 301
	    return
	 else if (alpha(4) .ge. 0.0d0) then
	    ierr = 309
            write(*,*) 'error 309 in swaptf'
            write(*,*) vcl(1,va), vcl(2,va), vcl(3,va)
            write(*,*) vcl(1,vb), vcl(2,vb), vcl(3,vb)
            write(*,*) vcl(1,vc), vcl(2,vc), vcl(3,vc)
            write(*,*) vcl(1,vd), vcl(2,vd), vcl(3,vd)
            write(*,*) vcl(1,ve), vcl(2,ve), vcl(3,ve)
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
   20    continue
	 if (kneg .eq. 1 .and. kzero .eq. 0) then
c           swap 2 tetrahedra for 3.
	    call updatf(a,b,d,c,e,0,npt,sizht,front,back,fc,ht)
	    call updatf(a,c,d,b,e,0,npt,sizht,front,back,fc,ht)
	    call updatf(b,c,d,a,e,0,npt,sizht,front,back,fc,ht)
            call updatf(a,b,e,c,d,0,npt,sizht,front,back,fc,ht)
     	    call updatf(a,c,e,b,d,0,npt,sizht,front,back,fc,ht)
     	    call updatf(b,c,e,a,d,0,npt,sizht,front,back,fc,ht)
            if (ierr .ne. 0) return
     	    call htdel(ind,npt,sizht,fc,ht)
     	    call htins(ind,a,d,e,b,c,npt,sizht,fc,ht)
     	    call availf(hdavfc,nfc,maxfc,fc,ind)
            if (ierr .ne. 0) return
     	    call htins(ind,b,d,e,a,c,npt,sizht,fc,ht)
     	    call availf(hdavfc,nfc,maxfc,fc,ind)
            if (ierr .ne. 0) return
	    call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)
	    ntetra = ntetra + 1
	    if (msglvl .eq. 4) write (iprt,610)
	 else if (kneg .eq. 2 .and. kzero .eq. 0) then
c           relabel so new tetrahedra are acde, bcde (edge ab deleted).
	    if (alpha(1) .lt. 0.0d0) then
	       j = a
	       a = c
	       c = j
	    else if (alpha(2) .lt. 0.0d0) then
	       j = b
	       b = c
	       c = j
	    endif
	    ind1 = htsrc(a,b,d,npt,sizht,fc,ht)
	    if (ind1 .le. 0) go to 50
	    if (fc(4,ind1) .eq. e .or. fc(5,ind1) .eq. e) then
c              swap 3 tetrahedra for 2.
	       call updatf(a,c,d,b,e,0,npt,sizht,front,back,fc,ht)
               call updatf(a,c,e,b,d,0,npt,sizht,front,back,fc,ht)
	       call updatf(a,d,e,b,c,0,npt,sizht,front,back,fc,ht)
	       call updatf(b,c,d,a,e,0,npt,sizht,front,back,fc,ht)
	       call updatf(b,c,e,a,d,0,npt,sizht,front,back,fc,ht)
	       call updatf(b,d,e,a,c,0,npt,sizht,front,back,fc,ht)
	       if (ierr .ne. 0) return
	       call htdel(ind,npt,sizht,fc,ht)
	       call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)
	       call htdel(ind1,npt,sizht,fc,ht)
	       if (fc(7,ind1) .ge. 0) then
		  fc(2,ind1) = 0
	       else
	          if (ind1 .eq. nfc) then
	             nfc = nfc - 1
	          else
	             fc(1,ind1) = -hdavfc
                     hdavfc = ind1
	          endif
	       endif
	       ind1 = htsrc(a,b,e,npt,sizht,fc,ht)
	       if (ind1 .le. 0) go to 50
	       call htdel(ind1,npt,sizht,fc,ht)
	       if (fc(7,ind1) .ge. 0) then
	          fc(2,ind1) = 0
	       else
	          if (ind1 .eq. nfc) then
	             nfc = nfc - 1
	          else
	             fc(1,ind1) = -hdavfc
                     hdavfc = ind1
	          endif
	       endif
	       ntetra = ntetra - 1
	       if (msglvl .eq. 4) write (iprt,620) c,d,e
	    else
	       ierr = 308
	       return
	    endif
	 else if (kneg .eq. 1 .and. kzero .eq. 1) then
c           relabel vertices so that de intersects ab.
	    if (alpha(1) .eq. 0.0d0) then
	       j = a
	       a = c
	       c = j
	    else if (alpha(2) .eq. 0.0d0) then
	       j = b
	       b = c
	       c = j
	    endif
	    ind1 = htsrc(a,b,d,npt,sizht,fc,ht)
	    ind2 = htsrc(a,b,e,npt,sizht,fc,ht)
	    if (ind1 .le. 0 .or. ind2 .le. 0) go to 50
	    if (fc(4,ind1) .eq. c) then
	       f = fc(5,ind1)
	    else
	       f = fc(4,ind1)
	    endif
	    if (fc(4,ind2) .eq. c) then
	       g = fc(5,ind2)
	    else
	       g = fc(4,ind2)
	    endif
	    if (f .le. 0 .or. g .le. 0 .or. f .ne. g) then
	       ierr = 308
	       return
	    else
c              swap 4 tetrahedra for 4.
	       call updatf(a,c,d,b,e,0,npt,sizht,front,back,fc,ht)
	       call updatf(a,c,e,b,d,0,npt,sizht,front,back,fc,ht)
	       call updatf(b,c,d,a,e,0,npt,sizht,front,back,fc,ht)
	       call updatf(b,c,e,a,d,0,npt,sizht,front,back,fc,ht)
	       call updatf(a,d,f,b,e,0,npt,sizht,front,back,fc,ht)
	       call updatf(a,e,f,b,d,0,npt,sizht,front,back,fc,ht)
	       call updatf(b,d,f,a,e,0,npt,sizht,front,back,fc,ht)
	       call updatf(b,e,f,a,d,0,npt,sizht,front,back,fc,ht)
	       if (ierr .ne. 0) return
	       call htdel(ind,npt,sizht,fc,ht)
	       call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)
	       ind = htsrc(a,b,f,npt,sizht,fc,ht)
	       if (ind .le. 0) go to 50
	       call htdel(ind,npt,sizht,fc,ht)
	       if (fc(7,ind) .ge. 0) then
		  fc(2,ind) = 0
		  call availf(hdavfc,nfc,maxfc,fc,ind)
		  if (ierr .ne. 0) return
	       endif
	       call htins(ind,d,e,f,a,b,npt,sizht,fc,ht)
	       call htdel(ind1,npt,sizht,fc,ht)
	       if (fc(7,ind1) .ge. 0) then
		  fc(2,ind1) = 0
		  call availf(hdavfc,nfc,maxfc,fc,ind1)
		  if (ierr .ne. 0) return
	       endif
	       call htins(ind1,a,d,e,c,f,npt,sizht,fc,ht)
	       call htdel(ind2,npt,sizht,fc,ht)
	       if (fc(7,ind2) .ge. 0) then
		  fc(2,ind2) = 0
		  call availf(hdavfc,nfc,maxfc,fc,ind2)
		  if (ierr .ne. 0) return
	       endif
	       call htins(ind2,b,d,e,c,f,npt,sizht,fc,ht)
	       if (front .eq. 0) then
		  front = ind1
	       else
		  fc(7,back) = ind1
	       endif
	       fc(7,ind1) = ind2
	       fc(7,ind2) = 0
	       back = ind2
	       if (msglvl .eq. 4) write (iprt,630) a,b,d,e,f
	    endif
	 else
	    ierr = 308
	    return
	 endif
      go to 10
c
   50 continue
      ierr = 300
c
  600 format (1x,'index =',i7,' : ',5i7)
  610 format (4x,'swap 2-3')
  620 format (4x,'swap 3-2 with new common face:',3i7)
  630 format (4x,'swap 4-4: edge ',2i7,' repl by ',2i7,'   f =',i7)
      end
