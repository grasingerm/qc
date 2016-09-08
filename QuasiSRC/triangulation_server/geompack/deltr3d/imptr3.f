      subroutine imptr3(bndcon,postlt,crit,npt,sizht,maxfc,vcl,vm,nfc,
     $   ntetra,bf,fc,ht,nface)
      implicit logical (a-z)
      logical bndcon,postlt
      integer crit,maxfc,nface,nfc,npt,ntetra,sizht
      integer bf(3,*),fc(7,maxfc),ht(0:sizht-1),vm(npt)
      double precision vcl(3,*)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: improve a given 3-d triangulation by applying local
c        transformations based on some local criterion.
c
c     input parameters:
c        bndcon - .true. iff boundary faces are constrained (i.e. not
c              swapped by local transformations)
c        postlt - .true. iff further local transformations applied by
c              postprocessing routine imptrf
c        crit - criterion code: 1 for (local max-min) solid angle
c              criterion, 2 for radius ratio criterion, 3 for mean ratio
c              criterion, < 1 or > 3 for empty circumsphere criterion
c	 npt - number of 3-d vertices (points)
c	 sizht - size of hash table ht; a good choice is a prime number
c              which is about 1/8 * nface (or 3/2 * npt for random
c              points from the uniform distribution)
c        maxfc - maximum size available for fc array
c        vcl(1:3,1:*) - vertex coordinate list
c        vm(1:npt) - indices of vertices of vcl being triangulated
c        nfc - number of positions used in fc array
c        ntetra - number of tetrahedra in triangulation
c        bf(1:3,1:*) -  array of boundary face records; see dtris3
c        fc(1:7,1:maxfc) - array of face records; see routine dtris3
c        ht(0:sizht-1) - hash table using direct chaining
c        [note: bf, fc, ht should be as output by dtris3 or dtriw3.]
c
c     updated parameters:
c        nfc,ntetra,bf,fc,ht
c
c     output parameters:
c        nface - number of faces in triangulation; nface <= nfc
c
c     abnormal return:
c        ierr is set to 11, 300, 301, 305, 308, or 309
c
c     routines called:
c        imptrd,imptrf,swapes,swapmu
c
      integer ierr,iprt,msglvl
      common /gerror/ ierr
      common /gprint/ iprt,msglvl
      save /gerror/,/gprint/
c
      integer back,front,hdavbf,hdavfc,i,ptr
c
c     create initial queue of interior faces.
c
      hdavbf = fc(7,1)
      hdavfc = fc(7,2)
      fc(7,1) = -1
      fc(7,2) = -1
      front = 0
      do 10 i = 1,nfc
	 if (fc(1,i) .gt. 0 .and. fc(5,i) .gt. 0) then
	    if (front .eq. 0) then
	       front = i
	    else
	       fc(7,back) = i
	    endif
	    back = i
	 endif
   10 continue
      if (front .ne. 0) fc(7,back) = 0
      if (msglvl .eq. 4) write (iprt,600) crit
c
      if (crit .ge. 1 .and. crit .le. 3) then
	 call swapmu(bndcon,crit,npt,sizht,nfc,maxfc,vcl,vm,bf,fc,ht,
     $      ntetra,hdavfc,front,back,i)
      else
	 call swapes(bndcon,0,npt,sizht,nfc,maxfc,vcl,vm,bf,fc,ht,
     $      ntetra,hdavfc,front,back,i)
      endif
      if (ierr .ne. 0) return
      if (crit .ge. 1 .and. crit .le. 3 .and. postlt) then
	 call imptrf(bndcon,crit,npt,sizht,maxfc,vcl,vm,nfc,ntetra,
     $      hdavfc,bf,fc,ht)
      else if (postlt) then
	 call imptrd(bndcon,npt,sizht,maxfc,vcl,vm,nfc,ntetra,hdavfc,
     $      bf,fc,ht)
      endif
      if (ierr .ne. 0) return
c
      nface = nfc
      ptr = hdavfc
   20 continue
      if (ptr .eq. 0) go to 30
	 nface = nface - 1
	 ptr = -fc(1,ptr)
	 go to 20
   30 continue
      fc(7,1) = hdavbf
      fc(7,2) = hdavfc
c
  600 format (/1x,'imptr3: criterion =',i3)
      end
