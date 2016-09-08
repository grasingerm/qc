      subroutine nwthou(i,npt,sizht,nbf,nfc,maxbf,maxfc,bf,fc,ht,ntetra,
     $   hdavbf,hdavfc,front,back,bfi)
      implicit logical (a-z)
      integer back,bfi,front,hdavbf,hdavfc,i,maxbf,maxfc,nbf,nfc,npt
      integer ntetra,sizht
      integer bf(3,maxbf),fc(7,maxfc),ht(0:sizht-1)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: create new tetrahedra in 3-d triangulation outside
c        convex hull by joining vertex i to visible boundary faces.
c
c     input parameters:
c        i - (local) index of next vertex inserted in triangulation
c        npt - number of 3-d points to be triangulated
c        sizht - size of hash table ht
c        nbf - number of positions used in bf array
c        nfc - number of positions used in fc array
c        maxbf - maximum size available for bf array
c        maxfc - maximum size available for fc array
c        bf(1:3,1:maxbf) -  array of boundary face records; see dtris3
c        fc(1:7,1:maxfc) - array of face records; see routine dtris3
c        ht(0:sizht-1) - hash table using direct chaining
c        ntetra - number of tetrahedra in triangulation
c        hdavbf - head pointer to available bf records
c        hdavfc - head pointer to available fc records
c        front - index of front of queue (or top of stack) of visible
c              boundary faces
c
c     updated parameters:
c        nbf,nfc,bf,fc,ht,ntetra,hdavbf,hdavfc
c
c     output parameters:
c        back - index of back of queue (or bottom of stack) of visible
c              boundary faces (which become interior faces)
c        bfi - index of fc of a boundary face containing vertex i
c
c     abnormal return:
c        ierr is set to 11, 12, or 300
c
c     routines called:
c        availf,htins,htsrc
c
      integer ierr,iprt,msglvl
      common /gerror/ ierr
      common /gprint/ iprt,msglvl
      save /gerror/,/gprint/
c
      integer a,b,bfnew,c,d,e,htsrc,ind,j,k,l,nbr,ptr
c
c     for abc in queue, form tetrahedron abci + add faces abi, aci, bci.
c     ptr, nbr, ind are indices of fc; k, l, bfnew indices of bf.
c
      bfi = 0
      ptr = front
   10 continue
	 back = ptr
	 a = fc(1,ptr)
	 b = fc(2,ptr)
	 c = fc(3,ptr)
	 k = -fc(5,ptr)
	 fc(5,ptr) = i
	 ntetra = ntetra + 1
	 if (msglvl .eq. 4) write (iprt,600) a,b,c,i
	 do 20 e = 1,3
	    if (e .eq. 2) then
	       d = b
	       b = a
	       a = d
	    else if (e .eq. 3) then
	       d = a
	       a = c
	       c = d
	    endif
	    nbr = bf(e,k)
	    if (fc(7,nbr) .ne. -1) then
	       if (fc(5,nbr) .eq. i) go to 20
	    endif
	    call availf(hdavfc,nfc,maxfc,fc,ind)
	    if (ierr .ne. 0) return
	    l = -fc(5,nbr)
	    if (bf(1,l) .eq. ptr) then
	       j = 1
	    else if (bf(2,l) .eq. ptr) then
	       j = 2
	    else
	       j = 3
	    endif
	    if (fc(7,nbr) .ne. -1) then
	       call htins(ind,b,c,i,a,fc(j,nbr),npt,sizht,fc,ht)
	    else
	       if (hdavbf .ne. 0) then
		  bfnew = hdavbf
		  hdavbf = -bf(1,hdavbf)
	       else
		  if (nbf .ge. maxbf) then
		     ierr = 12
		     return
		  endif
		  nbf = nbf + 1
		  bfnew = nbf
	       endif
	       if (bfi .eq. 0) bfi = ind
	       call htins(ind,b,c,i,a,-bfnew,npt,sizht,fc,ht)
	       bf(j,l) = ind
	       bf(3,bfnew) = nbr
	    endif
   20    continue
	 if (k .eq. nbf) then
	    nbf = nbf - 1
	 else
	    bf(1,k) = -hdavbf
	    hdavbf = k
	 endif
	 ptr = fc(7,ptr)
      if (ptr .ne. 0) go to 10
c
c     set bf(1:2,bfnew) fields for new boundary faces.
c
      ptr = bfi
      a = fc(1,ptr)
      j = 2
   30 continue
	 b = fc(j,ptr)
	 c = fc(4,ptr)
   40    continue
            nbr = htsrc(a,c,i,npt,sizht,fc,ht)
	    if (nbr .le. 0) then
	       ierr = 300
	       return
	    endif
	    if (fc(5,nbr) .gt. 0) then
	       if (fc(4,nbr) .eq. b) then
	          d = fc(5,nbr)
	       else
	          d = fc(4,nbr)
	       endif
	       b = c
	       c = d
	       go to 40
	    endif
	 k = -fc(5,ptr)
	 l = -fc(5,nbr)
	 if (fc(1,ptr) .eq. a) then
	    bf(2,k) = nbr
	 else
	    bf(1,k) = nbr
	 endif
	 if (fc(1,nbr) .eq. a) then
	    j = 1
	 else
	    j = 2
	 endif
	 bf(3-j,l) = ptr
	 a = fc(3-j,nbr)
	 ptr = nbr
      if (ptr .ne. bfi) go to 30
c
  600 format (1x,'new tetra: ',4i7)
      end
