      subroutine nwthfc(i,ifac,npt,sizht,nbf,nfc,maxbf,maxfc,bf,fc,ht,
     $   ntetra,hdavbf,hdavfc,front,back)
      implicit logical (a-z)
      integer back,front,hdavbf,hdavfc,i,ifac,maxbf,maxfc,nbf,nfc,npt
      integer ntetra,sizht
      integer bf(3,maxbf),fc(7,maxfc),ht(0:sizht-1)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: create new tetrahedra in 3-d triangulation from the
c        insertion of vertex i on face fc(*,ifac).
c
c     input parameters:
c        i - (local) index of next vertex inserted in triangulation;
c              it is assumed i is largest index so far
c        ifac - face containing vertex i is fc(*,ifac)
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
c
c     updated parameters:
c        nbf,nfc,bf,fc,ht,ntetra,hdavbf,hdavfc
c
c     output parameters:
c        front,back - indices of front and back of queue of interior
c              faces abc such that abci is a new tetrahedron
c
c     abnormal return:
c        ierr is set to 11, 12, or 300
c
c     routines called:
c        availf,htdel,htins,htsrc
c
      integer ierr,iprt,msglvl
      common /gerror/ ierr
      common /gprint/ iprt,msglvl
      save /gerror/,/gprint/
c
      integer a,aa,b,bb,bf1,bf2,c,cc,d,e,htsrc,ind,ind2,iv,j
      integer nbrac,nbrbc,nv
      logical bface
c
      front = 0
      back = 0
      nv = 5
      bface = (fc(5,ifac) .le. 0)
      if (bface) nv = 4
      a = fc(1,ifac)
      b = fc(2,ifac)
      c = fc(3,ifac)
      do 20 iv = 4,nv
	 ntetra = ntetra + 2
	 d = fc(iv,ifac)
	 do 10 j = 1,3
	    if (j .eq. 1) then
	       aa = a
	       bb = b
	       cc = c
	    else if (j .eq. 2) then
	       aa = b
	       bb = c
	       cc = a
	    else
	       aa = c
	       bb = a
	       cc = b
	    endif
	    call availf(hdavfc,nfc,maxfc,fc,ind)
	    if (ierr .ne. 0) return
	    call htins(ind,aa,d,i,bb,cc,npt,sizht,fc,ht)
	    ind = htsrc(aa,bb,d,npt,sizht,fc,ht)
	    if (ind .le. 0) then
	       ierr = 300
	       return
	    endif
	    if (fc(4,ind) .eq. cc) then
	       fc(4,ind) = i
	    else
	       fc(5,ind) = i
	    endif
	    if (fc(5,ind) .gt. 0) then
	       if (front .eq. 0) then
	          front = ind
	       else
		  fc(7,back) = ind
	       endif
	       back = ind
	    endif
	    if (msglvl .eq. 4) write (iprt,600) aa,bb,d,i
   10    continue
   20 continue
      if (front .ne. 0) fc(7,back) = 0
c
      call availf(hdavfc,nfc,maxfc,fc,ind)
      call availf(hdavfc,nfc,maxfc,fc,ind2)
      if (ierr .ne. 0) return
      if (bface) then
	 e = fc(5,ifac)
      else
         e = fc(4,ifac)
      endif
      call htdel(ifac,npt,sizht,fc,ht)
      call htins(ifac,a,b,i,d,e,npt,sizht,fc,ht)
      call htins(ind,a,c,i,d,e,npt,sizht,fc,ht)
      call htins(ind2,b,c,i,d,e,npt,sizht,fc,ht)
      if (bface) then
	 e = -e
	 if (hdavbf .ne. 0) then
	    bf1 = hdavbf
	    hdavbf = -bf(1,hdavbf)
	    if (hdavbf .ne. 0) then
	       bf2 = hdavbf
	       hdavbf = -bf(1,hdavbf)
	    else
	       nbf = nbf + 1
	       bf2 = nbf
	    endif
	 else
	    nbf = nbf + 2
	    bf1 = nbf - 1
	    bf2 = nbf
	 endif
	 if (nbf .gt. maxbf) then
	    nbf = maxbf
	    ierr = 12
	    return
	 endif
	 fc(5,ind) = -bf1
	 fc(5,ind2) = -bf2
	 nbrac = bf(2,e)
	 nbrbc = bf(1,e)
	 bf(1,e) = ind2
	 bf(2,e) = ind
	 bf(1,bf1) = ind2
	 bf(2,bf1) = ifac
	 bf(3,bf1) = nbrac
	 bf(1,bf2) = ind
	 bf(2,bf2) = ifac
	 bf(3,bf2) = nbrbc
	 j = -fc(5,nbrac)
	 if (bf(1,j) .eq. ifac) then
	    bf(1,j) = ind
	 else if (bf(2,j) .eq. ifac) then
	    bf(2,j) = ind
	 else
	    bf(3,j) = ind
	 endif
	 j = -fc(5,nbrbc)
	 if (bf(1,j) .eq. ifac) then
	    bf(1,j) = ind2
	 else if (bf(2,j) .eq. ifac) then
	    bf(2,j) = ind2
	 else
	    bf(3,j) = ind2
	 endif
      endif
c
  600 format (1x,'new tetra: ',4i7)
      end
