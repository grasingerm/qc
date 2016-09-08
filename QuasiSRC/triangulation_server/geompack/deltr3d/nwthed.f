      subroutine nwthed(i,ifac,iedg,npt,sizht,nbf,nfc,maxbf,maxfc,bf,fc,
     $   ht,ntetra,hdavbf,hdavfc,front,back)
      implicit logical (a-z)
      integer back,front,hdavbf,hdavfc,i,iedg,ifac,maxbf,maxfc,nbf,nfc
      integer npt,ntetra,sizht
      integer bf(3,maxbf),fc(7,maxfc),ht(0:sizht-1)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: create new tetrahedra in 3-d triangulation from the
c        insertion of vertex i on edge fc(iedg,ifac).
c
c     input parameters:
c        i - (local) index of next vertex inserted in triangulation;
c              it is assumed i is largest index so far
c        ifac,iedg - edge containing i is fc(iedg,ifac); 1 <= iedg <= 3
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
      integer a,aa,b,bb,bfn(2),bfo(2),c,cp,csav,d,e,fcn(2),fco(2),htsrc
      integer ind,inew,j,k,l,nbrac,nbrbc,nfcin
      logical bedge
c
      front = 0
      back = 0
      if (iedg .eq. 1) then
	 a = fc(1,ifac)
	 b = fc(2,ifac)
	 c = fc(3,ifac)
      else if (iedg .eq. 2) then
	 a = fc(2,ifac)
	 b = fc(3,ifac)
	 c = fc(1,ifac)
      else
	 a = fc(1,ifac)
	 b = fc(3,ifac)
	 c = fc(2,ifac)
      endif
      csav = c
      d = fc(4,ifac)
c
c     determine faces incident on edge ab in circular order, and store
c     their indices at end of fc array.
c
      bedge = .false.
      ind = ifac
      k = 0
   10 continue
	 k = k + 1
	 if (nfc + k .gt. maxfc) then
	    ierr = 11
	    return
	 endif
         fc(1,nfc+k) = ind
	 if (bedge) go to 20
	 if (d .eq. csav) go to 50
	 ind = htsrc(a,b,d,npt,sizht,fc,ht)
	 if (ind .le. 0) then
	    ierr = 300
	    return
	 endif
	 if (fc(5,ind) .le. 0) then
	    bedge = .true.
	    cp = d
	 else if (fc(4,ind) .eq. c) then
	    c = d
	    d = fc(5,ind)
	 else
	    c = d
	    d = fc(4,ind)
	 endif
      go to 10
c
c     determine further faces in case of ab being a boundary edge.
c
   20 continue
      if (fc(5,ifac) .le. 0) go to 50
      l = k
      do 30 j = 1,k/2
	 e = fc(1,nfc+j)
	 fc(1,nfc+j) = fc(1,nfc+l)
	 fc(1,nfc+l) = e
	 l = l - 1
   30 continue
      c = csav
      csav = cp
      d = fc(5,ifac)
   40 continue
	 ind = htsrc(a,b,d,npt,sizht,fc,ht)
	 if (ind .le. 0) then
	    ierr = 300
	    return
	 endif
	 k = k + 1
	 if (nfc + k .gt. maxfc) then
	    ierr = 11
	    return
	 endif
         fc(1,nfc+k) = ind
	 if (fc(5,ind) .le. 0) then
	    go to 50
	 else if (fc(4,ind) .eq. c) then
	    c = d
	    d = fc(5,ind)
	 else
	    c = d
	    d = fc(4,ind)
	 endif
      go to 40
c
c     create new faces and tetrahedra, and add faces to queue.
c
   50 continue
      nfcin = nfc
      nfc = nfc + k
      ntetra = ntetra + k
      if (bedge) then
	 ntetra = ntetra - 1
	 fcn(1) = nfcin + 1
	 fcn(2) = nfcin + k
	 fco(1) = fc(1,nfcin+1)
	 fco(2) = fc(1,nfcin+k)
	 bfo(1) = -fc(5,fco(1))
	 bfo(2) = -fc(5,fco(2))
      endif
      do 70 j = 1,k
	 inew = nfcin + j
	 ind = fc(1,inew)
	 if (fc(1,ind) .eq. a) then
	    if (fc(2,ind) .eq. b) then
	       c = fc(3,ind)
	    else
	       c = fc(2,ind)
	    endif
	 else
	    c = fc(1,ind)
	 endif
	 d = fc(4,ind)
	 e = fc(5,ind)
	 call htdel(ind,npt,sizht,fc,ht)
	 call htins(ind,a,c,i,d,e,npt,sizht,fc,ht)
	 call htins(inew,b,c,i,d,e,npt,sizht,fc,ht)
	 if (j .eq. k) then
	    if (bedge) go to 70
	    d = csav
	 else if (j .gt. 1) then
	    if (d .eq. cp) d = e
	 endif
	 call availf(hdavfc,nfc,maxfc,fc,ind)
	 if (ierr .ne. 0) return
	 call htins(ind,c,d,i,a,b,npt,sizht,fc,ht)
	 aa = a
	 bb = b
	 do 60 l = 1,2
	    if (l .eq. 2) then
	       aa = b
	       bb = a
	    endif
	    ind = htsrc(aa,c,d,npt,sizht,fc,ht)
	    if (ind .le. 0) then
	       ierr = 300
	       return
	    endif
	    if (fc(5,ind) .le. 0) then
	       fc(4,ind) = i
	    else
	       if (fc(4,ind) .eq. bb) then
		  fc(4,ind) = i
	       else
		  fc(5,ind) = i
	       endif
	       if (front .eq. 0) then
	          front = ind
	       else
	          fc(7,back) = ind
	       endif
	       back = ind
	    endif
	    if (msglvl .eq. 4) write (iprt,600) aa,c,d,i
   60    continue
	 cp = c
   70 continue
      if (front .ne. 0) fc(7,back) = 0
c
      if (bedge) then
	 d = c
	 c = csav
	 if (hdavbf .ne. 0) then
	    bfn(1) = hdavbf
	    hdavbf = -bf(1,hdavbf)
	    if (hdavbf .ne. 0) then
	       bfn(2) = hdavbf
	       hdavbf = -bf(1,hdavbf)
	    else
	       nbf = nbf + 1
	       bfn(2) = nbf
	    endif
	 else
	    nbf = nbf + 2
	    bfn(1) = nbf - 1
	    bfn(2) = nbf
	 endif
	 if (nbf .gt. maxbf) then
	    nbf = maxbf
	    ierr = 12
	    return
	 endif
	 fc(5,nfcin+1) = -bfn(1)
	 fc(5,nfcin+k) = -bfn(2)
	 do 80 j = 1,2
	    if (j .eq. 2) c = d
	    if (c .lt. a) then
	       nbrac = bf(3,bfo(j))
	       nbrbc = bf(2,bfo(j))
	       bf(1,bfo(j)) = fco(3-j)
	       bf(2,bfo(j)) = fcn(j)
	       bf(1,bfn(j)) = fcn(3-j)
	       bf(2,bfn(j)) = fco(j)
	    else if (c .lt. b) then
	       nbrac = bf(3,bfo(j))
	       nbrbc = bf(1,bfo(j))
	       bf(1,bfo(j)) = fcn(j)
	       bf(2,bfo(j)) = fco(3-j)
	       bf(1,bfn(j)) = fcn(3-j)
	       bf(2,bfn(j)) = fco(j)
	    else
	       nbrac = bf(2,bfo(j))
	       nbrbc = bf(1,bfo(j))
	       bf(1,bfo(j)) = fcn(j)
	       bf(2,bfo(j)) = fco(3-j)
	       bf(1,bfn(j)) = fco(j)
	       bf(2,bfn(j)) = fcn(3-j)
	    endif
	    bf(3,bfo(j)) = nbrac
	    bf(3,bfn(j)) = nbrbc
	    l = -fc(5,nbrbc)
	    if (bf(1,l) .eq. fco(j)) then
	       bf(1,l) = fcn(j)
	    else if (bf(2,l) .eq. fco(j)) then
	       bf(2,l) = fcn(j)
	    else
	       bf(3,l) = fcn(j)
	    endif
   80    continue
      endif
c
  600 format (1x,'new tetra: ',4i7)
      end
