      subroutine nwthin(i,ifac,ivrt,npt,sizht,nfc,maxfc,fc,ht,ntetra,
     $   hdavfc,front,back)
      implicit logical (a-z)
      integer back,front,hdavfc,i,ifac,ivrt,maxfc,nfc,npt,ntetra,sizht
      integer fc(7,maxfc),ht(0:sizht-1)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: create new tetrahedra in 3-d triangulation from the
c        insertion of vertex i in interior of tetrahedron.
c
c     input parameters:
c        i - (local) index of next vertex inserted in triangulation
c        ifac - face of tetrahedron containing vertex i is fc(*,ifac)
c        ivrt - 4 or 5 where 4th vertex of tetrahedron is fc(ivrt,ifac)
c        npt - number of 3-d points to be triangulated
c        sizht - size of hash table ht
c        nfc - number of positions used in fc array
c        maxfc - maximum size available for fc array
c        fc(1:7,1:maxfc) - array of face records; see routine dtris3
c        ht(0:sizht-1) - hash table using direct chaining
c        ntetra - number of tetrahedra in triangulation
c        hdavfc - head pointer to available fc records
c
c     updated parameters:
c        nfc,fc,ht,ntetra,hdavfc
c
c     output parameters:
c        front,back - indices of front and back of queue of interior
c              faces abc such that abci is a new tetrahedron
c
c     abnormal return:
c        ierr is set to 11 or 300
c
c     routines called:
c        availf,htins,htsrc
c
      integer ierr,iprt,msglvl
      common /gerror/ ierr
      common /gprint/ iprt,msglvl
      save /gerror/,/gprint/
c
      integer a,aa,b,bb,c,cc,d,dd,htsrc,ind,indx(6),j
c
      front = 0
      back = 0
      ntetra = ntetra + 3
      a = fc(1,ifac)
      b = fc(2,ifac)
      c = fc(3,ifac)
      d = fc(ivrt,ifac)
      do 10 j = 1,4
	 if (j .eq. 1) then
	    aa = a
	    bb = b
	    cc = c
	    dd = d
	    ind = ifac
	 else
	    if (j .eq. 2) then
	       cc = d
	       dd = c
	    else if (j .eq. 3) then
	       bb = c
	       dd = b
	    else
	       aa = b
	       dd = a
	    endif
	    ind = htsrc(aa,bb,cc,npt,sizht,fc,ht)
	    if (ind .le. 0) then
	       ierr = 300
	       return
	    endif
	 endif
	 if (fc(4,ind) .eq. dd) then
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
	 if (msglvl .eq. 4) write (iprt,600) aa,bb,cc,i
   10 continue
      if (front .ne. 0) fc(7,back) = 0
      do 20 j = 1,6
	 call availf(hdavfc,nfc,maxfc,fc,indx(j))
	 if (ierr .ne. 0) return
   20 continue
      call htins(indx(1),a,b,i,c,d,npt,sizht,fc,ht)
      call htins(indx(2),a,c,i,b,d,npt,sizht,fc,ht)
      call htins(indx(3),a,d,i,b,c,npt,sizht,fc,ht)
      call htins(indx(4),b,c,i,a,d,npt,sizht,fc,ht)
      call htins(indx(5),b,d,i,a,c,npt,sizht,fc,ht)
      call htins(indx(6),c,d,i,a,b,npt,sizht,fc,ht)
c
  600 format (1x,'new tetra: ',4i7)
      end
