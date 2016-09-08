      program drfcht
      implicit logical (a-z)
c
c     driver program for testing routines order3, availf, htins, htdel,
c        htsrc, updatf.
c     routines called:
c        availf, htdel, htins, htsrc, initcb, updatf
c
      integer maxfc,n,p
      parameter (maxfc = 10)
      parameter (n = 10)
      parameter (p = 7)
c
      integer ierr
      common /gerror/ ierr
c
      integer a,b,back,c,case,d,e,front,hdavfc,htsrc,i,ind,iprt,irdr,nfc
      integer fc(7,maxfc),ht(0:p-1)
c
      irdr = 5
      iprt = 6
      call initcb(0.0d0)
      do 10 i = 0,p-1
	 ht(i) = 0
   10 continue
      front = 0
      back = 0
      hdavfc = 0
      nfc = 0
   20 continue
	 write (iprt,*) 'enter case:'
	 read (irdr,*) case
	 if (case .lt. 1 .or. case .gt. 5) stop
	 if (case .eq. 1) then
 	    write (iprt,*) 'htins - enter a,b,c,d,e:'
 	    read (irdr,*) a,b,c,d,e
 	    call availf(hdavfc,nfc,maxfc,fc,ind)
	    if (ierr .ne. 0) then
	       write (iprt,*) '*** array fc is full'
	       ierr = 0
	    else
	       call htins(ind,a,b,c,d,e,n,p,fc,ht)
	    endif
	 else if (case .eq. 2) then
 	    write (iprt,*) 'htdel - enter ind:'
 	    read (irdr,*) ind
	    if (ind .le. 0 .or. ind .gt. nfc .or. fc(1,ind) .le. 0)
     $      then
	       write (iprt,*) '*** invalid index'
	    else
	       call htdel(ind,n,p,fc,ht)
	       fc(1,ind) = -hdavfc
	       hdavfc = ind
	    endif
	 else if (case .eq. 3) then
 	    write (iprt,*) 'htsrc - enter a,b,c:'
 	    read (irdr,*) a,b,c
	    ind = htsrc(a,b,c,n,p,fc,ht)
 	    write (iprt,*) 'ind=',ind
	 else if (case .eq. 4) then
 	    write (iprt,*) 'updatf - enter a,b,c,d,e,i:'
 	    read (irdr,*) a,b,c,d,e,i
	    call updatf(a,b,c,d,e,i,n,p,front,back,fc,ht)
	    if (ierr .ne. 0) then
	       write (iprt,*) '*** unsuccessful search'
	       ierr = 0
	    endif
	 else
	    write (iprt,600) hdavfc,nfc,front,back,(ht(i),i=0,p-1),
     $         (i,(fc(a,i),a=1,7),i=1,nfc)
	 endif
      go to 20
  600 format (' hdavfc=',i3,'   nfc=',i3,'   front=',i3,'   back=',i3/
     $   ' ht:  ',7i5,5x,'fc:'/(1x,8i5))
      end
