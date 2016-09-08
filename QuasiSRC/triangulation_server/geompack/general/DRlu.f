      program drlu
      implicit logical (a-z)
c
c     driver program for testing routines lufac, lusol.
c     routines called:
c        gtime, initcb, lufac, lusol, urand
c
      double precision pi,tol
      common /gconst/ pi,tol
c
      integer nmax
      parameter (nmax = 200)
c
      integer i,iprt,ipvt(nmax),irdr,j,n,seed
      double precision a(nmax,nmax),b(nmax),emax,esum,t,tolin
      real t0,tf,ts,urand
      logical singlr
c
      irdr = 5
      iprt = 6
      write (iprt,*) 'enter tolin,n,seed:'
      read (irdr,*) tolin,n,seed
      if (n .le. 0 .or. n .gt. nmax) stop
      call initcb(tolin)
      do 10 i = 1,n
	 b(i) = 0.0d0
   10 continue
      do 30 j = 1,n
	 do 20 i = 1,n
	    a(i,j) = dble(urand(seed)*2.0e0 - 1.0e0)
	    b(i) = b(i) + a(i,j)
   20    continue
   30 continue
      if (n .le. 4) then
	 write (iprt,600) 'a, b'
	 do 40 i = 1,n
	    write (iprt,610) (a(i,j),j=1,n),b(i)
   40    continue
      endif
c
      call gtime(t0)
      call lufac(a,nmax,n,tol,ipvt,singlr)
      call gtime(tf)
      tf = tf - t0
      if (singlr) then
	 write (iprt,600) '*** matrix is singular'
	 stop
      endif
      call gtime(t0)
      call lusol(a,nmax,n,ipvt,b)
      call gtime(ts)
      ts = ts - t0
      if (n .le. 4) then
	 write (iprt,600) 'ipvt, lu'
	 ipvt(n) = n
	 do 50 i = 1,n
	    write (iprt,620) ipvt(i),(a(i,j),j=1,n)
   50    continue
      endif
      emax = 0.0d0
      esum = 0.0d0
      do 60 i = 1,n
	 t = abs(b(i) - 1.0d0)
	 emax = max(emax,t)
	 esum = esum + t
   60 continue
      write (iprt,630) (b(i),i=1,min(4,n))
      write (iprt,640) emax,esum,tf,ts
c
  600 format (1x,a)
  610 format (1x,5f15.7)
  620 format (1x,i5,4f15.7)
  630 format (1x,'x = ',4f15.7)
  640 format (1x,'emax,esum = ',2e15.7/1x,'tf, ts = ',2f9.3)
      end
