      program drinit
      implicit logical (a-z)
c
c     driver program for testing routine initcb.
c     routines called:
c        gtime, initcb
c
      integer ierr,iprt,msglvl
      double precision pi,tol
      common /gerror/ ierr
      common /gconst/ pi,tol
      common /gprint/ iprt,msglvl
c
      integer i,irdr,j
      real t0,t1
      double precision tolin
c
      irdr = 5
      iprt = 6
      write (iprt,*) 'enter tolin:'
      read (irdr,*) tolin
      call initcb(tolin)
      write (iprt,*) 'ierr=',ierr,'   pi=',pi,'   tol=',tol
      write (iprt,*) 'iprt=',iprt,'   msglvl=',msglvl
      call gtime(t0)
      j = 0
      do 10 i = 1,10000000
	 j = j + 1
   10 continue
      call gtime(t1)
      write (iprt,*) 'time=',t1-t0
      end
