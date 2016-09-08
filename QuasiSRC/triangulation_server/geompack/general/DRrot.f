      program drrot
      implicit logical (a-z)
c
c     driver program for testing routine rotiar.
c     routines called:
c        rotiar
c
      integer maxn
      parameter (maxn = 50)
c
      integer a(0:maxn-1),i,iprt,irdr,n,shift
c
      irdr = 5
      iprt = 6
      write (iprt,*) 'enter n:'
      read (irdr,*) n
      if (n .le. 0 .or. n .gt. maxn) stop
      do 10 i = 0,n-1
	 a(i) = i
   10 continue
      write (iprt,600) (a(i),i=0,n-1)
   20 continue
         write (iprt,*) 'enter shift:'
         read (irdr,*) shift
	 if (shift .eq. 0) stop 
	 call rotiar(n,a,shift)
         write (iprt,600) (a(i),i=0,n-1)
      go to 20
c
  600 format (15i5)
      end
