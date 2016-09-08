      program drprim
      implicit logical (a-z)
c
c     driver program for testing routine prime.
c     routines called:
c        prime
c
      integer iprt,irdr,k,prime
c
      irdr = 5
      iprt = 6
   10 continue
         write (iprt,*) 'enter integer:'
         read (irdr,*) k
	 if (k .le. 0) stop
         write (iprt,*) 'prime=',prime(k)
      go to 10
      end
