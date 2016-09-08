      program drsort
      implicit logical (a-z)
c
c     driver program for testing routines randpt, dhpsrt, ihpsrt.
c     routines called:
c        dhpsrt, ihpsrt, initcb, randpt
c
      integer maxk,maxn
      parameter (maxk = 4)
      parameter (maxn = 100)
c
      integer axis,i,iprt,irdr,j,k,n,nptav,seed
      integer ia(maxk,maxn),map(maxn)
      double precision da(maxk,maxn),scale(maxk),trans(maxk)
      logical iflag
      data scale/1.0d0,1.0d0,1.0d0,1.0d0/
      data trans/0.0d0,0.0d0,0.0d0,0.0d0/
c
      irdr = 5
      iprt = 6
      call initcb(0.0d0)
      read (irdr,*) k,n,seed,axis,nptav
      iflag = (k .lt. 0)
      k = abs(k)
      if (k .lt. 1 .or. k .gt. maxk .or. n .lt. 1 .or. n .gt. maxn) stop
      do 10 j = 1,n
	 map(j) = j
   10 continue
      call randpt(k,n,seed,axis,nptav,scale,trans,maxk,da)
      if (iflag) then
	 do 30 j = 1,n
	    do 20 i = 1,k
	       ia(i,j) = int(n*da(i,j))
   20       continue
   30    continue
	 do 40 j = 1,n
	    write (iprt,610) j,(ia(i,j),i=1,k)
   40    continue
	 call ihpsrt(k,n,maxk,ia,map)
	 write (iprt,600)
	 do 50 j = 1,n
	    write (iprt,610) map(j),(ia(i,map(j)),i=1,k)
   50    continue
      else
	 do 60 j = 1,n
	    write (iprt,620) j,(da(i,j),i=1,k)
   60    continue
	 call dhpsrt(k,n,maxk,da,map)
	 write (iprt,600)
	 do 70 j = 1,n
	    write (iprt,620) map(j),(da(i,map(j)),i=1,k)
   70    continue
      endif
c
  600 format (1x)
  610 format (1x,5i5)
  620 format (1x,i5,4f15.7)
      end
