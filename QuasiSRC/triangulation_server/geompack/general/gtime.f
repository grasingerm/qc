      subroutine gtime(time)
      implicit logical (a-z)
      real time
c
c     purpose: get current cpu time in seconds from c routine clock_.
c
      integer clock
      time = clock()/60.0
      return
      end
