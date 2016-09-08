      real function urand(iy)
      implicit logical (a-z)
      integer iy
c
c     from forsythe, malcolm, moler, page 246.
c     urand is a uniform random number generator based on theory and
c     suggestions given in d. e. knuth (1969), vol. 2. the integer iy
c     should be initialized to an arbitrary integer prior to the first
c     call to urand. the calling program should not alter the value of
c     iy between subsequent calls to urand. values of urand will be
c     returned in the interval (0,1).
c
      integer ia,ic,itwo,m2,m,mic
      double precision halfm
      real s
      data m2/0/, itwo/2/
      data ia/0/, ic/0/, mic/0/, s/0.0/
      save ia,ic,itwo,m2,mic,s
c
      if (m2 .ne. 0) go to 20
c
c     if first entry, compute machine integer word length
c
      m = 1
   10 m2 = m
      m = itwo*m2
      if (m .gt. m2) go to 10
      halfm = m2
c
c     compute multiplier and increment for linear congruential method
c
      ia = 8*int(halfm*atan(1.0d0)/8.0d0) + 5
      ic = 2*int(halfm*(0.5d0-sqrt(3.0d0)/6.0d0)) + 1
      mic = (m2 - ic) + m2
c
c     s is the scale factor for converting to floating point
c
      s = 0.5/halfm
c
c     compute next random number
c
   20 iy = iy*ia
c
c     the following statement is for computers which do not allow
c     integer overflow on addition
c
      if (iy .gt. mic) iy = (iy - m2) - m2
c
      iy = iy + ic
c
c     the following statement is for computers where the word
c     length for addition is greater than for multiplication
c
      if (iy/2 .gt. m2) iy = (iy - m2) - m2
c
c     the following statement is for computers where integer
c     overflow affects the sign bit
c
      if (iy .lt. 0) iy = (iy + m2) + m2
      urand = real(iy)*s
      return
      end
