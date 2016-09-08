      double precision function tetmu(crit,a,b,c,d,s)
      implicit logical (a-z)
      integer crit
      double precision a(3),b(3),c(3),d(3),s(4)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: compute tetrahedron shape measure, scaled so that
c        largest possible value is 1.
c
c     input parameters:
c        crit - criterion number: 1 for sin(min solid angle / 2),
c              2 for radius ratio, 3 for mean ratio
c        a(1:3),b(1:3),c(1:3),d(1:3) - 4 vertices of tetrahedron
c
c     working parameter:
c        s(1:4) - needed for sangmn routine
c
c     returned function value:
c        tetmu - sangmn(abcd)*sf, radrth(abcd), or emnrth(abcd)
c              if crit = 1, 2, or 3 resp.; else 1.0
c
c     routines called:
c        emnrth,radrth,sangmn
c
      double precision emnrth,radrth,sangmn,sf
      data sf/3.674234614174767/
      save sf
c
      if (crit .eq. 1) then
	 tetmu = sangmn(a,b,c,d,s)*sf
      else if (crit .eq. 2) then
	 tetmu = radrth(a,b,c,d)
      else if (crit .eq. 3) then
	 tetmu = emnrth(a,b,c,d)
      else
         tetmu = 1.0d0
      endif
      end
