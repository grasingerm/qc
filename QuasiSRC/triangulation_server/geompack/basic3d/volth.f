      double precision function volth(a,b,c,d)
      implicit logical (a-z)
      double precision a(3),b(3),c(3),d(3)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: compute 6 times volume of tetrahedron.
c
c     input parameters:
c        a(1:3),b(1:3),c(1:3),d(1:3) - 4 vertices of tetrahedron
c
c     returned function value:
c        volth - 6 * volume of tetrahedron
c
      integer i
      double precision u(3),v(3),w(3)
c
      do 10 i = 1,3
	 u(i) = b(i) - a(i)
	 v(i) = c(i) - a(i)
	 w(i) = d(i) - a(i)
   10 continue
      volth = abs( u(1)*(v(2)*w(3) - v(3)*w(2)) + u(2)*(v(3)*w(1) -
     $   v(1)*w(3)) + u(3)*(v(1)*w(2) - v(2)*w(1)) )
      end
