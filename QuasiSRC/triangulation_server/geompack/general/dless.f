      logical function dless(k,p,q)
      implicit logical (a-z)
      integer k
      double precision p(k),q(k)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: determine whether p is lexicographically less than q in
c        floating point arithmetic?
c
c     input parameters:
c        k - dimension of points
c        p(1:k),q(1:k) - two k-dimensional double precision points
c
c     returned function value:
c        dless - .true. if p < q, .false. otherwise
c
      double precision pi,tol
      common /gconst/ pi,tol
      save /gconst/
c
      integer i
      double precision cmax
c
      do 10 i = 1,k
	 cmax = max(abs(p(i)),abs(q(i)))
	 if (abs(p(i) - q(i)) .le. tol*cmax .or. cmax .le. tol) go to 10
         if (p(i) .lt. q(i)) then
	    dless = .true.
         else
	    dless = .false.
	 endif
	 return
   10 continue
      dless = .false.
      end
