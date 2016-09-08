      logical function iless(k,p,q)
      implicit logical (a-z)
      integer k
      integer p(k),q(k)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: determine whether p is lexicographically less than q.
c
c     input parameters:
c        k - dimension of points
c        p(1:k),q(1:k) - two k-dimensional integer points
c
c     returned function value:
c        iless - .true. if p < q, .false. otherwise
c
      integer i
c
      do 10 i = 1,k
	 if (p(i) .eq. q(i)) go to 10
         if (p(i) .lt. q(i)) then
	    iless = .true.
         else
	    iless = .false.
	 endif
	 return
   10 continue
      iless = .false.
      end
