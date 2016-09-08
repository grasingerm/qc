      subroutine lusol(a,lda,n,ipvt,b)
      implicit logical (a-z)
      integer lda,n
      integer ipvt(n-1)
      double precision a(lda,n),b(n)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: solve linear system a*x = b given lu factorization of a.
c        it is assumed that a is nonsingular.
c
c     input parameters:
c        a(1:n,1:n) - contains factors l, u output from routine lufac
c        lda - leading dimension of array a in calling routine
c        n - order of matrix a
c        ipvt(1:n-1) - pivot indices from routine lufac
c        b(1:n) - right hand side vector
c
c     output parameters:
c        b(1:n) - solution vector x
c
      integer i,k,m
      double precision t
c
c     forward elimination
c
      if (n .le. 1) go to 50
      do 20 k = 1,n-1
	 m = ipvt(k)
	 t = b(m)
	 b(m) = b(k)
	 b(k) = t
	 do 10 i = k+1,n
	    b(i) = b(i) - a(i,k)*t
   10    continue
   20 continue
c
c     back substitution
c
      do 40 k = n,2,-1
	 t = b(k)/a(k,k)
	 b(k) = t
	 do 30 i = 1,k-1
	    b(i) = b(i) - a(i,k)*t
   30    continue
   40 continue
   50 continue
      b(1) = b(1)/a(1,1)
      end
