      subroutine lufac(a,lda,n,tol,ipvt,singlr)
      implicit logical (a-z)
      integer lda,n
      integer ipvt(n-1)
      double precision a(lda,n),tol
      logical singlr
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: obtain lu factorization of matrix a, i.e. apply gaussian
c        elimination with partial pivoting to a.
c
c     input parameters:
c        a(1:n,1:n) - n by n matrix to be factored
c        lda - leading dimension of array a in calling routine
c        n - order of matrix a
c        tol - relative tolerance for detecting singularity of a
c
c     output parameters:
c        a(1:n,1:n) - upper triangular matrix u and multipliers of unit
c           lower triangular matrix l (if matrix a is nonsingular)
c        ipvt(1:n-1) - pivot indices
c        singlr - .true. iff matrix is singular; this occurs when the
c              magnitude of a pivot element is <= tol*max(|a(i,j)|)
c
      integer i,j,k,kp1,m
      double precision t,tolabs
c
      if (n .lt. 1) return
      singlr = .true.
      t = 0.0d0
      do 20 j = 1,n
	 do 10 i = 1,n
	    t = max(t,abs(a(i,j)))
   10    continue
   20 continue
      tolabs = tol*t
c
      do 70 k = 1,n-1
	 kp1 = k + 1
	 m = k
	 do 30 i = kp1,n
	    if (abs(a(i,k)) .gt. abs(a(m,k))) m = i
   30    continue
	 ipvt(k) = m
	 t = a(m,k)
	 a(m,k) = a(k,k)
         a(k,k) = t
	 if (abs(t) .le. tolabs) return
	 do 40 i = kp1,n
	    a(i,k) = a(i,k)/t
   40    continue
	 do 60 j = kp1,n
	    t = a(m,j)
	    a(m,j) = a(k,j)
	    a(k,j) = t
	    if (t .eq. 0.0d0) go to 60
	    do 50 i = kp1,n
	       a(i,j) = a(i,j) - a(i,k)*t
   50       continue
   60    continue
   70 continue
      if (abs(a(n,n)) .gt. tolabs) singlr = .false.
      end
