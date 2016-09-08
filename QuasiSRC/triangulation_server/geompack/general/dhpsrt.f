      subroutine dhpsrt(k,n,lda,a,map)
      implicit logical (a-z)
      integer k,lda,n
      integer map(n)
      double precision a(lda,*)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: use heapsort to obtain the permutation of n k-dimensional
c        double precision points so that the points are in lexicographic
c        increasing order.
c
c     input parameters:
c        k - dimension of points
c        n - number of points
c        lda - leading dimension of array a in calling routine; should
c              be >= k
c        a(1:k,1:*) - array of >= n k-d double precision points
c        map(1:n) - the points of a with indices map(1), map(2), ...,
c              map(n) are to be sorted
c
c     updated parameters:
c        map(1:n) - elements are permuted so that a(*,map(1)) <=
c              a(*,map(2)) <= ... <= a(*,map(n))
c
c     routines called:
c        dsftdw
c
      integer i,t
c
      do 10 i = n/2,1,-1
         call dsftdw(i,n,k,lda,a,map)
   10 continue
      do 20 i = n,2,-1
	 t = map(1)
	 map(1) = map(i)
	 map(i) = t
         call dsftdw(1,i-1,k,lda,a,map)
   20 continue
      end
