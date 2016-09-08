      integer function prime(k)
      implicit logical (a-z)
      integer k
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: return a prime >= k (if possible) from internal array
c        of primes. more primes can be added if desired.
c
c     input parameters:
c        k - positive integer
c
c     returned function value:
c        prime - smallest prime >= k from internal array (or largest
c              in array)
c
      integer nprime
      parameter (nprime = 150)
c
      integer l,m,u
c
      integer primes(nprime)
      data primes/17,31,47,61,79,97,113,127,149,163,179,193,211,227,241,
     $   257,271,293,307,331,353,379,401,431,457,479,503,541,563,587,
     $   613,641,673,701,727,751,773,797,821,853,877,907,929,953,977,
     $   1009,1049,1087,1123,1163,1201,1237,1277,1319,1361,1399,1433,
     $   1471,1511,1543,1579,1613,1657,1699,1741,1783,1831,1873,1931,
     $   1973,2017,2069,2129,2203,2267,2333,2389,2441,2503,2557,2609,
     $   2663,2719,2789,2851,2917,2999,3061,3137,3209,3299,3371,3449,
     $   3527,3613,3697,3779,3863,3947,4049,4211,4421,4621,4813,5011,
     $   5227,5413,5623,5813,6011,6211,6421,6619,6823,7013,7211,7411,
     $   7621,7817,8011,8219,8419,8623,8819,9011,9221,9413,9613,9811,
     $   10037,10211,10427,10613,10831,11027,11213,11411,11617,11813,
     $   12011,12211,12413,12611,12821,13033,13217,13411,13613,13829,
     $   14011/
      save primes
c
      if (k .le. primes(1)) then
	 prime = primes(1)
	 return
      else if (k .ge. primes(nprime)) then
	 prime = primes(nprime)
	 return
      endif
c
c     use binary search to find prime >= k.
c
      l = 1
      u = nprime
   10 continue
	 m = (l + u)/2
	 if (k .lt. primes(m)) then
	    u = m - 1
	 else if (k .gt. primes(m)) then
	    l = m + 1
	 else
	    prime = primes(m)
	    return
	 endif
      if (l .le. u) go to 10
      prime = primes(u+1)
      end
