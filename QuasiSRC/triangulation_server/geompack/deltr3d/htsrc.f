      integer function htsrc(a,b,c,n,p,fc,ht)
      implicit logical (a-z)
      integer a,b,c,n,p
      integer fc(7,*),ht(0:p-1)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: search for record fc(1:7,ind) containing key a,b,c
c        in hash table ht.
c
c     input parameters:
c        a,b,c - first 3 fields of fc record (in any order)
c        n - upper bound on vertex indices
c        p - size of hash table
c        fc(1:7,1:*) - array of face records; see routine dtris3
c        ht(0:p-1) - hash table using direct chaining
c
c     returned function value:
c        htsrc - index of fc record with key a,b,c if found,
c              or 0 if not found
c
c     routines called:
c        order3
c
      integer htfun
      integer aa,bb,cc,ind,k
c
      aa = a
      bb = b
      cc = c
      call order3(aa,bb,cc)
c      k = mod(aa*n + bb, p)
c      k = mod(k*n + cc, p)
      k = htfun(aa, bb, cc, n, p)
      ind = ht(k)
   10 continue
      if (ind .ne. 0) then
	 if (fc(1,ind) .eq. aa .and. fc(2,ind) .eq. bb .and.
     $      fc(3,ind) .eq. cc) go to 20
	 ind = fc(6,ind)
	 go to 10
      endif
   20 continue
      htsrc = ind
      end
