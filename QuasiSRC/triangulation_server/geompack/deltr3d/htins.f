      subroutine htins(ind,a,b,c,d,e,n,p,fc,ht)
      implicit logical (a-z)
      integer a,b,c,d,e,ind,n,p
      integer fc(7,*),ht(0:p-1)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: insert record fc(1:7,ind) containing a,b,c,d,e,htlink,-1
c        into hash table ht.
c
c     input parameters:
c        ind - index of fc array
c        a,b,c,d,e - first 5 fields of fc record (or column)
c        n - upper bound on vertex indices
c        p - size of hash table
c        fc(1:7,1:*) - array of face records; see routine dtris3
c        ht(0:p-1) - hash table using direct chaining
c
c     updated parameters:
c        fc(1:7,ind) - fields of record are set
c        ht - 1 head ptr updated (record inserted at front of chain)
c
c     routines called:
c        order3
c
      integer htfun
      integer aa,bb,cc,k
c
      aa = a
      bb = b
      cc = c
      call order3(aa,bb,cc)
c      k = mod(aa*n + bb, p)
c      k = mod(k*n + cc, p)
      k = htfun(aa, bb, cc, n, p)
      fc(1,ind) = aa
      fc(2,ind) = bb
      fc(3,ind) = cc
      fc(4,ind) = d
      fc(5,ind) = e
      fc(6,ind) = ht(k)
      fc(7,ind) = -1
      ht(k) = ind
      end
