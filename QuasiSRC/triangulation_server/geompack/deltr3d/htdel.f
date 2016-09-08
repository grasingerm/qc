      subroutine htdel(ind,n,p,fc,ht)
      implicit logical (a-z)
      integer ind,n,p
      integer fc(7,*),ht(0:p-1)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: delete record fc(1:7,ind) from hash table ht.
c
c     input parameters:
c        ind - index of fc array
c        n - upper bound on vertex indices
c        p - size of hash table
c        fc(1:7,1:*) - array of face records; see routine dtris3
c        ht(0:p-1) - hash table using direct chaining
c
c     updated parameters:
c        fc,ht - one link is updated
c
      integer htfun
      integer k,ptr
c
c      k = mod(fc(1,ind)*n + fc(2,ind), p)
c      k = mod(k*n + fc(3,ind), p)
      k = htfun(fc(1,ind), fc(2,ind), fc(3,ind), n, p)
      ptr = ht(k)
      if (ptr .eq. ind) then
	 ht(k) = fc(6,ind)
      else
   10    continue
	 if (fc(6,ptr) .ne. ind) then
	    ptr = fc(6,ptr)
	    go to 10
	 endif
	 fc(6,ptr) = fc(6,ind)
      endif
      end
