      subroutine updatf(a,b,c,d,e,i,n,p,front,back,fc,ht)
      implicit logical (a-z)
      integer a,b,back,c,d,e,front,i,n,p
      integer fc(7,*),ht(0:p-1)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: update record in fc due to a local transformation.
c        tetrahedron abcd becomes abce. add face abc to queue if it is
c        interior face, not yet in queue, and its largest index isn't i.
c    
c     input parameters:
c        a,b,c - first 3 fields of fc record (in any order)
c        d,e - fourth vertex indices of old and new tetrahedrons
c        i - vertex index determining whether face put on queue
c        n - upper bound on vertex indices
c        p - size of hash table
c        front,back - front and back pointers of queue
c        fc(1:7,1:*) - array of face records; see routine dtris3
c        ht(0:p-1) - hash table using direct chaining
c
c     updated parameters:
c        front,back,fc - queue pointers are updated
c
c     abnormal return:
c        ierr is set to 300
c
c     routines called:
c        htsrc
c
      integer ierr
      common /gerror/ ierr
      save /gerror/
c
      integer htsrc,ind
c
      ind = htsrc(a,b,c,n,p,fc,ht)
      if (ind .le. 0) then
	 ierr = 300
	 return
      endif
      if (fc(4,ind) .eq. d) then
	 fc(4,ind) = e
      else
	 fc(5,ind) = e
      endif
      if (fc(7,ind) .eq. -1 .and. fc(3,ind) .ne. i .and.
     $   fc(5,ind) .gt. 0) then
	 fc(7,ind) = 0
	 if (front .eq. 0) then
	    front = ind
	 else
	    fc(7,back) = ind
	 endif
	 back = ind
      endif
      end
