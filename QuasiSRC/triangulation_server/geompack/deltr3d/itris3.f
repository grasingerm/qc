      subroutine itris3(npt,sizht,maxbf,maxfc,vcl,vm,nbf,nfc,nface,
     $   ntetra,bf,fc,ht)
      implicit logical (a-z)
      integer maxbf,maxfc,nbf,nface,nfc,npt,ntetra,sizht
      integer bf(3,maxbf),fc(7,maxfc),ht(0:sizht-1),vm(npt)
      double precision vcl(3,*)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: construct initial triangulation of 3-d vertices by first
c        sorting them in lexicographically increasing (x,y,z) order and
c        then inserting 1 vertex at a time from outside the convex hull.
c
c     input parameters:
c	 npt - number of 3-d vertices (points)
c	 sizht - size of hash table ht; a good choice is a prime number
c              which is about 1/8 * nface (or 3/2 * npt for random
c              points from the uniform distribution)
c        maxbf - maximum size available for bf array
c        maxfc - maximum size available for fc array
c        vcl(1:3,1:*) - vertex coordinate list
c        vm(1:npt) - indices of vertices of vcl being triangulated
c
c     updated parameters:
c        vm(1:npt) - indices are permuted, so that vcl(*,vm(1)), ... ,
c              vcl(*,vm(npt)) are in lexicographic increasing order,
c              with possible slight reordering so first 4 vertices are
c              non-coplanar
c
c     output parameters:
c        nbf - number of positions used in bf array; nbf <= maxbf
c        nfc - number of positions used in fc array; nfc <= maxfc
c        nface - number of faces in triangulation; nface <= nfc
c        ntetra - number of tetrahedra in triangulation
c        bf(1:3,1:nbf) -  array of boundary face records containing ptrs
c              (indices) to fc; if fc(5,i) = -j < 0 and fc(1:3,i) = abc,
c              then bf(1,j) points to other boundary face with edge bc,
c              bf(2,j) points to other boundary face with edge ac, and
c              bf(3,j) points to other boundary face with edge ab;
c              if bf(1,j) <= 0, record is not used and is in avail list
c        fc(1:7,1:nfc) - array of face records which are in linked lists
c              in hash table with direct chaining. fields are:
c              fc(1:3,*) - a,b,c with 1<=a<b<c<=npt; indices in vm of 3
c                 vertices of face; if a <= 0, record is not used (it is
c                 in linked list of avail records with indices <= nfc);
c                 internal use: if b <= 0, face in queue, not in triang
c              fc(4:5,*) - d,e; indices in vm of 4th vertex of 1 or 2
c                 tetrahedra with face abc; if abc is boundary face
c                 then e < 0 and |e| is an index of bf array
c              fc(6,*) - htlink; pointer (index in fc) of next element
c                 in linked list (or null = 0)
c              fc(7,*) - used internally for qlink (link for queues or
c                 stacks); pointer (index in fc) of next face in queue/
c                 stack (or null = 0); qlink = -1 indicates face is not
c                 in any queue/stack, and is output value (for records
c                 not in avail list), except:
c        fc(7,1:2) - hdavbf,hdavfc : head ptrs of avail list in bf, fc
c        ht(0:sizht-1) - hash table using direct chaining; entries are
c              head pointers of linked lists (indices of fc array)
c              containing the faces and tetrahedra of triangulation
c
c     abnormal return:
c        ierr is set to 11, 12, 300, 301, 302, 303, 304, or 306
c
c     routines called:
c        dhpsrt,frstet,htins,nwthou,opside,vbfac
c
      integer ierr,iprt,msglvl
      common /gerror/ ierr
      common /gprint/ iprt,msglvl
      save /gerror/,/gprint/
c
      integer a,b,back,bfi,e,front,hdavbf,hdavfc,i,i3,i4,ip,j,k,op
      integer opside,ptr,topnv,top,va,vb,vc,vi
      double precision ctr(3)
c
c     permute elements of vm so that vertices are in lexicographic
c     order, and initialize data structures.
c
      call dhpsrt(3,npt,3,vcl,vm)
      call frstet(.true.,npt,vcl,vm,i3,i4)
      if (ierr .ne. 0) return
      do 10 i = 1,3
	 ctr(i) = ( vcl(i,vm(1)) + vcl(i,vm(2)) + vcl(i,vm(3)) +
     $      vcl(i,vm(4)) )/4.0d0
   10 continue
      do 20 i = 0,sizht-1
	 ht(i) = 0
   20 continue
      hdavbf = 0
      hdavfc = 0
      nbf = 4
      nfc = 4
      ntetra = 1
      call htins(1,1,2,3,4,-1,npt,sizht,fc,ht)
      call htins(2,1,2,4,3,-2,npt,sizht,fc,ht)
      call htins(3,1,3,4,2,-3,npt,sizht,fc,ht)
      call htins(4,2,3,4,1,-4,npt,sizht,fc,ht)
      bf(1,1) = 4
      bf(2,1) = 3
      bf(3,1) = 2
      bf(1,2) = 4
      bf(2,2) = 3
      bf(3,2) = 1
      bf(1,3) = 4
      bf(2,3) = 2
      bf(3,3) = 1
      bf(1,4) = 3
      bf(2,4) = 2
      bf(3,4) = 1
      if (msglvl .eq. 4) write (iprt,600) (vm(i),i=1,4),i3,i4
c
c     insert ith vertex into triangulation of first i-1 vertices.
c
      do 120 i = 5,npt
	 vi = vm(i)
	 if (msglvl .eq. 4) write (iprt,610) i,vi
	 ip = i - 1
	 if (i .eq. 5) ip = 2
	 if (i .eq. i3 + 2) ip = 3
	 if (i .eq. i4 + 1) ip = 4
c
c        form stacks of boundary faces involving vertex ip.
c        top is for stack of boundary faces to be tested for visibility.
c        front is for stack of boundary faces visible from vertex i.
c        topnv is for stack of boundary faces not visible from i.
c
         front = 0
	 topnv = 0
	 if (i .eq. 5) then
	    top = 4
	    a = 3
	    b = 2
	    if (ip .eq. 2) a = 2
	    if (ip .le. 3) b = 1
	    fc(7,top) = a
	    fc(7,a) = b
	    fc(7,b) = 0
	 else if (ip .eq. i - 1) then
	    top = bfi
	    fc(7,bfi) = 0
	    b = fc(2,bfi)
	    ptr = bf(1,-fc(5,bfi))
   30       continue
	       if (fc(1,ptr) .eq. b) then
		  b = fc(2,ptr)
		  j = 1
	       else
		  b = fc(1,ptr)
		  j = 2
	       endif
	       fc(7,ptr) = top
	       top = ptr
	       ptr = bf(j,-fc(5,ptr))
	    if (ptr .ne. bfi) go to 30
	 else
	    do 50 k = 1,nbf
	       if (bf(1,k) .le. 0) go to 50
	       do 40 e = 1,3
		  ptr = bf(e,k)
		  if (fc(1,ptr) .eq. ip) then
		     b = fc(2,ptr)
		     j = 3
		     go to 60
		  else if (fc(2,ptr) .eq. ip) then
		     b = fc(1,ptr)
		     j = 3
		     go to 60
                  else if (fc(3,ptr) .eq. ip) then
		     b = fc(1,ptr)
		     j = 2
		     go to 60
		  endif
   40          continue
   50       continue
   60       continue
	    bfi = ptr
	    top = bfi
	    fc(7,bfi) = 0
	    ptr = bf(j,-fc(5,bfi))
   70       continue
	       if (fc(1,ptr) .eq. b) then
		  j = 1
		  if (fc(2,ptr) .eq. ip) then
		     b = fc(3,ptr)
		  else
		     b = fc(2,ptr)
		  endif
	       else if (fc(2,ptr) .eq. b) then
		  j = 2
		  if (fc(1,ptr) .eq. ip) then
		     b = fc(3,ptr)
		  else
		     b = fc(1,ptr)
		  endif
	       else
		  j = 3
		  if (fc(1,ptr) .eq. ip) then
		     b = fc(2,ptr)
		  else
		     b = fc(1,ptr)
		  endif
	       endif
	       fc(7,ptr) = top
	       top = ptr
	       ptr = bf(j,-fc(5,ptr))
	    if (ptr .ne. bfi) go to 70
	 endif
c
c        find a boundary face visible from vertex i.
c
   80    continue
	 if (top .eq. 0) go to 110
	    ptr = top
	    top = fc(7,ptr)
	    va = vm(fc(1,ptr))
	    vb = vm(fc(2,ptr))
	    vc = vm(fc(3,ptr))
	    op = opside(vcl(1,va),vcl(1,vb),vcl(1,vc),ctr,vcl(1,vi))
	    if (op .eq. 2) then
	       ierr = 301
	       return
	    endif
	    if (op .eq. 1) then
	       front = ptr
   90          continue
	       if (top .eq. 0) go to 110
		  ptr = top
		  top = fc(7,ptr)
		  fc(7,ptr) = -1
	       go to 90
	    else
	       fc(7,ptr) = topnv
	       topnv = ptr
	    endif
	 go to 80
  110    continue
	 if (front .eq. 0) then
	    ierr = 306
	    return
	 endif
c
c        find remaining vis boundary faces, add new tetra with vertex i.
c
	 call vbfac(vcl(1,vi),ctr,vcl,vm,bf,fc,front,topnv)
         if (ierr .ne. 0) return
	 call nwthou(i,npt,sizht,nbf,nfc,maxbf,maxfc,bf,fc,ht,ntetra,
     $      hdavbf,hdavfc,front,back,bfi)
         if (ierr .ne. 0) return
  120 continue
c
      nface = nfc
      ptr = hdavfc
  130 continue
      if (ptr .eq. 0) go to 140
	 nface = nface - 1
	 ptr = -fc(1,ptr)
	 go to 130
  140 continue
      fc(7,1) = hdavbf
      fc(7,2) = hdavfc
c
  600 format (/1x,'itris3: first tetrahedron: ',4i7/4x,'i3, i4 =',2i7)
  610 format (/1x,'step',i7,':   vertex i =',i7)
      end
