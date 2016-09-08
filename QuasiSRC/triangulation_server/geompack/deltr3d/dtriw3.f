      subroutine dtriw3(npt,sizht,maxbf,maxfc,vcl,vm,nbf,nfc,nface,
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
c     purpose: construct delaunay triangulation of 3-d vertices using
c        incremental approach and local transformations. vertices are
c        inserted one at a time in order given by vm array. the initial
c        tetrahedra created due to a new vertex are obtained by a walk
c        through the triangulation until location of vertex is known.
c
c     input parameters:
c	 npt - number of 3-d vertices (points)
c	 sizht - size of hash table ht; a good choice is a prime number
c              which is about 1/8 * nface (or 3/2 * npt for random
c              points from the uniform distribution)
c        maxbf - maximum size available for bf array
c        maxfc - maximum size available for fc array
c        vcl(1:3,1:*) - vertex coordinate list
c        vm(1:npt) - indices of vertices of vcl being triangulated;
c              vertices are inserted in order given by vm, except:
c
c     updated parameters:
c        vm(1:npt) - the third and fourth elements may be swapped so
c              that first 4 vertices are non-coplanar
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
c        ierr is set to 11, 12, 300, 301, 302, 303, 304, 305,307, or 309
c
c     routines called:
c        frstet,htins,nwthed,nwthfc,nwthin,nwthou,swapes,vbfac,walkt3
c
      integer ierr,iprt,msglvl
      common /gerror/ ierr
      common /gprint/ iprt,msglvl
      save /gerror/,/gprint/
c
      integer back,front,hdavbf,hdavfc,i,i3,i4,ifac,ind,ivrt,ptr,vi
      double precision ctr(3)
c
c     find initial valid tetrahedron, and initialize data structures.
c
      call frstet(.false.,npt,vcl,vm,i3,i4)
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
      ifac = 4
      if (msglvl .eq. 4) write (iprt,600) (vm(i),i=1,4),i3,i4
c
c     insert ith vertex into delaunay triang of first i-1 vertices.
c     walk through triang to find location of vertex i, create new
c     tetrahedra, apply local transf based on empty sphere criterion.
c
      do 30 i = 5,npt
	 vi = vm(i)
	 if (msglvl .eq. 4) write (iprt,610) i,vi
	 if (fc(5,ifac) .eq. i-1) then
	    ivrt = 5
	 else
	    ivrt = 4
	 endif
	 call walkt3(vcl(1,vi),npt,sizht,ntetra,vcl,vm,fc,ht,ifac,ivrt)
	 if (ierr .ne. 0) return
	 if (ivrt .eq. 6) then
	    call nwthfc(i,ifac,npt,sizht,nbf,nfc,maxbf,maxfc,bf,fc,ht,
     $         ntetra,hdavbf,hdavfc,front,back)
	 else if (ivrt .ge. 4) then
	    call nwthin(i,ifac,ivrt,npt,sizht,nfc,maxfc,fc,ht,ntetra,
     $         hdavfc,front,back)
	 else if (ivrt .eq. 0) then
	    front = ifac
	    call vbfac(vcl(1,vi),ctr,vcl,vm,bf,fc,front,0)
            if (ierr .ne. 0) return
	    call nwthou(i,npt,sizht,nbf,nfc,maxbf,maxfc,bf,fc,ht,ntetra,
     $         hdavbf,hdavfc,front,back,ind)
	 else if (ivrt .ge. 1) then
	    call nwthed(i,ifac,ivrt,npt,sizht,nbf,nfc,maxbf,maxfc,bf,fc,
     $         ht,ntetra,hdavbf,hdavfc,front,back)
	 else
	    ierr = 302
	 endif
         if (ierr .ne. 0) return
	 call swapes(.false.,i,npt,sizht,nfc,maxfc,vcl,vm,bf,fc,ht,
     $      ntetra,hdavfc,front,back,ind)
         if (ierr .ne. 0) return
	 if (ind .ne. 0) ifac = ind
   30 continue
c
      nface = nfc
      ptr = hdavfc
   40 continue
      if (ptr .eq. 0) go to 50
	 nface = nface - 1
	 ptr = -fc(1,ptr)
	 go to 40
   50 continue
      fc(7,1) = hdavbf
      fc(7,2) = hdavfc
c
  600 format (/1x,'dtriw3: first tetrahedron: ',4i7/4x,'i3, i4 =',2i7)
  610 format (/1x,'step',i7,':   vertex i =',i7)
      end
