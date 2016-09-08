      subroutine ifacty(ind,npt,sizht,vcl,vm,fc,ht,type,a,b,c)
      implicit logical (a-z)
      integer a,b,c,ind,npt,sizht
      integer fc(7,*),ht(0:sizht-1),vm(npt)
      double precision vcl(3,*)
      character type*3
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: determine type of interior face of 3-d triangulation.
c
c     input parameters:
c        ind - index of fc array, assumed to be an interior face
c        npt - size of vm array
c        sizht - size of ht array
c        vcl(1:3,1:*) - vertex coordinate list
c        vm(1:npt) - vertex mapping array, from local to global indices
c        fc(1:7,1:*) - array of face records; see routine dtris3
c        ht(0:sizht-1) - hash table using direct chaining
c
c     output parameters:
c        type - 't23','t32','t22','t44','n32','n44','n40','n30',or 'n20'
c        a,b,c - local indices of interior face; for t32, n32, t22, t44,
c              or n44 face, ab is edge that would get swapped out and c
c              is third vertex; for t23 face, a < b < c; for n40 face,
c              a is interior vertex; for n30 face, a is inside a face
c              with vertex b; for n20 face, a is on an edge
c
c     abnormal return:
c        ierr is set to 300, 301, or 309
c
c     routines called:
c        baryth,htsrc
c
      integer ierr
      common /gerror/ ierr
      save /gerror/
c
      integer d,e,f,htsrc,ind1,j,kneg,kzero
      double precision alpha(4)
      logical degen
c
      a = fc(1,ind)
      b = fc(2,ind)
      c = fc(3,ind)
      d = fc(4,ind)
      e = fc(5,ind)
      call baryth(vcl(1,vm(a)),vcl(1,vm(b)),vcl(1,vm(c)),vcl(1,vm(d)),
     $   vcl(1,vm(e)),alpha,degen)
      if (degen) then
         ierr = 301
         return
      else if (alpha(4) .ge. 0.0d0) then
	 ierr = 309
         write(*,*) 'error 309 in ifacty'
         write(*,*) vcl(1,vm(a)), vcl(2,vm(a)), vcl(3,vm(a))
         write(*,*) vcl(1,vm(b)), vcl(2,vm(b)), vcl(3,vm(b))
         write(*,*) vcl(1,vm(c)), vcl(2,vm(c)), vcl(3,vm(c))
         write(*,*) vcl(1,vm(d)), vcl(2,vm(d)), vcl(3,vm(d))
         write(*,*) vcl(1,vm(e)), vcl(2,vm(e)), vcl(3,vm(e))
         return
      endif
      kneg = 1
      kzero = 0
      do 10 j = 1,3
         if (alpha(j) .lt. 0.0d0) then
            kneg = kneg + 1
         else if (alpha(j) .eq. 0.0d0) then
            kzero = kzero + 1
         endif
   10 continue
      type = 'xxx'
      if (kneg .eq. 1 .and. kzero .eq. 0) then
         type = 't23'
      else if (kneg .eq. 2 .and. kzero .eq. 0) then
         if (alpha(1) .lt. 0.0d0) then
            j = a
            a = c
            c = j
         else if (alpha(2) .lt. 0.0d0) then
            j = b
            b = c
            c = j
         endif
         ind1 = htsrc(a,b,d,npt,sizht,fc,ht)
         if (ind1 .le. 0) then
            ierr = 300
            return
         else if (fc(4,ind1) .eq. e .or. fc(5,ind1) .eq. e) then
            type = 't32'
         else
            type = 'n32'
         endif
      else if (kneg .eq. 3 .and. kzero .eq. 0) then
         type = 'n40'
         if (alpha(2) .gt. 0.0d0) then
            j = a
            a = b
            b = j
         else if (alpha(3) .gt. 0.0d0) then
            j = a
            a = c
            c = j
         endif
      else if (kneg .eq. 1 .and. kzero .eq. 1) then
         if (alpha(1) .eq. 0.0d0) then
            j = a
            a = c
            c = j
         else if (alpha(2) .eq. 0.0d0) then
            j = b
            b = c
            c = j
         endif
         ind1 = htsrc(a,b,d,npt,sizht,fc,ht)
         if (ind1 .le. 0) then
            ierr = 300
            return
         else if (fc(5,ind1) .le. 0) then
            type = 't22'
         else
            if (fc(4,ind1) .eq. c) then
               f = fc(5,ind1)
            else
               f = fc(4,ind1)
            endif
            ind1 = htsrc(a,b,e,npt,sizht,fc,ht)
            if (ind1 .le. 0) then
               ierr = 300
               return
            else if (fc(4,ind1) .eq. f .or. fc(5,ind1) .eq. f) then
               type = 't44'
            else
               type = 'n44'
            endif
         endif
      else if (kneg .eq. 2 .and. kzero .eq. 1) then
         type = 'n30'
         if (alpha(1) .eq. 0.0d0) then
            j = a
            a = c
            c = j
         else if (alpha(2) .eq. 0.0d0) then
            j = b
            b = c
            c = j
	    alpha(2) = alpha(3)
         endif
         if (alpha(2) .gt. 0.0d0) then
            j = a
            a = b
            b = j
         endif
      else if (kneg .eq. 1 .and. kzero .eq. 2) then
         type = 'n20'
         if (alpha(2) .gt. 0.0d0) then
            j = a
            a = b
            b = j
         else if (alpha(3) .gt. 0.0d0) then
            j = a
            a = c
            c = j
         endif
      endif
      end
