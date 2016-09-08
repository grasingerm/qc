      subroutine vbfac(pt,ctr,vcl,vm,bf,fc,topv,topnv)
      implicit logical (a-z)
      integer bf(3,*),fc(7,*),topnv,topv,vm(*)
      double precision ctr(3),pt(3),vcl(3,*)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: determine boundary faces of 3-d triangulation visible
c        from point pt, given a starting visible boundary face.
c
c     input parameters:
c	 pt(1:3) - 3-d point
c        ctr(1:3) - 3-d point in interior of triangulation
c        vcl(1:3,1:*) - vertex coordinate list
c        vm(1:*) - indices of vertices of vcl being triangulated
c        bf(1:3,1:*) -  array of boundary face records; see dtris3
c        fc(1:7,1:*) - array of face records; see routine dtris3; row 7
c              is used for links of 3 stacks in this routine
c        topv - index of fc of visible boundary face
c        topnv - index of top of stack of boundary faces already found
c              to be not visible from pt, or 0 for empty stack
c
c     updated parameters:
c        fc(7,*) - gets updated; only stack of visible boundary faces
c              remains at end of routine
c        topv - index of top of stack of visible boundary faces
c
c     abnormal return:
c        ierr is set to 301
c
c     routines called:
c        opside
c
      integer ierr
      integer j,k,nbr,op,opside,ptr,topn,topt,va,vb,vc
c
c     topn is index of top of stack of non-visible boundary faces.
c     topt is index of top of stack of boundary faces to be tested.
c
      topn = topnv
      topt = 0
      fc(7,topv) = 0
      k = -fc(5,topv)
      do 10 j = 1,3
	 nbr = bf(j,k)
	 if (fc(7,nbr) .eq. -1) then
	    fc(7,nbr) = topt
	    topt = nbr
	 endif
   10 continue
c
   20 continue
      if (topt .eq. 0) go to 40
	 ptr = topt
	 topt = fc(7,ptr)
	 va = vm(fc(1,ptr))
	 vb = vm(fc(2,ptr))
	 vc = vm(fc(3,ptr))
	 op = opside(vcl(1,va),vcl(1,vb),vcl(1,vc),ctr,pt)
	 if (op .eq. 2) then
	    ierr = 301
	    return
	 endif
	 if (op .eq. 1) then
	    fc(7,ptr) = topv
	    topv = ptr
	    k = -fc(5,ptr)
	    do 30 j = 1,3
	       nbr = bf(j,k)
	       if (fc(7,nbr) .eq. -1) then
		  fc(7,nbr) = topt
		  topt = nbr
	       endif
   30       continue
	 else
	    fc(7,ptr) = topn
	    topn = ptr
	 endif
      go to 20
c
c     for boundary faces not visible from pt, set fc(7,*) = -1.
c
   40 continue
      if (topn .eq. 0) return
	 ptr = topn
	 topn = fc(7,ptr)
	 fc(7,ptr) = -1
      go to 40
      end
