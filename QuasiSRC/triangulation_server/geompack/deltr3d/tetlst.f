      subroutine tetlst(nfc,vm,fc,nt,tetra)
      implicit logical (a-z)
      integer nfc,nt
      integer fc(7,nfc),tetra(4,*),vm(*)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: construct list of tetrahedra from fc array. global vertex
c        indices from vm are produced. the vertex indices for each
c        tetrahedron are sorted in increasing order.
c
c     input parameters:
c        nfc - number of positions used in fc array
c        vm(1:*) - indices of vertices of vcl that are triangulated
c        fc(1:7,1:nfc) - array of face records; see routine dtris3
c
c     output parameters:
c        nt - number of tetrahedra
c        tetra(1:4,1:nt) - contains global tetrahedron indices; it is
c              assumed there is enough space
c
      integer a,i,j,k,l,t(4)
c
      nt = 0
      do 20 i = 1,nfc
	 if (fc(1,i) .le. 0) go to 20
	 do 10 k = 4,5
	    if (fc(k,i) .gt. fc(3,i)) then
	       nt = nt + 1
               tetra(1,nt) = fc(1,i)
               tetra(2,nt) = fc(2,i)
               tetra(3,nt) = fc(3,i)
               tetra(4,nt) = fc(k,i)
	    endif
   10    continue
   20 continue
c
      do 70 k = 1,nt
	 do 30 i = 1,4
	    t(i) = vm(tetra(i,k))
   30    continue
	 do 50 i = 1,3
	    l = i
	    do 40 j = i+1,4
	       if (t(j) .lt. t(l)) l = j
   40       continue
	    a = t(i)
	    t(i) = t(l)
	    t(l) = a
   50    continue
	 do 60 i = 1,4
	    tetra(i,k) = t(i)
   60    continue
   70 continue
      end
