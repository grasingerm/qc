      subroutine availf(hdavfc,nfc,maxfc,fc,ind)
      implicit logical (a-z)
      integer hdavfc,ind,maxfc,nfc
      integer fc(7,*)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: return index of next available record in fc array,
c        either hdavfc or nfc+1.
c
c     input parameters:
c        hdavfc - head pointer of available records in fc
c        nfc - current number of records used in fc
c        maxfc - maximum number of records available in fc
c        fc(1:7,1:*) - array of face records; see routine dtris3
c
c     updated parameters:
c        hdavfc,nfc - 1 of these may be updated
c
c     output parameters:
c        ind - index of available record (if fc not full)
c
c     abnormal return:
c        ierr is set to 11
c
      integer ierr
      common /gerror/ ierr
      save /gerror/
c
      if (hdavfc .ne. 0) then
	 ind = hdavfc
	 hdavfc = -fc(1,hdavfc)
      else
	 if (nfc .ge. maxfc) then
	    ierr = 11
	 else
	    nfc = nfc + 1
	    ind = nfc
	 endif
      endif
      end
