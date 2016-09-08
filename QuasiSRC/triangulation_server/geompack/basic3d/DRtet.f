      program drtet
      implicit logical (a-z)
c
c     driver program for testing routines ccsph, baryth, opside, volth,
c        angle3, sdang, insph, radrth, emnrth, sangmn, ccradi.
c     routines called:
c        baryth, ccradi, ccsph, emnrth, initcb, insph, opside, radrth,
c        sangmn, sdang, volth
c
      double precision pi,tol
      common /gconst/ pi,tol
c
      integer case,i,in,iprt,irdr,op,opside
      double precision a(3),alpha(4),b(3),c(3),centre(3),d(3),dang(6)
      double precision e(3),sang(4),rad,radsq,ratio,vol,volth
      double precision ccradi,emnrth,radrth,sangmn
      logical degen
c
      irdr = 5
      iprt = 6
      call initcb(0.0d0)
   10 continue
	 write (iprt,*) 'enter case:'
	 read (irdr,*) case
	 if (case .lt. 1 .or. case .gt. 7) stop
	 if (case .eq. 1) then
 	    write (iprt,*) 'ccsph - enter a,b,c,d,e:'
 	    read (irdr,*) a(1),a(2),a(3),b(1),b(2),b(3),c(1),c(2),c(3),
     $         d(1),d(2),d(3),e(1),e(2),e(3)
 	    call ccsph(.true.,a,b,c,d,e,centre,radsq,in)
 	    write (iprt,600) 'centre=',(centre(i),i=1,3)
	    if (radsq .ge. 0.0d0) radsq = sqrt(radsq)
 	    write (iprt,*) 'rad=',radsq,'   in=',in
	 else if (case .eq. 2) then
 	    write (iprt,*) 'baryth - enter a,b,c,d,e:'
 	    read (irdr,*) a(1),a(2),a(3),b(1),b(2),b(3),c(1),c(2),c(3),
     $         d(1),d(2),d(3),e(1),e(2),e(3)
 	    call baryth(a,b,c,d,e,alpha,degen)
 	    write (iprt,600) 'alpha=',(alpha(i),i=1,4)
 	    write (iprt,*) 'degen=',degen
	 else if (case .eq. 3) then
 	    write (iprt,*) 'opside - enter a,b,c,d,e:'
 	    read (irdr,*) a(1),a(2),a(3),b(1),b(2),b(3),c(1),c(2),c(3),
     $         d(1),d(2),d(3),e(1),e(2),e(3)
 	    op = opside(a,b,c,d,e)
 	    write (iprt,*) 'opside=',op
	 else if (case .eq. 4) then
 	    write (iprt,*) 'volth - enter a,b,c,d:'
 	    read (irdr,*) a(1),a(2),a(3),b(1),b(2),b(3),c(1),c(2),c(3),
     $         d(1),d(2),d(3)
 	    vol = volth(a,b,c,d)/6.0d0
 	    write (iprt,*) 'volth=',vol
	 else if (case .eq. 5) then
 	    write (iprt,*) 'sdang,sangmn - enter a,b,c,d:'
 	    read (irdr,*) a(1),a(2),a(3),b(1),b(2),b(3),c(1),c(2),c(3),
     $         d(1),d(2),d(3)
 	    call sdang(a,b,c,d,sang,dang)
 	    write (iprt,600) 'sang=',(sang(i) - pi,i=1,4)
 	    write (iprt,600) 'dang=',(dang(i),i=1,6)
 	    write (iprt,*) 'sin(sangmn/2)=',sangmn(a,b,c,d,sang)
 	    write (iprt,600) 'sang=',(2.0d0*asin(sang(i)),i=1,4)
	 else if (case .eq. 6) then
 	    write (iprt,*) 'insph,ccradi,radrth - enter a,b,c,d:'
 	    read (irdr,*) a(1),a(2),a(3),b(1),b(2),b(3),c(1),c(2),c(3),
     $         d(1),d(2),d(3)
 	    call insph(a,b,c,d,centre,rad)
 	    write (iprt,600) 'incent=',(centre(i),i=1,3)
 	    radsq = ccradi(a,b,c,d)
	    if (radsq .eq. 0.0d0) then
	       ratio = 0.0d0
	    else
	       ratio = rad*sqrt(radsq)*12.0d0
	    endif
 	    write (iprt,*) 'inrad=',rad,'   ratio=',ratio
	    write (iprt,*) 'radrth=',radrth(a,b,c,d)
	 else
 	    write (iprt,*) 'emnrth - enter a,b,c,d:'
 	    read (irdr,*) a(1),a(2),a(3),b(1),b(2),b(3),c(1),c(2),c(3),
     $         d(1),d(2),d(3)
	    write (iprt,*) 'emnrth=',emnrth(a,b,c,d)
	 endif
      go to 10
  600 format (a,6f12.6)
      end
