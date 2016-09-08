      program drshap
      implicit logical (a-z)
c
c     driver program for testing routines radrth, emnrth, sangmn.
c     routines called:
c        ccsph, emnrth, initcb, insph, radrth, sangmn, sdang
c
      double precision pi,tol
      common /gconst/ pi,tol
c
      integer in,iprt,irdr
      double precision a(3),b(3),c(3),centre(3),d(3),dang(6),sang(4)
      double precision e,emnrth,eta,rad,radrth,radsq,rho1,rho2,samin1
      double precision samin2,sangmn,xc,xc0,xc1,xcinc,xd,xd0,xd1,xdinc
      double precision yc,yc0,yc1,ycinc,yd,yd0,yd1,ydinc,zd,zd0,zd1
      double precision zdinc
      data a/0.0d0,0.0d0,0.0d0/
      data b/1.0d0,0.0d0,0.0d0/
      data c(3)/0.0d0/
c
      irdr = 5
      iprt = 6
      call initcb(0.0d0)
      read (irdr,*) xc0,xc1,xcinc,yc0,yc1,ycinc,xd0,xd1,xdinc,
     $   yd0,yd1,ydinc,zd0,zd1,zdinc
      write (iprt,600)
      do 50 xc = xc0,xc1,xcinc
	 c(1) = xc
         do 40 yc = yc0,yc1,ycinc
	    c(2) = yc
            do 30 xd = xd0,xd1,xdinc
	       d(1) = xd
               do 20 yd = yd0,yd1,ydinc
		  d(2) = yd
                  do 10 zd = zd0,zd1,zdinc
		     d(3) = zd
                     call sdang(a,b,c,d,sang,dang)
                     samin1 = min(sang(1),sang(2),sang(3),sang(4)) - pi
                     samin2 = 2.0d0*asin(sangmn(a,b,c,d,sang))
 	             call insph(a,b,c,d,centre,rad)
 	             call ccsph(.false.,a,b,c,d,e,centre,radsq,in)
	             rho1 = 3.0d0*rad/sqrt(abs(radsq))
	             rho2 = radrth(a,b,c,d)
	             eta = emnrth(a,b,c,d)
 	             write (iprt,610) c(1),c(2),d(1),d(2),d(3),
     $                  samin2,rho2,eta
		     if (abs(samin1 - samin2) .gt. 100.0d0*tol)
     $                  write (iprt,620) 'samin:',samin1,samin2 
		     if (abs(rho1 - rho2) .gt. 100.0d0*tol)
     $                  write (iprt,620) 'rho:',rho1,rho2 
   10             continue
   20          continue
   30       continue
   40    continue
   50 continue
  600 format (3x,'xc',6x,'yc',6x,'xd',6x,'yd',6x,'zd',7x,'samin',8x,
     $   'rho',9x,'eta')
  610 format (5f8.3,3f12.6)
  620 format (10x,a,2f20.14)
      end
