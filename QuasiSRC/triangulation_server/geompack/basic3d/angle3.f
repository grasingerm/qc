      double precision function angle3(u,v,rtolsq)
      implicit logical (a-z)
      double precision rtolsq,u(3),v(3)
c
c     written and copyright by:
c        barry joe, dept. of computing science, univ. of alberta
c        edmonton, alberta, canada  t6g 2h1
c        phone: (403) 492-5757      email: barry@cs.ualberta.ca
c
c     purpose: compute angle in range [0,pi] between 3-d vectors u, v.
c
c     input parameters:
c        u(1:3),v(1:3) - vectors of length 3
c        rtolsq - relative tolerance used to detect 0 vector based on
c              square of euclidean length
c
c     returned function value:
c        angle3 - angle between 2 vectors in range [0,pi]
c        [note: if u or v is the 0 vector, angle3 = pi is returned.]
c
      double precision pi,tol
      common /gconst/ pi,tol
      save /gconst/
c
      double precision dotp,lu,lv,t
c
      dotp = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
      lu = u(1)**2 + u(2)**2 + u(3)**2
      lv = v(1)**2 + v(2)**2 + v(3)**2
      if (lu .gt. rtolsq .and. lv .gt. rtolsq) then
	 t = dotp/sqrt(lu*lv)
	 if (abs(t) .gt. 1.0d0 - tol) t = sign(1.0d0,t)
         angle3 = acos(t)
      else
	 angle3 = pi
      endif
      end
