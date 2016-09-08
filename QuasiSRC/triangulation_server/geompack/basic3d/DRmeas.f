      program drmeas
      implicit logical (a-z)
c
c     driver program for studying relationship between tetrahedron
c     shape measures.
c     routines called:
c        emnrth, radrth, sangmn
c
      integer iprt,irdr
      double precision a(3),b(3),c(3),d(3),sang(4)
      double precision emnrth,eta,radrth,rho,sangmn,sigma
      double precision l1,l2,l3,u1,u2,u3,q1,q2,q3,r1,r2,r3
      data a/0.0d0,0.0d0,0.0d0/
      data b/1.0d0,0.0d0,0.0d0/
      data c(3)/0.0d0/
      data l1,l2,l3,u1,u2,u3/3*1.0d0,3*0.0d0/
c
      irdr = 5
      iprt = 6
   10 continue
      write (iprt,*) 'enter s,t,x,y,z:'
      read (irdr,*) c(1),c(2),d(1),d(2),d(3)
      if (c(2) .le. 0.0d0 .or. d(3) .le. 0.0d0) go to 20
      sigma = sangmn(a,b,c,d,sang)
      rho = radrth(a,b,c,d)
      eta = emnrth(a,b,c,d)
      q1 = rho/eta**3
      r1 = rho/eta**0.75d0
      q2 = sigma/(eta*sqrt(eta))
      r2 = sigma/eta**0.75d0
      q3 = sigma/rho**2
      r3 = sigma/sqrt(rho)
      l1 = min(l1,q1)
      l2 = min(l2,q2)
      l3 = min(l3,q3)
      u1 = max(u1,r1)
      u2 = max(u2,r2)
      u3 = max(u3,r3)
      write (iprt,600) sigma,rho,eta,q1,q2,q3,r1,r2,r3
      go to 10
   20 continue
      write (iprt,*) 'lower and upper bounds'
      write (iprt,600) l1,l2,l3,u1,u2,u3
  600 format (3e15.7)
      end
