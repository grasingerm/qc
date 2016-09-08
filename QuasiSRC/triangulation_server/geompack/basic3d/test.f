      program test
      implicit none

      external initcb
      external insph

      double precision a(3), b(3), c(3), d(3)
      double precision center(3)
      double precision radius

      integer i
      
      a(1)=-2.0
      a(2)=1.0
      a(3)=-3.0

      b(1)=10.0
      b(2)=2.0
      b(3)=1.0

      c(1)=2.0
      c(2)=8.0
      c(3)=-5.0

      d(1)=-3.0
      d(2)=-2.0
      d(3)=10.0
      
      call initcb(0.0d0)

      call insph(a,b,c,d,center,radius)

      write(6,*) 'center=',(center(i),i=1,3)
      write(6,*) 'radius=', radius

      end
      
