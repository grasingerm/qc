      integer function htfun(a, b, c, n, p)
      implicit none
      integer a, b, c, n, p
c
c hash function that prevents overflows
c
      integer*8 tmp_k
      integer   k
      tmp_k = kmod(kint(real(a))*kint(real(n)) + kint(real(b)), 
     1        kint(real(p)))
      k     = int(kmod(tmp_k*kint(real(n)) + kint(real(c)), 
     1        kint(real(p))))
      htfun = k
      end
