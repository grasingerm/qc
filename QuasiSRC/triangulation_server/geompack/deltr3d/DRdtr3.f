      program drdtr3
      implicit logical (a-z)
c
c     driver program for testing routines dtris3,dtriw3,imptr3,itris3,
c        bnsrt3,vbfac,walkt3,nwthou,nwthin,nwthfc,nwthed,frstet,swapes,
c        swapmu,tetmu,tetlst,imptrf,swaptf,fndmsw,imptrd,fnmswd
c     routines called:
c        bnsrt3, ccsph, dtris3, dtriw3, emnrth, gtime, imptr3, initcb,
c        itris3, prime, radrth, randpt, sangmn, tetlst
c
      integer ierr,iprt,msglvl
      double precision pi,tol
      common /gerror/ ierr
      common /gconst/ pi,tol
      common /gprint/ iprt,msglvl
c
      integer maxalg,maxbf,maxfc,maxht,maxth,maxvc
      parameter (maxalg = 10)
      parameter (maxbf = 100000)
      parameter (maxfc = 800000)
      parameter (maxht = 8011)
      parameter (maxth = 400000)
      parameter (maxvc = 100000)
c
      integer bf(3,maxbf),crit(maxalg),fc(7,maxfc),ht(0:maxht-1)
      integer tetra(4,maxth),vm(maxvc)
      integer a,alg,axis,b,c,d,e,i,ialg,imeas,in,inmode,irdr,j,mlvl,nalg
      integer nav,nbf,nface,nfc,nlo,npt,nt,ntetra,nx,nz,prime,seed,sizht
      real t0,tritim
      double precision sang(4),scal(3),tran(3),vcl(3,maxvc)
      double precision binexp,eta,mineta,minrho,minsig,radsq,rho,sigma
      double precision sx,sz,tolin,emnrth,radrth,sangmn
      logical bndcon
c
c     inmode = 1: read vertices; inmode = 2: generate random vertices;
c        inmode = 3: worst case set of vertices.
c     alg = 1: dtris3; alg = 2: dtriw3; alg = 3: dtriw3 with call to
c        bnsrt3 first; alg = 4: itris3.
c     nalg = 1: construct delaunay triang; 2 <= abs(nalg) <= maxalg:
c        improve based on local max-min empty circumsphere, solid angle,
c        radius ratio, mean ratio criterion as specified by crit(i) = 0,
c        1, 2, 3, resp, 2 <= i <= abs(nalg), i.e. there are abs(nalg)-1
c        improvements; nalg < -1: improvements are boundary-constrained.
c     msglvl = -1: print only measurements and allow multiple runs;
c        msglvl = 0: print arrays; msglvl = 2: print list of tetrahedra
c        in output units 1, 2, 3; msglvl = 4: print new tetrahedra,
c        local optimality tests, local transformations.
c
      irdr = 5
      imeas = 7
   10 continue
      read (irdr,*) inmode,alg,nalg,mlvl,npt,tolin,binexp
      if (inmode .le. 0 .or. alg .le. 0 .or. alg .gt. 4) stop
      call initcb(tolin)
      if (npt .gt. maxvc) then
	 write (iprt,600) 'maxvc',npt
	 stop
      endif
      msglvl = mlvl
      write (imeas,700) inmode,alg,nalg,npt,tol,binexp
      bndcon = (nalg .lt. 0)
      nalg = min(abs(nalg),maxalg)
      if (nalg .gt. 1) read (irdr,*) (crit(i),i=2,nalg)
      crit(1) = 0
      if (alg .eq. 4) crit(1) = -1
      if (inmode .eq. 1) then
	 read (irdr,*) (vcl(1,i),vcl(2,i),vcl(3,i),i=1,npt)
	 sizht = prime(npt*3/2)
      else if (inmode .eq. 2) then
	 read (irdr,*) seed,axis,nav,(scal(i),i=1,3),(tran(i),i=1,3)
         write (imeas,710) seed,axis,nav,(scal(i),i=1,3),(tran(i),i=1,3)
	 call randpt(3,npt,seed,axis,nav,scal,tran,3,vcl)
	 sizht = prime(npt*3/2)
      else 
	 nx = npt/2
	 nz = npt - nx
	 sx = 2.0d0*pi/nx
	 sz = 1.0d0/(nz - 1)
	 do 20 i = 1,nx
	    vcl(1,i) = cos(i*sx)
	    vcl(2,i) = sin(i*sx)
	    vcl(3,i) = 0.0d0
   20    continue
	 do 30 i = 1,nz
	    vcl(1,nx+i) = 0.0d0
	    vcl(2,nx+i) = 0.0d0
	    vcl(3,nx+i) = (i - 1)*sz
   30    continue
	 sizht = prime(npt**2/20)
      endif
      sizht = min(sizht,maxht)
      do 40 i = 1,npt
	 vm(i) = i
   40 continue
      if (msglvl .ge. 0) write (iprt,650) npt,(i,vcl(1,i),vcl(2,i),
     $   vcl(3,i),i=1,npt)
c
      do 70 ialg = 1,nalg
         call gtime(t0)
	 if (ialg .eq. 1) then
            if (alg .eq. 1) then
	       call dtris3(npt,sizht,maxbf,maxfc,vcl,vm,nbf,nfc,nface,
     $            ntetra,bf,fc,ht)
            else if (alg .le. 3) then
	       if (alg .eq. 3) call bnsrt3(binexp,npt,vcl,vm,fc,ht)
	       call dtriw3(npt,sizht,maxbf,maxfc,vcl,vm,nbf,nfc,nface,
     $            ntetra,bf,fc,ht)
	    else
	       call itris3(npt,sizht,maxbf,maxfc,vcl,vm,nbf,nfc,nface,
     $            ntetra,bf,fc,ht)
	    endif
	 else
	    call imptr3(bndcon,.true.,crit(ialg),npt,sizht,maxfc,vcl,vm,
     $         nfc,ntetra,bf,fc,ht,nface)
	 endif
         call gtime(tritim)
         tritim = tritim - t0
         if (ierr .ne. 0) then
	    write (iprt,610) ierr
	    write (imeas,610) ierr
	    stop
         endif
         if (ntetra .gt. maxth) then
	    write (iprt,600) 'maxth',ntetra
	    stop
         endif
	 call tetlst(nfc,vm,fc,nt,tetra)
	 if (nt .ne. ntetra) then
	    write (iprt,620) nt,ntetra
	    write (imeas,620) nt,ntetra
	    stop
	 endif
c
	 nlo = 0
	 do 50 i = 1,nfc
	    if (fc(1,i) .gt. 0 .and. fc(5,i) .gt. 0) then
	       a = vm(fc(1,i))
	       b = vm(fc(2,i))
	       c = vm(fc(3,i))
	       d = vm(fc(4,i))
	       e = vm(fc(5,i))
	       call ccsph(.true.,vcl(1,a),vcl(1,b),vcl(1,c),vcl(1,d),
     $            vcl(1,e),sang,radsq,in)
	       if (in .ge. 1) then
		  nlo = nlo + 1
		  if (msglvl .eq. 4) then
		     if (nlo .eq. 1) write (iprt,630)
		     write (iprt,640) i,a,b,c,d,e
		  endif
	       endif
	    endif
   50    continue
	 minsig = 2.0d0
	 minrho = 2.0d0
	 mineta = 2.0d0
	 do 60 i = 1,nt
	    a = tetra(1,i)
	    b = tetra(2,i)
	    c = tetra(3,i)
	    d = tetra(4,i)
	    sigma = sangmn(vcl(1,a),vcl(1,b),vcl(1,c),vcl(1,d),sang)
	    rho = radrth(vcl(1,a),vcl(1,b),vcl(1,c),vcl(1,d))
	    eta = emnrth(vcl(1,a),vcl(1,b),vcl(1,c),vcl(1,d))
	    minsig = min(minsig,sigma)
	    minrho = min(minrho,rho)
	    mineta = min(mineta,eta)
   60    continue
	 minsig = minsig*1.5d0*sqrt(6.0d0)
c
         write (imeas,720) crit(ialg),nbf,nfc,nface,ntetra,sizht,nlo,
     $      minsig,minrho,mineta,tritim
	 if (msglvl .eq. 2 .and. ialg .le. 4)
     $      write (ialg,660) ((tetra(j,i),j=1,4),i=1,nt)
         if (msglvl .ge. 0) then
	    if (ialg .eq. 1) write (iprt,670) 'vm',npt,(vm(i),i=1,npt)
	    write (iprt,670) 'ht',sizht,(ht(i),i=0,sizht-1)
	    write (iprt,680) nfc,(i,(fc(j,i),j=1,7),i=1,nfc)
	    if (ialg .eq. 1 .or. .not. bndcon)
     $         write (iprt,690) nbf,(i,(bf(j,i),j=1,3),i=1,nbf)
	 endif
         write(12, 730)
         write(12, 740)
         write(12, 750) npt, ntetra
         write(12, 760) (vcl(1,i),vcl(2,i),vcl(3,i),i=1,npt)
         write(12, *)
         write(12, 660) ((tetra(j,i),j=1,4),i=1,nt)
   70 continue
      if (msglvl .lt. 0) go to 10
c
  600 format (1x,'*** ',a,' must be increased to',i8)
  610 format (/1x,'ierr=',i5)
  620 format (/1x,'error from tetlst, nt != ntetra :',2i7)
  630 format (1x,'nonlocally optimal faces - sphere criterion')
  640 format (1x,'face',i7,': a,b,c,d,e=',5i7)
  650 format (1x,'vcl   ',i7/(1x,i7,3f23.15))
  660 format (1x,4i7)
  670 format (/1x,a,3x,i7/(1x,10i7))
  680 format (/1x,'fc   ',i7/(1x,8i7))
  690 format (/1x,'bf   ',i7/(1x,4i7))
  700 format (/1x,'inm=',i1,'   alg=',i1,'   nalg=',i2,'   npt=',i7,
     $   '   tol=',d14.7,'   binexp=',f9.7)
  710 format (1x,'seed=',i9,'   axis=',i2,'   nav=',i7/1x,
     $   'scal(1:3),tran(1:3)=',6f9.4)
  720 format (1x,'crit=',i2,'   nbf=',i5,'   nfc=',i6,'   nface=',i6,
     $   '   ntetra=',i6,'   sizht=',i5/1x,'nlo=',i5,'   msig=',f9.7,
     $   '   mrho=',f9.7,'   meta=',f9.7,'   tim=',f9.3)
 730  format (1x, 'title = "fem"')
 740  format (1x, 'variables = "x" "y" "z"')
 750  format (1x, 'zone n= ',i6,', e= ',i6,
     $   ', f=fepoint, et=tetrahedron')
 760  format (1x, 3f23.15 )
      end
