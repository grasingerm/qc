c
c this dummy FORTRAN function extracts error and feeds it back to C
c

      subroutine geterr( auxerr )
      integer auxerr

      common /GERROR/ ierr

      auxerr=ierr

      end

      subroutine setout( lvl )
      integer lvl
      integer IPRT,MSGLVL

      common /GPRINT/ IPRT,MSGLVL
      MSGLVL=lvl

      end
