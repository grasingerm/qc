dnl
dnl-----------------------------------------------------------------------------
dnl                                      
dnl              Jaroslaw Knap, Michael Ortiz and Anna Pandolfi
dnl                    California Institute of Technology
dnl                       (C) 2001 All Rights Reserved
dnl                                      
dnl-----------------------------------------------------------------------------
dnl $Id: compiler.m4,v 1.3 2002/12/06 01:49:49 fago Exp $
dnl
dnl $Log: compiler.m4,v $
dnl Revision 1.3  2002/12/06 01:49:49  fago
dnl Significant fixes to automake build especially on AIX and HPUX.
dnl
dnl Revision 1.2  2002/12/04 17:45:46  fago
dnl Ported to automake 1.7.1, autoconf 2.56. perror -> ERROR. Fixed build bugs,
dnl especially on HPUX.
dnl
dnl Revision 1.1  2002/03/08 00:02:47  knap
dnl Updated build procedure to use automake/autoconf.
dnl
dnl Revision 1.2  2001/10/02 22:43:57  knap
dnl Imported sources.
dnl
dnl Revision 1.6  2001/05/14 05:15:40  knap
dnl Ported to HP-UX.
dnl
dnl Revision 1.5  2001/05/14 04:19:06  knap
dnl Ported to Linux. Corrected various errors on AIX.
dnl
dnl Revision 1.4  2001/05/10 20:52:03  knap
dnl Corrected settting of 64-bit compiler mode options.
dnl
dnl

dnl
dnl Compiler definitions
dnl

dnl
dnl reset the order of checks for compilers to native, gcc
dnl

AC_PROG_F77(xlf95_r xlf90_r xlf_r f77 f90 g77 f2c g77)dnl

dnl
dnl arch dependend
dnl

FFLAGS_64_BIT=""

case "$host" in
  *-hp-hpux*)
  if test $debug_defined = "yes"; then
    FFLAGS="-g"
  else
    FFLAGS="-O"
  fi
  if test "$F77" = "f77"; then
    $FLAGS="$FLAGS +ppu"
  fi
  FFLAGS_64_BIT="+DA2.0W"
  ;;

  *-dec-osf*)
  if test $debug_defined = "yes"; then
    FFLAGS="-g"
  else
    FFLAGS="-O"
  fi
  ;;

  *-sun-solaris*)
  if test $debug_defined = "yes"; then
    FFLAGS="-g"
  else
    FFLAGS="-O"
  fi
  ;;

  *-sgi-irix*)
  if test $debug_defined = "yes"; then
    FFLAGS="-g"
  else
    FFLAGS="-O"
  fi
  FFLAGS_64_BIT="-64"
  FORTRAN_LIBS="-lfortran -lffio -lftn"
  ;;

  *-ibm-aix*)
  # if using f77 adjust flags
  if test $debug_defined = "yes"; then
    FFLAGS="-g -qmaxmem=8192"
  else
    FFLAGS="-O -qmaxmem=8192"
  fi
  
  if test "$F77" = "f77" -o "$F77" = "xlf" -o "$F77" = "xlf_r" ; then
    FFLAGS="$FFLAGS -qnosave"
  fi
 
  if test "$F77" = "xlf90_r" -o "$F77" = "xlf95_r" ; then
    FFLAGS="$FFLAGS -qfixed"
  fi
 
  FFLAGS_64_BIT="-q64"
  ;;

  *-linux-*)
  if test $debug_defined = "yes"; then
    FFLAGS="-g"
  else
    FFLAGS="-O"
  fi
  ;;

  *)
  ;;
esac

dnl
dnl add flags for 64-bit mode
dnl
if test "x$enable_64_bit" = "xyes"; then
  FFLAGS="$FFLAGS $FFLAGS_64_BIT"
fi

dnl
dnl
dnl
AC_F77_WRAPPERS
