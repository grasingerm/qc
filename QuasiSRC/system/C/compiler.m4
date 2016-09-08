dnl
dnl-----------------------------------------------------------------------------
dnl                                      
dnl              Jaroslaw Knap, Michael Ortiz and Anna Pandolfi
dnl                    California Institute of Technology
dnl                       (C) 2001 All Rights Reserved
dnl                                      
dnl-----------------------------------------------------------------------------
dnl $Id: compiler.m4,v 1.4 2002/12/06 01:49:48 fago Exp $
dnl
dnl $Log: compiler.m4,v $
dnl Revision 1.4  2002/12/06 01:49:48  fago
dnl Significant fixes to automake build especially on AIX and HPUX.
dnl
dnl Revision 1.3  2002/12/04 17:45:46  fago
dnl Ported to automake 1.7.1, autoconf 2.56. perror -> ERROR. Fixed build bugs,
dnl especially on HPUX.
dnl
dnl Revision 1.2  2002/10/21 17:19:45  fago
dnl Tweaking autotools build.
dnl
dnl Revision 1.1.1.1  2002/04/22 03:33:47  knap
dnl Initial import
dnl
dnl Revision 1.2  2002/04/05 18:28:27  fago
dnl Added posix4 library check and sched_yield() test.
dnl
dnl Revision 1.1  2002/03/08 00:02:44  knap
dnl Updated build procedure to use automake/autoconf.
dnl
dnl Revision 1.2  2001/10/02 22:43:56  knap
dnl Imported sources.
dnl
dnl Revision 1.2  2001/05/10 20:51:06  knap
dnl Corrected setting of 64-bit compiler options.
dnl
dnl
dnl AC_REQUIRE([AC_CANONICAL_HOST])

dnl
dnl Compiler definitions
dnl

dnl
dnl reset the order of checks for compilers to native, gcc
dnl
AC_PROG_CC(xlc_r cc gcc)dnl
AC_PROG_CPP
dnl
dnl arch dependend

CFLAGS_64_BIT=""

case "$host" in
  *-hp-hpux*)
  if test $debug_defined = "yes"; then
    CFLAGS="-g -D__QC_HPUX"
  else
    CFLAGS="-O -D__QC_HPUX"
  fi
  CFLAGS_64_BIT="+DA2.0W"
  ;;

  *-dec-osf*)
  if test $debug_defined = "yes"; then
    CFLAGS="-g -D__QC_OSF"
  else
    CFLAGS="-O -D__QC_OSF"
  fi
  ;;

  *-sun-solaris*)
  if test $debug_defined = "yes"; then
    CFLAGS="-g -D__QC_SUN"
  else
    CFLAGS="-O -D__QC_SUN"
  fi
  ;;

  *-sgi-irix*)
  if test $debug_defined = "yes"; then
    CFLAGS="-g -D__QC_SGI"
  else
    CFLAGS="-O -D__QC_SGI"
  fi
  CFLAGS_64_BIT="-64"
  ;;
  
  *-ibm-aix*)
  if test -z "$GCC"; then
    if test $debug_defined = "yes"; then
      CFLAGS="-g -qthreaded -qmaxmem=8192 -D__QC_AIX"
    else
      CFLAGS="-O -qthreaded -qmaxmem=8192 -D__QC_AIX"
    fi
    CFLAGS_64_BIT="-q64"
  else
    if test $debug_defined = "yes"; then
      CFLAGS="-g -D__QC_AIX"
    else
      CFLAGS="-O -D__QC_AIX"
    fi
  fi
  ;;

  *-linux-*)
  if test $debug_defined = "yes"; then
    CFLAGS="-g -Wall -D__QC_LINUX -D__USE_UNIX98"
  else
    CFLAGS="-O -D__QC_LINUX -D__USE_UNIX98"
dnl Need -D__USE_UNIX98 for pthread.h definitions?
  fi
  ;;
  *)
  ;;
esac

dnl
dnl add flags for 64-bit mode
dnl
if test "x$enable_64_bit" = "xyes"; then
    CFLAGS="$CFLAGS $CFLAGS_64_BIT"
fi
