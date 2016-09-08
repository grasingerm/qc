dnl
dnl-----------------------------------------------------------------------------
dnl                                      
dnl              Jaroslaw Knap, Michael Ortiz and Anna Pandolfi
dnl                    California Institute of Technology
dnl                       (C) 2001 All Rights Reserved
dnl                                      
dnl-----------------------------------------------------------------------------
dnl $Id: compiler.m4,v 1.3 2002/12/06 01:49:48 fago Exp $
dnl
dnl $Log: compiler.m4,v $
dnl Revision 1.3  2002/12/06 01:49:48  fago
dnl Significant fixes to automake build especially on AIX and HPUX.
dnl
dnl Revision 1.2  2002/12/04 17:45:46  fago
dnl Ported to automake 1.7.1, autoconf 2.56. perror -> ERROR. Fixed build bugs,
dnl especially on HPUX.
dnl
dnl Revision 1.3  2002/10/21 17:19:45  fago
dnl Tweaking autotools build.
dnl
dnl Revision 1.2  2002/10/21 14:20:23  fago
dnl Rewrite of mtwist library and touch-up to langevin_force.*
dnl
dnl Revision 1.1.1.1  2002/04/22 03:33:47  knap
dnl Initial import
dnl
dnl Revision 1.1  2002/03/08 00:02:46  knap
dnl Updated build procedure to use automake/autoconf.
dnl
dnl Revision 1.3  2001/12/26 18:44:58  knap
dnl Corrected support for GCC.
dnl
dnl Revision 1.2  2001/10/02 22:43:56  knap
dnl Imported sources.
dnl
dnl Revision 1.8  2001/05/14 18:33:22  knap
dnl Ported to Compaq's Alpha.
dnl
dnl Revision 1.7  2001/05/14 05:28:20  knap
dnl Added additional preprocessor flags on HP-UX.
dnl
dnl Revision 1.6  2001/05/14 05:15:40  knap
dnl Ported to HP-UX.
dnl
dnl Revision 1.5  2001/05/14 04:17:11  knap
dnl Added CPPFLAGS on irix.
dnl
dnl Revision 1.4  2001/05/10 20:50:17  knap
dnl Corrected setting of 64-bit compiler mode options.
dnl
dnl Revision 1.3  2001/05/08 23:35:51  knap
dnl Ported automake/autoconf files to IBM AIX.
dnl
dnl Revision 1.2  2001/05/06 06:06:18  knap
dnl Updated to include more platforms.
dnl
dnl Revision 1.1  2001/05/04 23:37:38  knap
dnl Original sources added.
dnl
dnl
dnl AC_REQUIRE([AC_CANONICAL_HOST])

dnl
dnl Compiler definitions
dnl

dnl
dnl Various macro definitions
dnl

dnl
dnl reset the order of checks for C++ to native then GNU
dnl
AC_PROG_CXX(CC aCC cxx xlC_r c++ gcc cc++ cl g++)dnl
AC_PROG_CXXCPP
dnl
dnl arch dependend
dnl
CXXFLAGS_64_BIT=""

case "$host" in
  *-hp-hpux*)
  if test -z "$GCC"; then
    if test $debug_defined = "yes"; then
      CXXFLAGS="-AA -g -D_HPUX_SOURCE -D_RWSTD_MULTI_THREAD -D_REENTRANT -D__QC_HPUX"
    else
      CXXFLAGS="-AA -O -D_HPUX_SOURCE -D_RWSTD_MULTI_THREAD -D_REENTRANT -D__QC_HPUX"
    fi
    CXXFLAGS_64_BIT="+DA2.0W"
    CXXCPP="$CXXCPP -AA -D_RWSTD_MULTI_THREAD -D_REENTRANT -D__QC_HPUX"
  else
    if test $debug_defined = "yes"; then
      CXXFLAGS="-g -D__QC_HPUX"
    else
      CXXFLAGS="-O -D__QC__HPUX"
    fi
  fi
  ;;

  *-dec-osf*)
  if test $debug_defined = "yes"; then
    CXXFLAGS="-g -D__QC_OSF"
  else
    CXXFLAGS="-O -D__QC_OSF"
  fi
  ;;

  *-sun-solaris*)
  if test $debug_defined = "yes"; then
    CXXFLAGS="-g -D__QC_SUN"
  else
    CXXFLAGS="-O -D__QC_SUN"
  fi
  ;;

  *-sgi-irix*)
  if test -z "$GCC"; then
    if test $debug_defined = "yes"; then
      CXXFLAGS="-g -LANG:std -D__QC_SGI"
    else
      CXXFLAGS="-O -LANG:std -D__QC_SGI"
    fi
    CPPFLAGS="-LANG:std"
    CXXFLAGS_64_BIT="-64"
    AR="$CXX -ar -o"
    AR_FLAGS=" "
    LD=$CXX
    CXX="$CXX -LANG:std"
  else
    if test $debug_defined = "yes"; then
      CXXFLAGS="-g"
    else
      CXXFLAGS="-O"
    fi
  fi
  ;;

  *-ibm-aix*)
  if test -z "$GCC"; then
    if test $debug_defined = "yes"; then
      CXXFLAGS="-g -qthreaded -qmaxmem=8192 -D__QC_AIX"
    else
      CXXFLAGS="-O -qthreaded -qmaxmem=8192 -D__QC_AIX"
    fi
    CPPFLAGS="-qlanglvl=extended -D__QC_AIX"
    CXXFLAGS_64_BIT="-q64"
    ac_cv_cxx_exception_flag=no
  else
    if test $debug_defined = "yes"; then
      CXXFLAGS="-g -D__QC_AIX"
    else
      CXXFLAGS="-O -D__QC_AIX"
    fi
  fi
  ;;

  *-linux-*)
  if test $debug_defined = "yes"; then
    CXXFLAGS="-g -Wall -D__QC_LINUX -D__USE_UNIX98"
  else
    CXXFLAGS="-O -D__QC_LINUX -D__USE_UNIX98"
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
  CXXFLAGS="$CXXFLAGS $CXXFLAGS_64_BIT"
fi
