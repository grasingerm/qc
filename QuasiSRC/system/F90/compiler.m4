dnl
dnl-----------------------------------------------------------------------------
dnl                                      
dnl              Jaroslaw Knap, Michael Ortiz and Anna Pandolfi
dnl                    California Institute of Technology
dnl                       (C) 2001 All Rights Reserved
dnl                                      
dnl-----------------------------------------------------------------------------
dnl Compiler definitions
dnl

dnl
dnl reset the order of checks for compilers to native, gcc
dnl
AC_CHECK_PROGS(CC, $CC cc gcc, gcc)dnl
AC_CHECK_PROGS(CXX, $CCC CC aCC c++ g++ gcc cxx cc++ cl, gcc)dnl
AC_CHECK_PROGS(F77, f77 g77 f2c)dnl

dnl
dnl arch dependend

FFLAGS_64=""
CFLAGS_64=""

case "$host" in
  *-hp-hpux*)
  # if we are using f77 add +ppu
  if test "$F77" = "f77"; then
    FFLAGS="-O +ppu"
  fi
  CFLAGS_64_BIT="+DA2.0W"
  ;;
  *-dec-osf*)
  # if using f77 adjust flags
  if test "$F77" = "f77"; then
    FFLAGS="-O"
  fi
  ;;
  *-sun-solaris*)
  # if using f77 adjust flags
  if test "$F77" = "f77"; then
    FFLAGS="-O"
  fi
  ;;
  *-sgi-irix*)
  # 
  # use MIPSPro compilers as default
  #
  
  ;;
  *-ibm-aix*)
  # if using f77 adjust flags
  if test "$F77" = "f77"; then
    FFLAGS="-O -qextname -qmaxmem=8192"
    FFLAGS_64="-q64"
  fi
  if test -z "$GCC"; then
    CFLAGS="-O"
    CFLAGS_64="-q64"
  fi
  ;;

  *)
  ;;
esac
