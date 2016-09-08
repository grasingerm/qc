dnl
dnl Configuration for general library
dnl



dnl
dnl Initialize
dnl
AC_INIT
AC_CONFIG_AUX_DIR(../../../bin)
AC_CONFIG_SRCDIR([dhpsrt.f])
AM_INIT_AUTOMAKE(general, [" "])
AM_MAINTAINER_MODE
AC_CANONICAL_HOST
AC_SRC_DIR_EXPAND
AM_CONFIG_HEADER(config.h)
AC_PROG_MAKE_SET
AM_DISABLE_SHARED

dnl
dnl
dnl
builtin(include, ../../../system/debug.m4)dnl
builtin(include, ../../../system/64_bit.m4)dnl
builtin(include, ../../../system/ar.m4)dnl

builtin(include, ../../../system/C/compiler.m4)dnl
builtin(include, ../../../system/F77/compiler.m4)dnl

dnl
dnl Checks for programs.
dnl
AC_PROG_CC
AC_PROG_F77
AC_PROG_LIBTOOL

dnl
dnl Checks for libraries.
dnl

dnl
dnl Checks for headers
dnl 

dnl
dnl Checks for typedefs, structures, and compiler characteristics.
dnl

dnl
dnl Checks for specific definitions in headers
dnl

dnl
dnl Check for specific functions
dnl

dnl
dnl Finish up
dnl
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
