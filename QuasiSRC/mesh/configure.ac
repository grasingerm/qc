dnl
dnl Configuration for mesh_generation util programs
dnl


dnl
dnl Initialize
dnl
AC_INIT
AC_CONFIG_SRCDIR([mesh.cc])
AC_CONFIG_AUX_DIR(../bin)
AM_INIT_AUTOMAKE(MESH_GENERATE, [" "])
AM_MAINTAINER_MODE
AC_CANONICAL_HOST
AC_SRC_DIR_EXPAND
AM_CONFIG_HEADER(config.h)
AC_PROG_MAKE_SET
AM_DISABLE_SHARED

dnl
dnl Jarek's system/compiler checks
dnl
builtin(include, ../system/debug.m4)dnl
builtin(include, ../system/64_bit.m4)dnl
builtin(include, ../system/ar.m4)dnl

builtin(include, ../system/C/compiler.m4)dnl
builtin(include, ../system/C++/compiler.m4)dnl
builtin(include, ../system/F77/compiler.m4)dnl
				   

dnl
dnl Checks for programs.
dnl
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CC_C_O
AC_PROG_LIBTOOL
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_LN_S

dnl
dnl Checks for libraries.
dnl
AC_CHECK_LIB(pthread, pthread_create)
AC_CHECK_LIB(m, sin)
AC_CHECK_LIB(nsl, xdrstdio_create)


dnl
dnl Checks for headers
dnl 
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS(pthread.h sched.h sys/mpctl.h errno.h unistd.h \
math.h assert.h search.h signal.h sys/wait.h ctype.h limits.h \
fcntl.h sys/mman.h libgen.h rpc/rpc.h rpc/xdr.h alloca.h sys/time.h \
sys/resource.h float.h synch.h values.h)
AC_CHECK_HEADER(cassert, ,CXXFLAGS="$CXXFLAGS -I${srcdir}/../system/C++")
AC_CHECK_HEADERS(algorithm algo.h fstream fstream.h vector vector.h limits)


dnl
dnl Checks for typedefs, structures, and compiler characteristics.
dnl
AC_C_CONST
AC_TYPE_SIZE_T
AC_HEADER_TIME


dnl
dnl Check for specific functions
dnl
dnl AC_FUNC_MALLOC
AC_CHECK_FUNCS([pow strstr rint sqrt strchr])


dnl
dnl Finish up
dnl
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
