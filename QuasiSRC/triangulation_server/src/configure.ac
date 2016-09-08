dnl
dnl Configuration for triangulation server
dnl



dnl
dnl Initialize
dnl
AC_INIT
AC_CONFIG_AUX_DIR(../../bin)
AC_CONFIG_SRCDIR([triangulation_server.c])
AM_INIT_AUTOMAKE(triangulation server, [" "])
AM_MAINTAINER_MODE
AC_CANONICAL_HOST
AC_SRC_DIR_EXPAND
AM_CONFIG_HEADER(config.h)
AC_PROG_MAKE_SET
AM_DISABLE_SHARED


dnl
dnl Jarek's system/compiler checks
dnl
builtin(include, ../../system/debug.m4)dnl
builtin(include, ../../system/64_bit.m4)dnl
builtin(include, ../../system/ar.m4)dnl

builtin(include, ../../system/C/compiler.m4)dnl
builtin(include, ../../system/C++/compiler.m4)dnl
builtin(include, ../../system/F77/compiler.m4)dnl


dnl
dnl Checks for programs.
dnl
AC_PROG_CC
AC_PROG_CPP
AC_PROG_F77
AC_PROG_CC_C_O
AC_PROG_LIBTOOL
AC_PROG_AWK
AC_PROG_INSTALL


dnl
dnl Fortran support.
dnl
AC_F77_WRAPPERS
AC_F77_LIBRARY_LDFLAGS


dnl
dnl Checks for libraries.
dnl
AC_CHECK_LIB(pthread, pthread_create)
AC_CHECK_LIB(pthread, __pthread_create_system)
AC_CHECK_LIB(m, sin)
AC_CHECK_LIB(socket, socket)


dnl
dnl Checks for headers
dnl 
AC_HEADER_STDC
AC_CHECK_HEADERS(pthread.h sched.h sys/mpctl.h errno.h unistd.h \
math.h assert.h search.h signal.h sys/wait.h ctype.h limits.h \
fcntl.h sys/mman.h libgen.h rpc/rpc.h rpc/xdr.h alloca.h sys/time.h \
sys/resource.h sys/utsname.h sys/socket.h sys/un.h)


dnl
dnl Checks for typedefs, structures, and compiler characteristics.
dnl
AC_C_CONST
AC_TYPE_SIZE_T
AC_TYPE_PID_T


dnl
dnl Checks for specific definitions in headers
dnl


dnl
dnl Check for specific functions
dnl
AC_FUNC_FORK
dnl AC_FUNC_MALLOC
AC_CHECK_FUNCS([atexit memset socket sqrt strchr uname])
AC_CHECK_FUNCS(kmod)

dnl
dnl Finish up
dnl
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
