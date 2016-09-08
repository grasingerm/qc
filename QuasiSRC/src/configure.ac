dnl
dnl Configuration for main Quasi sources
dnl


dnl
dnl Initialize
dnl
AC_INIT
AC_CONFIG_SRCDIR([monitor.c])
AC_CONFIG_AUX_DIR(../bin)
AM_INIT_AUTOMAKE(QC3D, [" "])
AM_MAINTAINER_MODE
AC_CANONICAL_HOST
AC_SRC_DIR_EXPAND
AM_CONFIG_HEADER(config.h)
AC_PROG_MAKE_SET
AM_DISABLE_SHARED

dnl
dnl
dnl
builtin(include, ../system/debug.m4)dnl
builtin(include, ../system/64_bit.m4)dnl
builtin(include, ../system/ar.m4)dnl
builtin(include, ../system/dmalloc.m4)dnl

builtin(include, ../system/C/compiler.m4)dnl
builtin(include, ../system/C++/compiler.m4)dnl
builtin(include, ../system/F77/compiler.m4)dnl
				   

dnl
dnl Checks for programs.
dnl
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CC_C_O
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_INSTALL
AC_PROG_AWK

AC_LANG([C++])
AM_PROG_LIBTOOL

dnl
dnl Checks for libraries.
dnl
AC_CHECK_LIB(pthread, pthread_create)
AC_CHECK_LIB(pthread, __pthread_create_system)
AC_CHECK_LIB(m, sin)
AC_CHECK_LIB(socket, socket)
AC_CHECK_LIB(nsl, xdrstdio_create)
AC_CHECK_LIB(gen, basename)
AC_CHECK_LIB(posix4, sched_yield)


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
AC_TYPE_OFF_T
AC_TYPE_PID_T
AC_HEADER_TIME

dnl The following don't work, but we want them in site_element_map_cache.c
dnl AC_CHECK_TYPES(rwlock_t)
dnl AC_CHECK_TYPES(pthread_rwlock_t)


dnl
dnl Checks for specific definitions in headers
dnl

AC_MSG_CHECKING(for PTHREAD_SCOPE_SYSTEM in pthread.h)
AC_EGREP_CPP(yes, 
	[#include <pthread.h>
         #ifdef PTHREAD_SCOPE_SYSTEM
           yes
         #endif
        ], AC_DEFINE(HAVE_PTHREAD_SCOPE_SYSTEM,1,
    [Define if have PTHREAD_SCOPE_SYSTEM in pthread.h])
    AC_MSG_RESULT(yes),AC_MSG_RESULT(no))

AC_MSG_CHECKING(for _SC_NPROC_ONLN in unistd.h)
AC_EGREP_CPP(yes,
        [#include <unistd.h>
         #ifdef _SC_NPROC_ONLN
           yes
         #endif
        ], AC_DEFINE(HAVE__SC_NPROC_ONLN,1,
    [Define if have _SC_NPROC_ONLN in unistd.h])
    AC_MSG_RESULT(yes),AC_MSG_RESULT(no))

AC_MSG_CHECKING(for _SC_NPROCESSORS_ONLN in unistd.h)
AC_EGREP_CPP(yes,
        [#include <unistd.h>
         #ifdef _SC_NPROCESSORS_ONLN
           yes
         #endif
        ], AC_DEFINE(HAVE__SC_NPROCESSORS_ONLN,1,
    [Define if have _SC_NPROCESSORS_ONLN in unistd.h])
    AC_MSG_RESULT(yes),AC_MSG_RESULT(no))



dnl
dnl Check for specific functions
dnl
AC_FUNC_ALLOCA
AC_FUNC_FORK
dnl AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_FUNC_MMAP
AC_FUNC_STAT
AC_CHECK_FUNCS([floor gettimeofday memset mkdir munmap pathconf pow \
rint sqrt strchr strerror strstr])
AC_CHECK_FUNCS(pthread_attr_getguardsize,,)
AC_CHECK_FUNCS(pthread_num_processors_np,,)
AC_CHECK_FUNCS(pthread_launch_policy_np,,)
AC_CHECK_FUNCS(pthread_setconcurrency,,)
AC_CHECK_FUNCS([rwlock_init pthread_rwlock_init],,)
AC_CHECK_FUNCS([rw_wrlock pthread_rwlock_wrlock],,)
AC_CHECK_FUNCS([rw_unlock pthread_rwlock_unlock],,)
AC_CHECK_FUNCS([rw_rdlock pthread_rwlock_rdlock],,)


dnl
dnl Finish up
dnl
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
