dnl @synopsis ARES_WITH_DMALLOC
dnl
dnl The dmalloc library is an useful tool to debug memory problems in
dnl your programs, but you don't really want to compile
dnl dmalloc-support into every binary you produce, because dmalloc
dnl brings performance loss.
dnl
dnl The ARES_WITH_DMALLOC macro defines a user switch '--with-dmalloc'
dnl which can be used likes this:
dnl
dnl      ./configure --with-dmalloc[=BASE-PATH]
dnl
dnl If no BASE-PATH has been provided, "/usr/local" will be used as
dnl default.
dnl
dnl The BASE-PATH is the place where autoconf expects to find the
dnl include- and link library files of dmalloc, specifically in:
dnl
dnl      $(BASE-PATH)/include
dnl      $(BASE-PATH)/lib
dnl
dnl If dmalloc-support has been enabled, the pre-processor defines
dnl "WITH_DMALLOC" will added to CPPFLAGS
dnl as well as the apropriate "-I" statement.
dnl
dnl Use the first define in your source codes to determine whether you
dnl have dmalloc support enabled or not. Usually something like this
dnl will suffice:
dnl
dnl      #ifdef WITH_DMALLOC
dnl      #  include <dmalloc.h>
dnl      #endif
dnl
dnl You will find dmalloc at <http://www.dmalloc.com/>.
dnl
dnl @version $Id: dmalloc.m4,v 1.2 2003/08/06 18:32:59 knap Exp $
dnl @author Peter Simons <simons@computer.org>
dnl
AC_MSG_CHECKING(whether to use dmalloc library)
AC_ARG_WITH(dmalloc,
[  --with-dmalloc[=ARG]      compile with dmalloc library],
if test "$withval" = "" -o "$withval" = "yes"; then
    ac_cv_dmalloc="/usr/local"
else
    ac_cv_dmalloc="$withval"
fi
AC_MSG_RESULT(yes)
AC_DEFINE(WITH_DMALLOC,1,
          [Define if using the dmalloc debugging malloc package])
CPPFLAGS="$CPPFLAGS -I$ac_cv_dmalloc/include"
CFLAGS="$CFLAGS -I$ac_cv_dmalloc/include"
CXXFLAGS="$CXXFLAGS -I$ac_cv_dmalloc/include"
LDFLAGS="$LDFLAGS -L$ac_cv_dmalloc/lib"
DMALLOC_C_LIBS="-ldmalloc"
DMALLOC_CXX_LIBS="-ldmallocxx"
AC_SUBST(DMALLOC_C_LIBS)
AC_SUBST(DMALLOC_CXX_LIBS)
,AC_MSG_RESULT(no))
