dnl
dnl-----------------------------------------------------------------------------
dnl                                      
dnl              Jaroslaw Knap, Michael Ortiz and Anna Pandolfi
dnl                    California Institute of Technology
dnl                       (C) 2001 All Rights Reserved
dnl                                      
dnl-----------------------------------------------------------------------------
dnl
dnl $Id: inventor.m4,v 1.2 2002/12/06 01:49:48 fago Exp $
dnl
dnl $Log: inventor.m4,v $
dnl Revision 1.2  2002/12/06 01:49:48  fago
dnl Significant fixes to automake build especially on AIX and HPUX.
dnl
dnl Revision 1.1  2002/03/08 00:02:44  knap
dnl Updated build procedure to use automake/autoconf.
dnl
dnl Revision 1.1  2002/01/02 01:53:06  knap
dnl Initial revision.
dnl
dnl
AC_DEFUN([AC_CHECK_INVENTOR],
[ AC_MSG_CHECKING([whether to add support for Open Inventor])
  AC_ARG_WITH([inventor],
  [  [--with-inventor[=PATH]]  use Open Inventor (default is no)],
[ case "$withval" in
  no)
    AC_MSG_RESULT(no)
    ;;
  yes)
    AC_MSG_RESULT(yes)
    if test '!' -d /usr/include/Inventor; then
      AC_ERROR(Inventor headers not found in /usr/include/Inventor: you must supply the path.)
    fi
    AC_MSG_RESULT(yes)
    AC_MSG_RESULT(Assuming Inventor headers  and libraries are in $withval/include/Inventor and $withval/lib.)
    AC_DEFINE(HAVE_INVENTOR)
    LIBS="$LIBS -lInventor -lInventorXt"
    ;;
  *)
    AC_MSG_RESULT(yes)
    if test '!' -d $withval/include/Inventor; then
      AC_ERROR(Inventor headers not found in $withval/include/Inventor: please supply the correct path.)
    fi
    AC_MSG_RESULT(Assuming Inventor headers and libraries are in $withval/include/Inventor and $withval/lib.)
    AC_DEFINE([HAVE_INVENTOR], 0, [Define if Open Inventor is available.])
    CXXFLAGS="$CXXFLAGS -I$withval/include"
    LIBS="$LIBS -L$withval/lib -lInventor -lInventorXt"
    ;;
  esac ],
  [AC_MSG_RESULT(no)]
)])
