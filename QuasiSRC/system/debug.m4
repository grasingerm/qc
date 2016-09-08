dnl
dnl-----------------------------------------------------------------------------
dnl                                      
dnl              Jaroslaw Knap, Michael Ortiz and Anna Pandolfi
dnl                    California Institute of Technology
dnl                       (C) 2001 All Rights Reserved
dnl                                      
dnl-----------------------------------------------------------------------------
dnl $Id: debug.m4,v 1.1 2002/03/08 00:02:43 knap Exp $
dnl
dnl $Log: debug.m4,v $
dnl Revision 1.1  2002/03/08 00:02:43  knap
dnl Updated build procedure to use automake/autoconf.
dnl
dnl Revision 1.2  2001/10/02 22:43:55  knap
dnl Imported sources.
dnl
dnl Revision 1.1  2001/05/06 06:07:33  knap
dnl Initial sources.
dnl
dnl
AH_TEMPLATE([DEBUG], [Define if compiling for debug])
dnl
AC_MSG_CHECKING(whether to enable debugging)
debug_defined=no
AC_ARG_ENABLE(debug,
[  --enable-debug          Enable debugging.],
[ case "$enableval" in
  yes)
    AC_MSG_RESULT(yes)
    AC_DEFINE(DEBUG)
    debug_defined=yes
    ;;
  *)
    AC_MSG_RESULT(no)
    ;;
  esac ], AC_MSG_RESULT(no)
)
