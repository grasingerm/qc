dnl
dnl-----------------------------------------------------------------------------
dnl                                      
dnl              Jaroslaw Knap, Michael Ortiz and Anna Pandolfi
dnl                    California Institute of Technology
dnl                       (C) 2001 All Rights Reserved
dnl                                      
dnl-----------------------------------------------------------------------------
dnl $Id: 64_bit.m4,v 1.3 2003/07/07 19:29:53 knap Exp $ 
dnl 
dnl
dnl $Log: 64_bit.m4,v $
dnl Revision 1.3  2003/07/07 19:29:53  knap
dnl Introduced WITH_64_BIT_ENABLED define set when compiling in the 64-bit
dnl mode.
dnl
dnl Revision 1.2  2002/12/06 01:49:47  fago
dnl Significant fixes to automake build especially on AIX and HPUX.
dnl
dnl Revision 1.1  2002/03/08 00:02:43  knap
dnl Updated build procedure to use automake/autoconf.
dnl
dnl Revision 1.2  2001/10/02 22:43:55  knap
dnl Imported sources.
dnl
dnl Revision 1.1  2001/05/10 21:18:13  knap
dnl Initial sources.
dnl
dnl
dnl AC_REQUIRE([AC_CANONICAL_HOST])

case "$host" in
  *-hp-hpux*)
  HAVE_64_BIT="yes"
  ;;

  *-sgi-irix*)
  HAVE_64_BIT="yes"
  ;;

  *-ibm-aix*)
  HAVE_64_BIT="yes"
  # export OBJECT_MODE=64 # this is used by AR 
  ;;

  *)
  HAVE_64_BIT="no"
  ;;
esac

AC_MSG_CHECKING(whether to use 64-bit compiler mode)
AC_ARG_ENABLE(64-bit,
  [  --enable-64-bit Use 64-bit compiler mode],
  [ if test "$HAVE_64_BIT" = "no"; then
      AC_ERROR("$host" does not support 64-bit executables.)
    fi
    AC_DEFINE(WITH_64_BIT_ENABLED,1,[Define if compiling for 64-bit])
    AC_MSG_RESULT(yes) ],
  [ AC_MSG_RESULT(no) ]
)
