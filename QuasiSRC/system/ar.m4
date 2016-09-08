dnl
dnl-----------------------------------------------------------------------------
dnl                                      
dnl              Jaroslaw Knap, Michael Ortiz and Anna Pandolfi
dnl                    California Institute of Technology
dnl                       (C) 2001 All Rights Reserved
dnl                                      
dnl-----------------------------------------------------------------------------
dnl $Id: ar.m4,v 1.3 2002/12/06 01:49:48 fago Exp $
dnl
dnl $Log: ar.m4,v $
dnl Revision 1.3  2002/12/06 01:49:48  fago
dnl Significant fixes to automake build especially on AIX and HPUX.
dnl
dnl Revision 1.6  2002/05/13 05:45:27  knap
dnl Synced automake/autoconf files.
dnl
dnl Revision 1.5  2002/05/09 19:09:10  knap
dnl Commented AR_FLAGS defs on AIX.
dnl
dnl Revision 1.4  2002/03/11 23:44:04  knap
dnl Changed handling of AR and AR_FLAGS.
dnl
dnl Revision 1.3  2002/03/11 20:59:24  knap
dnl Moved ar flags from AR to AR_FLAGS.
dnl
dnl Revision 1.2  2001/10/02 22:43:55  knap
dnl Imported sources.
dnl
dnl Revision 1.2  2001/05/14 04:21:40  knap
dnl Default flags were missing for AR-corrected.
dnl
dnl Revision 1.1  2001/05/10 21:18:14  knap
dnl Initial sources.
dnl
dnl
dnl AC_REQUIRE([AC_CANONICAL_HOST])

dnl
dnl default value
dnl
AR="ar"
AR_FLAGS="cru"
AR_FLAGS_64=""

dnl
dnl platform specific data
dnl

case "$host" in

  *-ibm-aix*)
    AR_FLAGS_64="-X64"
  ;;

  *)
  ;;

esac

dnl
dnl add flags for 64-bit mode
dnl
if test "x$enable_64_bit" = "xyes"; then
  AR_FLAGS="$AR_FLAGS_64 $AR_FLAGS"
fi


AC_SUBST(AR)
AC_SUBST(AR_FLAGS)

