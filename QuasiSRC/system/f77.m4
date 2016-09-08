#
# $Id: f77.m4,v 1.1 2002/03/08 00:02:43 knap Exp $
#
# $Log: f77.m4,v $
# Revision 1.1  2002/03/08 00:02:43  knap
# Updated build procedure to use automake/autoconf.
#
# Revision 1.1  2001/07/17 05:43:24  knap
# Added autoconf macros.
#
# Revision 1.1  2001/05/14 04:15:53  knap
# Initial sources.
#
#
#
# check if FORTRAN compiler appends underscore to symbol names
#
AC_MSG_CHECKING(whether $F77 appends underscore to symbol names)

if test "${ac_cv_have_f77_symbol_underscore+set}" = set; then
  if test "$ac_cv_have_f77_symbol_underscore" = "yes"; then
    echo "(cached) yes"
    AC_DEFINE(HAVE_F77_SYMBOL_UNDERSCORE)
  else
    echo "(cached) no"
  fi
else
  AC_LANG_SAVE
  AC_LANG_FORTRAN77
  AC_TRY_RUN_FORTRAN([
      program main
      call underscorefunction
      end
      subroutine underscorefunction
      return
      end],
      [ if test -n "`nm conftest 2>/dev/null | grep underscorefunction_ 2>/dev/null`"; then
          AC_MSG_RESULT(yes)
          AC_DEFINE(HAVE_F77_SYMBOL_UNDERSCORE)
          ac_cv_have_f77_symbol_underscore=yes
          elif test -n "`nm conftest 2>/dev/null | grep underscorefunction 2>/dev/null`"; then
          AC_MSG_RESULT(no)
          ac_cv_have_f77_symbol_underscore=no
          fi ],)dnl
  AC_LANG_RESTORE
fi

