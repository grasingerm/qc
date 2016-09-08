#
#-------------------------------------------------------------------------------
#                                        
#                Jaroslaw Knap, Michael Ortiz and Anna Pandolfi
#                      California Institute of Technology
#                         (C) 2001 All Rights Reserved
#                                        
#-------------------------------------------------------------------------------
# $Id: acinclude.m4,v 1.1 2002/03/08 00:02:43 knap Exp $
#
# $Log: acinclude.m4,v $
# Revision 1.1  2002/03/08 00:02:43  knap
# Updated build procedure to use automake/autoconf.
#
# Revision 1.2  2001/10/02 22:43:55  knap
# Imported sources.
#
# Revision 1.1  2001/05/21 18:48:24  knap
# Initial source of local m4 macro file.
#
#
AC_DEFUN([AC_SRC_DIR_EXPAND], [
# expand $srcdir to an absolute path
srcdir=`CDPATH=:; cd $srcdir && pwd`
])

#
# try run a fortran program
#
AC_DEFUN([AC_TRY_RUN_FORTRAN],
[
[cat > conftest.$ac_ext <<EOF]
[$1]
EOF
if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext} && (./conftest; exit) 2>/dev/null
then
dnl Don't remove the temporary files here, so they can be examined.
  ifelse([$2], , :, [$2])
else
  echo "configure: failed program was:" >&AC_FD_CC
  cat conftest.$ac_ext >&AC_FD_CC
ifelse([$3], , , [  rm -fr conftest*
  $3
])dnl
fi
rm -fr conftest*
])
