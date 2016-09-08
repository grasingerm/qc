dnl
dnl $Id: link.m4,v 1.1 2002/03/08 00:02:44 knap Exp $
dnl
dnl $log
dnl
AC_REQUIRE([AC_CANONICAL_HOST])

dnl
dnl
dnl

case "$host" in
  *-sgi-irix*)
  LINK='${CXX} ${AM_CXXFLAGS} ${CXXFLAGS} ${AM_LDFLAGS} ${LDFLAGS}  -o $@'
  ;;

  *-linux-*)
  LINK='${F77LINK} -lstdc++ -lm -lgcc -lc -lgcc'
  ;;

  *)
  LINK='${F77LINK}'
  ;;
esac

AC_SUBST(LINK)


