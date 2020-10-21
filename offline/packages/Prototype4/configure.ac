AC_INIT(prototype4,[1.00])
AC_CONFIG_SRCDIR([configure.ac])
AM_INIT_AUTOMAKE

AC_PROG_CXX(g++)
LT_INIT([disable-static])

case $CXX in
 clang++)
dnl boost accumulators
 CXXFLAGS="$CXXFLAGS -Wall -Werror -Wno-deprecated-declarations"
 ;;
esac

dnl test for gcc 8.3
if test $ac_cv_prog_gxx = yes; then
  if test `g++ -dumpversion | gawk '{print $1>=8.0?"1":"0"}'` = 1; then
   CXXFLAGS="$CXXFLAGS -Wall -Wno-deprecated-declarations -Werror"
  else
   CXXFLAGS="$CXXFLAGS -Wall -Werror"
  fi
fi

CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
