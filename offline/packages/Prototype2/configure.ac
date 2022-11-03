AC_INIT(prototype2,[1.00])
AC_CONFIG_SRCDIR([configure.ac])
AM_INIT_AUTOMAKE

AC_PROG_CXX(CC g++)
LT_INIT([disable-static])

CXXFLAGS="$CXXFLAGS -Wall -Werror -Wno-unused-parameter"

case $CXX in
 clang++)
  CXXFLAGS="$CXXFLAGS -Wno-bitwise-instead-of-logical -Wno-unknown-warning-option"
 ;;
esac

CINTDEFS=" -noIncludePaths  -inlineInputHeader "
AC_SUBST(CINTDEFS)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
