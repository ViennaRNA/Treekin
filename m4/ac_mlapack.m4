AC_DEFUN([AC_TREEKIN_MLAPACK],[

dnl ---------------------------------------------------------------------------
dnl Check for mlapack
dnl ---------------------------------------------------------------------------

AC_ARG_ENABLE([mlapack],
              AC_HELP_STRING(
                  [--disable-mlapack],
                  [disable support for mlapack models (default=enabled)]
              ),
              [enable_mlapack=$enableval],
              [enable_mlapack=yes]
)

# MLAPACK package library path support, if not installed in usual directories
AC_ARG_WITH([mlapack],
            AC_HELP_STRING(
                [--with-mlapack=PREFIX],
                [alternative prefix path to mlapack library]
            ),
            MLAPACKPATHSET=1,
            MLAPACKPATHSET=0
)

AS_IF([test $MLAPACKPATHSET = 1],
      [ AX_APPEND_FLAG(CXXFLAGS, [-I$with_mlapack/include])
        AX_APPEND_FLAG(LDFLAGS, [-L$with_mlapack/lib])
      ])

mlapack_LIBS=""
mlapack_wanted_but_failed=""
mlapack_data_types=""

NEED_LIB_QD=0
NEED_LIB_GMP=0
NEED_LIB_MPC=0
NEED_LIB_MPFR=0

if test "$enable_mlapack" = "yes"; then
  if test  $MLAPACKPATHSET = 1 ; then
        AC_LANG_PUSH([C++])
        AC_CHECK_HEADER(mlapack/mpack_config.h,
            [],
            [
              mlapack_wanted_but_failed="( can't find header file 'mlapack/mpack_config.h' )";
              enable_mlapack="no";
            ],
            [
#include <mlapack/mpack_config.h>
            ])
        AC_LANG_POP
  else
        PKG_CHECK_MODULES([mlapack],
                          [mlapack >= 0.8.1],
                          [],
                          [
                            mlapack_wanted_but_failed="( mlapack library missing or not of version 0.8 or higher )";
                            enable_mlapack="no";
                          ])
  fi
fi

if test "$enable_mlapack" = "yes"; then
  AC_DEFINE([WITH_MLAPACK], [1], [MLAPACK support])

  #test which libraries of the mlapack package are installed and set the flags.
  AC_LANG_PUSH([C++])
  ac_CPPFLAGS_bak=$CPPFLAGS
  CPPFLAGS="$mlapack_CFLAGS $CPPFLAGS"

  AC_CHECK_HEADER(mlapack/mlapack_qd.h, [
      AC_DEFINE([WITH_MLAPACK_QD], [1], [Quad Double support])
      NEED_LIB_QD=1
  ], [],
  [#include <mlapack/mpack_config.h>])

  AC_CHECK_HEADER(mlapack/mlapack_dd.h, [
      AC_DEFINE([WITH_MLAPACK_DD], [1], [Double Double support])
      AC_TREEKIN_APPEND_VAR_COMMA(mlapack_data_types, [DD])
      NEED_LIB_QD=1
  ], [],
  [#include <mlapack/mpack_config.h>])

  AC_CHECK_HEADER(mlapack/mlapack_double.h, [
      AC_DEFINE([WITH_MLAPACK_DOUBLE], [1], [Double support])
      AC_TREEKIN_APPEND_VAR_COMMA(mlapack_data_types, [double])
  ], [],
  [#include <mlapack/mpack_config.h>])

  AC_CHECK_HEADER(mlapack/mlapack_longdouble.h, [
      AC_DEFINE([WITH_MLAPACK_LD], [1], [Long Double support])
      AC_TREEKIN_APPEND_VAR_COMMA(mlapack_data_types, [long double])
  ], [],
  [#include <mlapack/mpack_config.h>])

  AC_CHECK_HEADER(mlapack/mlapack___float128.h, [
      AC_DEFINE([WITH_MLAPACK___FLOAT128], [1], [__float128 support])
      AC_TREEKIN_APPEND_VAR_COMMA(mlapack_data_types, [__float128])
  ], [],
  [#include <mlapack/mpack_config.h>])

  AC_CHECK_HEADER(mlapack/mlapack_mpfr.h, [
      AC_DEFINE([WITH_MLAPACK_MPFR], [1], [MPFR support])
      AC_TREEKIN_APPEND_VAR_COMMA(mlapack_data_types, [MPFR])
      NEED_LIB_MPFR=1
      NEED_LIB_MPC=1
  ], [],
  [#include <mlapack/mpack_config.h>])

  AC_CHECK_HEADER(mlapack/mlapack_gmp.h, [
      AC_DEFINE([WITH_MLAPACK_GMP], [1], [GMP support])
      AC_TREEKIN_APPEND_VAR_COMMA(mlapack_data_types, [GMP])
      NEED_LIB_GMP=1
  ], [],
  [#include <mlapack/mpack_config.h>])

  CPPFLAGS=$ac_CPPFLAGS_bak

  AC_LANG_POP
fi


## try to determine libs to link against, if mlapack is not provided through MLAPACK_LIBS or pkg-config
if test "$MLAPACKPATHSET" = "1"; then
  # check for further dependencies of MLAPACK
  AC_LANG_PUSH([C++])
  if test "$NEED_LIB_GMP" = "1"; then
      AC_CHECK_HEADER(gmp.h, HAVE_GMP_HEADER=yes)

      if test "x$HAVE_GMP_HEADER" = "xyes"; then
          AC_CHECK_LIB(gmp, __gmpz_init, [mlapack_LIBS="$mlapack_LIBS -lgmpxx -lgmp"; HAVE_GMP=yes; AC_TREEKIN_APPEND_VAR_COMMA(mlapack_data_types, [GMP])])
      fi

      if test "x$HAVE_GMP" != "xyes"; then
          AC_MSG_RESULT([No GMP library with C++ wrapper found! You could try to set the CPPFLAGS=-I/path/to/gmp/include and LDFLAGS=-L/path/to/gmp/lib])
      fi
  fi

  if test "$NEED_LIB_MPFR" = "1"; then
      AC_CHECK_HEADER(mpfr.h, HAVE_MPFR_HEADER=yes)

      if test "x$HAVE_MPFR_HEADER" = "xyes"; then
          AC_CHECK_LIB(mpfr, mpfr_init, [mlapack_LIBS="$mlapack_LIBS -lmpfr"; HAVE_MPFR=yes; AC_TREEKIN_APPEND_VAR_COMMA(mlapack_data_types, [MPFR])])
      fi

      if test "x$HAVE_MPFR" != "xyes"; then
          AC_MSG_RESULT([No MPFR library found! You could try to set the CPPFLAGS=-I/path/to/mpfr/include and LDFLAGS=-L/path/to/mpfr/lib])
      fi
  fi

  if test "$NEED_LIB_MPC" = "1"; then
      AC_CHECK_HEADER(mpc.h, HAVE_MPC_HEADER=yes)

      if test "x$HAVE_MPC_HEADER" = "xyes"; then
          AC_CHECK_LIB(mpc, mpc_init2, [mlapack_LIBS="$mlapack_LIBS -lmpc"; HAVE_MPC=yes;  AC_TREEKIN_APPEND_VAR_COMMA(mlapack_data_types, [MPC])])
      fi

      if test "x$HAVE_MPC" != "xyes"; then
          AC_MSG_RESULT([No MPC library found! You could try to set the CPPFLAGS=-I/path/to/mpc/include and LDFLAGS=-L/path/to/mpc/lib])
      fi
  fi

  if test "$NEED_LIB_QD" = "1"; then
      AC_CHECK_HEADER(qd/qd_real.h, HAVE_QD_HEADER=yes)

      if test "x$HAVE_QD_HEADER" = "xyes"; then
          AC_CHECK_LIB(qd, c_qd_sqrt, [mlapack_LIBS="$mlapack_LIBS -lqd"; HAVE_QD=yes;  AC_TREEKIN_APPEND_VAR_COMMA(mlapack_data_types, [QD])])
      fi

      if test "x$HAVE_QD" != "xyes"; then
          AC_MSG_RESULT([No QD library found! You could try to set the CPPFLAGS=-I/path/to/qd/include and LDFLAGS=-L/path/to/qd/lib])
      fi
  fi
  AC_LANG_POP
fi


# error output if MLAPACK not found
if test "x$mlapack_wanted_but_failed" != "x"; then
    if test "$MLAPACKPATHSET" = "1"; then
        AC_MSG_WARN([
**********************************************************************
Failed to setup linking the program against the MLAPACK library.

Couldn't find mlapack in specified path '$with_mlapack'.
If you haven't installed it yet, you can obtain the mlapack library from

https://github.com/RaumZeit/mlapack

**********************************************************************
        ])
    else
        AC_MSG_WARN([
**********************************************************************
Failed to setup linking the program against the MLAPACK library.

In case you have installed mlapack in a non-standard directory, consider
using the

  --with-mlapack=/path/to/mlapack

parameter to specify the location where to find MLAPACK.
Otherwise, you can obtain the mlapack library from

https://github.com/RaumZeit/mlapack

**********************************************************************
        ])
    fi

fi
])
