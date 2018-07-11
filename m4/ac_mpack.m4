AC_DEFUN([AC_TREEKIN_MPACK],[

dnl ---------------------------------------------------------------------------
dnl Check for mlapack
dnl ---------------------------------------------------------------------------

AC_ARG_ENABLE([MPACK],
              AC_HELP_STRING(
                  [--disable-MPACK],
                  [disable support for MPACK models (default=enabled)]
              ),
              [enable_mpack=$enableval],
              [enable_mpack=yes]
)

# MPACK package library path support, if not installed in usual directories
AC_ARG_WITH([MPACK],
            AC_HELP_STRING(
                [--with-MPACK=PREFIX],
                [alternative prefix path to MPACK library]
            ),
            MPACKPATHSET=1,
            MPACKPATHSET=0
)

AS_IF([test $MPACKPATHSET = 1],
      [ CXXFLAGS="-I$with_MPACK/include $CXXFLAGS";
        LDFLAGS="-L$with_MPACK/lib $LDFLAGS" ])

mpack_LIBS=""
mpack_wanted_but_failed=""
mpack_DIR=""

NEED_LIB_QD=0
NEED_LIB_GMP=0
NEED_LIB_MPC=0
NEED_LIB_MPFR=0

if test "$enable_mpack" = "yes"; then
  if test  $MPACKPATHSET = 1 ; then
        AC_LANG_PUSH([C++])
        AC_CHECK_HEADER(mpack/mpack_config.h,
            [],
            [
              mpack_wanted_but_failed="( can't find header file 'mpack/mpack_config.h' )";
              enable_mpack="no";
            ],
            [
#include <mpack/mpack_config.h>
            ])
        AC_LANG_POP
  else
        PKG_CHECK_MODULES([mpack],
                          [mpack >= 0.8],
                          [
                            CXXFLAGS="$mpack_CFLAGS $CXXFLAGS"
                            # PKG_CHECK_MODULES sets mpack_LIBS automatically
                            mpack_DIR=`pkg-config --variable includedir mpack`
                          ],
                          [
                            mpack_wanted_but_failed="( mpack library missing or not of version 0.8 or higher )";
                            enable_mpack="no";
                          ])
  fi
fi

if test "$enable_mpack" = "yes"; then
  ppFlags=" "
##  if test "x$mpack_DIR" != "x"; then
##    ac_save_CPPFLAGS="$CPPFLAGS"
##    CPPFLAGS="-I${mpack_DIR}/mpack $CPPFLAGS"
##  fi

  #test which libraries of the mpack package are installed and set the flags.
  AC_LANG_PUSH([C++])
  AC_CHECK_HEADER(mpack/mlapack_qd.h, [ppFlags="-DWITH_MPACK_QD $ppFlags"; NEED_LIB_QD=1 ],[],
  [ #include <mpack/mpack_config.h>
   ])
  AC_CHECK_HEADER(mpack/mlapack_dd.h, [ppFlags="-DWITH_MPACK_DD $ppFlags"; NEED_LIB_QD=1 ],[],
  [#include <mpack/mpack_config.h>
   ])
  AC_CHECK_HEADER(mpack/mlapack_double.h, [ppFlags="-DWITH_MPACK_DOUBLE $ppFlags" ],[],
  [#include <mpack/mpack_config.h>
   ])
  AC_CHECK_HEADER(mpack/mlapack_longdouble.h, [ppFlags="-DWITH_MPACK_LD $ppFlags" ],[],
  [#include <mpack/mpack_config.h>
   ])
  AC_CHECK_HEADER(mpack/mlapack___float128.h, [ppFlags="-DWITH_MPACK___FLOAT128 $ppFlags" ],[],
  [#include <mpack/mpack_config.h>
   ])
  AC_CHECK_HEADER(mpack/mlapack_mpfr.h, [ppFlags="-DWITH_MPACK_MPFR $ppFlags"; NEED_LIB_MPFR=1; NEED_LIB_MPC=1 ],[],
  [#include <mpack/mpack_config.h>
   ])
  AC_CHECK_HEADER(mpack/mlapack_gmp.h, [ppFlags="-DWITH_MPACK_GMP $ppFlags"; NEED_LIB_GMP=1 ],[],
  [#include <mpack/mpack_config.h>
   ])
  AC_LANG_POP
  
##  if test "x$mpack_DIR" != "x"; then
##    CPPFLAGS="$ac_save_CPPFLAGS"
##  fi

  CPPFLAGS=$CPPFLAGS$ppFlags
fi


if test "$enable_mpack" = "yes"; then
  CPPFLAGS=" -DWITH_MPACK $CPPFLAGS"

  # check for further dependencies of MPACK
  AC_LANG_PUSH([C++])
  if test "$NEED_LIB_GMP" = "1"; then
      AC_CHECK_HEADER(gmp.h, HAVE_GMP_HEADER=yes,[],
      [
#if HAVE_SUPPRESS_COMPILER_WARNING_H
#include <SUPPRESS_COMPILER_WARNING.h>
#endif
      ])

      if test "x$HAVE_GMP_HEADER" = "xyes"; then
          AC_CHECK_LIB(gmp, __gmpz_init, [mpack_LIBS="$mpack_LIBS -lgmpxx -lgmp"; HAVE_GMP=yes])
      fi

      if test "x$HAVE_GMP" != "xyes"; then
          AC_MSG_RESULT([No GMP library with C++ wrapper found! You could try to set the CPPFLAGS=-I/path/to/gmp/include and LDFLAGS=-L/path/to/gmp/lib])
      fi
  fi

  if test "$NEED_LIB_MPFR" = "1"; then
      AC_CHECK_HEADER(mpfr.h, HAVE_MPFR_HEADER=yes,[],
      [
#if HAVE_SUPPRESS_COMPILER_WARNING_H
#include <SUPPRESS_COMPILER_WARNING.h>
#endif
      ])

      if test "x$HAVE_MPFR_HEADER" = "xyes"; then
          AC_CHECK_LIB(mpfr, mpfr_init, [mpack_LIBS="$mpack_LIBS -lmpfr"; HAVE_MPFR=yes])
      fi

      if test "x$HAVE_MPFR" != "xyes"; then
          AC_MSG_RESULT([No MPFR library found! You could try to set the CPPFLAGS=-I/path/to/mpfr/include and LDFLAGS=-L/path/to/mpfr/lib])
      fi
  fi

  if test "$NEED_LIB_MPC" = "1"; then
      AC_CHECK_HEADER(mpc.h, HAVE_MPC_HEADER=yes,[],
      [
#if HAVE_SUPPRESS_COMPILER_WARNING_H
#include <SUPPRESS_COMPILER_WARNING.h>
#endif
      ])

      if test "x$HAVE_MPC_HEADER" = "xyes"; then
          AC_CHECK_LIB(mpc, mpc_init2, [mpack_LIBS="$mpack_LIBS -lmpc"; HAVE_MPC=yes])
      fi

      if test "x$HAVE_MPC" != "xyes"; then
          AC_MSG_RESULT([No MPC library found! You could try to set the CPPFLAGS=-I/path/to/mpc/include and LDFLAGS=-L/path/to/mpc/lib])
      fi
  fi

  if test "$NEED_LIB_QD" = "1"; then
      AC_CHECK_HEADER(qd/qd_real.h, HAVE_QD_HEADER=yes,[],
      [
#if HAVE_SUPPRESS_COMPILER_WARNING_H
#include <SUPPRESS_COMPILER_WARNING.h>
#endif
      ])

      if test "x$HAVE_QD_HEADER" = "xyes"; then
          AC_CHECK_LIB(qd, c_qd_sqrt, [mpack_LIBS="$mpack_LIBS -lqd"; HAVE_QD=yes])
      fi

      if test "x$HAVE_QD" != "xyes"; then
          AC_MSG_RESULT([No QD library found! You could try to set the CPPFLAGS=-I/path/to/qd/include and LDFLAGS=-L/path/to/qd/lib])
      fi
  fi
  AC_LANG_POP
fi


# error output if MPACK not found
if test "x$mpack_wanted_but_failed" != "x"; then
    if test "$MPACKPATHSET" = "1"; then
        AC_MSG_WARN([
**********************************************************************
Failed to setup linking the program against the MPACK library.

Couldn't find MPack in specified path '$with_MPACK'.
If you haven't installed it yet, you can obtain the MPack library from

https://github.com/nakatamaho/mpack

**********************************************************************
        ])
    else
        AC_MSG_WARN([
**********************************************************************
Failed to setup linking the program against the MPACK library.

In case you have installed MPack in a non-standard directory, consider
using the

  --with-MPACK=/path/to/mpack

parameter to specify the location where to find MPack.
Otherwise, you can obtain the MPack library from

https://github.com/nakatamaho/mpack

**********************************************************************
        ])
    fi

fi

AM_CONDITIONAL(with_MPACK, 
        [test "$enable_mpack" = "yes"])
        
])

