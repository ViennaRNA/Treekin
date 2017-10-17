dnl ---------------------------------------------------------------------------
dnl Autoconf initialization.
dnl
dnl Bug reports are sent to mtw@tbi.univie.ac.at
dnl ---------------------------------------------------------------------------

AC_INIT([Treekin], 0.4.1, [mtw@tbi.univie.ac.at], [Treekin], [http://www.tbi.univie.ac.at/RNA/Treekin])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(calc.cpp)
AC_CONFIG_HEADERS(config.h)

AM_INIT_AUTOMAKE

AC_CANONICAL_HOST

TestDataDirectory="\$(srcdir)"
dnl SHAREDLIBEXT="so"

case $host_os in
  *linux*)
    HOST_TYPE="linux"
    SHAREDLIBEXT="so"
    ;;

  *solaris*)
    HOST_TYPE="solaris"
    SHAREDLIBEXT="so"
    ;;

  *cygwin*)
    HOST_TYPE="cygwin"
    SHAREDLIBEXT="so"
    ;;

  *darwin*)
    HOST_TYPE="darwin"
    SHAREDLIBEXT="dylib"
    TestDataDirectory="."
    ;;
esac
AC_SUBST(TestDataDirectory)
AC_SUBST(HOST_TYPE)
AC_SUBST(SHAREDLIBEXT)
AC_SUBST(CFLAGS)

dnl ---------------------------------------------------------------------------
dnl Checks for various programs and packages.
dnl ---------------------------------------------------------------------------

AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_PROG_MAKE_SET

AC_PATH_PROG([AUTOCONF],[autoconf],[autoconf])
AC_PATH_PROG([ACLOCAL],[aclocal],[aclocal])
AC_PATH_PROG([MAKE],[make],[make])
AC_PATH_PROG([AR],[ar],[ar])

dnl ---------------------------------------------------------------------------
dnl Checks for gengetopt
dnl ---------------------------------------------------------------------------
AC_DEFUN([AX_TREEKIN_GENGETOPT],[
  AC_PATH_PROG([GENGETOPT],[gengetopt],[],[])
  # if gengetopt is unavailable, check whether commandline parser sources are present
  AS_IF([ test "x$GENGETOPT" = "x"],[
    create_commandline_parser="no"
    AS_IF([ ! test -f treekin_cmdline.c ], [
      AC_MSG_ERROR(gengetopt not found)
    ])
  ])
  AM_CONDITIONAL(MAKE_CMDLINE_PARSER, test "x$create_commandline_parser" != "xno")
])

AX_TREEKIN_GENGETOPT



dnl ---------------------------------------------------------------------------
dnl Checks for libraries.
dnl ---------------------------------------------------------------------------

AC_LANG(C)
#AC_CHECK_LIB(gfortran, exp)
dnl AC_SEARCH_LIBS([ddot_],[blas],[],[
dnl		AC_MSG_ERROR([unable to find the ddot_() function in lblas])
dnl		])
AC_CHECK_LIB(blas,ddot_,[],[
	AC_MSG_ERROR([unable to find the ddot_() function in lblas. Is blas installed?])
])				
AC_CHECK_LIB(lapack,dgeevx_,[],[
	AC_MSG_ERROR([unable to find the dgeevx_() function in lblas. Is lapack installed?])
])

#search for the correct include path in the default locations
AC_CHECK_HEADERS([/usr/include/lapacke.h], [
  CXXFLAGS="-I/usr/include $CXXFLAGS"
], [
  AC_CHECK_HEADERS([/usr/include/lapacke/lapacke.h], [
    CXXFLAGS="-I/usr/include/lapacke $CXXFLAGS"
  ], [
    AC_CHECK_HEADERS([/usr/local/include/lapacke.h], [
      CXXFLAGS="-I/usr/local/include $CXXFLAGS"
    ], [
      AC_CHECK_HEADERS([/usr/local/include/lapacke/lapacke.h], [
        CXXFLAGS="-I/usr/local/include/lapacke $CXXFLAGS"
      ], [])
    ])
  ])
])

dnl AC_CHECK_LIB(m, isnan)
dnl AC_CHECK_LIB(arpack,dseupd_)

dnl ---------------------------------------------------------------------------
dnl Checks for header files.
dnl ---------------------------------------------------------------------------

AC_HEADER_STDC
AC_CHECK_HEADERS(errno.h)

dnl Check for glib-2.0 using pkgconfig
dnl PKG_CHECK_MODULES([GLIB], [glib-2.0], [have_glib=true], [have_glib=false])
dnl if test "x${have_glib}" = "xfalse" ; then
dnl         AC_MSG_ERROR([No Glib package information found])
dnl fi
dnl AC_SUBST(GLIB_CFLAGS)
dnl AC_SUBST(GLIB_LIBS)

AC_CHECK_HEADERS(math.h)

dnl ---------------------------------------------------------------------------
dnl Check for mlapack
dnl ---------------------------------------------------------------------------

DEPENDENCYNOTFOUND=0;

AC_ARG_ENABLE(MPACK,
    AC_HELP_STRING(
    	[--disable-MPACK],
	    [disable support for MPACK models (default=enabled)]
	),
[enable_mpack=$enableval],
[enable_mpack=yes]
)

AM_CONDITIONAL(with_MPACK, 
	[test "$enable_mpack" = "yes"])
	
# MPACK package library path support, if not installed in usual directories
AC_ARG_WITH(MPACK,
	AC_HELP_STRING(
	    [--with-MPACK=PREFIX],
	    [alternative prefix path to MPACK library]
	  ),
	  MPACKPATHSET=1,
	  MPACKPATHSET=0
)
if test  $MPACKPATHSET = 1 ; then
		CXXFLAGS="-I$with_MPACK/include $CXXFLAGS"
		LDFLAGS="-L$with_MPACK/lib $LDFLAGS"
fi
mpack_LIBS=""

NEED_LIB_QD=0
NEED_LIB_GMP=0
NEED_LIB_MPC=0
NEED_LIB_MPFR=0

if test "$enable_mpack" = "yes"; then
  if test  $MPACKPATHSET = 1 ; then
        AC_LANG_PUSH([C++])
		AC_CHECK_HEADER(mpack/mpack_config.h,
          [
            MPACKNOTFOUND=0;
          ],
          [
            MPACKNOTFOUND=1;
          ],
           [#if HAVE_MPACK_MPACK_CONFIG_H
           #include <mpack/mpack_config.h>
           #endif
        ])
		AC_LANG_POP
  else
		PKG_CHECK_MODULES([mpack], [mpack >= 0.8], 
		[
		  MPACKNOTFOUND=0;
		  CXXFLAGS="$mpack_CFLAGS $CXXFLAGS"
		  # PKG_CHECK_MODULES sets mpack_LIBS automatically
		], [MPACKNOTFOUND=1;])
  fi
  
  ppFlags=" "
  #test which libraries of the mpack package are installed and set the flags.
  AC_LANG_PUSH([C++])
  AC_CHECK_HEADER(mpack/mlapack_qd.h, [ppFlags="-DWITH_MPACK_QD $ppFlags"; NEED_LIB_QD=1 ],[],
  [#if HAVE_MPACK_MPACK_CONFIG_H
   #include <mpack/mpack_config.h>
   #endif
   ])
  AC_CHECK_HEADER(mpack/mlapack_dd.h, [ppFlags="-DWITH_MPACK_DD $ppFlags"; NEED_LIB_QD=1 ],[],
  [#if HAVE_MPACK_MPACK_CONFIG_H
   #include <mpack/mpack_config.h>
   #endif
   ])
  AC_CHECK_HEADER(mpack/mlapack_double.h, [ppFlags="-DWITH_MPACK_DOUBLE $ppFlags" ],[],
  [#if HAVE_MPACK_MPACK_CONFIG_H
   #include <mpack/mpack_config.h>
   #endif
   ])
  AC_CHECK_HEADER(mpack/mlapack_longdouble.h, [ppFlags="-DWITH_MPACK_LD $ppFlags" ],[],
  [#if HAVE_MPACK_MPACK_CONFIG_H
   #include <mpack/mpack_config.h>
   #endif
   ])
  AC_CHECK_HEADER(mpack/mlapack___float128.h, [ppFlags="-DWITH_MPACK___FLOAT128 $ppFlags" ],[],
  [#if HAVE_MPACK_MPACK_CONFIG_H
   #include <mpack/mpack_config.h>
   #endif
   ])
  AC_CHECK_HEADER(mpack/mlapack_mpfr.h, [ppFlags="-DWITH_MPACK_MPFR $ppFlags"; NEED_LIB_MPFR=1; NEED_LIB_MPC=1 ],[],
  [#if HAVE_MPACK_MPACK_CONFIG_H
   #include <mpack/mpack_config.h>
   #endif
   ])
  AC_CHECK_HEADER(mpack/mlapack_gmp.h, [ppFlags="-DWITH_MPACK_GMP $ppFlags"; NEED_LIB_GMP=1 ],[],
  [#if HAVE_MPACK_MPACK_CONFIG_H
   #include <mpack/mpack_config.h>
   #endif
   ])
  AC_LANG_POP
  
  CPPFLAGS=$CPPFLAGS$ppFlags
fi

if test "$enable_mpack" = "yes"; then
	CPPFLAGS=" -DWITH_MPACK $CPPFLAGS"
fi

# error output if MPACK not found
if test "$MPACKNOTFOUND" = "1"; then
	AC_MSG_NOTICE()
	AC_MSG_NOTICE(the MPACK library is required for building Treekin.)
	AC_MSG_NOTICE(it can be obtained from https://github.com/nakatamaho/mpack/archive/master.zip.)
	AC_MSG_NOTICE()
 	if test "$MPACKPATHSET" = "1"; then
		AC_MSG_NOTICE(cant find MPACK library in given path '$with_MPACK'.)
 	else
		AC_MSG_NOTICE(if MPACK is installed elsewhere use --with-MPACK=PREFIX.)
 	fi
 	DEPENDENCYNOTFOUND=1;
fi
# error ABORT if on of the libraries was not found
if test "$DEPENDENCYNOTFOUND" = "1"; then
	AC_MSG_NOTICE()
	AC_MSG_ERROR(library not found! Try './configure --help'.)
fi


# check for further dependencies of MPACK
AC_LANG_PUSH([C++])
if test "$NEED_LIB_GMP" = "1"; then
    AC_CHECK_HEADER(gmp.h, HAVE_GMP_HEADER=yes,[],
    [#if HAVE_SUPPRESS_COMPILER_WARNING_H
     #include <SUPPRESS_COMPILER_WARNING.h>
     #endif
    ])
    if test x$HAVE_GMP_HEADER = xyes; then
        AC_CHECK_LIB(gmp, __gmpz_init, [mpack_LIBS="$mpack_LIBS -lgmpxx -lgmp"; HAVE_GMP=yes])
    fi
    if test x$HAVE_GMP != xyes; then
        AC_MSG_RESULT([No GMP library with C++ wrapper found! You could try to set the CPPFLAGS=-I/path/to/gmp/include and LDFLAGS=-L/path/to/gmp/lib])
    fi
fi
if test "$NEED_LIB_MPFR" = "1"; then
    AC_CHECK_HEADER(mpfr.h, HAVE_MPFR_HEADER=yes,[],
    [#if HAVE_SUPPRESS_COMPILER_WARNING_H
     #include <SUPPRESS_COMPILER_WARNING.h>
     #endif
    ])
    if test x$HAVE_MPFR_HEADER = xyes; then
        AC_CHECK_LIB(mpfr, mpfr_init, [mpack_LIBS="$mpack_LIBS -lmpfr"; HAVE_MPFR=yes])
    fi
    if test x$HAVE_MPFR != xyes; then
        AC_MSG_RESULT([No MPFR library found! You could try to set the CPPFLAGS=-I/path/to/mpfr/include and LDFLAGS=-L/path/to/mpfr/lib])
    fi
fi
if test "$NEED_LIB_MPC" = "1"; then
    AC_CHECK_HEADER(mpc.h, HAVE_MPC_HEADER=yes,[],
    [#if HAVE_SUPPRESS_COMPILER_WARNING_H
     #include <SUPPRESS_COMPILER_WARNING.h>
     #endif
    ])
    if test x$HAVE_MPC_HEADER = xyes; then
        AC_CHECK_LIB(mpc, mpc_init2, [mpack_LIBS="$mpack_LIBS -lmpc"; HAVE_MPC=yes])
    fi
    if test x$HAVE_MPC != xyes; then
        AC_MSG_RESULT([No MPC library found! You could try to set the CPPFLAGS=-I/path/to/mpc/include and LDFLAGS=-L/path/to/mpc/lib])
    fi
fi
if test "$NEED_LIB_QD" = "1"; then
    AC_CHECK_HEADER(qd/qd_real.h, HAVE_QD_HEADER=yes,[],
    [#if HAVE_SUPPRESS_COMPILER_WARNING_H
     #include <SUPPRESS_COMPILER_WARNING.h>
     #endif
    ])
    if test x$HAVE_QD_HEADER = xyes; then
        AC_CHECK_LIB(qd, c_qd_sqrt, [mpack_LIBS="$mpack_LIBS -lqd"; HAVE_QD=yes])
    fi
    if test x$HAVE_QD != xyes; then
        AC_MSG_RESULT([No QD library found! You could try to set the CPPFLAGS=-I/path/to/qd/include and LDFLAGS=-L/path/to/qd/lib])
    fi
fi
AC_LANG_POP


dnl ---------------------------------------------------------------------------
dnl Checks for typedefs, structures, and compiler characteristics.
dnl ---------------------------------------------------------------------------

AC_C_BIGENDIAN
AC_C_CONST
AC_TYPE_SIZE_T
AC_EXEEXT
AC_OBJEXT

dnl ---------------------------------------------------------------------------
dnl Miscellaneous
dnl ---------------------------------------------------------------------------

# Root directory
# Translate path for native Windows compilers for use with 'make check'

ROOT_DIR=`pwd`
case $host in
*-*-cygwin* | *-*-mingw*)
  if (cygpath --mixed $ROOT_DIR) >/dev/null 2>/dev/null; then
    ROOT_DIR=`cygpath --mixed $ROOT_DIR`
  fi
  # Extra files generated by some Windows compilers
  EXTRA_CLEAN="*.stackdump *.exp *.lib"
  ;;
esac

dnl ---------------------------------------------------------------------------
dnl Output
dnl ---------------------------------------------------------------------------

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([packaging/treekin.spec packaging/PKGBUILD])

AC_OUTPUT

dnl ---------------------------------------------------------------------------
dnl Inform users of settings
dnl ---------------------------------------------------------------------------

AC_MSG_NOTICE(
[

You can run make now.
])
