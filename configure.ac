AC_INIT(buff-em, 0.1, homer@homerreid.com)
AC_CONFIG_SRCDIR([README])

# Shared-library version.  This has a technical
# meaning related to binary compatibility and is NOT
# the same as the human version number.  See http://tinyurl.com/yegqrlz
SHARED_VERSION_INFO="0:0:0" # CURRENT:REVISION:AGE

AM_INIT_AUTOMAKE
LT_INIT
AM_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_MAINTAINER_MODE
AC_SUBST(SHARED_VERSION_INFO)
AC_DISABLE_SHARED dnl shared libs by default cause too much trouble

AC_PROG_CXX
AC_PROG_F77
AC_F77_WRAPPERS
AC_PROG_INSTALL

AC_LANG([C++])

##################################################
# compiler flags
##################################################

AC_ARG_ENABLE(debug, [AC_HELP_STRING([--enable-debug],[compile for debugging])], enable_debug=$enableval, enable_debug=no)

if test "x$enable_debug" = xyes; then
   CXXFLAGS=-g
   FFLAGS=-g
fi

# TODO: set optimization flags? Does it matter, or is the default -g -O2 okay?

# add gcc warnings in debug or maintainer mode
if test "x$enable_debug" = xyes || test "$USE_MAINTAINER_MODE" = yes; then
if test "$ac_test_CFLAGS" != "set"; then
   if test $ac_cv_prog_gcc = yes; then
      CXXFLAGS="$CXXFLAGS -Wall -W -Wcast-qual -Wpointer-arith -Wcast-align -pedantic -Wno-long-long -Wshadow -Wwrite-strings -Wredundant-decls -pedantic"
   fi
fi
fi

##################################################
# threading library
##################################################

AC_ARG_WITH(openmp, [AC_HELP_STRING([--without-openmp],[do not use OpenMP directives for parallelism])], with_openmp=$withval, with_openmp=yesdefault)

if test "x$with_openmp" = xyesdefault; then with_openmp=yes; fi
if test "x$with_openmp" = xyes; then
   AC_DEFINE([USE_OPENMP], [1], [Define to use OpenMP threading.])
   AX_OPENMP([], [AC_MSG_ERROR([Could not find OpenMP flag; configure with --without-openmp])])
   CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
fi

##################################################
# checks for blas/lapack
##################################################

ACX_BLAS([], [AC_MSG_ERROR([BLAS was not found!])])
ACX_LAPACK([], [AC_MSG_ERROR([LAPACK was not found!])])

LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

##################################################
# checks for readline 
# (which is used by some test programs)
##################################################
AC_CHECK_LIB(readline, readline)

##################################################
# checks for backtrace functionality, used to 
# output useful debugging info to log files in
# the event of code crashes
##################################################
AC_CHECK_HEADERS([execinfo.h])
AC_CHECK_FUNCS([backtrace])

##################################################
# Checks for miscellaneous math functions
##################################################
AC_CHECK_FUNCS([isnan isinf erf erfc copysign])

##################################################
# Checks for scuff-em
##################################################
AC_LANG_PUSH([C++])

AC_CHECK_LIB(scuff, cevaluator_create, [],
             [AC_MSG_ERROR([Could not find SCUFF-EM libraries; add -L/path/to/scuff-installation/libs to LDFLAGS.])])

AC_CHECK_HEADERS([libscuff.h], [],
             [AC_MSG_ERROR([Could not find SCUFF-EM headers; add -I/path/to/scuff-installation/include/scuff-em to CPPFLAGS.])])

AC_LANG_POP

##################################################
##################################################
##################################################
AC_CONFIG_FILES([
 Makefile
 m4/Makefile
 src/Makefile
 src/libs/Makefile
 src/libs/libbuff/Makefile
 src/libs/libbuff/UpsilonFiles/Makefile
 src/applications/Makefile
 src/applications/buff-analyze/Makefile
 src/applications/buff-scatter/Makefile
 src/unitTests/Makefile
])
AC_OUTPUT

##################################################
##################################################
##################################################
AC_MSG_NOTICE([**])
AC_MSG_NOTICE([** The configure process was successful.])
AC_MSG_NOTICE([** Now type 'make' to compile and sudo 'make install'.])
AC_MSG_NOTICE([**])     
