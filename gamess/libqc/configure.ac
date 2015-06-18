#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(libqc, 0.1, asadchev@gmail.com)
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([-Wall foreign])
AC_CONFIG_SRCDIR([src/core/molecule.hpp])
AC_CONFIG_HEADER([config/config.h])


LT_INIT([disable-shared])
AC_CONFIG_MACRO_DIR([config/m4])

# store current user given compiler flags to avoid default setup via AC_PROG_CXX
OLD_CXXFLAGS=$CXXFLAGS
OLD_CFLAGS=$CFLAGS

# check for C++ preprocessor and compiler and the library compiler
# (might change the compiler flags)
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_CXXCPP

# reset compiler flags to initial flags
CXXFLAGS=$OLD_CXXFLAGS
CFLAGS=$OLD_CFLAGS



# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdint.h stdlib.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_INT8_T
AC_TYPE_SIZE_T
AC_HEADER_TIME

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_CHECK_FUNCS([gettimeofday sqrt])

# VTK
sinclude(config/am_path_vtk.m4)
AM_OPTIONS_VTK
AM_PATH_VTK
AM_CONDITIONAL(HAVE_VTK, test "$with_vtk" = "yes")
if test "$with_rysq" != "no"; then
   AC_DEFINE([HAVE_VTK], [1], [with VTK])   
fi

# DOXYGEN SUPPORT
sinclude(config/ax_prog_doxygen.m4)

# doxygen features
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)

# initialize doxygen
DX_INIT_DOXYGEN($PACKAGE_NAME, [doc/Doxyfile], [doc])

# boost
sinclude(config/ax_boost_base.m4)
AX_BOOST_BASE([])
AM_CPPFLAGS="$AM_CPPFLAGS $BOOST_CPPFLAGS"

AC_CONFIG_FILES([Makefile
		src/Makefile
		src/core/Makefile
		src/hf/Makefile
		src/util/Makefile
		src/adapter/Makefile
		src/matrix/Makefile
		src/bindings/Makefile
		src/contrib/Makefile])

# debug compilation
AC_ARG_ENABLE(debug,
    AC_HELP_STRING(--enable-debug, [Debug compilation (Default = no)]),
    enable_debug=$enableval, enable_debug=no)

if test "$enable_debug" = "yes" ; then
    CFLAGS="$CFLAGS  -g -O0 -Wall -Wno-uninitialized"
    CXXFLAGS="$CXXFLAGS -g -O0 -Wall -Wno-uninitialized"
fi


## enable ga
AC_ARG_WITH(ga, AC_HELP_STRING([--with-ga=PATH], [where ga is installed]),
            [with_ga=$withval], [with_ga="no"])
AM_CONDITIONAL(_HAVE_GA, [test "$with_ga" != "no"])
if test "$with_ga" != "no"; then
   if test "$with_ga" != "yes"; then
       AM_CPPFLAGS="$AM_CPPFLAGS -I$with_ga/include"
       LDFLAGS="$LDFLAGS -L$with_ga/lib"
   fi
   #LIBS="$LIBS -lglobal -lga++"
   with_ga=1
   #AC_DEFINE(_HAVE_GA, 1, [enable ga])
else
   with_ga=0
fi
AC_SUBST(_HAVE_GA, ["$with_ga"])

## enable cuda
AC_ARG_WITH(cuda, AC_HELP_STRING([--with-cuda=PATH], [where cuda is installed]),
            [with_cuda=$withval], [with_cuda="no"])

# Rysq settings
AC_ARG_WITH(rysq, 
            AC_HELP_STRING([--with-rysq=PATH], [where rysq is installed]))

if test "$with_rysq" != "no"; then
   if test "$with_rysq" != "yes"; then
      AM_CPPFLAGS="$AM_CPPFLAGS -I$with_rysq/include"
      LDFLAGS="$LDFLAGS -L$with_rysq/lib"
   fi
   LIBS="$LIBS -lrysq"
fi  


# gamess bindings
AC_ARG_WITH(gamess,
	AS_HELP_STRING([--with-gamess],[Adds GAMESS support to libqc]),
	[if test "$withval" = "yes"; then
	    AC_DEFINE([LIBQC_WITH_GAMESS], [1], [Compile with GAMESS interface])
	 fi])
AM_CONDITIONAL(LIBQC_WITH_GAMESS, test "$with_gamess" = "yes")

# integer8 Fortran bindings
AC_ARG_WITH(integer8,
	AS_HELP_STRING([--with-integer8],[Adds Integer8 Fortran support]),
	[if test x$withval = xyes; then
	    AC_DEFINE([LIBQC_WITH_INTEGER8], [1], [8-byte Fortran integer])
	 fi])

		

dnl AC_ARG_ENABLE(lmax, AS_HELP_STRING([--enable-lmax[=lmax]],[Max angular momentum]), [RYSQ_LMAX=$enableval],[RYSQ_LMAX=3])
dnl AC_SUBST([RYSQ_LMAX])

dnl AC_ARG_ENABLE(sp, AS_HELP_STRING([--disable-sp],[disables SP orbitals]),[
dnl if test x$enableval = xno; then
dnl    RYSQ_SP=0
dnl else
dnl    RYSQ_SP=1
dnl fi],[RYSQ_SP=1])
dnl AC_SUBST([RYSQ_SP])

dnl AC_ARG_ENABLE(unroll-bra, AS_HELP_STRING([--enable-unroll-bra[=unroll_size]],[Sets the size of Unrolled Bra integral loop]),
dnl 			  [RYSQ_UNROLL_BRA=$enableval],[RYSQ_UNROLL_BRA=10])
dnl AC_SUBST([RYSQ_UNROLL_BRA])

dnl AC_ARG_ENABLE(cuda, AS_HELP_STRING([--enable-cuda],[Enables CUDA code to be compiled]), [
dnl if test x$enableval = xyes; then
dnl    CURYSQ=1
dnl else
dnl 	CURSYQ=0
dnl fi],[CURYSQ=0])
dnl AM_CONDITIONAL(CURYSQ, test x$CURYSQ = x1)

dnl AC_ARG_ENABLE(cu-lmax, AS_HELP_STRING([--enable-cu-lmax[=lmax]],[Max angular momentum for CUDA]), [CURYSQ_LMAX=$enableval],[CURYSQ_LMAX=3])
dnl AC_SUBST([CURYSQ_LMAX])

dnl AC_ARG_WITH(integer8, AS_HELP_STRING([--with-integer8],[With Fortran integer8.]), [
dnl  if test x$withval = xyes; then
dnl     libqc_with_integer8=1
dnl  else
dnl 	libqc_with_integer8=0
dnl  fi], [libqc_with_integer8=0])
dnl AM_CONDITIONAL(RYSQ_WITH_GAMESS, test x$rysq_with_gamess = x1)

#AM_CPPFLAGS="$AM_CPPFLAGS -I $top_srcdir"

     dnl AC_INIT([arm], [1.0])
     dnl AC_CONFIG_AUX_DIR([.])
     dnl AM_INIT_AUTOMAKE
     dnl AC_PROG_CC
     dnl AC_CONFIG_FILES([Makefile])
     # Call hand's ./configure script recursively.
AC_CONFIG_SUBDIRS([rysq])
AM_CPPFLAGS="$AM_CPPFLAGS -I $PWD/rysq/src"
     dnl AC_OUTPUT


AC_SUBST([LIBS])
AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AM_LDFLAGS])

AC_OUTPUT
