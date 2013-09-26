# SYNOPSIS
#
#   Test for Planet
#
#   AX_PATH_PLANET( <Minimum Required Version>, <package-required=yes/no> )
#
# DESCRIPTION
#
#   Provides a --with-planet=DIR option. Searches --with-planet,
#   $PLANET_DIR, and the usual places for PLANET headers and libraries.
#
#   On success, sets PLANET_CPPFLAGS, PLANET_LIBS, and #defines HAVE_PLANET.
#   Also defines automake conditional PLANET_ENABLED.  Assumes package
#   is optional unless overridden with $2=yes.
#
# LAST MODIFICATION
#
#   $Id: planet.m4 37037 2013-02-16 01:03:09Z pbauman $
#
# COPYLEFT
#
#   Copyright (c) 2013 Roy H. Stogner <roystgnr@ices.utexas.edu>
#   Copyright (c) 2013 Paul T. Bauman <pbauman@ices.utexas.edu>
#   Copyright (c) 2010 Karl W. Schulz <karl@ices.utexas.edu>
#   Copyright (c) 2009 Rhys Ulerich <rhys.ulerich@gmail.com>
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2008 Caolan McNamara <caolan@skynet.ie>
#   Copyright (c) 2008 Alexandre Duret-Lutz <adl@gnu.org>
#   Copyright (c) 2008 Matthew Mueller <donut@azstarnet.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_PATH_PLANET],
[

AC_ARG_VAR(PLANET_DIR,[root directory of Planet installation])

AC_ARG_WITH(planet,
  [AS_HELP_STRING([--with-planet[=DIR]],[root directory of Planet installation (default = PLANET_DIR)])],
  [with_planet=$withval
if test "${with_planet}" != yes; then
    PLANET_PREFIX=$withval
fi
],[
with_planet=$withval
if test "x${PLANET_DIR}" != "x"; then
   PLANET_PREFIX=${PLANET_DIR}
fi
])

# package requirement; if not specified, the default is to assume that
# the package is optional

is_package_required=ifelse([$2], ,no, $2 )

HAVE_PLANET=0

# logic change: if the user called the macro, check for package,
# decide what to do based on whether the package is required or not.

# if test "${with_planet}" != no ; then

    if test -d "${PLANET_PREFIX}/lib" ; then
       PLANET_LDFLAGS="-L${PLANET_PREFIX}/lib -Wl,-rpath,${PLANET_PREFIX}/lib"
       PLANET_LIBS="-lplanet"
    fi

    if test -d "${PLANET_PREFIX}/include" ; then
       PLANET_CPPFLAGS="-I${PLANET_PREFIX}/include -I${PLANET_PREFIX}/src"
    fi

    ac_PLANET_save_CPPFLAGS="$CPPFLAGS"
    ac_PLANET_save_LDFLAGS="$LDFLAGS"
    ac_PLANET_save_LIBS="$LIBS"

    CPPFLAGS="${PLANET_CPPFLAGS} ${CPPFLAGS}"
    LDFLAGS="${PLANET_LDFLAGS} ${LDFLAGS}"
    LIBS="${PLANET_LIBS} ${LIBS}"

    AC_LANG_PUSH([C++])
    AC_CHECK_HEADER([planet/planet_version.h],[found_header=yes],[found_header=no])
    AC_LANG_POP([C++])

    #-----------------------
    # Minimum version check
    #----------------------

    min_planet_version=ifelse([$1], ,0.0.0, $1)

    # looking for major.minor.micro style versioning

    MAJOR_VER=`echo $min_planet_version | sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${MAJOR_VER}" = "x" ; then
       MAJOR_VER=0
    fi

    MINOR_VER=`echo $min_planet_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${MINOR_VER}" = "x" ; then
       MINOR_VER=0
    fi

    MICRO_VER=`echo $min_planet_version | sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${MICRO_VER}" = "x" ; then
       MICRO_VER=0
    fi

    # begin additional test(s) if header if available

    if test "x${found_header}" = "xyes" ; then

        AC_MSG_CHECKING(for planet - version >= $min_planet_version)
        version_succeeded=no

	AC_LANG_PUSH([C++])
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
        @%:@include "planet/planet_version.h"
            ]], [[
            #if PLANET_MAJOR_VERSION > $MAJOR_VER
            /* Sweet nibblets */
            #elif (PLANET_MAJOR_VERSION >= $MAJOR_VER) && (PLANET_MINOR_VERSION > $MINOR_VER)
            /* Winner winner, chicken dinner */
            #elif (PLANET_MAJOR_VERSION >= $MAJOR_VER) && (PLANET_MINOR_VERSION >= $MINOR_VER) && (PLANET_MICRO_VERSION >= $MICRO_VER)
            /* I feel like chicken tonight, like chicken tonight? */
            #else
            #  error version is too old
            #endif
        ]])],[
            AC_MSG_RESULT(yes)
            version_succeeded=yes
        ],[
            AC_MSG_RESULT(no)
        ])
	AC_LANG_POP([C++])

    if test "$version_succeeded" != "yes";then
       if test "$is_package_required" = yes; then
          AC_MSG_ERROR([

   Your Planet version does not meet the minimum versioning
   requirements ($min_planet_version).  Please use --with-planet to specify the location
   of an updated installation or consider upgrading the system version.

          ])
       fi
    fi

    # Library availability

    AC_MSG_CHECKING([for -lplanet linkage])

    AC_LANG_PUSH([C++])

    AC_LINK_IFELSE(
                  [AC_LANG_PROGRAM([#include "planet/planet_version.h"],
                                   [Planet::get_planet_version()])],
                  [AC_MSG_RESULT(yes)
                   found_library=yes],
                  [AC_MSG_RESULT(no) 
                   found_library=no])

    fi   dnl end test if header if available

    AC_LANG_POP([C++])

    CPPFLAGS="$ac_PLANET_save_CPPFLAGS"
    LDFLAGS="$ac_PLANET_save_LDFLAGS"
    LIBS="$ac_PLANET_save_LIBS"

    succeeded=no
    if test "$found_header" = yes; then
        if test "$version_succeeded" = yes; then
           if test "$found_library" = yes; then
              succeeded=yes
           fi
        fi
    fi

    if test "$succeeded" = no; then
       if test "$is_package_required" = yes; then
          AC_MSG_ERROR([Planet not found.  Try either --with-planet or setting PLANET_DIR.])
       else
          AC_MSG_NOTICE([optional Planet library not found])
          PLANET_CPPFLAGS=""   # PLANET_CFLAGS empty on failure
          PLANET_LDFLAGS=""    # PLANET_LDFLAGS empty on failure
          PLANET_LIBS=""       # PLANET_LIBS empty on failure
       fi
    else
        HAVE_PLANET=1
        AC_DEFINE(HAVE_PLANET,1,[Define if Planet is available])
        AC_SUBST(PLANET_CPPFLAGS)
        AC_SUBST(PLANET_LDFLAGS)
        AC_SUBST(PLANET_LIBS)
        AC_SUBST(PLANET_PREFIX)
    fi

    AC_SUBST(HAVE_PLANET)

# fi

AM_CONDITIONAL(PLANET_ENABLED,test x$HAVE_PLANET = x1)

])
