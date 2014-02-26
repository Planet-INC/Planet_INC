 dnl                                               -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([planet], [0.0.0], [splessis@ices.utexas.edu], [planet])
PACKAGE_DESCRIPTION="A planet's atmospheric code using Antioch, adapted to Titan"
AC_SUBST([PACKAGE_DESCRIPTION])
PACKAGE_URL=""
AC_SUBST([PACKAGE_URL])
AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_HEADER(planet_config.h.tmp)
dnl Prefixes everything in the temporary file with FINS_
AX_PREFIX_CONFIG_H(planet_config.h,PLANET,planet_config.h.tmp)

AC_CONFIG_AUX_DIR([build-aux])
AC_CANONICAL_TARGET()
AM_INIT_AUTOMAKE([dist-bzip2 tar-ustar 1.11 no-define color-tests parallel-tests])
AM_SILENT_RULES(yes)  # use silent rules where available - automake 1.11

LT_INIT

dnl Release versioning
AX_SPLIT_VERSION
GENERIC_MAJOR_VERSION=$AX_MAJOR_VERSION
GENERIC_MINOR_VERSION=$AX_MINOR_VERSION
GENERIC_MICRO_VERSION=$AX_POINT_VERSION
GENERIC_LIBRARY_NAME=planet

dnl shared library versioning

PACKAGE=$GENERIC_LIBRARY_NAME
GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION
GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION

dnl AC_SUBST(GENERIC_LIBRARY_VERSION)
dnl AC_SUBST(GENERIC_LIBRARY_NAME)

AC_SUBST(GENERIC_VERSION)
AC_SUBST(GENERIC_RELEASE)
AC_SUBST(GENERIC_MAJOR_VERSION)
AC_SUBST(GENERIC_MINOR_VERSION)
AC_SUBST(GENERIC_MICRO_VERSION)

AC_DEFINE_UNQUOTED([MAJOR_VERSION],$GENERIC_MAJOR_VERSION,[Major version])
AC_DEFINE_UNQUOTED([MINOR_VERSION],$GENERIC_MINOR_VERSION,[Minor version])
AC_DEFINE_UNQUOTED([MICRO_VERSION],$GENERIC_MICRO_VERSION,[Micro version])

dnl------------------------------
dnl Checks for compilers
dnl------------------------------

AC_PROG_CC
AC_PROG_CXX

dnl Added for AX_CODE_COVERAGE macro to work correctly, even though there's no Fortran here.
AC_PROG_FC

dnl Allow users to make vectorized tests larger
AC_ARG_WITH([n_tuples],
	    AC_HELP_STRING([--with-n_tuples=<1...1000000000>],
			   [number of tuples used in each vector test]),
	    [n_tuples="$withval"],
	    [n_tuples=2])
AC_DEFINE_UNQUOTED(N_TUPLES, $n_tuples, [number of tuples used in each vector test])

dnl Allow users to make vectorized tests longer
AC_ARG_WITH([n_runs],
	    AC_HELP_STRING([--with-n_runs=<1...1000000000>],
			   [number of runs done in each vector test]),
	    [n_runs="$withval"],
	    [n_runs=1])
AC_DEFINE_UNQUOTED(N_RUNS, $n_runs, [number of runs done in each vector test])

dnl------------------------
dnl Checks dependencies
dnl--------------------------

# Check for Antioch
AX_PATH_ANTIOCH(0.0.0,yes)

AX_PATH_GRINS(0.4.0,yes)

AC_ARG_ENABLE([libmesh-flags],
               AC_HELP_STRING([--disable-libmesh-flags],[disable use of libmesh cxxflags for CXXFLAGS]),
               [case "${enableval}" in
                   yes) USE_LIBMESH_CXX_FLAGS=1 ;;
                    no) USE_LIBMESH_CXX_FLAGS=0 ;;
                     *) AC_MSG_ERROR(bad value ${enableval} for libmesh-flags) ;;
                esac],
               [USE_LIBMESH_CXX_FLAGS=1])

if test "$USE_LIBMESH_CXX_FLAGS" = "1"; then
   AC_DEFINE(USE_LIBMESH_CXX_FLAGS,1,[Define if libmesh cxxflags enabled])
   CXXFLAGS=`$LIBMESH_CONFIG --cppflags --cxxflags`
else
 dnl TODO: WE NEED TO DISTINGUISH BETWEEN COMPILERS IN A BETTER WAY
   AC_DEFINE(USE_LIBMESH_CXX_FLAGS,0,[Define if libmesh cxxflags enabled])
   CXXFLAGS="-O0 -g -fno-unsafe-math-optimizations"
fi

dnl -Wall warnings, -Wall the time.
AX_CXXFLAGS_WARN_ALL

dnl--------------------------
dnl Even if no dependency needed it,
dnl C++11 makes our asserts nicer,
dnl so let's try and enable it.
dnl--------------------------
AC_ARG_ENABLE(cxx11,
              AC_HELP_STRING([--enable-cxx11],
                             [build with C++11 support]),
	      [case "${enableval}" in
	        yes)  enablecxx11=yes ;;
	         no)  enablecxx11=no ;;
	          *)  AC_MSG_ERROR(bad value ${enableval} for --enable-cxx11) ;;
	       esac],
	       [enablecxx11=optional])

if (test x$enablecxx11 = xyes); then
  AX_CXX_COMPILE_STDCXX_11(noext, mandatory)
  AX_CXX_AUTO_THIS(optional)
elif (test x$enablecxx11 = xoptional); then
  AX_CXX_COMPILE_STDCXX_11(noext, optional)
  AX_CXX_AUTO_THIS(optional)
fi

CONFIGURE_EIGEN(2.3.0,yes)
AX_PATH_GRINS(0.4.0,yes)

dnl--------------------------
dnl Checks for code coverage
dnl--------------------------
AX_CODE_COVERAGE

dnl---------------------------------
dnl Query configuration environment
dnl Must be done *before* Doxygen to
dnl get some happy variables defined.
dnl---------------------------------
AX_SUMMARIZE_ENV

dnl------------------------------------------
dnl            Doxygen support
dnl------------------------------------------

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN(planet, doxygen/planet.dox, docs/doxygen)

AC_CONFIG_FILES(doxygen/txt_common/about_vpath.page)

dnl-----------------------------------------------
dnl Generate files for unit and regression testing
dnl-----------------------------------------------
AC_CONFIG_FILES(test/temperature_unit.sh,                     [chmod +x test/temperature_unit.sh])
AC_CONFIG_FILES(test/photon_opacity_unit.sh,                  [chmod +x test/photon_opacity_unit.sh])
AC_CONFIG_FILES(test/atmospheric_mixture_unit.sh,             [chmod +x test/atmospheric_mixture_unit.sh])
AC_CONFIG_FILES(test/photon_evaluator_unit.sh,                [chmod +x test/photon_evaluator_unit.sh])
AC_CONFIG_FILES(test/eddy_diffusion_evaluator_unit.sh,        [chmod +x test/eddy_diffusion_evaluator_unit.sh])
AC_CONFIG_FILES(test/molecular_diffusion_evaluator_unit.sh,   [chmod +x test/molecular_diffusion_evaluator_unit.sh])
AC_CONFIG_FILES(test/diffusion_evaluator_unit.sh,             [chmod +x test/diffusion_evaluator_unit.sh])
AC_CONFIG_FILES(test/physics_helper_unit.sh,                  [chmod +x test/physics_helper_unit.sh])
AC_CONFIG_FILES(test/solver_test.sh,                          [chmod +x test/solver_test.sh])
AC_CONFIG_FILES(test/ionospheric_test.sh,                     [chmod +x test/ionospheric_test.sh])

AC_CONFIG_FILES(test/input/solver_test.in)
AC_CONFIG_FILES(test/input/grins_input_physics_helper.in)

dnl-----------------------------------------------
dnl Generate header files
dnl-----------------------------------------------
AC_CONFIG_FILES(src/utilities/include/planet/planet_version.h)

dnl-----------------------------------------------
dnl Generate Makefiles and other templated things
dnl-----------------------------------------------
AC_CONFIG_FILES([
  planet.pc
  doxygen/Makefile
  doxygen/planet.dox
  Makefile
  src/Makefile
  test/Makefile
])

# Must still call AC_OUTPUT() after generating all the files
AC_OUTPUT()

dnl------------------------------
dnl Final summary
dnl------------------------------
AX_SUMMARIZE_CONFIG
