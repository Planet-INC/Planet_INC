BUILT_SOURCES   = .license.stamp

EXTRA_DIST = input_files

check_PROGRAMS  = 
check_PROGRAMS += atmosphere_tests
check_PROGRAMS += photon_tests
check_PROGRAMS += chapman_tests

AM_CPPFLAGS  = 
AM_CPPFLAGS += -I$(top_srcdir)/src/core/include
AM_CPPFLAGS += -I$(top_srcdir)/src/atmosphere/include
AM_CPPFLAGS += -I$(top_srcdir)/src/photon_flux/include
AM_CPPFLAGS += -I$(top_srcdir)/src/absorption/include
AM_CPPFLAGS += -I$(top_srcdir)/src/sampling/include
AM_CPPFLAGS += -I$(top_srcdir)/src/diffusion/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utilities/include
AM_CPPFLAGS += -I$(top_builddir)/src/utilities/include #planet_version.h
AM_CPPFLAGS +=  $(planet_optional_test_INCLUDES)

AM_LDFLAGS = $(planet_optional_test_LDFLAGS)
LIBS = $(planet_optional_test_LIBS)

LDADD = $(top_builddir)/src/libplanet.la

# Sources for these tests
atmosphere_tests_SOURCES = atmosphere_tests.C
photon_tests_SOURCES = photon_tests.C
chapman_tests_SOURCES = chapman_tests.C

#Define tests to actually be run
TESTS = 
TESTS  += atmosphere_tests
TESTS  += photon_tests
TESTS  += chapman_tests

CLEANFILES =
if CODE_COVERAGE_ENABLED
  CLEANFILES += *.gcda *.gcno
endif

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE \
	$(top_srcdir)/test/*.C
	echo 'updated source license headers' >$@

CLEANFILES += .license.stamp

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
