# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

AM_TESTSUITE_SUMMARY_HEADER = ' for $(PACKAGE_STRING) ($(subdir))'
if ENABLE_OPENMP
AM_TESTS_ENVIRONMENT = OMP_NUM_THREADS=$${OMP_NUM_THREADS-1}; export OMP_NUM_THREADS;
endif

TEST_EXTENSIONS = .x .sh

test_PROGRAMS_ = \
    test_grid2vtk.x \
    test_vtk_output.x

xfail_test_PROGRAMS_ =

# Build test programs:
check_PROGRAMS = $(test_PROGRAMS_) $(xfail_test_PROGRAMS_)

test_SCRIPTS_ = \
    test_duplicate_stencils_parallel.sh \
    test_duplicate_stencils_parallel_c.sh \
    test_generate_cubed_sphere.sh \
    test_read_cube_csv.sh \
    test_read_fesom.sh \
    test_read_exodus_parallel.sh \
    test_read_icon.sh \
    test_read_icon_c.sh \
    test_read_icon_parallel.sh \
    test_read_icon_parallel_c.sh \
    test_read_mpiom.sh \
    test_read_scrip.sh \
    test_read_scrip_c.sh \
    test_read_scrip_parallel.sh \
    test_read_scrip_parallel_c.sh \
    test_read_woa_data.sh

xfail_test_SCRIPTS_ =

# Build programs associated with test scripts:
check_PROGRAMS += \
    test_read_cube_csv.x

if HAVE_NETCDF

check_PROGRAMS += \
    test_generate_cubed_sphere.x \
    test_read_fesom.x \
    test_read_icon_c.x \
    test_read_mpiom.x \
    test_read_scrip_c.x \
    test_read_woa_data.x

if ENABLE_FORTRAN_BINDINGS
check_PROGRAMS += \
    test_read_icon.x \
    test_read_scrip.x
endif ENABLE_FORTRAN_BINDINGS

if TEST_MPI

check_PROGRAMS += \
    test_read_exodus_parallel.x \
    test_read_icon_parallel_c.x \
    test_read_scrip_parallel_c.x \
    test_duplicate_stencils_parallel_c.x

if ENABLE_FORTRAN_BINDINGS
check_PROGRAMS += \
    test_read_icon_parallel.x \
    test_read_scrip_parallel.x \
    test_duplicate_stencils_parallel.x
endif ENABLE_FORTRAN_BINDINGS

endif TEST_MPI

endif HAVE_NETCDF

TESTS = \
    $(test_PROGRAMS_) \
    $(test_SCRIPTS_) \
    $(xfail_test_PROGRAMS_) \
    $(xfail_test_SCRIPTS_)

XFAIL_TESTS = \
    $(xfail_test_PROGRAMS_) \
    $(xfail_test_SCRIPTS_)

if TEST_MPI
if !ENABLE_CONCURRENT_MPI_TESTS
# Run MPI tests one by one:
#   The following dependencies are generated with
#   $(top_srcdir)/config/mpi_test_deps.sh.
test_read_exodus_parallel.log: test_read_scrip_parallel.log
test_read_scrip_parallel.log: test_read_scrip_parallel_c.log
test_read_scrip_parallel_c.log: test_read_icon_parallel.log
test_read_icon_parallel.log: test_read_icon_parallel_c.log
test_read_icon_parallel_c.log: test_duplicate_stencils_parallel.log
test_duplicate_stencils_parallel.log: test_duplicate_stencils_parallel_c.log
endif
endif

AM_CPPFLAGS = \
    -I$(top_srcdir)/src/core \
    -I$(top_srcdir)/src/utils \
    -I$(top_srcdir)/tests/common

AM_CFLAGS = \
    $(OPENMP_CFLAGS) \
    $(MPI_CFLAGS) \
    $(YAXT_CFLAGS) \
    $(NETCDF_CFLAGS)

AM_FCFLAGS = \
    $(FCPPINC)$(top_srcdir)/tests/common \
    $(FCMODINC)$(top_builddir)/tests/common \
    $(FCMODINC)$(top_builddir)/src/core \
    $(FCMODINC)$(top_builddir)/src/utils \
    $(MPI_FCFLAGS)

LDADD_COMMON = \
    $(top_builddir)/tests/common/libtest.a \
    $(top_builddir)/src/utils/libyac_utils.a \
    $(top_builddir)/src/core/libyac_core.a \
    $(LAPACK_CLIBS) \
    $(NETCDF_CLIBS)

LDADD = $(LDADD_COMMON) $(YAXT_CLIBS) $(MPI_CLIBS) $(RPATH_CFLAGS)
FCLDADD = $(LDADD_COMMON) $(OPENMP_FCLIBS)
if TEST_YAXT_FC
FCLDADD += $(YAXT_FCLIBS)
else
FCLDADD += $(YAXT_CLIBS)
endif
FCLDADD += $(MPI_FCLIBS) $(RPATH_FCFLAGS)

if ENABLE_FORTRAN_BINDINGS

if FCMODUC
utest_mod = UTEST.$(FCMODEXT)
else
utest_mod = utest.$(FCMODEXT)
endif

utest_FCDEPS = ../common/$(utest_mod) ../common/test_macros.inc

endif ENABLE_FORTRAN_BINDINGS

test_duplicate_stencils_parallel_x_LDADD = $(FCLDADD)
test_duplicate_stencils_parallel_x_SOURCES = test_duplicate_stencils_parallel.F90
test_duplicate_stencils_parallel.$(OBJEXT): $(utest_FCDEPS)

test_read_icon_x_LDADD = $(FCLDADD)
test_read_icon_x_SOURCES = test_read_icon.F90
test_read_icon.$(OBJEXT): $(utest_FCDEPS)

test_read_icon_parallel_c_x_SOURCES = test_read_icon_parallel_c.c test_read_icon_common.c test_read_icon_common.h

test_read_icon_parallel_x_LDADD = $(FCLDADD)
test_read_icon_parallel_x_SOURCES = test_read_icon_parallel.F90 test_read_icon_common.c
test_read_icon_parallel.$(OBJEXT): $(utest_FCDEPS)

test_read_scrip_x_LDADD = $(FCLDADD)
test_read_scrip_x_SOURCES = test_read_scrip.F90
test_read_scrip.$(OBJEXT): $(utest_FCDEPS)

test_read_scrip_parallel_x_LDADD = $(FCLDADD)
test_read_scrip_parallel_x_SOURCES = test_read_scrip_parallel.F90
test_read_scrip_parallel.$(OBJEXT): $(utest_FCDEPS)

if ENABLE_FORTRAN_BINDINGS
mostlyclean-local:
	rm -f *.$(FCMODEXT)
endif ENABLE_FORTRAN_BINDINGS

# Compile *.F90 without CPPFLAGS, which are normally meant for the C compiler
# and might not be compatible with the Fortran compiler:
FCDEFS = $(FCPPDEF)HAVE_CONFIG_H
DEFAULT_FCINCLUDES = $(FCPPINC).@yac_fcisrc@ $(FCPPINC)$(top_builddir)/src
PPFCCOMPILE = $(FC) $(FCDEFS) $(DEFAULT_FCINCLUDES) $(FCINCLUDES) \
	$(AM_FCFLAGS) $(FCFLAGS)
