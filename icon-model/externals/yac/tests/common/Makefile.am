# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

check_LIBRARIES = libtest.a

libtest_a_SOURCES = \
    dist_grid_utils.c \
    dist_grid_utils.h \
    test_common.c \
    test_common.h \
    test_macros.inc \
    tests.c \
    tests.h

if ENABLE_FORTRAN_BINDINGS
libtest_a_SOURCES += \
    utest.F90
endif ENABLE_FORTRAN_BINDINGS

if HAVE_NETCDF
libtest_a_SOURCES += \
    grid_file_common.c \
    grid_file_common.h \
    weight_file_common.c \
    weight_file_common.h
endif

AM_CPPFLAGS = \
    -I$(top_srcdir)/src/core

AM_CFLAGS = \
    $(MPI_CFLAGS) \
    $(YAXT_CFLAGS) \
    $(NETCDF_CFLAGS)

EXTRA_DIST = utest.F90 test_macros.inc

if ENABLE_FORTRAN_BINDINGS

if FCMODUC
utest_mod = UTEST.$(FCMODEXT)
else
utest_mod = utest.$(FCMODEXT)
endif

$(utest_mod): utest.$(OBJEXT)
	@if test ! -f $@; then rm -f $<; $(MAKE) $(AM_MAKEFLAGS) $<; fi

mostlyclean-local:
	rm -f *.$(FCMODEXT)

endif ENABLE_FORTRAN_BINDINGS

# Compile *.F90 without CPPFLAGS, which are normally meant for the C compiler
# and might not be compatible with the Fortran compiler:
FCDEFS = $(FCPPDEF)HAVE_CONFIG_H
DEFAULT_FCINCLUDES = $(FCPPINC).@yac_fcisrc@ $(FCPPINC)$(top_builddir)/src
PPFCCOMPILE = $(FC) $(FCDEFS) $(DEFAULT_FCINCLUDES) $(FCINCLUDES) \
	$(AM_FCFLAGS) $(FCFLAGS)
