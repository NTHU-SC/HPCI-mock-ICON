# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#

pkgpyexec_PYTHON = mtime/__init__.py
nodist_pkgpyexec_PYTHON = mtime/_mtime.py
pkgpyexec_LTLIBRARIES = mtime/__mtime.la

mtime___mtime_la_SOURCES =
mtime___mtime_la_LDFLAGS = -shared -avoid-version -module
mtime___mtime_la_LIBADD = $(top_builddir)/src/libmtime_c.la

all-local: mtime/.symlinkstamp

mtime/.symlinkstamp: mtime/__mtime.la
	$(AM_V_at)($(am__cd) mtime && rm -f __mtime$(acx_modext) && $(LN_S) $(acx_objdir)/__mtime$(acx_modext) __mtime$(acx_modext)) && touch $@

CLEANFILES = mtime/.symlinkstamp mtime/__mtime$(acx_modext)

# overrides to make regular compilations also use libtool
COMPILE = $(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=compile $(CC) -static $(DEFS) \
	$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	$(AM_CFLAGS) $(CFLAGS)
FCCOMPILE = $(LIBTOOL) $(AM_V_lt) --tag=FC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=compile $(FC) -static $(AM_FCFLAGS) $(FCFLAGS)
