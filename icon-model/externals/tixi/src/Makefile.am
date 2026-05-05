# SPDX-FileCopyrightText: 2025 DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# SPDX-License-Identifier: BSD-3-Clause

lib_LIBRARIES = libicon_tixi.a

if FCMODUC
tixi_mod = TIXI.$(FCMODEXT)
else !FCMODUC
tixi_mod = tixi.$(FCMODEXT)
endif !FCMODUC

nodist_include_HEADERS = $(tixi_mod)

libicon_tixi_a_SOURCES = \
  tixi.f90 \
  tixi.h \
  tixiData.h \
  tixiImpl.c \
  tixiInternal.c \
  tixiInternal.h \
  tixiUtils.c \
  tixiUtils.h \
  uidHelper.c \
  uidHelper.h \
  xpathFunctions.c \
  xpathFunctions.h

nodist_libicon_tixi_a_SOURCES = \
  tixi_version.h

AM_CFLAGS = $(XML2_CFLAGS)

$(tixi_mod): tixi.$(OBJEXT)
	@if test ! -f $@; then rm -f $<; $(MAKE) $<; fi

MOSTLYCLEANFILES = $(tixi_mod)
