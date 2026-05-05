# SPDX-FileCopyrightText: 2025 DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# SPDX-License-Identifier: BSD-3-Clause

m4_define([tixi_version], [2, 2, 2])

AC_PREREQ([2.69])
AC_INIT(
  [icon-tixi],
  [m4_argn(1, tixi_version).m4_argn(2, tixi_version).m4_argn(3, tixi_version)])

AC_SUBST([TIXI_VERSION_MAJOR], ['m4_argn([1], tixi_version)'])
AC_SUBST([TIXI_VERSION_MINOR], ['m4_argn([2], tixi_version)'])
AC_SUBST([TIXI_VERSION_PATCH], ['m4_argn([3], tixi_version)'])
AC_SUBST([TIXI_VERSION], [$PACKAGE_VERSION])

AC_CONFIG_MACRO_DIR([m4])
dnl Make sure all required m4 files are in place:
m4_pattern_forbid([^(A[SC]X)_])dnl
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([src/tixi.h])

AM_INIT_AUTOMAKE([1.16.1 foreign])
dnl Enable silent rules by default with a workaround for the case when the
dnl version constraint above is relaxed:
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
AM_MAINTAINER_MODE

dnl Tests for the C compiler:
dnl Fix a bug in Autoconf 2.70 (see https://savannah.gnu.org/support/?110396):
m4_if(m4_version_compare(m4_defn([AC_AUTOCONF_VERSION]), [2.70]), [0],
  [m4_pushdef([_AC_PROG_CC_C99],
        m4_bpatsubst(m4_dquote(m4_defn([_AC_PROG_CC_C99])),
          [\[ac_c_conftest_c89_program\]], [[ac_c_conftest_c99_program]]))dnl
   AC_PROG_CC
   m4_popdef([_AC_PROG_CC_C99])],
  [AC_PROG_CC])
dnl Expand the legacy macro for older version of Autoconf (2.69 and older):
m4_if(m4_version_compare(m4_defn([AC_AUTOCONF_VERSION]), [2.70]), [-1],
  [AC_PROG_CC_C99])
AS_VAR_IF([ac_cv_prog_cc_c99], [no],
  [AC_MSG_FAILURE([C compiler does not support ISO C99])])

dnl Tests for the Fortran compiler:
AC_PROG_FC
AC_LANG_PUSH([Fortran])
AC_FC_SRCEXT([f90])
ACX_FC_MODULE_NAMING([AC_SUBST([FCMODEXT], ["$acx_cv_fc_module_naming_ext"])])
AM_CONDITIONAL([FCMODUC], [test x"$acx_cv_fc_module_naming_upper" = xyes])
ACX_FC_MODULE_CHECK([ISO_C_BINDING])
AC_LANG_POP([Fortran])

dnl Tests for static library generation tools:
AC_ARG_VAR([AR], [archiver command])
AC_CHECK_TOOL([AR], [ar], [:])
AS_VAR_IF([AR], [:], [AC_MSG_ERROR([the archiver command not found])])
AC_ARG_VAR([ARFLAGS], [archiver flags])
ARFLAGS=${ARFLAGS-cr}
AC_ARG_VAR([RANLIB], [archive indexer command])
AC_PROG_RANLIB

dnl Checks for XML library:
ACX_LANG_PACKAGE_INIT([xml2], [-I], [-L], [/include/libxml2])

tixi_save_CFLAGS=$CFLAGS
XML2_CFLAGS=${XML2_CFLAGS-$acx_c_xml2_inc_search_args}
CFLAGS="$XML2_CFLAGS $CFLAGS"
AC_CHECK_HEADER([libxml/xmlversion.h],
  [ACX_LANG_LIB_SEARCH([XML2_CLIBS], [xmlReadFile],
     [ASX_PREPEND_LDFLAGS([$acx_c_xml2_lib_search_args], [-lxml2])],
     [XML2_CLIBS=$acx_cv_c_lib_func_xmlReadFile],
     [AC_MSG_FAILURE([cannot link to libxml2 library])])],
  [AC_MSG_FAILURE([libxml/xmlversion.h not found])], [//])
CFLAGS=$tixi_save_CFLAGS

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  src/tixi_version.h
])

AC_OUTPUT
