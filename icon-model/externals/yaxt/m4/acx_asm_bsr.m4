dnl acx_asm_bsr.m4 --- check whether the compiler supports bsr
dnl                    instruction in inline assembly
dnl
dnl Copyright  (C)  2023  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords: inline assembly
dnl Author: Thomas Jahns <jahns@dkrz.de>
dnl Maintainer: Thomas Jahns <jahns@dkrz.de>
dnl URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
dnl
dnl Redistribution and use in source and binary forms, with or without
dnl modification, are  permitted provided that the following conditions are
dnl met:
dnl
dnl Redistributions of source code must retain the above copyright notice,
dnl this list of conditions and the following disclaimer.
dnl
dnl Redistributions in binary form must reproduce the above copyright
dnl notice, this list of conditions and the following disclaimer in the
dnl documentation and/or other materials provided with the distribution.
dnl
dnl Neither the name of the DKRZ GmbH nor the names of its contributors
dnl may be used to endorse or promote products derived from this software
dnl without specific prior written permission.
dnl
dnl THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
dnl IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
dnl TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
dnl PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
dnl OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
dnl EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
dnl PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
dnl PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
dnl LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
dnl NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
dnl SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
dnl
dnl
dnl ACX_CHECK_ASM_BSR([ACTION-IF-AVAILABLE],[ACTION-IF-NOT-FOUND])
dnl ---------------------------------------------------------------
AC_DEFUN([ACX_CHECK_ASM_BSR],
  [AC_CHECK_SIZEOF([long])
   AC_CACHE_CHECK([whether inline assembly can be used to access bsr instruction],
     [acx_cv_have_asm_bsr],
     [AC_LINK_IFELSE([AC_LANG_PROGRAM(,
       [  unsigned long v = 56, ms1bpos;
@%:@if SIZEOF_LONG == 8
  __asm__ ("bsrq %1, %0" : "=r" (ms1bpos) : "r" (v));
@%:@elif SIZEOF_LONG == 4
  __asm__ ("bsrl %1, %0" : "=r" (ms1bpos) : "r" (v));
@%:@else
@%:@error "unexpected size of size_t!"
@%:@endif]
  return ms1bpos == 6;)],
     [AS_VAR_SET([acx_cv_have_asm_bsr],[yes])],
     [AS_VAR_SET([acx_cv_have_asm_bsr],[no])])])
   AS_IF([test "$acx_cv_have_asm_bsr" = yes],
     [m4_default([$2],[AC_DEFINE([HAVE_ASM_BSR],[1],
       [defined to 1 if inline assembly can emit bsrl or bsrq instruction])])],
     [$3])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/"
dnl license-default: "bsd"
dnl End:
dnl
