/**
 * @file xt_ddt.c
 *
 * @copyright Copyright  (C)  2022 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdbool.h>
#include <string.h>
#include <mpi.h>

#ifdef _OPENACC
#define STR(s) #s
#define xt_Pragma(args) _Pragma(args)
#define XtPragmaACC(args) xt_Pragma(STR(acc args))
#else
#define XtPragmaACC(args)
#endif

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_ddt.h"
#include "xt_ddt_internal.h"

//static const char filename[] = "xt_ddt.c";


static void xt_ddt_pack_8(
  size_t count, ssize_t *restrict displs, const uint8_t *restrict src,
  uint8_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_pack_16(
  size_t count, ssize_t *restrict displs, const uint16_t *restrict src,
  uint16_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_pack_32(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_pack_32_2(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_pack_96(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_pack_64(
  size_t count, ssize_t *restrict displs, const uint64_t *restrict src,
  uint64_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_pack_128(
  size_t count, ssize_t *restrict displs, const uint64_t *restrict src,
  uint64_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_pack_160(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_pack_256(
  size_t count, ssize_t *restrict displs, const uint64_t *restrict src,
  uint64_t (*restrict dst)[4], enum xt_memtype memtype);

static void xt_ddt_unpack_8(
  size_t count, ssize_t *restrict displs, const uint8_t *restrict src,
  uint8_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_unpack_16(
  size_t count, ssize_t *restrict displs, const uint16_t *restrict src,
  uint16_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_unpack_32(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_unpack_32_2(
  size_t count, ssize_t *restrict displs, const uint32_t (*restrict src)[2],
  uint32_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_unpack_96(
  size_t count, ssize_t *restrict displs, const uint32_t (*restrict src)[3],
  uint32_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_unpack_64(
  size_t count, ssize_t *restrict displs, const uint64_t *restrict src,
  uint64_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_unpack_128(
  size_t count, ssize_t *restrict displs, const uint64_t (*restrict src)[2],
  uint64_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_unpack_160(
  size_t count, ssize_t *restrict displs, const uint32_t (*restrict src)[5],
  uint32_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_unpack_256(
  size_t count, ssize_t *restrict displs, const uint64_t (*restrict src)[4],
  uint64_t *restrict dst, enum xt_memtype memtype);

struct xt_ddt_kernels xt_ddt_valid_kernels[] = {
  {.base_pack_size = 1,
   .element_size = 1,
   .pack = (xt_ddt_kernel_func)xt_ddt_pack_8,
   .unpack = (xt_ddt_kernel_func)xt_ddt_unpack_8},
  {.base_pack_size = 2,
   .element_size = 2,
   .pack = (xt_ddt_kernel_func)xt_ddt_pack_16,
   .unpack = (xt_ddt_kernel_func)xt_ddt_unpack_16},
  {.base_pack_size = 4,
   .element_size = 4,
   .pack = (xt_ddt_kernel_func)xt_ddt_pack_32,
   .unpack = (xt_ddt_kernel_func)xt_ddt_unpack_32},
  {.base_pack_size = 4,
   .element_size = 8,
   .pack = (xt_ddt_kernel_func)xt_ddt_pack_32_2,
   .unpack = (xt_ddt_kernel_func)xt_ddt_unpack_32_2},
  {.base_pack_size = 8,
   .element_size = 8,
   .pack = (xt_ddt_kernel_func)xt_ddt_pack_64,
   .unpack = (xt_ddt_kernel_func)xt_ddt_unpack_64},
  {.base_pack_size = 4,
   .element_size = 12,
   .pack = (xt_ddt_kernel_func)xt_ddt_pack_96,
   .unpack = (xt_ddt_kernel_func)xt_ddt_unpack_96},
  {.base_pack_size = 8,
   .element_size = 16,
   .pack = (xt_ddt_kernel_func)xt_ddt_pack_128,
   .unpack = (xt_ddt_kernel_func)xt_ddt_unpack_128},
  {.base_pack_size = 4,
   .element_size = 20,
   .pack = (xt_ddt_kernel_func)xt_ddt_pack_160,
   .unpack = (xt_ddt_kernel_func)xt_ddt_unpack_160},
  {.base_pack_size = 8,
   .element_size = 32,
   .pack = (xt_ddt_kernel_func)xt_ddt_pack_256,
   .unpack = (xt_ddt_kernel_func)xt_ddt_unpack_256},
};

size_t xt_ddt_get_pack_size_internal(Xt_ddt ddt) {

  return (ddt == NULL)?0:(ddt->pack_size);
}

size_t xt_ddt_get_pack_size(MPI_Datatype mpi_ddt) {

  return xt_ddt_get_pack_size_internal(xt_ddt_from_mpi_ddt(mpi_ddt));
}

static void xt_ddt_copy_displs(Xt_ddt ddt, enum xt_memtype memtype) {

  // count total number of displacements
  size_t total_displs_size = 0, count = ddt->count;
  for (size_t i = 0; i < count; ++i)
    total_displs_size += ddt->data[i].displ_count;

  // allocate displacements in specified memory type
  ssize_t *displs;
  size_t buffer_size = total_displs_size * sizeof(*displs);
  displs = xt_gpu_malloc(buffer_size, memtype);

  // copy displacements from host to specified memory type
  xt_gpu_memcpy(
    displs, ddt->data[0].displs[XT_MEMTYPE_HOST],
    buffer_size, memtype, XT_MEMTYPE_HOST);

  // set displacements for all data entries
  for (size_t i = 0, offset = 0; i < count; ++i) {
    ddt->data[i].displs[memtype] = displs + offset;
    offset += ddt->data[i].displ_count;
  }

  ddt->displs_available[memtype] = 1;
}

#define add_rhs_byte_displ(rtype,ptr,disp) \
  ((const rtype *)(const void *)((const unsigned char *)(ptr) + (disp)))

static void xt_ddt_pack_8(
  size_t count, ssize_t *restrict displs, const uint8_t *restrict src,
  uint8_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i)
    dst[i] = *add_rhs_byte_displ(uint8_t, src, displs[i]);
}

static void xt_ddt_pack_16(
  size_t count, ssize_t *restrict displs, const uint16_t *restrict src,
  uint16_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i)
    dst[i] = *add_rhs_byte_displ(uint16_t, src, + displs[i]);
}

static void xt_ddt_pack_32(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i)
    dst[i] = *add_rhs_byte_displ(uint32_t, src, + displs[i]);
}

static void xt_ddt_pack_32_2(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype) {
  uint32_t (*restrict dst_)[2] = (uint32_t(*)[2])dst;
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst_, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    const uint32_t *src_32 = add_rhs_byte_displ(uint32_t, src, displs[i]);
XtPragmaACC(loop independent)
    for (int j = 0; j < 2; ++j) dst_[i][j] = src_32[j];
  }
}

static void xt_ddt_pack_96(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype) {
  uint32_t (*restrict dst_)[3] = (uint32_t(*)[3])dst;
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst_, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    const uint32_t *src_32 = add_rhs_byte_displ(uint32_t, src, displs[i]);
XtPragmaACC(loop independent)
    for (int j = 0; j < 3; ++j) dst_[i][j] = src_32[j];
  }
}

static void xt_ddt_pack_64(
  size_t count, ssize_t *restrict displs, const uint64_t *restrict src,
  uint64_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i)
    dst[i] = *add_rhs_byte_displ(uint64_t, src, displs[i]);
}

static void xt_ddt_pack_128(
  size_t count, ssize_t *restrict displs, const uint64_t *restrict src,
  uint64_t *restrict dst, enum xt_memtype memtype) {
  uint64_t (*restrict dst_)[2] = (uint64_t(*)[2])dst;
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst_, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    const uint64_t *src_64 = add_rhs_byte_displ(uint64_t, src, displs[i]);
XtPragmaACC(loop independent)
    for (int j = 0; j < 2; ++j) dst_[i][j] = src_64[j];
  }
}

static void xt_ddt_pack_160(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype) {
  uint32_t (*restrict dst_)[5] = (uint32_t(*)[5])dst;
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst_, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    const uint32_t *src_32 = add_rhs_byte_displ(uint32_t, src, displs[i]);
XtPragmaACC(loop independent)
    for (int j = 0; j < 5; ++j) dst_[i][j] = src_32[j];
  }
}

static void xt_ddt_pack_256(
  size_t count, ssize_t *restrict displs, const uint64_t *restrict src,
  uint64_t (*restrict dst)[4], enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    const uint64_t *src_64 = add_rhs_byte_displ(uint64_t, src, displs[i]);
XtPragmaACC(loop independent)
    for (int j = 0; j < 4; ++j) dst[i][j] = src_64[j];
  }
}

void xt_ddt_pack_internal(
  Xt_ddt ddt, const void *src, void *dst, enum xt_memtype memtype) {

  XT_GPU_INSTR_PUSH(xt_ddt_pack_internal);

  size_t dst_offset = 0;

  // if the displacements are not avaible in the required memory type
  if (!ddt->displs_available[memtype]) xt_ddt_copy_displs(ddt, memtype);

  size_t count = ddt->count;

  // for all sections with the same elemental datatype extent
  for (size_t i = 0; i < count; ++i) {

    struct xt_ddt_kernels * kernel =
      &xt_ddt_valid_kernels[ddt->data[i].kernel_idx];
    size_t displ_count = ddt->data[i].displ_count;
    kernel->pack(
      displ_count, ddt->data[i].displs[memtype], src,
      (unsigned char *)dst + dst_offset, memtype);

    dst_offset += displ_count * kernel->element_size;
  }
  XT_GPU_INSTR_POP;
}

void xt_ddt_pack(MPI_Datatype mpi_ddt, const void *src, void *dst) {

  XT_GPU_INSTR_PUSH(xt_ddt_pack);
  XT_GPU_INSTR_PUSH(xt_ddt_pack:initialise);

  enum xt_memtype src_memtype = xt_gpu_get_memtype(src);
  enum xt_memtype dst_memtype = xt_gpu_get_memtype(dst);

  size_t pack_size;
  void *orig_dst;
  /* pacify buggy -Wmaybe-uninitialized */
#if defined __GNUC__ && __GNUC__ <= 11
  pack_size = 0;
  orig_dst = NULL;
#endif

  // if the source and destination are in different memory types
  if (src_memtype != dst_memtype) {
    pack_size = xt_ddt_get_pack_size(mpi_ddt);
    orig_dst = dst;
    dst = xt_gpu_malloc(pack_size, src_memtype);
  }

  XT_GPU_INSTR_POP; //xt_ddt_pack:initialise

  xt_ddt_pack_internal(
    xt_ddt_from_mpi_ddt(mpi_ddt), src, dst, src_memtype);

  XT_GPU_INSTR_PUSH(xt_ddt_pack:finalise);

  // if the source and destination are in different memory types
  if (src_memtype != dst_memtype) {
    xt_gpu_memcpy(orig_dst, dst, pack_size, dst_memtype, src_memtype);
    xt_gpu_free(dst, src_memtype);
  }

  XT_GPU_INSTR_POP; // xt_ddt_pack:finalise
  XT_GPU_INSTR_POP; // xt_ddt_pack
}

static void xt_ddt_unpack_8(
  size_t count, ssize_t *restrict displs, const uint8_t *restrict src,
  uint8_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i)
    dst[displs[i]] = src[i];
}


static void xt_ddt_unpack_16(
  size_t count, ssize_t *restrict displs, const uint16_t *restrict src,
  uint16_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    uint16_t *dst_ = (void *)((unsigned char *)dst + displs[i]);
    dst_[0] = src[i];
  }
}

static void xt_ddt_unpack_32(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    uint32_t *dst_ = (void *)((unsigned char *)dst + displs[i]);
    dst_[0] = src[i];
  }
}

static void xt_ddt_unpack_32_2(
  size_t count, ssize_t *restrict displs, const uint32_t (*restrict src)[2],
  uint32_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    uint32_t *dst_32 = (void *)((unsigned char *)dst + displs[i]);
    dst_32[0] = src[i][0];
    dst_32[1] = src[i][1];
  }
}

static void xt_ddt_unpack_96(
  size_t count, ssize_t *restrict displs, const uint32_t (*restrict src)[3],
  uint32_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    uint32_t *dst_32 = (void *)((unsigned char *)dst + displs[i]);
    dst_32[0] = src[i][0];
    dst_32[1] = src[i][1];
    dst_32[2] = src[i][2];
  }
}

static void xt_ddt_unpack_64(
  size_t count, ssize_t *restrict displs, const uint64_t *restrict src,
  uint64_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    uint64_t *dst_ = (void *)((unsigned char *)dst + displs[i]);
    dst_[0] = src[i];
  }
}

static void xt_ddt_unpack_128(
  size_t count, ssize_t *restrict displs, const uint64_t (*restrict src)[2],
  uint64_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    uint64_t *dst_64 = (void *)((unsigned char *)dst + displs[i]);
    dst_64[0] = src[i][0];
    dst_64[1] = src[i][1];
  }
}

static void xt_ddt_unpack_160(
  size_t count, ssize_t *restrict displs, const uint32_t (*restrict src)[5],
  uint32_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    uint32_t *dst_32 = (void *)((unsigned char *)dst + displs[i]);
    dst_32[0] = src[i][0];
    dst_32[1] = src[i][1];
    dst_32[2] = src[i][2];
    dst_32[3] = src[i][3];
    dst_32[4] = src[i][4];
  }
}

static void xt_ddt_unpack_256(
  size_t count, ssize_t *restrict displs, const uint64_t (*restrict src)[4],
  uint64_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    uint64_t *dst_64 = (void *)((unsigned char *)dst + displs[i]);
    dst_64[0] = src[i][0];
    dst_64[1] = src[i][1];
    dst_64[2] = src[i][2];
    dst_64[3] = src[i][3];
  }
}

void xt_ddt_unpack_internal(
  Xt_ddt ddt, const void *src, void *dst, enum xt_memtype memtype) {

  XT_GPU_INSTR_PUSH(xt_ddt_unpack_internal);

  size_t src_offset = 0;

  // if the displacements are not avaible in the required memory type
  if (!ddt->displs_available[memtype]) xt_ddt_copy_displs(ddt, memtype);

  size_t count = ddt->count;

  // for all sections with the same elemental datatype extent
  for (size_t i = 0; i < count; ++i) {

    struct xt_ddt_kernels * kernel =
      &xt_ddt_valid_kernels[ddt->data[i].kernel_idx];
    size_t displ_count = ddt->data[i].displ_count;
    kernel->unpack(
      displ_count, ddt->data[i].displs[memtype],
      add_rhs_byte_displ(void, src, src_offset), dst, memtype);

    src_offset += displ_count * kernel->element_size;
  }
  XT_GPU_INSTR_POP;
}

void xt_ddt_unpack(MPI_Datatype mpi_ddt, const void *src, void *dst) {

  XT_GPU_INSTR_PUSH(xt_ddt_unpack);
  XT_GPU_INSTR_PUSH(xt_ddt_unpack:initialise);

  enum xt_memtype src_memtype = xt_gpu_get_memtype(src);
  enum xt_memtype dst_memtype = xt_gpu_get_memtype(dst);

  void *src__ = NULL;
  const void *src_;
  // if the source and destination are in different memory types
  if (src_memtype != dst_memtype) {
    size_t pack_size = xt_ddt_get_pack_size(mpi_ddt);
    src__ = xt_gpu_malloc(pack_size, dst_memtype);
    xt_gpu_memcpy(src__, src, pack_size, dst_memtype, src_memtype);
    src_ = src__;
  } else
    src_ = src;

  XT_GPU_INSTR_POP; // xt_ddt_unpack:initialise

  xt_ddt_unpack_internal(
    xt_ddt_from_mpi_ddt(mpi_ddt), src_, dst, dst_memtype);

  XT_GPU_INSTR_PUSH(xt_ddt_unpack:finalise);

  // if the source and destination are in different memory types
  if (src_memtype != dst_memtype) xt_gpu_free(src__, dst_memtype);

  XT_GPU_INSTR_POP; // xt_ddt_unpack:finalise
  XT_GPU_INSTR_POP; // xt_ddt_unpack
}

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
