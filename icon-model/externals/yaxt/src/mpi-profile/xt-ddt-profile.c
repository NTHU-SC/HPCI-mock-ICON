/**
 * @file xt-ddt-profile.c
 * @brief track creation and destruction of MPI derived datatypes
 *
 * @copyright Copyright  (C)  2023 Jörg Behrens <behrens@dkrz.de>
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
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

static FILE *dt_out = NULL;

static void init_instr(void) __attribute__((constructor));

static void init_instr(void)
{
  static const char rank_vars[][24]
    = { { "MPI_LOCALRANKID" } , { "PMI_RANK" }, { 0 } };
  enum {
    nrank_vars = sizeof (rank_vars) / sizeof (rank_vars[0]) - 1,
  };
  const char *rank_str;
  for (size_t i = 0; i < nrank_vars; ++i) {
    if ((rank_str = getenv(rank_vars[i])))
      goto found_rank;
  }
  rank_str = "unknown";
found_rank:
  ;
  static const char dt_pfx[] = "dt_log.",
    dt_sfx[] = ".txt";
  size_t pfx_size = sizeof (dt_pfx), sfx_size = sizeof (dt_sfx),
    rank_str_len = strlen(rank_str);
  char dt_log_fn[pfx_size + sfx_size - 1 + rank_str_len];
  memcpy(dt_log_fn, dt_pfx, pfx_size - 1);
  memcpy(dt_log_fn+pfx_size-1, rank_str, rank_str_len);
  memcpy(dt_log_fn+pfx_size-1+rank_str_len, dt_sfx, sfx_size);
  dt_out = fopen(dt_log_fn, "w");
  setlinebuf(dt_out);
}


int
MPI_Type_dup(MPI_Datatype oldtype, MPI_Datatype *newtype)
{
  fprintf(dt_out, "MPI_Type_dup(%jd);\n", (intmax_t)oldtype);
  int rc = PMPI_Type_dup(oldtype, newtype);
  fprintf(dt_out, "%jd\n", (intmax_t)(*newtype));
  return rc;
}

int
MPI_Type_contiguous(int count, MPI_Datatype oldtype,
                    MPI_Datatype *newtype)
{
  fprintf(dt_out, "MPI_Type_contiguous(%d, %jd);\n", count, (intmax_t)oldtype);
  int rc = PMPI_Type_contiguous(count, oldtype, newtype);
  fprintf(dt_out, "%jd\n", (intmax_t)(*newtype));
  return rc;
}

int
MPI_Type_vector(int count, int blocklength, int stride,
                MPI_Datatype oldtype, MPI_Datatype *newtype)
{
  fprintf(dt_out, "MPI_Type_vector(%d, %d, %d, %jd);\n",
          count, blocklength, stride, (intmax_t)oldtype);
  int rc = PMPI_Type_vector(count, blocklength, stride, oldtype, newtype);
  fprintf(dt_out, "%jd\n", (intmax_t)(*newtype));
  return rc;
}

int
MPI_Type_hvector(int count, int blocklength, MPI_Aint stride,
                 MPI_Datatype oldtype, MPI_Datatype *newtype)
{
  fprintf(dt_out, "MPI_Type_hvector(%d, %d, %jd, %jd);\n",
          count, blocklength, (intmax_t)stride, (intmax_t)oldtype);
  int rc = PMPI_Type_hvector(count, blocklength, stride, oldtype, newtype);
  fprintf(dt_out, "%jd\n", (intmax_t)(*newtype));
  return rc;
}

int
MPI_Type_create_hvector(int count, int blocklength,
                        MPI_Aint stride, MPI_Datatype oldtype,
                        MPI_Datatype *newtype)
{
  fprintf(dt_out, "MPI_Type_create_hvector(%d, %d, %jd, %jd);\n",
          count, blocklength, (intmax_t)stride, (intmax_t)oldtype);
  int rc = PMPI_Type_create_hvector(count, blocklength, stride, oldtype, newtype);
  fprintf(dt_out, "%jd\n", (intmax_t)(*newtype));
  return rc;
}


int
MPI_Type_indexed(int count,
                 XT_MPI2_CONST int array_of_blocklengths[],
                 XT_MPI2_CONST int array_of_displacements[],
                 MPI_Datatype oldtype,
                 MPI_Datatype *newtype)
{
  fprintf(dt_out, "MPI_Type_indexed(%d, (int){\n", count);
  for (int i = 0; i < count; ++i)
    fprintf(dt_out, "  %d,\n", array_of_blocklengths[i]);
  fputs("}, (int[]){\n", dt_out);
  for (int i = 0; i < count; ++i)
    fprintf(dt_out, "  %d,\n", array_of_displacements[i]);
  fprintf(dt_out, "}, %jd);\n", (intmax_t)oldtype);
  int rc = PMPI_Type_indexed(count, array_of_blocklengths,
                             array_of_displacements, oldtype, newtype);
  fprintf(dt_out, "%jd\n", (intmax_t)(*newtype));
  return rc;
}

int
MPI_Type_hindexed(int count,
                  int array_of_blocklengths[],
                  MPI_Aint array_of_displacements[],
                  MPI_Datatype oldtype,
                  MPI_Datatype *newtype)
{
  fprintf(dt_out, "MPI_Type_hindexed(%d, (int[]){\n", count);
  for (int i = 0; i < count; ++i)
    fprintf(dt_out, "  %d,\n", array_of_blocklengths[i]);
  fputs("}, (MPI_Aint[]){\n", dt_out);
  for (int i = 0; i < count; ++i)
    fprintf(dt_out, "  %jd,\n", (intmax_t)array_of_displacements[i]);
  fprintf(dt_out, "}, %jd);\n", (intmax_t)oldtype);
  int rc = PMPI_Type_hindexed(count, array_of_blocklengths,
                              array_of_displacements, oldtype, newtype);
  fprintf(dt_out, "%jd\n", (intmax_t)(*newtype));
  return rc;
}


int
MPI_Type_create_hindexed(int count,
                         XT_MPI2_CONST int array_of_blocklengths[],
                         XT_MPI2_CONST MPI_Aint array_of_displacements[],
                         MPI_Datatype oldtype, MPI_Datatype *newtype)
{
  fprintf(dt_out, "MPI_Type_create_hindexed(%d, (int){\n", count);
  for (int i = 0; i < count; ++i)
    fprintf(dt_out, "  %d,\n", array_of_blocklengths[i]);
  fputs("}, (MPI_Aint[]){\n", dt_out);
  for (int i = 0; i < count; ++i)
    fprintf(dt_out, "  %jd,\n", (intmax_t)array_of_displacements[i]);
  fprintf(dt_out, "}, %jd);\n", (intmax_t)oldtype);
  int rc = PMPI_Type_create_hindexed(count, array_of_blocklengths,
                                     array_of_displacements, oldtype,
                                     newtype);
  fprintf(dt_out, "%jd\n", (intmax_t)(*newtype));
  return rc;
}

int
MPI_Type_create_indexed_block(int count, int blocklength,
                              XT_MPI2_CONST int array_of_displacements[],
                              MPI_Datatype oldtype, MPI_Datatype *newtype)
{
  fprintf(dt_out, "MPI_Type_create_indexed_block(%d, %d, (int){\n", count, blocklength);
  for (int i = 0; i < count; ++i)
    fprintf(dt_out, "  %d,\n", array_of_displacements[i]);
  fprintf(dt_out, "}, %jd);\n", (intmax_t)oldtype);
  int rc = PMPI_Type_create_indexed_block(count, blocklength,
                                          array_of_displacements,
                                          oldtype, newtype);
  fprintf(dt_out, "%jd\n", (intmax_t)(*newtype));
  return rc;
}

#if MPI_VERSION >= 3
int
MPI_Type_create_hindexed_block(
  int count, int blocklength, XT_MPI2_CONST MPI_Aint array_of_displacements[],
  MPI_Datatype oldtype, MPI_Datatype *newtype)
{
  fprintf(dt_out, "MPI_Type_create_indexed_block(%d, %d, (MPI_Aint){\n", count, blocklength);
  for (int i = 0; i < count; ++i)
    fprintf(dt_out, "  %jd,\n", (intmax_t)array_of_displacements[i]);
  fprintf(dt_out, "}, %jd);\n", (intmax_t)oldtype);
  int rc = PMPI_Type_create_hindexed_block(count, blocklength,
                                           array_of_displacements,
                                           oldtype, newtype);
  fprintf(dt_out, "%jd\n", (intmax_t)(*newtype));
  return rc;
}
#endif

int
MPI_Type_create_resized(MPI_Datatype oldtype, MPI_Aint lb,
                        MPI_Aint extent, MPI_Datatype *newtype)
{
  fprintf(dt_out, "MPI_Type_create_resized(%jd, %jd, %jd);\n",
          (intmax_t)oldtype, (intmax_t)lb, (intmax_t)extent);
  int rc = PMPI_Type_create_resized(oldtype, lb, extent, newtype);
  fprintf(dt_out, "%jd\n", (intmax_t)(*newtype));
  return rc;
}


int
MPI_Type_struct(int count, int *array_of_blocklengths,
                MPI_Aint *array_of_displacements, MPI_Datatype *array_of_types,
                MPI_Datatype *newtype)
{
  fprintf(dt_out, "MPI_Type_struct(%d, (int){\n", count);
  for (int i = 0; i < count; ++i)
    fprintf(dt_out, "  %d,\n", array_of_blocklengths[i]);
  fputs("}, (MPI_Aint[]){\n", dt_out);
  for (int i = 0; i < count; ++i)
    fprintf(dt_out, "  %jd,\n", (intmax_t)array_of_displacements[i]);
  fputs("}, (MPI_Datatype[]){\n", dt_out);
  for (int i = 0; i < count; ++i)
    fprintf(dt_out, "  %jd,\n", (intmax_t)array_of_types[i]);
  fputs("});\n", dt_out);
  int rc = PMPI_Type_struct(
    count, array_of_blocklengths, array_of_displacements,
    array_of_types, newtype);
  fprintf(dt_out, "%jd\n", (intmax_t)(*newtype));
  return rc;
}

int
MPI_Type_create_struct(int count, XT_MPI2_CONST int array_of_block_lengths[],
                       XT_MPI2_CONST MPI_Aint array_of_displacements[],
                       XT_MPI2_CONST MPI_Datatype array_of_types[],
                       MPI_Datatype *newtype)
{
  fprintf(dt_out, "MPI_Type_create_struct(%d, (int){\n", count);
  for (int i = 0; i < count; ++i)
    fprintf(dt_out, "  %d,\n", array_of_block_lengths[i]);
  fputs("}, (MPI_Aint[]){\n", dt_out);
  for (int i = 0; i < count; ++i)
    fprintf(dt_out, "  %jd,\n", (intmax_t)array_of_displacements[i]);
  fputs("}, (MPI_Datatype[]){\n", dt_out);
  for (int i = 0; i < count; ++i)
    fprintf(dt_out, "  %jd,\n", (intmax_t)array_of_types[i]);
  fputs("});\n", dt_out);
  int rc = PMPI_Type_create_struct(
    count, array_of_block_lengths, array_of_displacements,
    array_of_types, newtype);
  fprintf(dt_out, "%jd\n", (intmax_t)(*newtype));
  return rc;
}

int
MPI_Type_create_darray(int size, int rank, int ndims,
                       XT_MPI2_CONST int gsize_array[], XT_MPI2_CONST int distrib_array[],
                       XT_MPI2_CONST int darg_array[], XT_MPI2_CONST int psize_array[],
                       int order, MPI_Datatype oldtype,
                       MPI_Datatype *newtype)
{
  int rc = PMPI_Type_create_darray(
    size, rank, ndims, gsize_array, distrib_array, darg_array, psize_array,
    order, oldtype, newtype);
  fprintf(dt_out, "%jd\n", (intmax_t)(*newtype));
  return rc;
}

int
MPI_Type_create_subarray(int ndims,
                         XT_MPI2_CONST int size_array[],
                         XT_MPI2_CONST int subsize_array[],
                         XT_MPI2_CONST int start_array[],
                         int order,
                         MPI_Datatype oldtype,
                         MPI_Datatype *newtype)
{
  int rc = PMPI_Type_create_subarray(ndims, size_array, subsize_array,
                                     start_array, order, oldtype, newtype);
  fprintf(dt_out, "%jd\n", (intmax_t)(*newtype));
  return rc;
}


int
MPI_Type_free(MPI_Datatype *datatype)
{
  fprintf(dt_out, "MPI_Type_free(&%jd);\n", (intmax_t)(*datatype));
  int rc = PMPI_Type_free(datatype);
  return rc;
}

int
MPI_Type_commit(MPI_Datatype *datatype)
{
  fprintf(dt_out, "MPI_Type_commit(&%jd);\n", (intmax_t)(*datatype));
  int rc = PMPI_Type_commit(datatype);
  return rc;
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
