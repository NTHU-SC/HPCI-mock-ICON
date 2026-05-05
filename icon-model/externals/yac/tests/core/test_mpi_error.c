// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <mpi.h>

#include "yac_mpi_common.h"

/** \file test_mpi_error.c
 *  \test
 * This example tests the yac_mpi_error interface.
 */

int main (void) {

  MPI_Init(NULL, NULL);

  yac_mpi_error(MPI_ERR_COMM, MPI_COMM_WORLD);

  MPI_Finalize();

  return EXIT_SUCCESS;
}

