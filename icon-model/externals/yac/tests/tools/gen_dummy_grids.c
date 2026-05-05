// Copyright (c) 2025 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grid_file_common.h"

int main (int argc, char* argv[]) {

  if (argc != 5) exit(EXIT_FAILURE);

  {
    size_t num_lon = 271, num_lat = 181;
    double lon_range[] = {-135.0, 135.0}, lat_range[] = {-90.0, 90.0};
    write_dummy_exodus_grid_file(
      argv[1], num_lon, num_lat, lon_range, lat_range);
  }

  {
    write_dummy_scrip_grid_file(
      argv[2], argv[3], argv[4], 1,
      360, 10, (double[]){0.0,360.0}, (double[]){0.0, 10.0});
  }

  return EXIT_SUCCESS;
}
