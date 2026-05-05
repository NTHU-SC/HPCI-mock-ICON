// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TEST_READ_ICON_COMMON_H

/* the grid file contains the following grid
 * (numbers == global cell/vertex indices, XX == masked cell)
 *
 *00-00-01-02-02-05-03-08-04   ------------------------
 * \ 00 /\ 11 /\ 12 /\ 13 /    \    /\ XX /\    /\    /
 * 01 03 04 06 07 09 10 11      \  /  \  /  \  /  \  /
 *   \/ 01 \/ 10 \/ 14 \/        \/    \/ XX \/    \/
 *   05-12-06-14-07-17-08         ------------------
 *    \ 02 /\ 09 /\ 15 /          \    /\ XX /\    /
 *    13 15 16 18 19 20            \  /  \  /  \  /
 *      \/ 03 \/ 08 \/              \/    \/ XX \/
 *      09-21-10-23-11               ------------
 *       \ 04 /\ 07 /                \    /\ XX /
 *       22 24 25 26                  \  /  \  /
 *         \/ 05 \/                    \/    \/
 *         12-27-13                     ------
 *          \ 06 /                      \    /
 *          28 29                        \  /
              \/                          \/
 *            14
 */

static double vlon[15] = {-2,-1,0,1,2, -1.5,-0.5,0.5,1.5, -1,0,1, -0.5,0.5, 0};
static double vlat[15] = {1.5,1.5,1.5,1.5,1.5, 0.5,0.5,0.5,0.5, -0.5,-0.5,-0.5,
                          -1.5,-1.5, -2.5};
static double clon[16] = {1,1,0,0,-1,-1,-2,-1,0,0,1,1,1,1,1,0};
static double clat[16] = {-1.5,-1,-1,-0.5,-0.5,0,0,0.5,0.5,0,0,-0.5,0.5,1,1,1.5};
static int mask[16] = {0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0};

enum coord_units {DEG = 0, RAD = 1};

void write_test_grid_file(
  char const * file_name, enum coord_units coord_unit);


#endif // TEST_READ_ICON_COMMON_H
