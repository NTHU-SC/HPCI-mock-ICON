// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdio.h>
#include <string.h>
#ifdef YAC_NETCDF_ENABLED
#include <netcdf.h>
#endif

#include "tests.h"
#include "test_common.h"
#include "io_utils.h"
#include "geometry.h"
#include "test_read_icon_common.h"

#ifdef YAC_NETCDF_ENABLED

int vertex_of_cell[3][16] = {{1,2,6,7,10,11,13,11,8,7,3,2,3,4,4,8},
                             {2,6,7,10,11,13,14,12,11,8,7,3,4,5,8,9},
                             {6,7,10,11,13,14,15,14,12,11,8,7,8,9,9,12}};
int edge_of_cell[3][16] = {{1,4,13,16,22,25,28,24,19,15,7,3,6,9,10,18},
                           {4,13,16,22,25,28,30,27,24,19,15,7,10,12,18,21},
                           {2,5,14,17,23,26,29,26,20,17,8,5,8,11,11,20}};
int cells_of_vertex[6][15] = {{1, 1,11,13,14,1, 2, 9,14,3, 4, 8,5,6,7},
                              {0, 2,12,14, 0,2, 3,10,15,4, 5, 9,6,7,0},
                              {0,12,13,15, 0,3, 4,11,16,5, 6,16,7,8,0},
                              {0, 0, 0, 0, 0,0,10,13, 0,0, 8, 0,0,0,0},
                              {0, 0, 0, 0, 0,0,11,15, 0,0, 9, 0,0,0,0},
                              {0, 0, 0, 0, 0,0,12,16, 0,0,10, 0,0,0,0}};
int vertex_of_edge[2][30] = {{1,1,2,2,2,3,3,3,4,4,4,5, 6,6,7,7,7,8,8,8,9, 10,10,11,11,11,12, 13,13,14},
                             {2,6,3,6,7,4,7,8,5,8,9,9, 7,10,8,10,11,9,11,12,12, 11,13,12,13,14,14, 14,15,15}};


static void def_coord_unit_att(
  int ncid, int var_id, enum coord_units coord_unit) {

  char const * units_str = (coord_unit==DEG)?"degree":"radian";
  YAC_HANDLE_ERROR(
    nc_put_att_text(ncid, var_id, "units", strlen(units_str), units_str));
}

static void write_coord(
  int ncid, int var_id, double * var, size_t var_size,
  enum coord_units coord_unit) {

  if (coord_unit == RAD) {
    double * temp_var = xmalloc(var_size * sizeof(*temp_var));
    for (size_t i = 0; i < var_size; ++i)
      temp_var[i] = var[i] * YAC_RAD;
    var = temp_var;
  }

  YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_id, var));

  if (coord_unit == RAD) free(var);
}

void write_test_grid_file(
  char const * file_name, enum coord_units coord_unit) {

  int ncid;

  // create file
  yac_nc_create(file_name, NC_CLOBBER, &ncid);

  int dim_cell_id, dim_vertex_id, dim_edge_id, dim_nv_id, dim_ne_id, dim_nc_id;

  // define dimensions
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "cell", 16, &dim_cell_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "vertex", 15, &dim_vertex_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "edge", 30, &dim_edge_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "nv", 3, &dim_nv_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "ne", 6, &dim_ne_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "nc", 2, &dim_nc_id));


  int var_vlon_id, var_vlat_id, var_clon_id, var_clat_id, var_mask_id,
      var_v2c_id, var_c2v_id, var_c2e_id, var_v2e_id;

  // define variables
  YAC_HANDLE_ERROR(nc_def_var(ncid, "vlon", NC_DOUBLE, 1, &dim_vertex_id,
                          &var_vlon_id));
  YAC_HANDLE_ERROR(nc_def_var(ncid, "vlat", NC_DOUBLE, 1, &dim_vertex_id,
                          &var_vlat_id));
  YAC_HANDLE_ERROR(nc_def_var(ncid, "clon", NC_DOUBLE, 1, &dim_cell_id,
                          &var_clon_id));
  YAC_HANDLE_ERROR(nc_def_var(ncid, "clat", NC_DOUBLE, 1, &dim_cell_id,
                          &var_clat_id));
  YAC_HANDLE_ERROR(nc_def_var(ncid, "cell_sea_land_mask", NC_INT, 1, &dim_cell_id,
                          &var_mask_id));
  def_coord_unit_att(ncid, var_vlon_id, coord_unit);
  def_coord_unit_att(ncid, var_vlat_id, coord_unit);
  def_coord_unit_att(ncid, var_clon_id, coord_unit);
  def_coord_unit_att(ncid, var_clat_id, coord_unit);
  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, "vertex_of_cell", NC_INT, 2, (int[]){dim_nv_id, dim_cell_id},
      &var_c2v_id));
  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, "edge_of_cell", NC_INT, 2, (int[]){dim_nv_id, dim_cell_id},
      &var_c2e_id));
  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, "cells_of_vertex", NC_INT, 2, (int[]){dim_ne_id, dim_vertex_id},
      &var_v2c_id));
  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, "edge_vertices", NC_INT, 2, (int[]){dim_nc_id, dim_edge_id},
      &var_v2e_id));

  // end definition
  YAC_HANDLE_ERROR(nc_enddef(ncid));

  // write grid data
  write_coord(ncid, var_vlon_id, vlon, sizeof(vlon)/sizeof(vlon[0]), coord_unit);
  write_coord(ncid, var_vlat_id, vlat, sizeof(vlat)/sizeof(vlat[0]), coord_unit);
  write_coord(ncid, var_clon_id, clon, sizeof(clon)/sizeof(clon[0]), coord_unit);
  write_coord(ncid, var_clat_id, clat, sizeof(clat)/sizeof(clat[0]), coord_unit);
  YAC_HANDLE_ERROR(nc_put_var_int(ncid, var_mask_id, mask));
  YAC_HANDLE_ERROR(nc_put_var_int(ncid, var_c2v_id, &(vertex_of_cell[0][0])));
  YAC_HANDLE_ERROR(nc_put_var_int(ncid, var_c2e_id, &(edge_of_cell[0][0])));
  YAC_HANDLE_ERROR(nc_put_var_int(ncid, var_v2c_id, &(cells_of_vertex[0][0])));
  YAC_HANDLE_ERROR(nc_put_var_int(ncid, var_v2e_id, &(vertex_of_edge[0][0])));

  YAC_HANDLE_ERROR(nc_close(ncid));
}

#else // YAC_NETCDF_ENABLED

void write_test_grid_file(
  char const * file_name, enum coord_units coord_unit) {

  UNUSED(file_name);
  UNUSED(coord_unit);
  die("ERROR(write_test_grid_file): YAC is built without the NetCDF support");
}
#endif // YAC_NETCDF_ENABLED

void write_dummy_grid_file(
  char * grid_name, char * grid_filename, char * mask_filename,
  int with_corners, size_t num_lon, size_t num_lat,
  double lon_range[2], double lat_range[2]) {

#ifndef YAC_NETCDF_ENABLED
  UNUSED(grid_name);
  UNUSED(grid_filename);
  UNUSED(mask_filename);
  UNUSED(with_corners);
  UNUSED(num_lon);
  UNUSED(num_lat);
  UNUSED(lon_range);
  UNUSED(lat_range);
  die("ERROR(write_dummy_grid_file): YAC is built without the NetCDF support");
#else

  { // grid file
    int ncid;

    // create file
    yac_nc_create(grid_filename, NC_CLOBBER, &ncid);

    char crn_dim_name[128];
    char x_dim_name[128];
    char y_dim_name[128];

    sprintf(crn_dim_name, "crn_%s", grid_name);
    sprintf(x_dim_name, "x_%s", grid_name);
    sprintf(y_dim_name, "y_%s", grid_name);

    int dim_crn_id = -1;
    int dim_x_id;
    int dim_y_id;

    // define dimensions
    if (with_corners)
      YAC_HANDLE_ERROR(nc_def_dim(ncid, crn_dim_name, 4, &dim_crn_id));
    YAC_HANDLE_ERROR(nc_def_dim(ncid, x_dim_name, num_lon, &dim_x_id));
    YAC_HANDLE_ERROR(nc_def_dim(ncid, y_dim_name, num_lat, &dim_y_id));

    char cla_var_name[128];
    char clo_var_name[128];
    char lat_var_name[128];
    char lon_var_name[128];

    sprintf(cla_var_name, "%s.cla", grid_name);
    sprintf(clo_var_name, "%s.clo", grid_name);
    sprintf(lat_var_name, "%s.lat", grid_name);
    sprintf(lon_var_name, "%s.lon", grid_name);

    int corner_dim_ids[3] = {dim_crn_id, dim_y_id, dim_x_id};
    int cell_dim_ids[2] = {dim_y_id, dim_x_id};

    int var_cla_id = -1;
    int var_clo_id = -1;
    int var_lat_id;
    int var_lon_id;

    char degree[] = "degree";
    char title[] = "This is a reg lon-lat dummy grid";

    // define variable
    if (with_corners) {
      YAC_HANDLE_ERROR(
        nc_def_var(
          ncid, cla_var_name, NC_DOUBLE, 3, corner_dim_ids, &var_cla_id));
      YAC_HANDLE_ERROR(
        nc_put_att_text(ncid, var_cla_id, "units", strlen(degree), degree));
      YAC_HANDLE_ERROR(
        nc_put_att_text(ncid, var_cla_id, "title", strlen(title), title));

      YAC_HANDLE_ERROR(
        nc_def_var(
          ncid, clo_var_name, NC_DOUBLE, 3, corner_dim_ids, &var_clo_id));
      YAC_HANDLE_ERROR(
        nc_put_att_text(ncid, var_clo_id, "units", strlen(degree), degree));
      YAC_HANDLE_ERROR(
        nc_put_att_text(ncid, var_clo_id, "title", strlen(title), title));
    }

    YAC_HANDLE_ERROR(
      nc_def_var(
        ncid, lat_var_name, NC_DOUBLE, 2, cell_dim_ids, &var_lat_id));
    YAC_HANDLE_ERROR(
      nc_put_att_text(ncid, var_lat_id, "units", strlen(degree), degree));
    YAC_HANDLE_ERROR(
      nc_put_att_text(ncid, var_lat_id, "title", strlen(title), title));

    YAC_HANDLE_ERROR(
      nc_def_var(
        ncid, lon_var_name, NC_DOUBLE, 2, cell_dim_ids, &var_lon_id));
    YAC_HANDLE_ERROR(
      nc_put_att_text(ncid, var_lon_id, "units", strlen(degree), degree));
    YAC_HANDLE_ERROR(
      nc_put_att_text(ncid, var_lon_id, "title", strlen(title), title));


    // end definition
    YAC_HANDLE_ERROR(nc_enddef(ncid));

    // write grid data

    double cla[4][num_lat][num_lon];
    double clo[4][num_lat][num_lon];
    double lat[num_lat][num_lon];
    double lon[num_lat][num_lon];

    for (size_t i = 0; i < num_lon; ++i) {
      double vertex_lon[2] =
        {((lon_range[1] - lon_range[0]) * (double)i)/(double)num_lon,
         ((lon_range[1] - lon_range[0]) * (double)(i+1))/(double)num_lon};
      for (size_t j = 0; j < num_lat; ++j) {
        double vertex_lat[2] =
            {((lat_range[1] - lat_range[0]) * (double)j)/(double)num_lat,
            ((lat_range[1] - lat_range[0]) * (double)(j+1))/(double)num_lat};
        cla[0][j][i] = lat_range[0] + vertex_lat[0];
        cla[1][j][i] = lat_range[0] + vertex_lat[0];
        cla[2][j][i] = lat_range[0] + vertex_lat[1];
        cla[3][j][i] = lat_range[0] + vertex_lat[1];
        clo[0][j][i] = lon_range[0] + vertex_lon[0];
        clo[1][j][i] = lon_range[0] + vertex_lon[1];
        clo[2][j][i] = lon_range[0] + vertex_lon[1];
        clo[3][j][i] = lon_range[0] + vertex_lon[0];
        lat[j][i] = lat_range[0] + (vertex_lat[0] + vertex_lat[1]) * 0.5;
        lon[j][i] = lon_range[0] + (vertex_lon[0] + vertex_lon[1]) * 0.5;
      }
    }


    if (with_corners) {
      YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_cla_id, &cla[0][0][0]));
      YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_clo_id, &clo[0][0][0]));
    }
    YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_lat_id, &lat[0][0]));
    YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_lon_id, &lon[0][0]));

    YAC_HANDLE_ERROR(nc_close(ncid));
  }

  { // mask file
    int ncid;

    // create file
    yac_nc_create(mask_filename, NC_CLOBBER, &ncid);

    char x_dim_name[128];
    char y_dim_name[128];

    sprintf(x_dim_name, "x_%s", grid_name);
    sprintf(y_dim_name, "y_%s", grid_name);

    int dim_x_id;
    int dim_y_id;

    // define dimensions
    YAC_HANDLE_ERROR(nc_def_dim(ncid, x_dim_name, num_lon, &dim_x_id));
    YAC_HANDLE_ERROR(nc_def_dim(ncid, y_dim_name, num_lat, &dim_y_id));

    char frc_var_name[128];
    char msk_var_name[128];

    sprintf(frc_var_name, "%s.frc", grid_name);
    sprintf(msk_var_name, "%s.msk", grid_name);

    int dim_ids[2] = {dim_y_id, dim_x_id};

    int var_frc_id;
    int var_msk_id;

    char adim[] = "adim";

    // define variable
    YAC_HANDLE_ERROR(
      nc_def_var(
        ncid, frc_var_name, NC_DOUBLE, 2, dim_ids, &var_frc_id));
    YAC_HANDLE_ERROR(
      nc_put_att_text(ncid, var_frc_id, "units", strlen(adim), adim));

    YAC_HANDLE_ERROR(
      nc_def_var(
        ncid, msk_var_name, NC_INT, 2, dim_ids, &var_msk_id));
    YAC_HANDLE_ERROR(
      nc_put_att_text(ncid, var_msk_id, "units", strlen(adim), adim));


    // end definition
    YAC_HANDLE_ERROR(nc_enddef(ncid));

    // write grid data

    double frc[num_lat][num_lon];
    int msk[num_lat][num_lon];

    for (size_t i = 0; i < num_lon; ++i) {
      for (size_t j = 0; j < num_lat; ++j) {
        frc[j][i] = 1;
        msk[j][i] = 0;
      }
    }

    YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_frc_id, &frc[0][0]));
    YAC_HANDLE_ERROR(nc_put_var_int(ncid, var_msk_id, &msk[0][0]));

    YAC_HANDLE_ERROR(nc_close(ncid));
  }
#endif // YAC_NETCDF_ENABLED
}

void write_test_grid_file_f2c(char const * file_name, int coord_unit) {

  write_test_grid_file(file_name, (enum coord_units)coord_unit);
}
