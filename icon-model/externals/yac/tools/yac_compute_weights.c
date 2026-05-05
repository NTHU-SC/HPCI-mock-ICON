// Copyright (c) 2025 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <unistd.h>
#include <yaxt.h>
#include <string.h>
#include <netcdf.h>
#include "yac.h"
#include "yac_utils.h"
#include "yac_mpi_common.h"
#include "yac_mpi_internal.h"

// redefine YAC assert macros
#undef YAC_ASSERT
#undef YAC_ASSERT_F
static char const * cmd;
#define DEFAULT_INTERP_STACK \
  "          [\n" \
  "            {\"conservative\":\n" \
  "              {\"order\": 1,\n" \
  "                \"enforced_conservation\": false,\n" \
  "                \"partial_coverage\": false,\n" \
  "                \"normalisation\": \"fracarea\"}\n" \
  "            },\n" \
  "            {\"fixed\":\n" \
  "              {\"user_value\": -1}\n" \
  "            }\n" \
  "          ]\n"
#define STR_USAGE \
  "YAC weight file generation tool\n" \
  "  Reads in or generates a source and target grid file, computes\n" \
  "  interpolation weights, and writes them to file.\n" \
  "  Program is parallelised using MPI and can be run with an arbitrary number\n" \
  "  of processes.\n" \
  "\n" \
  "Usage:\n" \
  "\n" \
  "  mpirun -n $N %s [OPTION]\n" \
  "\n" \
  "  Mandatory arguments:\n" \
  "    -s/-t GRID_TYPE,{CONFIG}\n" \
  "      type and configuraion of source/target grid\n" \
  "      GRID_TYPE:\n" \
  "        \"icon\":   ICON-formated NetCDF grid file\n" \
  "            CONFIG: {FILE,NAME}\n" \
  "              FILE: grid file name\n" \
  "              NAME: grid name\n" \
  "        \"scrip\":  OASIS-SCRIP-formated NetCDF grid file \n" \
  "            CONFIG: {FILE,NAME,EDGE_TYPE,MASK_FILE}\n" \
  "              FILE: grid file name\n" \
  "              NAME: grid name\n" \
  "              EDGE_TYPE: grid edge type\n" \
  "                \"gc\": great circles edge\n" \
  "                \"ll\": lon/lat circles edge\n" \
  "              MASK_FILE: mask file name\n" \
  "        \"exodus\": EXODUS-formated NetCDF grid file\n" \
  "            CONFIG: {FILE,NAME,EDGE_TYPE}\n" \
  "              FILE: grid file name\n" \
  "              NAME: grid name\n" \
  "              EDGE_TYPE: grid edge type\n" \
  "                \"gc\": great circles edge\n" \
  "                \"ll\": lon/lat circles edge\n" \
  "        \"reg2d\": Regular lon/lat grid\n" \
  "                   (edges follow either circles of constant longitude or\n" \
  "                    latitude)\n" \
  "            CONFIG: {NAME,NLON,NLAT,MIN_LON,MAX_LON,MIN_LAT,MAX_LAT}\n" \
  "              NAME: grid name\n" \
  "              NLON: number of grid vertices in longitude direction\n" \
  "              NLAT: number of grid vertices in latitude direction\n" \
  "              MIN_LON: minimum longitude in degree\n" \
  "              MAX_LON: maximum longitude in degree\n" \
  "              MIN_LAT: minimum latitude in degree\n" \
  "              MAX_LAT: maximum latitude in degree\n" \
  "        \"reg2drot\": Regular lon/lat grid with rotated north pole\n" \
  "                       (all edges follow great circles)\n" \
  "            CONFIG: {NAME,NLON,NLAT,MIN_LON,MAX_LON,MIN_LAT,MAX_LAT," \
                       "POL_LON,POL_LAT}\n" \
  "              NAME: grid name\n" \
  "              NLON: number of grid vertices in longitude direction\n" \
  "              NLAT: number of grid vertices in latitude direction\n" \
  "              MIN_LON: minimum longitude in degree\n" \
  "              MAX_LON: maximum longitude in degree\n" \
  "              MIN_LAT: minimum latitude in degree\n" \
  "              MAX_LAT: maximum latitude in degree\n" \
  "              POL_LON: longitude of north in degree\n" \
  "              POL_LAT: latitude of north in in degree\n" \
  "    -o FILE\n" \
  "      weight file name\n" \
  "\n" \
  "  Optional arguments:\n" \
  "    -i {INTERP_STACK}\n" \
  "      JSON-formated interpolation stack configuration " \
         "(see YAC documentation)\n" \
  "      default value:\n" \
  DEFAULT_INTERP_STACK  \
  "    -T\n" \
  "      enables writing of run-time performance measurments to stdout\n" \
  "\n" \
  "  Examples:\n" \
  "\n" \
  "    mpirun -n $N %s -s exodus,CSMesh129.nc,CSMesh129,gc\\\n" \
  "                    -t exodus,ICOMesh100.nc,ICOMesh100,gc\\\n" \
  "                    -o CSMesh129_to_ICOMesh100.nc\n"\
  "\n" \
  "    mpirun -n $N %s -s scrip,grids.nc,torc,gc,masks.nc\\\n" \
  "                    -t icon,icon_grid_R02B05.nc,icon_R02B05\\\n" \
  "                    -o torc_to_iconR02B05.nc -T -i \"[\\\"nnn\\\"]\"\n"

#define YAC_ASSERT(exp, msg) \
  { \
    if(!((exp))) { \
      fprintf(stderr, "ERROR: %s\n" STR_USAGE, msg, cmd, cmd, cmd); \
      exit(EXIT_FAILURE); \
    } \
  }

#define YAC_ASSERT_F(exp, format, ...) \
  { \
    if(!((exp))) { \
      fprintf( \
        stderr, "ERROR: " format "\n\n" STR_USAGE, \
        __VA_ARGS__, cmd, cmd, cmd); \
      exit(EXIT_FAILURE); \
    } \
  }

enum grid_edge_type {GC_EDGES, LL_EDGES};

struct grid_config {
  enum grid_type {
    EXODUS,
    ICON,
    SCRIP,
    REG2D,
    REG2DROT,
    UNDEFINED_GRID
  } type;
  union {
    struct {
      char const * grid_filename;
      enum grid_edge_type edge_type;
      size_t cell_coordinate_idx;
    } exodus;
    struct {
      char const * grid_filename;
      size_t cell_coordinate_idx;
    } icon;
    struct {
      char const * grid_filename;
      char const * mask_filename;
      enum grid_edge_type edge_type;
      size_t * duplicated_cell_idx;
      yac_int * orig_cell_global_ids;
      size_t nbr_duplicated_cells;
    } scrip;
    struct {
      size_t nlon, nlat;
      double min_lon, max_lon, min_lat, max_lat;
    } reg2d;
    struct {
      size_t nlon, nlat;
      double min_lon, max_lon, min_lat, max_lat;
      double pole_lon, pole_lat;
    } reg2drot;
  } data;
  size_t cell_coordinate_idx;
  char const * grid_name;
  size_t global_num_cells;
};

static void parse_arguments(int argc, char ** argv,
                            struct grid_config * src_grid_config,
                            struct grid_config * tgt_grid_config,
                            char const ** weight_filename,
                            char const ** interp_stack_string,
                            int * print_timer,
                            char const ** debug_grid_file);

extern int const YAC_YAML_PARSER_JSON_FORCE;
struct yac_interp_stack_config *
  yac_yaml_parse_interp_stack_config_string(
    char const * interp_stack_config, int parse_flags);

static struct yac_basic_grid * get_basic_grid_from_config(
  struct grid_config * grid_config, char const * debug_grid_file);

static void grid_config_delete(struct grid_config grid_config);

struct time_rank {
  double time;
  int rank;
};
static struct time_rank local_time_rank;
static void timer_start(int print_timer);
static void timer_stop(int print_timer, char const * timer_name);

int main (int argc, char *argv[]) {

  MPI_Init(&argc, &argv);
  xt_initialize(MPI_COMM_WORLD);

  yac_mpi_call(
    MPI_Comm_rank(MPI_COMM_WORLD, &local_time_rank.rank), MPI_COMM_WORLD);

  cmd = argv[0];

  // parse command line arguments
  char const * weight_filename;
  char const * interp_stack_string;
  struct grid_config src_grid_config;
  struct grid_config tgt_grid_config;
  int print_timer;
  char const * debug_grid_file;
  parse_arguments(
    argc, argv, &src_grid_config, &tgt_grid_config, &weight_filename,
    &interp_stack_string, &print_timer, &debug_grid_file);

  // read grid data
  timer_start(print_timer);
  struct yac_basic_grid * src_grid =
    get_basic_grid_from_config(&src_grid_config, debug_grid_file);
  timer_stop(print_timer, "read src_grid");
  timer_start(print_timer);
  struct yac_basic_grid * tgt_grid =
    get_basic_grid_from_config(&tgt_grid_config, debug_grid_file);
  timer_stop(print_timer, "read tgt_grid");

  // generate distributed grid pair
  timer_start(print_timer);
  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);
  timer_stop(print_timer, "dist_grid_pair generation");

  // setup field information
  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CELL,
      .coordinates_idx = src_grid_config.cell_coordinate_idx,
      .masks_idx = SIZE_MAX}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CELL,
     .coordinates_idx = src_grid_config.cell_coordinate_idx,
     .masks_idx = SIZE_MAX};

  // generate interpolation grid
  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(
      grid_pair, yac_basic_grid_get_name(src_grid),
      yac_basic_grid_get_name(tgt_grid),
      num_src_fields, src_fields, tgt_field);

  // generate interpolation stack configuration from JSON-formated string
  struct yac_interp_stack_config * interp_stack_config =
    yac_yaml_parse_interp_stack_config_string(
      interp_stack_string, YAC_YAML_PARSER_JSON_FORCE);

  // generate interpolation stack
  struct interp_method ** interp_stack =
    yac_interp_stack_config_generate(interp_stack_config);

  // compute interpolation weights
  timer_start(print_timer);
  struct yac_interp_weights * weights =
    yac_interp_method_do_search(interp_stack, interp_grid);
  timer_stop(print_timer, "weight computation");

  // OASIS SCRIP formated grid files may contain duplicated cell, these are
  // masked out in the basic grid, but still require a interpolation stencil
  if (tgt_grid_config.type == SCRIP) {
    yac_duplicate_stencils(
      weights, tgt_grid, tgt_grid_config.data.scrip.orig_cell_global_ids,
      tgt_grid_config.data.scrip.duplicated_cell_idx,
      tgt_grid_config.data.scrip.nbr_duplicated_cells, YAC_LOC_CELL);
  }

  // write weights to file
  timer_start(print_timer);
  enum yac_weight_file_on_existing on_existing = YAC_WEIGHT_FILE_OVERWRITE;
  yac_interp_weights_write_to_file(
    weights, weight_filename, yac_basic_grid_get_name(src_grid),
    yac_basic_grid_get_name(tgt_grid), src_grid_config.global_num_cells,
    tgt_grid_config.global_num_cells, on_existing);
  timer_stop(print_timer, "writing weight file");

  // cleanup
  grid_config_delete(tgt_grid_config);
  grid_config_delete(src_grid_config);
  yac_interp_weights_delete(weights);
  yac_interp_method_delete(interp_stack);
  free(interp_stack);
  yac_interp_stack_config_delete(interp_stack_config);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);
  yac_basic_grid_delete(tgt_grid);
  yac_basic_grid_delete(src_grid);

  xt_finalize();
  MPI_Finalize();

  return EXIT_SUCCESS;
}

static enum grid_edge_type parse_edge_type(char const * edge_type_string) {

  int is_gc = !strcmp(edge_type_string, "gc") ||
              !strcmp(edge_type_string, "GC");
  int is_ll = !strcmp(edge_type_string, "ll") ||
              !strcmp(edge_type_string, "LL");
  YAC_ASSERT_F(
    is_gc || is_ll, "invalid grid edge type (\"%s\")", edge_type_string);
  free((void*)edge_type_string);

  return is_gc?GC_EDGES:LL_EDGES;
}

static size_t parse_size_t(char const * size_t_string) {

  char * endptr;
  long int long_value = strtol(size_t_string, &endptr, 10);

  YAC_ASSERT_F(
    (endptr != size_t_string) && (*endptr == '\0') && (long_value >= 0),
    "\"%s\" is not a valid size_t value", size_t_string);
  free((void*)size_t_string);

  return (size_t)long_value;
}

static double parse_double(char const * double_string) {

  char * endptr;
  double dble_value = strtod(double_string, &endptr);

  YAC_ASSERT_F(
    (endptr != double_string) && (*endptr == '\0'),
    "\"%s\" is not a valid double value", double_string);
  free((void*)double_string);

  return dble_value;
}

static char const * get_next_token(char const * token_name) {
  char const * token = strtok(NULL, ",");
  YAC_ASSERT_F(token, "missing %s", token_name);
  return strdup(token);
}

size_t read_netcdf_dimension(char const * filename, char const * dim_name) {

  int rank;
  yac_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &rank), MPI_COMM_WORLD);

  size_t dimlen;

  if (rank == 0) {
    int ncid, dimid;
    yac_nc_open(filename, NC_NOWRITE, &ncid);
    yac_nc_inq_dimid(ncid, dim_name, &dimid);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, &dimlen));
  }

  yac_mpi_call(
    MPI_Bcast(&dimlen, 1, YAC_MPI_SIZE_T, 0, MPI_COMM_WORLD), MPI_COMM_WORLD);

  return dimlen;
}

struct grid_config parse_grid_config_exodus() {

  char const * grid_filename = get_next_token("exodus grid filename");
  char const * grid_name = get_next_token("exodus grid name");
  enum grid_edge_type edge_type =
    parse_edge_type(get_next_token("exodus edge type"));
  size_t global_num_cells = read_netcdf_dimension(grid_filename, "num_elem");
  return
    (struct grid_config) {
      .type = EXODUS,
      .grid_name = grid_name,
      .global_num_cells = global_num_cells,
      .data.exodus.grid_filename = grid_filename,
      .data.exodus.edge_type = edge_type
    };
}

struct grid_config parse_grid_config_icon() {

  char const * grid_filename = get_next_token("icon grid filename");
  char const * grid_name = get_next_token("icon grid name");
  size_t global_num_cells = read_netcdf_dimension(grid_filename, "cell");
  return
    (struct grid_config) {
      .type = ICON,
      .grid_name = grid_name,
      .global_num_cells = global_num_cells,
      .data.icon.grid_filename = grid_filename
    };
}

struct grid_config parse_grid_config_scrip() {

  char const * grid_filename = get_next_token("scrip grid filename");
  char const * grid_name = get_next_token("scrip grid name");
  enum grid_edge_type edge_type =
    parse_edge_type(get_next_token("scrip edge type"));
  char const * mask_filename = get_next_token("scrip mask filename");
  char * x_dim_name = malloc(strlen(grid_name) + 3);
  char * y_dim_name = malloc(strlen(grid_name) + 3);
  strcpy(x_dim_name, "x_"), strcat(x_dim_name, grid_name);
  strcpy(y_dim_name, "y_"), strcat(y_dim_name, grid_name);
  size_t global_num_cells =
    read_netcdf_dimension(grid_filename, x_dim_name) *
    read_netcdf_dimension(grid_filename, y_dim_name);
  free(y_dim_name);
  free(x_dim_name);
  return
    (struct grid_config) {
      .type = SCRIP,
      .grid_name = grid_name,
      .global_num_cells = global_num_cells,
      .data.scrip.grid_filename = grid_filename,
      .data.scrip.edge_type = edge_type,
      .data.scrip.mask_filename = mask_filename
    };
}

struct grid_config parse_grid_config_reg2d() {

  char const * grid_name = get_next_token("reg2d grid name");
  size_t nlon = parse_size_t(get_next_token("reg2d number of vertices in lon"));
  size_t nlat = parse_size_t(get_next_token("reg2d number of vertices in lat"));
  double min_lon = parse_double(get_next_token("reg2d minimum longitude"));
  double max_lon = parse_double(get_next_token("reg2d maximum longitude"));
  double min_lat = parse_double(get_next_token("reg2d minimum latitude"));
  double max_lat = parse_double(get_next_token("reg2d maximum latitude"));
  return
    (struct grid_config) {
      .type = REG2D,
      .grid_name = grid_name,
      .global_num_cells = (nlon - 1) * (nlat - 1),
      .data.reg2d.nlon = nlon,
      .data.reg2d.nlat = nlat,
      .data.reg2d.min_lon = min_lon,
      .data.reg2d.max_lon = max_lon,
      .data.reg2d.min_lat = min_lat,
      .data.reg2d.max_lat = max_lat
    };
}

struct grid_config parse_grid_config_reg2drot() {

  char const * grid_name = get_next_token("reg2d grid name");
  size_t nlon = parse_size_t(get_next_token("reg2d number of vertices in lon"));
  size_t nlat = parse_size_t(get_next_token("reg2d number of vertices in lat"));
  double min_lon = parse_double(get_next_token("reg2d minimum longitude"));
  double max_lon = parse_double(get_next_token("reg2d maximum longitude"));
  double min_lat = parse_double(get_next_token("reg2d minimum latitude"));
  double max_lat = parse_double(get_next_token("reg2d maximum latitude"));
  double pole_lon = parse_double(get_next_token("reg2d pole longitude"));
  double pole_lat = parse_double(get_next_token("reg2d pole latitude"));
  return
    (struct grid_config) {
      .type = REG2DROT,
      .grid_name = grid_name,
      .global_num_cells = (nlon - 1) * (nlat - 1),
      .data.reg2drot.nlon = nlon,
      .data.reg2drot.nlat = nlat,
      .data.reg2drot.min_lon = min_lon,
      .data.reg2drot.max_lon = max_lon,
      .data.reg2drot.min_lat = min_lat,
      .data.reg2drot.max_lat = max_lat,
      .data.reg2drot.pole_lon = pole_lon,
      .data.reg2drot.pole_lat = pole_lat
    };
}

static struct grid_config parse_grid_config(
  char const * grid_config_string_, char const * src_tgt) {

  struct grid_config grid_config = {.type = UNDEFINED_GRID};

  char * grid_config_string = strdup(grid_config_string_);
  char const * grid_type_string = strtok(grid_config_string, ",");
  if (!strcmp("exodus", grid_type_string))
    grid_config = parse_grid_config_exodus();
  if (!strcmp("icon", grid_type_string))
    grid_config = parse_grid_config_icon();
  if (!strcmp("scrip", grid_type_string))
    grid_config = parse_grid_config_scrip();
  if (!strcmp("reg2d", grid_type_string))
    grid_config = parse_grid_config_reg2d();
  if (!strcmp("reg2drot", grid_type_string))
    grid_config = parse_grid_config_reg2drot();

  free(grid_config_string);

  YAC_ASSERT_F(
    grid_config.type != UNDEFINED_GRID,
    "invalid %s grid type (\"%s\")", src_tgt, grid_type_string);

  return grid_config;
}

static void parse_arguments(int argc, char ** argv,
                            struct grid_config * src_grid_config,
                            struct grid_config * tgt_grid_config,
                            char const ** weight_filename,
                            char const ** interp_stack_string,
                            int * print_timer,
                            char const ** debug_grid_file) {

  src_grid_config->type = UNDEFINED_GRID;
  tgt_grid_config->type = UNDEFINED_GRID;
  *weight_filename = NULL;
  *interp_stack_string = DEFAULT_INTERP_STACK;
  *print_timer = 0;
  *debug_grid_file = NULL;

  int opt;
  while ((opt = getopt(argc, argv, "s:t:o:i:d:T")) != -1) {
    YAC_ASSERT(
      (opt == 's') ||
      (opt == 't') ||
      (opt == 'o') ||
      (opt == 'i') ||
      (opt == 'd') ||
      (opt == 'T'), "invalid command argument")
    switch (opt) {
      default:
      case 's':
        YAC_ASSERT(
          src_grid_config->type == UNDEFINED_GRID,
          "multiple source grid arguments")
        *src_grid_config = parse_grid_config(optarg, "source");
        break;
      case 't':
        YAC_ASSERT(
          tgt_grid_config->type == UNDEFINED_GRID,
          "multiple target grid arguments")
        *tgt_grid_config = parse_grid_config(optarg, "target");
        break;
      case 'o':
        *weight_filename = optarg;
        break;
      case 'i':
        *interp_stack_string = optarg;
        break;
      case 'd':
        *debug_grid_file = optarg;
        break;
      case 'T':
        *print_timer = 1;
        break;
    }
  }
  YAC_ASSERT_F(
    optind >= argc, "non-option ARGV-element: \"%s\"", argv[optind])
  YAC_ASSERT(argc != 1, "too few arguments")
  YAC_ASSERT(
    src_grid_config->type != UNDEFINED_GRID,  "source grid argument is missing")
  YAC_ASSERT(
    tgt_grid_config->type != UNDEFINED_GRID,  "target grid argument is missing")
  YAC_ASSERT(*weight_filename != NULL, "weight_filename argument is missing")
}

static inline void normalise_vector(double v[]) {

   double norm = 1.0 / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

   v[0] *= norm;
   v[1] *= norm;
   v[2] *= norm;
}

static size_t generate_cell_center_coordinates(
  struct yac_basic_grid * basic_grid) {

  struct yac_basic_grid_data * basic_grid_data =
    yac_basic_grid_get_data(basic_grid);

  size_t coordinates_idx;

  if (basic_grid_data->num_cells == 0) {

    coordinates_idx = SIZE_MAX;

  } else {

    yac_coordinate_pointer cell_center_coords =
      malloc(basic_grid_data->num_cells * sizeof(*cell_center_coords));

    for (size_t i = 0; i < basic_grid_data->num_cells; ++i) {

      double cell_center_coord[3] = {0.0, 0.0, 0.0};
      size_t * curr_cell_to_vertex =
        basic_grid_data->cell_to_vertex +
        basic_grid_data->cell_to_vertex_offsets[i];

      for (int j = 0; j < basic_grid_data->num_vertices_per_cell[i]; ++j) {

        double * curr_vertex =
          basic_grid_data->vertex_coordinates[curr_cell_to_vertex[j]];
        cell_center_coord[0] += curr_vertex[0];
        cell_center_coord[1] += curr_vertex[1];
        cell_center_coord[2] += curr_vertex[2];
      }
      normalise_vector(cell_center_coord);
      memcpy(cell_center_coords[i], cell_center_coord, 3 * sizeof(double));
    }
    coordinates_idx =
      yac_basic_grid_add_coordinates_nocpy(
        basic_grid, YAC_LOC_CELL, cell_center_coords);
  }

  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, &coordinates_idx, 1,
      YAC_MPI_SIZE_T, MPI_MIN, MPI_COMM_WORLD), MPI_COMM_WORLD);

  return coordinates_idx;
}

static void generate_reg2d_vertices(
  size_t * nbr_vertices, double min_lon, double max_lon,
  double min_lat, double max_lat,
  double ** lon_vertices, double ** lat_vertices) {

  *lon_vertices = malloc(nbr_vertices[0] * sizeof(**lon_vertices));
  *lat_vertices = malloc(nbr_vertices[1] * sizeof(**lat_vertices));

  double lon_diff = max_lon - min_lon;
  double lat_diff = max_lat - min_lat;

  for (size_t i = 0; i < nbr_vertices[0]; ++i)
    (*lon_vertices)[i] =
      min_lon + (lon_diff * (double)i) / (double)(nbr_vertices[0] - 1);
  (*lon_vertices)[nbr_vertices[0]-1] = max_lon;

  for (size_t i = 0; i < nbr_vertices[1]; ++i)
    (*lat_vertices)[i] =
      min_lat + (lat_diff * (double)i) / (double)(nbr_vertices[1] - 1);
  (*lat_vertices)[nbr_vertices[1]-1] = max_lat;
}

static struct yac_basic_grid * generate_reg2d_grid(
  struct grid_config grid_config) {

  int rank;
  yac_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &rank), MPI_COMM_WORLD);

  struct yac_basic_grid * basic_grid;

  if (rank == 0) {

    size_t nbr_vertices[2] =
      {grid_config.data.reg2d.nlon, grid_config.data.reg2d.nlat};
    double * lon_vertices, * lat_vertices;
    generate_reg2d_vertices(
      nbr_vertices,
      grid_config.data.reg2d.min_lon, grid_config.data.reg2d.max_lon,
      grid_config.data.reg2d.min_lat, grid_config.data.reg2d.max_lat,
      &lon_vertices, &lat_vertices);

    int cyclic[2] = {0,0};

    basic_grid =
      yac_basic_grid_reg_2d_deg_new(
        grid_config.grid_name, nbr_vertices, cyclic,
        lon_vertices, lat_vertices);

    free(lon_vertices);
    free(lat_vertices);

  } else {

    basic_grid = yac_basic_grid_empty_new(grid_config.grid_name);

  }

  return basic_grid;
}

static struct yac_basic_grid * generate_reg2drot_grid(
  struct grid_config grid_config) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  struct yac_basic_grid * basic_grid;

  if (rank == 0) {

    size_t nbr_vertices[2] =
      {grid_config.data.reg2drot.nlon, grid_config.data.reg2drot.nlat};
    double * lon_vertices, * lat_vertices;
    generate_reg2d_vertices(
      nbr_vertices,
      grid_config.data.reg2drot.min_lon, grid_config.data.reg2drot.max_lon,
      grid_config.data.reg2drot.min_lat, grid_config.data.reg2drot.max_lat,
      &lon_vertices, &lat_vertices);

    int cyclic[2] = {0,0};

    basic_grid =
      yac_basic_grid_reg_2d_rot_deg_new(
        grid_config.grid_name, nbr_vertices, cyclic, lon_vertices, lat_vertices,
        grid_config.data.reg2drot.pole_lon, grid_config.data.reg2drot.pole_lat);

    free(lon_vertices);
    free(lat_vertices);

  } else {

    basic_grid = yac_basic_grid_empty_new(grid_config.grid_name);

  }

  return basic_grid;
}

static struct yac_basic_grid * get_basic_grid_from_config(
  struct grid_config * grid_config, char const * debug_grid_file) {

  struct yac_basic_grid * basic_grid = NULL;

  YAC_ASSERT(
    (grid_config->type == EXODUS) ||
    (grid_config->type == ICON) ||
    (grid_config->type == SCRIP) ||
    (grid_config->type == REG2D) ||
    (grid_config->type == REG2DROT), "invalid grid type");

  switch (grid_config->type) {
    default:
    case (EXODUS): {
      int use_ll_edges = grid_config->data.exodus.edge_type == LL_EDGES;
      basic_grid =
        yac_read_exodus_basic_grid_parallel(
          grid_config->data.exodus.grid_filename, grid_config->grid_name,
          use_ll_edges, MPI_COMM_WORLD);
      grid_config->cell_coordinate_idx =
        generate_cell_center_coordinates(basic_grid);
      break;
    }
    case (ICON): {
      yac_read_icon_basic_grid_parallel_2(
        grid_config->data.icon.grid_filename, grid_config->grid_name,
        MPI_COMM_WORLD, &basic_grid, &grid_config->cell_coordinate_idx, NULL);
      break;
    }
    case (SCRIP): {
      int valid_mask_value = 0;
      int use_ll_edges = grid_config->data.scrip.edge_type == LL_EDGES;
      basic_grid =
        yac_read_scrip_basic_grid_parallel(
          grid_config->data.scrip.grid_filename,
          grid_config->data.scrip.mask_filename,
          MPI_COMM_WORLD, grid_config->grid_name, valid_mask_value,
          grid_config->grid_name, use_ll_edges,
          &grid_config->cell_coordinate_idx,
          &grid_config->data.scrip.duplicated_cell_idx,
          &grid_config->data.scrip.orig_cell_global_ids,
          &grid_config->data.scrip.nbr_duplicated_cells);
      break;
    }
    case (REG2D): {
      basic_grid = generate_reg2d_grid(*grid_config);
      grid_config->cell_coordinate_idx =
        generate_cell_center_coordinates(basic_grid);
      break;
    }
    case (REG2DROT): {
      basic_grid = generate_reg2drot_grid(*grid_config);
      grid_config->cell_coordinate_idx =
        generate_cell_center_coordinates(basic_grid);
      break;
    }
  }

  if (debug_grid_file != NULL)
    yac_basic_grid_to_file_parallel(
      basic_grid, debug_grid_file, MPI_COMM_WORLD);

  return basic_grid;
}

static void grid_config_delete(struct grid_config grid_config) {

  YAC_ASSERT(
    (grid_config.type == EXODUS) ||
    (grid_config.type == ICON) ||
    (grid_config.type == SCRIP) ||
    (grid_config.type == REG2D) ||
    (grid_config.type == REG2DROT), "invalid grid type");

  free((void*)grid_config.grid_name);

  switch (grid_config.type) {
    default:
    case (EXODUS): {
      free((void*)grid_config.data.exodus.grid_filename);
      break;
    }
    case (ICON): {
      free((void*)grid_config.data.icon.grid_filename);
      break;
    }
    case (SCRIP): {
      free((void*)grid_config.data.scrip.grid_filename);
      free((void*)grid_config.data.scrip.mask_filename);
      free(grid_config.data.scrip.duplicated_cell_idx);
      free(grid_config.data.scrip.orig_cell_global_ids);
      break;
    }
    case (REG2D): {
      break;
    }
    case (REG2DROT): {
      break;
    }
  }
}

static void timer_start(int print_timer) {

  if (!print_timer) return;

  yac_mpi_call(MPI_Barrier(MPI_COMM_WORLD), MPI_COMM_WORLD);
  local_time_rank.time = MPI_Wtime();
}

static void timer_stop(int print_timer, char const * timer_name) {

  if (!print_timer) return;

  local_time_rank.time = MPI_Wtime() - local_time_rank.time;
  struct time_rank time_rank_min, time_rank_max;
  double time_sum;
  yac_mpi_call(
    MPI_Reduce(
      &local_time_rank, &time_rank_min, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0,
      MPI_COMM_WORLD), MPI_COMM_WORLD);
  yac_mpi_call(
    MPI_Reduce(
      &local_time_rank, &time_rank_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0,
      MPI_COMM_WORLD), MPI_COMM_WORLD);
  yac_mpi_call(
    MPI_Reduce(
      &local_time_rank.time, &time_sum, 1, MPI_DOUBLE, MPI_SUM, 0,
      MPI_COMM_WORLD), MPI_COMM_WORLD);

  int comm_size;
  yac_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);
  if (local_time_rank.rank == 0)
    fprintf(stdout, "%s: min %.3lfs (%d) avg %.3lfs max %.3lfs (%d)\n",
            timer_name, time_rank_min.time, time_rank_min.rank,
            time_sum / (double)comm_size,
            time_rank_max.time, time_rank_max.rank);
}
