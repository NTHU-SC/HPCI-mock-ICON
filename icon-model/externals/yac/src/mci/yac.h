// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef YAC_H
#define YAC_H

#include <stddef.h>

#include <mpi.h>

extern int const YAC_LOCATION_CELL;
extern int const YAC_LOCATION_CORNER;
extern int const YAC_LOCATION_EDGE;

extern int const YAC_EXCHANGE_TYPE_NONE;
extern int const YAC_EXCHANGE_TYPE_SOURCE;
extern int const YAC_EXCHANGE_TYPE_TARGET;

extern int const YAC_ACTION_NONE;            //!< no data exchanges
extern int const YAC_ACTION_REDUCTION;       //!< data reduction, but data exchange
extern int const YAC_ACTION_COUPLING;        //!< data exchange
extern int const YAC_ACTION_GET_FOR_RESTART; //!< last valid get
extern int const YAC_ACTION_PUT_FOR_RESTART; //!< last valid put
extern int const YAC_ACTION_OUT_OF_BOUND;    //!< put/get is outside of the valid range

extern int const YAC_REDUCTION_TIME_NONE;
extern int const YAC_REDUCTION_TIME_ACCUMULATE;
extern int const YAC_REDUCTION_TIME_AVERAGE;
extern int const YAC_REDUCTION_TIME_MINIMUM;
extern int const YAC_REDUCTION_TIME_MAXIMUM;

extern int const YAC_TIME_UNIT_MILLISECOND;
extern int const YAC_TIME_UNIT_SECOND;
extern int const YAC_TIME_UNIT_MINUTE;
extern int const YAC_TIME_UNIT_HOUR;
extern int const YAC_TIME_UNIT_DAY;
extern int const YAC_TIME_UNIT_MONTH;
extern int const YAC_TIME_UNIT_YEAR;
extern int const YAC_TIME_UNIT_ISO_FORMAT;

extern int const YAC_CALENDAR_NOT_SET;
extern int const YAC_PROLEPTIC_GREGORIAN;
extern int const YAC_YEAR_OF_365_DAYS;
extern int const YAC_YEAR_OF_360_DAYS;

extern int const YAC_AVG_ARITHMETIC;
extern int const YAC_AVG_DIST;
extern int const YAC_AVG_BARY;

extern int const YAC_NCC_AVG;
extern int const YAC_NCC_DIST;

extern int const YAC_NNN_AVG;
extern int const YAC_NNN_DIST;
extern int const YAC_NNN_GAUSS;
extern int const YAC_NNN_RBF;
extern int const YAC_NNN_ZERO;

extern int const YAC_CONSERV_DESTAREA;
extern int const YAC_CONSERV_FRACAREA;

extern int const YAC_SPMAP_AVG;
extern int const YAC_SPMAP_DIST;

extern int const YAC_SPMAP_NONE;
extern int const YAC_SPMAP_SRCAREA;
extern int const YAC_SPMAP_INVTGTAREA;
extern int const YAC_SPMAP_FRACAREA;

extern int const YAC_FILE_MISSING_ERROR;
extern int const YAC_FILE_MISSING_CONT;

extern int const YAC_FILE_SUCCESS_STOP;
extern int const YAC_FILE_SUCCESS_CONT;

extern int const YAC_YAML_EMITTER_DEFAULT; //!<  emit to YAML format
extern int const YAC_YAML_EMITTER_JSON;    //!<  emit to JSON format

extern int const YAC_CONFIG_OUTPUT_FORMAT_YAML;
extern int const YAC_CONFIG_OUTPUT_FORMAT_JSON;

extern int const YAC_CONFIG_OUTPUT_SYNC_LOC_DEF_COMP; //!< after component definition
extern int const YAC_CONFIG_OUTPUT_SYNC_LOC_SYNC_DEF; //!< after synchronization of definition
extern int const YAC_CONFIG_OUTPUT_SYNC_LOC_ENDDEF;   //!< after end of definitions

extern double const YAC_FRAC_MASK_NO_VALUE;

extern int const YAC_WGT_ON_EXISTING_ERROR;     //!< generate error and abort, if
                                                //!< file with the same name already
                                                //!< existis
extern int const YAC_WGT_ON_EXISTING_KEEP;      //!< keep existing weight file
extern int const YAC_WGT_ON_EXISTING_OVERWRITE; //!< overwrite existing weight file

#define YAC_MAX_CHARLEN (132)

/* -----------------------------------------------------------------------------

    C and Fortran description of the user API

   ----------------------------------------------------------------------------- */

/** MPI Handshake (\ref mpi_handshake_detail)

    Splits the provided communicator into group communicators. Each group
    communicator contains all processes of the provided communicator that
    provided the same group name. The order of the group names can be
    arbitrary on each process. A process can be part of multiple groups.

    @param[in]  comm        MPI communicator
    @param[in]  n           number of group communicators to be generated
    @param[in]  group_names group names
    @param[out] group_comms group communicators

    @remark This call is collective for processes in comm.
    @remark \ref yac_cinit will call this routine for `MPI_COMM_WORLD`
            and use a fixed \ref mpi_handshake_group_name to build
            its internal communicator.
    @remark If a process does not wish to use YAC but is part of a
            communicator (e.g. `MPI_COMM_WORLD`) that is used by other
            processes that actually do use YAC, this routine can be
            called, while not providing \ref mpi_handshake_group_name.
            This will exclude a process from later collective YAC-calls.
    @remark \ref yac_cget_mpi_handshake_group_name allows to query the
            \ref mpi_handshake_group_name internally used by YAC.
*/
     void yac_cmpi_handshake( MPI_Comm comm,
                              size_t n,
                              char const** group_names,
                              MPI_Comm * group_comms );

/** Getter function for the default instance

    Returns the instance id of the default instance.

    @remark The default instance needs to be initialised before calling
            this function.
    @returns The instance id of the default instance

 */
     int yac_cget_default_instance_id();

/** Elementary initialisation of the whole system

     @remark This call initialises the default YAC instance.
     @remark A call to \ref yac_ccleanup or \ref yac_cfinalize is required to
             free all memory associated with the default YAC instance.
     @remark This call executes an MPI handshake (\ref mpi_handshake_detail) on
             `MPI_COMM_WORLD` with the group name \ref mpi_handshake_group_name
             and uses the created communicator for the initialization. If you
             do not want to execute an MPI Handshake on `MPI_COMM_WORLD`, call
             \ref yac_cinit_comm with `MPI_COMM_WORLD` directly.
*/
    void yac_cinit (void);

/** Elementary initialisation of the whole system

     @param[out] yac_instance_id id of the YAC instance initialised by
                                 this call
     @remark This call initialises a YAC instance.
     @remark A call to \ref yac_ccleanup_instance or
             \ref yac_cfinalize_instance is required to free all memory
             associated with this YAC instance.
     @remark This call executes an MPI handshake (\ref mpi_handshake_detail) on
             `MPI_COMM_WORLD` with the group name \ref mpi_handshake_group_name
             and uses the created communicator for the initialization. If you
             do not want to execute an MPI Handshake on `MPI_COMM_WORLD`, call
             \ref yac_cinit_comm with `MPI_COMM_WORLD` directly.
*/
    void yac_cinit_instance ( int * yac_instance_id );

/** Elementary initialisation of the whole system using a user-provided
 *  MPI communicator

     @param[in] comm            MPI communicator
     @remark This call initialises the default YAC instance.
     @remark A call to \ref yac_ccleanup or \ref yac_cfinalize is required to
             free all memory associated with the default YAC instance.
     @remark the MPI communicator provided to this routine has to contain
             all processes that will take part in the coupling
     @remark this call is collective for all processes in comm
*/
    void yac_cinit_comm ( MPI_Comm comm );

/** Elementary initialisation of the whole system using a user-provided
 *  MPI communicator

     @param[in]  comm            MPI communicator
     @param[out] yac_instance_id id of the YAC instance initialised by
                                 this call
     @remark This call initialises a YAC instance.
     @remark A call to \ref yac_ccleanup_instance or
             \ref yac_cfinalize_instance is required to free all memory
             associated with this YAC instance.
     @remark In case the user has multiple YAC instances in parallel, he has
             to initialise yaxt himself.
*/
    void yac_cinit_comm_instance ( MPI_Comm comm,
                                   int * yac_instance_id );

/* -------------------------------------------------------------------------------- */

/** Dummy for initialisation of the whole system
     @remark This routine can be called instead of \ref yac_cinit or
             \ref yac_cinit_instance, if the local process does not wish to
             initialise a YAC instance.
*/
    void yac_cinit_dummy(void);

/** Dummy initialisation of the whole system using a user-provided
 *  MPI world communicator

     @param[in] comm MPI communicator
     @remark This routine can be called instead of \ref yac_cinit_comm or
             \ref yac_cinit_comm_instance, if the local process does not wish to
             initialise a YAC instance.
*/
    void yac_cinit_comm_dummy ( MPI_Comm comm );

/* -------------------------------------------------------------------------------- */

/** Function for reading configuration from YAML configuration file for a specific
 *  YAC instance

     @param[in] yac_instance_id YAC instance id
     @param[in] yaml_file       YAML configuration file
*/
     void yac_cread_config_yaml_instance( int yac_instance_id,
                                          const char * yaml_file);

/** Function for reading configuration from YAML configuration file

     @param[in] yaml_file       YAML configuration file
*/
     void yac_cread_config_yaml(const char * yaml_file);

/** Function for reading configuration from JSON configuration file for a specific
 *  YAC instance

     @param[in] yac_instance_id YAC instance id
     @param[in] json_file       JSON configuration file
*/
     void yac_cread_config_json_instance( int yac_instance_id,
                                          const char * json_file);

/** Function for reading configuration from JSON configuration file

     @param[in] json_file       JSON configuration file
*/
     void yac_cread_config_json( const char * json_file );

/* -------------------------------------------------------------------------------- */

/** Activates writing out of the coupling configuration

     @param[in] yac_instance_id     YAC instance_id
     @param[in] filename            name of the file to be written
     @param[in] fileformat          file format (YAC_CONFIG_OUTPUT_FORMAT_YAML
                                    or YAC_CONFIG_OUTPUT_FORMAT_JSON)
     @param[in] sync_location       synchronisation point after which the file
                                    is to be written
                                    (YAC_CONFIG_OUTPUT_SYNC_LOC_DEF_COMP,
                                    YAC_CONFIG_OUTPUT_SYNC_LOC_SYNC_DEF, or
                                    YAC_CONFIG_OUTPUT_SYNC_LOC_ENDDEF)
     @param[in] include_definitions include user definitions (components, grids,
                                    and fields)
*/
     void yac_cset_config_output_file_instance( int yac_instance_id,
                                                const char * filename,
                                                int fileformat,
                                                int sync_location,
                                                int include_definitions);

/** Activates writing out of the coupling configuration

     @param[in] filename            name of the file to be written
     @param[in] fileformat          file format (YAC_CONFIG_OUTPUT_FORMAT_YAML
                                    or YAC_CONFIG_OUTPUT_FORMAT_JSON)
     @param[in] sync_location       synchronisation point after which the file
                                    is to be written
                                    (YAC_CONFIG_OUTPUT_SYNC_LOC_DEF_COMP,
                                    YAC_CONFIG_OUTPUT_SYNC_LOC_SYNC_DEF, or
                                    YAC_CONFIG_OUTPUT_SYNC_LOC_ENDDEF)
     @param[in] include_definitions include user definitions (components, grids,
                                    and fields)
*/
     void yac_cset_config_output_file( const char * filename,
                                       int fileformat,
                                       int sync_location,
                                       int include_definitions);

/* -------------------------------------------------------------------------------- */

/** Activates writing out of the grid data

     @param[in] yac_instance_id YAC instance_id
     @param[in] gridname        name of the grid to be written
     @param[in] filename        name of the file to be written
     @remark the writing is done in parallel (see \ref io_config_detail)
*/
     void yac_cset_grid_output_file_instance( int yac_instance_id,
                                              const char * gridname,
                                              const char * filename);

/** Activates writing out of the grid data

     @param[in] gridname        name of the grid to be written
     @param[in] filename        name of the file to be written
     @remark the writing is done in parallel (see \ref io_config_detail)
*/
     void yac_cset_grid_output_file( const char * gridname,
                                     const char * filename);

/* -------------------------------------------------------------------------------- */

/** \brief Clean-up default YAC instance (see \ref phase_restart)

     Frees all memory associated with the default YAC instance.

     @remark Call \ref yac_cfinalize for freeing all memory allocated by YAC.
*/
     void yac_ccleanup ();

/** \brief Clean-up a YAC instance (see \ref phase_restart)

     Frees all memory associated with the provided YAC instance.

     @param[in] yac_instance_id id of the YAC instance to be cleaned up
     @remark Call \ref yac_cfinalize for freeing all memory allocated by YAC.
*/
     void yac_ccleanup_instance ( int yac_instance_id );

/* -------------------------------------------------------------------------------- */

/** \brief Finalises YAC

     First cleans up the default YAC instance, if it was initialised and not
     yet cleaned up by a call to \ref yac_ccleanup.
     If there are no remaining YAC instances, MPI and yaxt will be finalised
     (if initialised by YAC) and instance independed memory allocated by YAC
     will be freed.
*/
     void yac_cfinalize ();

/** \brief Finalises YAC

     First cleans up the provided YAC instance, if it was initialised and not
     yet cleaned up by a call to \ref yac_ccleanup_instance.
     If there are no remaining YAC instances, MPI and yaxt will be finalised
     (if initialised by YAC) and instance independed memory allocated by YAC
     will be freed.
     @param[in] yac_instance_id id of the YAC instance to be finalised
*/
     void yac_cfinalize_instance (int yac_instance_id);

/* -------------------------------------------------------------------------------- */

/** Definition of job start and end datetime for the default YAC instance

     @param[in] start_datetime calendar job start datetime
     @param[in] end_datetime   calendar job end datetime
 */

     void yac_cdef_datetime ( const char * start_datetime,
                              const char * end_datetime );

/** Definition of job start and end datetime

     @param[in] yac_instance_id id of the YAC instance
     @param[in] start_datetime  calendar job start datetime
     @param[in] end_datetime    calendar job end datetime
 */

     void yac_cdef_datetime_instance ( int yac_instance_id,
                                       const char * start_datetime,
                                       const char * end_datetime );

/** Defines the calendar of the default instance

     @param[in] calendar The calendar @see YAC_PROLEPTIC_GREGORIAN etc.
*/

     void yac_cdef_calendar ( int calendar );

/** Gets the calendar of the default instance
 */

     int yac_cget_calendar ( );

/* -------------------------------------------------------------------------------- */

/** Providing a group communicator for the default YAC instance

     @param[out] group_comm communicator containing all processes that
                            passed the same group name to the init-routine
                            or MPI_COMM_NULL if not group name was provided.
 */

     void yac_cget_groupcomm( MPI_Comm * group_comm );

/** Providing a group communicator

     @param[in] yac_instance_id id of the YAC instance
     @param[out] group_comm communicator containing all processes that
                            passed the same group name to the init-routine
                            or MPI_COMM_NULL if not group name was provided.
 */

     void yac_cget_groupcomm_instance( int yac_instance_id,
                                       MPI_Comm * group_comm );

/* -------------------------------------------------------------------------------- */

/** Non-blocking pre-definition of a component for the default YAC instance.

     @param[in]  comp_name name of the component
     @param[out] comp_id   component Id
     @remark this call does not finish the component definition phase. I.e. the
             component communicators are not set up after this call.
*/

     void yac_cpredef_comp( char const * comp_name,
                            int * comp_id );

/** Non-blocking pre-definition of a component.

     @param[in]  yac_instance_id id of the YAC instance
     @param[in]  comp_name name of the component
     @param[out] comp_id   component Id
     @remark this call does not finish the component definition phase. I.e. the
             component communicators are not set up after this call.
*/

     void yac_cpredef_comp_instance( int yac_instance_id,
                                     char const * comp_name,
                                     int * comp_id );

/* -------------------------------------------------------------------------------- */

/** Elementary definition of the component for the default YAC instance

     @param[in]  comp_name name of the component
     @param[out] comp_id   component Id
     @remark this call is collective for all processes that initialised the
             default YAC instance
     @remark components can only be defined once in the initialisation phase of
             a YAC instance. If you want to predefine a component, please use `yac_cpredef_comp`.
*/

     void yac_cdef_comp ( const char * comp_name,
                          int * comp_id );

/** Elementary definition of the component.

     @param[in]  yac_instance_id id of the YAC instance
     @param[in]  comp_name name of the component
     @param[out] comp_id   component Id
     @remark this call is collective for all processes that initialised the
             provided YAC instance
     @remark components can only be defined once in the initialisation phase of
             a YAC instance. If you want to predefine a component, please use
             `yac_cpredef_comp_instance`.
*/

     void yac_cdef_comp_instance ( int yac_instance_id,
                                   const char * comp_name,
                                   int * comp_id );

/* -------------------------------------------------------------------------------- */

/** Elementary definition of the components for the default YAC instance

     @param[in]  comp_names names of the components
     @param[in]  num_comps  number of components
     @param[out] comp_ids   component Id's
     @remark this call is collective for all processes that initialised the
             default YAC instance
     @remark components can only be defined once in the initialisation phase of
             a YAC instance
*/

     void yac_cdef_comps ( const char ** comp_names,
                           int num_comps,
                           int * comp_ids );

/** Elementary definition of the components.

     @param[in]  yac_instance_id id of the YAC instance
     @param[in]  comp_names names of the components
     @param[in]  num_comps  number of components
     @param[out] comp_ids   component Id's
     @remark this call is collective for all processes that initialised the
             provided YAC instance
     @remark components can only be defined once in the initialisation phase of
             a YAC instance
*/

     void yac_cdef_comps_instance ( int yac_instance_id,
                                    const char ** comp_names,
                                    int num_comps,
                                    int * comp_ids );

/* -------------------------------------------------------------------------------- */

/** Providing the component MPI communicator

     @param[in]  comp_id   local component ID (from yac_cdef_comp),
                           handle to the component struct
     @param[out] comp_comm component communicator

*/
     void yac_cget_comp_comm ( int comp_id,
                               MPI_Comm *comp_comm );

/* -------------------------------------------------------------------------------- */

/** Generates an MPI communicator that contains all processes of the provided
 *  components

     @param[in]  comp_names name of components
     @param[in]  num_comps  number of components
     @param[out] comps_comm generated communicator
     @remark the components have to defined in the default YAC instance
     @remark the local process has to be in at least one of the provided
             components
     @remark this call is collective for all processes in the provided
             components
*/
     void yac_cget_comps_comm ( const char ** comp_names,
                                int num_comps,
                                MPI_Comm * comps_comm);

/** Generates an MPI communicator that contains all processes of the provided
 *  components

     @param[in]  yac_instance_id    id of the YAC instance
     @param[in]  comp_names name of components
     @param[in]  num_comps  number of components
     @param[out] comps_comm generated communicator
     @remark the local process has to be in at least one of the provided
             components
     @remark this call is collective for all processes in the provided
             components
*/
     void yac_cget_comps_comm_instance ( int yac_instance_id,
                                         const char ** comp_names,
                                         int num_comps,
                                         MPI_Comm * comps_comm );

/* -------------------------------------------------------------------------------- */

/** Definition of a set of points for 2D regular grids.

     @param[in]  grid_id    handle to the grid struct
     @param[in]  nbr_points number of points
     @param[in]  location   location of points
     @param[in]  x_points   array of point longitudes, in radians
     @param[in]  y_points   array of point latitudes, in radians
     @param[out] point_id   returned handle to the points struct
     @remark
       - This routine supports points defined on cells or vertices. For edges
         use \ref yac_cdef_points_unstruct instead.
       - The array x_points and y_points are expected to be of size
         nbr_points[0] and nbr_points[1] respectively.
       - See \ref yac_cdef_grid_reg2d on how the coordinates are being
         interpreted.
*/

     void yac_cdef_points_reg2d ( int const grid_id,
                                  int const *nbr_points,
                                  int const location,
                                  double const *x_points,
                                  double const *y_points,
                                  int *point_id );

/* -------------------------------------------------------------------------------- */

/** Definition of a set of points for 2D curvilinear grids.

     @param[in]  grid_id    handle to the grid struct
     @param[in]  nbr_points number of points
     @param[in]  location   location of points
     @param[in]  x_points   array of point longitudes, in radians
     @param[in]  y_points   array of point latitudes, in radians
     @param[out] point_id   returned handle to the points struct
     @remark
       - This routine supports points defined on cells or vertices. For edges
         use \ref yac_cdef_points_unstruct instead.
       - The array x_points and y_points are expected to be of size
         nbr_points[0] * nbr_points[1].
       - See \ref yac_cdef_grid_curve2d on how the coordinates are being
         interpreted.
*/

     void yac_cdef_points_curve2d ( int const grid_id,
                                    int const *nbr_points,
                                    int const location,
                                    double const *x_points,
                                    double const *y_points,
                                    int *point_id );

/* -------------------------------------------------------------------------------- */

/** Definition of a set of points for 2D unstructured grids.

     @param[in]  grid_id    handle to the grid struct
     @param[in]  nbr_points number of points
     @param[in]  location   location of points
     @param[in]  x_points   array of point longitudes, in radians
     @param[in]  y_points   array of point latitudes, in radians
     @param[out] point_id   returned handle to the points struct
     @remark The array x_points and y_points are expected to be of size
             nbr_points.
     @see \ref yac_cdef_grid_unstruct
*/

     void yac_cdef_points_unstruct ( int const grid_id,
                                     int const nbr_points,
                                     int const location,
                                     double const *x_points,
                                     double const *y_points,
                                     int *point_id );

/* -------------------------------------------------------------------------------- */

/** Definition of a set of points for 2D regular rotated grid.

     @param[in]  grid_id    handle to the grid struct
     @param[in]  nbr_points number of points
     @param[in]  location   location of points
     @param[in]  x_points   array of point longitudes, in radians
     @param[in]  y_points   array of point latitudes, in radians
     @param[in]  x_north_pole longitude of the north pole, in radians
     @param[in]  y_north_pole latitude of the north pole, in radians
     @param[out] point_id   returned handle to the points struct
     @remark
       - This routine supports points defined on cells or vertices. For edges
         use \ref yac_cdef_points_unstruct instead.
       - The array x_points and y_points are expected to be of size
         nbr_points[0] and nbr_points[1] respectively.
       - See \ref yac_cdef_grid_reg2d_rot on how the coordinates are being
         interpreted.
*/

     void yac_cdef_points_reg2d_rot ( int const grid_id,
                                      int const *nbr_points,
                                      int const location,
                                      double const *x_points,
                                      double const *y_points,
                                      double x_north_pole,
                                      double y_north_pole,
                                      int *point_id );

/* -------------------------------------------------------------------------------- */

/** Definition of a 2d regular grid

     @param[in]  grid_name    name of the grid
     @param[in]  nbr_vertices 2d array containing the number of vertices in each dimension
     @param[in]  cyclic       2d array containing information about cyclic behaviour in each dimension
     @param[in]  x_vertices   array of vertex longitudes, in radians
     @param[in]  y_vertices   array of vertex latitudes, in radians
     @param[out] grid_id      id of generated grid
     @remark
       - This call generate a grid with the following basic grid data:
         \code{.c}
         size_t internal_nbr_cells_2d[2];
         size_t internal_nbr_vertices_2d[2];

         if (cyclic[0]) {
           internal_nbr_cells_2d[0] = nbr_vertices[0];
           internal_nbr_cells_2d[1] = nbr_vertices[1] - 1;
           internal_nbr_vertices_2d[0] = nbr_vertices[0];
           internal_nbr_vertices_2d[1] = nbr_vertices[1];
         } else {
           internal_nbr_cells_2d[0] = nbr_vertices[0] - 1;
           internal_nbr_cells_2d[1] = nbr_vertices[1] - 1;
           internal_nbr_vertices_2d[0] = nbr_vertices[0];
           internal_nbr_vertices_2d[1] = nbr_vertices[1];
         }

         size_t internal_nbr_cells =
           internal_nbr_cells_2d[0] * internal_nbr_cells_2d[1];
         size_t internal_nbr_vertices =
           internal_nbr_vertices_2d[0] * internal_nbr_vertices_2d[1];
         size_t internal_nbr_edges =
           internal_nbr_cells_2d[1] * internal_nbr_vertices_2d[0] +
           internal_nbr_cells_2d[0] * internal_nbr_vertices_2d[1];

         double * internal_x_vertices =
           malloc(internal_nbr_vertices * sizeof(*internal_x_vertices));
         double * internal_y_vertices =
           malloc(internal_nbr_vertices * sizeof(*internal_y_vertices));

         for (size_t i = 0, k = 0; i < nbr_vertices[1]; ++i) {
           for (size_t j = 0; j < nbr_vertices[0]; ++j, ++k) {
             internal_x_vertices[k] = x_vertices[j];
             internal_y_vertices[k] = y_vertices[i];
           }
         }
         \endcode
       - The edges of the grid either follow circles of longitude or latitude.
       - See \ref phase_def_grid_edge_ordering
*/
     void yac_cdef_grid_reg2d ( const char * grid_name,
                                int nbr_vertices[2],
                                int cyclic[2],
                                double *x_vertices,
                                double *y_vertices,
                                int *grid_id);

/* -------------------------------------------------------------------------------- */

/** Definition of a 2d curvilinear grid

     @param[in]  grid_name    name of the grid
     @param[in]  nbr_vertices 2d array containing the number of vertices in each dimension
     @param[in]  cyclic       2d array containing information about cyclic behaviour in each dimension
     @param[in]  x_vertices   array of vertex longitudes, in radians
     @param[in]  y_vertices   array of vertex latitudes, in radians
     @param[out] grid_id      id of generated grid
     @remark
       - The array x_vertices and y_vertices are expected to be of size
         nbr_vertices[0] * nbr_vertices[1].
       - This call generate a grid with the following basic grid data:
         \code{.c}
         size_t internal_nbr_cells_2d[2];
         size_t internal_nbr_vertices_2d[2];

         if (cyclic[0]) {
           internal_nbr_cells_2d[0] = nbr_vertices[0];
           internal_nbr_cells_2d[1] = nbr_vertices[1] - 1;
           internal_nbr_vertices_2d[0] = nbr_vertices[0] + 1;
           internal_nbr_vertices_2d[1] = nbr_vertices[1];
         } else {
           internal_nbr_cells_2d[0] = nbr_vertices[0] - 1;
           internal_nbr_cells_2d[1] = nbr_vertices[1] - 1;
           internal_nbr_vertices_2d[0] = nbr_vertices[0];
           internal_nbr_vertices_2d[1] = nbr_vertices[1];
         }

         size_t internal_nbr_cells =
           internal_nbr_cells_2d[0] * internal_nbr_cells_2d[1];
         size_t internal_nbr_vertices =
           internal_nbr_vertices_2d[0] * internal_nbr_vertices_2d[1];
         size_t internal_nbr_edges =
           (internal_nbr_cells_2d[0] + 1) * internal_nbr_cells_2d[1] +
           internal_nbr_cells_2d[0] * (internal_nbr_cells_2d[1] + 1);

         double * internal_x_vertices =
           malloc(internal_nbr_vertices * sizeof(*internal_x_vertices));
         double * internal_y_vertices =
           malloc(internal_nbr_vertices * sizeof(*internal_y_vertices));

         for (size_t i = 0; i < internal_nbr_vertices; ++i) {
           internal_x_vertices[i] = x_vertices[i];
           internal_y_vertices[i] = y_vertices[i];
         }
         \endcode
       - The edges of the grid follow great circles.
       - See \ref phase_def_grid_edge_ordering
*/
     void yac_cdef_grid_curve2d ( const char * grid_name,
                                  int nbr_vertices[2],
                                  int cyclic[2],
                                  double *x_vertices,
                                  double *y_vertices,
                                  int *grid_id);

/* -------------------------------------------------------------------------------- */

/** Definition of an unstructured grid

     @param[in]  grid_name             name of the grid
     @param[in]  nbr_vertices          number of vertices
     @param[in]  nbr_cells             number of cells
     @param[in]  num_vertices_per_cell array containing the number of vertices for each cell
     @param[in]  x_vertices            array of vertex longitudes, in radians
     @param[in]  y_vertices            array of vertex latitudes, in radians
     @param[in]  cell_to_vertex        connectivity of vertices belonging to cells\n
                                       (the vertex indices per cell have to be in clockwise or counterclockwise ordering)
     @param[out] grid_id               id of generated grid
     @remark
       - The edges of the grid follow great circles.
       - See \ref phase_def_grid_edge_ordering
*/
     void yac_cdef_grid_unstruct ( const char * grid_name,
                                   int nbr_vertices,
                                   int nbr_cells,
                                   int *num_vertices_per_cell,
                                   double *x_vertices,
                                   double *y_vertices,
                                   int *cell_to_vertex,
                                   int *grid_id);

/* -------------------------------------------------------------------------------- */

/** Definition of an unstructured grid with lon-lat edges

     @param[in]  grid_name             name of the grid
     @param[in]  nbr_vertices          number of vertices
     @param[in]  nbr_cells             number of cells
     @param[in]  num_vertices_per_cell array containing the number of vertices for each cell
     @param[in]  x_vertices            array of vertex longitudes, in radians
     @param[in]  y_vertices            array of vertex latitudes, in radians
     @param[in]  cell_to_vertex        connectivity of vertices belonging to cells\n
                                       (the vertex indices per cell have to be in clockwise or counterclockwise ordering)
     @param[out] grid_id               id of generated grid
     @remark See \ref phase_def_grid_edge_ordering
     @remark YAC will check all edges of the grid and determine whether they are on
             circles of longitudes (same x coordinate) or latitudes (same y coordinate).
             An edge that does not fulfill this condition will cause an error.
*/
     void yac_cdef_grid_unstruct_ll ( const char * grid_name,
                                      int nbr_vertices,
                                      int nbr_cells,
                                      int *num_vertices_per_cell,
                                      double *x_vertices,
                                      double *y_vertices,
                                      int *cell_to_vertex,
                                      int *grid_id);

/* -------------------------------------------------------------------------------- */

/** Definition of an unstructured grid with explicit edge definition

     @param[in]  grid_name          name of the grid
     @param[in]  nbr_vertices       number of vertices
     @param[in]  nbr_cells          number of cells
     @param[in]  nbr_edges          number of edges
     @param[in]  num_edges_per_cell array containing the number of edges for each cell
     @param[in]  x_vertices         array of vertex longitudes, in radians
     @param[in]  y_vertices         array of vertex latitudes, in radians
     @param[in]  cell_to_edge       connectivity of edges belonging to cells\n
                                    (the edge indices per cell have to be in clockwise or counterclockwise ordering)
     @param[in]  edge_to_vertex     connectivity of vertices belonging to edges\n
     @param[out] grid_id            id of generated grid
     @remark
       - The edges of the grid follow great circles.
*/
     void yac_cdef_grid_unstruct_edge ( const char * grid_name,
                                        int nbr_vertices,
                                        int nbr_cells,
                                        int nbr_edges,
                                        int *num_edges_per_cell,
                                        double *x_vertices,
                                        double *y_vertices,
                                        int *cell_to_edge,
                                        int *edge_to_vertex,
                                        int *grid_id);

/* -------------------------------------------------------------------------------- */

/** Definition of an unstructured grid with lon-lat edges with explicit edge definition

     @param[in]  grid_name          name of the grid
     @param[in]  nbr_vertices       number of vertices
     @param[in]  nbr_cells          number of cells
     @param[in]  nbr_edges          number of edges
     @param[in]  num_edges_per_cell array containing the number of edges for each cell
     @param[in]  x_vertices         array of vertex longitudes, in radians
     @param[in]  y_vertices         array of vertex latitudes, in radians
     @param[in]  cell_to_edge       connectivity of edges belonging to cells\n
                                    (the edge indices per cell have to be in clockwise or counterclockwise ordering)
     @param[in]  edge_to_vertex     connectivity of vertices belonging to edges\n
     @param[out] grid_id            id of generated grid
     @remark YAC will check all edges of the grid and determine whether they are on
             circles of longitudes (same x coordinate) or latitudes (same y coordinate).
             An edge that does not fulfill this condition will cause an error.
*/
     void yac_cdef_grid_unstruct_edge_ll ( const char * grid_name,
                                           int nbr_vertices,
                                           int nbr_cells,
                                           int nbr_edges,
                                           int *num_edges_per_cell,
                                           double *x_vertices,
                                           double *y_vertices,
                                           int *cell_to_edge,
                                           int *edge_to_vertex,
                                           int *grid_id);

/* -------------------------------------------------------------------------------- */

/** Definition of grid consisting of a cloud of points

     @param[in]  grid_name  name of the grid
     @param[in]  nbr_points number of points
     @param[in]  x_points   array of point longitudes, in radians
     @param[in]  y_points   array of point latitudes, in radians
     @param[out] grid_id    id of generated grid
*/
     void yac_cdef_grid_cloud ( const char * grid_name,
                                int nbr_points,
                                double *x_points,
                                double *y_points,
                                int *grid_id);

/** Definition of a rotated 2d regular grid

     @param[in]  grid_name    name of the grid
     @param[in]  nbr_vertices 2d array containing the number of vertices in each dimension
     @param[in]  cyclic       2d array containing information about cyclic behaviour in each dimension
     @param[in]  x_vertices   array of vertex longitudes, in radians
     @param[in]  y_vertices   array of vertex latitudes, in radians
     @param[in]  x_north_pole longitude of the rotated north pole, in radians
     @param[in]  y_north_pole longitude of the rotated north pole, in radians
     @param[out] grid_id      id of generated grid
     @remark
       - This call generate a grid with the following basic grid data:
         \code{.c}
         size_t internal_nbr_cells_2d[2];
         size_t internal_nbr_vertices_2d[2];

         if (cyclic[0]) {
           internal_nbr_cells_2d[0] = nbr_vertices[0];
           internal_nbr_cells_2d[1] = nbr_vertices[1] - 1;
           internal_nbr_vertices_2d[0] = nbr_vertices[0];
           internal_nbr_vertices_2d[1] = nbr_vertices[1];
         } else {
           internal_nbr_cells_2d[0] = nbr_vertices[0] - 1;
           internal_nbr_cells_2d[1] = nbr_vertices[1] - 1;
           internal_nbr_vertices_2d[0] = nbr_vertices[0];
           internal_nbr_vertices_2d[1] = nbr_vertices[1];
         }

         size_t internal_nbr_cells =
           internal_nbr_cells_2d[0] * internal_nbr_cells_2d[1];
         size_t internal_nbr_vertices =
           internal_nbr_vertices_2d[0] * internal_nbr_vertices_2d[1];
         size_t internal_nbr_edges =
           internal_nbr_cells_2d[1] * internal_nbr_vertices_2d[0] +
           internal_nbr_cells_2d[0] * internal_nbr_vertices_2d[1];

         double * internal_x_vertices =
           malloc(internal_nbr_vertices * sizeof(*internal_x_vertices));
         double * internal_y_vertices =
           malloc(internal_nbr_vertices * sizeof(*internal_y_vertices));

         for (size_t i = 0, k = 0; i < nbr_vertices[1]; ++i) {
           for (size_t j = 0; j < nbr_vertices[0]; ++j, ++k) {
             double x_vertex_rot, y_vertex_rot;
             rotate(
               x_vertices[j], y_vertices[i], &x_vertex_rot, &y_vertex_rot);
             internal_x_vertices[k] = x_vertex_rot;
             internal_y_vertices[k] = y_vertex_rot;
           }
         }
         \endcode
       - The edges of the grid follow great circles.
       - See \ref phase_def_grid_edge_ordering
*/
     void yac_cdef_grid_reg2d_rot ( const char * grid_name,
                                    int nbr_vertices[2],
                                    int cyclic[2],
                                    double *x_vertices,
                                    double *y_vertices,
                                    double x_north_pole,
                                    double y_north_pole,
                                    int *grid_id);

/* -------------------------------------------------------------------------------- */

/** Set global ids for a grid

     @param[in] global_index array of global indices
     @param[in] location     cell/vertex/edge
     @param[in] grid_id      grid id
     @remark global indices are to be provided in the local order of their
             respective cells/vertices/edges
*/

     void yac_cset_global_index ( int const * global_index,
                                  int location,
                                  int grid_id);

/* -------------------------------------------------------------------------------- */

/** Set core mask for a grid

     @param[in] is_core   0 for cells/vertices/edges that are halos,
                          1 for cells/vertices/edges that are core
     @param[in] location  cell/vertex/edge
     @param[in] grid_id   grid id
     @remark cells/vertices/edges who are halos are not used or set
             in a put or get operation
     @remark core mask values are to be provided in the local order of their
             respective cells/vertices/edges
*/

     void yac_cset_core_mask ( int const * is_core,
                               int location,
                               int grid_id);

/* -------------------------------------------------------------------------------- */

/** Set the default mask for points

     @param[in] is_valid   0 for points that are masked out, 1 for valid points
     @param[in] points_id  id of points/cells
*/

     void yac_cset_mask ( int const * is_valid,
                          int points_id );

/* -------------------------------------------------------------------------------- */

/** define a mask for a grid

     @param[in]  grid_id    grid ID
     @param[in]  nbr_points number of points
     @param[in]  location   cell/vertex/edge
     @param[in]  is_valid   0 for points that are masked out, 1 for valid points
     @param[out] mask_id    returned handle to the mask struct
*/

     void yac_cdef_mask ( int const grid_id,
                          int const nbr_points,
                          int const location,
                          int const * is_valid,
                          int *mask_id );

/* -------------------------------------------------------------------------------- */

/** define a named mask for a grid

     @param[in]  grid_id    grid ID
     @param[in]  nbr_points number of points
     @param[in]  location   cell/vertex/edge
     @param[in]  is_valid   0 for points that are masked out, 1 for valid points
     @param[in]  name       name of the mask
     @param[out] mask_id    returned handle to the mask struct
*/

     void yac_cdef_mask_named ( int const grid_id,
                                int const nbr_points,
                                int const location,
                                int const * is_valid,
                                char const * name,
                                int *mask_id );

/* -------------------------------------------------------------------------------- */

/** Definition of the coupling field (using default mask, if defined)

     @param[in]  field_name      character string providing the name of the coupling field
     @param[in]  component_id    component ID
     @param[in]  point_ids       point IDs
     @param[in]  num_pointsets   number of pointsets per grid
     @param[in]  collection_size collection size
     @param[in]  timestep        timestep
     @param[in]  time_unit       time unit
     @param[out] field_id        returned field_id which has to be used to identify coupling fields
                                 in yac_cput and yac_cget.
     @remark A calendar has to be defined before this routine is called.
             This can be done by a call to \ref yac_cdef_calendar, by reading
             of configuration file that defines a calendar, or by
             synchronizing with other processes that already have defined
             a calendar using \ref yac_csync_def.
*/

       void yac_cdef_field ( char const * field_name,
                             int const component_id,
                             int const * point_ids,
                             int const num_pointsets,
                             int collection_size,
                             const char* timestep,
                             int time_unit,
                             int * field_id );

/* -------------------------------------------------------------------------------- */

/** Definition of the coupling field

     @param[in]  field_name      character string providing the name of the coupling field
     @param[in]  component_id    component ID
     @param[in]  point_ids       point IDs
     @param[in]  mask_ids        mask IDs
     @param[in]  num_pointsets   number of pointsets per grid
     @param[in]  collection_size collection size
     @param[in]  timestep        timestep
     @param[in]  time_unit
     @param[out] field_id        returned field_id which has to be used to identify coupling fields
                                 in yac_cput and yac_cget.
     @remark A calendar has to be defined before this routine is called.
             This can be done by a call to \ref yac_cdef_calendar, by reading
             of configuration file that defines a calendar, or by
             synchronizing with other processes that already have defined
             a calendar using \ref yac_csync_def.
*/

     void yac_cdef_field_mask ( char const * field_name,
                                int const component_id,
                                int const * point_ids,
                                int const * mask_ids,
                                int const num_pointsets,
                                int collection_size,
                                const char* timestep,
                                int time_unit,
                                int * field_id );

/* -------------------------------------------------------------------------- */

/** Enables fractional masking for a coupling field and sets the fallback value

     @param[in]  comp_name
     @param[in]  grid_name
     @param[in]  field_name
     @param[in]  frac_mask_fallback_value fractional mask fallback value
     @remark If this field is used as a source in couple, the user will have
             to provide a fractional mask for all source points along with the
             source field values in all put/exchange operations of this
             field. (see \ref frac_mask_desc)
*/

     void yac_cenable_field_frac_mask( const char* comp_name,
                                       const char* grid_name,
                                       const char* field_name,
                                       double frac_mask_fallback_value);

/** Enables fractional masking for a coupling field and sets the fallback value

     @param[in] yac_instance_id
     @param[in]  comp_name
     @param[in]  grid_name
     @param[in]  field_name
     @param[in]  frac_mask_fallback_value fractional mask fallback value
     @remark If this field is used as a source in couple, the user will have
             to provide a fractional mask for all source points along with the
             source field values in all put/exchange operations of this
             field. (see \ref frac_mask_desc)
*/

     void yac_cenable_field_frac_mask_instance(
       int yac_instance_id,
       const char* comp_name,
       const char* grid_name,
       const char* field_name,
       double frac_mask_fallback_value);

/* -------------------------------------------------------------------------- */

/** Define metadata for a component
     @param[in] comp_name
     @param[in] metadata
*/
     void yac_cdef_component_metadata( const char* comp_name,
                                       const char* metadata);

/** Define metadata for a component
     @param[in] yac_instance_id
     @param[in] comp_name
     @param[in] metadata
*/
     void yac_cdef_component_metadata_instance( int yac_instance_id,
                                                const char* comp_name,
                                                const char* metadata);

/** Define metadata for a grid for the default YAC instance
     @param[in] grid_name
     @param[in] metadata
*/
     void yac_cdef_grid_metadata( const char* grid_name,
                                  const char* metadata);

/** Define metadata for a grid
     @param[in] yac_instance_id
     @param[in] grid_name
     @param[in] metadata
*/
     void yac_cdef_grid_metadata_instance( int yac_instance_id,
                                           const char* grid_name,
                                           const char* metadata);

/** Define metadata for a field
     @param[in] comp_name
     @param[in] grid_name
     @param[in] field_name
     @param[in] metadata
*/
     void yac_cdef_field_metadata( const char* comp_name,
                                   const char* grid_name,
                                   const char* field_name,
                                   const char* metadata);

/** Define metadata for a field
     @param[in] yac_instance_id
     @param[in] comp_name
     @param[in] grid_name
     @param[in] field_name
     @param[in] metadata
*/
     void yac_cdef_field_metadata_instance( int yac_instance_id,
                                            const char* comp_name,
                                            const char* grid_name,
                                            const char* field_name,
                                            const char* metadata);

/* -------------------------------------------------------------------------- */

/** Returns the action a put/get would return for the current timestep

     @param[in]  field_id field id returned by \ref yac_cdef_field
     @param[out] action   action for the current timestep
                          (\ref YAC_ACTION_NONE,
                           \ref YAC_ACTION_REDUCTION,
                           \ref YAC_ACTION_COUPLING,
                           \ref YAC_ACTION_GET_FOR_RESTART,
                           \ref YAC_ACTION_PUT_FOR_RESTART, or
                           \ref YAC_ACTION_OUT_OF_BOUND)
     @remark If the exchange type of the field is \ref YAC_EXCHANGE_TYPE_NONE,
             this routine will return \ref YAC_ACTION_NONE.
     @remark If this routine is called before \ref yac_cenddef, the behaviour
             is undefined.
*/

     void yac_cget_action( int field_id,
                           int * action);

/** Returns the current datetime of the field

     @param[in] field_id
     @returns current datetime of the field
*/
     const char * yac_cget_field_datetime(int field_id);

/* -------------------------------------------------------------------------- */

/** If the action for the current timestep is \ref YAC_ACTION_NONE,
    this routine will advance the internal clock to next timestep.

     @param[in]  field_id field id returned by \ref yac_cdef_field
     @remark If this routine is called before \ref yac_cenddef, the behaviour
             is undefined.
     @remark This routine can be called instead of \ref yac_cput or \ref yac_cget,
             if these routines would not have taken any action for the current
             timestep, in order to advance the internal clock.
*/

     void yac_cupdate(int field_id);

/* -------------------------------------------------------------------------- */

/** Get an extended coupling configuration (all parameters are set
    to the default values)
     @param[out] ext_couple_config_id extended coupling configuration
 */
     void yac_cget_ext_couple_config(int * ext_couple_config_id);

/** Frees an extended coupling configuration
     @param[in] ext_couple_config_id
 */
     void yac_cfree_ext_couple_config(int ext_couple_config_id);

/** Sets the weight file name
     @param[in] ext_couple_config_id extended coupling configuration
     @param[in] weight_file name of a weight file
     @remark if (weight_file == NULL) no weight file will be written
     @remark YAC will write out the weight file in parallel.
             This parallel output can be configured as described here:
             \ref io_config_detail
 */
     void yac_cset_ext_couple_config_weight_file( int ext_couple_config_id,
                                                  char const * weight_file);

/** Gets the weight file name
     @param[in]  ext_couple_config_id extended coupling configuration
     @param[out] weight_file          name of the weight file, if it is set;
                                      NULL otherwise
 */
     void yac_cget_ext_couple_config_weight_file( int ext_couple_config_id,
                                                  char const ** weight_file);

/** Specify handling of already existing weight files
     @param[in] ext_couple_config_id    extended coupling configuration
     @param[in] weight_file_on_existing handling of existing weight files
     @remark this only applies, if a weight file name was specified
     @remark Avaiable options:
             - \ref YAC_WGT_ON_EXISTING_ERROR
             - \ref YAC_WGT_ON_EXISTING_KEEP
             - \ref YAC_WGT_ON_EXISTING_OVERWRITE
     @see yac_cset_ext_couple_config_weight_file
 */
     void yac_cset_ext_couple_config_weight_file_on_existing(
          int ext_couple_config_id, int weight_file_on_existing);

/** Gets how already existing weight files are handled
     @param[in]  ext_couple_config_id   extended coupling configuration
     @param[in] weight_file_on_existing handling of existing weight files
 */
     void yac_cget_ext_couple_config_weight_file_on_existing(
          int ext_couple_config_id, int * weight_file_on_existing);

/** Sets the mapping side
     @param[in] ext_couple_config_id extended coupling configuration
     @param[in] mapping_side mapping side
     @remark Valid values are: source = 1 and target = 0
 */
     void yac_cset_ext_couple_config_mapping_side( int ext_couple_config_id,
                                                   int mapping_side);

/** Gets the mapping side
     @param[in]  ext_couple_config_id extended coupling configuration
     @param[out] mapping_side         mapping side
 */
     void yac_cget_ext_couple_config_mapping_side( int ext_couple_config_id,
                                                   int * mapping_side);

/** Sets the scale factor
     @param[in] ext_couple_config_id extended coupling configuration
     @param[in] scale_factor         scale factor
 */
     void yac_cset_ext_couple_config_scale_factor( int ext_couple_config_id,
                                                   double scale_factor);

/** Gets the scale factor
     @param[in]  ext_couple_config_id extended coupling configuration
     @param[out] scale_factor         scale_factor
 */
     void yac_cget_ext_couple_config_scale_factor( int ext_couple_config_id,
                                                   double * scale_factor);

/** Sets the scale summand
     @param[in] ext_couple_config_id extended coupling configuration
     @param[in] scale_summand        scale summand
 */
     void yac_cset_ext_couple_config_scale_summand( int ext_couple_config_id,
                                                    double scale_summand);

/** Gets the scale summand
     @param[in]  ext_couple_config_id extended coupling configuration
     @param[out] scale_summand        scale_summand
 */
     void yac_cget_ext_couple_config_scale_summand( int ext_couple_config_id,
                                                    double * scale_summand);

/** Sets source mask names
     @param[in] ext_couple_config_id extended coupling configuration
     @param[in] num_src_mask_names   number of source mask names
     @param[in] src_mask_names       source mask names
 */
     void yac_cset_ext_couple_config_src_mask_names(
       int ext_couple_config_id,
       size_t num_src_mask_names,
       char const * const * src_mask_names);

/** Gets source mask names
     @param[in]  ext_couple_config_id extended coupling configuration
     @param[out] num_src_mask_names   number of source mask names
     @param[out] src_mask_names       source mask names
 */
     void yac_cget_ext_couple_config_src_mask_names(
       int ext_couple_config_id,
       size_t * num_src_mask_names,
       char const * const ** src_mask_names);

/** Sets target mask name
     @param[in] ext_couple_config_id extended coupling configuration
     @param[in] tgt_mask_name        target mask name
 */
     void yac_cset_ext_couple_config_tgt_mask_name( int ext_couple_config_id,
                                                    char const * tgt_mask_name);

/** Gets target mask name
     @param[in]  ext_couple_config_id extended coupling configuration
     @param[out] tgt_mask_name        target mask name
 */
     void yac_cget_ext_couple_config_tgt_mask_name( int ext_couple_config_id,
                                                    char const ** tgt_mask_name);

/** Sets yaxt exchanger name
     @param[in] ext_couple_config_id extended coupling configuration
     @param[in] yaxt_exchanger_name  yaxt exchange name
 */
     void yac_cset_ext_couple_config_yaxt_exchanger_name(
          int ext_couple_config_id, char const * yaxt_exchanger_name);

/** Gets yaxt exchanger name
     @param[in]  ext_couple_config_id extended coupling configuration
     @param[out] yaxt_exchanger_name  yaxt exchanger name
 */
     void yac_cget_ext_couple_config_yaxt_exchanger_name(
          int ext_couple_config_id, char const ** yaxt_exchanger_name);


/** Sets use of raw exchange
     @param[in] ext_couple_config_id extended coupling configuration
     @param[in] use_raw_exchange     use raw exchange
 */
     void yac_cset_ext_couple_config_use_raw_exchange(
          int ext_couple_config_id, int use_raw_exchange);

/** Gets use of raw exchange
     @param[in]  ext_couple_config_id extended coupling configuration
     @param[out] use_raw_exchange     use raw exchange
 */
     void yac_cget_ext_couple_config_use_raw_exchange(
          int ext_couple_config_id, int * use_raw_exchange);

/* -------------------------------------------------------------------------- */

/** Define couple using an extended coupling configuration
     @param[in] yac_instance_id id of the YAC instance for that the couple
                should be defined
     @param[in] src_comp_name component name of the source component
     @param[in] src_grid_name grid name of the source grid
     @param[in] src_field_name field name of the source field
     @param[in] tgt_comp_name component name of the target component
     @param[in] tgt_grid_name grid name of the target grid
     @param[in] tgt_field_name field name of the target field
     @param[in] coupling_timestep time step for the coupling
     @param[in] time_unit unit of coupling_timestep argument
     @param[in] time_reduction type for reducing multiple timesteps
                @see YAC_REDUCTION_TIME_NONE etc.
     @param[in] interp_stack_config_id id of the interpolation stack config to
                be used
     @param[in] src_lag lag for this couple on the source component
     @param[in] tgt_lag lag for this couple on the target component
     @param[in] ext_couple_config_id extended coupling configuration
 */

     void yac_cdef_couple_custom_instance( int yac_instance_id,
                                           char const * src_comp_name,
                                           char const * src_grid_name,
                                           char const * src_field_name,
                                           char const * tgt_comp_name,
                                           char const * tgt_grid_name,
                                           char const * tgt_field_name,
                                           char const * coupling_timestep,
                                           int time_unit,
                                           int time_reduction,
                                           int interp_stack_config_id,
                                           int src_lag,
                                           int tgt_lag,
                                           int ext_couple_config_id);

/** Define couple in the default YAC instance using an extended
    coupling configuration
     @param[in] src_comp_name component name of the source component
     @param[in] src_grid_name grid name of the source grid
     @param[in] src_field_name field name of the source field
     @param[in] tgt_comp_name component name of the target component
     @param[in] tgt_grid_name grid name of the target grid
     @param[in] tgt_field_name field name of the target field
     @param[in] coupling_timestep time step for the coupling given in the defined
                timestep unit
     @param[in] time_unit unit of coupling_timestep argument
     @param[in] time_reduction type for reducing multiple timesteps
                @see YAC_REDUCTION_TIME_NONE etc.
     @param[in] interp_stack_config_id id of the interpolation stack config to
                be used
     @param[in] src_lag lag for this couple on the source component
     @param[in] tgt_lag lag for this couple on the target component
     @param[in] ext_couple_config_id extended coupling configuration
 */

     void yac_cdef_couple_custom( char const * src_comp_name,
                                  char const * src_grid_name,
                                  char const * src_field_name,
                                  char const * tgt_comp_name,
                                  char const * tgt_grid_name,
                                  char const * tgt_field_name,
                                  char const * coupling_timestep,
                                  int time_unit,
                                  int time_reduction,
                                  int interp_stack_config_id,
                                  int src_lag,
                                  int tgt_lag,
                                  int ext_couple_config_id);

/* -------------------------------------------------------------------------- */

/** Define couple
     @param[in] yac_instance_id id of the YAC instance for that the couple
                should be defined
     @param[in] src_comp_name component name of the source component
     @param[in] src_grid_name grid name of the source grid
     @param[in] src_field_name field name of the source field
     @param[in] tgt_comp_name component name of the target component
     @param[in] tgt_grid_name grid name of the target grid
     @param[in] tgt_field_name field name of the target field
     @param[in] coupling_timestep time step for the coupling
     @param[in] time_unit unit of coupling_timestep argument
     @param[in] time_reduction type for reducing multiple timesteps
                @see YAC_REDUCTION_TIME_NONE etc.
     @param[in] interp_stack_config_id id of the interpolation stack config to
                be used
     @param[in] src_lag lag for this couple on the source component
     @param[in] tgt_lag lag for this couple on the target component
 */

     void yac_cdef_couple_instance( int yac_instance_id,
                                    char const * src_comp_name,
                                    char const * src_grid_name,
                                    char const * src_field_name,
                                    char const * tgt_comp_name,
                                    char const * tgt_grid_name,
                                    char const * tgt_field_name,
                                    char const * coupling_timestep,
                                    int time_unit,
                                    int time_reduction,
                                    int interp_stack_config_id,
                                    int src_lag,
                                    int tgt_lag);

/** Define couple in the default YAC instance
     @param[in] src_comp_name component name of the source component
     @param[in] src_grid_name grid name of the source grid
     @param[in] src_field_name field name of the source field
     @param[in] tgt_comp_name component name of the target component
     @param[in] tgt_grid_name grid name of the target grid
     @param[in] tgt_field_name field name of the target field
     @param[in] coupling_timestep time step for the coupling given in the defined
                timestep unit
     @param[in] time_unit unit of coupling_timestep argument
     @param[in] time_reduction type for reducing multiple timesteps
                @see YAC_REDUCTION_TIME_NONE etc.
     @param[in] interp_stack_config_id id of the interpolation stack config to
                be used
     @param[in] src_lag lag for this couple on the source component
     @param[in] tgt_lag lag for this couple on the target component
 */

     void yac_cdef_couple( char const * src_comp_name,
                           char const * src_grid_name,
                           char const * src_field_name,
                           char const * tgt_comp_name,
                           char const * tgt_grid_name,
                           char const * tgt_field_name,
                           char const * coupling_timestep,
                           int time_unit,
                           int time_reduction,
                           int interp_stack_config_id,
                           int src_lag,
                           int tgt_lag);

/* -------------------------------------------------------------------------------- */

/** Checks whether the field dimensions match with field definition and aborts if
 *  there is a mismatch

     @param[in] field_id
     @param[in] collection_size collection size
     @param[in] num_pointsets   number of pointsets
     @param[in] pointset_sizes  size of each pointset
     @remark if num_pointsets is -1, its value will not be checked
     @remark if pointset_sizes is NULL, its values will not be checked
*/
     void yac_ccheck_field_dimensions ( int field_id,
                                        int collection_size,
                                        int num_pointsets,
                                        int const * pointset_sizes );

/* -------------------------------------------------------------------------------- */

/** Checks whether the source field buffer size matches the one expected by
 *  YAC for this field and abort if there is a mismatch

     @param[in] field_id
     @param[in] collection_size       collection size
     @param[in] src_field_buffer_size source buffer size
                                      (SUM(src_field_buffer_sizes(:)))
*/
     void yac_ccheck_src_field_buffer_size ( int field_id,
                                             int collection_size,
                                             int src_field_buffer_size );

/** Checks whether the source field buffer sizes matches the one expected by
 *  YAC for this field and abort if there is a mismatch

     @param[in] field_id
     @param[in] num_src_fields         number of source fields
     @param[in] collection_size        collection size
     @param[in] src_field_buffer_sizes source field buffer sizes
*/
     void yac_ccheck_src_field_buffer_sizes ( int field_id,
                                              int num_src_fields,
                                              int collection_size,
                                              int * src_field_buffer_sizes );

/* -------------------------------------------------------------------------------- */

/** Get interpolation information for raw data exchange

     @param[in]  field_id
     @param[out] frac_mask_fallback_value fractional masking fallback value
     @param[out] scaling_factor           scaling factor
     @param[out] scaling_summand          scaling summand

     @param[out] num_fixed_values         number of fixed values
     @param[out] fixed_values             array of fixed values
     @param[out] num_tgt_per_fixed_value  number of target points per fixed
                                          value
     @param[out] tgt_idx_fixed            local target point indices that
                                          receive a fixed value

     @param[out] num_wgt_tgt              number of target points receiving a
                                          weighted sum of source points
     @param[out] wgt_tgt_idx              local target point indices that
                                          receive weighted sum of source points
     @param[out] num_src_per_tgt          number of source points per target
                                          point
     @param[out] weights                  weights
     @param[out] src_field_idx            source field indices
     @param[out] src_idx                  source field buffer indices
     @param[out] num_src_fields           number of source fields
     @param[out] src_field_buffer_sizes   buffer sizes of source field buffers
     @remark If a coupling is configured for raw data exchange, target
             processes can extract once all information required to execute the
             interpolation by hand.
     @remark The user is responsible for freeing the memory associated with
             returned arrays.
*/
     void yac_cget_raw_interp_weights_data ( int const field_id,
                                             double * frac_mask_fallback_value,
                                             double * scaling_factor,
                                             double * scaling_summand,
                                             size_t * num_fixed_values,
                                             double ** fixed_values,
                                             size_t ** num_tgt_per_fixed_value,
                                             size_t ** tgt_idx_fixed,
                                             size_t * num_wgt_tgt,
                                             size_t ** wgt_tgt_idx,
                                             size_t ** num_src_per_tgt,
                                             double ** weights,
                                             size_t ** src_field_idx,
                                             size_t ** src_idx,
                                             size_t * num_src_fields,
                                             size_t ** src_field_buffer_sizes );

/** Get interpolation information for raw data exchange
 *  (weight matrix is stored as a sparse matrix using CSR format)

     @param[in]  field_id
     @param[out] frac_mask_fallback_value fractional masking fallback value
     @param[out] scaling_factor           scaling factor
     @param[out] scaling_summand          scaling summand

     @param[out] num_fixed_values         number of fixed values
     @param[out] fixed_values             array of fixed values
     @param[out] num_tgt_per_fixed_value  number of target points per fixed
                                          value
     @param[out] tgt_idx_fixed            local target point indices that
                                          receive a fixed value
     @param[out] src_indptr               offset to weights/src_idx/src_field_idx
                                          for each local target point
     @param[out] weights                  weights
     @param[out] src_field_idx            source field indices
     @param[out] src_idx                  source field buffer indices
     @param[out] num_src_fields           number of source fields
     @param[out] src_field_buffer_sizes   buffer sizes of source field buffers
     @remark If a coupling is configured for raw data exchange, target
             processes can extract once all information required to execute the
             interpolation by hand.
     @remark The user is responsible for freeing the memory associated with
             returned arrays.
     @remark Number of weights per target is computed by
             \c "src_indptr[tgt_idx+1] - src_indptr[tgt_idx]".
*/
     void yac_cget_raw_interp_weights_data_csr ( int const field_id,
                                                 double * frac_mask_fallback_value,
                                                 double * scaling_factor,
                                                 double * scaling_summand,
                                                 size_t * num_fixed_values,
                                                 double ** fixed_values,
                                                 size_t ** num_tgt_per_fixed_value,
                                                 size_t ** tgt_idx_fixed,
                                                 size_t ** src_indptr,
                                                 double ** weights,
                                                 size_t ** src_field_idx,
                                                 size_t ** src_idx,
                                                 size_t * num_src_fields,
                                                 size_t ** src_field_buffer_sizes );

/* -------------------------------------------------------------------------------- */

/** Receiving of the coupling fields

     @param[in]  field_id           -
     @param[in]  collection_size    -
     @param[in]  recv_field         - receive buffer (all data is stored in one
                                                      contiguous part of the memory)
                                      dimensions: recv_field[collection_size][nbr_points]
     @param[in]  info               - returned info argument indicating the action performed
     @param[out] ierror             - returned error
*/
     void yac_cget_ ( int const field_id,
                      int const collection_size,
                      double *recv_field,
                      int    *info,
                      int    *ierror );

/** Receiving of the coupling fields

     @param[in]  field_id           -
     @param[in]  collection_size    -
     @param[in]  recv_field         - receive buffer
                                      dimensions: recv_field[collection_size][nbr_points]
     @param[out] info               - returned info argument indicating the action performed
     @param[out] ierror             - returned error
*/
     void yac_cget ( int const field_id,
                     int const collection_size,
                     double **recv_field,
                     int    *info,
                     int    *ierror );

/** Receiving of the coupling fields

     @param[in]  field_id           -
     @param[in]  collection_size    -
     @param[in]  recv_field         - receive buffer (all data is stored in one
                                                      contiguous part of the memory)
                                      dimensions: recv_field[collection_size][nbr_points]
     @param[in]  info               - returned info argument indicating the action performed
     @param[out] ierror             - returned error
     @remark This routine returns immediately without waiting for the receive
             operation to be completed. The user should not access the memory
             associated with recv_field until this is the case. This can be
             achieved by calling \ref yac_cwait.
*/
     void yac_cget_async_ ( int const field_id,
                            int const collection_size,
                            double *recv_field,
                            int    *info,
                            int    *ierror );

/** Receiving of the coupling fields

     @param[in]  field_id           -
     @param[in]  collection_size    -
     @param[in]  recv_field         - receive buffer
                                      dimensions: recv_field[collection_size][nbr_points]
     @param[out] info               - returned info argument indicating the action performed
     @param[out] ierror             - returned error
     @remark This routine returns immediately without waiting for the receive
             operation to be completed. The user should not access the memory
             associated with recv_field until this is the case. This can be
             achieved by calling \ref yac_cwait.
*/
     void yac_cget_async ( int const field_id,
                           int const collection_size,
                           double **recv_field,
                           int    *info,
                           int    *ierror );

/** Receiving source points required for the interpolation of the local
    target field

     @param[in]  field_id
     @param[in]  collection_size
     @param[in]  src_field_buffer source field buffer
                                  (all data is stored in one contiguous part of
                                  the memory)
                                  (dimensions:
                                     src_field_buffer
                                       [collection_size]
                                       [num_src_fields]
                                       [src_field_buffer_sizes[src_field_idx]])
     @param[out] info             returned info argument indicating the
                                  action performed
     @param[out] ierror           returned error
*/
     void yac_cget_raw_ ( int const field_id,
                          int const collection_size,
                          double *src_field_buffer,
                          int *info,
                          int *ierror );

/** Receiving source points required for the interpolation of the local
    target field

     @param[in]  field_id
     @param[in]  collection_size
     @param[in]  src_field_buffer source field buffer
                                  (dimensions:
                                     src_field_buffer
                                       [collection_size]
                                       [num_src_fields]
                                       [src_field_buffer_sizes[src_field_idx]])
     @param[out] info             returned info argument indicating the
                                  action performed
     @param[out] ierror           returned error
*/
     void yac_cget_raw ( int const field_id,
                         int const collection_size,
                         double ***src_field_buffer,
                         int *info,
                         int *ierror );

/** Receiving source points required for the interpolation of the local
    target field

     @param[in]  field_id
     @param[in]  collection_size
     @param[in]  src_field_buffer source field buffer
                                  (all data is stored in one contiguous part of
                                  the memory)
                                  (dimensions:
                                     src_field_buffer
                                       [collection_size *
                                        SUM(src_field_buffer_sizes)])
     @param[out] info             returned info argument indicating the
                                  action performed
     @param[out] ierror           returned error
     @remark This routine returns immediately without waiting for the receive
             operation to be completed. The user should not access the memory
             associated with recv_field until this is the case. This can be
             achieved by calling \ref yac_cwait.
*/
     void yac_cget_raw_async_ ( int const field_id,
                                int const collection_size,
                                double *src_field_buffer,
                                int *info,
                                int *ierror );

/** Receiving source points required for the interpolation of the local
    target field

     @param[in]  field_id
     @param[in]  collection_size
     @param[in]  src_field_buffer source field buffer
                                  (dimensions:
                                     src_field_buffer
                                       [collection_size]
                                       [num_src_fields]
                                       [src_field_buffer_sizes[src_field_idx]])
     @param[out] info             returned info argument indicating the
                                  action performed
     @param[out] ierror           returned error
     @remark This routine returns immediately without waiting for the receive
             operation to be completed. The user should not access the memory
             associated with recv_field until this is the case. This can be
             achieved by calling \ref yac_cwait.
*/
     void yac_cget_raw_async ( int const field_id,
                               int const collection_size,
                               double ***src_field_buffer,
                               int *info,
                               int *ierror );

/** Receiving source points required for the interpolation of the local
    target field

     @param[in]  field_id
     @param[in]  collection_size
     @param[in]  src_field_buffer source field buffer
                                  (all data is stored in one contiguous part of
                                   the memory)
                                  (dimensions:
                                     src_field_buffer
                                       [collection_size]
                                       [num_src_fields]
                                       [src_field_buffer_sizes[src_field_idx]])
     @param[in]  src_frac_mask_buffer source fractional mask buffer
                                      (all data is stored in one contiguous
                                       part of the memory)
                                      (dimensions:
                                         src_frac_mask_buffer
                                           [collection_size]
                                           [num_src_fields]
                                           [src_field_buffer_sizes[src_field_idx]])
     @param[out] info             returned info argument indicating the
                                  action performed
     @param[out] ierror           returned error
     @remark source field values in the src_field_buffer have already been
             weighted, by the fractional mask
*/
     void yac_cget_raw_frac_ ( int const field_id,
                               int const collection_size,
                               double *src_field_buffer,
                               double *src_frac_mask_buffer,
                               int *info,
                               int *ierror );

/** Receiving source points required for the interpolation of the local
    target field

     @param[in]  field_id
     @param[in]  collection_size
     @param[in]  src_field_buffer source field buffer
                                  (dimensions:
                                     src_field_buffer
                                       [collection_size]
                                       [num_src_fields]
                                       [src_field_buffer_sizes[src_field_idx]])
     @param[in]  src_frac_mask_buffer source fractional mask buffer
                                      (dimensions:
                                         src_frac_mask_buffer
                                           [collection_size]
                                           [num_src_fields]
                                           [src_field_buffer_sizes[src_field_idx]])
     @param[out] info             returned info argument indicating the
                                  action performed
     @param[out] ierror           returned error
     @remark source field values in the src_field_buffer have already been
             weighted, by the fractional mask
*/
     void yac_cget_raw_frac ( int const field_id,
                              int const collection_size,
                              double ***src_field_buffer,
                              double ***src_frac_mask_buffer,
                              int *info,
                              int *ierror );

/** Receiving source points required for the interpolation of the local
    target field

     @param[in]  field_id
     @param[in]  collection_size
     @param[in]  src_field_buffer source field buffer
                                  (all data is stored in one contiguous part of
                                   the memory)
                                  (dimensions:
                                     src_field_buffer
                                       [collection_size]
                                       [num_src_fields]
                                       [src_field_buffer_sizes[src_field_idx]])
     @param[in]  src_frac_mask_buffer source fractional mask buffer
                                      (all data is stored in one contiguous
                                       part of the memory)
                                      (dimensions:
                                         src_frac_mask_buffer
                                           [collection_size]
                                           [num_src_fields]
                                           [src_field_buffer_sizes[src_field_idx]])
     @param[out] info             returned info argument indicating the
                                  action performed
     @param[out] ierror           returned error
     @remark source field values in the src_field_buffer have already been
             weighted, by the fractional mask
     @remark This routine returns immediately without waiting for the receive
             operation to be completed. The user should not access the memory
             associated with recv_field until this is the case. This can be
             achieved by calling \ref yac_cwait.
*/
     void yac_cget_raw_frac_async_ ( int const field_id,
                                     int const collection_size,
                                     double *src_field_buffer,
                                     double *src_frac_mask_buffer,
                                     int *info,
                                     int *ierror );

/** Receiving source points required for the interpolation of the local
    target field

     @param[in]  field_id
     @param[in]  collection_size
     @param[in]  src_field_buffer source field buffer
                                  (dimensions:
                                     src_field_buffer
                                       [collection_size]
                                       [num_src_fields]
                                       [src_field_buffer_sizes[src_field_idx]])
     @param[in]  src_frac_mask_buffer source fractional mask buffer
                                      (dimensions:
                                         src_frac_mask_buffer
                                           [collection_size]
                                           [num_src_fields]
                                           [src_field_buffer_sizes[src_field_idx]])
     @param[out] info             returned info argument indicating the
                                  action performed
     @param[out] ierror           returned error
     @remark source field values in the src_field_buffer have already been
             weighted, by the fractional mask
     @remark This routine returns immediately without waiting for the receive
             operation to be completed. The user should not access the memory
             associated with recv_field until this is the case. This can be
             achieved by calling \ref yac_cwait.
*/
     void yac_cget_raw_frac_async ( int const field_id,
                                    int const collection_size,
                                    double ***src_field_buffer,
                                    double ***src_frac_mask_buffer,
                                    int *info,
                                    int *ierror );

/* -------------------------------------------------------------------------------- */

/** Sending of the coupling fields

     @param[in]  field_id           -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer (all data is stored in one
                                                   contiguous part of the memory)
                                      dimensions send_field[collection_size]
                                                           [nbr_fields]
                                                           [nbr_points]
     @param[out] info               - returned info
     @param[out] ierror             - returned error
*/

     void yac_cput_ ( int const field_id,
                      int const collection_size,
                      double *send_field,
                      int    *info,
                      int    *ierror );

/** Sending of the coupling fields

     @param[in]  field_id           -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer
                                      dimensions send_field[collection_size]
                                                           [nbr_fields]
                                                           [nbr_points]
     @param[out] info               - returned info
     @param[out] ierror             - returned error
*/

     void yac_cput ( int const field_id,
                     int const collection_size,
                     double *** const send_field,
                     int *info,
                     int *ierror );

/* -------------------------------------------------------------------------------- */

/** Sending of the coupling fields

     @param[in]  field_id           -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer (all data is stored in one
                                                   contiguous part of the memory)
                                      dimensions send_field[collection_size]
                                                           [nbr_fields]
                                                           [nbr_points]
     @param[in]  send_frac_mask    - send fractional mask (all data is stored in one
                                                           contiguous part of the memory)
                                      dimensions send_frac_mask[collection_size]
                                                               [nbr_fields]
                                                               [nbr_points]
     @param[out] info               - returned info
     @param[out] ierror             - returned error
*/

     void yac_cput_frac_ ( int const field_id,
                           int const collection_size,
                           double *send_field,
                           double *send_frac_mask,
                           int    *info,
                           int    *ierror );

/** Sending of the coupling fields

     @param[in]  field_id           -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer
                                      dimensions send_field[collection_size]
                                                           [nbr_fields]
                                                           [nbr_points]
     @param[in]  send_frac_mask     - send fractional mask
                                      dimensions send_frac_mask[collection_size]
                                                               [nbr_fields]
                                                               [nbr_points]
     @param[out] info               - returned info
     @param[out] ierror             - returned error
*/

     void yac_cput_frac ( int const field_id,
                          int const collection_size,
                          double *** const send_field,
                          double *** const send_frac_mask,
                          int *info,
                          int *ierror );

/* -------------------------------------------------------------------------------- */

/** Sending of the coupling fields

     @param[in]  field_id           -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer
                                      dimensions send_field[collection_size * nbr_fields]
                                                           [nbr_points]
     @param[out] info               - returned info
     @param[out] ierror             - returned error
*/

     void yac_cput_ptr_ ( int const field_id,
                          int const collection_size,
                          double ** send_field,
                          int    *info,
                          int    *ierror );

/* -------------------------------------------------------------------------------- */

/** Sending of the coupling fields

     @param[in]  field_id           -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer
                                      dimensions send_field[collection_size * nbr_fields]
                                                           [nbr_points]
     @param[in]  send_frac_mask     - send fractional mask
                                      dimensions send_frac_mask[collection_size * nbr_fields]
                                                               [nbr_points]
     @param[out] info               - returned info
     @param[out] ierror             - returned error
*/

     void yac_cput_frac_ptr_ ( int const field_id,
                               int const collection_size,
                               double ** send_field,
                               double ** send_frac_mask,
                               int    *info,
                               int    *ierror );

/* -------------------------------------------------------------------------------- */

/** Exchange of the coupling fields

     @param[in]  send_field_id      -
     @param[in]  recv_field_id      -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer (all data is stored in one
                                                   contiguous part of the memory)
                                      dimensions send_field[collection_size]
                                                           [nbr_fields]
                                                           [nbr_points]
     @param[in]  recv_field         - receive buffer (all data is stored in one
                                                      contiguous part of the memory)
                                      dimensions: recv_field[collection_size][nbr_points]
     @param[out] send_info          - returned send info
     @param[out] recv_info          - returned recv info
     @param[out] ierror             - returned error
*/

     void yac_cexchange_ ( int const send_field_id,
                           int const recv_field_id,
                           int const collection_size,
                           double *send_field ,
                           double *recv_field,
                           int    *send_info,
                           int    *recv_info,
                           int    *ierror );


/** Exchange of the coupling fields

     @param[in]  send_field_id      -
     @param[in]  recv_field_id      -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer
                                      dimensions send_field[collection_size]
                                                           [nbr_fields]
                                                           [nbr_points]
     @param[in]  recv_field         - receive buffer
                                      dimensions: recv_field[collection_size][nbr_points]
     @param[out] send_info          - returned send info
     @param[out] recv_info          - returned recv info
     @param[out] ierror             - returned error
*/

     void yac_cexchange ( int const send_field_id,
                          int const recv_field_id,
                          int const collection_size,
                          double *** const send_field,
                          double ** recv_field,
                          int *send_info,
                          int *recv_info,
                          int *ierror );

/** Exchange of the source points required for the interpolation of the local
    target field

     @param[in]  send_field_id
     @param[in]  recv_field_id
     @param[in]  collection_size
     @param[in]  send_field       send buffer (all data is stored in one
                                               contiguous part of the memory)
                                  dimensions send_field[collection_size]
                                                       [nbr_fields]
                                                       [nbr_points]
     @param[in]  src_field_buffer source field buffer
                                  (dimensions:
                                     src_field_buffer
                                       [collection_size]
                                       [num_src_fields]
                                       [src_field_buffer_sizes[src_field_idx]])
     @param[out] send_info        returned send info
     @param[out] recv_info        returned recv info
     @param[out] ierror           returned error
*/

     void yac_cexchange_raw ( int const send_field_id,
                              int const recv_field_id,
                              int const collection_size,
                              double *** const send_field,
                              double ***src_field_buffer,
                              int    *send_info,
                              int    *recv_info,
                              int    *ierror );

/** Exchange of the source points required for the interpolation of the local
    target field

     @param[in]  send_field_id
     @param[in]  recv_field_id
     @param[in]  collection_size
     @param[in]  send_field       send buffer (all data is stored in one
                                               contiguous part of the memory)
                                  dimensions send_field[collection_size]
                                                       [nbr_fields]
                                                       [nbr_points]
     @param[in]  src_field_buffer source field buffer
                                  (dimensions:
                                     src_field_buffer
                                       [collection_size]
                                       [num_src_fields]
                                       [src_field_buffer_sizes[src_field_idx]])
     @param[out] send_info        returned send info
     @param[out] recv_info        returned recv info
     @param[out] ierror           returned error
*/

     void yac_cexchange_raw_ ( int const send_field_id,
                               int const recv_field_id,
                               int const collection_size,
                               double *send_field,
                               double *src_field_buffer,
                               int    *send_info,
                               int    *recv_info,
                               int    *ierror );

/* -------------------------------------------------------------------------------- */

/** Exchange of the source points required for the interpolation of the local
    target field

     @param[in]  send_field_id
     @param[in]  recv_field_id
     @param[in]  collection_size
     @param[in]  send_field       send buffer
                                  dimensions send_field[collection_size]
                                                       [nbr_fields]
                                                       [nbr_points]
     @param[in]  send_frac_mask   send fractional mask
                                  dimensions send_frac_mask[collection_size]
                                                           [nbr_fields]
                                                           [nbr_points]
     @param[in]  src_field_buffer source field buffer
                                  (dimensions:
                                     src_field_buffer
                                       [collection_size]
                                       [num_src_fields]
                                       [src_field_buffer_sizes[src_field_idx]])
     @param[in]  src_frac_mask_buffer source fractional mask buffer
                                      (dimensions:
                                         src_frac_mask_buffer
                                           [collection_size]
                                           [num_src_fields]
                                           [src_field_buffer_sizes[src_field_idx]])
     @param[out] send_info        returned send info
     @param[out] recv_info        returned recv info
     @param[out] ierror           returned error
     @remark source field values in the src_field_buffer have already been
             weighted, by the fractional mask
*/

     void yac_cexchange_raw_frac ( int const send_field_id,
                                   int const recv_field_id,
                                   int const collection_size,
                                   double *** const send_field,
                                   double *** const send_frac_mask,
                                   double ***src_field_buffer,
                                   double ***src_frac_mask_buffer,
                                   int    *send_info,
                                   int    *recv_info,
                                   int    *ierror );

/** Exchange of the source points required for the interpolation of the local
    target field

     @param[in]  send_field_id
     @param[in]  recv_field_id
     @param[in]  collection_size
     @param[in]  send_field       send buffer (all data is stored in one
                                               contiguous part of the memory)
                                  dimensions send_field[collection_size]
                                                       [nbr_fields]
                                                       [nbr_points]
     @param[in]  send_frac_mask   send fractional mask (all data is stored in one
                                               contiguous part of the memory)
                                  dimensions send_frac_mask[collection_size]
                                                           [nbr_fields]
                                                           [nbr_points]
     @param[in]  src_field_buffer source field buffer (all data is stored in
                                                       one contiguous part of
                                                       the memory)
                                  (dimensions:
                                     src_field_buffer
                                       [collection_size]
                                       [num_src_fields]
                                       [src_field_buffer_sizes[src_field_idx]])
     @param[in]  src_frac_mask_buffer source fractional mask buffer
                                      (all data is stored in one contiguous
                                      part of the memory)
                                      (dimensions:
                                         src_field_buffer
                                           [collection_size]
                                           [num_src_fields]
                                           [src_field_buffer_sizes[src_field_idx]])
     @param[out] send_info        returned send info
     @param[out] recv_info        returned recv info
     @param[out] ierror           returned error
     @remark source field values in the src_field_buffer have already been
             weighted, by the fractional mask
*/

     void yac_cexchange_raw_frac_ ( int const send_field_id,
                                    int const recv_field_id,
                                    int const collection_size,
                                    double *send_field,
                                    double *send_frac_mask,
                                    double *src_field_buffer,
                                    double *src_frac_mask_buffer,
                                    int    *send_info,
                                    int    *recv_info,
                                    int    *ierror );

/* -------------------------------------------------------------------------------- */

/** Exchange of the coupling fields

     @param[in]  send_field_id      -
     @param[in]  recv_field_id      -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer (all data is stored in one
                                                   contiguous part of the memory)
                                      dimensions send_field[collection_size]
                                                           [nbr_fields]
                                                           [nbr_points]
     @param[in]  send_frac_mask     - send fractional mask (all data is stored in one
                                                   contiguous part of the memory)
                                      dimensions send_frac_mask[collection_size]
                                                               [nbr_fields]
                                                               [nbr_points]
     @param[in]  recv_field         - receive buffer (all data is stored in one
                                                      contiguous part of the memory)
                                      dimensions: recv_field[collection_size][nbr_points]
     @param[out] send_info          - returned send info
     @param[out] recv_info          - returned recv info
     @param[out] ierror             - returned error
*/

     void yac_cexchange_frac_ ( int const send_field_id,
                                int const recv_field_id,
                                int const collection_size,
                                double *send_field,
                                double *send_frac_mask,
                                double *recv_field,
                                int    *send_info,
                                int    *recv_info,
                                int    *ierror );


/** Exchange of the coupling fields

     @param[in]  send_field_id      -
     @param[in]  recv_field_id      -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer
                                      dimensions send_field[collection_size]
                                                           [nbr_fields]
                                                           [nbr_points]
     @param[in]  send_frac_mask     - send fractional mask
                                      dimensions send_frac_mask[collection_size]
                                                               [nbr_fields]
                                                               [nbr_points]
     @param[in]  recv_field         - receive buffer
                                      dimensions: recv_field[collection_size][nbr_points]
     @param[out] send_info          - returned send info
     @param[out] recv_info          - returned recv info
     @param[out] ierror             - returned error
*/

     void yac_cexchange_frac ( int const send_field_id,
                               int const recv_field_id,
                               int const collection_size,
                               double *** const send_field,
                               double *** const send_frac_mask,
                               double ** recv_field,
                               int *send_info,
                               int *recv_info,
                               int *ierror );

/* -------------------------------------------------------------------------------- */

/** Sending of the coupling fields

     @param[in]  send_field_id      -
     @param[in]  recv_field_id      -
     @param[in]  collection_size    -
     @param[in]  send_field         - send buffer
                                      dimensions send_field[collection_size * nbr_fields]
                                                           [nbr_points]
     @param[in]  recv_field         - receive buffer
                                      dimensions: recv_field[collection_size][nbr_points]
     @param[out] send_info          - returned send info
     @param[out] recv_info          - returned recv info
     @param[out] ierror             - returned error
*/

     void yac_cexchange_ptr_ ( int const send_field_id,
                               int const recv_field_id,
                               int const collection_size,
                               double ** send_field,
                               double ** recv_field,
                               int    *send_info,
                               int    *recv_info,
                               int    *ierror );

/** Sending of the coupling fields

     @param[in]  send_field_id
     @param[in]  recv_field_id
     @param[in]  collection_size
     @param[in]  send_field       send buffer
                                  dimensions send_field[collection_size * nbr_fields]
                                                       [nbr_points]
     @param[in]  src_field_buffer source field buffer (all data is stored in
                                                       one contiguous part of
                                                       the memory)
                                  (dimensions:
                                     src_field_buffer
                                       [collection_size * num_src_fields]
                                       [src_field_buffer_sizes[src_field_idx]])
     @param[out] send_info        returned send info
     @param[out] recv_info        returned recv info
     @param[out] ierror           returned error
*/

     void yac_cexchange_raw_ptr_ ( int const send_field_id,
                                   int const recv_field_id,
                                   int const collection_size,
                                   double ** send_field,
                                   double ** src_field_buffer,
                                   int    *send_info,
                                   int    *recv_info,
                                   int    *ierror );

/* -------------------------------------------------------------------------------- */

/** Sending of the coupling fields

     @param[in]  send_field_id
     @param[in]  recv_field_id
     @param[in]  collection_size
     @param[in]  send_field      send buffer
                                 dimensions send_field[collection_size * nbr_fields]
                                                      [nbr_points]
     @param[in]  send_frac_mask  send buffer
                                 dimensions send_frac_mask[collection_size * nbr_fields]
                                                          [nbr_points]
     @param[in]  recv_field      receive buffer
                                 dimensions: recv_field[collection_size][nbr_points]
     @param[out] send_info       returned send info
     @param[out] recv_info       returned recv info
     @param[out] ierror          returned error
*/

     void yac_cexchange_frac_ptr_ ( int const send_field_id,
                                    int const recv_field_id,
                                    int const collection_size,
                                    double ** send_field,
                                    double ** send_frac_mask,
                                    double ** recv_field,
                                    int    *send_info,
                                    int    *recv_info,
                                    int    *ierror );

/** Sending of the coupling fields

     @param[in]  send_field_id
     @param[in]  recv_field_id
     @param[in]  collection_size
     @param[in]  send_field      send buffer
                                 dimensions send_field[collection_size * nbr_fields]
                                                      [nbr_points]
     @param[in]  send_frac_mask  send buffer
                                 dimensions send_frac_mask[collection_size * nbr_fields]
                                                          [nbr_points]
     @param[in]  src_field_buffer source field buffer (all data is stored in
                                                       one contiguous part of
                                                       the memory)
                                  (dimensions:
                                     src_field_buffer
                                       [collection_size * num_src_fields]
                                       [src_field_buffer_sizes[src_field_idx]])
     @param[in]  src_frac_mask_buffer source fractional mask buffer
                                      (all data is stored in one contiguous
                                      part of the memory)
                                      (dimensions:
                                         src_field_buffer
                                           [collection_size * num_src_fields]
                                           [src_field_buffer_sizes[src_field_idx]])
     @param[out] send_info       returned send info
     @param[out] recv_info       returned recv info
     @param[out] ierror          returned error
*/

     void yac_cexchange_raw_frac_ptr_ ( int const send_field_id,
                                        int const recv_field_id,
                                        int const collection_size,
                                        double ** send_field,
                                        double ** send_frac_mask,
                                        double ** src_field_buffer,
                                        double ** src_frac_mask_buffer,
                                        int    *send_info,
                                        int    *recv_info,
                                        int    *ierror );

/* -------------------------------------------------------------------------------- */

/** Determines whether there is an asynchronous communication associated with
    a field, which is not yet completed (for example by a previous put)

     @param[in]   field_id
     @param[out]  flag     "0" if there is an uncompleted asynchronous
                           communication associated with the field,
                           "1" otherwise
*/

     void yac_ctest ( int field_id, int * flag );

/* -------------------------------------------------------------------------------- */

/** Waits until all previous asynchronous communication associated with
    a field is completed (for example by a previous put)

     @param[in]   field_id
*/

     void yac_cwait ( int field_id );

/* -------------------------------------------------------------------------------- */

/** synchronize grids and fields of the default instance
 */
     void yac_csync_def ( void );

/** synchronize grids and fields

     @param[in]  yac_instance_id id of the YAC instance
 */
     void yac_csync_def_instance ( int yac_instance_id );

/** End of the definition phase for the default YAC instance and
    setup of internal data structures for coupling (communication matrices
    and interpolation weights)
 */
     void yac_cenddef ( void );

/** End of the definition phase, invocation of the search and
    setup of internal data structures for coupling (communication matrices
    and interpolation weights)

     @param[in]  yac_instance_id id of the YAC instance
 */
     void yac_cenddef_instance ( int yac_instance_id );

/** End of the definition phase for the default YAC instance,
    invocation of the search and
    setup of internal data structures for coupling (communication matrices
    and interpolation weights). In addition, the collected coupling configuration
    from all processes is emitted.

     @param[in]  emit_flags flags for configuring the generated coupling
                            configuration output
                            (\ref YAC_YAML_EMITTER_DEFAULT or
                             \ref YAC_YAML_EMITTER_JSON)
     @param[out] config     coupling configuration
     @remark the user is responsible for freeing the returned coupling
             configuration
 */

     void yac_cenddef_and_emit_config ( int emit_flags,
                                        char ** config );

/** End of the definition phase, invocation of the search and
    setup of internal data structures for coupling (communication matrices
    and interpolation weights). In addition, the collected coupling configuration
    from all processes is emitted.

     @param[in]  yac_instance_id id of the YAC instance
     @param[in]  emit_flags flags for configuring the generated coupling
                            configuration output
                            (\ref YAC_YAML_EMITTER_DEFAULT or
                             \ref YAC_YAML_EMITTER_JSON)
     @param[out] config     coupling configuration
     @remark the user is responsible for freeing the returned coupling
             configuration
 */
     void yac_cenddef_and_emit_config_instance ( int yac_instance_id,
                                                 int emit_flags,
                                                 char ** config);

/* --------------------------------------------------------------------------------
           query routines
   -------------------------------------------------------------------------------- */

/** query routine for the start datetime of the job of the default YAC instance
*/

     char * yac_cget_start_datetime ( void );

/** query routine for the start datetime of the job

     @param[in] yac_instance_id id of the YAC instance
*/

     char * yac_cget_start_datetime_instance ( int yac_instance_id );

/** query routine for the end datetime of the job of the default YAC instance
*/

     char * yac_cget_end_datetime ( void );

/** query routine for the end datetime of the job

     @param[in] yac_instance_id id of the YAC instance
*/

     char * yac_cget_end_datetime_instance ( int yac_instance_id );

/** query routine for the YAC version
*/

     char * yac_cget_version ( void );

/**  @brief Retrieve the MPI handshake group name used by YAC.

     This function returns (via an output pointer)
     \ref mpi_handshake_group_name.

     @return name of the MPI group used by YAC in MPI handshake.

     @note Returns a null-terminated string containing the group name.
     The string is managed internally and must **not** be modified or freed.
 */
     const char * yac_cget_mpi_handshake_group_name( void );

/* -------------------------------------------------------------------------------- */

/** query routine for number of components defined in the default YAC
 *  instance

     @return   number of components

*/
     int yac_cget_nbr_comps ( void );

/** query routine for number of components

     @param[in] yac_instance_id   id of the YAC instance
     @return                      number of components
*/
     int yac_cget_nbr_comps_instance ( int yac_instance_id );

/** query routine for number of grids defined in the default YAC
 *  instance that are referenced by at least one defined field or have
 *  metadata assigned to it

     @return                  number of fields defined on the calling process

*/
     int yac_cget_nbr_grids ( );

/** query routine for number of grids that are referenced by at least
 *  one defined field or have metadata assigned to it

     @param[in] yac_instance_id   id of the YAC instance
     @return                      number of fields defined on the calling process
*/
     int yac_cget_nbr_grids_instance ( int yac_instance_id );

/** query routine for number of grids defined in the default YAC
 *  instance for a given components that is referenced by at least
 *  one defined field.

     @param[in] comp_name
     @return                      number of components
*/
     int yac_cget_comp_nbr_grids ( const char* comp_name );

/** query routine for number of grids defined for a given components
 *  that is referenced by at least one defined field.

     @param[in] yac_instance_id   id of the YAC instance
     @param[in] comp_name
     @return                      number of components
*/
     int yac_cget_comp_nbr_grids_instance ( int yac_instance_id,
                                            const char* comp_name );

/** query routine for number of coupling fields defined in the given component
    and grids of the default YAC instance

     @param[in] comp_name     name of the component
     @param[in] grid_name     name of the grid
     @return                  number of fields defined on the calling process

*/
     int yac_cget_nbr_fields ( const char* comp_name,
                               const char* grid_name );

/** query routine for number of coupling fields defined in the given component
    and grids

     @param[in] yac_instance_id   id of the YAC instance
     @param[in] comp_name         name of the component
     @param[in] grid_name         name of the grid
     @return                      number of fields defined on the calling process

*/
     int yac_cget_nbr_fields_instance ( int yac_instance_id,
                                        const char* comp_name,
                                        const char* grid_name);

/* -------------------------------------------------------------------------------- */

/** query routine to get the list of component names for defined the
 *  default YAC instance.

     @param[in]  nbr_comps     number of components
     @param[out] comp_names    list of component names
*/
     void yac_cget_comp_names ( int nbr_comps,
                                const char ** comp_names );

/** query routine to get the list of component names for defined

     @param[in] yac_instance_id id of the YAC instance
     @param[in]  nbr_comps     number of components
     @param[out] comp_names  list of component names

*/
     void yac_cget_comp_names_instance ( int yac_instance_id,
                                         int nbr_comps,
                                         const char ** comp_names );

/** query routine to get the list of grid names of
 *  the default YAC instance that are referenced by at least
 *  one defined field or have metadata assigned to it

     @param[in]  nbr_grids     number of grids
     @param[out] grid_names    list of grid names

*/
     void yac_cget_grid_names ( int nbr_grids,
                                const char ** grid_names );

/** query routine to get the list of grid names that are referenced by at least
 *  one defined field or have metadata assigned to it

     @param[in]  yac_instance_id id of the YAC instance
     @param[in]  nbr_grids     number of grids
     @param[out] grid_names      list of grid names

*/
     void yac_cget_grid_names_instance ( int yac_instance_id,
                                         int nbr_grids,
                                         const char ** grid_names );

/** query routine to get the list of grid names of a given component
 *  that are referenced by at least one defined field or have metadata
 *  assigned to it

     @param[in]  comp_name       name of the component
     @param[in]  nbr_grids       number of grids
     @param[out] grid_names      list of grid names

*/
     void yac_cget_comp_grid_names ( const char* comp_name,
                                     int nbr_grids,
                                     const char ** grid_names );

/** query routine to get the list of grid names of a given component
 *  that are referenced by at least one defined field or have metadata
 *  assigned to it

     @param[in]  yac_instance_id id of the YAC instance
     @param[in]  comp_name       name of the component
     @param[in]  nbr_grids       number of grids
     @param[out] grid_names      list of grid names

*/
     void yac_cget_comp_grid_names_instance ( int yac_instance_id,
                                              const char* comp_name,
                                              int nbr_grids,
                                              const char ** grid_names );

/** query routine to get the list of field names defined on the given component
    and grid for the default YAC instance.

     @param[in]  comp_name     component name
     @param[in]  grid_name     grid name
     @param[in]  nbr_fields    number of fields
     @param[out] field_names   list of field names

*/
     void yac_cget_field_names ( const char* comp_name,
                                 const char* grid_name,
                                 int nbr_fields,
                                 const char ** field_names );

/** query routine to get the list of field names on the given component and grid

     @param[in]  yac_instance_id id of the YAC instance
     @param[in]  comp_name       component name
     @param[in]  grid_name       grid name
     @param[in]  nbr_fields      number of fields
     @param[out] field_names      list of grid names

*/
     void yac_cget_field_names_instance ( int yac_instance_id,
                                          const char * comp_name,
                                          const char* grid_name,
                                          int nbr_fields,
                                          const char ** field_names );

/* ---------------------------------------------------------------------- */

/** query routine to get the component name for a given ID

     @param[in] field_id       ID as provided by yac_cdef_field
     @return comp_name     component name of the fields component.
                               (array of length YAC_MAX_CHARLEN)

*/

     const char* yac_cget_component_name_from_field_id ( int field_id );

/** query routine to get the grid name for a given ID

     @param[in] field_id       ID as provided by yac_cdef_field
     @return grid_name     grid name of the fields grid.
                               (array of length YAC_MAX_CHARLEN)

*/

     const char* yac_cget_grid_name_from_field_id ( int field_id );

/** query routine to get the field name for a given ID

     @param[in] field_id       ID as provided by yac_cdef_field
     @return field_name    field name that has been given to yac_cdef_field.
                               (array of length YAC_MAX_CHARLEN)

*/

     const char* yac_cget_field_name_from_field_id ( int field_id );

/** query routine to get the timestep of a field for a given ID

     @param[in] field_id   ID as provided by yac_cdef_field
     @return timestep  timestep that has been given to yac_cdef_field
                           iso8601 string.
                           (array of length YAC_MAX_CHARLEN)

*/

     const char* yac_cget_timestep_from_field_id ( int field_id );

/** query routine to get the collection_size of a field for a given ID

     @param[in] field_id   ID as provided by yac_cdef_field
     @return collection_size

*/

     int yac_cget_collection_size_from_field_id ( int field_id );

/** query routine to get the role of a field for a given ID

     @param[in] field_id   ID as provided by yac_cdef_field
     @return role          exchange type of the field\n
                           (\ref YAC_EXCHANGE_TYPE_NONE,
                            \ref YAC_EXCHANGE_TYPE_SOURCE, or
                            \ref YAC_EXCHANGE_TYPE_TARGET)
     @remark the role is set by the call to \ref yac_cenddef or
             \ref yac_cenddef_instance

*/

     int yac_cget_role_from_field_id ( int field_id );

/** determines whether a field with the provided component, grid, and field
    name is defined on the local process

     @param[in] comp_name   component name
     @param[in] grid_name   grid name
     @param[in] field_name  field name
     @return != 0, if the local process has defined a field with the provided
             component, grid, and field name
     @return 0, otherwise
*/
     int yac_cget_field_is_defined( const char* comp_name,
                                    const char* grid_name,
                                    const char* field_name);

/** determines whether a field with the provided component, grid, and field
    name is defined on the local process

     @param[in] yac_instance_id
     @param[in] comp_name   component name
     @param[in] grid_name   grid name
     @param[in] field_name  field name
     @return != 0, if the local process has defined a field with the provided
             component, grid, and field name
     @return 0, otherwise
*/
     int yac_cget_field_is_defined_instance( int yac_instance_id,
                                              const char* comp_name,
                                              const char* grid_name,
                                              const char* field_name);

/** query routine to get the field_id from component, grid and field
    name (if defined on this process)

     @param[in] comp_name   component name
     @param[in] grid_name   grid name
     @param[in] field_name  field name
     @return field_id
*/
     int yac_cget_field_id( const char* comp_name,
                            const char* grid_name,
                            const char* field_name);

/** query routine to get the field_id from component, grid and field
    name (if defined on this process)

     @param[in] yac_instance_id
     @param[in] comp_name   component name
     @param[in] grid_name   grid name
     @param[in] field_name  field name
     @return field_id
*/
     int yac_cget_field_id_instance( int yac_instance_id,
                                     const char* comp_name,
                                     const char* grid_name,
                                     const char* field_name);

/* ---------------------------------------------------------------------- */

/** Get metadata for a component
     @param[in] comp_name
     @return metadata (NULL if no metadata is defined)
*/
     const char* yac_cget_component_metadata(const char* comp_name);

/** Get metadata for a component
     @param[in] yac_instance_id
     @param[in] comp_name
     @return metadata (NULL if no metadata is defined)
*/
     const char* yac_cget_component_metadata_instance( int yac_instance_id,
                                                       const char* comp_name);

/** Get metadata for a grid
     @param[in] grid_name
     @return metadata (NULL if no metadata is defined)
*/
     const char* yac_cget_grid_metadata(const char* grid_name);

/** Get metadata for a grid
     @param[in] yac_instance_id
     @param[in] grid_name
     @return metadata (NULL if no metadata is defined)
*/
     const char* yac_cget_grid_metadata_instance( int yac_instance_id,
                                                  const char* grid_name);

/** Get metadata for a field
     @param[in] comp_name
     @param[in] grid_name
     @param[in] field_name
     @return metadata (NULL if no metadata is defined)
*/
     const char* yac_cget_field_metadata( const char* comp_name,
                                          const char* grid_name,
                                          const char* field_name);

/** Get metadata for a field
     @param[in] yac_instance_id
     @param[in] comp_name
     @param[in] grid_name
     @param[in] field_name
     @return metadata (NULL if no metadata is defined)
*/
     const char* yac_cget_field_metadata_instance( int yac_instance_id,
                                                   const char* comp_name,
                                                   const char* grid_name,
                                                   const char* field_name);

/* ---------------------------------------------------------------------- */

/** query routine to get the timestep of a field

     @param[in] comp_name
     @param[in] grid_name
     @param[in] field_name
     @return timestep timestep that has been given to yac_cdef_field
             iso8601 string (array of length YAC_MAX_CHARLEN)

*/

     const char* yac_cget_field_timestep ( const char* comp_name,
                                           const char* grid_name,
                                           const char* field_name );

/** query routine to get the timestep of a field

     @param[in] yac_instance_id
     @param[in] comp_name
     @param[in] grid_name
     @param[in] field_name
     @return timestep timestep that has been given to yac_cdef_field
             iso8601 string (array of length YAC_MAX_CHARLEN)

*/

     const char* yac_cget_field_timestep_instance ( int yac_instance_id,
                                                    const char* comp_name,
                                                    const char* grid_name,
                                                    const char* field_name );

/** query routine to get the fractional mask fallback value of a field

     @param[in] comp_name
     @param[in] grid_name
     @param[in] field_name
     @return frac_mask_fallback_value
     @remark if YAC_FRAC_MASK_NO_VALUE is returned, then no
             fractional mask fallback value is defined for this field

*/

     double yac_cget_field_frac_mask_fallback_value( const char* comp_name,
                                                     const char* grid_name,
                                                     const char* field_name);

/** query routine to get the collection_size of a field

     @param[in] comp_name
     @param[in] grid_name
     @param[in] field_name
     @return collection_size

*/

     int yac_cget_field_collection_size( const char* comp_name,
                                         const char* grid_name,
                                         const char* field_name);

/** query routine to get the fractional mask fallback value of a field

     @param[in] yac_instance_id
     @param[in] comp_name
     @param[in] grid_name
     @param[in] field_name
     @return frac_mask_fallback_value
     @remark if \ref YAC_FRAC_MASK_NO_VALUE is returned, then no
             fractional mask fallback value is defined for this field

*/

     double yac_cget_field_frac_mask_fallback_value_instance(
       int yac_instance_id,
       const char* comp_name,
       const char* grid_name,
       const char* field_name );

/** query routine to get the collection_size of a field

     @param[in] yac_instance_id
     @param[in] comp_name
     @param[in] grid_name
     @param[in] field_name
     @return collection_size

*/

     int yac_cget_field_collection_size_instance(
       int yac_instance_id,
       const char* comp_name,
       const char* grid_name,
       const char* field_name);

/** query routine to get the role of a field for a given ID

     @param[in] comp_name
     @param[in] grid_name
     @param[in] field_name
     @return role          exchange type of the field\n
                           (\ref YAC_EXCHANGE_TYPE_NONE,
                           \ref YAC_EXCHANGE_TYPE_SOURCE, or
                           \ref YAC_EXCHANGE_TYPE_TARGET)
     @remark the role is set by the call to \ref yac_cenddef or
                           \ref yac_cenddef_instance

*/

     int yac_cget_field_role ( const char* comp_name,
                               const char* grid_name,
                               const char* field_name );

/** query routine to get the role of a field for a given ID

     @param[in] yac_instance_id
     @param[in] comp_name
     @param[in] grid_name
     @param[in] field_name
     @return role          exchange type of the field\n
                           (\ref YAC_EXCHANGE_TYPE_NONE,
                           \ref YAC_EXCHANGE_TYPE_SOURCE, or
                           \ref YAC_EXCHANGE_TYPE_TARGET)
     @remark the role is set by the call to \ref yac_cenddef or
                           \ref yac_cenddef_instance

*/

     int yac_cget_field_role_instance ( int yac_instance_id,
                                        const char* comp_name,
                                        const char* grid_name,
                                        const char* field_name );

/** query routine to get the source field for a given target field

     @param[in] tgt_comp_name
     @param[in] tgt_grid_name
     @param[in] tgt_field_name
     @param[out] src_comp_name
     @param[out] src_grid_name
     @param[out] src_field_name

*/

     void yac_cget_field_source ( const char* tgt_comp_name,
                                  const char* tgt_grid_name,
                                  const char* tgt_field_name,
                                  const char** src_comp_name,
                                  const char** src_grid_name,
                                  const char** src_field_name);

/** query routine to get the source field for a given target field

     @param[in] yac_instance_id
     @param[in] tgt_comp_name
     @param[in] tgt_grid_name
     @param[in] tgt_field_name
     @param[out] src_comp_name
     @param[out] src_grid_name
     @param[out] src_field_name

*/

     void yac_cget_field_source_instance ( int yac_instance_id,
                                           const char* tgt_comp_name,
                                           const char* tgt_grid_name,
                                           const char* tgt_field_name,
                                           const char** src_comp_name,
                                           const char** src_grid_name,
                                           const char** src_field_name);

/* --------------------------------------------------------------------------------
           auxiliary routines
   -------------------------------------------------------------------------------- */

#ifndef __GNUC__
#  define  __attribute__(x)  /*NOTHING*/
#endif

// prevent warning if compile with cython
#ifdef YAC_CYTHON
#  define  __attribute__(x)  /*NOTHING*/
#endif

/**
 * functions used as error handler must conform to this interface
 */
typedef void (*yac_abort_func)(MPI_Comm comm, const char *msg,
                               const char *source, int line)
  __attribute__((noreturn));


#ifdef YAC_CYTHON
#  undef  __attribute__
#endif

/** Calls the currently set abort handler (\ref abort_default "yac_abort_default"
    by default)
     @param[in] comm   MPI communicator used to call MPI_Abort
     @param[in] msg    message text to print
     @param[in] source string describing source file name
     @param[in] line   line number of caller
 */
     void yac_abort( MPI_Comm comm,
                     const char *msg,
                     const char *source,
                     int line)
     __attribute__((noreturn));

/** Call the \ref yac_abort function (providing the default communicator for
    the comm argument).
     @param msg    message text to print
     @param source string describing source file name
     @param line   line number of caller
 */
     void yac_abort_message( char const *msg,
                             const char *source,
                             int line);

/** Restores default abort handler
 */
     void yac_restore_default_abort_handler(void);

/**  Sets custom abort handler
     @param[in] custom_abort custom abort handler
     @remark This abort handler should call MPI_Abort if possible. Once this
             handler was called, the internal state of YAC is undefined and no
             YAC functions should be called anymore.
 */
     void yac_set_abort_handler(yac_abort_func custom_abort);

/** Gets abort handler
     @return currently set abort handler
 */
     yac_abort_func yac_get_abort_handler(void);

/** Gets default abort handler
     @return default abort handler
 */
     yac_abort_func yac_get_default_abort_handler(void);

/** Sets default MPI communicator (MPI_COMM_WORLD by default)
     @param[in] comm default MPI communicator
 */
     void yac_set_default_comm(MPI_Comm comm);

/* -------------------------------------------------------------------------------- */

/** query routine to get the number of points in a grid

     @param[in] location
     @param[in] grid_id

*/

     size_t yac_cget_grid_size ( int location,
                                 int grid_id );

/* -------------------------------------------------------------------------------- */

/** computes the areas of all cells in a grid (an earth radius of 1.0 is assumed)

     @param[in]  grid_id
     @param[out] cell_areas

*/

     void yac_ccompute_grid_cell_areas ( int grid_id,
                                         double * cell_areas );


/* -------------------------------------------------------------------------------- */

/** query routine to get the number of points in a pointset

     @param[in] points_id

*/

     size_t yac_cget_points_size ( int points_id  );

/* -------------------------------------------------------------------------------- */

/** Get an extended source point mapping interpolation configuration
     @param[out] ext_spmap_config_id extended source point mapping
                                     interpolation configuration
     @see \ref spmap_ext_config
 */
     void yac_cget_ext_spmap_config(int * ext_spmap_config_id);

/** Frees an extended source point mapping interpolation configuration
     @param[in] ext_spmap_config_id
     @see \ref spmap_ext_config
 */
     void yac_cfree_ext_spmap_config(int ext_spmap_config_id);

/** Sets the spread distance
     @param[out] ext_spmap_config_id extended source point mapping
                                     interpolation configuration
     @param[in]  spread_distance     spread distance
     @see \ref spmap_ext_config
 */
     void yac_cset_ext_spmap_config_spread_distance(
          int ext_spmap_config_id, double spread_distance);

/** Sets the maximum search distance
     @param[out] ext_spmap_config_id extended source point mapping
                                     interpolation configuration
     @param[in]  max_search_distance maximum search distance
     @see \ref spmap_ext_config
 */
     void yac_cset_ext_spmap_config_max_search_distance(
          int ext_spmap_config_id, double max_search_distance);

/** Sets the weightening type
     @param[out] ext_spmap_config_id extended source point mapping
                                     interpolation configuration
     @param[in]  weight_type         weightening type
     @see \ref spmap_ext_config
 */
     void yac_cset_ext_spmap_config_weight_type(
          int ext_spmap_config_id, int weight_type);

/** Sets the scaling type
     @param[out] ext_spmap_config_id extended source point mapping
                                     interpolation configuration
     @param[in]  scale_type          scaling type
     @see \ref spmap_ext_config
 */
     void yac_cset_ext_spmap_config_scale_type(
          int ext_spmap_config_id, int scale_type);

/** enables the computation of the source cell areas by yac for the scaling
     @param[out] ext_spmap_config_id extended source point mapping
                                     interpolation configuration
     @param[in]  sphere_radius       sphere radius
     @see \ref spmap_ext_config
 */
     void yac_cset_ext_spmap_config_src_cell_area_config_yac(
          int ext_spmap_config_id, double sphere_radius);

/** enables the computation of the target cell areas by yac
     @param[out] ext_spmap_config_id extended source point mapping
                                     interpolation configuration
     @param[in]  sphere_radius       sphere radius
     @see \ref spmap_ext_config
 */
     void yac_cset_ext_spmap_config_tgt_cell_area_config_yac(
          int ext_spmap_config_id, double sphere_radius);

/** enables the reading of the source cell areas for the scaling
     @param[out] ext_spmap_config_id extended source point mapping
                                     interpolation configuration
     @param[in]  filename            netCDF file contain the cell areas
     @param[in]  varname             variable name of cell areas
     @param[in]  min_global_id       minimum global id
     @see \ref spmap_ext_config
 */
     void yac_cset_ext_spmap_config_src_cell_area_config_file(
          int ext_spmap_config_id, char const * filename,
          char const * varname, int min_global_id);

/** enables the reading of the target cell areas for the scaling
     @param[out] ext_spmap_config_id extended source point mapping
                                     interpolation configuration
     @param[in]  filename            netCDF file contain the cell areas
     @param[in]  varname             variable name of cell areas
     @param[in]  min_global_id       minimum global id
     @see \ref spmap_ext_config
 */
     void yac_cset_ext_spmap_config_tgt_cell_area_config_file(
          int ext_spmap_config_id, char const * filename,
          char const * varname, int min_global_id);


/** Get an source point mapping interpolation overwrite configuration
     @param[out] spmap_overwrite_config_id source point mapping interpolation
                                           overwrite configuration
     @see \ref spmap_ext_config
 */
     void yac_cget_spmap_overwrite_config(int * spmap_overwrite_config_id);

/** Frees an source point mapping interpolation overwrite configuration
     @param[in] spmap_overwrite_config_id source point mapping interpolation
                                          overwrite configuration
     @see \ref spmap_ext_config
 */
     void yac_cfree_spmap_overwrite_config(int spmap_overwrite_config_id);

/** Sets the source point selection method to a bounding circle
     @param[out] spmap_overwrite_config_id source point mapping interpolation
                                           overwrite configuration
     @param[in]  center_lon longitude coordinate of the center of the
                            bounding circle (in radians)
     @param[in]  center_lat latitude coordinate of the center of the
                            bounding circle (in radians)
     @param[in]  inc_angle  the angle between center vector and a vector
                            pointing to any point on the bounding circle
                            (in radians)
     @see \ref spmap_ext_config
 */
     void yac_cset_spmap_overwrite_config_src_point_selection_bnd_circle(
          int spmap_overwrite_config_id, double center_lon, double center_lat,
          double inc_angle);

/** Sets the spread distance
     @param[out] spmap_overwrite_config_id source point mapping interpolation
                                           overwrite configuration
     @param[in]  spread_distance           spread distance
     @see \ref spmap_ext_config
 */
     void yac_cset_spmap_overwrite_config_spread_distance(
          int spmap_overwrite_config_id, double spread_distance);

/** Sets the maximum search distance
     @param[out] spmap_overwrite_config_id source point mapping interpolation
                                           overwrite configuration
     @param[in]  max_search_distance       maximum search distance
     @see \ref spmap_ext_config
 */
     void yac_cset_spmap_overwrite_config_max_search_distance(
          int spmap_overwrite_config_id, double max_search_distance);

/** Sets the weightening type
     @param[out] spmap_overwrite_config_id source point mapping interpolation
                                           overwrite configuration
     @param[in]  weight_type               weightening type
     @see \ref spmap_ext_config
 */
     void yac_cset_spmap_overwrite_config_weight_type(
          int spmap_overwrite_config_id, int weight_type);

/* -------------------------------------------------------------------------------- */

/** gets an empty stack trace
     @param[out] interp_stack_config_id interpolation stack
*/
     void yac_cget_interp_stack_config(int * interp_stack_config_id);

/** frees a stack trace
     @param[in] interp_stack_config_id interpolation stack
*/
     void yac_cfree_interp_stack_config(int interp_stack_config_id);

/** adds average interpolation to the bottom of an interpolation stack
     @param[in] interp_stack_config_id interpolation stack
     @param[in] reduction_type         reduction type
     @param[in] partial_coverage       allow partial coverage
*/
     void yac_cadd_interp_stack_config_average( int interp_stack_config_id,
                                                int reduction_type,
                                                int partial_coverage);

/** adds nearest corner cells interpolation to the bottom of a
    interpolation stack
     @param[in] interp_stack_config_id interpolation stack
     @param[in] weight_type            reduction type
     @param[in] partial_coverage       allow partial coverage
*/
     void yac_cadd_interp_stack_config_ncc( int interp_stack_config_id,
                                            int weight_type,
                                            int partial_coverage);

/** adds n-nearest-neighbour interpolation to the bottom of a
    interpolation stack
     @param[in] interp_stack_config_id interpolation stack
     @param[in] type                   reduction type
     @param[in] n                      number of nearest neighbour points
     @param[in] max_search_distance    maximum search distance for each point
     @param[in] scale                  scale parameter required by some
                                       reduction types
     @remark a max_search_distance of 0.0 results in the search distance not
             being restricted
*/
     void yac_cadd_interp_stack_config_nnn( int interp_stack_config_id,
                                            int type,
                                            size_t n,
                                            double max_search_distance,
                                            double scale);

/** adds rbf interpolation to the bottom of an interpolation stack
     @param[in] interp_stack_config_id interpolation stack
     @param[in] n                      number of nearest neighbour points
     @param[in] max_search_distance    maximum search distance for each point
     @param[in] scale                  scale parameter
     @remark a max_search_distance of 0.0 results in the search distance not
             being restricted
*/
     void yac_cadd_interp_stack_config_rbf( int interp_stack_config_id,
                                            size_t n,
                                            double max_search_distance,
                                            double scale);

/** adds conservative interpolation to the bottom of an interpolation stack
     @param[in] interp_stack_config_id interpolation stack
     @param[in] order                  first or second order
     @param[in] enforced_conserv       enforce local conservation
     @param[in] partial_coverage       allow partial coverage
     @param[in] normalisation          normalisation type
*/
     void yac_cadd_interp_stack_config_conservative(
       int interp_stack_config_id,
       int order,
       int enforced_conserv,
       int partial_coverage,
       int normalisation);

/** adds source point mapping interpolation to the bottom
    of an interpolation stack
     @param[in] interp_stack_config_id interpolation stack
     @param[in] spread_distance        spread distance (in rad)
     @param[in] max_search_distance    maximum search distance (in rad)
     @param[in] weight_type            reduction type
     @param[in] scale_type             scaling type
     @param[in] src_sphere_radius      sphere radius used for
                                       source cell area computation
     @param[in] src_filename           netCDF file contain the source cell areas
     @param[in] src_varname            variable name of source cell areas
     @param[in] src_min_global_id      minimum global id used to map global cell
                                       ids to values in the netCDF file
     @param[in] tgt_sphere_radius      sphere radius used for
                                       target cell area computation
     @param[in] tgt_filename           netCDF file contain the target cell areas
     @param[in] tgt_varname            variable name of target cell areas
     @param[in] tgt_min_global_id      minimum global id used to map global cell
                                       ids to values in the netCDF file
     @remark cell areas for the scaling can either be computed by YAC or read
             from a netCDF file (see \ref interp_method_spmap)
     @remark for computing cell areas:<br>
             (sphere_radius > 0.0; filename == NULL; varname == NULL)
     @remark for reading cell areas:<br>
             (sphere_radius == 0.0; filename != NULL; varname != NULL)
*/
     void yac_cadd_interp_stack_config_spmap(
       int interp_stack_config_id,
       double spread_distance,
       double max_search_distance,
       int weight_type,
       int scale_type,
       double src_sphere_radius,
       char const * src_filename,
       char const * src_varname,
       int src_min_global_id,
       double tgt_sphere_radius,
       char const * tgt_filename,
       char const * tgt_varname,
       int tgt_min_global_id);

/** Adds source point mapping interpolation with extended configuration options
    to the bottom of an interpolation stack
     @param[in] interp_stack_config_id       interpolation stack
     @param[in] ext_spmap_config_id          main configuration
     @param[in] spmap_overwrite_config_ids   overwrite configurations
     @param[in] spmap_overwrite_config_count number of overwrite configurations
     @see \ref spmap_ext_config
*/
     void yac_cadd_interp_stack_config_spmap_ext(
       int interp_stack_config_id,
       int ext_spmap_config_id,
       int * spmap_overwrite_config_ids,
       int spmap_overwrite_config_count);

/** adds HCSBB interpolation to the bottom of an interpolation stack
     @param[in] interp_stack_config_id interpolation stack
*/
     void yac_cadd_interp_stack_config_hcsbb(int interp_stack_config_id);

/** adds user weight file interpolation to the bottom
    of an interpolation stack
     @param[in] interp_stack_config_id interpolation stack
     @param[in] filename               weight file name
     @remark YAC will read the weight file in parallel.
             This parallel output can be configured as described here:
             \ref io_config_detail
*/
     void yac_cadd_interp_stack_config_user_file(
       int interp_stack_config_id, char const * filename);

/** adds user weight file interpolation to the bottom
    of an interpolation stack
     @param[in] interp_stack_config_id interpolation stack
     @param[in] filename               weight file name
     @param[in] on_missing_file        specifies how YAC should behave if no
                                       file was found
     @param[in] on_success             specifies how YAC should behave in case
                                       a weight file was successfully read
     @remark YAC will read the weight file in parallel.
             This parallel output can be configured as described here:
             \ref io_config_detail
*/
     void yac_cadd_interp_stack_config_user_file_2(
       int interp_stack_config_id, char const * filename,
       int on_missing_file, int on_success);

/** adds fixed interpolation to the bottom of an interpolation stack
     @param[in] interp_stack_config_id interpolation stack
     @param[in] value                  fixed value
*/
     void yac_cadd_interp_stack_config_fixed( int interp_stack_config_id,
                                              double value);

/** adds fixed interpolation to the bottom of an interpolation stack
     @param[in] interp_stack_config_id interpolation stack
     @param[in] constructor_key        key provided to
                                       \ref yac_interp_method_check_add_constructor_callback
                                       for a constructor callback routine
     @param[in] do_search_key          key provided to
                                       \ref yac_interp_method_check_add_do_search_callback
                                       for a do_search callback routine
*/
     void yac_cadd_interp_stack_config_check( int interp_stack_config_id,
                                              char const * constructor_key,
                                              char const * do_search_key);

/** adds creep interpolation to the bottom of an interpolation stack
     @param[in] interp_stack_config_id interpolation stack
     @param[in] creep_distance         creep distance
*/
     void yac_cadd_interp_stack_config_creep( int interp_stack_config_id,
                                              int creep_distance);

/** adds user callback interpolation to the bottom of an interpolation stack
     @param[in] interp_stack_config_id   interpolation stack
     @param[in] func_compute_weights_key key provided to
                                         \ref yac_cadd_compute_weights_callback
                                         for a compute_weights callback routine
*/
     void yac_cadd_interp_stack_config_user_callback(
       int interp_stack_config_id,
       char const * func_compute_weights_key);

/** method signature for weight computation function used by
    \ref interp_method_callback

  If the field value of a target point, whose weights have been computed
  by this routine, is to be computed on the source process at which this
  routine was called, YAC would use the weights as follows (it is assumed
  that all required source points are available at this process):

  \code{.c}

  // this routine return an index to a source field value associated to
  // the provided global id
  size_t global_id_2_local_idx(int src_field_idx, int global_id);

  [...]

  double tgt_value[collection_size];
  for (int collection_idx = 0; collection_idx < collection_size;
       collection_idx)
    tgt_value[collection_idx] = 0.0;

  // for all source fields
  for (int src_field_idx = 0; src_field_idx < num_pointsets; ++src_field_idx) {

    for (size_t i = 0; i < result_count[src_field_idx]; ++i) {

      size_t src_field_point_idx =
        global_id_2_local_idx(
          src_field_idx, global_results_points[src_field_idx][i]);
      double weight = result_weights[src_field_idx][src_field_point_idx];

      for (int collection_idx = 0; collection_idx < collection_size;
           collection_idx) {

        tgt_value[collection_idx] +=
          src_field[collection_idx][src_field_idx][src_field_point_idx] *
          weight;
      }
    }
  }

  \endcode

  If one or more source points associated to the global ids return by this
  routine are not actually available at this process, YAC will ensure that
  the respective value will be available for the interpolation.

  @param[in]  tgt_coords            3D coordinates of the target point
  @param[in]  src_cell_id           global id of the source cell matching
                                    the target point
  @param[in]  src_cell_idx          (zero based) index of the source
                                    cell matching the target point
  @param[out] global_results_points global ids of source points to be used
                                    for the interpolation of the target point
  @param[out] result_weights        weights to be used for the interpolation
                                    of the target point
  @param[out] result_count          number of source points to be used
                                    for the interpolation of the target point
                                    per pointset
  @param[in]  user_data             user_data pointer provided to
                                    \ref yac_cadd_compute_weights_callback
  @remark The size of the first dimension of the arrays global_results_points,
          result_weights, and result_count has to be equal to num_pointsets
          provided to \ref yac_cdef_field or \ref yac_cdef_field_mask that was
          used to define the field involved in the coupling that is
          supposed to use this weight computation function.
  @remark 3D coordinates are for a point of a unit sphere
*/

#ifndef TYPEDEF_YAC_FUNC_COMPUTE_WEIGHTS
#define TYPEDEF_YAC_FUNC_COMPUTE_WEIGHTS

// Remark: make sure that this typedef is consistent with the one in interp_method_callback.h
typedef void (*yac_func_compute_weights)(
  double const tgt_coords[3], int src_cell_id, size_t src_cell_idx,
  int const ** global_results_points, double ** result_weights,
  size_t * result_count, void * user_data);
#endif

/** registers a callback routine for the computation of weights by the
 *  interpolation method user_callback

     @param[in] compute_weights_callback pointer to a weight computation routine
     @param[in] user_data                data pointer that will be passed to
                                         compute_weights_callback
     @param[in] key                      key for identifying the callback routine
     @remark the callback has to be set on all source processes taking part in
             respective coupling
     @remark the key has to match the one provided in coupling configuration
*/
     void yac_cadd_compute_weights_callback(
       yac_func_compute_weights compute_weights_callback,
       void * user_data,
       char const * key);

/** generate an interpolation stack from a "0" terminated string that contains a
 *  YAML formated description of the stack

     @param[in]  interp_stack_config    interpolation stack description
     @param[out] interp_stack_config_id interpolation stack
*/
     void yac_cget_interp_stack_config_from_string_yaml(
       char const * interp_stack_config, int * interp_stack_config_id);

/** generate an interpolation stack from a "0" terminated string that contains a
 *  JSON formated description of the stack

     @param[in]  interp_stack_config    interpolation stack description
     @param[out] interp_stack_config_id interpolation stack
*/
     void yac_cget_interp_stack_config_from_string_json(
       char const * interp_stack_config, int * interp_stack_config_id);

#endif // YAC_H
