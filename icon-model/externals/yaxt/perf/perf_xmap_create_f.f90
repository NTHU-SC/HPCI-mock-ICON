!>
!! @file perf_xmap_create_f.f90
!!
!! Demonstrate performance impact of various choices in xmap creation.
!!
!! @copyright Copyright  (C)  2024 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Thomas Jahns <jahns@dkrz.de>
!!
!
! Keywords:
! Maintainer: Thomas Jahns <jahns@dkrz.de>
! URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
#include "fc_feature_defs.inc"
#define XT_UNUSED(x) IF (SIZE( (/(x)/) ) < 0) CONTINUE
MODULE perf_xmap_create_config
  USE iso_c_binding
  USE xt_xmap_intersection, ONLY: xt_com_list
  USE yaxt, ONLY: xi => xt_int_kind, xt_xmdd_bucket_gen, &
       xt_mpi_fint_kind, xt_idxlist, xt_config, xt_stripe, &
       xt_idxlist_get_max_index, xt_idxlist_get_min_index, &
       xt_idxlist_get_num_indices, xt_idxlist_delete, &
       xt_idxlist_is_stripe_conversion_profitable, &
       xt_idxsection_new, xt_idxstripes_new, xt_idxvec_new, &
       xt_is_null, xt_xmdd_bucket_gen_new, &
       xt_config_set_xmdd_bucket_gen, xt_config_set_sort_algorithm_by_id, &
       xt_config_set_mem_saving, xt_config_set_xmap_stripe_align, &
       xt_xmdd_bucket_gen_define_interface, &
       xt_xmdd_bucket_gen_comms
  ! older PGI compilers do not handle generic interface correctly
#if defined __PGI && (__PGIC__ < 12 || (__PGIC__ ==  12 && __PGIC_MINOR__ <= 10))
  USE xt_idxlist_abstract, ONLY: xt_is_null
#endif
  USE mpi
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER :: pi8 =  14
  INTEGER, PUBLIC, PARAMETER :: i8 = SELECTED_INT_KIND(pi8)
  INTEGER, PUBLIC, PARAMETER :: &
       INDEX_GEN_COL_MAJOR = 0, &
       INDEX_GEN_ROW_MAJOR = 1, &
       INDEX_GEN_IDXVEC = 2, &
       INDEX_GEN_IDXSTRIPES = 3, &
       INDEX_GEN_IDXSECTION = 4
  TYPE config_settings
    TYPE(xt_xmdd_bucket_gen) :: bg
    INTEGER(xi) :: num_columns
    INTEGER(xi) :: nlev
    INTEGER :: comm_world_size, comm_world_rank, local_size
    INTEGER :: index_generation_sequence, &
         index_generation_method, &
         bucket_generation_sequence, &
         bucket_generation_method
    LOGICAL :: is_src
#ifdef HAVE_FC_ALLOCATABLE_CHARACTER
    CHARACTER(:), ALLOCATABLE :: fname_pat_src, fname_pat_dst
#else
    CHARACTER(len=1024) :: fname_pat_src, fname_pat_dst
#endif
  END TYPE config_settings
  TYPE(config_settings), SAVE :: run_config
  PUBLIC :: run_config, config_settings
  PUBLIC :: generate_3d_list
  PUBLIC :: set_custom_bucket_gen

  TYPE, BIND(c) :: custom_bucket_gen_state
    TYPE(xt_xmdd_bucket_gen_comms) :: comms
    INTEGER(xt_mpi_fint_kind) :: prev_rank_generated, &
         last_rank_generated, intra_comm_size
    TYPE(Xt_idxlist) :: prev_list
  END TYPE custom_bucket_gen_state

CONTAINS

  FUNCTION start_of_rank(rank, num_ranks, num_columns) RESULT(start_idx)
    INTEGER, INTENT(in) :: rank, num_ranks
    INTEGER(xi), INTENT(in) :: num_columns
    INTEGER(xi) :: start_idx
    start_idx = INT(INT(num_columns, i8) * INT(rank, i8) &
         &           / INT(num_ranks, i8), xi)
  END FUNCTION start_of_rank

  FUNCTION generate_3d_list(local_rank, local_size_, &
       index_generation_method, index_generation_sequence) &
       RESULT(idxlist_3d)
    INTEGER, INTENT(in) :: local_rank, local_size_, index_generation_method, &
         index_generation_sequence
    TYPE(xt_idxlist) :: idxlist_3d
    INTERFACE
      SUBROUTINE posix_abort() BIND(c, name='abort')
      END SUBROUTINE posix_abort
    END INTERFACE
    INTEGER(xi) :: start_idx, next_start_idx, local_num_columns, lev, idx
    INTEGER(xi) :: num_stripes, j
    INTEGER(xi), ALLOCATABLE :: idx3Dv(:)
    TYPE(xt_stripe), ALLOCATABLE :: idx3Ds(:)

    start_idx = start_of_rank(local_rank, local_size_, &
         INT(run_config%num_columns, xi))

    next_start_idx = start_of_rank(local_rank+1, local_size_, &
         INT(run_config%num_columns, xi))

    local_num_columns = next_start_idx - start_idx

    IF (index_generation_method == INDEX_GEN_IDXVEC) THEN
      ALLOCATE(idx3Dv(local_num_columns * run_config%nlev))

      IF (index_generation_sequence == INDEX_GEN_ROW_MAJOR) THEN
        j = 1_xi
        DO lev = 1_xi, run_config%nlev
          DO idx = 1_xi, local_num_columns
            idx3Dv(j) = start_idx + (idx-1_xi) + (lev-1_xi) * run_config%num_columns
            j = j + 1_xi
          END DO
        END DO
      ELSE IF (index_generation_sequence == INDEX_GEN_COL_MAJOR) THEN
        j = 1_xi
        DO idx = 1_xi, local_num_columns
          DO lev = 1_xi, run_config%nlev
            idx3Dv(j) = start_idx + (idx-1_xi) + (lev-1_xi) * run_config%num_columns
            j = j + 1_xi
          END DO
        END DO
      END IF
      idxlist_3d = xt_idxvec_new(idx3Dv)
    ELSE IF (index_generation_method == INDEX_GEN_IDXSTRIPES) THEN
      num_stripes = MERGE(run_config%nlev, local_num_columns, &
           index_generation_sequence == INDEX_GEN_ROW_MAJOR)
      ALLOCATE(idx3Ds(num_stripes))

      IF (index_generation_sequence == INDEX_GEN_ROW_MAJOR) THEN
        DO lev = 1_xi, run_config%nlev
          idx3Ds(lev) = xt_stripe(&
               start = start_idx + (lev-1_xi) * run_config%num_columns, &
               stride = 1_xi, &
               nstrides = INT(local_num_columns))
        END DO
      ELSE IF (index_generation_sequence == INDEX_GEN_COL_MAJOR) THEN
        DO idx = 1_xi, local_num_columns
          idx3Ds(idx) = xt_stripe(&
               start = start_idx + (idx-1_xi), &
               stride = run_config%num_columns, &
               nstrides = INT(run_config%nlev))
        END DO
      END IF
      idxlist_3d = xt_idxstripes_new(idx3Ds)
    ELSE ! IF (index_generation_method == INDEX_GEN_IDXSECTION)
      IF (index_generation_sequence == INDEX_GEN_ROW_MAJOR) THEN
        idxlist_3d = xt_idxsection_new(0_xi, &
             (/  run_config%nlev, run_config%num_columns /), &
             (/  INT(run_config%nlev), INT(local_num_columns) /), &
             (/  0_xi, start_idx /))
      ELSE
        ! currently unsupported
        CALL posix_abort()
      END IF
    END IF
  END FUNCTION generate_3d_list

  FUNCTION rank_of_column(column, num_ranks, num_columns)
    INTEGER, INTENT(in) :: num_ranks
    INTEGER(xi), INTENT(in) :: column, num_columns
    INTEGER :: rank_of_column
    rank_of_column = INT((INT(column, c_long_long) &
         &              * INT(num_ranks, c_long_long)) &
         &              / INT(num_columns, c_long_long))
  END FUNCTION rank_of_column

  FUNCTION custom_bucket_gen_init_state(gen_state_, src_idxlist, &
       dst_idxlist, config, comms, init_params) BIND(c) RESULT(stripify)

    INTEGER(c_int) :: stripify
    TYPE(c_ptr), VALUE :: gen_state_
    TYPE(xt_idxlist), VALUE, INTENT(in) :: src_idxlist, dst_idxlist
    TYPE(xt_config), VALUE, INTENT(in) :: config
    TYPE(xt_xmdd_bucket_gen_comms), INTENT(in) :: comms
    TYPE(c_ptr), VALUE, INTENT(in) :: init_params
    TYPE(xt_idxlist) :: init_prev_list

    TYPE(custom_bucket_gen_state), POINTER :: gen_state
    INTEGER :: ierror, start_rank, end_rank
    LOGICAL :: p
    INTEGER :: num_indices_src, num_indices_dst
    INTEGER(xi) :: first_column, last_column, start_src, end_src, &
         start_dst, end_dst

    XT_UNUSED(init_params)
    CALL C_F_POINTER(gen_state_, gen_state)
    gen_state%comms = comms
    CALL MPI_Comm_size(comms%intra_comm, gen_state%intra_comm_size, ierror)
    gen_state%prev_list = init_prev_list
    IF (run_config%index_generation_method == INDEX_GEN_IDXSTRIPES) THEN
      p = .TRUE.
    ELSE IF (xt_idxlist_is_stripe_conversion_profitable(&
         src_idxlist, config) > 0) THEN
      p = .TRUE.
    ELSE IF (xt_idxlist_is_stripe_conversion_profitable(&
         dst_idxlist, config) > 0) THEN
      p = .TRUE.
    ELSE
      p = .FALSE.
    END IF
    CALL mpi_allreduce(mpi_in_place, p, 1, mpi_logical, MPI_LOR, &
         comms%intra_comm, ierror)
    stripify = MERGE(1_c_int, 0_c_int, p)
    ! establish local ranges
    num_indices_src = xt_idxlist_get_num_indices(src_idxlist)
    num_indices_dst = xt_idxlist_get_num_indices(dst_idxlist)
    first_column = HUGE(first_column)
    last_column = -HUGE(last_column)
    IF (num_indices_src > 0) THEN
      start_src = xt_idxlist_get_min_index(src_idxlist)
      end_src = xt_idxlist_get_max_index(src_idxlist) &
           - (run_config%nlev - 1_xi) * INT(run_config%num_columns, xi)
      IF (start_src < first_column) first_column = start_src
      IF (end_src > last_column) last_column = end_src
    END IF
    IF (num_indices_dst > 0) THEN
      start_dst = xt_idxlist_get_min_index(dst_idxlist)
      end_dst = xt_idxlist_get_max_index(dst_idxlist) &
           - (run_config%nlev - 1_xi) * run_config%num_columns
      IF (start_dst < first_column) first_column = start_dst
      IF (end_dst > last_column) last_column = end_dst
    END IF
    start_rank = rank_of_column(first_column, gen_state%intra_comm_size, &
         INT(run_config%num_columns, xi))
    end_rank = rank_of_column(last_column, gen_state%intra_comm_size, &
         INT(run_config%num_columns, xi))
    gen_state%prev_rank_generated = start_rank-1
    gen_state%last_rank_generated = end_rank
  END FUNCTION custom_bucket_gen_init_state


  SUBROUTINE custom_bucket_gen_destroy_state(gen_state_) BIND(c)
    TYPE(c_ptr), VALUE, INTENT(in) :: gen_state_
    TYPE(custom_bucket_gen_state), POINTER :: gen_state

    CALL C_F_POINTER(gen_state_, gen_state)

    IF (.NOT. xt_is_null(gen_state%prev_list)) &
         CALL xt_idxlist_delete(gen_state%prev_list)
  END SUBROUTINE custom_bucket_gen_destroy_state

  FUNCTION custom_bucket_gen_get_intersect_max_num(&
       gen_state_, bucket_type) &
       RESULT(max_num) BIND(c)
    INTEGER(c_int) :: max_num
    TYPE(c_ptr), VALUE :: gen_state_
    INTEGER(c_int), VALUE :: bucket_type
    TYPE(custom_bucket_gen_state), POINTER :: gen_state

    XT_UNUSED(bucket_type)
    CALL C_F_POINTER(gen_state_, gen_state)
    max_num = gen_state%last_rank_generated - gen_state%prev_rank_generated
  END FUNCTION custom_bucket_gen_get_intersect_max_num

  FUNCTION custom_bucket_gen_next(gen_state_, bucket_type) &
       BIND(C) RESULT(bucket)
    TYPE(c_ptr), VALUE :: gen_state_
    INTEGER(c_int), VALUE :: bucket_type
    TYPE(xt_com_list) :: bucket
    TYPE(custom_bucket_gen_state), POINTER :: gen_state
    INTEGER :: next_rank

    XT_UNUSED(bucket_type)
    CALL C_F_POINTER(gen_state_, gen_state)
    IF (.NOT. xt_is_null(gen_state%prev_list)) &
         CALL xt_idxlist_delete(gen_state%prev_list)
    next_rank = gen_state%prev_rank_generated + 1
    gen_state%prev_rank_generated = next_rank
    bucket%rank = next_rank
    IF (next_rank <= gen_state%last_rank_generated) THEN
      bucket%list = generate_3d_list(next_rank, gen_state%intra_comm_size, &
           run_config%bucket_generation_method, &
           run_config%bucket_generation_sequence)
    END IF
    gen_state%prev_list = bucket%list
  END FUNCTION custom_bucket_gen_next

  SUBROUTINE set_custom_bucket_gen(config)
    TYPE(xt_config), INTENT(inout) :: config
    TYPE(custom_bucket_gen_state) :: dummy(2)
    INTEGER(mpi_address_kind) :: addr(2)
    INTEGER(c_size_t) :: sizeof_cbgs
    INTEGER :: ierror, i
    DO i = 1, 2
      CALL mpi_get_address(dummy(i), addr(i), ierror)
    END DO
    sizeof_cbgs = INT(addr(2) - addr(1), c_size_t)
    run_config%bg = xt_xmdd_bucket_gen_new()
    CALL xt_xmdd_bucket_gen_define_interface( &
         run_config%bg, &
         custom_bucket_gen_init_state, &
         custom_bucket_gen_destroy_state, &
         custom_bucket_gen_get_intersect_max_num, &
         custom_bucket_gen_next, &
         sizeof_cbgs, c_null_ptr)
    CALL xt_config_set_xmdd_bucket_gen(config, run_config%bg)
  END SUBROUTINE set_custom_bucket_gen

END MODULE perf_xmap_create_config

MODULE perf_xmap_create_init
  USE yaxt
  USE mpi
  USE iso_c_binding, ONLY: c_int, c_char, c_null_char
  USE perf_xmap_create_config, ONLY: config_settings, set_custom_bucket_gen, &
       index_gen_col_major, index_gen_row_major, index_gen_idxvec, &
       index_gen_idxstripes, index_gen_idxsection
  IMPLICIT NONE
  PRIVATE
  INTERFACE
    SUBROUTINE posix_exit(status) BIND(c, name='exit')
      IMPORT :: c_int
      INTEGER(c_int), VALUE, INTENT(in) :: status
    END SUBROUTINE posix_exit
  END INTERFACE
  PUBLIC :: parse_options
  PUBLIC :: read_idxlist
CONTAINS
  SUBROUTINE set_list_gen_methods(generation_sequence, generation_method, &
       arg)
    INTEGER, INTENT(inout) :: generation_sequence, generation_method
    CHARACTER(len=*), INTENT(in) :: arg

    INTEGER :: cur_parm_ofs, cur_parm_len, cur_parm_end

    cur_parm_ofs = 1
    ! split arg into comma-separated parts
    DO WHILE (cur_parm_ofs <= LEN(arg))
      cur_parm_len = SCAN(arg(cur_parm_ofs:), ",")
      IF (cur_parm_len /= 0) THEN
        cur_parm_end = cur_parm_ofs+cur_parm_len-2
      ELSE
        cur_parm_end = LEN(arg)
        cur_parm_len = cur_parm_end - cur_parm_ofs + 1
      END IF
      SELECT CASE (arg(cur_parm_ofs:cur_parm_end))
      CASE ("col-major")
        generation_sequence = INDEX_GEN_COL_MAJOR
      CASE ("row-major")
        generation_sequence = INDEX_GEN_ROW_MAJOR
      CASE ("idxvec")
        generation_method = INDEX_GEN_IDXVEC
      CASE ("idxstripes")
        generation_method = INDEX_GEN_IDXSTRIPES
      CASE ("idxsection")
        generation_method = INDEX_GEN_IDXSECTION
      CASE default
        WRITE (0, '(2a)') "error: unknown index generation method name: ", &
             arg(cur_parm_ofs:cur_parm_end)
        CALL posix_exit(1)
      END SELECT
      cur_parm_ofs = cur_parm_end &
           + MERGE(2, 1, cur_parm_end /= LEN(arg))
    END DO
  END SUBROUTINE set_list_gen_methods

  SUBROUTINE parse_options(run_config, config)
    TYPE(config_settings), INTENT(inout) :: run_config
    TYPE(xt_config), INTENT(inout) :: config

    INTERFACE
      FUNCTION xt_sort_algo_id_by_name(name) BIND(c) RESULT(rc)
        IMPORT :: c_int, c_char
        INTEGER(c_int) :: rc
        CHARACTER(kind=c_char, len=1), INTENT(in) :: name(*)
      END FUNCTION xt_sort_algo_id_by_name
    END INTERFACE
    INTEGER :: i, j, num_cmd_args, arg_len, eq_pos, opt_name_end, &
         algo, ialign, size_hamocc, ierror, optarg_len
    LOGICAL :: next_is_optarg, contradiction
    INTEGER, PARAMETER :: max_opt_arg_len = 256
    CHARACTER(max_opt_arg_len) :: opt, optarg
    CHARACTER(kind=c_char, len=1) :: optarg_c(max_opt_arg_len+1)

    num_cmd_args = COMMAND_ARGUMENT_COUNT()
    i = 1
    DO WHILE (i <= num_cmd_args)
      CALL GET_COMMAND_ARGUMENT(i, opt, arg_len, status=ierror)
      IF (ierror /= 0) THEN
        IF (arg_len > LEN(opt)) THEN
          WRITE (0, '(a,i0,a)') 'error: command-line argument (number ', &
               i, ') too long!'
        ELSE
          WRITE (0, '(2(a, i0))') 'error handling command-line argument ', i, &
               ', status=', ierror
        END IF
        CALL posix_exit(1)
      END IF
      IF (opt(1:2) == '--') THEN
        eq_pos = INDEX(opt, '=')
        IF (eq_pos /= 0) THEN
          opt_name_end = eq_pos-1
        ELSE
          opt_name_end = arg_len
        END IF
        IF (eq_pos /= 0) THEN
          optarg = opt(eq_pos+1:)
          optarg_len = arg_len-eq_pos
        ELSE IF (i < num_cmd_args) THEN
          CALL GET_COMMAND_ARGUMENT(i+1, optarg, optarg_len)
        ELSE
          optarg_len = 0
        END IF
        SELECT CASE (opt(3:opt_name_end))
        CASE ("index-generation")
          CALL set_list_gen_methods(run_config%index_generation_sequence, &
               run_config%index_generation_method, optarg(:optarg_len))
          next_is_optarg = eq_pos == 0
        CASE ("index-src-file-pat")
          IF (optarg_len == 0) THEN
            WRITE (0, '(a)') "error: invalid source list file pattern!"
            FLUSH(0)
            CALL posix_exit(1)
          END IF
#ifdef HAVE_FC_ALLOCATABLE_CHARACTER
          ALLOCATE(CHARACTER(optarg_len) :: run_config%fname_pat_src)
#endif
          run_config%fname_pat_src(:) = optarg(1:optarg_len)
          next_is_optarg = eq_pos == 0
        CASE ("index-dst-file-pat")
          IF (optarg_len == 0) THEN
            WRITE (0, '(a)') "error: invalid destination list file pattern!"
            FLUSH(0)
            CALL posix_exit(1)
          END IF
#ifdef HAVE_FC_ALLOCATABLE_CHARACTER
          ALLOCATE(CHARACTER(optarg_len) :: run_config%fname_pat_dst)
#endif
          run_config%fname_pat_dst(:) = optarg(1:optarg_len)
          next_is_optarg = eq_pos == 0
        CASE ("num-columns")
          READ (optarg, *) run_config%num_columns
          next_is_optarg = eq_pos == 0
        CASE ("nlev", "num-levels")
          READ (optarg, *) run_config%nlev
          next_is_optarg = eq_pos == 0
        CASE ("sort-algorithm")
          DO j = 1, optarg_len
            optarg_c(j) = optarg(j:j)
          END DO
          optarg_c(optarg_len+1) = c_null_char
          algo = xt_sort_algo_id_by_name(optarg_c)
          IF (algo == -1) THEN
            WRITE (0, '(2a)') "error: invalid sort algorithm name: ", &
                 optarg
            FLUSH(0)
            CALL posix_exit(1)
          END IF
          CALL xt_config_set_sort_algorithm_by_id(config, algo)
          next_is_optarg = eq_pos == 0
        CASE ("enable-mem-saving", "disable-mem-saving")
          IF (eq_pos /= 0) THEN
            WRITE (0, '(4a)') &
                 "error: invalid argument to ", TRIM(opt), ": ", &
                 optarg
            FLUSH(0)
            CALL posix_exit(1)
          END IF
          next_is_optarg = .FALSE.
          CALL xt_config_set_mem_saving(config, MERGE(1, 0, arg_len == 19))
        CASE ("enable-custom-bucket-generator")
          CALL set_custom_bucket_gen(config)
          IF (eq_pos /= 0) THEN
            CALL set_list_gen_methods(run_config%bucket_generation_sequence, &
                 run_config%bucket_generation_method, optarg)
          END IF
          next_is_optarg = .FALSE.
        CASE ("stripe-alignment-mode", &
             "disable-stripe-alignment", "enable-stripe-alignment")
          IF (opt(3:8) == 'stripe') THEN
            IF (eq_pos == 0 .AND. i == num_cmd_args) THEN
              WRITE (0, '(2a)') "error: option --stripe-alignment-mode ", &
                   "requires an argument"
              FLUSH(0)
              CALL posix_exit(1)
            END IF
            READ (optarg, *) ialign
            next_is_optarg = eq_pos == 0
          ELSE
            IF (eq_pos /= 0) THEN
              WRITE (0, '(4a)') &
                   "error: invalid argument to ", TRIM(opt),": ", &
                   optarg
              FLUSH(0)
              CALL posix_exit(1)
            END IF
            ialign = MERGE(1, 0, opt(3:4) == "en")
            next_is_optarg = .FALSE.
          END IF
          CALL xt_config_set_xmap_stripe_align(config, ialign)
        CASE DEFAULT
          WRITE (0, '(2a)') 'unrecognized command line argument: ', &
               optarg(1:arg_len)
          FLUSH(0)
          CALL posix_exit(1)
          STOP
        END SELECT
        i = i + MERGE(2, 1, next_is_optarg)
      ELSE
        EXIT
      END IF
    END DO
    IF (i <= num_cmd_args) THEN
      IF (i == num_cmd_args) THEN
        READ(opt, *, iostat=ierror) size_hamocc
        IF (ierror /= 0) THEN
          WRITE (0, '(3a, i0)') 'error reading group size ', TRIM(opt), &
               ', iostat=', ierror
          FLUSH(0)
          CALL posix_exit(1)
        END IF
      ELSE
        WRITE (0, '(a)') 'too many command line arguments:', opt(1:arg_len)
        DO i = i+1, num_cmd_args
          CALL GET_COMMAND_ARGUMENT(i, opt, arg_len)
          WRITE (0, '(a)') opt(1:arg_len)
        END DO
        CALL posix_exit(1)
      END IF
    ELSE
      size_hamocc = run_config%comm_world_size / 2
    END IF
    run_config%is_src = run_config%comm_world_rank &
         < (run_config%comm_world_size - size_hamocc)
    run_config%local_size = MERGE(run_config%comm_world_size - size_hamocc,  &
         size_hamocc, run_config%is_src)
    IF (run_config%bucket_generation_sequence == -1) &
         run_config%bucket_generation_sequence &
         = run_config%index_generation_sequence
    IF (run_config%bucket_generation_method == -1) &
         run_config%bucket_generation_method &
         = run_config%index_generation_sequence

#ifdef HAVE_FC_ALLOCATABLE_CHARACTER
    contradiction = ALLOCATED(run_config%fname_pat_src) &
         .NEQV. ALLOCATED(run_config%fname_pat_dst)
#else
    contradiction = run_config%fname_pat_src == '' &
         .NEQV. run_config%fname_pat_dst == ''
#endif
    IF (contradiction) THEN
      WRITE (0, '(3a)') "error: ", &
           "option --index-src-file-pat must be used together", &
           " with --index-dst-file-pat"
      CALL posix_exit(1)
    END IF
  END SUBROUTINE parse_options

  FUNCTION read_idxlist(idxgen_method, name_pat, subst_mpi_rank)
    INTEGER, INTENT(in) :: idxgen_method, subst_mpi_rank
    CHARACTER(len=*), INTENT(in) :: name_pat
    CHARACTER(len=LEN(name_pat)+5) :: fname
    TYPE(xt_idxlist) :: read_idxlist
    INTEGER :: ofs, num_indices, num_idx_per_line, i, last, &
         num_lines
    INTEGER :: ierror
    INTEGER, PARAMETER :: input_unit = 10
    INTEGER(xt_int_kind), ALLOCATABLE :: indices(:)
    IF (subst_mpi_rank < 0) THEN
      read_idxlist = xt_idxempty_new()
      RETURN
    END IF
    ofs = INDEX(name_pat, "%{rank}")
    IF (ofs == 0) THEN
      fname = name_pat
    ELSE
      WRITE (fname, "(a,i0,a)") name_pat(1:ofs-1), &
           subst_mpi_rank, name_pat(ofs+7:)
    END IF
    OPEN(unit=input_unit, file=fname, action="read", status="old")
    SELECT CASE (idxgen_method)
    CASE (index_gen_idxvec)
      READ(input_unit, *) num_indices, num_idx_per_line
      ALLOCATE(indices(num_indices))
      num_lines = (num_indices+num_idx_per_line-1)/num_idx_per_line
      DO i = 1, num_lines
        ofs = (i - 1) * num_idx_per_line
        last = MERGE(ofs + num_idx_per_line, num_indices, i < num_lines)
        READ(input_unit, *) indices(ofs + 1:last)
      END DO
      read_idxlist = xt_idxvec_new(indices)
    CASE default
      WRITE (0, '(a)') 'unsupported index list method selected for reading!'
      CALL mpi_abort(mpi_comm_world, 1, ierror)
    END SELECT
    CLOSE(input_unit)
  END FUNCTION read_idxlist

END MODULE perf_xmap_create_init

PROGRAM perf_xmap_create
  USE mpi
  USE yaxt, xi => xt_int_kind
  USE iso_c_binding
  USE perf_xmap_create_config
  USE perf_xmap_create_init
  IMPLICIT NONE
  TYPE(xt_config) :: conf
  TYPE(xt_idxlist) :: idxlist_3d, idxlist_empty, &
       idxlist_src, idxlist_dst
  TYPE(xt_xmap) :: xmap, reo
  DOUBLE PRECISION :: tic, dt, min, max, sum

  LOGICAL :: use_synthetic_indices
  INTEGER :: ierror

  IF (HUGE(1_xi) >= 1240000) THEN
    run_config%num_columns = INT(1240, xi)
    run_config%num_columns = run_config%num_columns * 1000_xi
  ELSE
    run_config%num_columns = 12400_xi
  END IF
  run_config%nlev = 50_xi
  run_config%index_generation_sequence = INDEX_GEN_COL_MAJOR
  run_config%index_generation_method = INDEX_GEN_IDXVEC
  run_config%bucket_generation_sequence = -1
  run_config%bucket_generation_method = -1
#ifndef HAVE_FC_ALLOCATABLE_CHARACTER
  run_config%fname_pat_src = ''
  run_config%fname_pat_dst = ''
#endif

  CALL mpi_init(ierror)
  CALL xt_initialize(mpi_comm_world)
  CALL mpi_comm_rank(mpi_comm_world, run_config%comm_world_rank, ierror)
  CALL mpi_comm_size(mpi_comm_world, run_config%comm_world_size, ierror)

  conf = xt_config_new()

  CALL parse_options(run_config, conf)

#ifdef HAVE_FC_ALLOCATABLE_CHARACTER
  use_synthetic_indices = .NOT. ALLOCATED(run_config%fname_pat_src)
#else
  use_synthetic_indices = run_config%fname_pat_src == ''
#endif
  IF (use_synthetic_indices) THEN
    idxlist_3d = generate_3d_list(&
         MOD(run_config%comm_world_rank, run_config%local_size), &
         run_config%local_size, run_config%index_generation_method, &
         run_config%index_generation_sequence)
    idxlist_empty = xt_idxempty_new()
    idxlist_src = MERGE(idxlist_3d, idxlist_empty, run_config%is_src)
    idxlist_dst = MERGE(idxlist_empty, idxlist_3d, run_config%is_src)
  ELSE
    idxlist_src = read_idxlist(run_config%index_generation_method, &
         run_config%fname_pat_src, run_config%comm_world_rank)
    idxlist_dst = read_idxlist(run_config%index_generation_method, &
         run_config%fname_pat_dst, run_config%comm_world_rank)
  END IF

  tic = mpi_wtime()

  xmap = xt_xmap_dist_dir_custom_new(&
       idxlist_src, idxlist_dst, mpi_comm_world, conf)

  dt = mpi_wtime() - tic

  reo = xt_xmap_reorder_custom(xmap, xt_reorder_recv_up, conf)

  CALL xt_xmap_delete(reo)
  CALL xt_xmap_delete(xmap)
  CALL xt_idxlist_delete(idxlist_src)
  CALL xt_idxlist_delete(idxlist_dst)

  CALL mpi_reduce(dt, min, 1, mpi_double_precision, &
       mpi_min, 0, mpi_comm_world, ierror)
  CALL mpi_reduce(dt, max, 1, mpi_double_precision, &
       mpi_max, 0, mpi_comm_world, ierror)
  CALL mpi_reduce(dt, sum, 1, mpi_double_precision, &
       mpi_sum, 0, mpi_comm_world, ierror)

  IF (run_config%comm_world_rank == 0) THEN
    WRITE (0, "(3(a,' ',es12.6,' '))") "min", min, "max", max, &
         "avg", sum / REAL(run_config%comm_world_size, KIND(1.0d0))
  END IF

  CALL xt_config_delete(conf);
  IF (.NOT. xt_is_null(run_config%bg)) &
    CALL xt_xmdd_bucket_gen_delete(run_config%bg)

  CALL xt_finalize()
  CALL MPI_Finalize(ierror)

END PROGRAM perf_xmap_create
!
! Local Variables:
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! license-project-url: "https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/"
! license-default: "bsd"
! End:
!
