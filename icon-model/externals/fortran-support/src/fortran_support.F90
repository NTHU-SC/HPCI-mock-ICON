! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2025, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! Public interface for the fortran-support library
MODULE fortran_support

  USE mo_exception, ONLY: init_logger, debug_on, debug_off, &
    & set_msg_timestamp, enable_logging, message, warning, finish, &
    & print_value, message_to_own_unit, message_text
  USE mo_expression, ONLY: expression, parse_expression_string
  USE mo_fortran_tools, ONLY: assign_if_present, assign_if_present_allocatable, if_associated, &
    & t_ptr_1d_dp, t_ptr_1d_sp, t_ptr_1d_int, &
    & t_ptr_2d_dp, t_ptr_2d_sp, t_ptr_2d_int, &
    & t_ptr_3d_dp, t_ptr_3d_sp, t_ptr_3d_int, t_ptr_i2d3d, &
    & t_ptr_4d_dp, t_ptr_4d_sp, t_ptr_4d_int, &
    & t_ptr_2d3d, t_ptr_2d3d_vp, t_ptr_tracer, &
    & copy, init, swap, negative2zero, var_scale, &
    & var_add, init_zero_contiguous, init_contiguous, &
    & init_zero_contiguous_dp, init_zero_contiguous_sp, &
    & init_contiguous_dp, init_contiguous_sp, init_contiguous_i4, &
    & init_contiguous_l, minval_1d, minval_2d, resize_arr_c1d, DO_DEALLOCATE, &
    & DO_PTR_DEALLOCATE, insert_dimension, assert_acc_host_only, &
    & assert_acc_device_only, set_acc_host_or_device, set_acc_async_queue
#ifdef __SINGLE_PRECISION
  USE mo_fortran_tools, ONLY: &
    & t_ptr_1d_wp => t_ptr_1d_sp, &
    & t_ptr_2d_wp => t_ptr_2d_sp, &
    & t_ptr_3d_wp => t_ptr_3d_sp, &
    & t_ptr_4d_wp => t_ptr_4d_sp
#else
  USE mo_fortran_tools, ONLY: &
    & t_ptr_1d_wp => t_ptr_1d_dp, &
    & t_ptr_2d_wp => t_ptr_2d_dp, &
    & t_ptr_3d_wp => t_ptr_3d_dp, &
    & t_ptr_4d_wp => t_ptr_4d_dp
#endif

  USE mo_hash_table, ONLY: t_HashTable, hashTable_make, t_HashIterator
  USE mo_io_units, ONLY: filename_max, nerr, nlog, nnml, nstat, ngmt, nin, &
    & nout, nnml_output, find_next_free_unit
  USE mo_namelist, ONLY: position_nml, open_nml, close_nml, POSITIONED, &
    & MISSING, LENGTH_ERROR, READ_ERROR, open_nml_output, close_nml_output
  USE mo_octree, ONLY: octree_init, octree_finalize, octree_count_point, &
    & octree_query_point, t_range_octree, OCTREE_DEPTH
  USE mo_simple_dump, ONLY: dump2text
  USE mo_util_backtrace, ONLY: ftn_util_backtrace
  USE mo_util_file, ONLY: util_symlink, util_unlink, util_islink, util_rename, &
    & util_tmpnam, util_filesize, util_file_is_writable, createSymlink, &
    & get_filename, get_filename_noext, get_path
  USE mo_util_libc, ONLY: memset_f, memcmp_f, memcpy_f, strerror
  USE mo_util_nml, ONLY: util_annotate_nml
  USE mo_util_rusage, ONLY: add_rss_list, add_rss_usage, print_rss_usage, &
    & close_rss_lists
  USE mo_util_sort, ONLY: quicksort, insertion_sort
#ifdef __SX__
  USE mo_util_sort, ONLY: radixsort, radixsort_int
#endif
  USE mo_util_stride, ONLY: util_stride_1d, util_stride_2d, util_get_ptrdiff
  USE mo_util_string, ONLY: tolower, lowcase, toupper, separator, int2string, &
    & real2string, logical2string, split_string, string_contains_word, &
    & tocompact, str_replace, t_keyword_list, associate_keyword, &
    & with_keywords, remove_duplicates, difference, add_to_list, one_of, &
    & insert_group, delete_keyword_list, sort_and_compress_list, tohex, &
    & remove_whitespace, pretty_print_string_list, find_trailing_number, &
    & toCharArray, toCharacter, c2f_char, charArray_dup, charArray_equal, &
    & charArray_toLower, normal, bold, fg_black, fg_red, fg_green, fg_yellow, &
    & fg_blue, fg_magenta, fg_cyan, fg_white, fg_default, bg_black, bg_red, &
    & bg_green, bg_yellow, bg_blue, bg_magenta, bg_cyan, bg_white, bg_default, &
    & new_list
  USE mo_util_string_parse, ONLY: util_do_parse_intlist
  USE mo_util_system, ONLY: util_exit, util_abort, util_system
#ifdef __XT3__
  USE mo_util_system, ONLY: util_base_iobuf
#endif
  USE mo_util_table, ONLY: initialize_table, finalize_table, add_table_column, &
    & set_table_entry, print_table, t_table, t_column
  USE mo_util_texthash, ONLY: text_hash, text_hash_c, text_isEqual, sel_char
#if defined(__PGI) || defined(__FLANG)
  USE mo_util_texthash, ONLY: t_char_workaround
#endif
  USE mo_util_timer, ONLY: util_cputime, util_walltime, util_gettimeofday, &
    & util_init_real_time, util_get_real_time_size, util_read_real_time, &
    & util_diff_real_time

  ! From mo_exception
  PUBLIC :: init_logger, debug_on, debug_off, set_msg_timestamp, &
    & enable_logging, message, warning, finish, print_value, &
    & message_to_own_unit, message_text

  ! From mo_expression
  PUBLIC :: expression, parse_expression_string

  ! From mo_fortran_tools
  PUBLIC :: assign_if_present, assign_if_present_allocatable, if_associated, &
    & t_ptr_1d_wp, t_ptr_1d_dp, t_ptr_1d_sp, t_ptr_1d_int, &
    & t_ptr_2d_wp, t_ptr_2d_dp, t_ptr_2d_sp, t_ptr_2d_int, &
    & t_ptr_3d_wp, t_ptr_3d_dp, t_ptr_3d_sp, t_ptr_3d_int, t_ptr_i2d3d, &
    & t_ptr_4d_wp, t_ptr_4d_dp, t_ptr_4d_sp, t_ptr_4d_int, &
    & t_ptr_2d3d, t_ptr_2d3d_vp, t_ptr_tracer, &
    & copy, init, swap, negative2zero, var_scale, &
    & var_add, init_zero_contiguous, init_contiguous, &
    & init_zero_contiguous_dp, init_zero_contiguous_sp, &
    & init_contiguous_dp, init_contiguous_sp, init_contiguous_i4, &
    & init_contiguous_l, minval_1d, minval_2d, resize_arr_c1d, DO_DEALLOCATE, &
    & DO_PTR_DEALLOCATE, insert_dimension, assert_acc_host_only, &
    & assert_acc_device_only, set_acc_host_or_device, set_acc_async_queue

  ! From mo_hash_table
  PUBLIC :: t_HashTable, hashTable_make, t_HashIterator

  ! From mo_io_units
  PUBLIC :: filename_max, nerr, nlog, nnml, nstat, ngmt, nin, nout, &
    & nnml_output, find_next_free_unit

  ! From mo_namelist
  PUBLIC :: position_nml, open_nml, close_nml, POSITIONED, MISSING, &
    & LENGTH_ERROR, READ_ERROR, open_nml_output, close_nml_output

  ! From mo_octree
  PUBLIC :: octree_init, octree_finalize, octree_count_point, &
    & octree_query_point, t_range_octree, OCTREE_DEPTH

  ! From mo_simple_dump
  PUBLIC :: dump2text

  ! From mo_util_backtrace
  PUBLIC :: ftn_util_backtrace

  ! From mo_util_file
  PUBLIC :: util_symlink, util_unlink, util_islink, util_rename, util_tmpnam, &
    & util_filesize, util_file_is_writable, createSymlink, get_filename, &
    & get_filename_noext, get_path

  ! From mo_util_libc
  PUBLIC :: memset_f, memcmp_f, memcpy_f, strerror

  ! From mo_util_nml
  PUBLIC :: util_annotate_nml

  ! From mo_util_rusage
  PUBLIC :: add_rss_list, add_rss_usage, print_rss_usage, close_rss_lists

  ! From mo_util_sort
  PUBLIC :: quicksort, insertion_sort
#ifdef __SX__
  PUBLIC :: radixsort, radixsort_int
#endif

  ! From mo_util_stride
  PUBLIC :: util_stride_1d, util_stride_2d, util_get_ptrdiff

  ! From mo_util_string
  PUBLIC :: tolower, lowcase, toupper, separator, int2string, real2string, &
    & logical2string, split_string, string_contains_word, tocompact, &
    & str_replace, t_keyword_list, associate_keyword, with_keywords, &
    & remove_duplicates, difference, add_to_list, one_of, insert_group, &
    & delete_keyword_list, sort_and_compress_list, tohex, remove_whitespace, &
    & pretty_print_string_list, find_trailing_number, toCharArray, &
    & toCharacter, c2f_char, charArray_dup, charArray_equal, &
    & charArray_toLower, normal, bold, fg_black, fg_red, fg_green, fg_yellow, &
    & fg_blue, fg_magenta, fg_cyan, fg_white, fg_default, bg_black, bg_red, &
    & bg_green, bg_yellow, bg_blue, bg_magenta, bg_cyan, bg_white, bg_default, &
    & new_list

  ! From mo_util_string_parse
  PUBLIC :: util_do_parse_intlist

  ! From mo_util_system
  PUBLIC :: util_exit, util_abort, util_system
#ifdef __XT3__
  PUBLIC :: util_base_iobuf
#endif

  ! From mo_util_table
  PUBLIC :: initialize_table, finalize_table, add_table_column, &
    & set_table_entry, print_table, t_table, t_column

  ! From mo_util_texthash
  PUBLIC :: text_hash, text_hash_c, text_isEqual, sel_char
#if defined(__PGI) || defined(__FLANG)
  PUBLIC :: t_char_workaround
#endif

  ! From mo_util_timer
  PUBLIC :: util_cputime, util_walltime, util_gettimeofday, &
    & util_init_real_time, util_get_real_time_size, util_read_real_time, &
    & util_diff_real_time

END MODULE
