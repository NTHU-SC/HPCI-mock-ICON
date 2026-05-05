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

MODULE TEST_mo_util_texthash
  USE FORTUTF
  USE fortran_support, ONLY: text_hash_c, text_hash, text_isEqual, sel_char

CONTAINS
  SUBROUTINE TEST_text_hash_c_short
    CALL TAG_TEST("TEST_text_c_hash")
    CALL ASSERT_EQUAL(text_hash_c('Unittest'), 345061529)
  END SUBROUTINE

  SUBROUTINE TEST_text_hash_c_long
    CALL TAG_TEST("TEST_text_hash_c_long")
    CALL ASSERT_EQUAL(text_hash_c('UnittestFrameworkFortran'), 1314345762)
  END SUBROUTINE

  SUBROUTINE TEST_text_hash_short
    CLASS(*), POINTER :: ptr_word
    CHARACTER(8), TARGET :: word
    CALL TAG_TEST("TEST_text_hash")
    word = 'Unittest'
    ptr_word => word

    CALL ASSERT_EQUAL(text_hash(ptr_word), 345061529)
  END SUBROUTINE

  SUBROUTINE TEST_text_hash_long
    CLASS(*), POINTER :: ptr_word
    CHARACTER(24), TARGET :: word
    CALL TAG_TEST("TEST_text_hash_long")
    word = 'UnittestFrameworkFortran'
    ptr_word => word

    CALL ASSERT_EQUAL(text_hash(ptr_word), 1314345762)
  END SUBROUTINE

  SUBROUTINE TEST_text_isEqual
    CLASS(*), POINTER :: ptr_word_1, ptr_word_2
    CHARACTER(17), TARGET :: word_1, word_2, word_3
    word_1 = 'UnittestFramework'
    word_2 = 'UnittestFramework'
    word_3 = 'UnittestWorkframe'

    ptr_word_1 => word_1
    ptr_word_2 => word_2
    CALL TAG_TEST("TEST_text_isEqual_word1")
    CALL ASSERT_EQUAL(text_isEqual(ptr_word_1, ptr_word_2), .TRUE.)

    ptr_word_2 => word_3
    CALL TAG_TEST("TEST_text_isEqual_word2")
    CALL ASSERT_EQUAL(text_isEqual(ptr_word_1, ptr_word_2), .FALSE.)
  END SUBROUTINE

  SUBROUTINE TEST_sel_char
    CLASS(*), POINTER :: ptr_word_1
    CHARACTER(17), TARGET :: word_1, word_2
    word_1 = 'UnittestFramework'
    word_2 = 'UnittestWorkframe'

    CALL TAG_TEST("TEST_sel_char_word1")
    CALL ASSERT_EQUAL(word_1 == word_2, .FALSE.)

    ptr_word_1 => word_1
    word_2 = sel_char(ptr_word_1, 'test_mo_util_text_hash', 'Internal Error')
    CALL TAG_TEST("TEST_sel_char_word2")
    CALL ASSERT_EQUAL(word_1 == word_2, .TRUE.)
  END SUBROUTINE

END MODULE TEST_mo_util_texthash
