MODULE global_index_lookup

  IMPLICIT NONE

  ! from mo_decomposition_tools.f90
  TYPE t_glb2loc_index_lookup
    ! sorted list of global indices.
    ! supposed to be used in a binary_search
    INTEGER, ALLOCATABLE :: glb_index(:)
    ! list of local indices. Supposed to lookup after the global index has been
    ! found in glb_index
    INTEGER, ALLOCATABLE :: glb_index_to_loc(:)
  END TYPE  t_glb2loc_index_lookup

  PUBLIC :: setup_glb2loc, glb2loc_lookup, t_glb2loc_index_lookup

CONTAINS

  SUBROUTINE setup_glb2loc(glb, glb2loc_index_lookup)
    INTEGER, INTENT(IN)                         :: glb(:)
    TYPE(t_glb2loc_index_lookup), INTENT(INOUT) :: glb2loc_index_lookup

    INTEGER :: local_size, i

    local_size = SIZE(glb)

    ! compute glb2loc lookup by sorting the global indices and storing the perm.
    ALLOCATE(glb2loc_index_lookup%glb_index(local_size))
    glb2loc_index_lookup%glb_index = glb
    ALLOCATE(glb2loc_index_lookup%glb_index_to_loc(local_size))
    glb2loc_index_lookup%glb_index_to_loc = (/(i, i=1,local_size)/)
    CALL quicksort_permutation_int(glb2loc_index_lookup%glb_index, &
                                   glb2loc_index_lookup%glb_index_to_loc)
  END SUBROUTINE setup_glb2loc

  ! taken from ICONs `mo_util_sort.f90`:
  SUBROUTINE swap_permutation_int(a, i,j, permutation)
    !> array for in-situ sorting
    INTEGER,  INTENT(INOUT)           :: a(:)
    !> indices to be exchanged
    INTEGER,  INTENT(IN)              :: i,j
    !> (optional) permutation of indices
    INTEGER,  INTENT(INOUT)           :: permutation(:)
    ! local variables
    INTEGER :: t, t_p

    t    = a(i)
    a(i) = a(j)
    a(j) = t
    t_p            = permutation(i)
    permutation(i) = permutation(j)
    permutation(j) = t_p
  END SUBROUTINE swap_permutation_int

  ! --------------------------------------------------------------------
  !> Simple recursive implementation of Hoare's QuickSort algorithm
  !  for a 1D array of INTEGER values.
  !
  !  Ordering after the sorting process: smallest...largest.
  !
  RECURSIVE SUBROUTINE quicksort_permutation_int(a, permutation, l_in, r_in)
    INTEGER,  INTENT(INOUT)           :: a(:)           !< array for in-situ sorting
    INTEGER,  INTENT(INOUT)           :: permutation(:) !< (optional) permutation of indices
    INTEGER,  INTENT(IN),    OPTIONAL :: l_in,r_in      !< left, right partition indices
    ! local variables
    INTEGER :: i,j,l,r,t_p,t,v,m

    IF (PRESENT(l_in)) THEN
      l = l_in
    ELSE
      l = 1
    END IF
    IF (PRESENT(r_in)) THEN
      r = r_in
    ELSE
      r = SIZE(a,1)
    END IF
    IF (r>l) THEN
      i = l-1
      j = r

      ! median-of-three selection of partitioning element
      IF ((r-l) > 3) THEN
        m = (l+r)/2
        IF (a(l)>a(m))  CALL swap_permutation_int(a, l,m, permutation)
        IF (a(l)>a(r)) THEN
          CALL swap_permutation_int(a, l,r, permutation)
        ELSE IF (a(r)>a(m)) THEN
          CALL swap_permutation_int(a, r,m, permutation)
        END IF
      END IF

      v = a(r)
      LOOP : DO
        CNTLOOP1 : DO
          i = i+1
          IF (a(i) >= v) EXIT CNTLOOP1
        END DO CNTLOOP1
        CNTLOOP2 : DO
          j = j-1
          IF ((a(j) <= v) .OR. (j==1)) EXIT CNTLOOP2
        END DO CNTLOOP2
        t    = a(i)
        a(i) = a(j)
        a(j) = t

        t_p            = permutation(i)
        permutation(i) = permutation(j)
        permutation(j) = t_p

        IF (j <= i) EXIT LOOP
      END DO LOOP
      a(j) = a(i)
      a(i) = a(r)
      a(r) = t

      permutation(j) = permutation(i)
      permutation(i) = permutation(r)
      permutation(r) = t_p
      CALL quicksort_permutation_int(a,permutation,l,i-1)
      CALL quicksort_permutation_int(a,permutation,i+1,r)
    END IF
  END SUBROUTINE quicksort_permutation_int

  ! taken from ICONs mo_decomposition_tools.f90
  PURE FUNCTION binary_search(array, key)

    INTEGER, INTENT(IN) :: array(:), key
    INTEGER :: binary_search

    INTEGER :: lb, ub, middle

    !$ACC ROUTINE SEQ

    lb = 1
    ub = SIZE(array)
    middle = ub / 2

    IF (ub == 0) THEN
      binary_search = 0
      RETURN
    END IF

    DO WHILE (ub >= lb)

      middle = (ub + lb) / 2
      IF (array(middle) < key) THEN
        lb = middle + 1
      ELSE IF (array(middle) > key) THEN
        ub = middle - 1
      ELSE
        EXIT
      END IF
    END DO

    IF (array(middle) == key) THEN
      binary_search = middle
    ELSE IF (array(middle) > key) THEN
      binary_search = -middle + 1
    ELSE
      binary_search = -middle
    END IF
  END FUNCTION binary_search

  INTEGER FUNCTION glb2loc_lookup(this, glb) RESULT(loc)
    CLASS(t_glb2loc_index_lookup), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: glb
    INTEGER :: pos
    pos = binary_search(this%glb_index, glb)
    loc = this%glb_index_to_loc(pos)
  END FUNCTION glb2loc_lookup

END MODULE global_index_lookup
