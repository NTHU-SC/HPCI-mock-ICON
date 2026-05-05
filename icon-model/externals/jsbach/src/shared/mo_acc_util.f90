!> A set of subroutines for comparing the values of matrix elements (1D-5D)
!> on CPU and GPU using OpenACC features
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>

#define _ICON_STYLE_

#ifdef _ICON_STYLE_
#define _REAL_TYPE_ wp
#define _ZERO_DP_ 0.0_wp
#else
#define _REAL_TYPE_ 8
#define _ZERO_DP_ 0.0
#endif

#ifdef _OPENACC

MODULE mo_acc_util

#ifdef _ICON_STYLE_
  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish
#endif

  USE OPENACC

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: accCompareRealMatrixHostDevice

  INTERFACE accCompareRealMatrixHostDevice
    MODULE PROCEDURE accCompareRealMatrixHostDevice1D
    MODULE PROCEDURE accCompareRealMatrixHostDevice2D
    MODULE PROCEDURE accCompareRealMatrixHostDevice3D
    MODULE PROCEDURE accCompareRealMatrixHostDevice4D
    MODULE PROCEDURE accCompareRealMatrixHostDevice5D
  END INTERFACE accCompareRealMatrixHostDevice

  CHARACTER(len=*), PARAMETER :: modname = 'mo_acc_util'

CONTAINS

  ! ====================================================================================================== !
  !
  SUBROUTINE accCompareRealMatrixHostDevice1D(matrix, matrixName, deviation, detailedOutput)

    REAL(_REAL_TYPE_), DIMENSION(:), POINTER, INTENT(IN) :: matrix
    CHARACTER(len=*), INTENT(IN) :: matrixName
    REAL(_REAL_TYPE_), OPTIONAL, INTENT(OUT) :: deviation
    LOGICAL, OPTIONAL, INTENT(IN) :: detailedOutput

    CHARACTER(len=*), PARAMETER :: routine = modname//':accCompareRealMatrixHostDevice1D'

    LOGICAL :: detailedOutputFlag
    INTEGER :: iL, iU, matrixSize, iLP
    INTEGER :: notPresentCellNumber
    REAL(_REAL_TYPE_), ALLOCATABLE :: matrixGPU(:)
    REAL(_REAL_TYPE_) :: minMatrix(1:2), maxMatrix(1:2), devLocal, devGlobal

    detailedOutputFlag = .FALSE.
    IF (PRESENT(detailedOutput)) detailedOutputFlag = detailedOutput

    iL = LBOUND(matrix, 1); iU = UBOUND(matrix, 1)
    matrixSize = iU - iL + 1;

    IF (.NOT. acc_is_present(matrix(iL:iU))) THEN

      notPresentCellNumber = 0
      DO iLP = iL, iU
        IF (.NOT. acc_is_present(matrix(iL:iU))) THEN
          notPresentCellNumber = notPresentCellNumber + 1
        END IF
      END DO

      WRITE(0,*)
      WRITE(0,*) '>>>>> accCompareRealArrayHostDevice1D'
      WRITE(0,*) '>>>>> GPU DATA PARTIALLY PRESENT'
      WRITE(0,*) '>>>>> MATRIX  NAME: ', TRIM(matrixName)
      WRITE(0,'(1x, A21, I6, A, I6, A, I6, A)') '>>>>> MATRIX SHAPE: (', iL, ':', iU, ')'
      WRITE(0,'(1x, A20, I14)') '>>>>> MATRIX  SIZE: ', matrixSize
      WRITE(0,'(1x, A20, I14)') '>>>>>  GPU PRESENT: ', matrixSize - notPresentCellNumber
      WRITE(0,*)

#ifndef _ICON_STYLE_
      STOP 'STOP PROGRAM'
#else
      CALL finish(routine, 'STOP PROGRAM')
#endif

    END IF

    ALLOCATE(matrixGPU(iL:iU))
    !$ACC ENTER DATA CREATE(matrixGPU(iL:iU))

    !$ACC PARALLEL LOOP DEFAULT(PRESENT)
    DO iLP = iL, iU
        matrixGPU(iLP) = matrix(iLP)
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC UPDATE HOST(matrixGPU(iL:iU))

    !$ACC EXIT DATA DELETE(matrixGPU(iL:iU))

    minMatrix(1) =    matrix(iL); maxMatrix(1) =    matrix(iL)
    minMatrix(2) = matrixGPU(iL); maxMatrix(2) = matrixGPU(iL)

    devGlobal = _ZERO_DP_

    DO iLP = iL, iU

      IF (minMatrix(1) > matrix(iLP)) THEN
        minMatrix(1) = matrix(iLP)
      END IF
      IF (maxMatrix(1) < matrix(iLP)) THEN
        maxMatrix(1) = matrix(iLP)
      END IF

      IF (minMatrix(2) > matrixGPU(iLP)) THEN
        minMatrix(2) = matrixGPU(iLP)
      END IF
      IF (maxMatrix(2) < matrixGPU(iLP)) THEN
        maxMatrix(2) = matrixGPU(iLP)
      END IF

      devLocal = ABS(matrix(iLP) - matrixGPU(iLP))
      IF (devGlobal < devLocal) THEN
        devGlobal = devLocal
      END IF

    END DO

    IF (PRESENT(deviation)) THEN
      deviation = devGlobal
    END IF

    DEALLOCATE(matrixGPU)

    IF (detailedOutputFlag ) THEN

      WRITE(0,*) '>>>>> accCompareRealArrayHostDevice1D'
      WRITE(0,*) '>>>>> MATRIX  NAME: ', TRIM(matrixName)
      WRITE(0,'(1x, A21, I6, A, I6, A, I6, A)') '>>>>> MATRIX SHAPE: (', iL, ':', iU, ')'
      WRITE(0,'(1x, A25, E16.6, A, E16.6, A, E16.6, A, E16.6, A, E16.6)') '>>>>>      MIN/MAX: CPU [', &
        & minMatrix(1), '; ', maxMatrix(1), ']  GPU [', minMatrix(2), '; ', maxMatrix(2), '] DEVIATION ', devGlobal

    ELSE

      WRITE(0,'(1x, A6, A, A7, E16.6, A, E16.6, A, E16.6, A, E16.6, A, E16.6)') '>>>>> ', TRIM(matrixName), ': CPU [', &
        & minMatrix(1), '; ', maxMatrix(1), ']  GPU [', minMatrix(2), '; ', maxMatrix(2), '] DEVIATION ', devGlobal

    END IF

  END SUBROUTINE accCompareRealMatrixHostDevice1D

  ! ====================================================================================================== !
  !
  SUBROUTINE accCompareRealMatrixHostDevice2D(matrix, matrixName, deviation, detailedOutput)

    REAL(_REAL_TYPE_), DIMENSION(:,:), POINTER, INTENT(IN) :: matrix
    CHARACTER(len=*), INTENT(IN) :: matrixName
    REAL(_REAL_TYPE_), OPTIONAL, INTENT(OUT) :: deviation
    LOGICAL, OPTIONAL, INTENT(IN) :: detailedOutput

    CHARACTER(len=*), PARAMETER :: routine = modname//':accCompareRealMatrixHostDevice2D'

    LOGICAL :: detailedOutputFlag
    INTEGER :: iL1, iL2, iU1, iU2, matrixSize, iLP1, iLP2
    INTEGER :: notPresentCellNumber
    REAL(_REAL_TYPE_), ALLOCATABLE :: matrixGPU(:,:)
    REAL(_REAL_TYPE_) :: minMatrix(1:2), maxMatrix(1:2), devLocal, devGlobal

    detailedOutputFlag = .FALSE.
    IF (PRESENT(detailedOutput)) detailedOutputFlag = detailedOutput

    iL1 = LBOUND(matrix, 1); iU1 = UBOUND(matrix, 1)
    iL2 = LBOUND(matrix, 2); iU2 = UBOUND(matrix, 2)
    matrixSize = (iU1 - iL1 + 1) * (iU2 - iL2 + 1);

    IF (.NOT. acc_is_present(matrix(iL1:iU1, iL2:iU2))) THEN

      notPresentCellNumber = 0
      DO iLP2 = iL2, iU2
        DO iLP1 = iL1, iU1
          IF (.NOT. acc_is_present(matrix(iLP1:iLP1, iLP2:iLP2))) THEN
            notPresentCellNumber = notPresentCellNumber + 1
          END IF
        END DO
      END DO

      WRITE(0,*)
      WRITE(0,*) '>>>>> accCompareRealArrayHostDevice2D'
      WRITE(0,*) '>>>>> GPU DATA PARTIALLY PRESENT'
      WRITE(0,*) '>>>>> MATRIX  NAME: ', TRIM(matrixName)
      WRITE(0,'(1x, A21, I6, A, I6, A, I6, A, I6, A)') '>>>>> MATRIX SHAPE: (', iL1, ':', iU1, ') (', iL2, ':', iU2, ')'
      WRITE(0,'(1x, A20, I14)') '>>>>> MATRIX  SIZE: ', matrixSize
      WRITE(0,'(1x, A20, I14)') '>>>>>  GPU PRESENT: ', matrixSize - notPresentCellNumber
      WRITE(0,*)

#ifndef _ICON_STYLE_
      STOP 'STOP PROGRAM'
#else
      CALL finish(routine, 'STOP PROGRAM')
#endif

    END IF

    ALLOCATE(matrixGPU(iL1:iU1, iL2:iU2))
    !$ACC ENTER DATA CREATE(matrixGPU(iL1:iU1, iL2:iU2))

    !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(PRESENT)
    DO iLP2 = iL2, iU2
      DO iLP1 = iL1, iU1
        matrixGPU(iLP1, iLP2) = matrix(iLP1, iLP2)
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC UPDATE HOST(matrixGPU(iL1:iU1, iL2:iU2))

    !$ACC EXIT DATA DELETE(matrixGPU(iL1:iU1, iL2:iU2))

    minMatrix(1) =    matrix(iL1, iL2); maxMatrix(1) =    matrix(iL1, iL2)
    minMatrix(2) = matrixGPU(iL1, iL2); maxMatrix(2) = matrixGPU(iL1, iL2)

    devGlobal = _ZERO_DP_

    DO iLP2 = iL2, iU2
      DO iLP1 = iL1, iU1

        IF (minMatrix(1) > matrix(iLP1, iLP2)) THEN
          minMatrix(1) = matrix(iLP1, iLP2)
        END IF
        IF (maxMatrix(1) < matrix(iLP1, iLP2)) THEN
          maxMatrix(1) = matrix(iLP1, iLP2)
        END IF

        IF (minMatrix(2) > matrixGPU(iLP1, iLP2)) THEN
          minMatrix(2) = matrixGPU(iLP1, iLP2)
        END IF
        IF (maxMatrix(2) < matrixGPU(iLP1, iLP2)) THEN
          maxMatrix(2) = matrixGPU(iLP1, iLP2)
        END IF

        devLocal = ABS(matrix(iLP1, iLP2) - matrixGPU(iLP1, iLP2))
        IF (devGlobal < devLocal) THEN
          devGlobal = devLocal
        END IF

      END DO
    END DO

    IF (PRESENT(deviation)) THEN
      deviation = devGlobal
    END IF

    DEALLOCATE(matrixGPU)

    IF (detailedOutputFlag) THEN

      WRITE(0,*) '>>>>> accCompareRealArrayHostDevice2D'
      WRITE(0,*) '>>>>> MATRIX  NAME: ', TRIM(matrixName)
      WRITE(0,'(1x, A21, I6, A, I6, A, I6, A, I6, A)') '>>>>> MATRIX SHAPE: (', iL1, ':', iU1, ') (', iL2, ':', iU2, ')'
      WRITE(0,'(1x, A25, E16.6, A, E16.6, A, E16.6, A, E16.6, A, E16.6)') '>>>>>      MIN/MAX: CPU [', &
        & minMatrix(1), '; ', maxMatrix(1), ']  GPU [', minMatrix(2), '; ', maxMatrix(2), '] DEVIATION ', devGlobal

    ELSE

      WRITE(0,'(1x, A6, A, A7, E16.6, A, E16.6, A, E16.6, A, E16.6, A, E16.6)') '>>>>> ', TRIM(matrixName), ': CPU [', &
        & minMatrix(1), '; ', maxMatrix(1), ']  GPU [', minMatrix(2), '; ', maxMatrix(2), '] DEVIATION ', devGlobal

    END IF

  END SUBROUTINE accCompareRealMatrixHostDevice2D

  ! ====================================================================================================== !
  !
  SUBROUTINE accCompareRealMatrixHostDevice3D(matrix, matrixName, deviation, detailedOutput)

    REAL(_REAL_TYPE_), DIMENSION(:,:,:), POINTER, INTENT(IN) :: matrix
    CHARACTER(len=*), INTENT(IN) :: matrixName
    REAL(_REAL_TYPE_), OPTIONAL, INTENT(OUT) :: deviation
    LOGICAL, OPTIONAL, INTENT(IN) :: detailedOutput

    CHARACTER(len=*), PARAMETER :: routine = modname//':accCompareRealMatrixHostDevice3D'

    LOGICAL :: detailedOutputFlag
    INTEGER :: iL1, iL2, iL3, iU1, iU2, iU3, matrixSize, iLP1, iLP2, iLP3
    INTEGER :: notPresentCellNumber
    REAL(_REAL_TYPE_), ALLOCATABLE :: matrixGPU(:,:,:)
    REAL(_REAL_TYPE_) :: minMatrix(1:2), maxMatrix(1:2), devLocal, devGlobal

    detailedOutputFlag = .FALSE.
    IF (PRESENT(detailedOutput)) detailedOutputFlag = detailedOutput

    iL1 = LBOUND(matrix, 1); iU1 = UBOUND(matrix, 1)
    iL2 = LBOUND(matrix, 2); iU2 = UBOUND(matrix, 2)
    iL3 = LBOUND(matrix, 3); iU3 = UBOUND(matrix, 3)
    matrixSize = (iU1 - iL1 + 1) * (iU2 - iL2 + 1) * (iU3 - iL3 + 1);

    IF (.NOT. acc_is_present(matrix(iL1:iU1, iL2:iU2, iL3:iU3))) THEN

      notPresentCellNumber = 0
      DO iLP3 = iL3, iU3
        DO iLP2 = iL2, iU2
          DO iLP1 = iL1, iU1
            IF (.NOT. acc_is_present(matrix(iLP1:iLP1, iLP2:iLP2, iLP3:iLP3))) THEN
              notPresentCellNumber = notPresentCellNumber + 1
            END IF
          END DO
        END DO
      END DO

      WRITE(0,*)
      WRITE(0,*) '>>>>> accCompareRealArrayHostDevice3D'
      WRITE(0,*) '>>>>> GPU DATA PARTIALLY PRESENT'
      WRITE(0,*) '>>>>> MATRIX  NAME: ', TRIM(matrixName)
      WRITE(0,'(1x, A21, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A)') '>>>>> MATRIX SHAPE: (', iL1, ':', iU1, ') (', &
        & iL2, ':', iU2, ') (', iL3, ':', iU3, ')'
      WRITE(0,'(1x, A20, I14)') '>>>>> MATRIX  SIZE: ', matrixSize
      WRITE(0,'(1x, A20, I14)') '>>>>>  GPU PRESENT: ', matrixSize - notPresentCellNumber
      WRITE(0,*)

#ifndef _ICON_STYLE_
      STOP 'STOP PROGRAM'
#else
      CALL finish(routine, 'STOP PROGRAM')
#endif

    END IF

    ALLOCATE(matrixGPU(iL1:iU1, iL2:iU2, iL3:iU3))

    !$ACC ENTER DATA CREATE(matrixGPU(iL1:iU1, iL2:iU2, iL3:iU3))
    !$ACC PARALLEL LOOP COLLAPSE(3) DEFAULT(PRESENT)
    DO iLP3 = iL3, iU3
      DO iLP2 = iL2, iU2
        DO iLP1 = iL1, iU1
          matrixGPU(iLP1, iLP2, iLP3) = matrix(iLP1, iLP2, iLP3)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC UPDATE HOST(matrixGPU(iL1:iU1, iL2:iU2, iL3:iU3))
    !$ACC EXIT DATA DELETE(matrixGPU(iL1:iU1, iL2:iU2, iL3:iU3))

    minMatrix(1) =    matrix(iL1, iL2, iL3); maxMatrix(1) =    matrix(iL1, iL2, iL3)
    minMatrix(2) = matrixGPU(iL1, iL2, iL3); maxMatrix(2) = matrixGPU(iL1, iL2, iL3)

    devGlobal = _ZERO_DP_

    DO iLP3 = iL3, iU3
      DO iLP2 = iL2, iU2
        DO iLP1 = iL1, iU1

          IF (minMatrix(1) > matrix(iLP1, iLP2, iLP3)) THEN
            minMatrix(1) = matrix(iLP1, iLP2, iLP3)
          END IF
          IF (maxMatrix(1) < matrix(iLP1, iLP2, iLP3)) THEN
            maxMatrix(1) = matrix(iLP1, iLP2, iLP3)
          END IF

          IF (minMatrix(2) > matrixGPU(iLP1, iLP2, iLP3)) THEN
            minMatrix(2) = matrixGPU(iLP1, iLP2, iLP3)
          END IF
          IF (maxMatrix(2) < matrixGPU(iLP1, iLP2, iLP3)) THEN
            maxMatrix(2) = matrixGPU(iLP1, iLP2, iLP3)
          END IF

          devLocal = ABS(matrix(iLP1, iLP2, iLP3) - matrixGPU(iLP1, iLP2, iLP3))
          IF (devGlobal < devLocal) THEN
            devGlobal = devLocal
          END IF

        END DO
      END DO
    END DO

    IF (PRESENT(deviation)) THEN
      deviation = devGlobal
    END IF

    DEALLOCATE(matrixGPU)

    IF (detailedOutputFlag) THEN

      WRITE(0,*) '>>>>> accCompareRealArrayHostDevice3D'
      WRITE(0,*) '>>>>> MATRIX  NAME: ', TRIM(matrixName)
      WRITE(0,'(1x, A21, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A)') '>>>>> MATRIX SHAPE: (', iL1, ':', iU1, ') (', &
        & iL2, ':', iU2, ') (', iL3, ':', iU3, ')'
      WRITE(0,'(1x, A25, E16.6, A, E16.6, A, E16.6, A, E16.6, A, E16.6)') '>>>>>      MIN/MAX: CPU [', &
        & minMatrix(1), '; ', maxMatrix(1), ']  GPU [', minMatrix(2), '; ', maxMatrix(2), '] DEVIATION ', devGlobal

    ELSE

      WRITE(0,'(1x, A6, A, A7, E16.6, A, E16.6, A, E16.6, A, E16.6, A, E16.6)') '>>>>> ', TRIM(matrixName), ': CPU [', &
        & minMatrix(1), '; ', maxMatrix(1), ']  GPU [', minMatrix(2), '; ', maxMatrix(2), '] DEVIATION ', devGlobal

    END IF

  END SUBROUTINE accCompareRealMatrixHostDevice3D

  ! ====================================================================================================== !
  !
  SUBROUTINE accCompareRealMatrixHostDevice4D(matrix, matrixName, deviation, detailedOutput)

    REAL(_REAL_TYPE_), DIMENSION(:,:,:,:), POINTER, INTENT(IN) :: matrix
    CHARACTER(len=*), INTENT(IN) :: matrixName
    REAL(_REAL_TYPE_), OPTIONAL, INTENT(OUT) :: deviation
    LOGICAL, OPTIONAL, INTENT(IN) :: detailedOutput

    CHARACTER(len=*), PARAMETER :: routine = modname//':accCompareRealMatrixHostDevice4D'

    LOGICAL :: detailedOutputFlag
    INTEGER :: iL1, iL2, iL3, iL4, iU1, iU2, iU3, iU4, matrixSize, iLP1, iLP2, iLP3, iLP4
    INTEGER :: notPresentCellNumber
    REAL(_REAL_TYPE_), ALLOCATABLE :: matrixGPU(:,:,:,:)
    REAL(_REAL_TYPE_) :: minMatrix(1:2), maxMatrix(1:2), devLocal, devGlobal

    detailedOutputFlag = .FALSE.
    IF (PRESENT(detailedOutput)) detailedOutputFlag = detailedOutput

    iL1 = LBOUND(matrix, 1); iU1 = UBOUND(matrix, 1)
    iL2 = LBOUND(matrix, 2); iU2 = UBOUND(matrix, 2)
    iL3 = LBOUND(matrix, 3); iU3 = UBOUND(matrix, 3)
    iL4 = LBOUND(matrix, 4); iU4 = UBOUND(matrix, 4)
    matrixSize = (iU1 - iL1 + 1) * (iU2 - iL2 + 1) * (iU3 - iL3 + 1) * (iU4 - iL4 + 1);

    IF (.NOT. acc_is_present(matrix(iL1:iU1, iL2:iU2, iL3:iU3, iL4:iU4))) THEN

      notPresentCellNumber = 0
      DO iLP4 = iL4, iU4
        DO iLP3 = iL3, iU3
          DO iLP2 = iL2, iU2
            DO iLP1 = iL1, iU1
              IF (.NOT. acc_is_present(matrix(iLP1:iLP1, iLP2:iLP2, iLP3:iLP3, iLP4:iLP4))) THEN
                notPresentCellNumber = notPresentCellNumber + 1
              END IF
            END DO
          END DO
        END DO
      END DO

      WRITE(0,*)
      WRITE(0,*) '>>>>> accCompareRealArrayHostDevice4D'
      WRITE(0,*) '>>>>> GPU DATA PARTIALLY PRESENT'
      WRITE(0,*) '>>>>> MATRIX  NAME: ', TRIM(matrixName)
      WRITE(0,'(1x, A21, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A)') '>>>>> MATRIX SHAPE: (', &
        & iL1, ':', iU1, ') (', iL2, ':', iU2, ') (', iL3, ':', iU3, ') (', iL4, ':', iU4, ')'
      WRITE(0,'(1x, A20, I14)') '>>>>> MATRIX  SIZE: ', matrixSize
      WRITE(0,'(1x, A20, I14)') '>>>>>  GPU PRESENT: ', matrixSize - notPresentCellNumber
      WRITE(0,*)

#ifndef _ICON_STYLE_
      STOP 'STOP PROGRAM'
#else
      CALL finish(routine, 'STOP PROGRAM')
#endif

    END IF

    ALLOCATE(matrixGPU(iL1:iU1, iL2:iU2, iL3:iU3, iL4:iU4))
    !$ACC ENTER DATA CREATE(matrixGPU(iL1:iU1, iL2:iU2, iL3:iU3, iL4:iU4))

    !$ACC PARALLEL LOOP COLLAPSE(4) DEFAULT(PRESENT)
    DO iLP4 = iL4, iU4
      DO iLP3 = iL3, iU3
        DO iLP2 = iL2, iU2
          DO iLP1 = iL1, iU1
            matrixGPU(iLP1, iLP2, iLP3, iLP4) = matrix(iLP1, iLP2, iLP3, iLP4)
          END DO
        END DO
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC UPDATE HOST(matrixGPU(iL1:iU1, iL2:iU2, iL3:iU3, iL4:iU4))

    !$ACC EXIT DATA DELETE(matrixGPU(iL1:iU1, iL2:iU2, iL3:iU3, iL4:iU4))

    minMatrix(1) =    matrix(iL1, iL2, iL3, iL4); maxMatrix(1) =    matrix(iL1, iL2, iL3, iL4)
    minMatrix(2) = matrixGPU(iL1, iL2, iL3, iL4); maxMatrix(2) = matrixGPU(iL1, iL2, iL3, iL4)

    devGlobal = _ZERO_DP_

    DO iLP4 = iL4, iU4
      DO iLP3 = iL3, iU3
        DO iLP2 = iL2, iU2
          DO iLP1 = iL1, iU1

            IF (minMatrix(1) > matrix(iLP1, iLP2, iLP3, iLP4)) THEN
              minMatrix(1) = matrix(iLP1, iLP2, iLP3, iLP4)
            END IF
            IF (maxMatrix(1) < matrix(iLP1, iLP2, iLP3, iLP4)) THEN
              maxMatrix(1) = matrix(iLP1, iLP2, iLP3, iLP4)
            END IF

            IF (minMatrix(2) > matrixGPU(iLP1, iLP2, iLP3, iLP4)) THEN
              minMatrix(2) = matrixGPU(iLP1, iLP2, iLP3, iLP4)
            END IF
            IF (maxMatrix(2) < matrixGPU(iLP1, iLP2, iLP3, iLP4)) THEN
              maxMatrix(2) = matrixGPU(iLP1, iLP2, iLP3, iLP4)
            END IF

            devLocal = ABS(matrix(iLP1, iLP2, iLP3, iLP4) - matrixGPU(iLP1, iLP2, iLP3, iLP4))
            IF (devGlobal < devLocal) THEN
              devGlobal = devLocal
            END IF

          END DO
        END DO
      END DO
    END DO

    IF (PRESENT(deviation)) THEN
      deviation = devGlobal
    END IF

    DEALLOCATE(matrixGPU)

    IF (detailedOutputFlag) THEN

      WRITE(0,*) '>>>>> accCompareRealArrayHostDevice4D'
      WRITE(0,*) '>>>>> MATRIX  NAME: ', TRIM(matrixName)
      WRITE(0,'(1x, A21, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A)') '>>>>> MATRIX SHAPE: (', &
        & iL1, ':', iU1, ') (', iL2, ':', iU2, ') (', iL3, ':', iU3, ') (', iL4, ':', iU4, ')'
      WRITE(0,'(1x, A25, E16.6, A, E16.6, A, E16.6, A, E16.6, A, E16.6)') '>>>>>      MIN/MAX: CPU [', &
        & minMatrix(1), '; ', maxMatrix(1), ']  GPU [', minMatrix(2), '; ', maxMatrix(2), '] DEVIATION ', devGlobal

    ELSE

      WRITE(0,'(1x, A6, A, A7, E16.6, A, E16.6, A, E16.6, A, E16.6, A, E16.6)') '>>>>> ', TRIM(matrixName), ': CPU [', &
        & minMatrix(1), '; ', maxMatrix(1), ']  GPU [', minMatrix(2), '; ', maxMatrix(2), '] DEVIATION ', devGlobal

    END IF

  END SUBROUTINE accCompareRealMatrixHostDevice4D

  ! ====================================================================================================== !
  !
  SUBROUTINE accCompareRealMatrixHostDevice5D(matrix, matrixName, deviation, detailedOutput)

    REAL(_REAL_TYPE_), DIMENSION(:,:,:,:,:), POINTER, INTENT(IN) :: matrix
    CHARACTER(len=*), INTENT(IN) :: matrixName
    REAL(_REAL_TYPE_), OPTIONAL, INTENT(OUT) :: deviation
    LOGICAL, OPTIONAL, INTENT(IN) :: detailedOutput

    CHARACTER(len=*), PARAMETER :: routine = modname//':accCompareRealMatrixHostDevice5D'

    LOGICAL :: detailedOutputFlag
    INTEGER :: iL1, iL2, iL3, iL4, iL5, iU1, iU2, iU3, iU4, iU5, matrixSize, iLP1, iLP2, iLP3, iLP4, iLP5
    INTEGER :: notPresentCellNumber
    REAL(_REAL_TYPE_), ALLOCATABLE :: matrixGPU(:,:,:,:,:)
    REAL(_REAL_TYPE_) :: minMatrix(1:2), maxMatrix(1:2), devLocal, devGlobal

    detailedOutputFlag = .FALSE.
    IF (PRESENT(detailedOutput)) detailedOutputFlag = detailedOutput

    iL1 = LBOUND(matrix, 1); iU1 = UBOUND(matrix, 1)
    iL2 = LBOUND(matrix, 2); iU2 = UBOUND(matrix, 2)
    iL3 = LBOUND(matrix, 3); iU3 = UBOUND(matrix, 3)
    iL4 = LBOUND(matrix, 4); iU4 = UBOUND(matrix, 4)
    iL5 = LBOUND(matrix, 5); iU5 = UBOUND(matrix, 5)
    matrixSize = (iU1 - iL1 + 1) * (iU2 - iL2 + 1) * (iU3 - iL3 + 1) * (iU4 - iL4 + 1) * (iU5 - iL5 + 1);

    IF (.NOT. acc_is_present(matrix(iL1:iU1, iL2:iU2, iL3:iU3, iL4:iU4, iL5:iU5))) THEN

      notPresentCellNumber = 0
      DO iLP5 = iL5, iU5
        DO iLP4 = iL4, iU4
          DO iLP3 = iL3, iU3
            DO iLP2 = iL2, iU2
              DO iLP1 = iL1, iU1
                IF (.NOT. acc_is_present(matrix(iLP1:iLP1, iLP2:iLP2, iLP3:iLP3, iLP4:iLP4, iLP5:iLP5))) THEN
                  notPresentCellNumber = notPresentCellNumber + 1
                END IF
              END DO
            END DO
          END DO
        END DO
      END DO

      WRITE(0,*)
      WRITE(0,*) '>>>>> accCompareRealArrayHostDevice5D'
      WRITE(0,*) '>>>>> GPU DATA PARTIALLY PRESENT'
      WRITE(0,*) '>>>>> MATRIX  NAME: ', TRIM(matrixName)
      WRITE(0,'(1x, A21, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A)') &
        & '>>>>> MATRIX SHAPE: (', iL1, ':', iU1, ') (', iL2, ':', iU2, ') &
        & (', iL3, ':', iU3, ') (', iL4, ':', iU4, ') (', iL5, ':', iU5, ')'
      WRITE(0,'(1x, A20, I14)') '>>>>> MATRIX  SIZE: ', matrixSize
      WRITE(0,'(1x, A20, I14)') '>>>>>  GPU PRESENT: ', matrixSize - notPresentCellNumber
      WRITE(0,*)

#ifndef _ICON_STYLE_
      STOP 'STOP PROGRAM'
#else
      CALL finish(routine, 'STOP PROGRAM')
#endif

    END IF

    ALLOCATE(matrixGPU(iL1:iU1, iL2:iU2, iL3:iU3, iL4:iU4, iL5:iU5))
    !$ACC ENTER DATA CREATE(matrixGPU(iL1:iU1, iL2:iU2, iL3:iU3, iL4:iU4, iL5:iU5))

    !$ACC PARALLEL LOOP COLLAPSE(5) DEFAULT(PRESENT)
    DO iLP5 = iL5, iU5
      DO iLP4 = iL4, iU4
        DO iLP3 = iL3, iU3
          DO iLP2 = iL2, iU2
            DO iLP1 = iL1, iU1
              matrixGPU(iLP1, iLP2, iLP3, iLP4, iLP5) = matrix(iLP1, iLP2, iLP3, iLP4, iLP5)
            END DO
          END DO
        END DO
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC UPDATE HOST(matrixGPU(iL1:iU1, iL2:iU2, iL3:iU3, iL4:iU4, iL5:iU5))

    !$ACC EXIT DATA DELETE(matrixGPU(iL1:iU1, iL2:iU2, iL3:iU3, iL4:iU4, iL5:iU5))

    minMatrix(1) =    matrix(iL1, iL2, iL3, iL4, iL5); maxMatrix(1) =    matrix(iL1, iL2, iL3, iL4, iL5)
    minMatrix(2) = matrixGPU(iL1, iL2, iL3, iL4, iL5); maxMatrix(2) = matrixGPU(iL1, iL2, iL3, iL4, iL5)

    devGlobal = _ZERO_DP_

    DO iLP5 = iL5, iU5
      DO iLP4 = iL4, iU4
        DO iLP3 = iL3, iU3
          DO iLP2 = iL2, iU2
            DO iLP1 = iL1, iU1

              IF (minMatrix(1) > matrix(iLP1, iLP2, iLP3, iLP4, iLP5)) THEN
                minMatrix(1) = matrix(iLP1, iLP2, iLP3, iLP4, iLP5)
              END IF
              IF (maxMatrix(1) < matrix(iLP1, iLP2, iLP3, iLP4, iLP5)) THEN
                maxMatrix(1) = matrix(iLP1, iLP2, iLP3, iLP4, iLP5)
              END IF

              IF (minMatrix(2) > matrixGPU(iLP1, iLP2, iLP3, iLP4, iLP5)) THEN
                minMatrix(2) = matrixGPU(iLP1, iLP2, iLP3, iLP4, iLP5)
              END IF
              IF (maxMatrix(2) < matrixGPU(iLP1, iLP2, iLP3, iLP4, iLP5)) THEN
                maxMatrix(2) = matrixGPU(iLP1, iLP2, iLP3, iLP4, iLP5)
              END IF

              devLocal = ABS(matrix(iLP1, iLP2, iLP3, iLP4, iLP5) - matrixGPU(iLP1, iLP2, iLP3, iLP4, iLP5))
              IF (devGlobal < devLocal) THEN
                devGlobal = devLocal
              END IF

            END DO
          END DO
        END DO
      END DO
    END DO

    IF (PRESENT(deviation)) THEN
      deviation = devGlobal
    END IF

    DEALLOCATE(matrixGPU)

    IF (detailedOutputFlag) THEN

      WRITE(0,*) '>>>>> accCompareRealArrayHostDevice5D'
      WRITE(0,*) '>>>>> MATRIX  NAME: ', TRIM(matrixName)
      WRITE(0,'(1x, A21, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A)') &
        & '>>>>> MATRIX SHAPE: (', iL1, ':', iU1, ') (', iL2, ':', iU2, ') &
        & (', iL3, ':', iU3, ') (', iL4, ':', iU4, ') (', iL5, ':', iU5, ')'
      WRITE(0,'(1x, A25, E16.6, A, E16.6, A, E16.6, A, E16.6, A, E16.6)') '>>>>>      MIN/MAX: CPU [', &
        & minMatrix(1), '; ', maxMatrix(1), ']  GPU [', minMatrix(2), '; ', maxMatrix(2), '] DEVIATION ', devGlobal

    ELSE

      WRITE(0,'(1x, A6, A, A7, E16.6, A, E16.6, A, E16.6, A, E16.6, A, E16.6)') '>>>>> ', TRIM(matrixName), ': CPU [', &
        & minMatrix(1), '; ', maxMatrix(1), ']  GPU [', minMatrix(2), '; ', maxMatrix(2), '] DEVIATION ', devGlobal

    END IF

  END SUBROUTINE accCompareRealMatrixHostDevice5D
  ! ====================================================================================================== !

END MODULE mo_acc_util

#endif

#ifndef _OPENACC

MODULE mo_acc_util

#ifdef _ICON_STYLE_
  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: accCompareRealMatrixHostDevice

  INTERFACE accCompareRealMatrixHostDevice
    MODULE PROCEDURE accCompareRealMatrixHostDevice1D
    MODULE PROCEDURE accCompareRealMatrixHostDevice2D
    MODULE PROCEDURE accCompareRealMatrixHostDevice3D
    MODULE PROCEDURE accCompareRealMatrixHostDevice4D
    MODULE PROCEDURE accCompareRealMatrixHostDevice5D
  END INTERFACE accCompareRealMatrixHostDevice

CONTAINS

  ! ====================================================================================================== !
  !
  SUBROUTINE accCompareRealMatrixHostDevice1D(matrix, matrixName, deviation, detailedOutput)

    REAL(_REAL_TYPE_), DIMENSION(:), POINTER, INTENT(IN) :: matrix
    CHARACTER(len=*), INTENT(IN) :: matrixName
    REAL(_REAL_TYPE_), OPTIONAL, INTENT(OUT) :: deviation
    LOGICAL, OPTIONAL, INTENT(IN) :: detailedOutput

    LOGICAL :: detailedOutputFlag
    INTEGER :: iL, iU
    REAL(_REAL_TYPE_) :: minMatrix, maxMatrix

    detailedOutputFlag = .FALSE.
    IF (PRESENT(detailedOutput)) detailedOutputFlag = detailedOutput
    IF (PRESENT(deviation)) deviation = _ZERO_DP_

    iL = LBOUND(matrix, 1); iU = UBOUND(matrix, 1)
    minMatrix = MINVAL(matrix)
    maxMatrix = MAXVAL(matrix)

    IF (detailedOutputFlag) THEN

      WRITE(0,*) '>>>>> accCompareRealArrayHostDevice1D'
      WRITE(0,*) '>>>>> MATRIX  NAME: ', TRIM(matrixName)
      WRITE(0,'(1x, A21, I6, A, I6, A, I6, A)') '>>>>> MATRIX SHAPE: (', iL, ':', iU, ')'
      WRITE(0,'(1x, A25, E16.6, A, E16.6, A)') '>>>>>      MIN/MAX: CPU [', &
        & minMatrix, '; ', maxMatrix, ']  GPU NOT USED'

    ELSE

      WRITE(0,'(1x, A6, A, A7, E16.6, A, E16.6, A)') '>>>>> ', TRIM(matrixName), ': CPU [', &
        & minMatrix, '; ', maxMatrix, ']  GPU NOT USED'

    END IF

  END SUBROUTINE accCompareRealMatrixHostDevice1D

  ! ====================================================================================================== !
  !
  SUBROUTINE accCompareRealMatrixHostDevice2D(matrix, matrixName, deviation, detailedOutput)

    REAL(_REAL_TYPE_), DIMENSION(:,:), POINTER, INTENT(IN) :: matrix
    CHARACTER(len=*), INTENT(IN) :: matrixName
    REAL(_REAL_TYPE_), OPTIONAL, INTENT(OUT) :: deviation
    LOGICAL, OPTIONAL, INTENT(IN) :: detailedOutput

    LOGICAL :: detailedOutputFlag
    INTEGER :: iL1, iL2, iU1, iU2
    REAL(_REAL_TYPE_) :: minMatrix, maxMatrix

    detailedOutputFlag = .FALSE.
    IF (PRESENT(detailedOutput)) detailedOutputFlag = detailedOutput
    IF (PRESENT(deviation)) deviation = _ZERO_DP_

    iL1 = LBOUND(matrix, 1); iU1 = UBOUND(matrix, 1)
    iL2 = LBOUND(matrix, 2); iU2 = UBOUND(matrix, 2)
    minMatrix = MINVAL(matrix)
    maxMatrix = MAXVAL(matrix)

    IF (detailedOutputFlag) THEN

      WRITE(0,*) '>>>>> accCompareRealArrayHostDevice2D'
      WRITE(0,*) '>>>>> MATRIX  NAME: ', TRIM(matrixName)
      WRITE(0,'(1x, A21, I6, A, I6, A, I6, A, I6, A)') '>>>>> MATRIX SHAPE: (', iL1, ':', iU1, ') (', iL2, ':', iU2, ')'
      WRITE(0,'(1x, A25, E16.6, A, E16.6, A)') '>>>>>      MIN/MAX: CPU [', minMatrix, '; ', maxMatrix, ']  GPU NOT USED'

    ELSE

      WRITE(0,'(1x, A6, A, A7, E16.6, A, E16.6, A)') '>>>>> ', TRIM(matrixName), ': CPU [', &
        & minMatrix, '; ', maxMatrix, ']  GPU NOT USED'

    END IF

  END SUBROUTINE accCompareRealMatrixHostDevice2D

  ! ====================================================================================================== !
  !
  SUBROUTINE accCompareRealMatrixHostDevice3D(matrix, matrixName, deviation, detailedOutput)

    REAL(_REAL_TYPE_), DIMENSION(:,:,:), POINTER, INTENT(IN) :: matrix
    CHARACTER(len=*), INTENT(IN) :: matrixName
    REAL(_REAL_TYPE_), OPTIONAL, INTENT(OUT) :: deviation
    LOGICAL, OPTIONAL, INTENT(IN) :: detailedOutput

    LOGICAL :: detailedOutputFlag
    INTEGER :: iL1, iL2, iL3, iU1, iU2, iU3
    REAL(_REAL_TYPE_) :: minMatrix, maxMatrix

    detailedOutputFlag = .FALSE.
    IF (PRESENT(detailedOutput)) detailedOutputFlag = detailedOutput
    IF (PRESENT(deviation)) deviation = _ZERO_DP_

    iL1 = LBOUND(matrix, 1); iU1 = UBOUND(matrix, 1)
    iL2 = LBOUND(matrix, 2); iU2 = UBOUND(matrix, 2)
    iL3 = LBOUND(matrix, 3); iU3 = UBOUND(matrix, 3)
    minMatrix = MINVAL(matrix)
    maxMatrix = MAXVAL(matrix)

    IF (detailedOutputFlag) THEN

      WRITE(0,*) '>>>>> accCompareRealArrayHostDevice3D'
      WRITE(0,*) '>>>>> MATRIX  NAME: ', TRIM(matrixName)
      WRITE(0,'(1x, A21, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A)') '>>>>> MATRIX SHAPE: (', iL1, ':', iU1, ') (', &
        & iL2, ':', iU2, ') (', iL3, ':', iU3, ')'
      WRITE(0,'(1x, A25, E16.6, A, E16.6, A)') '>>>>>      MIN/MAX: CPU [', minMatrix, '; ', maxMatrix, ']  GPU NOT USED'

    ELSE

      WRITE(0,'(1x, A6, A, A7, E16.6, A, E16.6, A)') '>>>>> ', TRIM(matrixName), ': CPU [', &
        & minMatrix, '; ', maxMatrix, ']  GPU NOT USED'

    END IF

  END SUBROUTINE accCompareRealMatrixHostDevice3D

  ! ====================================================================================================== !
  !
  SUBROUTINE accCompareRealMatrixHostDevice4D(matrix, matrixName, deviation, detailedOutput)

    REAL(_REAL_TYPE_), DIMENSION(:,:,:,:), POINTER, INTENT(IN) :: matrix
    CHARACTER(len=*), INTENT(IN) :: matrixName
    REAL(_REAL_TYPE_), OPTIONAL, INTENT(OUT) :: deviation
    LOGICAL, OPTIONAL, INTENT(IN) :: detailedOutput

    LOGICAL :: detailedOutputFlag
    INTEGER :: iL1, iL2, iL3, iL4, iU1, iU2, iU3, iU4
    REAL(_REAL_TYPE_) :: minMatrix, maxMatrix

    detailedOutputFlag = .FALSE.
    IF (PRESENT(detailedOutput)) detailedOutputFlag = detailedOutput
    IF (PRESENT(deviation)) deviation = _ZERO_DP_

    iL1 = LBOUND(matrix, 1); iU1 = UBOUND(matrix, 1)
    iL2 = LBOUND(matrix, 2); iU2 = UBOUND(matrix, 2)
    iL3 = LBOUND(matrix, 3); iU3 = UBOUND(matrix, 3)
    iL4 = LBOUND(matrix, 4); iU4 = UBOUND(matrix, 4)
    minMatrix = MINVAL(matrix)
    maxMatrix = MAXVAL(matrix)

    IF (detailedOutputFlag) THEN

      WRITE(0,*) '>>>>> accCompareRealArrayHostDevice4D'
      WRITE(0,*) '>>>>> MATRIX  NAME: ', TRIM(matrixName)
      WRITE(0,'(1x, A21, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A)') '>>>>> MATRIX SHAPE: (', &
       & iL1, ':', iU1, ') (', iL2, ':', iU2, ') (', iL3, ':', iU3, ') (', iL4, ':', iU4, ')'
      WRITE(0,'(1x, A25, E16.6, A, E16.6, A)') '>>>>>      MIN/MAX: CPU [', minMatrix, '; ', maxMatrix, ']  GPU NOT USED'

    ELSE

      WRITE(0,'(1x, A6, A, A7, E16.6, A, E16.6, A)') '>>>>> ', TRIM(matrixName), ': CPU [', &
        & minMatrix, '; ', maxMatrix, ']  GPU NOT USED'

    END IF

  END SUBROUTINE accCompareRealMatrixHostDevice4D

  ! ====================================================================================================== !
  !
  SUBROUTINE accCompareRealMatrixHostDevice5D(matrix, matrixName, deviation, detailedOutput)

    REAL(_REAL_TYPE_), DIMENSION(:,:,:,:,:), POINTER, INTENT(IN) :: matrix
    CHARACTER(len=*), INTENT(IN) :: matrixName
    REAL(_REAL_TYPE_), OPTIONAL, INTENT(OUT) :: deviation
    LOGICAL, OPTIONAL, INTENT(IN) :: detailedOutput

    LOGICAL :: detailedOutputFlag
    INTEGER :: iL1, iL2, iL3, iL4, iL5, iU1, iU2, iU3, iU4, iU5
    REAL(_REAL_TYPE_) :: minMatrix, maxMatrix

    detailedOutputFlag = .FALSE.
    IF (PRESENT(detailedOutput)) detailedOutputFlag = detailedOutput
    IF (PRESENT(deviation)) deviation = _ZERO_DP_

    iL1 = LBOUND(matrix, 1); iU1 = UBOUND(matrix, 1)
    iL2 = LBOUND(matrix, 2); iU2 = UBOUND(matrix, 2)
    iL3 = LBOUND(matrix, 3); iU3 = UBOUND(matrix, 3)
    iL4 = LBOUND(matrix, 4); iU4 = UBOUND(matrix, 4)
    iL5 = LBOUND(matrix, 5); iU5 = UBOUND(matrix, 5)
    minMatrix = MINVAL(matrix)
    maxMatrix = MAXVAL(matrix)

    IF (detailedOutputFlag) THEN

      WRITE(0,*) '>>>>> accCompareRealArrayHostDevice5D'
      WRITE(0,*) '>>>>> MATRIX  NAME: ', TRIM(matrixName)
      WRITE(0,'(1x, A21, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A, I6, A)') &
        & '>>>>> MATRIX SHAPE: (', iL1, ':', iU1, ') (', iL2, ':', iU2, ') &
        & (', iL3, ':', iU3, ') (', iL4, ':', iU4, ') (', iL5, ':', iU5, ')'
      WRITE(0,'(1x, A25, E16.6, A, E16.6, A)') '>>>>>      MIN/MAX: CPU [', minMatrix, '; ', maxMatrix, ']  GPU NOT USED'

    ELSE

      WRITE(0,'(1x, A6, A, A7, E16.6, A, E16.6, A)') '>>>>> ', TRIM(matrixName), ': CPU [', &
        & minMatrix, '; ', maxMatrix, ']  GPU NOT USED'

    END IF

  END SUBROUTINE accCompareRealMatrixHostDevice5D
  ! ====================================================================================================== !

END MODULE mo_acc_util

#endif
