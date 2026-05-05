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
!------------------------------------------------------------------------------
!
! program reading downstream grid cells (FDIR) from HD parameter file
! and generates matrix with upstream grid cells
!
! # Note that Veronikas original program needed some more grid info
!         --> excluded with switch IGRID = 0.
!
!
!------------------------------------------------------------------------------
PROGRAM upstream_cell

  IMPLICIT NONE
  INCLUDE 'netcdf.inc'

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
  INTEGER, PARAMETER :: IGRID = 0
  LOGICAL :: allow_bifurcations = .FALSE.

  INTEGER :: c, cds, ncells, ncells_up, nfilled, nlevels, num_args, i
  INTEGER :: stat, ncid, ncin, ncin2, varid, varndims
  INTEGER :: clonid, clatid, cvlatid, cvlonid, idnc, idni, idnv
  INTEGER :: clonin, clatin, cvlatin, cvlonin
  INTEGER, ALLOCATABLE :: vardimids(:)
  INTEGER, ALLOCATABLE :: num_cells_up(:)
  INTEGER :: flag

  REAL(dp), ALLOCATABLE :: clon(:)
  REAL(dp), ALLOCATABLE :: clon_vertices(:,:)
  REAL(dp), ALLOCATABLE :: clat(:)
  REAL(dp), ALLOCATABLE :: clat_vertices(:,:)
  REAL(dp), ALLOCATABLE :: cell_ds(:)
  REAL(dp), ALLOCATABLE :: cell_ds_bifurcated(:,:)
  REAL(dp), ALLOCATABLE :: cells_up(:,:)

  CHARACTER*100 :: filename,bifurcated_rdirs_filename
  CHARACTER*3  :: output_str,flag_str

  INTEGER, PARAMETER :: nvert = 3   ! number of vertices per grid cell
  INTEGER, PARAMETER :: ninmax = 12  ! 12 neighbors --> 12 = maximum number of inflow cells

  num_args = command_argument_count()
  IF (num_args == 0) THEN
    allow_bifurcations = .FALSE.
  ELSE IF (num_args == 1) THEN
    call get_command_argument(1,value=flag_str)
    READ(flag_str,*) flag
    IF (flag == 0) THEN
      allow_bifurcations = .FALSE.
    ELSE
      allow_bifurcations = .TRUE.
    END IF
  ELSE
    write(*,*) "Wrong number of command line arguments given"
    stop
  END IF
  write(*,*) "allow_bifurcations= ", allow_bifurcations

  filename='hdpara_icon.nc'
  bifurcated_rdirs_filename='bifurcated_next_cell_index_for_upstream_cell.nc'

  !-- read indices of grid cell downstream
  stat = nf_open(TRIM(filename),NF_NOWRITE, ncin)
  CALL hdlerr(stat)
  stat = nf_inq_varid(ncin,'FDIR',varid)
  CALL hdlerr(stat)

  stat = nf_inq_varndims(ncin,varid,varndims)
  CALL hdlerr(stat)
  ALLOCATE(vardimids(varndims))
  stat = nf_inq_vardimid(ncin,varid,vardimids)
  CALL hdlerr(stat)
  stat = nf_inq_dimlen(ncin,vardimids(1),ncells)
  CALL hdlerr(stat)
  DEALLOCATE(vardimids)

  ALLOCATE(cell_ds(ncells))
  stat = nf_get_var_double(ncin,varid,cell_ds)
  CALL hdlerr(stat)

  IF(allow_bifurcations) THEN
      stat = nf_open(TRIM(bifurcated_rdirs_filename),NF_NOWRITE, ncin2)
    CALL hdlerr(stat)
    stat = nf_inq_varid(ncin2,'bifurcated_next_cell_index',varid)
    CALL hdlerr(stat)

    stat = nf_inq_varndims(ncin2,varid,varndims)
    CALL hdlerr(stat)
    ALLOCATE(vardimids(varndims))
    stat = nf_inq_vardimid(ncin2,varid,vardimids)
    CALL hdlerr(stat)
    stat = nf_inq_dimlen(ncin2,vardimids(1),ncells)
    CALL hdlerr(stat)
    stat = nf_inq_dimlen(ncin2,vardimids(2),nlevels)
    CALL hdlerr(stat)
    DEALLOCATE(vardimids)
    ALLOCATE(cell_ds_bifurcated(ncells,nlevels))
    stat = nf_get_var_double(ncin2,varid,cell_ds_bifurcated)
    CALL hdlerr(stat)
  END IF

  ALLOCATE(cells_up(ncells,ninmax))
  cells_up(:,:) = -1

  ALLOCATE(num_cells_up(ncells))
  num_cells_up(:) = 0

  !-- find upstream cells
  print*, 'Finding upstream cells: START '

  DO c = 1, ncells
    cds = cell_ds(c)
    if(cds < 1) cycle
    ncells_up = num_cells_up(cds) + 1
    cells_up(cds,ncells_up) = c
    num_cells_up(cds) =  ncells_up
  END DO
  IF (allow_bifurcations) THEN
    DO i = 1, 11
      IF (ANY(cell_ds_bifurcated(:,i) > 0)) THEN
        DO c = 1, ncells
          cds = cell_ds_bifurcated(c,i)
          if(cds < 1) cycle
          ncells_up = num_cells_up(cds) + 1
          cells_up(cds,ncells_up) = c
          num_cells_up(cds) =  ncells_up
        END DO
      END IF
    END DO
  END IF
  nfilled = maxval(num_cells_up)
  write (output_str,'(I2)') nfilled
  print*, 'Highest number of upstream cells: '//output_str
  print*, 'Finding upstream cells: END'

  !-- define the output file

  stat = nf_inq_varid(ncin,'clon',clonin)
  CALL hdlerr(stat)
  ALLOCATE(clon(ncells))
  stat = nf_get_var_double(ncin,clonin,clon)
  CALL hdlerr(stat)

  IF (IGRID.GT.0.5) THEN
    print*, 'get clon_vertices '
    stat = nf_inq_varid(ncin,'clon_vertices',cvlonin)
    CALL hdlerr(stat)
    ALLOCATE(clon_vertices(ncells,nvert))
    stat = nf_get_var_double(ncin,cvlonin,clon_vertices)
    CALL hdlerr(stat)
  ENDIF

  print*, 'get clat '
  stat = nf_inq_varid(ncin,'clat',clatin)
  CALL hdlerr(stat)
  ALLOCATE(clat(ncells))
  stat = nf_get_var_double(ncin,clatin,clat)
  CALL hdlerr(stat)

  IF (IGRID.GT.0.5) THEN
    stat = nf_inq_varid(ncin,'clat_vertices',cvlatin)
    CALL hdlerr(stat)
    ALLOCATE(clat_vertices(ncells,nvert))
    stat = nf_get_var_double(ncin,cvlatin,clat_vertices)
  ENDIF

  stat = nf_create('upstream_cells.nc',OR(NF_CLOBBER,NF_NETCDF4),ncid)
  CALL hdlerr(stat)
  stat = nf_def_dim(ncid, 'ncells', ncells, idnc)
  CALL hdlerr(stat)
  stat = nf_def_dim(ncid, 'vertices', nvert, idnv)
  CALL hdlerr(stat)
  stat = nf_def_dim(ncid, 'nneigh', ninmax, idni)
  CALL hdlerr(stat)

  print*, 'get clon '
  stat = nf_def_var(ncid, 'clon', NF_DOUBLE, 1, (/idnc/), clonid)
  CALL hdlerr(stat)
  stat = nf_copy_att(ncin, clonin, 'standard_name', ncid, clonid)
  CALL hdlerr(stat)
  stat = nf_copy_att(ncin, clonin, 'long_name', ncid, clonid)
  CALL hdlerr(stat)
  stat = nf_copy_att(ncin, clonin, 'units', ncid, clonid)
  CALL hdlerr(stat)
  stat = nf_copy_att(ncin, clonin, 'bounds', ncid, clonid)
  CALL hdlerr(stat)

  IF (IGRID.GT.0.5) THEN
    stat = nf_def_var(ncid, 'clon_vertices', NF_DOUBLE, 2, (/idnv,idnc/), cvlonid)
    CALL hdlerr(stat)
  ENDIF
  stat = nf_def_var(ncid, 'clat', NF_DOUBLE, 1, (/idnc/), clatid)
  CALL hdlerr(stat)
  stat = nf_copy_att(ncin, clatin, 'standard_name', ncid, clatid)
  CALL hdlerr(stat)
  stat = nf_copy_att(ncin, clatin, 'long_name', ncid, clatid)
  CALL hdlerr(stat)
  stat = nf_copy_att(ncin, clatin, 'units', ncid, clatid)
  CALL hdlerr(stat)
  stat = nf_copy_att(ncin, clatin, 'bounds', ncid, clatid)
  CALL hdlerr(stat)
  IF (IGRID.GT.0.5) THEN
    stat = nf_def_var(ncid, 'clat_vertices', NF_DOUBLE, 2, (/idnv,idnc/), cvlatid)
    CALL hdlerr(stat)
  ENDIF
  stat = nf_def_var(ncid, 'CELLS_UP', NF_DOUBLE, 2, (/idnc,idni/), varid)
  CALL hdlerr(stat)
  stat = nf_put_att_text(ncid, varid, 'long_name', 19, 'grid cells upstream')
  CALL hdlerr(stat)
  stat = nf_put_att_text(ncid, varid, 'grid_type', 12, 'unstructured')
  CALL hdlerr(stat)
  stat = nf_put_att_text(ncid, varid, 'coordinates', 9, 'clon clat')
  CALL hdlerr(stat)
  stat = nf_enddef(ncid)
  CALL hdlerr(stat)

  !-- write output file variables

  print*, 'write clon '
  stat = nf_put_var_double(ncid,clonid,clon)
  CALL hdlerr(stat)

  IF (IGRID.GT.0.5) THEN
    stat = nf_put_var_double(ncid,cvlonid,clon_vertices)
    CALL hdlerr(stat)
  ENDIF
  stat = nf_put_var_double(ncid,clatid,clat)
  CALL hdlerr(stat)
  IF (IGRID.GT.0.5) THEN
    stat = nf_put_var_double(ncid,cvlatid,clat_vertices)
    CALL hdlerr(stat)
  ENDIF
  stat = nf_put_var_double(ncid,varid,cells_up)
  CALL hdlerr(stat)

  stat = nf_close(ncid)
  CALL hdlerr(stat)

  stat = nf_close(ncin)
  CALL hdlerr(stat)

  IF (allow_bifurcations) THEN
    stat = nf_close(ncin2)
    CALL hdlerr(stat)
  END IF

END PROGRAM upstream_cell

!------------------------------------------------------------------------------
!
!  Routine to handle netcdf errors
!
!------------------------------------------------------------------------------
SUBROUTINE hdlerr(stat)

  IMPLICIT NONE

  include 'netcdf.inc'

! INTENT(in)
  INTEGER,       INTENT(in) :: stat

  IF (stat /= NF_NOERR) THEN
     WRITE (6,*) '--------'
     WRITE (6,*) ' ERROR:  ', nf_strerror(stat)
     WRITE (6,*) '--------'
     STOP
  END IF

END SUBROUTINE hdlerr
!------------------------------------------------------------------------------
