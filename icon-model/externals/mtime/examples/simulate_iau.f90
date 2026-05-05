!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
PROGRAM simulate_iau

  USE mtime

  IMPLICIT NONE

  TYPE(datetime) :: start_date, stop_date, current_date, previous_date
  TYPE(timedelta), POINTER :: time_step, iau_time_shift

  LOGICAL, PARAMETER :: iterate_iau = .TRUE.
  INTEGER :: iau_iter

  REAL :: dt_shift, dtime

  INTEGER :: jstep, jstep0, jstep_shift

  CHARACTER(len=max_datetime_str_len) :: dstring

  WRITE (0, *) "Start execution, set calendar ..."

  CALL setCalendar(proleptic_gregorian)

  WRITE (0, *) "Assign values ..."

  start_date = newDatetime("2016-01-01T00:00:00")
  stop_date = newDatetime("2016-01-02T00:00:00")

  time_step => newTimedelta("PT15M")
  dtime = 900.0

  iau_time_shift => newTimedelta("-PT1H30M")
  dt_shift = -5400.0

  WRITE (0, *) "Prepare time loop ..."

  current_date = start_date
  current_date = current_date + iau_time_shift

  CALL datetimeToString(start_date, dstring)
  WRITE (0, *) '           start date ', dstring
  CALL datetimeToString(current_date, dstring)
  WRITE (0, *) '   shifted start date ', dstring

  iau_iter = MERGE(1, 0, iterate_iau)

  IF (iterate_iau) THEN
    jstep_shift = NINT(dt_shift/dtime)
  ELSE
    jstep_shift = 0
  END IF

  previous_date = current_date

  jstep0 = 0
  jstep = (jstep0 + 1) + jstep_shift

  WRITE (0, *) "Start time loop ..."

  time_loop: DO

    current_date = current_date + time_step

    CALL datetimeToString(current_date, dstring)
    WRITE (0, *) "   Time loop ", dstring, jstep

    IF ((current_date%date%day /= previous_date%date%day) .AND. .NOT. (jstep == 0 .AND. iau_iter == 1)) THEN
      previous_date = current_date
    END IF

    WRITE (0, *) '   --- integrate nh input: ', dstring, jstep - jstep_shift, iau_iter

    IF (current_date >= stop_date) THEN
      EXIT time_loop
    END IF

    IF (jstep == 0 .AND. iau_iter == 1) THEN
      iau_iter = 2
      jstep = (jstep0 + 1) + jstep_shift
      current_date = current_date + iau_time_shift
    ELSE
      jstep = jstep + 1
    END IF

  END DO time_loop

  CALL deallocateTimeDelta(time_step)
  CALL deallocateTimeDelta(iau_time_shift)

END PROGRAM simulate_iau
