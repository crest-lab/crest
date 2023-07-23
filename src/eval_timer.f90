!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021-2023 Philipp Pracht
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! crest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

!c  Transform a time in seconds to a sting in the format hh:mm:ss
subroutine eval_time(time,printout)
  use iso_fortran_env,dp => int64,wp => real64
  implicit none
  real(wp) :: time
  character(len=*) :: printout
  integer(dp) :: hours,minutes
  real(wp) ::  seconds

  hours = floor(time/3600.0d0)
  time = time-(real(hours)*3600.0d0)
  minutes = floor(time/60.0d0)
  time = time-(real(minutes)*60.0d0)
!  seconds = floor(time)
  seconds = time

  write (printout,'(i0,a,i2,a,f6.3,a)') &
  &     hours,'h :',minutes,'m :',seconds,'s'
end subroutine eval_time

subroutine eval_sub_timer(tim)
  use iso_fortran_env,wp => real64,dp => int64
  use crest_data
  implicit none
  type(timer) :: tim
  character(len=64) :: ftime
  integer(dp) ::  t1
  real(wp) :: t2
  integer :: i,j
  j = tim%times
  if (j .lt. 1) return
  do i = 1,j
    t1 = tim%t(i,3)
    if (t1 .ne. 0) then
      t2 = real(t1)/real(tim%rate)
      call eval_time(t2,ftime)
      if (trim(tim%names(i)) .ne. '') then
        write (*,'(a20,'' wall time : '',a20)') trim(tim%names(i)),trim(ftime)
      end if
    else
      cycle
    end if
  end do
  return
end subroutine eval_sub_timer

subroutine eval_timer(tim)
  use iso_fortran_env,wp => real64,dp => int64
  use crest_data
  implicit none
  type(timer) :: tim
  character(len=64) :: ftime
  integer(dp) ::  ttot1
  real(wp) :: ttot2
  integer :: j
  write (*,*)
  call smallhead('Wall Time Summary')
  call eval_sub_timer(tim)
  j = tim%times
  ttot1 = sum(tim%t(1:j,3))
  ttot2 = real(ttot1)/real(tim%rate)
  call eval_time(ttot2,ftime)
  write (*,'(''--------------------'')')
  write (*,'(''Overall wall time  : '',a)') trim(ftime)

  call tim%clear
end subroutine eval_timer

subroutine propquit(tim)
  use crest_data
  implicit none
  type(timer) :: tim
  call eval_timer(tim)
  write (*,*)
  write (*,*) 'CREST terminated normally.'
  stop
end subroutine propquit

