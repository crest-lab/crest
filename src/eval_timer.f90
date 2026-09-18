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
subroutine eval_timer(tim)
!********************************
!* The final timer evaluation to
!* be called at the end of CREST
!********************************
  use crest_parameters
  use crest_data
  use crest_calculator,only:engrad_total
  use crest_restartlog
  use iomod,only:get_peak_rss_kb,to_str
  implicit none
  type(timer) :: tim
  real(wp) :: time_total,time_avg,mem
  character(len=40) :: atmp
  write (stdout,*)
  call smallhead('Wall Time Summary')
  call tim%write(stdout,'CREST runtime',verbose=.true.)
  time_total = tim%get()
  call tim%clear
  mem = real(get_peak_rss_kb(),wp)
  write (stdout,'(" * Peak RSS: ",f8.2, " MiB")') mem/1024.0_wp
  if (engrad_total > 0.0_wp) then
    write (atmp,'(f30.3)') time_total/real(engrad_total,wp)
    if (engrad_total < 10.0_wp**5) then
      write (stdout,'(" * Total number of energy+grad calls: ",i0)') &
      &  nint(engrad_total)
    else
      write (stdout,'(" * Total number of energy+grad calls: ",es11.4)') &
      &  engrad_total
    end if
    write (stdout,*)
    !call dump_restart()
  end if
end subroutine eval_timer

subroutine propquit(tim)
  use crest_parameters,only:stdout
  use crest_data
  implicit none
  type(timer) :: tim
  call eval_timer(tim)
  call creststop(status_normal)
end subroutine propquit
