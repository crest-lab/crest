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
  use crest_calculator,only: engrad_total
  use crest_restartlog
  implicit none
  type(timer) :: tim
  real(wp) :: time_total,time_avg
  character(len=40) :: atmp
  write (stdout,*)
  call smallhead('Wall Time Summary')
  call tim%write(stdout,'CREST runtime',verbose=.true.)
  time_total = tim%get()
  call tim%clear
  if(engrad_total > 0)then
  write(atmp,'(f30.3)') time_total/real(engrad_total,wp)
  write(stdout,'(" * Total number of energy+grad calls: ",i0)') & !,a,1x,a,a)') & 
  &  engrad_total!,' (avg. wall-time',trim(adjustl(atmp)),' sec)'
  write(stdout,*)
  call dump_restart() 
  endif
end subroutine eval_timer

subroutine propquit(tim)
  use crest_parameters, only: stdout
  use crest_data
  implicit none
  type(timer) :: tim
  call eval_timer(tim)
  write (stdout,*) 'CREST terminated normally.'
  stop
end subroutine propquit
