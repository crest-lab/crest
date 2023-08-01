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
  use crest_parameters
  use crest_data
  implicit none
  type(timer) :: tim
  write (stdout,*)
  call smallhead('Wall Time Summary')
  call tim%write(stdout,'CREST runtime',verbose=.true.)
  call tim%clear
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
