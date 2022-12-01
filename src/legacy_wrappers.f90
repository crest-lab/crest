!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022 Philipp Pracht
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

subroutine confscript2i(env,tim)
  use iso_fortran_env,only:wp => real64
  use crest_data
  implicit none
  type(systemdata) :: env
  type(timer)   :: tim
  if(env%legacy)then
      call confscript2i_legacy(env,tim)
  else
      call crest_search_imtdgc(env,tim)
  endif
end subroutine confscript2i
