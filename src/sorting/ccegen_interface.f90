!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020 Philipp Pracht, Stefan Grimme
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

module ccegen_interface
!**********************************************************
!* module to load an interface to the CCEGEN routine.    *
!**********************************************************
  implicit none
  interface
    subroutine CCEGEN(env,pr,fname)
      use crest_parameters
      use crest_data
      implicit none
      type(systemdata),intent(inout) :: env
      logical,intent(in) :: pr
      character(len=*),intent(in) :: fname
    end subroutine CCEGEN
  end interface
end module ccegen_interface
