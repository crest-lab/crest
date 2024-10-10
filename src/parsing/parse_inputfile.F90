!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Philipp Pracht
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

module parse_inputfile
  use crest_parameters
  use crest_data
  use parse_datastruct
  use parse_toml
  implicit none
  private

  public :: parse_test
  public :: parse_input

  external creststop

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine parse_test(fname)
    implicit none
    character(len=*) :: fname !> name of the input file
    type(root_object) :: dict
    call parse_input(fname,dict)
    call dict%print()
    return
  end subroutine parse_test

!========================================================================================!

  subroutine parse_input(fname,dict)
    implicit none

    character(len=*) :: fname !> name of the input file
    type(root_object),intent(out) :: dict

    character(len=:),allocatable :: extension
    integer :: p
    logical :: ex

!>-- check for existance of file
    inquire (file=fname,exist=ex)
    if (.not.ex) then
      return
    end if
    dict%filename = trim(fname)

!>-- get file-type extension
    p = index(fname,'.',.true.)
    if (p .ne. 0) then
      extension = trim(fname(p:))
    else
      extension = 'none'
    end if

    select case (extension)
    case default
      write (stdout,'(a,a)') 'Unknown file format of input file ',trim(fname)
      call creststop(status_input)
    case ('.toml')
#ifdef WITH_TOMLF
!>--- parse .toml file via the toml-f library (the DEFAULT setting)
      call parse_tomlf(fname,dict)
#else
!>--- or, as a fallback only, via the intrinisc routine
      call parse_toml_input_fallback(fname,dict)
#endif
    end select

  end subroutine parse_input

!========================================================================================!
end module parse_inputfile
