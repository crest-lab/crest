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
  public :: find_input_file

  external creststop

  integer,parameter :: nf = 2
  character(len=*),parameter :: ftypes(nf) = [&
          & '.toml','.TOML']
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
    case ('.toml','.TOML')
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

  subroutine find_input_file(args,nra,pos)
    !*******************************************
    !* A routine to look up an input file from
    !* the list of command line arguments 
    !*******************************************
    implicit none
    integer,intent(in) :: nra
    character(len=*) :: args(nra)
    integer,intent(out) :: pos

    character(len=:),allocatable :: argument
    logical,allocatable :: isinputfile(:)
    logical,allocatable :: inputprio(:)
    integer :: i,j,k,l
    logical :: ex

    pos = 0
    allocate (isinputfile(nra),source=.false.)
    allocate (inputprio(nra),source=.false.)

    do i = 1,nra
      argument = args(i)
      inquire (file=argument,exist=ex)
      if (ex) then
        do j = 1,nf
          if (index(argument,ftypes(j)) .ne. 0) then
            isinputfile(i) = .true.
            if (i > 1) then
              if (args(i-1) == '--input'.or.args(i-1) == '-i') then
                inputprio(i) = .true.
              end if
            end if
          end if
        end do
      end if
    end do

    !> if there are multiple inputs given, we select the last one,
    !> except if it was provided with --input/-i
    !> (same logic applies for multiple --input/-i)
    if (any(inputprio(:))) then
      do i = 1,nra
        if (inputprio(i)) pos = i
      end do
    elseif (any(isinputfile(:))) then
      do i = 1,nra
        if (isinputfile(i)) pos = i
      end do
    else
      pos = 0
    end if

  end subroutine find_input_file
!========================================================================================!
end module parse_inputfile
