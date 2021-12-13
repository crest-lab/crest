!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2022 Philipp Pracht
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

module calc_type
  use iso_fortran_env,only:wp => real64

  implicit none

  character(len=1),public,parameter :: sep = '/'
  character(len=12),public,parameter :: dev0 = ' 2>/dev/null'

  public :: calcdata
  !=========================================================================================!
  !data object that contains settings for a QM/SQM/FF calculation.
  type :: calcdata

    integer :: id

    integer :: molcharge = 0
    integer :: uhf = 0

    character(len=:),allocatable :: calcspace
    character(len=:),allocatable :: calcfile
    character(len=:),allocatable :: gradfile
    character(len=:),allocatable :: path
    character(len=:),allocatable :: other
    character(len=:),allocatable :: binary
    character(len=:),allocatable :: systemcall

    !>--- property requests
    logical :: rdwbo = .false.
    real(wp),allocatable :: wbo(:,:)

    logical :: rddip = .false.
    real(wp) :: dip(3)
    logical :: rddipgrad = .false.
    real(wp),allocatable :: dipgrad(:,:)

  contains
    procedure :: reset => calculation_reset

  end type calcdata

contains
!==========================================================================================!

  subroutine calculation_reset(self)
    implicit none
    class(calcdata) :: self

    if (allocated(self%calcspace)) deallocate (self%calcspace)
    if (allocated(self%calcfile)) deallocate (self%calcfile)
    if (allocated(self%gradfile)) deallocate (self%gradfile)
    if (allocated(self%path)) deallocate (self%path)
    if (allocated(self%other)) deallocate (self%other)
    if (allocated(self%binary)) deallocate (self%binary)
    if (allocated(self%systemcall)) deallocate (self%systemcall)
    if (allocated(self%wbo)) deallocate (self%wbo)
    if (allocated(self%dipgrad)) deallocate (self%dipgrad)

    self%rdwbo = .false.
    self%rddip = .false.

    return
  end subroutine calculation_reset

!==========================================================================================!
end module calc_type
