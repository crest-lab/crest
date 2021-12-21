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
  use constraints
  implicit none

  character(len=1),public,parameter :: sep = '/'
  character(len=12),public,parameter :: dev0 = ' 2>/dev/null'


  !=========================================================================================!
  type :: calculation_settings

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

    !>--- results/property requests
    real(wp) :: epot
    real(wp) :: efix
    real(wp) :: etot

    logical :: rdwbo = .false.
    real(wp),allocatable :: wbo(:,:)

    logical :: rddip = .false.
    real(wp) :: dip(3)
    logical :: rddipgrad = .false.

  end type calculation_settings
  !=========================================================================================!


  public :: calcdata
  !=========================================================================================!
  !data object that contains settings for a calculation.
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

    integer :: type 

    !>--- calculations
    integer :: ncalculations = 0
    type(calculation_settings),allocatable :: calcs(:)

    !>--- constraints
    integer :: nconstraints = 0
    type(constraint),allocatable :: cons(:)

    !>--- results/property requests
    real(wp) :: epot
    real(wp) :: efix
    real(wp) :: etot

    logical :: rdwbo = .false.
    real(wp),allocatable :: wbo(:,:)

    logical :: rddip = .false.
    real(wp) :: dip(3)
    logical :: rddipgrad = .false.
    real(wp),allocatable :: dipgrad(:,:)

  contains
    procedure :: reset => calculation_reset
    generic,public :: add => calculation_add_constraint,calculation_add_settings
    procedure,private :: calculation_add_constraint,calculation_add_settings

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

    if (allocated(self%calcs)) deallocate (self%calcs)
    self%ncalculations = 0

    if (allocated(self%cons)) deallocate (self%cons)
    self%nconstraints = 0

    self%rdwbo = .false.
    self%rddip = .false.

    return
  end subroutine calculation_reset

!==========================================================================================!

  subroutine calculation_add_settings(self,cal)
    implicit none
    class(calcdata) :: self
    type(calculation_settings) :: cal
    type(calculation_settings),allocatable :: callist(:)
    integer :: i,j

    if (self%ncalculations < 1) then
      allocate (self%calcs(1))
      self%ncalculations = 1
      self%calcs(1) = cal
    else
      i = self%ncalculations + 1
      j = self%ncalculations
      allocate (callist(i))
      callist(1:j) = self%calcs(1:j)
      callist(i) = cal
      call move_alloc(callist,self%calcs)
      self%ncalculations = i
    end if

    return
  end subroutine calculation_add_settings

!==========================================================================================!

  subroutine calculation_add_constraint(self,constr)
    implicit none
    class(calcdata) :: self
    type(constraint) :: constr
    type(constraint),allocatable :: conslist(:)
    integer :: i,j

    if (self%nconstraints < 1) then
      allocate (self%cons(1))
      self%nconstraints = 1
      self%cons(1) = constr
    else
      i = self%nconstraints + 1
      j = self%nconstraints
      allocate (conslist(i))
      conslist(1:j) = self%cons(1:j)
      conslist(i) = constr
      call move_alloc(conslist,self%cons)
      self%nconstraints = i
    end if

    return
  end subroutine calculation_add_constraint

!==========================================================================================!
end module calc_type
