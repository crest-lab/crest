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
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use constraints
  implicit none

  character(len=1),public,parameter :: sep = '/'
  character(len=12),public,parameter :: dev0 = ' 2>/dev/null'


   !> job type enumerator
   type ,private:: enum_jobtype
      integer :: unknown   = 0
      integer :: xtbsys    = 1
      integer :: generic   = 2
      integer :: turbomole = 3
      integer :: orca      = 4
      integer :: terachem  = 5
   end type enum_jobtype
   type(enum_jobtype), parameter,public :: jobtype = enum_jobtype()

!=========================================================================================!
  !> data object that contains the data for a *SINGLE* calculation
  public :: calculation_settings
  type :: calculation_settings

    integer :: id

    integer :: chrg = 0
    integer :: uhf = 0

    character(len=:),allocatable :: calcspace
    character(len=:),allocatable :: calcfile
    character(len=:),allocatable :: gradfile
    character(len=:),allocatable :: path
    character(len=:),allocatable :: other
    character(len=:),allocatable :: binary
    character(len=:),allocatable :: systemcall

    !>--- gradient format specifications
    integer :: gradtype = 0
    integer :: gradfmt = 0
    character(len=:),allocatable :: gradkey
    character(len=:),allocatable :: efile

    !>--- results/property requests
    real(wp) :: epot = 0.0_wp
    real(wp) :: efix = 0.0_wp
    real(wp) :: etot = 0.0_wp

    logical :: rdwbo = .false.
    real(wp),allocatable :: wbo(:,:)

    logical :: rddip = .false.
    real(wp) :: dip(3) = 0
    logical :: rddipgrad = .false.
    real(wp),allocatable :: dipgrad(:,:,:)

  contains
    procedure :: deallocate => calculation_settings_deallocate

  end type calculation_settings
!=========================================================================================!

!=========================================================================================!
  !> data object that contains settings for *ALL* calculations and constraints.
  public :: calcdata
  type :: calcdata
    integer :: id = 0
    integer :: which = 0

    !>--- calculations
    integer :: ncalculations = 0
    type(calculation_settings),allocatable :: calcs(:)
    real(wp),allocatable :: etmp(:)
    real(wp),allocatable :: grdtmp(:,:,:)

    !>--- constraints
    integer :: nconstraints = 0
    type(constraint),allocatable :: cons(:)

    !>--- results/property requests
    real(wp) :: epot
    real(wp) :: efix
    real(wp) :: etot

    !>--- optimization settings
    integer  :: optlev = 0
    integer  :: micro_opt = 20
    integer  :: maxcycle = 0
    real(wp) :: maxdispl_opt = 1.0_wp
    real(wp) :: hlow_opt = 0.010_wp
    real(wp) :: hmax_opt = 5.0_wp
    real(wp) :: acc_opt = 1.0_wp
    logical  :: exact_rf = .false.
    logical  :: average_conv = .false.
    logical  :: tsopt = .false.
    integer  :: iupdat = 0  !> 0=BFGS, 1=Powell, 2=SR1, 3=Bofill, 4=Schlegel

    !>--- printouts and io
    logical :: pr_energies = .false.
    integer :: eout_unit = stdout
    character(len=:),allocatable :: elog
  contains
    procedure :: reset => calculation_reset
    procedure :: init => calculation_init
    generic,public :: add => calculation_add_constraint,calculation_add_settings
    procedure,private :: calculation_add_constraint,calculation_add_settings
    procedure :: copy => calculation_copy
    procedure :: printconstraints => calculation_print_constraints
    procedure :: removeconstraint => calculation_remove_constraint
  end type calcdata
!=========================================================================================!

contains
!=========================================================================================!

  subroutine calculation_reset(self)
    implicit none
    class(calcdata) :: self

    self%id = 0
    self%which = 0

    if (allocated(self%calcs)) deallocate (self%calcs)
    self%ncalculations = 0

    if (allocated(self%cons)) deallocate (self%cons)
    self%nconstraints = 0

    self%optlev = 0
    self%micro_opt = 20
    self%maxcycle = 0
    self%maxdispl_opt = 1.0_wp
    self%hlow_opt = 0.010_wp
    self%hmax_opt = 5.0_wp
    self%acc_opt = 1.0_wp
    self%exact_rf = .false.
    self%average_conv = .false.
    self%tsopt = .false.
    self%iupdat = 0

    self%pr_energies = .false.
    self%eout_unit = stdout
    if (allocated(self%elog)) deallocate (self%elog)

    return
  end subroutine calculation_reset

  subroutine calculation_settings_deallocate(self)
    implicit none
    class(calculation_settings) :: self

    if (allocated(self%calcspace)) deallocate (self%calcspace)
    if (allocated(self%calcfile)) deallocate (self%calcfile)
    if (allocated(self%gradfile)) deallocate (self%gradfile)
    if (allocated(self%path)) deallocate (self%path)
    if (allocated(self%other)) deallocate (self%other)
    if (allocated(self%binary)) deallocate (self%binary)
    if (allocated(self%systemcall)) deallocate (self%systemcall)
    if (allocated(self%wbo)) deallocate (self%wbo)
    if (allocated(self%dipgrad)) deallocate (self%dipgrad)
    if (allocated(self%gradkey)) deallocate (self%gradkey)
    if (allocated(self%efile)) deallocate (self%efile)

    self%id = 0
    self%chrg = 0
    self%uhf = 0
    self%epot = 0.0_wp
    self%efix = 0.0_wp
    self%etot = 0.0_wp

    self%rdwbo = .false.
    self%rddip = .false.
    self%dip = 0.0_wp
    self%rddipgrad = .false.
    self%gradtype = 0
    self%gradfmt = 0

    return
  end subroutine calculation_settings_deallocate

!=========================================================================================!

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

!=========================================================================================!

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

!=========================================================================================!

  subroutine calculation_remove_constraint(self,d)
    implicit none
    class(calcdata) :: self
    type(constraint) :: constr
    type(constraint),allocatable :: conslist(:)
    integer :: i,j,d,d1,d2

    if (self%nconstraints < d) return

    i = self%nconstraints - 1
    j = self%nconstraints
    allocate (conslist(i))
    if(d == 1)then
    conslist(1:i) = self%cons(2:j)
    else if( d == j)then
    conslist(1:i) = self%cons(1:i)    
    else
    d1=d-1
    d2=d+1
    conslist(1:d1) = self%cons(1:d1)
    conslist(d:i) = self%cons(d2:j)
    endif
    call move_alloc(conslist,self%cons)
    self%nconstraints = i

    return
  end subroutine calculation_remove_constraint



!=========================================================================================!
  subroutine calculation_print_constraints(self,chnl)
    implicit none
    class(calcdata) :: self
    integer :: i,j
    integer,optional :: chnl

    if (self%nconstraints < 1) then
      return
    else
      if(present(chnl))then
      i = chnl
      else
      i = stdout
      endif
      !write(i,*)
      do j =1,self%nconstraints
        call self%cons(j)%print(i)
      enddo
    end if

    return
  end subroutine calculation_print_constraints

!=========================================================================================!

  subroutine calculation_init(self)
    class(calcdata) :: self

    if (allocated(self%elog)) then
      self%pr_energies = .true.
      open (newunit=self%eout_unit,file=self%elog)
    end if

  end subroutine calculation_init

!=========================================================================================!
!> copy a calcdata object from src to self
  subroutine calculation_copy(self,src)
    class(calcdata) :: self
    type(calcdata) :: src
    integer :: i
    !call self%reset()

    self%id = src%id
    self%which = src%which

    self%ncalculations = src%ncalculations
    if (allocated(self%calcs)) deallocate (self%calcs)
    !self%calcs = src%calcs
    do i = 1,self%ncalculations
      call self%add(src%calcs(i))
    enddo

    self%nconstraints = src%nconstraints
    if (allocated(self%cons)) deallocate (self%cons)
    !self%cons = src%cons
    do i = 1,self%nconstraints
      call self%add(src%cons(i))
    enddo

    self%optlev = src%optlev
    self%micro_opt = src%micro_opt
    self%maxcycle = src%maxcycle
    self%maxdispl_opt = src%maxdispl_opt
    self%hlow_opt = src%hlow_opt
    self%hmax_opt = src%hmax_opt
    self%acc_opt = src%acc_opt
    self%exact_rf = src%exact_rf
    self%average_conv = src%average_conv
    self%tsopt = src%tsopt
    self%iupdat = src%iupdat

    self%pr_energies = src%pr_energies 
    self%eout_unit = src%eout_unit
    self%elog = src%elog

    return
  end subroutine calculation_copy


!=========================================================================================!
end module calc_type
