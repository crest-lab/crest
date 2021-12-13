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

module calc_module

  use iso_fortran_env,only:wp => real64
  use strucrd
  use calc_type
  use xtb_sc
  use lj

  implicit none

  !=========================================================================================!
  !--- private module variables and parameters
  private
  integer :: i,j,k,l,ich,och,io
  logical :: ex

  !--- some constants and name mappings
  real(wp),parameter :: bohr = 0.52917726_wp
  real(wp),parameter :: autokcal = 627.509541_wp

  public :: engrad
  interface engrad
    module procedure :: engrad_mol
    module procedure :: engrad_xyz
  end interface engrad

  public :: test_engrad

contains

  subroutine test_engrad(fname)
    implicit none

    character(len=*) :: fname

    type(coord) :: mol
    type(coord) :: molo
    type(calcdata) :: dat

    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)

    logical :: pr

    pr = .false.

    call mol%open(fname)

    allocate (grad(3,mol%nat),source=0.0_wp)

    !dat%id = 99
    !dat%other = '15.0   1.5'
    dat%id = 10
    dat%calcspace = 'testdir'
    call engrad(mol,dat,energy,grad,io)

    write (*,*) 'Energy: ',energy
    write (*,*) 'Gradient:'
    do i = 1,mol%nat
      write (*,'(3f18.8)') grad(1:3,i)
    end do

    write (*,*)
    molo = mol
    if (pr) then
      call numgrad(molo,dat)
    end if
    deallocate (grad)

    return
  end subroutine test_engrad

!========================================================================================!
! subroutine engrad
! main routine to perform some energy and gradient calculation
  subroutine engrad_mol(mol,dat,energy,gradient,iostatus)
    implicit none
    type(coord) :: mol
    type(calcdata) :: dat
    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: gradient(3,mol%nat)
    integer,intent(out) :: iostatus

    real(wp) :: dum1,dum2

    call initsignal()

    iostatus = 0
    dum1 = 1.0_wp
    dum2 = 1.0_wp

    select case (dat%id)
    case (10) !-- xtb system call
      call xtb_engrad(mol,dat,energy,gradient,iostatus)
    case (99) !-- Lennard-Jones dummy calculation
      if (allocated(dat%other)) then
        read (dat%other,*) dum1,dum2
      end if
      call lj_engrad(mol%nat,mol%xyz,dum1,dum2,energy,gradient)
    case default
      write (*,*) 'Nothing selected for energy and gradient calculation.'
    end select

    return
  end subroutine engrad_mol

  subroutine engrad_xyz(n,xyz,at,dat,energy,gradient,iostatus)
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n) ! coordinates should be in Bohr (a.u.)
    integer,intent(in) :: at(n)
    type(calcdata) :: dat
    type(coord) :: mol

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: gradient(3,n)
    integer,intent(out) :: iostatus

    real(wp) :: dum1,dum2

    mol%nat = n
    mol%at = at
    mol%xyz = xyz
    call engrad_mol(mol,dat,energy,gradient,iostatus)

    return
  end subroutine engrad_xyz

!========================================================================================!
! subroutine numgrad
! routine to perform a numerical gradient calculation
  subroutine numgrad(mol,dat)
    implicit none

    type(coord) :: mol
    type(calcdata) :: dat

    real(wp) :: energy,el,er
    real(wp),allocatable :: grad(:,:)
    real(wp),allocatable :: numgrd(:,:)
    real(wp),parameter :: step = 0.00001_wp,step2 = 0.5_wp / step

    allocate (grad(3,mol%nat),source=0.0_wp)
    allocate (numgrd(3,mol%nat),source=0.0_wp)

    do i = 1,mol%nat
      do j = 1,3
        mol%xyz(j,i) = mol%xyz(j,i) + step
        call engrad(mol%nat,mol%xyz,mol%at,dat,er,grad,io)

        mol%xyz(j,i) = mol%xyz(j,i) - 2 * step
        call engrad(mol%nat,mol%xyz,mol%at,dat,el,grad,io)

        mol%xyz(j,i) = mol%xyz(j,i) + step
        numgrd(j,i) = step2 * (er - el)
      end do
    end do

    write (*,*) 'Numerical Gradient:'
    do i = 1,mol%nat
      write (*,'(3f18.8)') numgrd(1:3,i)
    end do

    deallocate (numgrd,grad)

    return
  end subroutine numgrad

!==========================================================================================!
end module calc_module
