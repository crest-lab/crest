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
  use constraints
  use calc_type
  use xtb_sc
  use lj
  use nonadiabatic_module
  implicit none

  !=========================================================================================!
  !--- private module variables and parameters
  private
  !integer :: i,j,k,l,ich,och,io
  !logical :: ex

  !--- some constants and name mappings
  real(wp),parameter :: bohr = 0.52917726_wp
  real(wp),parameter :: autokcal = 627.509541_wp

  public :: engrad
  interface engrad
    module procedure :: engrad_mol
    module procedure :: engrad_xyz
  end interface engrad

  public :: test_engrad
  public :: numhess
  public :: constrhess

contains
!========================================================================================!
! subroutine test_engrad
! test engrad routine
  subroutine test_engrad(fname)
    implicit none

    character(len=*) :: fname

    type(coord) :: mol
    type(coord) :: molo
    type(calcdata) :: calc
    type(constraint) :: constr
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    type(calculation_settings) :: job
    integer :: i,j,k,l,ich,och,io

    logical :: pr

    pr = .true.

    call mol%open(fname)

    allocate (grad(3,mol%nat),source=0.0_wp)

    !job%id = 99
    !job%other = '15.0   1.5'
    !job%id = 10
    job%id = 0
    job%calcspace = 'testdir'
    call calc%add(job)

    !call constr%bondconstraint(1,2,2.0_wp,0.1_wp)
    !call constr%print()
    !call constr%sphereconstraint(mol%nat,4.0_wp,1.0_wp,6.0_wp,.true.)
    !call constr%angleconstraint(2,1,3,1387.0_wp,0.02_wp)
    call constr%dihedralconstraint(1,2,3,9,27.0_wp,0.02_wp)
    !call constr%print()
    call calc%add(constr)
    call calc%cons(1)%print()

    call engrad(mol,calc,energy,grad,io)

    write (*,*) 'Energy: ',energy
    write (*,*) 'Gradient:'
    do i = 1,mol%nat
      write (*,'(3f18.8)') grad(1:3,i)
    end do

    write (*,*)
    molo = mol
    if (pr) then
      call numgrad(molo,calc,grad)
    end if
    deallocate (grad)

    return
  end subroutine test_engrad

!========================================================================================!
! subroutine engrad
! main routine to perform some energy and gradient calculation
  subroutine engrad_mol(mol,calc,energy,gradient,iostatus)
    implicit none
    type(coord) :: mol
    type(calcdata) :: calc
    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: gradient(3,mol%nat)
    integer,intent(out) :: iostatus

    integer :: i,j,k,l,n,io

    real(wp) :: dum1,dum2
    real(wp) :: efix
    real(wp),allocatable :: grdfix(:,:)
    real(wp),allocatable :: grdtmp(:,:,:)
    real(wp),allocatable :: etmp(:)

    call initsignal()

    iostatus = 0
    dum1 = 1.0_wp
    dum2 = 1.0_wp
    etmp = 0.0_wp

    !>--- Calculation
    n = calc%ncalculations
    if (n > 0) then
      allocate (grdtmp(3,mol%nat,n),etmp(n),source=0.0_wp)
      !>--- loop over all calculations to be done
      do i = 1,calc%ncalculations
        if (calc%which > 0) then
          if (i < calc%which) cycle
          if (i > calc%which) exit
        end if
        select case (calc%calcs(i)%id)
        case (10) !-- xtb system call
          call xtb_engrad(mol,calc%calcs(i),etmp(i),grdtmp(:,:,i),iostatus)
        case (99) !-- Lennard-Jones dummy calculation
          if (allocated(calc%calcs(i)%other)) then
            read (calc%calcs(i)%other,*) dum1,dum2
          end if
          call lj_engrad(mol%nat,mol%xyz,dum1,dum2,etmp(i),grdtmp(:,:,i))
        case default
          !write (*,*) 'Nothing selected for energy and gradient calculation.'
          etmp(i) = 0.0_wp
          grdtmp(:,:,i) = 0.0_wp
        end select
      end do
      !>--- switch case for what to to with the energies
      select case (calc%id)
      case default !> take e+grd only from first level
        energy = etmp(1)
        gradient = grdtmp(:,:,1)
      case (2:)
        j = calc%id
        if (j <= calc%ncalculations) then
          energy = etmp(j)
          gradient = grdtmp(:,:,j)
        end if
      case( -1 ) !> non-adiabatic arithmetic mean
        if(calc%ncalculations > 1)then
          call engrad_mean(mol%nat,etmp(1),etmp(2),grdtmp(:,:,1), &
          &                grdtmp(:,:,2),energy,gradient)
        endif
      end select
      !>--- printout (to file or stdout)
      call calc_print_energies(calc,energy,etmp)
      !>--- deallocate
      deallocate (etmp,grdtmp)
    end if

    !>--- Constraints
    if (calc%nconstraints > 0) then
      allocate (grdfix(3,mol%nat),source=0.0_wp)
      do i = 1,calc%nconstraints
        efix = 0.0_wp
        grdfix = 0.0_wp
        call calc_constraint(mol%nat,mol%xyz,calc%cons(i),efix,grdfix)
        energy = energy + efix
        gradient = gradient + grdfix
      end do
      deallocate (grdfix)
    end if

    return
  end subroutine engrad_mol

  subroutine engrad_xyz(n,xyz,at,calc,energy,gradient,iostatus)
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n) ! coordinates should be in Bohr (a.u.)
    integer,intent(in) :: at(n)
    type(calcdata) :: calc
    type(coord) :: mol

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: gradient(3,n)
    integer,intent(out) :: iostatus

    real(wp) :: dum1,dum2

    mol%nat = n
    mol%at = at
    mol%xyz = xyz
    call engrad_mol(mol,calc,energy,gradient,iostatus)

    return
  end subroutine engrad_xyz

!========================================================================================!
! subroutine numgrad
! routine to perform a numerical gradient calculation
  subroutine numgrad(mol,calc,angrad)
    implicit none

    type(coord) :: mol
    type(calcdata) :: calc
    real(wp) :: angrad(3,mol%nat)
    integer :: i,j,k,l,ich,och,io

    real(wp) :: energy,el,er
    real(wp),allocatable :: grad(:,:)
    real(wp),allocatable :: numgrd(:,:)
    real(wp),parameter :: step = 0.00001_wp,step2 = 0.5_wp / step

    allocate (grad(3,mol%nat),source=0.0_wp)
    allocate (numgrd(3,mol%nat),source=0.0_wp)

    do i = 1,mol%nat
      do j = 1,3
        mol%xyz(j,i) = mol%xyz(j,i) + step
        call engrad(mol%nat,mol%xyz,mol%at,calc,er,grad,io)

        mol%xyz(j,i) = mol%xyz(j,i) - 2 * step
        call engrad(mol%nat,mol%xyz,mol%at,calc,el,grad,io)

        mol%xyz(j,i) = mol%xyz(j,i) + step
        numgrd(j,i) = step2 * (er - el)
      end do
    end do

    write (*,*) 'Numerical Gradient:'
    do i = 1,mol%nat
      write (*,'(3f18.8)') numgrd(1:3,i)
    end do

    write (*,*)
    write (*,*) 'Gradient Difference:'
    do i = 1,mol%nat
      write (*,'(3f18.8)') numgrd(1:3,i) - angrad(1:3,i)
    end do

    deallocate (numgrd,grad)

    return
  end subroutine numgrad

!========================================================================================!
!> subroutine numhess
!> routine to perform a numerical hessian calculation
  subroutine numhess(nat,at,xyz,calc,hess)
    implicit none

    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(calcdata) :: calc
    real(wp),intent(out) :: hess(nat * 3,nat * 3)

    type(coord) :: mol !> coord type, so that the input remains unchanged
    real(wp) :: energy,el,er
    real(wp),allocatable :: gradr(:,:),gradl(:,:)
    real(wp),parameter :: step = 0.00001_wp,step2 = 0.5_wp / step
    integer :: i,j,k,l,ii,jj,io

    hess = 0.0_wp
    mol%nat = nat
    mol%at = at
    mol%xyz = xyz
    allocate (gradr(3,mol%nat),source=0.0_wp)
    allocate (gradl(3,mol%nat),source=0.0_wp)

    do i = 1,mol%nat
      do j = 1,3
        ii = (i - 1) * 3 + j
        gradr = 0.0_wp
        mol%xyz(j,i) = mol%xyz(j,i) + step
        call engrad(mol%nat,mol%xyz,mol%at,calc,er,gradr,io)

        gradl = 0.0_wp
        mol%xyz(j,i) = mol%xyz(j,i) - 2.0_wp * step
        call engrad(mol%nat,mol%xyz,mol%at,calc,el,gradl,io)

        mol%xyz(j,i) = mol%xyz(j,i) + step
        do k = 1,mol%nat
          do l = 1,3
            jj = (k - 1) * 3 + l
            hess(jj,ii) = (gradr(l,k) - gradl(l,k)) * step2
          end do
        end do
      end do
    end do

    deallocate (gradl,gradr)
    call mol%deallocate()
    return
  end subroutine numhess

!========================================================================================!
!> subroutine constrhess
!> routine to perform a numerical Hessian calculation
!> but ONLY include contributions of the constraints.
!>
!> phess is the packed Hessian on which the constraint
!> contributions are added.
!>--------------------------------------------------------
  subroutine constrhess(nat,at,xyz,calc,phess)
    implicit none

    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(calcdata),intent(in) :: calc
    real(wp),intent(inout) :: phess((nat * 3) * ((nat * 3) + 1) / 2)

    type(calcdata) :: dummycalc
    integer :: tc_backup
    real(wp) :: energy,el,er
    real(wp),allocatable :: hess(:,:)

    integer :: i,j,k,n3

    !phess = 0.0_wp
    if (calc%nconstraints <= 0) return
    dummycalc = calc !> new dummy calculation
    dummycalc%id = 0  !> set to zero so that only constraints are considered
    dummycalc%pr_energies = .false.
    n3 = nat * 3
    allocate (hess(n3,n3),source=0.0_wp)

    call numhess(nat,at,xyz,dummycalc,hess)

    k = 0
    do i = 1,n3
      do j = 1,i
        k = k + 1
        phess(k) = phess(k) + 0.5_wp * (hess(j,i) + hess(i,j))
      end do
    end do

    deallocate (hess)
    return
  end subroutine constrhess
!==========================================================================================!

  subroutine calc_print_energies(calc,energy,energies)
    implicit none
    type(calcdata) :: calc
    real(wp) :: energy
    real(wp) :: energies(calc%ncalculations)
    integer :: i,j
    character(len=20) :: atmp
    character(len=:),allocatable :: btmp

    if (.not. calc%pr_energies) return
    btmp = ''
    write (atmp,'(f20.12)') energy
    btmp = btmp//atmp
    do i = 1,calc%ncalculations
      write (atmp,'(f20.12)') energies(i)
      btmp = btmp//atmp
    end do
    write (calc%eout_unit,'(a)') btmp
    deallocate (btmp)
    return
  end subroutine calc_print_energies

!==========================================================================================!
end module calc_module
