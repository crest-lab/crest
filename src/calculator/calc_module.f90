!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2023 Philipp Pracht
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
!>--- types and readers
  use iso_fortran_env,only:wp => real64
  use strucrd
  use calc_type
!>--- potentials and API's
  use xtb_sc
  use generic_sc
  use lj
  use api_engrad
!>--- other
  use constraints
  use nonadiabatic_module
  implicit none
!=========================================================================================!
!>--- private module variables and parameters
  private

!>--- some constants and name mappings
  real(wp),parameter :: bohr = 0.52917726_wp
  real(wp),parameter :: autokcal = 627.509541_wp

!>--- public module routines
  public :: engrad
  interface engrad
    module procedure :: engrad_mol
    module procedure :: engrad_xyz
  end interface engrad

  public :: calc_eprint
  interface calc_eprint
    module procedure :: calc_print_energies
    module procedure :: calc_print_energies2
  end interface calc_eprint

  public :: numhess1,numhess2
  public :: constrhess

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

!========================================================================================!
!> subroutine engrad
!> main routine to perform some energy and gradient calculation
  subroutine engrad_mol(mol,calc,energy,gradient,iostatus)
    implicit none
    type(coord) :: mol
    type(calcdata) :: calc
    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: gradient(3,mol%nat)
    integer,intent(out) :: iostatus

    integer :: i,j,k,l,n,io,nocc

    real(wp) :: dum1,dum2
    real(wp) :: efix
    real(wp),allocatable :: grdfix(:,:)
    real(wp),allocatable :: grdtmp(:,:,:)
    real(wp),allocatable :: etmp(:)

!==========================================================!
    call initsignal()

    !>--- Calculation setup
    n = calc%ncalculations

    !$omp critical
    if (n > 0) then
      if (.not. allocated(calc%etmp)) allocate (calc%etmp(n),source=0.0_wp)
      if (.not. allocated(calc%grdtmp)) allocate (calc%grdtmp(3,mol%nat,n),source=0.0_wp)
      if (.not. allocated(calc%eweight)) then
        allocate (calc%eweight(n),source=0.0_wp)
        do i = 1,n
          calc%eweight(i) = calc%calcs(i)%weight
        end do
      end if
    end if
    !$omp end critical

    iostatus = 0
    dum1 = 1.0_wp
    dum2 = 1.0_wp
    calc%etmp = 0.0_wp
    !calc%grdtmp = 0.0_wp

!==========================================================!
    !>--- Calculation
    if (n > 0) then

      !==================================================================================!
      !>--- loop over all calculations to be done
      do i = 1,calc%ncalculations
        select case (calc%calcs(i)%id)
        case (jobtype%xtbsys)  !>-- xtb system call
          call xtb_engrad(mol,calc%calcs(i),calc%etmp(i),calc%grdtmp(:,:,i),iostatus)

        case (jobtype%generic) !>-- generic script/program call
          call generic_engrad(mol,calc%calcs(i),calc%etmp(i),calc%grdtmp(:,:,i),iostatus)

        case (jobtype%tblite)  !>-- tblite api call
          call tblite_engrad(mol,calc%calcs(i),calc%etmp(i),calc%grdtmp(:,:,i),iostatus)

        case (jobtype%gfn0) !>-- GFN0-xTB api
          call gfn0_engrad(mol,calc%calcs(i),calc%calcs(i)%g0calc,calc%etmp(i), &
          &                calc%grdtmp(:,:,i),iostatus)

        case (jobtype%gfn0occ) !>-- Special GFN0-xTB api given orbital population
          !> note the use of calc%g0calc instead of calc%calcs(i)%g0calc !
          call gfn0occ_engrad(mol,calc%calcs(i),calc%g0calc,calc%etmp(i), &
          &                   calc%grdtmp(:,:,i),iostatus)
 
        case (jobtype%gfnff) !>-- GFN-FF api
          call gfnff_engrad(mol,calc%calcs(i),calc%etmp(i),calc%grdtmp(:,:,i),iostatus)


        case (99) !-- Lennard-Jones dummy calculation
          if (allocated(calc%calcs(i)%other)) then
            read (calc%calcs(i)%other,*) dum1,dum2
          end if
          call lj_engrad(mol%nat,mol%xyz,dum1,dum2,calc%etmp(i),calc%grdtmp(:,:,i))
        case default
          !write (*,*) 'Nothing selected for energy and gradient calculation.'
          calc%etmp(i) = 0.0_wp
          calc%grdtmp(:,:,i) = 0.0_wp
        end select
        if (iostatus /= 0) then
          return
        end if
      end do

      !==================================================================================!
      !>--- switch case for what to to with the energies
      select case (calc%id)
      case default
        !> take e+grd only from first level
        energy = calc%etmp(1)
        gradient = calc%grdtmp(:,:,1)
      case (2:)
        !> take e+grd from one of the specified calculations
        j = calc%id
        if (j <= calc%ncalculations) then
          energy = calc%etmp(j)
          gradient = calc%grdtmp(:,:,j)
        end if
      case (-1)
        !> take non-adiabatic arithmetic mean (of first two calculations)
        call engrad_mean(mol%nat,calc%ncalculations,calc%etmp,calc%grdtmp, &
        &                energy,gradient)
      end select
      !>--- printout (to file or stdout)
      call calc_eprint(calc,energy,calc%etmp)
    end if

!==========================================================!
    !>--- Constraints
    if (calc%nconstraints > 0) then
      !$omp critical
      allocate (grdfix(3,mol%nat),source=0.0_wp)
      !$omp end critical
      do i = 1,calc%nconstraints
        efix = 0.0_wp
        grdfix = 0.0_wp
        if (calc%cons(i)%type > 0) then
          !>--- structural constraints
          call calc_constraint(mol%nat,mol%xyz,calc%cons(i),efix,grdfix)
        else if (allocated(calc%etmp) .and. allocated(calc%grdtmp)) then
          !>--- non-adiabatic constraints
          if (n > 1) then
            call calc_nonadiabatic_constraint(mol%nat,calc%cons(i),n,calc%etmp,calc%grdtmp,efix,grdfix)
          else !> this "else" is necessary for constrained model hessians
            efix = 0.0_wp
            grdfix = 0.0_wp
          end if
        end if
        energy = energy + efix
        gradient = gradient + grdfix
      end do
      !$omp critical
      deallocate (grdfix)
      !$omp end critical
    end if

    return
  end subroutine engrad_mol

  !>--- wrapper for the engrad_mol routine
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
!========================================================================================!
!========================================================================================!

  subroutine numgrad(mol,calc,angrad)
!*******************************************************
!* subroutine numgrad
!* routine to perform a numerical gradient calculation
!*******************************************************
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
  subroutine numhess1(nat,at,xyz,calc,hess,io)
!***************************************************
!* Calculate and return the Hessian for the
!* COMBINED energy and gradient from engrad,
!* including constraints
!***************************************************
    implicit none

    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(calcdata) :: calc
    real(wp),intent(out) :: hess(nat * 3, nat*3)
    integer,intent(out)  :: io

    type(coord) :: mol !> coord type, so that the input remains unchanged
    real(wp) :: energy,el,er,hij
    real(wp),allocatable :: gradr(:,:),gradl(:,:)
    real(wp),allocatable :: gradr_tmp(:,:,:), gradl_tmp(:,:,:)
    real(wp),parameter :: step = 0.005_wp,step2 = 0.5_wp / step !0.00001_wp
    integer :: i,j,k,l,m,ii,jj

    hess = 0.0_wp
    io = 0
    mol%nat = nat
    mol%at = at
    mol%xyz = xyz
    
    allocate (gradr(3,mol%nat),source=0.0_wp) !dummy
    allocate (gradl(3,mol%nat),source=0.0_wp) !dummy

    do i = 1,mol%nat
      do j = 1,3
        ii = (i - 1) * 3 + j
        gradr = 0.0_wp
        mol%xyz(j,i) = mol%xyz(j,i) + step
        call engrad(mol, calc,er,gradr,io)

        gradl = 0.0_wp
        mol%xyz(j,i) = mol%xyz(j,i) - 2.0_wp * step
        call engrad(mol, calc,el,gradl,io)

        mol%xyz(j,i) = mol%xyz(j,i) + step

          do k = 1,mol%nat
            do l = 1,3
              jj = (k - 1) * 3 + l
              hess(jj,ii) = (gradr(l,k) - gradl(l,k)) * step2
            end do
          end do

      end do
    end do

    !Symmetrize Hessian
      do i = 1,nat*3
        do j = i,nat*3
          hij = (hess(i,j) + hess(j,i))*0.5_wp 
          hess(i,j) = hij
          hess(j,i) = hij
        end do
      end do

    call engrad(mol,calc,el,gradl,io) !>- to get the gradient of the non-displaced structure

    deallocate (gradl_tmp,gradr_tmp)
    deallocate (gradl,gradr)
    call mol%deallocate()
    return
  end subroutine numhess1

!========================================================================================!
  subroutine numhess2(nat,at,xyz,calc,hess,io)
!***************************************************
!* Calculate and return the Hessian for EACH of the
!* calculation levels, but contributions from
!* constraints are ignored.
!***************************************************
    implicit none

    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(calcdata) :: calc
    real(wp),intent(out) :: hess(nat * 3,nat * 3,calc%ncalculations)
    integer,intent(out)  :: io

    type(coord) :: mol !> coord type, so that the input remains unchanged
    real(wp) :: energy,el,er,hij
    real(wp),allocatable :: gradr(:,:),gradl(:,:)
    real(wp),allocatable :: gradr_tmp(:,:,:), gradl_tmp(:,:,:)
    real(wp),parameter :: step = 0.005_wp,step2 = 0.5_wp / step !0.00001_wp
    integer :: i,j,k,l,m,ii,jj

    hess = 0.0_wp
    io = 0
    mol%nat = nat
    mol%at = at
    mol%xyz = xyz
    
    allocate (gradr(3,mol%nat),source=0.0_wp) !dummy
    allocate (gradl(3,mol%nat),source=0.0_wp) !dummy

    allocate (gradr_tmp(3,mol%nat,calc%ncalculations),source=0.0_wp)
    allocate (gradl_tmp(3,mol%nat,calc%ncalculations),source=0.0_wp)

    do i = 1,mol%nat
      do j = 1,3
        ii = (i - 1) * 3 + j
        !gradr = 0.0_wp
        mol%xyz(j,i) = mol%xyz(j,i) + step
        call engrad(mol%nat,mol%xyz,mol%at,calc,er,gradr,io)

        gradr_tmp = calc%grdtmp

        !gradl = 0.0_wp
        mol%xyz(j,i) = mol%xyz(j,i) - 2.0_wp * step
        call engrad(mol%nat,mol%xyz,mol%at,calc,el,gradl,io)

        gradl_tmp = calc%grdtmp

        mol%xyz(j,i) = mol%xyz(j,i) + step

        do m = 1,calc%ncalculations
          do k = 1,mol%nat
            do l = 1,3
              jj = (k - 1) * 3 + l
              hess(jj,ii,m) = (gradr_tmp(l,k,m) - gradl_tmp(l,k,m)) * step2
            end do
          end do
        end do
      end do
    end do

    !Symmetrize Hessian
    do m = 1,calc%ncalculations
      do i = 1,nat*3
        do j = i,nat*3
          hij = (hess(i,j,m) + hess(j,i,m))*0.5_wp 
          hess(i,j,m) = hij
          hess(j,i,m) = hij
        end do
      end do
    end do

    call engrad(mol%nat,mol%xyz,mol%at,calc,el,gradl,io) !>- to get the gradient of the non-displaced s

    deallocate (gradl_tmp,gradr_tmp)
    deallocate (gradl,gradr)
    call mol%deallocate()
    return
  end subroutine numhess2


!========================================================================================!
  subroutine constrhess(nat,at,xyz,calc,phess)
!*********************************************************
!* subroutine constrhess
!* routine to perform a numerical Hessian calculation
!* but ONLY include contributions of the constraints.
!*
!* phess is the packed Hessian on which the constraint
!* contributions are added.
!*********************************************************
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
    logical :: consgeo
    integer :: i,j,k,n3,io

    if (calc%nconstraints <= 0) return
    !>--- skip if only nonadiabatic constraints

    dummycalc = calc !> new dummy calculation
    dummycalc%id = 0  !> set to zero so that only constraints are considered
    dummycalc%ncalculations = 0
    dummycalc%pr_energies = .false.
    n3 = nat * 3
    allocate (hess(n3,n3),source=0.0_wp)

    call numhess1(nat,at,xyz,dummycalc,hess,io)

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
    real(wp) :: energies(:)
    integer :: i,j,l
    character(len=20) :: atmp
    character(len=:),allocatable :: btmp

    if (.not. calc%pr_energies) return
    btmp = ''
    write (atmp,'(f20.12)') energy
    btmp = btmp//atmp
    l = size(energies,1)
    do i = 1,l
      write (atmp,'(f20.12)') energies(i)
      btmp = btmp//atmp
    end do
    write (calc%eout_unit,'(a)') btmp
    deallocate (btmp)
    return
  end subroutine calc_print_energies

!==========================================================================================!
  subroutine calc_print_energies2(calc,energy,energies,chnl)
    implicit none
    type(calcdata) :: calc
    real(wp) :: energy
    real(wp) :: energies(:)
    integer :: chnl
    integer :: i,j,l
    character(len=20) :: atmp
    character(len=:),allocatable :: btmp
    btmp = ''
    write (atmp,'(f20.12)') energy
    btmp = btmp//atmp
    l = size(energies,1)
    do i = 1,l
      write (atmp,'(f20.12)') energies(i)
      btmp = btmp//atmp
    end do
    write (chnl,'(a)') btmp
    deallocate (btmp)
    return
  end subroutine calc_print_energies2

!==========================================================================================!
!==========================================================================================!
end module calc_module
