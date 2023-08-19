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

!> module api_engrad
!> a collection of engrad calls for different APIs
!> this builds the communication between CRESTs
!> "calculation_settings" and the respective API setups

!=========================================================================================!
module api_engrad

  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use strucrd
  use calc_type
  use iomod,only:makedir,directory_exist,remove
  !> API modules
  use api_helpers
  use tblite_api
  use gfn0_api
  use gfnff_api
  use xhcff_api
  implicit none
  !--- private module variables and parameters
  private

  public :: tblite_engrad
  public :: gfn0_engrad,gfn0occ_engrad
  public :: gfnff_engrad
  public :: xhcff_engrad

!=========================================================================================!
!=========================================================================================!
contains    !> MODULE PROCEDURES START HERE
!=========================================================================================!
!=========================================================================================!

  subroutine tblite_engrad(mol,calc,energy,grad,iostatus)
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus

    character(len=:),allocatable :: cpath
    logical :: loadnew,pr

    integer :: i,j,k,l,ich,och,io
    logical :: ex
    iostatus = 0
    pr = .false.

!>--- setup system call information
    !$omp critical
    call tblite_init(calc,loadnew)
!>--- tblite printout handling
    call api_handle_output(calc,'tblite.out',mol,pr)
    if(pr)then
       !> tblite uses its context (ctx)( type, rather than calc%prch
       calc%ctx%unit = calc%prch
       calc%ctx%verbosity = 1 
    else
       calc%ctx%verbosity = 0
    endif

!>-- populate parameters and wavefunction
    if (loadnew) then
      call tblite_setup(mol,calc%chrg,calc%uhf,calc%tblitelvl,calc%etemp, &
      &    calc%ctx,calc%wfn,calc%tbcalc)
      call tblite_addsettings(calc%tbcalc,calc%maxscc,calc%rdwbo,calc%saveint)
      call tblite_add_solv(mol,calc%chrg,calc%uhf,calc%ctx,calc%wfn,calc%tbcalc, &
      &    calc%solvmodel,calc%solvent)
    end if
    !$omp end critical

!>--- do the engrad call
    call initsignal()
    call tblite_singlepoint(mol,calc%chrg,calc%uhf,calc%accuracy, &
    & calc%ctx,calc%wfn,calc%tbcalc,energy,grad,calc%tbres,iostatus)
    if (iostatus /= 0) return
    call api_print_e_grd(pr,calc%ctx%unit,mol,energy,grad)

!>--- postprocessing, getting other data
    !$omp critical
    call tblite_wbos(calc,mol,iostatus)
    !$omp end critical

    return
  end subroutine tblite_engrad

!========================================================================================!

  subroutine gfn0_engrad(mol,calc,g0calc,energy,grad,iostatus)
!> This is the GFN0 engrad call that uses the standard implementation
    implicit none
    !> INPUT
    type(coord) :: mol
    type(calculation_settings) :: calc
    type(gfn0_data),intent(inout),allocatable  :: g0calc
    !> OUTPUT
    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus
    !> LOCAL
    type(gfn0_results) :: res
    character(len=:),allocatable :: cpath
    logical :: loadnew
    logical :: pr

    integer :: i,j,k,l,ich,och,io
    logical :: ex
    iostatus = 0
    pr = .false.
!>--- setup system call information
    !$omp critical
    call gfn0_init(calc,g0calc,loadnew)
!>--- printout handling
    call api_handle_output(calc,'gfn0.out',mol,pr)
!>-- populate parameters and wavefunction
    if (loadnew) then
      call gfn0_setup(mol,calc%chrg,calc%uhf,g0calc)
      call gfn0_init2(mol,calc,g0calc)
    end if
    call gfn0_init3(mol,calc,g0calc)
    !$omp end critical

!>--- do the engrad call
    call initsignal()
    call gfn0_sp(mol,calc%chrg,calc%uhf,g0calc,energy,grad,iostatus,res)
    if (iostatus /= 0) return
    if (pr) then
      call gfn0_print(calc%prch,g0calc,res)
      call api_print_e_grd(pr,calc%prch,mol,energy,grad)
    end if

!>--- postprocessing, getting other data
    !$omp critical
    call gfn0_wbos(calc,calc%g0calc,mol,iostatus)
    !$omp end critical

    return
  end subroutine gfn0_engrad

!========================================================================================!

  subroutine gfn0occ_engrad(mol,calc,g0calc,energy,grad,iostatus)
!> This is the GFN0 engrad call in which a config can be specified
    implicit none
    !> INPUT
    type(coord) :: mol
    type(calculation_settings) :: calc
    type(gfn0_data),intent(inout),allocatable  :: g0calc
    !> OUTPUT
    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(:,:)
    integer,intent(out) :: iostatus
    !> LOCAL
    type(gfn0_results) :: res
    character(len=:),allocatable :: cpath
    logical :: loadnew,pr
    integer :: i,j,k,l,ich,och,io
    logical :: ex
    iostatus = 0
    pr = .false.
!>--- setup system call information
    !$omp critical
    call gfn0occ_init(calc,g0calc,loadnew)
!>--- printout handling
    call api_handle_output(calc,'gfn0.out',mol,pr)
!>--- populate parameters and wavefunction
    if (loadnew) then
      call gfn0_setup(mol,calc%chrg,calc%uhf,g0calc)
      call gfn0occ_init2(mol,calc,g0calc)
    end if
    call gfn0occ_init3(mol,calc,g0calc)
    !$omp end critical

!>--- do the engrad call
    call initsignal()
    call gfn0_sp_occ(mol,calc%chrg,calc%uhf,calc%occ,g0calc, &
    &    energy,grad,iostatus,res)
    if (iostatus /= 0) return
    if (pr) then
      call gfn0_print(calc%prch,g0calc,res)
      call api_print_e_grd(pr,calc%prch,mol,energy,grad)
    end if

!>--- postprocessing, getting other data
    !$omp critical
    call gfn0_wbos(calc,g0calc,mol,iostatus)
    !$omp end critical

    return
  end subroutine gfn0occ_engrad

!========================================================================================!
  subroutine gfnff_engrad(mol,calc,energy,grad,iostatus)
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus

    character(len=:),allocatable :: cpath
    logical :: loadnew,pr
    integer :: i,j,k,l,ich,och,io
    logical :: ex
    iostatus = 0
    pr = .false.
!>--- setup system call information
    !$omp critical
    call gfnff_init(calc,loadnew)
!>--- printout handling
    call api_handle_output(calc,'gfnff.out',mol,pr) 

!>--- populate parameters and neighbourlists
    if (loadnew) then
      call gfnff_api_setup(mol,calc%chrg,calc%ff_dat,iostatus,pr,calc%prch)
    end if
    !$omp end critical
    if (iostatus /= 0) return

!>--- do the engrad call
    call initsignal()
    call gfnff_sp(mol,calc%ff_dat,energy,grad,iostatus)
    if (iostatus /= 0) return

!>--- printout
    if (pr) then
      call gfnff_printout(calc%prch,calc%ff_dat)
      call api_print_e_grd(pr,calc%prch,mol,energy,grad)
    end if

!>--- postprocessing, getting other data
    !$omp critical
    call gfnff_wbos(calc,mol,iostatus)
    !$omp end critical

    return
  end subroutine gfnff_engrad

!========================================================================================!
  subroutine xhcff_engrad(mol,calc,energy,grad,iostatus)
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus

    character(len=:),allocatable :: cpath
    logical :: loadnew,pr
    integer :: i,j,k,l,ich,och,io
    logical :: ex
    iostatus = 0
    pr = .false.
!>--- setup system call information
    !$omp critical
    call xhcff_initcheck(calc,loadnew)
!>--- printout handling
    call api_handle_output(calc,'xhcff.out',mol,pr)
!>--- populate parameters
    if (loadnew) then
      !> call xhcff with verbosity turned off
      call xhcff_setup(mol,calc%xhcff,calc%extpressure,calc%ngrid,calc%proberad, &
      &                calc%vdwset,pr,calc%prch,iostatus)
    end if
    !$omp end critical
    if (iostatus /= 0) return

!>--- do the engrad call
    call initsignal()
    call xhcff_sp(mol,calc%xhcff,energy,grad,iostatus)
    if (iostatus /= 0) return

!>--- printout
    if (pr) then
      !> the xhcff_sp call includes the printout within xhcff-lib
      call api_print_e_grd(pr,calc%prch,mol,energy,grad)
    end if

!>--- postprocessing, getting other data

    return
  end subroutine xhcff_engrad

!========================================================================================!
end module api_engrad
