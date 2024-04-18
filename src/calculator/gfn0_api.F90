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

!====================================================!
! module gfn0_api
! An interface to gfn0-xtb standalone calculations
!====================================================!

module gfn0_api
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use strucrd
#ifdef WITH_GFN0
  use gfn0_interface
  use gfn0_module,only:gfn0_gbsa_init,generate_config
#endif
  use wiberg_mayer,only:get_wbo,get_wbo_rhf,density_matrix
  implicit none
  private

#ifndef WITH_GFN0
  !> these are placeholders if no gfn0 module is used!
  type :: gfn0_results
    integer :: id = 0
  end type gfn0_results
  type :: gfn0_data
    integer :: id = 0
  end type gfn0_data
#endif

  public :: gfn0_results,gfn0_data
  public :: gfn0_setup,gfn0_addsettings
  public :: gfn0_sp,gfn0_sp_occ
  public :: gfn0_getwbos
  public :: gfn0_gen_occ
  public :: gfn0_print
  public :: gfn0_getdipole

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine gfn0_setup(mol,chrg,uhf,g0calc)
    implicit none
    type(coord),intent(in)  :: mol
    integer,intent(in)      :: chrg
    integer,intent(in)      :: uhf
    type(gfn0_data),intent(inout) :: g0calc
#ifdef WITH_GFN0

    !> initialize parametrization of GFN0
    call gfn0_init(mol%nat,mol%at,mol%xyz,chrg,uhf,g0calc)

#else /* WITH_GFN0 */
    write (stdout,*) 'Error: Compiled without GFN0-xTB support!'
    write (stdout,*) 'Use -DWITH_GFN0=true in the setup to enable this function'
    error stop
#endif
  end subroutine gfn0_setup

!========================================================================================!
!> gfn0_addsettings is used to add other settings from
!> CRESTs calculation object to the gfn0_data
  subroutine gfn0_addsettings(mol,g0calc,solv,model,etemp,loadwbo)
    implicit none
    type(coord),intent(in)  :: mol
    type(gfn0_data),intent(inout) :: g0calc
    character(len=*),intent(in),optional :: solv
    character(len=*),intent(in),optional :: model
    real(wp),intent(in),optional :: etemp
    logical,intent(in),optional  :: loadwbo

    integer :: nao
#ifdef WITH_GFN0
    !> add solvation?
    if (present(solv)) then
      if (allocated(g0calc%gbsa)) deallocate (g0calc%gbsa)
      allocate (g0calc%gbsa)
      if (present(model)) then
        select case (model)
        case ('alpb')
          call gfn0_gbsa_init(mol%nat,mol%at,.true.,solv,g0calc%gbsa)
        case default !> default GBSA
          call gfn0_gbsa_init(mol%nat,mol%at,.false.,solv,g0calc%gbsa)
        end select
      else
        call gfn0_gbsa_init(mol%nat,mol%at,.false.,solv,g0calc%gbsa)
      end if
    end if
    if (present(etemp)) then
      g0calc%xtbData%etemp = max(0.0_wp,etemp)
    end if
    if (present(loadwbo)) then
      if (loadwbo) then
        nao = g0calc%basis%nao
        if (.not.allocated(g0calc%wfn%S)) allocate (g0calc%wfn%S(nao,nao),source=0.0_wp)
      end if
    end if
#else
    write (stdout,*) 'Error: Compiled without GFN0-xTB support!'
    write (stdout,*) 'Use -DWITH_GFN0=true in the setup to enable this function'
    error stop
#endif
  end subroutine gfn0_addsettings

!========================================================================================!

  subroutine gfn0_sp(mol,chrg,uhf,g0calc,energy,gradient,iostatus,res)
    implicit none
    !> INPUT
    type(coord),intent(in)  :: mol
    integer,intent(in) :: chrg
    integer,intent(in) :: uhf
    type(gfn0_data),intent(inout) :: g0calc
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,mol%nat)
    integer,intent(out) :: iostatus
    type(gfn0_results),intent(inout),optional :: res
    !> LOCAL
    logical :: fail
    energy = 0.0_wp
    gradient = 0.0_wp
    iostatus = 0
    fail = .false.
#ifdef WITH_GFN0
    !if (present(res)) then
    call gfn0_singlepoint(mol%nat,mol%at,mol%xyz,chrg,uhf,g0calc, &
    &          energy,gradient,fail,res=res)
    ! else
    !   call gfn0_singlepoint(mol%nat,mol%at,mol%xyz,chrg,uhf,g0calc, &
    !   &          energy,gradient,fail)
    ! end if
    if (fail) then
      iostatus = -1
    end if
#else
    write (stdout,*) 'Error: Compiled without GFN0-xTB support!'
    write (stdout,*) 'Use -DWITH_GFN0=true in the setup to enable this function'
    error stop
#endif
  end subroutine gfn0_sp

!========================================================================================!

  subroutine gfn0_sp_occ(mol,chrg,uhf,occ,g0calc,energy,gradient,iostatus,res)
    implicit none
    !> INPUT
    type(coord),intent(in)  :: mol
    integer,intent(in) :: chrg
    integer,intent(in) :: uhf
    type(gfn0_data),intent(inout) :: g0calc
    real(wp),intent(in) :: occ(:)
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,mol%nat)
    integer,intent(out) :: iostatus
    type(gfn0_results),intent(inout),optional :: res
    !> LOCAL
    logical :: fail
    energy = 0.0_wp
    gradient = 0.0_wp
    iostatus = 0
    fail = .false.
#ifdef WITH_GFN0
    !if (present(res)) then
    call gfn0_occ_singlepoint(mol%nat,mol%at,mol%xyz,chrg,uhf,occ,g0calc, &
    &          energy,gradient,fail,res=res)
    !else
    !  call gfn0_occ_singlepoint(mol%nat,mol%at,mol%xyz,chrg,uhf,occ,g0calc, &
    !  &          energy,gradient,fail)
    !end if
    if (fail) then
      iostatus = -1
    end if
#else
    write (stdout,*) 'Error: Compiled without GFN0-xTB support!'
    write (stdout,*) 'Use -DWITH_GFN0=true in the setup to enable this function'
    error stop
#endif
  end subroutine gfn0_sp_occ

!========================================================================================!

  subroutine gfn0_gen_occ(nel,nao,active,occ)
    implicit none
    integer,intent(in) :: nel
    integer,intent(in) :: nao
    integer,intent(in) :: active(:)
    real(wp),intent(out) :: occ(nao)
    integer :: i
    occ = 0.0_wp
#ifdef WITH_GFN0
    call generate_config(nel,nao,occ,active)
#endif
  end subroutine gfn0_gen_occ

!========================================================================================!

  subroutine gfn0_print(iunit,g0calc,res)
    implicit none
    integer,intent(in) :: iunit
    type(gfn0_data),intent(in) :: g0calc
    type(gfn0_results),intent(in) :: res
#ifdef WITH_GFN0
    call gfn0_print_summary(iunit,g0calc,res)
#endif
    return
  end subroutine gfn0_print

!========================================================================================!

  subroutine gfn0_getwbos(g0calc,nat,wbo)
!*************************
!* obtain wbos from gfn0
!*************************
    implicit none
    type(gfn0_data),intent(in) :: g0calc
    integer,intent(in) :: nat
    real(wp),intent(out) :: wbo(nat,nat)
    real(wp),allocatable :: Pa(:,:),Pb(:,:)
    integer :: ndim
    wbo = 0.0_wp
#ifdef WITH_GFN0

    ndim = g0calc%basis%nao
    allocate (Pa(ndim,ndim),Pb(ndim,ndim))
    call density_matrix(ndim,g0calc%wfn%focca,g0calc%wfn%C,Pa)
    call density_matrix(ndim,g0calc%wfn%foccb,g0calc%wfn%C,Pb)
    call get_wbo(nat,g0calc%basis%nao,Pa,Pb, &
    &         g0calc%wfn%S,g0calc%basis%aoat2,wbo)

    deallocate (Pa,Pb)
#endif
  end subroutine gfn0_getwbos

!========================================================================================!

  subroutine gfn0_getdipole(g0calc,mol,dipole)
!*****************************************
!* obtain molecular dipole from gfn0 wfn
!* Note, these come directly from an EEQ!
!*****************************************
    implicit none
    type(gfn0_data),intent(in) :: g0calc
    type(coord),intent(in) :: mol
    real(wp),intent(out) :: dipole(3)
    dipole = 0.0_wp
#ifdef WITH_GFN0
    dipole = matmul(mol%xyz,g0calc%wfn%q)
#endif
  end subroutine gfn0_getdipole

!========================================================================================!
!========================================================================================!
end module gfn0_api

