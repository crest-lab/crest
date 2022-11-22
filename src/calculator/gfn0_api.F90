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
  use wiberg_mayer,only:get_wbo
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
  public :: gfn0_gen_occ

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

#else /* WITH_TBLITE */
    write (stdout,*) 'Error: Compiled without GFN0-xTB support!'
    write (stdout,*) 'Use -DWITH_GFN0=true in the setup to enable this function'
    error stop
#endif
  end subroutine gfn0_setup

!========================================================================================!
!> gfn0_addsettings is used to add other settings from
!> CRESTs calculation object to the gfn0_data
  subroutine gfn0_addsettings(mol,g0calc,solv,model)
    implicit none
    type(coord),intent(in)  :: mol
    type(gfn0_data),intent(inout) :: g0calc
    character(len=*),intent(in),optional :: solv
    character(len=*),intent(in),optional :: model
#ifdef WITH_GFN0
    !> add solvation?
    if (present(solv)) then
      if (allocated(g0calc%gbsa)) deallocate (g0calc%gbsa)
      allocate (g0calc%gbsa)
      if (present(model)) then
        select case(model)
        case('alpb')
        call gfn0_gbsa_init(mol%nat,mol%at,.true.,solv,g0calc%gbsa)
        case default !> default GBSA
        call gfn0_gbsa_init(mol%nat,mol%at,.false.,solv,g0calc%gbsa)
        end select
      else
        call gfn0_gbsa_init(mol%nat,mol%at,.false.,solv,g0calc%gbsa)
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
    if (present(res)) then
      call gfn0_singlepoint(mol%nat,mol%at,mol%xyz,chrg,uhf,g0calc, &
      &          energy,gradient,fail,res)
    else
      call gfn0_singlepoint(mol%nat,mol%at,mol%xyz,chrg,uhf,g0calc, &
      &          energy,gradient,fail)
    end if
    if(fail)then
      iostatus = -1
    endif
#else
    write (stdout,*) 'Error: Compiled without GFN0-xTB support!'
    write (stdout,*) 'Use -DWITH_GFN0=true in the setup to enable this function'
    error stop
#endif
  end subroutine gfn0_sp

!========================================================================================!

  subroutine gfn0_sp_occ(mol,chrg,uhf,nlev,occ,g0calc,energies,gradients,iostatus,res)
    implicit none
    !> INPUT
    type(coord),intent(in)  :: mol
    integer,intent(in) :: chrg
    integer,intent(in) :: uhf
    integer,intent(in) :: nlev
    type(gfn0_data),intent(inout) :: g0calc
!    real(wp),intent(in) :: occ(g0calc%basis%nao, nlev)
    real(wp),intent(in) :: occ(:,:)
    !> OUTPUT
    real(wp),intent(out) :: energies(nlev)
    real(wp),intent(out) :: gradients(3,mol%nat,nlev)
    integer,intent(out) :: iostatus
    type(gfn0_results),intent(inout),optional :: res
    !> LOCAL
    logical :: fail
    energies = 0.0_wp
    gradients = 0.0_wp
    iostatus = 0 
    fail = .false.
#ifdef WITH_GFN0
    if (present(res)) then
      call gfn0_occ_singlepoint(mol%nat,mol%at,mol%xyz,chrg,uhf,nlev,occ,g0calc, &
      &          energies,gradients,fail,res)
    else
      call gfn0_occ_singlepoint(mol%nat,mol%at,mol%xyz,chrg,uhf,nlev,occ,g0calc, &
      &          energies,gradients,fail)
    end if
    if(fail)then
      iostatus = -1
    endif
#else
    write (stdout,*) 'Error: Compiled without GFN0-xTB support!'
    write (stdout,*) 'Use -DWITH_GFN0=true in the setup to enable this function'
    error stop
#endif
  end subroutine gfn0_sp_occ

!========================================================================================!

  subroutine gfn0_gen_occ(nlev,nel,nao,active,occ)
    implicit none
    integer,intent(in) :: nlev
    integer,intent(in) :: nel
    integer,intent(in) :: nao
    integer,intent(in) :: active(:,:)
    real(wp),intent(out) :: occ(nao,nlev)
    integer :: i
    occ = 0.0_wp
    do i=1,nlev
     call generate_config(nel,nao,occ(:,i),active(:,i))
    enddo
  end subroutine gfn0_gen_occ

!========================================================================================!
!========================================================================================!
end module gfn0_api

