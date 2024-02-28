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
! module tblite_api
! An interface of CREST to tblite
!====================================================!

module tblite_api
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use strucrd
#ifdef WITH_TBLITE
  use mctc_env,only:error_type
  use mctc_io,only:structure_type,new
  use tblite_context_type,only:tblite_ctx => context_type
  use tblite_wavefunction_type,only:wavefunction_type,new_wavefunction
  use tblite_wavefunction,only:sad_guess,eeq_guess
  use tblite_xtb, tblite_calculator => xtb_calculator 
  use tblite_results,only:tblite_resultstype => results_type
#endif
  use wiberg_mayer,only:get_wbo_rhf
  implicit none
  private

#ifndef WITH_TBLITE
  !> these are placeholders if no tblite is used!
  type :: wavefunction_type
    integer :: id = 0
  end type wavefunction_type
  type :: tblite_calculator
    integer :: id = 0
  end type tblite_calculator
  type :: tblite_ctx
    integer :: unit = stdout
    integer :: verbosity = 0
  end type tblite_ctx
  type :: tblite_resultstype
    integer :: id = 0
  end type tblite_resultstype
  type :: tblite_solvation_type
    integer :: id = 0
  end type tblite_solvation_type
#endif

!>--- tblite calculator bundle
  type :: tblite_data
    type(wavefunction_type)     :: wfn
    type(tblite_calculator)     :: calc
    type(tblite_ctx)            :: ctx
    type(tblite_resultstype)    :: res
  end type tblite_data
  public :: tblite_data

  !> Type enumerator
  type :: enum_tblite_method
    integer :: unknown = 0
    integer :: gfn1 = 1
    integer :: gfn2 = 2
    integer :: ipea1 = 3
  end type enum_tblite_method
  type(enum_tblite_method),parameter,public :: xtblvl = enum_tblite_method()

  !> Conversion factor from Kelvin to Hartree
  real(wp),parameter :: ktoau = 3.166808578545117e-06_wp

  integer :: verbosity = 0 
  !> IMPORTANT: tblite is not entirely thread-safe
  !> if verbosity is >0. We'll have to turn it off.
  !> At least for statically compiled binaries

  public :: wavefunction_type,tblite_calculator
  public :: tblite_ctx,tblite_resultstype
  public :: tblite_setup,tblite_singlepoint,tblite_addsettings
  public :: tblite_getwbos
  public :: tblite_add_solv

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine tblite_setup(mol,chrg,uhf,lvl,etemp,ctx,wfn,tbcalc)
    implicit none
    type(coord),intent(in)  :: mol
    integer,intent(in)      :: chrg
    integer,intent(in)      :: uhf
    integer,intent(in)      :: lvl
    real(wp),intent(in)     :: etemp
    type(tblite_ctx)        :: ctx
    type(wavefunction_type) :: wfn
    type(tblite_calculator) :: tbcalc
#ifdef WITH_TBLITE
    type(structure_type) :: mctcmol
    type(error_type),allocatable :: error

    real(wp) :: etemp_au,energy
    real(wp),allocatable :: grad(:,:)
    logical :: pr

    pr = (ctx%verbosity > 0)

    !>--- make an mctcmol object from mol
    call tblite_mol2mol(mol,chrg,uhf,mctcmol)

    !>--- select parametrization and set up calculator
    select case (lvl)
    case (xtblvl%gfn1)
      if(pr) call ctx%message("tblite> setting up GFN1-xTB calculation")
      call new_gfn1_calculator(tbcalc,mctcmol)
    case (xtblvl%gfn2)
      if(pr) call ctx%message("tblite> setting up GFN2-xTB calculation")
      call new_gfn2_calculator(tbcalc,mctcmol)
    case (xtblvl%ipea1)
      if(pr) call ctx%message("tblite> setting up IPEA1-xTB calculation")
      call new_ipea1_calculator(tbcalc,mctcmol)
    case default
      call ctx%message("Error: Unknown method in tblite!")
      error stop
    end select

    !>-- setup wavefunction object
    etemp_au = etemp*ktoau
    call new_wavefunction(wfn,mol%nat,tbcalc%bas%nsh,  &
    &              tbcalc%bas%nao,1,etemp_au)

#else /* WITH_TBLITE */
    write (stdout,*) 'Error: Compiled without tblite support!'
    write (stdout,*) 'Use -DWITH_TBLITE=true in the setup to enable this function'
    error stop
#endif
  end subroutine tblite_setup

!========================================================================================!

  subroutine tblite_add_solv(mol,chrg,uhf,ctx,wfn,tbcalc,smodel,solvent)
#ifdef WITH_TBLITE
    use tblite_container,only:container_type
    use tblite_solvation,only:new_solvation,tblite_solvation_type => solvation_type, &
    &                            solvent_data,get_solvent_data,solvation_input,  &
    &                            cpcm_input,alpb_input,alpb_solvation
#endif
    implicit none
    type(coord),intent(in)  :: mol
    integer,intent(in)      :: chrg
    integer,intent(in)      :: uhf
    type(tblite_ctx)        :: ctx
    type(wavefunction_type) :: wfn
    type(tblite_calculator) :: tbcalc
    character(len=:),allocatable,intent(in) :: smodel
    character(len=:),allocatable,intent(in) :: solvent
#ifdef WITH_TBLITE
    type(structure_type) :: mctcmol
    type(error_type),allocatable :: error

    class(container_type),allocatable :: cont
    class(tblite_solvation_type),allocatable :: solv
    type(solvation_input),allocatable :: solv_inp
    type(solvent_data) :: solv_data
    character(len=:),allocatable :: str,solvdum
    logical :: pr

    real(wp) :: etemp_au,energy

    if (.not.allocated(smodel).or..not.allocated(solvent)) then
      return
    end if

    !> special case: tblite dosen't know 'h2o', only 'water'

    pr = (ctx%verbosity > 0)

    !>--- make an mctcmol object from mol
    call tblite_mol2mol(mol,chrg,uhf,mctcmol)

    if(pr) call ctx%message("tblite> setting up tblite implicit solvation")
    !>--- generat solvation parametrization
    if(solvent == 'h2o')then !> special case: tblite doesn't know 'h2o', only 'water' ...
      solvdum = 'water'
    else
      solvdum = solvent
    endif
    solv_data = get_solvent_data(solvdum)
    if (solv_data%eps <= 0.0_wp) then
      if(pr) call ctx%message("tblite> Unknown solvent!")
      return
    end if
    allocate (solv_inp)
    select case (trim(smodel))
    case ('gbsa')
      if(pr) call ctx%message("tblite> using GBSA/"//solvdum)
      allocate (solv_inp%alpb)
      solv_inp%alpb = alpb_input(solv_data%eps,alpb=.false.)
    case ('cpcm')
      if(pr) call ctx%message("tblite> using CPCM/"//solvdum)
      allocate (solv_inp%cpcm)
      solv_inp%cpcm = cpcm_input(solv_data%eps)
    case ('alpb')
      if(pr) call ctx%message("tblite> using ALPB/"//solvdum)
      allocate (solv_inp%alpb)
      solv_inp%alpb = alpb_input(solv_data%eps,alpb=.true.)
    case default
      if(pr) call ctx%message("tblite> Unknown tblite implicit solvation model!")
      return
    end select
    str = 'tblite> WARNING: implicit solvation energies are not entirely '// &
    &'consistent with the xtb implementation.'
    if(pr) call ctx%message(str)

!>--- add to calculator
    call new_solvation(solv,mctcmol,solv_inp,error)
    if (allocated(error)) then
      if(pr) call ctx%message("tblite> failed to set up tblite implicit solvation!")
      return
    end if
    call move_alloc(solv,cont)
    call tbcalc%push_back(cont)

    deallocate (solv_inp)

#else /* WITH_TBLITE */
    write (stdout,*) 'Error: Compiled without tblite support!'
    write (stdout,*) 'Use -DWITH_TBLITE=true in the setup to enable this function'
    error stop
#endif
  end subroutine tblite_add_solv

!========================================================================================!

  subroutine tblite_singlepoint(mol,chrg,uhf,accuracy,ctx,wfn,tbcalc,energy,gradient,tbres,iostatus)
    implicit none
    type(coord),intent(in)   :: mol
    integer,intent(in)       :: chrg
    integer,intent(in)       :: uhf
    real(wp),intent(in)      :: accuracy
    type(tblite_ctx)         :: ctx
    type(wavefunction_type)  :: wfn
    type(tblite_calculator)  :: tbcalc
    real(wp),intent(out)     :: energy
    real(wp),intent(out)     :: gradient(3,mol%nat)
    type(tblite_resultstype) :: tbres
    integer,intent(out)      :: iostatus
#ifdef WITH_TBLITE
    type(structure_type) :: mctcmol
    type(error_type),allocatable :: error
    real(wp) :: sigma(3,3)
    logical :: pr

    iostatus = 0
    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
    pr = (ctx%verbosity > 0)

    !>--- make an mctcmol object from mol
    call tblite_mol2mol(mol,chrg,uhf,mctcmol)

    !>--- call the singlepoint routine
    call xtb_singlepoint(ctx,mctcmol,tbcalc,wfn,accuracy,energy,gradient, &
    & sigma,verbosity,results=tbres)

    if (ctx%failed()) then
      ! Tear down the error stack to send the actual error messages back
      if(pr) call ctx%message("tblite> Singlepoint calculation failed")
      iostatus = 1
    end if

#else /* WITH_TBLITE */
    iostatus = 0
    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
    write (stdout,*) 'Error: Compiled without tblite support!'
    write (stdout,*) 'Use -DWITH_TBLITE=true in the setup to enable this function'
    error stop
#endif
  end subroutine tblite_singlepoint

!========================================================================================!
#ifdef WITH_TBLITE
  subroutine tblite_mol2mol(mol,chrg,uhf,mctcmol)
    implicit none
    !> input & output
    type(coord) :: mol
    integer,intent(in) :: chrg
    integer,intent(in) :: uhf
    type(structure_type),intent(out) :: mctcmol
    !> locals
    real(wp) :: fchrg

    fchrg = real(chrg,wp)

    !>--- make an mctcmol object from mol
    if (.not.allocated(mol%lat)) then
      call new(mctcmol,mol%at,mol%xyz,charge=fchrg,uhf=uhf)
    else
      call new(mctcmol,mol%at,mol%xyz,charge=fchrg,uhf=uhf,lattice=mol%lat)
    end if

  end subroutine tblite_mol2mol
#endif
!========================================================================================!
!> tblite_addsettings is used to add other settings from
!> CRESTs calculation object to the tblite_calculator
  subroutine tblite_addsettings(tbcalc,maxscc,rdwbo,saveint)
    implicit none
    type(tblite_calculator),intent(inout) :: tbcalc
    integer,intent(in) :: maxscc
    logical,intent(in) :: rdwbo
    logical,intent(in) :: saveint
#ifdef WITH_TBLITE
    tbcalc%max_iter = maxscc
    tbcalc%save_integrals = (rdwbo.or.saveint)
#endif
  end subroutine tblite_addsettings
!========================================================================================!
!> obtain wbos from tblite
  subroutine tblite_getwbos(tbcalc,wfn,tbres,nat,wbo)
    implicit none
    type(tblite_calculator),intent(in) :: tbcalc
    type(wavefunction_type)  :: wfn
    type(tblite_resultstype),intent(in) :: tbres
    integer,intent(in) :: nat
    real(wp),intent(out) :: wbo(nat,nat)

    wbo = 0.0_wp
#ifdef WITH_TBLITE
    call get_wbo_rhf(nat,tbcalc%bas%nao,wfn%density, &
    &         tbres%overlap,tbcalc%bas%ao2at,wbo)
#endif
  end subroutine tblite_getwbos

!========================================================================================!
end module tblite_api
