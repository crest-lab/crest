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
  use tblite_xtb,xtb_calculator => xtb_calculator
  use tblite_results,only:tblite_resultstype => results_type
  use tblite_wavefunction_mulliken,only:get_molecular_dipole_moment
  use tblite_ceh_singlepoint,only:ceh_guess
  use tblite_ceh_ceh,only:new_ceh_calculator
#endif
  use wiberg_mayer
  implicit none
  private

#ifndef WITH_TBLITE
  !> these are placeholders if no tblite is used!
  type :: wavefunction_type
    integer :: id = 0
  end type wavefunction_type
  type :: xtb_calculator
    integer :: id = 0
  end type xtb_calculator
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
    integer  :: lvl = 0
    real(wp) :: accuracy = 1.0_wp
    type(wavefunction_type)     :: wfn
    type(xtb_calculator)        :: calc
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
    !> the guesses can be used for charges, but NOT for e+grd!
    integer :: eeq = 4
    integer :: ceh = 5
  end type enum_tblite_method
  type(enum_tblite_method),parameter,public :: xtblvl = enum_tblite_method()

  !> Conversion factor from Kelvin to Hartree
  real(wp),parameter :: ktoau = 3.166808578545117e-06_wp

  integer :: verbosity = 0
  !> IMPORTANT: tblite is not entirely thread-safe
  !> if verbosity is >0. We'll have to turn it off.
  !> At least for statically compiled binaries

  public :: wavefunction_type,xtb_calculator
  public :: tblite_ctx,tblite_resultstype
  public :: tblite_setup,tblite_singlepoint,tblite_addsettings
  public :: tblite_getwbos
  public :: tblite_add_solv
  public :: tblite_getcharges
  public :: tblite_getdipole

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine tblite_setup(mol,chrg,uhf,lvl,etemp,tblite)
!*****************************************************************
!* subroutine tblite_setup initializes the tblite object which is
!* passed between the CREST calculators and this module
!*****************************************************************
    implicit none
    type(coord),intent(in)  :: mol
    integer,intent(in)      :: chrg
    integer,intent(in)      :: uhf
    type(tblite_data),intent(inout) :: tblite
    integer,intent(in)      :: lvl
    real(wp),intent(in)     :: etemp
#ifdef WITH_TBLITE
    type(structure_type) :: mctcmol
    type(error_type),allocatable :: error

    real(wp) :: etemp_au,energy
    real(wp),allocatable :: grad(:,:)
    logical :: pr

    pr = (tblite%ctx%verbosity > 0)

!>--- make an mctcmol object from mol
    call tblite_mol2mol(mol,chrg,uhf,mctcmol)

!>--- select parametrization and set up calculator
    tblite%lvl = lvl
    select case (tblite%lvl)
    case (xtblvl%gfn1)
      if (pr) call tblite%ctx%message("tblite> setting up GFN1-xTB calculation")
      call new_gfn1_calculator(tblite%calc,mctcmol)
    case (xtblvl%gfn2)
      if (pr) call tblite%ctx%message("tblite> setting up GFN2-xTB calculation")
      call new_gfn2_calculator(tblite%calc,mctcmol)
    case (xtblvl%ipea1)
      if (pr) call tblite%ctx%message("tblite> setting up IPEA1-xTB calculation")
      call new_ipea1_calculator(tblite%calc,mctcmol)
    case (xtblvl%ceh)
      if (pr) call tblite%ctx%message("tblite> setting up CEH calculation")
      call new_ceh_calculator(tblite%calc,mctcmol)
    case (xtblvl%eeq)
      if (pr) call tblite%ctx%message("tblite> setting up D4 EEQ charges calculation")
      call new_ceh_calculator(tblite%calc,mctcmol) !> doesn't matter but needs initialization
    case default
      call tblite%ctx%message("Error: Unknown method in tblite!")
      error stop
    end select

!>-- setup wavefunction object
    etemp_au = etemp*ktoau
    call new_wavefunction(tblite%wfn,mol%nat,tblite%calc%bas%nsh,  &
    &              tblite%calc%bas%nao,1,etemp_au)

#else /* WITH_TBLITE */
    write (stdout,*) 'Error: Compiled without tblite support!'
    write (stdout,*) 'Use -DWITH_TBLITE=true in the setup to enable this function'
    error stop
#endif
  end subroutine tblite_setup

!========================================================================================!

  subroutine tblite_add_solv(mol,chrg,uhf,tblite,smodel,solvent)
!***********************************************************************
!* This subroutine sets up the implicit solvation container for tblite
!***********************************************************************
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
    type(tblite_data),intent(inout) :: tblite
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

    if (.not.allocated(smodel).or..not.allocated(solvent)) then
      return
    end if
    pr = (tblite%ctx%verbosity > 0)

!>--- some tblite calculators have nothing to do with implicit solvation
    if (tblite%lvl > 3) then
      if (pr) call tblite%ctx%message("tblite> skipping implicit solvation setup for this potential")
      return
    end if

!>--- make an mctcmol object from mol
    call tblite_mol2mol(mol,chrg,uhf,mctcmol)

    if (pr) call tblite%ctx%message("tblite> setting up tblite implicit solvation")
!>--- generat solvation parametrization
    if (solvent == 'h2o') then !> special case: tblite doesn't know 'h2o', only 'water' ...
      solvdum = 'water'
    else
      solvdum = solvent
    end if
    solv_data = get_solvent_data(solvdum)
    if (solv_data%eps <= 0.0_wp) then
      if (pr) call tblite%ctx%message("tblite> Unknown solvent!")
      return
    end if
    allocate (solv_inp)
    select case (trim(smodel))
    case ('gbsa')
      if (pr) call tblite%ctx%message("tblite> using GBSA/"//solvdum)
      allocate (solv_inp%alpb)
      solv_inp%alpb = alpb_input(solv_data%eps,alpb=.false.)
    case ('cpcm')
      if (pr) call tblite%ctx%message("tblite> using CPCM/"//solvdum)
      allocate (solv_inp%cpcm)
      solv_inp%cpcm = cpcm_input(solv_data%eps)
    case ('alpb')
      if (pr) call tblite%ctx%message("tblite> using ALPB/"//solvdum)
      allocate (solv_inp%alpb)
      solv_inp%alpb = alpb_input(solv_data%eps,alpb=.true.)
    case default
      if (pr) call tblite%ctx%message("tblite> Unknown tblite implicit solvation model!")
      return
    end select
    str = 'tblite> WARNING: implicit solvation energies are not entirely '// &
    &'consistent with the xtb implementation.'
    if (pr) call tblite%ctx%message(str)

!>--- add to calculator
    call new_solvation(solv,mctcmol,solv_inp,error)
    if (allocated(error)) then
      if (pr) call tblite%ctx%message("tblite> failed to set up tblite implicit solvation!")
      return
    end if
    call move_alloc(solv,cont)
    call tblite%calc%push_back(cont)

    deallocate (solv_inp)

#else /* WITH_TBLITE */
    write (stdout,*) 'Error: Compiled without tblite support!'
    write (stdout,*) 'Use -DWITH_TBLITE=true in the setup to enable this function'
    error stop
#endif
  end subroutine tblite_add_solv

!========================================================================================!

  subroutine tblite_singlepoint(mol,chrg,uhf,tblite,energy,gradient,iostatus)
!**************************************************
!* The actual calculator call.
!* The tblite object must be set up at this point
!**************************************************
    implicit none
    type(coord),intent(in)   :: mol
    integer,intent(in)       :: chrg
    integer,intent(in)       :: uhf
    type(tblite_data),intent(inout) :: tblite
    real(wp),intent(out)     :: energy
    real(wp),intent(out)     :: gradient(3,mol%nat)
    integer,intent(out)      :: iostatus
#ifdef WITH_TBLITE
    type(structure_type) :: mctcmol
    type(error_type),allocatable :: error
    real(wp) :: sigma(3,3)
    logical :: pr

    iostatus = 0
    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
    pr = (tblite%ctx%verbosity > 0)

!>--- make an mctcmol object from mol
    call tblite_mol2mol(mol,chrg,uhf,mctcmol)

!>--- call the singlepoint routine
    select case (tblite%lvl)
    case default
      call xtb_singlepoint(tblite%ctx,mctcmol,tblite%calc,tblite%wfn,tblite%accuracy, &
     &                    energy,gradient, &
     &                    sigma,verbosity,results=tblite%res)
    case (xtblvl%ceh)
      call ceh_guess(tblite%ctx,tblite%calc,mctcmol,error,tblite%wfn, &
      &              tblite%accuracy,verbosity)
    case(xtblvl%eeq)
      call eeq_guess(mctcmol, tblite%calc, tblite%wfn)
    end select

    if (tblite%ctx%failed()) then
      !> Tear down the error stack to send the actual error messages back
      if (pr) call tblite%ctx%message("tblite> Singlepoint calculation failed")
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
!*************************************************************************
!* tblite uses its own molecule type thats different from our coord type
!* This routine does the minimal conversion
!*************************************************************************
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

  subroutine tblite_addsettings(tblite,maxscc,rdwbo,saveint,accuracy)
!**********************************************************
!* tblite_addsettings is used to add other settings from
!* CRESTs calculation object to the xtb_calculator
!**********************************************************
    implicit none
    type(tblite_data),intent(inout) :: tblite
    integer,intent(in) :: maxscc
    logical,intent(in) :: rdwbo
    logical,intent(in) :: saveint
    real(wp),intent(in) :: accuracy
#ifdef WITH_TBLITE
    tblite%calc%max_iter = maxscc
    tblite%calc%save_integrals = (rdwbo.or.saveint)
    tblite%accuracy = accuracy
#endif
  end subroutine tblite_addsettings

!========================================================================================!

  subroutine tblite_getwbos(tblite,nat,wbo)
!***************************
!* obtain wbos from tblite
!***************************
    implicit none
    type(tblite_data),intent(inout) :: tblite
    integer,intent(in) :: nat
    real(wp),intent(out) :: wbo(nat,nat)
    real(wp),allocatable :: S(:,:)
    integer :: nao,i
    real(wp),allocatable :: focca(:),foccb(:)
    real(wp),allocatable :: Pa(:,:),Pb(:,:)
    wbo = 0.0_wp
#ifdef WITH_TBLITE
    select case (tblite%lvl)
    case default

      nao = tblite%calc%bas%nao
      allocate(Pa(nao,nao),Pb(nao,nao))
      call split_foccab(nao,tblite%wfn%focc, tblite%wfn%nel(1), tblite%wfn%nel(2), &
      & focca, foccb)
      call density_matrix(nao,focca,tblite%wfn%coeff(:,:,1),Pa)
      call density_matrix(nao,foccb,tblite%wfn%coeff(:,:,1),Pb)
      call get_wbo(nat, nao, Pa,Pb, tblite%res%overlap, tblite%calc%bas%ao2at, wbo)

    case (xtblvl%ceh)
    !> no external access to the overlap in CEH, hence use the Wiberg BO with S=I
      nao = tblite%calc%bas%nao
      allocate(S(nao,nao), source=0.0_wp) 
      do i=1,nao
        S(i,i) = 1.0_wp
      enddo
      call get_wbo_rhf(nat,tblite%calc%bas%nao,tblite%wfn%density, &
      &                S,tblite%calc%bas%ao2at,wbo)
      wbo = wbo*2.0_wp !> somehow this is much better

    case( xtblvl%eeq )
      wbo = 0.0_wp
    end select
#endif
  end subroutine tblite_getwbos

!========================================================================================!

  subroutine tblite_getcharges(mol,tblite,q)
!**************************************
!* obtain molecular dipole from tblite
!**************************************
    implicit none
    type(coord) :: mol
    type(tblite_data),intent(inout) :: tblite
    real(wp),intent(out) :: q(mol%nat)
#ifdef WITH_TBLITE
    q = 0.0_wp
    q(:) = tblite%wfn%qat(:,1)
#else
    q = 0.0_wp
#endif
  end subroutine tblite_getcharges

!========================================================================================!

  subroutine tblite_getdipole(mol,chrg,uhf,tblite,dipole)
!**************************************
!* obtain molecular dipole from tblite
!**************************************
    implicit none
    type(coord) :: mol
    integer :: chrg,uhf
    type(tblite_data),intent(inout) :: tblite
    real(wp),intent(out) :: dipole(3)
#ifdef WITH_TBLITE
    type(structure_type) :: mctcmol
    dipole = 0.0_wp
    call tblite_mol2mol(mol,chrg,uhf,mctcmol)
    !> dipole moment is obtained from molecular charges and atomic dipole moments
    call get_molecular_dipole_moment(mctcmol,tblite%wfn%qat(:,1), &
    &    tblite%wfn%dpat(:,:,1),dipole)
#else
    dipole = 0.0_wp
#endif
  end subroutine tblite_getdipole

!========================================================================================!
!========================================================================================!
end module tblite_api
