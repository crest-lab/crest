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
  use iso_fortran_env,only:wp => real64,stdout=>output_unit
  use strucrd
#ifdef WITH_TBLITE
  use mctc_env,only:error_type
  use mctc_io,only:structure_type,new
  use tblite_context_type,only:tblite_ctx => context_type
  use tblite_wavefunction_type,only:wavefunction_type,new_wavefunction
  use tblite_xtb_calculator,only:tblite_calculator => xtb_calculator
  use tblite_xtb_gfn2,only:new_gfn2_calculator
  use tblite_xtb_gfn1,only:new_gfn1_calculator
  use tblite_xtb_ipea1,only:new_ipea1_calculator
  use tblite_xtb_singlepoint,only:xtb_singlepoint
  use tblite_results, only : tblite_resultstype => results_type
#endif
  use wiberg_mayer, only: get_wbo
  implicit none
  private
 
#ifdef WITH_TBLITE
#else
  !> these are placeholders if no tblite is used!
  type :: wavefunction_type  
    integer :: id = 0
  end type wavefunction_type
  type :: tblite_calculator
    integer :: id = 0 
  end type tblite_calculator
  type :: tblite_ctx
    integer :: unit = stdout
  end type tblite_ctx
  type :: tblite_resultstype
    integer :: id = 0
  end type tblite_resultstype
#endif

   !> Type enumerator
   type :: enum_tblite_method
      integer :: unknown   = 0
      integer :: gfn1      = 1
      integer :: gfn2      = 2
      integer :: ipea1     = 3
   end type enum_tblite_method
   type(enum_tblite_method), parameter, public :: xtblvl = enum_tblite_method()

  !> Conversion factor from Kelvin to Hartree
  real(wp),parameter :: ktoau = 3.166808578545117e-06_wp

  integer :: verbosity = 1

  public :: wavefunction_type,tblite_calculator
  public :: tblite_ctx,tblite_resultstype
  public :: tblite_setup,tblite_singlepoint,tblite_addsettings
  public :: tblite_getwbos

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
    type(error_type), allocatable :: error
    
    real(wp) :: etemp_au,energy
    real(wp),allocatable :: grad(:,:)

    !>--- make an mctcmol object from mol
    call tblite_mol2mol(mol,chrg,uhf,mctcmol)

    !>--- select parametrization and set up calculator 
    select case(lvl)
    case( xtblvl%gfn1 )
     call ctx%message("tblite> setting up GFN1-xTB calculation")
     call new_gfn1_calculator(tbcalc, mctcmol)
    case( xtblvl%gfn2 )
     call ctx%message("tblite> setting up GFN2-xTB calculation")
     call new_gfn2_calculator(tbcalc, mctcmol)
    case( xtblvl%ipea1 )
     call ctx%message("tblite> setting up IPEA1-xTB calculation")
     call new_ipea1_calculator(tbcalc, mctcmol) 
    case default 
     call ctx%message("Error: Unknown method in tblite!")
     error stop 
    end select
 
    !>-- setup wavefunction object
    etemp_au = etemp*ktoau
    call new_wavefunction(wfn, mol%nat, tbcalc%bas%nsh,  &
    &              tbcalc%bas%nao, 1, etemp_au) 

    

#else /* WITH_TBLITE */
    write(stdout,*) 'Error: Compiled without tblite support!'
    write(stdout,*) 'Use -DWITH_TBLITE=true in the setup to enable this function'
    error stop 
#endif
  end subroutine tblite_setup

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
    type(error_type), allocatable :: error    
    real(wp) :: sigma(3,3) 

    iostatus = 0 
    energy = 0.0_wp
    gradient(:,:) = 0.0_wp

    !>--- make an mctcmol object from mol
    call tblite_mol2mol(mol,chrg,uhf,mctcmol)

    !>--- call the singlepoint routine
    call xtb_singlepoint(ctx, mctcmol, tbcalc, wfn, accuracy, energy, gradient, & 
    & sigma, verbosity, results=tbres)
 
    if (ctx%failed()) then
      ! Tear down the error stack to send the actual error messages back
      call ctx%message("tblite> Singlepoint calculation failed")
      iostatus = 1
    end if

#else /* WITH_TBLITE */
    iostatus = 0
    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
    write(stdout,*) 'Error: Compiled without tblite support!'
    write(stdout,*) 'Use -DWITH_TBLITE=true in the setup to enable this function'
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
       
    fchrg = real(chrg, wp)  

    !>--- make an mctcmol object from mol
    if(.not.allocated(mol%lat))then
    call new(mctcmol, mol%at, mol%xyz, charge=fchrg, uhf=uhf)
    else
    call new(mctcmol, mol%at, mol%xyz, charge=fchrg, uhf=uhf, lattice=mol%lat)
    endif

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
     tbcalc%save_integrals = (rdwbo .or. saveint)
#endif    
  end subroutine tblite_addsettings
!========================================================================================!
!> tblite_addsettings is used to add other settings from
!> CRESTs calculation object to the tblite_calculator 
  subroutine tblite_getwbos(tbcalc,wfn,tbres,nat,wbo)
    implicit none
    type(tblite_calculator),intent(in) :: tbcalc
    type(wavefunction_type)  :: wfn
    type(tblite_resultstype),intent(in) :: tbres
    integer,intent(in) :: nat
    real(wp),intent(out) :: wbo(nat,nat)

    wbo = 0.0_wp
#ifdef WITH_TBLITE
    call get_wbo(nat, tbcalc%bas%nao, wfn%density, &
    &         tbres%overlap, tbcalc%bas%ao2at, wbo)
#endif    
  end subroutine tblite_getwbos



!========================================================================================!
end module tblite_api
