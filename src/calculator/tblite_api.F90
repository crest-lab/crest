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
  use tblite_context_type,only:context_type
  use tblite_wavefunction_type,only:wavefunction_type,new_wavefunction
  use tblite_xtb_calculator,only:tblite_calculator => xtb_calculator
  use tblite_xtb_gfn2,only:new_gfn2_calculator
  use tblite_xtb_gfn1,only:new_gfn1_calculator
  use tblite_xtb_ipea1,only:new_ipea1_calculator
  use tblite_xtb_singlepoint,only:xtb_singlepoint
#endif
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
#endif

   !> Type enumerator
   type :: enum_tblite_method
      integer :: unknown   = 0
      integer :: gfn1      = 1
      integer :: gfn2      = 2
      integer :: ipea1     = 3
   end type enum_tblite_method
   type(enum_tblite_method), parameter :: xtblvl = enum_tblite_method()

  !> Conversion factor from Kelvin to Hartree
  real(wp),parameter :: ktoau = 3.166808578545117e-06_wp
  !> Calculation accuracy
  real(wp),parameter :: accuracy = 1.0_wp

  integer :: verbosity = 1

  public :: wavefunction_type,tblite_calculator,tblite_setup

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!


  subroutine tblite_setup(mol,chrg,uhf,lvl,wfn,tbcalc)
    implicit none
    type(coord),intent(in)  :: mol
    integer,intent(in)      :: chrg
    integer,intent(in)      :: uhf
    integer,intent(in)      :: lvl
    type(wavefunction_type) :: wfn
    type(tblite_calculator) :: tbcalc
#ifdef WITH_TBLITE
    type(structure_type) :: mctcmol
    type(context_type)   :: ctx
    type(error_type), allocatable :: error
    
    ctx%unit = stdout

    !>--- select parametrization and set up calculator 
    select case(lvl)
    case( xtblvl%gfn1 )
     call ctx%message("tblite> setting up GFN1-xTB calculation")
     call new_gfn2_calculator(tbcalc, mctcmol)
    case( xtblvl%gfn2 )
     call ctx%message("tblite> setting up GFN2-xTB calculation")
     call new_gfn1_calculator(tbcalc, mctcmol)
    case( xtblvl%ipea1 )
     call ctx%message("tblite> setting up IPEA1-xTB calculation")
     call new_ipea1_calculator(tbcalc, mctcmol) 
    case default 
     call ctx%message("Error: Unknown method in tblite!")
     error stop 
    end select
  

#else /* WITH_TBLITE */
    write(stdout,*) 'Error: Compiled without tblite support!'
    write(stdout,*) 'Use -DWITH_TBLITE=true in the setup to enable this function'
    error stop 
#endif
  end subroutine tblite_setup



!========================================================================================!
#ifdef WITH_TBLITE
  subroutine tblite_mol2mol(mol,mctcmol)
    implicit none
    type(coord) :: mol
    type(structure_type) :: mctcmol

  end subroutine tblite_mol2mol
#endif
!========================================================================================!
end module tblite_api
