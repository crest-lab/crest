!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2026 Philipp Pracht
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

module molecule_parameters
  use iso_fortran_env,only:wp => real64,stdout=>output_unit
  use iso_c_binding
  implicit none

  public :: wp,stdout  !> RE-EXPORTS
!&<
!>--- some constants and name mappings
  real(wp),parameter,public :: bohr     = 0.52917726_wp
  real(wp),parameter,public :: aatoau   = 1.0_wp/bohr
  real(wp),parameter,public :: autoaa   = bohr
  real(wp),parameter,public :: autokcal = 627.509541_wp
  real(wp),parameter,public :: autoeV   = 27.211324570273_wp

!>--- global extxyz output unit preference (mutable, set at runtime via TOML)
  character(len=32),public :: extxyz_units_global = 'hartree'

!>-- filetypes as integers
  type ,private:: enum_coordtype
    integer :: unknown    = 0
    integer :: turbomole  = 1
    integer :: xyz        = 2
    integer :: extxyz     = 22
    integer :: sdf        = 3
    integer :: sdfV2000   = 31
    integer :: sdfV3000   = 32
    integer :: PDB        = 4
  end type enum_coordtype
  type(enum_coordtype), parameter,public :: coordtype = enum_coordtype()

  !> Element symbols
  character(len=2),parameter :: PSE(118) = [ &
   & 'H ',                                                                                'He', &
   & 'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne', &
   & 'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar', &
   & 'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
   & 'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe', &
   & 'Cs','Ba','La',                                                                            &
   &                'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',      &
   &                'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn', &
   & 'Fr','Ra','Ac',                                                                            &
   &                'Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',      &
   &                'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og' ]
!&>

  public :: ncore

! ══════════════════════════════════════════════════════════════════════════════
contains  !> MODULE PROCEDURES START HERE
! ══════════════════════════════════════════════════════════════════════════════

  pure elemental integer function ncore(at)
    integer,intent(in) :: at
    if (at .le. 2) then
      ncore = 0
    elseif (at .le. 10) then
      ncore = 2
    elseif (at .le. 18) then
      ncore = 10
    elseif (at .le. 29) then   !zn
      ncore = 18
    elseif (at .le. 36) then
      ncore = 28
    elseif (at .le. 47) then
      ncore = 36
    elseif (at .le. 54) then
      ncore = 46
    elseif (at .le. 71) then
      ncore = 54
    elseif (at .le. 79) then
      ncore = 68
    elseif (at .le. 86) then
      ncore = 78
    elseif (at .le. 103) then !> Rn core
      ncore = 86
    elseif (at .le. 118) then !> Og core
      ncore = 102
    end if
  end function ncore

end module molecule_parameters
