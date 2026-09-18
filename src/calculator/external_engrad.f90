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

module crest_external_engrad
!****************************************************************************
!* Interface for host-supplied (externally defined) energy+gradient
!* routines that can be injected into the CREST calculator via the
!* jobtype%external mechanism.
!*
!* This module bundles the contract that a host program (which links CREST
!* as a library) must satisfy, kept separate from the calculator data
!* structures so it can be used standalone:
!*   - engrad_interface : native Fortran callback
!*
!* A host written in C/C++/Python can still be plugged in by wrapping its
!* routine in a thin Fortran shim matching engrad_interface.
!****************************************************************************
  use iso_fortran_env,only:wp => real64
  implicit none
  public

!=========================================================================================!

  abstract interface
    subroutine engrad_interface(nat,at,xyz,chrg,uhf,energy,gradient,iostatus,userdata)
    !*********************************************************************
    !* Interface for an externally supplied energy+gradient routine.
    !* A host program that links CREST as a library can implement a
    !* routine matching this signature and register it on a
    !* calculation_settings object (see %set_external). CREST will then
    !* call it through the regular engrad/potential_core dispatch without
    !* any compile-time knowledge of the implementation.
    !*
    !* In/output:
    !*  nat      : number of atoms                                 (in)
    !*  at(nat)  : atomic numbers                                  (in)
    !*  xyz(3,*) : Cartesian coordinates in Bohr (CREST internal)  (in)
    !*  chrg     : molecular charge                                (in)
    !*  uhf      : Nalpha-Nbeta = 2S                               (in)
    !*  energy   : total energy in Hartree                         (out)
    !*  gradient : gradient dE/dxyz in Hartree/Bohr                (out)
    !*  iostatus : 0 on success, /=0 signals a failed calculation  (out)
    !*  userdata : optional opaque host context, passed verbatim   (inout)
    !*********************************************************************
      import :: wp
      integer,intent(in)              :: nat
      integer,intent(in)              :: at(nat)
      real(wp),intent(in)             :: xyz(3,nat)
      integer,intent(in)              :: chrg
      integer,intent(in)              :: uhf
      real(wp),intent(out)            :: energy
      real(wp),intent(out)            :: gradient(3,nat)
      integer,intent(out)             :: iostatus
      class(*),intent(inout),optional :: userdata
    end subroutine engrad_interface
  end interface

!=========================================================================================!
end module crest_external_engrad
