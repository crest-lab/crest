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
! a small module containing Lennard-Jones potentials
! (mainly for testing other implementations)
!====================================================!

module nonadiabatic_module

  use iso_fortran_env,only:wp => real64

  implicit none

!=========================================================================================!
  !--- private module variables and parameters
  private
  !--- some constants and name mappings
  real(wp),parameter :: bohr = 0.52917726_wp
  real(wp),parameter :: autokcal = 627.509541_wp

  public :: engrad_mean

contains

!========================================================================================!
!> subrotuine engrad_mean
!> given two energies and the respective gradients,
!> the routine returns the arithmetic mean between
!> those two.
!>---------------------------------------------------
  subroutine engrad_mean(nat,e1,e2,grd1,grd2,em,grdm)
    implicit none
    integer,intent(in) :: nat
    real(wp),intent(in) :: e1,e2
    real(wp),intent(in) :: grd1(3,nat)
    real(wp),intent(in) :: grd2(3,nat)
    real(wp),intent(out) :: em
    real(wp),intent(out) :: grdm(3,nat)
    
    em = 0.5_wp * (e1 + e2)

    grdm = 0.5_wp * (grd1 + grd2)

    return
  end subroutine engrad_mean  


end module nonadiabatic_module
