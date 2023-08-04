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
!
! Routines were adapted from the xtb code (github.com/grimme-lab/xtb)
! under the Open-source software LGPL-3.0 Licencse.
!================================================================================!

module optimize_module
  use iso_fortran_env,only:wp => real64
  use crest_calculator
  use strucrd
  use ancopt_module,only:ancopt
  implicit none

  real(wp),private,parameter :: autoaa = 0.52917726_wp
  real(wp),private,parameter :: aatoau = 1.0_wp/autoaa

  public :: optimize_geometry

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine optimize_geometry(mol,molnew,calc,etot,grd,pr,wr,iostatus)
    implicit none
    !> Input
    type(coord)    :: mol
    type(calcdata) :: calc
    logical,intent(in)        :: pr
    logical,intent(in)        :: wr
    !> Output
    type(coord)   :: molnew
    integer,intent(out)       :: iostatus
    real(wp),intent(inout)    :: etot
    real(wp),intent(inout)    :: grd(3,mol%nat)

    iostatus = -1
    !$omp critical
    molnew%at = mol%at  !> do not overwrite original geometry
    molnew%xyz = mol%xyz
    molnew%nat = mol%nat
    !$omp end critical

    !> initial singlepoint
    call engrad(molnew,calc,etot,grd,iostatus)

    !> optimization
    call ancopt(molnew,calc,etot,grd,pr,wr,iostatus)

    return
  end subroutine optimize_geometry

!========================================================================================!
end module optimize_module
