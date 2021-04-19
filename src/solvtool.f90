!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 Sebastian Spicher, Philipp Pracht
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
!===================================================================!
! This file contains routines related to QCG and microsolvation
!===================================================================!
!======================================================!
! main routine 
!======================================================!
subroutine crest_solvtool(env, tim)
  use iso_fortran_env, wp => real64
  use crest_data
  use iomod
  use zdata
  use strucrd
  implicit none

  type(systemdata):: env    ! MAIN STORAGE OS SYSTEM DATA
  type(timer):: tim


  !>-----------------------------------
  call qcg_head()
  call tim%start(2,'QCG') !start a timer
  !>-----------------------------------



  !<----------------------------------
  call tim%stop(2) !stop a timer

  return
end subroutine crest_solvtool
