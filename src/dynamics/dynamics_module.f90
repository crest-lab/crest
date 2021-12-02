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

module dynamics_module

   use iso_fortran_env, only: wp => real64
   use strucrd

   implicit none

   !=========================================================================================! 
   !--- private module variables and parameters
   private
     integer :: i,j,k,l,ich,och,io
     logical :: ex

   !--- some constants and name mappings
     real(wp),parameter :: bohr  = 0.52917726_wp
     real(wp),parameter :: autokcal = 627.509541_wp
     !-- filetypes as integers
     integer,parameter :: type_md     = 1
     integer,parameter :: type_mtd    = 2
  

     public :: mddata
   !=========================================================================================!
   !data object that contains settings for a molecular dynamics simulation.
   type :: mddata
       
       integer  :: md_index = 0         ! some index for parallelization
       integer  :: simtype = type_md    ! type of the molecular dynamics simulation
    !--- data
       real(wp) :: length_ps = 0.0_wp   ! total simulation length in ps
       integer  :: length_steps = 0     ! total simulation length in steps
       real(wp) :: tstep = 0.0_wp       ! timestep in fs
       integer  :: sdump = 0            ! snapshot dump to trajectory every x steps

       real(wp) :: tsoll = 298.15_wp    ! wanted temperature
       integer  :: shake = 0            ! SHAKE algorithm selection
       logical  :: thermostat = .true.  ! apply thermostat?

   end type mddata


contains


end module dynamics_module
