!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020 Philipp Pracht
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
!======================================================!
! metadynamics based sampling with logfermi potential
!======================================================!
subroutine compress(env, tim)
  use iso_fortran_env, wp => real64
  use crest_data
  use iomod
  use zdata
  use strucrd, only: rdensembleparam,rdensemble
  implicit none

  !type(options) :: opt
  type(systemdata) :: env
  type(timer) :: tim


  !--- start timer
  call tim%start(1, 'Setup')

  !--- 1. set up reactor metadynamics
  !    wall-potential, mass density, simulation length
  call V2mdlength(env)  !set the MD length according to a flexibility measure
  env%preactormtd=.false. !MTD settings already obtained
  env%preactorpot=.true.  !Generate logfermi potential 
  call reactor_setup(env)
  call tim%stop(1)

  call rdrcontrol('rcontrol',env)

  call tim%start(2,'test MD')
  call trialmd(env)    !calculate a short 1ps test MTD to check settings
  call tim%stop(2)



  !--- 2. run the metadynamics
     write(*,*)
     write(*,'(''========================================'')')
     write(*,'(''            MTD Simulations '',15x)')
     write(*,'(''========================================'')')
     call tim%start(3,'MTD')
     call MetaMD_para_OMP(env)
     call tim%stop(3)

  !--- 3. loose threshold optimization
     call tim%start(4,'geometry opt.')
     call append_INPUT_to('coord','input')   !include the input structure into the last optimization
     call multilevel_opt(env,2)
     call tim%stop(4)

     write(*,*)
     write(*,'(3x,''================================================'')')
     write(*,'(3x,''|           Final Geometry Optimization        |'')')
     write(*,'(3x,''================================================'')')
     if(env%doNMR) env%cgf(3)=.true. !--- if NMR equivalencies are requested, turn them on here
     call tim%start(4,'geometry opt.')
     call multilevel_opt(env,99)   !--- the last CREGEN is done within this subroutine
     call tim%stop(4)

  !--- 4. read the final ensemble and analyze



!---- print CREGEN results and clean up Directory a bit
      call V2terminating()

  return
end subroutine compress

