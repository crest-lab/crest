!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Philipp Pracht
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

!================================================================================!
subroutine env2calc(env,calc,molin)
!******************************************
!* This piece of code generates a calcdata
!* object from the global settings in env
!******************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use iomod
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env
  type(coord),intent(in),optional :: molin
  !> OUTPUT
  type(calcdata) :: calc
  !> LOCAL
  type(calculation_settings) :: cal,cal2
  type(coord) :: mol

!>--- Calculator level

!>-- defaults to whatever env has selected or gfn0
  call cal%create(trim(env%gfnver))
  if (present(molin)) then
    mol = molin
  end if

  cal%uhf = env%uhf
  cal%chrg = env%chrg
!>-- obtain WBOs OFF by default
  cal%rdwbo = .false.
  cal%rddip = .true.
  !> except for SP runtype (from command line!)
  if( env%crestver == crest_sp )then
    cal%rdwbo = .true.
    cal%rddip = .true.
  endif

  !> implicit solvation
  if (env%gbsa) then
    if (index(env%solv,'gbsa') .ne. 0) then
      cal%solvmodel = 'gbsa'
    else if (index(env%solv,'alpb') .ne. 0) then
      cal%solvmodel = 'alpb'
    else
      cal%solvmodel = 'unknown'
    end if
    cal%solvent = trim(env%solvent)
  end if

  !> do not reset parameters between calculations (opt for speed)
  cal%apiclean = .false.

  call cal%autocomplete(1)
  call calc%add(cal)

!>--- Refinement level
  if (trim(env%gfnver2) .ne. '') then
    env%gfnver2 = lowercase(env%gfnver2)
    call cal2%create(trim(env%gfnver2))

    cal2%chrg = cal%chrg
    cal2%uhf = cal%uhf
    if (env%gbsa) then
      cal2%solvmodel = cal%solvmodel
      cal2%solvent = cal%solvent
    end if

    call cal2%autocomplete(2)
        
    cal2%refine_lvl = refine%singlepoint
    call calc%add(cal2)
    if(allocated(env%refine_queue)) deallocate(env%refine_queue)
    call env%addrefine( refine%singlepoint )  
  end if

  return
end subroutine env2calc

subroutine env2calc_setup(env)
!***********************************
!* Setup the calc object within env
!* (wrapper to get the mol object)
!***********************************
  use crest_data
  use crest_calculator
  use strucrd
  implicit none
  !> INOUT
  type(systemdata),intent(inout) :: env
  !> LOCAL
  type(calcdata) :: calc
  type(coord) :: mol
  interface
    subroutine env2calc(env,calc,molin)
      use crest_parameters
      use crest_data
      use crest_calculator
      use strucrd
      implicit none
      type(systemdata),intent(inout) :: env
      type(coord),intent(in),optional :: molin
      type(calcdata) :: calc
    end subroutine env2calc

  end interface

  call env%ref%to(mol) 

  call env2calc(env,env%calc,mol)

 ! env%calc = calc
end subroutine env2calc_setup

!================================================================================!
subroutine confscript2i(env,tim)
  use iso_fortran_env,only:wp => real64
  use crest_data
  implicit none
  type(systemdata) :: env
  type(timer)   :: tim
  if (env%legacy) then
    call confscript2i_legacy(env,tim)
  else
    if (.not.env%entropic) then
      call crest_search_imtdgc(env,tim)
    else
      call crest_search_entropy(env,tim)
    end if
  end if
end subroutine confscript2i

!================================================================================!
subroutine mdopt(env,tim)
  use iso_fortran_env,only:wp => real64
  use crest_data
  implicit none
  type(systemdata) :: env
  type(timer)   :: tim
  if (env%legacy) then
    call mdopt_legacy(env,tim)
  else
    call crest_ensemble_optimization(env,tim)
  end if
end subroutine mdopt

!================================================================================!
subroutine screen(env,tim)
  use iso_fortran_env,only:wp => real64
  use crest_data
  implicit none
  type(systemdata) :: env
  type(timer)   :: tim
  if (env%legacy) then
    call screen_legacy(env,tim)
  else
    call crest_ensemble_screening(env,tim)
  end if
end subroutine screen

!=================================================================================!

subroutine xtbsp(env,xtblevel)
  use iso_fortran_env,only:wp => real64
  use crest_data
  use strucrd,only:coord
  implicit none
  type(systemdata) :: env
  integer,intent(in),optional :: xtblevel
  interface
    subroutine crest_xtbsp(env,xtblevel,molin)
      import :: systemdata,coord
      type(systemdata) :: env
      integer,intent(in),optional :: xtblevel
      type(coord),intent(in),optional :: molin
    end subroutine crest_xtbsp
  end interface
  if (env%legacy) then
    call xtbsp_legacy(env,xtblevel)
  else
    call crest_xtbsp(env,xtblevel)
  end if
end subroutine xtbsp
subroutine xtbsp2(fname,env)
  use iso_fortran_env,only:wp => real64
  use crest_data
  use strucrd
  implicit none
  type(systemdata) :: env
  character(len=*),intent(in) :: fname
  type(coord) :: mol
  interface
    subroutine crest_xtbsp(env,xtblevel,molin)
      import :: systemdata,coord
      type(systemdata) :: env
      integer,intent(in),optional :: xtblevel
      type(coord),intent(in),optional :: molin
    end subroutine crest_xtbsp
  end interface
  if (env%legacy) then
    call xtbsp2_legacy(fname,env)
  else
    call mol%open(trim(fname))
    call crest_xtbsp(env,xtblevel=-1,molin=mol)
  end if
end subroutine xtbsp2

!=================================================================================!

subroutine confscript1(env,tim)
  use crest_parameters
  use crest_data
  implicit none
  type(systemdata) :: env
  type(timer)   :: tim
  write (stdout,*)
  write (stdout,*) 'This runtype has been entirely deprecated.'
  write (stdout,*) 'You may try an older version of the program if you want to use it.'
  stop
end subroutine confscript1

!=================================================================================!

subroutine nciflexi(env,flexval)
  use crest_parameters
  use crest_data
  use strucrd
  implicit none
  type(systemdata) :: env
  type(coord) ::mol
  real(wp) :: flexval
  if (env%legacy) then
    call nciflexi_legacy(env,flexval)
  else
    call env%ref%to(mol)
    call nciflexi_gfnff(mol,flexval)
  end if
end subroutine nciflexi

!================================================================================!

subroutine thermo_wrap(env,pr,nat,at,xyz,dirname, &
        &  nt,temps,et,ht,gt,stot,bhess)
!******************************************
!* Wrapper for a Hessian calculation
!* to get thermodynamics of the molecule
!*****************************************
  use crest_parameters,only:wp
  use crest_data
  implicit none
  !> INPUT
  type(systemdata) :: env
  logical,intent(in) :: pr
  integer,intent(in) :: nat
  integer,intent(inout) :: at(nat)
  real(wp),intent(inout) :: xyz(3,nat)  !> in Angstroem!
  character(len=*) :: dirname
  integer,intent(in)  :: nt
  real(wp),intent(in)  :: temps(nt)
  logical,intent(in) :: bhess       !> calculate bhess instead?
  !> OUTPUT
  real(wp),intent(out) :: et(nt)    !> enthalpy in Eh
  real(wp),intent(out) :: ht(nt)    !> enthalpy in Eh
  real(wp),intent(out) :: gt(nt)    !> free energy in Eh
  real(wp),intent(out) :: stot(nt)  !> entropy in cal/molK

  if (env%legacy) then
    call thermo_wrap_legacy(env,pr,nat,at,xyz,dirname, &
    &                    nt,temps,et,ht,gt,stot,bhess)
  else
    call thermo_wrap_new(env,pr,nat,at,xyz,dirname, &
    &                    nt,temps,et,ht,gt,stot,bhess)
  end if
end subroutine thermo_wrap

!================================================================================!

subroutine trialMD(env)
!***********************************************
!* subroutine trialMD
!* Takes the global metadynamics settings
!* And performs a short 1 ps simulation to
!* check if the molecular dynamis/metadynamics
!* will run, or if the timestep is too large
!***********************************************
  use crest_parameters,only:wp
  use crest_data
  implicit none
  !> INPUT
  type(systemdata) :: env

  if (env%legacy) then
    !> old xtb subprocess version
    call trialMD_legacy(env)
  else
    !> new calculator implementation
    call trialMD_calculator(env)
  end if

end subroutine trialMD

!================================================================================!

subroutine trialOPT(env)
!**********************************************************
!* subroutine trialOPT
!* Performs a geometry optimization of the structure
!* saved to env%ref and checks for changes in the topology
!**********************************************************
  use crest_data
  use crest_parameters, only: stdout
  implicit none
  !> INPUT
  type(systemdata) :: env

  if (env%legacy) then
    call xtbopt_legacy(env)
  else
    call trialOPT_calculator(env)
  end if

  if(env%crestver == crest_trialopt)then
!>-- if we reach this point in the standalone trialopt the geometry is ok!
   write(stdout,*)
   stop 'Geometry ok!'
  endif
end subroutine trialOPT

!================================================================================!

