!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022 Philipp Pracht
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
!> This piece of code generates a calcdata object from the global settings in env
subroutine env2calc(env,calc,molin)
   use crest_data 
   use calc_type
   use strucrd
   implicit none 
   !> INPUT
   type(systemdata),intent(in) :: env
   type(coord),intent(in),optional :: molin
   !> OUTPUT 
   type(calcdata) :: calc
   !> LOCAL
   type(calculation_settings) :: cal
   type(coord) :: mol

   call calc%reset()

   cal%uhf  = env%uhf
   cal%chrg = env%chrg
   !>-- obtain WBOs OFF by default
   cal%rdwbo = .false.

   !>-- defaults to whatever env has selected or gfn0
   select case(trim(env%gfnver))
     case( '--gfn0' )
       cal%id = jobtype%gfn0
     case( '--gfn1' )
       cal%id = jobtype%tblite
       cal%tblitelvl = 1
     case( '--gfn2' )
       cal%id = jobtype%tblite
       cal%tblitelvl = 2
     case default
       cal%id = jobtype%gfn0
   end select
   if(present(molin))then
     mol = molin
   else
     call mol%open('coord')
   endif

   call calc%add( cal )   

   return
end subroutine env2calc

!================================================================================!
subroutine confscript2i(env,tim)
  use iso_fortran_env,only:wp => real64
  use crest_data
  implicit none
  type(systemdata) :: env
  type(timer)   :: tim
  if(env%legacy)then
      call confscript2i_legacy(env,tim)
  else
      call crest_search_imtdgc(env,tim)
  endif
end subroutine confscript2i

!=================================================================================!
subroutine xtbsp(env,xtblevel)
  use iso_fortran_env,only:wp => real64
  use crest_data
  use strucrd, only: coord
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
  if(env%legacy)then
      call xtbsp_legacy(env,xtblevel)
  else
      call crest_xtbsp(env,xtblevel)
  endif
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
  if(env%legacy)then
      call xtbsp2_legacy(fname,env)
  else
      call mol%open(trim(fname))
      call crest_xtbsp(env,xtblevel=-1,molin=mol)
  endif
end subroutine xtbsp2

!=================================================================================!

subroutine confscript1(env,tim)
  use crest_parameters
  use crest_data
  implicit none
  type(systemdata) :: env
  type(timer)   :: tim
  write(stdout,*)
  write(stdout,*) 'This runtype has been entirely deprecated.'
  write(stdout,*) 'You may try an older version of the program if you want to use it.'
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
  if(env%legacy)then
    call nciflexi_legacy(env,flexval)
  else
    call env%ref%to(mol)
    call nciflexi_gfnff(mol,flexval)
  endif  
end subroutine nciflexi

!================================================================================!

subroutine thermo_wrap(env,pr,nat,at,xyz,dirname, &
        &  nt,temps,et,ht,gt,stot,bhess)
!******************************************
!* Wrapper for a Hessian calculation
!* to get thermodynamics of the molecule
!*****************************************
  use crest_parameters, only: wp
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

  if(env%legacy)then
    call thermo_wrap_legacy(env,pr,nat,at,xyz,dirname, &
    &                    nt,temps,et,ht,gt,stot,bhess)
  else
    error stop 'thermo_wrap not yet implemented for API version'
  endif
end subroutine thermo_wrap

!================================================================================!
