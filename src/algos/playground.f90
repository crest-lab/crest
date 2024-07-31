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

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> Implementation of whatever, for testing implementations
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
!> Input/Output:
!>  env  -  crest's systemdata object
!>  tim  -  timer object
!>-----------------------------------------------
subroutine crest_playground(env,tim)
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd 
  use canonical_mod
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich 
  logical :: pr,wr
!========================================================================================!
  type(calcdata) :: calc
  real(wp) :: accuracy,etemp
   
  integer :: V,maxgen
  integer,allocatable :: A(:,:)
  logical,allocatable :: rings(:,:)
  integer,allocatable :: tmp(:)
  logical :: connected,fail

  real(wp) :: energy
  real(wp),allocatable :: grad(:,:),geo(:,:),csv(:,:)

  type(canonical_sorter) ::  can
!========================================================================================!
  call tim%start(14,'Test implementation') 
!========================================================================================!
  !call system('figlet welcome')
  write(*,*) "              _                          "
  write(*,*) "__      _____| | ___ ___  _ __ ___   ___ "
  write(*,*) "\ \ /\ / / _ \ |/ __/ _ \| '_ ` _ \ / _ \"
  write(*,*) " \ V  V /  __/ | (_| (_) | | | | | |  __/"
  write(*,*) "  \_/\_/ \___|_|\___\___/|_| |_| |_|\___|"
  write(*,*) 
!========================================================================================!
  call env%ref%to(mol)
  write(*,*)
  write(*,*) 'Input structure:'
  call mol%append(stdout)
  write(*,*) 
!========================================================================================!

  allocate(grad(3,mol%nat), source=0.0_wp)
  call env2calc(env,calc,mol)
  calc%calcs(1)%rdwbo=.true.
  call calc%info(stdout)

  call engrad(mol,calc,energy,grad,io)
  call calculation_summary(calc,mol,energy,grad)
  

  write(stdout,*)
  write(stdout,*) 'CANGEN algorithm' 
  call can%init(mol,calc%calcs(1)%wbo)
  call can%rankprint(mol) 
  call can%stereo(mol)

!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_playground
