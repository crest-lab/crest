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
  !use iso_fortran_env,only: wp =>real64,stdout=>output_unit
  use crest_parameters
  use crest_data
  use strucrd 
  use calc_type
  use calc_module
  use optimize_module
  use tblite_api
  use wiberg_mayer, only: write_wbo 
  use adjacency
  use zdata
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
  type(zmolecule) :: zmol

  real(wp) :: energy
  real(wp),allocatable :: grad(:,:),geo(:,:)
  integer,allocatable :: na(:),nb(:),nc(:),at2(:)
  integer :: nat2
!========================================================================================!
  call tim%start(14,'test implementation') 
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

  call test_vonMises(env)


!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_playground
