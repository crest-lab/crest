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
  type(wavefunction_type) :: wfn  
  type(tblite_calculator) :: tbcalc  
  type(tblite_ctx)        :: ctx
  real(wp) :: accuracy,etemp
   
  integer :: V
  integer,allocatable :: A(:,:)
  logical,allocatable :: rings(:,:)
  integer,allocatable :: tmp(:)
  logical :: connected,fail
  type(zmolecule) :: zmol

  real(wp) :: energy
  real(wp),allocatable :: grad(:,:),geo(:,:)
  integer,allocatable :: na(:),nb(:),nc(:)
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

  allocate(grad(3,mol%nat),source=0.0_wp)
  calc = env%calc
  

  write(*,*) 'job type',calc%calcs(1)%id
  write(*,*) 'etemp',calc%calcs(1)%etemp
  write(*,*) 'chrg',calc%calcs(1)%chrg,'uhf', calc%calcs(i)%uhf
  write(*,*) 'accuracy',calc%calcs(1)%accuracy
  write(*,*) 'maxscc',calc%calcs(1)%maxscc
  call engrad(mol,calc,energy,grad,io)
  write(*,*) 'iostatus',io  

  call write_wbo(calc%calcs(1)%wbo, cutoff=0.05_wp)

   if(calc%calcs(1)%rdwbo .and. allocated(calc%calcs(1)%wbo))then
   write(*,*)
   write(*,*) 'WBOs:'
   do i=1,mol%nat
     do j=i+1,mol%nat
       if(calc%calcs(1)%wbo(i,j) .gt. 0.05_wp)then
         write(*,*),i,j,calc%calcs(1)%wbo(i,j)
       endif
     enddo
   enddo
   endif
   write(*,*)
   write (*,*) 'Energy: ',energy
   write (*,*) 'Gradient:'
   do i = 1,mol%nat
      write (*,'(3f18.8)') grad(1:3,i)
   end do
   write(*,*)


  deallocate(grad)

!========================================================================================!

  V = mol%nat
  allocate(A(V,V),na(V),nb(V),nc(V), source = 0)
  call wbo2adjacency(V,calc%calcs(1)%wbo,A,0.02_wp)
  allocate(geo(3,V), source = 0.0_wp)
  call BETTER_XYZINT(mol%nat,mol%xyz,A,NA,NB,NC,geo)
  !do i=1,mol%nat 
  !   write(stdout,*) i,'=>',na(i),nb(i),nc(i)
  !enddo

  !call XYZGEO2(mol%nat,mol%xyz,NA,NB,NC,radtodeg,GEO)
  call print_zmat(6,mol%nat,geo,NA,NB,NC) 

  call GMETRY2(mol%nat,geo,mol%xyz,na,nb,nc,fail)

  call mol%write('test.xyz')

  call XYZGEO2(mol%nat,mol%xyz,NA,NB,NC,1.0_wp,GEO)
  call print_zmat(6,mol%nat,geo,NA,NB,NC)


!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_playground
