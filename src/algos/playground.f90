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
  use probabilities_module
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
  real(wp),allocatable :: mu(:)
  real(wp) :: kappa,rrad
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

   
!  allocate(mu(30000))
!  do i=1,30000 
!    call random_number(rrad)
!    rrad = (rrad-0.5_wp) * 10.0_wp * degtorad
!    mu(i) = rrad
!  enddo 
!  mu(1:10000) = mu(1:10000) !+ 120.0_wp*degtorad
!  mu(10001:20000) = mu(10001:20000) + 120.0_wp*degtorad
!  mu(20001:30000) = mu(20001:30000) + 240.0_wp*degtorad
!  !mu(1:1500) = mu(1:1500) + 210.0_wp*degtorad
!  !mu(1501:3000) = mu(1501:3000) + 240.0_wp*degtorad
!  kappa = 1.5_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','111-1.5-vonmises.txt')
!
!  kappa = 1.0_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','111-1.0-vonmises.txt')
!  
!  kappa = 0.5_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','111-0.5-vonmises.txt')
!
!  kappa = 2.0_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','111-2.0-vonmises.txt')
!
!  kappa = 5.0_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','111-5.0-vonmises.txt')
!
!  deallocate(mu)
!  allocate(mu(3000))
!  do i=1,3000
!    call random_number(rrad)
!    rrad = (rrad-0.5_wp) * 10.0_wp * degtorad
!    mu(i) = rrad
!  enddo
!  mu(1:1000) = mu(1:1000) !+ 120.0_wp*degtorad
!  mu(1001:2000) = mu(1001:2000) + 120.0_wp*degtorad
!  mu(2001:3000) = mu(2001:3000) + 240.0_wp*degtorad
!  !mu(1:1500) = mu(1:1500) + 210.0_wp*degtorad
!  !mu(1501:3000) = mu(1501:3000) + 240.0_wp*degtorad
!  kappa = 1.5_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','3000-111-1.5-vonmises.txt')
!
!  kappa = 0.5_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','3000-111-0.5-vonmises.txt')
!
!  kappa = 5.0_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','3000-111-5.0-vonmises.txt')
!
!
!  deallocate(mu)
!  allocate(mu(300))
!  do i=1,300
!    call random_number(rrad)
!    rrad = (rrad-0.5_wp) * 10.0_wp * degtorad
!    mu(i) = rrad
!  enddo
!  mu(1:100) = mu(1:100) !+ 120.0_wp*degtorad
!  mu(101:200) = mu(101:200) + 120.0_wp*degtorad
!  mu(201:300) = mu(201:300) + 240.0_wp*degtorad
!  !mu(1:1500) = mu(1:1500) + 210.0_wp*degtorad
!  !mu(1501:3000) = mu(1501:3000) + 240.0_wp*degtorad
!  kappa = 1.5_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','300-111-1.5-vonmises.txt')
!
!  kappa = 0.5_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','300-111-0.5-vonmises.txt')
!
!  kappa = 5.0_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','300-111-5.0-vonmises.txt')
!
!
!
!  deallocate(mu)
!  allocate(mu(30000))
!  do i=1,30000 
!    call random_number(rrad)
!    rrad = (rrad-0.5_wp) * 10.0_wp * degtorad
!    mu(i) = rrad
!  enddo 
!
!  mu(1:5000) = mu(1:5000) !+ 120.0_wp*degtorad
!  mu(5001:10000) = mu(5001:10000) + 60.0_wp*degtorad
!  mu(10001:15000) = mu(10001:15000) + 120.0_wp*degtorad
!  mu(15001:20000) = mu(15001:20000) + 180.0_wp*degtorad
!  mu(20001:25000) = mu(20001:25000) + 240.0_wp*degtorad
!  mu(25001:30000) = mu(25001:30000) + 300.0_wp*degtorad
!
!  kappa = 1.5_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','30k-6-1.5-vonmises.txt')
!
!  kappa = 0.5_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','30k-6-0.5-vonmises.txt')
!
!  kappa = 5.0_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','30k-6-5.0-vonmises.txt')
!
!  deallocate(mu)
!  allocate(mu(3000))
!  do i=1,3000
!    call random_number(rrad)
!    rrad = (rrad-0.5_wp) * 10.0_wp * degtorad
!    mu(i) = rrad
!  enddo
!  mu(1:500) = mu(1:500) !+ 120.0_wp*degtorad
!  mu(501:1000) = mu(501:1000) + 60.0_wp*degtorad
!  mu(1001:1500) = mu(1001:1500) + 120.0_wp*degtorad
!  mu(1501:2000) = mu(1501:2000) + 180.0_wp*degtorad
!  mu(2001:2500) = mu(2001:2500) + 240.0_wp*degtorad
!  mu(2501:3000) = mu(2501:3000) + 300.0_wp*degtorad
!
!  kappa = 1.5_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','3000-6-1.5-vonmises.txt')
!
!  kappa = 0.5_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','3000-6-0.5-vonmises.txt')
!
!  kappa = 5.0_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','3000-6-5.0-vonmises.txt')
!
!
!  deallocate(mu)
!  allocate(mu(300))
!  do i=1,300
!    call random_number(rrad)
!    rrad = (rrad-0.5_wp) * 10.0_wp * degtorad
!    mu(i) = rrad
!  enddo
!  mu(1:50) = mu(1:50) !+ 120.0_wp*degtorad
!  mu(51:100) = mu(51:100) + 60.0_wp*degtorad
!  mu(101:150) = mu(101:150) + 120.0_wp*degtorad
!  mu(151:200) = mu(151:200) + 180.0_wp*degtorad
!  mu(201:250) = mu(201:250) + 240.0_wp*degtorad
!  mu(251:300) = mu(251:300) + 300.0_wp*degtorad
!  kappa = 1.5_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','300-6-1.5-vonmises.txt')
!
!  kappa = 0.5_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','300-6-0.5-vonmises.txt')
!
!  kappa = 5.0_wp
!  call vonMises_plot(kappa,mu,3600)
!  call rename('vonmises.txt','300-6-5.0-vonmises.txt')




  if(allocated(mu)) deallocate(mu)
  allocate(mu(3000))
  do i=1,3000
    call random_number(rrad)
    rrad = (rrad-0.5_wp) * 10.0_wp * degtorad
    mu(i) = rrad
  enddo
  mu(1:2000) = mu(1:2000) 
  mu(2001:3000) = mu(2001:3000) + 30.0_wp*degtorad

  kappa = 1.5_wp
  call vonMises_plot(kappa,mu,3600)
  call rename('vonmises.txt','3000-21-1.5-vonmises.txt')

  kappa = 0.5_wp
  call vonMises_plot(kappa,mu,3600)
  call rename('vonmises.txt','3000-21-0.5-vonmises.txt')

  kappa = 5.0_wp
  call vonMises_plot(kappa,mu,3600)
  call rename('vonmises.txt','3000-21-5.0-vonmises.txt')


  if(allocated(mu)) deallocate(mu)
  allocate(mu(30000))
  do i=1,30000
    call random_number(rrad)
    rrad = (rrad-0.5_wp) * 10.0_wp * degtorad
    mu(i) = rrad
  enddo
  mu(1:20000) = mu(1:20000) 
  mu(20001:30000) = mu(20001:30000) + 30.0_wp*degtorad

  kappa = 1.5_wp
  call vonMises_plot(kappa,mu,3600)
  call rename('vonmises.txt','30k-21-1.5-vonmises.txt')

  kappa = 0.5_wp
  call vonMises_plot(kappa,mu,3600)
  call rename('vonmises.txt','30k-21-0.5-vonmises.txt')

  kappa = 5.0_wp
  call vonMises_plot(kappa,mu,3600)
  call rename('vonmises.txt','30k-21-5.0-vonmises.txt')




!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_playground
