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
!> Implementation of a singlepoint calculation
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
!> Input/Output:
!>  env  -  crest's systemdata object
!>  tim  -  timer object
!>-----------------------------------------------
subroutine crest_singlepoint(env,tim)
  use crest_parameters
  use crest_data
  use strucrd 
  use calc_type
  use calc_module
  use optimize_module
  use tblite_api
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
   

  real(wp) :: energy
  real(wp),allocatable :: grad(:,:)
!========================================================================================!
  call tim%start(14,'singlepoint calc.') 
!========================================================================================!
  write(*,*)
  !call system('figlet singlepoint')
  write(*,*) "     _             _                  _       _   "
  write(*,*) " ___(_)_ __   __ _| | ___ _ __   ___ (_)_ __ | |_ "
  write(*,*) "/ __| | '_ \ / _` | |/ _ \ '_ \ / _ \| | '_ \| __|"
  write(*,*) "\__ \ | | | | (_| | |  __/ |_) | (_) | | | | | |_ "
  write(*,*) "|___/_|_| |_|\__, |_|\___| .__/ \___/|_|_| |_|\__|"
  write(*,*) "             |___/       |_|                      "
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
  write(*,*) 'etemp   ',calc%calcs(1)%etemp
  write(*,*) 'chrg    ',calc%calcs(1)%chrg
  write(*,*) 'uhf     ', calc%calcs(i)%uhf
  write(*,*) 'accuracy',calc%calcs(1)%accuracy
  write(*,*) 'maxscc  ',calc%calcs(1)%maxscc
  call engrad(mol,calc,energy,grad,io)
  write(*,*) 'iostatus',io  

   write(*,*)
   write (*,*) 'Energy: ',energy
   write (*,*) 'Gradient:'
   do i = 1,mol%nat
      write (*,'(3f18.8)') grad(1:3,i)
   end do
   write(*,*)
 
  deallocate(grad)
!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_singlepoint

!========================================================================================!
!========================================================================================!
!> replacement for the legacy xtbsp routine, makes use of gfn0 or tblite
!> the purpose of this routine is usually to generate WBOs
subroutine crest_xtbsp(env,xtblevel,molin)
   use crest_data
   use strucrd
   use calc_type
   use calc_module
   implicit none
   !> INPUT 
   type(systemdata) :: env
   integer,intent(in),optional :: xtblevel
   type(coord),intent(in),optional :: molin
   !> LOCAL
   integer :: lv,io
   type(calcdata) :: tmpcalc
   type(calculation_settings) :: cal
   type(coord) :: mol
   real(wp) :: energy
   real(wp),allocatable :: grad(:,:)

   call env2calc(env,tmpcalc,molin)

   tmpcalc%calcs(1)%rdwbo = .true. !> obtain WBOs
   if(present(xtblevel))then
     select case(xtblevel) !> no default
     case( 0 )
       tmpcalc%calcs(1)%id = jobtype%gfn0
     case( 1 )
       tmpcalc%calcs(1)%id = jobtype%tblite
       tmpcalc%calcs(1)%tblitelvl = 1
     case( 2 ) 
       tmpcalc%calcs(1)%id = jobtype%tblite
       tmpcalc%calcs(1)%tblitelvl = 2
     end select
   endif
   if(present(molin))then
     mol = molin
   else
     call mol%open('coord')
   endif
   allocate(grad(3,mol%nat), source = 0.0_wp)

   call engrad(mol,tmpcalc,energy,grad,io)
   if(io .ne. 0)then
     error stop 'crest_xtbsp failed'
   endif

   deallocate(grad)
   call tmpcalc%reset()
   call mol%deallocate()
end subroutine crest_xtbsp

