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

module calc_module

   use iso_fortran_env, only: wp => real64
   use strucrd
   use lj

   implicit none

   !=========================================================================================!
   !--- private module variables and parameters
   private
     integer :: i,j,k,l,ich,och,io
     logical :: ex

   !--- some constants and name mappings
     real(wp),parameter :: bohr  = 0.52917726_wp
     real(wp),parameter :: autokcal = 627.509541_wp


   public :: calcdata
   !=========================================================================================!
   !data object that contains settings for a molecular dynamics simulation.
   type :: calcdata

       integer :: id

       integer  :: molcharge    
       integer  :: uhf

       character(len=:),allocatable :: path
       character(len=:),allocatable :: other    
       character(len=:),allocatable :: systemcall

   end type calcdata


   public :: engrad
   public :: test_engrad

contains

subroutine engrad(n,xyz,at,dat,energy,gradient)
        implicit none
        integer,intent(in)  :: n
        real(wp),intent(in) :: xyz(3,n)
        integer,intent(in)  :: at(n)
        type(calcdata) :: dat

        real(wp),intent(inout) :: energy
        real(wp),intent(inout) :: gradient(3,n)
       
        real(wp) :: dum1,dum2


        dum1 = 1.0_wp
        dum2 = 1.0_wp

        select case (dat%id)
         case( 99 ) !-- Lennard-Jones dummy calculation
            call lj_engrad(n,xyz,dum1,dum2,energy,gradient)     
         case default
            write(*,*) 'Nothing selected for energy and gradient calculation.'
        end select


        return
end subroutine engrad        


subroutine test_engrad(fname)
        implicit none

        character(len=*) :: fname

        type(coord)    :: mol
        type(calcdata) :: dat

        real(wp) :: energy
        real(wp),allocatable :: grad(:,:)


        call mol%open(fname)

        allocate(grad(3,mol%nat), source = 0.0_wp)


        dat%id = 99
        call engrad(mol%nat,mol%xyz,mol%at,dat,energy,grad)
        
        write(*,'(3x,a,f10.8)') 'Energy: ',energy
        write(*,*) 'Gradient:'
        do i=1,mol%nat
           write(*,'(3f18.8)') grad(1:3,i)
        enddo


        deallocate(grad)

        return
end subroutine test_engrad


!==========================================================================================!        
end module calc_module
