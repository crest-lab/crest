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

subroutine crest_refine(env,input,output)
!*******************************************************************
!* subroutine crest_refine
!* This subroutine will process the ensemble specified by "input"
!* and either overwrite it or write a new file "output".
!* The routine is intended to be called after geometry optimization
!* to re-rank structures by singlepoint energies, or add other
!* contributions to the energy.
!*******************************************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use cregen_interface
  use parallel_interface
  implicit none
  type(systemdata),intent(inout) :: env
  character(len=*),intent(in) :: input
  character(len=*),intent(in),optional :: output
!===========================================================!
  integer :: i,j,k,l,io,ich,m,t1,t2
  logical :: pr,wr,ex
!===========================================================!
  character(len=:),allocatable :: outname
  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)
  character(len=:),allocatable :: ensnam
  integer :: nat,nall
  real(wp),allocatable :: eread(:),etmp(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  integer :: nrefine,refine_stage
  type(coord),allocatable :: structures(:)
!===========================================================!
!>--- setup
  if (present(output)) then
    outname = output !> new file
  else
    outname = input  !> overwrite
  end if

  if(.not.allocated(env%refine_queue))then
    call rename(trim(input),trim(output))
    return
  endif

!>--- presorting step, if necessary
  if (env%refine_presort) then
    call newcregen(env,0,input)
    call rename('crest_ensemble.xyz',input)
  end if

!>--- read in
!  call rdensemble(input,nat,nall,at,xyz,eread)
!  allocate (etmp(nall),source=0.0_wp)
!!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!!>--- Important: crest_sploop requires coordinates in Bohrs
!  xyz = xyz/bohr
!!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
  call rdensembleparam(input,nat,nall)
  allocate (xyz(3,nat,nall),at(nat),etmp(nall),eread(nall))
  call rdensemble(input,nall,structures)
  at(:) = structures(1)%at(:)
  do j = 1,nall
    eread(j) = structures(j)%energy
    xyz(1:3,1:nat,j) = structures(j)%xyz(1:3,1:nat)
  end do

!===========================================================!
  DO_REFINE: if (allocated(env%refine_queue)) then
!===========================================================!

    call smallhead('ensemble refinement')

    nrefine = size(env%refine_queue,1)

    do i = 1,nrefine
      refine_stage = env%refine_queue(i)
      !> set the calculator to the correct stage
      env%calc%refine_stage = refine_stage

      select case (refine_stage)
      case (refine%singlepoint)
        write (stdout,'("> Singlepoint re-ranking for ",i0," structures")') nall
        do j = 1,nall
          structures(j)%xyz(1:3,1:nat) = xyz(1:3,1:nat,j)
        end do
        call crest_sploop(env,nall,structures,eread)

      case (refine%correction)
        write (stdout,'("> Additive correction for ",i0," structures")') nall
        do j = 1,nall
          structures(j)%xyz(1:3,1:nat) = xyz(1:3,1:nat,j)
        end do
        call crest_sploop(env,nall,structures,etmp)
        eread(:) = eread(:)+etmp(:)

      case (refine%geoopt)
        write (stdout,'("> Geometry optimization of ",i0," structures")') nall
        do j = 1,nall
          structures(j)%xyz(1:3,1:nat) = xyz(1:3,1:nat,j)
        end do
        call crest_oloop(env,nall,structures,.false.)
        do j = 1,nall
          xyz(1:3,1:nat,j) = structures(j)%xyz(1:3,1:nat)
          eread(j) = structures(j)%energy
        end do

      case (refine%deltaG)
        write (stdout,'("> Free energy correction (δG) for ",i0," structures")') nall
        call crest_hessloop(env,nat,nall,at,xyz,etmp)
        eread(:) = eread(:)+etmp(:)

      end select
      write (stdout,*)
    end do

    !> reset the refinement stage of the calculator
    env%calc%refine_stage = 0

!===========================================================!
  end if DO_REFINE
!===========================================================!

!>--- sync refined energies and coordinates back into coord structures
  do j = 1,nall
    structures(j)%energy = eread(j)
    structures(j)%xyz(1:3,1:nat) = xyz(1:3,1:nat,j)
    structures(j)%wrextxyz = .true.
  end do
!>--- write output ensemble in extxyz format
  call wrensemble(outname,nall,structures)

!===========================================================!
  deallocate (etmp,eread,xyz,at,structures)
  return
end subroutine crest_refine
!========================================================================================!
!========================================================================================!
