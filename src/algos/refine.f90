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
  implicit none
  type(systemdata),intent(inout) :: env
  character(len=*),intent(in) :: input
  character(len=*),intent(in),optional :: output
!===========================================================!
  integer :: i,j,k,l,io,ich,m
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
!===========================================================!
!>--- setup
  if (present(output)) then
    outname = output !> new file
  else
    outname = input  !> overwrite
  end if
 
!>--- presorting step, if necessary
  if(env%refine_presort)then
    call newcregen(env,0,input)
    call rename('crest_ensemble.xyz',input)
  endif

!>--- read in
  call rdensemble(input,nat,nall,at,xyz,eread)
  allocate (etmp(nall),source=0.0_wp)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: crest_sploop requires coordinates in Bohrs
    xyz = xyz / bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

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
        call crest_sploop(env,nat,nall,at,xyz,eread)

      case (refine%correction)
        write (stdout,'("> Additive correction for ",i0," structures")') nall
        call crest_sploop(env,nat,nall,at,xyz,etmp)
        eread(:) = eread(:)+etmp(:)

      case (refine%geoopt)
        write (stdout,'("> Geometry optimization of ",i0," structures")') nall
        call crest_oloop(env,nat,nall,at,xyz,eread,.false.)

      case(refine%confsolv)
        write (stdout,'("> ConfSolv: ΔΔGsoln estimation from 3D directed message passing neural networks (D-MPNN)")')
        call confsolv_request( input, nall, env%threads, etmp, io)
        if(io == 0)then
        eread(:) = etmp(:)*kcaltoau  !> since CREGEN deals with Eh energies
        endif   
      end select
      write(stdout,*) 
    end do

    !> reset the refinement stage of the calculator
    env%calc%refine_stage = 0

!===========================================================!
  end if DO_REFINE
!===========================================================!

!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: ensemble file must be written in AA
  xyz = xyz / angstrom
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- write output ensemble
  call wrensemble(outname,nat,nall,at,xyz,eread)

!===========================================================!
  deallocate (etmp,eread,xyz,at)
  return
end subroutine crest_refine
!========================================================================================!
!========================================================================================!
