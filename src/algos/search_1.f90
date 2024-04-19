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

!> This is an implementation of a generic sampling workflow
!> consisting of one batch of metadynamics with snapshot optimization.

subroutine crest_search_1(env,tim)
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol
!===========================================================!
  type(mddata) :: mddat

  type(mddata),allocatable :: mddats(:)
  integer :: nsim

  character(len=:),allocatable :: ensnam
  integer :: nat,nall,T,Tn
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  logical :: dump

!===========================================================!
!>--- printout header
  write (stdout,*)
  write (stdout,'(10x,"┍",49("━"),"┑")')
  write (stdout,'(10x,"│",13x,a,12x,"│")') "CREST SAMPLING ALGORITHM"
  write (stdout,'(10x,"┕",49("━"),"┙")')
  write (stdout,*)

!===========================================================!
!>--- setup
  call env%ref%to(mol)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)

!===========================================================!
!>--- Dynamics

  write(stdout,*)
  write(stdout,'(1x,a)') '------------------------------'
  write(stdout,'(1x,a)') 'Molecular Dynamics Simulations'
  write(stdout,'(1x,a)') '------------------------------'

  nsim = -1 !>--- enambles automatic MTD setup in init routines
  call crest_search_multimd_init(env,mol,mddat,nsim)
  allocate (mddats(nsim), source=mddat)
  call crest_search_multimd_init2(env,mddats,nsim)

  call tim%start(2,'Molecular dynamics (MD)')
  call crest_search_multimd(env,mol,mddats,nsim)
  call tim%stop(2)
!>--- a file called crest_dynamics.trj should have been written
  ensnam = 'crest_dynamics.trj'

!==========================================================!
!>--- Reoptimization of trajectories

  write(stdout,*)
  write(stdout,'(1x,a)') '---------------------'
  write(stdout,'(1x,a)') 'Ensemble Optimization'
  write(stdout,'(1x,a)') '---------------------'

!>--- read ensemble
  call rdensembleparam(ensnam,nat,nall)
  if (nall .lt. 1) then
    write(stdout,*) 'empty ensemble file'
    return
  endif
  allocate (xyz(3,nat,nall),at(nat),eread(nall))
  call rdensemble(ensnam,nat,nall,at,xyz,eread)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: crest_oloop requires coordinates in Bohrs
  xyz = xyz / bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
 
  write(stdout,'(1x,a,i0,a,a,a)')'Optimizing all ',nall, &
  & ' structures from file "',trim(ensnam),'" ...'

!>--- set threads
  call new_ompautoset(env,'auto',nall,T,Tn) 

!>--- optimize
  call tim%start(3,'Geometry optimization')
  dump = .true.
  call crest_oloop(env,nat,nall,at,xyz,eread,dump)
  call tim%stop(3)

!==========================================================!
  return
end subroutine crest_search_1

