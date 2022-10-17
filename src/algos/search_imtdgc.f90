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

!> This is the re-implementation of a 
!> CREST's iMTD-GC default workflow

subroutine crest_search_imtdgc(env,tim)
  use crest_parameters, only: wp,stdout
  use crest_data
  use strucrd
  use calc_type
  use calc_module
  use dynamics_module
  use shake_module
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,m
  logical :: pr,wr
!===========================================================!
  type(calcdata) :: calc
  type(mddata) :: mddat
  type(shakedata) :: shk

  type(mddata),allocatable :: mddats(:)
  integer :: nsim

  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)
  character(len=:),allocatable :: ensnam
  integer :: nat,nall
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  logical :: dump,ex
  character(len=80) :: atmp
  logical :: multilevel(6)
!===========================================================!
!>--- printout header
  write (stdout,*)
  write (stdout,'(10x,"┍",49("━"),"┑")')
  write (stdout,'(10x,"│",14x,a,13x,"│")') "CREST iMTD-GC SAMPLING"
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

  call crest_search_multimd_init(env,mol,mddat,nsim)
  allocate (mddats(nsim), source=mddat)
  call crest_search_multimd_init2(env,mddats,nsim)

  call tim%start(2,'MD simulations')
  call crest_search_multimd(env,mol,mddats,nsim)
  call tim%stop(2)
  !>--- a file called crest_dynamics.trj should have been written
  ensnam = 'crest_dynamics.trj'

!==========================================================!
!>--- Reoptimization of trajectories
  call tim%start(3,'geom. optimization')
  multilevel = .true.
  multilevel(5) = .true.
  multilevel(2) = .true.
  call crest_multilevel_oloop(env,ensnam,multilevel)
  call tim%stop(3)


!==========================================================!
!>--- rename ensemble and sort
  call rename(ensemblefile,mecpensemble)
  call newcregen(env,12)
 

!==========================================================!
  return
end subroutine crest_search_imtdgc

!========================================================================================!

subroutine crest_multilevel_oloop(env,ensnam,multilevel)
  use crest_parameters, only: wp,stdout
  use crest_data
  use strucrd
  use calc_type
  use calc_module
  implicit none
  type(systemdata) :: env 
  character(len=*),intent(in) :: ensnam
  logical,intent(in) :: multilevel(6)
  integer :: nat,nall
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  logical :: dump
  character(len=128) :: inpnam,outnam
  integer :: i

  write(stdout,*)
  write(stdout,'(1x,a)') '--------------------------------'
  write(stdout,'(1x,a)') 'Multilevel Ensemble Optimization'
  write(stdout,'(1x,a)') '--------------------------------'

  !>--- read ensemble
  call rdensembleparam(ensnam,nat,nall)
  if (nall .lt. 1) then
    write(stdout,*) 'empty ensemble file',trim(ensnam)
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

!>--- sequential optimizations of ensembles
  dump = .true. !> optimized structures will be written to crest_ensemble.xyz
  do i=6,1,-1
    if(multilevel(i))then
     !>--- set threads
       call ompautoset(env%threads,7,env%omp,env%MAXRUN,nall)
     !>--- run parallel optimizations
       call crest_oloop(env,nat,nall,at,xyz,eread,dump)
     !>--- rename ensemble and sort
       call checkname_xyz(crefile,inpnam,outnam)
       call rename(ensemblefile,trim(inpnam))
       call sort_and_check(env,trim(inpnam))
       call checkname_xyz(crefile,inpnam,outnam)
     !>--- read new ensemble for next iteration
       deallocate(eread,at,xyz)
       call rdensembleparam(trim(inpnam),nat,nall)
       if (nall .lt. 1) then
         write(stdout,*) 'empty ensemble file',trim(inpnam)
         stop
       endif
       allocate (xyz(3,nat,nall),at(nat),eread(nall))
       call rdensemble(trim(inpnam),nat,nall,at,xyz,eread)
    endif
  enddo




  if(allocated(eread)) deallocate(eread)
  if(allocated(at))  deallocate(at)
  if(allocated(xyz)) deallocate(xyz)
  return
end subroutine crest_multilevel_oloop

