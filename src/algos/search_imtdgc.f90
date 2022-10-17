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
  character(len=80) :: atmp,btmp,str
  logical :: multilevel(6)
  logical :: start,lower
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
!>--- Start mainloop 
  env%nreset = 0
  start = .true.
  MAINLOOP : do
    call printiter
    if (.not. start) then
      call clean_V2i  !--- clean Dir for new iterations
      env%nreset = env%nreset + 1
    else !--at the beginning clean existing backup-ensembles
      call rmrfw('.cre_')
    end if
!===========================================================!
!>--- Meta-dynamics loop
  mtdloop: do i = 1,env%Maxrestart

    write(stdout,*)
    write(stdout,'(1x,a)') '------------------------------'
    write(stdout,'(1x,a,i0)') 'Meta-Dynamics Iteration ',i
    write(stdout,'(1x,a)') '------------------------------'

    nsim = -1 !>--- enambles automatic MTD setup in init routines
    call crest_search_multimd_init(env,mol,mddat,nsim)
    allocate (mddats(nsim), source=mddat)
    call crest_search_multimd_init2(env,mddats,nsim)

    call tim%start(2,'MTD simulations')
    call crest_search_multimd(env,mol,mddats,nsim)
    call tim%stop(2)
    !>--- a file called crest_dynamics.trj should have been written
    ensnam = 'crest_dynamics.trj'
    !>--- deallocate for next iteration
    if(allocated(mddats))deallocate(mddats)

!==========================================================!
!>--- Reoptimization of trajectories
    call tim%start(3,'geom. optimization')
    multilevel = .false.
    multilevel(5) = .true.
    multilevel(2) = .true.
    call crest_multilevel_oloop(env,ensnam,multilevel)
    call tim%stop(3)

    !>--- save the CRE under a backup name
    call checkname_xyz(crefile,atmp,str)
    call checkname_xyz('.cre',str,btmp)
    call rename(atmp,btmp)
    !>--- save cregen output
    call checkname_tmp('cregen',atmp,btmp)
    call rename('cregen.out.tmp',btmp)

!=========================================================!
!>--- cleanup after first iteration and prepare next
    if (i .eq. 1 .and. start) then
      start = .false.
      !>-- obtain a first lowest energy as reference
      env%eprivious = env%elowest
      !>-- remove the two extreme-value MTDs
      if (.not. env%readbias .and.  env%runver .ne. 33 .and. &
      &   env%runver .ne. 787878 ) then
        env%nmetadyn = env%nmetadyn - 2
      end if
      !>-- the cleanup 
      call clean_V2i   
      !>-- and always do two cycles of MTDs
      cycle mtdloop 
    endif
!=========================================================!
!>--- Check for lowest energy
    call elowcheck(lower,env)
    if (.not. lower) then
      exit mtdloop
    end if
  enddo mtdloop

!==========================================================!
!>--- exit mainloop
   exit MAINLOOP
  enddo MAINLOOP

!==========================================================!
!>--- rename ensemble and sort
  call rename(ensemblefile,mecpensemble)
  call newcregen(env,12)
 

!==========================================================!
!>--- print CREGEN results and clean up Directory a bit
    call V2terminating()

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
       !>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
       !>--- Important: crest_oloop requires coordinates in Bohrs
       xyz = xyz / bohr
       !>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
    endif
  enddo




  if(allocated(eread)) deallocate(eread)
  if(allocated(at))  deallocate(at)
  if(allocated(xyz)) deallocate(xyz)
  return
end subroutine crest_multilevel_oloop

