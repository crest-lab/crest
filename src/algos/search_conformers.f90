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

subroutine crest_search_imtdgc(env,tim)
!*******************************************************************
!* This is the re-implementation of CREST's iMTD-GC default workflow
!*
!* Compared to the legacy implementation, this version
!* is separated from the entropy algo to keep things clean
!* The entropy algo (sMTD-iMTD) can be found in search_entropy.f90
!*******************************************************************
  use crest_parameters,only:wp,stdout
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  use iomod
  use utilities
  use cregen_interface
  use crest_restartlog
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,m
  logical :: pr,wr,doreturn
!===========================================================!
  type(calcdata) :: calc
  type(mddata) :: mddat
  type(shakedata) :: shk

  type(mddata),allocatable :: mddats(:)
  integer :: nsim,nallout

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
  type(restart_data) :: rdat
  logical :: do_restart,skip_mtdloop,skip_collect,firstiter,fex
!===========================================================!
!>--- printout header
  write (stdout,*)
  write (stdout,'(10x,"┍",49("━"),"┑")')
  write (stdout,'(10x,"│",14x,a,13x,"│")') "CREST iMTD-GC SAMPLING"
  write (stdout,'(10x,"┕",49("━"),"┙")')
  write (stdout,*)

! ── restart detection ─────────────────────────────────────────────
  do_restart = .false.
  skip_mtdloop = .false.
  skip_collect = .false.
  if (env%allowrestart .and. restart_file_exists()) then
    call read_restart_log(rdat)
    if (rdat%runtype == crest_imtd .and. rdat%stage /= 'done') then
      do_restart = .true.
      call print_restart_info(rdat)
      !> skip entire mtdloop and collectcre only when past the MTD loop
      skip_mtdloop = (rdat%stage == 'post_collect')
      skip_collect = (rdat%stage == 'post_collect')
    end if
  end if

!===========================================================!
!>--- setup
  call env%ref%to(mol)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)

!>--- saftey terminations
  call crest_sampling_skip(env,doreturn)
  if (doreturn) return

!>--- sets the MD length according to a flexibility measure
  call md_length_setup(env)
!>--- create the MD calculator saved to env
  call env_to_mddat(env)

  if (env%performMTD) then
!>--- (optional) calculate a short 1ps test MTD to check settings
    call tim%start(1,'Trial metadynamics (MTD)')
    call trialmd(env)
    call tim%stop(1)
    if (env%iostatus_meta .ne. 0) return
  end if

!===========================================================!
!>--- Start mainloop
  env%nreset = 0
  start = .true.
! ── apply restart state ───────────────────────────────────────────
  if (do_restart) then
    env%nreset   = rdat%main_iter
    env%elowest  = rdat%elowest
    env%eprivious = rdat%eprivious
    env%nmetadyn = rdat%nmetadyn
    start = .false.
! ── restore lowest structure as reference geometry ─────────────────
    call restart_restore_reference(env,rdat%last_file)
    call env%ref%to(mol)
  end if
  MAINLOOP: do
    call printiter
    if (do_restart) then
!>--- restart: preserve .cre_*.xyz files, skip cleanup
      continue
    else if (.not.start) then
!>--- clean Dir for new iterations, but leave iteration backup files
      call clean_V2i
      env%nreset = env%nreset+1
    else
!>--- at the beginning, wipe directory clean
      call V2cleanup(.false.)
    end if
!===========================================================!
!>--- Meta-dynamics loop (skipped on restart to use existing .cre_*.xyz)
    if (.not.skip_mtdloop) then
    mtdloop: do i = 1,env%Maxrestart

! ── restart: skip based on stage ──────────────────────────────────
      if (do_restart) then
        if (rdat%stage == 'mtd_loop' .and. i <= rdat%mtd_iter) cycle mtdloop
        if (rdat%stage == 'mtd_trj'  .and. i <  rdat%mtd_iter) cycle mtdloop
      end if

      write (stdout,*)
      write (stdout,'(1x,a)') '------------------------------'
      write (stdout,'(1x,a,i0)') 'Meta-Dynamics Iteration ',i
      write (stdout,'(1x,a)') '------------------------------'

!==========================================================!
!>--- MTD run (skipped for mtd_trj restart: trajectory already exists)
      if (do_restart .and. i == rdat%mtd_iter .and. &
        &  rdat%stage == 'mtd_trj') then
        write (stdout,'(1x,a,i0,a)') 'Restarting iteration ',i, &
          & ' from existing trajectory/ensemble'
        ensnam = trim(rdat%last_file)
      else
        nsim = -1 !>--- enambles automatic MTD setup in init routines
        call crest_search_multimd_init(env,mol,mddat,nsim)
        allocate (mddats(nsim),source=mddat)
        call crest_search_multimd_init2(env,mddats,nsim)

        call tim%start(2,'Metadynamics (MTD)')
        call crest_search_multimd(env,mol,mddats,nsim)
        call tim%stop(2)
!>--- a file called crest_dynamics.trj.xyz should have been written
        ensnam = 'crest_dynamics.trj.xyz'
        if (allocated(mddats)) deallocate (mddats)
!>--- checkpoint: trajectory ready, optimization about to start
        call write_restart_log(crest_imtd,'mtd_trj',env%nreset,i, &
          &  env%nmetadyn,env%elowest,env%eprivious,ensnam)
      end if

!==========================================================!
!>--- Reoptimization of trajectories
      call tim%start(3,'Geometry optimization')
      call optlev_to_multilev(env%optlev,multilevel)
      call crest_multilevel_oloop(env,ensnam,multilevel,i)
      call tim%stop(3)
      if (env%iostatus_meta .ne. 0) return

!>--- save the CRE under a backup name
      call checkname_xyz(crefile,atmp,str)
      call checkname_xyz('.cre',str,btmp)
      call rename(atmp,btmp)
!>--- save cregen output
      call checkname_tmp('cregen',atmp,btmp)
      call rename('cregen.out.tmp',btmp)

!=========================================================!
!>--- cleanup and state update after first iteration (before checkpoint)
      firstiter = (i .eq. 1 .and. start)
      if (firstiter) then
        start = .false.
!>-- obtain a first lowest energy as reference
        env%eprivious = env%elowest
!>-- remove the two extreme-value MTDs
        if (.not.env%readbias.and.env%runver .ne. 33.and. &
        &   env%runver .ne. 787878) then
          env%nmetadyn = env%nmetadyn-2
        end if
!>-- the cleanup
        call clean_V2i
      end if
!>--- checkpoint after this MTD iteration (nmetadyn already updated above)
      call write_restart_log(crest_imtd,'mtd_loop',env%nreset,i, &
        &  env%nmetadyn,env%elowest,env%eprivious,trim(str))
!>-- always do two cycles of MTDs
      if (firstiter) cycle mtdloop
!=========================================================!
!>--- Check for lowest energy
      call elowcheck(lower,env)
      if (.not.lower) then
        exit mtdloop
      end if
!>--- a lower conformer was found: seed the next MTD round from it
      call env%ref%to(mol)
    end do mtdloop
    end if !> end skip_mtdloop guard
    skip_mtdloop = .false.
    do_restart = .false.
!=========================================================!
!>--- collect all ensembles from mtdloop and merge
    if (skip_collect) then
!>--- post_collect restart: collectcre already ran, reuse last file
      inquire(file=trim(rdat%last_file),exist=fex)
      if (.not.fex) then
        write (stdout,'(/,a)') '**ERROR** restart ensemble not found: ' &
          &  //trim(rdat%last_file)
        write (stdout,'(a,/)') ' Delete crest.restart and rerun from scratch.'
        call creststop(status_safety)
      end if
      atmp = trim(rdat%last_file)
      write (stdout,'(1x,a,a)') 'Restarting from ensemble: ',trim(atmp)
      skip_collect = .false.
    else
      write (stdout,*)
      write (stdout,'(''========================================'')')
      write (stdout,'(''           MTD Simulations done         '')')
      write (stdout,'(''========================================'')')
      write (stdout,'(1x,''Collecting ensmbles.'')')
!>-- collecting all ensembles saved as ".cre_*.xyz"
      call collectcre(env)
      call newcregen(env,0)
      call checkname_xyz(crefile,atmp,btmp)
!>--- checkpoint after collection and CREGEN
      call write_restart_log(crest_imtd,'post_collect',env%nreset,0, &
        &  env%nmetadyn,env%elowest,env%eprivious,trim(atmp))
    end if
!>--- remaining number of structures
    call remaining_in(atmp,env%ewin,nallout)

!=========================================================!
!>--- (optional) Perform additional MDs on the lowest conformers
    if (env%rotamermds) then
      call tim%start(4,'Molecular dynamics (MD)')
      call crest_rotamermds(env,conformerfile)
      call tim%stop(4)
      if (env%iostatus_meta .ne. 0) return

!>--- Reoptimization of trajectories
      call checkname_xyz(crefile,atmp,btmp)
      write (stdout,'('' Appending file '',a,'' with new structures'')') trim(atmp)
      ensnam = 'crest_dynamics.trj.xyz'
      call appendto(ensnam,trim(atmp))
      call tim%start(3,'Geometry optimization')
      call crest_multilevel_wrap(env,trim(atmp),-1)
      call tim%stop(3)
      if (env%iostatus_meta .ne. 0) return

      call elowcheck(lower,env)
      if (lower) then
        call checkname_xyz(crefile,atmp,str)
        call checkname_xyz('.cre',str,btmp)
        call rename(atmp,btmp)
!>--- restart sampling from the new lowest structure
        call env%ref%to(mol)
        cycle MAINLOOP
      end if
    end if

!=========================================================!
!>--- (optional) Perform GC step
    if (env%performCross) then
      call tim%start(5,'Genetic crossing (GC)')
      call crest_newcross3(env)
      call tim%stop(5)
      if (env%iostatus_meta .ne. 0) return

      call confg_chk3(env)
      call elowcheck(lower,env)
      if (lower) then
        call checkname_xyz(crefile,atmp,str)
        call checkname_xyz('.cre',str,btmp)
        call rename(atmp,btmp)
        if (env%iterativeV2) then
!>--- restart sampling from the new lowest structure
          call env%ref%to(mol)
          cycle MAINLOOP
        end if
      end if
    end if

!==========================================================!
!>--- exit mainloop
    exit MAINLOOP
  end do MAINLOOP

!==========================================================!
!>--- final ensemble optimization
  write (stdout,'(/)')
  write (stdout,'(3x,''================================================'')')
  write (stdout,'(3x,''|           Final Geometry Optimization        |'')')
  write (stdout,'(3x,''================================================'')')
  call tim%start(3,'Geometry optimization')
  call checkname_xyz(crefile,atmp,str)
  call crest_multilevel_wrap(env,trim(atmp),0)
  call tim%stop(3)
  if (env%iostatus_meta .ne. 0) return

!==========================================================!
!>--- run is complete: drop the restart checkpoint
  call delete_restart_log()

!==========================================================!
!>--- final ensemble sorting
!  call newcregen(env,0)
!> this is actually done within the last crest_multilevel_
!> call, so I comment it out here

!==========================================================!
!>--- print CREGEN results and clean up Directory a bit
  write (stdout,'(/)')
  call smallhead('Final Ensemble Information')
  call V2terminating()

!==========================================================!
  return
end subroutine crest_search_imtdgc

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!

subroutine restart_restore_reference(env,fname)
!****************************************************************
!* Reload the lowest (first) structure of the last CREGEN-sorted
!* ensemble into env%ref%xyz on restart.
!*
!* During an uninterrupted run CREGEN sets env%ref%xyz to the
!* lowest conformer (saveelow branch). The crest.restart
!* checkpoint stores no geometry, so after a restart env%ref%xyz
!* still holds the input structure. Re-reading it here from the
!* recorded ensemble restores the reference point a continuous
!* run would have had. Units match: both env%ref%xyz and the
!* coord ensemble are in Bohr. Silently returns if the file is
!* missing/empty; the working 'mol' geometry is not modified.
!*
!* Arguments:
!*   env   - system data (env%ref%xyz is updated in place)
!*   fname - last ensemble file recorded in the checkpoint
!****************************************************************
  use crest_parameters,only:wp,stdout
  use crest_data
  use strucrd
  implicit none
  type(systemdata),intent(inout) :: env
  character(len=*),intent(in)    :: fname
  type(coord),allocatable :: structures(:)
  integer :: nall
  logical :: ex

  if (len_trim(fname) == 0) return
  inquire (file=trim(fname),exist=ex)
  if (.not.ex) return

  call rdensemble(trim(fname),nall,structures)
  if (nall < 1) return
  if (.not.allocated(structures(1)%xyz)) return

!>--- only overwrite when the atom count matches the reference
  if (allocated(env%ref%xyz)) then
    if (size(env%ref%xyz,2) == structures(1)%nat) then
      env%ref%xyz = structures(1)%xyz
      write (stdout,'(1x,a,a)') &
        & 'Restart: reference geometry restored from ',trim(fname)
    end if
  end if
end subroutine restart_restore_reference

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!

subroutine crest_multilevel_wrap(env,ensnam,level)
!*************************************************
!* wrapper for the multilevel_oloop to select
!* only a single optimization level
!*************************************************
  use crest_parameters,only:wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  implicit none
  type(systemdata) :: env
  character(len=*),intent(in) :: ensnam
  integer,intent(in) :: level
  logical :: multilevel(6)
  integer :: k
  multilevel = .false.
  select case (level)
  case (1:) !> explicit selection (level is a positie integer)
    k = min(level,6)
    multilevel(k) = .true.
  case default
    !>-- map global variable to multilevel selection (level is 0 or negative)
    k = optlevmap_alt(env%optlev)+level
    k = max(1,k)
    multilevel(k) = .true.
  end select
  call crest_multilevel_oloop(env,ensnam,multilevel,0)
end subroutine crest_multilevel_wrap

!========================================================================================!
subroutine crest_multilevel_oloop(env,ensnam,multilevel_in,mtd_iter_in)
!*******************************************************
!* multilevel optimization loop.
!* construct consecutive optimizations starting with
!* crude thresholds to very tight ones.
!* mtd_iter_in: when > 0, writes a mtd_trj restart
!* checkpoint after each CREGEN step (pass 0 to skip).
!*******************************************************
  use crest_parameters,only:wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
  use utilities
  use parallel_interface
  use crest_restartlog
  implicit none
  type(systemdata) :: env
  character(len=*),intent(in) :: ensnam
  logical,intent(in) :: multilevel_in(6)
  integer,intent(in) :: mtd_iter_in
  integer :: nat,nall
  logical :: dump,pr
  character(len=128) :: inpnam,outnam
  integer :: i,l,k,T,Tn,j
  real(wp) :: ewinbackup,rthrbackup
  real(wp) :: hlowbackup
  integer :: microbackup
  integer :: optlevelbackup
  logical :: multilevel(6)
  type(coord),allocatable :: structures(:)

  interface
    subroutine crest_refine(env,input,output)
      use crest_data
      implicit none
      type(systemdata),intent(inout) :: env
      character(len=*),intent(in) :: input
      character(len=*),intent(in),optional :: output
    end subroutine crest_refine
  end interface

!>--- save backup thresholds
  ewinbackup = env%ewin
  rthrbackup = env%rthr
  optlevelbackup = env%calc%optlev
  hlowbackup = env%calc%hlow_opt
  microbackup = env%calc%micro_opt

!>--- set multilevels, or enforce just one
  multilevel(:) = .false.
  if (env%multilevelopt) then
    multilevel(:) = multilevel_in(:)
  else
    k = optlevmap_alt(env%optlev)
    multilevel(k) = .true.
  end if

  pr = .false.
  l = count(multilevel)
  if (l > 1) then
    pr = .true.
    write (stdout,*)
    write (stdout,'(1x,a)') '======================================'
    write (stdout,'(1x,a)') '|  Multilevel Ensemble Optimization  |'
    write (stdout,'(1x,a)') '======================================'
  end if

!>--- read ensemble
  call rdensembleparam(ensnam,nat,nall)
  if (nall .lt. 1) then
    write (stdout,*) '**ERROR** empty ensemble file ',trim(ensnam)
    env%iostatus_meta = status_failed
    return
  end if
  call rdensemble(ensnam,nall,structures)

  write (stdout,'(1x,a,i0,a,a,a)') 'Optimizing all ',nall, &
  & ' structures from file "',trim(ensnam),'" ...'

!>--- sequential optimizations of ensembles
  dump = .true. !> optimized structures will be written to crest_ensemble.xyz
  do i = 1,6
    if (multilevel(i)) then
      !>--- set threads
      call new_ompautoset(env,'auto',nall,T,Tn)
      !>--- set optimization parameters
      call set_multilevel_options(env,i,.true.)
      !>--- run parallel optimizations
      call crest_oloop(env,nall,structures,dump)
      deallocate (structures)
      !>--- rename ensemble and sort
      call checkname_xyz(crefile,inpnam,outnam)
      call rename(ensemblefile,trim(inpnam))
      !>--- check for empty ensemble content
      call rdensembleparam(trim(inpnam),nat,nall)
      if (nall .lt. 1) then
        write (stdout,*) '**ERROR** empty ensemble file',trim(inpnam)
        env%iostatus_meta = status_failed
        return
      end if

      write (stdout,*)
      !==========================================================!
      !>-- dedicated ensemble refinement step (overwrites inpnam)
      call crest_refine(env,trim(inpnam))
      !==========================================================!

      !>--- CREGEN sorting
      call sort_and_check(env,trim(inpnam))
      call checkname_xyz(crefile,inpnam,outnam)
! ── restart checkpoint: intermediate ensemble after this opt. level ──
      if (mtd_iter_in > 0) then
        call write_restart_log(env%crestver,'mtd_trj',env%nreset, &
          &  mtd_iter_in,env%nmetadyn,env%elowest,env%eprivious,trim(inpnam))
      end if
      !>--- check for empty ensemble content (again)
      call rdensembleparam(trim(inpnam),nat,nall)
      if (nall .lt. 1) then
        write (stdout,*) '**ERROR** empty ensemble file',trim(inpnam)
        env%iostatus_meta = status_failed
        return
      end if
      !>--- read new ensemble for next iteration
      call rdensemble(trim(inpnam),nall,structures)

      !>--- restore default sorting thresholds
      env%ewin = ewinbackup
      env%rthr = rthrbackup
      env%calc%optlev = optlevelbackup
      env%calc%hlow_opt = hlowbackup
      env%calc%micro_opt = microbackup
    end if
  end do

  if (allocated(structures)) deallocate (structures)
  return
contains
  subroutine set_multilevel_options(env,i,pr)
    implicit none
    type(systemdata) :: env
    integer,intent(in) :: i
    logical,intent(in) :: pr

    env%calc%hlow_opt = env%hlowopt
    env%calc%micro_opt = nint(env%microopt)

    select case (i)
    case (1)
      if (pr) call smallhead('crude pre-optimization')
      env%calc%optlev = -3
      !> larger thresholds
      env%rthr = env%rthr*2.0d0
      env%ewin = aint(env%ewin*2.0d0)
    case (2)
      if (pr) call smallhead('optimization with very loose thresholds')
      env%calc%optlev = -2
      env%rthr = env%rthr*1.5d0
      env%ewin = aint(env%ewin*2.0d0)
    case (3)
      if (pr) call smallhead('optimization with loose thresholds')
      env%calc%optlev = -1
      env%ewin = aint(env%ewin*(10.0d0/6.0d0))
    case (4)
      if (pr) call smallhead('optimization with regular thresholds')
      env%calc%optlev = 0
    case (5)
      if (pr) call smallhead('optimization with tight thresholds')
      env%calc%optlev = 1
    case (6)
      if (pr) call smallhead('optimization with very tight thresholds')
      env%calc%optlev = 2
    case default
      if (pr) call smallhead('optimization with default thresholds')
      env%ewin = 6.0_wp
      env%rthr = 0.125_wp
      env%calc%optlev = 0
    end select

    call print_opt_data(env%calc,stdout,natoms=env%ref%nat)

  end subroutine set_multilevel_options
end subroutine crest_multilevel_oloop

!========================================================================================!
subroutine crest_rotamermds(env,ensnam)
!***********************************************************
!* set up and perform several MDs at different temperatures
!* on the lowest few conformers
!***********************************************************
  use crest_parameters,only:wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  implicit none
  type(systemdata),intent(inout) :: env
  character(len=*),intent(in) :: ensnam

  integer :: nsim
  type(mddata) :: mddat
  type(mddata),allocatable :: mddats(:)
  type(coord) :: mol
  type(coord),allocatable :: mols(:)
  integer :: nat,nall
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  integer :: nstrucs,i,j,k,io
  real(wp) :: temp,newtemp
  character(len=80) :: atmp

!>--- coord setup
  call env%ref%to(mol)
  call rdensembleparam(ensnam,nat,nall)
  if (nall .lt. 1) then
    write (stdout,*) '**ERROR** empty ensemble file',trim(ensnam)
    env%iostatus_meta = status_failed
    return
  end if

!>--- determine how many MDs need to be run and setup
  call adjustnormmd(env)
  nstrucs = min(nall,env%nrotammds)
  nsim = nstrucs*env%temps
  call crest_search_multimd_init(env,mol,mddat,nsim)
  allocate (mddats(nsim),source=mddat)
  call crest_search_multimd_init2(env,mddats,nsim)
!>--- adjust T's and runtimes
  k = 0
  do i = 1,env%temps
    !> each T block 100K higher
    temp = env%nmdtemp+(i-1)*100.0_wp
    do j = 1,nstrucs
      k = k+1
      mddats(k)%tsoll = temp
      !> reduce runtime by 50% compared to MTDs
      mddats(k)%length_ps = mddats(k)%length_ps*0.5_wp
      call mdautoset(mddats(k),io)
    end do
  end do

!>--- read ensemble and prepare mols
  allocate (xyz(3,nat,nall),at(nat),eread(nall))
  call rdensemble(ensnam,nat,nall,at,xyz,eread)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: mols must be in Bohrs
  xyz = xyz/bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
  allocate (mols(nsim),source=mol)
  k = 0
  do i = 1,env%temps
    do j = 1,nstrucs
      k = k+1
      mols(k)%at = at
      mols(k)%xyz(:,:) = xyz(:,:,j)
    end do
  end do
  deallocate (eread,at,xyz)

!>--- print what we are doing
  write (stdout,*)
  write (atmp,'(''Additional regular MDs on lowest '',i0,'' conformer(s)'')') nstrucs
  call smallheadline(trim(atmp))

!>--- and finally, run the MDs
  call crest_search_multimd2(env,mols,mddats,nsim)

  if (allocated(mols)) deallocate (mols)
  if (allocated(mddats)) deallocate (mddats)
  return
end subroutine crest_rotamermds

!========================================================================================!
subroutine crest_newcross3(env)
!*********************************************************
!* wrapper for the conformational crossing
!* takes the latest crest_rotamers_*, crosses structures
!* and writes the optimized ones back to the file
!*********************************************************
  use crest_parameters
  use crest_data
  use iomod
  use utilities
  implicit none
  type(systemdata) :: env
  real(wp) :: ewinbackup
  integer  :: i,imax,tmpconf,nremain
  character(len=128) :: inpnam,outnam,refnam
  character(len=512) :: thispath,tmppath
  logical :: multilevel(6)
  real(wp),allocatable :: backupthr(:)

  multilevel = .false.
  call getcwd(thispath)

  do i = 1,1  !>-- technically it would be possible to repeat the crossing
!>-- determine max number of new structures
    imax = min(nint(env%mdtime*50.0d0),5000)
    if (env%setgcmax) then
      imax = nint(env%gcmax)
    else if (imax < 0) then
      imax = 5000
    end if
    if (env%quick) then
      imax = nint(float(imax)*0.5d0)
    end if

!>-- call the crossing routine
    call checkname_xyz(crefile,refnam,tmppath)
    call touch(trim(tmppath))
    call crest_crossing(env,imax,trim(refnam),env%gcmaxparent)
    if (imax .lt. 1) then
      call remove(trim(tmppath))
      return
      exit
    end if

!>-- optimize ensemble
    if (env%gcmultiopt) then !>-- optionally split into two steps
      call optlev_to_multilev(env%optlev,multilevel)
    else
      multilevel(4) = .true.
    end if
    call crest_multilevel_oloop(env,'confcross.xyz',multilevel,0)
    if (env%iostatus_meta .ne. 0) return

!>-- append optimized crossed structures and original to a single file
    call checkname_xyz(crefile,inpnam,outnam)
    write (stdout,'(a,a)') 'appending new structures to ',trim(refnam)
    call appendto(trim(inpnam),trim(refnam))
    do while (trim(inpnam) .ne. trim(refnam))
      call remove(trim(inpnam))
      call checkname_xyz(crefile,inpnam,outnam)
    end do
  end do
end subroutine crest_newcross3

!=============================================================================!
!#############################################################################!
!=============================================================================!

subroutine crest_sampling_skip(env,doreturn)
  use crest_data
  implicit none
  type(systemdata),intent(inout) :: env
  logical,intent(out) :: doreturn

  doreturn = .false.

  !> Diatomic molecules need no sampling
  if (env%ref%nat .le. 2) then
    call catchdiatomic(env)
    doreturn = .true.
    return
  end if

  !> "alkylize" special runtype
  call crest_proxy_nalkane(env,doreturn)

end subroutine crest_sampling_skip
