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

subroutine crest_search_entropy(env,tim)
!*******************************************************************
!* This is the re-implementation of CREST's sMTD-iMTD workflow
!* from https://doi.org/10.1039/d1sc00621e
!* with calculation of conformational entropy
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
!> Entropy algo variables
  logical :: stopiter,fail
  integer :: bref,dum,eit,eit2
!===========================================================!
  type(restart_data) :: rdat
  logical :: do_restart,skip_mtdloop,skip_collect,skip_emtdcopy0,firstiter,fex
!===========================================================!
!>--- printout header
  write (stdout,*)
  write (stdout,'(10x,"┍",49("━"),"┑")')
  write (stdout,'(10x,"│",14x,a,13x,"│")') "CREST ENTROPY SAMPLING"
  write (stdout,'(10x,"┕",49("━"),"┙")')
  write (stdout,*)
  write (stdout,'(1x,a)') 'please cite:'
  write (stdout,'(1x,a)') '• P.Pracht, S.Grimme, Chem. Sci., 2021, 12, 6551-6568.'
  write (stdout,'(1x,a)') '• J.Gorges, S.Grimme, A.Hansen, P.Pracht, PCCP, 2022,24, 12249-12259.'
  write (stdout,*)

! ── restart detection ─────────────────────────────────────────────
  do_restart = .false.
  skip_mtdloop = .false.
  skip_collect = .false.
  skip_emtdcopy0 = .false.
  if (env%allowrestart .and. restart_file_exists()) then
    call read_restart_log(rdat)
    if (rdat%runtype == env%crestver .and. rdat%stage /= 'done') then
      do_restart = .true.
      call print_restart_info(rdat)
      !> skip entire mtdloop and collectcre when past the MTD loop
      skip_mtdloop = (rdat%stage == 'post_collect' .or. &
        &              rdat%stage == 'entropy_smtd')
      skip_collect = (rdat%stage == 'post_collect' .or. &
        &              rdat%stage == 'entropy_smtd')
      !> additionally skip emtdcopy(iter=0) when that call already ran
      skip_emtdcopy0 = (rdat%stage == 'entropy_smtd')
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
   if(env%iostatus_meta .ne. 0) return
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
        call write_restart_log(env%crestver,'mtd_trj',env%nreset,i, &
          &  env%nmetadyn,env%elowest,env%eprivious,ensnam)
      end if

!==========================================================!
!>--- Reoptimization of trajectories
      call tim%start(3,'Geometry optimization')
      call optlev_to_multilev(env%optlev,multilevel)
      call crest_multilevel_oloop(env,ensnam,multilevel,i)
      call tim%stop(3)
      if(env%iostatus_meta .ne. 0 ) return

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
      call write_restart_log(env%crestver,'mtd_loop',env%nreset,i, &
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
!>--- (env%ref is kept at the current lowest by CREGEN)
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
      call write_restart_log(env%crestver,'post_collect',env%nreset,0, &
        &  env%nmetadyn,env%elowest,env%eprivious,trim(atmp))
    end if
!>--- remaining number of structures
    call remaining_in(atmp,env%ewin,nallout)

!=========================================================!
!>---- Entropy mode iterative statically biased MDs
    if (env%entropymd) then
!>--- determine how many MDs need to be run and setup
!>--- and other entropy mode parameters
      call adjustnormmd(env)
      call mtdatoms(env)
      if (.not.skip_emtdcopy0) then
        call emtdcopy(env,0,stopiter,fail)
! ── checkpoint: entropy rotamer file written, sMTD iterations about to start ──
        if (env%crestver == crest_imtd2) then
          write (btmp,'(a,i0,a)') 'crest_smtd_',0,'.xyz'
        else
          write (btmp,'(a,i0,a)') 'crest_entropy_rotamer_',0,'.xyz'
        end if
        call write_restart_log(env%crestver,'entropy_smtd',env%nreset,0, &
          &  env%nmetadyn,env%elowest,env%eprivious,trim(btmp))
      end if
      bref = env%emtd%nbias

!>--- sMTD iterations, done until max iterations or convergence
      ENTROPYITER: do eit = 1,env%emtd%iter
        !> Modify bias
        dum = nint(float(env%emtd%nbias)*env%emtd%nbiasgrow)
        env%emtd%nbias = max(env%emtd%nbias+1,dum)
        fail = .false.

!>--- Loop handling fallbacks
        EFALLBACK: do k = 1,env%emtd%maxfallback
          call printiter2(eit)
          call tim%start(6,'Static metadynamics (sMTD)')
          !>-- start from the current crest_conformers.xyz
          call crest_smtd_mds(env,conformerfile)
          call tim%stop(6)
          if(env%iostatus_meta .ne. 0) return
          call emtdcheckempty(env,fail,env%emtd%nbias)

          if (fail) then
            if (k == env%emtd%maxfallback) then
              stopiter = .true.
            else
              cycle EFALLBACK
            end if
          else

!!>--- Reoptimization of trajectories
            call checkname_xyz(crefile,atmp,btmp)
            call tim%start(3,'Geometry optimization')
            multilevel = (/.true.,.false.,.false.,.false.,.false.,.true./)
            call crest_multilevel_oloop(env,trim(atmp),multilevel,0)
            call tim%stop(3)
            if(env%iostatus_meta .ne. 0 ) return

!>--- if in the entropy mode a lower structure was found -> cycle (required for extrapolation)
            call elowcheck(lower,env)
            if (lower.and.env%entropic) then
              env%emtd%nbias = bref  !> IMPORTANT, reset for restart
!>--- restart sampling from the new lowest structure
              call env%ref%to(mol)
              cycle MAINLOOP
            end if

!>--- otherwise, handle files andfile handling
            eit2 = eit
            call emtdcopy(env,eit2,stopiter,fail)
            env%emtd%iterlast = eit2
! ── checkpoint: update last_file to current entropy rotamer file ──────
            if (env%crestver == crest_imtd2) then
              write (btmp,'(a,i0,a)') 'crest_smtd_',eit2,'.xyz'
            else
              write (btmp,'(a,i0,a)') 'crest_entropy_rotamer_',eit2,'.xyz'
            end if
            call write_restart_log(env%crestver,'entropy_smtd',env%nreset,0, &
              &  env%nmetadyn,env%elowest,env%eprivious,trim(btmp))
          end if

          if (.not.lower.and.fail.and..not.stopiter) then
            cycle EFALLBACK
          end if

          exit EFALLBACK  !> fallback loop is exited on first opportuinity
        end do EFALLBACK

        if (stopiter) then
          exit ENTROPYITER
        end if

      end do ENTROPYITER
    end if

!==========================================================!
!>--- exit mainloop
    exit MAINLOOP
  end do MAINLOOP

!==========================================================!
!>--- run is complete: drop the restart checkpoint
  call delete_restart_log()

!==========================================================!
!>--- print CREGEN results and clean up Directory a bit
  write (stdout,'(/)')
  call smallhead('Final Ensemble Information')
  call V2terminating()

!==========================================================!
  return
end subroutine crest_search_entropy

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!

subroutine crest_smtd_mds(env,ensnam)
!***********************************************************
!* set up and perform several sMTD's on a number of
!* conformers obtained from clustering.
!* The input ensemble (read from ensnam) is typically the
!* conformer file.
!***********************************************************
  use crest_parameters,only:wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use iomod
  use utilities
  use dynamics_module
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
  character(len=128) :: atmp,btmp
!============================================================!
  integer :: nclustbackup
  integer :: TOTAL
!============================================================!
!>--- coord setup
  call env%ref%to(mol)
  call rdensembleparam(ensnam,nat,nall)
  if (nall .lt. 1) then
    write (stdout,*) '**ERROR** empty ensemble file',trim(ensnam)
    env%iostatus_meta = status_failed
    return
  end if

!============================================================!
!>--- PCA/k-Means Cluster setup
!============================================================!
  nclustbackup = env%maxcluster

  !>--- first, generate the structures will be used as bias
  env%nclust = env%emtd%nbias  !> this determines how many clusters will be build
  call create_anmr_dummy(nat)
  call smallhead('determining bias structures via PCA/k-Means')
  call CCEGEN(env,.false.,ensnam)  !> this routine does PCA/k-Means
  call rdensembleparam(clusterfile,nat,TOTAL)
  if (TOTAL < 1) then
    call copy('crest_best.xyz',clusterfile)
    TOTAL = 1
  end if
  write (*,'(1x,i0,a)') TOTAL,' structures were selected'
  write (*,'(1x,a,/)') 'done.'
  env%mtdstaticfile = "crest_bias.xyz"
  env%nstatic = TOTAL
  call rename(clusterfile,env%mtdstaticfile)

  !>--- then, get the input structures
  env%nclust = env%emtd%nMDs  !> this determines how many clusters will be build
  call smallhead('determining MTD seed structures via PCA/k-Means')
  call CCEGEN(env,.false.,ensnam)  !> this routine does PCA/k-Means
  call rdensembleparam(clusterfile,nat,TOTAL)
  write (stdout,'(1x,i0,a)') TOTAL,' structures were selected'
  write (stdout,'(1x,a,/)') 'done.'

  !>--- and cleanup
  call remove('anmr_nucinfo')
  env%nclust = nclustbackup
!============================================================!
!============================================================!

!>--- Generate the required number of static MD calculators
  nsim = min(TOTAL,env%emtd%nMDs) !> from the generated cluster, but limited to env%emtd%nMDs
  call crest_search_multimd_init(env,mol,mddat,nsim) !> general mddat setup
  allocate (mddats(nsim),source=mddat)
!>--- adjust T's and runtimes, and load the bias
  call crest_init_multimd_smtd(env,mddats,nsim,env%mtdstaticfile)

!>--- read cluster ensemble and prepare mols to start MTDs from
  call rdensemble(clusterfile,nall,mols)

!>--- print what we are doing
  write (atmp,'(''Static MTDs (umbrella sampling) on '',i0,'' selected conformer(s)'')') nsim
  call smallheadline(trim(atmp))
  write (stdout,'("> Using ",i0," constant RMSD bias potentials per MTD")') env%nstatic

!===================================================================!
!>--- and finally, run the sMTDs on the different starting structures
  call crest_search_multimd2(env,mols,mddats,nsim)
!>--- output will be collected in crest_dynamics.trj.xyz
!>--- but the entropy routines look for the crest_rotamers_ files
  call checkname_xyz(crefile,atmp,btmp)
  call rename('crest_dynamics.trj.xyz',atmp)
!===================================================================!
!>--- by default, clean up the directory
  if (.not.env%keepModef) call cleanMTD

!>--- deallocate molecule and MD containers
  if (allocated(mols)) deallocate (mols)
  if (allocated(mddats)) deallocate (mddats)
  return
end subroutine crest_smtd_mds

!=========================================================================================!
subroutine crest_init_multimd_smtd(env,mddats,nsim,biasfile)
!**************************************************************
!* Append a list of MD calculators (mddats),
!* change them to static metadynamics and
!* adujst otherwise needed parameter such as the temperature.
!* Bias structures will be read from biasfile
!*
!* The routines adjustnormmd() and mtdatoms() must have been
!* called before calling this routine so all the required data
!* is initialized!
!**************************************************************
  use crest_parameters,only:wp,stdout,bohr,sep
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use iomod,only:makedir,directory_exist,remove
  use utilities 
!$ use omp_lib
  implicit none
  type(systemdata),intent(inout) :: env
  type(mddata),intent(inout) :: mddats(nsim)
  integer,intent(in) :: nsim
  character(len=*),intent(in) :: biasfile
  integer :: i,io
  integer :: nat,nall
  integer,allocatable :: at(:)
  real(wp),allocatable :: xyz(:,:,:)
  real(wp),allocatable :: eread(:)
  logical :: ex
  integer :: idum1
  real(wp) :: dum1
  type(mtdpot),allocatable :: mtds(:)
  type(mtdpot) :: mtdtmp
  character(len=80) :: atmp
  character(len=*),parameter :: mdir = 'STATICMTD'

!>--- parallel MD setup, prepare files
  ex = directory_exist(mdir)
  if (ex) then
    call rmrf(mdir)
  end if
  io = makedir(mdir)
  do i = 1,nsim
    mddats(i)%md_index = i
    write (atmp,'(a,i0,a)') 'crest_',i,'.trj'
    mddats(i)%trajectoryfile = mdir//sep//trim(atmp)
    write (atmp,'(a,i0,a)') 'crest_',i,'.mdrestart'
    mddats(i)%restartfile = mdir//sep//trim(atmp)
!>--- append settings
    mddats(i)%simtype = type_mtd  !> set to MTD runtype (includes the static version)
    mddats(i)%tsoll = env%emtd%temperature !> temperature
    mddats(i)%length_ps = env%mdtime*env%emtd%lenfac  !> simulation length
!>--- complete real-time settings to steps again
    call mdautoset(mddats(i),io)
  end do

!>--- kpush & alpha, and parameters that will be the same for all
  mtdtmp%kpush = env%emtd%katoms*env%emtd%kpush
  mtdtmp%alpha = env%emtd%alpha
  mtdtmp%cvdump_fs = huge(dum1)   !> set to large to avoid new structure dumps
  mtdtmp%cvdumpstep = huge(idum1) !> same
  mtdtmp%mtdtype = cv_rmsd_static  !> set the correct bias type

!>--- load static bias stuctures
  inquire (file=biasfile,exist=ex)
  if (.not.ex)then
    write(stdout,'(a,a)') 'Could not initialize static metadynamics: missing ',trim(biasfile)
    call creststop(status_input)
  endif
  call rdensembleparam(biasfile,nat,nall)
  allocate (xyz(3,nat,nall),at(nat),eread(nall))
  call rdensemble(biasfile,nat,nall,at,xyz,eread)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: bias structures must be in Bohrs
  xyz = xyz/bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

!>--- transfer a copy of mtdtmp to each MD container
  do i = 1,nsim
    if (allocated(mddats(i)%mtd)) deallocate (mddats(i)%mtd)
    if (allocated(mddats(i)%cvtype)) deallocate (mddats(i)%cvtype)
    mddats(i)%npot = 1
    allocate (mddats(i)%mtd(1),source=mtdtmp)
    allocate (mddats(i)%cvtype(1),source=cv_rmsd_static)

    !>--- a ramp parameter depending on timestep (relative to old GFN2-xTB default)
    mddats(i)%mtd(1)%ramp = (mddats(i)%tstep/5.0_wp)*env%emtd%mtdramp

    !>--- the bias structures are transferred here
    allocate (mddats(i)%mtd(1)%cvxyz(3,nat,nall))
    mddats(i)%mtd(1)%cvxyz(:,:,:) = xyz(:,:,:)
    mddats(i)%mtd(1)%ncur = nall    !> will not change
    mddats(i)%mtd(1)%maxsave = nall !> won't change either

    !>--- transfer the atomlist (ther sMTD pot is only acting on the heavy atoms)
    allocate (mddats(i)%mtd(1)%atinclude(nat),source=.false.)
    mddats(i)%mtd(1)%atinclude(:) = env%emtd%atomlist2(:)
  end do

  deallocate (eread,at,xyz)
  return
end subroutine crest_init_multimd_smtd

