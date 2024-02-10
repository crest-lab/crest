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
!* This is a TODO
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

!===========================================================!
!>--- setup
  call env%ref%to(mol)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)

!>--- sets the MD length according to a flexibility measure
  call md_length_setup(env)
!>--- create the MD calculator saved to env
  call env_to_mddat(env)

  if (env%performMTD) then
!>--- (optional) calculate a short 1ps test MTD to check settings
   call tim%start(1,'Trial metadynamics (MTD)')
   call trialmd(env)    
   call tim%stop(1)
  end if

!===========================================================!
!>--- Start mainloop
  env%nreset = 0
  start = .true.
  MAINLOOP: do
    call printiter
    if (.not.start) then
!>--- clean Dir for new iterations, but leave iteration backup files
      call clean_V2i
      env%nreset = env%nreset+1
    else
!>--- at the beginning, wipe directory clean
      call V2cleanup(.false.)
    end if
!===========================================================!
!>--- Meta-dynamics loop
    mtdloop: do i = 1,env%Maxrestart

      write (stdout,*)
      write (stdout,'(1x,a)') '------------------------------'
      write (stdout,'(1x,a,i0)') 'Meta-Dynamics Iteration ',i
      write (stdout,'(1x,a)') '------------------------------'

      nsim = -1 !>--- enambles automatic MTD setup in init routines
      call crest_search_multimd_init(env,mol,mddat,nsim)
      allocate (mddats(nsim),source=mddat)
      call crest_search_multimd_init2(env,mddats,nsim)

      call tim%start(2,'Metadynamics (MTD)')
      call crest_search_multimd(env,mol,mddats,nsim)
      call tim%stop(2)
!>--- a file called crest_dynamics.trj should have been written
      ensnam = 'crest_dynamics.trj'
!>--- deallocate for next iteration
      if (allocated(mddats)) deallocate (mddats)

!==========================================================!
!>--- Reoptimization of trajectories
      call tim%start(3,'Geometry optimization')
      call optlev_to_multilev(env%optlev,multilevel)
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
      if (i .eq. 1.and.start) then
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
!>-- and always do two cycles of MTDs
        cycle mtdloop
      end if
!=========================================================!
!>--- Check for lowest energy
      call elowcheck(lower,env)
      if (.not.lower) then
        exit mtdloop
      end if
    end do mtdloop
!=========================================================!
!>--- collect all ensembles from mtdloop and merge
    write (stdout,*)
    write (stdout,'(''========================================'')')
    write (stdout,'(''           MTD Simulations done         '')')
    write (stdout,'(''========================================'')')
    write (stdout,'(1x,''Collecting ensmbles.'')')
!>-- collecting all ensembles saved as ".cre_*.xyz"
    call collectcre(env)
    call newcregen(env,0)
    call checkname_xyz(crefile,atmp,btmp)
!>--- remaining number of structures
    call remaining_in(atmp,env%ewin,nallout)

!=========================================================!
!>---- Entropy mode iterative statically biased MDs
    if (env%entropymd) then
!>--- determine how many MDs need to be run and setup
!>--- and other entropy mode parameters
      call adjustnormmd(env)
      call mtdatoms(env)
      call emtdcopy(env,0,stopiter,fail)
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
            call crest_multilevel_oloop(env,trim(atmp),multilevel)
            call tim%stop(3)

!>--- if in the entropy mode a lower structure was found -> cycle (required for extrapolation)
            call elowcheck(lower,env)
            if (lower.and.env%entropic) then
              env%emtd%nbias = bref  !> IMPORTANT, reset for restart
              cycle MAINLOOP
            end if

!>--- otherwise, handle files andfile handling
            eit2 = eit
            call emtdcopy(env,eit2,stopiter,fail)
            env%emtd%iterlast = eit2
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
    write (stdout,*) 'empty ensemble file',trim(ensnam)
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
!>--- output will be collected in crest_dynamics.trj
!>--- but the entropy routines look for the crest_rotamers_ files
  call checkname_xyz(crefile,atmp,btmp)
  call rename('crest_dynamics.trj',atmp)
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
  if (.not.ex) error stop 'Could not initialize static metadynamics'
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

