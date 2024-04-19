!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2023 Philipp Pracht
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

!=========================================================================================!
!  This is the code of the Conformer-Rotamer Ensemble Sampling Tool (CREST).
!=========================================================================================!
program CREST
  use iso_fortran_env,wp => real64
  !> module for the main data storage
  use crest_data
  use crest_restartlog
  implicit none
  type(systemdata) :: env  !> MAIN STORAGE OF SYSTEM DATA
  type(timer)   :: tim     !> timer object

  integer :: i,j,l,args,io
  character(len=:),allocatable :: arg(:)
  character(len=:),allocatable :: infile
  character(len=512) :: thisdir
  character(len=1024) :: cmd
  real(wp) :: dumfloat,dumfloat2,d3,d4,d5,d6,d7,d8
  logical :: ex,ex1,ex2

  intrinsic :: iargc,getarg

  call initsignal() !SIGTERM catcher

!=========================================================================================!
!>  Initialize system clock time
  call tim%init(20)

!=========================================================================================!
!> set defaults and pars flags
  args = iargc()
  l = len_trim(cmd)
  allocate (arg(args),source=repeat(' ',l))
  do i = 1,args
    call getarg(i,arg(i))
  end do
  call parseflags(env,arg,args)
  deallocate (arg)
  call restart_save_env(env)
!=========================================================================================!
!> scratch dir handling

  if (env%scratch) then
    call getcwd(thisdir)
    call scrdir(env)
  end if

!=========================================================================================!
!>   OMP_NUM_THREAD handling
!=========================================================================================!
  if (.not.env%autothreads) then
    call ompenvset(env%omp)
  else
    if (.not.env%threadssetmanual) then
      call ompgetauto(env%threads,env%omp,env%MAXRUN)
    end if
    call new_ompautoset(env,'max',0,i,j)
  end if

!=========================================================================================!
!>   DRY run stop
!=========================================================================================!
  if (env%dryrun) then
    call crest_dry(env)
  end if

!=========================================================================================!
!> SOME I/O STUFF
!=========================================================================================!
!>--- check for the coord file in the working directory
  if (env%crestver /= crest_solv) then
    inquire (file='coord',exist=ex)
    if (.not.ex) then
      error stop 'No coord file found. Exit.'
    end if
  end if

!=========================================================================================!
!>        PRE-CONFSEARCH PROPERTY CALCS
!=========================================================================================!
  select case (env%properties)
!>--- only CREGEN routine
  case (p_cregen)
    call tim%start(1,'CREGEN')
    write (*,*) 'Using only the cregen sorting routine.'
    env%cgf(6) = .true.   !write confg output to file <tmp>
    if (env%doNMR) then
      env%cgf(3) = .true.
      if (.not.env%fullcre) then
        env%cgf(2) = .false.
      end if
    end if
    if (env%newcregen) then
      block
      use cregen_interface
      call newcregen(env,0)
      end block
    else
      call cregen2(env)
    end if
    if (env%doNMR.and.env%fullcre) then
      call entropic(env,.true.,.false.,.false.,env%ensemblename, &
      &    env%tboltz,dumfloat,dumfloat2)
    end if
    if (env%cluster) then
      call ccegen(env,.true.,ensemblefile)
    end if
    call tim%stop(1)
    call propquit(tim)
!>--- zsort routine
  case(p_zsort) 
    call zsort
    write (*,*)
    write (*,*) 'The z-matrix of the input coord file has been sorted.'
    write (*,*) 'The sorted file in TM format is called "zcoord"'
    write (*,*)
    write (*,*) 'exit.'
    call propquit(tim)

!>--- only ensemble comparison
  case (p_compare)
    call compare_ensembles(env)     
    call propquit(tim)
 !>--- protonation tool
  case (p_protonate)
    call protonate(env,tim)
    call propquit(tim)
!>--- deprotonation
  case (p_deprotonate)
    call deprotonate(env,tim)
    call propquit(tim)
!>--- tautomerization
  case (p_tautomerize)
    call tautomerize(env,tim)
    call propquit(tim)
!>--- extended tautomerization
  case (p_tautomerize2)
    call tautomerize_ext(env%ensemblename,env,tim)
    call propquit(tim)
!>--- stereoisomerization
  case (p_isomerize)
    call stereoisomerize(env,tim)
    call propquit(tim)

!>--- reactor setup
  case (p_reactorset)
    call reactor_setup(env)

    stop
!>--- enhanched ensemble entropy
  case (p_CREentropy)
    call entropic(env,.true.,.true.,.false.,env%ensemblename, &
    &    env%tboltz,dumfloat,dumfloat2)
    call propquit(tim)
!>--- calculate hessians and average thermo. contrib
  case (p_rrhoaverage)
    call tim%start(4,'freq+thermo')
    call calcSrrhoav(env,env%ensemblename)
    call tim%stop(4)
    call propquit(tim)
!>--- to PCA and k-Means clustering for given file
  case (p_cluster)
    call ccegen(env,.true.,env%ensemblename)
    call propquit(tim)
!>--- properties for enesemble file
  case (p_propcalc)
    call propcalc(env%ensemblename,env%properties2,env,tim)
    call propquit(tim)
!>--- calculate potential correction for acid/base reaction
  case (p_acidbase)
    call tim%start(4,'acid/base')
    if (env%ptb%pka_mode == 0) then
      call acidbase(env,env%ptb%pka_acidensemble,env%ptb%pka_baseensemble,env%chrg,.true., &
          & .false.,dumfloat,.false.,d3,d4,d5,d6,d7,d8)
    else
      call rewrite_AB_ensemble(env,env%ptb%pka_acidensemble,env%ptb%pka_baseensemble)
    end if
    call tim%stop(4)
    call propquit(tim)
!>--- calculate potential correction for acid/base reaction
  case (p_ligand)
    call tim%start(4,'')
    call ligandtool(env%ptb%infile,env%ptb%newligand, &
    &    env%ptb%centeratom,env%ptb%ligand)
    call tim%stop(4)
    call propquit(tim)
!>--- wrapper for the thermo routine
  case (p_thermo)
    call tim%start(4,'')
    !call thermo_mini(env)
    call thermo_standalone(env)
    call tim%stop(4)
    call propquit(tim)
!>--- ensemble merging tool
  case (p_gesc1,p_gesc2)
    call tim%start(9,'')
    call biasmerge(env)
    call tim%stop(9)
    if (env%properties == -9224) call propquit(tim)
!>--- do nothing here
  case default
    continue
  end select
!=========================================================================================!
!>         PRE-OPTIMIZATION OF THE GEOMETRY
!=========================================================================================!
  if (env%preopt) then
    call trialOPT(env)
  else if (env%presp) then
    call xtbsp(env)
  end if
!=========================================================================================!
!>         MAIN WORKFLOW CALLS START HERE
!=========================================================================================!
!> many of these routine calls take a detour through legacy_wrappers.f90 !
  select case (env%crestver)
  case (crest_mfmdgc)           !> MF-MD-GC algo (deprecated)
    call confscript1(env,tim)

  case (crest_imtd,crest_imtd2) !> MTD-GC algo
    call confscript2i(env,tim)

  case (crest_mdopt, crest_mdopt2)
    call mdopt(env,tim)        !> MDOPT

  case (crest_screen)
    call screen(env,tim)       !> SCREEN

  case (crest_nano)
    call reactor(env,tim)      !> NANO-REACTOR

  case (crest_compr)
    call compress(env,tim)     !> MTD COMPRESS mode

  case (crest_msreac)
    call msreact_handler(env,tim) !> MSREACT sub-program

  case (crest_pka)
    call pkaquick(env,tim)

  case (crest_solv)             !> microsolvation tools
    call crest_solvtool(env,tim)

  case (crest_sp)
    call crest_singlepoint(env,tim)

  case (crest_optimize)
    call crest_optimization(env,tim)

  case (crest_moldyn)
    call crest_moleculardynamics(env,tim)

  case (crest_s1)
    call crest_search_1(env,tim)

  case (crest_mecp)
    call crest_search_mecp(env,tim)

  case (crest_numhessian)
    call crest_numhess(env,tim)

  case (crest_scanning)
    call crest_scan(env,tim)

  case (crest_rigcon) !> rule-based conformer generation
    call crest_rigidconf(env,tim)

  case (crest_trialopt) !> test optimization standalone
    call trialOPT(env)

  case (crest_ensemblesp) !> singlepoints along ensemble
    call crest_ensemble_singlepoints(env,tim)    

  case (crest_test)
    call crest_playground(env,tim)

  case default
    continue
  end select

  if (env%outputsdf.or.env%sdfformat) then
    if (any((/crest_mfmdgc,crest_imtd,crest_imtd2/) == env%crestver)) then
      call new_wrsdfens(env,conformerfile,conformerfilebase//'.sdf',.false.)
    end if
    if (any((/crest_screen,crest_mdopt/) == env%crestver)) then
      call new_wrsdfens(env,'crest_ensemble.xyz','crest_ensemble.sdf',.false.)
    end if
  end if

!=========================================================================================!
!>        POST-CONFSEARCH PROPERTY CALCS
!=========================================================================================!
  if (env%npq .gt. 0) then
    infile = "crest_rotamers.xyz"
    do i = 1,env%npq
      j = env%pqueue(i)
      select case (j)
      case (1:8,10,20,100,998)
        call propcalc(conformerfile,j,env,tim)
      case (45)
        call tim%start(15,'Conf. entropy evaluation')
        call newentropyextrapol(env)
        call tim%stop(15)
      case (50:59)  !hybrid reoptimization (e.g. gfn2@gff)
        call propcalc(infile,j,env,tim)
        infile = 'crest_reopt.xyz'
      case (70) !PCA and clustering
        call ccegen(env,.true.,conformerfile)
      case (555)
        call tautomerize_ext(infile,env,tim)
      case default
        continue
      end select
    end do
  end if

!=========================================================================================!
!> go back from scratch directory
  if (env%scratch) then
    call chdir(thisdir)
    call scrend(env)
  end if

!=========================================================================================!
!> shout down hosted subprocesses
  block
  use ConfSolv_module
  call cs_shutdown(io)
  end block

!=========================================================================================!
!> Evaluate and print timings
  call eval_timer(tim)
  write (*,*) 'CREST terminated normally.'
!> end of main program
end program CREST

!=========================================================================================!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!=========================================================================================!
