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
!  use iso_fortran_env,wp => real64
  use crest_parameters !> Datatypes and constants
  use crest_data !> module for the main data storage (imports systemdata and timer)
  use crest_calculator
  USE,INTRINSIC :: IEEE_EXCEPTIONS
  implicit none
  type(systemdata) :: env  !> MAIN STORAGE OF SYSTEM DATA
  type(calcdata),target :: calc_origin
  type(timer)   :: tim     !> timer object

  integer :: i,j,l,args,io
  character(len=:),allocatable :: arg(:)
  character(len=:),allocatable :: infile
  character(len=512) :: thisdir
  character(len=1024) :: cmd
  real(wp) :: dumfloat,dumfloat2,d3,d4,d5,d6,d7,d8
  logical :: ex,ex1,ex2,iterate

  intrinsic :: iargc,getarg
  LOGICAL :: overflow,division_by_zero,invalid_operation

  call initsignal() !SIGTERM catcher

!=========================================================================================!
!>  Initialize system clock time
  call tim%init(20)

  env%calc => calc_origin
!=========================================================================================!
!> set defaults and pars flags
  args = iargc()
  allocate (arg(args),source=repeat(' ',1024))
  do i = 1,args
    call getarg(i,arg(i))
  end do
  call parseflags(env,arg,args)
  deallocate (arg)

!=========================================================================================!
!> post-parsing sanity checks on the assembled env%calc object
!=========================================================================================!
!>--- g-xTB (05/2026) has no implicit solvation; stop safely if requested
  call gxtb_solvation_guard(env)

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
    call crest_dry_run(env,tim)
  end if

!=========================================================================================!
!> SOME I/O STUFF
!=========================================================================================!
!>--- check for the coord file in the working directory
  if (env%crestver /= crest_solv.and.env%crestver /= crest_sorting) then
    inquire (file='coord',exist=ex)
    if (.not.ex) then
      error stop 'No coord file found. Exit.'
    end if
  end if

!=========================================================================================!
!>        PRE-CONFSEARCH PROPERTY CALCS
!=========================================================================================!
  select case (env%properties)
!>--- zsort routine
  case (p_zsort)
    call zsort
    write (*,*)
    write (*,*) 'The z-matrix of the input coord file has been sorted.'
    write (*,*) 'The sorted file in TM format is called "zcoord"'
    write (*,*)
    write (*,*) 'exit.'
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
!>--- properties for enesemble file
  case (p_propcalc)
    call propcalc(env%ensemblename,env%properties2,env,tim)
    call propquit(tim)
!>--- calculate potential correction for acid/base reaction
  case (p_acidbase)
    call tim%start(4,'acid/base')
    if (env%protb%pka_mode == 0) then
      call acidbase(env,env%protb%pka_acidensemble,env%protb%pka_baseensemble,env%chrg,.true., &
          & .false.,dumfloat,.false.,d3,d4,d5,d6,d7,d8)
    else
      call rewrite_AB_ensemble(env,env%protb%pka_acidensemble,env%protb%pka_baseensemble)
    end if
    call tim%stop(4)
    call propquit(tim)
!>--- calculate potential correction for acid/base reaction
  case (p_ligand)
    call tim%start(4,'')
    call ligandtool(env%protb%infile,env%protb%newligand, &
    &    env%protb%centeratom,env%protb%ligand)
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

!>--- alkylation prep
  if (env%alkylize) then
    call crest_setup_alkylize(env)
  end if

!=========================================================================================!
!>         PRE-OPTIMIZATION OF THE GEOMETRY
!=========================================================================================!
  if (env%crestver /= crest_none) then
    if (env%preopt) then
      call trialOPT(env)
    else if (env%presp) then
      call xtbsp(env)
    end if
  end if
!=========================================================================================!
!>         SET UP QUEUES, IF REQUIRED
!=========================================================================================!
  call crest_queue_setup(env,iterate)
!=========================================================================================!
!>         MAIN WORKFLOW CALLS START HERE (including iterator)
!=========================================================================================!

  ITERATOR: do while (iterate)
    call crest_queue_iter(env,iterate)

!> NOTE: many of these routine calls take a detour through legacy_wrappers.f90 !
    select case (env%crestver)
    case (crest_mfmdgc)           !> MF-MD-GC algo (deprecated)
      call confscript1(env,tim)

    case (crest_imtd,crest_imtd2) !> MTD-GC algo
      call confscript2i(env,tim)

    case (crest_mdopt,crest_mdopt2)
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

    case (crest_ttc) !> tensor-train (light) conformer search
      call crest_ttconf(env,tim)

    case (crest_trialopt) !> test optimization standalone
      call trialOPT(env)

    case (crest_ensemblesp) !> singlepoints along ensemble
      call crest_ensemble_singlepoints(env,tim)

    case (crest_ensemblehess) !> Hessians + thermochemistry along ensemble
      call crest_ensemble_hessians(env,tim)

    case (crest_protonate)
      call protonate(env,tim)

    case (crest_deprotonate)
      call deprotonate(env,tim)

    case (crest_tautomerize)
      call tautomerize(env,tim)

    case (crest_sorting) !> interface to standalone ensemble sorting
      call crest_sort(env,tim)

    case (crest_bh) !> Standard basin-hopping
      call crest_basinhopping(env,tim)

    case (crest_test)
      call crest_playground(env,tim)

    case (crest_none)
      call crest_no_runtype_selected()

    case default
      continue
    end select

    !> additional processing
    call crest_queue_iter_resort(env,iterate)

  end do ITERATOR

  env%calc => calc_origin
  call crest_queue_reconstruct(env,tim)

!=========================================================================================!
!>        ADDITIONAL OUTPUT FORMATTING
!=========================================================================================!
  call crest_ensemble_reformat(env)

!=========================================================================================!
!>        POST-CONFSEARCH PROPERTY CALCS
!=========================================================================================!
  if (env%npq .gt. 0) then
    infile = "crest_rotamers.xyz"
    do i = 1,env%npq
      j = env%pqueue(i)
      select case (j)
      case (p_prop_hess,p_prop_autoir,p_prop_ohess,p_prop_reopt,p_prop_dipole,p_prop_finalhess)
        call propcalc(conformerfile,j,env,tim)

      case (abs(p_CREentropy))
        call tim%start(15,'Conf. entropy evaluation')
        call newentropyextrapol(env)
        call tim%stop(15)

      case (p_prop_multilevel:p_prop_multilevel+9)  !hybrid reoptimization (e.g. gfn2@gff)
        call propcalc(infile,j,env,tim)

      case (abs(p_cluster)) !PCA and clustering
        call ccegen(env,.true.,conformerfile)

      case (abs(p_tautomerize2))
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
!> one final cleanup
  call custom_cleanup(env)

!=========================================================================================!
!> Print a summary of output files written in this run
  call crest_output_summary(env)

!=========================================================================================!
!> Repeat the g-xTB syscall note (non-WITH_GXTB path) if g-xTB was requested
  if (env%gfnver == '--gxtb') call gxtb_syscall_warning()

!=========================================================================================!
!> Evaluate and print timings, then stop the program
  call eval_timer(tim)
  call creststop(env%iostatus_meta)
!> end of main program
end program CREST

!=========================================================================================!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!=========================================================================================!
