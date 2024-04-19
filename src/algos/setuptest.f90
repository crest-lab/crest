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

!> The following contains routines that are run in the automatic
!> setup of CREST to check for consistent settings.
!> First is the trialMD calculator that can automatically adjust the timestep
!> Second is the initial geometry optimization with topology check

!========================================================================================!
!========================================================================================!
!> Routines for checking metadynamics settings (timestep & shake, mostly)
!========================================================================================!
!========================================================================================!
subroutine trialMD_calculator(env)
!*****************************************************
!* subroutine trialMD_calculator
!* Runs a short metadynamics simulation and, in case
!* the simulation terminates with error, tries to
!* adjust the time step and tries again. Up to a
!* maximum of 6 tries.
!*****************************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  use dynamics_module
  use strucrd
  use iomod
  implicit none

  type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA

  integer :: io,Vdumpfreq,counter
  real(wp) :: prefac,alpha,length
  real(wp) :: rtime,tstep
  integer :: shakemode,maxiter

  character(len=512) :: thispath
  character(len=80)  :: trajectory,atmp

  logical :: shakefail,mdfail,ex,pr
  type(coord) :: mol,molstart
  type(mddata) :: MD,MDSTART
  type(mtdpot) :: MTD
  type(timer) :: profiler

  type(calcdata) :: tmpcalc
  real(wp) :: energy
  real(wp),allocatable :: grd(:,:)
  integer :: T,Tn
  character(len=*),parameter :: dirnam = 'TRIALMD'

!>--- OMP settings (should be set to 1 to simulate max parallelization)
  call new_ompautoset(env,'min',0,T,Tn)

  call getcwd(thispath)

!>--- Make a dir
  call rmrf(dirnam)
  io = makedir(trim(dirnam))

!>--- set up mol and MD calculator
  prefac = 0.003d0*env%rednat !> Vbias k
  alpha = 0.5d0         !> Vbias alpha

  call env%ref%to(molstart)
  MDSTART = env%mddat         !> env%mddat should already be set up at this point
  MDSTART%length_ps = 1.0_wp  !> we want only a short simulaton
  MDSTART%length_steps = 0    !> to zero so it will be calculated automatically
  MDSTART%dumpstep = 20.0_wp  !> fs dump step to trajectory file, so we end up with 50 structures
  MDSTART%sdump = 0           !> also zero to reset
  if (allocated(env%ref%wbo)) then  !> should be allocated from main program
    MDSTART%shk%wbo = env%ref%wbo
  else !> otherwise, obtain from scratch
    tmpcalc = env%calc
    mol = molstart
    tmpcalc%calcs(1)%rdwbo = .true. !> obtain WBOs
    allocate(grd(3,mol%nat))
    call engrad(mol,tmpcalc,energy,grd,io)
    call move_alloc(tmpcalc%calcs(1)%wbo, env%ref%wbo)
    deallocate(grd)
    call tmpcalc%reset()
    MDSTART%shk%wbo = env%ref%wbo
  end if
  !MDSTART%printstep = 10

  MTD%kpush = prefac
  MTD%alpha = alpha
  MTD%mtdtype = cv_rmsd
  MTD%cvdump_fs = 550.0_wp
  call MDSTART%add(MTD)

  pr = .false. !> supress stdout printout of MD

!>--- Header
  write (stdout,*)
  call smallhead('Starting trial MTD to test settings')

!>--- Iterative loop, since it is also tested if the MD runs at all
  counter = 1
  maxiter = 6
  tstep = MDSTART%tstep
  shakemode = MDSTART%shk%shake_mode
  call profiler%init(maxiter)
  iterativ: do

!>--- Restore initial starting geometry
    mol = molstart
!>--- Modify MD output trajectory
    MD = MDSTART
    MD%tstep = tstep
    MD%shk%shake_mode = shakemode
    MD%shake = (shakemode > 0)
    write (trajectory,'("crest_trial_md_",i0,".trj")') counter
    MD%trajectoryfile = trim(dirnam)//sep//trim(trajectory)

!>--- And run it
    io = 1
    !================================!
    call profiler%start(counter)
    call dynamics(mol,MD,env%calc,pr,io)
    call profiler%stop(counter)
    !================================!

    if (io == 0) then
      write (atmp,'(1x,"Trial MTD ",i0," runtime (",f3.1," ps)")') counter,MD%length_ps
      call profiler%write_timing(stdout,counter,trim(atmp))
      exit iterativ

    else

      write (stdout,'(1x,"Trial MTD ",i0," did not converge!")') counter
      if (counter >= maxiter) then
        write (stdout,'(1x,"Automatic MD restart failed ",i0," times!")') counter
        write (stdout,'(1x,"Please try other settings manually.")')
        write (stdout,*)
        error stop
      end if
      counter = counter+1

      if (tstep <= 1.0d0.and.shakemode > 0) then
        write (stdout,'(1x,"Time step can´t be reduced further!")')
        write (stdout,'(1x,"Trying to turn SHAKE off ...")')
        write (stdout,*)
        shakemode = 0
        cycle iterativ

      else if (tstep <= 1.0d0.and.shakemode == 0) then
        write (stdout,'(1x,"Automatic MTD settings check failed!")')
        write (stdout,'(1x,"Please try other settings manually.")')
        error stop
      end if

      !> don't reduce the timestep below 1 fs automatically
      tstep = max(1.0d0, (tstep-1.0d0))
      write (stdout,'(1x,"Reducing the time step to ",f4.1," fs and trying again ...")') tstep
      write (stdout,*)
    end if
!>--- End loop
  end do iterativ

!>--- transfer final settings to global settings
  env%mdstep = MD%tstep
  env%mddat%tstep = MD%tstep
  env%shake = shakemode
  env%mddat%shk%shake_mode = shakemode
  env%mddat%shake = (shakemode > 0)

!>--- If we succeded and exited the above loop, estimate runtime for multiple MTDs
  rtime = profiler%get(counter)
  call timeestimate(rtime,env%mdtime,env%nmetadyn,env%threads)

  return
!=========================================================================================!
contains
!=========================================================================================!
  subroutine timeestimate(rtime,mdtime,nmetadyn,threads)
!******************************************
!* calaculate and print the time estimate
!******************************************
    implicit none
    real(wp)  :: mdtime
    integer :: nmetadyn
    integer :: threads
    real(wp)  :: rtime,rtime2
    real(wp)  :: w,t
    character(len=80) :: stime,stime2

    rtime2 = rtime

    rtime = rtime*mdtime             !scale the runtime to the MTD length
    call time2string(rtime,stime)
    write (stime2,'(f14.1)') mdtime
    stime2 = adjustl(stime2)
    write (stdout,'(1x,a,a,a,a)') 'Estimated runtime for one MTD (',trim(stime2), &
    &                      ' ps) on a single thread: ',trim(stime)

    if (threads < nmetadyn) then
      w = float(nmetadyn)/float(threads)
      t = float(ceiling(w))
      rtime = rtime*t
      call time2string(rtime,stime)
    else if (threads == nmetadyn) then
      call time2string(rtime,stime)
    else if (threads > nmetadyn) then
      w = float(threads)/float(nmetadyn)
      t = float(floor(w))
      rtime = rtime/t
      call time2string(rtime,stime)
    end if

    write (stdout,'(1x,a,i0,a,i0,a,a)') 'Estimated runtime for a batch of ',nmetadyn, &
    &                      ' MTDs on ',threads,' threads: ',trim(stime)

  end subroutine timeestimate

  subroutine time2string(rtime,str)
    implicit none
    real(wp) :: rtime
    real(wp) :: hours,minutes,seconds
    character(len=*) :: str

    hours = aint(rtime/3600d0)
    minutes = aint((rtime-(3600d0*hours))/60d0)
    seconds = anint(rtime-(3600d0*hours+minutes*60d0))

    if (hours > 0.0d0) then
      write (str,'(i0,1x,a,1x,i0,1x,a,1x,i0,1x,a)') &
      & nint(hours),'h',nint(minutes),'min',nint(seconds),'sec'
    else if (minutes > 0.0d0) then
      write (str,'(i0,1x,a,1x,i0,1x,a)') &
      & nint(minutes),'min',nint(seconds),'sec'
    else
      write (str,'(i0,1x,a)') nint(seconds),'sec'
    end if

    return
  end subroutine time2string
!========================================================================================!
end subroutine trialMD_calculator

!========================================================================================!
!========================================================================================!
!> Routines to check topology behaviour upon optimization of input
!========================================================================================!
!========================================================================================!
subroutine trialOPT_calculator(env)
!**********************************************************
!* subroutine trialOPT_calculator
!* A new calculator-based implementation of xtbopt_legacy
!* Performs a geometry optimization of the structure
!* saved to env%ref and checks for changes in the topology
!*
!* requires to have env%ref and env%calc initialized
!**********************************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  use optimize_module
  use strucrd
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env
  !> LOCAL
  type(coord) :: mol,molopt
  type(calcdata) :: tmpcalc
  integer :: io,T,Tn
  real(wp) :: energy
  real(wp),allocatable :: grd(:,:)
  logical :: success,pr,wr

!>--- get all available threads
  call new_ompautoset(env,'max',0,T,Tn)

!>--- small header
  write (stdout,*)
  call smallhead('Initial Geometry Optimization')

!>--- setup
  call env%ref%to(mol)
  call env%ref%to(molopt)
  allocate(grd(3,mol%nat), source=0.0_wp)
  tmpcalc = env%calc  !> create copy of calculator
  tmpcalc%optlev = -1 !> set loose convergence thresholds 

!>--- perform geometry optimization
  pr = .false. !> stdout printout
  wr = .true.  !> write crestopt.log
  call optimize_geometry(mol,molopt,tmpcalc,energy,grd,pr,wr,io)

!>--- check success 
  success = (io == 0)
  call trialOPT_warning(env,molopt,success)

  deallocate(grd) 
end subroutine trialOPT_calculator

!========================================================================================!
subroutine trialOPT_warning(env,mol,success)
!*******************************************************
!* subroutine trialOPT_warning
!* Processes the trialOPT status and, if successfull,
!* overwrites env%ref with the opt. structure
!* If the checks fail, CREST is stopped by this routine
!*******************************************************
  use crest_parameters
  use crest_data
  use strucrd
  use utilities
  implicit none
  type(systemdata),intent(inout) :: env
  type(coord),intent(in) :: mol  !> optimized geometry to be checked
  logical,intent(in) :: success  !> status of geometry optimization
  integer :: ntopo
  integer,allocatable :: topo(:)
  logical,allocatable :: changed(:)
  integer(dp) :: i  !> revlin expects int64 input
  integer :: j,k,l
  logical :: tchange

  if (.not.success) then
    write (stdout,*)
    write (stdout,*) ' Initial geometry optimization failed!'
    write (stdout,*) ' Please check your input.'
    error stop
  end if
  write (stdout,*) 'Geometry successfully optimized.'
!---- if necessary, check if the topology has changed!
  if (allocated(env%ref%topo)) then
    ntopo = mol%nat*(mol%nat+1)/2
    allocate (topo(ntopo))
    allocate (changed(mol%nat),source=.false.)
    call quicktopo(mol%nat,mol%at,mol%xyz,ntopo,topo)
    tchange = .false.
    do i = 1,ntopo
      if (topo(i) .ne. env%ref%topo(i)) then
        tchange = .true.
        call revlin(i,k,l)
        changed(k) = .true.
        changed(l) = .true.
      end if
    end do
    if (tchange) then
      write (stdout,*)
      write (stdout,'(1x,a)') '*WARNING* Change in topology detected!'
      write (stdout,'(1x,a)') 'Topology change compared to the input affects atoms:'
      do i = 1,mol%nat
        if (changed(i)) then
          write (stdout,'(1x,i0,"(",a,")")',advance='no') i,trim(i2e(mol%at(i),'nc'))
        end if
      end do
      write (stdout,*)

!>--- either update the topology (see option B below)
      if (.not.env%reftopo) then
        write (stdout,'(1x,a)') 'Taking new topology as reference and continue ...'
        env%ref%topo = topo

!>--- or abort the run
      else
        write (stdout,*)
        call smallhead('* READ THE FOLLOWING CAREFULLY *')
        write (stdout,'(1x,a)') 'A topology change was seen in the initial geometry optimization.'
        write (stdout,'(1x,a,a,a)') 'This could be an artifact of the chosen theory level (e.g. xTB).'
        if (env%legacy) then
          write (stdout,'(1x,a)') 'You can check the optimization trajectory in the "xtbopt.log" file.'
        else
          write (stdout,'(1x,a)') 'You can check the optimization trajectory in the "crestopt.log" file.'
        end if
        write (stdout,'(1x,a)') 'Try either of these options:'
        write (stdout,'(/,4x,a)') 'A) Pre-optimize your input seperately and use the optimized'
        write (stdout,'(4x,a)') '   structure as input for CREST. (Only recommended if structure is intact)'
        write (stdout,'(/,4x,a)') 'B) Restart the same CREST call as before, but ignore the topology change'
        write (stdout,'(4x,a)') '   by using the "--noreftopo" keyword. (May produce artifacts)'
        write (stdout,'(/,4x,a)') 'C) Fix the initial input geometry by introducing bond length constraints'
        write (stdout,'(4x,a)') '   or by using a method with fixed topology (e.g. GFN-FF).'
        write (stdout,*)
        error stop 'safety termination of CREST'
      end if
    end if
  end if

!>--- If all checks succeded, update reference with optimized geometry
  env%ref%nat = mol%nat
  env%ref%at = mol%at
  env%ref%xyz = mol%xyz

end subroutine trialOPT_warning

!========================================================================================!
!========================================================================================!
