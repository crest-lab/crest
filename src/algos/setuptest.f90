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
subroutine trialMD_calculator(env)
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

  character(len=*),parameter :: dirnam = 'TRIALMD'

!>--- OMP settings (should be set to 1 to simulate max parallelization)
  if (env%autothreads) then
    call ompautoset(env%threads,8,env%omp,env%MAXRUN,8)
  end if

  call getcwd(thispath)

!>--- Make a dir
  call rmrf(dirnam)
  io = makedir(trim(dirnam))

!>--- set up mol and MD calculator
  prefac = 0.003d0*env%rednat !> Vbias k
  alpha = 0.5d0         !> Vbias alpha
  !Vdumpfreq = 10000     !> Vbias dumpfrequency, i.e., never updated !SG

  call env%ref%to(molstart)
  MDSTART = env%mddat         !> env%mddat should already be set up at this point
  MDSTART%length_ps = 1.0_wp  !> we want only a short simulaton
  MDSTART%length_steps = 0    !> to zero so it will be calculated automatically
  MDSTART%dumpstep = 20.0_wp  !> fs dump step to trajectory file, so we end up with 50 structures
  MDSTART%sdump = 0           !> also zero to reset
  if(allocated(env%ref%wbo))then  !> should be allocated from main program
    MDSTART%shk%wbo = env%ref%wbo
  endif
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
    write(trajectory,'("crest_trial_md_",i0,".trj")') counter
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
