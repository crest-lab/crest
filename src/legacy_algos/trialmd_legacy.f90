!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Philipp Pracht
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
subroutine trialMD_legacy(env)
  use crest_parameters
  use crest_data
  use ls_rmsd
  use iomod
  use utilities
  implicit none

  type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA

  integer :: io,Vdumpfreq,counter

  real(wp) :: prefac,alpha,length
  real(wp) :: perc
  real(wp) :: rtime

  character(len=512) :: thispath,jobcall
  character(len=64)  :: dirnam,fname
  character(len=80)  :: pipe

  real(wp),allocatable :: gdum(:,:),Udum(:,:),xdum(:),ydum(:)  !rmsd dummy stuff

  logical :: shakefail,mdfail,ex
  integer :: sysio,T,Tn

  allocate (gdum(3,3),Udum(3,3),xdum(3),ydum(3))

!>--- some settings
  call new_ompautoset(env,'auto',1,T,Tn)

  call getcwd(thispath)

  dirnam = 'TRIALMD'
  fname = 'coord'

  prefac = 0.003d0*env%rednat   !Vbias k
  alpha = 0.5d0             !Vbias alpha
  Vdumpfreq = 10000         !Vbias dumpfrequency, i.e., never updated !SG
  length = 1.0d0   !MTD length in ps
  perc = 0.50d0    !which last percentage of the xtb.trj should be included in the mRMSD calculation?
  ! 50% = the last half

  pipe = ' > xtb.out 2>/dev/null'

  write (jobcall,'(a,1x,a,1x,a,'' --md '',a,a)') &
  &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),pipe

  !--- slightly different jobcall for qmdff usage, and check for solvent file
  if (env%useqmdff) then
    write (jobcall,'(a,1x,a,1x,a,'' --md --qmdff '',a,a)') &
    &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),pipe

    inquire (file='solvent',exist=ex)
    if (.not.ex) then
      write (stdout,*)
      write (stdout,*) "Use of QMDFF was requested by '-qmdff', but there is no solvent file!"
      write (stdout,*) "exit."
      error stop
    end if
  end if

!>--- Header
  write (stdout,*)
  write (stdout,'(''-------------------------------------'')')
  write (stdout,'(''Starting a trial MTD to test settings'')')
  write (stdout,'(''-------------------------------------'')')

!>--- Iterative loop, since it is also tested if the MD runs at all
  counter = 1
  iterativ: do
!>--- Make a dir
    call rmrf(dirnam)
    io = makedir(trim(dirnam))

    call copysub('coord',trim(dirnam))
    call env%wrtCHRG(trim(dirnam))
    call copysub(env%fixfile,trim(dirnam))
    if (env%useqmdff) then
      call copysub('solvent',trim(dirnam))
    end if

    call chdir(trim(dirnam))

!>--- set up MD
    call setMDrun2(fname,env%hmass,length,env%mdtemp,env%mdstep,env%shake,env%mddumpxyz, &
    &             env%mdskip,Vdumpfreq,-1,env%cts)
    call setMetadyn2(fname,prefac,alpha,100)  ! Metadynamic settings

    !================================!
    call command(trim(jobcall),sysio)
    !================================!

    inquire (file='xtbmdok',exist=ex) !this file should only exist if the MTD finished flawlessly
    call minigrep('xtb.out','SHAKE did not converge',shakefail)
    call minigrep('xtb.out','MD is unstable, emergency exit',mdfail)
    if (ex.and..not.shakefail.and..not.mdfail) then
      call gettime('xtb.out',rtime)
      call timeestimate(rtime,env%mdtime,env%nmetadyn,env%threads)

!>--- change Dir back
      call chdir(thispath)
      exit iterativ
    else
      call chdir(thispath)
      write (stdout,'(1x,''Trial MTD '',i0,'' did not converge!'')') counter
      if (counter .ge. 6) then
        write (stdout,'(1x,''Automatic xtb restart failed '',i0,'' times!'')') counter
        write (stdout,'(1x,''Please try other settings by hand.'')')
        write (stdout,*)
        error stop
      end if
      counter = counter+1
      if (env%mdstep == 1.0d0) then
        write (stdout,'(1x,''Time step can`t be reduced further!'')')
        write (stdout,'(1x,''Trying to turn SHAKE off...'')')
        write (stdout,*)
        env%shake = 0
        cycle iterativ
      end if
      env%mdstep = max(1.0d0, (env%mdstep-1.0d0)) !don't reduce the timestep below 1 fs
      write (stdout,'(1x,''Reducing the time step to '',f4.1,'' fs and trying again...'')') env%mdstep
      write (stdout,*)
    end if
!>--- End loop
  end do iterativ

  return
!=========================================================================================!
contains
!=========================================================================================!
  subroutine gettime(fil,secs)
!*********************************
!* Obtain runtime from xtb output
!*********************************
    implicit none
    real(wp) :: secs
    real(wp) :: floats(10)
    character(len=*) :: fil
    character(len=512) :: tmp,tmp2
    integer :: io,ich,n
    secs = 0.0d0
    open (newunit=ich,file=fil)
    do
      read (ich,'(a)',iostat=io) tmp
      if (io < 0) exit
      if (index(tmp,'finished run on') .ne. 0) then
        read (ich,'(a)') tmp
        read (ich,'(a)') tmp
        read (ich,'(a)') tmp
        read (ich,'(a)') tmp
        call rdarg(tmp,'time:',tmp2)
        call readl(tmp2,floats,n)
      end if
    end do
    secs = secs+floats(2)*86400.0d0  !days to seconds
    secs = secs+floats(3)*3600.0d0   !hours to seconds
    secs = secs+floats(4)*60.0d0     !minutes to seconds
    secs = secs+floats(5)
    return
  end subroutine gettime

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
end subroutine trialMD_legacy
