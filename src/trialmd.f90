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
subroutine trialMD(env)
      use iso_fortran_env, only : wp => real64
      use crest_data
      use ls_rmsd
      use iomod
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
      integer :: sysio

      allocate(gdum(3,3),Udum(3,3),xdum(3),ydum(3))

      settingOMP: associate( autothreads => env%autothreads, threads => env%threads, &
      &           omp => env%omp, MAXRUN => env%MAXRUN)

      settingData: associate( nat => env%nat, rednat => env%rednat, ProgName => env%ProgName,         &
      &            gfnver => env%gfnver, gbsa => env%gbsa, solvent => env%solvent,                    &
      &            constraints => env%constraints, fixfile => env%fixfile, useqmdff => env%useqmdff,  &
      &            hmass => env%hmass, mdtemp => env%mdtemp, mdstep => env%mdstep, shake => env%shake,&
      &            mddumpxyz => env%mddumpxyz, mdskip => env%mdskip)

!---- some settings
      if(autothreads)then
         call ompautoset(threads,8,omp,MAXRUN,8) !use maximum number of threads for the xtb job, with a maximum of 8
      endif

      call getcwd(thispath)

      dirnam='MRMSD'
      fname='coord'

      prefac=0.003d0*rednat   !Vbias k
      alpha=0.5d0             !Vbias alpha
      Vdumpfreq=10000         !Vbias dumpfrequency, i.e., never updated !SG
      length=1.0d0   !MTD length in ps
      perc=0.50d0    !which last percentage of the xtb.trj should be included in the mRMSD calculation?
                     ! 50% = the last half


      pipe=' > xtb.out 2>/dev/null'

      write(jobcall,'(a,1x,a,1x,a,'' --md '',a,a)') &
      &     trim(ProgName),trim(fname),trim(gfnver),trim(env%solv),pipe 

      !--- slightly different jobcall for qmdff usage, and check for solvent file
      if(useqmdff)then
       write(jobcall,'(a,1x,a,1x,a,'' --md --qmdff '',a,a)') &
       &     trim(ProgName),trim(fname),trim(gfnver),trim(env%solv),pipe

       inquire(file='solvent',exist=ex)
       if(.not.ex)then
          write(*,*)
          write(*,*) "Use of QMDFF was requested by '-qmdff', but there is no solvent file!"
          write(*,*) "exit."
          error stop
       endif
      endif

!---- Header
       write(*,*)
       write(*,'(''-------------------------------------'')')
       write(*,'(''Starting a trial MTD to test settings'')')
       write(*,'(''-------------------------------------'')')



!---- Iterative loop, since it is also tested if the MD runs at all
      counter = 1
      iterativ: do
!---- Make a dir
      call rmrf(dirnam)
      io = makedir(trim(dirnam))

      call copysub('coord',trim(dirnam))
      call env%wrtCHRG(trim(dirnam))   
      call copysub(env%fixfile,trim(dirnam))
      if(useqmdff)then
        call copysub('solvent',trim(dirnam))
      endif
      
      call chdir(trim(dirnam))

!---- set up MD
      call setMDrun2(fname,hmass,length,mdtemp,mdstep,shake,mddumpxyz, &
      &             mdskip,Vdumpfreq,-1,env%cts)
      call setMetadyn2(fname,prefac,alpha,100)  ! Metadynamic settings
 
      call command(trim(jobcall), sysio)

      inquire(file='xtbmdok',exist=ex) !this file should only exist if the MTD finished flawlessly
      call minigrep('xtb.out','SHAKE did not converge',shakefail)
      call minigrep('xtb.out','MD is unstable, emergency exit',mdfail)
      if(ex.and. .not.shakefail .and. .not.mdfail)then
      call gettime('xtb.out',rtime)
      call timeestimate(rtime,env%mdtime,env%nmetadyn,env%threads)

      
!---- change Dir back
         call chdir(thispath)
         exit iterativ
      else
         call chdir(thispath) 
         write(*,'(1x,''Trial MTD '',i0,'' did not converge!'')')counter
         if(counter.ge.6)then
           write(*,'(1x,''Automatic xtb restart failed '',i0,'' times!'')')counter
           write(*,'(1x,''Please try other settings by hand.'')')
           write(*,*)
           error stop 
         endif
         counter = counter + 1
         if(mdstep==1.0d0)then
            write(*,'(1x,''Time step can`t be reduced further!'')')
            write(*,'(1x,''Trying to turn SHAKE off...'')')
            write(*,*)
            shake=0
            cycle iterativ
         endif
         mdstep = max(1.0d0,(mdstep - 1.0d0)) !don't reduce the timestep below 1 fs
         write(*,'(1x,''Reducing the time step to '',i0,'' fs and trying again...'')')nint(mdstep)
         write(*,*)
      endif
!---- End loop
      enddo iterativ


      end associate settingData
      end associate settingOMP

end subroutine trialMD

!--------------------------------------------------------------------------------------------------------------
! calculate mRMSD
!--------------------------------------------------------------------------------------------------------------
subroutine c_mRMSD(rmsds,n,mRMSD)
      implicit none
      integer :: n
      real*8 ::  rmsds(n)     !mean RMSD
      real*8 :: mRMSD
     
      mRMSD=sum(rmsds)/float(n)

      return
end subroutine c_mRMSD

!--------------------------------------------------------------------------------------------------------------
! calculate f(mRMSD) with a Logistic funciton
! f(mRMSD) = ul/(1+exp(-steep*(mRMSD-shift))) + lb    with ul=(ub-lb)
!--------------------------------------------------------------------------------------------------------------
function f_mRMSD(mRMSD,lb,ub,steep,shift)
      implicit none
      real*8 ::  f_mRMSD  !f(mRMSD) retrun value
      real*8 ::  mRMSD    !mean RMSD
      real*8 ::  lb       !lower limit
      real*8 ::  ub       !upper limit
      real*8 ::  steep    !steepness of the function
      real*8 ::  shift    !mRMSD shift of the function
      real*8 :: ul
      ul=ub-lb
      f_mRMSD = ul / (1+ exp(-steep*(mRMSD - shift))) + lb
      return
end function f_mRMSD

!--------------------------------------------------------------------------------------------------------------
! calaculate and print the time estimate
!--------------------------------------------------------------------------------------------------------------
subroutine timeestimate(rtime,mdtime,nmetadyn,threads)
     implicit none
     real*8  :: mdtime
     integer :: nmetadyn
     integer :: threads
     real*8  :: rtime,rtime2
     real*8  :: w,t
     character(len=80) :: stime,stime2
 
     rtime2=rtime

     rtime=rtime*mdtime             !scale the runtime to the MTD length
     call time2string(rtime,stime)
     !write(*,*) trim(stime)
     write(stime2,'(f14.1)')mdtime
     stime2=adjustl(stime2)
     write(*,'(1x,a,a,a,a)')'Estimated runtime for one MTD (',trim(stime2), &
     &                      ' ps) on a single thread: ',trim(stime)

     if(threads < nmetadyn)then
        w=float(nmetadyn)/float(threads)
        t=float(ceiling(w))
        rtime=rtime*t
        call time2string(rtime,stime)
     else if(threads == nmetadyn)then
        call time2string(rtime,stime)
     else if(threads > nmetadyn)then
        w=float(threads)/float(nmetadyn)
        t=float(floor(w))
        rtime=rtime/t
        call time2string(rtime,stime)
     endif    

     write(*,'(1x,a,i0,a,i0,a,a)')'Estimated runtime for a batch of ',nmetadyn, &
     &                      ' MTDs on ',threads ,' threads: ',trim(stime)


end subroutine timeestimate

subroutine time2string(rtime,str)
     implicit none
     real*8 :: rtime
     real*8 :: hours,minutes,seconds
     character(len=*) :: str

     hours   =  aint( rtime / 3600d0)
     minutes =  aint( (rtime-(3600d0*hours))/60d0 )
     seconds = anint( rtime-(3600d0*hours + minutes*60d0) )

     if(hours>0.0d0)then
       write(str,'(i0,1x,a,1x,i0,1x,a,1x,i0,1x,a)') &
       & nint(hours),'h',nint(minutes),'min',nint(seconds),'sec'  
     else if(minutes>0.0d0)then
       write(str,'(i0,1x,a,1x,i0,1x,a)') &
       & nint(minutes),'min',nint(seconds),'sec'
     else
       write(str,'(i0,1x,a)')nint(seconds),'sec'
     endif

     return
end subroutine time2string


