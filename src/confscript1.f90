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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Standalone Version of the confscript1 code (OUTDATED ALGO)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine confscript1(env,tim)
      use iso_fortran_env, wp => real64
      use iomod
      use crest_data
      implicit none

      type(systemdata) :: env
      type(timer)      :: tim
      integer :: nmodes
      character(len=512) :: str

      logical :: better

      associate( performModef => env%performModef, performMD => env%performMD, &
      & performCross => env%performCross, autothreads => env%autothreads, doNMR => env%doNMR, &
      & quick => env%quick, keepModef => env%keepModef)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         CONFORMATION SEARCH CALLS START HERE
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      env%icount = 0 
      write(*,*)'starting with hessian calculation...'

!c main loop                  
!c the logical better indicates if a better structure has been found
! 9999 continue

      do
      call clean  ! remove dirs and confg files

      call tim%start(1,'hessian')
      call ohess(env)
      call tim%stop(1)

      env%icount = env%icount + 1

      nmodes=min(nmodes,200)    ! limit to 200 modes

      if(performModef)then
         call tim%start(2,'MF')
         call modef(env,better)
         call tim%stop(2)
         if(better.and.env%icount.lt.env%Maxrestart) cycle ! goto 9999  ! avoid infinite loop
      endif
      if(performMD)then
         call tim%start(3,'MD')
         call md(env,better)
         call tim%stop(3)
         if(better.and.env%icount.lt.env%Maxrestart) cycle !goto 9999  ! avoid infinite loop
      endif
      if(performCross)then
         if(autothreads)then
           call ompautoset(env%threads,2,env%omp,env%MAXRUN,0) !max. number of Jobs
         endif
         call tim%start(4,'GC')
         if(.not.performMD)call confg_chk(env,0,better)
         call cross(env,better)
         call tim%stop(4)
         if(better.and.env%icount.lt.env%Maxrestart) cycle !goto 9999  ! avoid infinite loop
      endif
      exit
      enddo

      call tim%start(5,'')
      if(.not.performCross)then
         !if(doNMR) cgf(3)=.true.
         call confg_chk(env,0,better)
      endif

      if(doNMR)then
         env%cgf(3)=.true.
         env%cgf(2)=.false.
         call confg_chk(env,0,better)
      endif
      call tim%stop(5)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!C last bit of formal stuff
      str='cregen.out.tmp'
      call catdel(trim(str))

      if(.not.keepModef)then
        call clean()
      endif

      end associate
!c end of confscript1 routine
end subroutine confscript1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c RUN THE MODE FOLLOWING 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine modef(env,better)
      use iso_fortran_env, wp => real64
      use crest_data
      implicit none

      type(systemdata) :: env
      logical better
      integer start_mode,end_mode,end_model

      start_mode=7
      end_mode  = start_mode+env%nmodes-1
      end_model = start_mode+2*env%nmodes/3-1  ! 2/3 number of local modes
      !if(end_model<start_mode)end_model=start_mode
      end_mode = max(start_mode,end_mode)
      end_model= max(start_mode,end_model)

      write(*,*)
      !write(*,'(5x,''========================================'')')
      !write(*,'(5x,''|         Modefollowing (MF)           |'')')
      !write(*,'(5x,''========================================'')')
   
      write(*,'(''-------------------'')')
      write(*,'(''MODEFOLLOWING (MF)'')')
      write(*,'(''-------------------'')')



      if(.not.env%quick)then
!     cheaper cycles first

      ! 3x normal mode following
        call xtbmodef(env,start_mode,end_mode,21,0.5d0,2.0d0,5,1)      
        call confg_chk(env,1,better)
        if(better) return

        call xtbmodef(env,start_mode,end_mode,31,0.0d0,0.5d0,5,1)
        call confg_chk(env,1,better)
        if(better) return

        call xtbmodef(env,start_mode,end_mode,41,0.1d0,3d0,15,1)
        call confg_chk(env,1,better)
        if(better) return


      ! 2x local mode following
        call xtbmodef(env,start_mode,end_model,21,0.5d0,1.0d0,5,2)
        call confg_chk(env,1,better)
        if(better) return

        call xtbmodef(env,start_mode,end_model,41,0.1d0,2.0d0,5,2)
        call confg_chk(env,1,better)
        if(better) return

      else !if quick search is activated do only one singe NMF

        call xtbmodef(env,start_mode,end_mode,41,0.1d0,3.0d0,15,1)
        call confg_chk(env,1,better)
        if(better) return

      endif
end subroutine modef

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c RUN THE MD                
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine md(env,better)
      use iso_fortran_env, wp => real64
      use crest_data
      implicit none

      type(systemdata) :: env

      logical :: better
      real(wp) :: time
      character(len=512) :: atmp,btmp
      character(len=1024) :: str
      integer :: io
      
      write(*,*)
      write(*,'(''------------'')')
      write(*,'(''MD SAMPLING'')')
      write(*,'(''------------'')')
      
      time = env%mdtime

      !if(env%mdmode.eq.1)then
      !  idump=nint(5000.0d0*env%level)
      !else
        env%mddumpxyz=nint( (time*1000.0d0) / real(env%snapshots))
      !endif
      !write(*,*) idump
!c set MD parameters in subroutine <xtbconfmd.f>
      !call confmd(env%mdmode,time,idump)
      call normalMD_para_OMP(env,2,env%temps)
      call multilevel_opt(env,6)

      call checkname_xyz(crefile,atmp,btmp)
      write(str,'(''cat MODEF*/xtb_modemin_*.xyz >> '',a)')trim(atmp)
      !call system(trim(str))
      call execute_command_line(trim(str), exitstat=io)

      call confg_chk(env,0,better)

      return
end subroutine md

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c RUN THE CROSSING        
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine cross(env,better)     
      use iso_fortran_env, wp => real64
      use crest_data
      use iomod
      implicit none

      type(systemdata) :: env

      logical :: better
      integer :: i,imax,tmpconf
      character(len=128) :: inpnam,outnam
      character(len=512) :: thispath

      write(*,*)
      write(*,'(''------------------------'')')
      write(*,'(''STRUCTURE CROSSING (GC)'')')
      write(*,'(''------------------------'')')

      call getcwd(thispath)

      imax=min(30*env%nat,2000)  ! not too many

      do i=1,1               
          call confcross(env,imax,tmpconf)
          if(tmpconf.lt.1)then
             return
             exit
          endif
          call chdir('OPTIM')
          call confopt(env,'confcross.xyz',tmpconf,.true.)
          call chdir(trim(thispath))
          call checkname_xyz(crefile,inpnam,outnam)
          call appendto('OPTIM'//'/'//'opt.xyz',trim(inpnam))
          call rmrf('TMPCONF* OPTIM')

          call confg_chk(env,0,better)
          if(better) return                               
      enddo
end subroutine cross

!ccccccccccccccccccccccccccccccccccccccc
!c RUN THE HESSIAN    
!ccccccccccccccccccccccccccccccccccccccc

subroutine ohess(env)
      use crest_data
      use iomod
      use strucrd, only: coord2xyz
      implicit none

      type(systemdata) :: env

      integer :: io

      integer :: iost,ich
      character(len=512) :: str,atmp
      character(len=868) :: tmp
      character(len=1024):: jobcall
      real*8 tmpx
      logical :: ex
      character(len=:),allocatable :: pipe

      if(env%autothreads)then
        call ompautoset(env%threads,1,env%omp,env%MAXRUN,0) !max. number of OMP Threads
      endif

      open(newunit=ich,file='coord')
      do
        read(ich,'(a)',iostat=iost)atmp
        if(iost < 0)exit
        if(index(atmp,'$coord').ne.0) cycle
        if(index(atmp,'$').ne.0)then
          exit
        endif
      enddo
      write(ich,'(a)')'$modef'
      write(ich,'(2x,a)') 'vthr=300.00'
      write(ich,'(a)') '$opt'
      write(ich,'(2x,''optlevel='',i6)')nint(env%optlev)
      write(ich,'(a)')'$write'
      write(ich,'(2x,a)')'modef=true'
      call write_cts(ich,env%cts)
      write(ich,'(a)') '$end'
      close(ich)

!c Initial geometry optimization and freq. calculation      
      str=''
      tmp=''
      pipe=' > xtb.out 2>xtb.out'
      write(tmp,'(a,'' coord '',a,'' --ohess'')') & 
      & trim(env%ProgName),trim(env%gfnver)
      if(env%gbsa) then
         write(jobcall,'(a,'' --gbsa '',a)') &
         & trim(tmp),trim(env%solvent)
      else
         jobcall=trim(tmp)
      endif
      !call system(trim(jobcall)//pipe)
      call execute_command_line(trim(jobcall)//pipe, exitstat=io)

      call remove("g98_canmode.out")
      call remove("g98_locmode.out")
      call remove("g98.out")

      call rename('xtbopt.coord','coord')
      call ConstrainsToSET('coord',.false.,env%constraints)
      if(env%gbsa) then  ! makes charges for qmdff
         jobcall=''
         tmp=trim(env%ProgName)//' coord --scc '
         write(jobcall,'(a,a,'' --gbsa '',a,'' > tmp 2>/dev/null'')') &
         & trim(tmp),trim(env%gfnver),trim(env%solv)
      else
         write(jobcall,'(a,'' coord --scc '',a,'' > tmp 2>/dev/null'')') &
         & trim(env%ProgName),trim(env%gfnver)
      endif
      !call system(trim(jobcall))
      call execute_command_line(trim(jobcall), exitstat=io)


      if(env%mdmode.eq.1)then
        !call system('qmdff coord > qmdff.out')
        call execute_command_line('qmdff coord > qmdff.out', exitstat=io)
      endif

      ex=.false.
      inquire(file='tmpxx',exist=ex) 
      if(.not.ex)then
         call grepval('xtb.out','recommended # of modes for mode following         ',ex,tmpx)
         if(ex)then
            open(newunit=ich,file='tmpxx')
            write(ich,*) nint(tmpx)
            close(ich)
         endif
      endif

      open(unit=1,file='tmpxx')
      read(1,*) env%nmodes
      close(1,status='delete')
      if(env%nmodes.lt.1) error stop 'somethings wrong with Hessian calculation'

      env%nmodes=nint(float(env%nmodes)/env%level) ! search depth


      write(*,*) '####### (RE)START #######'
      write(*,*) 'nmodes : ', env%nmodes
      write(*,*) '#atoms : ', env%nat   
      
      if(env%icount.eq.0) then
         call coord2xyz('coord','xtbopt.xyz')
         call appendto('xtbopt.xyz','.history.xyz')
         call remove('xtbopt.xyz')
      endif
end subroutine ohess

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c RUN THE CONFG CHECK    
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine confg_chk(env,m,better)
      use crest_data                      
      use iomod
      implicit none

      type(systemdata) :: env
      !type(options)    :: opt

      integer :: m,io
      logical :: better
      character(len=512) :: str,fname
      character(len=512) :: atmp,btmp

      fname='cregen.out.tmp'
      better=.false.
      if(env%autothreads)then
         call ompautoset(env%threads,4,env%omp,env%MAXRUN,0) !mode=4 --> Program intern Threads max
      endif

      if(m.eq.1)then
         call rmrfw(crefile)
         call checkname_xyz(crefile,atmp,btmp)
         write(str,'(''cat MODEF*/xtb_modemin_*.xyz >> '',a)')trim(atmp)
         !call system(trim(str))
         call execute_command_line(trim(str), exitstat=io)
         env%cgf(2)=.false. !turn off NEWFILE for confg call
         if(.not.env%newcregen)then
             call cregen2(env)
         else
             call newcregen(env,0)
         endif
         env%cgf(2)=.true.  !turn it on again
      else
         if(.not.env%newcregen)then
             call cregen2(env)
         else
             call newcregen(env,0)
         endif
      endif

      inquire(file='LOWER_FOUND',exist=better)

      if(better) then
         write(*,*) 'LOWER ENERGY STRUCTURE FOUND AND RESET'
         call getcregenwarning(trim(fname))
         call rename('scoord.1','coord')
         call ConstrainsToSET('coord',.false.,env%constraints)
         call appendto('crest_best.xyz','.history.xyz')
         call remove('LOWER_FOUND')
      endif

      if(env%autothreads)then
         call ompautoset(env%threads,5,env%omp,env%MAXRUN,0) !mode=5 --> Program intern Threads min
      endif

end subroutine confg_chk

