!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Stefan Grimme
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
subroutine xtbmodef(env,s,e,n,u,l,o,x)
      use iso_fortran_env, only : wp => real64    
      use setdgmod
      use iomod
      use crest_data
      implicit none
      type(systemdata) :: env
      !type(options)    :: opt
      integer :: x  !1=normal mode, 2=local mode
      integer :: s  !start mode
      integer :: e  !end mode
      integer :: n  !number of points on mode
      integer :: o  !opt cycles per point per mode
      integer :: i,k
      integer :: TOTAL,MINM,MAXM
      integer :: DONE,dir,ACTIVE,MODE,CHECK
      integer :: io,r
      real(wp)  :: u !mode update factor
      real(wp)  :: l !steplength between points on mode

      character(len=20):: modefile1,modefile2
      character(len=20)::extn
      character(len=80):: solv,fname,outf,pipe
      character(len=512):: str,thispath,tmppath,dg
      character(len=1024):: jobcall
      !character(len=:), allocatable :: ProgName, gfnver

      logical:: l2,update,ex,clo
      logical :: niceprint

      real(wp) :: percent
      character(len=52) :: bar


      call getcwd(thispath)  !get current working directory

      associate( gbsa => env%gbsa, solvent => env%solvent,  &
      & fixfile => env%fixfile, constraints => env%constraints, optlev => env%optlev, &
      & trackorigin => env%trackorigin, gfnver => env%gfnver, ProgName => env%ProgName)
      niceprint = env%niceprint

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      solv=''
      extn=''
      modefile1='xtb_normalmodes'
      modefile2='xtb_localmodes'
      fname='coord'
      clo=.false.
      outf='xtb.out'
      update=.true. !print the current mode


      if(gbsa)then
         write(solv,'(''--gbsa '',a)')trim(solvent)
      else
         solv=''
      endif

      pipe='2>/dev/null'


      TOTAL= e - s +1
      MINM=1
      MAXM=500

      if(env%autothreads)then
        call ompautoset(env%threads,7,env%omp,env%MAXRUN,TOTAL) !set the global OMP/MKL variables for the xtb jobs
      endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !error section
      
      if(e < s)then
         error stop 'Stopped. Someting is wrong with the mode recommendation.'
      end if

      inquire(file=modefile1,exist=l2)
      if(.not.l2)then
         error stop 'can not find file <xtb_*modes>'
      end if

      if(TOTAL > MAXM)then
         error stop 'too many modes. EXIT.'
      end if

      if(TOTAL < MINM)then
         error stop '< 1 modes. EXIT.'
      end if

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !setting up the coord file correctly
      
      if(x.eq.0)then
         write(*,*)'setting defaults for PES scan for anharmonicity'
         call setdg(fname,'mode_local','-1')
      end if
      
      inquire(file=constraints,exist=ex)
      if(ex)then
         call clear_setblock(fname)
      endif

      if(x.eq.1)then
         write(*,*)'setting defaults for canonical NC search'
         call setdg(fname,'mode_updat',u)
         call setdg(fname,'mode_n',n)
         call setdg(fname,'mode_step',l)
         call setdg_optlev(fname,optlev)
         call setdg(fname,'mode_local','0')
      end if

      if(x.eq.2)then
         write(*,*)'setting defaults for local NC search'
         call setdg(fname,'mode_updat',u)
         call setdg(fname,'mode_n',n)
         call setdg(fname,'mode_step',l)
         call setdg_optlev(fname,optlev)
         call setdg(fname,'mode_local','1')
      end if

      call setdg(fname,'maxopt',o)

      call ConstrainsToSET(fname,clo,env%constraints)

      write(*,'(''  number of points on each mode  :'',i6)')n
      write(*,'(''  step length between points     :'',f9.4)')l
      write(*,'(''  structure update factor (0-1)  :'',f9.4)')u
      write(*,'(''  max. opt cycles per point      :'',i6)')o

      write(*,'(1x,''# Modes: '',i0,1x,''first: '',i0,1x,''last: '',i0)')TOTAL,s,e
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !set up MODEF* directories

      DONE=0
      do while(DONE.ne.TOTAL)
         DONE=DONE + 1
         dir=DONE + s - 1
         write(tmppath,'(''MODEF'',i0)')dir
         inquire(file=tmppath,exist=l2)
         if(x.lt.2)then !don't delete for local mode run
          call rmrf(tmppath)
          r = makedir(trim(tmppath))
         endif
         call copysub('coord',trim(tmppath))
         !call copysub('.CHRG',trim(tmppath))
         !call copysub('.UHF',trim(tmppath))
         call env%wrtCHRG(trim(tmppath))   
         call copysub(fixfile,trim(tmppath))
         call copysub(constraints,trim(tmppath))

         call chdir(trim(tmppath))
         call getcwd(tmppath)
         call clearmodef()

          str=trim(thispath)//'/'//'xtb_normalmodes'
          io = sylnk(trim(str),'xtb_normalmodes')

          str=trim(thispath)//'/'//'xtb_localmodes'
          io = sylnk(trim(str),'xtb_localmodes')

         call chdir(thispath)
      end do
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !main loop
      CHECK=0
      DONE=0
      ACTIVE=0
      
      k=0
      if(niceprint)then
         call printemptybar()
      else
          write(6,'(1x,i0)',advance='no') 0
          flush(6)
      endif
!$omp parallel &
!$omp shared( MODE,TOTAL,percent,k,niceprint,bar ) &
!$omp shared( env, solv,outf,pipe )
!$omp single
      do i=1,TOTAL
         MODE= i + s - 1
      !$omp task firstprivate( MODE ) private( tmppath,io )
         write(tmppath,'(a,i0)')'MODEF',MODE
         write(jobcall,'(a,'' coord '',a,'' -modef '',i0,1x,a,'' > '',a,1x,a)') &
         & trim(env%ProgName),trim(env%gfnver),MODE,trim(solv),trim(outf),trim(pipe)
         !call system('cd '//trim(tmppath)//' && '//trim(jobcall))
         call execute_command_line('cd '//trim(tmppath)//' && '//trim(jobcall), exitstat=io)
      !$omp critical
        k=k+1
        if(niceprint)then
          percent=float(k)/float(TOTAL)*100
          call  progbar(percent,bar)
          call printprogbar(percent,bar)
        else
          write(6,'(1x,i0)',advance='no')MODE
          flush(6)
        endif
      !$omp end critical
      !$omp end task
      enddo
!$omp taskwait
!$omp end single
!$omp end parallel

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(trackorigin)then
      DONE=0
      do while(DONE.ne.TOTAL)
         DONE=DONE + 1
         dir=DONE + s - 1
         write(tmppath,'(''MODEF'',i0)')dir
         call chdir(tmppath)
         do i=1,n
          if(x.lt.2)then
             write(str,'(a,i0,a,i0,a)')'xtb_modemin_',i,'_',dir,'.xyz'
             inquire(file=trim(str),exist=ex)
             if(ex)then
               write(dg,'(a,i0)')'mode',dir
               call addorigin(trim(str),trim(dg))
             endif
          else
             write(str,'(a,i0,a,i0,a)')'xtb_modemin_',i,'_',dir+1000,'.xyz'
             inquire(file=trim(str),exist=ex)
             if(ex)then
               write(dg,'(a,i0)')'locmode',dir
               call addorigin(trim(str),trim(dg))
             endif
          endif
         enddo
         call chdir(thispath)
      end do
      endif


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call chdir(thispath)

      write(*,*)
      write(*,*)'done.'

      end associate

end subroutine xtbmodef



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine clearmodef

      use iomod
      implicit none

      call remove('xtbmodefok')
      call remove('RUNNING')
      call remove('FINISHED')
end subroutine clearmodef
