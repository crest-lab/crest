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

!================================================================================!
! some very simple drivers to execute turbomole calculations via xtb
!================================================================================!

!--------------------------------------------------------------------------------!
! read the required remote control (rc) file.
subroutine dftrc_reader(env,b973c)
      use iso_fortran_env, wp => real64
      use crest_data
      use iomod
      use filemod
      implicit none

      type(systemdata) :: env
      !type(options)    :: opt
      type(filetype)   :: rc
      logical,intent(in) :: b973c


      integer :: i,j,k,l
      integer :: n

      logical :: ex,ex2
      character(len=:),allocatable :: rcfile
      character(len=:),allocatable :: atmp

      if(allocated(env%dftrcfile))then
      rcfile=env%dftrcfile
      else
      rcfile='~/.dftrc'
      endif 
      inquire(file=rcfile,exist=ex)
      inquire(file='~/.dftrc',exist=ex2)
      if(.not.ex .and. ex2) rcfile='~/.dftrc' !fall back to .dftrc if the given rc file is non existent

!----- defaults
      env%dftprog=1 ! =TM
      env%dftsetup='cefine -func b973c -sym c1'
      env%dftruntype=1 !=OPT
      if(b973c)return !for a quick setup of B97-3c exit the routine

!---- read RC file
      inquire(file=rcfile,exist=ex)
      if( ex )then
        env%dftrcfile=rcfile
        call rc%open(rcfile)    
        n=rc%nlines
!---- first the QM program
        PROG : do i=1,n
           if(index(rc%line(i),'TURBOMOLE').ne.0)then
              env%dftprog=1 ! =TM
              exit PROG
           endif
           if(index(rc%line(i),'ORCA').ne.0)then
              write(*,'(1x,a)')'No ORCA support yet. exit'
              error stop
              env%dftprog=2 ! =ORCA
              exit PROG
           endif
        enddo PROG
!---- then, depending on the program a setup input line
        SETUP :select case( env%dftprog )
         case( 1 ) !TURBOMOLE
            TMSET : do i=1,n
               if(index(rc%line(i),'SETUP:').ne.0)then
                  env%dftsetup=rc%line(i)
                  call clinex(env%dftsetup,1)
                  env%dftsetup=adjustl(env%dftsetup)
               endif
               if(index(rc%line(i),'SETUP2:').ne.0)then
                  env%dftsetup2=rc%line(i)
                  call clinex(env%dftsetup2,1)
                  env%dftsetup2=adjustl(env%dftsetup2)
                  env%resetsetup=.true.
               endif
            enddo TMSET
         case( 2 ) !ORCA
          write(*,'(1x,a)')'No ORCA support yet. exit'
          error stop
         case default
           continue
        end select SETUP
!---- and a jobtype
        do i=1,n
         if(index(rc%line(i),'JOBTYPE:').ne.0 .or. &
         &  index(rc%line(i),'RUNTYPE:').ne.0)then
             if(index(rc%line(i),'OPT'   ).ne.0)env%dftruntype=1
             if(index(rc%line(i),'SINGLE').ne.0)env%dftruntype=2
             if(index(rc%line(i),'SP'    ).ne.0)env%dftruntype=2
             if(index(rc%line(i),'FREQ'  ).ne.0)env%dftruntype=3
             if(index(rc%line(i),'AUTOIR').ne.0)env%dftruntype=4
             k=i
         endif
        enddo
!---- cross-check with property runtypes
      call DFTcrosscheck(env)
!--- possible opt level, only used for OPT,FREQ,AUTOIR
        if( any( (/1,3,4/)==env%dftruntype) )then
             if(index(rc%line(k),'CRUDE'    ).ne.0)env%dftoptlev='crude'
             if(index(rc%line(k),'LOOSE'    ).ne.0)env%dftoptlev='loose'
             if(index(rc%line(k),'NORMAL'   ).ne.0)env%dftoptlev='normal'
             if(index(rc%line(k),'TIGHT'    ).ne.0)env%dftoptlev='tight'
             if(index(rc%line(k),'VERYTIGHT').ne.0)env%dftoptlev='verytight'
        endif
!--- specialized mass scaling for automated IR averaging
        if( env%dftruntype==4 )then
          do i=1,n
           if(index(rc%line(i),'MASSSCAL'    ).ne.0)then
              env%dftmasscal =.true.
              atmp = getlarg(rc%line(i),2)
              if(trim(atmp).ne.'')then
                env%dftmasspar=atmp
              else
                env%dftmasspar='b3lyp-3c'
              endif
           endif
          enddo
          !write(*,*) env%dftmasspar
        endif

        call rc%close  !close file
      endif

!---- printout 
      call printDFTsettings(env)
     
      return
end subroutine dftrc_reader
!----------------------------------------------------------------------------------------!
subroutine DFTcrosscheck(env)
      use iso_fortran_env, wp => real64
      use crest_data
      implicit none
      type(systemdata) :: env
      select case(env%properties)
        case( 3,5 )
          if(env%dftruntype.ne.1)then
             env%dftruntype=1
          endif
        case( 4,6 )
          if(env%dftruntype.ne.4)then
             env%dftruntype=4
          endif
        case( 7 )
          if(env%dftruntype.ne.4)then
             env%dftruntype=2
          endif
        case( 8 )
          if(env%dftruntype.ne.4)then
             env%dftruntype=3
          endif
        case default
          continue
      end select
      return
end subroutine DFTcrosscheck
!---------------------------------------------------------------------------------------!
subroutine printDFTsettings(env)
      use iso_fortran_env, wp => real64
      use crest_data
      implicit none
      type(systemdata) :: env
      character(len=:),allocatable :: run
      write(*,'(1x,a)')'=============================================='
      write(*,'(1x,a)')'The following DFT settings have been chosen:'
      write(*,'(1x,a)')'----------------------------------------------'
      if(allocated(env%dftrcfile))then
      write(*,'(1x,a,a)')   'File   : ',env%dftrcfile
      endif
      PROG : select case( env%dftprog )
       case( 1 )
         write(*,'(1x,a)')  'Program: TURBOMOLE'
       case( 2 )
         write(*,'(1x,a)')  'Program: ORCA'
       case default
         write(*,'(1x,a)')'error: no DFT settings.'
         return
      end select PROG

         write(*,'(1x,a,a)')'Setup  : ',env%dftsetup
      if( allocated(env%dftsetup2) .and. any((/3,4/)==env%dftruntype))then
         write(*,'(1x,a,a)')'Setup2 : ',env%dftsetup2
      endif

      RUNTYPE : select case( env%dftruntype )
       case( 1 )
         run='Runtype: geometry optimization'
         if( allocated(env%dftoptlev) )then
           run=run//', '//env%dftoptlev//' opt.'
         endif
       case( 2 )
         run='Runtype: singlepoint calculation'
       case( 3 )
         run='Runtype: frequency calculation'
         if( allocated(env%dftoptlev) )then
           run=run//', '//env%dftoptlev//' opt.'
         endif
       case( 4 )
         run='Runtype: autoIR procedure'
         if( allocated(env%dftoptlev) )then
           run=run//', '//env%dftoptlev//' opt.'
         endif

       case default
         run='Runtype: none'
      end select RUNTYPE
      write(*,'(1x,a)') run

      if( env%dftmasscal )then
         write(*,'(1x,a)')'IR mass param: ON'
      endif


      write(*,'(1x,a)')'=============================================='
      write(*,*)

      WARN : select case( env%dftprog )
       case( 1 )
         call dftTMwarning
       case( 2 )
         call dftORCAwarning
       case default
         return
      end select WARN
      return
end subroutine printDFTsettings
!--------------------------------------------------------------------------------------!
subroutine xtbDFTdriver(env,xname,jobcall)
      use iso_fortran_env, wp => real64
      use crest_data
      implicit none
      type(systemdata) :: env
      character(len=*),intent(inout) :: jobcall
      character(len=*),intent(in)  :: xname
      character(len=:),allocatable :: runtype
      character(len=:),allocatable :: qmprog

      RUN : select case( env%dftruntype )
       !case( 1 )
       !  runtype='--opt'
       case( 2 )
         runtype='--sp'
       case( 1,3:4 )
         !runtype='--ohess'
         runtype='--opt'  !not '--ohess' since only the optimization needs to be dirven by xtb
         if(allocated(env%dftoptlev))then
            runtype='--opt '//env%dftoptlev
         endif
       case default
         runtype=''
      end select RUN

      QM : select case( env%dftprog )
       case( 1 )
         qmprog='--tm'
       case( 2 )
         qmprog='--orca'
       case default
         qmprog=''
      end select QM

      write(jobcall,'(a,1x,a,1x,a,1x,a,'' >xtb.out 2>/dev/null'')') &
      &    trim(env%ProgName),trim(xname),runtype,qmprog
      return
end subroutine xtbDFTdriver
!------------------------------------------------------------------------------------------------------
subroutine DFTprocessing(env,TMPCONF,nat,at)
      use iso_fortran_env, wp => real64
      use crest_data
      implicit none
      type(systemdata) :: env
      !type(options)    :: opt
      integer :: TMPCONF
      integer :: at(nat)
      integer :: nat
      real(wp),allocatable :: pop(:)
      real(wp),allocatable :: xyz(:,:,:),eread(:)

      allocate(xyz(3,nat,TMPCONF),eread(TMPCONF))

      select case(env%dftruntype)
        case( 1,2 )
         allocate(pop(TMPCONF))
         call etotprop(TMPCONF,env,pop,.true.)
         deallocate(pop)
         call rdpropens(TMPCONF,nat,xyz) !get updated geometries
         call wrpropens(TMPCONF,nat,xyz,at,eread)
         if(env%hardcutDFT)then
            call wrpropens_pop(env,TMPCONF,nat,xyz,at,eread,env%pthr)
         endif
       case( 3 )
         write(*,*)
         if(env%dftprog == 1)then !TM
           call aoforce_turbomole(env,TMPCONF)
         else if(env%dftprog == 2)then !ORCA
           !ORCA goes here
         endif
         call rdpropens(TMPCONF,nat,xyz) !get updated geometries
         call wrpropens(TMPCONF,nat,xyz,at,eread)
       case( 4 )
         write(*,*)
         if(env%dftprog == 1)then !TM
           call replaceGridsize(env,TMPCONF,3)
           call aoforce_turbomole(env,TMPCONF)
         else if(env%dftprog == 2)then !ORCA
           !ORCA goes here
         endif
         call autoir(TMPCONF,0,env)
         call rdpropens(TMPCONF,nat,xyz) !get updated geometries
         call wrpropens(TMPCONF,nat,xyz,at,eread)
       case default
          continue
      end select

      deallocate(eread,xyz)
 
!--- possible cleanup     
      if(env%pclean)then
       call cleanDFT(TMPCONF) !remove large files such as "mos" and "dh"
      endif

      return
end subroutine DFTprocessing
!------------------------------------------------------------------------------------------------------
! cut DFT populations HARD
subroutine cutDFTpop(env,pop,TMPCONF)!,eread)
      use iso_fortran_env, wp => real64
      use crest_data
      use iomod
      implicit none

      type(systemdata) :: env
      !type(options)    :: opt
      integer          :: TMPCONF
      real(wp)         :: pop(TMPCONF)
      !real(wp)         :: eread(TMPCONF)

      integer  :: maxs
      real(wp) :: popthr,popsum
      real(wp),allocatable :: tmppop(:)
      integer :: tmpmax,i
      integer,allocatable :: trac(:)

      allocate(tmppop(TMPCONF),trac(TMPCONF)) 

       maxs=max(1,env%hardcutnst)    !max nr. of structures to be considered. at least one.
       popthr=env%harcutpthr  !sum of populations that must be included.

      tmppop=pop
      do i=1,TMPCONF
         trac(i)=maxloc(tmppop(:),1)
         tmppop(trac(i))= -1.0_wp 
      enddo 

      popsum=0.0_wp
      tmpmax=0
      tmppop=pop
      pop=0.0_wp
      do i=1,TMPCONF
         pop(trac(i))=1.0_wp
         tmpmax=tmpmax+1
         if(tmpmax.ge.maxs)exit
         popsum=popsum+tmppop(trac(i))  
         if(popsum.ge.popthr)exit
      enddo 

      deallocate(tmppop)
      return
end subroutine cutDFTpop

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine dftTMwarning
      write(*,'(1x,a)')'! Warning: this feature requires a working TURBOMOLE infrastructure! !'
      write(*,'(1x,a)')'! Additionally the "cefine" tool is required.                        !'
      write(*,'(1x,a)')'! Tested with TURBOMOLE 7.3.1, older version might have troubles.    !'
      write(*,*)
end subroutine dftTMwarning
!-----------------------------------------------------------------------------------------------------
subroutine cefine_setup(env,TMPCONF)
      use iso_fortran_env, wp => real64
      use crest_data
      use iomod
      use strucrd, only: xyz2coord
      implicit none

      type(systemdata) :: env
      !type(options)    :: opt
      integer          :: TMPCONF

      integer :: i,j,k,l
      integer :: vz,io

      character(len=52) :: bar
      real(wp) :: percent

      character(len=512) :: thispath
      character(len=64) :: tmppath,val
      character(len=:),allocatable  :: cefine

      if(env%autothreads)then
        call ompautoset(env%threads,7,env%omp,env%MAXRUN,TMPCONF) !to enforce serial performance replace TMPCONF by 1
      endif
      !write(*,*) env%omp,env%MAXRUN
      !call ompprint_intern()
      !call ompprint()
      call TMparnodes(env%omp)
      !call getenv('PARNODES',val)
      !write(*,*)'Parnodes=',trim(val)

      call getcwd(thispath)

      write(*,'(1x,a)')'Setting up DFT calculations using the following command:'
      write(*,'(1x,a,1x,a)')'>',env%dftsetup
      k=0
      do i=1,TMPCONF
        write(tmppath,'(a,i0)')'TMPCONF',i
        k=k+1
        call chdir(trim(tmppath))

        call xyz2coord('struc.xyz','coord')

        cefine=trim(env%dftsetup)//' >cefine.out 2>>cefine.out'
        !call system(cefine)
        call execute_command_line(cefine, exitstat=io)


        if(env%dftmasscal)then
          call rdcontrol('control','control',env%dftmasspar,.false.)
        endif

        if(env%niceprint)then
          percent=float(k)/float(TMPCONF)*100
          call  progbar(percent,bar)
          call printprogbar(percent,bar)
        else
          write(6,'(1x,i0)',advance='no')k
          flush(6)
        endif
        call chdir(trim(thispath))
      enddo
      write(*,*)
      write(*,'(1x,a)')'done.'

      return
end subroutine cefine_setup

!---- quick ridft
subroutine ridft_turbomole(env,TMPCONF)
      use iso_fortran_env, wp => real64
      use crest_data
      implicit none
      type(systemdata) :: env
      !type(options)    :: opt
      integer          :: TMPCONF
      character(len=:),allocatable :: jobcall
      real(wp),allocatable :: pop(:)
!--- setting the threads for correct parallelization
      if(env%autothreads)then
        call ompautoset(env%threads,7,env%omp,env%MAXRUN,TMPCONF) !set global OMP/MKL variable for xtb jobs
      endif
!--- construct jobcall and run in parallel
      jobcall='ridft > ridft.out 2>>ridft.out'
      allocate(pop(TMPCONF), source=1.0_wp)
      write(*,'(1x,a)')'Running ridft with new settings ...'
      call prop_OMP_loop(env,TMPCONF,jobcall,pop)
      write(*,'(/,1x,a)')'done.'
      deallocate(pop)
      return
end subroutine ridft_turbomole

!---- aoforce calculation
subroutine aoforce_turbomole(env,TMPCONF)
      use iso_fortran_env, wp => real64
      use crest_data
      use iomod
      implicit none

      type(systemdata) :: env
      !type(options)    :: opt
      integer          :: TMPCONF
      
      real(wp),allocatable :: pop(:)
      real(wp) :: pthr
      integer :: i,j,k,l
      integer :: npop
      integer :: maxpop
      character(len=1024) :: jobcall
      
      pthr=env%pthr ! population threshold
 
   !--- thread handling
      if(env%autothreads)then
        call ompautoset(env%threads,7,env%omp,env%MAXRUN,TMPCONF) !to enforce serial performance replace TMPCONF with 1
      endif
      call TMparnodes(env%omp)

   !--- construct the aoforce call, incl. parallelization
      write(jobcall,'(a,1x,i0,1x,a)')'aoforce -smpcpus',env%omp,' >aoforce.out 2>>aoforce.out'

   !--- obtain populations
      allocate(pop(TMPCONF))
      call etotprop(TMPCONF,env,pop,.false.)
      if(env%hardcutDFT)then
         call cutDFTpop(env,pop,TMPCONF)
      endif

      maxpop=maxloc(pop(:),1)
      npop=0
      do i=1,TMPCONF
         if(pop(i).ge.pthr .or. i.eq.maxpop )npop=npop+1
      enddo

      call smallhead('Calculating DFT frequencies with TURBOMOLE (aoforce)')
      write(*,'(1x,a,i0,a,i0,a,i0,a)') &
      & 'For ',npop,' out of ',TMPCONF,' structures populations >',nint(pthr*100.0_wp), &
      & '% (-pthr) were obtained.'
      write(*,'(1x,a)')'Running aoforce for those structures and skipping the rest...'
      call prop_OMP_loop(env,TMPCONF,jobcall,pop)
      write(*,'(1x,a)')'done.'


      deallocate(pop)
      return
end subroutine aoforce_turbomole
!-----------------------------------------------------------------------------------------------------
! cut DFT populations HARD
subroutine replaceGridsize(env,TMPCONF,g)!,eread)
      use iso_fortran_env, wp => real64
      use crest_data
      use iomod
      use filemod
      implicit none

      type(systemdata) :: env
      !type(options)    :: opt
      type(filetype)   :: control
      integer          :: TMPCONF
      integer          :: g
      integer :: i,j,k,l
      character(len=:),allocatable :: tmppath
      character(len=64) :: dum
      character(len=64) :: gridsize
      integer,allocatable :: greps(:)
      logical :: ex

      write(gridsize,'(a,i0)')'   gridsize   m',g

      tmppath=''

      do i=1,TMPCONF
        write(dum,'(a,i0)')'TMPCONF',i
        tmppath=trim(dum)//'/control'
        inquire(file=tmppath,exist=ex)
        
        if(ex)then
           call control%open(tmppath)
           greps = control%findall('gridsize',k)
           !l = control%grephead('gridsize')
           do j=1,k
             l =greps(j)
             call control%replace(l,trim(gridsize))
           enddo
           call control%flush
           call control%close
        endif
      enddo
      if(allocated(greps))deallocate(greps)
      return
end subroutine replaceGridsize


!---- parallelization settings for Turbomole (SMD enviroment)
subroutine TMparnodes(omp)
      use iomod
      implicit none
      integer,intent(inout) :: omp
      integer :: io

      io = setenv('PARA_ARCH','SMP')
      io = setenv('PARNODES',omp)
      io = setenv('TM_PAR_FORK','on')

end subroutine TMparnodes
!------------------------------------------------------------------------------------------------------
subroutine dftORCAwarning
      write(*,'(1x,a)')'! Warning: this feature requires a working ORCA infrastructure! !'
      write(*,'(1x,a)')'! Tested with  ORCA 4.1, older version might have troubles.    !'
      write(*,*)
end subroutine dftORCAwarning 

