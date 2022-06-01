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

!==========================================================!
! subroutine confopt
! Optimize all stuctures of an ensemle within previously
! set up directories.
! The optimizations are done within an OMP loop.
! On Input: opt - main storage of program boolean data
!           sys - main storage of general program data
!           xyz - the coord file name in each directory
!           TMPCONF - total number of structures
!           confcross - boolean, is this the opt. after GC?
! On Output: none
!            structures are written to an ensemble file 
!            called 'opt.xyz'
!==========================================================!
subroutine confopt(env,xyz,TMPCONF,confcross)
      use iso_fortran_env, only : wp => real64
      use crest_data
      use iomod
      use strucrd
      implicit none
 
      type(systemdata) :: env

      character(len=*),intent(in)  :: xyz   !file base name
      integer,intent(in) :: TMPCONF  !number of structures to be optimized
      logical,intent(in) :: confcross ! used after confcross?

      integer :: i
      integer :: vz
      integer :: fileid

      character(len=20) :: pipe
      character(len=80) :: solv
      character(len=256) :: tmpname,oname         
      character(len=512) :: str,thispath,tmppath 
      character(len=1024):: jobcall,jobcall2       
      
      logical :: ex
      logical :: l1,l2,update

      logical :: niceprint

      real(wp),allocatable :: c0(:,:)
      integer,allocatable  :: at(:)
      real(wp) :: energy
      character(len=128) :: comment

! setting the threads for correct parallelization
      if(env%autothreads)then
        call ompautoset(env%threads,7,env%omp,env%MAXRUN,TMPCONF) !set global OMP/MKL variable for xtb jobs
      endif

      write(*,'(1x,''Starting optimization of generated structures'')')
      write(*,'(1x,i0,'' jobs to do.'')')TMPCONF

      solv=''
      pipe='2>/dev/null'
      call getcwd(thispath)
      update=.true.

      niceprint=env%niceprint

      !create the system call (it is the same for every optimization)
      write(jobcall,'(a,1x,a,1x,a,'' --opt '',a,1x,a,'' >xtb.out'')') &
      &    trim(env%ProgName),trim(xyz),trim(env%gfnver),trim(env%solv),trim(pipe)
      if(env%useqmdff)then
      write(jobcall,'(a,1x,a,1x,a,'' --opt --qmdff '',a,1x,a,'' >xtb.out'')') &
      &    trim(env%ProgName),trim(xyz),trim(env%gfnver),trim(env%solv),trim(pipe)
      endif
      if(env%altopt)then
      write(jobcall,'(a,1x,a,1x,a,'' --ceasefiles --opt '',a,1x,a,'' >xtb.out'')') &
      &    trim(env%ProgName),trim(xyz),trim(env%gfnver2),trim(env%solv),trim(pipe)
      env%reweight=.false.
      endif

 
      if(env%reweight)then
        write(jobcall2,'(a,1x,a,1x,a,'' --sp '',a,1x,a,'' >xtb.out'')') &
      &    trim(env%ProgName),'xtbopt.xyz',trim(env%gfnver2),trim(env%solv),trim(pipe)
        jobcall = trim(jobcall)//' ; '//trim(jobcall2)
      endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(TMPCONF < 1 ) then
         !write(*,*)'No new structures to be optimized. Exiting.'
         return
         goto 667
      endif
 
      !-- Do the optimizations (parallel system calls) 
      call opt_OMP_loop(TMPCONF,'TMPCONF',jobcall,niceprint)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      oname='opt.xyz'
      do i=1,20                                                  
         write(tmpname,'(a,''_'',i0,''.xyz'')')crefile,i
                 
         inquire(file=tmpname, exist=l2)
         if(l2) then
           oname=tmpname
           !exit
         endif
      enddo
      
      write(*,*) ''
      write(*,*) 'done.'
      write(str,'(''Now appending '',a,'' file with new structures'')')trim(oname)
      write(*,*) trim(str)

      allocate(c0(3,env%nat),at(env%nat))
      open(newunit=fileid,file=trim(oname))
      if(.not.confcross)then
       do vz=1,TMPCONF
         write(tmppath,'(''TMPCONF'',i0)')vz
         str=trim(tmppath)//'/'//'.xtboptok'
         inquire(file=trim(str),exist=ex)
         if(ex)then
         str=trim(tmppath)//'/'//'xtbopt.xyz'
         call rdxmol(trim(str),env%nat,at,c0,comment)
         if(env%reweight)then
           str=trim(tmppath)//'/'//'xtb.out'
           call grepval(trim(str),'total energy',l1,energy)
           write(comment,'(2x,f18.8)') energy
         endif
         call wrxyz(fileid,env%nat,at,c0,comment)
         endif
       enddo
      else
       do vz=1,TMPCONF
         write(tmppath,'(''TMPCONF'',i0)')vz
         str=trim(tmppath)//'/'//'.xtboptok'
         inquire(file=trim(str),exist=ex)
         if(ex)then
         str=trim(tmppath)//'/'//'xtbopt.xyz'
         call rdxmol(trim(str),env%nat,at,c0,comment)
         if(env%reweight)then
           str=trim(tmppath)//'/'//'xtb.out'
           call grepval(trim(str),'total energy',l1,energy)
           write(comment,'(2x,f18.8)') energy
         endif
         if(env%trackorigin)then
           comment = trim(comment)//' !gc'
         endif
         call wrxyz(fileid,env%nat,at,c0,comment)
         endif
       enddo
      endif
      close(fileid)
      deallocate(at,c0)
      
667   continue
      return
end subroutine confopt
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine clearconfopt
      use iomod
      implicit none
      call remove('.xtboptok')
      call remove('RUNNING')
      call remove('FINISHED')
      return
end subroutine clearconfopt


!=====================================================================!
! a general version of an omp-parallel optimization loop
! 
! On Input: TMPCONF   - number of dirs in total
!           base      - base name of the dirs, e.g., 'TMPCONF'
!           jobcall   - the system call to be executed in all dirs
!           niceprint - boolean for nicer printout
!=====================================================================!
subroutine opt_OMP_loop(TMPCONF,base,jobcall,niceprint)
      use iso_fortran_env, wp => real64
      use crest_data
      use iomod
      implicit none

      !type(systemdata) :: env
      !!type(options)    :: opt

      integer :: TMPCONF
      character(len=1024) :: jobcall
      character(len=*) :: base
      character(len=:),allocatable :: bdir

      logical :: niceprint

      character(len=52) :: bar
      real(wp) :: percent

      integer :: vz,k,i,maxpop
      integer :: io

      character(len=512) :: tmppath


      !niceprint=env%niceprint
      if(niceprint)then
        call printemptybar()
      endif
      k=0        ! count finished jobs

!$omp parallel &
!$omp shared( vz,jobcall,bdir,TMPCONF,percent,k,niceprint,bar,maxpop )
!$omp single
      do i=1,TMPCONF
         vz=i
      !$omp task firstprivate( vz ) private( tmppath,io )
         call initsignal()
         write(tmppath,'(a,i0)')trim(base),vz
         !write(str,'("cd ",a," && ",a)')trim(tmppath),trim(jobcall)
         !write(*,*) trim(str)
         call execute_command_line('cd '//trim(tmppath)//' && '//trim(jobcall), exitstat=io)
      !$omp critical
        k=k+1
        if(niceprint)then
          percent=float(k)/float(TMPCONF)*100
          call  progbar(percent,bar)
          call printprogbar(percent,bar)
        else
          if(gui)then
             call wrGUIpercent(k,TMPCONF,100)
          else
            write(6,'(1x,i0)',advance='no')k
            flush(6)
          endif
        endif
      !$omp end critical
      !$omp end task
      enddo
!$omp taskwait
!$omp end single
!$omp end parallel

      return
end subroutine opt_OMP_loop
