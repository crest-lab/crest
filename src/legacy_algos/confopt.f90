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
      integer :: vz, T,Tn
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
      call new_ompautoset(env,'auto',TMPCONF,T,Tn)

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
      integer :: TMPCONF
      character(len=1024) :: jobcall
      character(len=*) :: base
      character(len=:),allocatable :: bdir
      logical :: niceprint
      real(wp) :: percent
      integer :: vz,k,i,maxpop
      integer :: io
      character(len=512) :: tmppath

      if(niceprint)then
        call printprogbar(0.0_wp)
      endif
      k=0        ! count finished jobs

!$omp parallel &
!$omp shared( vz,jobcall,bdir,TMPCONF,percent,k,niceprint,maxpop )
!$omp single
      do i=1,TMPCONF
         vz=i
      !$omp task firstprivate( vz ) private( tmppath,io )
         call initsignal()
      !$omp critical
         write(tmppath,'(a,i0)')trim(base),vz
      !$omp end critical
         call command('cd '//trim(tmppath)//' && '//trim(jobcall), io)
      !$omp critical
        k=k+1
        if(niceprint)then
          percent=float(k)/float(TMPCONF)*100
          call printprogbar(percent)
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

!---------------------------------------------------------------------
! Parallel Optimization along trajectory
! Directories are created as needed.
!---------------------------------------------------------------------
subroutine MDopt_para_inplace(env,ensnam,multilev)
         use iso_c_binding
         use iso_fortran_env, only : wp => real64
         use crest_data
         use iomod
         use strucrd
      !$ use omp_lib
         implicit none

         type(systemdata) :: env

         character(len=*),intent(in)  :: ensnam
         integer,intent(in) :: multilev

         character(len=80)  :: fname,pipe
         character(len=256) :: btmp,ctmp
         character(len=512) :: thispath,optpath,tmppath
         character(len=40),allocatable :: origin(:)
         character(len=1024) :: jobcall,jobcall2

         integer :: nat
         real(wp),allocatable :: eread(:)
         real(wp),allocatable :: xyz(:,:,:)
         integer,allocatable  :: at(:)
         
         integer :: i,k,l
         integer :: ich,r,io
         integer :: nall
         integer :: vz,TMPCONF,T,Tn
         integer,allocatable  :: at0(:)
         logical :: xo

         logical :: niceprint   
         real(wp) :: percent

         logical :: verbose,l1
       
         nat = env%nat

!---- get current path
         call getcwd(thispath)

!--- some general settings
         fname='struc.xyz'

         l1=.true.

         verbose=.true.
         if(verbose)then
           pipe=' 2>/dev/null'
         else
           pipe=' > xtb.out 2>/dev/null'
         endif


!---- read the input ensemble
        call rdensembleparam(ensnam,nat,nall)
        if(nall.lt.1) return
        allocate(xyz(3,nat,nall),at(nat),eread(nall))
        if(.not.env%trackorigin)then
        call rdensemble(ensnam,nat,nall,at,xyz,eread)
        else
        allocate(origin(nall))
        call rdensemble_origin(ensnam,nat,nall,at,xyz,eread,origin)
        endif

!---- create meta-directory for the optimizations
         optpath='OPTIM'
         call rmrf(optpath)
         r = makedir(trim(optpath))

         !call copysub('.CHRG',trim(optpath))
         !call copysub('.UHF',trim(optpath))
         call env%wrtCHRG(trim(optpath))
         call copysub(env%fixfile,trim(optpath))
         !call copysub(env%constraints,trim(optpath))
         if(env%useqmdff)then
         call copysub('solvent',trim(optpath))
         endif
         if(env%gfnver=='--gff')then
            l = sylnk(trim(thispath)//'/'//'gfnff_topo',trim(optpath)//'/'//'gfnff_topo')
         endif

         call chdir(trim(optpath))
         call getcwd(optpath)


!---- the jobcall
      !create the system call (it is the same for every optimization)
      write(jobcall,'(a,1x,a,1x,a,'' --ceasefiles --opt '',a,1x,a,'' >xtb.out'')') &
      &    trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),trim(pipe)
      if(env%useqmdff)then
      write(jobcall,'(a,1x,a,1x,a,'' --ceasefiles --opt --qmdff '',a,1x,a,'' >xtb.out'')') &
      &    trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),trim(pipe)
      endif
      if(env%altopt)then
      write(jobcall,'(a,1x,a,1x,a,'' --ceasefiles --opt '',a,1x,a,'' >xtb.out'')') &
      &    trim(env%ProgName),trim(fname),trim(env%gfnver2),trim(env%solv),trim(pipe)
      env%reweight=.false.
      endif

      if(env%reweight)then
        write(jobcall2,'(a,1x,a,1x,a,'' --ceasefiles --sp '',a,1x,a,'' >xtb.out'')') &
      &    trim(env%ProgName),'xtbopt.xyz',trim(env%gfnver2),trim(env%solv),trim(pipe)
        jobcall = trim(jobcall)//' ; '//trim(jobcall2)
      endif

!---- parallel loop

      TMPCONF = nall
      niceprint = env%niceprint
      k=0

! setting the threads for correct parallelization
      call new_ompautoset(env,'auto',TMPCONF,T,Tn)

      write(*,'(1x,a,i0,a,a,a)')'Optimizing all ',nall,' structures from file "',trim(ensnam),'" ...'

!$omp parallel &
!$omp shared( vz,jobcall,optpath,tmppath,ctmp,TMPCONF,percent,k,niceprint) &
!$omp shared( nat,at,xyz,multilev,xo,btmp,env ) 
!$omp single
      do i=1,TMPCONF
         call initsignal()
         vz=i
         !$omp task firstprivate( vz ) private( tmppath,ctmp,io,ich,xo,btmp,at0 )
         call initsignal()
         !$omp critical
         tmppath=""
         write(tmppath,'(a,i0)')"TMPCONF",vz
         io = makedir(trim(tmppath))
         ich = vz + 1000
         !!$omp critical
         ctmp=trim(tmppath)//'/'//trim(fname)
         open(unit=ich,file=ctmp)
              call wrxyz(ich,nat,at,xyz(:,:,vz))
              if(multilev.ge.0)then
                  l1=.false.
                  write(ich,'(a)')'$end'
                  write(ich,'(a)')'$opt'
                  write(ich,'(2x,a,f14.4)')'hlow=',env%hlowopt
                  write(ich,'(2x,a,i0)')'microcycle=',nint(env%microopt)
                  write(ich,'(2x,a,f14.4)')'s6=',env%s6opt
              endif
              select case(multilev)
                case( 0 ) !--- from global variable
                  write(ich,'(2x,''optlevel='',i6)')nint(env%optlev)
                  write(ich,'(''$end'')')
                case( 1 ) !--- crude
                  write(ich,'(2x,''optlevel='',i6)')-3
                  write(ich,'(''$end'')')
                case( 2 ) !--- loose
                  write(ich,'(2x,''optlevel='',i6)')-1
                  write(ich,'(''$end'')')
                case( 3 ) !--- normal
                  write(ich,'(2x,''optlevel='',i6)')0
                  write(ich,'(''$end'')')
                case( 4 ) !--- tight
                  write(ich,'(2x,''optlevel='',i6)')1
                  write(ich,'(''$end'')')
                case( 5 ) !--- vtight 
                  write(ich,'(2x,''optlevel='',i6)')2
                  write(ich,'(''$end'')')
              end select
              call write_cts(ich,env%cts)
              call write_cts_biasext(ich,env%cts)
              if (env%crestver .EQ. crest_solv) then
                call write_cts_NCI(ich,env%cts)
                write(ich,'(''$end'')')
              end if

          close(ich)

          if(index(env%fixfile,'none selected').eq.0)then
          io = sylnk(trim(optpath)//'/'//env%fixfile,trim(tmppath)//'/'//env%fixfile)
          endif

          !$omp end critical

         call env%wrtCHRG(trim(tmppath))
         if(env%useqmdff)then
         call copysub('solvent',trim(tmppath))
         endif
         if(env%gfnver=='--gff')then
            io = sylnk(trim(optpath)//'/'//'gfnff_topo',trim(tmppath)//'/'//'gfnff_topo')
         endif

       !-- run the optimization
       call command('cd '//trim(tmppath)//' && '//trim(jobcall), io)  

       !$omp critical
        !-- write structure to ensemble
        ctmp=trim(tmppath)//'/'//'xtbopt.xyz'
        inquire(file=trim(ctmp),exist=xo)
        if(xo)then
           allocate(at0(nat))
           call rdxmol(trim(ctmp),nat,at0,xyz(:,:,vz),btmp)
           eread(vz) = grepenergy(btmp)
           deallocate(at0)
           if(env%reweight)then
             ctmp=trim(tmppath)//'/'//'xtb.out'
             call grepval(trim(ctmp),'total energy',l1,eread(vz))
           endif
        else   
           eread(vz) = 1.0d0 
        endif
        !-- clean up
        if(env%crestver .ne. crest_solv) call rmrf(trim(tmppath))
        k=k+1
        if(niceprint)then
          percent=float(k)/float(TMPCONF)*100
          call printprogbar(percent)
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
 
      if(.not.niceprint)then
        write(*,'(/,1x,a)') 'done.'
      else
        write(*,*)
      endif  

!---- write file and clear some memory
      if(env%trackorigin)then
          do i=1,nall
             origin(i) = '!'//trim(origin(i))
          enddo
          call wrensemble('opt.xyz',nat,nall,at,xyz,eread,origin)
          deallocate(origin,eread,at,xyz)
      else
          call wrensemble('opt.xyz',nat,nall,at,xyz,eread)
          deallocate(eread,at,xyz)
      endif


!---- go back to original directory
      call chdir(thispath)

end subroutine MDopt_para_inplace
