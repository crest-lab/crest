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

!--------------------------------------------------------------------------------------------
! A quick single point xtb calculation
!--------------------------------------------------------------------------------------------
subroutine xtbsp(env,xtblevel)
         use iso_fortran_env, only : wp => real64
         use iomod
         use crest_data
         implicit none
         type(systemdata) :: env
         !type(options)    :: opt
         integer,optional :: xtblevel
         character(len=80) :: fname,pipe,xtbflag
         character(len=512) :: jobcall
         integer :: io
!---- some options
         pipe=' > xtb.out 2>/dev/null'
         call remove('gfnff_topo')
         call remove('energy')
         if(.not.env%chargesfile)call remove('charges')
         call remove('xtbrestart')
!---- (OPTIONAL) select xtb level and set flag
         if(present(xtblevel))then
             select case(xtblevel)
               case( 0 )
                   xtbflag = '--gfn0'
               case( 1 )
                   xtbflag = '--gfn1'
               case( 2 )
                   xtbflag = '--gfn2'
               case( 3 )
                   xtbflag = '--gfnff'    
               case default
                   xtbflag = trim(env%gfnver)
             end select
         else
             xtbflag = trim(env%gfnver)
         endif
!---- setting threads
         if(env%autothreads)then
            call ompautoset(env%threads,7,env%omp,env%MAXRUN,1) !set the global OMP/MKL variables for the xtb jobs
         endif
!---- new plain coord file
         fname='tmpcoord'
         call copy('coord',fname)
         call clear_setblock(fname)
!---- jobcall
         write(jobcall,'(a,1x,a,1x,a,'' --sp '',a,1x,a,a)') &
         &     trim(env%ProgName),trim(fname),trim(xtbflag),trim(env%solv),trim(pipe)
         !call system(trim(jobcall))
         call execute_command_line(trim(jobcall), exitstat=io)
!---- cleanup
         call remove(fname)
         call remove('xtb.out')
         call remove('energy')
         if(.not.env%chargesfile)call remove('charges')
         call remove('xtbrestart')
         call remove('xtbtopo.mol')
         call remove('gfnff_topo')
end subroutine xtbsp

!--------------------------------------------------------------------------------------------
! A quick single point xtb calculation
!--------------------------------------------------------------------------------------------
subroutine xtbsp2(fname,env)
         use iso_fortran_env, only : wp => real64
         use iomod
         use crest_data
         implicit none
         character(len=*) :: fname
         type(systemdata) :: env
         !type(options)    :: opt
         character(len=512) :: jobcall
         character(*),parameter :: pipe=' > xtbcalc.out 2>/dev/null'
         integer :: io
         call remove('gfnff_topo')
         call remove('energy')
         if(.not.env%chargesfile)call remove('charges')
         call remove('xtbrestart')
!---- setting threads
         if(env%autothreads)then
            call ompautoset(env%threads,7,env%omp,env%MAXRUN,1) !set the global OMP/MKL variables for the xtb jobs
         endif
!---- jobcall
         write(jobcall,'(a,1x,a,1x,a,'' --sp --wbo '',a,1x,a,a)') &
         &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),trim(pipe)
         !call system(trim(jobcall))
         call execute_command_line(trim(jobcall), exitstat=io)
!---- cleanup
         call remove('xtbcalc.out')
         call remove('energy')
         if(.not.env%chargesfile)call remove('charges')
         call remove('xtbrestart')
         call remove('xtbtopo.mol')
         call remove('gfnff_topo')
end subroutine xtbsp2

!--------------------------------------------------------------------------------------------
! A quick xtb geometry optimization at the beginning of the program
!--------------------------------------------------------------------------------------------
subroutine xtbopt(env)
         use iso_fortran_env, only : wp => real64
         use iomod
         use crest_data
         use strucrd
         implicit none
         type(systemdata) :: env
         !type(options)    :: opt
         character(len=80) :: fname,pipe,solv
         character(len=512) :: jobcall
         logical :: fin
         character(len=256) :: atmp
         integer :: ich,iost,io,i
         type(coord) :: mol
         integer :: ntopo
         integer,allocatable :: topo(:)
         logical :: tchange = .false.

!---- small header
         write(*,*)
         call smallhead('xTB Geometry Optimization')
!---- some options
         pipe=' > xtb.out 2>/dev/null'
         call remove('gfnff_topo')
         if(.not.env%chargesfile)call remove('charges')
         call remove('grad')
         call remove('mos')
         call remove('xtbopt.log')
         call remove('xtbrestart')

!---- setting threads
         if(env%autothreads)then
            call ompautoset(env%threads,7,env%omp,env%MAXRUN,1) !set the global OMP/MKL variables for the xtb jobs
         endif
!---- new plain coord file
         fname='tmpcoord'
         !call copy('coord',fname)
         call wrc0(fname,env%ref%nat,env%ref%at,env%ref%xyz)
         call clear_setblock(fname)
!---- coord setup section
         open(newunit=ich,file=fname)
         do
              read(ich,'(a)',iostat=iost)atmp
              if(iost < 0)exit
              if(index(atmp,'$coord').ne.0) cycle
              if(index(atmp(1:1),'$').ne.0)then
                !write(ich,'(a)')'$end'
                exit
              endif
         enddo
         !add constraints (only if given, else the routine returns)
         call write_cts(ich,env%cts)
         call write_cts_biasext(ich,env%cts)
         write(ich,'(a)') '$end'
         close(ich)

!---- jobcall
         write(jobcall,'(a,1x,a,1x,a,'' --opt '',a,1x,a,a)') &
         &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),trim(pipe)
         call execute_command_line(trim(jobcall), exitstat=io)

         call minigrep('xtb.out','optimized geometry written to:',fin)
         if(.not.fin)then
             write(*,*)
             write(*,*) ' Initial geometry optimization failed!'
             write(*,*) ' Please check your input.'
             error stop
         endif
         write(*,*) 'Geometry successfully optimized.'
!---- if necessary, check if the topology has changed!
         call mol%open('xtbopt.coord')
         if(allocated(env%ref%topo))then
          ntopo = mol%nat*(mol%nat+1)/2
          allocate(topo(ntopo))
          call quicktopo(mol%nat,mol%at,mol%xyz,ntopo,topo)
          do i=1,ntopo
            if(topo(i).ne.env%ref%topo(i)) tchange=.true.
          enddo
          if(tchange)then
           write(*,'(1x,a)') 'WARNING! Change in topology detected!'   
           !--- either update the topology (see option B below)
           if(.not.env%reftopo)then
             env%ref%topo = topo
           !--- or abort the run
           else
           write(*,'(1x,a)') 'The topology change was seen in the initial geometry optimization.'  
           write(*,'(1x,a,a,a)') 'This could be an artifact of the chosen theory level (',trim(env%gfnver),').'
           write(*,'(1x,a)') 'You can check the optimization trajectory in the "xtbopt.log" file.'
           write(*,'(1x,a)') 'Try either of these options:'
           write(*,'(/,4x,a)') 'A) Pre-optimize your input seperately with xtb and use the optimized'
           write(*,'(4x,a)') '   structure as input for CREST. (Only recommended if structure is intact)'
           write(*,'(/,4x,a)') 'B) Restart the same CREST call as before, but ignore the topology change'
           write(*,'(4x,a)') '   by using the "--noreftopo" keyword. (May produce artifacts)'
           write(*,'(/,4x,a)') 'C) Fix the initial input geometry by introducing bond length constraints'
           write(*,'(4x,a)') '   or by using a method with fixed topology (GFN-FF).'
           write(*,*)
           error stop 'abnormal termination of crest'
           endif
          endif
         endif
!---- update reference with optimized geometry
         env%ref%nat = mol%nat
         env%ref%at = mol%at
         env%ref%xyz = mol%xyz
         call mol%deallocate()
         call rename('xtbopt.coord','coord')

!---- cleanup
         call remove(fname)
         call remove('xtb.out')
         call remove('energy')
         if(.not.env%chargesfile)call remove('charges')
         call remove('grad')
         call remove('mos')
         call remove('xtbopt.log')
         call remove('xtbrestart')
         call remove('gfnff_topo')
end subroutine xtbopt


!--------------------------------------------------------------------------------------------
! A single METADYN run (resp. its setup)
!--------------------------------------------------------------------------------------------
subroutine MetaMD(env,nr,mdtime,fac,expo,dumplist)
         use iso_c_binding
         use iso_fortran_env, only : wp => real64
         use iomod
         use crest_data
         use strucrd, only: wrc0,rdensembleparam,rdensemble
         implicit none
         type(systemdata) :: env
         !type(options)    :: opt

         real(wp) :: mdtime

         real(wp)  :: fac
         real(wp)  :: expo
         integer   :: dumplist,dumplist2
         integer   :: nr

         character(len=20)  :: fname
         character(len=80)  :: solv,pipe
         character(len=256) :: basename,dirname
         character(len=512) :: tmppath,thispath,str

         integer :: r
         logical :: verbose,parallel,ex

         real(wp),allocatable :: confs(:,:,:),eread(:)
         integer,allocatable ::  ats(:)
         integer :: nall2,iz1,iz2
         integer :: i,j,k,l
 
         associate( nat => env%nat, hmass => env%hmass, mdtemp => env%mdtemp,    &
         & mdstep => env%mdstep, shake => env%shake, mddumpxyz => env%mddumpxyz, &
         & mdskip => env%mdskip, mddump => env%mddump )

!---- some settings
         basename='METADYN'  !base name of the directories
         fname='coord'

         dumplist2=dumplist

!----

         call getcwd(thispath)              !current dir= thispath
         !call getname_dir(trim(basename),dirname)
         write(dirname,'(a,i0)') trim(basename),nr
         tmppath=trim(dirname)
         call rmrf(tmppath)         !clear old directory
         r = makedir(trim(tmppath)) !create new directory

         call copysub('coord',trim(tmppath))
         !call copysub('.CHRG',trim(tmppath))
         !call copysub('.UHF',trim(tmppath))
         call env%wrtCHRG(trim(tmppath))
         call copysub(env%fixfile,trim(tmppath))
         call copysub(env%constraints,trim(tmppath))
         if(env%useqmdff)then
            call copysub('solvent',trim(tmppath))
         endif
         if(env%gfnver=='--gff')then
!            l = sylnk(trim(thispath)//'/'//'gfnff_topo',trim(tmppath)//'/'//'gfnff_topo')
         endif

         call chdir(trim(tmppath))  !switch to working directory

!---- do stuff here
         call setMDrun2(fname,hmass,mdtime,mdtemp,mdstep,shake,mddumpxyz, &
         &             mdskip,mddump,-1,env%cts)
         call setMetadyn2(fname,fac,expo,dumplist2)  ! Metadynamic settings

!----
         call chdir(thispath)  !go back to orgiginal directory  
         end associate
end subroutine MetaMD

!--------------------------------------------------------------------------------------------
! Run several METADYN in parallel, OMP VERSION!
!--------------------------------------------------------------------------------------------
subroutine MetaMD_para_OMP(env)
      use iso_fortran_env, only : wp => real64
      use iomod
      use crest_data
      implicit none

      type(systemdata) :: env
      !type(options)    :: opt

      integer :: i,j,k,l
      integer :: ACTIVE,DONE,TOTAL,vz
      real(wp) :: time
      character(len=512) :: thispath,tmppath
      character(len=512) :: jobcall
      character(len=80)  :: fname,pipe,solv,atmp,btmp
      integer :: dum,io
      logical :: ex

      if(env%autothreads)then
         dum = env%nmetadyn
         call ompautoset(env%threads,7,env%omp,env%MAXRUN,dum) !set the global OMP/MKL variables for the xtb jobs
      endif

      time = env%mdtime !for some reason this is necessary

      call getcwd(thispath)

      fname='coord'
      pipe=' > xtb.out 2>/dev/null'

      write(jobcall,'(a,1x,a,1x,a,1x,''--md'',1x,a,1x,a,a)') &
      &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),pipe
    !--- slightly different jobcall for QMDFF usage
      if(env%useqmdff)then
      write(jobcall,'(a,1x,a,1x,a,1x,''--md --qmdff'',1x,a,1x,a,a)') &
      &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),pipe
      endif

!---- Small Header
      write(*,*)
      write(*,'(5x,''========================================'')')
      write(*,'(5x,''|         Meta-MD (MTD) Sampling       |'')')
      write(*,'(5x,''========================================'')')
      write(*,*)

!---- set up directories
      do i=1,env%nmetadyn
           call MetaMD(env,i,time,env%metadfac(i),env%metadexp(i), &
       &               env%metadlist(i))
      enddo

!$omp parallel &
!$omp shared( vz,jobcall,time,env )
!$omp single
      do i=1,env%nmetadyn
         call initsignal()
         vz=i
       !$omp task firstprivate( vz ) private( tmppath,io,ex )
         call initsignal()
       !$omp critical
             write(*,'(a,i4,a)') 'Starting Meta-MD',vz,' with the settings:'
             write(*,'(''     MD time /ps        :'',f8.1)')time
             write(*,'(''     dt /fs             :'',f8.1)')env%mdstep
             write(*,'(''     dumpstep(trj) /fs  :'',i8)')env%mddumpxyz
             write(*,'(''     dumpstep(Vbias)/ps :'',f8.1)')float(env%mddump)/1000d0
             write(*,'(''     Vbias factor k /Eh :'',f8.4)')env%metadfac(vz)
             write(*,'(''     Vbias exp α /bohr⁻²:'',f8.2)')env%metadexp(vz)
       !$omp end critical
         write(tmppath,'(a,i0)')'METADYN',vz
         call execute_command_line('cd '//trim(tmppath)//' && '//trim(jobcall), exitstat=io)
         inquire(file=trim(tmppath)//'/'//'xtb.trj',exist=ex)
         if(.not.ex .or. io.ne.0)then
         write(*,'(a,i0,a)')'*Warning: Meta-MTD ',vz,' seemingly failed (no xtb.trj)*'
         call system('cp -r '//trim(tmppath)//' FAILED-MTD')
         else
         write(*,'(a,i0,a)')'*Meta-MTD ',vz,' finished*'
         endif
      !$omp end task
      enddo
!$omp taskwait
!$omp end single
!$omp end parallel

     if(env%trackorigin)then
        call set_trj_origins('METADYN','mtd')
     endif
     call checkname_xyz(crefile,atmp,btmp)
     call collect_trj_skipfirst('METADYN',trim(atmp))  !collect all 'xtb.trj' from the METADYN directories,skip the first point on .trj

     if(.not.env%keepModef)then
       call cleanMTD
     endif

end subroutine MetaMD_para_OMP
 
!---------------------------------------------------------------------
! Parallel Optimization along trajectory
!---------------------------------------------------------------------
subroutine MDopt_para(env,ensnam,multilev)
         use iso_c_binding
         use iso_fortran_env, only : wp => real64
         use crest_data
         use iomod
         use strucrd, only: rdensembleparam,rdensemble,wrxyz
         implicit none
         type(systemdata) :: env
         !type(options)    :: opt
         character(len=*),intent(in)  :: ensnam
         integer,intent(in) :: multilev
         character(len=80)  :: solv,fname,pipe
         character(len=256) :: btmp,ctmp
         character(len=512) :: str,dg
         character(len=512) :: thispath,optpath
         character(len=40),allocatable :: origin(:)
         real(wp) :: lev

         real(wp),allocatable :: eread(:)
         real(wp),allocatable :: xyz(:,:,:)
         integer,allocatable  :: at(:)
         
         integer :: i,j,k,l
         integer :: ich,r
         integer :: iz1,iz2,nall

         logical :: verbose,l1
       
         associate( nat => env%nat )

!------- IMPORTANT: wrapper if in-place mode is active ------------!
! The in-place version does the same as this routine, but is
! less heavy on the disk-space...
         if(env%inplaceMode)then
            call MDopt_para_inplace(env,ensnam,multilev)
            return
         endif
!------------------------------------------------------------------!

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

!---- create directory for the optimizations
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

!---- set up directories
      write(6,'(1x,a,a,a)',advance='no')'writing TMPCONF* Dirs from file "',trim(ensnam),'" ...'
      do i=1,nall
         write(ctmp,'(''TMPCONF'',i0)')i
         r = makedir(trim(ctmp))
         call chdir(ctmp)
              open(newunit=ich,file=fname)
              call wrxyz(ich,nat,at,xyz(:,:,i))
              if(multilev.ge.0)then
                  l1=.false.
                  write(ich,'(a)')'$end'
                  write(ich,'(a)')'$opt'
                  write(ich,'(1x,a,f14.4)')'hlow=',env%hlowopt
                  write(ich,'(1x,a,i0)')'microcycle=',nint(env%microopt)
                  write(ich,'(1x,a,f14.4)')'s6=',env%s6opt
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
              close(ich)
         call chdir(optpath)

         !call copysub('.CHRG',trim(ctmp))
         !call copysub('.UHF',trim(ctmp))
         call env%wrtCHRG(trim(ctmp))
         call copysub(env%fixfile,trim(ctmp))
         if(env%useqmdff)then
         call copysub('solvent',trim(ctmp))
         endif
         if(env%gfnver=='--gff')then
            l = sylnk(trim(optpath)//'/'//'gfnff_topo',trim(ctmp)//'/'//'gfnff_topo')
         endif
      enddo
      write(*,'(1x,a)') 'done.'

!---- clear some memory
      deallocate(eread,at,xyz)


!---- perform the optimizations using the confopt routine
      call confopt(env,fname,nall,.false.)


      if(env%trackorigin)then
        call addorigin2('opt.xyz',origin,nall)
        deallocate(origin)
      endif

!---- go back to original directory
      call chdir(thispath)

      end associate
end subroutine MDopt_para

!--------------------------------------------------------------------------------------------
! Different handling of optimizations in multi-step-filtering
!--------------------------------------------------------------------------------------------
subroutine multilevel_opt(env,modus)
     use iso_fortran_env, only : wp => real64
     use iomod
     use crest_data
     implicit none

     type(systemdata) :: env
     !type(options)    :: opt
 

     character(len=128) :: inpnam,outnam
     character(len=512) :: thispath,filename
     character(len=:),allocatable :: olev
     character(len=:),allocatable :: headder
     real(wp)  :: time
     real(wp)  :: newthr
     real(wp)  :: ewinbackup,rthrbackup
     real(wp),allocatable :: backupthr(:)
 
     integer :: modus,nremain

     call getcwd(thispath)

     allocate(backupthr(8))
     backupthr=env%thresholds    !default thresholds
     ewinbackup=env%ewin
     rthrbackup=env%rthr

!     if(multilevel)then
      select case( modus )
!--- first optimization with "maxopt"
      case( 1 )
        call smallhead('1. crude pre-optimization')
        call checkname_xyz(crefile,inpnam,outnam)
        call MDopt_para(env,trim(inpnam),modus)
        filename=trim(thispath)//'/'//trim(outnam) 
        call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
    !---using cregen to sort the optimized structures
        call checkname_xyz(crefile,inpnam,outnam)
        !--- get other thresholds
        newthr=aint(ewinbackup*( 15.0d0 / 6.0d0 ))
        env%ewin=newthr
        env%rthr=rthrbackup*3.0d0 !larger RTHR
        call sort_and_check(env,trim(filename))
      !--- only if modus=1 the input file is overwritten with the sorted file
        call remove(inpnam)
        call rename(outnam,inpnam)
    !-----
        call rmoptim()
        write(*,*)
!--- then vloose optimization
      case( 2,22 )
        call checkname_xyz(crefile,inpnam,outnam)
        if(modus==2)then  
        call smallhead('Ensemble optimization with loose thresholds')
        call MDopt_para(env,trim(inpnam),modus)
        else
        call smallhead('Ensemble optimization with crude thresholds')
        call MDopt_para(env,trim(inpnam),1)
        endif
        filename=trim(thispath)//'/'//trim(outnam)
        call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
        !--- get other thresholds
        newthr=aint(ewinbackup*(10.0d0/6.0d0))
        env%ewin=newthr   !
        call sort_and_check(env,trim(filename))
        call rmoptim()
        write(*,*)
!--- lastly vtight optimization
      case( 3 )
        call smallhead('3. optimization with very tight thresholds')
        call checkname_xyz(crefile,inpnam,outnam)
        call MDopt_para(env,trim(inpnam),0) !optlev is set from the module variable
        filename=trim(thispath)//'/'//trim(outnam)
        call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
        call sort_and_check(env,trim(filename))
        call rmoptim()
        write(*,*)
!--- crude optimization
      case( 4 )
        call smallhead('1. crude pre-optimization')
        call checkname_xyz(crefile,inpnam,outnam)
        call MDopt_para(env,trim(inpnam),1)
        filename=trim(thispath)//'/'//trim(outnam)
        call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
    !---using cregen to sort the optimized structures
        call checkname_xyz(crefile,inpnam,outnam)
        !--- get other thresholds
        newthr=aint(ewinbackup* 2.0d0)
        env%ewin=newthr
        env%rthr=rthrbackup*2.0d0 !larger RTHR
        call sort_and_check(env,trim(filename))
      !--- only if modus=1 the input file is overwritten with the sorted file
        call remove(inpnam)
        call rename(outnam,inpnam)
    !-----
        call rmoptim()
        write(*,*)
!--- tight optimization
      case( 5:6 )
        if(modus.eq.5)then
        call smallhead('2. optimization with tight thresholds')
        else if(modus.eq.6)then
        call smallhead('Ensemble optimization with tight thresholds')
        endif
        call checkname_xyz(crefile,inpnam,outnam)
        call MDopt_para(env,trim(inpnam),4) !optlev is set from the module variable
        filename=trim(thispath)//'/'//trim(outnam)
        call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
        call sort_and_check(env,trim(filename))
        call rmoptim()
        write(*,*)
!--- default optimization
      case default
        olev=optlevflag(env%optlev)
        headder='Ensemble optimization with '//olev//' thresholds'  
        call smallhead(headder)
        !call smallhead('Ensemble optimization')
        call checkname_xyz(crefile,inpnam,outnam)
        call MDopt_para(env,trim(inpnam),0) !optlev is set from the module variable
        filename=trim(thispath)//'/'//trim(outnam)
        call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
!        if (env%crestver .ne. crest_solv) then
        if(.not. env%QCG) then
            call sort_and_check(env,trim(filename))
            call rmoptim()
        end if
        write(*,*)
      end select

     env%thresholds=backupthr
     env%ewin=ewinbackup
     env%rthr=rthrbackup
     deallocate(backupthr)
 
end subroutine multilevel_opt
!--------------------------------------------------------------------------------------------
! Construct MD $set-block version 2 (updated)
!--------------------------------------------------------------------------------------------
subroutine  setMDrun2(fname,hmass,mdtime,mdtemp,mdstep,shake,mddumpxyz, &
 &                   mdskip,mddump,nvt,cts)
    use iso_fortran_env, only : wp => real64
    use crest_data
    implicit none
    type(constra) :: cts
    character(len=*) :: fname
    real(wp) :: hmass
    real(wp) :: mdtime
    real(wp) :: mdtemp
    real(wp) :: mdstep
    integer  :: shake
    integer  :: mddumpxyz
    integer  :: mdskip
    integer  :: mddump
    integer  :: nvt
    integer :: ich,iost
    character(len=256):: atmp
    logical :: ex

    !--- coord setup section
    open(newunit=ich,file=fname)
    !--- don't modify coords:
    do
      read(ich,'(a)',iostat=iost)atmp
      if(iost < 0)exit
      if(index(atmp,'$coord').ne.0) cycle
      if(index(atmp(1:1),'$').ne.0)exit
    enddo
    !--- then write the MD settings 
    write(ich,'(a)') '$md'
         if(hmass>0.0_wp)then      !set H-atom mass
            write(ich,'(2x,a,i0)') 'hmass=',nint(hmass)
         endif
         if(mdtime>0.0_wp)then     !set MD simulation time in ps
            write(ich,'(2x,a,f10.2)') 'time=',mdtime
         endif
         if(mdtemp>0.0_wp)then     !set MD Temperature
            write(ich,'(2x,a,f10.2)') 'temp=',mdtemp
         endif
         if(mdstep.gt.0)then     !set MD timestep in fs
            write(ich,'(2x,a,f10.2)') 'step=',mdstep
         endif
         if(shake>=0)then     !set MD shake mode
            write(ich,'(2x,a,i0)')'shake=',shake
         endif
         !the order of setting mddump and mddumpxyz is important!!!
         if(mddumpxyz>0)then    ! Frequency of structure dump into xtb.trj in fs
            write(ich,'(2x,a,i0)')'dump=',mddumpxyz
         endif
         if(mdskip>=0)then    !skipping optimizations in -mdopt
             write(ich,'(2x,a,i0)')'skip=',mdskip
         endif
         if(nvt>=0)then      ! NVT ensemble
           if(nvt==0)then
             write(ich,'(2x,a)')'nvt=false'
           else
             write(ich,'(2x,a)')'nvt=true'
           endif
         endif
         !the order of setting mddump and mddumpxyz is important!!!
         if(mddump>0)then     ! Vbias dump for metadyn in fs                       
            write(ich,'(a)')'$set'
            write(ich,'(2x,a,2x,i0)')'mddump',mddump
         endif
    !--- add constraints (only if given, else the routines return)
    call write_cts(ich,cts)
    call write_cts_NCI(ich,cts)
    if(cts%cbonds_md)then
    call write_cts_CBONDS(ich,cts)   
    endif
    if(cts%dispscal_md)then
    call write_cts_DISP(ich,cts)
    endif
    call write_cts_rcontrol(ich,cts)
    write(ich,'(a)') '$end'
    close(ich)
    return
end subroutine setMDrun2

!--------------------------------------------------------------------------------------------
subroutine setMetadyn2(fname,fac,expo,dumplist)
         use iso_fortran_env, only : wp => real64
         implicit none
         character(len=*) :: fname
         real(wp)  :: fac
         real(wp)  :: expo
         integer   :: dumplist
         character(len=256) :: atmp
         character(len=256),allocatable :: mdyn(:)
         integer :: i,j
         integer :: ich,ich2,iost,lz

         allocate(mdyn(10))
         mdyn(1:10)=''
         lz=1

         open(newunit=ich,file=fname)
         open(newunit=ich2,file='tmpcoordfile')
         do
              read(ich,'(a)',iostat=iost)atmp
              if(iost < 0)exit
              if(index(atmp,'$metadyn').ne.0)then
                 
                 do
                   read(ich,'(a)',iostat=iost)atmp
                   if(iost < 0)exit
                   if(index(atmp,'$').ne.0)exit
                   mdyn(lz)=adjustl(atmp)
                   lz=lz+1
                 enddo
              endif
              write(ich2,'(a)')trim(atmp)
         enddo

         write(ich2,'(a)') '$metadyn'
         do i=1,lz
            if(index(mdyn(i),'save').ne.0) cycle
            if(index(mdyn(i),'kpush').ne.0)cycle
            if(index(mdyn(i),'alp').ne.0)  cycle
            if(trim(mdyn(i)).eq.'')cycle
            write(ich2,'(2x,a)') trim(mdyn(i))
         enddo
         write(ich2,'(2x,a,i0)') 'save=',dumplist
         write(atmp,'(f12.6)') fac
         write(ich2,'(2x,a,a)') 'kpush=',adjustl(trim(atmp))
         write(atmp,'(f12.6)') expo
         write(ich2,'(2x,a,a)') 'alp=',adjustl(trim(atmp))
         write(ich2,'(a)') '$end'
         close(ich2)
         close(ich)
         deallocate(mdyn)
         call rename('tmpcoordfile',fname)
         return
end subroutine setMetadyn2

!---------------------------------------------------------------------
! Modified verison of the GC (confscript.v.2)
!---------------------------------------------------------------------
subroutine cross2(env)
      use iso_fortran_env, only : wp => real64
      use crest_data
      use iomod
      implicit none
      type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
      !type(options) :: opt       ! MAIN STORAGE OF BOOLEAN SETTINGS

      logical  :: better
      real(wp) :: ewin,ewinbackup
      integer  :: i,imax,tmpconf,nremain
      logical  :: ex
      character(len=80)  :: str
      character(len=128) :: inpnam,outnam
      character(len=512) :: thispath,tmppath

      real(wp),allocatable :: backupthr(:)

      allocate(backupthr(8))
      backupthr=env%thresholds
      ewinbackup=env%ewin

      call getcwd(thispath)
      tmppath='OPTIM'
      
      write(*,*)
      write(*,'(5x,''========================================'')')
      write(*,'(5x,''|        Structure Crossing (GC)       |'')')
      write(*,'(5x,''========================================'')')

      imax=min(nint(env%mdtime*50.0d0),5000) ! long setting
      if(env%setgcmax)then
        imax=nint(env%gcmax)
      endif

      if(env%quick)then
        imax=nint(float(imax)*0.5d0)
      endif

      
      do i=1,1
          call confcross(env,imax,tmpconf)
          if(tmpconf.lt.1) then
              return
              exit
          endif
          if(env%gcmultiopt)then !---printout
            call smallhead('GC: loose pre-optimization')
          endif

          call chdir('OPTIM')
          call confopt(env,'confcross.xyz',tmpconf,.true.)
              
          if(env%gcmultiopt)then
          !--- stay in dir OPTIM
            call rmrfw('TMPCONF')
            call checkname_xyz(crefile,inpnam,outnam)
            call rename('opt.xyz',trim(inpnam))
            !env%thresholds(1)=aint(backupthr(1)*(10.0d0/6.0d0))
            env%ewin=aint(ewinbackup*(10.0d0/6.0d0))
            call confg_chk3(env)
            call checkname_xyz(crefile,inpnam,outnam)
            call remaining_in(inpnam,env%ewin,nremain)
            env%thresholds=backupthr
            env%ewin=ewinbackup
            call smallhead('GC: optimization with tight thresholds')
            if(env%iterativeV2)then
               call MDopt_para(env,inpnam,4) !<--- creates another dir OPTIM within the current dir OPTIM (fixed optlev tight)
            else
               call MDopt_para(env,inpnam,0) !<--- creates another dir OPTIM within the current dir OPTIM
            endif
            call rename('OPTIM'//'/'//'opt.xyz','opt.xyz') !move the optimized GC ensemble into current direcotry
          endif

          !--- exit dir OPTIM and get file
          call chdir(trim(thispath))  
          call checkname_xyz(crefile,inpnam,outnam)
          call appendto('OPTIM'//'/'//'opt.xyz',trim(inpnam))
          call rmrfw('TMPCONF')
          call rmrf('OPTIM')
          call rmrf('zmat*.zmat')
      enddo
      deallocate(backupthr)
end subroutine cross2


!---------------------------------------------------------------------
! run confg check (confscript.v.2)
!---------------------------------------------------------------------
subroutine confg_chk3(env)
      use crest_data
      implicit none
      type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
      !type(options) :: opt       ! MAIN STORAGE OF BOOLEAN SETTINGS

      associate(aut => env%autothreads, th => env%threads, omp => env%omp, &
      &         MX => env%MAXRUN )

      if(aut)then
         call ompautoset(th,4,omp,MX,0) !mode=4 --> Program intern Threads max
         if(.not.env%newcregen)then
             call cregen2(env)
         else
             !--- Special handling qcg, no RMSD, because a CMA transformed structure would cause wrong wall pot.
             if(env%crestver .eq. crest_solv) then
                call newcregen(env,6)
             else
                call newcregen(env,0)
             end if
         endif
         call ompautoset(th,5,omp,MX,0) !mode=5 --> Program intern Threads min
      else
         if(.not.env%newcregen)then
             call cregen2(env)
         else
             if(env%crestver .eq. crest_solv) then
                call newcregen(env,6)
             else
                call newcregen(env,0)
            end if
         endif
      endif
      end associate
end subroutine confg_chk3

!---------------------------------------------------------------------
! collect xtb.trj
!---------------------------------------------------------------------
subroutine collect_trj(base,whichfile)
      use iomod
      implicit none
      character(len=*) :: base,whichfile
      character(len=256) :: str,dir
      integer :: i
      logical :: ex
      i=1
      do
        write(dir,'(a,i0)')trim(base),i
        ex = directory_exist(trim(dir))
        if(.not.ex)then
          exit
        else
          write(str,'(a,a,''xtb.trj'')')trim(dir),'/'
          call appendto(trim(str),whichfile)
          i=i+1
        endif
      enddo
      return
end subroutine collect_trj

!---------------------------------------------------------------------
! collect xtb.trj, but skip the first point on the trj
!---------------------------------------------------------------------
subroutine collect_trj_skipfirst(base,whichfile)
      use iomod
      implicit none
      character(len=*) :: base,whichfile
      character(len=256) :: str,dir
      integer :: i
      logical :: ex
      i=1
      do
        write(dir,'(a,i0)')trim(base),i
        ex = directory_exist(trim(dir))
        if(.not.ex)then
          exit
        else
          write(str,'(a,a,''xtb.trj'')')trim(dir),'/'
          call TRJappendto_skipfirst(trim(str),whichfile)
          i=i+1
        endif
      enddo
      return
end subroutine collect_trj_skipfirst

!------------------------------------------------------------------------
! writes the Input coord to the first position of crest_rotamers_*.xyz
!------------------------------------------------------------------------
subroutine append_INPUT_to(fnam,nat,newtag)
      use iomod
      use crest_data
      use strucrd, only: coord2xyz
      implicit none
      integer :: nat
      character(len=*)  :: fnam
      character(len=80) :: tmpnam,inam,onam    
      character(len=*) :: newtag

      tmpnam='ensemble.xyz.tmp'
      call checkname_xyz(crefile,inam,onam)

      call coord2xyz(fnam,trim(tmpnam))
      call addoriginXYZ(trim(tmpnam),trim(newtag))

      call appendto(inam,tmpnam)
      call remove(inam)
      call rename(tmpnam,inam)
      return
end subroutine append_INPUT_to

!-------------------------------------------------------------------------
! print the number of remaining files in an ensemble file for a given energy window
!-------------------------------------------------------------------------
subroutine remaining_in(filename,ewin,nall)
      use strucrd, only: rdensembleparam,rdensemble
      implicit none
      integer :: iz1,iz2,nall,remain
      real*8 :: ewin
      character(len=*) :: filename
      integer :: i,j,k,nat

      open(newunit=k,file=trim(filename))
      read(k,*)nat
      close(k)

      call rdensembleparam(trim(filename),nat,nall)

      if(nall.lt.1)then
        write(*,*) 'No conformer was left. Something must be something seriously wrong.'
        write(*,*) 'Terminating the run.'
        error stop
      endif

      write(*,'(1x,i0,'' structures remain within '',f8.2,'' kcal/mol window'')') &
      &        nall,ewin
      return
end subroutine remaining_in

!-------------------------------------------------------------------------
! print the number of remaining files in an ensemble file for a given energy window
!-------------------------------------------------------------------------
subroutine sort_and_check(env,filename)
      use iso_fortran_env, only : wp => real64
      use crest_data
      use iomod
      use strucrd, only: rdensembleparam,rdensemble
      implicit none

      type(systemdata) :: env
      !type(options)    :: opt

      integer :: iz1,iz2,nall,remain
      real(wp) :: ewin,thresh(8)
      character(len=*) :: filename
      integer :: i,j,k,nallin,nallout,nallthr
      real(wp) :: nthr,increase
     
      character(len=80) :: inpnam,outnam

      real(wp),allocatable :: dummythr(:)

      increase=1.5d0  ! factor 1.5 increase 
      nthr=0.05d0     ! 5%

      ewin=env%ewin

      call rdensembleparam(trim(filename),env%nat,nallin) !--initial number of structures
      call checkname_xyz(crefile,inpnam,outnam)
      call confg_chk3(env)
      call remaining_in(outnam,ewin,nallout) !--- remaining number of structures
      
      if(.not.env%entropic)then !don't do this for entropy mode
      nallthr=nint(float(nallin)*nthr)
      if(nallout.lt.nallthr)then
       write(*,'(1x,''This is less than '',i0,''% of the initial '',i0,'' structures.'')') &
       & nint(nthr*100.0d0),nallin
       write(*,'(1x,''Increasing energy window to include more...'')')
       call remove(outnam)                    !--- delete this file since it contains to few structures
       env%ewin=ewin*increase
       call confg_chk3(env)
       call remaining_in(outnam,ewin*increase,nallout) !--- remaining number of structures
      endif
      endif

      return    
end subroutine sort_and_check

!--------------------------------------------------------------------------------------------
! collect all the backup ensembles from the different runs into a new file.
!--------------------------------------------------------------------------------------------
subroutine collectcre(env)
       use iomod
       use crest_data
       implicit none

       type(systemdata) :: env

       character(len=80) :: atmp,btmp
       character(len=80) :: crename
       integer :: i
       logical :: ex
     
     !---- remove crest_romtamer and crest_conformer files
       call rmcres()

     !---- get the crename
       call checkname_xyz(crefile,crename,btmp)
      
       do i=1,env%Maxrestart
          call checkname_xyz('.cre',atmp,btmp)
          inquire(file=atmp,exist=ex)
          if(.not.ex)exit
          call appendto(trim(atmp),trim(crename))
          call remove(atmp)
       enddo
       call rmrfw('.cre_')
       return
end subroutine collectcre

!--------------------------------------------------------------------------------------------
! check lower found
!--------------------------------------------------------------------------------------------
subroutine elowcheck(lower,env)
       use iso_fortran_env, only : wp => real64
       use iomod
       use crest_data
       use strucrd, only: xyz2coord
       implicit none
       type(systemdata) :: env
       !type(options) :: opt
       real(wp) :: ediff,ethr,ewin
       real(wp),parameter :: autokcal = 627.509541d0
       logical :: lower

!--- some defaults
     ethr = env%ethr  !ETHR threshold that is also used in CREGEN (=0.05 kcal/mol)
     if(env%entropic)then !except for entropy mode
         ethr = 0.1_wp
     endif
     ewin=env%ewin
     lower=.false.

         ediff=(env%eprivious-env%elowest)*autokcal
         if(ediff.ge.ethr)then
           write(*,'(1x,a)') '...............................................'
           write(*,'(1x,a)') 'A new lower conformer was found!'
           write(*,'(1x,a,f10.5,a,f10.5,a)') 'Improved by ', ediff/autokcal,' Eh or ',ediff, 'kcal/mol'
           write(*,'(1x,a)') '...............................................'
           env%eprivious=env%elowest
           lower=.true.
           if(.not.env%allowrestart)then
               lower=.false.
               return
           endif
           if(env%entropic)then
             write(*,'(1x,a)')'Restart based on new lowest structure.' 
             call clean_V2i
             call rmrfw('crest_entropy_rotamer_')
             return
           endif    
           !if(env%crestver==22)then ! "--v4"
           !  write(*,'(1x,a)') 'Energy window (ewin) will be increased to include the old ensemble'
           !  env%ewin = env%ewin + ediff
           !  return
           !endif
         !---- clean the dir
           call clean_V2i
         !---- save the new best conformer
           call XYZappendto('crest_best.xyz','.history.xyz')
           call xyz2coord('crest_best.xyz','coord')                     !new reference coord to start the MTDs with
         else
           lower=.false.
         endif

        return
end subroutine elowcheck

!-----------------------------------------------------------------------
! diatomic "error" catcher ---> we will have no conformers
!-----------------------------------------------------------------------
subroutine catchdiatomic(env)
    use iso_fortran_env, wp => real64
    use crest_data
    use strucrd
    use iomod
    implicit none
    !type(options) :: opt
    type(systemdata) :: env
    type(coord) :: strc
    integer :: ich
    character(len=1024) :: jobcall

    call strc%open('coord')

    open(file=conformerfile,newunit=ich)
    strc%xyz = strc%xyz*bohr !to ang
    call wrxyz(ich,strc%nat,strc%at,strc%xyz)
    close(ich)
    !create the system call (it is the same for every optimization)
    write(jobcall,'(a,1x,a,1x,a,'' --opt '',a,1x,a,'' --ceasefiles  >xtb.out'')') &
   &    trim(env%ProgName),conformerfile,trim(env%gfnver),trim(env%solv),' 2>/dev/null'
    call execute_command_line(trim(jobcall), exitstat=ich)
    call copy('xtbopt.xyz',conformerfile)
    call copy(conformerfile,'crest_rotamers.xyz')
    call copy(conformerfile,'crest_best.xyz')

    call strc%deallocate
    return
end subroutine catchdiatomic

!-------------------------------------------------------------------------------------!
! entropy ensemble file copy routine
!-------------------------------------------------------------------------------------!
subroutine emtdcopy(env,iter,stopiter,broken)
    use iso_fortran_env, wp => real64
    use crest_data
    use iomod
    use strucrd
    implicit none
    !type(options) :: opt
    type(systemdata) :: env
    integer :: iter,iter2
    logical :: stopiter
    logical,intent(inout) :: broken
    integer :: i,nall
    real(wp) :: nallfrac
    real(wp) :: sdiff
    logical :: ex,conv1,conv2
    character(len=80) :: atmp,btmp
    character(len=80) :: crename
    character(len=22),parameter :: sfile='crest_entropy_rotamer_'
    character(len=11),parameter :: bfile='crest_smtd_'
    character(len=:),allocatable :: filname
    real(wp) :: T,S,Cp
    integer :: nt


    stopiter = .false.
    broken = .false.
    T = 298.15d0
    !--- determine temperature dependence
    if(.not.allocated(env%thermo%temps))then
      call env%thermo%get_temps()
    endif
    nt=env%thermo%ntemps

    if(.not.allocated(env%emtd%soft))then
      allocate(env%emtd%soft(nt), source=0.0d0)
    endif
    if(.not.allocated(env%emtd%cpoft))then
      allocate(env%emtd%cpoft(nt), source=0.0d0)
    endif
    if(.not.allocated(env%emtd%hoft))then
      allocate(env%emtd%hoft(nt), source=0.0d0)
    endif

    if(env%crestver==22)then
      filname=trim(bfile)
    else
      filname=trim(sfile)
    endif

    if(iter == 0)then
       call checkname_xyz(crefile,atmp,crename)
       write(btmp,'(a,i0,a)') filname,iter,'.xyz'
       call rename(trim(atmp),trim(btmp))
       call rmrfw('crest_rotamers_')
       call rmrfw('.data')
       env%emtd%sapproxlast = env%emtd%sapprox
       call rdensembleparam(conformerfile,i,nall)
       env%emtd%nconflast=nall
       call entropic(env,.false.,.true.,.true.,trim(btmp),T,S,Cp)
       call writesdata(env,nall,iter)
       return
    endif

    if(iter >= 1)then
       call checkname_xyz(crefile,atmp,crename)
       iter2 = iter - 1
       write(btmp,'(a,i0,a)') filname,iter2,'.xyz' !-- file from last iter.
       inquire(file=trim(btmp),exist=ex)
       if(ex)then
        call appendto(trim(btmp),trim(atmp)) !-- cat together
        call newcregen(env,2)
        sdiff = (env%emtd%sapproxlast/env%emtd%sapprox)
        sdiff = 1.0d0 - sdiff
        sdiff = abs(sdiff)     !can dS be negative?
        if(sdiff < 0.0d0)then
            broken =.true.
        else
            broken = .false.
            env%emtd%sapproxlast = env%emtd%sapprox
        endif
       endif
    endif

    if(.not.broken)then
    call checkname_xyz(crefile,crename,btmp)
    write(btmp,'(a,i0,a)') filname,iter,'.xyz'
    call rename(trim(crename),trim(btmp))
    endif

    call rdensembleparam(conformerfile,i,nall)
    if(nall > 50000)then !safety fallback for extremely large ensembles (e.g. C18)
      env%emtd%confthr=0.1d0
    endif
    nallfrac = 1.0d0 - (float(env%emtd%nconflast)/float(nall))
    conv1 = (sdiff < env%emtd%sconvthr) !.and.(sdiff > 0.0d0)
    conv2 = (nallfrac < env%emtd%confthr).and.(nallfrac >= 0.0d0)
    if(nallfrac < 0.0d0)then !if we for some reason got less conformers in this iteration
        broken = .true.
    endif
    write(*,*)
    write(*,'(1x,a,i0,a)') 'Static MTD simulation iteration ',iter,' finished'
    write(*,'(1x,a,i0,a,f5.2,a,a)') 'A total of ',nall,' conformers (+',  &
    &    nallfrac*100.0d0, '%) in ',conformerfile
    write(*,'(1x,a,f12.6,a,f5.2,a)') 'S(approx.) =',env%emtd%sapprox, &
    &     ' cal/molK (dS =',sdiff*100.0d0,'%)'
    write(*,'(1x,a,l)') 'Convergence w.r.t. conformers: ',conv2
    write(*,'(1x,a,l)') 'Convergence w.r.t. entropy   : ',conv1
    write(*,*)

    if(.not.broken)then !-- NORMAL case
      if( conv1 .and. conv2 )then
         stopiter = .true.
         if( iter < env%emtd%iter)then
             call copy(trim(btmp),trim(crename))
         endif
      endif
      if(env%entropic .and. env%crestver.ne.22)then
      call entropic(env,.false.,.true.,.true.,trim(btmp),T,S,Cp)
      call writesdata(env,nall,iter)
      endif
      env%emtd%nconflast = nall
      env%emtd%sapproxlast = env%emtd%sapprox
    else                !-- ROLLBACK  
      !--- NOTE: instead of exiting one could also try to re-do the iteration,
      !          but this is commented out. Relevant comments are marked by "!*"
      write(*,'(1x,a)') 'Warning: Nconf DECREASED, which is unphysical.'  
      !write(*,'(1x,a)') 'Trying to re-do last iteration ...'   !*
      write(*,'(1x,a)') 'Rewind ensemble and exit iterations ...'    !*
      call checkname_xyz(crefile,crename,btmp)
      call remove(crename)
      write(btmp,'(a,i0,a)') filname,iter-1,'.xyz'
      call copy(trim(btmp),trim(crename))
      call xyz2coord(trim(crename),'coord')  !get the right coord again
      call newcregen(env,2) !get the right crest_conformers again
      iter=iter-1  !*
      stopiter = .true.
    endif

    if(iter < env%emtd%iter .and. .not.stopiter)then
    call rmrfw('crest_rotamers_')
    endif
    return
contains
subroutine writesdata(env,nall,inum)
    implicit none
    !type(options) :: opt
    type(systemdata) :: env
    integer :: nall
    integer :: inum
    character(len=64) :: adum
    integer :: ich,nt,i
    nt=env%thermo%ntemps
    write(adum,'(a,i0)')'Sdata',inum
    open(newunit=ich,file=trim(adum))
    write(ich,*) nall
    do i=1,nt
     write(ich,'(3F18.10)')env%emtd%soft(i), &
     &   env%emtd%Cpoft(i),env%emtd%hoft(i)
    enddo
    close(ich)
    return
end subroutine writesdata   
end subroutine emtdcopy

subroutine emtdcheckempty(env,empty,nbias)
    use crest_data
    use strucrd
    use iomod
    implicit none
    !type(options) :: opt
    type(systemdata) :: env
    logical :: empty
    integer :: nbias
    integer :: i,j,k,nall
    character(len=128) :: atmp,btmp
    call checkname_xyz(crefile,atmp,btmp)
    call rdensembleparam(trim(atmp),i,nall)
    empty = .false.

    if(nall<1)then
        call remove(trim(atmp))
        if(env%emtd%maxfallback>1)then
        write(*,'(1x,a)') 'Empty ensemble, retry iteration with less bias ...'
        else
        write(*,'(1x,a)') 'Empty ensemble, skipping ahead ...'
        endif
        call remove('crest_clustered.xyz')
        call remove('crest_bias.xyz')
        nbias = nbias - 1
        nbias = max(0,nbias)
        empty = .true.
    endif
    return
end subroutine emtdcheckempty

!========================================================================!
! Sample additional OH orientations
!========================================================================!
subroutine XHorient(env,infile)
    use iso_fortran_env, only: wp=>real64
    use crest_data
    use strucrd
    use zdata
    use iomod
    implicit none
    type(systemdata) :: env
    character(len=*) :: infile
    integer :: nat,nall
    character(len=:),allocatable :: newfile
    newfile='oh_ensemble.xyz'
    !--- shouldn't cost much
      call ohflip_ensemble(env,infile,env%maxflip)
      call rdensembleparam(newfile,nat,nall)
    !--- only proceed if there are potential new structures
      if(nall>0)then
        call smallhead('Additional orientation sampling')
        call MDopt_para(env,newfile,0)
      !---- printout and copy
        call rename('OPTIM'//'/'//'opt.xyz',trim(newfile))
        call rmrf('OPTIM')
      else
        call remove(newfile)  
      endif
    return
end subroutine XHorient
