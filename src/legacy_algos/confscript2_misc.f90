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
subroutine xtbsp_legacy(env,xtblevel)
  use crest_parameters
  use iomod
  use utilities
  use crest_data
  implicit none
  type(systemdata) :: env
  integer,optional :: xtblevel
  character(len=80) :: fname,xtbflag
  character(len=512) :: jobcall
  integer :: io
  character(*),parameter :: pipe = ' >xtb.out 2>/dev/null'
!---- some options
  !pipe=' > xtb.out 2>/dev/null'
  call remove('gfnff_topo')
  call remove('energy')
  call remove('wbo')
  if (.not.env%chargesfile) call remove('charges')
!---- (OPTIONAL) select xtb level and set flag
  if (present(xtblevel)) then
    select case (xtblevel)
    case (0)
      xtbflag = '--gfn0'
    case (1)
      xtbflag = '--gfn1'
    case (2)
      xtbflag = '--gfn2'
    case (3)
      xtbflag = '--gfnff'
    case default
      xtbflag = trim(env%gfnver)
    end select
  else
    xtbflag = trim(env%gfnver)
  end if
!---- setting threads
  if (env%autothreads) then
    call ompautoset(env%threads,7,env%omp,env%MAXRUN,1) !set the global OMP/MKL variables for the xtb jobs
  end if
!---- new plain coord file
  fname = 'tmpcoord'
  call copy('coord',fname)

!---- jobcall
  jobcall = ""
  jobcall = trim(jobcall)//trim(env%ProgName)
  jobcall = trim(jobcall)//" "//trim(fname)//" --sp"
  jobcall = trim(jobcall)//" "//trim(xtbflag)
  jobcall = trim(jobcall)//" "//trim(env%solv)
  jobcall = trim(jobcall)//pipe
  call command(trim(jobcall), io)
!---- cleanup
  call remove(fname)
  call remove('xtb.out')
  call remove('energy')
  if (.not.env%chargesfile) call remove('charges')
  call remove('xtbrestart')
  call remove('xtbtopo.mol')
  call remove('gfnff_topo')
end subroutine xtbsp_legacy

!--------------------------------------------------------------------------------------------
! A quick single point xtb calculation
!--------------------------------------------------------------------------------------------
subroutine xtbsp2_legacy(fname,env)
  use crest_parameters
  use iomod
  use crest_data
  implicit none
  character(len=*) :: fname
  type(systemdata) :: env
  character(len=512) :: jobcall
  character(*),parameter :: pipe = ' > xtbcalc.out 2>/dev/null'
  character(len=4) :: chrgstr
  integer :: io
  call remove('gfnff_topo')
  call remove('energy')
  if (.not.env%chargesfile) call remove('charges')
  call remove('xtbrestart')
!---- setting threads
  if (env%autothreads) then
    call ompautoset(env%threads,7,env%omp,env%MAXRUN,1) !set the global OMP/MKL variables for the xtb jobs
  end if
!---- jobcall
  jobcall = ""
  jobcall = trim(jobcall)//trim(env%ProgName)
  jobcall = trim(jobcall)//" "//trim(fname)//" --sp --wbo"
  jobcall = trim(jobcall)//" "//trim(env%gfnver)
  jobcall = trim(jobcall)//" "//trim(env%solv)
  if (env%chrg /= 0) then
    write (chrgstr,'(i0)') env%chrg
    jobcall = trim(jobcall)//" --chrg "//trim(chrgstr)
  end if
  if (env%uhf /= 0) then
    write (chrgstr,'(i0)') env%uhf
    jobcall = trim(jobcall)//" --uhf "//trim(chrgstr)
  end if
  jobcall = trim(jobcall)//pipe
  call command(trim(jobcall), io)
!---- cleanup
  call remove('xtbcalc.out')
  call remove('energy')
  if (.not.env%chargesfile) call remove('charges')
  call remove('xtbrestart')
  call remove('xtbtopo.mol')
  call remove('gfnff_topo')
end subroutine xtbsp2_legacy

!--------------------------------------------------------------------------------------------
! A quick xtb geometry optimization at the beginning of the program
!--------------------------------------------------------------------------------------------
subroutine xtbopt_legacy(env)
  use crest_parameters
  use iomod
  use crest_data
  use strucrd
  implicit none
  type(systemdata) :: env
  character(len=80) :: fname,pipe
  character(len=512) :: jobcall
  logical :: fin
  character(len=256) :: atmp
  character(len=4) :: chrgstr
  integer :: ich,iost,io,i
  type(coord) :: mol
  integer :: ntopo
  integer,allocatable :: topo(:)
  logical :: tchange = .false.
  logical :: ex 

!---- small header
  write (*,*)
  call smallhead('xTB Geometry Optimization')
!---- some options
  pipe = ' > xtb.out 2>/dev/null'
  call remove('gfnff_topo')
  if (.not.env%chargesfile) call remove('charges')
  call remove('grad')
  call remove('mos')
  call remove('xtbopt.log')
  call remove('xtbrestart')

!---- setting threads
  if (env%autothreads) then
    call ompautoset(env%threads,7,env%omp,env%MAXRUN,1) !set the global OMP/MKL variables for the xtb jobs
  end if
!---- new plain coord file
  fname = 'tmpcoord'
  call wrc0(fname,env%ref%nat,env%ref%at,env%ref%xyz)

!---- coord setup section
  open (newunit=ich,file=fname)
  do
    read (ich,'(a)',iostat=iost) atmp
    if (iost < 0) exit
    if (index(atmp,'$coord') .ne. 0) cycle
    if (index(atmp(1:1),'$') .ne. 0) then
      !write(ich,'(a)')'$end'
      exit
    end if
  end do
  !add constraints (only if given, else the routine returns)
  call write_cts(ich,env%cts)
  call write_cts_biasext(ich,env%cts)
  write (ich,'(a)') '$end'
  close (ich)

!>---- jobcall
  jobcall = ""
  jobcall = trim(jobcall)//trim(env%ProgName)
  jobcall = trim(jobcall)//" "//trim(fname)//' --opt'
  jobcall = trim(jobcall)//" "//trim(env%gfnver)
  jobcall = trim(jobcall)//" "//trim(env%solv)
  if (env%chrg /= 0) then
    jobcall = trim(jobcall)//" --chrg "//to_str(env%chrg)
  end if
  if (env%uhf /= 0) then
    jobcall = trim(jobcall)//" --uhf "//to_str(env%uhf)
  end if
  jobcall = trim(jobcall)//pipe
  call command(trim(jobcall), io)

  call minigrep('xtb.out','optimized geometry written to:',fin)
  inquire(file='xtbopt.coord',exist=ex)
  if(ex)then
    call mol%open('xtbopt.coord')
  endif

!>--- process the optimization status
  call trialOPT_warning(env,mol,fin)
  call mol%deallocate()
  call rename('xtbopt.coord','coord')

!>---- cleanup
  call remove(fname)
  call remove('xtb.out')
  call remove('energy')
  if (.not.env%chargesfile) call remove('charges')
  call remove('grad')
  call remove('mos')
  call remove('xtbopt.log')
  call remove('xtbrestart')
  call remove('gfnff_topo')
end subroutine xtbopt_legacy

!------------------------------------------------------------------------------
! A single METADYN run (resp. its setup)
!------------------------------------------------------------------------------
subroutine MetaMD(env,nr,mdtime,fac,expo,dumplist)
  use iso_c_binding
  use crest_parameters
  use iomod
  use crest_data
  use strucrd,only:wrc0,rdensembleparam,rdensemble
  use utilities
  implicit none
  type(systemdata) :: env

  real(wp) :: mdtime

  real(wp)  :: fac
  real(wp)  :: expo
  integer   :: dumplist,dumplist2
  integer   :: nr

  character(len=20)  :: fname
  character(len=256) :: basename,dirname
  character(len=512) :: tmppath,thispath
  integer :: r

  associate (nat => env%nat,hmass => env%hmass,mdtemp => env%mdtemp,    &
  & mdstep => env%mdstep,shake => env%shake,mddumpxyz => env%mddumpxyz, &
  & mdskip => env%mdskip,mddump => env%mddump)

!---- some settings
    basename = 'METADYN'  !base name of the directories
    fname = 'coord'

    dumplist2 = dumplist

!----

    call getcwd(thispath)              !current dir= thispath
    !call getname_dir(trim(basename),dirname)
    write (dirname,'(a,i0)') trim(basename),nr
    tmppath = trim(dirname)
    call rmrf(tmppath)         !clear old directory
    r = makedir(trim(tmppath)) !create new directory

    call copysub('coord',trim(tmppath))
    !call copysub('.CHRG',trim(tmppath))
    !call copysub('.UHF',trim(tmppath))
    call env%wrtCHRG(trim(tmppath))
    call copysub(env%fixfile,trim(tmppath))
    call copysub(env%constraints,trim(tmppath))
    if (env%useqmdff) then
      call copysub('solvent',trim(tmppath))
    end if
    if (env%gfnver == '--gff') then
!            l = sylnk(trim(thispath)//'/'//'gfnff_topo',trim(tmppath)//'/'//'gfnff_topo')
    end if

    call chdir(trim(tmppath))  !switch to working directory

!---- do stuff here
    call setMDrun2(fname,hmass,mdtime,mdtemp,mdstep,shake,mddumpxyz, &
    &             mdskip,mddump,-1,env%cts)
    call setMetadyn2(fname,fac,expo,dumplist2)  ! Metadynamic settings

!----
    call chdir(thispath)  !go back to orgiginal directory
  end associate
end subroutine MetaMD

!-----------------------------------------------------------------------------
! Run several METADYN in parallel, OMP VERSION!
!-----------------------------------------------------------------------------
subroutine MetaMD_para_OMP(env)
  use crest_parameters
  use iomod
  use crest_data
  use utilities
  implicit none

  type(systemdata) :: env

  integer :: i,vz
  real(wp) :: time
  character(len=512) :: thispath,tmppath
  character(len=512) :: jobcall
  character(len=80)  :: fname,pipe,atmp,btmp
  integer :: dum,io
  logical :: ex

  if (env%autothreads) then
    dum = env%nmetadyn
    call ompautoset(env%threads,7,env%omp,env%MAXRUN,dum) !set the global OMP/MKL variables for the xtb jobs
  end if

  time = env%mdtime !for some reason this is necessary

  call getcwd(thispath)

  fname = 'coord'
  pipe = ' > xtb.out 2>/dev/null'

  write (jobcall,'(a,1x,a,1x,a,1x,''--md'',1x,a,1x,a,a)') &
  &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),pipe
  !--- slightly different jobcall for QMDFF usage
  if (env%useqmdff) then
    write (jobcall,'(a,1x,a,1x,a,1x,''--md --qmdff'',1x,a,1x,a,a)') &
    &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),pipe
  end if

!---- Small Header
  write (*,*)
  write (*,'(5x,''========================================'')')
  write (*,'(5x,''|         Meta-MD (MTD) Sampling       |'')')
  write (*,'(5x,''========================================'')')
  write (*,*)

!---- set up directories
  do i = 1,env%nmetadyn
    call MetaMD(env,i,time,env%metadfac(i),env%metadexp(i), &
&               env%metadlist(i))
  end do

!$omp parallel &
!$omp shared( vz,jobcall,time,env )
!$omp single
  do i = 1,env%nmetadyn
    call initsignal()
    vz = i
    !$omp task firstprivate( vz ) private( tmppath,io,ex )
    call initsignal()
    !$omp critical
    write (*,'(a,i4,a)') 'Starting Meta-MD',vz,' with the settings:'
    write (*,'(''     MD time /ps        :'',f8.1)') time
    write (*,'(''     dt /fs             :'',f8.1)') env%mdstep
    write (*,'(''     dumpstep(trj) /fs  :'',i8)') env%mddumpxyz
    write (*,'(''     dumpstep(Vbias)/ps :'',f8.1)') float(env%mddump)/1000d0
    write (*,'(''     Vbias factor k /Eh :'',f8.4)') env%metadfac(vz)
    write (*,'(''     Vbias exp α /bohr⁻²:'',f8.2)') env%metadexp(vz)
    !$omp end critical
    write (tmppath,'(a,i0)') 'METADYN',vz
    call command('cd '//trim(tmppath)//' && '//trim(jobcall), io)
    inquire (file=trim(tmppath)//'/'//'xtb.trj',exist=ex)
    if (.not.ex.or.io .ne. 0) then
      write (*,'(a,i0,a)') '*Warning: Meta-MTD ',vz,' seemingly failed (no xtb.trj)*'
      call system('cp -r '//trim(tmppath)//' FAILED-MTD')
    else
      write (*,'(a,i0,a)') '*Meta-MTD ',vz,' finished*'
    end if
    !$omp end task
  end do
!$omp taskwait
!$omp end single
!$omp end parallel

  if (env%trackorigin) then
    call set_trj_origins('METADYN','mtd')
  end if
  call checkname_xyz(crefile,atmp,btmp)
  call collect_trj_skipfirst('METADYN',trim(atmp))  !collect all 'xtb.trj' from the METADYN directories,skip the first point on .trj

  if (.not.env%keepModef) then
    call cleanMTD
  end if

end subroutine MetaMD_para_OMP

!---------------------------------------------------------------------
! Parallel Optimization along trajectory
!---------------------------------------------------------------------
subroutine MDopt_para(env,ensnam,multilev)
  use iso_c_binding
  use crest_parameters
  use crest_data
  use iomod
  use strucrd,only:rdensembleparam,rdensemble,wrxyz
  use utilities
  implicit none
  type(systemdata) :: env
  character(len=*),intent(in)  :: ensnam
  integer,intent(in) :: multilev
  character(len=80)  :: fname,pipe
  character(len=256) :: ctmp
  character(len=512) :: thispath,optpath
  character(len=40),allocatable :: origin(:)

  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)

  integer :: i,l
  integer :: ich,r
  integer :: nall

  logical :: verbose,l1

  associate (nat => env%nat)

!------- IMPORTANT: wrapper if in-place mode is active ------------!
! The in-place version does the same as this routine, but is
! less heavy on the disk-space...
    if (env%inplaceMode) then
      call MDopt_para_inplace(env,ensnam,multilev)
      return
    end if
!------------------------------------------------------------------!

!---- get current path
    call getcwd(thispath)

!--- some general settings
    fname = 'struc.xyz'
    l1 = .true.
    verbose = .true.
    if (verbose) then
      pipe = ' 2>/dev/null'
    else
      pipe = ' > xtb.out 2>/dev/null'
    end if

!---- read the input ensemble
    call rdensembleparam(ensnam,nat,nall)
    if (nall .lt. 1) return
    allocate (xyz(3,nat,nall),at(nat),eread(nall))
    if (.not.env%trackorigin) then
      call rdensemble(ensnam,nat,nall,at,xyz,eread)
    else
      allocate (origin(nall))
      call rdensemble_origin(ensnam,nat,nall,at,xyz,eread,origin)
    end if

!---- create directory for the optimizations
    optpath = 'OPTIM'
    call rmrf(optpath)
    r = makedir(trim(optpath))

    !call copysub('.CHRG',trim(optpath))
    !call copysub('.UHF',trim(optpath))
    call env%wrtCHRG(trim(optpath))
    call copysub(env%fixfile,trim(optpath))
    !call copysub(env%constraints,trim(optpath))
    if (env%useqmdff) then
      call copysub('solvent',trim(optpath))
    end if
    if (env%gfnver == '--gff') then
      l = sylnk(trim(thispath)//'/'//'gfnff_topo',trim(optpath)//'/'//'gfnff_topo')
    end if

    call chdir(trim(optpath))
    call getcwd(optpath)

!---- set up directories
    write (6,'(1x,a,a,a)',advance='no') 'writing TMPCONF* Dirs from file "',trim(ensnam),'" ...'
    do i = 1,nall
      write (ctmp,'(''TMPCONF'',i0)') i
      r = makedir(trim(ctmp))
      call chdir(ctmp)
      open (newunit=ich,file=fname)
      call wrxyz(ich,nat,at,xyz(:,:,i))
      if (multilev .ge. 0) then
        l1 = .false.
        write (ich,'(a)') '$end'
        write (ich,'(a)') '$opt'
        write (ich,'(1x,a,f14.4)') 'hlow=',env%hlowopt
        write (ich,'(1x,a,i0)') 'microcycle=',nint(env%microopt)
        write (ich,'(1x,a,f14.4)') 's6=',env%s6opt
      end if
      select case (multilev)
      case (0) !--- from global variable
        write (ich,'(2x,''optlevel='',i6)') nint(env%optlev)
        write (ich,'(''$end'')')
      case (1) !--- crude
        write (ich,'(2x,''optlevel='',i6)')-3
        write (ich,'(''$end'')')
      case (2) !--- loose
        write (ich,'(2x,''optlevel='',i6)')-1
        write (ich,'(''$end'')')
      case (3) !--- normal
        write (ich,'(2x,''optlevel='',i6)') 0
        write (ich,'(''$end'')')
      case (4) !--- tight
        write (ich,'(2x,''optlevel='',i6)') 1
        write (ich,'(''$end'')')
      case (5) !--- vtight
        write (ich,'(2x,''optlevel='',i6)') 2
        write (ich,'(''$end'')')
      end select
      call write_cts(ich,env%cts)
      call write_cts_biasext(ich,env%cts)
      close (ich)
      call chdir(optpath)

      !call copysub('.CHRG',trim(ctmp))
      !call copysub('.UHF',trim(ctmp))
      call env%wrtCHRG(trim(ctmp))
      call copysub(env%fixfile,trim(ctmp))
      if (env%useqmdff) then
        call copysub('solvent',trim(ctmp))
      end if
      if (env%gfnver == '--gff') then
        l = sylnk(trim(optpath)//'/'//'gfnff_topo',trim(ctmp)//'/'//'gfnff_topo')
      end if
    end do
    write (*,'(1x,a)') 'done.'

!---- clear some memory
    deallocate (eread,at,xyz)

!---- perform the optimizations using the confopt routine
    call confopt(env,fname,nall,.false.)

    if (env%trackorigin) then
      call addorigin2('opt.xyz',origin,nall)
      deallocate (origin)
    end if

!---- go back to original directory
    call chdir(thispath)

  end associate
end subroutine MDopt_para

!--------------------------------------------------------------------------------------------
! Different handling of optimizations in multi-step-filtering
!--------------------------------------------------------------------------------------------
subroutine multilevel_opt(env,modus)
  use crest_parameters
  use iomod
  use crest_data
  use utilities
  implicit none

  type(systemdata) :: env

  character(len=128) :: inpnam,outnam
  character(len=512) :: thispath,filename
  character(len=:),allocatable :: olev
  character(len=:),allocatable :: headder
  real(wp)  :: newthr
  real(wp)  :: ewinbackup,rthrbackup
  real(wp),allocatable :: backupthr(:)
  integer :: modus

  call getcwd(thispath)

  allocate (backupthr(8))
  backupthr = env%thresholds    !default thresholds
  ewinbackup = env%ewin
  rthrbackup = env%rthr

!     if(multilevel)then
  select case (modus)
!--- first optimization with "maxopt"
  case (1)
    call smallhead('1. crude pre-optimization')
    call checkname_xyz(crefile,inpnam,outnam)
    call MDopt_para(env,trim(inpnam),modus)
    filename = trim(thispath)//'/'//trim(outnam)
    call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
    !---using cregen to sort the optimized structures
    call checkname_xyz(crefile,inpnam,outnam)
    !--- get other thresholds
    newthr = aint(ewinbackup*(15.0d0/6.0d0))
    env%ewin = newthr
    env%rthr = rthrbackup*3.0d0 !larger RTHR
    call sort_and_check(env,trim(filename))
    !--- only if modus=1 the input file is overwritten with the sorted file
    call remove(inpnam)
    call rename(outnam,inpnam)
    !-----
    call rmoptim()
    write (*,*)
!--- then vloose optimization
  case (2,22)
    call checkname_xyz(crefile,inpnam,outnam)
    if (modus == 2) then
      call smallhead('Ensemble optimization with loose thresholds')
      call MDopt_para(env,trim(inpnam),modus)
    else
      call smallhead('Ensemble optimization with crude thresholds')
      call MDopt_para(env,trim(inpnam),1)
    end if
    filename = trim(thispath)//'/'//trim(outnam)
    call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
    !--- get other thresholds
    newthr = aint(ewinbackup*(10.0d0/6.0d0))
    env%ewin = newthr   !
    call sort_and_check(env,trim(filename))
    call rmoptim()
    write (*,*)
!--- lastly vtight optimization
  case (3)
    call smallhead('3. optimization with very tight thresholds')
    call checkname_xyz(crefile,inpnam,outnam)
    call MDopt_para(env,trim(inpnam),0) !optlev is set from the module variable
    filename = trim(thispath)//'/'//trim(outnam)
    call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
    call sort_and_check(env,trim(filename))
    call rmoptim()
    write (*,*)
!--- crude optimization
  case (4)
    call smallhead('1. crude pre-optimization')
    call checkname_xyz(crefile,inpnam,outnam)
    call MDopt_para(env,trim(inpnam),1)
    filename = trim(thispath)//'/'//trim(outnam)
    call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
    !---using cregen to sort the optimized structures
    call checkname_xyz(crefile,inpnam,outnam)
    !--- get other thresholds
    newthr = aint(ewinbackup*2.0d0)
    env%ewin = newthr
    env%rthr = rthrbackup*2.0d0 !larger RTHR
    call sort_and_check(env,trim(filename))
    !--- only if modus=1 the input file is overwritten with the sorted file
    call remove(inpnam)
    call rename(outnam,inpnam)
    !-----
    call rmoptim()
    write (*,*)
!--- tight optimization
  case (5:6)
    if (modus .eq. 5) then
      call smallhead('2. optimization with tight thresholds')
    else if (modus .eq. 6) then
      call smallhead('Ensemble optimization with tight thresholds')
    end if
    call checkname_xyz(crefile,inpnam,outnam)
    call MDopt_para(env,trim(inpnam),4) !optlev is set from the module variable
    filename = trim(thispath)//'/'//trim(outnam)
    call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
    call sort_and_check(env,trim(filename))
    call rmoptim()
    write (*,*)
!--- default optimization
  case default
    olev = optlevflag(env%optlev)
    headder = 'Ensemble optimization with '//olev//' thresholds'
    call smallhead(headder)
    !call smallhead('Ensemble optimization')
    call checkname_xyz(crefile,inpnam,outnam)
    call MDopt_para(env,trim(inpnam),0) !optlev is set from the module variable
    filename = trim(thispath)//'/'//trim(outnam)
    call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
!        if (env%crestver .ne. crest_solv) then
    if (.not.env%QCG) then
      call sort_and_check(env,trim(filename))
      call rmoptim()
    end if
    write (*,*)
  end select

  env%thresholds = backupthr
  env%ewin = ewinbackup
  env%rthr = rthrbackup
  deallocate (backupthr)

end subroutine multilevel_opt
!--------------------------------------------------------------------------------------------
! Construct MD $set-block version 2 (updated)
!--------------------------------------------------------------------------------------------
subroutine setMDrun2(fname,hmass,mdtime,mdtemp,mdstep,shake,mddumpxyz, &
 &                   mdskip,mddump,nvt,cts)
  use crest_parameters
  use crest_data
  use utilities
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

  !--- coord setup section
  open (newunit=ich,file=fname)
  !--- don't modify coords:
  do
    read (ich,'(a)',iostat=iost) atmp
    if (iost < 0) exit
    if (index(atmp,'$coord') .ne. 0) cycle
    if (index(atmp(1:1),'$') .ne. 0) exit
  end do
  !--- then write the MD settings
  write (ich,'(a)') '$md'
  if (hmass > 0.0_wp) then      !set H-atom mass
    write (ich,'(2x,a,i0)') 'hmass=',nint(hmass)
  end if
  if (mdtime > 0.0_wp) then     !set MD simulation time in ps
    write (ich,'(2x,a,f10.2)') 'time=',mdtime
  end if
  if (mdtemp > 0.0_wp) then     !set MD Temperature
    write (ich,'(2x,a,f10.2)') 'temp=',mdtemp
  end if
  if (mdstep .gt. 0) then     !set MD timestep in fs
    write (ich,'(2x,a,f10.2)') 'step=',mdstep
  end if
  if (shake >= 0) then     !set MD shake mode
    write (ich,'(2x,a,i0)') 'shake=',shake
  end if
  !the order of setting mddump and mddumpxyz is important!!!
  if (mddumpxyz > 0) then    ! Frequency of structure dump into xtb.trj in fs
    write (ich,'(2x,a,i0)') 'dump=',mddumpxyz
  end if
  if (mdskip >= 0) then    !skipping optimizations in -mdopt
    write (ich,'(2x,a,i0)') 'skip=',mdskip
  end if
  if (nvt >= 0) then      ! NVT ensemble
    if (nvt == 0) then
      write (ich,'(2x,a)') 'nvt=false'
    else
      write (ich,'(2x,a)') 'nvt=true'
    end if
  end if
  !the order of setting mddump and mddumpxyz is important!!!
  if (mddump > 0) then     ! Vbias dump for metadyn in fs
    write (ich,'(a)') '$set'
    write (ich,'(2x,a,2x,i0)') 'mddump',mddump
  end if
  !--- add constraints (only if given, else the routines return)
  call write_cts(ich,cts)
  call write_cts_NCI(ich,cts)
  if (cts%cbonds_md) then
    call write_cts_CBONDS(ich,cts)
  end if
  if (cts%dispscal_md) then
    call write_cts_DISP(ich,cts)
  end if
  call write_cts_rcontrol(ich,cts)
  write (ich,'(a)') '$end'
  close (ich)
  return
end subroutine setMDrun2

!--------------------------------------------------------------------------------------------
subroutine setMetadyn2(fname,fac,expo,dumplist)
  use crest_parameters
  implicit none
  character(len=*) :: fname
  real(wp)  :: fac
  real(wp)  :: expo
  integer   :: dumplist
  character(len=256) :: atmp
  character(len=256),allocatable :: mdyn(:)
  integer :: i
  integer :: ich,ich2,iost,lz

  allocate (mdyn(10))
  mdyn(1:10) = ''
  lz = 1

  open (newunit=ich,file=fname)
  open (newunit=ich2,file='tmpcoordfile')
  do
    read (ich,'(a)',iostat=iost) atmp
    if (iost < 0) exit
    if (index(atmp,'$metadyn') .ne. 0) then

      do
        read (ich,'(a)',iostat=iost) atmp
        if (iost < 0) exit
        if (index(atmp,'$') .ne. 0) exit
        mdyn(lz) = adjustl(atmp)
        lz = lz+1
      end do
    end if
    write (ich2,'(a)') trim(atmp)
  end do

  write (ich2,'(a)') '$metadyn'
  do i = 1,lz
    if (index(mdyn(i),'save') .ne. 0) cycle
    if (index(mdyn(i),'kpush') .ne. 0) cycle
    if (index(mdyn(i),'alp') .ne. 0) cycle
    if (trim(mdyn(i)) .eq. '') cycle
    write (ich2,'(2x,a)') trim(mdyn(i))
  end do
  write (ich2,'(2x,a,i0)') 'save=',dumplist
  write (atmp,'(f12.6)') fac
  write (ich2,'(2x,a,a)') 'kpush=',adjustl(trim(atmp))
  write (atmp,'(f12.6)') expo
  write (ich2,'(2x,a,a)') 'alp=',adjustl(trim(atmp))
  write (ich2,'(a)') '$end'
  close (ich2)
  close (ich)
  deallocate (mdyn)
  call rename('tmpcoordfile',fname)
  return
end subroutine setMetadyn2

!---------------------------------------------------------------------
! Modified verison of the GC (confscript.v.2)
!---------------------------------------------------------------------
subroutine cross3(env)
  use crest_parameters
  use crest_data
  use iomod
  use utilities
  use cregen_interface
  implicit none
  type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
  real(wp) :: ewinbackup
  integer  :: i,imax,tmpconf,nremain
  character(len=128) :: inpnam,outnam
  character(len=512) :: thispath,tmppath

  real(wp),allocatable :: backupthr(:)

  allocate (backupthr(8))
  backupthr = env%thresholds
  ewinbackup = env%ewin

  call getcwd(thispath)

  do i = 1,1  !>-- technically it would be possible to repeat the crossing
    !>-- determine max number of new structures
    imax = min(nint(env%mdtime*50.0d0),5000)
    if (env%setgcmax) then
      imax = nint(env%gcmax)
    end if
    if (env%quick) then
      imax = nint(float(imax)*0.5d0)
    end if

    !>-- call the crossing routine
    call checkname_xyz(crefile,inpnam,outnam)
    call crest_crossing(env,imax,trim(inpnam),env%gcmaxparent)
    if (imax .lt. 1) then
      return
      exit
    end if

    !>-- optimize ensemble
    if (env%gcmultiopt) then !> for printout
      call smallhead('GC: loose pre-optimization')
    end if
    call MDopt_para(env,'confcross.xyz',2)

    !>-- optional multilevel opt
    if (env%gcmultiopt) then
      call checkname_xyz('confcross',inpnam,outnam)
      call rename('OPTIM'//'/'//'opt.xyz',trim(inpnam))
      env%ewin = aint(ewinbackup*(10.0d0/6.0d0))
      call newcregen(env,15)
      call checkname_xyz('confcross',inpnam,outnam)
      call remaining_in(inpnam,env%ewin,nremain)
      env%thresholds = backupthr
      env%ewin = ewinbackup
      call smallhead('GC: optimization with tight thresholds')
      if (env%iterativeV2) then
        call MDopt_para(env,inpnam,4)
      else
        call MDopt_para(env,inpnam,0)
      end if
    end if

    !>-- append optimized crossed structures to original crest_rotamers input
    call checkname_xyz(crefile,inpnam,outnam)
    call appendto('OPTIM'//'/'//'opt.xyz',trim(inpnam))
  end do
  call rmrfw('confcross_')
  deallocate (backupthr)
end subroutine cross3

!---------------------------------------------------------------------
! run confg check (confscript.v.2)
!---------------------------------------------------------------------
subroutine confg_chk3(env)
  use crest_data
  use cregen_interface 
  implicit none
  type(systemdata) :: env    !> MAIN SYSTEM DATA

  call ompautoset(env%threads,4,env%omp,env%MAXRUN,0) !mode=4 --> Program intern Threads max
  !>-- Special handling qcg, no RMSD,
  !    because a CMA transformed structure would cause wrong wall pot.
  if (env%crestver .eq. crest_solv) then
    call newcregen(env,6)
  else
    call newcregen(env,0)
  end if
  call ompautoset(env%threads,5,env%omp,env%MAXRUN,0) !mode=5 --> Program intern Threads min
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
  i = 1
  do
    write (dir,'(a,i0)') trim(base),i
    ex = directory_exist(trim(dir))
    if (.not.ex) then
      exit
    else
      write (str,'(a,a,''xtb.trj'')') trim(dir),'/'
      call appendto(trim(str),whichfile)
      i = i+1
    end if
  end do
  return
end subroutine collect_trj

!---------------------------------------------------------------------
! collect xtb.trj, but skip the first point on the trj
!---------------------------------------------------------------------
subroutine collect_trj_skipfirst(base,whichfile)
  use iomod
  use utilities
  implicit none
  character(len=*) :: base,whichfile
  character(len=256) :: str,dir
  integer :: i
  logical :: ex
  i = 1
  do
    write (dir,'(a,i0)') trim(base),i
    ex = directory_exist(trim(dir))
    if (.not.ex) then
      exit
    else
      write (str,'(a,a,''xtb.trj'')') trim(dir),'/'
      call TRJappendto_skipfirst(trim(str),whichfile)
      i = i+1
    end if
  end do
  return
end subroutine collect_trj_skipfirst

!------------------------------------------------------------------------
! writes the Input coord to the first position of crest_rotamers_*.xyz
!------------------------------------------------------------------------
subroutine append_INPUT_to(fnam,newtag)
  use iomod
  use crest_data
  use strucrd,only:coord2xyz
  use utilities
  implicit none
  character(len=*)  :: fnam
  character(len=80) :: tmpnam,inam,onam
  character(len=*) :: newtag

  tmpnam = 'ensemble.xyz.tmp'
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
  use crest_parameters
  use crest_restartlog, only: restart_write_dummy
  use strucrd,only:rdensembleparam,rdensemble
  implicit none
  integer :: nall
  real(wp) :: ewin
  character(len=*) :: filename
  integer :: k,nat

  call restart_write_dummy(trim(filename))

  open (newunit=k,file=trim(filename))
  read (k,*) nat
  close (k)

  call rdensembleparam(trim(filename),nat,nall)

  if (nall .lt. 1) then
    write (*,*) 'No conformer was left. Something must be seriously wrong.'
    write (*,*) 'Terminating the run.'
    error stop
  end if

  write (*,'(1x,i0,'' structures remain within '',f8.2,'' kcal/mol window'')') &
  &        nall,ewin
  return
end subroutine remaining_in

!-------------------------------------------------------------------------
! print the number of remaining files in an ensemble file for a given energy window
!-------------------------------------------------------------------------
subroutine sort_and_check(env,filename)
  use crest_parameters
  use crest_data
  use iomod
  use strucrd,only:rdensembleparam,rdensemble
  use utilities
  implicit none

  type(systemdata) :: env
  real(wp) :: ewin
  character(len=*) :: filename
  integer :: nallin,nallout,nallthr
  real(wp) :: nthr,increase

  character(len=80) :: inpnam,outnam

  increase = 1.5d0  ! factor 1.5 increase
  nthr = 0.05d0     ! 5%

  ewin = env%ewin

  call rdensembleparam(trim(filename),env%nat,nallin) !--initial number of structures
  call checkname_xyz(crefile,inpnam,outnam)
  call confg_chk3(env)
  call remaining_in(outnam,ewin,nallout) !--- remaining number of structures

  if (.not.env%entropic) then !don't do this for entropy mode
    nallthr = nint(float(nallin)*nthr)
    if (nallout .lt. nallthr) then
      write (*,'(1x,''This is less than '',i0,''% of the initial '',i0,'' structures.'')') &
      & nint(nthr*100.0d0),nallin
      write (*,'(1x,''Increasing energy window to include more...'')')
      call remove(outnam)                    !--- delete this file since it contains to few structures
      env%ewin = ewin*increase
      call confg_chk3(env)
      call remaining_in(outnam,ewin*increase,nallout) !--- remaining number of structures
    end if
  end if

  return
end subroutine sort_and_check

!--------------------------------------------------------------------------------------------
! collect all the backup ensembles from the different runs into a new file.
!--------------------------------------------------------------------------------------------
subroutine collectcre(env)
  use iomod
  use crest_data
  use utilities
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

  do i = 1,env%Maxrestart
    call checkname_xyz('.cre',atmp,btmp)
    inquire (file=atmp,exist=ex)
    if (.not.ex) exit
    call appendto(trim(atmp),trim(crename))
    call remove(atmp)
  end do
  call rmrfw('.cre_')
  return
end subroutine collectcre

!--------------------------------------------------------------------------------------------
! check lower found
!--------------------------------------------------------------------------------------------
subroutine elowcheck(lower,env)
  use crest_parameters
  use iomod
  use crest_data
  use strucrd,only:xyz2coord
  use utilities 
  implicit none
  type(systemdata) :: env
  real(wp) :: ediff,ethr,ewin
  logical :: lower

!--- some defaults
  ethr = env%ethr  !ETHR threshold that is also used in CREGEN (=0.05 kcal/mol)
  if (env%entropic) then !except for entropy mode
    ethr = 0.1_wp
  end if
  ewin = env%ewin
  lower = .false.

  ediff = (env%eprivious-env%elowest)*autokcal
  if (ediff .ge. ethr) then
    write (*,'(1x,a)') '...............................................'
    write (*,'(1x,a)') 'A new lower conformer was found!'
    write (*,'(1x,a,f10.5,a,f10.5,a)') 'Improved by ',ediff/autokcal,' Eh or ',ediff,'kcal/mol'
    write (*,'(1x,a)') '...............................................'
    env%eprivious = env%elowest
    lower = .true.
    if (.not.env%allowrestart) then
      lower = .false.
      return
    end if
    if (env%entropic) then
      write (*,'(1x,a)') 'Restart based on new lowest structure.'
      call clean_V2i
      call rmrfw('crest_entropy_rotamer_')
      return
    end if
    !---- clean the dir
    call clean_V2i
    !---- save the new best conformer
    call XYZappendto('crest_best.xyz','.history.xyz')
    call xyz2coord('crest_best.xyz','coord') !new reference coord to start the MTDs with
  else
    lower = .false.
  end if

  return
end subroutine elowcheck

!-----------------------------------------------------------------------
! diatomic "error" catcher ---> we will have no conformers
!-----------------------------------------------------------------------
subroutine catchdiatomic(env)
  use crest_parameters
  use crest_data
  use strucrd
  use iomod
  implicit none
  type(systemdata) :: env
  type(coord) :: mol
  integer :: ich
  character(len=1024) :: jobcall

  call mol%open('coord')

  open (file=conformerfile,newunit=ich)
  mol%xyz = mol%xyz*bohr !to ang
  call wrxyz(ich,mol%nat,mol%at,mol%xyz)
  close (ich)
  !create the system call (it is the same for every optimization)
  write (jobcall,'(a,1x,a,1x,a,'' --opt '',a,1x,a,'' --ceasefiles  >xtb.out'')') &
 &    trim(env%ProgName),conformerfile,trim(env%gfnver),trim(env%solv),' 2>/dev/null'
  call command(trim(jobcall), ich)
  call copy('xtbopt.xyz',conformerfile)
  call copy(conformerfile,'crest_rotamers.xyz')
  call copy(conformerfile,'crest_best.xyz')

  call mol%deallocate
  return
end subroutine catchdiatomic

!-------------------------------------------------------------------------------------!
! entropy ensemble file copy routine
!-------------------------------------------------------------------------------------!
subroutine emtdcopy(env,iter,stopiter,broken)
  use crest_parameters
  use crest_data
  use iomod
  use strucrd
  use utilities 
  use cregen_interface
  implicit none
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
  character(len=22),parameter :: sfile = 'crest_entropy_rotamer_'
  character(len=11),parameter :: bfile = 'crest_smtd_'
  character(len=:),allocatable :: filname
  real(wp) :: T,S,Cp
  integer :: nt

  stopiter = .false.
  broken = .false.
  T = 298.15d0
!>--- determine temperature dependence
  if (.not.allocated(env%thermo%temps)) then
    call env%thermo%get_temps()
  end if
  nt = env%thermo%ntemps
!>--- space to save S,Cp,H(T) at different temperatures
  if (.not.allocated(env%emtd%soft)) then
    allocate (env%emtd%soft(nt),source=0.0d0)
  end if
  if (.not.allocated(env%emtd%cpoft)) then
    allocate (env%emtd%cpoft(nt),source=0.0d0)
  end if
  if (.not.allocated(env%emtd%hoft)) then
    allocate (env%emtd%hoft(nt),source=0.0d0)
  end if

!>--- output file name
  if (env%crestver == crest_imtd2) then
    filname = trim(bfile)
  else
    filname = trim(sfile)
  end if

!>--- Setup in the very first call
  if (iter == 0) then
    call checkname_xyz(crefile,atmp,crename)
    write (btmp,'(a,i0,a)') filname,iter,'.xyz'
    call rename(trim(atmp),trim(btmp))
    call rmrfw('crest_rotamers_')
    call rmrfw('.data')
    env%emtd%sapproxlast = env%emtd%sapprox
    call rdensembleparam(conformerfile,i,nall)
    env%emtd%nconflast = nall
    call entropic(env,.false.,.true.,.true.,trim(btmp),T,S,Cp)
    call writesdata(env,nall,iter)
    return
  end if

!>--- Setup in the iterations
  if (iter >= 1) then
    call checkname_xyz(crefile,atmp,crename)
    iter2 = iter-1
    write (btmp,'(a,i0,a)') filname,iter2,'.xyz' !-- file from last iter.
    inquire (file=trim(btmp),exist=ex)
    if (ex) then
      call appendto(trim(btmp),trim(atmp)) !-- cat together
      call newcregen(env,2)
      sdiff = (env%emtd%sapproxlast/env%emtd%sapprox)
      sdiff = 1.0d0-sdiff
      sdiff = abs(sdiff)     !can dS be negative?
      if (sdiff < 0.0d0) then
        broken = .true.
      else
        broken = .false.
        env%emtd%sapproxlast = env%emtd%sapprox
      end if
    end if
  end if

!>--- run checks
  if (.not.broken) then
    call checkname_xyz(crefile,crename,btmp)
    write (btmp,'(a,i0,a)') filname,iter,'.xyz'
    call rename(trim(crename),trim(btmp))
  end if

  call rdensembleparam(conformerfile,i,nall)
  if (nall > 50000) then !safety fallback for extremely large ensembles (e.g. C18)
    env%emtd%confthr = 0.1d0
  end if
!>---get convergence criteria
  nallfrac = 1.0d0-(float(env%emtd%nconflast)/float(nall))
  conv1 = (sdiff < env%emtd%sconvthr) !.and.(sdiff > 0.0d0)
  conv2 = (nallfrac < env%emtd%confthr).and.(nallfrac >= 0.0d0)
  if (nallfrac < 0.0d0) then !> if we for some reason got less conformers in this iteration
    broken = .true.
  end if
  write (*,*)
  write (*,'(1x,a,i0,a)') 'Static MTD simulation iteration ',iter,' finished'
  write (*,'(1x,a,i0,a,f5.2,a,a)') 'A total of ',nall,' conformers (+',  &
  &    nallfrac*100.0d0,'%) in ',conformerfile
  write (*,'(1x,a,f12.6,a,f5.2,a)') 'S(approx.) =',env%emtd%sapprox, &
  &     ' cal/molK (dS =',sdiff*100.0d0,'%)'
  write (*,'(1x,a,l)') 'Convergence w.r.t. conformers: ',conv2
  write (*,'(1x,a,l)') 'Convergence w.r.t. entropy   : ',conv1
  write (*,*)

  if (.not.broken) then !-- NORMAL case
    if (conv1.and.conv2) then
      stopiter = .true.
      if (iter < env%emtd%iter) then
        call copy(trim(btmp),trim(crename))
      end if
    end if
    if (env%entropic.and.env%crestver .ne. crest_imtd2) then
      call entropic(env,.false.,.true.,.true.,trim(btmp),T,S,Cp)
      call writesdata(env,nall,iter)
    end if
!>--- save data for next iteraton
    env%emtd%nconflast = nall
    env%emtd%sapproxlast = env%emtd%sapprox
  else                !-- ROLLBACK
    !--- NOTE: instead of exiting one could also try to re-do the iteration,
    !          but this is commented out. Relevant comments are marked by "!*"
    write (*,'(1x,a)') 'Warning: Nconf DECREASED, which is unphysical.'
    !write(*,'(1x,a)') 'Trying to re-do last iteration ...'   !*
    write (*,'(1x,a)') 'Rewind ensemble and exit iterations ...'    !*
    call checkname_xyz(crefile,crename,btmp)
    call remove(crename)
    write (btmp,'(a,i0,a)') filname,iter-1,'.xyz'
    call copy(trim(btmp),trim(crename))
    call xyz2coord(trim(crename),'coord')  !get the right coord again
    call newcregen(env,2) !get the right crest_conformers again
    iter = iter-1  !*
    stopiter = .true.
  end if

  if (iter < env%emtd%iter.and..not.stopiter) then
    call rmrfw('crest_rotamers_')
  end if
  return
contains
  subroutine writesdata(env,nall,inum)
    implicit none
    type(systemdata) :: env
    integer :: nall
    integer :: inum
    character(len=64) :: adum
    integer :: ich,nt,i
    nt = env%thermo%ntemps
    write (adum,'(a,i0)') 'Sdata',inum
    open (newunit=ich,file=trim(adum))
    write (ich,*) nall
    do i = 1,nt
      write (ich,'(3F18.10)') env%emtd%soft(i), &
      &   env%emtd%Cpoft(i),env%emtd%hoft(i)
    end do
    close (ich)
    return
  end subroutine writesdata
end subroutine emtdcopy

subroutine emtdcheckempty(env,empty,nbias)
  use crest_data
  use strucrd
  use iomod
  use utilities 
  implicit none
  type(systemdata) :: env
  logical :: empty
  integer :: nbias
  integer :: i,nall
  character(len=128) :: atmp,btmp
  call checkname_xyz(crefile,atmp,btmp)
  call rdensembleparam(trim(atmp),i,nall)
  empty = .false.

  if (nall < 1) then
    call remove(trim(atmp))
    if (env%emtd%maxfallback > 1) then
      write (*,'(1x,a)') 'Empty ensemble, retry iteration with less bias ...'
    else
      write (*,'(1x,a)') 'Empty ensemble, skipping ahead ...'
    end if
    call remove('crest_clustered.xyz')
    call remove('crest_bias.xyz')
    nbias = nbias-1
    nbias = max(0,nbias)
    empty = .true.
  end if
  return
end subroutine emtdcheckempty

!========================================================================!
! Sample additional OH orientations
!========================================================================!
subroutine XHorient(env,infile)
  use crest_parameters
  use crest_data
  use strucrd
  use zdata
  use iomod
  implicit none
  type(systemdata) :: env
  character(len=*) :: infile
  integer :: nat,nall
  character(len=:),allocatable :: newfile
  newfile = 'oh_ensemble.xyz'
  !--- shouldn't cost much
  call ohflip_ensemble(infile,env%maxflip)
  call rdensembleparam(newfile,nat,nall)
  !--- only proceed if there are potential new structures
  if (nall > 0) then
    call smallhead('Additional orientation sampling')
    call MDopt_para(env,newfile,0)
    !---- printout and copy
    call rename('OPTIM'//'/'//'opt.xyz',trim(newfile))
    call rmrf('OPTIM')
  else
    call remove(newfile)
  end if
  return
end subroutine XHorient
