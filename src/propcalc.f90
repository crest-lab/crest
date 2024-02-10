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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Routines related to additonal property calculations
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine protreffrag(env)
  use iso_fortran_env,wp => real64
  use crest_data
  use strucrd,only:rdnat,rdcoord
  implicit none
  type(systemdata) :: env
  integer,allocatable  :: molvec(:)
  integer,allocatable  :: at(:)
  real(wp),allocatable :: xyz(:,:)

  associate (nat => env%nat)
!------ get number of fragments for original structure
    allocate (xyz(3,nat),at(nat),molvec(nat))
    call rdcoord('coord',nat,at,xyz)
    call mrec(env%ptb%nfrag,xyz,nat,at,molvec) !requires xyz in bohr
    deallocate (molvec,at,xyz)
  end associate
end subroutine protreffrag

!---------------------------------------------------------------------------------------
! perform a property calculation for a given ensemble file
!---------------------------------------------------------------------------------------
subroutine propcalc(iname,imode,env,tim)
  use iso_fortran_env,wp => real64
  use crest_data
  use iomod
  use strucrd,only:rdensembleparam,rdensemble,wrxyz
  use utilities,only:boltz2
  use cregen_interface
  implicit none

  type(systemdata) :: env
  type(timer)      :: tim
  integer :: imode

  character(len=*),intent(in)  :: iname   !file name

  integer :: i,k,r,ii
  integer :: TMPCONF
  integer :: P
  integer :: ich,ich2

  interface
    subroutine prop_OMP_loop(env,TMPCONF,jobcall,pop)
      import :: systemdata,wp
      implicit none

      type(systemdata) :: env
      integer :: TMPCONF
      character(len=1024) :: jobcall
      real(wp),intent(in),optional :: pop(TMPCONF)
    end subroutine
  end interface

  character(len=20) :: xname
  character(len=20) :: pipe
  character(len=80) :: solv
  character(len=256) :: ctmp
  character(len=512) :: str,thispath,tmppath,optpath
  character(len=1024):: jobcall
  character(len=:),allocatable :: largejobcall

  real(wp) :: pthr,sumpop
  integer :: maxpop
  integer :: nat,nall,ng
  logical :: ex,update
  logical :: niceprint

  character(len=40),allocatable :: origin(:)
  real(wp),allocatable :: eread(:),popul(:),dumm(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  logical,allocatable :: mask(:)
  integer,allocatable :: degen(:,:)

!---
  call largehead('P R O P E R T Y   C A L C U L A T I O N')

!--- some settings
  solv = ''
  pipe = '2>/dev/null'
  xname = 'struc.xyz'
  call getcwd(thispath)
  update = .true.
  maxpop = 1

  if (env%properties == -666) then
    imode = env%properties2
  end if

  niceprint = env%niceprint

  pthr = env%pthr

!---- read the input ensemble
  call rdensembleparam(iname,nat,nall)
  env%nat = nat
  if (nall .lt. 1) return
  allocate (xyz(3,nat,nall),at(nat),eread(nall),mask(nall))
  mask = .true.
  if (.not.env%trackorigin) then
    call rdensemble(iname,nat,nall,at,xyz,eread)
  else
    allocate (origin(nall))
    call rdensemble_origin(iname,nat,nall,at,xyz,eread,origin)
  end if
  TMPCONF = nall

!---- Use only populated conformers for special applications
  if (any((/2,111/) == imode)) then

    allocate (popul(nall))
    call boltz2(nall,eread,popul)
    !write(*,*) eread
    !write(*,*) popul
    k = 0
    maxpop = maxloc(popul(:),1) !locate max. populated structure (this one is always taken)
    sumpop = 0.0_wp
    do i = 1,nall
      if (popul(i) .ge. pthr.or.i .eq. maxpop) then
        k = k+1
      else
        exit
      end if
    end do
    write (*,'(1x,a,i0,a)') 'Population threshold (-pthr) : ',nint(pthr*100.0_wp),' %'

    write (*,'(1x,i0,a,i0,a)') k,' populated structure(s) (out of a total ', &
    & nall,') will be considered.'

    mask(:) = popul(:) .ge. pthr
    mask(maxpop) = .true.
  end if

!---- for multilevel reoptimization don't use all rotamers
  if (imode .ge. 50.and.imode .lt. 60) then
    inquire (file='cre_members',exist=ex)
    if (ex) then
      open (newunit=ich,file='cre_members')
      read (ich,*) ng
      allocate (degen(3,ng))
      do i = 1,ng
        read (ich,*) degen(1:3,i)
      end do
      close (ich)
      !-- always include lowest rotamer for all conf.
      do i = 1,ng
        mask(degen(2,i)) = .true.
      end do
    end if
  end if

!---- create the PROP directory
  !---- create directory for the optimizations
  optpath = 'PROP'
  call rmrf(optpath)
  r = makedir(trim(optpath))

  call copysub('coord',trim(optpath))

  call env%wrtCHRG(trim(optpath))
  call copysub(env%fixfile,trim(optpath))
  !call copysub(env%constraints,trim(optpath))
  if (env%useqmdff) then
    call copysub('solvent',trim(optpath))
  end if
  if (env%gfnver == '--gff') then
    r = sylnk(trim(thispath)//'/'//'gfnff_topo',trim(optpath)//'/'//'gfnff_topo')
  end if

  call chdir(trim(optpath))
  call getcwd(optpath)
!---- set up sub-directories
  write (*,'(1x,a,a,a)',advance='no') 'writing TMPCONF* Dirs from file "',trim(iname),'" ...'
  ii = 1
  do i = 1,nall
    if ((imode .eq. 2).and.(allocated(popul))) then
      if ((popul(i) .lt. pthr).and.i .ne. maxpop) cycle !skip unpopulated for '-prop autoir'
    end if
    if (.not.mask(i)) cycle

    write (ctmp,'(''TMPCONF'',i0)') ii
    r = makedir(trim(ctmp))
    call chdir(ctmp)
    open (newunit=ich,file=xname)
    call wrxyz(ich,nat,at,xyz(:,:,i))

    call write_cts(ich,env%cts)
    close (ich)

    if (imode .eq. 2.and.trim(env%gfnver) .eq. '--gfn2') then
      call add_mass_xtb(xname)
    end if

    call chdir(optpath)

    call env%wrtCHRG(trim(ctmp))
    call copysub(env%fixfile,trim(ctmp))
    if (env%useqmdff) then
      call copysub('solvent',trim(ctmp))
    end if
    if (env%gfnver == '--gff') then
      r = sylnk(trim(optpath)//'/'//'gfnff_topo',trim(ctmp)//'/'//'gfnff_topo')
    end if

    ii = ii+1
  end do
  write (*,'(1x,a)') 'done.'

  if (any(mask.eqv..false.)) then
    TMPCONF = ii-1
    nall = ii-1
    if (allocated(popul)) deallocate (popul)
  end if

!--- setting the threads for correct parallelization
  if (env%autothreads) then
    call ompautoset(env%threads,7,env%omp,env%MAXRUN,TMPCONF) !set global OMP/MKL variable for xtb jobs
  end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  call chdir(thispath)
!--- select what to do
  P = imode
  select case (P)
  case (1)
    call smallhead('Hessian calculations for all conformers')
    write (jobcall,'(a,1x,a,1x,a,'' --hess '',a,1x,a,'' >xtb.out'')') &
    &    trim(env%ProgName),trim(xname),trim(env%gfnver),trim(env%solv),trim(pipe)
  case (10)
    call smallhead('Optimization + Hessian calculations for all conformers')
    write (jobcall,'(a,1x,a,1x,a,'' --ohess '',a,1x,a,'' >xtb.out'')') &
    &    trim(env%ProgName),trim(xname),trim(env%gfnver),trim(env%solv),trim(pipe)
  case (13)
    call smallhead('Free energy calculation in solvation')
    write (jobcall,'(a,1x,a,1x,a,'' --sp '',a,'' >sp.out'')') &
    &    trim(env%ProgName),trim(xname),trim(env%gfnver),trim(pipe) !E_gas(Solv_geom) singlepoint
    largejobcall = trim(jobcall)//' ; '
    write (jobcall,'(a,1x,a,1x,a,'' --ohess '',a,1x,a,'' >xtb.out'')') &
    &    trim(env%ProgName),trim(xname),trim(env%gfnver),trim(env%solv),trim(pipe)
    jobcall = largejobcall//trim(jobcall)
  case (2)
    call smallhead('IR calculation for populated conformers')
    write (jobcall,'(a,1x,a,1x,a,'' --ohess '',a,1x,a,'' >xtb.out'')') &
    &    trim(env%ProgName),trim(xname),trim(env%gfnver),trim(env%solv),trim(pipe)
!       case( 3:6,7,8,100 ) ! unspecific case DFT
!         call smallhead('DFT calculation using xtb as driver')
!         if( any((/3,4/)==P) )then
!           call dftrc_reader(env,.true.)  !B97-3c OPT default
!           call dftTMwarning
!         else
!           call dftrc_reader(env,.false.) !read DFT settings
!         endif
!         call chdir(optpath)
!         call cefine_setup(env,TMPCONF)
!         call xtbDFTdriver(env,xname,jobcall) !create jobcall
  case (20)
    call smallhead('Reoptimization for all conformers')
    write (jobcall,'(a,1x,a,1x,a,'' --opt vtight '',a,1x,a,'' >xtb.out'')') &
    &    trim(env%ProgName),trim(xname),trim(env%gfnver),trim(env%solv),trim(pipe)
  case (50:59)
    call smallhead('Reoptimization of entire CRE')
    write (jobcall,'(a,1x,a,1x,a,'' --opt vtight '',a,1x,a,'' >xtb.out'')') &
    &    trim(env%ProgName),trim(xname),trim(env%gfnver2),trim(env%solv),trim(pipe)
  case default
    write (jobcall,'(a,1x,a,1x,a,'' --sp '',a,1x,a,'' >xtb.out'')') &
    &    trim(env%ProgName),trim(xname),trim(env%gfnver),trim(env%solv),trim(pipe)
  end select
  call chdir(optpath)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  write (*,'(1x,a,i0,a)') 'Performing calculations for ', &
  & TMPCONF,' structures ...'
  call sleep(1)
  call tim%start(10,'PROPERTY calc.')
  allocate (dumm(TMPCONF),source=1.0_wp)
  call prop_OMP_loop(env,TMPCONF,jobcall,dumm)  !<------- this is where the "magic" happens
  deallocate (dumm)
  write (*,*)
  call tim%stop(10)
  write (*,*) 'done.'

  call tim%start(10,'PROPERTY calc.')

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!--- select what to do with the output
  select case (imode)

!------ generate an ensemble with free energies from hessian calculations
  case (1,10)
    call rdpropens(TMPCONF,nat,xyz) !get updated geometries
    open (newunit=ich,file='crest_property.xyz')
    do i = 1,TMPCONF
      write (tmppath,'(a,i0,a,a)') 'TMPCONF',i,'/','xtb.out'
      call grepval(tmppath,'TOTAL FREE ENERGY',ex,eread(i))
      if (.not.env%trackorigin) then
        write (str,'(3x,f16.10)') eread(i)
      else
        write (str,'(3x,f16.10,2x,a,a)') eread(i),'!',trim(origin(i))
      end if
      call wrxyz(ich,nat,at,xyz(:,:,i),trim(str))
    end do
    close (ich)
    env%ensemblename = 'crest_property.xyz'
    env%confgo = .true.

    call newcregen(env,0)

    call rename('crest_property.xyz.sorted', &
    & trim(thispath)//'/crest_property.xyz')
    call remove(env%ensemblename)
!------ IR averaging
  case (2)
    call autoir(TMPCONF,imode,env)
    call rdpropens(TMPCONF,nat,xyz) !get updated geometries
    call wrpropens(TMPCONF,nat,xyz,at,eread)
!!----- DFT handling
!       case( 3:8,100 )
!         call DFTprocessing(env,TMPCONF,nat,at)
!------ vtight reoptimization only for conformers!
  case (20)
    call rdpropens(TMPCONF,nat,xyz) !get updated geometries
    call wrpropens(TMPCONF,nat,xyz,at,eread)
    env%ensemblename = 'crest_property.xyz'
    env%confgo = .true.

    call newcregen(env,0)

    call rename('crest_ensemble.xyz', &
    & trim(thispath)//'/crest_conformers.xyz')
  case (50:59)
    call rdpropens(TMPCONF,nat,xyz) !get updated geometries
    env%ensemblename = 'crest_reopt.xyz'
    open (newunit=ich,file=env%ensemblename)
    do i = 1,TMPCONF
      write (tmppath,'(a,i0,a,a)') 'TMPCONF',i,'/','xtb.out'
      call grepval(tmppath,'total energy',ex,eread(i))
      if (.not.env%trackorigin) then
        write (str,'(3x,f16.10)') eread(i)
      else
        write (str,'(3x,f16.10,2x,a,a)') eread(i),'!',trim(origin(i))
      end if
      call wrxyz(ich,nat,at,xyz(:,:,i),trim(str))
    end do
    close (ich)
    call rename(env%ensemblename,trim(thispath)//'/'//env%ensemblename)
    call chdir(thispath)
    if (imode .lt. 59) then !TODO temporary skip for some testing
      env%confgo = .true.

      call newcregen(env,0)

      env%confgo = .false.
      call rename(trim(env%ensemblename)//'.sorted', &
      & env%ensemblename)
    end if

  case (998) !singlpoint (no reranking) + dipoles
    open (newunit=ich,file='crest_property.xyz')
    open (newunit=ich2,file='crest.dipoles')
    do i = 1,TMPCONF
      write (tmppath,'(a,i0,a,a)') 'TMPCONF',i,'/','xtb.out'
      call grepval(tmppath,'| TOTAL ENERGY',ex,eread(i))
      write (str,'(3x,f16.10)') eread(i)
      call wrxyz(ich,nat,at,xyz(:,:,i),trim(str))

      call grepcntxt(tmppath,'molecular dipole:',ex,ctmp,3)
      if (ex) then
        write (ich2,'(a)') trim(ctmp(10:))
      else
        write (ich2,'(a)') ''
      end if
    end do
    close (ich)
    close (ich2)

    call rename('crest.dipoles', &
    & trim(thispath)//'/'//'crest.dipoles')

    write (*,*)
    write (*,*) 'Dipole moments for each conformer (x,y,z,total) written to crest.dipoles'

  case (999) !singlpoint reranking
    open (newunit=ich,file='crest_property.xyz')
    do i = 1,TMPCONF
      write (tmppath,'(a,i0,a,a)') 'TMPCONF',i,'/','xtb.out'
      call grepval(tmppath,'| TOTAL ENERGY',ex,eread(i))
      write (str,'(3x,f16.10)') eread(i)
      call wrxyz(ich,nat,at,xyz(:,:,i),trim(str))
    end do
    close (ich)
    env%ensemblename = 'crest_property.xyz'
    env%confgo = .true.

    call newcregen(env,0)

    call rename('crest_property.xyz.sorted', &
    & trim(thispath)//'/crest_property.xyz')
    call remove(env%ensemblename)
    call rename('crest_ensemble.xyz', &
    & trim(thispath)//'/crest_ensemble.xyz')

  case default
    continue
  end select

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!>--- move around some files

  inquire (file='cregen.out.tmp',exist=ex)
  if (ex) then
    call cat('cregen.out.tmp')
  end if

  inquire (file='crest.vibspectrum',exist=ex)
  if (ex) then
    call rename('crest.vibspectrum', &
    & trim(thispath)//'/'//'crest.vibspectrum')
  end if

  inquire (file='crest_property.xyz',exist=ex)
  if (ex) then
    call rename('crest_property.xyz', &
    & trim(thispath)//'/'//'crest_property.xyz')
  end if

  inquire (file='crest_populated.xyz',exist=ex)
  if (ex) then
    call rename('crest_populated.xyz', &
    & trim(thispath)//'/'//'crest_populated.xyz')
  end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!>--- stop timer

  call tim%stop(10)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!>--- change back to original directory (just to be sure)

  call chdir(thispath)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!>--- decide if something has to be cleaned up here

  if (env%pclean) then
    select case (imode)
    case (20,50:59,998,999)
      call rmrf('PROP')
    case default
      continue
    end select
  end if

  if (allocated(xyz)) deallocate (xyz)
  if (allocated(at)) deallocate (at)
  if (allocated(eread)) deallocate (eread)
  if (allocated(origin)) deallocate (origin)
  if (allocated(mask)) deallocate (mask)
  if (allocated(popul)) deallocate (popul)

end subroutine propcalc

!----------------------------------------------------------------------------------------------------
! THE OMP-PARALLEL LOOP
!----------------------------------------------------------------------------------------------------
subroutine prop_OMP_loop(env,TMPCONF,jobcall,pop)
  use iso_fortran_env,wp => real64
  use crest_data
  use iomod
  implicit none

  type(systemdata) :: env
  integer :: TMPCONF
  character(len=1024) :: jobcall
  real(wp),intent(in),optional :: pop(TMPCONF)

  real(wp) :: pthr
  logical :: niceprint
  real(wp) :: percent
  integer :: vz,k,i,maxpop,io
  character(len=512) :: tmppath

!----- quick settings
  niceprint = env%niceprint
  pthr = env%pthr
  maxpop = maxloc(pop(:),1)
  if (niceprint) then
    call printprogbar(0.0_wp)
  end if
  k = 0        ! count finished jobs

!$omp parallel &
!$omp shared( vz,jobcall,TMPCONF,pop,pthr,percent,k,niceprint,maxpop )
!$omp single
  do i = 1,TMPCONF
    vz = i
    !$omp task firstprivate( vz ) private( tmppath,io )
    call initsignal()
    !$omp critical
    write (tmppath,'(a,i0)') 'TMPCONF',vz
    !$omp end critical
    if (pop(vz) .ge. pthr.or.vz .eq. maxpop) then
      call command('cd '//trim(tmppath)//' && '//trim(jobcall),io)
    end if
    !$omp critical
    k = k+1
    if (niceprint) then
      percent = float(k)/float(TMPCONF)*100
      call printprogbar(percent)
    else
      write (6,'(1x,i0)',advance='no') k
      flush (6)
    end if
    !$omp end critical
    !$omp end task
  end do
!$omp taskwait
!$omp end single
!$omp end parallel

  return
end subroutine prop_OMP_loop

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! grep total energies and printout energy list
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine etotprop(TMPCONF,pop,pr)
  use iso_fortran_env,wp => real64
  use crest_data
  use iomod
  use utilities,only:boltz2
  implicit none
  integer,intent(in) :: TMPCONF
  real(wp),intent(inout) :: pop(TMPCONF)
  logical :: pr
  character(len=512) :: tmppath
  logical :: ex
  integer :: i
  real(wp) :: dE
  real(wp),allocatable :: eread(:)

  real(wp),parameter :: kcal = 627.5095_wp

  write (*,'(1x,a)') "Calculating populations from total energies ..."
  allocate (eread(TMPCONF))
  do i = 1,TMPCONF
    write (tmppath,'(a,i0,a,a)') 'TMPCONF',i,'/','xtb.out'
    call grepval(tmppath,'total energy',ex,eread(i))
  end do

!--- convert to populations
  call boltz2(TMPCONF,eread,pop)
  if (pr) then
    write (*,'(a)') '========================================================='
    write (*,'(a)') '============= total energies & populations  ============='
    write (*,'(a)') '========================================================='
    write (*,'('' structure    ΔE(kcal/mol)    Etot(Eh)        weight'')')
    do i = 1,TMPCONF
      dE = (eread(i)-eread(1))*kcal
      write (*,'(i5,6x,F10.2,4x,F14.6,F13.4)') i,dE,eread(i),pop(i)
    end do
    write (*,*)
  end if
  deallocate (eread)
  return
end subroutine etotprop

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! grep optimized geometries (ONLY that, no energies)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine rdpropens(TMPCONF,n,xyz)
  use iso_fortran_env,wp => real64
  use crest_data
  use iomod
  implicit none

  integer,intent(in)  :: TMPCONF
  integer,intent(in)  :: n
  real(wp),intent(inout) :: xyz(3,n,TMPCONF)

  integer :: ich
  character(len=512) :: tmppath,atmp
  character(len=64)  :: dum
  logical :: ex
  integer :: i,j

  real(wp),parameter :: kcal = 627.5095_wp

  do i = 1,TMPCONF
    write (tmppath,'(a,i0,a,a)') 'TMPCONF',i,'/','xtbopt.xyz'
    inquire (file=tmppath,exist=ex)
    if (.not.ex) then
      write (tmppath,'(a,i0,a,a)') 'TMPCONF',i,'/','struc.xyz'
    end if
    open (newunit=ich,file=tmppath)
    read (ich,'(a)') atmp
    read (ich,'(a)') atmp
    do j = 1,n
      read (ich,*) dum,xyz(1:3,j,i)
    end do
    close (ich)
  end do

  return
end subroutine rdpropens
subroutine wrpropens(TMPCONF,n,xyz,at,eread)
  use iso_fortran_env,wp => real64
  use crest_data
  use iomod
  use strucrd,only:wrxyz
  implicit none

  integer,intent(in)    :: TMPCONF
  integer,intent(in)    :: n
  real(wp),intent(in)   :: xyz(3,n,TMPCONF)
  real(wp),intent(inout):: eread(TMPCONF)
  integer,intent(in)    :: at(n)
  integer :: ich
  character(len=512) :: tmppath
  logical :: ex
  integer :: i

  open (newunit=ich,file='crest_property.xyz')
  do i = 1,TMPCONF
    write (tmppath,'(a,i0,a,a)') 'TMPCONF',i,'/','xtb.out'
    call grepval(tmppath,'total energy',ex,eread(i))
    call wrxyz(ich,n,at,xyz(:,:,i),eread(i))
  end do
  return
end subroutine wrpropens
subroutine wrpropens_pop(env,TMPCONF,n,xyz,at,eread,pthr)
  use iso_fortran_env,wp => real64
  use crest_data
  use iomod
  use strucrd,only:wrxyz
  use utilities,only:boltz2
  implicit none

  type(systemdata) :: env
  integer,intent(in)    :: TMPCONF
  integer,intent(in)    :: n
  real(wp),intent(in)   :: xyz(3,n,TMPCONF)
  real(wp),intent(in)   :: eread(TMPCONF)
  integer,intent(in)    :: at(n)
  real(wp),intent(in)   :: pthr
  integer :: ich,ich5
  integer :: i
  integer :: pmax
  real(wp),allocatable :: pop(:),dpop(:),edum(:)

  allocate (pop(TMPCONF),dpop(TMPCONF))
  call boltz2(TMPCONF,eread,pop)
  pmax = maxloc(pop,1)

  dpop = pop

!      if(env%hardcutDFT)then
!         call cutDFTpop(env,dpop,TMPCONF)
!         allocate(edum(TMPCONF))
!         edum=eread
!         do i=1,TMPCONF
!            if(dpop(i).le.0.00001_wp)then
!              edum(i)=0.0_wp
!            endif
!         enddo
!         call boltz2(TMPCONF,edum,pop)
!         deallocate(edum)
!      endif

  open (newunit=ich5,file='autoir.pop')
  open (newunit=ich,file='crest_populated.xyz')
  do i = 1,TMPCONF
    if (dpop(i) .lt. pthr.and.i .ne. pmax) cycle
    write (ich5,'(1x,i0,1x,f6.4)') i,pop(i)
    call wrxyz(ich,n,at,xyz(:,:,i),eread(i))
  end do
  close (ich)
  close (ich5)
  deallocate (dpop,pop)
  return
end subroutine wrpropens_pop

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! obtain an averaged IR spectrum for the ensemble
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine autoir(TMPCONF,imode,env)
  use iso_fortran_env,wp => real64
  use crest_data
  use iomod
  use utilities,only:boltz2
  implicit none

  type(systemdata) :: env

  integer,intent(in) :: TMPCONF
  integer :: imode

  integer :: i,j,k,l
  integer :: nall,nmodes,mnew,npop
  integer :: ich,ich5
  character(len=512) :: tmppath
  real(wp) :: dE
  real(wp) :: pthr
  logical :: ex
  integer  :: minl,minll(1),maxl

  real(wp),allocatable :: eread(:),pop(:),pop2(:),edum(:)
  real(wp),allocatable :: vibspec(:,:,:)

  real(wp),allocatable :: freq(:),tmpfreq(:),dum(:)
  real(wp),allocatable :: inten(:),tmpint(:)

  real(wp),parameter :: kcal = 627.5095_wp

  nall = TMPCONF
  pthr = env%pthr

  write (*,*)
  write (*,'(1x,a)') "Obtaining calculated vibspectra ..."

  select case (imode)
  case (2)
!--- get the free energies for better weights
    write (*,'(1x,a)') "Calculating populations from free energies ..."
    allocate (eread(TMPCONF),pop(TMPCONF))
    do i = 1,TMPCONF
      write (tmppath,'(a,i0,a,a)') 'TMPCONF',i,'/','xtb.out'
      call grepval(tmppath,'TOTAL FREE ENERGY',ex,eread(i))
    end do
  case default
!--- just energies
    write (*,'(1x,a)') "Calculating populations from total energies ..."
    allocate (eread(TMPCONF),pop(TMPCONF))
    do i = 1,TMPCONF
      write (tmppath,'(a,i0,a,a)') 'TMPCONF',i,'/','xtb.out'
      call grepval(tmppath,'total energy',ex,eread(i))
    end do
  end select
!--- convert to populations
  call boltz2(nall,eread,pop)
!      if(env%hardcutDFT .and. imode.ne.2)then !special sort mode
!         call cutDFTpop(env,pop,nall)
!         allocate(edum(nall))
!         edum(1:nall)=eread(1:nall)*pop(1:nall) !set energies of removed structures to 0
!         call boltz2(nall,edum,pop) !convert to new populations, only with used structures
!      endif
  maxl = maxloc(pop(:),1)

!--- short energy printout
  write (*,*)
  select case (imode)
  case (2)
    write (*,'(a)') '========================================================='
    write (*,'(a)') '============= free energies & populations  =============='
    write (*,'(a)') '========================================================='
    write (*,'('' structure    ΔG(kcal/mol)    Gtot(Eh)        weight'')')
  case default
    write (*,'(a)') '========================================================='
    write (*,'(a)') '============= total energies & populations  ============='
    write (*,'(a)') '========================================================='
    write (*,'('' structure    ΔE(kcal/mol)    Etot(Eh)        weight'')')
  end select
  do i = 1,nall
    dE = (eread(i)-eread(1))*kcal
    write (*,'(i5,6x,F10.2,4x,F14.6,F13.4)') i,dE,eread(i),pop(i)
  end do
  write (*,*)

  !--- if only populated structures are required.
  npop = 0
  do i = 1,TMPCONF
    if (pop(i) .ge. pthr.or.i .eq. maxl) npop = npop+1
  end do
  nall = npop
  write (*,'(1x,i0,a,i0,a)') nall, &
  & ' structures are above the population threshold of ',nint(pthr*100.0_wp),'%.'

!--- read-in spectra
  nmodes = env%nat*3
  allocate (vibspec(2,nmodes,nall),pop2(nall))  !vibspec(1,:,i)=frequencies of mol i
  !vibspec(2,:,i)=intensities of mol i

  open (newunit=ich5,file='autoir.pop')
  npop = 1
  do i = 1,TMPCONF
    if (pop(i) .ge. pthr.or.i .eq. maxl) then
      write (ich5,'(1x,i0,1x,f6.4)') i,pop(i)
      write (tmppath,'(a,i0,a,a)') 'TMPCONF',i,'/','vibspectrum'
      call rdvibs(tmppath,nmodes,vibspec(1,:,npop),vibspec(2,:,npop))
      pop2(npop) = pop(i)
      npop = npop+1
    end if
  end do
  close (ich5)

  write (*,'(1x,a,i0,a)',advance='no') 'Weighting ',nall,' vibspectra ...'
!--- weight spectra and sort frequencies
  mnew = (nmodes-6)*nall+6
  allocate (freq(mnew),inten(mnew))
!--- write to new arrays
  !--- translation first
  do i = 1,6
    freq(i) = 0.0_wp
    inten(i) = 0.0_wp
  end do
  !--- then the vibspectra
  l = 7
  do k = 1,nall
    do j = 7,nmodes
      freq(l) = vibspec(1,j,k)
      inten(l) = vibspec(2,j,k)*pop2(k) !scaled intensity
      l = l+1
    end do
  end do

!--- sort in ascending order
  allocate (tmpfreq(mnew),tmpint(mnew),dum(mnew-6))
  tmpfreq = 0.0_wp
  tmpint = 0.0_wp
  dum(:) = freq(7:mnew)
  do i = 7,mnew
    minll = minloc(dum)
    minl = minll(1)
    tmpfreq(i) = dum(minl)
    tmpint(i) = inten(minl+6)
    dum(minl) = 100000.0_wp
  end do
  freq = tmpfreq
  inten = tmpint
  deallocate (dum,tmpint,tmpfreq)

  write (*,'(1x,a)') "done."
  write (*,*)

!--- write new file
  open (file='crest.vibspectrum',newunit=ich)
  call write_tm_vibspectrum(ich,mnew,freq,inten)
  write (*,'(1x,a)') 'Written to file <crest.vibspectrum>'

  deallocate (inten,freq,pop2,vibspec,pop,eread)

  return
end subroutine autoir

!--- read vibspectrum file in TM format
subroutine rdvibs(fname,nmodes,freq,inten)
  use iso_fortran_env,wp => real64
  use crest_data
  use iomod
  implicit none

  character(len=*),intent(in) :: fname
  integer,intent(in)   :: nmodes
  real(wp),intent(out) :: freq(nmodes)    !frequencies
  real(wp),intent(out) :: inten(nmodes)   !intensities

  integer :: k
  integer :: ich,io,n
  character(len=256) :: atmp
  real(wp) :: floats(10)
  logical :: ex

  freq = 0.0_wp
  inten = 0.0_wp

  inquire (file=fname,exist=ex)
  if (.not.ex) return

  k = 1 !modes
  open (file=fname,newunit=ich)
  rdfile: do
    read (ich,'(a)',iostat=io) atmp
    if (io < 0) exit
    if (index(atmp,'$vibrational spectrum') .ne. 0) then
      rdblock: do
        read (ich,'(a)',iostat=io) atmp
        if (io < 0) exit rdfile
        if (index(atmp,'$end') .ne. 0) exit rdfile
        if (index(atmp,'#') .ne. 0) cycle rdblock !skip comment lines
        call readl(atmp,floats,n)
        freq(k) = floats(2)
        inten(k) = floats(3)
        k = k+1
      end do rdblock
    end if
  end do rdfile

  return
end subroutine rdvibs

subroutine write_tm_vibspectrum(ich,n3,freq,ir_int)
  use iso_fortran_env,wp => real64
  integer,intent(in)  :: ich ! file handle
  integer,intent(in)  :: n3
  real(wp),intent(in)  :: freq(n3)
  real(wp),intent(in)  :: ir_int(n3)
  integer  :: i
  real(wp) :: thr = 0.01_wp
  write (ich,'("$vibrational spectrum")')
  write (ich,'("#  mode     symmetry     wave number   IR intensity    selection rules")')
  write (ich,'("#                         cm**(-1)        (amu)          IR     RAMAN")')
  do i = 1,n3
    if (abs(freq(i)) .lt. thr) then
      write (ich,'(i6,9x,    f18.2,f16.5,7x," - ",5x," - ")') &
        i,freq(i),0.0_wp
    else
      write (ich,'(i6,8x,"a",f18.2,f16.5,7x,"YES",5x,"YES")') &
        i,freq(i),ir_int(i)
    end if
  end do
  write (ich,'("$end")')
end subroutine

