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

subroutine tauthead
  implicit none
  write (*,*) '       __________________________________________'
  write (*,*) '      |                                          |'
  write (*,*) '      |     automated tautomerization script     |'
  write (*,*) '      |__________________________________________|'
  write (*,*) ' Universitaet Bonn, MCTC'
  write (*,*) ' P.Pracht, Wed 28. Nov 13:11:52 CEST 2018'
  write (*,*)
  write (*,*) 'Cite as:'
  write (*,*) 'P.Pracht, R.Wilcken, A.Udvarhelyi, S.Rodde, S.Grimme'
  write (*,*) 'JCAMD, 2018, 32, 1139-1149.'
  write (*,*)
end subroutine tauthead

!--------------------------------------------------------------------------------------------
! Tautomerization workflow with GFNn-xTB
!--------------------------------------------------------------------------------------------
subroutine tautomerize(env,tim)
  use crest_parameters
  use crest_data
  use iomod
  use strucrd,only:coord2xyz
  use utilities
  implicit none
  type(systemdata) :: env
  type(timer)      :: tim
  type(protobj)    :: taut

  character(len=32)  :: dirn
  character(len=64)  :: tautname
  character(len=256) :: thispath
  character(len=256) :: filename
  character(len=128) :: inpnam,outnam
  character(len=128) :: dummy

  integer :: ich,i
  integer :: natp,nallout,refchrg

  logical :: ex

!--- printout & clean directory
  call tautclean
  call tauthead

  if (.not.allocated(env%ptb%atmap)) allocate (env%ptb%atmap(env%nat))
  if (.not.env%ptb%strictPDT.and..not.env%ptb%fixPDT) then
!--- sort the input file (H atoms to the bottom)
    call htothebottom('coord',env%chrg,env%nat,env%ptb%atmap)
  else
!--- or sort AND apply heavy atom bond constraints
    call PDT_constraints(env)
  end if

!--- get some settings
  call getcwd(thispath)
  dirn = 'PROT'
  tautname = 'tautomerize_0.xyz'
  taut = env%ptb
  natp = env%nat     !backup Nat
  refchrg = env%chrg

  !-- for the regular mode we have to start with a modified protonation cycle
  if (.not.taut%deprotprot) then
!--- do Protonation
    write (dummy,'(a,i0,a,i0)') 'P R O T O N A T I O N   C Y C L E     ', &
    & 1,' of ',taut%iter+1
    call largehead(trim(dummy))
    call protsmall(env,taut,tim)
    inquire (file='coordprot.0',exist=ex)
    if (.not.ex) return
    call checkname_xyz('tautomerize',inpnam,outnam)
    call rename('protonated.xyz',outnam)
!--- do Deprotonation
    write (dummy,'(a,i0,a,i0)') 'D E P R O T O N A T I O N   C Y C L E     ', &
    & 1,' of ',taut%iter+1
    call largehead(trim(dummy))
    call checkname_xyz('tautomerize',inpnam,outnam)
    call deprotens(inpnam,env,taut,tim)
    call rename('deprotonated.xyz',outnam)
  else
    !-- for the reversed deprotonation/protonation mode we can use the iteration loop
    !   below instead, but we have to set up some minor things
    call checkname_xyz('tautomerize',inpnam,outnam)
    call coord2xyz('coord',outnam) !get input in the xyz format
    taut%iter = taut%iter+1
  end if

  !taut%ewin=taut%ewin/2.0d0
!--- further iterations
  do i = 1,taut%iter
    call tautclean2
    if (.not.taut%deprotprot) then
      !--- do Protonation
      write (dummy,'(a,i0,a,i0)') 'P R O T O N A T I O N   C Y C L E     ', &
      & i+1,' of ',taut%iter+1
      call largehead(trim(dummy))
      call checkname_xyz('tautomerize',inpnam,outnam)
      call protens(inpnam,env,taut,tim)
      call rename('protonated.xyz',outnam)
      !--- do Deprotonation
      write (dummy,'(a,i0,a,i0)') 'D E P R O T O N A T I O N   C Y C L E     ', &
      & i+1,' of ',taut%iter+1
      call largehead(trim(dummy))

      call checkname_xyz('tautomerize',inpnam,outnam)
      call deprotens(inpnam,env,taut,tim)
      call rename('deprotonated.xyz',outnam)
    else
      !--- do Deprotonation
      write (dummy,'(a,i0,a,i0)') 'D E P R O T O N A T I O N   C Y C L E     ', &
      & i,' of ',taut%iter
      call largehead(trim(dummy))

      call checkname_xyz('tautomerize',inpnam,outnam)
      call deprotens(inpnam,env,taut,tim)
      call rename('deprotonated.xyz',outnam)
      !--- do Protonation
      write (dummy,'(a,i0,a,i0)') 'P R O T O N A T I O N   C Y C L E     ', &
      & i,' of ',taut%iter
      call largehead(trim(dummy))
      call checkname_xyz('tautomerize',inpnam,outnam)
      call protens(inpnam,env,taut,tim)
      call rename('protonated.xyz',outnam)
    end if
  end do

!--- clean all temporary ensemble files
  call tautclean2

!--- reverse (Deprotonation - Protonation)
  !call tautclean

  !call coord2xyz('coord',"struc_0.xyz")
  !call deprotens("struc_0.xyz",env,taut,tim)

!--- reset data for main dir
  env%chrg = refchrg
  if (env%chrg .eq. 0) then
    call remove('.CHRG')
  else
    open (newunit=ich,file='.CHRG')
    write (ich,*) env%chrg
    close (ich)
  end if
  env%nat = natp    !reset Nat

!--- get the new charge and set up the calculations
  !call printtaut
  call largehead('T A U T O M E R I Z E')
  call tim%start(2,'multilevel OPT')

  call smallhead('Final Geometry Optimization')
  call checkname_xyz('tautomerize',inpnam,outnam)
  call MDopt_para(env,inpnam,0)
  filename = trim(thispath)//'/'//trim(outnam)
  call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
  call rmrf('OPTIM')
  call sort_ens(taut,outnam,.false.)
  call remaining_in(outnam,taut%ewin,nallout) !--- remaining number of structures

  call cosort(outnam,'tautomers.xyz',.false.,.false.)
  call sort_ens(taut,'tautomers.xyz',.true.)
  call tim%stop(2)

!>--- (optional) post-processing
  if (env%relax) then
    call relaxensemble('tautomers.xyz',env,tim)
  end if

  if (env%outputsdf) then
    call new_wrsdfens(env,'tautomers.xyz','tautomers.sdf',.true.)
  end if

end subroutine tautomerize

!--------------------------------------------------------------------------------------------
! small verion of the protonate routine
!--------------------------------------------------------------------------------------------
subroutine protsmall(env,prot,tim)
  use crest_parameters
  use iomod
  use crest_data
  use strucrd,only:coord2xyz
  use utilities
  implicit none
  type(systemdata) :: env
  type(protobj)    :: prot
  type(timer)      :: tim

  integer :: ich,natp,nallout

  character(len=64)  :: protname
  character(len=256) :: thispath
  character(len=256) :: filename
  character(len=128) :: inpnam,outnam

  logical :: ex
  integer :: refchrg

  !call printprotcy

  call getcwd(thispath)
!--- do the xTB calculation for the LMOs
  call tim%start(1,'LMO calc.')
  call xtblmo(env)
  call tim%stop(1)
  inquire (file='coordprot.0',exist=ex)
  if (.not.ex) then
    write (*,*)
    write (*,*) '***Warning***'
    write (*,*) 'No "coordprot.0" file was written, it is possible that'
    write (*,*) 'there are no suitable LP- or π-centers in the molecule.'
    write (*,*) 'Hence the procedure could not be automatized. (sorry)'
    write (*,*) '***Warning***'
    return
  end if

!--- get the new charge and set up the calculations
  natp = env%nat+1
  prot%newchrg = env%chrg+1  !increase chrg by one
  refchrg = env%chrg
  env%chrg = env%chrg+1

  protname = 'protonate_0.xyz'
  call tim%start(2,'multilevel OPT')
  open (newunit=ich,file='.CHRG')
  write (ich,*) prot%newchrg     !new charge written here
  close (ich)
  write (*,'(''-----------------------'')')
  write (*,'(''Multilevel Optimization'')')
  write (*,'(''-----------------------'')')

  call coord2xyz('coordprot.0',trim(protname))
  call appendto('xtbscreen.xyz',protname)
  env%nat = natp

  call smallhead('1. crude pre-optimization')
  call checkname_xyz('protonate',inpnam,outnam)
  call MDopt_para(env,protname,1)
  filename = trim(thispath)//'/'//trim(outnam)
  call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
  call rmrf('OPTIM')
  prot%ewin = prot%ewin*2.0d0
  call sort_ens(prot,outnam,.false.)
  call remaining_in(outnam,prot%ewin,nallout) !--- remaining number of structures
  write (*,*)

  call smallhead('2. loose optimization')
  call checkname_xyz('protonate',inpnam,outnam)
  call MDopt_para(env,inpnam,2)
  filename = trim(thispath)//'/'//trim(outnam)
  call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
  call rmrf('OPTIM')
  prot%ewin = prot%ewin/2.0d0
  call sort_ens(prot,outnam,.false.)
  call remaining_in(outnam,prot%ewin,nallout) !--- remaining number of structures
  call cosort(outnam,'protonated.xyz',.false.,.true.)

  call tim%stop(2)

  env%nat = natp-1   !reset Nat
  env%chrg = refchrg   !reset chrg
  return
end subroutine protsmall

!--------------------------------------------------------------------------------------------
! deprotonate every structure of an ensemble file
!--------------------------------------------------------------------------------------------
subroutine deprotens(ens,env,prot,tim)
  use crest_parameters
  use iomod
  use crest_data
  use strucrd,only:rdensembleparam,rdensemble,i2e
  use utilities
  implicit none
  type(systemdata) :: env
  type(protobj)    :: prot
  type(timer)      :: tim

  integer :: ich,nallout
  integer :: nat,nall
  integer :: i,j,k,l

  character(len=*)   :: ens
  character(len=64)  :: protname
  character(len=256) :: thispath
  character(len=256) :: filename
  character(len=128) :: inpnam,outnam

  real(wp),allocatable :: xyz(:,:,:),eread(:)
  integer,allocatable  :: at(:)
  integer :: refchrg

  !call printdeprotcy

!--- settings
  call getcwd(thispath)
  protname = 'deprotonate_0.xyz'
  prot%newchrg = prot%newchrg-1
  refchrg = env%chrg

!--- read the file
  call rdensembleparam(ens,nat,nall)
  allocate (xyz(3,nat,nall),at(nat),eread(nall))
  call rdensemble(ens,nat,nall,at,xyz,eread)

  open (newunit=ich,file=protname)

  l = 0
  do k = 1,nall
    do i = 1,nat
      if (at(i) .ne. 1) cycle
      write (ich,'(1x,i6)') nat-1
      write (ich,*)
      do j = 1,nat
        if (i .eq. j) then
          cycle
        else
          write (ich,'(a2,2x,3F16.10)') i2e(at(j),'nc'),xyz(1:3,j,k)
        end if
      end do
      l = l+1
    end do
  end do

  close (ich)
  deallocate (eread,at,xyz)
  if (l .lt. 1) then
    error stop 'no new structures written in deprotonation setup!'
  end if

!--- get the new charge and set up the calculations
  call tim%start(2,'multilevel OPT')
  open (newunit=ich,file='.CHRG')
  write (ich,*) prot%newchrg     !new charge written here
  close (ich)
  write (*,'(''-----------------------'')')
  write (*,'(''Multilevel Optimization'')')
  write (*,'(''-----------------------'')')

!--- update Nat for optimization
  env%nat = nat-1
  env%chrg = prot%newchrg

  call smallhead('1. crude pre-optimization')
  call checkname_xyz('deprotonate',inpnam,outnam)
  call MDopt_para(env,protname,1)
  filename = trim(thispath)//'/'//trim(outnam)
  call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
  call rmrf('OPTIM')
  prot%ewin = prot%ewin*2.0d0
  call sort_ens(prot,outnam,.false.)
  call remaining_in(outnam,prot%ewin,nallout) !--- remaining number of structures
  write (*,*)

  call smallhead('2. loose optimization')
  call checkname_xyz('deprotonate',inpnam,outnam)
  call MDopt_para(env,inpnam,2)
  filename = trim(thispath)//'/'//trim(outnam)
  call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
  call rmrf('OPTIM')
  prot%ewin = prot%ewin/2.0d0
  call sort_ens(prot,outnam,.false.)
  call remaining_in(outnam,prot%ewin,nallout) !--- remaining number of structures
  call cosort(outnam,'deprotonated.xyz',.false.,.true.)

  call tim%stop(2)

!--- reset Nat
  env%nat = nat
  env%chrg = refchrg

end subroutine deprotens

!--------------------------------------------------------------------------------------------
! protonate routine to be used on an ensemble (quite lengthy and a lot of bookkeeping)
!--------------------------------------------------------------------------------------------
subroutine protens(ens,env,prot,tim)
  use crest_parameters
  use iomod
  use crest_data
  use strucrd,only:coord2xyz,wrc0,rdensembleparam,rdensemble
  use utilities
  implicit none
  type(systemdata) :: env
  type(protobj)    :: prot
  type(timer)      :: tim

  integer :: ich,natp,nallout
  integer :: nat,nall
  integer :: i,k,r
  integer :: vz,io,refchrg

  real(wp) :: percent

  character(len=*)   :: ens
  character(len=256) :: thispath,tmppath
  character(len=256) :: filename
  character(len=128) :: inpnam,outnam
  character(len=:),allocatable :: jobcall
  
  logical :: niceprint

  real(wp),allocatable :: xyz(:,:,:),eread(:)
  integer,allocatable  :: at(:)
  character(len=*),parameter :: dirn='PROT'
  !call printprotcy

!--- some settings
  call getcwd(thispath)
  niceprint = env%niceprint
  refchrg = env%chrg

  r = makedir(dirn)

!--- read the file
  call rdensembleparam(ens,nat,nall)
  allocate (xyz(3,nat,nall),at(nat),eread(nall))
  call rdensemble(ens,nat,nall,at,xyz,eread)
  xyz = xyz/bohr  !--- Ang to bohr
  natp = nat+1

!--- change dir
  call chdir(dirn)
  call getcwd(tmppath)
!--- make new dirs
  do i = 1,nall
    write (filename,'(a,i0)') dirn,i
    r = makedir(trim(filename))
    call chdir(trim(filename))
    call wrc0('coord',nat,at,xyz(:,:,i))
    open (newunit=ich,file='.CHRG')
    write (ich,*) prot%newchrg           !--- not yet updated; for LMO calculation
    close (ich)
    call chdir(tmppath)
  end do
  deallocate (eread,at,xyz)

!--- thread stuff
  if (env%autothreads) then
    call ompautoset(env%threads,7,env%omp,env%MAXRUN,nall) !set the global OMP/MKL variables for the xtb jobs
  end if

!--- creating the job
  jobcall = trim(env%ProgName)
  jobcall = trim(jobcall)//' '//'coord'
  jobcall = trim(jobcall)//' '//trim(env%gfnver)
  jobcall = trim(jobcall)//' --sp --lmo '//trim(env%solv)
  jobcall = trim(jobcall)//' > xtb.out 2>/dev/null'
  

!--- calculation loop for LMOs
  call tim%start(1,'LMO calc.')
  write (*,'(a,a,a)') 'Calculating LMOs for all structures in file <',trim(ens),'>'

  k = 0 !counting the finished jobs
  if (niceprint) then
    call printprogbar(0.0_wp)
  end if

!!$omp parallel &
!!$omp shared( vz,jobcall,nall,percent,k,niceprint ) &
!!$omp private(filename)
!!$omp single
  do i = 1,nall
    vz = i
!    !$omp task firstprivate( vz ) private( io )
    call initsignal()
    write (filename,'(a,i0)') dirn,vz
    call command('cd '//trim(filename)//' && '//trim(jobcall),io)
!   !$omp critical
    k = k+1
    if (niceprint) then
      percent = float(k)/float(nall)*100.0d0
      call printprogbar(percent)
    else
      write (stdout,'(1x,i0)',advance='no') k
      flush (stdout)
    end if
!    !$omp end critical
!    !$omp end task
  end do
!!$omp taskwait
!!$omp end single
!!$omp end parallel
  !--- this is a test for BASF
  if (env%threads > 8) then
    call sleep(5)
  end if
  write (*,*)
  write (*,'(a)',advance='no') 'Collecting generated protomers ...'
  jobcall = trim(tmppath)//'/'//'protomers.xyz'
  do i = 1,nall
    write (filename,'(a,i0)') dirn,i
    call chdir(trim(filename))
    call coord2xyz('coordprot.0','struc_0.xyz')
    call appendto('struc_0.xyz',jobcall)
    call appendto('xtbscreen.xyz',jobcall)
    call chdir(tmppath)
  end do
  write (*,*) 'done.'
  call tim%stop(1)

!--- clean up the sub-dirs
  call rmrfw('PROT')

!--- optimize
  call tim%start(2,'multilevel OPT')
  prot%newchrg = prot%newchrg+1
  open (newunit=ich,file='.CHRG')
  write (ich,*) prot%newchrg     !new charge written here
  close (ich)
  write (*,*)
  write (*,'(''-----------------------'')')
  write (*,'(''Multilevel Optimization'')')
  write (*,'(''-----------------------'')')

!--- update nat
  env%nat = nat+1
  env%chrg = prot%newchrg

  call smallhead('1. crude pre-optimization')
  call checkname_xyz('protonate',inpnam,outnam)
  call MDopt_para(env,'protomers.xyz',1)
  filename = trim(tmppath)//'/'//trim(outnam)
  call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
  call rmrf('OPTIM')
  prot%ewin = prot%ewin*2.0d0
  call sort_ens(prot,outnam,.false.)
  call remaining_in(outnam,prot%ewin,nallout) !--- remaining number of structures
  write (*,*)

  call smallhead('2. loose optimization')
  call checkname_xyz('protonate',inpnam,outnam)
  call MDopt_para(env,inpnam,2)
  filename = trim(tmppath)//'/'//trim(outnam)
  call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
  call rmrf('OPTIM')
  prot%ewin = prot%ewin/2.0d0
  call sort_ens(prot,outnam,.false.)
  call remaining_in(outnam,prot%ewin,nallout) !--- remaining number of structures
  call cosort(outnam,'protonated.xyz',.false.,.true.)

  call tim%stop(2)

!--- reset Nat
  env%nat = nat
  env%chrg = refchrg

!--- change back to original dir and copy the file with optimized protomers
  call chdir(thispath)
  jobcall = trim(tmppath)//'/'//'protonated.xyz' 
  call rename(trim(jobcall),'protonated.xyz')
  call rmrf(dirn)
  return
end subroutine protens

!==============================================================================!
! build the blacklist to decide which atoms are valid deprotonation candidates
!==============================================================================!
subroutine tautomerize_blacklist(env,fname,nat,atlist)
  use crest_parameters
  use iomod
  use crest_data
  use strucrd,only:rdcoord
  use utilities
  implicit none
  type(systemdata) :: env
  integer :: nat
  character(len=*) :: atlist
  integer,allocatable :: unconstrained(:)
  character(len=*) :: fname
  integer :: ncon,i,j,k
  integer,allocatable :: at(:)
  real(wp),allocatable :: xyz(:,:)
  real(wp),allocatable :: dist(:,:)

  write (*,'(1x,a,a)') 'Input list of atoms: ',trim(atlist)

  allocate (unconstrained(nat))
  call parse_atlist(atlist,ncon,nat,unconstrained) !"unconstrained" contains all the selected atoms

  if (ncon .lt. 1) then
    deallocate (unconstrained)
    return
  end if

  if (.not.allocated(env%ptb%blacklist)) then
    allocate (env%ptb%blacklist(nat))
    env%ptb%blacklist = .false. !none of the atoms is initially blacklisted
  end if

  allocate (xyz(3,nat),dist(nat,nat))
  call rdcoord(fname,nat,at,xyz)

  !setup distmat
  do i = 1,nat
    dist(i,i) = 9999.9_wp
    do j = 1,i-1
      dist(i,j) = (xyz(1,i)-xyz(1,j))**2+ &
     &            (xyz(2,i)-xyz(2,j))**2+ &
     &            (xyz(3,i)-xyz(3,j))**2
      dist(i,j) = sqrt(dist(i,j))
      dist(j,i) = dist(i,j)
    end do
  end do

  !loop over the atoms and decide which to blacklist
  do i = 1,nat
    if (unconstrained(i) == 1) then
      if (at(i) == 1) then
        k = minloc(dist(i,:),1)
        env%ptb%blacklist(k) = .true.
      else
        env%ptb%blacklist(i) = .true.
      end if
    end if
  end do

  deallocate (dist,xyz,unconstrained)
  return
end subroutine tautomerize_blacklist

!========================================================!
! check if the input atom order is kompatible with
! the tautomerization tools, i.e., if all hydrogen
! atoms are at the bottom of the list.
!========================================================!
function taut_check_atomorder(n,at) result(bool)
  implicit none
  integer :: n
  integer :: at(n)
  logical :: bool
  integer :: i
  bool = .true.
  do i = 2,n
    if (at(i-1) == 1) then
      if (at(i) .ne. 1) then
        bool = .false.
        exit
      end if
    end if
  end do
  return
end function taut_check_atomorder

!============================================================================================!
! Tautomerization workflow with GFNn-xTB (*extended version)
!============================================================================================!
subroutine tautomerize_ext(ensemb,env,tim)
  use crest_parameters
  use crest_data
  use iomod
  use strucrd
  use utilities
  implicit none
  character(len=*) :: ensemb
  type(systemdata) :: env
  type(timer)      :: tim
  type(protobj)    :: taut

  character(len=32)  :: dirn
  character(len=64)  :: tautname
  character(len=256) :: thispath
  character(len=256) :: filename
  character(len=128) :: inpnam,outnam
  character(len=128) :: dummy,atmp
  character(len=:),allocatable :: btmp

  type(ensemble) :: efile
  integer,allocatable :: slist(:)
  integer :: smax,s,sc
  real(wp),allocatable :: c0(:,:)

  integer :: ich,i,io
  integer :: natp,nallout,refchrg
  integer,allocatable :: atmaps(:,:)
  logical :: ex

!--- printout & clean directory
  call tautclean
  call tauthead

!--- get some settings
  call getcwd(thispath)
  dirn = 'PROT'
  tautname = 'tautomerize_0.xyz'
  taut = env%ptb
  natp = env%nat     !backup Nat
  refchrg = env%chrg

!--- read the given ensemble
  call efile%open(ensemb)
  allocate (c0(3,efile%nat))
  env%rednat = efile%nat
  env%nat = efile%nat

  if (taut%alldivers) then
    smax = efile%nall
    allocate (slist(smax))
    do i = 1,smax
      slist(i) = i
    end do
  else
    smax = taut%divers
    allocate (slist(smax))
    !ANALYZE DIVERSITY AND SELECT STRUCTURES
    do i = 1,smax   !<-- temporary workaround
      slist(i) = i
    end do
  end if
  allocate (atmaps(smax,natp))

!---------------------------
  SLOOP: do s = 1,smax
!---------------------------
    sc = slist(s) !sc is the current structure
    write (atmp,'(a,i0)') 'TP',sc
    btmp = trim(atmp)
    io = makedir(btmp)

    call chdir(btmp)
    c0(:,:) = efile%xyz(:,:,sc)/bohr
    !---  write .CHRG and .UHF files
    env%chrg = refchrg
    if (env%chrg .ne. 0) then
      open (newunit=ich,file='.CHRG')
      write (ich,*) env%chrg
      close (ich)
    end if
    if (env%uhf .ne. 0) then
      open (newunit=ich,file='.UHF')
      write (ich,*) env%uhf
      close (ich)
    end if
    call wrc0('coord',efile%nat,efile%at,c0)

!--- sort the input file (H atoms to the bottom)
    call htothebottom('coord',env%chrg,natp,atmaps(s,1:natp))

    !-- for the regular mode we have to start with a modified protonation cycle
    if (.not.taut%deprotprot) then
!--- do Protonation
      write (dummy,'(a,i0,a,i0)') 'P R O T O N A T I O N   C Y C L E     ', &
      & 1,' of ',taut%iter+1
      call largehead(trim(dummy))
      call protsmall(env,taut,tim)
      inquire (file='coordprot.0',exist=ex)
      if (.not.ex) return
      call checkname_xyz('tautomerize',inpnam,outnam)
      call rename('protonated.xyz',outnam)
!--- do Deprotonation
      write (dummy,'(a,i0,a,i0)') 'D E P R O T O N A T I O N   C Y C L E     ', &
      & 1,' of ',taut%iter+1
      call largehead(trim(dummy))
      call checkname_xyz('tautomerize',inpnam,outnam)
      call deprotens(inpnam,env,taut,tim)
      call rename('deprotonated.xyz',outnam)
!--- Relax structures by performing a small conformational search
      env%nat = natp
      call relaxensemble(outnam,env,tim)
      call rename('relax.'//outnam,outnam)
    else
      !-- for the reversed deprotonation/protonation mode we can use the iteration loop
      !   below instead, but we have to set up some minor things
      call checkname_xyz('tautomerize',inpnam,outnam)
      call coord2xyz('coord',outnam) !get input in the xyz format
      taut%iter = taut%iter+1
    end if

    !taut%ewin=taut%ewin/2.0d0
!--- further iterations
    do i = 1,taut%iter
      call tautclean2
      if (.not.taut%deprotprot) then
        !--- do Protonation
        write (dummy,'(a,i0,a,i0)') 'P R O T O N A T I O N   C Y C L E     ', &
        & i+1,' of ',taut%iter+1
        call largehead(trim(dummy))
        call checkname_xyz('tautomerize',inpnam,outnam)
        call protens(inpnam,env,taut,tim)
        call rename('protonated.xyz',outnam)
        !--- do Deprotonation
        write (dummy,'(a,i0,a,i0)') 'D E P R O T O N A T I O N   C Y C L E     ', &
        & i+1,' of ',taut%iter+1
        call largehead(trim(dummy))

        call checkname_xyz('tautomerize',inpnam,outnam)
        call deprotens(inpnam,env,taut,tim)
        call rename('deprotonated.xyz',outnam)
      else
        !--- do Deprotonation
        write (dummy,'(a,i0,a,i0)') 'D E P R O T O N A T I O N   C Y C L E     ', &
        & i,' of ',taut%iter
        call largehead(trim(dummy))

        call checkname_xyz('tautomerize',inpnam,outnam)
        call deprotens(inpnam,env,taut,tim)
        call rename('deprotonated.xyz',outnam)
        !--- do Protonation
        write (dummy,'(a,i0,a,i0)') 'P R O T O N A T I O N   C Y C L E     ', &
        & i,' of ',taut%iter
        call largehead(trim(dummy))
        call checkname_xyz('tautomerize',inpnam,outnam)
        call protens(inpnam,env,taut,tim)
        call rename('protonated.xyz',outnam)
      end if
!--- relax structures after each cycle
      env%nat = natp
      call relaxensemble(outnam,env,tim)
      call rename('relax.'//outnam,outnam)
    end do

!--- clean all temporary ensemble files
    call tautclean2

    call chdir(thispath)

    btmp = trim(btmp)//'/'//trim(outnam)
    call appendto(btmp,'collected.xyz')
!-------------------------
  end do SLOOP
!-------------------------
!      return

!--- reset data for main dir
  env%chrg = refchrg
  if (env%chrg .eq. 0) then
    call remove('.CHRG')
  else
    open (newunit=ich,file='.CHRG')
    write (ich,*) env%chrg
    close (ich)
  end if
  env%nat = natp    !reset Nat

!--- get the new charge and set up the calculations
  !call printtaut
  call largehead('T A U T O M E R I Z E')
  call tim%start(2,'multilevel OPT')

  call smallhead('Final Geometry Optimization')
  call checkname_xyz('tautomerize',inpnam,outnam)

  call rename('collected.xyz',inpnam)

  call MDopt_para(env,inpnam,0)
  filename = trim(thispath)//'/'//trim(outnam)
  call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
  call rmrf('OPTIM')
  call sort_ens(taut,outnam,.false.)
  call remaining_in(outnam,taut%ewin,nallout) !--- remaining number of structures

  call cosort(outnam,'tautomers.xyz',.false.,.false.)
  call sort_ens(taut,'tautomers.xyz',.true.)
  call tim%stop(2)

end subroutine tautomerize_ext

!--------------------------------------------------------------------------------------------!
! first sort file and then
! check if strict mode is active to define bond constraints (only for the heavy atoms)
!--------------------------------------------------------------------------------------------!
subroutine PDT_constraints(env)
  use crest_parameters
  use crest_data
  use iomod
  use zdata
  use utilities
  implicit none
  type(systemdata) :: env
  type(zmolecule) :: zmol
  !-- a default for the force constant (quite strong already)

  logical,parameter :: vverbose = .false.
  integer :: i,h,nh

  !-- sort all H atoms to the end of the file
  call htothebottom('coord',env%chrg,env%nat,env%ptb%atmap)

  !-- decide between atom fixing and bond constraints
  if (env%ptb%strictPDT.or.env%ptb%fixPDT) then

    !-- if -for some reason- there is a constraint already specified, remove it
    if (env%cts%used) then
      call env%cts%deallocate()
    end if

    env%cts%used = .true.
    !-- if heavy-atom bond constraints shall be specified
    if (.not.env%ptb%fixPDT) then

      write (*,*) 'Strict mode active. Heavy atom bond constraints will be applied.'
      write (*,'(1x,a,f8.4,a)') 'Selected force constant:',env%forceconst,' Eh'

      !-- analyse sorted coord and write bondlength constraint
      call autoHeavyConstraint('coord',env%forceconst)

      !-- read (only) the bondlength file as constraint
      call read_constrainbuffer('bondlengths',env%cts)

      if (vverbose) then
        write (*,*) 'The following bonds are affected:'
        do i = 1,env%cts%ndim
          if (trim(env%cts%sett(i)) .ne. '') then
            write (*,'(''>'',1x,a)') trim(env%cts%sett(i))
          end if
        end do
      end if

      !--- otherwise do exact atom fixing (of heavy atoms only)
    else
      write (*,*) 'Very strict mode active. Heavy atom positions will be constrained.'
      write (*,'(1x,a,f8.4,a)') 'Selected force constant:',env%forceconst,' Eh'

      call simpletopo_file('coord',zmol,.false.,.false.,'')
      h = zmol%hydrogen() !-- count hydrogen
      nh = zmol%nat-h
      call fix_first_X_atoms(nh,env%forceconst,'fixpositions')
      if (env%ptb%strictPDT) then
        call getbmat(zmol,3,env%forceconst)
        call appendto('bondlengths','fixpositions')
      end if
      call zmol%deallocate()

      call read_constrainbuffer('fixpositions',env%cts)
    end if
  end if

  return
end subroutine PDT_constraints
