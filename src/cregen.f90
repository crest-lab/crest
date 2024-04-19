!===============================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2024 Philipp Pracht
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
!===============================================================================!

!=========================================================================================!
!=========================================================================================!
!> CREGEN is the universal ensemble sorting routine of CREST.
!> This is a rewrite of the original routines since the old ones
!> got a bit messy over time.
!> the quickset variable can be used for some special runtypes:
!>   quickset:  2 - do symmetry analysis
!>              3 - switch off equivalency analysis
!>              6 - energy sorting only
!>              9 - no sorting, only check groups
!>             12 - no topology check, turn ewin to infty
!>             13 - no topology check, ewin and rmsd checking (msreact settings)
!=========================================================================================!
!=========================================================================================!

module cregen_interface
!*******************************************************
!* module to load an interface to the newcregen routine
!* mandatory to handle the optional input arguments
!*******************************************************
  implicit none
  interface
    subroutine newcregen(env,quickset,infile)
      use crest_parameters
      use crest_data
      use crest_restartlog
      use strucrd
      implicit none
      type(systemdata),intent(inout) :: env
      integer,intent(in),optional :: quickset
      character(len=*),intent(in),optional :: infile
    end subroutine newcregen
  end interface
end module cregen_interface

subroutine newcregen(env,quickset,infile)
!****************************
!* The main CREGEN routine
!****************************
  use crest_parameters
  use crest_data
  use crest_restartlog
  use strucrd
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env   !> MAIN STORAGE OS SYSTEM DATA
  integer,intent(in),optional :: quickset !> quick access to predefined CREGEN modes
  character(len=*),intent(in),optional :: infile
  !> LOCAL
  integer :: simpleset
  character(len=258) :: fname  !> input file
  character(len=258) :: oname  !> sorted output file
  character(len=258) :: cname  !> unique structure file
!>--- ensemble arguments
  integer :: nat                      !> number of atoms
  integer :: nall                     !> number of structures
  integer,allocatable  :: at(:)       !> atom numbers
  real(wp),allocatable :: xyz(:,:,:)  !> Cartesian coordinates
  character(len=128),allocatable :: comments(:)
  character(len=128),allocatable :: comref(:)
  real(wp),allocatable :: er(:)       !> energies
!>--- dummy ensemble arguments
  integer :: nallref
  integer :: nallnew
  real(wp),allocatable :: xyzref(:,:,:)
!>--- sorting arguments
  integer,allocatable :: gref(:),group(:)
  integer :: ng
  integer,allocatable :: degen(:,:)

!>--- float data
  real(wp) :: ewin,rthr,bthr,pthr,ethr,athr
  real(wp) :: T,couthr

!>--- boolean data
  logical :: checkbroken
  logical :: topocheck
  logical :: checkez
  logical :: sortE
  logical :: sortRMSD
  logical :: sortRMSD2
  logical :: newfile
  logical :: repairord
  logical :: conffile
  logical :: bonusfiles
  logical :: anal
  logical :: saveelow = .true.
  logical :: userinput

!>--- printout directions
  integer :: prch  !> the main printout channel
  logical :: pr1,pr2,pr3,pr4

!>--- restart skip & tracking
  if (trackrestart(env)) return

!====================================================================!
!>  S E T T I N G S
!====================================================================!
  if (present(quickset)) then
    simpleset = quickset
  else
    simpleset = 0
  end if

!>-- determine filenames and output channel
  if (present(infile)) then
    fname = trim(infile)
    userinput = .true.
  else
    fname = trim(env%ensemblename)
    userinput = .false.
  end if
  call cregen_files(env,fname,oname,cname,simpleset,userinput,prch)

!>-- determine which printouts are required
  call cregen_prout(env,simpleset,pr1,pr2,pr3,pr4)

!>-- determine which subroutines are required
  call cregen_director(env,simpleset,checkbroken,sortE,sortRMSD,sortRMSD2, &
  &  repairord,newfile,conffile,bonusfiles,anal,topocheck,checkez)

!>--- DATA SECTION
  call cregen_filldata1(env,ewin,rthr,ethr,bthr,athr,pthr,T,couthr)
  call cregen_filldata2(simpleset,ewin)

!>--- setting the threads for OMP parallel usage
  call cregen_setthreads(prch,env,.false.)

!=====================================================================!
!>  E N S E M B L E   P R O C E S S I N G
!=====================================================================!

!>--- read in the ensemble parameters
  call rdensembleparam(fname,nat,nallref)

!>--- print a summary about the ensemble and thresholds
  if (pr1) call cregen_pr1(prch,env,nat,nallref,rthr,bthr,pthr,ewin)

!>--- allocate space and read in the ensemble
  allocate (at(nat),comments(nallref),xyz(3,nat,nallref))
  call rdensemble(fname,nat,nallref,at,xyz,comments)
!>--- track ensemble for restart
  call trackensemble(fname,nat,nallref,at,xyz,comments)

!>--- check if the ensemble contains broken structures? i.e., fusion or dissociation
  if (checkbroken) then
    call discardbroken(prch,env,topocheck,nat,nallref,at,xyz,comments,nall)
!>--- if structures were discarded, resize xyz
    if (nall .lt. nallref) then
      xyzref = xyz(:,:,1:nall)
      call move_alloc(xyzref,xyz)
      comref = comments(1:nall)
      call move_alloc(comref,comments)
    end if
  else
    nall = nallref
  end if

!>--- compare neighbourlists to sort out chemically transformed structures
  if (topocheck) then
    call cregen_topocheck(prch,env,checkez,nat,nall,at,xyz,comments,nallnew)
!>--- if structures were discarded, resize xyz
    if (nallnew .lt. nall) then
!>-- special fallback if all are discared
      if (nallnew == 0) then
        call rdcoord('coord',nat,at,xyz(:,:,1))
        xyz = xyz*bohr
        write (comments(1),'(f18.8)') env%elowest
        nallnew = 1
      end if
      nall = nallnew
      xyzref = xyz(:,:,1:nall)
      call move_alloc(xyzref,xyz)
      comref = comments(1:nall)
      call move_alloc(comref,comments)
    end if
  end if
  if (topocheck.or.checkbroken) then
    write (prch,'('' number of reliable points      :'',i6)') nall
  end if

!>--- sort the ensemble by its energies and make a cut (EWIN)
  if (sortE) then
    call cregen_esort(prch,nat,nall,xyz,comments,nallnew,ewin)
    !>--- if structures were discarded, resize xyz
    if (nallnew .lt. nall) then
      nall = nallnew
      xyzref = xyz(:,:,1:nall)
      call move_alloc(xyzref,xyz)
      comref = comments(1:nall)
      call move_alloc(comref,comments)
    end if
  end if

!>--- do the rotational constants and RMSD check
  if (sortRMSD) then
    allocate (group(0:nall))
    call cregen_CRE(prch,env,nat,nall,at,xyz,comments,nallnew,group)
!>--- if structures were discarded, resize xyz
    if (nallnew .lt. nall) then
      nall = nallnew
      xyzref = xyz(:,:,1:nall)
      call move_alloc(xyzref,xyz)
      comref = comments(1:nall)
      call move_alloc(comref,comments)
      allocate (gref(0:nallnew))
      gref(0:nallnew) = group(0:nallnew)
      call move_alloc(gref,group)
      nall = nallnew
    end if
!>--- repair the order
    if (repairord) then
      call cregen_repairorder(nat,nall,xyz,comments,group)
    end if
!>--- get group info to degen
    ng = group(0)
    allocate (degen(3,ng))
    call cregen_groupinfo(nall,ng,group,degen)
  end if
  if (sortRMSD2) then
    allocate (group(0:nall))
    call cregen_CRE_2(prch,env,nat,nall,at,xyz,comments,nallnew,group)
  end if

!=====================================================================!
!>  E N S E M B L E   O U T P U T
!=====================================================================!

!>--- align all structures to the first structure using the RMSD
  call cregen_rmsdalign(nat,nall,at,xyz)

!>--- write new file with ALL remaining structures
  if (newfile) then
    call cregen_file_wr(env,oname,nat,nall,at,xyz,comments)
!>--- track ensemble for restart
    call trackensemble(oname,nat,nall,at,xyz,comments)
  end if
!>--- write a file containing only conformers (no rotamers)
  if (conffile) then
    call cregen_conffile(env,cname,nat,nall,at,xyz,comments,ng,degen)
  end if
  if (saveelow) then
    env%elowest = grepenergy(comments(1))
!>-- and update reference geometry (in Bohr)
    env%ref%xyz = xyz(:,:,1)/bohr
  end if

!>--- additional files for entropy mode
  if (bonusfiles) then
    call cregen_bonusfiles(ng,degen)
  end if

!>--- several printouts
  if (pr2) then
    allocate (er(nall))
    call cregen_pr2(prch,env,nall,comments,ng,degen,er)
    call cregen_econf_list(prch,nall,er,ng,degen)
    deallocate (er)
  end if
  if (pr3) then !> alternative to pr2
    call cregen_pr3(prch,oname,nall,comments)
  end if
  if (pr4) then !> group dara printout
    call cregen_pr4(prch,fname,nall,group)
  end if

!>--- analyze nuclear equivalencies, e.g. for NMR and Entropy
  if (anal) then
    call cregen_EQUAL(prch,nat,nall,at,xyz,group,athr,.not.env%entropic)
  end if

!>--- deallocate data
  if (prch .ne. stdout) then
    close (prch)
  end if
  if (allocated(er)) deallocate (er)
  if (allocated(degen)) deallocate (degen)
  if (allocated(group)) deallocate (group)
  deallocate (xyz,comments,at)
  return
end subroutine newcregen

!=========================================================================================!
!=========================================================================================!
!>  CREGEN DATA SECTION
!=========================================================================================!
!=========================================================================================!

subroutine cregen_files(env,fname,oname,cname,simpleset,userinput,iounit)
!*************************************************************
!* subroutine cregen_files
!* handle all settings regarding input and output file names
!* including where to print the cregen output
!*************************************************************
  use crest_parameters
  use crest_data
  use iomod
  use utilities
  implicit none
  type(systemdata) :: env    !> MAIN STORAGE OS SYSTEM DATA
  character(len=*) :: fname  !> name of the ensemble to be read
  character(len=*) :: oname  !> output ensemble name (including rotamers)
  character(len=*) :: cname  !> output ensemble name (only conformers)
  integer,intent(in) :: simpleset
  logical,intent(in) :: userinput !> was an input file given via the optional subroutine arg?
  integer,intent(out) :: iounit
  character(len=:),allocatable :: outfile
  logical :: ex
  !>--------------------------------------------------------------------
  outfile = 'cregen.out.tmp'
  if (env%cgf(6)) outfile = 'tmp'

  !>-- the entire cregen output can be printed printed to a seperate file
  !>   or to the terminal
  call remove(outfile)
  if (simpleset > 0) then
    select case (simpleset)
    case (6,9,12)
      iounit = stdout
    case default
      open (newunit=iounit,file=outfile)
    end select
  else if (env%confgo.and..not. (env%properties .eq. -2).and..not.env%relax) then
    iounit = stdout
  else
    open (newunit=iounit,file=outfile)
  end if

  if ((env%confgo.and.(index(trim(fname),'none selected') .eq. 0)) &
  &    .OR.userinput) then
    if (.not.userinput) then
      fname = trim(env%ensemblename)
    end if
    oname = trim(fname)//'.sorted'
    cname = 'crest_ensemble.xyz'
    if (env%fullcre) then
      env%ensemblename = trim(oname)
    end if
  else !> internal mode for conformational search
    call checkname_xyz(crefile,fname,oname)
    cname = conformerfile
  end if
  if (simpleset == 12) then !> MECP files
    fname = "crest_mecp_search.xyz"
    oname = "crest_mecp_search.xyz.sorted"
    cname = "crest_ensemble.xyz"
  end if
  if (simpleset == 13) then !> MSREACT files
    fname = "crest_unique_products.xyz"
    oname = "crest_unique_products.sorted"
    cname = "crest_msreact_products.xyz"
  end if
  if (simpleset == 15) then !> crossing files
    call checkname_xyz('confcross',fname,oname)
    cname = trim(fname)//'.unique'
  end if

  write (iounit,*) 'input  file name : ',trim(fname)
  select case (simpleset)
  case (9)
    continue
  case default
    write (iounit,*) 'output file name : ',trim(oname)
  end select

  inquire (file=fname,exist=ex)
  if (.not.ex) then
    write (0,*) 'Warning, file ',trim(fname),' does not exist!'
    error stop
  end if

  return
end subroutine cregen_files

!=========================================================================================!

subroutine cregen_prout(env,simpleset,pr1,pr2,pr3,pr4)
!***********************************************************
!* subroutine cregen_prout
!* handle all settings regarding which printouts are active
!* (currently only those for default cregen runs)
!***********************************************************
  use crest_parameters
  use crest_data
  use iomod
  implicit none
  type(systemdata) :: env !> MAIN STORAGE OS SYSTEM DATA
  integer,intent(in) :: simpleset
  logical,intent(out) :: pr1,pr2,pr3,pr4

  pr1 = .true.  !> threshold summary
  pr2 = .true.  !> detailed energy/group list
  pr3 = .false. !> plain energy list
  pr4 = .false. !> group list printout

  if (simpleset == 6) then
    pr1 = .false.
    pr2 = .false.
    if (env%crestver .ne. crest_solv) pr3 = .true.
  end if

  if (simpleset == 9) then
    pr1 = .true.
    pr2 = .false.
    pr3 = .false.
    pr4 = .true.
  end if

  if (simpleset == 13) then
    pr1 = .false.
    pr2 = .false.
    pr3 = .false.
    pr4 = .false.
  end if

  return
end subroutine cregen_prout

!=========================================================================================!

subroutine cregen_director(env,simpleset,checkbroken,sortE,sortRMSD,sortRMSD2, &
        &  repairord,newfile,conffile,bonusfiles,anal,topocheck,checkez)
!**************************************************************
!* subroutine cregen_director !IMPORTANT!
!* handle which comparisons are required and which files shall
!* be written
!**************************************************************
  use crest_parameters
  use crest_data
  use iomod
  implicit none
  type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
  integer,intent(in) :: simpleset
  logical,intent(out) :: checkbroken
  logical,intent(out) :: sortE,sortRMSD,sortRMSD2
  logical,intent(out) :: repairord
  logical,intent(out) :: newfile,conffile
  logical,intent(out) :: bonusfiles
  logical,intent(out) :: anal
  logical,intent(out) :: topocheck
  logical,intent(out) :: checkez

  checkbroken = .true. !> fragmentized structures are sorted out
  sortE = .true.       !> sort based on energy
  sortRMSD = .true.    !> sort based on RMSD
  sortRMSD2 = .false.  !> check groups for whole ensemble
  repairord = .true.   !> double-check the sorted Ensemble

  newfile = .true.  !> sorted input file

  conffile = .true. !> sorted unique structure file

  topocheck = env%checktopo !> topology is compared to reference structure
  checkez = env%checkiso    !> check for C=C cis/trans isomerizations
  if (env%relax) then
    topocheck = .false.
  end if

  bonusfiles = .false.
  if (env%entropic.or.env%doNMR) then
    bonusfiles = .true.
  end if

  anal = .false.
  if (env%doNMR.or.env%cgf(3).or.simpleset == 2) then
    anal = .true.
  end if
  if (simpleset == 3) then
    anal = .false.
  end if

  if (simpleset == 6) then  !energy sorting only
    checkbroken = .false.
    sortE = .true.
    sortRMSD = .false.
    repairord = .false.
    newfile = .true.
    if ((env%crestver .eq. crest_solv).and.(.not.env%QCG)) then
      conffile = .true. !Conffile is needed for confscript in QCG
    else
      conffile = .false.
    end if
    topocheck = .false.
    checkez = .false.
    bonusfiles = .false.
    anal = .false.
  end if

  if (simpleset == 9) then  !optpurge mode
    checkbroken = .false.
    sortE = .false.
    sortRMSD = .false.
    sortRMSD2 = .true.
    repairord = .false.
    newfile = .false.
    conffile = .false.
    topocheck = .false.
    checkez = .false.
    bonusfiles = .false.
    anal = .false.
  end if

  !> MECP search final sorting
  if (simpleset == 12) then
    topocheck = .false.
    checkez = .false.
    bonusfiles = .false.
    anal = .false.
  end if

  if(simpleset == 13)then  !msreact mode 
    checkbroken = .false.
    sorte = .true.
    sortRMSD = .true.
    sortRMSD2 = .false.
    repairord = .false.
    newfile = .true.
    conffile = .true.
    topocheck = .false.
    checkez = .false.
    bonusfiles= .false.
    anal = .false.

  endif    

  return
end subroutine cregen_director

!=========================================================================================!

subroutine cregen_filldata1(env,ewin,rthr,ethr,bthr,athr,pthr,T,couthr)
!*******************************************************
!* subroutine cregen_filldata1
!* get important threshold from "opt" and "sys" objects
!*******************************************************
  use crest_parameters
  use crest_data
  implicit none
  type(systemdata) :: env  !> MAIN STORAGE OS SYSTEM DATA
  real(wp),intent(out) :: ewin,rthr,ethr,bthr,athr,pthr,T,couthr
  !>--------------------------------------------------------------------
  ewin = env%ewin      !> ensemble energy window in kcal/mol
  rthr = env%rthr      !> RMSD thr in Ang
  ethr = env%ethr      !> E threshold in kcal
  bthr = env%bthr2     !> rot const thr (lower bound)
  athr = env%athr      !> to det. int. rotation. equal atoms for NMR, CRITICAL!
  pthr = env%pthr      !> population thr
  T = env%tboltz       !> Temperature
  couthr = env%couthr  !> coulomb sorting threshold
  return
end subroutine cregen_filldata1
subroutine cregen_filldata2(simpleset,ewin)
  use crest_parameters
  use crest_data
  implicit none
  integer,intent(in) :: simpleset
  real(wp),intent(out) :: ewin
  if (simpleset == 6.or.simpleset == 12) then
    ewin = huge(ewin)
  end if
  return
end subroutine cregen_filldata2

!=========================================================================================!

subroutine cregen_groupinfo(nall,ng,group,degen)
!*************************************************************
!* subroutine cregen_groupinfo
!* get info about each conformer group and save it to "degen"
!*************************************************************
  use crest_parameters
  implicit none
  integer :: nall,ng
  integer :: group(0:nall)
  integer :: degen(3,ng)
  integer :: i,j,k
  integer :: a,b
  do i = 1,ng
    a = 0; b = 0; k = 0
    do j = 1,nall
      if (group(j) .eq. i) then
        k = k+1
        if (a == 0) a = j
        b = j
      end if
    end do
    degen(1,i) = k !>-- number of members in group i
    degen(2,i) = a !>-- first member of group i
    degen(3,i) = b !>-- last member of group i
  end do
  return
end subroutine cregen_groupinfo

!=========================================================================================!
!=========================================================================================!
!>  CREGEN SUBROUTINES
!=========================================================================================!
!=========================================================================================!

subroutine discardbroken(ch,env,topocheck,nat,nall,at,xyz,comments,newnall)
!**************************************************
!* subroutine discardbroken
!* analyze an ensemble and track broken structures
!* to be discarded.
!**************************************************
  use crest_parameters
  use crest_data
  use strucrd
  use miscdata,only:rcov
  implicit none
  !> INPUT
  type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
  integer,intent(in) :: ch ! printout channel
  integer,intent(in) :: nat,nall
  integer,intent(in) :: at(nat)
  logical,intent(in) :: topocheck
  !> OUTPUT
  real(wp),intent(inout) :: xyz(3,nat,nall)
  character(len=*),intent(inout) :: comments(nall)
  integer,intent(out) :: newnall
  !> LOCAL
  integer :: llan
  integer,allocatable :: order(:),orderref(:)
  integer :: nat0
  real(wp),allocatable :: cref(:,:),c0(:,:),c1(:,:)
  integer,allocatable  :: at0(:),atdum(:)
  real(wp),allocatable :: cn(:),bond(:,:)
  integer :: frag,frag0
  real(wp) :: erj
  integer :: j
  logical :: substruc
  logical :: distok,distcheck
  real(wp) :: cnorm
  logical :: dissoc

  !>--- if we don't wish to include all atoms:
  substruc = (nat .ne. env%rednat.and.env%subRMSD)

  !>--- read the reference structure
  allocate (cref(3,nat),atdum(nat))
  call rdcoord('coord',nat,atdum,cref)
  !>--- check fragements
  allocate (bond(nat,nat),cn(nat),source=0.0_wp)
  call mreclm(frag0,nat,at,cref,atdum,bond,rcov,cn)
  deallocate (bond,cn)

  write (ch,'('' # fragment in coord            :'',i6)') frag0
  deallocate (atdum)
  if (substruc) then
    nat0 = env%rednat
    allocate (c0(3,nat0),at0(nat0),c1(3,nat0),atdum(nat0))
    call maskedxyz(nat,nat0,cref,c0,at,at0,env%includeRMSD)
  else
    allocate (c0(3,nat),at0(nat),c1(3,nat),atdum(nat))
    c0 = cref
    at0 = at
    nat0 = nat
  end if
  !>--- fragments for actual reference
  allocate (bond(nat0,nat0),cn(nat0))
  call mreclm(frag0,nat0,at0,c0,atdum,bond,rcov,cn)

  allocate (order(nall),orderref(nall))
  !>--- loop over the structures
  newnall = 0
  llan = nall
  do j = 1,nall
    erj = grepenergy(comments(j)) !> get energy of structure j
    if (.not.substruc) then
      c1(:,:) = xyz(:,:,j)/bohr
    else
      call maskedxyz(nat,nat0,xyz(:,:,j),c1,at,at0,env%includeRMSD)
      c1 = c1/bohr
    end if
    distok = distcheck(nat0,c1) !> distance check
    cnorm = sum(abs(c1))        !> clash check

    !>--- further checks: dissociation?
    dissoc = .false.
    if (abs(erj) .gt. 1.0d-6.and.cnorm .gt. 1.0d-6 &
    &   .and.distok.and.topocheck) then
      dissoc = .false.
      call mreclm(frag,nat0,at0,c1,atdum,bond,rcov,cn)
      if (frag .gt. frag0) then
        dissoc = .true.
      end if
    end if

    if (dissoc.or.(cnorm .lt. 1.0d-6).or.(.not.distok)) then
      !>--- move broken structures to the end of the matrix
      orderref(j) = llan
      llan = llan-1
      !write(ch,*) 'removing structure',j
    else
      newnall = newnall+1
      orderref(j) = newnall
    end if
  end do

  !>--- sort the xyz array (only if structures have been discarded)
  if (newnall .lt. nall) then
    order = orderref
    call xyzqsort(nat,nall,xyz,c0,order,1,nall)
    order = orderref
    call stringqsort(nall,comments,1,nall,order)

    llan = nall-newnall
    write (ch,'('' number of removed clashes      :'',i6)') llan
  end if
  !>--- otherwise the ensemble is ok

  if (allocated(orderref)) deallocate (orderref)
  if (allocated(order)) deallocate (order)
  if (allocated(cn)) deallocate (cn)
  if (allocated(bond)) deallocate (bond)
  if (allocated(atdum)) deallocate (atdum)
  if (allocated(c1)) deallocate (c1)
  if (allocated(at0)) deallocate (at0)
  if (allocated(c0)) deallocate (c0)
  if (allocated(cref)) deallocate (cref)
  return
end subroutine discardbroken

!=========================================================================================!

subroutine cregen_topocheck(ch,env,checkez,nat,nall,at,xyz,comments,newnall)
!*************************************************************
!* subroutine cregen_topocheck
!* analyze an ensemble and compare topology (neighbourlist)
!* to the reference structure
!*************************************************************
  use crest_parameters
  use crest_data
  use strucrd
  use miscdata,only:rcov
  use utilities
  implicit none
  type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
  integer,intent(in) :: ch ! printout channel
  logical,intent(in) :: checkez
  integer,intent(in) :: nat,nall
  integer,intent(in) :: at(nat)
  real(wp),intent(inout) :: xyz(3,nat,nall)
  character(len=*),intent(inout) :: comments(nall)
  integer,intent(out) :: newnall
  integer :: llan
  integer,allocatable :: order(:),orderref(:)
  real(wp),allocatable :: cref(:,:),c1(:,:)
  integer,allocatable  :: atdum(:)
  real(wp),allocatable :: cn(:),bond(:,:)
  integer,allocatable :: toporef(:)
  integer,allocatable :: topo(:)
  logical,allocatable :: neighmat(:,:)
  integer :: nbonds
  integer :: j,l
  integer :: ntopo,ncc,ccfail
  logical :: discard
  integer,allocatable :: ezat(:,:)
  real(wp),allocatable :: ezdihedref(:)
  real(wp),allocatable :: ezdihed(:)
  real(wp) :: winkeldiff

  !>--- read the reference structure
  allocate (cref(3,nat),atdum(nat))
  call rdcoord('coord',nat,atdum,cref)

  !>--- get the reference topology matrix (bonds)
  ntopo = nat*(nat+1)/2
  allocate (toporef(ntopo),topo(ntopo))
  allocate (neighmat(nat,nat),source=.false.)
  allocate (bond(nat,nat),cn(nat),source=0.0_wp)
  cn = 0.0d0
  bond = 0.0d0
  call xcoord2(nat,atdum,cref,rcov,cn,400.0_wp,bond)

  if (allocated(env%excludeTOPO)) then
    call bondtotopo_excl(nat,at,bond,cn,ntopo,toporef,neighmat,env%excludeTOPO)
  else
    call bondtotopo(nat,at,bond,cn,ntopo,toporef,neighmat)
  end if

  nbonds = sum(toporef)
  write (ch,'('' # bonds in reference structure :'',i6)') nbonds
  !>--- if required, check for C=C bonds (based only on structure!)
  if (checkez) then
    cref = cref*bohr
    call nezcc(nat,atdum,cref,cn,ntopo,toporef,ncc)
    if (ncc > 0) then
      write (ch,'(''   => # of C=C bonds :'',i6)') ncc
      allocate (ezat(4,ncc))
      allocate (ezdihedref(ncc),ezdihed(ncc),source=0.0d0)
      call ezccat(nat,atdum,cref,cn,ntopo,toporef,ncc,ezat)
      call ezccdihed(nat,cref,ncc,ezat,ezdihedref)
      !do i=1,ncc
      !  write(*,'(1x,a,4i4,a,f6.2)') 'C=C bond atoms:',ezat(1:4,i)," angle: ",ezdihedref(i)
      !enddo
    end if
  end if

  allocate (order(nall),orderref(nall))
  allocate (c1(3,nat))
  !>--- loop over the structures
  ccfail = 0
  newnall = 0
  llan = nall
  do j = 1,nall
    c1(1:3,1:nat) = xyz(1:3,1:nat,j)/bohr
    !>--- generate topo and compare
    discard = .false.
    cn = 0.0d0
    bond = 0.0d0
    call xcoord2(nat,at,c1,rcov,cn,400.0_wp,bond)
    if (allocated(env%excludeTOPO)) then
      call bondtotopo_excl(nat,at,bond,cn,ntopo,topo,neighmat,env%excludeTOPO)
    else
      call bondtotopo(nat,at,bond,cn,ntopo,topo,neighmat)
    end if
    do l = 1,ntopo
      if (toporef(l) .ne. topo(l)) then
        discard = .true.   !> if there is any mismatch in neighbor lists
        exit
      end if
    end do
    !>--- get E/Z info of C=C, discard isomers
    if (checkez.and..not.discard.and.ncc > 0) then
      c1 = c1*bohr
      call ezccdihed(nat,c1,ncc,ezat,ezdihed)
      do l = 1,ncc
        winkeldiff = ezdihedref(l)-ezdihed(l)
        winkeldiff = abs(winkeldiff)
        if (winkeldiff > 90.0_wp) then
          discard = .true.
          ccfail = ccfail+1
          exit
        end if
      end do
    end if

    if (discard) then
      !>-- move broken structures to the end of the matrix
      orderref(j) = llan
      llan = llan-1
    else
      newnall = newnall+1
      orderref(j) = newnall
    end if
  end do

  !>--- sort the xyz array (only if structures have been discarded)
  if (newnall .lt. nall) then
    order = orderref
    call xyzqsort(nat,nall,xyz,c1,order,1,nall)
    order = orderref
    call stringqsort(nall,comments,1,nall,order)

    llan = nall-newnall
    write (ch,'('' number of topology mismatches  :'',i6)') llan
    !>--- report the removals during a run
    if (ch .ne. stdout) then
      write (stdout,'("CREGEN> number of topology-based structure removals: ",i0)') llan
    end if
    if (checkez.and.ccfail > 0) then
      write (ch,'(''  => discared due to E/Z isom.  :'',i6)') ccfail
    end if
  end if
  !>--- otherwise the ensemble is ok

  deallocate (c1)
  deallocate (orderref,order)
  if (allocated(ezdihedref)) deallocate (ezdihedref)
  if (allocated(ezdihed)) deallocate (ezdihed)
  if (allocated(ezat)) deallocate (ezat)
  deallocate (cn,bond)
  deallocate (neighmat)
  deallocate (topo,toporef)
  deallocate (atdum)
  deallocate (cref)
  return

contains
  subroutine nezcc(nat,at,xyz,cn,ntopo,topo,ncc)
    !***************************************************
    !* Check how many (potential) C=C bonds are present
    !***************************************************
    use crest_parameters
    integer,intent(in)  :: nat
    integer,intent(in)  :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),intent(in) :: cn(nat)
    integer,intent(in)  :: ntopo
    integer,intent(in)  :: topo(ntopo)
    integer,intent(out) :: ncc
    real(wp) :: dist
    integer :: l
    integer :: ci,cj
    real(wp),parameter :: distcc = 1.384_wp
    ncc = 0
    do ci = 1,nat
      do cj = 1,ci-1
        if (ci == cj) cycle
        l = lin(ci,cj)
        if (topo(l) == 0) cycle
        if (at(ci) == 6.and.at(cj) == 6.and. &
        &  nint(cn(ci)) == 3.and.nint(cn(cj)) == 3) then
          dist = (xyz(1,ci)-xyz(1,cj))**2+ &
          &    (xyz(2,ci)-xyz(2,cj))**2+ &
          &    (xyz(3,ci)-xyz(3,cj))**2
          dist = sqrt(dist)
          if (dist < distcc) then
            ncc = ncc+1
          end if
        end if
      end do
    end do
    return
  end subroutine nezcc
  subroutine ezccat(nat,at,xyz,cn,ntopo,topo,ncc,ezat)
    !********************************************************
    !* Check which atoms can be used for C=C dihedral angles
    !********************************************************
    use crest_parameters
    integer,intent(in)  :: nat
    integer,intent(in)  :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),intent(in) :: cn(nat)
    integer,intent(in)  :: ntopo
    integer,intent(in)  :: topo(ntopo)
    integer,intent(in)  :: ncc
    integer,intent(out) :: ezat(4,ncc)
    real(wp) :: dist
    integer :: i,j,k,l
    integer :: ci,cj
    real(wp),parameter :: distcc = 1.384_wp
    if (ncc < 1) return
    k = 0
    do ci = 1,nat
      do cj = 1,ci-1
        if (ci == cj) cycle
        l = lin(ci,cj)
        if (topo(l) == 0) cycle
        if (at(ci) == 6.and.at(cj) == 6.and. &
        &  nint(cn(ci)) == 3.and.nint(cn(cj)) == 3) then
          dist = (xyz(1,ci)-xyz(1,cj))**2+ &
          &    (xyz(2,ci)-xyz(2,cj))**2+ &
          &    (xyz(3,ci)-xyz(3,cj))**2
          dist = sqrt(dist)
          if (dist < distcc) then
            k = k+1
            ezat(2,k) = ci
            ezat(3,k) = cj
            !>-- get a neighbour for ci
            do i = 1,nat
              if (i == cj.or.i == ci) cycle
              l = lin(ci,i)
              if (topo(l) == 1) then
                ezat(1,k) = i
                exit
              end if
            end do
            !>-- get a neighbour for cj
            do j = 1,nat
              if (j == cj.or.j == ci) cycle
              l = lin(cj,j)
              if (topo(l) == 1) then
                ezat(4,k) = j
                exit
              end if
            end do
          end if
        end if
      end do
    end do
    return
  end subroutine ezccat
  subroutine ezccdihed(nat,xyz,ncc,ezat,ezdihed)
    !********************************************************
    !* Check which atoms can be used for C=C dihedral angles
    !********************************************************
    use crest_parameters
    integer,intent(in)  :: nat
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in)  :: ncc
    integer,intent(in) :: ezat(4,ncc)
    real(wp),intent(out) :: ezdihed(ncc)
    integer :: i,k
    integer :: a,b,c,d
    real(wp) :: winkel
    if (ncc < 1) return
    k = 0
    do i = 1,ncc
      a = ezat(1,i)
      b = ezat(2,i)
      c = ezat(3,i)
      d = ezat(4,i)
      call DIHED(xyz,a,b,c,d,winkel) !>-- from intmodes.f
      winkel = abs(winkel*(180.0_wp/pi))
      if (winkel > 180.0_wp) then
        winkel = 360.0_wp-winkel
      end if
      ezdihed(i) = winkel
    end do
    return
  end subroutine ezccdihed
end subroutine cregen_topocheck

!=========================================================================================!

subroutine cregen_esort(ch,nat,nall,xyz,comments,nallout,ewin)
!**************************************************************
!* subroutine cregen_esort
!* sort the ensemble by energy and determine the new
!* ensemble size within the energy threshold.
!* On Input: ch - printout channel
!*           nat - number of atoms
!*           nall - number of structure in ensemble
!*           xyz  - Cartesian coordinates
!*           comments - commentary lines containing the energy
!*           ewin - energy window in kcal/mol
!* On Output: nallout - number of strucutres after cutoff
!**************************************************************
  use crest_parameters
  use strucrd
  implicit none
  integer,intent(in) :: ch
  integer,intent(in) :: nat
  integer,intent(in) :: nall
  real(wp),intent(inout) :: xyz(3,nat,nall)
  character(len=*) :: comments(nall)
  integer :: nallout
  real(wp),intent(in) :: ewin

  real(wp),allocatable :: energies(:)
  integer,allocatable :: orderref(:)
  integer,allocatable :: order(:)
  real(wp),allocatable :: c0(:,:)
  integer :: i
  real(wp) :: de,emax

  allocate (energies(nall))
  allocate (orderref(nall),order(nall))
  do i = 1,nall
    energies(i) = grepenergy(comments(i))
    orderref(i) = i
    order(i) = i
  end do
  !>-- sort the energies and obtain the order
  call qsort(energies,1,nall,orderref)
  !>-- after the sorting orderref contains information:
  !>   before the sorting element "i" WAS at position "orderref(i)"
  !>   but to use it as a mask, we need to invert it,
  !>   so that it is: element "i" IS NOW at position "orderref(i)"
  call maskinvert(nall,orderref)

  !>-- sort structures and comments based on the order
  order = orderref
  allocate (c0(3,nat))
  call xyzqsort(nat,nall,xyz,c0,order,1,nall)
  deallocate (c0)
  order = orderref

  call stringqsort(nall,comments,1,nall,order)

  !>-- determine cut-off of energies
  emax = maxval(energies(:),1)
  de = (emax-energies(1))*autokcal
  if (de .gt. ewin) then
    nallout = 1 !> lowest is always taken
    do i = 2,nall
      de = (energies(i)-energies(1))*autokcal
      if (de .lt. ewin) then
        nallout = nallout+1
      else
        exit
      end if
    end do
    write (ch,'('' number of removed by energy    :'',i6)') (nall-nallout)
    write (ch,'('' number of remaining points     :'',i6)') nallout
  else
    nallout = nall
  end if
  write (ch,*) 'reference state Etot :',energies(1)

  deallocate (order,orderref)
  deallocate (energies)
  return
end subroutine cregen_esort

!=========================================================================================!

subroutine cregen_CRE(ch,env,nat,nall,at,xyz,comments,nallout,group)
!*************************************************************
!* subroutine cregen_CRE
!* sort the ensemble based on rotational constants,RMSD and
!* energy to determine rotamers and duplicates.
!* On Input: ch - printout channel
!*           nat - number of atoms
!*           nall - number of structure in ensemble
!*           at   - atom types
!*           xyz  - Cartesian coordinates
!*           comments - commentary lines containing the energy
!* On Output: nallout - new total number of structures
!*            group   - to which group every structure belongs
!**************************************************************
  use crest_parameters,id => dp
  use crest_data
  use strucrd
  use ls_rmsd
  use axis_module
  use utilities
  implicit none
  type(systemdata) :: env
  integer,intent(in) :: ch
  integer,intent(in) :: nat
  integer,intent(in) :: nall
  integer,intent(in) :: at(nat)
  real(wp),intent(inout) :: xyz(3,nat,nall)
  integer,intent(inout) :: group(0:nall)
  character(len=*) :: comments(nall)
  integer :: nallout
  logical :: enantio = .true.   !check for enantiomers?

  !>--- float data
  real(wp) :: ewin,rthr,bthr,pthr,ethr,athr,T
  !>--- energy data
  real(wp),allocatable :: er(:)
  integer,allocatable :: orderref(:),order(:)
  real(wp) :: de
  !>--- dummy structure data
  integer,allocatable :: includeRMSD(:)
  real(wp),allocatable :: c0(:,:),c1(:,:),cdum(:,:)
  real(wp),allocatable :: c0h(:,:),c1h(:,:)
  integer,allocatable  :: maskheavy(:)
  integer,allocatable :: at0(:)
  logical :: substruc
  integer :: nat0
  real(wp),allocatable :: rot(:,:)
  real(wp) :: rotdum(3),bdum
  !>--- RMSD data
  real(sp),allocatable :: rmat(:) !SINGLE PRECISION
  real(wp) :: rdum,dr,rdum2
  real(wp),allocatable :: gdum(:,:),Udum(:,:),xdum(:),ydum(:) !dummy tensors
  integer(id) :: klong
  integer(id),allocatable :: rmap1(:)
  integer,allocatable :: rmap2(:)
  logical :: l1,l2,l3
  !>--- CRE comparison data
  integer,allocatable :: double(:)
  logical :: equalrotaniso !> this is a function
  logical,allocatable :: mask(:)
  real(wp) :: couthr
  real(wp),allocatable :: enuc(:)
  real(wp),allocatable :: ecoul(:)
  real(wp) :: r
  integer :: i,j,k,l,natnoh
  logical :: heavy

!>--- set parameters
  call cregen_filldata1(env,ewin,rthr,ethr,bthr,athr,pthr,T,couthr)
  if (env%entropic) enantio = .false.
  heavy = env%heavyrmsd
  substruc = (nat .ne. env%rednat.and.env%subRMSD)
  if (substruc) then
    nat0 = env%rednat
    includeRMSD = env%includeRMSD
  end if
  allocate (rot(3,nall))

!>--- get energies from the comment line
  allocate (er(nall))
  allocate (orderref(nall),order(nall))
  do i = 1,nall
    er(i) = grepenergy(comments(i))
    orderref(i) = i
  end do

!>--- get dummy structure memory space
  if (substruc) then
    allocate (c0(3,nat0),c1(3,nat0),at0(nat0))
  else
    allocate (c0(3,nat),c1(3,nat),at0(nat))
    at0 = at
    nat0 = nat
  end if

!>--- transform the coordinates to CMA and get rot.constants
  do i = 1,nall
    call axis(nat,at,xyz(:,:,i)) !>-- all coordinates to CMA
    if (substruc) then
      call maskedxyz(nat,nat0,xyz(:,:,i),c1,at,at0,includeRMSD)
    else
      c1(:,:) = xyz(:,:,i)
    end if
    call axis(nat0,at0,c1,rot(1:3,i),bdum)  !>-- B0 in MHz
  end do

!>--- RMSD part
  allocate (double(nall),source=0)
  !========================================================!
  !>-- crucial point: rmat is huge. VERY huge, potentially.
  !>-- for large ensembles the size of rmat can be strongly reduced
  !>   but this requires additional tracking and counting.
  !>   It will pay off, however.
  allocate (rmap1(nall),rmap2(nall))
  klong = 0
  do i = 1,nall
    rmap1(i) = klong
    l1 = .true.
    do j = 1,i-1
      !> ekcal(j) should always be smaller than ekcal(i) because i>j
      de = (er(i)-er(j))*autokcal
      if (de .lt. ethr) then !>-- we only need RMSDs for structures close in energy
        klong = klong+1
        if (l1) then
          rmap2(i) = j
          l1 = .false.
        end if
      end if
    end do
  end do

  !>-- now klong is the size of rmat with only the minimum required RMSDs
  !>   rmat itself can be single precision.
  allocate (rmat(klong),source=0.0_sp)
  allocate (gdum(3,3),Udum(3,3),xdum(3),ydum(3))
  !>-- begin calculation of RMSDs
  klong = 0
  write (stdout,'(a)',advance='no') 'CREGEN> running RMSDs ...'
  flush (stdout)
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  if (.not.substruc) then !regular case, all atoms included in RMSD
    if (.not.heavy) then    !really, the regular case
      do i = 1,nall
        c0(1:3,1:nat) = xyz(1:3,1:nat,i)
!$OMP PARALLEL PRIVATE ( j,klong,c1,xdum,ydum,Udum,gdum,rdum,rdum2,de) &
!$OMP SHARED ( i,c0,rmat,nat,xyz,rmap1,rmap2,er,ethr,enantio)
!$OMP DO
        do j = 1,i-1
          de = (er(i)-er(j))*autokcal
          if (de .lt. ethr) then
            c1(1:3,1:nat) = xyz(1:3,1:nat,j)
            call rmsd(nat,c0,c1,0,Udum,xdum,ydum,rdum,.false.,gdum) ! all atoms
            if (enantio) then  !also check for enantiomer by inverting a coordinate
              c1(1,:) = -c1(1,:)
              call rmsd(nat,c0,c1,0,Udum,xdum,ydum,rdum2,.false.,gdum) ! all atoms
            else
              rdum2 = rdum
            end if
            klong = linr(rmap1(i),rmap2(i),j)
            rmat(klong) = real(min(rdum,rdum2),sp)
          end if
        end do
!$OMP END DO
!$OMP END PARALLEL
      end do
    else  !> heavy atom case
      natnoh = nat-counth(nat,at)
      allocate (c0h(3,natnoh),c1h(3,natnoh),source=0.0_wp)
      allocate (maskheavy(nat),source=0)
      call heavymask(nat,at,maskheavy)
      write (*,*) 'doing heavy atom rmsds with ',natnoh,' atoms'
      do i = 1,nall
        call maskedxyz2(nat,natnoh,xyz(:,:,i),c0h,maskheavy)
!$OMP PARALLEL PRIVATE ( j,klong,c1h,xdum,ydum,Udum,gdum,rdum,rdum2,de) &
!$OMP SHARED ( i,c0h,rmat,nat,xyz,rmap1,rmap2,er,ethr,enantio,maskheavy)
!$OMP DO
        do j = 1,i-1
          de = (er(i)-er(j))*autokcal
          if (de .lt. ethr) then
            call maskedxyz2(nat,natnoh,xyz(:,:,j),c1h,maskheavy)
            call rmsd(natnoh,c0h,c1h,0,Udum,xdum,ydum,rdum,.false.,gdum) ! all atoms
            if (enantio) then  !also check for enantiomer by inverting a coordinate
              c1h(1,:) = -c1h(1,:)
              call rmsd(natnoh,c0h,c1h,0,Udum,xdum,ydum,rdum2,.false.,gdum) ! all atoms
            else
              rdum2 = rdum
            end if
            klong = linr(rmap1(i),rmap2(i),j)
            rmat(klong) = real(min(rdum,rdum2),sp)
          end if
        end do
!$OMP END DO
!$OMP END PARALLEL
      end do
      deallocate (maskheavy,c1h,c0h)
    end if
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  else !substruc == .true., RMSDs only on a part of the structure
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    do i = 1,nall
      call maskedxyz2(nat,nat0,xyz(:,:,i),c0,includeRMSD)
!$OMP PARALLEL PRIVATE ( j,klong,c1,xdum,ydum,Udum,gdum,rdum,rdum2,de) &
!$OMP SHARED ( i,c0,rmat,nat,nat0,xyz,rmap1,rmap2,er,ethr,includeRMSD,enantio )
!$OMP DO
      do j = 1,i-1
        de = (er(i)-er(j))*autokcal
        if (de .lt. ethr) then
          call maskedxyz2(nat,nat0,xyz(:,:,j),c1,includeRMSD)
          call rmsd(nat0,c0,c1,0,Udum,xdum,ydum,rdum,.false.,gdum) ! all atoms
          if (enantio) then  !also check for enantiomer by inverting a coordinate
            c1(1,:) = -c1(1,:)
            call rmsd(nat0,c0,c1,0,Udum,xdum,ydum,rdum2,.false.,gdum) ! all atoms
          else
            rdum2 = rdum
          end if
          klong = linr(rmap1(i),rmap2(i),j)
          rmat(klong) = real(min(rdum,rdum2),sp)
        end if
      end do
!$OMP END DO
!$OMP END PARALLEL
    end do
  end if
  write (stdout,'(1x,a)') 'done.'
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!>-- Now, with the RMSDs and rotational constants we can kick out duplicates
  do i = 1,nall
    do j = 1,i-1
      !>-- only check for structures in energy range
      de = (er(i)-er(j))*autokcal
      if (de .lt. ethr) then
        klong = linr(rmap1(i),rmap2(i),j)
        dr = rmat(klong)
      else
        cycle
      end if
      !>-- very small RMSD --> same structure
      if (dr .lt. rthr) then
        double(i) = j !>-- "i" is the same structure as "j"
        !>-- slightly larger RMSD, but same rot. constants --> same structure
      elseif (dr .lt. 2.0_wp*rthr) then
        l1 = equalrotaniso(i,j,nall,rot,0.5d0*bthr,env%bthrmax,env%bthrshift)
        if (l1) then
          double(i) = j  !>-- "i" is the same structure as "j"
        end if
      end if
    end do
    !>-- find the original reference. k is a dummy variable
    call backtrack(double,i,k)
  end do
  deallocate (c0,c1,at0)
  deallocate (ydum,xdum,Udum,gdum)
  deallocate (rmat,rmap2,rmap1)

  !>-- for ENSO write a file with duplicate info (if required)
  call enso_duplicates(env,nall,double)

!>-- count how many duplicates we have found
  allocate (mask(nall))
  mask(:) = double(:) .ne. 0
  k = count(mask,1)
  nallout = nall-k
  deallocate (mask)
  write (ch,*) 'number of doubles removed by rot/RMSD         :',k
  write (ch,*) 'total number unique points considered further :',nallout
!>-- sort structures, energies, rot const. and comments
  j = 0
  l = nall+1
  do i = 1,nall
    if (double(i) .eq. 0) then
      j = j+1
      orderref(i) = j
    else
      l = l-1
      orderref(i) = l
    end if
  end do
  order = orderref
  allocate (cdum(3,nat))
  call xyzqsort(nat,nall,xyz,cdum,order,1,nall)
  deallocate (cdum)
  order = orderref
  call maskqsort(er,1,nall,order)
  order = orderref
  call stringqsort(nall,comments,1,nall,order)
  order = orderref
  call matqsort(3,nall,rot,rotdum,1,nall,order)

!>-- finally, determine conformer groups and their rotamers
  allocate (c1(3,nat))
  allocate (enuc(nallout))
  do k = 1,nallout
    c1(1:3,1:nat) = xyz(1:3,1:nat,k)
    enuc(k) = 0.0_wp
    do i = 1,nat-1
      do j = i+1,nat
        r = (c1(1,i)-c1(1,j))**2 &
    &    +(c1(2,i)-c1(2,j))**2 &
    &    +(c1(3,i)-c1(3,j))**2+1.d-12
        enuc(k) = enuc(k)+at(i)*at(j)/r
      end do
    end do
  end do
!>-- check energy, rot. const. and nuclear permutation
  double = 0 !>-- re-use "double"
  SORTI: do i = 1,nallout
    SORTJ: do j = 1,i-1
      !>-- energy difference
      de = (er(i)-er(j))*autokcal
      l3 = double(j) .eq. 0
      if (.not.l3) cycle
      if (abs(de) .lt. ethr) then
        !>-- rotational constant difference
        l1 = equalrotaniso(i,j,nall,rot,bthr,env%bthrmax,env%bthrshift)
        !>-- nuclear permutation
        l2 = 2.0d0*abs(enuc(i)-enuc(j))/(enuc(i)+enuc(j)) .lt. 1.d-3
        if (l1.and.l2.and.l3) then
          double(i) = j   !>-- "i" is a rotamer of "j"
          call backtrack(double,i,k)
          cycle SORTI
        end if
      end if
    end do SORTJ
  end do SORTI
!>-- assign conformer groups
  k = 0
  group = 0
  do i = 1,nallout
    if (double(i) .eq. 0) then
      k = k+1
      group(i) = k
    else
      j = double(i)
      group(i) = group(j)
    end if
  end do
  group(0) = k !>-- total number of groups

  deallocate (enuc,c1,double)
  deallocate (order,orderref)
  if (allocated(ecoul)) deallocate (ecoul)
  deallocate (er)
  deallocate (rot)
  return
contains
  function counth(nat,at) result(nh)
    implicit none
    integer :: nat
    integer :: at(nat)
    integer :: nh,i
    nh = 0
    do i = 1,nat
      if (at(i) == 1) nh = nh+1
    end do
    return
  end function counth
  subroutine heavymask(nat,at,mask)
    implicit none
    integer :: nat
    integer :: at(nat)
    integer :: mask(nat)
    integer :: i
    mask = 0
    do i = 1,nat
      if (at(i) .ne. 1) mask(i) = 1
    end do
    return
  end subroutine heavymask
end subroutine cregen_CRE

!=========================================================================================!

subroutine calc_ecoul(nat,at,xyz,cn,atommask,ecoulomb)
!*******************************************************
!* subroutine calc_ecoul
!* calculate a (artificial) coulomb energy for a given
!* structure. Some atoms can be ignored via a mask
!*******************************************************
  use crest_parameters
  implicit none
  integer :: nat
  integer :: at(nat)
  real(wp) :: xyz(3,nat) !should be in Bohrs
  real(wp) :: cn(nat)
  integer :: atommask(nat) !contains 1 or 0
  real(wp) :: ecoulomb
  integer :: i,j
  real(wp) :: cnexp,rexp,natexp
  real(wp) :: r,dum,zi,zj
  ecoulomb = 0.0_wp
  !cn=1.0_wp  !comment in to ignore cn
  cnexp = 1.0_wp
  rexp = 1.0_wp
  natexp = 3.5_wp
  do i = 1,nat
    if (atommask(i) == 0) cycle
    zi = ((cn(i)**cnexp)*float(at(i)))
    do j = 1,i-1
      if (atommask(j) == 0) cycle
      r = (xyz(1,i)-xyz(1,j))**2 &
      &  +(xyz(2,i)-xyz(2,j))**2 &
      &  +(xyz(3,i)-xyz(3,j))**2
      r = sqrt(r)
      zj = ((cn(j)**cnexp)*float(at(j)))
      dum = (zi*zj)
      dum = dum/(r**rexp)
      !dum = dum * exp( -0.5_wp * (r-5)**2)
      ecoulomb = ecoulomb+dum
    end do
  end do
  !>-- normalization to the number of included(!) atoms
  dum = float(sum(atommask))
  ecoulomb = ecoulomb/(dum**natexp)
  ecoulomb = ecoulomb*627.5095
  return
end subroutine calc_ecoul

!=========================================================================================!

subroutine cregen_CRE_2(ch,env,nat,nall,at,xyz,comments,nallout,group)
!**************************************************************
!* subroutine cregen_CRE_2
!* an alternaitve version to the above routine to be used
!* with an unsorted (small) ensemble to determine
!* duplicate groups, but nothing is sorted!
!* On Input: ch - printout channel
!*           nat - number of atoms
!*           nall - number of structure in ensemble
!*           at   - atom types
!*           xyz  - Cartesian coordinates
!*           comments - commentary lines containing the energy
!* On Output: nallout - new total number of structures
!*            group   - to which group every structure belongs
!**************************************************************
  use crest_parameters,id => dp
  use crest_data
  use strucrd
  use ls_rmsd
  use axis_module
  use utilities
  implicit none
  type(systemdata) :: env
  integer,intent(in) :: ch
  integer,intent(in) :: nat
  integer,intent(in) :: nall
  integer,intent(in) :: at(nat)
  real(wp),intent(inout) :: xyz(3,nat,nall)
  integer,intent(inout) :: group(0:nall)
  character(len=*) :: comments(nall)
  integer :: nallout
  logical :: enantio = .true.   !> check for enantiomers?

  !>--- float data
  real(wp) :: ewin,rthr,bthr,pthr,ethr,athr,T
  !>--- energy data
  real(wp),allocatable :: er(:)
  integer,allocatable :: orderref(:),order(:)
  real(wp) :: de
  !>--- dummy structure data
  real(wp),allocatable :: c0(:,:),c1(:,:)
  real(wp),allocatable :: c0h(:,:),c1h(:,:)
  integer,allocatable  :: maskheavy(:)
  integer,allocatable :: at0(:)
  integer :: nat0
  real(wp),allocatable :: rot(:,:)
  real(wp) :: bdum
  !>--- RMSD data
  real(sp),allocatable :: rmat(:) !> SINGLE PRECISION
  real(wp) :: rdum,dr,rdum2
  real(wp),allocatable :: gdum(:,:),Udum(:,:),xdum(:),ydum(:) !> dummy tensors
  integer(id) :: klong
  integer(id),allocatable :: rmap1(:)
  integer,allocatable :: rmap2(:)
  logical :: l1,l2,l3
  !>--- CRE comparison data
  integer,allocatable :: double(:)
  logical :: equalrotaniso !> this is a function
  logical,allocatable :: mask(:)
  real(wp),allocatable :: enuc(:)
  real(wp) :: couthr
  real(wp),allocatable :: ecoul(:)
  real(wp) :: r
  integer :: i,j,k,natnoh
  logical :: heavy

!>--- set parameters
  call cregen_filldata1(env,ewin,rthr,ethr,bthr,athr,pthr,T,couthr)
  if (env%entropic) enantio = .false.
  heavy = env%heavyrmsd
  allocate (rot(3,nall))

!>--- get energies from the comment line
  allocate (er(nall))
  allocate (orderref(nall),order(nall))
  do i = 1,nall
    er(i) = grepenergy(comments(i))
    orderref(i) = i
  end do

!>--- get dummy structure memory space
  allocate (c0(3,nat),c1(3,nat),at0(nat))
  at0 = at
  nat0 = nat

!>--- transform the coordinates to CMA and get rot.constants
  do i = 1,nall
    call axis(nat,at,xyz(:,:,i)) !>-- all coordinates to CMA
    c1(:,:) = xyz(:,:,i)
    call axis(nat0,at0,c1,rot(1:3,i),bdum)  !>-- B0 in MHz
  end do

!>--- RMSD part
  allocate (double(nall),source=0)
  !========================================================!
  !>-- crucial point: rmat is huge. VERY huge, potentially.
  !                  use only for small ensembles!
  klong = nall
  klong = klong*(nall+1)
  klong = klong/2
  !write(*,*) nall,klong
  allocate (rmat(klong),source=0.0_sp)
  allocate (gdum(3,3),Udum(3,3),xdum(3),ydum(3))
  !>-- begin calculation of RMSDs
  klong = 0
  write (stdout,'(a)',advance='no') 'CREGEN> running RMSDs ...'
  flush (stdout)
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  if (.not.heavy) then    !really, the regular case
    do i = 1,nall
      c0(1:3,1:nat) = xyz(1:3,1:nat,i)
!$OMP PARALLEL PRIVATE ( j,klong,c1,xdum,ydum,Udum,gdum,rdum,rdum2,de) &
!$OMP SHARED ( i,c0,rmat,nat,xyz,rmap1,rmap2,er,ethr,enantio)
!$OMP DO
      do j = 1,i-1
        de = (er(i)-er(j))*autokcal
        if (de .lt. ethr) then
          c1(1:3,1:nat) = xyz(1:3,1:nat,j)
          call rmsd(nat,c0,c1,0,Udum,xdum,ydum,rdum,.false.,gdum) ! all atoms
          if (enantio) then  !also check for enantiomer by inverting a coordinate
            c1(1,:) = -c1(1,:)
            call rmsd(nat,c0,c1,0,Udum,xdum,ydum,rdum2,.false.,gdum) ! all atoms
          else
            rdum2 = rdum
          end if
          !klong=linr(rmap1(i),rmap2(i),j)
          klong = lina(i,j)
          rmat(klong) = real(min(rdum,rdum2),sp)
        end if
      end do
!$OMP END DO
!$OMP END PARALLEL
    end do
  else  !heavy atom case
    natnoh = nat-counth(nat,at)
    allocate (c0h(3,natnoh),c1h(3,natnoh),source=0.0_wp)
    allocate (maskheavy(nat),source=0)
    call heavymask(nat,at,maskheavy)
    write (*,*) 'doing heavy atom rmsds with ',natnoh,' atoms'
    do i = 1,nall
      call maskedxyz2(nat,natnoh,xyz(:,:,i),c0h,maskheavy)
!$OMP PARALLEL PRIVATE ( j,klong,c1h,xdum,ydum,Udum,gdum,rdum,rdum2,de) &
!$OMP SHARED ( i,c0h,rmat,nat,xyz,rmap1,rmap2,er,ethr,enantio,maskheavy)
!$OMP DO
      do j = 1,i-1
        de = (er(i)-er(j))*autokcal
        if (de .lt. ethr) then
          call maskedxyz2(nat,natnoh,xyz(:,:,j),c1h,maskheavy)
          call rmsd(natnoh,c0h,c1h,0,Udum,xdum,ydum,rdum,.false.,gdum) ! all atoms
          if (enantio) then  !also check for enantiomer by inverting a coordinate
            c1h(1,:) = -c1h(1,:)
            call rmsd(natnoh,c0h,c1h,0,Udum,xdum,ydum,rdum2,.false.,gdum) ! all atoms
          else
            rdum2 = rdum
          end if
          !klong=linr(rmap1(i),rmap2(i),j)
          klong = lina(i,j)
          rmat(klong) = real(min(rdum,rdum2),sp)
        end if
      end do
!$OMP END DO
!$OMP END PARALLEL
    end do
    deallocate (maskheavy,c1h,c0h)
  end if
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  write (stdout,'(1x,a)') 'done.'
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!>-- Now, with the RMSDs and rotational constants we can kick out duplicates
  do i = 1,nall
    do j = 1,i-1
      !>-- only check for structures in energy range
      de = (er(i)-er(j))*autokcal
      if (de .lt. ethr) then
        !klong= linr(rmap1(i),rmap2(i),j)
        klong = lina(i,j)
        dr = rmat(klong)
      else
        cycle
      end if
      !>-- very small RMSD --> same structure
      if (dr .lt. rthr) then
        double(i) = j !>-- "i" is the same structure as "j"
        !>-- slightly larger RMSD, but same rot. constants --> same structure
      elseif (dr .lt. 2.0_wp*rthr) then
        l1 = equalrotaniso(i,j,nall,rot,0.5d0*bthr,env%bthrmax,env%bthrshift)
        if (l1) then
          double(i) = j  !>-- "i" is the same structure as "j"
        end if
      end if
    end do
    !>-- find the original reference. k is a dummy variable
    call backtrack(double,i,k)
  end do
  deallocate (c0,c1,at0)
  deallocate (ydum,xdum,Udum,gdum)
  deallocate (rmat)

!>-- for ENSO write a file with duplicate info (if required)
  call enso_duplicates(env,nall,double)

!>-- count how many duplicates we have found
  allocate (mask(nall))
  mask(:) = double(:) .ne. 0
  k = count(mask,1)
  nallout = nall-k
  deallocate (mask)
  write (ch,'(1x,a,i10)') 'number of doubles removed by rot/RMSD:',k
  !write(ch,*)'total number unique points remaining :',nallout

!>-- finally, determine conformer groups and their rotamers
  allocate (c1(3,nat))
  allocate (enuc(nall))
  do k = 1,nall
    c1(1:3,1:nat) = xyz(1:3,1:nat,k)
    enuc(k) = 0.0_wp
    do i = 1,nat-1
      do j = i+1,nat
        r = (c1(1,i)-c1(1,j))**2 &
    &    +(c1(2,i)-c1(2,j))**2 &
    &    +(c1(3,i)-c1(3,j))**2+1.d-12
        enuc(k) = enuc(k)+at(i)*at(j)/r
      end do
    end do
  end do
!>-- check energy, rot. const. and nuclear permutation
  SORTI: do i = 1,nall
    SORTJ: do j = 1,i-1
      !>-- energy difference
      de = (er(i)-er(j))*autokcal
      l3 = double(j) .eq. 0
      if (.not.l3) cycle
      if (abs(de) .lt. ethr) then
        !>-- rotational constant difference
        l1 = equalrotaniso(i,j,nall,rot,bthr,env%bthrmax,env%bthrshift)
        !>-- nuclear permutation
        l2 = 2.0d0*abs(enuc(i)-enuc(j))/(enuc(i)+enuc(j)) .lt. 1.d-3
        if (l1.and.l2.and.l3) then
          double(i) = j   !>-- "i" is a rotamer of "j"
          call backtrack(double,i,k)
          cycle SORTI
        end if
      end if
    end do SORTJ
  end do SORTI
!>-- assign conformer groups
  k = 0
  group = 0
  do i = 1,nall
    if (double(i) .eq. 0) then
      k = k+1
      group(i) = k
    else
      j = double(i)
      group(i) = group(j)
    end if
  end do
  group(0) = k !>-- total number of groups
  write (ch,'(1x,a,i10)') 'number of removed rotamers           :', (nallout-k)
  nallout = k
  write (ch,'(1x,a,i10)') 'total number unique points remaining :',nallout

  deallocate (enuc,c1,double)
  deallocate (order,orderref)
  if (allocated(ecoul)) deallocate (ecoul)
  deallocate (er)
  deallocate (rot)
  return
contains
  function counth(nat,at) result(nh)
    implicit none
    integer :: nat
    integer :: at(nat)
    integer :: nh,i
    nh = 0
    do i = 1,nat
      if (at(i) == 1) nh = nh+1
    end do
    return
  end function counth
  subroutine heavymask(nat,at,mask)
    implicit none
    integer :: nat
    integer :: at(nat)
    integer :: mask(nat)
    integer :: i
    mask = 0
    do i = 1,nat
      if (at(i) .ne. 1) mask(i) = 1
    end do
    return
  end subroutine heavymask
end subroutine cregen_CRE_2

!=========================================================================================!

subroutine cregen_EQUAL(ch,nat,nall,at,xyz,group,athr,rotfil)
!****************************************************************
!* subroutine cregen_EQUAL
!* subroutine for the generation of nuclear equivalencies
!* On Input: ch - printout channel
!*           nat - number of atoms
!*           nall - number of structure in ensemble
!*           at   - atom types
!*           xyz  - Cartesian coordinates
!*           group - to which group does every strucutre belong
!*           athr  - threshold for equivalency comparison
!*           rotfil - wirte anmr_rotamer file?
!* On Output: resorted xyz and comments
!****************************************************************
  use crest_parameters,id => dp
  use crest_data
  use strucrd
  use miscdata,only:rcov
  use utilities
  implicit none
  integer,intent(in) :: ch
  integer,intent(in) :: nat
  integer,intent(in) :: nall
  integer,intent(in) :: at(nat)
  real(wp),intent(inout) :: xyz(3,nat,nall)
  integer,intent(inout) :: group(0:nall)
  real(wp),intent(in) :: athr
  logical,intent(in) :: rotfil
  real(wp),allocatable :: cdum(:,:)
  integer :: ng,n
  integer :: i,j,k,l
  logical :: ex

  !>--- arrays and variable for the analysis
  integer :: gmax
  integer,allocatable :: glist(:,:)
  integer :: current
  real(wp) :: dum
  real(wp) :: shortest_distance
  integer,allocatable  :: equiv(:,:,:)
  integer,allocatable  :: pair(:),pre(:),nb(:,:)
  logical,allocatable  :: vis(:)
  real(wp),allocatable :: metric(:,:)
  real(wp),allocatable :: dist(:,:,:)
  integer,allocatable  :: relat(:,:)
  real(wp),allocatable :: tmp2(:)
  integer :: m,m1,m2,s1,s2,iat,j1,k2

  !>-- further NMR-mode related data
  integer,allocatable :: nmract(:)
  integer,allocatable :: elist(:,:),flist(:,:)
  integer,allocatable :: jnd(:)
  real(wp),allocatable :: sd(:,:),jfake(:),cn(:)

  character(len=:),allocatable :: atmp
  integer :: ig,ir,irr,nr

!>--- variable declarations
  n = nat       !> other variable name for number of atoms
  ng = group(0) !> number of different groups (conformers)
  gmax = 0      !> max number of
  do i = 1,ng
    k = 0
    do j = 1,nall
      if (group(j) == i) k = k+1
    end do
    if (k .gt. gmax) gmax = k
  end do
  allocate (glist(0:gmax,ng),source=0)
  do i = 1,ng
    k = 0
    do j = 1,nall
      if (group(j) == i) then
        k = k+1
        glist(k,i) = j !> the k-th member of group i is structure j
      end if
    end do
    glist(0,i) = k !> number of members in group i
  end do

!>---distance neighbor list
  allocate (cdum(3,nat))

!>--- set up the "pair" array --> how many bonds are between two nuclei?
  allocate (pair(n*(n+1)/2),metric(n,n),vis(n),pre(n),nb(200,n))
  cdum(1:3,1:n) = xyz(1:3,1:n,1)/bohr
  call neighdist(n,at,cdum,nb,metric)
  k = 0
  pair = 0
  do i = 1,n-1
    do j = i+1,n
!>---the shortest bond path
      current = j
      dum = shortest_distance(n,i,j,nb,metric,vis,pre)
      k = 0
      do while (pre(current) /= 0)
        current = pre(current)
        k = k+1
      end do !> End loop: while precessor(current) /= 0
      pair(lin(j,i)) = k  !> # of bonds between i and j
    end do
  end do
  deallocate (nb,pre,vis,metric)

  allocate (tmp2(n),relat(0:n,n))
  allocate (equiv(0:n,n,0:nall),dist(n,n,nall))
  equiv = 0
!>-- (costly) symmetry analyis of all rotamers for NMR. this is complicated stuff also
!>   and the end of the program where this is completed...
  do i = 1,nall
    call distance(n,xyz(:,:,i),dist(:,:,i)) !> distance matrix
    do j = 1,n
      do k = 1,n
        tmp2(k) = dist(k,j,i)*dble(at(k))  !> the distance of j to all atoms * Z to distinguish
      end do
      call qqsort(tmp2,1,n)
      dist(1:n,j,i) = tmp2(1:n)
    end do
  end do
  write (ch,*) 'compare nuclear equivalencies ...'
  do i = 1,ng
    m = glist(0,i)
    if (m .lt. 2) cycle  !> det equivalent atoms in each group
    relat = 0
    do m1 = 1,m
!$OMP PARALLEL PRIVATE ( m2, s1, s2 ) SHARED ( relat )
!$OMP DO
      do m2 = 1,m1-1        !> compare all members
        s1 = glist(m1,i)    !> struc 1
        s2 = glist(m2,i)    !> struc 2
        call compare(n,nall,s1,s2,dist,athr,relat) !> athr is distance vector equivalence threshold
      end do
!$OMP END DO
!$OMP END PARALLEL
    end do
    equiv(0:n,1:n,i) = relat(0:n,1:n)
  end do
  deallocate (dist)
!>-- symmetrize result i.e. if iat is in list of jat, jat must be in list of iat
!>   done again at the end of this part
  do i = 1,ng
    do j1 = 1,n
      m1 = equiv(0,j1,i)
      do k = 1,m1
        iat = equiv(k,j1,i)
!>-- is atom j1 in the list of atom iat?
        ex = .false.
        m2 = equiv(0,iat,i)
        do k2 = 1,m2
          if (j1 .eq. equiv(k2,iat,i)) ex = .true.
        end do
        if (.not.ex) then
          equiv(0,iat,i) = equiv(0,iat,i)+1
          equiv(equiv(0,iat,i),iat,i) = j1
        end if
      end do
    end do
  end do

!>-- inlcude equivalence info from the other conformers as well i.e.
!>   assume that all conformers have the same chemical equivalencies
!>   the result is put into equiv(:,:,0)
  equiv(0:n,1:n,0) = equiv(0:n,1:n,1)
  ILOOP: do i = 2,ng
    JLOOP: do j = 1,n
      m2 = equiv(0,j,0)      !> end of list of lowest
      MLOOP: do m = 1,equiv(0,j,i)  !> list of higher
        k = equiv(m,j,i)    !> in the one in the higher list
        M1LOOP: do m1 = 1,m2
          if (equiv(m1,j,0) .eq. k) then !> already there?
            cycle MLOOP
          end if
        end do M1LOOP
        equiv(0,j,0) = equiv(0,j,0)+1 !> append
        equiv(equiv(0,j,0),j,0) = k
      end do MLOOP
    end do JLOOP
  end do ILOOP

!>--- NMR part and writeout
!> get NMR-active nuclei
  allocate (nmract(86))
  call cregen_nmract(ch,nmract)

  allocate (elist(n,n),flist(n,n))
  ig = 0
  atmp = 'anmr_nucinfo'
  open (unit=3,file=atmp)
  write (ch,'(''::::::::::: conformer group all :::::::::::'')')
  write (3,*) n
!cccccccccccccccccc
!> chem eq. first
!cccccccccccccccccc
  elist = 0
  do i = 1,n
    m = equiv(0,i,ig)
    do k = 1,m
      l = equiv(k,i,ig)
      elist(l,i) = 1
    end do
    elist(i,i) = 1
  end do
  do i = 1,n
    do j = 1,equiv(0,i,ig)
      k = equiv(j,i,ig)
      elist(1:n,k) = elist(1:n,k)+elist(1:n,i)
    end do
  end do
!>---  prepare write out
  do i = 1,n
    k = 1
    equiv(1,i,ig) = i
    elist(i,i) = 0
    do j = 1,n
      if (elist(j,i) .ne. 0) then
        k = k+1
        equiv(k,i,ig) = j
      end if
    end do
    equiv(0,i,ig) = k
  end do
  write (ch,*) 'chemical equivalencies (mag.active nuclei):'

  allocate (jnd(n))
  jnd = 1
  do j = 1,n
    m = equiv(0,j,ig)
    write (3,'(3x,i0,3x,i0)') j,m
    do l = 1,m
      if (l .ne. m) then
        write (3,'(1x,i0)',advance='no') equiv(l,j,ig)  ! include the atom ie if there are no equiv.
      else
        write (3,'(1x,i0)',advance='yes') equiv(l,j,ig)
      end if
    end do
    if (nmract(at(j)) .eq. 0) cycle
    if (m .gt. 1.and.jnd(j) .eq. 1) then  ! just print
      write (ch,'(''reference atom'',i4,'' # :'',i2)') equiv(1,j,ig),m
      do k = 1,m
        jnd(equiv(k,j,ig)) = 0
      end do
    end if
  end do
!cccccccccccccccccc
!> mag eq.
!cccccccccccccccccc
!> make a check list of atoms for the mag. eq.
  elist = 0
  flist = 1
  do i = 1,n
    m = equiv(0,i,ig) !> the following lines fill the equiv list
    do k = 1,m
      l = equiv(k,i,ig)
      elist(l,i) = 1
    end do
    elist(i,i) = 1
  end do
  flist = elist

  do i = 1,n
    m = equiv(0,i,ig)
    do k = 1,m
      l = equiv(k,i,ig)
      if (l .eq. i) cycle
      do j = 1,n
        if (flist(j,i) .eq. 1.or.nmract(at(j)) .eq. 0) cycle !> don't check non-magnetic nuclei
!c              write(*,*) l,j,pair(lin(i,j)),pair(lin(l,j)) !> and chem. equiv. ones (ie in the same
        if (pair(lin(i,j)) .ne. pair(lin(l,j))) elist(l,i) = 0 !> group
      end do
    end do
  end do
!>---  symmetrize
  do i = 1,n
    k = 1
    equiv(1,i,ig) = i
    elist(i,i) = 0
    do j = 1,n
      if (elist(j,i) .ne. 0) then
        k = k+1
        equiv(k,i,ig) = j
      end if
    end do
    equiv(0,i,ig) = k
  end do
  do i = 1,n
    do j = 1,equiv(0,i,ig)
      k = equiv(j,i,ig)
      elist(1:n,k) = elist(1:n,k)+elist(1:n,i)
    end do
  end do
!>---  prepare write out
  do i = 1,n
    k = 1
    equiv(1,i,ig) = i
    elist(i,i) = 0
    do j = 1,n
      if (elist(j,i) .ne. 0) then
        k = k+1
        equiv(k,i,ig) = j
      end if
    end do
    if (k .gt. 2) then
      equiv(0,i,ig) = k    !> CH3 etc
    else
      equiv(0,i,ig) = 1    !> this makes CH2-CH2 not mag. equiv.
    end if
  end do
  jnd = 1
  write (ch,*) 'magnetic equivalencies:'
  do j = 1,n
    m = equiv(0,j,ig)
    write (3,*) j,m
    write (3,'(20i5)') (equiv(l,j,ig),l=1,m)  !> include the atom ie if there are no equiv.
    if (nmract(at(j)) .eq. 0) cycle
    if (m .gt. 1.and.jnd(j) .eq. 1) then  !> just print
      write (ch,'(''reference atom'',i4,'' # :'',i2)') equiv(1,j,ig),m
      do k = 1,m
        jnd(equiv(k,j,ig)) = 0
      end do
    end if
  end do
  close (3)

!ccccccccccccccccccccc
!c J averaging matrix
!ccccccccccccccccccccc
  if (rotfil) then
    allocate (jfake(n*(n+1)/2),sd(n,n),cn(n))
    atmp = 'anmr_rotamer'
    open (unit=112,file=atmp,form='unformatted')
    write (112) ng
    jfake = 0
    do ig = 1,ng       !> all conf groups
      nr = glist(0,ig) !> how many rotamers?
      write (112) nr
      do ir = 1,nr
        irr = glist(ir,ig)
        call distance(n,xyz(:,:,irr),sd)   !> distance matrix
        cdum(1:3,1:n) = xyz(1:3,1:n,irr)
        call ncoord(n,rcov,at,cdum,cn,500.0d0)
        do i = 1,n-1
          do j = i+1,n
            jfake(lin(j,i)) = cn(i)*cn(j)*sqrt(dble(at(i)*at(j))) &
       &    /(dble(pair(lin(j,i)))*sd(j,i)**5) !> the approx. "J" is topologically equivalent to J
            !> R^3 was wrong in one case because Hs were artificially paired
            !> R^5 seems to be save
          end do
        end do
        write (112) jfake(1:n*(n+1)/2)  !> read by anmr
      end do
    end do
    close (112)
    deallocate (cn)
    deallocate (sd,jfake)
  end if

  if (allocated(jnd)) deallocate (jnd)
  if (allocated(elist)) deallocate (elist)
  if (allocated(flist)) deallocate (flist)
  if (allocated(equiv)) deallocate (equiv)
  if (allocated(dist)) deallocate (dist)
  if (allocated(relat)) deallocate (relat)
  if (allocated(tmp2)) deallocate (tmp2)
  if (allocated(pair)) deallocate (pair)
  if (allocated(glist)) deallocate (glist)

  return
end subroutine cregen_EQUAL

!=========================================================================================!

subroutine cregen_nmract(ch,nmract)
!***************************************
!* utility routine to fill nmract array
!***************************************
  use crest_parameters
  use utilities
  implicit none
  integer :: nmract(86)
  character(len=:),allocatable :: atmp
  logical :: fail
  integer :: io,i
  real(wp) :: xx(10)
  integer :: ch,ich2

  nmract = 0 !reset
  !>--- get NMR active nuclei
  atmp = '.anmrrc'  ! <--- name of the .anmrrc written by ENSO
  call getanmrrc(atmp,fail)
  if (fail) then  !>--- there is no .anmrrc from ENSO
    !write(ch,*)'NMR mode.'
    nmract = 0       ! all nuclei inactive
    nmract(1) = 1  ! H active
    !nmract(6) = 1  ! C active
    nmract(9) = 1  ! F active
    nmract(15) = 1  ! P active
  else          !>--- there IS a .anmrrc, and it is used.
    write (ch,*) 'NMR mode. Reading <',trim(atmp),'> for atomic NMR data'
    open (newunit=ich2,file=atmp)
    read (ich2,'(a)') atmp
    read (ich2,'(a)') atmp
    if (index(atmp,'ENSO') .ne. 0) then
      read (ich2,'(a)') atmp
    end if
    do
      read (ich2,*,iostat=io) i,xx(1:2),nmract(i)
      if (io < 0) exit
    end do
    close (ich2)
  end if

  return
end subroutine cregen_nmract

!=========================================================================================!

subroutine cregen_repairorder(nat,nall,xyz,comments,group)
!***************************************************************
!* subroutine cregen_repairorder
!* resort the ensemble to have groups grouped together
!* (can be important for small energy diff. between conformers)
!* On Input: ch - printout channel
!*           nat - number of atoms
!*           nall - number of structure in ensemble
!*           xyz  - Cartesian coordinates
!*           comments - commentary lines containing the energy
!*           group - to which group does every strucutre belong
!* On Output: resorted xyz and comments
!***************************************************************
  use crest_parameters,id => dp
  use crest_data
  use strucrd
  use utilities
  implicit none
  integer,intent(in) :: nat
  integer,intent(in) :: nall
  real(wp),intent(inout) :: xyz(3,nat,nall)
  integer,intent(inout) :: group(0:nall)
  character(len=*) :: comments(nall)
  real(wp),allocatable :: cdum(:,:)
  integer,allocatable :: order(:),orderref(:)
  character(len=128) :: btmp
  real(wp) :: edum
  logical :: ttag
  integer,allocatable :: timetag(:)
  integer :: ng,tmax
  integer :: i,j,k,l

  !>-- check if timetag info is present?
  ttag = .false.

  ng = group(0)
  allocate (order(nall),orderref(nall))
  !>-- determine new order
  k = 0
  if (ttag) then
    do i = 1,ng
      do l = 1,tmax  !>-- with timetag info
        do j = 1,nall
          if (group(j) .eq. i.and.timetag(j) .eq. l) then
            k = k+1
            orderref(k) = j
            order(k) = i
          end if
        end do
      end do
    end do
  else
    do i = 1,ng
      do j = 1,nall !>-- without timetag info
        if (group(j) .eq. i) then
          k = k+1
          orderref(k) = j
          order(k) = i
        end if
      end do
    end do
  end if

  !>-- sort xyz and comments
  group(1:nall) = order(1:nall)
  order = orderref
  allocate (cdum(3,nat))
  call xyzqsort(nat,nall,xyz,cdum,order,1,nall)
  deallocate (cdum)
  order = orderref
  call stringqsort(nall,comments,1,nall,order)
  if (ttag) then
    edum = grepenergy(comments(1))
    write (btmp,*) edum,'!t1'
    comments(1) = trim(btmp)
  end if

  deallocate (orderref,order)
  if (allocated(timetag)) deallocate (timetag)
  return
end subroutine cregen_repairorder

!=========================================================================================!

recursive subroutine xyzqsort(nat,nall,xyz,c0,ord,first,last)
!*************************************************************
!* recursive subroutine xyzqsort
!* A quicksort derivative for sorting an ensemble.
!* On Input: nat  - number of atoms
!*           nall - number of structures
!*           xyz  - the ensemble ( xyz(3,nat,nall) )
!*           c0   - a dummy coord field for sorting
!*           ord  - order of the ensemble ( ord(nall) )
!*           first - lower limit of sorting (nall dimension)
!*           last  - upperl limit of sorting (nall dimension)
!*************************************************************
  use crest_parameters
  implicit none
  integer :: nat,nall
  real(wp) :: xyz(3,nat,nall)
  real(wp) :: c0(3,nat)
  integer :: ord(nall)
  integer :: first,last
  integer :: x,t
  integer :: i,j
  x = ord((first+last)/2)
  i = first
  j = last
  do
    do while (ord(i) < x)
      i = i+1
    end do
    do while (x < ord(j))
      j = j-1
    end do
    if (i >= j) exit
    t = ord(i); ord(i) = ord(j); ord(j) = t
    c0(:,:) = xyz(:,:,i)
    xyz(:,:,i) = xyz(:,:,j)
    xyz(:,:,j) = c0(:,:)
    i = i+1
    j = j-1
  end do
  if (first < i-1) call xyzqsort(nat,nall,xyz,c0,ord,first,i-1)
  if (j+1 < last) call xyzqsort(nat,nall,xyz,c0,ord,j+1,last)
end subroutine xyzqsort

!=========================================================================================!

subroutine maskedxyz(n,nm,c,cm,at,atm,mask)
!************************************************
!* a small routine to get masked xyz coordinates
!************************************************
  use crest_parameters,only:wp
  implicit none
  integer,intent(in) :: n
  integer,intent(in) :: nm
  real(wp),intent(in) :: c(3,n)
  real(wp),intent(out) :: cm(3,nm)
  integer,intent(in) :: at(n)
  integer,intent(out) :: atm(nm)
  integer,intent(in) :: mask(n)
  integer :: i,k
  k = 1
  do i = 1,n
    if (mask(i) .gt. 0) then
      cm(1:3,k) = c(1:3,i)
      atm(k) = at(i)
      k = k+1
    end if
  end do
  return
end subroutine maskedxyz
subroutine maskedxyz2(n,nm,c,cm,mask)
!************************************************
!* a small routine to get masked xyz coordinates
!* version without at array
!************************************************
  use crest_parameters,only:wp
  implicit none
  integer,intent(in) :: n
  integer,intent(in) :: nm
  real(wp),intent(in) :: c(3,n)
  real(wp),intent(out) :: cm(3,nm)
  integer,intent(in) :: mask(n)
  integer :: i,k
  k = 1
  do i = 1,n
    if (mask(i) .gt. 0) then
      cm(1:3,k) = c(1:3,i)
      k = k+1
    end if
  end do
  return
end subroutine maskedxyz2

!=========================================================================================!

subroutine cregen_file_wr(env,fname,nat,nall,at,xyz,comments)
!*********************************
!* write the output ensemble file
!*********************************
  use crest_parameters,only:wp
  use crest_data
  use strucrd
  use utilities,only:boltz
  implicit none
  type(systemdata) :: env
  character(len=*) :: fname
  integer :: nat,nall
  integer :: at(nat)
  real(wp) :: xyz(3,nat,nall)
  character(len=*) :: comments(nall)
  character(len=128) :: newcomment

  integer :: ich,i
  real(wp),allocatable :: c0(:,:),xdum(:)
  real(wp) :: eref,T
  real(wp),allocatable :: er(:),erel(:),p(:)
  character(len=40),allocatable :: origin(:)
  real(wp),parameter :: autokcal = 627.509541_wp

  allocate (er(nall),erel(nall),p(nall),origin(nall))
  eref = grepenergy(comments(1))
  do i = 1,nall
    er(i) = grepenergy(comments(i))
    erel(i) = (er(i)-eref)*autokcal
    if (env%trackorigin) then
      call getorigin(comments(i),origin(i))
    end if
  end do
  T = env%tboltz
  call boltz(nall,T,erel,p)

  allocate (c0(3,nat),xdum(3))
  open (newunit=ich,file=fname)
  do i = 1,nall
    c0(:,:) = xyz(:,:,i)
    if (env%trackorigin) then
      write (newcomment,*) er(i),p(i),'!'//trim(origin(i))
    else
      write (newcomment,*) er(i),p(i)
    end if
    call wrxyz(ich,nat,at,c0,newcomment)
  end do
  close (ich)
  deallocate (xdum,c0)
  deallocate (origin,p,erel,er)
  return
end subroutine cregen_file_wr

!=========================================================================================!

subroutine cregen_conffile(env,cname,nat,nall,at,xyz,comments,ng,degen)
!*********************************
!* write the output ensemble file
!*********************************
  use crest_parameters,only:wp,bohr
  use crest_data
  use strucrd
  use iomod
  use utilities
  implicit none
  type(systemdata) :: env
  character(len=*) :: cname
  integer :: nat,nall
  integer :: at(nat)
  real(wp) :: xyz(3,nat,nall)
  character(len=*) :: comments(nall)
  integer :: ng
  integer :: degen(3,ng)
  character(len=128) :: newcomment
  integer :: ich,ich3,ichenso
  integer :: i,k
  real(wp),allocatable :: c0(:,:)
  real(wp),allocatable :: er(:)
  allocate (er(nall))
  do i = 1,nall
    er(i) = grepenergy(comments(i))
  end do
  allocate (c0(3,nat))
  if (env%enso) then
    open (newunit=ichenso,file='enso.tags')
  end if
  c0(:,:) = xyz(:,:,1)
  write (newcomment,'(2x,f18.8)') er(1)
  call wrxyz('crest_best.xyz',nat,at,xyz(:,:,1),newcomment)
  open (newunit=ich,file=trim(cname))
  do i = 1,ng
    k = degen(2,i)
    if (i .eq. 1.or.env%printscoords) then   !write a scoord.* for each conformer? scoord.1 is always written
      call getname1(i,newcomment)
      c0(:,:) = xyz(:,:,k)/bohr
      call wrc0(newcomment,nat,at,c0)
    end if
    write (newcomment,'(2x,f18.8)') er(k)
    call wrxyz(ich,nat,at,xyz(:,:,k),newcomment)
    if (env%enso) write (ichenso,*) trim(comments(k))
  end do
  close (ich)
  deallocate (c0)
  deallocate (er)

  if (env%enso) then
    close (ichenso)
  end if

  call remove('cre_members')
  open (newunit=ich3,file='cre_members')
  write (ich3,'(3x,i0)') ng
  do i = 1,ng
    k = degen(1,i)
    write (ich3,'(3x,i8,1x,i10,1x,i10)') &
    &   k,degen(2,i),degen(3,i)
  end do
  close (ich3)

  return
end subroutine cregen_conffile

!=========================================================================================!

subroutine cregen_rmsdalign(nat,nall,at,xyz)
!*****************************************************
!* Algin all structures in xyz to the first structure
!* in the ensemble based on the heavy-atom RMSD
!*****************************************************
  use crest_parameters,only:wp
  use crest_data
  use ls_rmsd
  use iomod
  implicit none
  integer :: nat,nall
  integer :: at(nat)
  real(wp),intent(inout) :: xyz(3,nat,nall)
  integer :: nath
  integer :: i,j,k
  real(wp),allocatable :: c0(:,:)
  real(wp),allocatable :: c1(:,:)
  real(wp),allocatable :: c2(:,:)
  real(wp) ::  g(3,3),U(3,3),x_center(3),y_center(3),rmsdval

  nath = 0
  do j = 1,nat
    if (at(j) > 2) nath = nath+1
  end do
  allocate (c0(3,nath),c1(3,nath),source=0.0d0)

  allocate (c2(3,nat))
  !>--- get the reference structure (the first one)
  i = 0
  do j = 1,nat
    if (at(j) > 2) then
      i = i+1
      c0(1:3,i) = xyz(1:3,j,1)
    end if
  end do

  do k = 2,nall
    !>--- and the other strucutres into c1
    i = 0
    do j = 1,nat
      if (at(j) > 2) then
        i = i+1
        c1(1:3,i) = xyz(1:3,j,k)
      end if
    end do
    call rmsd(i,c1,c0,1,U,x_center,y_center,rmsdval,.false.,g)

    c2 = matmul(U(1:3,1:3),xyz(1:3,1:nat,k))
    xyz(1:3,1:nat,k) = c2
  end do

  deallocate (c2)
  deallocate (c1,c0)
  return
end subroutine cregen_rmsdalign

!=========================================================================================!

subroutine cregen_bonusfiles(ng,degen)
!*****************************************
!* write the time tag and degeneracy file
!*****************************************
  use crest_parameters,only:wp,bohr
  use crest_data
  implicit none
  integer :: ng
  integer :: degen(3,ng)
  integer :: i,ich

  !>--- how many rotamers per conformer
  open (newunit=ich,file='cre_degen')
  write (ich,'(3x,i0)') ng
  do i = 1,ng
    write (ich,'(3x,i0,2x,i0)') i,degen(1,i)
  end do
  close (ich)

  return
end subroutine cregen_bonusfiles

!=========================================================================================!
!=========================================================================================!
!>  CREGEN PRINTOUTS
!=========================================================================================!
!=========================================================================================!
subroutine cregen_setthreads(ch,env,pr)
  use crest_parameters
  use crest_data
  implicit none
  type(systemdata) :: env
  integer :: ch
  logical :: pr
  integer :: TID,OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM,nproc,T,Tn
!>---- setting the threads for OMP parallel usage
  if (env%autothreads) then
    call new_ompautoset(env,'max',0,T,Tn)
!$OMP PARALLEL PRIVATE(TID)
    TID = OMP_GET_THREAD_NUM()
    IF (TID .EQ. 0.and.pr) THEN
      nproc = OMP_GET_NUM_THREADS()
      write (ch,*) '============================='
      write (ch,*) ' # threads =',nproc
      write (ch,*) '============================='
    END IF
!$OMP END PARALLEL
  end if
  return
end subroutine cregen_setthreads

subroutine cregen_pr1(ch,env,nat,nall,rthr,bthr,pthr,ewin)
  use crest_parameters
  use crest_data
  implicit none
  integer :: ch
  type(systemdata) :: env
  integer :: nat
  integer :: nall
  real(wp) :: rthr,bthr,pthr,ewin
  logical :: substruc
  substruc = (nat .ne. env%rednat.and.env%subRMSD)
  write (ch,'('' number of atoms                :'',3x,i0)') nat
  if (substruc) then
    write (ch,'('' atoms included in RMSD         :'',3x,i0)') env%rednat
  end if
  write (ch,'('' number of points on xyz files  :'',3x,i0)') nall
  write (ch,'('' RMSD threshold                 :'',f9.4)') rthr
  write (ch,'('' Bconst threshold               :'',f9.4)') bthr
  write (ch,'('' population threshold           :'',f9.4)') pthr
  if (ewin < 9999.9_wp) then
    write (ch,'('' conformer energy window  /kcal :'',f9.4)') ewin
  else
    write (ch,'('' conformer energy window  /kcal :'',6x,a)') '+∞'
  end if

  return
end subroutine cregen_pr1

subroutine enso_duplicates(env,nall,double)
  use crest_parameters
  use crest_data
  implicit none
  type(systemdata) :: env
  integer :: nall
  integer :: double(nall)
  integer :: i,j,ich

  if (.not.env%ENSO.or..not.env%confgo) return

  j = sum(double)
  open (newunit=ich,file='cregen.enso')
  if (j .gt. 0) then
    do i = 1,nall
      if (double(i) .gt. 0) then
        write (ich,*) i,double(i)
      end if
    end do
  else
    write (ich,*) "ALL UNIQUE"
  end if

  return
end subroutine enso_duplicates

subroutine create_anmr_dummy(nat)
  implicit none
  integer :: nat
  integer :: i,ich

  open (newunit=ich,file='anmr_nucinfo')
  write (ich,*) nat
  do i = 1,nat
    write (ich,'(3x,i0,3x,i0)') i,1
    write (ich,'(3x,i0)') i
  end do
  do i = 1,nat
    write (ich,*) i,1
    write (ich,'(i5)') i
  end do
  close (ich)
  return
end subroutine create_anmr_dummy

subroutine cregen_pr2(ch,env,nall,comments,ng,degen,er)
  use crest_parameters
  use crest_data
  use strucrd
  use iomod,only:touch
  use utilities,only:boltz
  implicit none
  integer :: ch
  type(systemdata) :: env
  integer :: nall
  character(len=*) :: comments(nall)
  integer :: ng
  integer :: degen(3,ng)
  real(wp) :: er(nall)
  integer :: ich,chref
  integer :: i,j,k
  real(wp),allocatable :: erel(:),egrp(:)
  real(wp),allocatable :: p(:),pg(:)
  real(wp) :: eref,T
  character(len=40),allocatable :: origin(:)
  integer :: a,b
  logical :: ex
  real(wp) :: A0,eav,g,s,ss,beta,elow

  allocate (origin(nall),erel(nall),p(nall))
  eref = grepenergy(comments(1))
  env%elowest = eref
  do i = 1,nall
    er(i) = grepenergy(comments(i))
    erel(i) = (er(i)-eref)*autokcal
    if (env%trackorigin) then
      call getorigin(comments(i),origin(i))
    else
      origin(i) = ''
    end if
  end do
  T = env%tboltz
  call boltz(nall,T,erel,p)
  allocate (pg(ng),source=0.0_wp)
  do i = 1,ng
    a = degen(2,i)
    b = degen(3,i)
    do j = a,b
      pg(i) = pg(i)+p(j)
    end do
  end do

  !>-- really long energy list
  write (ch,'(7x,a,8x,a,1x,a,2x,a,5x,a,3x,a,5x,a)') 'Erel/kcal','Etot', &
  &    'weight/tot','conformer','set','degen','origin'
  if (env%entropic.and.ng > 50000) then
    write (ch,'(1x,a)') '<skipped due to lenght and written to seperate file>'
    chref = ch
    open (newunit=ch,file='crest.conformerlist')
  end if
  k = 0
  do i = 1,ng
    k = k+1
    a = degen(2,i)
    b = degen(3,i)
    write (ch,'(i8,f8.3,1x,3f11.5,2i8,5x,a)') &
  & k,erel(a),er(a),p(a),pg(i),i,degen(1,i),trim(origin(a))
    if (.not.env%entropic) then
      do j = a+1,b
        k = k+1
        write (ch,'(i8,f8.3,1x,2f11.5,32x,a)') &
    &   k,erel(j),er(j),p(j),trim(origin(j))
      end do
    end if
  end do
  if (env%entropic.and.ng > 50000) then
    close (ch)
    ch = chref
  end if

  !>-- file for the '-compare' mode
  if (env%compareens) then
    open (newunit=ich,file='.cretrack')
    write (ich,'(5x,i0)') ng
    do i = 1,ng
      write (ich,'(1x,i8,1x,i7,1x,i7)') i,degen(2,i),degen(3,i)
    end do
    close (ich)
  end if

  !>-- some ensemble data, entropy and G (including all structures)
  A0 = 0
  eav = 0
  do i = 1,nall
    A0 = A0+p(i)*log(p(i)+1.d-12)
    eav = eav+p(i)*erel(i)
  end do
  beta = 1.0d0/(T*8.314510/4.184/1000.+1.d-14)
  g = (1.0d0/beta)*A0
  s = -1000.0d0*4.184*g/T
  ss = -1000.0d0*g/T

  write (ch,'(''T /K                                  :'', F9.2)') T
  write (ch,'(''E lowest                              :'',f12.5)') eref
  !>---- elow printout in between routines
  if (.not.env%confgo) then
    write (stdout,'("CREGEN> E lowest :",f12.5)') eref
  end if
  if (env%QCG) then
    write (ch,'(''ensemble average energy (kcal)        :'', F14.8)') eav
    write (ch,'(''ensemble entropy (cal/mol K)          :'',F14.8)') ss
    write (ch,'(''ensemble free energy (kcal/mol)       : '',F14.8)') g
  else
    write (ch,'(''ensemble average energy (kcal)        :'', F9.3)') eav
    write (ch,'(''ensemble entropy (J/mol K, cal/mol K) :'',2F9.3)') s,ss
    write (ch,'(''ensemble free energy (kcal/mol)       : '',F8.3)') g
  end if
  write (ch,'(''population of lowest in %             : '',F8.3)') pg(1)*100.d0

  !>-- some ensemble data, entropy and G (including only unique conformers)
  allocate (egrp(ng),source=0.0_wp)
  do i = 1,ng
    a = degen(2,i)
    egrp(i) = (er(a)-eref)*autokcal
  end do
  call boltz(ng,T,egrp,pg)
  A0 = 0
  do i = 1,ng
    A0 = A0+pg(i)*log(pg(i)+1.d-12)
  end do
  deallocate (egrp)
  beta = 1.0d0/(T*8.314510/4.184/1000.+1.d-14)
  g = (1.0d0/beta)*A0
  ss = -1000.0d0*g/T
  env%emtd%sapprox = ss  !> save for entropy mode

  !>-- MF-MD-GC legacy option
  if ((env%crestver .eq. 1).and.(.not.env%confgo)) then
    inquire (file='.tmpxtbmodef',exist=ex)
    if (ex) then
      open (unit=66,file='.tmpxtbmodef')
      read (66,*) i
      read (66,*) elow
      close (66)
    else
      elow = er(1)
    end if
    if ((elow-eref)*autokcal .lt. -0.2) then
      write (ch,*) '...............................................'
      write (ch,*) 'WARNING: new (best) energy less than that from '
      write (ch,*) 'WARNING: preceding Hessian calculation:  '
      write (ch,*) 'Improved by ',elow-eref,' Eh or ', (elow-eref)*autokcal,'kcal'
      write (ch,*) '...............................................'
      call touch('LOWER_FOUND')
    end if
  end if

  deallocate (pg)
  deallocate (p,erel,origin)
  return
end subroutine cregen_pr2

subroutine cregen_econf_list(ch,nall,er,ng,degen)
  use crest_parameters
  implicit none
  integer :: nall
  real(wp) :: er(nall)
  integer :: ng
  integer :: degen(3,ng)
  integer :: ch,ich2,i,j
  real(wp) :: eref,ewrt

  write (ch,*) 'number of unique conformers for further calc ',ng
  write (ch,*) 'list of relative energies saved as "crest.energies"'
  open (newunit=ich2,file='crest.energies')
  eref = er(1)
  do i = 1,ng
    j = degen(2,i)
    ewrt = er(j)-eref
    ewrt = ewrt*autokcal
    write (ich2,'(2x,i0,2x,f12.3)') i,ewrt
  end do
  close (ich2)

  return
end subroutine cregen_econf_list

subroutine cregen_pr3(ch,infile,nall,comments)
  use crest_parameters
  use strucrd
  implicit none
  integer :: ch
  character(len=*) :: infile
  integer :: nall
  character(len=*) :: comments(nall)
  real(wp),allocatable :: er(:)
  real(wp) :: dE
  integer :: i
  allocate (er(nall))

  do i = 1,nall
    er(i) = grepenergy(comments(i))
  end do

  write (ch,*)
  write (ch,'(a)') '==================================================='
  write (ch,'(a)') '============= ordered structure list =============='
  write (ch,'(a)') '==================================================='
  write (ch,'(a,a,a)') ' written to file <',trim(infile),'>'
  write (ch,*)
  write (ch,'(''   structure    ΔE(kcal/mol)    Etot(Eh)'')')
  do i = 1,nall
    dE = (er(i)-er(1))*autokcal
    write (ch,'(i10,3x,F12.2,2x,F14.6)') i,dE,er(i)
  end do
  write (ch,*)

  deallocate (er)
  return
end subroutine cregen_pr3

subroutine cregen_pr4(ch,infile,nall,group)
  use crest_parameters
  use strucrd
  implicit none
  integer :: ch
  character(len=*) :: infile
  integer :: nall
  integer :: group(0:nall)
  integer :: i,ich
  integer :: maxgroup
  !write(ch,*) group(1:nall)
  maxgroup = group(0)
  write (ch,'(1x,i0,a,i0,a,a,a)') maxgroup,' unique groups for ', &
  &    nall,' structures in file <',trim(infile),'>'

  open (newunit=ich,file='.groups')
  write (ich,'(5x,i0,3x,i0)') nall,maxgroup
  do i = 1,nall
    write (ich,'(2x,i10,2x,i10)') i,group(i)
  end do
  close (ich)
  return
end subroutine cregen_pr4

!=========================================================================================!
!=========================================================================================!
!> END OF CREGEN FILE
!=========================================================================================!
!=========================================================================================!
