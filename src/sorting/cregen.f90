!===============================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2026 Philipp Pracht
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
!> CREGEN - also see cregen_interfaces.f90 for importable interfaces
!=========================================================================================!
!=========================================================================================!

subroutine newcregen(env,quickset,infile,structurelist)
!****************************************************************************************
!* The main CREGEN routine
!*
!* CREGEN is the universal ensemble sorting routine of CREST.
!* This is a rewrite of the original routines since the old ones
!* got a bit messy over time.
!* The quickset variable can be used for some special runtypes:
!*   quickset:  2   - do symmetry analysis
!*              3   - switch off equivalency analysis
!*              6,7 - energy sorting only with (7) or without (6) ewin energy cut-off
!*              9   - no sorting, only check groups
!*             12   - no topology check, turn ewin to infty
!*             13   - no topology check, ewin and rmsd checking (msreact settings)
!****************************************************************************************
  use crest_parameters
  use crest_data
  use strucrd
  use cregen_subroutines
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env   !> MAIN STORAGE OS SYSTEM DATA
  integer,intent(in),optional :: quickset !> quick access to predefined CREGEN modes
  character(len=*),intent(in),optional :: infile
  type(coord),allocatable,intent(inout),optional :: structurelist(:)
  !> LOCAL
  integer :: simpleset
  character(len=:),allocatable :: fname  !> input file
  character(len=:),allocatable :: oname  !> sorted output file
  character(len=:),allocatable :: cname  !> unique structure file
!>--- ensemble arguments
  integer :: nat                      !> number of atoms
  integer :: nall                     !> number of structures
  character(len=128),allocatable :: comments(:)
  real(wp),allocatable :: er(:)       !> energies
  type(coord),allocatable :: structures(:)  !> a list of structures using the coord type
!>--- dummy ensemble arguments
  integer :: nallref
  integer :: nallnew
!>--- sorting arguments
  integer,allocatable :: gref(:),group(:)
  integer :: ng
  integer :: i,ii
  integer,allocatable :: degen(:,:)

!>--- float data
  real(wp) :: ewin,rthr,bthr,pthr,ethr,athr
  real(wp) :: T,couthr

!>--- boolean data
  logical :: ensembleinput = .false.
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
  logical :: saveelow
  logical :: userinput
  logical :: periodic   !> periodic (fixed-cell) ensemble?

!>--- printout directions
  integer :: prch  !> the main printout channel
  logical :: pr1,pr2,pr3,pr4

!====================================================================!
!>  S E T T I N G S
!====================================================================!
  if (present(quickset)) then
    simpleset = quickset
  else
    simpleset = 0
  end if

!>-- was an actual list of structures (rather than a file name) provided?
  if (present(structurelist)) then
    if (size(structurelist,1) > 0) ensembleinput = .true.
  end if

!>-- determine filenames and output channel
  if (present(infile)) then
    fname = trim(infile)
    userinput = .true.
  else
    fname = trim(env%ensemblename)
    userinput = .false.
  end if
  call cregen_files(env,fname,oname,cname,simpleset,userinput,ensembleinput,prch)

!>-- determine which printouts are required
  call cregen_prout(env,simpleset,pr1,pr2,pr3,pr4)

!>-- determine which subroutines are required
  call cregen_director(env,simpleset,checkbroken,sortE,sortRMSD,sortRMSD2, &
  &  repairord,newfile,conffile,bonusfiles,anal,topocheck,checkez,saveelow)

!>--- DATA SECTION
  call cregen_filldata1(env,ewin,rthr,ethr,bthr,athr,pthr,T,couthr)
  call cregen_filldata2(simpleset,ewin)

!>--- setting the threads for OMP parallel usage
  call cregen_setthreads(prch,env,.false.)

!=====================================================================!
!>  E N S E M B L E   P R O C E S S I N G
!=====================================================================!

!>--- read in the ensemble parameters
  if (.not.ensembleinput) then
    call rdensembleparam(fname,nat,nallref)
  else
    nat = structurelist(1)%nat
    nallref = size(structurelist,1)
  end if

!>--- print a summary about the ensemble and thresholds
  if (pr1) call cregen_pr1(prch,env,nat,nallref,rthr,bthr,pthr,ewin)

!>--- allocate space and read in the ensemble
  if (.not.ensembleinput) then
    call rdensemble(fname,nallref,structures)
  else
    call move_alloc(structurelist,structures)
  end if

!>--- decide periodic vs molecular ensemble (majority vote on lattices).
!>    Under PBC the alignment-based comparisons (axis/rotational constants,
!>    superposition RMSD, rotamers) are ill-defined, so we switch to the
!>    alignment-free fingerprint path. A gas-phase ensemble, conversely,
!>    discards any stray frames that carry a lattice.
  periodic = cregen_majority_periodic(structures)
  if (periodic) then
    write (prch,'(" CREGEN periodic (fixed-cell) mode",t35,":",a10)') ' enabled'
    checkez = .false.   !> C=C E/Z check is not used under PBC
    anal = .false.      !> nuclear-equivalency analysis relies on alignment
  else
    call cregen_discard_periodic(prch,structures)
    nallref = size(structures,1)
  end if

!>--- track ensemble for restart
  !call trackensemble(fname,nat,nallref,at,xyz,comments)

!> NOTE: We check topology and broken structures FIRST before
!> sorting by the energy an making a cut, because chemical changes
!> may produce isomers that are lower in energy at the given
!> level of theory. We do not want that when looking for conformers specifically.

!>--- check if the ensemble contains broken structures? i.e., fusion or dissociation
  if (checkbroken) then
    call cregen_discardbroken(prch,env,topocheck,structures,nall)
  else
    nall = nallref
  end if

!>--- compare neighbourlists to sort out chemically transformed structures
  if (topocheck) then
    call cregen_topocheck(prch,env,checkez,structures,nallnew)
    nall = nallnew !> update
!>--- if structures were discarded, resize xyz
  end if
  if (topocheck.or.checkbroken) then
    write (prch,'(" number of reliable points",t35,":",i10)') nall
  end if

!>--- sort the ensemble by its energies and make a cut (EWIN)
  if (sortE) then
    call cregen_esort(prch,structures,nallnew,ewin=ewin)
    nall = nallnew !> update
  end if

!>--- classify duplicates: periodic fingerprint OR rot.const/RMSD
  if (periodic) then
!>--- alignment-free duplicate detection (energy + cell volume + SVD fingerprint)
    call cregen_CRE_periodic(env,nall,structures,group,rthr, &
    &                        ethr/autokcal,printlvl=2,ch=prch)
    ng = group(0)
    allocate (degen(3,ng))
    call cregen_groupinfo(nall,ng,group,degen)
  else if (sortRMSD) then
    call cregen_CRE_new(env,nall,structures,group,rthr, &
    &                   ethr/autokcal,bthr,printlvl=2,ch=prch)
!>--- get group info to degen
    ng = group(0)
    allocate (degen(3,ng))
    call cregen_groupinfo(nall,ng,group,degen)
  else
    ng = nall
    if (ng > 0) then
      allocate (degen(3,ng))
      do i = 1,ng
        degen(1,i) = 1
        degen(2,i) = i
        degen(3,i) = i
      end do
    else
      allocate (degen(3,1))
      degen = 0
    end if
  end if

!=====================================================================!
!>  E N S E M B L E   O U T P U T
!=====================================================================!

!>--- align all structures to the first structure using the RMSD
!>    (skipped under PBC, where superposition is ill-defined)
  if (.not.periodic) call cregen_rmsdalign(nall,structures)

!>--- write new file with ALL remaining structures
  if (newfile) then
    call cregen_file_wr(env,oname,structures)
!>--- track ensemble for restart
!    call trackensemble(oname,nat,nall,at,xyz,comments)
  end if
!>--- write a file containing only conformers (no rotamers)
  if (conffile) then
    call cregen_conffile(env,cname,structures,ng,degen)
  end if
  if (saveelow) then
    env%elowest = structures(1)%energy
!>-- and update reference geometry (in Bohr)
    env%ref%xyz = structures(1)%xyz
  end if

!>--- additional files for entropy mode
  if (bonusfiles) then
    call cregen_bonusfiles(ng,degen)
  end if

!>--- several printouts
  if (pr2.or.pr3.or.pr4) then
    allocate (er(nall))
    do ii = 1,nall
      er(ii) = structures(ii)%energy
    end do
  end if

  if (pr2) then
    call cregen_pr2(prch,env,nall,ng,degen,er)
    call cregen_econf_list(prch,nall,er,ng,degen)
  end if
  if (pr3) then !> alternative to pr2
    call cregen_pr3(prch,oname,nall,er)
  end if
  if (pr4) then !> group data printout
    call cregen_pr4(prch,fname,nall,group)
  end if

!>--- analyze nuclear equivalencies, e.g. for NMR and Entropy
  if (anal) then
    call cregen_EQUAL(prch,nall,structures,group,athr,.not.env%entropic)
  end if

!>-- in case we had a structurelist given, move the (sorted) memory space back there
  if (ensembleinput) then
    call move_alloc(structures,structurelist)
  end if

  if (newfile) then
    write (prch,'(a,a)') 'Full ensemble file written to:    ',trim(oname)
  end if
  if (conffile) then
    write (prch,'(a,a)') 'Unique-structure file written to: ',trim(cname)
  end if

!>--- deallocate data
  if (prch .ne. stdout) then
    close (prch)
  end if
  if (allocated(er)) deallocate (er)
  if (allocated(degen)) deallocate (degen)
  if (allocated(group)) deallocate (group)
  return
end subroutine newcregen

!=========================================================================================!
!=========================================================================================!
!>  CREGEN DATA SECTION
!=========================================================================================!
!=========================================================================================!

subroutine cregen_files(env,fname,oname,cname,simpleset,userinput,ensembleinput,iounit)
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
  type(systemdata),intent(inout) :: env    !> MAIN STORAGE OS SYSTEM DATA
  character(len=:),allocatable,intent(inout) :: fname  !> name of the ensemble to be read
  character(len=:),allocatable,intent(inout) :: oname  !> output ensemble name (including rotamers)
  character(len=:),allocatable,intent(inout) :: cname  !> output ensemble name (only conformers)
  integer,intent(in) :: simpleset
  logical,intent(in) :: userinput !> was an input file given via the optional subroutine arg?
  logical,intent(in) :: ensembleinput !> was an structure list provided?
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
    case (6,7,9,12)
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
  &    .OR.userinput.OR.ensembleinput) then
    if (.not.userinput.and..not.ensembleinput) then
      fname = trim(env%ensemblename)
    end if
    cname = 'crest_ensemble.xyz'
    oname = trim(fname)//'.sorted'
    if (env%fullcre) then
      env%ensemblename = trim(oname)
    end if
  else !> internal mode for conformational search
    fname = repeat(' ',256)  !> need initialization because checkname_xyz
    oname = repeat(' ',256)  !> can't handle allocatable names
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

  write (iounit,'(1x,a,a)') 'input  file name : ',trim(fname)
  select case (simpleset)
  case (9)
    continue
  case default
    write (iounit,'(1x,a,a)') 'output file name : ',trim(oname)
  end select

  inquire (file=fname,exist=ex)
  if (.not.ex.and..not.ensembleinput) then
    write (stdout,'(a,a,a)') 'CREGEN> **WARNING** file ',trim(fname),' does not exist!'
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

  if (any(simpleset == (/6,7/)).or.env%esort) then
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
        &  repairord,newfile,conffile,bonusfiles,anal,topocheck,checkez,saveelow)
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
  logical,intent(out) :: saveelow

  checkbroken = .true. !> fragmentized structures are sorted out
  sortE = .true.       !> sort based on energy
  sortRMSD = .true.    !> sort based on RMSD
  sortRMSD2 = .false.  !> check groups for whole ensemble
  repairord = .true.   !> double-check the sorted Ensemble

  newfile = .true.  !> sorted input file

  conffile = .true. !> sorted unique structure file

  saveelow = .true. !> save (overwrite) lowest structure to env%ref

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

  if (any(simpleset == (/6,7/)).or.env%esort) then  !energy sorting only
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
    saveelow = .false.
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

  if (simpleset == 13) then  !msreact mode
    checkbroken = .false.
    sorte = .true.
    sortRMSD = .true.
    sortRMSD2 = .false.
    repairord = .false.
    newfile = .true.
    conffile = .true.
    topocheck = .false.
    checkez = .false.
    bonusfiles = .false.
    anal = .false.

  end if

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
  if (any(simpleset == (/6,12/))) then
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
  implicit none
  integer :: nall,ng
  integer :: group(0:nall)
  integer :: degen(3,ng)
  integer :: i,j,k,a,b
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

subroutine cregen_discardbroken(ch,env,topocheck,structures,newnall)
!**************************************************
!* subroutine cregen_discardbroken
!* analyze an ensemble and track broken structures
!* to be discarded.
!**************************************************
  use crest_parameters
  use crest_data
  use strucrd
  use adjacency
  use cregen_utils
  implicit none
  !> INPUT
  type(systemdata),intent(in) :: env    ! MAIN STORAGE OS SYSTEM DATA
  integer,intent(in) :: ch ! printout channel
  logical,intent(in) :: topocheck
  type(coord),intent(inout),allocatable,target :: structures(:)
  integer,intent(out) :: newnall
  !> LOCAL
  integer :: llan,nall,frag,frag0
  real(wp) :: erj,cnorm
  integer :: ii,jj
  logical :: substruc
  logical :: dissoc,distok,distok2
  logical,allocatable :: broke(:)
  type(coord) :: mol0
  type(coord),pointer :: mol
  type(coord),allocatable :: tmpstructures(:)

  !>--- if we don't wish to include all atoms:
  substruc = (structures(1)%nat .ne. env%rednat.and.env%subRMSD)
  nall = size(structures,1)
  !> Check fragments
  call env%ref%to(mol0)
  call cregen_calculate_fragments(mol0,nfrag=frag0)
  write (ch,'(" # fragment in coord",t35,":",i10)') frag0

  !>--- loop over the structures
  allocate (broke(nall),source=.false.)
  newnall = 0
  llan = nall
  do ii = 1,nall
    mol => structures(ii)
    erj = mol%energy
    !if (substruc) then
    !  !...
    !end if

    !>--- close contact checks
    cnorm = sum(abs(mol%xyz))           !> clash check
    distok = distcheck(mol%nat,mol%xyz) !> distance check
    distok2 = (cnorm .gt. 1.0d-6)

    !>--- further checks: dissociation?
    dissoc = .false.
    if (abs(erj) .gt. 1.0d-6.and. &
    &   distok.and.distok2.and.topocheck) then
      call cregen_calculate_fragments(mol,nfrag=frag)
      dissoc = (frag .gt. frag0)
    end if

    if (dissoc.or.(.not.distok).or.(.not.distok2)) then
      !>--- move broken structures to the end of the matrix
      broke(ii) = .true.
      !write(ch,*) 'removing structure',ii
    else
      newnall = newnall+1
    end if
  end do

  !>--- sort the xyz array (only if structures have been discarded)
  if (newnall .lt. nall) then
    allocate (tmpstructures(newnall))
    jj = 0
    do ii = 1,nall
      if (.not.broke(ii)) then
        jj = jj+1
        tmpstructures(jj) = structures(ii)
      end if
    end do
    call move_alloc(tmpstructures,structures)
    llan = nall-newnall
    write (ch,'(" number of removed clashes",t35,":",i10)') llan
  end if
  !>--- otherwise the ensemble is ok
  if (allocated(broke)) deallocate (broke)
  return
end subroutine cregen_discardbroken

!=========================================================================================!

subroutine cregen_topocheck(ch,env,checkez,structures,newnall)
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
  use crest_cn_module
  use quicksort_interface
  use cregen_utils
  implicit none
  type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
  integer,intent(in) :: ch ! printout channel
  logical,intent(in) :: checkez
  type(coord),intent(inout),allocatable,target :: structures(:)
  integer,intent(out) :: newnall
  integer :: nat,nall,llan
  real(wp),allocatable :: cn(:),bond(:,:)
  integer,allocatable :: toporef(:)
  integer,allocatable :: topo(:)
  logical,allocatable :: neighmat(:,:)
  integer :: nbonds
  integer :: ii,jj,l
  integer :: ntopo,ncc,ccfail
  logical :: discard
  integer,allocatable :: ezat(:,:)
  real(wp),allocatable :: ezdihedref(:)
  real(wp),allocatable :: ezdihed(:)
  real(wp) :: winkeldiff

  type(coord) :: mol0
  type(coord),pointer :: mol
  type(coord),allocatable :: tmpstructures(:)
  logical,allocatable :: broke(:)

  !>--- read the reference structure
  call env%ref%to(mol0)
  nat = mol0%nat
  nall = size(structures,1)
  call mol0%cn_to_bond(cn,bond)
  !>--- calculate reference "topology"
  if (allocated(env%excludeTOPO)) then
    call bondtotopo(nat,mol0%at,bond,cn,ntopo,toporef,neighmat,excl=env%excludeTOPO)
  else
    call bondtotopo(nat,mol0%at,bond,cn,ntopo,toporef,neighmat)
  end if

  nbonds = sum(toporef)
  write (ch,'(" # bonds in reference structure",t35,":",i10)') nbonds
  !>--- if required, check for C=C bonds (based only on structure!)
  if (checkez) then
    call nezcc(nat,mol0%at,mol0%xyz,cn,ntopo,toporef,ncc)
    if (ncc > 0) then
      write (ch,'("  => # of C=C bonds : ",i0)') ncc
      allocate (ezat(4,ncc))
      allocate (ezdihedref(ncc),ezdihed(ncc),source=0.0d0)
      call ezccat(nat,mol0%at,mol0%xyz,cn,ntopo,toporef,ncc,ezat)
      call ezccdihed(nat,mol0%xyz,ncc,ezat,ezdihedref)
      !do i=1,ncc
      !  write(*,'(1x,a,4i4,a,f6.2)') 'C=C bond atoms:',ezat(1:4,i)," angle: ",ezdihedref(i)
      !enddo
    end if
  end if

  allocate (broke(nall),source=.false.)
  !>--- loop over the structures
  ccfail = 0
  newnall = 0
  llan = nall
  do jj = 1,nall
    !>--- generate topo and compare
    discard = .false.
    mol => structures(jj)
    call mol%cn_to_bond(cn,bond)
    if (allocated(env%excludeTOPO)) then
      call bondtotopo(mol%nat,mol%at,bond,cn,ntopo,topo,neighmat,excl=env%excludeTOPO)
    else
      call bondtotopo(mol%nat,mol%at,bond,cn,ntopo,topo,neighmat)
    end if
    do l = 1,ntopo
      if (toporef(l) .ne. topo(l)) then
        discard = .true.   !> if there is any mismatch in neighbor lists
        exit
      end if
    end do
    !>--- get E/Z info of C=C, discard isomers
    if (checkez.and..not.discard.and.ncc > 0) then
      call ezccdihed(mol%nat,mol%xyz,ncc,ezat,ezdihed)
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
      broke(jj) = .true.
    else
      newnall = newnall+1
    end if
  end do

  !>--- sort the xyz array (only if structures have been discarded)
  if (newnall .lt. nall) then
    llan = nall-newnall
    write (ch,'(" number of topology mismatches",t35,":",i10)') llan
    !>--- report the removals during a run
    if (ch .ne. stdout) then
      write (stdout,'("CREGEN> number of topology-based structure removals: ",i0)') llan
    end if
    if (checkez.and.ccfail > 0) then
      write (ch,'(''  => discared due to E/Z isom.  : '',i0)') ccfail
    end if
    if (newnall >= 1) then
      allocate (tmpstructures(newnall))
      jj = 0
      do ii = 1,nall
        if (.not.broke(ii)) then
          jj = jj+1
          tmpstructures(jj) = structures(ii)
        end if
      end do
      call move_alloc(tmpstructures,structures)
    else
      if (ch .ne. stdout) then
        write (stdout,'("CREGEN> ** WARNING ** Full removal of ensemble! Falling back to reference structure.")')
      end if
      allocate (tmpstructures(1),source=mol0)
      call move_alloc(tmpstructures,structures)
    end if
  end if
  !>--- otherwise the ensemble is ok
  if (allocated(ezdihedref)) deallocate (ezdihedref)
  if (allocated(ezdihed)) deallocate (ezdihed)
  if (allocated(ezat)) deallocate (ezat)
  deallocate (cn,bond)
  deallocate (neighmat)
  deallocate (topo,toporef)
  return
end subroutine cregen_topocheck

!=========================================================================================!

subroutine cregen_esort(ch,structures,nallout,ewin)
!**************************************************************
!* subroutine cregen_esort
!* sort the ensemble by energy and determine the new
!* ensemble size within the energy threshold.
!* On Input: ch - printout channel
!*           structures - the list of structures
!*           nallout - number of surviving structures
!*           ewin - energy window in kcal/mol
!* On Output: nallout - number of strucutres after cutoff
!**************************************************************
  use crest_parameters
  use strucrd
  use quicksort_interface
  implicit none
  integer,intent(in) :: ch
  type(coord),intent(inout),allocatable :: structures(:)
  integer,intent(out) :: nallout
  real(wp),intent(in),optional :: ewin
  integer :: nall,nat

  real(wp),allocatable :: energies(:)
  type(coord),allocatable :: tmpstructures(:)
  integer :: ii,jj
  real(wp) :: de,emax,frac

  nall = size(structures,1)
  nallout = nall
  call ensemble_qsort(nall,structures,1,nall)

  !>-- determine cut-off of energies (optional)
  if (present(ewin)) then
    write (ch,'(80("*"))')
    allocate (energies(nall))
    do ii = 1,nall
      energies(ii) = structures(ii)%energy
    end do

    if (ewin < 9999.9_wp) then
      write (ch,'(" sorting energy window (EWIN)",t32,":",1x,f9.4,a)') ewin,' / kcal/mol'
    else
      write (ch,'(" sorting energy window (EWIN)",t32,":",3x,a,a)') '+∞',' / kcal/mol'
    end if
    emax = maxval(energies(:),1)
    de = (emax-energies(1))*autokcal
    if (de .gt. ewin) then
      nallout = 1 !> lowest is always taken
      do ii = 2,nall
        de = (energies(ii)-energies(1))*autokcal
        if (de .lt. ewin) then
          nallout = nallout+1
        else
          exit
        end if
      end do
      frac = real(nall-nallout,wp)/real(nall,wp)
      write (ch,'(" number of removed by energy",t32,":",3x,i10,a,f6.2,a)') &
      &       (nall-nallout),' (',frac*100.d0,'%)'
      write (ch,'(" number of remaining points",t32,":",3x,i10,a,f6.2,a)') &
      &       nallout,' (', (1.0d0-frac)*100.d0,'%)'

      allocate (tmpstructures(nallout))
      do ii = 1,nallout
        tmpstructures(ii) = structures(ii)
      end do
      call move_alloc(tmpstructures,structures)
    else
      nallout = nall
    end if
    write (ch,'(" reference state Etot",t32,":",2x,es14.6)') energies(1)
    deallocate (energies)
  end if

  return
end subroutine cregen_esort

!=========================================================================================!
!> The actual core of CREGEN: Conformer/Rotamer Classification
!=========================================================================================!
subroutine cregen_CRE_new(env,nall,structures,groups,rthresh,ethr,bthr, &
    &                         printlvl,ch)
!**************************************************************************************
!* Re-ractored implementaiton of the original CREGEN workflow, classifying conformers
!* according to their quaternion RMSD, rotational constants, energy,
!* and pair-distance sum
!* "structures" will be updated so that all true duplicates are pruned and
!* molecules are clustered by their group. "group" and "nall" are also updated.
!*
!* Input arguments:
!*         env - CREST systemdata
!*        nall - total number of structures
!*  structures - the structures
!*      groups - group assignment for each structure (dimension 0:nall), group(0) = maxval(group(1:nall))
!*      rthresh - RMSD threshold (in ANGSTRÖM) for conformer distinction
!*        ethr - inter-conformer energy threshold (in HARTREE) for pre-sorting
!*        bthr - rotational constant similarity threshold (percentage based)
!*
!* Optionals:
!*    printlvl - integer to direct the print verbosity. (0=minimal, 1=verbose)
!*          ch - integer for print channel
!*
!* Output:
!*      groups - integer array assigning each structure to a group
!*************************************************************************************
  use crest_parameters
  use crest_data
  use rotaniso_mod
  use axis_module
  use strucrd
  use canonical_mod
  use irmsd_module
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env
  integer,intent(inout) :: nall
  type(coord),intent(inout),allocatable,target :: structures(:)
  integer,intent(out),allocatable :: groups(:)
  real(wp),intent(in) :: RTHRESH
  real(wp),intent(in) :: ETHR
  real(wp),intent(in) :: BTHR
  integer,intent(in),optional :: printlvl
  integer,intent(in),optional :: ch

  !> LOCAL
  integer :: i,ii,jj,kk,T,cc,nat,io,gg
  integer :: gcount,ggcount,nallnew
  integer :: prlvl,prch
  type(rmsd_cache),allocatable :: rcaches(:)
  type(coord),allocatable,target :: workmols(:)
  type(canonical_sorter),allocatable :: sorters(:)
  type(coord),pointer :: ref,mol
  real(wp) :: rmsdval,RTHR,ediff,eii,avmom,rsq,frac
  real(wp),allocatable :: rot(:,:)
  real(wp),allocatable :: xyz_aa(:,:)
  integer,allocatable :: prune_table(:)
  real(wp),allocatable :: enuc(:)
  logical :: l1,l2
  character(len=:),allocatable :: tmpstr
  logical :: heavy,substruc
  logical,allocatable :: mask(:)
  integer,allocatable :: tmpgroups(:),double(:)
  type(coord),allocatable :: tmpstructures(:)
  !type(progress_state) :: ps

  !> defaults that are practically never touched
  real(wp),parameter :: bthrmax = 0.025_wp
  real(wp),parameter :: bthrshift = 0.5_wp
  real(wp),parameter :: enuc_thr = 1.0d-3

  logical,parameter :: debug = .false.

!>--- handle optional arguments
  if (present(printlvl)) then
    prlvl = printlvl
  else
    prlvl = 1
  end if

  if (present(ch)) then
    prch = ch
  else
    prch = stdout
  end if

!>--- set up parallelization
!     ...
  T = 1 !> doing it serial for now

!>--- set up parameters (NOTE, we are working with BOHR internally)
  RTHR = RTHRESH*aatoau

!>--- reference structure (the first one) for some setup
  ref => structures(1)
  nat = ref%nat

!>--- print some sorting data
  if (prlvl > 0) then
    tmpstr = 'Info for CREGEN sorting:'
    if (prlvl > 1.and.prch == stdout) then
      !  call printc(style(S_BOLD)//fg(YELLOW,bright=.true.)//trim(tmpstr)//reset())
    else
      write (prch,'(a)') 'Info for CREGEN sorting:'
    end if
    !write (prch,'(2x,a,i10)') 'number of structures     :',nall
    write (prch,'(2x,a,t32,a,f10.5,a)') 'RTHR (RMSD threshold)',':',RTHR*autoaa,' Å'
    write (prch,'(2x,a,t32,a,es10.2,a)') 'ETHR (energy threshold)',':',ETHR,' Ha'
    write (prch,'(2x,a,t32,a,f10.2,a)') 'BTHR (rot. threshold)',':',BTHR*100,' %'
    !write (prch,'(2x,a,i9)') 'OpenMP threads           :',T
  end if

!>--- mask setup: We may not include all atoms in the checks
  heavy = env%heavyrmsd
  substruc = (nat .ne. env%rednat.and.env%subRMSD.and.allocated(env%includeRMSD))
  if (heavy.or.substruc) then
    allocate (mask(nat),source=.false.)
  end if
  if (heavy) then
    do ii = 1,nat
      if (structures(1)%at(ii) .ne. 1) mask(ii) = .true.
    end do
  end if
  if (substruc) then
    do ii = 1,nat
      mask(ii) = (env%includeRMSD(ii) .eq. 1)
    end do
  end if
  if ((heavy.or.substruc).and.(prlvl > 0)) then
    write (prch,'(" Heavy/masked atoms",t32,":",i10," / ",i0)') count(mask),nat
  end if

  if (prlvl > 0) then
    tmpstr = "Starting calculations..."
    if (prlvl > 1.and.prch == stdout) then
      !  call printc(style(S_BOLD)//fg(YELLOW,bright=.true.)//trim(tmpstr)//reset())
    else
      write (stdout,'(a)') trim(tmpstr)
    end if
  end if
!>--- allocate work cache
  if (prlvl > 0) then
    write (prch,'(a)',advance='no') 'Allocating RMSD work cache ... '
    flush (prch)
  end if
  allocate (rcaches(T))
  allocate (workmols(T))
  do i = 1,T
    mol => workmols(i)
    allocate (mol%at(ref%nat))
    allocate (mol%xyz(3,ref%nat))
    nullify (mol)
    call rcaches(i)%allocate(ref%nat)
  end do
  if (prlvl > 0) then
    write (prch,'(a)') 'done.'
  end if

!> ----------------------------------------------
!> PRE-PROCESSING for more efficient sorting
!> ----------------------------------------------
  !> prune_table keeps track of which structure to compare to
  !> so for a  list of structures (1...j...k...nall), the entry
  !> prune_table(k) = j, tells us structure k is compared to all
  !> structures j up to k-1. The table is initialized to 1, so
  !> the full comparison list is used.
  allocate (prune_table(nall),source=1)
  !> conveniently, we can use the energy threshold to set a better
  !> comparison table, as in the original CREGEN routine.
  do ii = 1,nall
    eii = structures(ii)%energy
    do jj = 1,ii
      ediff = abs(eii-structures(jj)%energy)
      if (ediff <= ETHR) then
        prune_table(ii) = jj
        exit
      end if
    end do
  end do

  !> Prepare axis comparison
  !> axis alignment and rotational constant calculation
  allocate (rot(3,nall),source=0.0_wp)
  allocate (xyz_aa(3,nat),source=0.0_wp)
  do ii = 1,nall
    mol => structures(ii)
    call axis(mol%nat,mol%at,mol%xyz) !> all coordinates to CMA
    xyz_aa = mol%xyz*autoaa
    call axis(mol%nat,mol%at,xyz_aa,rot(1:3,ii),avmom) !> B_0 in MHz
  end do

  !> Scaled sum of atom-atom-distances (empirical measure)
  allocate (enuc(nall),source=0.0_wp)
  do ii = 1,nall
    mol => structures(ii)
    do jj = 1,mol%nat-1
      do kk = jj+1,mol%nat
        rsq = (mol%xyz(1,jj)-mol%xyz(1,kk))**2 &
          &  +(mol%xyz(2,jj)-mol%xyz(2,kk))**2 &
          &  +(mol%xyz(3,jj)-mol%xyz(3,kk))**2+1.d-12
        enuc(ii) = enuc(ii)+real(mol%at(jj)*mol%at(kk),wp)/rsq
      end do
    end do
  end do

!> --------------------------------------------
!> pre-processing end
!> --------------------------------------------

!>--- run the checks
  if (prlvl > 0) then
    write (prch,'(a,6x,a)',advance='no') 'Running CREGEN checks','... '
    flush (prch)
    if (prlvl > 1.and.prch == stdout) then
      !  write (stdout,*)
      !  call progress_init(ps,width=50,prefix=" ↳", &
      !    & suffix="",show_time=.true.,show_eta=.false.)
      !  call progress_update(ps,0,nall)
    else
      write (stdout,'(a)',advance='no') 'CREGEN> running RMSDs ...'
      flush (stdout)
    end if
  end if
  allocate (groups(nall),source=0)
  gcount = maxval(groups(:))
  do ii = 1,nall
!>--- find next unassigned conformer and assign a new group
    if (groups(ii) .ne. 0) cycle
    gcount = gcount+1
    groups(ii) = gcount

!>--- Then, cross-check all other unassigned conformers
    cc = 1  !> again, serial implementation for now
    ! !$omp parallel &
    ! !$omp shared(nall, nat, groups, sorters, rcaches, rot) &
    ! !$omp shared(workmols, structures, ii, prune_table,heavy,substruc,mask) &
    ! !$omp private(jj,rmsdval,cc,io, l1, l2)
    ! !$omp do schedule(dynamic)
    do jj = ii+1,nall
      !cc = omp_get_thread_num()+1
      if (groups(jj) .ne. 0) cycle
      if (ii < prune_table(jj)) cycle
      workmols(cc)%nat = structures(jj)%nat
      workmols(cc)%at(:) = structures(jj)%at(:)
      workmols(cc)%xyz(:,:) = structures(jj)%xyz(:,:)
      !if (heavy.or.substruc) then
      !  rmsdval = rmsd(structures(ii),workmols(cc),mask=mask,&
      !    &       scratch=rcaches(cc)%xyzscratch,ccache=rcaches(cc)%ccache)
      !else
      !  rmsdval = rmsd(structures(ii),workmols(cc), &
      !    &       scratch=rcaches(cc)%xyzscratch,ccache=rcaches(cc)%ccache)
      !end if
      !if (rmsdval < RTHR) then
      !> only "true" duplicates will have tiny RMSD, assign negative gcount for pruning
      !  groups(jj) = -gcount
      !else
      l1 = equalrotaniso(ii,jj,nall,rot,BTHR,bthrmax,bthrshift)
      l2 = (2.0_wp*abs(enuc(ii)-enuc(jj))/(enuc(ii)+enuc(jj))) .lt. enuc_thr
      if (l1.and.l2) groups(jj) = gcount
      !end if
    end do
    if (prlvl > 1) then
      !  call progress_update(ps,ii,nall)
    end if
    ! !$omp end do
    ! !$omp end parallel
  end do

!> for all groups run RMSD checks
  gcount = maxval(groups(1:nall))
  do gg = 1,gcount
    do ii = 1,nall
      if (groups(ii) .ne. gg) cycle
      do jj = ii+1,nall
        kk = groups(jj)
        if (kk .ne. gg.or.kk < 0) cycle

        workmols(cc)%nat = structures(jj)%nat
        workmols(cc)%at(:) = structures(jj)%at(:)
        workmols(cc)%xyz(:,:) = structures(jj)%xyz(:,:)
        if (heavy.or.substruc) then
          rmsdval = rmsd(structures(ii),workmols(cc),mask=mask,&
            &       scratch=rcaches(cc)%xyzscratch,ccache=rcaches(cc)%ccache)
        else
          rmsdval = rmsd(structures(ii),workmols(cc), &
            &       scratch=rcaches(cc)%xyzscratch,ccache=rcaches(cc)%ccache)
        end if
        if (rmsdval < RTHR) then
          !> only "true" duplicates will have tiny RMSD, assign negative gcount for pruning
          groups(jj) = -gg
        end if
      end do
    end do
  end do

  if (prlvl > 0) then
    !if (prlvl > 1 .and.prch == stdout) then
    if (prlvl > 1 .and.prch == stdout) then
      !  call progress_update(ps,nall,nall)
      !  call progress_finish(ps)
      write (prch,'(a)') 'done.'
    else
      write (stdout,'(a)') ' done.'
    end if
  end if

!> finally, resizing the ensemble with remaining unique conformers+rotamers
!> Note, structures are grouped by assigned group and within group are ordered
!> with increasing energy (because the initial ensemble was energy-sorted)
  if (prlvl > 0) then
    write (prch,'(a,6x,a)',advance='no') 'Discarding duplicates','...'
    flush (prch)
  end if
  gcount = maxval(groups(1:nall))
  nallnew = count(groups(1:nall) > 0)
  allocate (tmpstructures(nallnew))
  allocate (tmpgroups(0:nallnew),source=0)
  allocate (double(nallnew),source=0)
  cc = 0
  do ii = 1,gcount
    do jj = 1,nall
      ggcount = groups(jj)
      if (ggcount .eq. ii.and.ggcount > 0) then
        cc = cc+1
        tmpstructures(cc) = structures(jj)
        tmpgroups(cc) = ggcount
        do kk = 1,cc-1
          if (tmpgroups(kk) .eq. ggcount) then
            double(cc) = kk
            exit
          end if
        end do
      end if
    end do
  end do
  tmpgroups(0) = gcount
  call move_alloc(tmpgroups,groups)
  call move_alloc(tmpstructures,structures)
  if (prlvl > 0) then
    write (prch,'(a)') ' done.'
    frac = real(nall-nallnew,wp)/real(nall,wp)
    write (prch,'(1x,a,t40,a,i10,a,f6.2,a)') &
    &      "number of doubles removed by rot/RMSD",":",nall-nallnew,' (',frac*100.d0,'%)'
    write (prch,'(1x,a,t40,a,i10,a,f6.2,a)') &
    &      "number of unique structures remaining",":",nallnew,' (', (1.0d0-frac)*100.d0,'%)'
    frac = real(gcount,wp)/real(nallnew,wp)
    write (prch,'(1x,a,t40,a,i10,a,f6.2,a,i0,a)') &
    &      "number of unique conformers identified",":",gcount,' (', (frac)*100.d0,'% of ',nallnew,')'
  end if
  nall = nallnew

  !>-- for ENSO write a file with duplicate info (if required)
  call enso_duplicates(env,nall,double)

  if (allocated(prune_table)) deallocate (prune_table)
  if (allocated(mask)) deallocate (mask)
  if (allocated(enuc)) deallocate (enuc)
  if (allocated(rot)) deallocate (rot)
  if (allocated(prune_table)) deallocate (prune_table)
end subroutine cregen_CRE_new

!=========================================================================================!

function cregen_majority_periodic(structures) result(periodic)
!**************************************************************
!* Decide whether an ensemble is periodic (fixed-cell) or
!* molecular (gas-phase) by a simple majority vote: a structure
!* counts as periodic if it carries an allocated, non-zero
!* lattice. Returns .true. if at least half of the structures
!* are periodic.
!**************************************************************
  use crest_parameters
  use strucrd
  implicit none
  type(coord),intent(in) :: structures(:)
  logical :: periodic
  integer :: ii,nall,nper
  nall = size(structures,1)
  nper = 0
  do ii = 1,nall
    if (allocated(structures(ii)%lat)) then
      if (any(abs(structures(ii)%lat) > 1.0e-8_wp)) nper = nper+1
    end if
  end do
  periodic = (nall > 0).and. (2*nper >= nall)
end function cregen_majority_periodic

!=========================================================================================!

subroutine cregen_discard_periodic(ch,structures)
!**************************************************************
!* Remove any structures that carry a (non-zero) lattice from
!* a molecular (gas-phase) ensemble. Stray PBC frames would
!* otherwise contaminate the alignment-based comparisons.
!**************************************************************
  use crest_parameters
  use strucrd
  implicit none
  integer,intent(in) :: ch
  type(coord),intent(inout),allocatable :: structures(:)
  integer :: ii,jj,nall,nkeep
  logical,allocatable :: drop(:)
  type(coord),allocatable :: tmpstructures(:)
  nall = size(structures,1)
  allocate (drop(nall),source=.false.)
  do ii = 1,nall
    if (allocated(structures(ii)%lat)) then
      if (any(abs(structures(ii)%lat) > 1.0e-8_wp)) drop(ii) = .true.
    end if
  end do
  nkeep = count(.not.drop)
  if (nkeep < nall) then
    write (ch,'(" number of discarded PBC frames",t35,":",i10)') nall-nkeep
    allocate (tmpstructures(nkeep))
    jj = 0
    do ii = 1,nall
      if (.not.drop(ii)) then
        jj = jj+1
        tmpstructures(jj) = structures(ii)
      end if
    end do
    call move_alloc(tmpstructures,structures)
  end if
  deallocate (drop)
end subroutine cregen_discard_periodic

!=========================================================================================!

subroutine cregen_CRE_periodic(env,nall,structures,groups,rthresh,ethr,printlvl,ch)
!**************************************************************************************
!* Duplicate/rotamer classification for periodic (fixed-cell) ensembles.
!*
!* Two-stage scheme mirroring the molecular cregen_CRE_new, but with periodic-safe
!* descriptors:
!*   Stage 1 (grouping): candidate structures are pre-grouped by an alignment-free
!*     SVD "fingerprint" of the symmetric minimum-image distance matrix (Σ, the
!*     singular values, cf. Pracht/Morgan/Wales, J. Chem. Phys. 159, 064801 (2023))
!*     together with the cell volume. The fingerprint replaces the rotational
!*     constants of the molecular code as the cheap rotamer pre-filter, using a
!*     *loose* threshold (SIGTHR_GROUP) so that all rotamers of a conformer land
!*     in one group.
!*   Stage 2 (pruning): within each group, true duplicates are removed by the
!*     minimum-image RMSD (MIC-unwrap + Kabsch, via irmsd_module::rmsd, which is
!*     automatically periodic when the structures carry a lattice). RTHR is the
!*     authoritative duplicate threshold, exactly as in the molecular case.
!*
!* Input arguments:
!*         env - CREST systemdata
!*        nall - total number of structures (updated on output)
!*  structures - the (energy-sorted) structures; pruned/clustered on output
!*     rthresh - RMSD duplicate threshold (ANGSTRÖM)
!*        ethr - inter-structure energy threshold (HARTREE) for the comparison window
!* Optionals:
!*    printlvl - print verbosity (0=minimal, 1=verbose)
!*          ch - print channel
!* Output:
!*      groups - group assignment (dimension 0:nall), groups(0) = number of groups
!*************************************************************************************
  use crest_parameters
  use crest_data
  use strucrd
  use pbc_fingerprint_module
  use irmsd_module,only:rmsd,rmsd_cache
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env
  integer,intent(inout) :: nall
  type(coord),intent(inout),allocatable,target :: structures(:)
  integer,intent(out),allocatable :: groups(:)
  real(wp),intent(in) :: RTHRESH
  real(wp),intent(in) :: ETHR
  integer,intent(in),optional :: printlvl
  integer,intent(in),optional :: ch
  !> LOCAL
  integer :: ii,jj,kk,nat,gcount,ggcount,nallnew,cc
  integer :: prlvl,prch
  real(wp) :: eii,ediff,frac,RTHR,rmsdval
  integer,allocatable :: prune_table(:)
  real(wp),allocatable :: sig(:,:),sigtmp(:),vol(:)
  integer,allocatable :: tmpgroups(:),double(:)
  type(coord),allocatable :: tmpstructures(:)
  logical :: heavy,substruc
  logical,allocatable :: mask(:)
  type(rmsd_cache) :: rcache

!>--- handle optional arguments
  if (present(printlvl)) then
    prlvl = printlvl
  else
    prlvl = 1
  end if
  if (present(ch)) then
    prch = ch
  else
    prch = stdout
  end if

  nat = structures(1)%nat
  RTHR = RTHRESH*aatoau  !> RMSD threshold to Bohr (internal units)

  if (prlvl > 0) then
    write (prch,'(a)') 'Info for periodic CREGEN sorting:'
    write (prch,'(2x,a,t32,a,f10.5,a)') 'RTHR (MIC-RMSD threshold)',':',RTHR*autoaa,' Å'
    write (prch,'(2x,a,t32,a,es10.2,a)') 'ETHR (energy threshold)',':',ETHR,' Ha'
    write (prch,'(2x,a,t32,a,es10.2,a)') 'VOLTHR (cell volume)',':',volthr*autoaa**3,' Å³'
    write (prch,'(2x,a,t32,a,es10.2,a)') 'SIGTHR (fingerprint group)',':',sigthr_group,' (loose)'
  end if

!>--- mask setup: we may not include all atoms in the RMSD check
  heavy = env%heavyrmsd
  substruc = (nat .ne. env%rednat.and.env%subRMSD.and.allocated(env%includeRMSD))
  if (heavy.or.substruc) allocate (mask(nat),source=.false.)
  if (heavy) then
    do ii = 1,nat
      if (structures(1)%at(ii) .ne. 1) mask(ii) = .true.
    end do
  end if
  if (substruc) then
    do ii = 1,nat
      mask(ii) = (env%includeRMSD(ii) .eq. 1)
    end do
  end if
  call rcache%allocate(nat)

!> ----------------------------------------------
!> PRE-PROCESSING
!> ----------------------------------------------
!>--- energy-based comparison window (structures are energy-sorted already)
  allocate (prune_table(nall),source=1)
  do ii = 1,nall
    eii = structures(ii)%energy
    do jj = 1,ii
      ediff = abs(eii-structures(jj)%energy)
      if (ediff <= ETHR) then
        prune_table(ii) = jj
        exit
      end if
    end do
  end do

!>--- precompute fingerprints (singular values) and cell volumes
  if (prlvl > 0) then
    write (prch,'(a)',advance='no') 'Computing periodic fingerprints ... '
    flush (prch)
  end if
  allocate (sig(nat,nall),source=0.0_wp)
  allocate (vol(nall),source=0.0_wp)
  do ii = 1,nall
    call pbc_fingerprint(structures(ii),sigtmp)
    sig(1:nat,ii) = sigtmp(1:nat)
    deallocate (sigtmp)
    vol(ii) = structures(ii)%cellvol()
  end do
  if (prlvl > 0) write (prch,'(a)') 'done.'

!> ----------------------------------------------
!> CLASSIFICATION: greedy seed loop. Each unassigned structure becomes a new
!> unique representative; every later in-window structure is tested against it
!> by minimum-image RMSD and pruned if within RTHR. The fingerprint is a *loose*
!> cheap skip-filter (energy + cell volume + Σ): a pair whose fingerprints differ
!> by more than SIGTHR_GROUP cannot be a duplicate, so the (more expensive) RMSD
!> is skipped. It restricts nothing — unlike a hard grouping, it never prevents
!> a genuine duplicate from being compared.
!> ----------------------------------------------
  if (prlvl > 0) then
    write (prch,'(a)',advance='no') 'Running minimum-image RMSDs ... '
    flush (prch)
  end if
  allocate (groups(nall),source=0)
  gcount = 0
  do ii = 1,nall
    if (groups(ii) .ne. 0) cycle
    gcount = gcount+1
    groups(ii) = gcount             !> new unique representative
    do jj = ii+1,nall
      if (groups(jj) .ne. 0) cycle
      if (ii < prune_table(jj)) cycle                       !> energy window
      if (abs(vol(ii)-vol(jj)) >= volthr) cycle             !> cell volume
      if (fp_distance(sig(:,ii),sig(:,jj)) >= sigthr_group) cycle  !> fingerprint skip
      if (heavy.or.substruc) then
        rmsdval = rmsd(structures(ii),structures(jj),mask=mask, &
        &              scratch=rcache%xyzscratch,ccache=rcache%ccache)
      else
        rmsdval = rmsd(structures(ii),structures(jj), &
        &              scratch=rcache%xyzscratch,ccache=rcache%ccache)
      end if
      if (rmsdval < RTHR) groups(jj) = -gcount  !> duplicate of seed ii -> prune
    end do
  end do
  if (prlvl > 0) write (prch,'(a)') 'done.'

!> ----------------------------------------------
!> resize the ensemble keeping only unique structures (groups > 0)
!> ----------------------------------------------
  if (prlvl > 0) then
    write (prch,'(a,6x,a)',advance='no') 'Discarding duplicates','...'
    flush (prch)
  end if
  gcount = maxval(groups(1:nall))
  nallnew = count(groups(1:nall) > 0)
  allocate (tmpstructures(nallnew))
  allocate (tmpgroups(0:nallnew),source=0)
  allocate (double(nallnew),source=0)
  cc = 0
  do ii = 1,gcount
    do jj = 1,nall
      ggcount = groups(jj)
      if (ggcount .eq. ii.and.ggcount > 0) then
        cc = cc+1
        tmpstructures(cc) = structures(jj)
        tmpgroups(cc) = ggcount
        do kk = 1,cc-1
          if (tmpgroups(kk) .eq. ggcount) then
            double(cc) = kk
            exit
          end if
        end do
      end if
    end do
  end do
  tmpgroups(0) = gcount
  call move_alloc(tmpgroups,groups)
  call move_alloc(tmpstructures,structures)
  if (prlvl > 0) then
    write (prch,'(a)') ' done.'
    frac = real(nall-nallnew,wp)/real(nall,wp)
    write (prch,'(1x,a,t40,a,i10,a,f6.2,a)') &
    &      "number of doubles removed by MIC-RMSD",":",nall-nallnew,' (',frac*100.d0,'%)'
    write (prch,'(1x,a,t40,a,i10,a,f6.2,a)') &
    &      "number of unique structures remaining",":",nallnew,' (', (1.0d0-frac)*100.d0,'%)'
    frac = real(gcount,wp)/real(max(nallnew,1),wp)
    write (prch,'(1x,a,t40,a,i10,a,f6.2,a,i0,a)') &
    &      "number of unique conformers identified",":",gcount,' (',frac*100.d0,'% of ',nallnew,')'
  end if
  nall = nallnew

  !>-- for ENSO write a file with duplicate info (if required)
  call enso_duplicates(env,nall,double)

  if (allocated(prune_table)) deallocate (prune_table)
  if (allocated(sig)) deallocate (sig)
  if (allocated(vol)) deallocate (vol)
  if (allocated(mask)) deallocate (mask)
end subroutine cregen_CRE_periodic

!=========================================================================================!

subroutine cregen_irmsd_all(nall,structures,printlvl,iinversion)
!********************************************
!* Proof-of-concept routine to run all
!* pairs of RMSD for an array of structures
!********************************************
  use crest_parameters
  use crest_data
  use strucrd
  use axis_module
  use canonical_mod
  use irmsd_module
  use utilities,only:lin
  implicit none
  !> INPUT
  integer,intent(in) :: nall
  type(coord),intent(inout),target :: structures(nall)
  integer,intent(in),optional :: printlvl
  integer,intent(in),optional :: iinversion
  !> LOCAL
  integer :: i,j,ii,jj,T,nallpairs,cc,nat
  integer :: prlvl,iunit
  type(rmsd_cache),allocatable :: rcaches(:)
  type(coord),allocatable,target :: workmols(:)
  type(canonical_sorter),allocatable :: sorters(:)
  real(wp),allocatable :: rmsds(:)
  type(coord),pointer :: ref,mol
  type(coord) :: molloc
  real(wp) :: rmsdval,runtime
  logical :: stereocheck
  type(timer) :: profiler

  logical,parameter :: debug = .false.
  real(wp),allocatable :: debugrmsds(:)

  !> for implementing OpenMP parallelism
  T = 1

  !> print level
  if (present(printlvl)) then
    prlvl = printlvl
  else
    prlvl = 0
  end if

  !> set up timer
  call profiler%init(3)

  !> prepare workspace
  nallpairs = (nall*(nall+1))/2
  allocate (rmsds(nallpairs),source=0.0_wp)
  if (debug) then
    allocate (debugrmsds(nallpairs),source=0.0_wp)
  end if

  allocate (rcaches(T))
  ref => structures(1)
  nat = ref%nat
  allocate (workmols(T))
  do i = 1,T
    mol => workmols(i)
    allocate (mol%at(ref%nat))
    allocate (mol%xyz(3,ref%nat))
    nullify (mol)
    call rcaches(i)%allocate(ref%nat)
  end do

  !> set up ranks for each structure
  call profiler%start(1)
  allocate (sorters(nall))
  if (prlvl > 0) then
    write (stdout,'(a)',advance='no') 'CREGEN> Setting up canonical atom ranks ... '
    flush (stdout)
  end if
  do ii = 1,nall
    mol => structures(ii)
    call axis(mol%nat,mol%at,mol%xyz)
    call sorters(ii)%init(mol,invtype='apsp+',heavy=.false.)
    !call sorters(ii)%add_h_ranks(mol)
    if (ii == 1) then
      stereocheck = .not. (sorters(ii)%hasstereo(ref))
    end if
    call sorters(ii)%shrink()
  end do
  call profiler%stop(1)
  if (prlvl > 0) then
    call profiler%write_timing(stdout,1,'done.',.true.)
    runtime = (profiler%get(1)/real(nall,wp))*1000.0_wp
    write (stdout,'(a,f0.3,a)') 'CREGEN> Corresponding to approximately ',runtime, &
    &                       ' ms per processed structure'
  end if

  !> allow user to set inversion check (false rotamers)
  if (present(iinversion)) then
    select case (iinversion)
    case (0)
      continue
    case (1)
      stereocheck = .true.
    case (2)
      stereocheck = .false.
    end select
    if (prlvl > 1) then
      write (stdout,'(a,l2)') 'CREGEN> Check for false rotamers (geometry inversion)? -->',stereocheck
    end if
  end if

  !> And finally, run the RMSD checks
  call profiler%start(2)
  if (prlvl > 0) then
    write (stdout,*)
    write (stdout,'(a)',advance='no') 'CREGEN> Running all pair RMSDs ... '
    flush (stdout)
  end if
  cc = 1
  do ii = 1,nall
    rcaches(cc)%stereocheck = stereocheck
    rcaches(cc)%rank(1:nat,1) = sorters(ii)%rank(1:nat)
    do jj = ii+1,nall
      workmols(cc)%nat = structures(jj)%nat
      workmols(cc)%at(:) = structures(jj)%at(:)
      workmols(cc)%xyz(:,:) = structures(jj)%xyz(:,:)
      !molloc = structures(jj)
      rcaches(cc)%rank(:,2) = sorters(jj)%rank(:)
      call min_rmsd(structures(ii),workmols(cc), &
      &        rcache=rcaches(cc),rmsdout=rmsdval)
      rmsds(lin(ii,jj)) = rmsdval
    end do
  end do
  call profiler%stop(2)
  if (prlvl > 0) then
    call profiler%write_timing(stdout,2,'done.',.true.)
    !write (stdout,'(a)',advance='yes') 'done.'
    runtime = (profiler%get(2)/real(nallpairs,wp))*1000.0_wp
    write (stdout,'(a,f0.3,a)') 'CREGEN> Corresponding to approximately ',runtime, &
    &                       ' ms per processed RMSD'

  end if

  if (debug) then
    !> RMSD without permutation
    do ii = 1,nall
      do jj = ii+1,nall
        rmsdval = rmsd(structures(ii),structures(jj))
        debugrmsds(lin(ii,jj)) = rmsdval
      end do
    end do
  end if

  if (prlvl > 1) then
    write (stdout,'(a)') 'CREGEN> Writing cregen_rmsds.csv with RMSDs in Angström'
    open (newunit=iunit,file='cregen_rmsds.csv')
    if (debug) then
      write (iunit,'(a,3(",",a))') 'A','B','rmsd','rmsdref'
      do ii = 1,nall
        do jj = ii+1,nall
          write (iunit,'(i0,",",i0,2(",",f0.7))') &
          & min(ii,jj),max(ii,jj),rmsds(lin(ii,jj))*autoaa,debugrmsds(lin(ii,jj))*autoaa
        end do
      end do
    else
      write (iunit,'(a,",",a,",",a)') 'A','B','rmsd'
      do ii = 1,nall
        do jj = ii+1,nall
          write (iunit,'(i0,",",i0,",",f0.7)') min(ii,jj),max(ii,jj),rmsds(lin(ii,jj))*autoaa
        end do
      end do
    end if
    close (iunit)
  end if

  deallocate (sorters)
  deallocate (workmols)
  deallocate (rcaches)
  deallocate (rmsds)
end subroutine cregen_irmsd_all

!=========================================================================================!

subroutine cregen_irmsd_sort(env,nall,structures,groups,allcanon,printlvl)
!*******************************************************
!* Proof-of-concept routine to analyze an
!* ensemble only via the iRMSD procedure.
!* Conformers are identified by the rthr threshold only
!*******************************************************
  use crest_parameters
  use crest_data
  use iomod,only:to_str
  use strucrd
  use axis_module
  use canonical_mod
  use irmsd_module
  use utilities,only:lin
  use quicksort_interface
  use omp_lib
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env
  integer,intent(in) :: nall
  type(coord),intent(inout),target :: structures(nall)
  integer,intent(inout) :: groups(nall)
  logical,intent(in),optional :: allcanon
  integer,intent(in),optional :: printlvl

  !> LOCAL
  integer :: i,j,ii,jj,T,Tn,nallpairs,cc,nat,k
  integer :: gcount
  integer :: prlvl,iunit
  type(rmsd_cache),allocatable :: rcaches(:)
  type(coord),allocatable,target :: workmols(:)
  type(canonical_sorter),allocatable :: sorters(:)
  real(wp),allocatable :: rmsds(:)
  type(coord),pointer :: ref,mol
  type(coord) :: molloc
  real(wp) :: rmsdval,runtime,RTHR,ETHR,ediff
  logical :: stereocheck,individual_IDs
  type(timer) :: profiler
  integer :: ng
  integer,allocatable :: group(:),degen(:,:)
  real(wp),allocatable :: er(:)
  type(coord),allocatable :: structures_new(:)

  logical,parameter :: debug = .false.

!>--- handle optional arguments
  if (present(allcanon)) then
    individual_IDs = .not.allcanon
  else
    individual_IDs = .false.
  end if
  if (present(printlvl)) then
    prlvl = printlvl
  else
    prlvl = 1
  end if

!>--- set up parallelization
  call new_ompautoset(env,'max',nall,T,Tn)

!>--- set up timer
  call profiler%init(3)

!>--- set up parameters (note we are working with BOHR internally)
  RTHR = env%rthr*aatoau
  ETHR = env%ethr/autokcal

!>--- print some sorting data
  if (prlvl > 0) then
    write (stdout,'(a)') 'CREGEN> Info for iRMSD sorting:'
    write (stdout,'(2x,a,t32,a,i10)') 'number of structures',':',nall
    write (stdout,'(2x,a,t32,a,f10.5,a)') 'RTHR (RMSD threshold)',':',RTHR*autoaa,' Å'
    write (stdout,'(2x,a,t32,a,es10.2,a)') 'ETHR (energy threshold)',':',ETHR,' Ha'
    write (stdout,'(2x,a,t32,a,i10)') 'OpenMP threads',':',T
    write (stdout,'(2x,a,t32,a,a10)') 'Individual atom IDs?',':',to_str(individual_IDs)
    write (stdout,'(2x,a,t32,a)',advance='no') 'False rotamer check?',':'
    select case (env%iinversion)
    case (0)
      write (stdout,'(a10)') 'auto'
    case (1)
      write (stdout,'(a10)') 'on'
    case (2)
      write (stdout,'(a10)') 'off'
    end select
    write (stdout,*)
  end if

!>--- sorting by energy
  write (stdout,'(a)',advance='no') 'CREGEN> Sorting ensemble by energy ... '
  flush (stdout)
  call ensemble_qsort(nall,structures,1,nall)
  write (stdout,'(a)') 'done.'

!>--- Set up atom identities (either for all, or just the first structure)
  if (individual_IDs) then
    allocate (sorters(nall))
  else
    allocate (sorters(1))
  end if
  if (prlvl > 0) then
    write (stdout,'(a)',advance='no') 'CREGEN> Setting up canonical atom ranks ... '
    flush (stdout)
    call profiler%start(1)
  end if
  ref => structures(1)
  !$omp parallel &
  !$omp shared(sorters, structures, stereocheck) &
  !$omp private(mol,ii)
  !$omp do schedule(dynamic)
  do ii = 1,nall
    mol => structures(ii)
    call axis(mol%nat,mol%at,mol%xyz)
    if (individual_IDs.or.ii == 1) then
      call sorters(ii)%init(mol,invtype='apsp+',heavy=.false.)
    end if
    if (ii == 1) then
      stereocheck = .not. (sorters(ii)%hasstereo(ref))
    end if
    if (individual_IDs.or.ii == 1) then
      call sorters(ii)%shrink()
    end if
  end do
  !$omp end do
  !$omp end parallel
  if (prlvl > 0) then
    call profiler%stop(1)
    call profiler%write_timing(stdout,1,'done.',.true.)
    runtime = (profiler%get(1)/real(nall,wp))*1000.0_wp
    write (stdout,'(1x,a,f0.3,a)') '* Corresponding to approximately ',runtime, &
    &                       ' ms per processed RMSD'
    write (stdout,*)
  end if

  !>--- allow user to set inversion check (false rotamers)
  select case (env%iinversion)
  case (0)
    continue
  case (1)
    stereocheck = .true.
  case (2)
    stereocheck = .false.
  end select
  if (prlvl > 1) then
    write (stdout,'(a,l2)') 'CREGEN> Check for false rotamers (geometry inversion)? -->',stereocheck
  end if

!>--- allocate work cache
  if (prlvl > 0) then
    write (stdout,'(a)',advance='no') 'CREGEN> Allocating iRMSD work cache ... '
    flush (stdout)
  end if
  allocate (rcaches(T))
  ref => structures(1)
  nat = ref%nat
  allocate (workmols(T))
  do i = 1,T
    mol => workmols(i)
    allocate (mol%at(ref%nat))
    allocate (mol%xyz(3,ref%nat))
    nullify (mol)
    call rcaches(i)%allocate(ref%nat)
    rcaches(i)%stereocheck = stereocheck
  end do
  if (prlvl > 0) then
    write (stdout,'(a)') 'done.'
    write (stdout,*)
  end if

!>--- run the checks
  if (prlvl > 0) then
    write (stdout,'(a)',advance='no') 'CREGEN> Running all-pair iRMSDs ... '
    flush (stdout)
    call profiler%start(2)
  end if
  gcount = maxval(groups(:))
  do ii = 1,nall
!>--- find next unassigned conformer and assign a new group
    if (groups(ii) .ne. 0) cycle
    gcount = gcount+1
    groups(ii) = gcount

!>--- Then, cross-check all other unassigned conformers
    !$omp parallel &
    !$omp shared(nall, nat, groups, individual_IDs, sorters, rcaches) &
    !$omp shared(workmols, structures, ii, ETHR) &
    !$omp private(jj,rmsdval,cc,ediff)
    !$omp do schedule(dynamic)
    do jj = ii+1,nall
      cc = omp_get_thread_num()+1
      if (groups(jj) .ne. 0) cycle
      ediff = abs(structures(ii)%energy-structures(jj)%energy)
      if (ediff > ETHR) cycle
      if (individual_IDs) then
        rcaches(cc)%rank(1:nat,1) = sorters(ii)%rank(1:nat)
        rcaches(cc)%rank(1:nat,2) = sorters(jj)%rank(1:nat)
      else
        rcaches(cc)%rank(1:nat,1) = sorters(1)%rank(1:nat)
        rcaches(cc)%rank(1:nat,2) = sorters(1)%rank(1:nat)
      end if
      workmols(cc)%nat = structures(jj)%nat
      workmols(cc)%at(:) = structures(jj)%at(:)
      workmols(cc)%xyz(:,:) = structures(jj)%xyz(:,:)
      call min_rmsd(structures(ii),workmols(cc), &
      &        rcache=rcaches(cc),rmsdout=rmsdval)
      if (rmsdval < RTHR) groups(jj) = gcount
    end do
    !$omp end do
    !$omp end parallel
  end do
  if (prlvl > 0) then
    call profiler%stop(2)
    call profiler%write_timing(stdout,2,'done.',.true.)
    write (stdout,*)
  end if

  if (debug) then
    write (*,*) 'assigned groups, and count'
    do ii = 1,maxval(groups(:))
      write (*,*) ii,count(groups(:) == ii)
    end do
  end if

  allocate (group(0:nall),source=0)
  allocate (structures_new(nall))
  ng = maxval(groups(:))
  k = 0
  do ii = 1,ng
    do jj = 1,nall
      if (groups(jj) == ii) then
        k = k+1
        group(k) = ii
        structures_new(k) = structures(jj)
      end if
    end do
  end do
  group(0) = maxval(groups(:))
  allocate (degen(3,ng))
  call cregen_groupinfo(nall,ng,group,degen)
  allocate (er(nall))
  do ii = 1,nall
    er(ii) = structures_new(ii)%energy
    structures(ii) = structures_new(ii)
  end do
  if (prlvl > 0) then
    call cregen_pr2(stdout,env,nall,ng,degen,er)
    call cregen_econf_list(stdout,nall,er,ng,degen)
  end if
  if (prlvl > 1) then
    write (stdout,'(a,a)') 'Unique-structure file written to: ',ensemblefile
    block
      use cregen_subroutines,only:cregen_conffile
      call cregen_conffile(env,ensemblefile,structures,ng,degen)
    end block
  end if

end subroutine cregen_irmsd_sort

!=========================================================================================!

subroutine cregen_EQUAL(ch,nall,structures,group,athr,rotfil)
!****************************************************************
!* subroutine cregen_EQUAL
!* subroutine for the generation of nuclear equivalencies
!* On Input: ch - printout channel
!*           nall - number of structures
!*           structures  - the list of structures
!*           group - to which group does every strucutre belong
!*           athr  - threshold for equivalency comparison
!*           rotfil - wirte anmr_rotamer file?
!* On Output: resorted xyz and comments
!****************************************************************
  use crest_parameters,id => dp
  use crest_data
  use strucrd
  use miscdata,only:rcov
  use crest_cn_module
  use utilities
  implicit none
  integer,intent(in) :: ch
  integer,intent(in) :: nall
  type(coord),intent(in) :: structures(nall)
  integer,intent(inout) :: group(0:nall)
  real(wp),intent(in) :: athr
  logical,intent(in) :: rotfil
  integer,allocatable :: at(:)
  real(wp),allocatable :: cdum(:,:)
  integer :: ng,n,nat
  integer :: i,j,k,l
  logical :: ex

  !>--- arrays and variable for the analysis
  integer :: gmax
  integer,allocatable :: glist(:,:)
  integer :: current
  real(wp) :: dum
  real(wp) :: shortest_distance
  integer,allocatable :: equiv(:,:,:)
  integer,allocatable :: pair(:),pre(:),nb(:,:)
  logical,allocatable :: vis(:)
  real(wp),allocatable :: metric(:,:)
  real(wp),allocatable :: dist(:,:,:)
  integer,allocatable :: relat(:,:)
  real(wp),allocatable :: tmp2(:)
  integer :: m,m1,m2,s1,s2,iat,j1,k2

  !>-- further NMR-mode related data
  integer,allocatable :: nmract(:)
  integer,allocatable :: elist(:,:),flist(:,:)
  integer,allocatable :: jnd(:)
  real(wp),allocatable :: sd(:,:),jfake(:),cn(:)

  character(len=:),allocatable :: atmp
  integer :: ig,ir,irr,nr

!>--- infer from structure list
  nat = structures(1)%nat
  allocate (at(nat))
  at(:) = structures(1)%at(:)

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
  !cdum(1:3,1:n) = xyz(1:3,1:n,1) / bohr
  cdum(1:3,1:n) = structures(1)%xyz(1:3,1:n)
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
    !call distance(n,xyz(:,:,i),dist(:,:,i)) !> distance matrix
    call distance(n,structures(i)%xyz(:,:),dist(:,:,i)) !> distance matrix
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
        call distance(n,structures(irr)%xyz(:,:),sd)   !> distance matrix
        cdum(1:3,1:n) = structures(irr)%xyz(1:3,1:n)
        call calculate_CN(n,at,cdum,cn)
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
!=========================================================================================!

subroutine cregen_file_wr(env,fname,structures)
!*********************************************************************
!* write the output ensemble file with all structures (rotamer file)
!*********************************************************************
  use crest_parameters
  use crest_data
  use strucrd
  use utilities,only:boltz
  implicit none
  type(systemdata),intent(inout) :: env
  character(len=*),intent(in) :: fname
  type(coord),intent(inout) :: structures(:)
  integer :: nat,nall
  character(len=128) :: newcomment

  integer :: ich,ii
  real(wp) :: eref,T
  real(wp),allocatable :: er(:),erel(:),p(:)
  character(len=40),allocatable :: origin(:)

  nall = size(structures,1)
  allocate (er(nall),erel(nall),p(nall))!,origin(nall))
  eref = structures(1)%energy
  do ii = 1,nall
    er(ii) = structures(ii)%energy
    erel(ii) = (er(ii)-eref)*autokcal
    !if (env%trackorigin) then
    !  call getorigin(comments(i),origin(i))
    !end if
  end do
  T = env%tboltz
  call boltz(nall,T,erel,p)

  open (newunit=ich,file=fname)
  do ii = 1,nall
    !if (env%trackorigin) then
    !write (newcomment,'(a,f10.8,1x,a)') 'population=',p(ii),trim(origin(ii))
    !else
    write (newcomment,'(a,f10.8)') 'population=',p(ii)
    !end if
    structures(ii)%comment = trim(newcomment)
    call structures(ii)%append(ich)
  end do
  close (ich)
  !deallocate (origin,p,erel,er)
  if (allocated(origin)) deallocate (origin)
  deallocate (p,erel,er)
  return
end subroutine cregen_file_wr

!=========================================================================================!

subroutine cregen_conffile(env,cname,structures,ng,degen)
!*********************************
!* write the output ensemble file
!*********************************
  use crest_parameters,only:wp,bohr
  use crest_data
  use strucrd
  use iomod
  use utilities
  implicit none
  type(systemdata),intent(inout) :: env
  character(len=*),intent(in) :: cname
  type(coord),intent(inout) :: structures(:)
  integer,intent(in) :: ng
  integer,intent(in) :: degen(3,ng)
  integer :: nat,nall
  integer :: ich,ich3,ichenso
  integer :: i,k,ii
  real(wp),allocatable :: er(:)

  nall = size(structures,1)
  allocate (er(nall))
  do ii = 1,nall
    er(ii) = structures(ii)%energy
    if (allocated(structures(ii)%comment)) &
    &  deallocate (structures(ii)%comment)
  end do
  if (env%enso) then
    open (newunit=ichenso,file='enso.tags')
  end if
  call structures(1)%write('crest_best.xyz')
  open (newunit=ich,file=trim(cname))
  do ii = 1,ng
    k = degen(2,ii)
    if (k <= 0.or.k > nall) cycle
    call structures(k)%append(ich)
    if (env%enso) write (ichenso,'(2x,f18.8)') er(k)
  end do
  close (ich)
  if (env%enso) then
    close (ichenso)
  end if
  deallocate (er)

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

subroutine cregen_rmsdalign(nall,structures)
!*****************************************************
!* Algin all structures in an array to the first one
!* in the ensemble, based on the heavy-atom RMSD
!*****************************************************
  use crest_parameters
  use irmsd_module
  use strucrd
  implicit none
  integer,intent(in) :: nall
  type(coord),intent(inout) :: structures(nall)
  integer :: ii,nat
  logical,allocatable :: mask(:)

  nat = structures(1)%nat
  allocate (mask(nat),source=.false.)
  do ii = 1,nat
    if (structures(1)%at(ii) > 1) mask(ii) = .true.
  end do

  do ii = 2,nall
    call rmsd_align(structures(1),structures(ii),mask=mask)
  end do

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
  use omp_lib
  implicit none
  type(systemdata) :: env
  integer :: ch
  logical :: pr
  !integer :: TID,OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM,
  integer :: TID,nproc,T,Tn
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
  write (ch,'(80("*"))')
  write (ch,'(" number of atoms",t35,":",i10)') nat
  if (substruc) then
    write (ch,'(" atoms included in RMSD",t35,":",i10)') env%rednat
  end if
  write (ch,'(" number of points on xyz file",t35,":",i10)') nall
  !write (ch,'('' RMSD threshold                 :'',f9.4)') rthr
  !write (ch,'('' Bconst threshold               :'',f9.4)') bthr
  !write (ch,'('' population threshold           :'',f9.4)') pthr
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

subroutine cregen_pr2(ch,env,nall,ng,degen,er)
  use crest_parameters
  use crest_data
  use strucrd
  use iomod,only:touch,remove
  use utilities,only:boltz
  implicit none
  integer,intent(in) :: ch
  type(systemdata),intent(inout) :: env
  integer,intent(in) :: nall
  integer,intent(in) :: ng
  integer,intent(in) :: degen(3,ng)
  real(wp),intent(in) :: er(nall)
  integer :: ich,chref,och,och2
  integer :: i,j,k
  real(wp),allocatable :: erel(:),egrp(:)
  real(wp),allocatable :: p(:),pg(:),paccu(:)
  real(wp) :: eref,T
  character(len=40),allocatable :: origin(:)
  integer :: a,b
  logical :: ex,abbrev,print_placeholder
  real(wp) :: A0,eav,g,s,ss,beta,elow
  integer,parameter :: printlimit = 100

  allocate (origin(nall),erel(nall),p(nall))
  eref = minval(er,1)
  env%elowest = eref
  do i = 1,nall
    !er(i) = grepenergy(comments(i))
    erel(i) = (er(i)-eref)*autokcal
    !if (env%trackorigin) then
    !  call getorigin(comments(i),origin(i))
    !else
    origin(i) = ''
    !end if
  end do
  T = env%tboltz
  call boltz(nall,T,erel,p)
  allocate (pg(ng),source=0.0_wp)
  allocate (paccu(0:nall),source=0.0_wp)
  do i = 1,ng
    a = degen(2,i)
    b = degen(3,i)
    do j = a,b
      pg(i) = pg(i)+p(j)
      paccu(j) = paccu(j-1)+p(j)
    end do
  end do

  och = ch
  abbrev = nall > printlimit

  !>-- really long energy list
  write (och,'(80("*"))')
  write (och,'(1x,a8,1x,a8,(1x,a16),(1x,a8),(1x,a12),1x,a9,1x,a5)') &
  &      '  ','ΔE','Etot','weight','conf.weight','conformer',''
  write (och,'(a8,1x,a8,(1x,a16),(1x,a8),(1x,a12),1x,a9,1x,a5,1x,a6)') &
  &        'id ','kcal/mol','hartree','p(i)','p(group)','group','degen','origin'
  write (och,'(4x,4("-"),1x,8("-"),(1x,16("-")),(1x,8("-")),(1x,12("-")),1x,9("-"),1x,5("-"),1x,6("-"))')
  if (abbrev) then
    call remove('cregen.full')
    open (newunit=och2,file='cregen.full',status='replace')
    write (och2,'(1x,a8,1x,a8,(1x,a16),(1x,a8),(1x,a12),1x,a9,1x,a5)') &
    &      '  ','ΔE','Etot','weight','conf.weight','conformer',''
    write (och2,'(a8,1x,a8,(1x,a16),(1x,a8),(1x,a12),1x,a9,1x,a5,1x,a6)') &
    &        'id ','kcal/mol','hartree','p(i)','p(group)','group','degen','origin'
    write (och2,'(4x,4("-"),1x,8("-"),(1x,16("-")),(1x,8("-")),(1x,12("-")),1x,9("-"),1x,5("-"),1x,6("-"))')
  else
    call remove('cregen.full')
  end if

  print_placeholder = .true.
  k = 0
  do i = 1,ng
    k = k+1
    a = degen(2,i)
    b = degen(3,i)
    if (k <= printlimit.or.k > nall-10) then
      write (och,'(i8,1x,f8.4,1x,f16.6,(1x,f8.5),(1x,f12.5),1x,i9,1x,i5,a)') &
      &     a,erel(a),er(a),p(a),pg(i),i,degen(1,i),trim(origin(a))
    else if (print_placeholder) then
      print_placeholder = .false.
      write (och,'(5x,"...",1x,"<skipped due to ensemble length, see file cregen.full> ...")')
    end if
    if (abbrev) then
      write (och2,'(i8,1x,f8.4,1x,f16.6,(1x,f8.5),(1x,f12.5),1x,i9,1x,i5,a)') &
      &     a,erel(a),er(a),p(a),pg(i),i,degen(1,i),trim(origin(a))
    end if
    do j = a+1,b
      k = k+1
      if (k <= printlimit.or.k > nall-10) then
        write (och,'(i8,1x,f8.4,1x,f16.6,1x,f8.5,1x,a12,1x,a9,1x,a5,1x,a)') &
               & k,erel(j),er(j),p(j),'.','.','.',trim(origin(j))
      else if (print_placeholder) then
        print_placeholder = .false.
        write (och,'(5x,"...",1x,"<skipped due to ensemble length, see file cregen.full> ...")')
      end if
      if (abbrev) then
        write (och2,'(i8,1x,f8.4,1x,f16.6,1x,f8.5,1x,a12,1x,a9,1x,a5,1x,a)') &
               & k,erel(j),er(j),p(j),'.','.','.',trim(origin(j))
      end if
    end do
  end do

  if (abbrev) close (och2)

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

  write (och,'(80("*"))')
  write (och,'("Statistics for *THIS* ensemble:")')
  write (och,'(35("-"))')
  write (och,'(" Number of groups & total",t42,":",2x, i0,", ",i0)') ng,nall
  write (och,'(" Temperature used for populations",t42,":",2x,F9.2," K")') T
  write (och,'(" Energy of lowest structure",t42,":",2x,es14.6)') eref
  !>---- elow printout in between routines
  if (.not.env%confgo) then
    write (stdout,'("CREGEN> E lowest :",f20.10,a)') eref,' Ha'
  end if
  write (och,'(" Ensemble average energy (kcal/mol)",t42,":",2x,F14.8)') eav
  if (env%QCG) then
    write (och,'(" Ensemble entropy (cal/mol K)",t42,":",2x,F14.8)') ss
  else
    write (och,'(" Ensemble entropy (J/mol K, cal/mol K)",t42,":",2x,2F9.3)') s,ss
  end if
  write (och,'(" Ensemble free energy (kcal/mol)",t42,":",2x,F14.8)') g
  write (och,'(" Population of lowest strucure",t42,":",2x,F9.3," %")') pg(1)*100.d0
  write (och,'(" Highest population & group",t42,":",2x,F9.3," %, ",i0)') maxval(pg,1)*100.d0,maxloc(pg,1)

  j = min(10,ng)
  i = degen(3,j)
  write (och,'(" Accum.population of lowest 10 groups",t42,":",2x,F9.3," %")') paccu(i)*100.d0
  do i = 1,ng
    j = degen(3,i)
    if (paccu(j) >= 0.5_wp) exit
  end do
  write (och,'(" 50% accum.population for groups",t42,":",6x,"1 - ",i0)') i
  do i = 1,ng
    j = degen(3,i)
    if (paccu(j) >= 0.95_wp) exit
  end do
  write (och,'(" 95% accum.population for groups",t42,":",6x,"1 - ",i0)') i

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

  write (och,'(80("*"))')

  deallocate (paccu,pg)
  deallocate (p,erel,origin)
  return
end subroutine cregen_pr2

subroutine cregen_econf_list(ch,nall,er,ng,degen)
  use crest_parameters
  implicit none
  integer,intent(in) :: ch
  integer,intent(in) :: nall
  real(wp),intent(in) :: er(nall)
  integer,intent(in) :: ng
  integer,intent(in) :: degen(3,ng)
  integer :: ich2,i,j
  real(wp) :: eref,ewrt

  write (ch,'(a,i0)') 'Number of unique conformers for further calculation: ',ng
  write (ch,'(a)') 'List of relative energies (kcal/mol) saved as "crest.energies"'
  open (newunit=ich2,file='crest.energies')
  eref = minval(er,1)
  do i = 1,ng
    j = degen(2,i)
    ewrt = er(j)-eref
    ewrt = ewrt*autokcal
    write (ich2,'(i10,1x,f12.4,es20.10)') i,ewrt,er(i)
  end do
  close (ich2)

  return
end subroutine cregen_econf_list

subroutine cregen_pr3(ch,infile,nall,er)
  use crest_parameters
  use strucrd
  implicit none
  integer,intent(in) :: ch
  character(len=*),intent(in) :: infile
  integer,intent(in) :: nall
  real(wp),intent(in) :: er(nall)
  real(wp) :: dE,eref
  integer :: i
  write (ch,*)
  write (ch,'(a)') '====================================================='
  write (ch,'(a)') '============== ordered structure list ==============='
  write (ch,'(a)') '====================================================='
  write (ch,'(a,a,a)') ' written to file <',trim(infile),'>'
  write (ch,*)
  write (ch,'(a10,4x,a15,a25)') 'structure','ΔE(kcal/mol)','Etot(Eh)'
  eref = minval(er,1)
  !write (ch,'(''   structure    ΔE(kcal/mol)    Etot(Eh)'')')
  do i = 1,nall
    dE = (er(i)-eref)*autokcal
    write (ch,'(i10,3x,F15.4,F25.10)') i,dE,er(i)
  end do
  write (ch,*)
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
