!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020-2024 Philipp Pracht, Johannes Gorges
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

module msmod
  use crest_parameters
  use crest_data
  use strucrd
  use iomod
  implicit none

  !>-- storage for a single mol
  type msmol
    integer :: nat
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    real(wp) :: etot
    integer :: chrg
    integer :: nfrag
    !> index of fragatoms 1 first fragment 2 second fragment
    integer,allocatable :: fragi(:) 
    !> number of fragments in product 
    !> 1: isomer 1< : fragmentation; 0, not determined yet, -1 means we sort this out
    integer :: fragcount
    real(wp) :: molmass
  contains
    procedure :: newmol => msmod_newmol
    procedure :: dist => msmoldist
  end type msmol

  public :: msmol

  !-- product list
  type plist
    integer :: nmol
    type(msmol),allocatable :: mol(:)
    logical,allocatable :: new(:)
  contains
    procedure :: dealloc => deallocate_plist
    procedure :: append => plist_append
    procedure :: remove => plist_remove
  end type plist

  !>-- global data for msreact
  type msobj
    real(wp) :: rcut = 1.3_wp        !> cutoff for factor of sum of covalent radii for detection of fragmentation of fragment_structures
    real(wp) :: T = 5000.0_wp        !> electronic temperature! favours open-shell fragment generation
    real(wp) :: fc = 0.05_wp         !> repulsive force constant in constrained optimizations
    real(wp) :: fc_attr = -0.05_wp   !> attractive force constant in constrained optimizations
    real(wp) :: cdist = 1.5_wp       !> constraing distance scaling factor of rcov
    real(wp) :: cdist_att = 0.5_wp   !> constraing distance scaling factor of rcov for attractive part
    real(wp) :: distthr_att = 4.0_wp !> distance threshold in Angstrom for attractive part of constraint
    real(wp) :: atomshift = 0.75_wp  !> shift of atoms in random atom displacement
    real(wp) :: fragdist = 0.0_wp    !> increase distance by to separate fragments (useful for subsequent transition state search)
    integer  :: maxc = 15            !>max optimization cycle in constrained optimizations

    type(plist) :: pl
  end type msobj

  public :: msobj

!=======================================================================================!
!=======================================================================================!
contains  !> MODULE PROCEDURES START HERE
!=======================================================================================!
!=======================================================================================!

  subroutine deallocate_plist(self)
    class(plist) :: self
    if (allocated(self%mol)) deallocate (self%mol)
    if (allocated(self%new)) deallocate (self%new)
    return
  end subroutine deallocate_plist

!=======================================================================================!

  subroutine msmod_newmol(self,nat,xyz,at,etot,chrg,fragcount,molmass)
    implicit none
    class(msmol) :: self
    integer :: nat
    real(wp) :: xyz(3,nat)
    integer :: at(nat)
    real(wp) :: etot
    integer :: chrg
    integer :: fragcount
    real(wp) :: molmass
    self%nat = nat
    allocate (self%xyz(3,nat))
    self%xyz = xyz
    allocate (self%at(nat))
    self%at = at
    self%etot = etot
    self%chrg = chrg
    self%fragcount = fragcount
    self%molmass = molmass
    return
  end subroutine msmod_newmol

!=======================================================================================!

  function msmoldist(self,i,j) result(dist)
!***************************************************
!* calculate the distance between two atoms i and j
!***************************************************
    implicit none
    class(msmol) :: self
    integer :: i,j
    real(wp) :: dist
    dist = 0.0_wp
    if (allocated(self%xyz)) then
      dist = (self%xyz(1,i)-self%xyz(1,j))**2+ &
     &     (self%xyz(2,i)-self%xyz(2,j))**2+ &
     &     (self%xyz(3,i)-self%xyz(3,j))**2
      dist = sqrt(dist)
    end if
    return
  end function msmoldist

!=======================================================================================!

  subroutine plist_append(self,nat,at,xyz,etot,chrg,fragcount,molmass)
    implicit none
    class(plist) :: self
    type(msmol) :: mol
    integer :: nat
    real(wp) :: xyz(3,nat)
    integer :: at(nat)
    real(wp) :: etot
    integer :: chrg
    integer :: fragcount
    real(wp) :: molmass

    type(msmol),allocatable :: dummy(:)
    logical,allocatable     :: btmp(:)
    integer :: n
    call mol%newmol(nat,xyz,at,etot,chrg,fragcount,molmass)
    if (.not.allocated(self%mol)) then
      self%nmol = 1
      allocate (self%mol(1))
      allocate (self%new(1))
      self%mol(1) = mol
      self%new(1) = .true.
    else
      n = self%nmol+1
      allocate (dummy(n))
      allocate (btmp(n))
      dummy(1:self%nmol) = self%mol(1:self%nmol)
      dummy(n) = mol
      btmp(1:self%nmol) = self%new(1:self%nmol)
      btmp(n) = .true.
      self%nmol = n
      call move_alloc(btmp,self%new)
      call move_alloc(dummy,self%mol)
    end if
    return
  end subroutine plist_append

!=======================================================================================!

  subroutine plist_remove(self,index)
!****************************************************
!* remove entry from product with index i from  list
!****************************************************
    implicit none
    class(plist) :: self
    type(msmol) :: mol

    integer :: n,i
    integer,intent(in) :: index

    if (index < 1.or.index > self%nmol) then
      write (*,*) "Invalid index. Cannot remove molecule."
      return
    end if

    do i = index,self%nmol-1
      self%mol(i) = self%mol(i+1) ! replace the entry with the next one
    end do

    self%nmol = self%nmol-1 ! ignore last entry

  end subroutine plist_remove

!=======================================================================================!

  subroutine get_input_energy(env,etemp,etot)
!******************************************************************************
!* A quick xtb geometry optimization in xyz coordinates to get starting energy
!******************************************************************************
    implicit none
    type(systemdata) :: env
    character(len=80) :: fname,pipe
    character(len=:),allocatable :: jobcall
    logical :: fin
    character(len=256) :: atmp
    integer :: ich,iost,io,i
    type(coord) :: mol
    integer :: ntopo
    integer,allocatable :: topo(:)
    real(wp),intent(out) :: etot
    real(wp) :: etemp ! electronic temperature in K
    logical :: tchange = .false.
    logical :: ldum

!---- setting threads
    if (env%autothreads) then
      call ompautoset(env%threads,7,env%omp,env%MAXRUN,1) !set the global OMP/MKL variables for the xtb jobs
    end if

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

    fname = '.CHRG'
    open (newunit=ich,file=fname)
    write (ich,'(i0)') env%chrg  ! EI +1, DEA -1, CID 0, give in crest call
    close (ich)

!---- input xyz file
    fname = env%inputcoords

!    write (jobcall,'(a,1x,a,f10.4,1x,a,1x,a)') &
!    &     trim(env%ProgName),trim(fname)//" --sp --etemp ",etemp,trim(env%gfnver),trim(pipe)
    jobcall = trim(env%ProgName)
    jobcall = trim(jobcall)//' '//trim(fname)
    write(atmp,'(f10.4)') etemp
    jobcall = trim(jobcall)//' --sp --etemp '//trim(atmp)
    jobcall = trim(jobcall)//trim(env%gfnver)//trim(pipe)

    call execute_command_line(trim(jobcall),exitstat=io)

    call minigrep('xtb.out','finished run',fin)
    if (.not.fin) then
      write (*,*)
      write (*,*) ' Initial singlepoint calculation failed!'
      write (*,*) ' Please check your input.'
      error stop
    end if
    write (*,*) 'Input energy successfully calculated.'

    call grepval('xtb.out',"| TOTAL ENERGY",ldum,etot)
!---- cleanup
    call remove('xtb.out')
    call remove('energy')
    if (.not.env%chargesfile) call remove('charges')
    call remove('grad')
    call remove('mos')
    call remove('xtbopt.log')
    call remove('xtbtopo.mol')
    call remove('.xtbopttok')
    call remove('xtbrestart')
    call remove('gfnff_topo')
  end subroutine get_input_energy

!=======================================================================================!

  subroutine readbasicpos(env,nbaseat,basicatlist)
!********************************************************************
!* this subroutine reads the lmo.out file of an xtb -lmo calculation
!* and identifies the atoms which have a LP or pi-Orbital
!********************************************************************
    implicit none
    type(systemdata) :: env
    integer :: i,io,ich,at1,at2,j
    integer,intent(out) ::  nbaseat
    integer,intent(out),allocatable :: basicatlist(:)
    integer,allocatable :: dumlist(:)
    character(len=512) :: tmp
    character(len=64) :: fname
    character(len=64) :: type,dumc
    real(wp) :: dumr
    logical :: ex
    fname = 'xtb.out'
    inquire (file=fname,exist=ex)
    if (.not.ex) then
      write (*,*) "xtb.out not found"
      stop
    end if
    nbaseat = 0
    open (newunit=ich,file=fname)
    do
      read (ich,'(a)',iostat=io) tmp
      if (index(tmp,'files:') .ne. 0) exit
      if (index(tmp,'pi ') .ne. 0) nbaseat = nbaseat+2 ! count first two  atoms of pi or delpi bond with highest participation and ignore the rest
      if (index(tmp,'LP ') .ne. 0) nbaseat = nbaseat+1 ! count LP
    end do
    allocate (dumlist(nbaseat))
    dumlist = 0
    j = 0
    rewind (ich)
    do
      read (ich,'(a)',iostat=io) tmp
      if (index(tmp,'starting deloc pi') .ne. 0) exit
      if (index(tmp,'files:') .ne. 0) exit
      if (index(tmp,'pi ') .ne. 0) then
        backspace (ich)
        if (tmp(64:64) == ' ') then ! if element symbol has one character, there is a space and we use this routine
          if (tmp(80:80) == ' ') then
            read (ich,*) dumr,type,dumr,dumr,dumr,dumr,dumr,at1,dumc,dumc,dumr,at2,dumc,dumc,dumr ! first element has one and second has one character
          else
            read (ich,*) dumr,type,dumr,dumr,dumr,dumr,dumr,at1,dumc,dumc,dumr,at2,dumc,dumr ! first element has one and second has two character
          end if
        else
          if (tmp(80:80) == ' ') then
            read (ich,*) dumr,type,dumr,dumr,dumr,dumr,dumr,at1,dumc,dumr,at2,dumc,dumc,dumr  ! first element has two and second has pm character
          else
            read (ich,*) dumr,type,dumr,dumr,dumr,dumr,dumr,at1,dumc,dumr,at2,dumc,dumr ! first element has two and second has two character
          end if
        end if
        if (findloc(dumlist,at1,1) .eq. 0) then
          j = j+1
          dumlist(j) = at1
        end if
        if (findloc(dumlist,at2,1) .eq. 0) then
          j = j+1
          dumlist(j) = at2
        end if
      end if

      if (index(tmp,'LP ') .ne. 0) then
        backspace (ich)
        if (tmp(64:64) == ' ') then ! if element symbol has one character, there is a space and we use this routine
          read (ich,*) dumr,type,dumr,dumr,dumr,dumr,dumr,at1,dumc,dumc,dumr
        else ! if element symbol has two characters, there is no space and we use this routine
          read (ich,*) dumr,type,dumr,dumr,dumr,dumr,dumr,at1,dumc,dumr
        end if
        if (findloc(dumlist,at1,1) .eq. 0) then
          j = j+1
          dumlist(j) = at1
        end if
      end if
    end do
    close (ich)
    basicatlist = pack(dumlist,dumlist .ne. 0) ! sort out zeroes
    nbaseat = size(basicatlist)
    write (*,"('Protonation sites at atom positions:', *(i5))") (basicatlist(i),i=1,nbaseat)
    call remove('xtb.out')
  end subroutine readbasicpos

!=======================================================================================!

  subroutine fragment_structure(nat,oz,xyz,rcut,at1,at2,frag,fragcount)
!*****************************************************************************************
!* Taken and slightly modified from QCxMS mass spectra code https://github.com/qcxms/QCxMS
!*  Subroutine for definition of two or more fragments
!*  if at1 = 0 :  global separation in (nonbonded parts), beginning with atom at2
!*  if at1 and at2 <> 0 : define fragments only if a at1-at2 bond (and no other) exists
!*  if at1 and at2 = 0 : delete all fragment assignments
!*  no bond if rcut times the cov.dist.
!*  works better than mrec
!*****************************************************************************************
    use miscdata,only:Rad => QCxMS_Rad
    implicit none
    integer  :: at1,at2,nat
    integer  :: i,j
    integer  :: attotal,currentfrag
    integer  :: oz(nat),frag(nat)
    integer  :: fragcount
    real(wp),intent(in) ::  xyz(3,nat)
    real(wp) :: rcov,r
    real(wp) :: rcut

    logical  :: finish
    logical,allocatable  :: connect(:,:)

    allocate (connect(nat,nat))
    connect(1:nat,1:nat) = .false.

    do i = 1,nat-1
      do j = i+1,nat
        r = sqrt((xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2 &
        & +(xyz(3,i)-xyz(3,j))**2)
        rcov = rcut*0.5_wp*(Rad(oz(i))+Rad(oz(j)))
        if (r .lt. rcov) then
          connect(i,j) = .true.
          connect(j,i) = .true.
        end if
      end do
    end do
    if ((at1 .eq. 0).and.(at2 .eq. 0)) then
      do i = 1,nat
        frag(i) = 1
      end do
      return
    else

      do i = 1,nat
        frag(i) = 0
      end do

      frag(at1) = 1
      attotal = 1

      if (at2 .ne. 0) then
        connect(at1,at2) = .false.
        connect(at2,at1) = .false.
      end if

      finish = .false.
      currentfrag = 0

      do while (attotal .ne. nat)

        currentfrag = currentfrag+1

        ! cycle through atoms and find connected ones

        do while (.not. (finish))
          finish = .true.
          do i = 1,nat
            if (frag(i) .eq. currentfrag) then
              do j = 1,nat
                if (connect(i,j)) then
                  if (frag(j) .eq. 0) then
                    frag(j) = currentfrag
                    attotal = attotal+1
                    finish = .false.
                  elseif (frag(j) .eq. currentfrag) then
                    cycle
                  end if
                end if
              end do
            end if
          end do
        end do

        ! find the first atom in the next fragment

        do i = 1,nat
          if (frag(i) .eq. 0) then
            frag(i) = currentfrag+1
            attotal = attotal+1
            exit
          end if
        end do
        finish = .false.
      end do

    end if

    do i = 1,3 ! is enough, we only need to know we have 1, 2 or more than 2 fragments
      if (count(frag == i) .gt. 0) then
        fragcount = i
      end if
    end do

    deallocate (connect)
    return
  end subroutine fragment_structure

!=======================================================================================!

  subroutine sortoutduplicates(env,mso)
!******************************************************************************
!* sort out duplicates with molbar or inchi by openbabel (needs to be sourced)
!******************************************************************************
    implicit none
    type(msobj) :: mso
    integer ::  np
    integer :: i,j
    integer :: rm !removecounter
    character(len=1024),allocatable :: barcodes(:)
    character(len=2048) :: duplicates
    logical,allocatable :: double(:)
    type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
    logical :: lprint
    logical :: stat

    lprint = env%mslargeprint
    np = mso%pl%nmol
    allocate (double(np),barcodes(np))

    ! compute barcodes
    if (env%msmolbar) then
      call calcmolbar(env,mso,np,barcodes,stat)
      if (.not.stat) return
    elseif (env%msinchi) then
      call calcinchi(env,mso,np,barcodes,stat)
      if (.not.stat) return
    else
      return ! no topology check for duplicated structures selected
    end if

    ! identify duplicates
    double = .false.
    do i = 1,np
      if (double(i)) cycle
      duplicates = ""
      do j = i+1,np
        if (.not.double(j)) then
          if (trim(barcodes(i)) == trim(barcodes(j))) then
            double(j) = .true.
            if (lprint) write (duplicates,'(a,i0)') trim(duplicates)//" ",j
          else
            cycle
          end if
        end if
      end do
      if (lprint.and.len(trim(duplicates)) .gt. 0) then
        write (*,'(a,i0,a)') "removed duplicates of ",i,": "//trim(duplicates)
      end if
    end do

    ! remove duplicates
    rm = 0
    do i = 1,np
      rm = rm+1
      if (double(i)) then
        call mso%pl%remove(rm)
        rm = rm-1
      end if
    end do

    write (*,*) "sorted out ",np-mso%pl%nmol," fragments accoding to topology check"
    write (*,*) "number of products is now:",mso%pl%nmol-1

  end subroutine sortoutduplicates

!=======================================================================================!

  subroutine calcmolbar(env,mso,np,barcodes,stat)
!******************************************************
!* compute molbar codes for all structures with molbar
!******************************************************
    implicit none
    type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
    type(msobj) :: mso
    integer ::  np
    integer :: i,ich,iocheck,io
    character(len=:),allocatable :: prog
    character(len=:),allocatable :: subdir
    character(len=1024) :: barcodes(np)
    character(len=512) :: jobcall
    character(len=512) :: thisdir
    character(len=80) :: fname
    logical :: sourced,stat

    stat = .true.

    prog = 'molbar'
    call checkprog_silent(prog,.true.,io)
    if (io .ne. 0) then
      write (*,*) trim(prog)," not found, no topology check for duplicated structures possible"
      stat = .false.
      return
    end if

    call getcwd(thisdir)
    subdir = 'topodir'
    io = makedir(subdir)
    call chdir(subdir)

    do i = 1,np
      write (fname,'(i0,a)') i,"temp.xyz"
      open (newunit=ich,file=fname)
      call wrxyz(ich,mso%pl%mol(i)%nat,mso%pl%mol(i)%at,mso%pl%mol(i)%xyz(:,:),mso%pl%mol(i)%etot)
      close (ich)
    end do

    ! molbar on multiple core is faster but was unstable in some cases, single core more stable
    !write(jobcall,'(a)') 'molbar *temp.xyz -s  > molbar.out'

    write (jobcall,'(a,i0,a)') 'molbar *temp.xyz -s -T ',env%threads,' > molbar.out'
    write (*,*) "Calling molbar for sorting out duplicates by molbar"

    call execute_command_line(trim(jobcall))

    do i = 1,np
      write (fname,'(i0,a)') i,'temp.mb'
      open (newunit=ich,file=fname,status="old",action="read")
      read (ich,'(a)',iostat=iocheck) barcodes(i)
      close (ich)
    end do

    call chdir(thisdir)
    call rmrf('topodir')
  end subroutine calcmolbar

!=======================================================================================!

  subroutine calcinchi(env,mso,np,barcodes,stat)
!********************************************************
!* compute inchi codes for all structures with openbabel
!********************************************************
    implicit none
    type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
    type(msobj) :: mso
    integer ::  np
    integer :: i,ich,iocheck,io
    character(len=:),allocatable :: prog
    character(len=:),allocatable :: subdir
    character(len=1024) :: barcodes(np)
    character(len=512) :: jobcall
    character(len=512) :: thisdir
    character(len=80) :: fname
    logical :: sourced,stat

    stat = .true.

    prog = 'obabel'
    call checkprog_silent(prog,.true.,io)
    if (io .ne. 0) then
      write (*,*) trim(prog)," not found, no topology check for duplicated structures possible"
      stat = .false.
      return
    end if

    call getcwd(thisdir)
    subdir = 'topodir'
    io = makedir(subdir)
    call chdir(subdir)

    open (newunit=ich,file='temp.xyz')
    do i = 1,np
      call wrxyz(ich,mso%pl%mol(i)%nat,mso%pl%mol(i)%at,mso%pl%mol(i)%xyz(:,:),mso%pl%mol(i)%etot)
    end do
    close (ich)

    write (jobcall,'(a)') 'obabel -i xyz temp.xyz -o inchi > inchicodes 2>obabelout'
    write (*,*) "Calling obabel for sorting out duplicates by InChi codes"
    call execute_command_line(trim(jobcall))

    open (newunit=ich,file='inchicodes',status="old",action="read")
    do i = 1,np
      read (ich,"(a)") barcodes(i)
    end do
    close (ich)

    call chdir(thisdir)
    call rmrf('topodir')
  end subroutine calcinchi

!=======================================================================================!

  subroutine write_fragments(env,mso,estart,nisomers,nfragpairs,fname,lprint)
!*********************************************************
!* write fragments to files and optionally to directories
!*********************************************************
    use atmasses
    implicit none
    character(len=*) :: fname
    real(wp),intent(in) :: estart
    integer,intent(in) :: nisomers,nfragpairs
    type(systemdata) :: env
    type(msobj) :: mso
    integer :: prstruc ! number of structures to print
    integer :: incr
    integer :: np0,np,np2
    integer :: i,j,k,ich1,ich2,ich,r
    integer :: nat
    integer :: natf(2) ! number of atoms of fragments 1 and 2
    integer,allocatable :: atf(:,:) ! atomtypes of fragment
    character(len=1024) :: thisdir,isodir,fragdir
    character(len=40) :: strucname
    character(len=80) :: comment
    logical :: ex
    logical :: lprint
    character(len=40) :: sumform,sumformula
    real(wp) :: mass,erel

    nat = mso%pl%mol(1)%nat
    call getcwd(thisdir)
    np = mso%pl%nmol
    np0 = np
    ! first write to isomers.xyz and fragmentpairs.xyz and filter out
    write (*,*) "writing isomers to <isomers.xyz> and fragmentpairs to <fragmentpairs.xyz"
    open (newunit=ich1,file='isomers.xyz',status='replace')
    open (newunit=ich2,file='fragmentpairs.xyz',status='replace')
    i = 0
    do
      i = i+1
      if (i .gt. np) exit ! np changes in the loop
      if (mso%pl%mol(i)%fragcount .eq. 2) then ! fragmentpairs
        call wrxyz(ich2,mso%pl%mol(i)%nat,mso%pl%mol(i)%at,mso%pl%mol(i)%xyz(:,:),mso%pl%mol(i)%etot)
        if (env%msiso) then
          call mso%pl%remove(i)
          i = i-1
          np = np-1
        end if
      elseif (mso%pl%mol(i)%fragcount .eq. 1) then ! isomers
        call wrxyz(ich1,mso%pl%mol(i)%nat,mso%pl%mol(i)%at,mso%pl%mol(i)%xyz(:,:),mso%pl%mol(i)%etot)
        if (env%msnoiso) then
          call mso%pl%remove(i)
          i = i-1
          np = np-1
        end if
      end if
    end do
    close (ich1)
    close (ich2)
    allocate (atf(2,mso%pl%mol(1)%nat))
    ! print to output, directories and products.xyz
    write (*,*) "writing products to <"//trim(fname)//">"
    if (lprint) write (*,*) "writing product structures and molecular masses to directories"

    prstruc = env%msnfrag
    if (prstruc .gt. 0) then
      write (*,*) "Printing only ",prstruc," selected structures to products.xyz"
      incr = np/prstruc
    else
      incr = 1
      prstruc = np
    end if

    write (*,'(a)') '========================================================'
    write (*,'(a,i0)') " number of printed structures: ",prstruc ! todo changes depending on print settings
    write (*,*) " directory | fragment type | rel. energy [kcal/mol]"
    if (.not.env%msiso) write (*,*) " fragmentpair: | sumformula | molecular mass"

    np2 = 0 ! number of printed sturctures
    open (newunit=ich1,file=fname,status='replace')
    do i = 1,np,incr
      if (np2 .ge. prstruc) exit
      erel = (mso%pl%mol(i)%etot-estart)*autokcal
      if (mso%pl%mol(i)%fragcount .lt. 1) write (*,*) "no fragment count for ",i
      if (mso%pl%mol(i)%fragcount .eq. 1) then
        np2 = np2+1
        write (isodir,'(a,i0)') "p",np2
        write (*,'(a4,a20,f18.8)') trim(isodir),"isomer",erel
        strucname = 'isomer.xyz'
      elseif (mso%pl%mol(i)%fragcount .gt. 1) then
        np2 = np2+1
        write (isodir,'(a,i0)') "p",np2
        write (*,'(a4,a20,f18.8)') trim(isodir),"fragmentpair",erel
        strucname = 'pair.xyz'
      end if
      write (comment,'(f9.5,1x,a)') mso%pl%mol(i)%etot,trim(isodir)
      call wrxyz(ich1,mso%pl%mol(i)%nat,mso%pl%mol(i)%at,mso%pl%mol(i)%xyz(:,:),trim(comment))

      ! write to directories
      if (lprint) then
        r = makedir(trim(isodir))
        call chdir(trim(isodir))
        open (newunit=ich,file=strucname,status='replace')
        call wrxyz(ich,mso%pl%mol(i)%nat,mso%pl%mol(i)%at,mso%pl%mol(i)%xyz(:,:),mso%pl%mol(i)%etot)
        close (ich)
        ! just write mass for pairs too
        mass = 0
        do k = 1,nat
          mass = mass+ams(mso%pl%mol(i)%at(k))
        end do
        call wrshort_real("molmass",mass)
        call chdir(trim(thisdir))
      end if
      ! no write separated fragment structures
      if (mso%pl%mol(i)%fragcount .gt. 1) then
        do j = 1,2 ! loop over fragment 1 and 2
          natf(j) = count(mso%pl%mol(i)%fragi == j)
          mass = 0
          atf(j,:) = 0
          do k = 1,nat
            if (mso%pl%mol(i)%fragi(k) == j) then
              mass = mass+ams(mso%pl%mol(i)%at(k))
              atf(j,k) = mso%pl%mol(i)%at(k)
            end if
          end do
          sumformula = sumform(nat,atf(j,:))
          ! write to directories
          if (lprint) then
            write (fragdir,'(a,i0)') trim(isodir)//"f",j
            r = makedir(trim(fragdir))
            call chdir(fragdir)
            strucname = 'fragment.xyz'
            open (newunit=ich,file=strucname,status='replace')
            write (ich,*) natf(j)
            write (ich,*)
            do k = 1,nat
              if (mso%pl%mol(i)%fragi(k) == j) then !only write xyz if fragment really exists
                write (ich,'(a2,5x,3F18.8)') i2e(mso%pl%mol(i)%at(k)),mso%pl%mol(i)%xyz(1,k),mso%pl%mol(i)%xyz(2,k),mso%pl%mol(i)%xyz(3,k)
              end if
            end do
            call wrshort_real("mass",mass)
            call chdir(trim(thisdir))
          end if
          !write(*,*) " directory | sumformula | mass"
          write (*,'(a6,i0,5x,a20,9x,f9.5)') trim(isodir)//"f",j,trim(sumformula),mass
        end do
      end if
    end do

    write (*,'(a)') '========================================================'
    write (*,*)
    call wrshort_int('npairs',np2)
    if (env%msnoiso) then
      write (*,*) "sorted out ",np0-np2,"non-dissociated structures"
    elseif (env%msiso) then
      write (*,*) "sorted out ",np0-np2,"dissociated structure pairs"
    end if
    write (*,*) "Number of generated isomers: ",nisomers
    write (*,*) "Number of generated fragmentpairs: ",nfragpairs
  end subroutine write_fragments

!=======================================================================================!

  subroutine msinputreader(mso,msinput)
!**********************************************************
!* a simple input reader to read in additional input file
!* to set more advanced options
!**********************************************************
    implicit none
    type(msobj) :: mso
    integer :: ich,io,iocheck
    character(len=40) :: line
    character(len=80) :: msinput
    real(wp) :: xx(10)
    integer :: j
    logical :: ex

    write (*,*) "msinput is:,",trim(msinput)
    !---- read input file
    inquire (file=trim(msinput),exist=ex)
    if (.not.ex) then
      return
    else
      write (*,*) 'Reading <'//trim(msinput)//'> file'
    end if

    mso%rcut = 1.3_wp        !> cutoff for factor of sum of covalent radii for detection of fragmentation of fragment_structures
    mso%T = 5000.0_wp        !> electronic temperature! favours open-shell fragment generation
    mso%fc = 0.05_wp         !> repulsive force constant in constrained optimizations
    mso%fc_attr = -0.05_wp   !> attractive force constant in constrained optimizations
    mso%cdist = 1.5_wp       !> constraing distance scaling factor of rcov
    mso%cdist_att = 0.5_wp   !> constraing distance scaling factor of rcov for attractive part
    mso%distthr_att = 4.0_wp !> distance threshold in Angstrom for attractive part of constraint
    mso%atomshift = 0.75_wp  !> shift of atoms in random atom displacement
    mso%fragdist = 0.0_wp    !> increase distance (in Angstrom) to separate fragments (useful for subsequent transition state search)
    mso%maxc = 15

    open (newunit=ich,file=trim(msinput),status='old')

    ! Read line-by-line
    do
      read (ich,'(a)',iostat=iocheck) line

      ! Check for errors in the input
      if (iocheck > 0) then     !Fail
        write (*,*) 'Something is wrong in the input. Exiting...'
        stop

        ! End-of-file
      elseif (iocheck < 0) then !EOF
        exit

        ! Keywords
      else
        write (*,'(''>'',a)') line
        ! capatilize and trimthe line
        line = uppercase(line)
        line = trim(line)

        ! read keywords
        if (index(line,'FRAGDIST') /= 0) then
          call readl(trim(line),xx,j)
          mso%fragdist = xx(1)
        end if
        if (index(line,'ATOMSHIFT') /= 0) then
          call readl(trim(line),xx,j)
          mso%atomshift = xx(1)
        end if
        if (index(line,'DISTTHR_ATTR') /= 0) then
          call readl(trim(line),xx,j)
          mso%distthr_att = xx(1)
        end if
        if (index(line,'FC_REP') /= 0) then
          call readl(trim(line),xx,j)
          mso%fc = xx(1)
        end if
        if (index(line,'FC_ATTR') /= 0) then
          call readl(trim(line),xx,j)
          mso%fc_attr = xx(1)
        end if
        if (index(line,'ETEMP') /= 0) then
          call readl(trim(line),xx,j)
          mso%T = xx(1)
        end if
      end if

    end do

  end subroutine msinputreader

!=======================================================================================!

  subroutine detect_fragments(mso,lprint,nisomers,nfragpairs)
!*******************************************************************
! detection of fragments and removal of multiple fragmented species
!*******************************************************************
    implicit none
    type(msobj) :: mso
    integer :: nat
    integer :: i
    integer :: nfragpairs,nisomers
    integer :: npoly ! number of multiple fragmented species
    integer :: np ! number of products
    character(len=2048) :: multfrags
    logical,intent(in) :: lprint

    write (*,*) "Detect fragments"

    nat = mso%pl%mol(1)%nat
    np = mso%pl%nmol

    nisomers = 0
    nfragpairs = 0
    npoly = 0

    ! determine fragments for every structure
    do i = 1,np
      allocate (mso%pl%mol(i)%fragi(nat))
      call fragment_structure(mso%pl%mol(i)%nat,mso%pl%mol(i)%at,mso%pl%mol(i)%xyz(:,:),mso%rcut,1,0,mso%pl%mol(i)%fragi,mso%pl%mol(i)%fragcount) ! works better than mrec
    end do
    ! remove multiple fragmented species
    multfrags = ''
    i = 0
    do
      i = i+1
      if (i .gt. np) exit ! np changes within the loop
      if (mso%pl%mol(i)%fragcount .gt. 2) then
        npoly = npoly+1
        if (lprint) write (multfrags,'(a,i0)') trim(multfrags)//" ",i
        call mso%pl%remove(i)
        i = i-1
        np = np-1
        cycle
      elseif (mso%pl%mol(i)%fragcount .eq. 2) then
        nfragpairs = nfragpairs+1
      else ! for isomers do nothing
        nisomers = nisomers+1
        cycle
      end if
    end do
    if (lprint.and.len(trim(multfrags)) .gt. 0) then
      write (*,'(a,i0,a)') "removed multiple fragmented structures: "//trim(multfrags)
    end if
    return
  end subroutine detect_fragments

!=======================================================================================!

  subroutine increase_fragdist(mso,lprint)
!****************************************************************************
!* increase distance between fragments, better for transistion state search
!* we remove here also multiple fragmented species,
!* and unreasonable fragments (too much distance between fragments)
!****************************************************************************
    use axis_module
    implicit none
    type(msobj) :: mso
    integer :: nat
    integer,allocatable :: at1(:),at2(:)
    integer,allocatable :: fragi(:) !fragment index
    integer :: frag1,frag2,i,j
    integer :: npoly ! number of multiple fragmented species
    integer :: np ! number of products
    real(wp) :: norm
    real(wp),allocatable :: xyz1(:,:),xyz2(:,:) ! xyz of fragment 1 and 2
    real(wp) :: cmass1(3),cmass2(3) ! center of mass of fragment 1 and 2
    character(len=2048) :: distfrags
    logical,intent(in) :: lprint

    frag1 = 0
    frag2 = 0
    if (mso%fragdist .gt. 0.0_wp) then
      write (*,'(a,f10.8,a)') " Increase distance of fragments by ",mso%fragdist," Angstrom"
    else
      return
    end if

    nat = mso%pl%mol(1)%nat
    np = mso%pl%nmol

    ! increase distance between fragments (based on center of mass) and remove unreasonable fragments
    ! that are too far away from each other, which cause problems in transition state search

    distfrags = ''
    i = 0
    do
      i = i+1
      if (i .gt. np) exit ! np changes within the loop
      ! increase distance between fragments
      if (mso%pl%mol(i)%fragcount .eq. 2) then
        do j = 1,nat
          if (mso%pl%mol(i)%fragi(j) == 1) then
            frag1 = frag1+1
          end if
          if (mso%pl%mol(i)%fragi(j) == 2) then
            frag2 = frag2+1
          end if
        end do
        allocate (at1(frag1),at2(frag2))
        allocate (xyz1(3,frag1),xyz2(3,frag2))
        frag1 = 0
        frag2 = 0
        do j = 1,nat
          if (mso%pl%mol(i)%fragi(j) == 1) then
            frag1 = frag1+1
            xyz1(:,frag1) = mso%pl%mol(i)%xyz(:,j)
            at1(frag1) = mso%pl%mol(i)%at(j)
          end if
          if (mso%pl%mol(i)%fragi(j) == 2) then
            frag2 = frag2+1
            xyz2(:,frag2) = mso%pl%mol(i)%xyz(:,j)
            at2(frag2) = mso%pl%mol(i)%at(j)
          end if
        end do
        call CMAv(frag1,at1,xyz1,cmass1)
        call CMAv(frag2,at2,xyz2,cmass2)
        norm = sqrt((cmass2(1)-cmass1(1))**2+(cmass2(2)-cmass1(2))**2+(cmass2(3)-cmass1(3))**2)
        ! remove structures with unreasonable distance between fragments
        if (norm .gt. 15) then !15 angstrom worked well in tests
          if (lprint) write (distfrags,'(a,i0)') trim(distfrags)//" ",i
          call mso%pl%remove(i)
          i = i-1
          np = np-1
          cycle
        end if
        do j = 1,nat
          ! move fragment 2 away from fragment 1
          if (mso%pl%mol(i)%fragi(j) == 2) then
            mso%pl%mol(i)%xyz(:,j) = mso%pl%mol(i)%xyz(:,j)+(cmass2-cmass1)/norm*mso%fragdist
          end if
        end do
        deallocate (at1,at2,xyz1,xyz2)
      else ! for isomers do nothing
        cycle
      end if
    end do

    if (lprint.and.len(trim(distfrags)) .gt. 0) then
      write (*,'(a,i0,a)') "removed unreasonable structures:"//trim(distfrags)
    end if

    return
  end subroutine increase_fragdist

!=======================================================================================!

  subroutine rdplist(mso,estart,fname)
!*******************************************************
!* collect structures of ensemble fil into product list
!*******************************************************
    use atmasses
    implicit none
    integer :: nat
    integer :: nall
    integer :: chrg ! charge of the molecule
    type(msobj) :: mso
    character(len=*) :: fname
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    real(wp) :: eread
    real(wp) :: molmass
    real(wp) :: estart
    character(len=6) :: sym
    character(len=512) :: line
    integer :: i,j,k,dum,io,ich
    logical :: ex,found

    found = .false.
    write (*,*) "reading product structures file ",fname

    call rdensembleparam(fname,nat,nall)
    allocate (at(nat),xyz(3,nat))
    call rdshort('.CHRG',chrg)
    eread = 0.0_wp
    xyz = 0.0_wp
    open (newunit=ich,file=fname)
    do i = 1,nall
      read (ich,*,iostat=io) dum
      if (io < 0) exit
      if (io > 0) cycle
      read (ich,'(a)',iostat=io) line
      if (io < 0) exit
      eread = grepenergy(line)
      if (abs(eread-estart) .lt. 0.00001_wp) then ! skip input structure, detect by energy
        found = .true.
        do k = 1,dum
          read (ich,'(a)',iostat=io) line
        end do
        cycle
        if (io < 0) exit
      end if
      do j = 1,dum
        read (ich,'(a)',iostat=io) line
        if (io < 0) exit
        call coordline(line,sym,xyz(1:3,j),io)
        if (io .ne. 0) then
          backspace (ich)
          exit
        end if
        at(j) = e2i(sym)
      end do
      molmass = molweight(nat,at)
      call mso%pl%append(nat,at,xyz,eread,chrg,0,molmass)
    end do
    close (ich)

    if (io < 0) then
      error stop 'error while reading product file.'
    end if

    deallocate (xyz,at)
    if (.not.found) then
      write (*,*) "Warning, input structure not found in ensemble and could not be sorted out"
    end if
    return
  end subroutine rdplist

!=======================================================================================!

  subroutine wrplist(mso,fname)
!*************************
!* write products to file
!*************************
    implicit none
    integer :: i,ich,io
    character(len=*) :: fname
    type(msobj) :: mso
    open (newunit=ich,file=fname)
    do i = 1,mso%pl%nmol
      call wrxyz(ich,mso%pl%mol(i)%nat,mso%pl%mol(i)%at,mso%pl%mol(i)%xyz(:,:),mso%pl%mol(i)%etot)
    end do
    close (ich)
    return
  end subroutine wrplist

!=======================================================================================!
!=======================================================================================!
end module msmod  !> END OF MODULE
!=======================================================================================!
!=======================================================================================!
