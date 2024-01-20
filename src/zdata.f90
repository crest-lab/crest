!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2022 Philipp Pracht
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
! The zdata module contains several polymorphic datatypes that
! are related to topology analysis routines within crest.
! they are used, e.g., in the Entropy calculation or in the
! identification of isomers
!================================================================================!
module zdata
  use crest_parameters
  implicit none

  public :: zatom
  public :: zring
  public :: zgrp
  public :: zequal
  public :: zmolecule
  public :: zensemble

  public :: readwbo

  public :: maxringident

  private

!=========================================================================================!
!=========================================================================================!

  integer :: maxringident = 100  !limit for ring sizes identified in the topology routines
  !this is necessary for large intertwined systems, e.g. proteins

!=========================================================================================!

  !>--- a single atom and all the data involved, the atom tracks its own neighbours
  type :: zatom
    integer :: atype        !> integer element
    integer :: pos = 0      !> position in the list of atoms
    integer :: nei = 0      !> # of neighbouring (bonded) atoms
    logical :: ring = .false.       !> is atom part of a ring?
    logical :: multiring = .false.  !> is the atom in a polycyclic ring system (e.g. steroids)?
    logical :: stereo = .false.     !> is atom a stereocenter?
    character(len=3) :: rs = ''     !> R or S
    real(wp) :: cn = 0.0_wp         !> the coordination number
    real(wp) :: p = 0.0_wp          !> partial charge
    !>--- arrays
    real(wp) :: cart(3)                !> cartesian coordinates x,y,z
    integer,allocatable  :: ngh(:)     !> list of neighbouring atoms
    integer,allocatable  :: ngt(:)     !> atom types of the neighbouring atoms
    integer,allocatable  :: prio(:)    !> priority of neighbor chains for Cahn-Ingold rules
    logical :: prioset = .false.       !> is the priority determined
    real(wp),allocatable :: nweig(:)   !> weight of the neighbouring chains
    character(len=:),allocatable :: el !> element symbol
  contains
    procedure :: allocate => allocate_zatom     !> allocate memory for atom
    procedure :: deallocate => deallocate_zatom !> deallocate memory
    procedure :: wrtmline => wrtmline !> write a single line of TM coordinates
    procedure :: wrneigh => wrneigh   !> write a list of neighbours
    procedure :: wrcn => wrcn         !> write CN
  end type zatom

!=========================================================================================!

  !>--- ring subtype. contains ring(-system) size and involved atoms
  type :: zring
    integer :: rs                   !> # atoms in ring
    integer,allocatable :: rlist(:) !> atoms in ring
  contains
    procedure :: deallocate => deallocate_zring  !deallocate (if allocated)
    procedure :: print => print_zring !print to screen
  end type zring

!=========================================================================================!

  !>--- object for a single group
  type :: zgrp
    integer :: nm !> # of members in group
    integer,allocatable :: mem(:)
  contains
    procedure :: allocate => allocate_zgrp
    procedure :: deallocate => deallocate_zgrp
    procedure :: append => append_zgrp
    procedure :: prgrp => print_zgrp
  end type zgrp

  !>--- equivalencies-class. contains information about equivalent nuclei
  type :: zequal
    integer :: nat
    integer :: ng = 0  !> number of groups of equivalent atoms (initialized to 0)
    integer :: eng = 0 !> number of groups with more than 1 nuclei
    integer,allocatable :: ord(:)    !> group order. which nuclei belongs to which group
    type(zgrp),allocatable :: grp(:) !> all groups (dimension <= nat)

  contains
    procedure :: allocate => allocate_zequal !> allocate nat, ord(:) and grp(:)
    procedure :: member => is_x_member       !> check if atom x is member of any group
    procedure :: prsum => prsummary_zequal   !> print a summary of all the groups
    procedure :: geteng => zequal_geteng     !> count number of groups with more than 1 member
  end type zequal

!=========================================================================================!

  !>--- molecule class. contains Nat, all atoms as "zatom" type, the number of rings
  type :: zmolecule
    !>--- regular atom
    integer :: nat                       !> number of atoms
    real(wp),allocatable :: xyz(:,:)     !> coordinates (in a.u.)
    type(zatom),allocatable :: zat(:)    !> atoms, each atom is an object, see type(atom)
    integer,allocatable :: at(:)         !> for quick access of atom type
    integer :: maxneigh = 0              !> maximum number of neighbors
    integer :: limitnei = 6              !> limit for highly coordinated atoms
    integer :: nonterm = 0               !> number of non-terminal atoms
    integer :: nfrag = 1                 !> number of fragments in the system (for NCI complexes)
    integer,allocatable :: molvec(:)     !> fragment tracking
    real(wp),allocatable :: zmat(:,:)    !> internal coordinates of the molecule,
    !> requires tracking zna,znb,znc
    integer,allocatable :: zna(:)
    integer,allocatable :: znb(:)
    integer,allocatable :: znc(:)
    real(wp) :: chrg = 0.0_wp   !> molecular charge
    real(wp) :: uhf = 0.0_wp   !> nα-nβ electrons
    real(wp) :: energy = 0.0_wp   !> energy

    !>--- bonds
    integer :: nb = 0
    integer,allocatable :: bondpairs(:,:)

    !>--- data concerning molecular rings
    integer :: nri = 0                   !> number of rings
    type(zring),allocatable :: zri(:)    !> rings
    integer :: maxring = 0               !> maximum ring size
    logical :: polycycle = .false.       !> is there any polycycle in the molecule?

    !>--- matirces for atom pairs
    real(wp),allocatable :: distmat(:,:) !> atom-atom distances in Cartesian space
    real(wp),allocatable :: wbo(:,:)     !> atom-atom bond orders

    !>--- stereoinformation for the molecule
    integer :: nstereo = 0                !> number of stereocenters
    integer,allocatable :: stereotrac(:)  !> track the order of stereocenters
    integer,allocatable :: inverter(:,:)  !> book keeping of inverting atoms at the resp. stereocenter
    integer,allocatable :: invector(:,:)  !> inversion vector

    !>--- atom maps (required if the object is a fragment)
    integer,allocatable :: map(:)
    integer,allocatable :: revmap(:)

    !>--- QCG information
    integer   :: nmol          !> number of molecules
    real(wp)  :: cma(3)        !> center of mass
    real(wp)  :: aniso         !> anisotropy factor
    real(wp)  :: ell_abc(3)    !> ellipsoid axis
    real(wp)  :: atot          !> surface area
    real(wp)  :: vtot          !> volume
    real(wp)  :: rtot          !> radius
    real(wp)  :: mass          !> mass
    real(wp)  :: gt            !> gibbs free energy
    real(wp)  :: ht            !> enthalpy
    real(wp)  :: svib          !> vibrational entropy
    real(wp)  :: srot          !> rotational entropy
    real(wp)  :: stra          !> translational entropy
    real(wp)  :: eax(3)        !> molecular axis

    !>--- procedures to be used with the zmol type
  contains
    procedure :: wrtable => wrtable          !> write CNs and neighbours
    procedure :: prsym => prsym              !> print element and number of atom i
    procedure :: prrings => printrings_zmol  !> print a list of all unique rings
    procedure :: deallocate => deallocate_zmol  !> deallocate everything
    procedure :: dist => zatdist             !> distance between two atoms
    procedure :: mrec => zmol_mrec           !> mrec routine directly acting upon neighbour list
    procedure :: getxyz => zmol_getxyz
    procedure :: hydrogen => count_hydrogen  !> function to retun number of hydrogen atoms
    procedure :: countbonds => count_bonds   !> cound bonds and track which atoms from them
    procedure :: methyl => is_methyl         !> check if a given atom is a methyl C
    procedure :: fragment => get_fragment    !> create a new zmolecule object for the i-th fragment
    procedure :: addring => zmol_add_ring    !> add a ring to the object
    procedure :: adjacency => zmol_adjacency !> generate adjacency matrix from zmol
  end type zmolecule

!=========================================================================================!

  !>--- ensemble class. contains all structures of an ensemble
  type :: zensemble
    !>--- data
    integer :: nat = 0
    integer :: nall = 0
    !> vnat is used instead of nat if not all structures have the same number of atoms
    !> in this case nat will be  =maxval(vnat,1)
    integer,allocatable :: vnat(:)

    integer,allocatable  :: at(:)      !> atom types as integer, dimension will be at(nat)
    real(wp),allocatable :: xyz(:,:,:) !> coordinates, dimension will be xyz(3,nat,nall)
    real(wp),allocatable :: eread(:)   !> energy of each structure, dimension will be eread(nall)

    integer :: nconf                 !> #conformers
    integer,allocatable :: nrot(:)   !> #rotamers per conformer, dimension will be nrot(nconf)
    integer :: nrotmax               !> maximum number of rotamers in the enesemble
    integer,allocatable :: nmat(:,:) !> track which rotamer belongs to which conformer, dimension will be nmat(nconf,nrotmax)

    !>--- procedures to be used with the zensemble
  contains
    procedure :: allocate => allocate_zensemble
    procedure :: deallocate => deallocate_zensemble
  end type zensemble

!=========================================================================================!
!=========================================================================================!
contains !> MODULE PROCEDURES START HERE
!=========================================================================================!
!=========================================================================================!

  subroutine allocate_zatom(self)
    implicit none
    class(zatom) :: self
    integer :: a
    a = self%nei
    allocate (self%ngh(a),source=0)
    return
  end subroutine allocate_zatom
!=========================================================================================!
  subroutine deallocate_zatom(self)
    implicit none
    class(zatom) :: self
    if (allocated(self%ngh)) deallocate (self%ngh)
    if (allocated(self%ngt)) deallocate (self%ngt)
    if (allocated(self%prio)) deallocate (self%prio)
    if (allocated(self%nweig)) deallocate (self%nweig)
    return
  end subroutine deallocate_zatom
!=========================================================================================!
  subroutine deallocate_zmol(self)
    implicit none
    class(zmolecule) :: self
    integer :: i
    if (allocated(self%at)) deallocate (self%at)
    if (allocated(self%xyz)) deallocate (self%xyz)
    if (allocated(self%distmat)) deallocate (self%distmat)
    if (allocated(self%wbo)) deallocate (self%wbo)
    if (allocated(self%stereotrac)) deallocate (self%stereotrac)
    if (allocated(self%inverter)) deallocate (self%inverter)
    if (allocated(self%invector)) deallocate (self%invector)
    if (allocated(self%map)) deallocate (self%map)
    if (allocated(self%revmap)) deallocate (self%revmap)

    if (allocated(self%zat)) then
      do i = 1,self%nat
        call self%zat(i)%deallocate()
      end do
      deallocate (self%zat)
    end if

    if (allocated(self%zri)) then
      do i = 1,self%nri
        call self%zri(i)%deallocate()
      end do
      deallocate (self%zri)
    end if

    return
  end subroutine deallocate_zmol
!=========================================================================================!
!> write a single line of coordinates
  subroutine wrtmline(self,ch)
    implicit none
    class(zatom) :: self
    integer :: ch
    write (ch,'(3F24.12,5x,a)') self%cart(1:3),self%el
    return
  end subroutine wrtmline
!=========================================================================================!
!> write a single line of coordinates
  subroutine wrneigh(self,ch)
    implicit none
    class(zatom) :: self
    integer :: ch
    integer :: i
    write (ch,'(1x,a,a,i0,a,5x)',advance='no') self%el,'(',self%pos,')'
    if (self%nei .lt. 1) then
      write (ch,'(a)') 'no neighbours!'
    else
      do i = 1,self%nei
        write (ch,'(i0,1x)',advance='no') self%ngh(i)
      end do
      write (ch,*)
    end if
    return
  end subroutine wrneigh
!=========================================================================================!
! wcheck if a specific atom is a emthyl (or similar) group
  function is_methyl(self,atm) result(itis)
    implicit none
    class(zmolecule) :: self
    integer :: atm
    logical :: itis
    integer :: i,j,c
    integer :: a,b
    itis = .false.
    if (self%zat(atm)%nei .lt. 4) then
      itis = .false.
      return
    else
      do i = 1,self%zat(atm)%nei
        c = 0
        do j = 1,self%zat(atm)%nei
          if (i == j) cycle
          a = self%zat(atm)%ngh(i)
          b = self%zat(atm)%ngh(j)
          if (self%at(a) == self%at(b)) then
            if (self%zat(a)%nei == 1.and. &
            &   self%zat(b)%nei == 1) c = c+1
          end if
        end do
        if (c == 2) itis = .true.
      end do
    end if
    return
  end function is_methyl
!=========================================================================================!
! write a single line
  subroutine wrcn(self,ch)
    implicit none
    class(zatom) :: self
    integer :: ch
    write (ch,'(1x,a,a,i0,a,5x)',advance='no') self%el,'(',self%pos,')'
    write (ch,'(f8.4)') self%cn
    return
  end subroutine wrcn
!=========================================================================================!
!write a table of CNs and neighbouring atoms for ZMOLECULE class
  subroutine wrtable(self,ch)
    implicit none
    class(zmolecule) :: self
    integer :: ch
    integer :: i,j
    character(len=10) :: dum
    write (ch,'(13x,a,9x,a)') 'CN','neighbours'
    do i = 1,self%nat
      write (dum,'(a,a,i0,a)') self%zat(i)%el,'(',i,')'
      write (ch,'(1x,a10,f8.4,5x)',advance='no') dum,self%zat(i)%cn
      if (self%zat(i)%nei .lt. 1) then
        write (ch,'(a)') 'no neighbours!'
      else
        do j = 1,self%zat(i)%nei
          write (ch,'(i0,1x)',advance='no') self%zat(i)%ngh(j)
        end do
        write (ch,*)
      end if
    end do
    return
  end subroutine wrtable
!=========================================================================================!
!> calculate the distance between two atoms i and j
  function zatdist(self,i,j) result(dist)
    implicit none
    class(zmolecule) :: self
    integer :: i,j
    real(wp) :: dist
    dist = 0.0_wp
    if (allocated(self%distmat)) then
      dist = self%distmat(i,j)
    else
      dist = (self%zat(i)%cart(1)-self%zat(j)%cart(1))**2+ &
     &     (self%zat(i)%cart(2)-self%zat(j)%cart(2))**2+ &
     &     (self%zat(i)%cart(3)-self%zat(j)%cart(3))**2
      dist = sqrt(dist)
    end if
    return
  end function zatdist
!=========================================================================================!
!> write a table of CNs and neighbouring atoms for ZMOLECULE class
  function prsym(self,i)
    implicit none
    class(zmolecule) :: self
    integer :: i
    character(len=20) :: prsym
    write (prsym,'(a,a,i0,a)') trim(self%zat(i)%el),'(',i,')'
    return
  end function prsym
!=========================================================================================!
  subroutine printrings_zmol(self,ch)
    implicit none
    class(zmolecule) :: self
    integer :: ch
    integer :: i
    if (allocated(self%zri)) then
      do i = 1,self%nri
        call self%zri(i)%print(ch)
      end do
    end if
    return
  end subroutine printrings_zmol
!=========================================================================================!
!> create a new zmolecule object for the i-th fragment
  subroutine get_fragment(self,i,znew)
    implicit none
    class(zmolecule) :: self
    type(zmolecule) :: znew
    integer :: i,natnew
    integer :: j,k,l,m
    integer,allocatable :: map(:)
    real(wp),allocatable :: xyz(:,:)
    if (i > self%nfrag) return
    call znew%deallocate()
    allocate (map(self%nat),source=0)
    natnew = 0
    k = 0
    do j = 1,self%nat
      if (self%molvec(j) == i) then
        k = k+1
        map(j) = k
        natnew = natnew+1
      end if
    end do
    znew%nat = natnew
    allocate (znew%zat(znew%nat),znew%at(znew%nat))
    allocate (znew%map(self%nat),znew%revmap(znew%nat))
    znew%map = map
    do j = 1,self%nat
      if (map(j) > 0) znew%revmap(map(j)) = j
    end do
    k = 0
    do j = 1,self%nat
      if (self%molvec(j) == i) then
        k = k+1
        znew%zat(k) = self%zat(j)
        znew%at(k) = self%at(j)
        do l = 1,znew%zat(k)%nei
          m = znew%zat(k)%ngh(l)
          znew%zat(k)%ngh(l) = map(m)
        end do
      end if
    end do
    allocate (znew%distmat(znew%nat,znew%nat),source=0.0_wp)
    allocate (xyz(3,znew%nat))
    call znew%getxyz(xyz)
    do j = 1,znew%nat
      do k = 1,znew%nat-1
        znew%distmat(j,k) = (xyz(1,j)-xyz(1,k))**2+ &
      &                   (xyz(2,j)-xyz(2,k))**2+ &
      &                   (xyz(3,j)-xyz(3,k))**2
        znew%distmat(j,k) = sqrt(znew%distmat(j,k))
        znew%distmat(k,j) = znew%distmat(j,k)
      end do
    end do
    deallocate (xyz)
    deallocate (map)
    return
  end subroutine get_fragment
!=========================================================================================!
  subroutine zmol_add_ring(self,ring)
    implicit none
    class(zmolecule) :: self
    type(zring) :: ring
    type(zring),allocatable :: rnew(:)
    integer :: i,k
    k = self%nri+1
    allocate (rnew(k))
    do i = 1,self%nri
      rnew(i) = self%zri(i)
    end do
    rnew(k) = ring
    if (ring%rs > self%maxring) self%maxring = ring%rs
    call move_alloc(rnew,self%zri)
    self%nri = k
    return
  end subroutine zmol_add_ring
!=========================================================================================!
  subroutine deallocate_zring(self)
    implicit none
    class(zring) :: self
    if (allocated(self%rlist)) deallocate (self%rlist)
    return
  end subroutine deallocate_zring
!=========================================================================================!
  subroutine print_zring(self,ch)
    implicit none
    class(zring) :: self
    integer :: ch
    character(len=:),allocatable :: atmp
    character(len=30) :: btmp
    integer :: i
    atmp = ' atoms:'
    if (allocated(self%rlist)) then
      write (ch,'(1x,a,i0)') 'ring size: ',self%rs
      do i = 1,self%rs
        write (btmp,'(1x,i0)') self%rlist(i)
        atmp = trim(atmp)//trim(btmp)
      end do
      write (ch,*) trim(atmp)
    end if
    deallocate (atmp)
    return
  end subroutine print_zring
!=========================================================================================!
  subroutine allocate_zgrp(self,nm)
    implicit none
    class(zgrp) :: self
    integer :: nm
    if (allocated(self%mem)) deallocate (self%mem)
    self%nm = nm
    allocate (self%mem(nm),source=0)
    return
  end subroutine allocate_zgrp
!=========================================================================================!
  subroutine deallocate_zgrp(self)
    implicit none
    class(zgrp) :: self
    if (allocated(self%mem)) deallocate (self%mem)
    self%nm = 0
    return
  end subroutine deallocate_zgrp
!=========================================================================================!
  subroutine append_zgrp(self,obj)
    implicit none
    class(zgrp) :: self
    integer :: obj
    integer,allocatable :: tmp(:)
    if (allocated(self%mem)) then
      if (any(self%mem(:) .eq. obj)) return !return if already present
      allocate (tmp(self%nm+1))
      tmp(1:self%nm) = self%mem(1:self%nm)
      deallocate (self%mem)
      self%nm = self%nm+1
      tmp(self%nm) = obj
      call move_alloc(tmp,self%mem)
    else
      call self%allocate(1)
      self%mem(1) = obj
      return
    end if
    return
  end subroutine append_zgrp
!=========================================================================================!
  function print_zgrp(self)
    implicit none
    class(zgrp) :: self
    character(len=:),allocatable :: print_zgrp
    character(len=16) :: dum
    integer :: i
    print_zgrp = 'group :'
    if (allocated(self%mem)) then
      do i = 1,self%nm
        write (dum,'(1x,i0)') self%mem(i)
        print_zgrp = print_zgrp//trim(dum)
      end do
    end if
    return
  end function print_zgrp
!=========================================================================================!
  subroutine allocate_zequal(self,nat)
    implicit none
    class(zequal) :: self
    integer :: nat
    self%nat = nat
    if (allocated(self%ord)) deallocate (self%ord)
    allocate (self%ord(nat),source=0)
    if (allocated(self%grp)) deallocate (self%grp)
    allocate (self%grp(nat))
    return
  end subroutine allocate_zequal
!=========================================================================================!
  function is_x_member(self,x)
    implicit none
    logical :: is_x_member
    class(zequal) :: self
    integer :: x
    is_x_member = .false.
    if (self%ord(x) .ne. 0) then
      is_x_member = .true.
    end if
    return
  end function is_x_member
!=========================================================================================!
  subroutine zequal_geteng(self)
    implicit none
    class(zequal) :: self
    integer :: x,i
    x = 0
    if (allocated(self%grp)) then
      do i = 1,self%ng
        if (self%grp(i)%nm .gt. 1) x = x+1
      end do
    end if
    self%eng = x
    return
  end subroutine zequal_geteng
!=========================================================================================!
  subroutine prsummary_zequal(self,ch)
    implicit none
    class(zequal) :: self
    integer :: ch
    integer :: i,j
    write (ch,'(1x,a,i0)') 'Number of atoms in molecule: ',self%nat
    write (ch,'(1x,a,i0)') 'Number of groups :           ',self%ng
    if (allocated(self%ord)) then
      do i = 1,self%ng
        if (self%grp(i)%nm .gt. 1) then
          write (ch,'(3x,a,i0,a)',advance='no') 'group ',i,': '
          do j = 1,self%grp(i)%nm
            write (ch,'(i0,1x)',advance='no') self%grp(i)%mem(j)
          end do
          write (ch,*)
        end if
      end do
    end if
    return
  end subroutine prsummary_zequal
!=========================================================================================!
  subroutine allocate_zensemble(self,nat,nall)
    implicit none
    class(zensemble) :: self
    integer :: nat,nall
    self%nat = nat
    self%nall = nall
    if (allocated(self%at)) deallocate (self%at)
    if (allocated(self%xyz)) deallocate (self%xyz)
    if (allocated(self%eread)) deallocate (self%eread)
    allocate (self%at(nat),source=0)
    allocate (self%xyz(3,nat,nall),source=0.0_wp)
    allocate (self%eread(nall),source=0.0_wp)
    return
  end subroutine allocate_zensemble
!=========================================================================================!
  subroutine deallocate_zensemble(self)
    implicit none
    class(zensemble) :: self
    if (allocated(self%at)) deallocate (self%at)
    if (allocated(self%xyz)) deallocate (self%xyz)
    if (allocated(self%eread)) deallocate (self%eread)
    if (allocated(self%nrot)) deallocate (self%nrot)
    if (allocated(self%nmat)) deallocate (self%nmat)
    return
  end subroutine deallocate_zensemble
!=========================================================================================!
!>  WBO reader
!=========================================================================================!
  subroutine readwbo(fname,nat,wbo)
    implicit none
    character(len=*) :: fname
    character(len=256) :: atmp
    integer,intent(in) :: nat
    real(wp),intent(inout) :: wbo(nat,nat)
    logical :: ex
    integer :: i,j,ich,io
    wbo = 0 ! if it does not exist, single and double bonds are not distinguished
    inquire (file=fname,exist=ex)
    if (ex) then
      open (newunit=ich,file=fname)
      do
        read (ich,'(a)',iostat=io) atmp
        if (io < 0) exit
        read (atmp,*,iostat=io) i,j,wbo(i,j)
        if (io .ne. 0) cycle
        wbo(j,i) = wbo(i,j)
      end do
      close (ich)
    end if
    return
  end subroutine readwbo

!=========================================================================================!
!> MREC routine for the ZMOLECULE class, i.e., parse the neighbourlist recursively
!> and determine the total number of fragments in the structure
!=========================================================================================!
  subroutine zmol_mrec(self)
    implicit none
    class(zmolecule) :: self
    logical,allocatable :: taken(:)
    integer :: i
    if (allocated(self%molvec)) deallocate (self%molvec)
    allocate (self%molvec(self%nat))
    allocate (taken(self%nat))
    self%molvec = 0
    self%nfrag = 1
    taken = .false.
    do i = 1,self%nat
      if (.not.taken(i)) then
        self%molvec(i) = self%nfrag
        taken(i) = .true.
        call zmol_frag(i,self,self%nat,taken,self%nfrag,self%molvec)
        self%nfrag = self%nfrag+1
      end if
    end do
    self%nfrag = self%nfrag-1
    deallocate (taken)
    return
  end subroutine zmol_mrec
  recursive subroutine zmol_frag(i,self,nat,taken,nfrag,molvec)
    implicit none
    integer :: i,nat
    class(zmolecule) :: self
    logical :: taken(nat)
    integer :: nfrag
    integer :: molvec(nat)
    integer :: j,k
    do j = 1,self%zat(i)%nei
      k = self%zat(i)%ngh(j)
      if (.not.taken(k)) then
        molvec(k) = nfrag
        taken(k) = .true.
        call zmol_frag(k,self,nat,taken,nfrag,molvec)
      end if
    end do
    return
  end subroutine zmol_frag

  subroutine zmol_getxyz(self,xyz)
    implicit none
    class(zmolecule) :: self
    real(wp) :: xyz(3,self%nat)
    integer :: i
    do i = 1,self%nat
      xyz(1:3,i) = self%zat(i)%cart(1:3)
    end do
    return
  end subroutine zmol_getxyz

!=========================================================================================!
!> count number of hydrogen atoms
  function count_hydrogen(self) result(hnum)
    implicit none
    class(zmolecule) :: self
    integer :: hnum,i
    hnum = 0
    do i = 1,self%nat
      if (self%at(i) == 1) then
        hnum = hnum+1
      end if
    end do
    return
  end function count_hydrogen

!=========================================================================================!
!> count number of bonds in the molecule
  subroutine count_bonds(self)
    implicit none
    class(zmolecule) :: self
    integer :: i,j
    integer :: nb
    if (allocated(self%bondpairs)) deallocate (self%bondpairs)
    nb = 0
    do i = 1,self%nat
      do j = 1,i
        if (any(self%zat(i)%ngh(:) .eq. j)) then  !only include bonds from the neighbour lists
          nb = nb+1
        else
          cycle
        end if
      end do
    end do
    self%nb = nb
    allocate (self%bondpairs(2,nb),source=0)
    nb = 0
    do i = 1,self%nat
      do j = 1,i
        if (any(self%zat(i)%ngh(:) .eq. j)) then  !only include bonds from the neighbour lists
          nb = nb+1
          self%bondpairs(1,nb) = j
          self%bondpairs(2,nb) = i
        else
          cycle
        end if
      end do
    end do
    return

  end subroutine count_bonds

!=========================================================================================!
!>  get the adjacency and edgelength matrices from the neighbour lists
  subroutine zmol_adjacency(zmol,A,E)
    implicit none
    class(zmolecule) :: zmol
    integer,intent(out) :: A(zmol%nat,zmol%nat)
    real(wp),intent(out),optional :: E(zmol%nat,zmol%nat)
    integer :: i,j
    integer :: nat
    nat = zmol%nat
    A = 0
    if (present(E)) E = 0.0_wp
    do i = 1,nat
      do j = 1,nat
        !> only include bonds from the neighbour lists
        if (any(zmol%zat(i)%ngh(:) .eq. j)) then
          if (present(E)) E(i,j) = zmol%distmat(i,j)
          A(i,j) = 1
        else
          cycle
        end if
      end do
    end do
    return
  end subroutine zmol_adjacency

!=========================================================================================!
!=========================================================================================!
end module zdata
