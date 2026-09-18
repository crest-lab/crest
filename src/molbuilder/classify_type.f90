!=============================================================================!
! This file is part of crest.
!
! Copyright (C) 2026 Philipp Pracht
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
!=============================================================================!

module molbuilder_classify_type
  use crest_parameters,only:wp,stdout
  use strucrd,only:coord,i2e,sumform
  use adjacency
  use canonical_mod
  use molbuilder_rigidconf_analyze
  use INTERNALS_mod
  implicit none
  private

  type :: functional_group
    character(len=:),allocatable :: name
    integer :: natms = 0
    integer,allocatable :: ids(:)
    integer :: attached_to = 0
    logical :: seeded = .false.
  contains
    procedure :: clear => clear_func_group
    procedure :: copy => copy_func_group
  end type functional_group

  type :: mol_ring
    !> one ring (cycle) of the molecular graph. The member atom
    !> indices are stored sorted ascending so that two rings can be
    !> compared for identity by a simple element-wise comparison.
    integer :: size = 0             !> number of member atoms
    integer,allocatable :: atoms(:) !> member atom indices, sorted ascending
  contains
    procedure :: set    => mol_ring_set
    procedure :: equals => mol_ring_equals
  end type mol_ring

  type,private:: dihedral_types
    integer :: unknown = 0
    integer :: single = 1
    integer :: improper = 2
    integer :: stiff = 3
    integer :: macrocycle = 4
  end type dihedral_types
  type(dihedral_types),parameter,public :: dtypes = dihedral_types()

  type,extends(coord) :: coord_classify
    !> new components that are added to the coord type:
    !integer,allocatable :: A(:,:)  !> molecular graph/adjacency matrix
    integer,allocatable :: Ah(:,:) !> heavy-atom molecular graph/adjacency

    !> per-atom properties/information
    real(wp),allocatable :: CN(:)               !> coordination number
    integer,allocatable :: hyb(:)               !> hybridization/neighbours count
    integer,allocatable :: nhn(:)               !> non-H-neighbours count
    integer,allocatable :: prio(:)              !> "invariants"/atom priorities
    logical,allocatable :: inring(:)            !> atom part of ring?
    logical,allocatable :: term(:)              !> terminal atom (H,F,Cl,...,=O,etc.)
    character(len=10),allocatable :: atinfo(:)  !> atom info

    !> functional groups
    integer :: nfuncs = 0
    type(functional_group),allocatable :: funcgroups(:)

    !> covalent fragments (connected components of the molecular graph)
    integer :: nfrag = 0                !> number of disconnected fragments
    integer,allocatable :: fragment(:) !> per-atom fragment id (1..nfrag)

    !> ring library (unique cycles of the molecular graph)
    integer :: nrings = 0
    type(mol_ring),allocatable :: ringlist(:)

    !> internal coordinates
    integer :: ndieder = 0
    real(wp),allocatable :: zmat(:,:)
    integer,allocatable  :: zmap(:,:) !> na,nb,nc
    integer,allocatable  :: ztod(:)
    integer,allocatable  :: hatsort(:,:)
    integer,allocatable  :: dtype(:)

    !> utility storage
    logical,allocatable :: lwork(:)
    integer,allocatable :: iwork(:)

  contains
    procedure :: as_coord
    procedure :: from_coord
    generic,public :: add => coord_classify_add_fg
    procedure,private :: coord_classify_add_fg
    procedure :: get_zmat => coord_classify_calculate_zmat
    procedure :: from_zmat => coord_classify_reconstruct_from_zmat
    procedure :: update_zmat => coord_classify_update_zmat
    procedure :: check_dihedrals => coord_classify_check_dihedrals
    procedure :: collect_rings => coord_classify_collect_rings
    procedure :: print_funcgroups => coord_classify_print_functional
    procedure :: print_zmat => coord_classify_print_zmat
    procedure :: print_rings => coord_classify_print_rings
  end type coord_classify

  public :: coord_classify   !> the extended coord type
  public :: functional_group !> subtype of coord_classify
  public :: mol_ring         !> a single molecular-graph ring
  public :: setup_classify   !> setup a coord_classify from coord
  public :: atinfo_classify  !> add atinfo string to a coord_classify

!==============================================================================!
contains  !> MODULE PROCEDURES START HERE
!==============================================================================!

!> BASIC TYPE PROCEDURES

  function as_coord(this) result(mol)
    class(coord_classify),intent(in) :: this
    type(coord) :: mol

    mol%nat = this%nat
    if (allocated(this%at)) mol%at = this%at
    if (allocated(this%xyz)) mol%xyz = this%xyz

    mol%energy = this%energy
    if (allocated(this%comment)) mol%comment = this%comment
    mol%chrg = this%chrg
    mol%uhf = this%uhf
    mol%nbd = this%nbd
    if (allocated(this%bond)) mol%bond = this%bond
    if (allocated(this%lat)) mol%lat = this%lat
    if (allocated(this%qat)) mol%qat = this%qat
    mol%pdb = this%pdb

  end function as_coord

  subroutine from_coord(this,mol)
    class(coord_classify),intent(inout) :: this
    type(coord),intent(in) :: mol

    this%nat = mol%nat
    if (allocated(mol%at)) this%at = mol%at
    if (allocated(mol%xyz)) this%xyz = mol%xyz

    this%energy = mol%energy
    if (allocated(mol%comment)) this%comment = mol%comment
    this%chrg = mol%chrg
    this%uhf = mol%uhf
    this%nbd = mol%nbd
    if (allocated(mol%bond)) this%bond = mol%bond
    if (allocated(mol%lat)) this%lat = mol%lat
    if (allocated(mol%qat)) this%qat = mol%qat
    this%pdb = mol%pdb
  end subroutine from_coord

  subroutine clear_func_group(self)
    implicit none
    class(functional_group) :: self
    if (allocated(self%name)) deallocate (self%name)
    if (allocated(self%ids)) deallocate (self%ids)
    self%attached_to = 0
    self%natms = 0
    self%seeded = .false.
  end subroutine clear_func_group

  subroutine copy_func_group(self,fg)
    implicit none
    class(functional_group) :: self
    type(functional_group) :: fg
    if (allocated(fg%name)) self%name = fg%name
    if (allocated(fg%ids)) self%ids = fg%ids
    self%attached_to = fg%attached_to
    self%natms = fg%natms
    self%seeded = fg%seeded
  end subroutine copy_func_group

  subroutine coord_classify_add_fg(self,fg)
    implicit none
    class(coord_classify) :: self
    type(functional_group) :: fg
    type(functional_group),allocatable :: fg_list(:)
    integer :: ii,jj
    if (.not.allocated(self%funcgroups)) then
      allocate (self%funcgroups(1))
      call self%funcgroups(1)%copy(fg)
    else
      ii = size(self%funcgroups,1)
      allocate (fg_list(ii+1))
      do jj = 1,ii
        call fg_list(jj)%copy(self%funcgroups(jj))
      end do
      call fg_list(ii+1)%copy(fg)
      call move_alloc(fg_list,self%funcgroups)
    end if
    self%nfuncs = size(self%funcgroups,1)
  end subroutine coord_classify_add_fg

!=============================================================================!
!#############################################################################!
!=============================================================================!

!> CLASSIFICATION ROUTINES

  subroutine setup_classify(mol,molc,wbo)
    !***************************************************
    !* set up the derived coord_classify object "molc"
    !* from a standard coord object "mol".
    !* in particular, adjacency graphs, CN, etc.
    !***************************************************
    implicit none
    type(coord),intent(in) :: mol
    type(coord_classify),intent(out) :: molc
    real(wp),intent(in),optional :: wbo(:,:)

    real(wp),allocatable :: Bmat(:,:)
    logical,allocatable :: rings(:,:)
    type(canonical_sorter),allocatable :: tmpcan
    integer :: nat
    integer :: ii,jj

    !> Initialize
    call molc%from_coord(mol)
    nat = molc%nat

    !> set up CN, and from that topology
    call mol%cn_to_bond(molc%CN,Bmat,'cov')
    if (present(wbo)) then
      Bmat(:,:) = wbo(:,:)
    end if
    call wbo2adjacency(molc%nat,Bmat,molc%bond,0.02_wp)
    deallocate (Bmat)

    !> set up other parameters
    allocate (molc%hyb(nat),source=0)
    allocate (molc%inring(nat),source=.false.)
    allocate (molc%term(nat),source=.false.)
    call check_rings_min(nat,molc%bond,rings)
    do ii = 1,nat
      molc%hyb(ii) = sum(molc%bond(:,ii))
      if (any(rings(:,ii))) molc%inring(ii) = .true.
      if (molc%hyb(ii) .eq. 1) molc%term(ii) = .true.
    end do

    allocate (tmpcan)
    call tmpcan%init(mol,invtype='apsp+',heavy=.false.)
    call move_alloc(tmpcan%rank,molc%prio)
    call move_alloc(tmpcan%hadjac,molc%Ah)
    deallocate (tmpcan)

    allocate (molc%nhn(nat),source=0)
    do ii = 1,nat
      molc%nhn(ii) = sum(molc%Ah(:,ii))
    end do

    !> label covalent fragments (connected components of the graph)
    call coord_classify_fragments(molc)

    !> collect the unique rings of the molecular graph
    call molc%collect_rings()

  end subroutine setup_classify

  subroutine coord_classify_fragments(molc)
    !***************************************************************
    !* Label the covalent fragments (connected components) of the
    !* molecular graph molc%bond, filling molc%fragment(1:nat) with
    !* a fragment id per atom and molc%nfrag with the component
    !* count. An iterative flood fill over the adjacency matrix;
    !* nfrag = 1 for a single connected molecule.
    !***************************************************************
    implicit none
    type(coord_classify),intent(inout) :: molc
    integer :: nat,i,j,head,tail,a
    integer,allocatable :: stack(:)

    nat = molc%nat
    molc%nfrag = 0
    if (allocated(molc%fragment)) deallocate (molc%fragment)
    if (nat < 1.or..not.allocated(molc%bond)) return
    allocate (molc%fragment(nat),source=0)
    allocate (stack(nat),source=0)

    do i = 1,nat
      if (molc%fragment(i) /= 0) cycle           !> already in a fragment
      molc%nfrag = molc%nfrag+1
! ── flood fill the component reachable from atom i ───────────────────────────
      head = 1; tail = 1; stack(1) = i
      molc%fragment(i) = molc%nfrag
      do while (head <= tail)
        a = stack(head); head = head+1
        do j = 1,nat
          if (molc%bond(j,a) > 0.and.molc%fragment(j) == 0) then
            molc%fragment(j) = molc%nfrag
            tail = tail+1; stack(tail) = j
          end if
        end do
      end do
    end do
    deallocate (stack)
  end subroutine coord_classify_fragments

  subroutine atinfo_classify(molc)
    !*****************************************
    !* Update a coord_classify object "molc"
    !* and fill in its atinfo strings based
    !* on some basic chemoinformatics.
    !****************************************
    implicit none
    type(coord_classify),intent(inout) :: molc
    integer :: nat
    integer :: ii,jj

    if (molc%nat <= 0) then
      write (stdout,*) 'molc not allocated in atinfo_classify()'
      return
    end if
    if (allocated(molc%atinfo)) deallocate (molc%atinfo)
    nat = molc%nat
    allocate (molc%atinfo(nat),source=repeat(' ',10))

    do ii = 1,molc%nat
      associate (str => molc%atinfo(ii),ati => molc%at(ii))
        str = trim(i2e(ati))
        select case (ati)

        case (6) !> carbon
          if (molc%hyb(ii) == 3) then !> sp2
            str = trim(i2e(ati,'lc'))
          else if (molc%hyb(ii) == 4) then !> sp3
            if (molc%nhn(ii) == 1) then
              str = trim(str)//'H3'
            else if (molc%nhn(ii) == 2) then
              str = trim(str)//'H2'
            else if (molc%nhn(ii) == 0) then
              str = 'methane'
            end if
          end if

        case (7) !> nitrogen
          if (molc%hyb(ii) == 3) then
            if (molc%nhn(ii) == 1) then
              str = trim(str)//'H2'
            else if (molc%nhn(ii) == 2) then
              str = trim(str)//'R2'
            else if (molc%nhn(ii) == 3) then
              str = trim(str)//'R3'
            else if (molc%nhn(ii) == 0) then
              str = 'ammonia'
            end if
          else if (molc%hyb(ii) == 4) then
            str = trim(str)//'4+'
          end if

        case (8) !> oxygen
          if (molc%hyb(ii) == 1) then
            str = 'o'
          else if (molc%hyb(ii) == 2) then
            if (molc%nhn(ii) == 1) then
              str = trim(str)//'H'
            else if (molc%nhn(ii) == 0) then
              str = 'water'
            end if
          end if

        case (16) !> sulfur
          if (molc%hyb(ii) == 2) then
            if (molc%nhn(ii) == 1) then
              str = trim(str)//'H'
            end if
          end if

        end select
      end associate
    end do
  end subroutine atinfo_classify

  subroutine coord_classify_calculate_zmat(molc,natural)
    implicit none
    class(coord_classify),intent(inout) :: molc
    logical,intent(in),optional :: natural

    if (.not.allocated(molc%xyz)) return

    if (allocated(molc%zmat)) deallocate (molc%zmat)
    if (allocated(molc%zmap)) deallocate (molc%zmap)
    if (allocated(molc%ztod)) deallocate (molc%ztod)

    if (present(natural)) then
      if (natural) then
        !write (stdout,'(/,a)') 'NOTE: atom order will temporarily be changed!'
        call coord_classify_hatsort(molc)
      end if
    end if

    allocate (molc%zmap(molc%nat,3),source=0)
    allocate (molc%zmat(3,molc%nat),source=0.0_wp)
    call BETTER_XYZINT(molc%nat,molc%xyz,molc%bond, &
    &    molc%zmap(:,1),molc%zmap(:,2),molc%zmap(:,3),molc%zmat)

    allocate (molc%ztod(molc%nat),source=0)
    call rigidconf_count_fallback(molc%nat, &
   & molc%zmap(:,1),molc%zmap(:,2),molc%zmap(:,3), &
   & molc%bond,molc%ndieder,molc%ztod)
    if (present(natural)) then
      if (natural) then
        call prune_zmat_dihedrals(molc,molc%zmat, &
        & molc%zmap(:,1),molc%zmap(:,2),molc%zmap(:,3),molc%ztod, &
        hpyrad=.true.,bond=molc%bond)
        !call molc%print_zmat(stdout)
        call coord_classify_hatsort_restore(molc)
        deallocate (molc%hatsort)
      end if
    end if

  end subroutine coord_classify_calculate_zmat

  subroutine coord_classify_reconstruct_from_zmat(molc,mol)
    implicit none
    class(coord_classify),intent(inout) :: molc
    type(coord),intent(out),optional :: mol

    if (.not.allocated(molc%zmat)) then
      write (stdout,*) '** ERROR ** in coord_classify_reconstruct_from_zmat(): zmat not allocated!'
      return
    end if
    call GMETRY2(molc%nat,molc%zmat, &
      &              molc%xyz,        &
      &  molc%zmap(:,1),molc%zmap(:,2),molc%zmap(:,3))

    if (present(mol)) then
      mol = molc%as_coord()
    end if
  end subroutine coord_classify_reconstruct_from_zmat

  subroutine coord_classify_update_zmat(molc,mol)
    !************************************************************
    !* Update the Z-matrix with fresh coords from Cartesian ones
    !* (The mapping must exist at this point)
    !************************************************************
    implicit none
    class(coord_classify),intent(inout) :: molc
    type(coord),intent(in),optional :: mol
    integer :: ii,jj,a,b,c,d

    if (.not.allocated(molc%zmat)) then
      write (stdout,*) '** ERROR ** in coord_classify_update_zmat(): zmat not allocated!'
      return
    end if
    if (.not.allocated(molc%zmap)) then
      write (stdout,*) '** ERROR ** in coord_classify_update_zmat(): zmapping not allocated!'
      return
    end if
    if (present(mol)) then
      if (.not.all(mol%at .eq. molc%at)) then
        write (stdout,*) '** ERROR ** in coord_classify_update_zmat(): mismatch in atom order'
        return
      end if
      molc%xyz = mol%xyz
    end if
    do a = 1,molc%nat
      b = molc%zmap(a,1)
      if (b == 0) cycle
      molc%zmat(1,a) = molc%dist(a,b)
      c = molc%zmap(a,2)
      if (c == 0) cycle
      molc%zmat(2,a) = molc%angle(a,b,c)
      d = molc%zmap(a,3)
      if (d == 0) cycle
      molc%zmat(3,a) = molc%dihedral(a,b,c,d)
    end do
  end subroutine coord_classify_update_zmat

  subroutine coord_classify_check_dihedrals(molc)
    !************************************************************
    !* Attempt to assign dihedral angles to a type of dihedral
    !************************************************************
    implicit none
    class(coord_classify),intent(inout) :: molc
    integer :: ii,jj,a,b,c,d

    if (.not.allocated(molc%zmap)) then
      write (stdout,*) '** ERROR ** in coord_classify_update_zmat(): zmapping not allocated!'
      return
    end if

    if (allocated(molc%dtype)) deallocate (molc%dtype)
    allocate (molc%dtype(molc%nat),source=dtypes%unknown)

    do ii = 1,molc%nat
      if (molc%zmap(3,ii) .eq. 0) cycle
      a = ii
      b = molc%zmap(ii,1)
      c = molc%zmap(ii,2)
      d = molc%zmap(ii,3)

      if (molc%bond(a,b) > 0.and. &
      &  molc%bond(b,c) > 0.and. &
      &  molc%bond(c,d) > 0) then
        molc%dtype(ii) = dtypes%single
      else if(molc%bond(a,b) > 0.and. & 
      &       molc%bond(b,c) > 0.and. &  
      &       molc%bond(b,d) > 0) then   
        molc%dtype(ii) = dtypes%improper
      end if
      !write(*,*) ii, molc%dtype(ii)
    end do

  end subroutine coord_classify_check_dihedrals

!=============================================================================!

  subroutine coord_classify_hatsort(molc)
    !**************************************************************
    !* a routine that resorts the atomorder in molc so that
    !* hydrogen atoms come last. required for natural z-mat setup
    !* Also mapps the order to restore it later on
    !**************************************************************
    implicit none
    class(coord_classify),intent(inout) :: molc

    real(wp),allocatable :: xyztmp(:,:)
    integer,allocatable :: attmp(:),bondtmp(:,:)

    integer :: ii,kk,jj
    if (allocated(molc%hatsort)) deallocate (molc%hatsort)
    allocate (molc%hatsort(molc%nat,2),source=0)
    allocate (attmp(molc%nat),source=0)
    allocate (xyztmp(3,molc%nat),source=0.0_wp)

    kk = 0
    !> heavy atoms
    do ii = 1,molc%nat
      if (molc%at(ii) .ne. 1) then
        kk = kk+1
        molc%hatsort(kk,1) = ii
        molc%hatsort(ii,2) = kk
        xyztmp(1:3,kk) = molc%xyz(1:3,ii)
        attmp(kk) = molc%at(ii)
      end if
    end do
    !> hydrogen
    do ii = 1,molc%nat
      if (molc%at(ii) .eq. 1) then
        kk = kk+1
        molc%hatsort(kk,1) = ii
        molc%hatsort(ii,2) = kk
        xyztmp(1:3,kk) = molc%xyz(1:3,ii)
        attmp(kk) = molc%at(ii)
      end if
    end do

    call move_alloc(xyztmp,molc%xyz)
    call move_alloc(attmp,molc%at)

    if (allocated(molc%bond)) then
      allocate (bondtmp(molc%nat,molc%nat),source=0)
      do ii = 1,molc%nat
        do jj = 1,molc%nat
          bondtmp(molc%hatsort(jj,2),molc%hatsort(ii,2)) = molc%bond(jj,ii)
        end do
      end do
      call move_alloc(bondtmp,molc%bond)
    end if
  end subroutine coord_classify_hatsort

  subroutine coord_classify_hatsort_restore(molc)
    !**********************************************
    !* Restore original order from h-atom sorting
    !*********************************************
    implicit none
    class(coord_classify),intent(inout) :: molc

    real(wp),allocatable :: xyztmp(:,:)
    integer,allocatable :: attmp(:),bondtmp(:,:)
    integer,allocatable :: ztodtmp(:),zmaptmp(:,:)
    real(wp),allocatable :: zmattmp(:,:)
    integer :: ii,kk,jj
    if (.not.allocated(molc%hatsort)) return

    allocate (attmp(molc%nat),source=0)
    allocate (xyztmp(3,molc%nat),source=0.0_wp)
    do ii = 1,molc%nat
      kk = molc%hatsort(ii,1)
      xyztmp(1:3,kk) = molc%xyz(1:3,ii)
      attmp(kk) = molc%at(ii)
    end do
    call move_alloc(xyztmp,molc%xyz)
    call move_alloc(attmp,molc%at)
    if (allocated(molc%bond)) then
      allocate (bondtmp(molc%nat,molc%nat),source=0)
      do ii = 1,molc%nat
        do jj = 1,molc%nat
          bondtmp(molc%hatsort(jj,1),molc%hatsort(ii,1)) = molc%bond(jj,ii)
        end do
      end do
      call move_alloc(bondtmp,molc%bond)
    end if

    allocate (ztodtmp(molc%nat),source=0)
    allocate (zmaptmp(molc%nat,3),source=0)
    allocate (zmattmp(3,molc%nat),source=0.0_wp)
    do ii = 1,molc%nat
      ztodtmp(molc%hatsort(ii,1)) = molc%ztod(ii)
      zmattmp(1:3,molc%hatsort(ii,1)) = molc%zmat(1:3,ii)
      do jj = 1,3
        if (molc%zmap(ii,jj) > 0) then
          zmaptmp(molc%hatsort(ii,1),jj) = molc%hatsort(molc%zmap(ii,jj),1)
        else
          zmaptmp(molc%hatsort(ii,1),jj) = 0
        end if
      end do
    end do
    call move_alloc(ztodtmp,molc%ztod)
    call move_alloc(zmaptmp,molc%zmap)
    call move_alloc(zmattmp,molc%zmat)
  end subroutine coord_classify_hatsort_restore

!=============================================================================!
!#############################################################################!
!=============================================================================!

!> PRINTOUTS and naming

  subroutine coord_classify_print_functional(self,prch)
    implicit none
    class(coord_classify) :: self
    integer,intent(in) :: prch

    integer,allocatable :: at(:)
    integer :: ii,jj,nn

    if (.not.allocated(self%funcgroups)) return

    do ii = 1,size(self%funcgroups,1)
      nn = self%funcgroups(ii)%natms
      allocate (at(nn),source=0)
      do jj = 1,nn
        at(jj) = self%at(self%funcgroups(ii)%ids(jj))
      end do
      write (prch,'(3(1x,a))') 'functional group:', &
        & self%funcgroups(ii)%name,sumform(nn,at)
      deallocate (at)
    end do
  end subroutine coord_classify_print_functional

  subroutine coord_classify_print_zmat(self,prch)
    implicit none
    class(coord_classify) :: self
    integer,intent(in) :: prch
    if (.not.allocated(self%zmat)) then
      write (prch,*) 'zmat not allocated!'
      return
    end if

    write (prch,'(/,a)') 'Internal coordinates:'
    call print_zmat(prch,self%nat,self%at,self%zmat, &
    &    self%zmap(:,1),self%zmap(:,2),self%zmap(:,3),.true.)
  end subroutine coord_classify_print_zmat

!=============================================================================!
!#############################################################################!
!=============================================================================!

!> RING-LIBRARY PROCEDURES

  subroutine mol_ring_set(self,members)
    !************************************************************
    !* Store a ring from its member atom indices, sorting them
    !* ascending so that rings can be compared element-wise.
    !************************************************************
    implicit none
    class(mol_ring),intent(inout) :: self
    integer,intent(in) :: members(:)
    integer :: n,i,j,key
    n = size(members)
    if (allocated(self%atoms)) deallocate (self%atoms)
    allocate (self%atoms(n))
    self%atoms(:) = members(:)
    self%size = n
    !> insertion sort (rings are small, so this is plenty fast)
    do i = 2,n
      key = self%atoms(i)
      j = i-1
      do while (j >= 1)
        if (self%atoms(j) <= key) exit
        self%atoms(j+1) = self%atoms(j)
        j = j-1
      end do
      self%atoms(j+1) = key
    end do
  end subroutine mol_ring_set

  logical function mol_ring_equals(self,other) result(eq)
    !************************************************************
    !* Two rings are identical iff they have the same size and
    !* the same (sorted) member atoms.
    !************************************************************
    implicit none
    class(mol_ring),intent(in) :: self
    type(mol_ring),intent(in) :: other
    integer :: i
    eq = .false.
    if (self%size /= other%size) return
    if (.not.allocated(self%atoms).or..not.allocated(other%atoms)) return
    do i = 1,self%size
      if (self%atoms(i) /= other%atoms(i)) return
    end do
    eq = .true.
  end function mol_ring_equals

  subroutine ring_add_unique(self,newring)
    !************************************************************
    !* Append "newring" to the molecule's ring library unless an
    !* identical ring is already stored (duplicate suppression).
    !************************************************************
    implicit none
    type(coord_classify),intent(inout) :: self
    type(mol_ring),intent(in) :: newring
    type(mol_ring),allocatable :: tmp(:)
    integer :: k
    do k = 1,self%nrings
      if (self%ringlist(k)%equals(newring)) return   !> already known
    end do
    if (.not.allocated(self%ringlist)) then
      allocate (self%ringlist(1))
      self%ringlist(1) = newring
    else
      allocate (tmp(self%nrings+1))
      tmp(1:self%nrings) = self%ringlist(1:self%nrings)
      tmp(self%nrings+1) = newring
      call move_alloc(tmp,self%ringlist)
    end if
    self%nrings = self%nrings+1
  end subroutine ring_add_unique

  subroutine coord_classify_collect_rings(self)
    !************************************************************
    !* Build the molecule's ring library from its adjacency
    !* graph. check_rings_min flags every directly bonded vertex
    !* pair (M,N) that still shares a path once their bond is
    !* removed (i.e. lies on a ring); get_ring_min then returns
    !* the smallest such ring. Iterating over all ring-bearing
    !* bonds and discarding duplicates yields the set of unique
    !* smallest rings (a smallest-set-of-smallest-rings flavour).
    !************************************************************
    implicit none
    class(coord_classify),intent(inout) :: self
    logical,allocatable :: ringmask(:,:)
    integer,allocatable :: path(:)
    integer :: nat,i,j,nring
    type(mol_ring) :: newring

    self%nrings = 0
    if (allocated(self%ringlist)) deallocate (self%ringlist)
    nat = self%nat
    if (nat < 3.or..not.allocated(self%bond)) return

    call check_rings_min(nat,self%bond,ringmask)
    allocate (path(nat),source=0)
    do i = 1,nat
      do j = 1,i-1
        if (.not.ringmask(i,j)) cycle           !> bond (i,j) not on a ring
        call get_ring_min(nat,self%bond,i,j,path,nring)
        if (nring < 3) cycle                     !> not a valid ring
        call newring%set(path(1:nring))
        call ring_add_unique(self,newring)
      end do
    end do
    deallocate (path)
    if (allocated(ringmask)) deallocate (ringmask)
  end subroutine coord_classify_collect_rings

  subroutine coord_classify_print_rings(self,prch)
    !************************************************************
    !* Print a short summary of the detected ring library.
    !************************************************************
    implicit none
    class(coord_classify) :: self
    integer,intent(in) :: prch
    integer :: ii,jj,a
    character(len=:),allocatable :: line
    character(len=16) :: tok

    if (self%nrings < 1) then
      write (prch,'(/,1x,a)') 'Ring perception: no rings detected.'
      return
    end if
    write (prch,'(/,1x,a,i0,a)') 'Ring perception: ',self%nrings, &
    &  ' unique ring(s) detected'
    do ii = 1,self%nrings
      line = ''
      do jj = 1,self%ringlist(ii)%size
        a = self%ringlist(ii)%atoms(jj)
        write (tok,'(a,i0)') trim(i2e(self%at(a))),a
        line = trim(line)//' '//trim(tok)
      end do
      write (prch,'(3x,a,i0,a,i0,a,a)') 'ring ',ii,'  (',self%ringlist(ii)%size, &
      &  '-membered):',trim(line)
    end do
  end subroutine coord_classify_print_rings

!=============================================================================!
!#############################################################################!
!=============================================================================!
end module molbuilder_classify_type

