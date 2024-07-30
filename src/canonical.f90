
module canonical_mod
!****************************************************************
!* Implementation of the CANGEN algorithm by Weininger et al.
!* D.Weininger et al., J. Chem. Inf. Comput. Sci., 1989, 29, 97-101.
!* doi.org/10.1021/ci00062a008
!*
!* The algorithm is used as the backend implementation to
!* determine the unqiue atom sequence in canonical SMILES.
!* It is pretty useful for all kinds of stuff, e.g.
!* for determining tautomers, or trying to restore atom order
!****************************************************************
  use crest_parameters
  use strucrd
  use adjacency
  implicit none
  private

  public :: canonical_sorter
  type :: canonical_sorter
    !> system mapping
    integer :: nat = 0
    integer :: hatms = 0
    integer,allocatable :: nmap(:)
    integer,allocatable :: hmap(:)
    integer,allocatable :: neigh(:,:) !> heavy-atom adjacency matrix

    !> the important bit for the algorithm
    integer(int64),allocatable :: invariants(:)
    integer,allocatable :: invariants0(:)
    integer,allocatable :: prime(:)
    integer,allocatable :: rank(:)        !> That's what we are after

    !> workspace helpers
    integer,allocatable :: newrank(:)
    integer(int64),allocatable :: newinv(:)

  contains
    procedure :: deallocate => deallocate_canonical_sorter
    procedure :: init => init_canonical_sorter
    procedure :: update_ranks
    procedure :: update_invariants
    procedure :: iterate
  end type canonical_sorter

  logical,parameter :: debug = .false.

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine deallocate_canonical_sorter(self)
    implicit none
    class(canonical_sorter),intent(inout) :: self
    self%nat = 0
    self%hatms = 0
    if (allocated(self%nmap)) deallocate (self%nmap)
    if (allocated(self%hmap)) deallocate (self%hmap)
    if (allocated(self%invariants)) deallocate (self%invariants)
    if (allocated(self%invariants0)) deallocate (self%invariants0)
    if (allocated(self%rank)) deallocate (self%rank)
    if (allocated(self%prime)) deallocate (self%prime)
    if (allocated(self%neigh)) deallocate (self%neigh)
    if (allocated(self%newrank)) deallocate (self%newrank)
    if (allocated(self%newinv)) deallocate (self%newinv) 
  end subroutine deallocate_canonical_sorter

!========================================================================================!

  subroutine init_canonical_sorter(self,mol,wbo)
    implicit none
    class(canonical_sorter),intent(inout) :: self
    type(coord),intent(in) :: mol
    real(wp),intent(in) :: wbo(mol%nat,mol%nat)
    integer :: nodes
    integer,allocatable :: Amat(:,:) !> adjacency matrix for FULL molecule
    integer :: counth,countb,countbo
    real(wp) :: countbo2
    integer :: i,j,k,l,ii,ati,atj

!>--- all atoms of the full mol. graph are nodes
    nodes = mol%nat

!>--- map to heavy atom-only representation
    k = 0
    do i = 1,mol%nat
      if (mol%at(i) .ne. 1) k = k+1
    end do
    self%nat = nodes
    self%hatms = k
    allocate (self%nmap(nodes))
    allocate (self%hmap(k))
    allocate (self%invariants(k),source=0_int64)
    allocate (self%invariants0(k),source=0)
    allocate (self%prime(k),source=2)
    allocate (self%rank(k),source=1)
    allocate (self%neigh(k,k),source=0)

!>--- get connectivity. Easiest is just via WBO
    call wbo2adjacency(nodes,wbo,Amat,0.02_wp)

    k = 0
    do i = 1,nodes
      if (mol%at(i) .ne. 1) then
        k = k+1
        self%nmap(i) = k
        self%hmap(k) = i
      else
        self%nmap(i) = 0
      end if
    end do
    do i = 1,k
      do j = 1,i-1
        self%neigh(j,i) = Amat(self%hmap(j),self%hmap(i))
        self%neigh(i,j) = self%neigh(j,i)
      end do
    end do

!>--- get the first invatiants
    do i = 1,k
      ii = self%hmap(i)
      ati = mol%at(ii)
      counth = 0
      countbo2 = 0.0_wp
      countb = 0
      do j = 1,nodes
        if (Amat(j,ii) .ne. 0) then
          if (mol%at(j) .eq. 1) counth = counth+1 !> count H neighbours
          countb = countb+1 !> count all neighbours
        end if
        countbo2 = countbo2+wbo(j,ii)  !> sum the total bond order
      end do
      countb = countb-counth          !> only heavy atom connections
      countbo = nint(countbo2)-counth !> same for number of bonds
      self%invariants(i) = get_invariant0(ati,countb,countbo,0,counth)
      self%invariants0(i) = self%invariants(i) !> back up for later tasks
    end do
   
    deallocate(Amat)

    if (debug) then
      call debugprint(self,mol)
    end if
!>--- start assignment
    allocate (self%newrank(k), source=0) !> workspace
    allocate (self%newinv(k), source=0_int64) !>workspace
    call self%update_ranks()
    self%rank(:) = self%newrank(:)
    if (debug) then
      call debugprint(self,mol)
    end if
    call self%update_invariants()
    if (debug) then
      call debugprint(self,mol)
    end if
    call self%iterate(mol) !> iterate recursively until ranking doesn't change
  end subroutine init_canonical_sorter

!========================================================================================!

  function get_invariant0(ati,nneigh,nbonds,chrg,hneigh) result(inv)
    implicit none
    integer :: inv
    integer,intent(in) :: ati !> atomic number
    integer,intent(in) :: nneigh !> # neighbours
    integer,intent(in) :: nbonds !> # bonds (sum of wbo)
    integer,intent(in) :: chrg   !> charge on atom
    integer,intent(in) :: hneigh !> # H atoms bound to atom
    inv = 0
    inv = inv+nneigh*10**5
    inv = inv+nbonds*10**4
    inv = inv+ati*10**2
    inv = inv+chrg*10**1
    inv = inv+hneigh
  end function get_invariant0

!=========================================================================================!

  subroutine update_ranks(self)
!>---update ranks and primes
    implicit none
    class(canonical_sorter) :: self
    integer :: maxrank,i,j,k,ii
    integer :: newrank,ngroup
    integer(int64) :: mincurr
    maxrank = maxval(self%rank,1)
    newrank = 0
    !>-- use newinv as workingspace
    self%newinv(:) = self%invariants(:)
    RANKLOOP: do i = 1,maxrank
      ngroup = count(self%rank(:) .eq. i)
      if (debug) then
        write (stdout,*) ngroup,' for rank ',i
      end if
      !>--- loop over all in the group
      GROUPLOOP: do ii = 1,ngroup
        mincurr = minval(self%newinv(:),self%rank(:) .eq. i) !> assign current minimum value within group
        if (mincurr == huge(mincurr)) cycle RANKLOOP !> cycle if all in group have been assigned
        newrank = newrank+1
        ASSIGNLOOP: do j = 1,self%hatms
          if (self%rank(j) == i.and.self%newinv(j) == mincurr) then
            self%newrank(j) = newrank
            if (debug) then
              write (stdout,*) 'new rank',newrank,'assigned to atom',j
            end if
            self%newinv(j) = huge(mincurr) !> remove this minimum value
          end if
        end do ASSIGNLOOP
      end do GROUPLOOP
    end do RANKLOOP
!>--- assign primes
    do i = 1,self%hatms
      !self%rank(i) = self%newrank(i) !> do this outside of routine
      self%prime(i) = nth_prime(self%newrank(i))
    end do
  end subroutine update_ranks

!=========================================================================================!

  subroutine update_invariants(self)
!>---update invariants
    implicit none
    class(canonical_sorter) :: self
    integer :: i,j,k,ii
    integer(int64) :: invprod
    do i = 1,self%hatms
      invprod = 1
      !>--- loop over all heavy-atom neighbours of i
      do j = 1,self%hatms
        if (self%neigh(j,i) > 0) then
          invprod = invprod*self%prime(j)
        end if
      end do
      self%newinv(i) = invprod
    end do
    if (debug) write (stdout,*) 'update invariants'
    self%invariants(:) = self%newinv(:)
  end subroutine update_invariants


!========================================================================================!

   recursive subroutine iterate(self,mol)
     implicit none
     class(canonical_sorter) :: self
     type(coord) :: mol
     call self%update_ranks()
     if (debug) then
      call debugprint(self,mol)
     end if
     if(all(self%rank(:).eq.self%newrank(:)))then
       return  !> termination condition
     else
       self%rank(:) = self%newrank(:)
     endif
     call self%update_invariants()
     if (debug) then
      call debugprint(self,mol)
     end if
     call self%iterate(mol)
   end subroutine iterate


!========================================================================================!

  function nth_prime(x) result(prime)
    implicit none
    integer,intent(in) :: x
    integer :: prime
    integer :: c,num,i
    logical :: is_prime
    integer,parameter :: prime_numbers(100) = (/2,3,5,7,11,13,17,19,23,29, &
    & 31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109, &
    & 113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197, &
    & 199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283, &
    & 293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389, &
    & 397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487, &
    & 491,499,503,509,521,523,541/)
    if (x <= 100) then
      prime = prime_numbers(x)
      return
    end if
    c = 0
    num = 1
    do while (c < x)
      num = num+1
      is_prime = .true.
      do i = 2,int(sqrt(real(num)))
        if (mod(num,i) == 0) then
          is_prime = .false.
          exit
        end if
      end do
      if (is_prime) then
        c = c+1
      end if
    end do
    prime = num
  end function nth_prime

!========================================================================================!

  subroutine debugprint(can,mol)
    implicit none
    type(canonical_sorter) :: can
    type(coord)  :: mol
    integer :: i,k,ii,ati
    write (stdout,'(a10,a5,a15,a10,a10)') 'heavy-atom','type','invariant','rank','prime'
    do i = 1,can%hatms
      ii = can%hmap(i)
      ati = mol%at(ii)
      write (stdout,'(i10,a5,i15,i10,i10)') i,i2e(ati,'nc'),can%invariants(i),can%rank(i),can%prime(i)
    end do
  end subroutine debugprint

!========================================================================================!
!========================================================================================!
end module canonical_mod
