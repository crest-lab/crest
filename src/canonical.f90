
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
  use geo
  implicit none
  private

  public :: canonical_sorter
  type :: canonical_sorter
    !> system mapping
    integer :: nat = 0    !> total number of atoms
    integer :: hatms = 0  !> number of heavy atoms (>H)
    integer :: nfrag = 0  !> number of molecules/subgraphs in structure
    integer,allocatable :: nmap(:)  !> map atom to heavy-atom list (H->0)
    integer,allocatable :: hmap(:)  !> map heavy-atom to full atom order
    integer :: maxnei = 0
    integer,allocatable :: neigh(:,:)  !> neighbour list (neigh(j,i) = j-th neighbour of atom i)
    integer,allocatable :: hadjac(:,:) !> heavy-atom adjacency matrix

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
    procedure :: shrink => shrink_canonical_sorter
    procedure :: init => init_canonical_sorter
    procedure :: update_ranks
    procedure :: update_invariants
    procedure :: iterate
    procedure :: rankprint
    procedure :: stereo => analyze_stereo
    procedure :: compare => compare_canonical_sorter
  end type canonical_sorter

  logical,parameter :: debug = .false.

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine deallocate_canonical_sorter(self)
!********************************************
!* reset canonical_sorter entirely
!********************************************
    implicit none
    class(canonical_sorter),intent(inout) :: self
    self%nat = 0
    self%hatms = 0
    self%nfrag = 0
    if (allocated(self%nmap)) deallocate (self%nmap)
    if (allocated(self%hmap)) deallocate (self%hmap)
    if (allocated(self%invariants)) deallocate (self%invariants)
    if (allocated(self%invariants0)) deallocate (self%invariants0)
    if (allocated(self%rank)) deallocate (self%rank)
    if (allocated(self%prime)) deallocate (self%prime)
    if (allocated(self%hadjac)) deallocate (self%hadjac)
    if (allocated(self%newrank)) deallocate (self%newrank)
    if (allocated(self%newinv)) deallocate (self%newinv)
    if (allocated(self%neigh)) deallocate(self%neigh)
  end subroutine deallocate_canonical_sorter

!========================================================================================!

  subroutine shrink_canonical_sorter(self)
!***************************************************************
!* Reduce memory consumption of canconical_sorter
!* by deallocating everything except the (original) invariants
!* mappings and the determined ranks
!***************************************************************
    implicit none
    class(canonical_sorter),intent(inout) :: self
    if (allocated(self%invariants)) deallocate (self%invariants)
    if (allocated(self%prime)) deallocate (self%prime)
    if (allocated(self%hadjac)) deallocate (self%hadjac)
    if (allocated(self%newrank)) deallocate (self%newrank)
    if (allocated(self%newinv)) deallocate (self%newinv)
    if (allocated(self%neigh)) deallocate(self%neigh)
  end subroutine shrink_canonical_sorter

!========================================================================================!

  subroutine init_canonical_sorter(self,mol,wbo)
!*****************************************************************
!* Initializes the canonical_sorter and runs the CANGEN algorithm
!*****************************************************************
    implicit none
    class(canonical_sorter),intent(inout) :: self
    type(coord),intent(in) :: mol
    real(wp),intent(in) :: wbo(mol%nat,mol%nat)
    integer :: nodes
    integer,allocatable :: Amat(:,:) !> adjacency matrix for FULL molecule
    integer :: counth,countb,countbo
    real(wp) :: countbo2
    integer :: i,j,k,l,ii,ati,atj,maxnei
    integer,allocatable :: ichrgs(:),frag(:)
    logical :: use_icharges

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
    allocate (self%hadjac(k,k),source=0)

!>--- get connectivity. Easiest is just via WBO (allocates Amat)
    call wbo2adjacency(nodes,wbo,Amat,0.02_wp)

!>--- determine number of subgraphs
    allocate(frag(nodes),source=0)
    call setup_fragments(nodes,Amat,frag)
    self%nfrag=maxval(frag(:),1)

!>--- documment neighbour list
    maxnei = 0
    do i = 1,mol%nat
      k = count(Amat(:,i) > 0)
      if (k > maxnei) maxnei = k
    end do
    if (debug) write (stdout,*) 'maximum number of neighbours',maxnei
    self%maxnei = maxnei
    allocate (self%neigh(maxnei,mol%nat),source=0)

!>--- fill rest of self
    k = 0
    do i = 1,nodes
      l = 0
      if (mol%at(i) .ne. 1) then
        k = k+1
        self%nmap(i) = k
        self%hmap(k) = i
      else
        self%nmap(i) = 0
      end if
      do j = 1,nodes
        if (Amat(j,i) > 0) then
          l = l+1
          self%neigh(l,i) = j
        end if
      end do
    end do
    do i = 1,k
      do j = 1,i-1
        self%hadjac(j,i) = Amat(self%hmap(j),self%hmap(i))
        self%hadjac(i,j) = self%hadjac(j,i)
      end do
    end do

!>--- get the first invatiants
    if(allocated(mol%qat))then
       use_icharges = .true.
       allocate(ichrgs(mol%nat), source=0)
       ichrgs(:) = nint(mol%qat(:))
    else
       use_icharges = .false.
    endif
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
      if(use_icharges)then
        self%invariants(i) = get_invariant0(ati,countb,countbo,ichrgs(ii),counth) 
      else
        self%invariants(i) = get_invariant0(ati,countb,countbo,0,counth)
      endif
      self%invariants0(i) = int(self%invariants(i)) !> back up for later tasks
    end do

    deallocate (Amat)

    if (debug) then
      call debugprint(self,mol)
    end if
!>--- start assignment
    allocate (self%newrank(k),source=0) !> workspace
    allocate (self%newinv(k),source=0_int64) !>workspace
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
        if (self%hadjac(j,i) > 0) then
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
    if (all(self%rank(:) .eq. self%newrank(:))) then
      return  !> termination condition
    else
      self%rank(:) = self%newrank(:)
    end if
    call self%update_invariants()
    if (debug) then
      call debugprint(self,mol)
    end if
    call self%iterate(mol)
  end subroutine iterate

!========================================================================================!

  subroutine analyze_stereo(self,mol)
    implicit none
    class(canonical_sorter) :: self
    type(coord),intent(in) :: mol
    integer :: i,ii,zero,nei,j,jj,maxrank
    integer :: k,l,rs
    integer,allocatable :: neiranks(:,:)
    real(wp) :: coords(3,4)
    logical,allocatable :: isstereo(:)
    allocate (isstereo(mol%nat),source=.false.)
    allocate (neiranks(4,mol%nat),source=0)
    maxrank = maxval(self%rank(:))
    do i = 1,self%hatms
      ii = self%hmap(i)
      zero = count(self%neigh(:,ii) == 0)
      nei = self%maxnei-zero
!>--- consider only atoms with 4 unique (in terms of CANGEN ranks) neighbours as stereocenter
      if (nei == 4) then
        do j = 1,4
          jj = self%neigh(j,ii)
          if (mol%at(jj) == 1) then !> one hydrogen allowed
            neiranks(j,ii) = maxrank
          else
            neiranks(j,ii) = self%rank(jj)
          end if
        end do
        !isstereo(ii) = unique_neighbours(nei,self%neigh(1:nei,ii))
        isstereo(ii) = unique_neighbours(4,neiranks(:,ii))
      end if
    end do
!>--- do some actual geometry checks to determine "R" and "S"
    do i = 1,mol%nat
      if (isstereo(i)) then
        !>--- transfer to dummy coords
        do j = 1,4
          k = minloc(neiranks(:,i),1)
          jj = self%neigh(k,i)
          neiranks(k,i) = huge(k)
          coords(1:3,j) = mol%xyz(1:3,jj)-mol%xyz(1:3,i) !> shift i to center of coordinates
        end do
        RS = determineRS(coords)
        ii = self%nmap(i)
        self%invariants0(ii) = (self%invariants0(ii)+abs(RS)*10**6)*RS
      end if
    end do

    if (debug) then
      do ii = 1,self%hatms
        self%invariants(ii) = self%invariants0(ii)
      end do
      call debugprint(self,mol)
    end if
    deallocate (neiranks,isstereo)
  end subroutine analyze_stereo

!========================================================================================!

  function compare_canonical_sorter(self, other) result(yesno)
!*****************************************************
!* compare two canonical_sorter objects to determine
!* if both correspond to the same molecule
!*****************************************************
    implicit none
    class(canonical_sorter) :: self
    type(canonical_sorter) :: other
    logical :: yesno
    integer,allocatable :: sorted_invariants(:,:)
    integer :: maxrank_a,maxrank_b,hatms
    integer :: i,j,jj,k,kk

    yesno=.false.
!>--- the obvious cases first
    if(self%nat .ne. other%nat) return
    if(self%hatms .ne. other%hatms) return  
    maxrank_a = maxval(self%rank(:),1)
    maxrank_b = maxval(other%rank(:),1)     
    if(maxrank_a .ne. maxrank_b) return

!>--- if the easy checks passed, compare the actual invariants!
    hatms=self%hatms
    allocate(sorted_invariants(hatms,2), source=0)
    jj=0
    kk=0
    do i=1,maxrank_a  !> maxrank_a == maxrank_b, see above
      do j=1,hatms
       if(self%rank(j)==i)then
          jj=jj+1
           sorted_invariants(jj,1) = self%invariants0(j)
       endif 
      enddo
      do k=1,hatms
       if(other%rank(k)==i)then
          kk=kk+1
          sorted_invariants(kk,2) = other%invariants0(k)
       endif
      enddo
    enddo

    if(debug)then
       do i=1,hatms
          write(stdout,'(i10,i10,i10)') i,sorted_invariants(i,:)
       enddo
    endif    

!>--- assign identity
    yesno = all(sorted_invariants(:,1).eq.sorted_invariants(:,2))

    deallocate(sorted_invariants)
    return
  end function compare_canonical_sorter


!========================================================================================!
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

  subroutine rankprint(can,mol)
    implicit none
    class(canonical_sorter) :: can
    type(coord)  :: mol
    integer :: i,k,ii,ati
    write (stdout,'(a10,a10,a10,a10,2x,a)') 'heavy-atom','type','invariant0','rank','neighbours'
    do i = 1,can%hatms
      ii = can%hmap(i)
      ati = mol%at(ii)
      write (stdout,'(i10,a10,i10,i10,2x,a)') i,i2e(ati,'nc'),can%invariants0(i),can%rank(i),print_neighbours(mol,can%neigh(:,ii))
    end do
  end subroutine rankprint

  function print_neighbours(mol,neigh) result(btmp)
    implicit none
    type(coord) :: mol
    integer,intent(in) :: neigh(:)
    character(len=:),allocatable :: btmp
    character(len=20) :: atmp
    integer :: i,j,k
    btmp = ''
    if (neigh(1) == 0) then
      btmp = ' ---'
      return
    end if
    do i = 1,size(neigh,1)
      if (neigh(i) .ne. 0) then
        write (atmp,'(a,a,a,i0,a)') ' ',trim(i2e(mol%at(neigh(i)))),'(',neigh(i),')'
        btmp = trim(btmp)//trim(atmp)
      end if
    end do
  end function print_neighbours

!========================================================================================!

  function unique_neighbours(n,arr) result(yesno)
    implicit none
    integer,intent(in) :: n
    integer,intent(in) :: arr(n)
    logical :: yesno
    integer :: i,j
    yesno = .true.
    do i = 1,n-1
      do j = i+1,n
        if (arr(i) == arr(j)) then
          yesno = .false.
          exit
        end if
      end do
      if (.not.yesno) exit
    end do
  end function unique_neighbours

!========================================================================================!

  function determineRS(coords) result(RS)
    implicit none
    integer :: RS !> 1 for "R", -1 for "S"
    real(wp),intent(inout) :: coords(3,4)
    real(wp) :: theta
    real(wp) :: vec(3),uec(3)
    integer :: k,l,m,n

    k = 4
    !> rotate the highest prio atom onto z axis (0,0,1)
    !> first into into xz plane, therefore get the -x unit vector first
    !> and the projection of the neighbour in the xy-plane
    vec = (/-1.0d0,0.0d0,0.0d0/)
    uec(1) = coords(1,1)
    uec(2) = coords(2,1)
    uec(3) = 0.0_wp
    theta = tangle(vec,uec)
    !> then rotate around z-axis
    vec = (/0.0d0,0.0d0,1.0d0/)
    if (uec(2) .lt. 0.0_wp) theta = -theta
    do l = 1,k
      call rodrot(coords(:,l),vec,theta)
    end do

    !> afterwards angle to the z-axis
    vec = (/0.0d0,0.0d0,1.0d0/)
    uec = coords(:,1)
    theta = tangle(vec,uec)
    !> then take this angle and rotate around y
    vec = (/0.0d0,1.0d0,0.0d0/)
    do l = 1,k
      call rodrot(coords(:,l),vec,theta)
    end do

    !> as a last step tak the lowest-prio neighbour and rotate it into the xy-plane
    vec = (/-1.0d0,0.0d0,0.0d0/)
    uec(1) = coords(1,k)
    uec(2) = coords(2,k)
    uec(3) = 0.0_wp
    theta = tangle(vec,uec)
    vec = (/0.0,0.0,1.0/)
    if (uec(2) .lt. 0.0_wp) theta = -theta
    do l = 1,k
      call rodrot(coords(:,l),vec,theta)
    end do

    !determine orientation
    if (coords(2,2) .gt. coords(2,3)) then
      RS = 1
    else
      RS = -1
    end if

  end function determineRS

!========================================================================================!
!========================================================================================!
end module canonical_mod
