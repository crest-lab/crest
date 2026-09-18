module molbuilder_construct_mod
  !***********************************************
  !* A module for constructing molecules
  !* E.g. replacing functional groups and so on
  !***********************************************
  use crest_parameters
  use strucrd
  use irmsd_module
  use adjacency
  use miscdata,only:rcov
  implicit none
  private

  public :: attach

  public :: split
  interface split
    module procedure split_onbond
    module procedure split_onshared
  end interface split

  public :: binarysplit
  interface binarysplit
    module procedure binarysplit_onshared
  end interface binarysplit

!=============================================================================!
contains  !> MODULE PROCEDURES START HERE
!=============================================================================!

  subroutine attach(base,side,alignmap,new,clash, &
      & original_map,remove_base,remove_side,remove_lastx,reficn)
    !***********************************************************************
    !* This routine attaches a side-molecule to a base-molecule
    !* The assumption is that we have (at least) 3 proxy atoms
    !* that are used to match the side to the base, specified
    !* via the alignmap.
    !* Args:
    !*   base  - base molecule
    !*   side  - side chain molecule to attach to base
    !*  alignmap - (x,2) specification list of x proxy atoms in side/base
    !*    new  - newly constructed molecule
    !*
    !* Optionl args:
    !*         clash - logical, were clashes produced?
    !*  original_map - original atom position mapping
    !*   remove_base - list of atoms to remove from base upon
    !*                 constructing the new mol
    !*   remove_side - list of atoms to remove from side upon
    !*                 constructing the new mol
    !*   remove_lastx - integer (one for base and side) to remove
    !*                  final x atoms in constructing new mol
    !*   reficn       - CNs of the reference structure as integers
    !***********************************************************************
    implicit none
    !> IN/OUTPUTS
    type(coord),intent(in) :: base
    type(coord),intent(in) :: side
    integer,intent(in) :: alignmap(:,:)
    type(coord),intent(out) :: new

    logical,intent(out),optional :: clash
    integer,intent(in),optional :: original_map(:,:)
    integer,intent(in),optional :: remove_base(:)
    integer,intent(in),optional :: remove_side(:)
    integer,intent(in),optional :: remove_lastx(:)
    integer,intent(in),optional :: reficn(:)

    !> LOCAL
    type(coord) :: cutout_base,cutout_side,side_tmp
    logical,allocatable :: cutlist_base(:),cutlist_side(:)
    integer,allocatable :: current_order(:),target_order(:),idx(:)
    integer,allocatable :: revorder_base(:),revorder_side(:)
    real(wp) :: rms,Umat(3,3),shift(3),center_base(3),center_side(3)
    real(wp),allocatable :: cn(:)
    integer :: nalign,nat_new
    logical :: closcontact,failconstruct
    integer :: ii,jj,kk,ll

    character(len=*),parameter :: source = "attach()"
    real(wp),parameter :: closedist = 0.95_wp

    !> defaults/checks of the alignmap
    nalign = size(alignmap,1)
    if (nalign < 3) then
      error stop source//": alignmap needs at least 3 atoms"
    end if
    if (size(alignmap,2) .ne. 2) then
      error stop source//": alignmap has wrong dimension"
    end if
    !> alignmap(:,1) --> base atoms
    if (any(alignmap(:,1) > base%nat).or.any(alignmap(:,1) < 1)) then
      error stop source//": alignmap(:,1) has invalid values"
    end if
    !> alignmap(:,2) --> side atoms
    if (any(alignmap(:,2) > side%nat).or.any(alignmap(:,2) < 1)) then
      error stop source//": alignmap(:,2) has invalid values"
    end if

    !> helper arrays
    allocate (target_order(nalign),current_order(nalign),idx(nalign),source=0)

    !> generate cutout of base
    allocate (cutlist_base(base%nat),source=.false.)
    kk = 0
    do ii = 1,nalign
      cutlist_base(alignmap(ii,1)) = .true.
      current_order(ii) = ii
    end do
    kk = 0
    do ii = 1,base%nat
      do jj = 1,nalign
        if (alignmap(jj,1) == ii) then
          kk = kk+1
          target_order(kk) = jj
        end if
      end do
    end do
    cutout_base = base%cutout(cutlist_base)
    call molatomsort(cutout_base,nalign,current_order,target_order,idx)

    !> generate cutout of side
    allocate (cutlist_side(side%nat),source=.false.)
    do ii = 1,nalign
      cutlist_side(alignmap(ii,2)) = .true.
      current_order(ii) = ii
    end do
    kk = 0
    do ii = 1,side%nat
      do jj = 1,nalign
        if (alignmap(jj,2) == ii) then
          kk = kk+1
          target_order(kk) = jj
        end if
      end do
    end do
    cutout_side = side%cutout(cutlist_side)
    call molatomsort(cutout_side,nalign,current_order,target_order,idx)

    !> determine offset and translate
    do ii = 1,3
      center_base(ii) = sum(cutout_base%xyz(ii,:))/real(cutout_base%nat)
      center_side(ii) = sum(cutout_side%xyz(ii,:))/real(cutout_side%nat)
    end do
    shift(1:3) = center_base(1:3)-center_side(1:3)

    !> determine the rotation matrix and translation
    rms = rmsd(cutout_base,cutout_side,rotmat=Umat)

    !> build the new rotated and translated side molecule
    side_tmp = side
    side_tmp%xyz = matmul(Umat,side_tmp%xyz)
    !> determine offset and translate
    shift = 0.0_wp
    do jj = 1,nalign
      shift(:) = shift(:)+base%xyz(:,alignmap(jj,1))-side_tmp%xyz(:,alignmap(jj,2))
    end do
    shift(:) = shift(:)/real(nalign)
    do ii = 1,side_tmp%nat
      side_tmp%xyz(1:3,ii) = side_tmp%xyz(1:3,ii)+shift(1:3)
    end do

    !> construct the new molecule, re-use cutlists
    cutlist_base(:) = .false.
    cutlist_side(:) = .false.
    do ii = 1,nalign
      !> always remove the align atoms on the side chain
      cutlist_side(alignmap(ii,2)) = .true.
    end do
    if (present(remove_base)) then
      !> optional removal on the base side
      kk = size(remove_base,1)
      do ii = 1,kk
        cutlist_base(remove_base(ii)) = .true.
      end do
    end if
    if (present(remove_side)) then
      !> optional removal on the base side
      kk = size(remove_side,1)
      do ii = 1,kk
        cutlist_side(remove_side(ii)) = .true.
      end do
    end if
    if (present(remove_lastx)) then
      !> useful for not transfering appended capping atoms
      do ii = base%nat-(remove_lastx(1)-1),base%nat
        cutlist_base(ii) = .true.
      end do
      do ii = side%nat-(remove_lastx(2)-1),side%nat
        cutlist_side(ii) = .true.
      end do
    end if
    allocate (revorder_base(base%nat),source=0)
    allocate (revorder_side(side%nat),source=0)
    if (present(original_map)) then
      !> if we have info on the original order, avoid all atoms
      !> present in all fragments
      closcontact = .false.
      do ii = 1,size(original_map,1)
        if (all(original_map(ii,:) .ne. 0)) then
          if (.not.closcontact) then
            cutlist_side(original_map(ii,2)) = .true.
          else
            cutlist_base(original_map(ii,1)) = .true.
          end if
          jj = original_map(ii,1)
          if (.not.closcontact) then
            do kk = 1,size(original_map,1)
              ll = original_map(kk,2)
              if (ll .ne. 0.and.original_map(kk,1) .eq. 0) then
                !write (*,*) ii,jj,ll,sqrt(sum((side_tmp%xyz(:,ll)-base%xyz(:,jj))**2))
                if ((sum((side_tmp%xyz(:,ll)-base%xyz(:,jj))**2)) < closedist) then
                  closcontact = .true.
                  cutlist_side(original_map(ii,2)) = .false.
                  cutlist_base(original_map(ii,1)) = .true.
                end if
              end if
            end do
          end if

        end if
        if (original_map(ii,1) .ne. 0) then
          revorder_base(original_map(ii,1)) = ii
        end if
        if (original_map(ii,2) .ne. 0) then
          revorder_side(original_map(ii,2)) = ii
        end if
      end do
      if (closcontact) then
        do kk = 1,size(original_map,1)
          if (all(original_map(kk,:) .ne. 0).and. &
              .not.any(original_map(kk,1) .eq. alignmap(:,1))) then
            cutlist_side(original_map(kk,2)) = .false.
            cutlist_base(original_map(kk,1)) = .true.
          end if
        end do
      end if
    end if
    kk = max(maxval(revorder_base),maxval(revorder_side))
    new%nat = 0
    do ii = 1,base%nat
      if (.not.cutlist_base(ii)) then
        new%nat = new%nat+1
        if (revorder_base(ii) .eq. 0) then
          kk = kk+1
          revorder_base(ii) = kk
        end if
      end if
    end do
    do ii = 1,side%nat
      if (.not.cutlist_side(ii)) then
        new%nat = new%nat+1
        if (revorder_side(ii) .eq. 0) then
          kk = kk+1
          revorder_side(ii) = kk
        end if
      end if
    end do
    allocate (new%at(new%nat),source=0)
    allocate (new%xyz(3,new%nat),source=0.0_wp)
    !kk = 0
    failconstruct = .false.
    do ii = 1,base%nat
      if (.not.cutlist_base(ii)) then
        !kk = kk+1
        kk = revorder_base(ii)
        if (kk > new%nat) then
          failconstruct = .true.
          exit
        end if
        new%at(kk) = base%at(ii)
        new%xyz(1:3,kk) = base%xyz(1:3,ii)
      end if
    end do
    do ii = 1,side%nat
      if (.not.cutlist_side(ii)) then
        !kk = kk+1
        kk = revorder_side(ii)
        if (kk > new%nat) then
          failconstruct = .true.
          exit
        end if
        new%at(kk) = side_tmp%at(ii)
        new%xyz(1:3,kk) = side_tmp%xyz(1:3,ii)
      end if
    end do

    if (present(clash)) then
      clash = .false.
      !> TODO implement geometric clash check
      if (failconstruct) then
        clash = .true.
      else if (present(reficn)) then
        call new%get_cn(cn)
        do ii = 1,new%nat
          if (nint(cn(ii)) .ne. reficn(ii)) then
            clash = .true.
            !write(*,*) "removing clash"
            exit
          end if
        end do
      end if
    end if
  end subroutine attach

!==============================================================================!

  subroutine split_onbond(input,bond,base,side,wbo)
    implicit none
    !> IN/OUTPUTS
    type(coord),intent(in) :: input
    integer,intent(in) :: bond(2)
    type(coord),intent(out) :: base
    type(coord),intent(out) :: side
    !> OPTIONAL
    real(wp),intent(in),optional :: wbo(input%nat,input%nat)
    !> LOCAL
    real(wp),allocatable :: cn(:),wbofake(:,:)
    integer :: V,fbase,fside
    integer,allocatable :: A(:,:),Anew(:,:)
    integer,allocatable :: frag(:),fragnew(:)
    logical,allocatable :: in_ring(:,:)
    integer,allocatable :: path_tmp(:)
    integer :: npath,nbase,nside
    logical :: needs_capping
    real(wp) :: distcap
    integer :: ii,jj,kk

    character(len=*),parameter :: source = "split_onbond()"

    !> checks
    if (any(bond(:) > input%nat).or.(any(bond(:) < 1))) then
      error stop source//": bond() has invalid atom specification"
    end if

    !> we will be working with graphs. define number of vertices = #atoms
    V = input%nat

    !> set up adjacency matrix
    if (present(wbo)) then
      call wbo2adjacency(V,wbo,A,0.01_wp)
    else
      call input%cn_to_bond(cn,wbofake)
      call wbo2adjacency(V,wbofake,A,0.01_wp)
      deallocate (wbofake,cn)
    end if

    !> get fragment array (indicates which atom is on which fragment)
    call setup_fragments(V,A,frag)

    !> if the two atoms are already on different fragments we
    !> can immediatly proceed, otherwise we need to make some cuts
    if (frag(bond(1)) .eq. frag(bond(2))) then
      !> check if the two atoms actually correspond to an bond
      if (A(bond(1),bond(2)) .eq. 0) then
        error stop source//": specified atoms are not directly connected!"
      end if

      needs_capping = .true.
      !> safety: check for rings (need more cuts)
      call check_rings_min(V,A,in_ring)
      !allocate (path_tmp(V),source=0)
      !call get_ring_min(V,A,bond(1),bond(2),path_tmp,npath)
      !if (npath > 0) then
      if (in_ring(bond(1),bond(2))) then
        !> TODO - Implement actual fallback
        error stop "Bond is in a ring! Can not split."
      else
        !> delete the bond and set up new fragment array
        allocate (Anew(V,V),source=A)
        Anew(bond(1),bond(2)) = 0
        Anew(bond(2),bond(1)) = 0
        call setup_fragments(V,Anew,frag)
      end if
    else
      needs_capping = .false.
    end if

    !> the logic here is: the fragment containing atom "bond(1)"
    !> will be base, and includes everything except the fragment
    !> containing atom "bond(2)"; even further fragments containing neither
    fbase = frag(bond(1))
    fside = frag(bond(2))

    nbase = 0
    nside = 0
    do ii = 1,input%nat
      if (frag(ii) == fside) then
        nside = nside+1
      else
        nbase = nbase+1
      end if
    end do

    !> capping off the cut bond (with H or smthg else)
    if (needs_capping) then
      nside = nside+1
      nbase = nbase+1
    end if

    !> construct output mols
    base%nat = nbase
    allocate (base%at(nbase),source=0)
    allocate (base%xyz(3,nbase),source=0.0_wp)
    side%nat = nside
    allocate (side%at(nside),source=0)
    allocate (side%xyz(3,nside),source=0.0_wp)

    jj = 0
    kk = 0
    do ii = 1,input%nat
      if (frag(ii) == fside) then
        jj = jj+1
        side%at(jj) = input%at(ii)
        side%xyz(1:3,jj) = input%xyz(1:3,ii)
      else
        kk = kk+1
        base%at(kk) = input%at(ii)
        base%xyz(1:3,kk) = input%xyz(1:3,ii)
      end if
    end do

    if (needs_capping) then
      jj = jj+1
      side%at(jj) = 1 !> hydrogen
      side%xyz(1:3,jj) = input%xyz(1:3,bond(1))
      distcap = (rcov(1)+rcov(input%at(bond(2))))*(3.0_wp/4.0_wp)
      call place_at_distance(input%xyz(1:3,bond(2)),side%xyz(1:3,jj),distcap)
      kk = kk+1
      base%at(kk) = 1   !> hydrogen
      base%xyz(1:3,kk) = input%xyz(1:3,bond(2))
      distcap = (rcov(1)+rcov(input%at(bond(1))))*(3.0_wp/4.0_wp)
      call place_at_distance(input%xyz(1:3,bond(1)),base%xyz(1:3,kk),distcap)
    end if

    if (allocated(Anew)) deallocate (Anew)
    if (allocated(frag)) deallocate (frag)
  end subroutine split_onbond

!==============================================================================!
  subroutine split_onshared(input,sharedlist,structures,sharedmap,&
      & wbo,ncap,position_mapping)
    implicit none
    !> IN/OUTPUTS
    type(coord),intent(in) :: input
    integer,intent(in) :: sharedlist(:)
    type(coord),intent(out),allocatable :: structures(:)
    integer,intent(out),allocatable :: sharedmap(:,:)
    !> OPTIONAL
    real(wp),intent(in),optional :: wbo(input%nat,input%nat)
    integer,intent(out),allocatable,optional :: ncap(:)
    integer,intent(out),allocatable,optional :: position_mapping(:,:)
    !> LOCAL
    type(coord) :: shared
    real(wp),allocatable :: cn(:),wbofake(:,:)
    integer :: V,fbase,fside,nshared,ftmp
    integer,allocatable :: A(:,:),Anew(:,:)
    integer,allocatable :: frag(:),fragnew(:),molassign(:)
    integer,allocatable :: ncapped(:)
    integer,allocatable :: pos_map(:,:)
    logical,allocatable :: in_ring(:,:)
    integer,allocatable :: path_tmp(:)
    integer,allocatable :: number_of_neighbours(:)
    logical,allocatable :: terminal_atom(:)
    logical,allocatable :: connected_to_share(:)
    logical,allocatable :: assign_to_mols(:,:)
    logical,allocatable :: unassigned_fragments(:)
    logical,allocatable :: capping_mapping(:,:)
    logical,allocatable :: methylizemapping(:,:)

    integer :: npath,nfrag,nfragnew,nbonds
    real(wp) :: distcap,methylproxy(3,3)
    integer :: ii,jj,jjj,kk,sii,sjj,M,mm,nn,ll,lll

    integer :: bond(2)

    character(len=*),parameter :: source = "split_onshared()"

    !> we will be working with graphs. define number of vertices = #atoms
    V = input%nat
    !> how many atoms are shared
    nshared = size(sharedlist,1)

    !> checks
    if (any(sharedlist(:) > input%nat).or.(nshared < 1)) then

      error stop source//": sharedlist() has invalid atom specification"
    end if

    !> set up adjacency matrix
    if (present(wbo)) then
      call wbo2adjacency(V,wbo,A,0.01_wp)
    else
      call input%cn_to_bond(cn,wbofake)
      call wbo2adjacency(V,wbofake,A,0.01_wp)
      deallocate (wbofake,cn)
    end if

!> ----------------------------------------------------------------------------
!> BOOKKEEPING START
!> ----------------------------------------------------------------------------
    !> get fragment array (indicates which atom is on which fragment)
    call setup_fragments(V,A,frag)

    !> some other mappings
    allocate (number_of_neighbours(V),source=0)
    do ii = 1,V
      number_of_neighbours(ii) = sum(A(:,ii))
    end do
    allocate (terminal_atom(V),source=.false.)
    do ii = 1,V
      terminal_atom(ii) = (number_of_neighbours(ii) == 1)
    end do

    !> The cutting logic starts here.
    !> First, we need to identify all atoms actually sharing fragments with the shared atoms
    allocate (connected_to_share(V),source=.false.)
    do ii = 1,V
      ftmp = frag(ii)
      do jj = 1,nshared
        if (ftmp == frag(sharedlist(jj))) then
          connected_to_share(ii) = .true.
          exit
        end if
      end do
    end do
    !> all atoms NOT part of that list will be present in both output fragments

    !> Then, we take the graph and construct a new one with detachted "shared" atoms
    allocate (Anew(V,V),source=A)
    do ii = 1,nshared
      sii = sharedlist(ii)
      do jj = 1,V
        if (jj == sii) cycle
        if ((A(jj,sii) == 1).and. &  !> look at existing bonds to shared section
         & .not.any(sharedlist(:) == jj).and. & !> exctept to other shared section atoms
         & .not.terminal_atom(jj)) then  !> and except terminal atoms (directly bound to shared section)
          Anew(jj,sii) = 0
          Anew(sii,jj) = 0
        end if
      end do
    end do
    !> get new fragments
    call setup_fragments(V,Anew,fragnew)
    !> due to the cutting we should have at least two more fragments in total, compared to before.
    !> for now we only handle exactly that +2 fragment case
    nfrag = maxval(frag)
    nfragnew = maxval(fragnew)
    M = 0
    if (nfragnew < (nfrag+2)) then
      error stop source//": system fragmentation yields currently unhandled edge-case"
    else
      M = nfragnew-nfrag
      !> now we can check asignment to new, split-up fragments
      !> we distinguish between the shared secion (:,1), and all M others (:,2:M)
      allocate (assign_to_mols(V,M+1),source=.false.)
      allocate (unassigned_fragments(nfragnew),source=.true.)
      !> first, get the separated shared-region
      mm = M+1
      do ii = 1,nshared
        sii = fragnew(sharedlist(ii))
        unassigned_fragments(sii) = .false.
        do jj = 1,V
          if (fragnew(jj) == sii) then
            assign_to_mols(jj,1:mm) = .true.
          end if
        end do
      end do
      !> then, all atoms that had no connection to the shared region
      !> and hence are present everywhere (except the shard region)
      do ii = 1,V
        if (.not.connected_to_share(ii)) then
          sii = fragnew(ii)
          unassigned_fragments(sii) = .false.
          do jj = 1,V
            if (fragnew(jj) == sii) then
              assign_to_mols(jj,2:mm) = .true.
            end if
          end do
        end if
      end do

      !!> exactly M fragments should be remaining
      if (count(unassigned_fragments) .ne. M) then
        error stop source//": wrong number of unassigned_fragments"
      end if
      mm = 1
      do ii = 1,V
        if (all(assign_to_mols(ii,:).eqv..false.)) then
          mm = mm+1 !> molecule new assignment number (except mm=1 because that is the shared region)
          sii = fragnew(ii) !> the reference fragment
          do jj = 1,V
            if (fragnew(jj) .eq. sii) assign_to_mols(jj,mm) = .true.
          end do
        end if
      end do
    end if

    !> prepare mapping and capping.
    allocate (capping_mapping(V,M),source=.false.)
    allocate (methylizemapping(nshared,M),source=.false.)
    !> First, simple, chemoinformatic rules
    do ii = 1,M
      mm = ii+1
      do jj = 1,nshared
        nbonds = 0
        kk = sharedlist(jj)
        do ll = 1,V
          if ((A(ll,kk) == 1).and.assign_to_mols(ll,mm)) then
            nbonds = nbonds+1
          end if
        end do
        if (nbonds == 1.and.input%at(kk) == 6) then
          methylizemapping(jj,ii) = .true.
        end if
      end do
    end do
    !> then "regular" capping, we determine original atoms as proxy for
    !> the cap (and later adjust the bondlength)
    !> Entries will be .true. for atoms that need to be added to fragment
    do ii = 1,M
      mm = ii+1
      do jj = 1,nshared
        if (methylizemapping(jj,ii)) cycle
        do kk = 1,V
          if (A(kk,sharedlist(jj)) == 1.and..not.assign_to_mols(kk,mm)) then
            capping_mapping(kk,ii) = .true.
          end if
        end do
      end do
    end do

!> ----------------------------------------------------------------------------
!> BOOKKEEPING END
!> ----------------------------------------------------------------------------
!> MOLECULE CONSTRUCTION START
!> ----------------------------------------------------------------------------

    allocate (structures(M)) !> we know that splitting produces M fragment
    allocate (sharedmap(nshared,M),source=0)
    allocate (ncapped(M),source=0)
    allocate (pos_map(V,M),source=0)

    do ii = 1,M
      mm = ii+1
      !> count atoms, allocate
      nn = count(assign_to_mols(:,mm),1)+ &
        &  count(capping_mapping(:,ii),1)+ &
        &  count(methylizemapping(:,ii),1)*3
      structures(ii)%nat = nn
      allocate (structures(ii)%at(nn),source=2)
      allocate (structures(ii)%xyz(3,nn),source=0.0_wp)
      kk = 0
      jjj = 0
      !> directly assigned atoms
      do jj = 1,V
        if (assign_to_mols(jj,mm)) then
          kk = kk+1
          structures(ii)%at(kk) = input%at(jj)
          structures(ii)%xyz(1:3,kk) = input%xyz(1:3,jj)
          if (any(sharedlist(:) == jj)) then
            jjj = jjj+1
            sharedmap(jjj,ii) = kk
          end if
          pos_map(jj,ii) = kk
        end if
      end do
      !> capping atoms
      do jj = 1,V
        if (capping_mapping(jj,ii)) then
          kk = kk+1
          ncapped(ii) = ncapped(ii)+1
          structures(ii)%at(kk) = 1 !input%at(jj)
          structures(ii)%xyz(1:3,kk) = input%xyz(1:3,jj)
          !> repair distance for capping atoms
          do ll = 1,V
            if ((A(ll,jj) == 1).and.assign_to_mols(ll,mm)) then
              distcap = (rcov(1)+rcov(input%at(ll)))*(3.0_wp/4.0_wp)
              call place_at_distance(input%xyz(1:3,ll),structures(ii)%xyz(1:3,kk),distcap)
            end if
          end do
        end if
      end do
      !> methylation atoms
      do jj = 1,nshared
        if (methylizemapping(jj,ii)) then
          sjj = sharedlist(jj)
          do ll = 1,V
            if (A(ll,sjj) == 1.and.assign_to_mols(ll,mm)) then
              call methylize(input%xyz(1:3,sjj),input%xyz(1:3,ll),methylproxy)
              do lll = 1,3
                kk = kk+1
                ncapped(ii) = ncapped(ii)+1
                structures(ii)%at(kk) = 1
                structures(ii)%xyz(1:3,kk) = methylproxy(1:3,lll)
              end do
              exit
            end if
          end do
        end if
      end do
    end do

!> ----------------------------------------------------------------------------
!> MOLECULE CONSTRUCTION END
!> ----------------------------------------------------------------------------

    if (present(ncap)) then
      call move_alloc(ncapped,ncap)
    end if

    if (present(position_mapping)) then
      call move_alloc(pos_map,position_mapping)
    end if

    if (allocated(pos_map)) deallocate (pos_map)
    if (allocated(ncapped)) deallocate (ncapped)
    if (allocated(assign_to_mols)) deallocate (assign_to_mols)
    if (allocated(capping_mapping)) deallocate (capping_mapping)
    if (allocated(unassigned_fragments)) deallocate (unassigned_fragments)
    if (allocated(number_of_neighbours)) deallocate (number_of_neighbours)
    if (allocated(terminal_atom)) deallocate (terminal_atom)
    if (allocated(connected_to_share)) deallocate (connected_to_share)
    if (allocated(Anew)) deallocate (Anew)
    if (allocated(frag)) deallocate (frag)
  end subroutine split_onshared
!==============================================================================!
  subroutine binarysplit_onshared(input,sharedlist,structures,sharedmap,&
      & wbo,ncap,position_mapping)
    implicit none
    !> IN/OUTPUTS
    type(coord),intent(in) :: input
    integer,intent(in) :: sharedlist(:)
    type(coord),intent(out),allocatable :: structures(:)
    integer,intent(out),allocatable :: sharedmap(:,:)
    !> OPTIONAL
    real(wp),intent(in),optional :: wbo(input%nat,input%nat)
    integer,intent(out),allocatable,optional :: ncap(:)
    integer,intent(out),allocatable,optional :: position_mapping(:,:)
    !> LOCAL
    type(coord) :: shared
    real(wp),allocatable :: cn(:),wbofake(:,:)
    integer :: V,fbase,fside,nshared,ftmp
    integer,allocatable :: A(:,:),Anew(:,:)
    integer,allocatable :: frag(:),fragnew(:),molassign(:)
    integer,allocatable :: ncapped(:)
    integer,allocatable :: pos_map(:,:)
    logical,allocatable :: in_ring(:,:)
    integer,allocatable :: path_tmp(:)
    integer,allocatable :: number_of_neighbours(:)
    logical,allocatable :: terminal_atom(:)
    logical,allocatable :: connected_to_share(:)
    logical,allocatable :: assign_to_mols(:,:)
    logical,allocatable :: unassigned_fragments(:)
    logical,allocatable :: capping_mapping(:,:)
    logical,allocatable :: methylizemapping(:,:)

    integer :: npath,nfrag,nfragnew,nbonds
    real(wp) :: distcap,methylproxy(3,3)
    integer :: ii,jj,jjj,kk,sii,sjj,M,mm,nn,ll,lll

    integer :: bond(2)

    character(len=*),parameter :: source = "binarysplit_onshared()"

    !> we will be working with graphs. define number of vertices = #atoms
    V = input%nat
    !> how many atoms are shared
    nshared = size(sharedlist,1)

    !> checks
    if (any(sharedlist(:) > input%nat).or.(nshared < 1)) then

      error stop source//": sharedlist() has invalid atom specification"
    end if

    !> set up adjacency matrix
    if (present(wbo)) then
      call wbo2adjacency(V,wbo,A,0.01_wp)
    else
      call input%cn_to_bond(cn,wbofake)
      call wbo2adjacency(V,wbofake,A,0.01_wp)
      deallocate (wbofake,cn)
    end if

!> ----------------------------------------------------------------------------
!> BOOKKEEPING START
!> ----------------------------------------------------------------------------
    !> get fragment array (indicates which atom is on which fragment)
    call setup_fragments(V,A,frag)

    !> some other mappings
    allocate (number_of_neighbours(V),source=0)
    do ii = 1,V
      number_of_neighbours(ii) = sum(A(:,ii))
    end do
    allocate (terminal_atom(V),source=.false.)
    do ii = 1,V
      terminal_atom(ii) = (number_of_neighbours(ii) == 1)
    end do

    !> The cutting logic starts here.
    !> First, we need to identify all atoms actually sharing fragments with the shared atoms
    allocate (connected_to_share(V),source=.false.)
    do ii = 1,V
      ftmp = frag(ii)
      do jj = 1,nshared
        if (ftmp == frag(sharedlist(jj))) then
          connected_to_share(ii) = .true.
          exit
        end if
      end do
    end do
    !> all atoms NOT part of that list will be present in both output fragments

    !> Then, we take the graph and construct a new one with detachted "shared" atoms
    nfrag = maxval(frag)
    allocate (Anew(V,V),source=A)
    iiloop1: do ii = 1,nshared
      sii = sharedlist(ii)
      do jj = 1,V
        if (jj == sii) cycle
        if ((A(jj,sii) == 1).and. &  !> look at existing bonds to shared section
         & any(sharedlist(:) == jj).and. & !> must be connected to other shared section atoms
         & .not.terminal_atom(jj)) then  !> and except terminal atoms (directly bound to shared section)
          Anew(jj,sii) = 0
          Anew(sii,jj) = 0

          !> get new fragments
          call setup_fragments(V,Anew,fragnew)
          nfragnew = maxval(fragnew)
          if ((nfragnew-nfrag) == 1) then
            exit iiloop1
          end if
        end if
      end do
    end do iiloop1

    M = 2
    if (nfragnew == nfrag) then
      error stop source//": system fragmentation yields currently unhandled edge-case"
    else
      !> now we can check asignment to new, split-up fragments
      !> we distinguish between the shared secion (:,1), and all M others (:,2:M)
      allocate (assign_to_mols(V,M+1),source=.false.)
      allocate (unassigned_fragments(nfragnew),source=.true.)
      !> distribute the separated fragments. Also assign the shared atoms themselves
      mm = 1
      do ii = 1,nshared
        sii = fragnew(sharedlist(ii))
        assign_to_mols(sharedlist(ii),:) = .true.
        if (.not.unassigned_fragments(sii)) cycle
        mm = mm+1
        unassigned_fragments(sii) = .false.
        do jj = 1,V
          if (fragnew(jj) == sii) then
            assign_to_mols(jj,mm) = .true.
          end if
        end do
      end do
      !> then, all atoms that had no connection to the shared region
      !> and hence are present everywhere (except the shard region)
      do ii = 1,V
        if (.not.connected_to_share(ii)) then
          sii = fragnew(ii)
          unassigned_fragments(sii) = .false.
          do jj = 1,V
            if (fragnew(jj) == sii) then
              assign_to_mols(jj,2:mm) = .true.
            end if
          end do
        end if
      end do

      !!> exactly 0 fragments should be remaining because we had the even split
      if (count(unassigned_fragments) .ne. 0) then
        error stop source//": wrong number of unassigned_fragments"
      end if
    end if

    !> prepare mapping and capping.
    allocate (capping_mapping(V,M),source=.false.)
    allocate (methylizemapping(nshared,M),source=.false.)
    !> First, simple, chemoinformatic rules
    do ii = 1,M
      mm = ii+1
      do jj = 1,nshared
        nbonds = 0
        sii = sharedlist(jj)
        do ll = 1,V
          if ((A(ll,sii) == 1).and.assign_to_mols(ll,mm)) then
            nbonds = nbonds+1
          end if
        end do
        if (nbonds == 1.and.input%at(sii) == 6) then
          methylizemapping(jj,ii) = .true.
        end if
      end do
    end do
    !> then "regular" capping, we determine original atoms as proxy for
    !> the cap (and later adjust the bondlength)
    !> Entries will be .true. for atoms that need to be added to fragment
    do ii = 1,M
      mm = ii+1
      do jj = 1,nshared
        if (methylizemapping(jj,ii)) cycle
        do kk = 1,V
          if (A(kk,sharedlist(jj)) == 1.and..not.assign_to_mols(kk,mm)) then
            capping_mapping(kk,ii) = .true.
          end if
        end do
      end do
    end do

!> ----------------------------------------------------------------------------
!> BOOKKEEPING END
!> ----------------------------------------------------------------------------
!> MOLECULE CONSTRUCTION START
!> ----------------------------------------------------------------------------

    allocate (structures(M)) !> we know that splitting produces M fragment
    allocate (sharedmap(nshared,M),source=0)
    allocate (ncapped(M),source=0)
    allocate (pos_map(V,M),source=0)

    do ii = 1,M
      mm = ii+1
      !> count atoms, allocate
      nn = count(assign_to_mols(:,mm),1)+ &
        &  count(capping_mapping(:,ii),1)+ &
        &  count(methylizemapping(:,ii),1)*3
      structures(ii)%nat = nn
      allocate (structures(ii)%at(nn),source=2)
      allocate (structures(ii)%xyz(3,nn),source=0.0_wp)
      kk = 0
      jjj = 0
      !> directly assigned atoms
      do jj = 1,V
        if (assign_to_mols(jj,mm)) then
          kk = kk+1
          structures(ii)%at(kk) = input%at(jj)
          structures(ii)%xyz(1:3,kk) = input%xyz(1:3,jj)
          if (any(sharedlist(:) == jj)) then
            jjj = jjj+1
            sharedmap(jjj,ii) = kk
          end if
          pos_map(jj,ii) = kk
        end if
      end do
      !> capping atoms
      do jj = 1,V
        if (capping_mapping(jj,ii)) then
          kk = kk+1
          ncapped(ii) = ncapped(ii)+1
          structures(ii)%at(kk) = 1 !input%at(jj)
          structures(ii)%xyz(1:3,kk) = input%xyz(1:3,jj)
          !> repair distance for capping atoms
          do ll = 1,V
            if ((A(ll,jj) == 1).and.assign_to_mols(ll,mm)) then
              distcap = (rcov(1)+rcov(input%at(ll)))*(3.0_wp/4.0_wp)
              call place_at_distance(input%xyz(1:3,ll),structures(ii)%xyz(1:3,kk),distcap)
            end if
          end do
        end if
      end do
      !> methylation atoms
      do jj = 1,nshared
        if (methylizemapping(jj,ii)) then
          sjj = sharedlist(jj)
          do ll = 1,V
            if (A(ll,sjj) == 1.and.assign_to_mols(ll,mm)) then
              call methylize(input%xyz(1:3,sjj),input%xyz(1:3,ll),methylproxy)
              do lll = 1,3
                kk = kk+1
                ncapped(ii) = ncapped(ii)+1
                structures(ii)%at(kk) = 1
                structures(ii)%xyz(1:3,kk) = methylproxy(1:3,lll)
              end do
              exit
            end if
          end do
        end if
      end do
    end do

!> ----------------------------------------------------------------------------
!> MOLECULE CONSTRUCTION END
!> ----------------------------------------------------------------------------

    if (present(ncap)) then
      call move_alloc(ncapped,ncap)
    end if

    if (present(position_mapping)) then
      call move_alloc(pos_map,position_mapping)
    end if

    if (allocated(pos_map)) deallocate (pos_map)
    if (allocated(ncapped)) deallocate (ncapped)
    if (allocated(assign_to_mols)) deallocate (assign_to_mols)
    if (allocated(capping_mapping)) deallocate (capping_mapping)
    if (allocated(unassigned_fragments)) deallocate (unassigned_fragments)
    if (allocated(number_of_neighbours)) deallocate (number_of_neighbours)
    if (allocated(terminal_atom)) deallocate (terminal_atom)
    if (allocated(connected_to_share)) deallocate (connected_to_share)
    if (allocated(Anew)) deallocate (Anew)
    if (allocated(frag)) deallocate (frag)
  end subroutine binarysplit_onshared

!=============================================================================!
  pure subroutine place_at_distance(ref,moving,dist,tol,ierr)
    implicit none
    real(wp),intent(in)    :: ref(3)
    real(wp),intent(inout) :: moving(3)
    real(wp),intent(in)    :: dist
    real(wp),intent(in),optional :: tol
    integer,intent(out),optional :: ierr

    real(wp) :: v(3),r,t
    real(wp),parameter :: default_tol = 1.0d-12

    t = default_tol
    if (present(tol)) t = tol

    v = moving-ref
    r = sqrt(dot_product(v,v))

    if (present(ierr)) ierr = 0
    if (r <= t) then
      !> Direction is undefined (moving ~= ref), so we cannot keep the axis.
      if (present(ierr)) ierr = 1
      return
    end if

    moving = ref+(dist/r)*v
  end subroutine place_at_distance

!==============================================================================!

  subroutine methylize(x,n,h_new,h_aligned,r_ch,ok)
    use geo
    implicit none
    real(wp),intent(in)  :: x(3),n(3)
    real(wp),intent(out) :: h_new(3,3)
    real(wp),intent(out),optional :: h_aligned(3)
    real(wp),intent(in),optional :: r_ch
    logical,intent(out),optional :: ok

    real(wp) :: rch,invs3,pi
    real(wp) :: Htmpl(3,4)
    real(wp) :: dir(3),a(3),axis(3),tmp(3)
    real(wp) :: angle,d
    integer  :: k
    real(wp),parameter :: eps = 1.0e-12_wp
    logical :: success

    pi = acos(-1.0_wp)

    rch = 1.09_wp*aatoau
    if (present(r_ch)) rch = r_ch

    !>--- methane proxy (C at origin), perfect tetrahedral directions
    invs3 = 1.0_wp/sqrt(3.0_wp)
    Htmpl(:,1) = rch*(/1.0_wp,1.0_wp,1.0_wp/)*invs3
    Htmpl(:,2) = rch*(/1.0_wp,-1.0_wp,-1.0_wp/)*invs3
    Htmpl(:,3) = rch*(/-1.0_wp,1.0_wp,-1.0_wp/)*invs3
    Htmpl(:,4) = rch*(/-1.0_wp,-1.0_wp,1.0_wp/)*invs3

    !>--- desired direction: outward from x away from neighbor n
    dir = n-x
    if (vec_len(dir) < eps) then
      h_new = 0.0_wp
      if (present(h_aligned)) h_aligned = 0.0_wp
      if (present(ok)) ok = .false.
      return
    end if
    call unitv(dir)

    !>--- template direction: use H1 bond direction (unit)
    a = Htmpl(:,1)
    call unitv(a)

    !>--- compute rotation: rotate a
    call crosp(a,dir,axis)
    d = dotp(a,dir,3)

    if (vec_len(axis) < eps) then
      !> a parallel or anti-parallel to dir
      if (d > 0.0_wp) then
        angle = 0.0_wp
        axis = (/1.0_wp,0.0_wp,0.0_wp/)  ! dummy axis (unused)
      else
        !> 180° rotation: choose any axis perpendicular to a
        tmp = (/1.0_wp,0.0_wp,0.0_wp/)
        if (abs(dotp(a,tmp,3)) > 0.9_wp) tmp = (/0.0_wp,1.0_wp,0.0_wp/)
        call crosp(a,tmp,axis)
        call unitv(axis)
        angle = pi
      end if
    else
      call unitv(axis)
      angle = tangle(a,dir)   !> radians, [0..pi]
    end if

    !> --- rotate & translate: H = x + R * Htmpl
    !> Return the other three (2..4) as "new"
    do k = 2,4
      tmp = Htmpl(:,k)
      if (abs(angle) > 0.0_wp) call rodrot(tmp,axis,angle)
      h_new(:,k-1) = x+tmp
    end do

    if (present(h_aligned)) then
      tmp = Htmpl(:,1)
      if (abs(angle) > 0.0_wp) call rodrot(tmp,axis,angle)
      h_aligned = x+tmp
    end if

    success = .true.
    if (present(ok)) ok = success
  end subroutine methylize

!=============================================================================!
!#############################################################################!
!=============================================================================!
end module molbuilder_construct_mod

