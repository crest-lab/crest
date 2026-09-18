module molbuilder_construct_list
  !**************************************************
  !* Bookkeeping module for reconstructing molecules
  !**************************************************
  use crest_parameters
  use crest_calculator
  use strucrd
  use quicksort_interface,only:qqsorti
  implicit none
  private

  type :: split_atms
    integer :: natms = 0
    integer,allocatable :: atms(:)
  end type split_atms

  type :: origin_map
    integer :: natms = 0
    integer,allocatable :: map(:)
  end type origin_map

  type :: construct_layer
    integer :: id = 0
    integer :: parent = 0
    integer :: parentnode = 0
    integer :: nnodes = 0
    real(wp) :: inverse_depth = 0.0_wp
    type(coord),allocatable :: node(:)
    integer,allocatable :: childlayer(:)
    integer,allocatable :: alignmap(:,:)
    integer,allocatable :: ncapped(:)
    integer,allocatable :: position_mapping(:,:)
    type(origin_map),allocatable :: origin(:)
    !> initial/reconstructed molecule
    integer :: nmols = 0
    type(coord),allocatable :: mols(:)
    !> reference
    type(coord) :: refmol
    real(wp),allocatable :: refcn(:)
    integer,allocatable :: reficn(:)
  end type construct_layer

  type :: construct_queue
    integer :: layer = 0
    integer :: node = 0
    integer :: duplicate_of_queue = 0
    character(len=:),allocatable :: workdir
    character(len=:),allocatable :: file
    type(calcdata) :: calc
  end type construct_queue

  type :: construct_heap
    !> The gerneral mapping
    integer :: nlayer = 0
    type(construct_layer),allocatable :: layer(:)
    !> the queue to treat endpoints
    integer :: nqueue = 0
    type(construct_queue),allocatable :: queue(:)
    !> originial data
    type(coord) :: originmol
    character(len=:),allocatable :: origindir
    type(calcdata),pointer :: origincalc
  contains
    procedure map_origins_for_layer
    procedure find_current_position
    procedure count_endpoints
    procedure setup_queue
    procedure fill_inverse_depth
  end type construct_heap

  !> exported types
  public :: split_atms
  public :: construct_heap
  !> exported helper functions
  public :: add_to_splitqueue

!=============================================================================!
contains  !> MODULE PROCEDURES START HERE
!=============================================================================!

  subroutine add_to_splitqueue(splitqueue,raw_split)
    implicit none
    type(split_atms),intent(inout),allocatable :: splitqueue(:)
    integer,intent(in) :: raw_split(:)
    type(split_atms) :: tmpsplit
    type(split_atms),allocatable :: tmpsplitqueue(:)
    integer :: ii,n,nall,nallnew
    logical :: duplicate

    n = size(raw_split,1)
    allocate (tmpsplit%atms(n))
    tmpsplit%natms = n
    tmpsplit%atms(:) = raw_split(:)
    call qqsorti(tmpsplit%atms,1,n)

    if (.not.allocated(splitqueue)) then
      allocate (splitqueue(1))
      splitqueue(1) = tmpsplit
    else
      nall = size(splitqueue,1)
      nallnew = nall+1
      duplicate = .false.
      do ii = 1,nall
        if (all(splitqueue(ii)%atms(:) .eq. tmpsplit%atms(:))) then
          duplicate = .true.
          exit
        end if
      end do
      if (.not.duplicate) then
        allocate (tmpsplitqueue(nallnew))
        do ii = 1,nall
          tmpsplitqueue(ii) = splitqueue(ii)
        end do
        tmpsplitqueue(nallnew) = tmpsplit
        call move_alloc(tmpsplitqueue,splitqueue)
      end if
    end if

  end subroutine add_to_splitqueue

  recursive subroutine find_original_atoms(heap,targetlayer,targetnode,atoms)
    implicit none
    type(construct_heap),intent(in) :: heap
    integer,intent(in) :: targetlayer,targetnode
    integer,intent(out),allocatable :: atoms(:)
    integer :: this
    integer :: ii,jj,kk,nat,dim1,dim2
    integer,allocatable :: tmpatoms(:)

    associate (layer => heap%layer(targetlayer))
      nat = layer%node(targetnode)%nat
      allocate (atoms(nat),source=0)
      if (layer%parent == 0) then
        call position_mapping_reverse( &
          & layer%position_mapping,targetnode,nat,atoms)
      else
        ii = layer%parent
        jj = layer%parentnode
        call find_original_atoms(heap,ii,jj,tmpatoms)
        dim1 = size(tmpatoms,1)
        dim2 = size(layer%position_mapping,1)
        if (dim1 .ne. dim2) then
          stop "something went wrong in find_original_atoms()"
        end if
        do ii = 1,dim1
          kk = layer%position_mapping(ii,targetnode)
          do jj = 1,nat
            if (kk == jj) then
              atoms(jj) = tmpatoms(ii)
            end if
          end do
        end do
        deallocate (tmpatoms)
      end if
    end associate
  end subroutine find_original_atoms

  subroutine position_mapping_reverse(position_mapping,node,nat,revats)
    !* get the original (one layer up) atom positions for ati
    !* of a given node from the saved position_mapping  of the associated layer
    !* Invalid setups return 0
    implicit none
    integer,intent(in) :: position_mapping(:,:)
    integer,intent(in) :: node,nat
    integer,intent(out),allocatable :: revats(:)
    integer :: ii,jj,dim1,dim2
    dim1 = size(position_mapping,1)
    dim2 = size(position_mapping,2)
    allocate (revats(nat),source=0)
    if (node < 1.or.node > dim2) return
    do ii = 1,dim1
      jj = position_mapping(ii,node)
      if (jj > 0) revats(jj) = ii
    end do
  end subroutine position_mapping_reverse

  subroutine map_origins_for_layer(heap,targetlayer)
    implicit none
    class(construct_heap),intent(inout) :: heap
    integer,intent(in) :: targetlayer
    integer :: ii,jj
    logical,parameter :: debug = .false.
    if (targetlayer < 1.or.targetlayer > heap%nlayer) return
    if (.not.allocated(heap%layer(targetlayer)%node)) return
    associate (layer => heap%layer(targetlayer))
      if (allocated(layer%origin)) deallocate (layer%origin)
      allocate (layer%origin(layer%nnodes))
      do ii = 1,layer%nnodes
        layer%origin(ii)%natms = layer%node(ii)%nat
        call find_original_atoms(heap,targetlayer,ii,layer%origin(ii)%map)
      end do
    end associate
  end subroutine map_origins_for_layer

  subroutine find_current_position(heap,atom,targetlayer,targetnode,current)
    implicit none
    class(construct_heap),intent(in) :: heap
    integer,intent(in) :: atom,targetlayer,targetnode
    integer,intent(out) :: current
    integer :: ii,jj
    current = 0
    associate (layer => heap%layer(targetlayer))
      do ii = 1,layer%node(targetnode)%nat
        jj = layer%origin(targetnode)%map(ii)
        if (jj == atom) then
          current = ii
          exit
        end if
      end do
    end associate
  end subroutine find_current_position

  function count_endpoints(heap) result(nendpoints)
    implicit none
    class(construct_heap) :: heap
    integer :: nendpoints
    integer :: ii,jj,kk
    nendpoints = 0
    do ii = 1,heap%nlayer
      if (.not.allocated(heap%layer(ii)%childlayer)) then
        nendpoints = nendpoints+heap%layer(ii)%nnodes
      else
        do jj = 1,heap%layer(ii)%nnodes
          if (heap%layer(ii)%childlayer(jj) == 0) then
            nendpoints = nendpoints+1
          end if
        end do
      end if
    end do
  end function count_endpoints

  subroutine setup_queue(heap)
    implicit none
    class(construct_heap) :: heap
    integer :: nqueue,kk,ii,jj

    nqueue = heap%count_endpoints()
    heap%nqueue = nqueue
    allocate (heap%queue(nqueue))

    kk = 0
    do ii = 1,heap%nlayer
      if (.not.allocated(heap%layer(ii)%childlayer)) then
        do jj = 1,heap%layer(ii)%nnodes
          kk = kk+1
          heap%queue(kk)%layer = ii
          heap%queue(kk)%node = jj
        end do
      else
        do jj = 1,heap%layer(ii)%nnodes
          if (heap%layer(ii)%childlayer(jj) == 0) then
            kk = kk+1
            heap%queue(kk)%layer = ii
            heap%queue(kk)%node = jj
          end if
        end do
      end if
    end do

  end subroutine setup_queue

  subroutine fill_inverse_depth(heap)
    !*****************************************************
    !* Precompute inverse_depth for all layers by        *
    !* recursing from each root layer (parent == 0).     *
    !*****************************************************
    implicit none
    class(construct_heap),intent(inout) :: heap
    integer :: ii
    do ii = 1,heap%nlayer
      if (heap%layer(ii)%parent == 0) then
        call fill_inverse_depth_layer(heap,ii)
      end if
    end do
  end subroutine fill_inverse_depth

  recursive subroutine fill_inverse_depth_layer(heap,layerid)
    !*************************************************************
    !* Recursively assign inverse_depth to layer layerid and all *
    !* its descendants. Leaves get 1.0; internal layers get the  *
    !* sum of their children's inverse_depth values.             *
    !*************************************************************
    implicit none
    type(construct_heap),intent(inout) :: heap
    integer,intent(in) :: layerid
    integer :: jj,childid
    real(wp) :: s
    associate (lyr => heap%layer(layerid))
      if (.not.allocated(lyr%childlayer)) then
        lyr%inverse_depth = 1.0_wp
        return
      end if
      s = 0.0_wp
      do jj = 1,lyr%nnodes
        childid = lyr%childlayer(jj)
        if (childid == 0) then
          s = s+1.0_wp
        else
          call fill_inverse_depth_layer(heap,childid)
          s = s+heap%layer(childid)%inverse_depth
        end if
      end do
      lyr%inverse_depth = s
    end associate
  end subroutine fill_inverse_depth_layer

end module molbuilder_construct_list
