
module irmsd_module
!*****************************************
!* Module that implements a more
!* modern interface to calculating RMSDs
!*****************************************
  use crest_parameters
  use strucrd
  use hungarian_module
  use axis_module
  implicit none
  private

  public :: rmsd,irmsd
  public :: min_rmsd

  public :: rmsd_align
  public :: checkranks,fallbackranks
  public :: setup_irmsd_ranks
  public :: molatomsort

  real(wp),parameter :: bigval = huge(bigval)

  public :: rmsd_core_cache
  type :: rmsd_core_cache
!*************************************
!* Memory cache for rmsd_core routine
!*************************************
    real(wp),allocatable :: x(:,:)
    real(wp),allocatable :: y(:,:)
    real(wp),allocatable :: xi(:)
    real(wp),allocatable :: yi(:)
    real(wp),allocatable :: xyzscratch(:,:,:)
  contains
    procedure :: allocate => allocate_rmsd_core_cache
  end type rmsd_core_cache

  public :: rmsd_cache
  type :: rmsd_cache
!****************************************************
!* cache implementation to avoid repeated allocation
!* and enable shared-memory parallelism
!****************************************************
    logical :: initialized = .false.
    real(wp),allocatable :: xyzscratch(:,:,:)
    integer,allocatable :: rank(:,:)
    integer,allocatable :: best_order(:,:)
    integer,allocatable :: current_order(:)
    integer,allocatable :: target_order(:)
    integer,allocatable :: order_bkup(:,:)
    integer,allocatable :: iwork(:)
    integer,allocatable :: iwork2(:,:)
    logical,allocatable :: assigned(:)  !> atom-wise
    logical,allocatable :: rassigned(:) !> rank-wise
    logical,allocatable :: lwork(:)

    integer :: nranks = 0
    integer,allocatable :: ngroup(:)
    logical :: stereocheck = .false.
    integer,allocatable :: proxy_topo_ref(:,:)
    integer,allocatable :: proxy_topo(:,:)

    type(rmsd_core_cache),allocatable :: ccache
    type(assignment_cache),allocatable :: acache
  contains
    procedure :: allocate => allocate_rmsd_cache
    procedure :: check_proxy_topo
    procedure :: initialize => initialize_rmsd_cache
  end type rmsd_cache

  real(wp),parameter :: inf = huge(1.0_wp)
  real(wp),parameter :: imat(3,3) = reshape([1.0_wp,0.0_wp,0.0_wp,  &
                                    &        0.0_wp,1.0_wp,0.0_wp,  &
                                    &        0.0_wp,0.0_wp,1.0_wp], &
                                    &        [3,3])

  real(wp),parameter :: Rx180(3,3) = reshape([1.0_wp,0.0_wp,0.0_wp,     &
                                     &          0.0_wp,-1.0_wp,0.0_wp,  &
                                     &          0.0_wp,0.0_wp,-1.0_wp], &
                                     &          [3,3])

  real(wp),parameter :: Ry180(3,3) = reshape([-1.0_wp,0.0_wp,0.0_wp,    &
                                     &          0.0_wp,1.0_wp,0.0_wp,   &
                                     &          0.0_wp,0.0_wp,-1.0_wp], &
                                     &          [3,3])

  !real(wp),parameter :: Rz180(3,3) = reshape([-1.0_wp,0.0_wp,0.0_wp,   &
  !                                   &          0.0_wp,-1.0_wp,0.0_wp, &
  !                                   &          0.0_wp,0.0_wp,1.0_wp], &
  !                                   &          [3,3])

  real(wp),parameter :: Rx90(3,3) = reshape([ &
                                     &    1.0_wp,0.0_wp,0.0_wp, &
                                     &    0.0_wp,0.0_wp,1.0_wp, &
                                     &    0.0_wp,-1.0_wp,0.0_wp &
                                     &    ], [3,3])
  real(wp),parameter :: Rx90T(3,3) = transpose(Rx90)

  real(wp),parameter :: Ry90(3,3) = reshape([ &
                                     &    0.0_wp,0.0_wp,-1.0_wp, &
                                     &    0.0_wp,1.0_wp,0.0_wp,  &
                                     &    1.0_wp,0.0_wp,0.0_wp   &
                                     &    ], [3,3])
  real(wp),parameter :: Ry90T(3,3) = transpose(Ry90)

  real(wp),parameter :: Rz90(3,3) = reshape([ &
                                     &    0.0_wp,1.0_wp,0.0_wp,  &
                                     &   -1.0_wp,0.0_wp,0.0_wp,  &
                                     &    0.0_wp,0.0_wp,1.0_wp   &
                                     &    ], [3,3])
  real(wp),parameter :: Rz90T(3,3) = transpose(Rz90)

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine allocate_rmsd_core_cache(self,nat,scratch)
    implicit none
    class(rmsd_core_cache),intent(inout) :: self
    integer,intent(in) :: nat
    logical,intent(in),optional :: scratch
    if (allocated(self%x)) deallocate (self%x)
    if (allocated(self%y)) deallocate (self%y)
    if (allocated(self%xi)) deallocate (self%xi)
    if (allocated(self%yi)) deallocate (self%yi)
    allocate (self%xi(nat),source=0.0_wp)
    allocate (self%yi(nat),source=0.0_wp)
    allocate (self%x(3,nat),source=0.0_wp)
    allocate (self%y(3,nat),source=0.0_wp)
    if(present(scratch))then
      if(scratch)then
        if(allocated(self%xyzscratch)) deallocate(self%xyzscratch)
        allocate(self%xyzscratch(3,nat,2),source=0.0_wp)
      endif
    endif
  end subroutine allocate_rmsd_core_cache

  subroutine allocate_rmsd_cache(self,nat)
    implicit none
    class(rmsd_cache),intent(inout) :: self
    integer,intent(in) :: nat
    if (allocated(self%xyzscratch)) deallocate (self%xyzscratch)
    if (allocated(self%rank)) deallocate (self%rank)
    if (allocated(self%best_order)) deallocate (self%best_order)
    if (allocated(self%current_order)) deallocate (self%current_order)
    if (allocated(self%target_order)) deallocate (self%target_order)
    if (allocated(self%order_bkup)) deallocate (self%order_bkup)
    if (allocated(self%iwork)) deallocate (self%iwork)
    if (allocated(self%iwork2)) deallocate (self%iwork2)
    if (allocated(self%assigned)) deallocate (self%assigned)
    if (allocated(self%rassigned)) deallocate (self%rassigned)
    if (allocated(self%lwork)) deallocate (self%lwork)
    if (allocated(self%ngroup)) deallocate (self%ngroup)
    if (allocated(self%proxy_topo_ref)) deallocate (self%proxy_topo_ref)
    if (allocated(self%proxy_topo)) deallocate (self%proxy_topo)
    if (allocated(self%ccache)) deallocate (self%ccache)
    if (allocated(self%acache)) deallocate (self%acache)
    allocate (self%assigned(nat),source=.false.)
    allocate (self%rassigned(nat),source=.false.)
    allocate (self%best_order(nat,3),source=0)
    allocate (self%current_order(nat),source=0)
    allocate (self%target_order(nat),source=0)
    allocate (self%order_bkup(nat,32),source=0)
    allocate (self%iwork(nat),source=0)
    allocate (self%iwork2(nat,2),source=0)
    allocate (self%rank(nat,2),source=0)
    self%nranks = 0
    allocate (self%ngroup(nat),source=0)
    allocate (self%proxy_topo(nat,2),source=0)
    allocate (self%proxy_topo_ref(nat,2),source=0)
    allocate (self%xyzscratch(3,nat,2),source=0.0_wp)
    allocate (self%ccache)
    allocate (self%acache)
    call self%ccache%allocate(nat)
    call self%acache%allocate(nat,nat,.true.) !> assume we are only using the LSAP implementation
  end subroutine allocate_rmsd_cache

  subroutine initialize_rmsd_cache(self,nat)
    implicit none
    class(rmsd_cache) :: self
    integer,intent(in) :: nat

    if (.not.self%initialized) then
      call self%allocate(nat)
      self%initialized = .true.
    else if (size(self%xyzscratch,2) .ne. nat) then
      call self%allocate(nat)
      self%initialized = .true.
    end if
  end subroutine initialize_rmsd_cache

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!

  function rmsd(ref,mol,mask,scratch,rotmat,gradient,ccache) result(rmsdval)
    !************************************************************************
    !* function rmsd
    !* Calculate the molecular RMSD via a quaternion algorithm
    !*
    !* Optional arguments are
    !*   mask - boolean array to select a substructure for RMSD calculation
    !*   scratch - workspace to create the substructures
    !*   rotmat  - rotation matrix as return argument
    !*   gradient - Cartesian gradient of the RMSD
    !************************************************************************
    implicit none
    real(wp) :: rmsdval
    type(coord),intent(in) :: ref
    type(coord),intent(in) :: mol
    !> OPTIONAL arguments
    logical,intent(in),optional :: mask(ref%nat)
    real(wp),intent(inout),target,optional :: scratch(3,ref%nat,2)
    real(wp),intent(out),optional :: rotmat(3,3)
    real(wp),intent(out),target,optional :: gradient(3,ref%nat)
    type(rmsd_core_cache),intent(inout),optional,target :: ccache
    !> variables
    type(rmsd_core_cache),allocatable,target :: ccachetmp
    type(rmsd_core_cache),pointer :: ccptr
    real(wp) :: Udum(3,3)
    real(wp),target :: gdum(3,3)
    integer :: nat,getrotmat
    logical :: calc_u
    real(wp),allocatable,target :: tmpscratch(:,:,:)
    logical :: getgrad
    real(wp),pointer :: grdptr(:,:)
    real(wp),pointer :: scratchptr(:,:,:)
    integer :: ic,k
    logical :: periodic
    real(wp) :: latuse(3,3)

    !> initialize to large value
    rmsdval = bigval
    !> check structure consistency
    if (mol%nat .ne. ref%nat) return

    !> periodic? if either structure carries a lattice, evaluate the RMSD under
    !> the minimum-image convention (no rotation fit). Prefer the probe's cell.
    periodic = .false.
    if (allocated(mol%lat)) then
      periodic = .true.
      latuse = mol%lat
    else if (allocated(ref%lat)) then
      periodic = .true.
      latuse = ref%lat
    end if

    !> get rotation matrix?
    getrotmat = 0
    calc_u = .false.
    if (present(rotmat)) then
      getrotmat = 1
      calc_u = .true.
    end if

    !> get gradient?
    if (present(gradient)) then
      getgrad = .true.
      gradient(:,:) = 0.0_wp
      grdptr => gradient
    else
      getgrad = .false.
      grdptr => gdum
    end if

    !> use present cache?
    if (present(ccache)) then
      ccptr => ccache
    else
      allocate (ccachetmp)
      call ccachetmp%allocate(ref%nat)
      ccptr => ccachetmp
    end if

!>--- substructure?
    if (present(mask)) then
      nat = count(mask(:))
      !> scratch workspace to use?
      if (present(scratch)) then
        scratchptr => scratch
      else if(allocated(ccptr%xyzscratch))then
        scratchptr => ccptr%xyzscratch
      else
        allocate (tmpscratch(3,nat,2))
        scratchptr => tmpscratch
      end if

      !> do the mapping
      k = 0
      do ic = 1,ref%nat
        if (mask(ic)) then
          k = k+1
          scratchptr(1:3,k,1) = mol%xyz(1:3,ic)
          scratchptr(1:3,k,2) = ref%xyz(1:3,ic)
        end if
      end do

      !> calculate
      if (periodic) then
        call rmsd_core(nat,scratchptr(1:3,1:nat,1),scratchptr(1:3,1:nat,2), &
        &          calc_u,Udum,rmsdval,getgrad,grdptr(1:3,:),ccptr,lat=latuse)
      else
        call rmsd_core(nat,scratchptr(1:3,1:nat,1),scratchptr(1:3,1:nat,2), &
        &          calc_u,Udum,rmsdval,getgrad,grdptr(1:3,:),ccptr)
      end if

      !> scatter the compact gradient (columns 1:nat) back onto the full atom
      !> order: the k-th selected atom sits at column k and must move to its
      !> true position ic. Iterate high->low so no value is overwritten before
      !> use; only clear a vacated source column (ic /= k). Unselected columns
      !> stay zero from the gradient(:,:) = 0 initialization above.
      if (getgrad) then
        k = nat
        do ic = ref%nat,1,-1
          if (mask(ic)) then
            grdptr(1:3,ic) = grdptr(1:3,k)
            if (ic /= k) grdptr(1:3,k) = 0.0_wp
            k = k-1
          end if
        end do
      end if

      nullify (scratchptr)
      if (allocated(tmpscratch)) deallocate (tmpscratch)

    else
!>--- standard calculation (no mask)
      if (periodic) then
        call rmsd_core(ref%nat,mol%xyz,ref%xyz, &
        &          calc_u,Udum,rmsdval,getgrad,grdptr(1:3,:),ccptr,lat=latuse)
      else
        call rmsd_core(ref%nat,mol%xyz,ref%xyz, &
        &          calc_u,Udum,rmsdval,getgrad,grdptr(1:3,:),ccptr)
      end if
    end if

    !> pass on rotation matrix if asked for
    if (calc_u) rotmat = Udum

  end function rmsd

!========================================================================================!

  subroutine rmsd_core(nat,xyz1,xyz2,calc_u,U,error,calc_g,grad,ccache,lat)
    !**********************************************************
    !* Rewrite or RMSD code with modified memory management
    !* Adapted from ls_rmsd, and using some of its subroutines
    !* The goal is to offload memory allocation to outside
    !* the routine in case it is repeadetly called
    !*
    !* If the optional lattice "lat" is given, the RMSD is
    !* evaluated under periodic boundary conditions: the probe
    !* (xyz1) is first minimum-image "unwrapped" relative to the
    !* reference (xyz2), x'_i = xyz2_i + MIC(xyz1_i - xyz2_i),
    !* which reassembles a molecule that wrapped around a cell
    !* face (and is a no-op for a compact molecule in a large
    !* box). The usual Kabsch superposition then removes BOTH
    !* translation and rotation, so the metric stays invariant to
    !* rigid-body motion (essential for the MTD conformer bias —
    !* otherwise the dynamics would just spin the molecule in the
    !* box). The integer image shift is piecewise-constant, so the
    !* gradient w.r.t xyz1 equals the standard Kabsch gradient.
    !**********************************************************
    use ls_rmsd,only:dstmev,rotation_matrix
    implicit none
    integer,intent(in) :: nat
    real(wp),intent(in) :: xyz1(3,nat)
    real(wp),intent(in) :: xyz2(3,nat)
    logical,intent(in) :: calc_u
    real(wp),dimension(3,3),intent(out) :: U
    real(wp),intent(out) :: error
    logical,intent(in) :: calc_g
    real(wp),intent(inout) :: grad(:,:)
    type(rmsd_core_cache),intent(inout) :: ccache
    real(wp),intent(in),optional :: lat(3,3) !> cell -> periodic (MIC) RMSD

    !> LOCAL
    integer :: i,j
    real(wp) :: x_center(3)
    real(wp) :: y_center(3)
    real(wp) :: x_norm,y_norm,lambda
    real(wp) :: Rmatrix(3,3)
    real(wp) :: S(4,4)
    real(wp) :: q(4)
    real(wp) :: tmp(3),rnat
    real(wp) :: latinv(3,3),sfrac(3)
    integer :: io

    !> associate
    associate (x => ccache%x,y => ccache%y,xi => ccache%xi,yi => ccache%yi)

      !> make copies of the original coordinates; under PBC the probe is
      !> minimum-image unwrapped onto the reference before the Kabsch fit
      if (present(lat)) then
        latinv = inv3x3(lat)
        do i = 1,nat
          sfrac = matmul(latinv,xyz1(1:3,i)-xyz2(1:3,i))
          sfrac = sfrac-anint(sfrac)
          x(1:3,i) = xyz2(1:3,i)+matmul(lat,sfrac)
        end do
      else
        x(1:3,1:nat) = xyz1(1:3,1:nat)
      end if
      y(1:3,1:nat) = xyz2(1:3,1:nat)

      !> calculate the barycenters, centroidal coordinates, and the norms
      x_norm = 0.0_wp
      y_norm = 0.0_wp
      rnat = 1.0_wp/real(nat,wp)
      do i = 1,3
        xi(:nat) = x(i,1:nat)
        yi(:nat) = y(i,1:nat)
        x_center(i) = sum(xi(1:nat))*rnat
        y_center(i) = sum(yi(1:nat))*rnat
        xi(1:nat) = xi(1:nat)-x_center(i)
        yi(1:nat) = yi(1:nat)-y_center(i)
        x(i,1:nat) = xi(1:nat)
        y(i,1:nat) = yi(1:nat)
        x_norm = x_norm+dot_product(xi,xi)
        y_norm = y_norm+dot_product(yi,yi)
      end do

      !> calculate the R matrix
      do i = 1,3
        do j = 1,3
          Rmatrix(i,j) = dot_product(x(i,1:nat),y(j,1:nat))
        end do
      end do

      !> S matrix
      S(1,1) = Rmatrix(1,1)+Rmatrix(2,2)+Rmatrix(3,3)
      S(2,1) = Rmatrix(2,3)-Rmatrix(3,2)
      S(3,1) = Rmatrix(3,1)-Rmatrix(1,3)
      S(4,1) = Rmatrix(1,2)-Rmatrix(2,1)

      S(1,2) = S(2,1)
      S(2,2) = Rmatrix(1,1)-Rmatrix(2,2)-Rmatrix(3,3)
      S(3,2) = Rmatrix(1,2)+Rmatrix(2,1)
      S(4,2) = Rmatrix(1,3)+Rmatrix(3,1)

      S(1,3) = S(3,1)
      S(2,3) = S(3,2)
      S(3,3) = -Rmatrix(1,1)+Rmatrix(2,2)-Rmatrix(3,3)
      S(4,3) = Rmatrix(2,3)+Rmatrix(3,2)

      S(1,4) = S(4,1)
      S(2,4) = S(4,2)
      S(3,4) = S(4,3)
      S(4,4) = -Rmatrix(1,1)-Rmatrix(2,2)+Rmatrix(3,3)

      !> Calculate eigenvalues and eigenvectors, and
      !> take the maximum eigenvalue lambda and the corresponding eigenvector q.
      call dstmev(S,lambda,q,io)
      if (io /= 0) then
        error = -1.0_wp
        return
      end if

      if (calc_u.or.calc_g) then
        !> reset
        U(:,:) = Imat(:,:)
        !> convert quaternion q to rotation matrix U
        call rotation_matrix(q,U)
      end if

      !> RMS Deviation
      error = sqrt(max(0.0_wp, ((x_norm+y_norm)-2.0_wp*lambda))*rnat)

      if (calc_g) then
        !> Gradient of the error of xyz1 w.r.t xyz2
        do i = 1,nat
          do j = 1,3
            tmp(:) = matmul(transpose(U(:,:)),y(:,i))
            grad(j,i) = ((x(j,i)-tmp(j))/error)*rnat
          end do
        end do
      end if

    end associate
  end subroutine rmsd_core

!========================================================================================!

  pure function inv3x3(a) result(ainv)
    !**********************************************************
    !* Analytic inverse of a 3x3 matrix (cofactor / adjugate).
    !**********************************************************
    implicit none
    real(wp),intent(in) :: a(3,3)
    real(wp) :: ainv(3,3)
    real(wp) :: det,detinv
    ainv(1,1) = a(2,2)*a(3,3)-a(2,3)*a(3,2)
    ainv(2,1) = a(2,3)*a(3,1)-a(2,1)*a(3,3)
    ainv(3,1) = a(2,1)*a(3,2)-a(2,2)*a(3,1)
    ainv(1,2) = a(1,3)*a(3,2)-a(1,2)*a(3,3)
    ainv(2,2) = a(1,1)*a(3,3)-a(1,3)*a(3,1)
    ainv(3,2) = a(1,2)*a(3,1)-a(1,1)*a(3,2)
    ainv(1,3) = a(1,2)*a(2,3)-a(1,3)*a(2,2)
    ainv(2,3) = a(1,3)*a(2,1)-a(1,1)*a(2,3)
    ainv(3,3) = a(1,1)*a(2,2)-a(1,2)*a(2,1)
    det = a(1,1)*ainv(1,1)+a(1,2)*ainv(2,1)+a(1,3)*ainv(3,1)
    detinv = 1.0_wp/det
    ainv = ainv*detinv
  end function inv3x3

!========================================================================================!
  subroutine rmsd_align(ref,mol,mask)
    !********************************************************
    !* Routine to align the structure "mol"
    !* to a reference structure "ref"
    !* via the quternion RMSD
    !*  Input:
    !*     ref - reference structure
    !*    mask - (optional) only consider specific atoms
    !*  In/Output:
    !*     mol - structure to have it's coordinates aligned
    !********************************************************

    implicit none
    type(coord),intent(in) :: ref
    type(coord),intent(inout) :: mol
    logical,intent(in),optional :: mask(:)

    real(wp) :: Umat(3,3),tmp
    real(wp) :: cref(3),cmol(3),shift(3)
    integer :: ii,jj,nn

    !> barycenter shift
    cref(:) = 0.0_wp
    cmol(:) = 0.0_wp
    shift(:) = 0.0_wp
    do ii = 1,ref%nat
      if (present(mask)) then
        if (mask(ii)) then
          cref(:) = cref(:)+ref%xyz(:,ii)
          cmol(:) = cmol(:)+mol%xyz(:,ii)
        end if
      else
        cref(:) = cref(:)+ref%xyz(:,ii)
        cmol(:) = cmol(:)+mol%xyz(:,ii)
      end if
    end do
    nn = ref%nat
    if (present(mask)) nn = count(mask)
    shift = cref-cmol
    shift = shift/nn
    do ii = 1,mol%nat
      mol%xyz(:,ii) = mol%xyz(:,ii)+shift(:)
    end do

    Umat(:,:) = 0.0_wp
    tmp = rmsd(ref,mol,mask=mask,rotmat=Umat)

    mol%xyz = matmul(Umat,mol%xyz)

  end subroutine rmsd_align

!========================================================================================!

  subroutine min_rmsd(ref,mol,rcache,rmsdout,align,topocheck,io)
    !****************************************************************************
    !* Main routine to determine minium RMSD considering atom permutation
    !* Input
    !*   ref  - the reference structure
    !*   mol  - the structure to be matched to ref
    !* Optinal arguments
    !*   rcache    - memory cache
    !*   rmsdout   - the calculated RMSD scalar
    !*   align     - quarternion-align mol in the last stage
    !*   topocheck - check molecule topology? if absent, doing check is default
    !*   io        - return status
    !****************************************************************************
    implicit none
    !> IN & OUTPUT
    type(coord),intent(in) :: ref
    type(coord),intent(inout) :: mol
    type(rmsd_cache),intent(inout),optional,target :: rcache
    real(wp),intent(out),optional :: rmsdout
    logical,intent(in),optional :: align
    logical,intent(in),optional :: topocheck
    integer,intent(out),optional :: io

    !> LOCAL
    type(rmsd_cache),pointer :: cptr
    type(rmsd_cache),allocatable,target :: local_rcache
    integer :: nat,ii,rnk,dumpunit,uniquenesscase,ioloc
    integer :: nunique
    real(wp) :: calc_rmsd
    real(wp) :: tmprmsd_sym(32)
    real(wp) :: rotmat(3,3),rotconst(3),shift(3)
    logical :: topocheck_l = .true.
    logical,parameter :: debug = .false.

!>--- defaults
    ioloc = 0

!>--- Initialization
    if (present(rcache)) then
      cptr => rcache
    else
      write (stdout,*) "WARNING: No iRMSD-cache provided. Attempting to fall back to atom types for sorting ranks."
      allocate (local_rcache)
      if (ref%nat .ne. mol%nat) then
        error stop 'Unequal molecule size in min_rmsd()'
      end if
      nat = max(ref%nat,mol%nat)
      call local_rcache%allocate(nat)
      call fallbackranks(ref,mol,nat,local_rcache%rank)
      cptr => local_rcache
    end if
    if (present(topocheck)) then
      topocheck_l = topocheck
    end if
    cptr%nranks = maxval(cptr%rank(:,1))

!>-- Consistency (topology) check
    if (topocheck_l) then
      ioloc = cptr%check_proxy_topo(ref,mol)
      if (ioloc > 0) then
        write (stdout,'(1x,a)') "WARNING: Different atom topologies detected in min_rmsd(), can't restore an atom order!"
        if (present(rmsdout)) then
          if (ioloc > 2) then !> topo check identified at least the same system size and maxrank --> quaternion RMSD may be feasible
            write (stdout,'(10x,a)') "Falling back to quaternion RMSD without reordering atoms. Values may be nonsensical."
            rmsdout = rmsd(ref,mol,ccache=cptr%ccache)
          else
            rmsdout = huge(rmsdout)
          end if
        end if
        if (present(io)) io = ioloc
        return
      end if
    end if

!>--- First sorting, to at least restore rank order (only if that's not the case!)
    if (.not.all(cptr%rank(:,1) .eq. cptr%rank(:,2))) then
      call rank_2_order(ref%nat,cptr%rank(:,1),cptr%target_order)
      call rank_2_order(mol%nat,cptr%rank(:,2),cptr%current_order)
      if (debug) then
        write (*,*) 'current order & rank & target order'
        do ii = 1,mol%nat
          write (*,*) cptr%current_order(ii),cptr%rank(ii,2),cptr%target_order(ii)
        end do
      end if
      call molatomsort(mol,mol%nat,cptr%current_order,cptr%target_order,cptr%iwork)
      cptr%rank(:,2) = cptr%rank(:,1) !> since the ranks must be equal now!
      if (debug) then
        write (*,*) 'sorted order & rank'
        do ii = 1,mol%nat
          write (*,*) cptr%current_order(ii),cptr%rank(ii,2)
        end do
      end if
    end if

!>--- Count symmetry equivalent groups and assign all unique atoms immediately
!     Note, the rank can be zero if we only are looking at heavy atoms
    if (all(cptr%ngroup(:) .eq. 0)) then
      do ii = 1,ref%nat
        rnk = cptr%rank(ii,1)
        if (rnk > 0) then
          cptr%ngroup(rnk) = cptr%ngroup(rnk)+1
        end if
      end do
    end if
    !> assignment reset
    cptr%assigned(:) = .false.
    cptr%rassigned(:) = .false.
    cptr%rassigned(cptr%nranks+1:) = .true. !> skip unneeded allocation space
    do ii = 1,ref%nat
      cptr%iwork(ii) = ii         !> also init iwork
      cptr%target_order(ii) = ii  !> also init target_order
      rnk = cptr%rank(ii,2)
      if (rnk < 1) then
        cptr%assigned(ii) = .true.
        cycle
      end if
      if (cptr%ngroup(rnk) .eq. 1) then
        cptr%assigned(ii) = .true.
        cptr%rassigned(rnk) = .true.
      end if
    end do
    if (debug) then
      write (*,*) 'rank & # members'
      do ii = 1,mol%nat
        if (cptr%ngroup(ii) > 0) then
          write (*,*) ii,cptr%ngroup(ii)
        end if
      end do
    end if

!>--- Perform the desired symmetry operations, align with rotational axis, run LSAP algo
!>    Since the rotational axis alignment can be a bit arbitrary w.r.t 180° rotations
!>    we need to check these as well.
    if (debug) then
      open (newunit=dumpunit,file='debugirmsd.xyz')
      call ref%append(dumpunit)
    end if

!>--- Check how many indices are unique
    cptr%lwork = unique_rank_mask(cptr%rank(:,1))
    nunique = count(cptr%lwork)

!> --------------------------------------------------------
!> SUBSTRUCTURE-BASED ALIGNMENT with enough unique indices
!> --------------------------------------------------------

    !> The logic here is: if we have enough unique atoms
    !> we can align the molecule with them and identify
    !> symmetry equivalent atoms via LSAP in those thereafter
    IF (nunique >= 3) then
      !> mol still needs a first alignment and CMA shift
      call CMAtrf(mol%nat,mol%nat,mol%at,mol%xyz)

      tmprmsd_sym(:) = inf
      tmprmsd_sym(1) = rmsd(ref,mol,cptr%lwork, &
        &                   cptr%xyzscratch,rotmat=rotmat, &
        &                   ccache=cptr%ccache)
      mol%xyz = matmul(rotmat,mol%xyz)
      call min_rmsd_iterate_through_groups(ref,mol,cptr,tmprmsd_sym(1))
      cptr%order_bkup(:,1) = cptr%iwork(:)
      if (cptr%stereocheck) then
        mol%xyz(3,:) = -mol%xyz(3,:)

        tmprmsd_sym(2) = rmsd(ref,mol,cptr%lwork, &
          &                   cptr%xyzscratch,rotmat=rotmat, &
          &                   ccache=cptr%ccache)
        mol%xyz = matmul(rotmat,mol%xyz)
        call min_rmsd_iterate_through_groups(ref,mol,cptr,tmprmsd_sym(2))
        cptr%order_bkup(:,2) = cptr%iwork(:)
        mol%xyz(3,:) = -mol%xyz(3,:)
      end if

      ii = minloc(tmprmsd_sym,1)
      if (ii == 2) then
        !> if the non-mirrored check was lower, revert the mirroring
        mol%xyz(3,:) = -mol%xyz(3,:)
      end if
!> ----------------------------------------------------
    ELSE
!> ----------------------------------------------------
!> ROTATIONAL AXIS ALIGNMENT AND LSAP CHECKS - START
!> ----------------------------------------------------

      !> initialize to huge
      tmprmsd_sym(:) = inf
      !> initial alignment of mol
      call axis(mol%nat,mol%at,mol%xyz,rotconst)
      call min_rmsd_rotcheck_unique(rotconst,uniquenesscase)

      !> Running the checks and check of uniqueness of rotational axes
      call min_rmsd_rotcheck_permute(ref,mol,cptr,tmprmsd_sym,1,uniquenesscase)
      if (debug) then
        write (*,*) 'Total LSAP cost:',minval(tmprmsd_sym(1:16))
        call mol%append(dumpunit)
      end if

      !> mirror z and re-run the same checks (i.e. the false rotamer inversion)
      if (cptr%stereocheck) then
        mol%xyz(3,:) = -mol%xyz(3,:)  !> mirror z
        call axis(mol%nat,mol%at,mol%xyz) !> align

        !> Running the checks
        call min_rmsd_rotcheck_permute(ref,mol,cptr,tmprmsd_sym,2,uniquenesscase)
        if (debug) then
          write (*,*) 'Total LSAP cost (inverted):',minval(tmprmsd_sym(17:32))
          call mol%append(dumpunit)
        end if
        mol%xyz(3,:) = -mol%xyz(3,:)  !> restore z
      end if

!>--- select the best match among the ones after symmetry operations and use its ordering
      ii = minloc(tmprmsd_sym(1:32),1)
      if (debug) then
        write (*,*) 'final alignment:',ii,"/ 32"
      end if
      if (ii > 16) then
        mol%xyz(3,:) = -mol%xyz(3,:)
        if (debug) write (*,*) 'inverting'
      end if
      if ((ii > 4.and.ii < 9).or.(ii > 20.and.ii < 25)) then
        if (uniquenesscase == 1) mol%xyz = matmul(Rx90,mol%xyz)
        if (uniquenesscase == 2) mol%xyz = matmul(Rz90,mol%xyz)
        if (uniquenesscase == 3) mol%xyz = matmul(Rz90,mol%xyz)
        if (debug) write (*,*) '90° tilt'
      else if ((ii > 8.and.ii < 13).or.(ii > 24.and.ii < 29)) then
        mol%xyz = matmul(Ry90,mol%xyz)
      else if ((ii > 12.and.ii < 17).or.(ii > 28)) then
        mol%xyz = matmul(Rx90,mol%xyz)
      end if
      select case (ii) !> 180° rotations
      case (1,5,9,13,17,21,25,29)
        continue
      case (2,6,10,14,18,22,26,30)
        mol%xyz = matmul(Rx180,mol%xyz)
        if (debug) write (*,*) '180°x'
      case (3,7,11,15,19,23,27,31)
        mol%xyz = matmul(Rx180,mol%xyz)
        mol%xyz = matmul(Ry180,mol%xyz)
        if (debug) write (*,*) '180°x, 180°y'
      case (4,8,12,16,20,24,28,32)
        mol%xyz = matmul(Ry180,mol%xyz)
        if (debug) write (*,*) '180°y'
      end select
!> ----------------------------------------------------
!> rotational axis alignment and LSAP checks - END
!> ----------------------------------------------------
    END IF
    cptr%current_order(:) = cptr%order_bkup(:,ii)
!> ----------------------------------------------------

    if (debug) then
      write (*,*) 'Determined remapping'
      do ii = 1,mol%nat
        write (*,*) cptr%current_order(ii),'-->',cptr%target_order(ii)
      end do
    end if

    call molatomsort(mol,mol%nat,cptr%current_order,cptr%target_order,cptr%iwork)
    if (debug) then
      call mol%append(dumpunit)
      close (dumpunit)
    end if

!>--- final RMSD with fully restored atom order
    if (present(align)) then
      calc_rmsd = rmsd(ref,mol,scratch=cptr%xyzscratch,ccache=cptr%ccache,rotmat=rotmat)
      if (align) then
        mol%xyz = matmul(rotmat,mol%xyz)
      end if
    else
      calc_rmsd = rmsd(ref,mol,scratch=cptr%xyzscratch,ccache=cptr%ccache)
    end if

    if (present(rmsdout)) rmsdout = calc_rmsd
    if (present(io)) io = ioloc
  end subroutine min_rmsd

  function irmsd(ref,mol,rcache, &
  & iinversion,align,topocheck,allcanon,io) result(rmsdval)
!***************************************************************
!* irmsd function
!* Standalone implementation to compare two structures
!* with the iRMSD method analog to the rmsd function
!* the optional rcache will get allocated if it not already is.
!*
!* Externally supplied per-atom ranks (ref%id / mol%id, allocated only
!* when a full valid set was provided) are used directly via
!* setup_irmsd_ranks; otherwise the canonical ranking is recomputed.
!***************************************************************
    implicit none

    type(coord),intent(inout) :: mol,ref
    integer,intent(in),optional :: iinversion
    type(rmsd_cache),intent(inout),optional,target :: rcache
    logical,intent(in),optional :: align
    logical,intent(in),optional :: topocheck
    logical,intent(in),optional :: allcanon
    integer,intent(out),optional :: io
    real(wp) :: rmsdval
    !> LOCAL
    type(rmsd_cache),pointer :: cptr
    type(rmsd_cache),allocatable,target :: local_rcache
    logical :: align_l = .true.
    logical :: topocheck_l = .true.
    logical :: allcanon_l = .false.
    integer :: io_l = 0
    real(wp) ::tmpd(3),tmpdist
    integer :: i,ich
    logical :: mirror
    logical,parameter :: debug = .false.

    !> move ref to CMA and align rotational axes
    call axis(ref%nat,ref%at,ref%xyz)

    !> optional args
    if (present(align)) align_l = align
    if (present(topocheck)) topocheck_l = topocheck
    if (present(allcanon)) allcanon_l = allcanon
    if (present(rcache)) then
      cptr => rcache
    else
      allocate (local_rcache)
      cptr => local_rcache
    end if
    call cptr%initialize(ref%nat)
    !call rcache%allocate(ref%nat)

    !> determine the per-atom ranks (provided ids or recomputed) and the
    !> false-enantiomer flag, then apply the inversion override
    call setup_irmsd_ranks(ref,mol,cptr%rank,cptr%stereocheck,allcanon=allcanon_l)
    if (present(iinversion)) then
      select case (iinversion)
      case (0)
        mirror = .true.
      case (1)
        mirror = .true.
        cptr%stereocheck = .true.
      case (2)
        mirror = .false.
        cptr%stereocheck = .false.
      end select
    end if

    call min_rmsd(ref,mol,rcache=cptr,rmsdout=rmsdval, &
    &    align=align_l,topocheck=topocheck_l,io=io_l)

    if (present(io)) io = io_l

    return
  end function irmsd

!========================================================================================!

  subroutine setup_irmsd_ranks(ref,mol,ranks,stereocheck,allcanon)
!**************************************************************
!* Determine the per-atom ranks used for iRMSD matching of the
!* ref/mol pair, plus the false-enantiomer (stereocheck) flag.
!*
!* Externally supplied ranks (coord%id, allocated only when a full
!* valid set was provided) are used directly; any side lacking them
!* is canonicalized on the fly. If the resulting pair fails the
!* consistency check, any side that came from provided ids is
!* recomputed and the check repeated; if it still fails, the ranks
!* fall back to atom types.
!*
!* ref,mol     : the two structures (may carry %id)
!* ranks       : (out) (nat,2) ranks for ref (col 1) and mol (col 2)
!* stereocheck : (out) .true. if false enantiomers should be considered
!* allcanon    : (optional) if .true., skip canonical ranks entirely and
!*               use atom types (legacy 'allcanon' behavior)
!**************************************************************
    use canonical_mod
    implicit none
    type(coord),intent(in) :: ref,mol
    integer,intent(out) :: ranks(:,:)
    logical,intent(out) :: stereocheck
    logical,intent(in),optional :: allcanon
    !> LOCAL
    type(canonical_sorter) :: canref,canmol
    integer,allocatable :: r1(:),r2(:)
    logical :: allcanon_l,ref_provided,mol_provided,refstereo

    allcanon_l = .false.
    if (present(allcanon)) allcanon_l = allcanon

    if (allcanon_l) then
      call fallbackranks(ref,mol,ref%nat,ranks)
      stereocheck = .false.
      return
    end if

    ref_provided = allocated(ref%id)
    mol_provided = allocated(mol%id)

    !> reference ranks (and the stereo flag derived from them)
    allocate (r1(ref%nat))
    if (ref_provided) then
      r1(:) = ref%id(:)
      canref%rank = ref%id              !> seed ranks; hasstereo self-fills the graph
      refstereo = canref%hasstereo(ref)
      call canref%deallocate()
    else
      call canref%init(ref,invtype='apsp+',heavy=.false.)
      r1(:) = canref%rank(:)
      refstereo = canref%hasstereo(ref)
      call canref%shrink()
    end if

    !> mobile ranks
    allocate (r2(mol%nat))
    if (mol_provided) then
      r2(:) = mol%id(:)
    else
      call canmol%init(mol,invtype='apsp+',heavy=.false.)
      r2(:) = canmol%rank(:)
      call canmol%shrink()
    end if

    !> if provided ranks do not reconcile, recompute the provided side(s)
    if (.not.checkranks(ref%nat,r1,r2).and. (ref_provided.or.mol_provided)) then
      if (ref_provided) then
        call canref%deallocate()
        call canref%init(ref,invtype='apsp+',heavy=.false.)
        r1(:) = canref%rank(:)
        refstereo = canref%hasstereo(ref)
        call canref%shrink()
      end if
      if (mol_provided) then
        call canmol%deallocate()
        call canmol%init(mol,invtype='apsp+',heavy=.false.)
        r2(:) = canmol%rank(:)
        call canmol%shrink()
      end if
    end if

    !> use the (reconciled) ranks, else fall back to atom types
    if (checkranks(ref%nat,r1,r2)) then
      ranks(:,1) = r1(:)
      ranks(:,2) = r2(:)
    else
      call fallbackranks(ref,mol,ref%nat,ranks)
    end if
    stereocheck = .not.refstereo

    deallocate (r1,r2)
  end subroutine setup_irmsd_ranks

!========================================================================================!

  subroutine min_rmsd_iterate_through_groups(ref,mol,rcache,val)
    implicit none
    type(coord),intent(in) :: ref
    type(coord),intent(inout) :: mol
    type(rmsd_cache),intent(inout),target :: rcache
    real(wp),intent(out) :: val
    integer :: rr,ii
    real(wp) :: val0
    type(assignment_cache),pointer :: aptr
    logical,parameter :: debug = .false.

    !> reset val
    val = 0.0_wp

    if (debug) then
      write (*,*) '# ranks:',rcache%nranks
    end if
    aptr => rcache%acache
    do rr = 1,rcache%nranks
      if (rcache%rassigned(rr)) cycle

      !> LSAP wrapper that computes the relevant Cost matrix for the atoms of rank rr
      call compute_linear_sum_assignment( &
      &        ref,mol,rcache%rank,rcache%ngroup,rr, &
      &        rcache%iwork2,aptr,val0)

      do ii = 1,rcache%ngroup(rr)
        rcache%iwork(rcache%iwork2(ii,1)) = rcache%iwork2(ii,2)
      end do

      !> add up the total LSAP cost (of considered ranks)
      !> we need this if we have to decide on a mapping in case of false enantiomers
      val = val+val0
    end do

  end subroutine min_rmsd_iterate_through_groups

!========================================================================================!

  subroutine min_rmsd_rotcheck_unique(rot,uniquenesscase,thr)
    !*******************************************************
    !* Based on the rotational constants, determine what we
    !* need to do with the molecule in the following
    !*******************************************************
    implicit none
    real(wp),intent(in) :: rot(3)
    integer,intent(out) :: uniquenesscase
    real(wp),intent(in),optional :: thr
    logical :: unique(3)
    integer :: nunique

    uniquenesscase = 0
    call uniqueax(rot,unique,thr)

    nunique = count(unique,1)
    select case (nunique)
    case (3) !> 3 unique principal axes
      uniquenesscase = 0
    case (1) !> one unique principal axis
      if (unique(1)) uniquenesscase = 1 !> A unique (long axis)
      if (unique(3)) uniquenesscase = 2 !> C unique (short axis)
    case (0) !> rotationally ambiguous system
      uniquenesscase = 3
    end select
  end subroutine min_rmsd_rotcheck_unique

!=======================================================================================!

  subroutine min_rmsd_rotcheck_permute(ref,mol,cptr,values,step,uniquenesscase)
    implicit none
    type(coord),intent(in) :: ref
    type(coord),intent(inout) :: mol
    type(rmsd_cache),intent(inout),target :: cptr
    real(wp),intent(inout) :: values(:)
    integer,intent(in) :: step,uniquenesscase
    integer :: ii,debugunit2
    real(wp) :: vals(16),dum
    logical,parameter :: debug = .false.

    !> reset val
    vals(:) = inf

    if (debug) then
      open (newunit=debugunit2,file='rotdebug.xyz')
      call ref%append(debugunit2)
    end if

    ALIGNLOOP: do ii = 1,4
      call min_rmsd_iterate_through_groups(ref,mol,cptr,dum)
      vals(1+4*(ii-1)) = dum
      if (debug) call mol%append(debugunit2)
      cptr%order_bkup(:,1+4*(ii-1)+16*(step-1)) = cptr%iwork(:)

      mol%xyz = matmul(Rx180,mol%xyz)
      call min_rmsd_iterate_through_groups(ref,mol,cptr,dum)
      vals(2+4*(ii-1)) = dum
      if (debug) call mol%append(debugunit2)
      cptr%order_bkup(:,2+4*(ii-1)+16*(step-1)) = cptr%iwork(:)

      mol%xyz = matmul(Ry180,mol%xyz)
      call min_rmsd_iterate_through_groups(ref,mol,cptr,dum)
      vals(3+4*(ii-1)) = dum
      if (debug) call mol%append(debugunit2)
      cptr%order_bkup(:,3+4*(ii-1)+16*(step-1)) = cptr%iwork(:)

      mol%xyz = matmul(Rx180,mol%xyz)
      call min_rmsd_iterate_through_groups(ref,mol,cptr,dum)
      vals(4+4*(ii-1)) = dum
      if (debug) call mol%append(debugunit2)
      cptr%order_bkup(:,4+4*(ii-1)+16*(step-1)) = cptr%iwork(:)

      mol%xyz = matmul(Ry180,mol%xyz) !> restore

      !exit ALIGNLOOP
      select case (uniquenesscase)
      case (0) !> 3 Unique moments of inertia
        exit ALIGNLOOP
      case (1) !> only one unique moment of inertia (A)
        if (ii == 2) then
          mol%xyz = matmul(Rx90T,mol%xyz)
          exit ALIGNLOOP
        end if
        mol%xyz = matmul(Rx90,mol%xyz)
      case (2) !> only one unique moment of inertia (C)
        if (ii == 2) then
          mol%xyz = matmul(Rz90T,mol%xyz)
          exit ALIGNLOOP
        end if
        mol%xyz = matmul(Rz90,mol%xyz)
      case (3)
        if (ii == 1) then
          mol%xyz = matmul(Rz90,mol%xyz)
        else if (ii == 2) then
          mol%xyz = matmul(Rz90T,mol%xyz)
          mol%xyz = matmul(Ry90,mol%xyz)
        else if (ii == 3) then
          mol%xyz = matmul(Ry90T,mol%xyz)
          mol%xyz = matmul(Rx90,mol%xyz)
        else
          mol%xyz = matmul(Rx90T,mol%xyz)
          exit ALIGNLOOP
        end if
      end select

    end do ALIGNLOOP

    if (debug) then
      close (debugunit2)
      write (*,*) 'vals:',vals(:)
    end if

    do ii = 1,16
      values(ii+16*(step-1)) = vals(ii)
    end do
  end subroutine min_rmsd_rotcheck_permute

!========================================================================================!

  subroutine fallbackranks(ref,mol,nat,ranks)
    !*****************************************************************
    !* If we are doing ranks on-the-fly (i.e. without canonical algo)
    !* we can fall back to just using the atom types
    !*****************************************************************
    implicit none
    type(coord),intent(in) :: ref,mol
    integer,intent(in) :: nat
    integer,intent(inout) :: ranks(nat,2)

    integer,allocatable :: typemap(:),rtypemap(:)
    integer :: k,ii
    allocate (typemap(nat),source=0)
    k = 0
    do ii = 1,ref%nat
      if (.not.any(typemap(:) .eq. ref%at(ii))) then
        k = k+1
        typemap(k) = ref%at(ii)
      end if
    end do
    do ii = 1,mol%nat
      if (.not.any(typemap(:) .eq. mol%at(ii))) then
        k = k+1
        typemap(k) = mol%at(ii)
      end if
    end do
    k = maxval(typemap(:))
    allocate (rtypemap(k),source=0)
    do ii = 1,nat
      if (typemap(ii) == 0) cycle
      rtypemap(typemap(ii)) = ii
    end do
    !> assign
    do ii = 1,ref%nat
      ranks(ii,1) = rtypemap(ref%at(ii))
    end do
    do ii = 1,mol%nat
      ranks(ii,2) = rtypemap(mol%at(ii))
    end do
    deallocate (rtypemap)
    deallocate (typemap)
  end subroutine fallbackranks

!========================================================================================!

  subroutine compute_linear_sum_assignment(ref,mol,ranks, &
                        & ngroups,targetrank,iwork2,acache,val0)
    !**************************************************************
    !* Run the linear assignment algorithm on the desired subset
    !* of atoms (via rank and targetrank)
    !**************************************************************
    implicit none
    !> IN & OUTPUT
    type(coord),intent(in) :: ref
    type(coord),intent(inout) :: mol
    integer,intent(in) :: ranks(:,:)
    integer,intent(in) :: ngroups(:)
    integer,intent(in) :: targetrank
    integer,intent(inout) :: iwork2(:,:)
    type(assignment_cache),intent(inout),optional,target :: acache
    real(wp),intent(out) :: val0

    !> LOCAL
    type(assignment_cache),pointer :: aptr
    type(assignment_cache),allocatable,target :: local_acache
    integer :: nat,i,j,ii,jj,rnknat,iostatus
    real(sp) :: dists(3)

    logical,parameter :: debug = .false.

    val0 = 0.0_wp

    if (present(acache)) then
      aptr => acache
    else
      allocate (local_acache)
      if (ref%nat .ne. mol%nat) then
        error stop 'Unequal molecule size in compute_linear_sum_assignment()'
      end if
      nat = max(ref%nat,mol%nat)
      call local_acache%allocate(nat,nat,.true.)
      aptr => local_acache
    end if

    !> Compute the cost matrix, which is simply the distance matrix
    !> between the two molecules.
    !> To avoid computational overhead we can skip the square root.
    !> It won't affect the result
    !> Also, since aptr%Cost is a flattened matrix, we only fill
    !> the first rnknat**2 entries
    rnknat = ngroups(targetrank)
    ii = 0
    do i = 1,ref%nat
      if (ranks(i,1) .ne. targetrank) cycle
      ii = ii+1
      iwork2(ii,1) = i !> mapping using the first column of iwork2
      jj = 0
      do j = 1,mol%nat
        if (ranks(j,2) .ne. targetrank) cycle
        jj = jj+1
        dists(:) = real((ref%xyz(:,i)-mol%xyz(:,j))**2,sp) !> use i and j
        aptr%Cost(jj+(ii-1)*rnknat) = sum(dists)
      end do
    end do

    if (debug) then
      write (*,*) 'target rank',targetrank,'# atoms',rnknat
    end if

    call lsap(aptr,rnknat,rnknat,.false.,iostatus)

    !> paasing back the determined order as second column of iwork2
    if (iostatus == 0) then
      if (debug) then
        do i = 1,rnknat
          write (*,*) iwork2(aptr%a(i),1),'-->',iwork2(aptr%b(i),1)
        end do
      end if
      do i = 1,rnknat
        jj = aptr%a(i)
        ii = aptr%b(i)
        if (ii == -1.or.jj == -1) cycle  !> cycle bad assignments
        val0 = val0+aptr%Cost(jj+(ii-1)*rnknat)
        iwork2(i,2) = iwork2(aptr%b(i),1)
      end do
    else
      !> in the unlikely case we have a failure of the LSAP
      !> we do just a 1:1 mapping, just so that the algo doesn't crash
      iwork2(1:rnknat,2) = iwork2(1:rnknat,1)
    end if

  end subroutine compute_linear_sum_assignment

!========================================================================================!

  subroutine rank_2_order(nat,rank,order)
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: rank(nat)
    integer,intent(out) :: order(nat)
    integer :: ii,jj,k,maxrank
    order(:) = 0
    maxrank = maxval(rank(:))
    k = 0
    do ii = 1,maxrank
      do jj = 1,nat
        if (rank(jj) == ii) then
          k = k+1
          order(jj) = k
        end if
      end do
    end do
  end subroutine rank_2_order

!========================================================================================!

  function checkranks(nat,ranks1,ranks2) result(yesno)
    !***********************************************************************
    !* Check two rank arrays to see if we have the same amount of
    !* atoms in the same ranks (a condition to bein able to work with them)
    !***********************************************************************
    implicit none
    logical :: yesno
    integer,intent(in) :: nat
    integer,intent(in) :: ranks1(nat)
    integer,intent(in) :: ranks2(nat)
    integer :: ii,jj,maxrank1,maxrank2
    integer :: count1,count2
    yesno = .false.

    maxrank1 = maxval(ranks1)
    maxrank2 = maxval(ranks2)
    !> different maxranks, so we can't have the same and return
    if (maxrank1 .ne. maxrank2) return

    do ii = 1,maxrank1
      count1 = 0
      count2 = 0
      do jj = 1,nat
        if (ranks1(jj) .eq. ii) count1 = count1+1
        if (ranks2(jj) .eq. ii) count2 = count2+1
      end do
      !> not the same amount of atoms in rank ii, return from function
      if (count1 .ne. count2) return
    end do

    !> if we reach this point we can assume the given ranks are o.k.
    yesno = .true.
  end function checkranks

!========================================================================================!
  function unique_rank_mask(ranks) result(mask)
    !*********************************************
    !* Takes a rank array and creates a mask that
    !* contains .true. if the respective rank
    !* appears only a single time
    !*********************************************
    implicit none
    integer,intent(in) :: ranks(:)
    logical,allocatable :: mask(:)
    integer :: ii,jj,n,k,rii
    n = size(ranks,1)
    allocate (mask(n),source=.false.)
    do ii = 1,n
      rii = ranks(ii)
      k = count(ranks == rii)
      if (k == 1) mask(ii) = .true.
    end do
  end function unique_rank_mask

!========================================================================================!

  subroutine molatomsort(mol,n,current_order,target_order,index_map)
    implicit none
    type(coord),intent(inout) :: mol
    integer,intent(in) :: n
    integer,intent(inout) :: current_order(n)
    integer,intent(in) :: target_order(n)
    integer,intent(inout) :: index_map(n)
    integer :: i,j,correct_atom,current_position

    !> Step 1: Create a mapping from target_order to current_order positions
    do i = 1,n
      index_map(current_order(i)) = i
    end do

    !> Step 2: Restore the target order
    do i = 1,n
      correct_atom = target_order(i)
      current_position = index_map(correct_atom)

      if (i /= current_position) then
        !> Swap atoms i and current_position in molecule
        call mol%swap(i,current_position)

        !> Update the index map since the atoms have been swapped
        index_map(current_order(i)) = current_position
        index_map(current_order(current_position)) = i

        !> Update the current_order array to reflect the swap
        j = current_order(i)
        current_order(i) = current_order(current_position)
        current_order(current_position) = j
      end if
    end do
  end subroutine molatomsort

!==========================================================================================!

  function check_proxy_topo(self,ref,mol) result(io)
    !******************************************************
    !* Attempt to compare the "topology" for the molecules ref and mol
    !* Assumes that ranks have been computed already.
    !* Checks are in order (cheap to expensive):
    !*   1) system size
    !*   2) max rank
    !*   3) joint sorted ranks and atom types
    !* Returns "io" with value 0 if successfull, or a number indatinc faliure condition
    implicit none
    class(rmsd_cache) :: self
    type(coord),intent(in) :: ref
    type(coord),intent(in) :: mol
    integer :: io
    integer :: n1,n2,m1,m2

    io = 0

    !> Check 1
    n1 = ref%nat
    n2 = mol%nat
    if (n1 .ne. n2) then
      io = 1
      return
    end if

    !> Check 2
    m1 = maxval(self%rank(:,1))
    m2 = maxval(self%rank(:,2))
    if (m1 .ne. m2) then
      io = 2; return
    end if

    !> Check 3
    self%proxy_topo_ref(:,1) = ref%at(:)
    self%proxy_topo_ref(:,2) = self%rank(:,1)
    call qsortm(self%proxy_topo_ref,2,self%iwork)

    self%proxy_topo(:,1) = mol%at(:)
    self%proxy_topo(:,2) = self%rank(:,2)
    call qsortm(self%proxy_topo,2,self%iwork)
    if (.not.all(self%proxy_topo .eq. self%proxy_topo_ref)) then
      io = 3
      return !> some difference in the sorting, return before setting passing to true
    end if

    !> All checks passed, io should still be 0
  end function check_proxy_topo

  recursive subroutine qsorti(v,ix,l,r)
    !*********************
    !* idx'ed quicksort
    !*********************
    integer,intent(in) :: v(:)
    integer,intent(inout) :: ix(:)
    integer,intent(in) :: l,r
    integer :: i,j,p,t,n
    if (l >= r) return
    p = v(ix((l+r)/2))
    n = size(v,1)
    i = l; j = r
    do
      do while (v(ix(i)) < p); i = i+1; end do
      do while (v(ix(j)) > p); j = j-1; end do
      if (i <= j) then
        t = ix(i); ix(i) = ix(j); ix(j) = t
        i = min(i+1,n); j = max(j-1,1)
      else
        exit
      end if
    end do
    if (l < j) call qsorti(v,ix,l,j)
    if (i < r) call qsorti(v,ix,i,r)
  end subroutine qsorti

  subroutine qsortm(a,k,ix)
    !************************************
    !* matrix wrapper to qsorti
    !* order is reflected to all columns
    !************************************
    integer,intent(inout) :: a(:,:)
    integer,intent(in) :: k
    integer :: n,i
    integer,intent(inout) :: ix(size(a,1))
    n = size(a,1)
    do i = 1,n
      ix(i) = i
    end do
    call qsorti(a(:,k),ix,1,n)
    a = a(ix,:)
  end subroutine qsortm

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!
end module irmsd_module
