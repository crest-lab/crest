
module cregen_utils
!*****************************************
!* Module that implements a utility routines
!* mainly used in CREGEN
!*****************************************
  use crest_parameters
  use strucrd
  use adjacency
  use axis_module
  use internals_mod
  use utilities
  implicit none
  public

  real(wp),parameter :: bigval = huge(bigval)

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine cregen_calculate_fragments(mol,frag,nfrag)
!*******************************************************
!* Assign each atom in a molecule to a fragment
!* based on CN connectivity
!* The fragment assignment array and the total number
!* of fragments are both optional return arguments
!*******************************************************
    class(coord),intent(in) :: mol
    integer,intent(out),allocatable,optional :: frag(:)
    integer,intent(out),optional :: nfrag
    integer,allocatable :: tmp(:),A(:,:)
    real(wp),allocatable :: cn(:),bond(:,:)
    integer :: nat

    call mol%cn_to_bond(cn,bond)
    nat = mol%nat
    allocate (A(nat,nat),source=0)
    call wbo2adjacency(mol%nat,bond,A,0.01_wp)
    call setup_fragments(mol%nat,A,tmp)
    if (present(nfrag)) then
      nfrag = maxval(tmp,1)
    end if
    if (present(frag)) then
      call move_alloc(tmp,frag)
    end if
  end subroutine cregen_calculate_fragments

  logical function distcheck(n,xyz)
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    real(wp) :: rij(3)
    integer :: i,j
    distcheck = .true.
    do i = 1,n-1
      do j = i+1,n
        rij = xyz(:,j)-xyz(:,i)
        if (sum(rij*rij) .lt. 1.d-3) then
          distcheck = .false.
          return
        end if
      end do
    end do
    return
  end function distcheck

!================================================================================!
!================================================================================!
!> Simplified "topology"-related routines
!================================================================================!
!================================================================================!

  subroutine bondtotopo(nat,at,bond,cn,ntopo,topo,neighbourmat,excl)
    !********************************************************************
    !* generate the topo array for a given structure
    !* This includes some empirical hacks for use in CREGEN
    !********************************************************************
    integer,intent(in)  :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(inout) :: bond(nat,nat)
    real(wp),intent(in) :: cn(nat)
    integer,intent(inout)  :: ntopo
    integer,intent(inout),allocatable :: topo(:)
    real(wp),allocatable :: cn2(:)
    logical,intent(inout),allocatable :: neighbourmat(:,:)
    logical,intent(in),optional :: excl(nat)
    integer :: i,j,k,l
    integer :: icn
    real(wp) :: rcn

    ntopo = nat*(nat+1)/2
    if (.not.allocated(topo)) allocate (topo(ntopo),source=0)
    if (.not.allocated(neighbourmat)) allocate (neighbourmat(nat,nat))
    allocate (cn2(nat),source=0.0_wp)
    topo(1:ntopo) = 0
    neighbourmat(:,:) = .false.

    !--- some heuristic rules and CN array setup
    do i = 1,nat
      cn2(i) = cn(i)
      rcn = real(floor(cn(i)),wp)
      select case (at(i)) !additional empirical topology rules
        ! case( 5 ) !B
        !   if( nint(cn(i)) > 4) cn2(i)=4.0_wp
        ! case( 9,17,35,53 ) !F,Cl,Br,I
        !   cn2(i) = min(cn(i),1.0_wp)
      case (6) !C
        if ((cn(i)-rcn) < 0.7_wp) then
          cn2(i) = rcn
        end if
      end select
      !-- extreme CN cases
      if (nint(cn(i)) > 8) cn2(i) = 8.0_wp
      !empirical: rounding down up to .6 is better for topo setup
      if ((cn(i)-rcn) < 0.6_wp) then
        cn2(i) = rcn
      end if
    end do
    !--- build the topology
    do i = 1,nat
      icn = nint(cn2(i))
      do k = 1,icn
        j = maxloc(bond(:,i),1)
        bond(j,i) = 0.0d0
        if (i .eq. j) cycle
        neighbourmat(i,j) = .true. !--important: not automatically (i,j)=(j,i)
        if (present(excl)) then
          if (excl(i).or.excl(j)) neighbourmat(i,j) = .false.
        end if
      end do
    end do
    do i = 1,nat
      do j = 1,nat
        if (i == j) cycle
        l = lin(i,j)
        !-- only save matching topology --> prevent high CN failures
        if (neighbourmat(i,j).and.neighbourmat(j,i)) then
          topo(l) = 1
        else
          ! special case for carbon (because the carbon CN is typically correct)
          ! this helps, e.g. with eta-coordination in ferrocene
          ! (used, except if both are carbon)
          if (.not. (at(i) == 6.and.at(j) == 6)) then
            if (at(i) == 6.and.neighbourmat(i,j)) topo(l) = 1
            if (at(j) == 6.and.neighbourmat(j,i)) topo(l) = 1
          end if
        end if
      end do
    end do
    deallocate (cn2)
    return
  end subroutine bondtotopo

  subroutine nezcc(nat,at,xyz,cn,ntopo,topo,ncc)
    !***************************************************
    !* Check how many (potential) C=C bonds are present
    !* Expecting xyz in BOHR
    !***************************************************
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
    real(wp),parameter :: distcc = 1.384_wp*aatoau
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
    !* Expecting xyz in BOHR
    !********************************************************
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
    real(wp),parameter :: distcc = 1.384_wp*aatoau
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
      call DIHED2(xyz,a,b,c,d,winkel) !>-- from intmodes.f
      winkel = abs(winkel*(180.0_wp/pi))
      if (winkel > 180.0_wp) then
        winkel = 360.0_wp-winkel
      end if
      ezdihed(i) = winkel
    end do
    return
  end subroutine ezccdihed

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!
end module cregen_utils
