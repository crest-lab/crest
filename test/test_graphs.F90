module test_graphs
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use crest_parameters,only:wp
  use adjacency
  implicit none
  private

  public :: collect_graphs

  real(wp),parameter :: thr = 1.0e-2_wp

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for the graph routines in adjacency.f90
!========================================================================================!
!========================================================================================!

!> Collect all exported unit tests
  subroutine collect_graphs(testsuite)
    !> Collection of tests
    type(unittest_type),allocatable,intent(out) :: testsuite(:)

!&<
    testsuite = [ &
    new_unittest("wbo2adjacency thresholding    ",test_wbo2adjacency), &
    new_unittest("setup_fragments components    ",test_setup_fragments), &
    new_unittest("check_adjacent connectivity   ",test_check_adjacent), &
    new_unittest("check_rings_min ring edges    ",test_check_rings_min), &
    new_unittest("get_ring_min smallest ring    ",test_get_ring_min)  &
    ]
!&>
  end subroutine collect_graphs

!========================================================================================!

!> A linear chain 1-2-3-4 encoded as a WBO matrix. wbo2adjacency must produce
!> the symmetric 0/1 adjacency for bonds above the threshold only.
  subroutine test_wbo2adjacency(error)
    type(error_type),allocatable,intent(out) :: error
    integer,parameter :: V = 4
    real(wp) :: wbo(V,V)
    integer,allocatable :: A(:,:)
    integer :: i,j

    wbo = 0.0_wp
    !> chain bonds with bond orders ~1
    call setbond_r(wbo,1,2,0.98_wp)
    call setbond_r(wbo,2,3,1.05_wp)
    call setbond_r(wbo,3,4,0.97_wp)
    !> a weak through-space contribution that must be ignored
    call setbond_r(wbo,1,4,0.01_wp)

    call wbo2adjacency(V,wbo,A,0.02_wp)

    call check(error,allocated(A))
    if (allocated(error)) return
    call check(error,size(A,1),V)
    if (allocated(error)) return

    !> expected adjacency of the chain
    call check(error,A(1,2),1)
    call check(error,A(2,3),1)
    call check(error,A(3,4),1)
    if (allocated(error)) return
    !> non-bonds (incl. the sub-threshold 1-4 pair) must be zero
    call check(error,A(1,3),0)
    call check(error,A(2,4),0)
    call check(error,A(1,4),0)
    if (allocated(error)) return

    !> matrix must be symmetric with a zero diagonal
    do i = 1,V
      call check(error,A(i,i),0)
      if (allocated(error)) return
      do j = 1,V
        call check(error,A(i,j),A(j,i))
        if (allocated(error)) return
      end do
    end do
  end subroutine test_wbo2adjacency

!========================================================================================!

!> Two disconnected triangles {1,2,3} and {4,5,6}: setup_fragments must label
!> them as exactly two connected components.
  subroutine test_setup_fragments(error)
    type(error_type),allocatable,intent(out) :: error
    integer,parameter :: V = 6
    integer :: A(V,V)
    integer,allocatable :: frag(:)

    call two_triangles(A)
    call setup_fragments(V,A,frag)

    call check(error,allocated(frag))
    if (allocated(error)) return
    call check(error,maxval(frag),2)
    if (allocated(error)) return

    !> all members of a triangle share one label, the two triangles differ
    call check(error,frag(1),frag(2))
    call check(error,frag(2),frag(3))
    call check(error,frag(4),frag(5))
    call check(error,frag(5),frag(6))
    if (allocated(error)) return
    if (frag(1) == frag(4)) then
      call test_failed(error,"disconnected fragments share a label")
    end if
  end subroutine test_setup_fragments

!========================================================================================!

!> check_adjacent must report connectivity within a fragment and the lack of it
!> across fragments of the two-triangle graph.
  subroutine test_check_adjacent(error)
    type(error_type),allocatable,intent(out) :: error
    integer,parameter :: V = 6
    integer :: A(V,V)
    integer :: tmp(V)
    logical :: conn

    call two_triangles(A)

    !> 1 and 3 are in the same triangle -> connected
    tmp = 0
    conn = check_adjacent(1,3,V,A,tmp)
    call check(error,conn)
    if (allocated(error)) return

    !> 1 and 4 are in different triangles -> not connected
    tmp = 0
    conn = check_adjacent(1,4,V,A,tmp)
    call check(error,.not.conn)
  end subroutine test_check_adjacent

!========================================================================================!

!> A six-membered ring (1-2-3-4-5-6-1) with one dangling substituent (7 on 1).
!> Every ring bond must be flagged, the dangling bond must not.
  subroutine test_check_rings_min(error)
    type(error_type),allocatable,intent(out) :: error
    integer,parameter :: V = 7
    integer :: A(V,V)
    logical,allocatable :: rings(:,:)
    integer :: i

    call ring6_with_tail(A)
    call check_rings_min(V,A,rings)

    call check(error,allocated(rings))
    if (allocated(error)) return

    !> all six ring edges are part of a ring
    do i = 1,5
      call check(error,rings(i,i+1))
      if (allocated(error)) return
    end do
    call check(error,rings(6,1))
    if (allocated(error)) return

    !> the dangling bond 1-7 is not part of any ring
    call check(error,.not.rings(1,7))
  end subroutine test_check_rings_min

!========================================================================================!

!> get_ring_min must recover the six-membered ring through the edge 1-2.
  subroutine test_get_ring_min(error)
    type(error_type),allocatable,intent(out) :: error
    integer,parameter :: V = 7
    integer :: A(V,V)
    integer :: path(V),nring
    integer :: i

    call ring6_with_tail(A)
    call get_ring_min(V,A,1,2,path,nring)

    !> the smallest ring through 1-2 contains all six ring atoms
    call check(error,nring,6)
    if (allocated(error)) return

    !> the dangling atom 7 must not be on the ring path
    do i = 1,nring
      if (path(i) == 7) then
        call test_failed(error,"dangling atom included in ring path")
        return
      end if
    end do
  end subroutine test_get_ring_min

!========================================================================================!
!========================================================================================!
!> shared helpers for building small test graphs
!========================================================================================!
!========================================================================================!

!> set a symmetric real bond entry
  subroutine setbond_r(M,i,j,val)
    real(wp),intent(inout) :: M(:,:)
    integer,intent(in) :: i,j
    real(wp),intent(in) :: val
    M(i,j) = val
    M(j,i) = val
  end subroutine setbond_r

!> set a symmetric integer adjacency entry
  subroutine setbond_i(A,i,j)
    integer,intent(inout) :: A(:,:)
    integer,intent(in) :: i,j
    A(i,j) = 1
    A(j,i) = 1
  end subroutine setbond_i

!> two disconnected triangles {1,2,3} and {4,5,6}
  subroutine two_triangles(A)
    integer,intent(out) :: A(6,6)
    A = 0
    call setbond_i(A,1,2)
    call setbond_i(A,2,3)
    call setbond_i(A,3,1)
    call setbond_i(A,4,5)
    call setbond_i(A,5,6)
    call setbond_i(A,6,4)
  end subroutine two_triangles

!> six-membered ring 1-2-3-4-5-6-1 with a dangling atom 7 attached to 1
  subroutine ring6_with_tail(A)
    integer,intent(out) :: A(7,7)
    A = 0
    call setbond_i(A,1,2)
    call setbond_i(A,2,3)
    call setbond_i(A,3,4)
    call setbond_i(A,4,5)
    call setbond_i(A,5,6)
    call setbond_i(A,6,1)
    call setbond_i(A,1,7)
  end subroutine ring6_with_tail

!========================================================================================!
!========================================================================================!
end module test_graphs
