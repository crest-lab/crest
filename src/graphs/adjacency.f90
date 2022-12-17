!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022 Philipp Pracht
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

!> routines to help handling graph adjacency
module adjacency
  use iso_fortran_env,wp => real64
  use dijkstra_mod
  use floydwarshall_mod
  implicit none

  !============================================!
  !>--- interfacing of shortest path algorithms
  public :: Dijkstra,getPathD
  public :: FloydWarshall,getPathFW
  !============================================!

  public

!================================================================================!
!================================================================================!
contains  !> MODULE PROCEDURES START HERE
!================================================================================!
!================================================================================!

  subroutine wbo2adjacency(V,wbo,A,thr)
    implicit none
    !> INPUT
    integer,intent(in)  :: V
    real(wp),intent(in) :: wbo(V,V)
    real(wp),intent(in) :: thr
    !> OUTPUT
    integer,intent(out),allocatable :: A(:,:)
    integer :: i,j
    if (.not. allocated(A)) allocate (A(V,V))
    A = 0
    do i = 1,V
      do j = 1,i-1
        if (wbo(j,i) > thr) A(j,i) = min(nint(wbo(j,i)),1)
        A(i,j) = A(j,i)
      end do
    end do
  end subroutine wbo2adjacency

!================================================================================!

!>--- recursive routine to visit all vertices connected to a given one
  recursive subroutine visit_adjacent(node,V,A,tmp,k)
    implicit none
    !> INPUT
    integer,intent(in) :: node,V
    integer,intent(in) :: A(V,V)
    !> OUTPUT
    integer,intent(inout) :: k
    integer,intent(inout) :: tmp(V)
    !> LOCAL
    integer :: i,j
    if (any(tmp(:) == node)) return
    tmp(k) = node
    k = k + 1
    do i = 1,V
      if (A(node,i) > 0) then
        j = i
        call visit_adjacent(j,V,A,tmp,k)
      end if
    end do
  end subroutine visit_adjacent

!================================================================================!

!>--- set up tracking array frag: vertex->fragment
  subroutine setup_fragments(V,A,frag)
    implicit none
    !> INPUT
    integer,intent(in) :: V
    integer,intent(in) :: A(V,V)
    !> OUTPUT
    integer,intent(out) :: frag(V)
    !> LOCAL
    integer :: i,j,k,nfrag
    integer,allocatable :: tmp(:)
    frag = 0
    nfrag = 0
    allocate (tmp(V),source=0)
    !> sweeps
    do i = 1,V
      if (frag(i) == 0) then
        nfrag = nfrag + 1
        j = i
        k = 1
        tmp = 0
        call visit_adjacent(j,V,A,tmp,k)
        do j = 1,V
          if (any(tmp(:) == j)) frag(j) = nfrag
        end do
      end if
    end do
  end subroutine setup_fragments

!================================================================================!

  function check_adjacent(vertex1,vertex2,V,A,tmp) result(adjacent)
    implicit none
    !> INPUT
    integer,intent(in) :: vertex1,vertex2,V
    integer,intent(in) :: A(V,V)
    !> OUTPUT
    integer,intent(inout) :: tmp(V)
    logical :: adjacent
    !> LOCAL
    integer :: i,j,k
    adjacent = .false.
    !> check connection if an empty list was provided
    if (sum(tmp) == 0) then
      k = 1
      call visit_adjacent(vertex1,V,A,tmp,k)
    end if
    !> check the provided fragment
    if (any(tmp(:) == vertex1) .and. &
    &  any(tmp(:) == vertex2)) then
      adjacent = .true.
      return
    end if
  end function check_adjacent

!================================================================================!

  !> check all directly connected(!) vertex pairs to see if they share a ring
  subroutine check_rings_min(V,A,rings)
    implicit none
    !> INPUT
    integer,intent(in) :: V
    integer,intent(in) :: A(V,V)
    !> OUTPUT
    logical,intent(out) :: rings(V,V)
    !> LOCAL
    integer,allocatable  :: Adum(:,:)
    integer,allocatable  :: tmp(:)
    integer :: i,j
    rings = .false.
    allocate(Adum(V,V), source = 0)
    allocate(tmp(V), source = 0)
    !>-- delete the adjacency between two neighbours
    !    and then check if they are still connected via
    !    another path --> a ring must exist if so
    do i=1,V
      if(sum(A(:,i)) == 1) cycle    !> skip dangling vertices (e.g. H atoms)
      do j=1,i-1
         if(sum(A(:,j)) == 1) cycle !> skip dangling vertices
         if(A(i,j) > 0)then
           Adum(:,:) = A(:,:)
           Adum(i,j) = 0
           Adum(j,i) = 0
           tmp = 0
           rings(i,j) = check_adjacent(i,j,V,Adum,tmp) 
           rings(j,i) = rings(i,j)
         endif
      enddo
    enddo
    deallocate(tmp)
    deallocate(Adum)
  end subroutine check_rings_min

!================================================================================!

  !> get the smallest ring connecting adjacenent vertices M and N
  subroutine get_ring_min(V,A,M,N,path,nring)
    implicit none
    !> INPUT
    integer,intent(in) :: M,N
    integer,intent(in) :: V
    integer,intent(in) :: A(V,V)
    !> OUTPUT
    integer,intent(out) :: nring   !> number of members in ring
    integer,intent(out) :: path(V) !> vertices in ring
    !> LOCAL
    real(wp),allocatable  :: dist(:)
    integer,allocatable   :: tmp(:)
    integer,allocatable   :: Atmp(:,:)
    integer :: i,j

    path  = 0
    nring = 0
    allocate(Atmp(V,V), source=0)
    Atmp(:,:) = A(:,:)
    Atmp(M,N) = 0
    Atmp(N,M) = 0
    allocate(dist(V), source=0.0_wp)
    allocate(tmp(V), source=0)
    call Dijkstra(V,Atmp,M,dist,tmp)
    call getPathD(V,nring,M,N,tmp,path)
    deallocate(tmp,dist)
  end subroutine get_ring_min

!================================================================================!
!================================================================================!
end module adjacency
