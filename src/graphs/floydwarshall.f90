!=======================================================================!
! Taken from github.com/pprcht/shortestpaths
! under the MIT license and slightly modified
!=======================================================================!

module floydwarshall_mod
  use iso_fortran_env,wp => real64
  implicit none

  public :: FloydWarshall,getPathFW
  interface FloydWarshall
    module procedure :: FloydWarshall_classic
    module procedure :: FloydWarshall_simple
  end interface FloydWarshall


!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

!> To run the Floyd-Warshall algorithm call the FloydWarshall() routine providing
!> an adjacency matrix, distance matrix, and starting point,
!> and then run getPathFW() to reconstruct the shortest path between two points.
!> Compared to Dijkstra's algorithm this one provides the shortest paths between
!> all vertices in one go. Though, for repeated executions it is more expensive.

!========================================================================================!

!>-- implementation of the algorithm including setup of "dist" and "prev"
  subroutine FloydWarshall_classic(V,A,E,dist,prev)
    implicit none
    !> INPUT
    integer,intent(in)  :: V         !> number of vertices (atoms)
    integer,intent(in)  :: A(V,V)    !> adjacency matrix
    real(wp),intent(in) :: E(V,V)    !> edge length/distance matrix
    !> OUTPUT
    real(wp),intent(inout) :: dist(V,V)
    integer,intent(inout) :: prev(V,V)
    !> LOCAL
    real(wp) :: inf
    real(wp) :: kdist
    integer :: i,j,k

    inf = 2 * sum(E)
    !>-- set distances and previously visited nodes
    dist = inf
    prev = -1
    do i = 1,V
      do j = 1,V
        if (A(i,j) == 1) then
          dist(i,j) = E(i,j)
          prev(i,j) = j
        end if
      end do
    end do

!>-- The algorithm is based on the following assumption:
!>   If a shortest path from vertex u to vertex v runns through a thrid
!>   vertex w, then the paths u-to-w and w-to-v are already minimal.
!>   Hence, the shorest paths are constructed by searching all path
!>   that run over an additional intermediate point k
!>   The following loop is the actual algorithm
    do k = 1,V
      do i = 1,V
        do j = 1,V
          kdist = dist(i,k) + dist(k,j)
          !>-- if the path ij runs over k, update
          if (dist(i,j) > kdist) then
            dist(i,j) = kdist
            prev(i,j) = prev(i,k)
          end if
        end do
      end do
    end do

    return
  end subroutine FloydWarshall_classic

!>-- implementation of the algorithm without edge length matrix
  subroutine FloydWarshall_simple(V,A,dist,prev)
    implicit none
    !> INPUT
    integer,intent(in)  :: V         !> number of vertices (atoms)
    integer,intent(in)  :: A(V,V)    !> adjacency matrix
    !> OUTPUT
    real(wp),intent(inout) :: dist(V,V)
    integer,intent(inout) :: prev(V,V)
    !> LOCAL
    real(wp) :: inf
    real(wp) :: kdist
    integer :: i,j,k

    inf = 2 * float(sum(A))
    !>-- set distances and previously visited nodes
    dist = inf
    prev = -1
    do i = 1,V
      do j = 1,V
        if (A(i,j) == 1) then
          dist(i,j) = float(A(i,j))
          prev(i,j) = j
        end if
      end do
    end do

!>-- The algorithm is based on the following assumption:
!>   If a shortest path from vertex u to vertex v runns through a thrid
!>   vertex w, then the paths u-to-w and w-to-v are already minimal.
!>   Hence, the shorest paths are constructed by searching all path
!>   that run over an additional intermediate point k
!>   The following loop is the actual algorithm
    do k = 1,V
      do i = 1,V
        do j = 1,V
          kdist = dist(i,k) + dist(k,j)
          !>-- if the path ij runs over k, update
          if (dist(i,j) > kdist) then
            dist(i,j) = kdist
            prev(i,j) = prev(i,k)
          end if
        end do
      end do
    end do

    return
  end subroutine FloydWarshall_simple


!========================================================================================!

!>-- reconstruct the path start to end for the Floyd-Warshall algorithm
  subroutine getPathFW(V,prev,start,end,path,lpath)
    implicit none
    integer :: V
    integer :: start
    integer :: end
    integer :: prev(V,V)
    integer :: path(V)
    integer :: lpath
    integer :: i,k
    if (prev(start,end) == -1) then
      path(1) = -1
      lpath = 0
      return
    end if
    k = 1
    path(k) = start
    i = start
    do while (i .ne. end)
      i = prev(i,end)
      k = k + 1
      path(k) = i
    end do
    lpath = k
    return
  end subroutine getPathFW

!========================================================================================!

end module floydwarshall_mod
