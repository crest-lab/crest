!=======================================================================!
! Taken from github.com/pprcht/shortestpaths
! under the MIT license and slightly modified
!=======================================================================!

module dijkstra_mod
  use iso_fortran_env,only:wp => real64,error_unit
  implicit none

  public :: getPathD
  public :: Dijkstra
  interface Dijkstra
    module procedure :: Dijkstra_classic
    module procedure :: Dijkstra_simple
  end interface Dijkstra 

!========================================================================================!
!========================================================================================!
contains
!========================================================================================!
!========================================================================================!

!> To run dijkstra's algorithm call the Dijkstra() routine providing
!> an adjacency matrix, distance matrix, and starting point,
!> and then run getPathD() to reconstruct the shortest path

!========================================================================================!
  subroutine Dijkstra_classic(V,A,E,start,dist,prev)
    implicit none
    !> INPUT
    integer,intent(in)  :: V       !> number of vertices
    integer,intent(in)  :: A(V,V)  !> adjacency matrix
    real(wp),intent(in) :: E(V,V)  !> vertex distance matrix
    integer,intent(in)  :: start   !> start vertex
    !> OUTPUT
    real(wp),intent(inout) :: dist(V) !> start -> vertex distances
    integer,intent(inout)  :: prev(V)
    !> LOCAL
    real(wp) :: inf
    real(wp) :: newdist
    integer,allocatable :: Q(:)
    integer :: Qcount
    integer :: i,j,u,vc

    inf = huge(inf)
    !>-- allocate array of unvisited vertices,
    !    distances and previously visited nodes
    allocate (Q(V))
    Qcount = V
    do i = 1,V
      Q(i) = i
      prev(i) = -1
      dist(i) = inf
    end do
    prev(start) = start
    dist(start) = 0.0_wp

    !>-- as long as there are unvisited vertices in Q
    do while (Qcount >= 1)
      !>-- get vertix with smallest reference value in dist
      j = minQ(Q,dist,Qcount)
      vc = Q(j)             !> track the vertex k
      Q(j) = Q(Qcount)      !> overwrite Q(j)
      Qcount = Qcount - 1   !> "resize" Q
      !>-- loop over the neighbours of vertex k
      do u = 1,V
        if (A(vc,u) == 0) cycle
        newdist = dist(vc) + E(vc,u)
        !>-- is the new dist value of u better than the previous one?
        if (newdist < dist(u)) then
          dist(u) = newdist
          prev(u) = vc
        end if
      end do
    end do
    deallocate (Q)
    return
  end subroutine Dijkstra_classic
  
  subroutine Dijkstra_simple(V,A,start,dist,prev)
  !> a variant that does not need the edge lenths E
    implicit none
    !> INPUT
    integer,intent(in)  :: V       !> number of vertices
    integer,intent(in)  :: A(V,V)  !> adjacency matrix
    integer,intent(in)  :: start   !> start vertex
    !> OUTPUT
    real(wp),intent(inout) :: dist(V) !> start -> vertex distances
    integer,intent(inout)  :: prev(V)
    !> LOCAL
    real(wp) :: inf
    real(wp) :: newdist
    integer,allocatable :: Q(:)
    integer :: Qcount
    integer :: i,j,u,vc

    inf = huge(inf)
    !>-- allocate array of unvisited vertices,
    !    distances and previously visited nodes
    allocate (Q(V))
    Qcount = V
    do i = 1,V
      Q(i) = i
      prev(i) = -1
      dist(i) = inf
    end do
    prev(start) = start
    dist(start) = 0.0_wp

    !>-- as long as there are unvisited vertices in Q
    do while (Qcount >= 1)
      !>-- get vertix with smallest reference value in dist
      j = minQ(Q,dist,Qcount)
      vc = Q(j)             !> track the vertex k
      Q(j) = Q(Qcount)      !> overwrite Q(j)
      Qcount = Qcount - 1   !> "resize" Q
      !>-- loop over the neighbours of vertex k
      do u = 1,V
        if (A(vc,u) == 0) cycle
        newdist = dist(vc) + float(A(vc,u))
        !>-- is the new dist value of u better than the previous one?
        if (newdist < dist(u)) then
          dist(u) = newdist
          prev(u) = vc
        end if
      end do
    end do
    deallocate (Q)
    return
  end subroutine Dijkstra_simple

  function minQ(Q,dist,Qcount)
    implicit none
    integer :: minQ
    integer :: Q(*)
    real(wp) :: dist(*)
    real(wp) :: dref,d
    integer :: Qcount
    integer :: i,j

    dref = huge(dref)
    do i = 1,Qcount
      j = Q(i)
      d = dist(j)
      if (d < dref) then
        minQ = i !> the position within Q has to be returned, not the vertex!
        dref = d
      end if
    end do
!    if(dref == huge(d))then
!        write(error_unit,*) "warning: disconnected graph!"
!    endif
    return
  end function minQ
!========================================================================================!

!> reconstruct path
  recursive subroutine getPathD(V,lpath,start,pos,prev,path)
    implicit none
    integer :: V
    integer :: lpath
    integer :: start
    integer :: pos
    integer :: prev(V)
    integer :: path(V)
    integer :: i,j,k
    lpath = lpath + 1
    if (pos == start) then !exit recursion
      path(lpath) = pos
      do i = 1,lpath / 2
        j = lpath - i + 1
        k = path(i)
        path(i) = path(j)
        path(j) = k
      end do
    else
      path(lpath) = pos
      j = prev(pos)
      call getPathD(V,lpath,start,j,prev,path)
    end if
    return
  end subroutine getPathD

!========================================================================================!
end module dijkstra_mod
