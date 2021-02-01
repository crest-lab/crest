!=======================================================================!
! Taken from github.com/pprcht/shortestpaths
! under the MIT license and slightly modified
!=======================================================================!
!-- implementation of the algorithm including setup of "dist" and "prev"
subroutine FloydWarshall(V,nmat,emat,dist, prev)
    use iso_fortran_env, wp => real64
    implicit none
    integer :: V !number of vertices (atoms)
    integer :: nmat(V,V)    ! adjacency matrix
    real(wp) :: emat(V,V)   ! edge lenght matrix
    real(wp),intent(inout) :: dist(V,V)
    integer,intent(inout) :: prev(V,V)
    real(wp) :: inf
    real(wp) :: kdist
    integer :: i,j,k

    inf = 2*sum(emat)
    !-- set distances and previously visited nodes
    dist=inf
    prev=-1
    do i=1,V
        do j=1,V
            if(nmat(i,j)==1)then
               dist(i,j) = emat(i,j)
               prev(i,j) = j
            endif
        enddo
    enddo

    !-- The algorithm is based on the following assumption:
!   If a shortest path from vertex u to vertex v runns through a thrid
!   vertex w, then the paths u-to-w and w-to-v are already minimal.
!   Hence, the shorest paths are constructed by searching all path
!   that run over an additional intermediate point k
!   The following loop is the actual algorithm
    do k=1,V
        do i=1,V
            do j=1,V
                kdist = dist(i,k) + dist(k,j)
                !-- if the path ij runs over k, update
                if(dist(i,j) > kdist)then
                    dist(i,j) = kdist
                    prev(i,j) = prev(i,k)
                endif
            enddo
        enddo
    enddo
    
    return
end subroutine FloydWarshall

!-- reconstruct the path start to end for the Floyd-Warshall algorithm
subroutine getPathFW(nmax,prev,start,end,path,lpath)
    implicit none
    integer :: nmax
    integer :: start
    integer :: end
    integer :: prev(nmax,nmax)
    integer :: path(nmax)
    integer :: lpath
    integer :: i,k
    if(prev(start,end) == -1)then
        path(1)=-1
        return
    endif
    k=1
    path(k) = start
    i=start
    do while (  i .ne. end  )
        i = prev(i,end)
        k=k+1
        path(k) = i
    enddo
    lpath=k
    return
end subroutine getPathFW
