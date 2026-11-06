!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020 Philipp Pracht, Stefan Grimme
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

!> Internal utility routines for CCEGEN (PCA/k-means clustering of conformer ensembles)
!
module ccegen_utils
!*************************************************************
!* Internal utilities for the CCEGEN routine.                *
!* Contains PCA helpers, k-means clustering routines,        *
!* and statistical evaluation procedures.                    *
!*************************************************************
  implicit none
  private

  public :: clustleveval,PCA_grpwrite,excludeLight,excludeSelected,svd_to_pc, &
  &         eucdist,kmeans,kmeans_seeds,kmeans_assign,kmeans_recenter, &
  &         cluststat,statanal,statwarning,getdiederatoms,calc_dieders

! ══════════════════════════════════════════════════════════════════════════════
contains  !> MODULE PROCEDURES START HERE
! ══════════════════════════════════════════════════════════════════════════════

subroutine clustleveval(env,maxclust,csthr,SSRSSTthr,pcthr)
  !*********************************************************
  !* Set clustering level defaults for maxclust, csthr,    *
  !* SSRSSTthr, and pcthr based on the clustering level.   *
  !*********************************************************
  use crest_parameters,idp => dp
  use crest_data
  implicit none
  type(systemdata) :: env
  integer :: clev
  integer :: maxclust
  real(wp) :: csthr
  real(wp) :: SSRSSTthr
  real(wp) :: pcthr

  SSRSSTthr = 0.90d0

  clev = env%clustlev
  if (env%clustlev >= 10) then  !> incremental clustering mode
    clev = env%clustlev-10
  end if

  select case (clev)
  case (-1)  !-- loose
    maxclust = 25
    csthr = 0.80d0
    pcthr = 0.80d0
  case (1)   !-- tight
    maxclust = 400
    if (env%clustlev >= 10) maxclust = 50
    csthr = 0.85d0
    pcthr = 0.90d0
  case (2)   !-- vtight
    maxclust = 400
    if (env%clustlev >= 10) maxclust = 100
    csthr = 0.9d0
    pcthr = 0.95d0
    SSRSSTthr = 0.92d0
  case default  !-- normal
    maxclust = 100
    if (env%clustlev >= 10) maxclust = 25
    csthr = 0.80d0
    pcthr = 0.85d0
  end select

  return
end subroutine clustleveval

subroutine PCA_grpwrite(nclust,npc,mm,pcvec,member)
  !*************************************************************
  !* Write the first two principal component projections and   *
  !* cluster membership for each structure to cluster.order.   *
  !*************************************************************
  use crest_parameters,idp => dp
  implicit none
  integer,intent(in) :: nclust
  integer,intent(in) :: npc,mm
  real(wp),intent(in) :: pcvec(npc,mm)
  integer,intent(in) :: member(mm)
  integer :: ich,i
  open (newunit=ich,file='cluster.order')
  write (ich,'(4x,i0,4x,i0,4x,i0)') mm,nclust,npc
  if (npc > 1) then
    do i = 1,mm
      write (ich,'(i8,1x,f16.8,1x,f16.8,1x,i8)') i,pcvec(1,i),pcvec(2,i),member(i)
    end do
  else
    do i = 1,mm
      write (ich,'(i8,1x,f16.8,1x,i8)') i,pcvec(1,i),member(i)
    end do
  end if
  close (ich)
  return
end subroutine PCA_grpwrite

subroutine excludeLight(zmol,inc)
  !*******************************************
  !* Zero out hydrogen atoms in the inc      *
  !* array to exclude them from the PCA.     *
  !*******************************************
  use crest_parameters,idp => dp
  use zdata
  implicit none
  type(zmolecule) :: zmol
  integer :: inc(zmol%nat)
  integer :: i
  do i = 1,zmol%nat
    if (zmol%at(i) == 1) then
      inc(i) = 0
    end if
  end do
  return
end subroutine excludeLight

subroutine excludeSelected(zmol,inc,atlist)
  !*******************************************
  !* Zero out user-specified atoms in inc    *
  !* to exclude them from the PCA.           *
  !*******************************************
  use crest_parameters,idp => dp
  use zdata
  implicit none
  type(zmolecule) :: zmol
  integer :: inc(zmol%nat)
  character(len=*) :: atlist
  integer :: i,ncon
  integer,allocatable :: inc2(:)
  allocate (inc2(zmol%nat),source=0)
  call parse_atlist_new(atlist,ncon,zmol%nat,zmol%at,inc2)
  do i = 1,zmol%nat
    if (inc2(i) == 1) inc(i) = 0
  end do
  deallocate (inc2)
  return
end subroutine excludeSelected

subroutine svd_to_pc(measure,m,n,sig,U,pr)
  !*************************************************************
  !* Singular value decomposition via LAPACK DGEJSV.           *
  !* Factorises X = U * diag(sig) * V^T after mean-centering.  *
  !* The singular values (sig) are the principal components.   *
  !* Requires M >= N.                                          *
  !*                                                           *
  !* Input:  measure(n,m) - descriptor matrix                  *
  !* Output: sig(n)       - singular values (eigenvalues)      *
  !*         U(m,n)       - left singular vectors              *
  !*************************************************************
  use crest_parameters,idp => dp
  implicit none
  integer :: n,m
  real(wp) :: measure(n,m)
  real(wp) :: sig(n)
  real(wp) :: U(m,n)
  integer :: i,j,info,lwork
  real(wp),allocatable :: mean(:),tmp(:)
  integer,allocatable :: ind(:)
  real(wp),allocatable :: X(:,:),V(:,:),work(:)
  integer,allocatable :: iwork(:)
  logical :: pr
  if (pr) then
    write (stdout,*) m,' mesaurements'
    write (stdout,*) n,' props'
  end if
  allocate (mean(n),ind(m),tmp(m))
  lwork = max(2*M+N,6*N+2*N*N)
  allocate (X(m,n),V(n,n),iwork(m+3*n),work(lwork))
  mean = 0.0d0
  do i = 1,m
    do j = 1,n
      mean(j) = mean(j)+measure(j,i)
    end do
  end do
  mean = mean/float(m)
  if (pr) write (stdout,*) mean
  do i = 1,m
    do j = 1,n
      X(i,j) = (mean(j)-measure(j,i))
    end do
  end do
  if (pr) then
    call PRMAT(6,X,m,n,'X')
  end if
! ── LAPACKs' DGEJSV ──────────────────────────────────────────────────────────
  call DGEJSV('C','U','V','N','N','N',  &
 &              m,n,X,m,sig,U,m,V,n,    &
 &              WORK,LWORK,IWORK,INFO)
  if (pr) then
    write (stdout,*) info
    write (stdout,*) sig
    call PRMAT(6,U**2,M,N,'U')
    call PRMAT(6,V,N,N,'V')
  end if
  deallocate (work,iwork,V,X,tmp,ind,mean)
  return
end subroutine svd_to_pc

function eucdist(ndim,p,q) result(dist)
  !*******************************************
  !* Euclidean distance between points p     *
  !* and q in ndim-dimensional space.        *
  !*******************************************
  use crest_parameters,idp => dp
  implicit none
  real(ap) :: dist
  integer :: ndim
  real(ap) :: p(ndim)
  real(ap) :: q(ndim)
  integer :: i
  dist = 0.0d0
  do i = 1,ndim
    dist = dist+(q(i)-p(i))**2
  end do
  dist = sqrt(dist)
  return
end function eucdist

subroutine kmeans(nclust,npc,mm,centroid,pcvec,ndist,dist,member)
  !*************************************************************
  !* K-means clustering: iteratively assigns structures to the *
  !* nearest centroid and recenters until convergence or       *
  !* maxiter iterations are reached.                           *
  !*************************************************************
  use crest_parameters,idp => dp
  implicit none
  integer,intent(in) :: nclust
  integer,intent(in) :: npc,mm
  real(wp),intent(in) :: pcvec(npc,mm)
  integer(idp),intent(in) :: ndist
  real(sp),intent(in)   :: dist(ndist)
  integer,intent(inout) :: member(mm)
  real(ap),intent(inout):: centroid(npc,nclust)
  integer,allocatable :: refmember(:)
  integer :: iter
  integer,parameter :: maxiter = 300

  if (nclust .le. 1) return

  allocate (refmember(mm),source=0)

  call kmeans_seeds(nclust,npc,mm,centroid,pcvec,ndist,dist)

  do iter = 1,maxiter
    ! ── assign each structure to its nearest centroid ────────────────────────
    member = 0
    call kmeans_assign(nclust,npc,mm,centroid,pcvec,member)

    if (all(member == refmember)) then
      exit
    else
      refmember = member
    end if
    call kmeans_recenter(nclust,npc,mm,centroid,pcvec,member)
  end do

  deallocate (refmember)
  return
end subroutine kmeans

subroutine kmeans_seeds(nclust,npc,mm,centroid,pcvec,ndist,dist)
  !*************************************************************
  !* Initialise k-means centroids by finding maximally         *
  !* separated points in PC space (greedy farthest-first).     *
  !*************************************************************
  use crest_parameters,idp => dp
  use utilities
  implicit none
  integer :: nclust,npc,mm
  real(ap) :: centroid(npc,nclust)
  real(wp) :: pcvec(npc,mm)
  integer(idp) :: ndist
  real(sp) :: dist(ndist)
  real(sp) :: ddum
  integer(idp) :: k,kiter
  integer :: i,j,l,c
  real(wp) :: distsum,maxdistsum
  real(ap),allocatable :: p(:),q(:)
  integer,allocatable :: taken(:)

  ! ── seed 1 & 2: the two most distant structures ────────────────────────────
  ddum = 0.0_sp
  do kiter = 1,ndist
    if (dist(kiter) > ddum) then
      ddum = dist(kiter)
      k = kiter
    end if
  end do
  call revlin(k,j,i)  !> reverse of the lin index to get (i,j)

  centroid(1:npc,1) = pcvec(1:npc,i)
  centroid(1:npc,2) = pcvec(1:npc,j)

  if (nclust .le. 2) return

  ! ── seeds 3+: farthest point from all existing centroids ───────────────────
  allocate (p(npc),q(npc),taken(nclust))
  taken = 0
  taken(1) = i
  taken(2) = j
  do i = 3,nclust
    maxdistsum = 0.0d0
    c = 0
!$OMP PARALLEL PRIVATE ( l, q, p, j, distsum ) &
!$OMP SHARED ( i, centroid, npc, mm, maxdistsum, pcvec, c, taken )
!$OMP DO
    do j = 1,mm
      distsum = 0.0d0
      p(1:npc) = pcvec(1:npc,j)
      do l = 1,i-1
        q(1:npc) = centroid(1:npc,l)
        distsum = distsum+eucdist(npc,p,q)
      end do
      !$OMP CRITICAL
      if (distsum .gt. maxdistsum) then
        if (.not.any(taken == j)) then
          maxdistsum = distsum
          c = j
          taken(i) = c
        end if
      end if
      !$OMP END CRITICAL
    end do
!$OMP END DO
!$OMP END PARALLEL
    if (c == 0) then
      exit
    else
      centroid(1:npc,i) = pcvec(1:npc,c)
    end if
  end do

  deallocate (taken,q,p)

  return
end subroutine kmeans_seeds

subroutine kmeans_assign(nclust,npc,mm,centroid,pcvec,member)
  !*************************************************************
  !* Assign each structure to the nearest centroid.            *
  !*************************************************************
  use crest_parameters,idp => dp
  implicit none
  integer :: nclust,npc,mm
  real(ap) :: centroid(npc,nclust)
  real(wp) :: pcvec(npc,mm)
  integer :: member(mm)
  integer :: i,j,c
  real(ap),allocatable :: centdist(:)
  real(ap),allocatable :: p(:),q(:)

  allocate (centdist(nclust),source=0.0_ap)
  allocate (p(npc),q(npc))
!$OMP PARALLEL PRIVATE ( i, j, p, q, c, centdist ) &
!$OMP SHARED ( mm, nclust, member, centroid, npc, pcvec )
!$OMP DO
  do i = 1,mm
    p(1:npc) = pcvec(1:npc,i)
    do j = 1,nclust
      q(1:npc) = centroid(1:npc,j)
      centdist(j) = eucdist(npc,p,q)
    end do
    c = minloc(centdist,1)
    member(i) = c
  end do
!$OMP END DO
!$OMP END PARALLEL
  deallocate (q,p)
  deallocate (centdist)
  return
end subroutine kmeans_assign

subroutine kmeans_recenter(nclust,npc,mm,centroid,pcvec,member)
  !*************************************************************
  !* Recompute each centroid as the mean of its member         *
  !* structures in PC space.                                   *
  !*************************************************************
  use crest_parameters,idp => dp
  implicit none
  integer :: nclust,npc,mm
  real(ap) :: centroid(npc,nclust)
  real(wp) :: pcvec(npc,mm)
  integer :: member(mm)
  integer :: i,j,c
  real(ap),allocatable :: p(:),q(:)

  allocate (p(npc),q(npc))
  do i = 1,nclust
    c = 0
    p = 0.0d0
    do j = 1,mm
      if (member(j) == i) then
        c = c+1
        p(1:npc) = p(1:npc)+pcvec(1:npc,j)
      end if
    end do
    if (c > 0) then
      p = p/float(c)
      centroid(1:npc,i) = p(1:npc)
    else
      p = 999.9d0
    end if
  end do
  deallocate (q,p)

  return
end subroutine kmeans_recenter

subroutine cluststat(nclust,npc,mm,centroid,pcvec,member,DBI,pSF,SSRSST)
  !*************************************************************
  !* Compute cluster quality metrics for a given partition:    *
  !*   DBI    = Davies-Bouldin index (lower is better)         *
  !*   pSF    = pseudo-F statistic (higher is better)          *
  !*   SSRSST = SSR/SST ratio (higher means more distinct)     *
  !*************************************************************
  use crest_parameters,idp => dp
  implicit none
  integer,intent(in) :: nclust
  integer,intent(in) :: npc,mm
  real(wp),intent(in) :: pcvec(npc,mm)
  integer,intent(in) :: member(mm)
  real(ap),intent(in):: centroid(npc,nclust)
  real(wp),intent(out) :: DBI,pSF,SSRSST
  real(wp) :: SSE,SSR,SST
  real(ap),allocatable :: p(:),q(:)
  real(wp),allocatable :: compact(:)
  real(wp),allocatable :: DBmat(:,:)
  real(wp) :: d,Rij,maxDB,weight
  integer :: i,c,k,c2

  DBI = 0.0d0
  pSF = 0.0d0
  SSRSST = 0.0d0

  if (nclust < 2) return

  allocate (p(npc),q(npc))

  ! ── sum of squares error (within-cluster) ──────────────────────────────────
  SSE = 0.0d0
  do c = 1,nclust
    p(1:npc) = centroid(1:npc,c)
    do i = 1,mm
      if (member(i) == c) then
        q(1:npc) = pcvec(1:npc,i)
        d = eucdist(npc,p,q)
        SSE = SSE+d**2
      end if
    end do
  end do

  ! ── total sum of squares (weighted centroid as global mean) ────────────────
  SST = 0.0d0
  p = 0.0d0
  do c = 1,nclust
    weight = real(count(member(:) == c,1),wp)/real(mm,wp)
    p(1:npc) = p(1:npc)+centroid(1:npc,c)*weight
  end do
  do i = 1,mm
    q(1:npc) = pcvec(1:npc,i)
    d = eucdist(npc,p,q)
    SST = SST+d**2
  end do

  SSR = SST-SSE
  SSRSST = SSR/SST

  ! ── pseudo-F statistic ─────────────────────────────────────────────────────
  if (nclust > 1) then
    pSF = (SSR/(float(nclust)-1.0d0))
    if (mm == nclust) then
      pSF = 0.0d0
    else
      pSF = pSF/(SSE/(float(mm)-float(nclust)))
    end if
  else
    pSF = 0.0d0
  end if

  ! ── Davies-Bouldin index ───────────────────────────────────────────────────
  allocate (compact(nclust),source=0.0d0)
  do c = 1,nclust
    p(1:npc) = centroid(1:npc,c)
    k = 0
    do i = 1,mm
      if (member(i) == c) then
        k = k+1
        q(1:npc) = pcvec(1:npc,i)
        d = eucdist(npc,p,q)
        compact(c) = compact(c)+d
      end if
    end do
    if (k > 0) then
      compact(c) = compact(c)/float(k)
    else
      compact(c) = 0
    end if
  end do
  allocate (DBmat(nclust,nclust),source=0.0d0)
  do c = 1,nclust
    p(1:npc) = centroid(1:npc,c)
    do c2 = 1,nclust
      if (c2 == c) cycle
      q(1:npc) = centroid(1:npc,c2)
      d = eucdist(npc,p,q)
      Rij = (compact(c)+compact(c2))/d
      DBmat(c,c2) = Rij
    end do
  end do
  do c = 1,nclust
    maxDB = maxval(DBmat(:,c),1)
    DBI = DBI+maxDB
  end do
  DBI = DBI/float(nclust)
  deallocate (DBmat)
  deallocate (compact)

  deallocate (q,p)
  return
end subroutine cluststat

subroutine statanal(n,nmax,statistics,extrema,pr,clust_sizes)
  !*************************************************************
  !* Identify local minima of DBI and local maxima of pSF      *
  !* across the tested cluster counts.                         *
  !*************************************************************
  use crest_parameters
  implicit none
  integer :: n,nmax
  real(wp) :: statistics(3,nmax)
  logical,intent(inout) :: extrema(2,n)
  logical :: pr
  integer,intent(in),optional :: clust_sizes(n)
  real(wp) :: last,next,current
  integer :: i,csize

  extrema = .false.
! ── identify local extrema of the DBI and pSF ────────────────────────────────
  do i = 2,n-1
    last = statistics(1,i-1)
    next = statistics(1,i+1)
    current = statistics(1,i)
    if ((current < last).and.(current < next)) then
      extrema(1,i) = .true.
    end if
    last = statistics(2,i-1)
    next = statistics(2,i+1)
    current = statistics(2,i)
    if ((current > last).and.(current > next)) then
      extrema(2,i) = .true.
    end if
  end do
  !>-- boundary check: one-sided comparison for the last cluster count
  if (n >= 2) then
    if (statistics(1,n) < statistics(1,n-1)) extrema(1,n) = .true.
    if (statistics(2,n) > statistics(2,n-1)) extrema(2,n) = .true.
  end if

  if (pr) then
    write (stdout,*)
    write (stdout,'(1x,a,/)') 'Suggestions for cluster sizes:'
    do i = 1,n
      if (extrema(1,i).or.extrema(2,i)) then
        csize = i
        if (present(clust_sizes)) csize = clust_sizes(i)
        if (extrema(1,i).and.extrema(2,i)) then
          write (stdout,'(1x,i8,''*'',3x,a,f8.4)') csize,'SSR/SST',statistics(3,i)
        else
          write (stdout,'(1x,i8,4x,a,f8.4)') csize,'SSR/SST',statistics(3,i)
        end if
      end if
    end do
    write (stdout,'(/,1x,a)') 'Cluster counts marked with a star (*) are reasonable'
    write (stdout,'(1x,a)') 'suggestions according to BOTH the DBI and pSF.'
  end if

  return
end subroutine statanal

subroutine statwarning(fname)
  !*************************************************************
  !* Print a note about the arbitrary nature of clustering.    *
  !*************************************************************
  use crest_parameters
  implicit none
  character(len=*) :: fname
  write (stdout,*)
  write (stdout,'(1x,a)') '!---------------------------- NOTE ----------------------------!'
  write (stdout,'(2x,a)') 'The partitioning of data (the ensemble) into clusters'
  write (stdout,'(2x,a)') 'of similar characteristics (structures) is ARBITRARY'
  write (stdout,'(2x,a)') 'and depends on many criteria (e.g. choice of PCs).'
  write (stdout,'(2x,a)') 'The selected cluster count is the smallest reasonable'
  write (stdout,'(2x,a)') 'number of clusters that can be formed according to'
  write (stdout,'(2x,a)') 'the DBI and pSF values for the given data.'
  write (stdout,*)
  write (stdout,'(2x,a)') 'If other cluster sizes are desired, rerun CREST with'
  write (stdout,'(2x,3a)') '"crest --sort ',trim(fname),' --cluster <number of clusters>"'
  write (stdout,*)
  write (stdout,'(2x,a)') 'Other default evaluation settings can be chosen with the'
  write (stdout,'(2x,a)') 'keywords "loose","normal", and "tight" as <level> via'
  write (stdout,'(2x,3a)') '"crest --sort ',trim(fname),' --cluster <level>"'
  write (stdout,'(1x,a)') '!--------------------------------------------------------------!'
end subroutine statwarning

subroutine getdiederatoms(zmol,nat,inc,nb,diedat,ndied)
  !*************************************************************
  !* Extract dihedral angle atom quartets from the molecular   *
  !* topology, skipping terminal, ignored, and methyl atoms.   *
  !*************************************************************
  use crest_parameters,idp => dp
  use zdata
  use strucrd
  implicit none
  type(zmolecule) :: zmol
  integer :: nat
  integer :: inc(nat)  !> 1 = include, 0 = ignore
  integer :: nb
  integer :: diedat(4,nb)
  integer,intent(out) :: ndied
  integer :: a,b,c,d
  integer :: i,j,k

  ndied = 0
  do i = 1,nb
    a = zmol%bondpairs(1,i)
    b = zmol%bondpairs(2,i)
    if (inc(a) == 0) cycle          !> ignored by user
    if (inc(b) == 0) cycle          !> ignored by user
    if (zmol%zat(a)%nei == 1) cycle !> terminal atom
    if (zmol%zat(b)%nei == 1) cycle !> terminal atom
    if (zmol%methyl(a)) cycle       !> methyl carbon
    if (zmol%methyl(b)) cycle       !> methyl carbon
    ! ── get one neighbour of a and one of b to form the quartet ──────────────
    do j = 1,zmol%zat(a)%nei
      c = zmol%zat(a)%ngh(j)
      if (c == b) then
        cycle
      else
        exit
      end if
    end do
    do k = 1,zmol%zat(b)%nei
      d = zmol%zat(b)%ngh(k)
      if (d == a) then
        cycle
      else
        exit
      end if
    end do
    ndied = ndied+1
    !> quartet layout: (1)=neighbour of a, (2)=a, (3)=b, (4)=neighbour of b
    diedat(2,ndied) = a
    diedat(3,ndied) = b
    diedat(1,ndied) = c
    diedat(4,ndied) = d
  end do

  return
end subroutine getdiederatoms

subroutine calc_dieders(mol,ndied,diedat,diedr)
  !*****************************************************
  !* Calculate dihedral angles for selected atom       *
  !* quartets. Results are in radians (-pi, pi).       *
  !*****************************************************
  use crest_parameters,idp => dp
  use strucrd
  implicit none
  type(coord),intent(in) :: mol
  integer,intent(in) :: ndied
  integer,intent(in) :: diedat(4,ndied)
  real(wp),intent(out) :: diedr(ndied)
  integer :: i

  diedr = 0.0_wp
  do i = 1,ndied
    !> quartet: (1)=neighbour of a, (2)=a, (3)=b, (4)=neighbour of b
    diedr(i) = mol%dihedral(diedat(1,i),diedat(2,i), &
    &                       diedat(3,i),diedat(4,i))
  end do

  return
end subroutine calc_dieders

! ══════════════════════════════════════════════════════════════════════════════
! ══════════════════════════════════════════════════════════════════════════════
end module ccegen_utils
