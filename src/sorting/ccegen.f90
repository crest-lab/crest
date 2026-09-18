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

subroutine CCEGEN(env,pr,fname)
  !*****************************************************************************
  !* PCA-based clustering of a conformational ensemble.                        *
  !* Performs SVD to extract principal components, then partitions structures  *
  !* into representative clusters via k-means.                                 *
  !*                                                                           *
  !* Input:  env   - system data (method settings, thresholds)                 *
  !*         pr    - printout flag                                             *
  !*         fname - ensemble file name                                        *
  !*****************************************************************************
  use ccegen_utils
  use crest_parameters,idp => dp
  use crest_data
  use zdata
  use strucrd
  use utilities
  implicit none
  type(systemdata) :: env
  type(timer) :: ctimer
  logical,intent(in)   :: pr
  character(len=*),intent(in) :: fname
  type(zmolecule) :: zmol
  type(zequal) :: groups
  type(zequal) :: subgroups
  integer,allocatable :: inc(:)
  logical :: heavyonly
  integer :: i,j,k,l,ich,c
  integer :: nat,nall
  real(wp) :: dum,dum2
  type(coord),allocatable :: mols(:)

  character(len=:),allocatable :: measuretype

  !>--- SVD params
  integer :: ntaken
  integer :: nallnew
  real(wp),allocatable :: xyznew(:,:,:)
  real(wp),allocatable :: measure(:,:)
  integer :: mn,mm
  real(wp),allocatable :: pc(:)
  real(wp),allocatable :: pcvec(:,:)
  real(wp),allocatable :: pcdum(:,:)
  integer :: nbnd,ndied
  integer,allocatable :: diedat(:,:)
  real(wp),allocatable :: diedr(:)
  real(wp) :: pcsum
  real(wp) :: pcthr
  real(wp) :: pcmin
  integer :: pccap
  integer :: npc
  real(wp),allocatable :: geo(:,:)
  integer,allocatable :: na(:),nb(:),nc(:)

  !>--- clustering params
  character(len=:),allocatable :: clusteralgo
  integer :: nclust
  integer :: nclustiter
  integer :: nclustmin,nclustmax
  integer,allocatable :: member(:)
  real(ap),allocatable :: p(:),q(:)
  real(sp),allocatable :: dist(:)
  real(ap),allocatable :: centroid(:,:)
  integer(idp) :: ndist,klong
  real(wp) :: DBI,pSF,SSRSST,SSRSSTthr
  real(wp) :: csthr
  integer :: ncb,ancb
  real(wp),allocatable :: eclust(:)
  integer,allocatable :: clustbest(:),ind(:)
  real(wp),allocatable :: statistics(:,:)
  integer,allocatable :: clust_sizes(:)
  logical,allocatable :: extrema(:,:)
  logical :: autolimit
  real(wp) :: fraclimit

  real(wp) :: emin,erel

  call ctimer%init(20)
  if (pr) then
    call largehead('Principal Component Analysis (PCA) and Clustering')
    write (stdout,'(1x,a,a)') 'Input file: ',trim(fname)
  end if

! ── set threads ─────────────────────────────────────────────────────────────────
  call cregen_setthreads(stdout,env,pr)

! ── read ensemble ────────────────────────────────────────────────────────────────
  call rdensemble(fname,nall,mols)
  if (nall < 1) then
    error stop "Ensemble is empty! must stop"
  end if
  nat = mols(1)%nat
  if (nall == 1) then
    if (pr) then
      write (stdout,*) 'Only one structure in ensemble!'
      write (stdout,*) 'Write structure to ',clusterfile,' and skip PCA parts'
    end if
    open (newunit=ich,file=clusterfile)
    call mols(1)%append(ich)
    close (ich)
    deallocate (mols)
    return
  end if

  heavyonly = .true.
  measuretype = env%pcmeasure
  clusteralgo = 'kmeans'
  pcthr = env%pcthr
  pcmin = env%pcmin
  csthr = env%csthr
  pccap = env%pccap
  autolimit = .true.
  fraclimit = 0.25d0

  if (env%maxcluster == 0) then
    call clustleveval(env,nclustmax,csthr,SSRSSTthr,pcthr)
  end if

! ── topology for reference structure ─────────────────────────────────────────────
  if (env%wbotopo) then
    env%wbofile = 'wbo'
  else
    env%wbofile = 'none given'
  end if
  call simpletopo(nat,mols(1)%at,mols(1)%xyz,zmol,pr,.false.,env%wbofile)
  allocate (inc(zmol%nat),source=0)

!===========================================================!
  if (measuretype .ne. 'dihedral') then
!===========================================================!
    if (pr) then
      write (stdout,*)
      call smallhead('READING NUCLEAR EQUIVALENCIES')
    end if
    call readequals('anmr_nucinfo',zmol,groups)
    if (pr) then
      call groups%prsum(6)
      write (stdout,'(1x,a)') 'Unlisted nuclei (groups) are unique.'
    end if

    if (pr) then
      write (stdout,*)
      call smallhead('ANALYZING EQUIVALENCIES')
    end if
    call distsubgr(zmol,groups,subgroups,inc,pr)

    if (pr) then
      write (stdout,*)
      call smallhead('DETERMINE ATOMS TO INCLUDE IN PCA')
    end if
    call excludeFromRMSD(zmol,inc)
    if (sum(inc) == 0) then
      if (pr) then
        write (stdout,*) 'WARNING: No atoms included in PCA'
        write (stdout,*) 'Including more atoms ...'
      end if
      inc = 1
      do i = 1,groups%ng
        if (groups%grp(i)%nm > 1) then
          write (stdout,*) groups%grp(i)%mem
          do j = 1,groups%grp(i)%nm
            k = groups%grp(i)%mem(j)
            inc(k) = 0
          end do
        end if
      end do
    end if
    if (env%pcaexclude) then
      call excludeSelected(zmol,inc,env%atlist)
    end if
    if (heavyonly) then
      call excludeLight(zmol,inc)
    end if
    if (pr) then
      do i = 1,zmol%nat
        if (inc(i) == 1) then
          write (stdout,'(1x,a,a,i0,a,5x,a)') zmol%zat(i)%el,'(',i,')','taken'
        end if
      end do
    end if
    ntaken = sum(inc)
    ! ── fallback: include all heavy atoms if too few were selected ──────────────
    if (ntaken <= 3) then
      do i = 1,zmol%nat
        if (zmol%at(i) /= 1) then
          inc(i) = 1
          if (pr) write (stdout,'(1x,a,a,i0,a,5x,a)') zmol%zat(i)%el,'(',i,')','taken'
        end if
      end do
    end if
    ntaken = sum(inc)

    call zmol%deallocate

    ! ── for large ensembles, limit the number of structures considered ──────────
    if (autolimit) then
      if ((env%nclust /= 0).and.(env%nclust*100 < nall)) then
        dum = float(nall)*fraclimit
        dum2 = float(env%nclust)
        nallnew = nint(max(dum,dum2))
      else
        nallnew = nall
      end if
    else
      nallnew = nall
    end if

    allocate (xyznew(3,ntaken,nallnew))
    do i = 1,nallnew
      k = 0
      do j = 1,nat
        if (inc(j) == 1) then
          k = k+1
          xyznew(:,k,i) = mols(i)%xyz(:,j)
        end if
      end do
    end do

!===================================================!
  else !measuretype=='dihedral'
!===================================================!

    if (autolimit) then
      if ((env%nclust /= 0).and.(env%nclust*100 < nall)) then
        dum = float(nall)*fraclimit
        dum2 = float(env%nclust)
        nallnew = nint(max(dum,dum2))
      else
        nallnew = nall
      end if
    else
      nallnew = nall
    end if

    inc = 1
    ntaken = sum(inc)

    call zmol%countbonds()
    nbnd = zmol%nb
    allocate (diedat(4,zmol%nb),source=0)
    call getdiederatoms(zmol,zmol%nat,inc,nbnd,diedat,ndied)
    ntaken = ndied
!==================================================!
  end if
!==================================================!

! ── SVD: principal component analysis ────────────────────────────────────────────
  if (ntaken > 3) then  !> requires at least 4 descriptors
    call ctimer%start(1,'PCA')
    if (pr) then
      write (stdout,*)
      call smallhead('PRINCIPAL COMPONENT ANALYSIS')
    end if
    mm = nallnew
    select case (measuretype)
    case ('cma','CMA','cmadist')
      if (pr) then
        write (stdout,'(1x,a)') 'Using CMA DISTANCES as descriptors:'
      end if
      !>-- all structures should have been shifted to the CMA by CREGEN;
      !>   somewhat robust but provides less information than zmatrix
      mn = min(ntaken,mm)
      allocate (measure(mn,mm),pc(mn),pcvec(mm,mn))
      do i = 1,mm
        do j = 1,mn
          measure(j,i) = xyznew(1,j,i)**2+ &
       &                 xyznew(2,j,i)**2+ &
       &                 xyznew(3,j,i)**2
          measure(j,i) = sqrt(measure(j,i))
        end do
      end do
    case ('cartesian','coords')
      if (pr) then
        write (stdout,'(1x,a)') 'Using CARTESIAN COORDINATES as descriptors:'
      end if
      !>-- REQUIRES PERFECT ALIGNMENT; not robust for flexible molecules
      mn = min(ntaken*3,mm)
      allocate (measure(mn,mm),pc(mn),pcvec(mm,mn))
      do i = 1,mm
        l = 0
        do j = 1,ntaken
          do k = 1,3
            l = l+1
            if (l > mn) exit
            measure(l,i) = xyznew(k,j,i)
          end do
          if (l > mn) exit
        end do
      end do
    case default !case( 'zmat','zmatrix' )
      if (pr) then
        write (stdout,'(1x,a)') 'Using ZMATRIX as descriptors (sin/cos of dihedrals):'
      end if
      !>-- dihedral angles, sin/cos transformed for periodicity
      l = ntaken-3  !>-- first three dihedral angles are zero by convention
      mn = min(mm,2*l)  !>-- two descriptors per dihedral, capped at structure count
      if (mn < 2) then  !>-- need at least one dihedral angle (two descriptors)
        if (pr) then
          write (stdout,*) "Not enough descriptors for PCA!"
          return
        end if
      end if
      allocate (measure(mn,mm),pc(mn),pcvec(mm,mn))
      allocate (geo(3,ntaken),source=0.0d0)
      allocate (na(ntaken),nb(ntaken),nc(ntaken))
      do i = 1,mm
        na = 0; nb = 0; nc = 0
        geo = 0.0d0
        call xyzint(xyznew(1:3,1:ntaken,i),ntaken,na,nb,nc,radtodeg,geo)
        do j = 1,mn/2
          k = j+3
          dum = geo(3,k)*degtorad
          measure(2*j-1,i) = sin(dum)
          measure(2*j,i) = cos(dum)
        end do
      end do
      deallocate (nc,nb,na,geo)
    case ('dihedral')
      mn = min(mm,2*ntaken)
      allocate (measure(mn,mm),diedr(ndied))
      if (pr) then
        write (stdout,'(1x,a)') 'Using DIHEDRAL ANGLES as descriptors (sin/cos transformed):'
        do i = 1,ntaken
          write (stdout,'(1x,a,4i6)') 'Atoms: ',diedat(1:4,i)
        end do
        write (stdout,*)
      end if
      do i = 1,mm
        call calc_dieders(mols(i),ndied,diedat,diedr)
        do j = 1,min(ntaken,mn/2)
          measure(2*j-1,i) = sin(diedr(j))
          measure(2*j,i) = cos(diedr(j))
        end do
      end do
      if (allocated(diedat)) deallocate (diedat)
      if (allocated(diedr)) deallocate (diedr)
      allocate (pc(mn),pcvec(mm,mn))
    end select
    if (pr) then
      write (stdout,*)
      write (stdout,'(1x,a,i0,a,i0,a)') 'Performing SVD for ', &
&          mm,' structures and ',mn,' props'
    end if
    call SVD_to_PC(measure,mm,mn,pc,pcvec,.false.)  !> MM >= MN required
    call ctimer%stop(1)
  else
    write (stdout,*) 'There are not enough descriptors for a PCA!'
    write (stdout,*) 'Taking all structures as representative and writing ',clusterfile
    open (newunit=ich,file=clusterfile)
    do i = 1,nall
      call mols(i)%append(ich)
    end do
    close (ich)
    return
  end if

  if (allocated(measure)) deallocate (measure)
  if (allocated(xyznew)) deallocate (xyznew)
  if (allocated(inc)) deallocate (inc)

  ! ── normalize eigenvalues and select contributing PCs ──────────────────────────
  pcsum = sum(pc)
  pc = pc/pcsum
  pcsum = 0.0d0
  npc = 0
  do i = 1,mn
    if (pc(i) < pcmin) exit
    pcsum = pcsum+pc(i)
    npc = npc+1
    if (pcsum .ge. pcthr) exit
  end do
  ! ── fall back if the pcmin cutoff discarded every PC ────────────────────────────
  !> For floppy systems with many descriptors the variance can be spread so evenly
  !> that even the largest normalized eigenvalue lies below pcmin, leaving npc=0.
  !> Clustering on zero dimensions yields a NaN DBI/pSF table and merges everything
  !> into a single cluster. Redo the selection ignoring pcmin (as with -nopcmin) so
  !> that at least the pcthr-relevant PCs are kept, and warn the user.
  if (npc == 0) then
    write (stdout,'(1x,a,f5.3,a)') 'WARNING: no principal component eigenvalue exceeds pcmin (', &
    &  pcmin,'); the'
    write (stdout,'(1x,a)') '         variance is spread too evenly across descriptors. Ignoring the'
    write (stdout,'(1x,a)') '         pcmin cutoff for this ensemble (equivalent to -nopcmin).'
    pcsum = 0.0d0
    do i = 1,mn
      pcsum = pcsum+pc(i)
      npc = npc+1
      if (pcsum .ge. pcthr) exit
    end do
  end if
  npc = min(npc,pccap)
  pcsum = 0.0d0
  do i = 1,npc
    pcsum = pcsum+pc(i)
  end do

  if (pr) then
    i = min(100,MM)
    k = min(npc,6)
    write (stdout,*)
    call smallhead('EIGENVECTORS AND NORMALIZED EIGENVALUES OF SVD ANALYIS')
    call PRMAT(6,pcvec,i,k,'Eigenvectors of principal components')
    write (stdout,'(1x,a,i0,a)') 'NOTE: eigenvectors are only shown for the first ',i,' structures'
    write (stdout,'(1x,a,i0,a)') '      and the first ',k,' contributing principal components.'
    write (stdout,*)

    write (stdout,*) mn,'principal component eigenvalues (normalized)'
    write (stdout,*) pc
    write (stdout,*)
    write (stdout,'(1x,a,i0,a,f6.2,a)') 'The first ',npc,' components account for a total of ',100.d0*pcsum,'% of the'
    write (stdout,'(1x,a)') 'ensembles unique structural features and are used for the clustering'
  end if

  !>-- rearrange eigenvectors: drop unused PCs and swap to (npc,mm) layout
  allocate (pcdum(npc,mm))
  do i = 1,mm
    pcdum(1:npc,i) = pcvec(i,1:npc)
  end do
  call move_alloc(pcdum,pcvec)  !> pcvec shape changes from (mm,mn) to (npc,mm)

! ── k-means clustering ───────────────────────────────────────────────────────────
  if (pr) then
    write (stdout,*)
    call smallhead('CLUSTERING ANALYSIS OF PRINCIPAL COMPONENTS')
  end if

  allocate (member(mm),source=0)

  !>-- packed distance matrix; split to avoid integer overflow for large mm
  ndist = mm
  ndist = ndist*(mm+1)
  ndist = ndist/2
  allocate (dist(ndist),source=0.0_sp)
  allocate (p(npc),q(npc),source=0.0_ap)

  do i = 1,mm
    p(1:npc) = pcvec(1:npc,i)
!$OMP PARALLEL PRIVATE ( j, klong, q, dum ) &
!$OMP SHARED ( i, dist, npc, p, pcvec )
!$OMP DO
    do j = 1,i
      q(1:npc) = pcvec(1:npc,j)
      dum = eucdist(npc,p,q)
      klong = lina(i,j)
      dist(klong) = real(dum,sp)
    end do
!$OMP END DO
!$OMP END PARALLEL
  end do

  if (pr) then
    select case (clusteralgo)
    case ('means','kmeans')
      write (stdout,'(1x,a)') 'Using a MEANS cluster algorithm.'
    end select
    write (stdout,'(1x,a)') 'For a good review of cluster algorithms see'
    write (stdout,'(1x,a)') 'JCTC, 2007, 3, 2312 (doi.org/10.1021/ct700119m)'
    write (stdout,*)
    write (stdout,'(1x,a)') 'DBI = Davies-Bouldin index'
    write (stdout,'(1x,a)') 'pSF = pseudo F-statistic'
    write (stdout,'(1x,a)') 'SSR/SST = ratio of explained and unexplained variation'
    write (stdout,*)
    write (stdout,'(1x,a8,4x,a14,4x,a14,4x,a14)') 'Nclust','DBI','pSF','SSR/SST'
    write (stdout,'(1x,a8,4x,a14,4x,a14,4x,a14)') '------','-------------','-------------','-------------'
  end if

! ── cluster evaluation settings ──────────────────────────────────────────────────
  if (env%maxcluster == 0) then
    call clustleveval(env,nclustmax,csthr,SSRSSTthr,pcthr)
    nclustmax = min(mm,nclustmax)
  else
    nclustmax = max(2,env%maxcluster)
    nclustmax = min(mm,env%maxcluster)
  end if
  if (env%nclust == 0) then
    nclustmin = 1
  else
    nclust = min(mm,env%nclust)
    nclustmin = nclust
    nclustmax = nclust
  end if

  allocate (statistics(3,nclustmax),source=0.0d0)
  allocate (clust_sizes(nclustmax),source=0)
  CLUSTERSIZES: do nclustiter = nclustmin,nclustmax

    nclust = nclustiter
    if (env%clustlev >= 10) then
      dum = float(mm)/float(nclustmax)
      dum2 = dum*float(nclustiter)
      nclust = nint(dum2)
    end if
    clust_sizes(nclustiter) = nclust

    allocate (centroid(npc,nclust),source=0.0_ap)
    centroid = 0.0_ap

    select case (clusteralgo)
    case ('means','kmeans')
      call ctimer%start(2,'k-Means clustering')
      call kmeans(nclust,npc,mm,centroid,pcvec,ndist,dist,member)
      call ctimer%stop(2)
    end select

    call ctimer%start(3,'statistics')
    call cluststat(nclust,npc,mm,centroid,pcvec,member,DBI,pSF,SSRSST)
    if (pr) then
      write (stdout,'(1x,i8,4x,f14.6,4x,f14.6,4x,f14.6)') nclust,DBI,pSF,SSRSST
    end if
    call ctimer%stop(3)
    deallocate (centroid)

    statistics(1,nclustiter) = DBI
    statistics(2,nclustiter) = pSF
    statistics(3,nclustiter) = SSRSST

    if (nclust == env%nclust) exit
    if (SSRSST > SSRSSTthr) exit
  end do CLUSTERSIZES
  if (allocated(centroid)) deallocate (centroid)

  write (stdout,*)
  if (env%nclust == 0) then
    if (pr) then
      write (stdout,'(1x,a,i0,a)') 'Ensemble checked up to a partitioning into ',nclust,' clusters.'
      write (stdout,'(1x,a)') 'Local MINIMA of the DBI indicate adequate cluster counts.'
      write (stdout,'(1x,a)') 'Local MAXIMA of the pSF indicate adequate cluster counts.'
      write (stdout,'(1x,a)') 'Higher SSR/SST vaules indicate more distinct clusters.'
      write (stdout,'(1x,a)') 'Analyzing statistical values ...'
    end if
    k = min(nclustiter,nclustmax)  !> last completed iteration index
    allocate (extrema(2,k))
    call ctimer%start(3,'statistics')
    call statanal(k,nclustmax,statistics,extrema,pr,clust_sizes)
    if (pr) call statwarning(fname)
    ! ── pick smallest cluster count with adequate SSR/SST ──────────────────────
    do i = 2,k
      if ((extrema(1,i).or.extrema(2,i)).and.(statistics(3,i) > csthr)) then
        nclust = clust_sizes(i)
        exit
      end if
    end do
    call ctimer%stop(3)
    deallocate (extrema)
    if (pr) then
      write (stdout,*)
      write (stdout,'(1x,a,f4.2,a,i0)') 'Suggested (SSR/SST >',csthr,') cluster count: ',nclust
    end if
    allocate (centroid(npc,nclust),source=0.0_ap)
    select case (clusteralgo)
    case ('means','kmeans')
      call ctimer%start(2,'k-Means clustering')
      call kmeans(nclust,npc,mm,centroid,pcvec,ndist,dist,member)
      call ctimer%stop(2)
    end select
    deallocate (centroid)
  else
    if (pr) then
      write (stdout,'(1x,a,i0,a)') 'Ensemble partitioning into ',nclust,' clsuters.'
    end if
  end if
  deallocate (statistics,clust_sizes)

  deallocate (q,p,dist)
  call PCA_grpwrite(nclust,npc,mm,pcvec,member)

  ncb = nclust
  ancb = ncb

  if (ancb .le. 1) return

  if (pr) then
    write (stdout,*)
    write (stdout,'(1x,a)') 'Representative structures'
    write (stdout,'(1x,a6,1x,a6,3x,a6,1x,a16,1x,a16)') 'Nr.','conf.','clust.','Etot/Eh','Erel/ kcal/mol'
    write (stdout,'(1x,a6,1x,a6,3x,a6,1x,a16,1x,a16)') '---','-----','------','-------','--------------'
  end if
  allocate (eclust(ncb),source=0.0d0)
  allocate (clustbest(ncb),ind(ncb),source=0)
  iiincb: do i = 1,ncb
    do j = 1,mm
      if (member(j) == i) then
        eclust(i) = mols(j)%energy
        clustbest(i) = j
        cycle iiincb
      end if
    end do
  end do iiincb
  do i = 1,ncb
    ind(i) = i
    c = 0
    do j = 1,mm
      if (member(j) == i) then
        c = c+1
        if (mols(j)%energy < eclust(i)) then
          eclust(i) = mols(j)%energy
          clustbest(i) = j
        end if
      end if
    end do
    if (c == 0) then
      clustbest(i) = -1  !> empty cluster, excluded from output
      ancb = ancb-1
    end if
  end do
  call qsort(eclust,1,ncb,ind)
  emin = minval(eclust,1)
  open (newunit=ich,file=clusterfile)
  do i = 1,ncb
    k = clustbest(ind(i))
    if (k > 0) then
      dum = mols(k)%energy
      call mols(k)%append(ich)
      if (pr) then
        erel = (dum-emin)*autokcal
        write (stdout,'(1x,i6,1x,i6,3x,i6,1x,f16.8,1x,f16.4)') i,k,member(k),dum,erel
      end if
    end if
  end do
  close (ich)
  if (pr) then
    write (stdout,'(/,1x,a)') '(The "clust." column refers to the cluster "ID")'
    write (stdout,*)
    write (stdout,'(1x,a,a,a,i0,a)') 'File ',clusterfile,' written with ',ancb,' representative structures.'
    if (ancb < ncb) then
      write (stdout,'(1x,a,i0,a)') '(',ncb-ancb,' clusters discarded due to cluster merge)'
    end if
  end if
  if (allocated(mols)) deallocate (mols)

  if (pr) then
    write (stdout,*)
    call ctimer%write(stdout,'PCA/k-Means clustering')
  end if
  call ctimer%clear()
  return
end subroutine CCEGEN
