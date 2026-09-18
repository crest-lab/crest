!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Philipp Pracht
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

!================================================================================!
! ENSEMBLE COMPARISON FUNCTION.
! To use execute:
!   crest --compare <ensemble1> <ensemble2>
!================================================================================!
subroutine compare_ensembles(env)
  !***********************************************************************
  !* Compare two molecular ensembles via iRMSD.                          *
  !* Reads both ensembles into coord arrays, sorts by energy, and        *
  !* computes iRMSD comparisons for the lowest structures using          *
  !* OMP-parallel permutation-invariant RMSD.                            *
  !* Per-structure tracking: group ID + best RMSD, O(ncomp1+ncomp2).     *
  !***********************************************************************
  use crest_parameters
  use crest_data
  use strucrd
  use quicksort_interface
  use axis_module
  use canonical_mod
  use irmsd_module
  use iomod
  use omp_lib
  use term_ui
  implicit none

  type(systemdata),intent(inout) :: env

  type(coord),allocatable :: strucs1(:),strucs2(:)
  type(coord),allocatable :: workmols(:)
  type(rmsd_cache),allocatable :: rcaches(:)
  type(canonical_sorter),allocatable :: sorters(:)
  real(wp),allocatable :: rmat(:,:),rmat_1d(:)

  !> per-structure best match and RMSD
  real(wp),allocatable :: best_rmsd_a(:),best_rmsd_b(:)
  integer,allocatable :: best_partner_a(:),best_partner_b(:)

  integer :: nall1,nall2,ncomp1,ncomp2,nat
  integer :: i,j,k,cc,T,Tn,npairs,ich
  integer :: nmatch_a,nmatch_b,pn1,pn2,pcount
  real(wp) :: RTHR,rmsdval,erel1,erel2,min_1,min_2,min_tot
  logical :: ex,stereocheck,first,store_rmat

  type(progress_state) :: pbar

  external :: PRMAT

! ── cleanup previous output files ───────────────────────────────────
  call compens_cleanup()

! ── validate input files ────────────────────────────────────────────
  inquire (file=env%ensemblename,exist=ex)
  if (.not.ex) then
    write (stdout,'(2x,a,a,a)') 'File <',trim(env%ensemblename),'> does not exist!'
    error stop
  end if
  inquire (file=env%ensemblename2,exist=ex)
  if (.not.ex) then
    write (stdout,'(2x,a,a,a)') 'File <',trim(env%ensemblename2),'> does not exist!'
    error stop
  end if

! ── read ensembles into coord arrays ────────────────────────────────
  write (stdout,'(1x,a,a,a)',advance='no') 'Reading ensemble <',trim(env%ensemblename),'> ...'
  flush (stdout)
  call rdensemble(env%ensemblename,nall1,strucs1)
  write (stdout,'(1x,i0,a)') nall1,' structures.'

  write (stdout,'(1x,a,a,a)',advance='no') 'Reading ensemble <',trim(env%ensemblename2),'> ...'
  flush (stdout)
  call rdensemble(env%ensemblename2,nall2,strucs2)
  write (stdout,'(1x,i0,a)') nall2,' structures.'
  write (stdout,*)

! ── validate compatibility ──────────────────────────────────────────
  nat = strucs1(1)%nat
  if (strucs1(1)%nat /= strucs2(1)%nat) then
    write (stdout,'(a)') "Nat1 /= Nat2 : Number of atoms of the two ensembles don't match!"
    write (stdout,'(a)') "You are trying to compare two different molecules!"
    error stop "exit."
  end if
  do i = 1,nat
    if (strucs1(1)%at(i) /= strucs2(1)%at(i)) then
      write (stdout,'(a)') "The ordering of atoms apparently is different between the two ensembles!"
      write (stdout,'(a)') "This way it is impossible to calculate RMSDs!"
      error stop "exit."
    end if
  end do

! ── sort each ensemble by energy ────────────────────────────────────
  call ensemble_qsort(nall1,strucs1,1,nall1)
  call ensemble_qsort(nall2,strucs2,1,nall2)

! ── select lowest structures ────────────────────────────────────────
  ncomp1 = min(nall1,env%maxcompare)
  ncomp2 = min(nall2,env%maxcompare)
  npairs = ncomp1*ncomp2

  call smallhead('Comparing the Ensembles')
  write (stdout,'(2x,a,a,a,i0,a,i0,a)') &
    'Ensemble A <',trim(env%ensemblename),'> : ', &
    nall1,' structures, using ',ncomp1,' lowest'
  write (stdout,'(2x,a,a,a,i0,a,i0,a)') &
    'Ensemble B <',trim(env%ensemblename2),'> : ', &
    nall2,' structures, using ',ncomp2,' lowest'

! ── set up OMP parallelization ──────────────────────────────────────
  call new_ompautoset(env,'max',0,T,Tn)
  write (stdout,'(2x,a,i0)') 'OpenMP threads: ',T

! ── set up RMSD threshold (in Bohr) ────────────────────────────────
  RTHR = env%rthr*aatoau
  write (stdout,'(2x,a,f8.4,a)') 'RMSD threshold: ',env%rthr,' Å'
  write (stdout,'(2x,a,i0)') 'Total comparisons: ',npairs
  write (stdout,*)

! ── set up canonical sorter and axis-align structures ───────────────
  write (stdout,'(1x,a)',advance='no') 'Setting up iRMSD infrastructure ...'
  flush (stdout)

  allocate (sorters(1))
  call axis(nat,strucs1(1)%at,strucs1(1)%xyz)
  call sorters(1)%init(strucs1(1),invtype='apsp+',heavy=.false.)
  stereocheck = .not.(sorters(1)%hasstereo(strucs1(1)))
  call sorters(1)%shrink()

  select case (env%iinversion)
  case (1)
    stereocheck = .true.
  case (2)
    stereocheck = .false.
  end select

  do i = 1,ncomp1
    call axis(nat,strucs1(i)%at,strucs1(i)%xyz)
  end do
  do i = 1,ncomp2
    call axis(nat,strucs2(i)%at,strucs2(i)%xyz)
  end do

! ── allocate per-thread work caches ─────────────────────────────────
  allocate (rcaches(T))
  allocate (workmols(T))
  do i = 1,T
    allocate (workmols(i)%at(nat))
    allocate (workmols(i)%xyz(3,nat))
    call rcaches(i)%allocate(nat)
    rcaches(i)%stereocheck = stereocheck
  end do
  write (stdout,'(1x,a)') 'done.'

! ── allocate per-structure tracking (O(ncomp1+ncomp2)) ──────────────
  store_rmat = (ncomp1 <= 20.and.ncomp2 <= 20)
  allocate (best_rmsd_a(ncomp1),source=huge(1.0_wp))
  allocate (best_rmsd_b(ncomp2),source=huge(1.0_wp))
  allocate (best_partner_a(ncomp1),source=0)
  allocate (best_partner_b(ncomp2),source=0)
  if (store_rmat) then
    allocate (rmat(ncomp1,ncomp2),source=0.0_wp)
  end if

! ── compute iRMSDs (OMP parallel) with progress bar ────────────────
  pcount = 0
  call progress_init(pbar,total=npairs,prefix=' iRMSD ')

  !$omp parallel &
  !$omp shared(strucs1,strucs2,rmat,store_rmat) &
  !$omp shared(best_rmsd_a,best_rmsd_b,best_partner_a,best_partner_b) &
  !$omp shared(sorters,rcaches,workmols,npairs,ncomp2,nat,RTHR,pbar) &
  !$omp private(k,i,j,cc,rmsdval)
  !$omp do schedule(dynamic)
  do k = 1,npairs
    cc = omp_get_thread_num()+1
    i = (k-1)/ncomp2+1
    j = mod(k-1,ncomp2)+1
    rcaches(cc)%rank(1:nat,1) = sorters(1)%rank(1:nat)
    rcaches(cc)%rank(1:nat,2) = sorters(1)%rank(1:nat)
    workmols(cc)%nat = nat
    workmols(cc)%at(:) = strucs2(j)%at(:)
    workmols(cc)%xyz(:,:) = strucs2(j)%xyz(:,:)
    call min_rmsd(strucs1(i),workmols(cc),rcache=rcaches(cc),rmsdout=rmsdval)
    if (store_rmat) rmat(i,j) = rmsdval
    if (rmsdval < RTHR) then
      !$omp critical
      if (rmsdval < best_rmsd_a(i)) then
        best_rmsd_a(i) = rmsdval
        best_partner_a(i) = j
      end if
      if (rmsdval < best_rmsd_b(j)) then
        best_rmsd_b(j) = rmsdval
        best_partner_b(j) = i
      end if
      !$omp end critical
    end if
    !$omp atomic
    pcount = pcount+1
    !$omp end atomic
    if (cc == 1) call progress_update(pbar,pcount,npairs)
  end do
  !$omp end do
  !$omp end parallel

  call progress_update(pbar,npairs,npairs,force=.true.)
  call progress_finish(pbar)
  write (stdout,*)

! ── print RMSD matrix (only for small matrices) ────────────────────
  if (store_rmat) then
    pn1 = ncomp1
    pn2 = ncomp2
    allocate (rmat_1d(pn1*pn2))
    do j = 1,pn2
      do i = 1,pn1
        rmat_1d((j-1)*pn1+i) = rmat(i,j)*autoaa
      end do
    end do
    call PRMAT(stdout,rmat_1d,pn1,pn2,'RMSD (Angstrom)')
    deallocate (rmat_1d,rmat)
  end if

! ── correlation printout ────────────────────────────────────────────
  min_1 = strucs1(1)%energy
  min_2 = strucs2(1)%energy
  min_tot = min(min_1,min_2)

  call smallhead('Correlation between Structures')
  write (stdout,'(2x,a,30x,a)') 'Ensemble A','Ensemble B'
  write (stdout,'(2x,a4,2x,a14,10x,a4,2x,a14,3x,a8)') &
    '#','Erel/kcal','#','Erel/kcal','RMSD/Å'
  write (stdout,'(2x,a)') repeat('-',62)

  do i = 1,ncomp1
    erel1 = (strucs1(i)%energy-min_1)*autokcal
    if (best_partner_a(i) > 0) then
      j = best_partner_a(i)
      erel2 = (strucs2(j)%energy-min_2)*autokcal
      write (stdout,'(2x,i4,2x,f14.5,2x,a5,2x,i4,2x,f14.5,3x,f8.4)') &
        i,erel1,'<--->',j,erel2,best_rmsd_a(i)*autoaa
    else
      write (stdout,'(2x,i4,2x,f14.5)') i,erel1
    end if
  end do

  write (stdout,'(2x,a)') repeat('-',62)

  nmatch_a = count(best_partner_a > 0)
  nmatch_b = count(best_partner_b > 0)
  write (stdout,'(2x,i0,a,i0,a)') nmatch_a,' of ',ncomp1, &
    ' structures in A have a match in B'
  write (stdout,'(2x,i0,a,i0,a)') nmatch_b,' of ',ncomp2, &
    ' structures in B have a match in A'

  first = .true.
  do i = 1,ncomp1
    if (best_partner_a(i) == 0) then
      if (first) then
        write (stdout,'(2x,a)',advance='no') 'Unmatched in A:'
        first = .false.
      end if
      write (stdout,'(1x,i0)',advance='no') i
    end if
  end do
  if (.not.first) write (stdout,*)

  first = .true.
  do j = 1,ncomp2
    if (best_partner_b(j) == 0) then
      if (first) then
        write (stdout,'(2x,a)',advance='no') 'Unmatched in B:'
        first = .false.
      end if
      write (stdout,'(1x,i0)',advance='no') j
    end if
  end do
  if (.not.first) write (stdout,*)

! ── write output files ──────────────────────────────────────────────
  open (newunit=ich,file='energy_1.dat')
  do i = 1,ncomp1
    write (ich,'(2x,f10.5,2x,f14.8)') &
      (strucs1(i)%energy-min_tot)*autokcal,strucs1(i)%energy
  end do
  close (ich)

  open (newunit=ich,file='energy_2.dat')
  do i = 1,ncomp2
    write (ich,'(2x,f10.5,2x,f14.8)') &
      (strucs2(i)%energy-min_tot)*autokcal,strucs2(i)%energy
  end do
  close (ich)

  open (newunit=ich,file='rmsdmatch.dat')
  do i = 1,ncomp1
    if (best_partner_a(i) > 0) then
      write (ich,'(2x,i6,i6,f10.4)') i,best_partner_a(i),best_rmsd_a(i)*autoaa
    end if
  end do
  close (ich)

! ── cleanup ─────────────────────────────────────────────────────────
  if (allocated(rmat)) deallocate (rmat)
  if (allocated(best_rmsd_a)) deallocate (best_rmsd_a)
  if (allocated(best_rmsd_b)) deallocate (best_rmsd_b)
  if (allocated(best_partner_a)) deallocate (best_partner_a)
  if (allocated(best_partner_b)) deallocate (best_partner_b)
  if (allocated(rcaches)) deallocate (rcaches)
  if (allocated(workmols)) deallocate (workmols)
  if (allocated(sorters)) deallocate (sorters)
  if (allocated(strucs1)) deallocate (strucs1)
  if (allocated(strucs2)) deallocate (strucs2)

end subroutine compare_ensembles

!---------------------------------------------------------------------------------------
subroutine compens_cleanup()
  !***********************************************
  !* Remove output files from a previous run.    *
  !***********************************************
  use iomod
  implicit none
  call remove('energy_1.dat')
  call remove('energy_2.dat')
  call remove('rmsdmatch.dat')
end subroutine compens_cleanup

