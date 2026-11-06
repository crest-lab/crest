!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2025 Philipp Pracht
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

module unionize_module
  use crest_parameters
  use strucrd,only:coord
  use axis_module
  use irmsd_module
  use canonical_mod
  use quicksort_interface,only:ensemble_qsort
  implicit none
  private

!  logical,parameter :: debug = .true.
  logical,parameter :: debug = .false.

  public :: unionizeEnsembles

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine unionizeEnsembles(nin,inputs,nmerge,newmols,rthr,ethr)
!***************************************************************************
!* Merge two ensembles into a new one
!* We asume "inputs" is our reference into which "newmols" shall be merged.
!* "inputs" must be an allocatable list of structures and is OVERWRITTEN
!* by the new list.
!* Setting rthr and ethr to zero (or omitting the arguments)
!* will lead to every structure being identified as unique and
!* append it to the output ensemble.
!***************************************************************************
    implicit none
    integer,intent(inout) :: nin
    type(coord),allocatable,intent(inout) :: inputs(:)
    integer,intent(in) :: nmerge
    type(coord),intent(in),target :: newmols(nmerge)
    real(wp),intent(in),optional :: rthr,ethr
    !> LOCAL
    integer :: nout
    type(coord),allocatable :: structures(:)
    logical :: dupe,broken
    integer :: i,j,k,l,nat,ntaken,first,last
    type(canonical_sorter) :: newsort,refsort
    real(wp) :: rthr_ref,ethr_ref
    real(wp) :: rmsdval,deltaE
    type(coord),pointer :: mol
    logical :: topocheck
    type(rmsd_cache),allocatable :: rcache   !> similarity check cache (iRMSD)
    integer,allocatable :: similarto(:)

    nout = 0
    dupe = .false.
    broken = .false.
    topocheck = .true.
    nat = newmols(1)%nat
    if (present(ethr)) then
      ethr_ref = ethr
    else
      ethr_ref = 0.0_wp
    end if
    if (present(rthr)) then
      rthr_ref = rthr
    else
      rthr_ref = 0.0_wp
    end if
!>--- allocate mapping
    allocate (similarto(nmerge),source=0)

    !>--- we can skip the sorting if "inputs" is empty
    if (nin .ne. 0) then
!>--- Prepare comparison data storage
      if (debug) write (*,*)
      if (.not.allocated(rcache)) then
        if (debug) write (*,*) "allocating RCACHE"
        !$omp critical
        allocate (rcache)
        call rcache%allocate(nat)
        !$omp end critical
      end if
      !$omp critical
      call refsort%init(inputs(1),invtype='apsp+',heavy=.false.)
      call newsort%init(newmols(1),invtype='apsp+',heavy=.false.)
      !$omp end critical

!>--- double loop to count duplicates
      COMPAREOUTER: do i = 1,nin
        COMPAREINNER: do j = 1,nmerge
          if (similarto(j) .ne. 0) cycle COMPAREINNER
          mol => newmols(j)
          !> Energy difference
          deltaE = (mol%energy-inputs(i)%energy)*autokcal
          !> we can skip some comparisons if the energy difference is too large
          if (abs(deltaE) .gt. ethr) cycle COMPAREINNER

          !> Geometry difference (permutation-invariant RMSD)
          if (topocheck) then
            rcache%rank(1:nat,1) = newsort%rank(1:nat)
            rcache%rank(1:nat,2) = refsort%rank(1:nat)
          end if
          call min_rmsd(mol,inputs(i), &
          &        rcache=rcache,rmsdout=rmsdval)

          if (debug) write (*,'(a,es15.4,a,es15.4,a)') 'RMSD=',rmsdval*autoaa, &
          &          ' Å, delta E=',deltaE,' kcal/mol'

          !> Check
          if (abs(deltaE) .lt. ethr_ref.and.rmsdval*autoaa .lt. rthr_ref) then
            dupe = .true.
            similarto(j) = i
            if (deltaE < 0.0_wp) then
              !> if the energy is lower, we replace the molecule (better conformation)
              inputs(i) = mol
            end if
            exit COMPAREINNER
          end if
        end do COMPAREINNER
      end do COMPAREOUTER
      nullify (mol)
      !$omp critical
      call newsort%deallocate()
      call refsort%deallocate()
      !$omp end critical
    end if

    ntaken = count(similarto(:) .eq. 0)
    nout = nin+ntaken

!>--- after having checked the molecules, allocate new (output) space
    allocate (structures(nout))
    k = 0
    if (nin .ne. 0) then
      do i = 1,nin
        k = k+1
        structures(k) = inputs(i)
      end do
    end if
    do i = 1,nmerge
      if (similarto(i) .eq. 0) then
        k = k+1
        structures(k) = newmols(i)
      end if
    end do

!>--- for good measure, sort by energy again
    first = 1
    last = nout
    call ensemble_qsort(nout,structures,first,last)

!>-- overwrite "inputs"
    nin = nout
    if(allocated(inputs)) deallocate(inputs)
    call move_alloc(structures,inputs)

  end subroutine unionizeEnsembles

!=========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!=========================================================================================!
end module unionize_module
