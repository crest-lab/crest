!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Patryk Wesołowski, Philipp Pracht
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

module oniom_hessian
  use crest_parameters
  use strucrd
  use crest_calculator
#ifdef WITH_LWONIOM
  use lwoniom_interface
#endif
  use lwoniom_module
  use crest_type_timer
  implicit none
  private

  public :: ONIOM_calc_hessians

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine ONIOM_calc_hessians(mol,calc,hessian)
!**********************************************************************
!* Calculate the ONIOM Hessian: the finite-difference loops for the
!* individual model system Hessians, their projection into the full
!* system basis and the ONIOM reconstruction are all handled by
!* lwONIOM's numerical Hessian driver; crest supplies the potentials
!* through the crest_oniom_calc adapter (see calculator.F90).
!**********************************************************************
    implicit none
    !> INPUT
    type(coord),target :: mol
    type(calcdata),target :: calc
    !> OUTPUT
    real(wp),intent(out) :: hessian(mol%nat*3,mol%nat*3)
    !> LOCAL
    integer :: iostatus
    integer :: i,j,n
    type(timer) :: profiler
#ifdef WITH_LWONIOM
    type(crest_oniom_calc) :: OCLC
    type(lwoniom_job),allocatable :: jobs(:)
#endif

    !> Hessian
    hessian(:,:) = 0.0_wp

    !> allocate energy and gradient storage
    n = calc%ncalculations
    if (n > 0) then
      if (.not.allocated(calc%etmp)) allocate (calc%etmp(n),source=0.0_wp)
      if (.not.allocated(calc%grdtmp)) allocate (calc%grdtmp(3,mol%nat,n),source=0.0_wp)
      if (.not.allocated(calc%eweight)) then
        allocate (calc%eweight(n),source=0.0_wp)
        do i = 1,n
          calc%eweight(i) = calc%calcs(i)%weight
        end do
      end if
    else
      error stop '***ERROR*** no calculations allocated'
    end if

#ifdef WITH_LWONIOM

! ── persistent fragment buffers (refreshed by the adapter callback) ───────────
    if (.not.allocated(calc%ONIOMmols)) then
      allocate (calc%ONIOMmols(calc%ONIOM%ncalcs))
    end if

! ── some printout on the upcoming work ────────────────────────────────────────
    call lwoniom_get_jobs(calc%ONIOM,jobs)
    write (stdout,'(/,1x,a,i0,a)') 'Calculating ',size(jobs),' ONIOM subsystem Hessians ...'
    do j = 1,size(jobs)
      write (stdout,'(" : ",a,i3,a,i5,a,i0,a,i0)') 'subsystem ',j,' with ', &
      &  jobs(j)%nat,' atoms, Hessian dimension ',3*jobs(j)%nat,' x ',3*jobs(j)%nat
    end do

! ── FD Hessians, projection and ONIOM reconstruction, all via lwONIOM ─────────
    call profiler%init(1)
    call profiler%start(1)
    OCLC%calc => calc
    call lwoniom_numhess_driver(calc%ONIOM,OCLC,mol%nat,mol%xyz,hessian, &
    &                           iostat=iostatus,verbose=.true.)

    !> count the finite-difference engrad calls in the global counter
    do j = 1,size(jobs)
      engrad_total = engrad_total+real(6*jobs(j)%nat,wp)
    end do

    call profiler%stop(1)
    call profiler%write_timing(stdout,1,'ONIOM Hessian construction done,')
    call profiler%deallocate()

    if (iostatus /= 0) then
      write (stdout,'(a,i0)') '**ERROR** ONIOM Hessian construction failed with code ',iostatus
    end if

#else
    call ONIOM_compile_error()
#endif
  end subroutine ONIOM_calc_hessians

!========================================================================================!
!========================================================================================!
end module oniom_hessian
