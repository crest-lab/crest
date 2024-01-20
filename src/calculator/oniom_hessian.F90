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
!********************************************************
!* Calculate individual Hessians for each ONIOM level
!********************************************************
    implicit none
    !> INPUT
    type(coord),target :: mol
    type(calcdata) :: calc
    !> OUTPUT
    real(wp),intent(out) :: hessian(mol%nat*3,mol%nat*3)
    !> LOCAL
    real(wp) :: energy,ep,em
    real(wp),allocatable :: gradient(:,:)
    real(wp),allocatable :: grad_plus(:,:),grad_minus(:,:)
    real(wp),allocatable :: hess(:,:)
    real(wp) :: hij
    integer :: iostatus
    integer :: i,j,k,l,ii,jj,i0,jid,n,ich
    type(coord),pointer :: molptr
    integer :: pnat,p3,highlowroot,Oid
    integer :: nhess
    character(len=100) :: atmp

    type(timer) :: profiler
    real(wp),parameter :: step = 0.005_wp,step2 = 0.5_wp/step

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

    !> allocate timer
    call profiler%init(n)

#ifdef WITH_LWONIOM

!>--- update ONIOMmols
    if (.not.allocated(calc%ONIOMmols)) then
      allocate (calc%ONIOMmols(calc%ONIOM%ncalcs))
    end if
    call ONIOM_update_geo(calc%ONIOM,mol,calc%ONIOMmols,calc%ONIOMmap)

!>--- Check # Hessians and allocate
    nhess = calc%ONIOM%ncalcs

!>--- loop over all calculations to be done, skip the ones
    do i = 1,calc%ncalculations
      call initsignal()
      !> Assign the molecule or skip
      if (calc%calcs(i)%ONIOM_id /= 0) then
        jid = calc%ONIOMrevmap(i)
        call ONIOM_associate_mol(calc%ONIOMmols(jid),molptr)
      else
        cycle
      end if
      pnat = molptr%nat
      p3 = 3*pnat
      highlowroot = calc%calcs(i)%ONIOM_highlowroot
      Oid = calc%calcs(i)%ONIOM_id

      !> some printpout
      write (atmp,'(1x,a,i3,a)') 'Calculating Hessian of ONIOM subsystem ',jid,' ...'
      write (stdout,'(/,a)') trim(atmp)
      write (stdout,'(1x,a)') repeat('-',len_trim(atmp))
      call calc%calcs(i)%info(stdout)
      write (stdout,'(" : ",a, i0,a,i0)') 'Hessian dimension ',3*pnat,' x ',3*pnat
      write (stdout,'(" : ",a,i0)') 'Number of e+grd calls (6N) ',6*pnat

      !> start timer
      call profiler%start(i)

      !> allocate gradients for plus and minus side
      if (allocated(grad_plus)) deallocate (grad_plus)
      if (allocated(grad_minus)) deallocate (grad_minus)
      allocate (grad_plus(3,pnat),source=0.0_wp)
      allocate (grad_minus(3,pnat),source=0.0_wp)
      if (allocated(hess)) deallocate (hess)
      allocate (hess(p3,p3),source=0.0_wp)

      !> calculate hessian for the fragment
      do i0 = 1,pnat
        do j = 1,3
          ii = (i0-1)*3+j

          !> positive direction
          grad_plus = 0.0_wp
          molptr%xyz(j,i0) = molptr%xyz(j,i0)+step
          !==========================================!
          call potential_core(molptr,calc,i,iostatus)
          !==========================================!
          grad_plus = calc%grdtmp(:,1:pnat,i)

          !> negative direction
          grad_minus = 0.0_wp
          molptr%xyz(j,i0) = molptr%xyz(j,i0)-2.0_wp*step
          !==========================================!
          call potential_core(molptr,calc,i,iostatus)
          !==========================================!
          grad_minus = calc%grdtmp(:,1:pnat,i)

          !> reset
          molptr%xyz(j,i0) = molptr%xyz(j,i0)+step

          !> count both calls in global variable
          engrad_total = engrad_total+2 !> global

          do k = 1,pnat
            do l = 1,3
              jj = (k-1)*3+l
              hess(jj,ii) = (grad_plus(l,k)-grad_minus(l,k))*step2
            end do
          end do

        end do
      end do

      !> Symmetrize Hessian
      do ii = 1,p3
        do j = ii,p3
          hij = (hess(ii,j)+hess(j,ii))*0.5_wp
          hess(ii,j) = hij
          hess(j,ii) = hij
        end do
      end do

      !write (atmp,'(a,i0,a,i0)') 'numhess.',Oid,'.',highlowroot
      !open (newunit=ich,file=trim(atmp))
      !write (ich,'(1x,a)') '$hessian'
      !do ii = 1,p3
      !  k = 0
      !  do j = 1,p3
      !    k = k+1
      !    if (k .le. 4) then
      !      write (ich,'(f16.8)',advance='no') hess(ii,j)
      !    else
      !      write (ich,'(f16.8)') hess(ii,j)
      !      k = 0
      !    end if
      !  end do
      !  if (k .ne. 0) then
      !    write (ich,*)
      !  end if
      !end do
      !write (ich,'(1x,a)') '$end'
      !close (ich)

      !> Pack the Hessian and put it into storage
      !write(*,*) Oid,highlowroot
      call lwonion_placehess(calc%ONIOM,mol%nat,pnat,hess, &
      & Oid,highlowroot)

      write (atmp,'(" : ",a,i3,a)') 'Hessian for subsytem ',jid,' done'
      call profiler%stop(i)
      call profiler%write_timing(stdout,i,trim(atmp))

      if (iostatus /= 0) then
        return
      end if
    end do

    call profiler%deallocate()

    write (stdout,'(/,a)',advance='no') ' Reconstructing ONIOM Hessian ...'
    flush (stdout)
    call lwoniom_gethess(mol%nat,calc%ONIOM,Hessian,.true.,iostatus)
    write (stdout,'(a)') ' done.'

#else
    call ONIOM_compile_error()
#endif
  end subroutine ONIOM_calc_hessians

!========================================================================================!

  subroutine ONIOM_compile_error()
    write (stdout,*) 'Error: Compiled without lwONIOM support!'
    write (stdout,*) 'Use -DWITH_LWONIOM=true in the setup to enable this function'
    error stop
  end subroutine ONIOM_compile_error

!========================================================================================!

!========================================================================================!
!========================================================================================!
end module oniom_hessian

