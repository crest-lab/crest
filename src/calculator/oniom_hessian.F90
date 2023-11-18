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
  implicit none
  private

  !> Hessian storage object for ONIOM projections
  !> NOTE: Hessians are stored as packed matrices
  type :: hessian_storage
    integer :: dimM
    real(wp),allocatable :: hessian_calculated(:)
    integer :: dimN
    real(wp),allocatable :: hessian_projected(:)
    real(wp),allocatable :: JacobianT(:,:)
  end type hessian_storage

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine ONIOM_calc_hessians(mol,calc)
!*******************************************************
!* Calculate individual Hessians for each ONIOM level
!********************************************************
    implicit none
    type(coord),target :: mol
    type(calcdata) :: calc
    real(wp) :: energy,ep,em
    real(wp),allocatable :: gradient(:,:)
    real(wp),allocatable :: grad_plus(:,:),grad_minus(:,:)
    real(wp),allocatable :: hess(:,:)
    real(wp) :: hij
    integer :: iostatus
    integer :: i,j,k,l,ii,jj,i0,jid
    type(coord),pointer :: molptr
    integer :: pnat,p3
    type(hessian_storage),allocatable :: ONIOMhess(:)
    integer :: nhess
    real(wp),parameter :: step = 0.005_wp,step2 = 0.5_wp/step

    !> INPUT
#ifdef WITH_LWONIOM
!>--- Check # Hessians and allocate
    nhess = calc%ONIOM%ncalcs  
    allocate(ONIOMhess(nhess))

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

      write(stdout,'(a,i3,a)',advance='no') 'Calculating Hessian of ONIOM subsystem ',jid,' ...'
      flush(stdout)

      if (allocated(grad_plus)) deallocate (grad_plus)
      if (allocated(grad_minus)) deallocate (grad_minus)
      allocate (grad_plus(3,pnat),source=0.0_wp)
      allocate (grad_minus(3,pnat),source=0.0_wp)
      if (allocated(hess)) deallocate (hess)
      allocate (hess(p3,p3),source=0.0_wp)

      !!==========================================!
      !!> the actual potential call
      !!==========================================!
      !call potential_core(molptr,calc,i,iostatus)
      !!==========================================!
      !!==========================================!
      !energy = calc%etmp(i)
      !gradient = calc%grdtmp(:,1:pnat,i)

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
       
      write(stdout,'(a)') 'done.'


      if (iostatus /= 0) then
        return
      end if
    end do

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
end module oniom_hessian

