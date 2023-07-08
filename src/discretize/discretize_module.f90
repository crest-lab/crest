!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Philipp Pracht
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

module discretize_module
  use iso_fortran_env,only:wp => real64
  use probabilities_module
  implicit none

  private
  real(wp),parameter :: pi = acos(0.0_wp)*2.0_wp

!******************************************************************************!
  public :: struc_info
  type :: struc_info
    !> data from ensemble/simulation
    integer :: ndata
    real(wp),allocatable :: data(:)

    !> discretized data
    integer :: nmax
    real(wp),allocatable :: maxpos(:)
    real(wp),allocatable :: maxintv(:,:)
    real(wp),allocatable :: maxprob(:)

  contains
    procedure :: deallocate => struc_info_deallocate
    procedure :: setup => struc_info_setup

  end type struc_info
!******************************************************************************!

  public :: probability_count_minima

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine probability_count_minima(R1,R2,kappa,mu,N,ngrid)
!****************************************************************
!* The subroutine identifies the number of (potential) minima N
!* for a probability density within the interval [R1,R2] via
!* a grid-based approach
!****************************************************************
    implicit none
    !> INPUT
    real(wp),intent(in) :: R1,R2
    real(wp),intent(in) :: kappa
    real(wp),intent(in) :: mu(:)
    integer,intent(in),optional :: ngrid
    !> OUTPUT
    integer,intent(out) :: N
    !> LOCAL
    integer :: ng,i
    real(wp) :: diff,d,theta,dl,dl2
    real(wp) :: p,p_old,dpdt,dpdt_old,dpdt2,dpdt2_old
    N = 0
    if (present(ngrid)) then
      ng = ngrid
    else
      ng = 360
    end if
    diff = abs(R1-R2)
    d = diff/float(ng)
    theta = min(R1,R2)
    call vonMises(theta,kappa,mu,p_old,dpdt_old,dpdt2_old)
    do i = 2,ng
      theta = theta+d
      call vonMises(theta,kappa,mu,p,dpdt,dpdt2)
      dl = dpdt*dpdt_old
      dl2 = (dpdt2+dpdt2_old)/2.0_wp
      if (dl < 0.0_wp.and.dl2 >= 0.0_wp) then
        N = N+1
        write (*,*) theta*(180_wp/pi)
      end if
      p_old = p
      dpdt_old = dpdt
      dpdt2_old = dpdt2
    end do

  end subroutine probability_count_minima

!========================================================================================!

  subroutine probability_locate_minima(R1,R2,kappa,mu,N,minpos,ngrid)
!****************************************************************
!* The subroutine locates the N minima expected for a probability 
!* density within the interval [R1,R2] via a grid-based approach
!****************************************************************
    implicit none
    !> INPUT
    real(wp),intent(in) :: R1,R2
    real(wp),intent(in) :: kappa
    real(wp),intent(in) :: mu(:)
    integer,intent(in)  :: N
    integer,intent(in),optional :: ngrid
    !> OUTPUT
    real(wp),intent(out) :: minpos(N)
    !> LOCAL
    integer :: ng,i,k
    real(wp) :: diff,d,theta,dl,dl2,tmin
    real(wp) :: p,p_old,dpdt,dpdt_old,dpdt2,dpdt2_old
    if (present(ngrid)) then
      ng = ngrid
    else
      ng = 360
    end if
    minpos(:) = 0.0_wp
    k = 0
    diff = abs(R1-R2)
    d = diff/float(ng)
    theta = min(R1,R2)
    call vonMises(theta,kappa,mu,p_old,dpdt_old,dpdt2_old)
    do i = 2,ng
      theta = theta+d
      call vonMises(theta,kappa,mu,p,dpdt,dpdt2)
      dl = dpdt*dpdt_old
      dl2 = (dpdt2+dpdt2_old)/2.0_wp
      if (dl < 0.0_wp.and.dl2 >= 0.0_wp) then
        k = k+1
        !> TODO implement refinement function here
        tmin = theta - d*0.5_wp !> for now just place it between theta and theta+d
        minpos(k) = tmin     
      end if
      p_old = p
      dpdt_old = dpdt
      dpdt2_old = dpdt2
    end do

  end subroutine probability_locate_minima

!========================================================================================!
  subroutine probability_num_int(R1,R2,kappa,mu,norm,integral,ngrid)
!****************************************************************
!* The subroutine performs a numerical integration of the
!* probability density within the interval [R1,R2]
!****************************************************************
    implicit none
    !> INPUT
    real(wp),intent(in) :: R1,R2
    real(wp),intent(in) :: kappa
    real(wp),intent(in) :: mu(:)
    real(wp),intent(in) :: norm
    integer,intent(in),optional :: ngrid
    !> OUTPUT
    real(wp),intent(out) :: integral
    !> LOCAL
    integer :: ng,i
    real(wp) :: diff,d,theta,dl,dl2
    real(wp) :: p,p_old,dpdt,dpdt_old,dpdt2,dpdt2_old
    integral = 1.0_wp
    if (present(ngrid)) then
      ng = ngrid
    else
      ng = 36000
    end if
    diff = abs(R1-R2)
    d = diff/float(ng)
    theta = min(R1,R2)
    integral = 0.0_wp
    do i = 1,ng
      call vonMises(theta,kappa,mu,p)
      integral = integral+p*d
      theta = theta+d
    end do
    integral = norm*integral
  end subroutine probability_num_int

!========================================================================================!
  subroutine probability_num_norm(R1,R2,kappa,mu,norm,ngrid)
!****************************************************************
!* The subroutine performs a numerical integration of the
!* probability density within the interval [R1,R2] and
!* determines the normalization constant
!****************************************************************
    implicit none
    !> INPUT
    real(wp),intent(in) :: R1,R2
    real(wp),intent(in) :: kappa
    real(wp),intent(in) :: mu(:)
    integer,intent(in),optional :: ngrid
    !> OUTPUT
    real(wp),intent(out) :: norm
    !> LOCAL
    integer :: ng,i
    real(wp) :: diff,d,theta,dl,dl2
    real(wp) :: p,p_old,dpdt,dpdt_old,dpdt2,dpdt2_old
    norm = 1.0_wp
    if (present(ngrid)) then
      ng = ngrid
    else
      ng = 36000
    end if
    call probability_num_int(R1,R2,kappa,mu,1.0_wp,dl,ng)
    norm = 1.0_wp/dl
  end subroutine probability_num_norm

!========================================================================================!

subroutine struc_info_deallocate(self)
   implicit none
   class(struc_info) :: self

   self%ndata = 0
   if(allocated(self%data)) deallocate(self%data)

   self%nmax = 0
   if(allocated(self%maxintv)) deallocate(self%maxintv)
   if(allocated(self%maxprob)) deallocate(self%maxprob)
end subroutine struc_info_deallocate

subroutine struc_info_setup(self,kappa,mu)
   implicit none
   class(struc_info) :: self
   real(wp),intent(in) :: kappa
   real(wp),intent(in) :: mu(:)
   integer :: N,nmin
   real(wp) :: norm,R1,R2
   real(wp),allocatable :: minpos(:)
   integer :: i,j,k

   call self%deallocate()

   !> set data
   N = size(mu,1)
   self%ndata = N
   allocate(self%data(N))
   self%data(:) = mu(:)

   !> find minima
   R1 = 0.0_wp
   R2 = 2.0_wp*pi
   call probability_count_minima(R1,R2,kappa,mu,nmin)
   allocate(minpos(nmin))
   call probability_locate_minima(R1,R2,kappa,mu,nmin,minpos) 

   !> allocate intervals
   self%nmax = nmin !> nmin = nmax due to periodicity
   allocate(self%maxintv(2,nmin))
   allocate(self%maxprob(nmin))
   do i=1,nmin
 
   enddo  

   !> get normalization constant
   call probability_num_norm(R1,R2,kappa,mu,norm)


   !> integrate to get probabilities#
   do i = 1,nmin

   enddo
   

end subroutine struc_info_setup



!========================================================================================!
!========================================================================================!
end module discretize_module
