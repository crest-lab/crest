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
  use crest_parameters
  use probabilities_module
  implicit none

  private
!  real(wp),parameter :: pi = acos(0.0_wp)*2.0_wp

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

    !> marginal entropy
    real(wp) :: S

  contains
    procedure :: deallocate => struc_info_deallocate
    procedure :: setup => struc_info_setup
    procedure :: entropy => struc_marginal_entropy
    procedure :: info => struc_prinfo
    procedure :: discrete => struc_discrete
  end type struc_info
!******************************************************************************!

  public :: probability_count_minima
  public :: write_data_matrix_bin,read_data_matrix_bin

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
        !write (*,*) theta*(180_wp/pi)
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
        tmin = theta-d*0.5_wp !> for now just place it between theta and theta+d
        minpos(k) = tmin
      end if
      p_old = p
      dpdt_old = dpdt
      dpdt2_old = dpdt2
    end do

  end subroutine probability_locate_minima

!========================================================================================!
  subroutine probability_locate_maximum(R1,R2,kappa,mu,N,maxpos,ngrid)
!****************************************************************
!* The subroutine locates the maximum for a probability
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
    real(wp),intent(out) :: maxpos
    !> LOCAL
    integer :: ng,i,k
    real(wp) :: diff,d,theta,dl,dl2,tmin
    real(wp) :: p,p_old,dpdt,dpdt_old,dpdt2,dpdt2_old
    real(wp) :: maxref
    if (present(ngrid)) then
      ng = ngrid
    else
      ng = 1000
    end if
    maxpos = 0.0_wp
    k = 0
    diff = abs(R1-R2)
    d = diff/float(ng)
    theta = min(R1,R2)
    maxref = 0.0_wp
    do i = 1,ng
      call vonMises(theta,kappa,mu,p)
      if (p > maxref) then
        maxref = p
        maxpos = theta
      end if
      theta = theta+d
    end do

  end subroutine probability_locate_maximum

!========================================================================================!
  subroutine probability_num_int(R1,R2,kappa,mu,norm,integral,ngrid,maxpos)
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
    real(wp),intent(out),optional :: maxpos
    !> LOCAL
    integer :: ng,i
    real(wp) :: diff,d,theta,dl,dl2
    real(wp) :: p,p_old,dpdt,dpdt_old,dpdt2,dpdt2_old
    real(wp) :: maxref,maxp
    integral = 1.0_wp
    maxref = 0.0_wp
    maxp = R1
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
      if (p > maxref) then
        maxref = p
        maxp = theta
      end if
      integral = integral+p*d
      theta = theta+d
    end do
    integral = norm*integral
    if(present(maxpos)) maxpos = maxp
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
    if (allocated(self%data)) deallocate (self%data)

    self%nmax = 0
    if (allocated(self%maxintv)) deallocate (self%maxintv)
    if (allocated(self%maxprob)) deallocate (self%maxprob)
  end subroutine struc_info_deallocate

  subroutine struc_info_setup(self,kappa,mu)
!************************************************************************
!* This routine does the entire setup for a struc_info object
!* The input mu(:) is the data input (an array of angular values in rad)
!* and kappa is the von-Mises discretization parameter (empirical).
!* A probability density is set up from mu, minima and maxima are
!* identified and ranges are determined. Also calculates the probability
!* for each range and a mariginal (first order) entropy for the struc_info
!***********************************************************************
    implicit none
    class(struc_info) :: self
    real(wp),intent(in) :: kappa
    real(wp),intent(in) :: mu(:)
    integer :: N,nmin
    real(wp) :: norm,R1,R2
    real(wp),allocatable :: minpos(:)
    integer :: i,j,k,ngrid
    real(wp) :: integral,maxpos

    call self%deallocate()

    !> set data
    N = size(mu,1)
    self%ndata = N
    allocate (self%data(N))
    self%data(:) = mu(:)

    !> find minima
    R1 = 0.0_wp
    R2 = 2.0_wp*pi
    call probability_count_minima(R1,R2,kappa,mu,nmin)
    allocate (minpos(nmin))
    call probability_locate_minima(R1,R2,kappa,mu,nmin,minpos)

    !> allocate intervals
    self%nmax = nmin !> nmin = nmax due to periodicity
    allocate (self%maxintv(2,nmin))
    allocate (self%maxprob(nmin),self%maxpos(nmin))
    do i = 1,nmin
      self%maxintv(2,i) = minpos(i)
      if (i == 1) then
        self%maxintv(1,i) = minpos(nmin)-2.0_wp*pi
      else
        self%maxintv(1,i) = minpos(i-1)
      end if
    end do

    !> get normalization constant
    ngrid = 3600
    call probability_num_norm(R1,R2,kappa,mu,norm,ngrid)

    !> integrate to get probabilities and get maxima positions
    ngrid = nint(float(ngrid)/float(nmin))
    do i = 1,nmin
      R1 = self%maxintv(1,i)
      R2 = self%maxintv(2,i)
      call probability_num_int(R1,R2,kappa,mu,norm,integral,ngrid,maxpos)
      self%maxprob(i) = integral
!      call probability_locate_maximum(R1,R2,kappa,mu,N,maxpos)
      self%maxpos(i) = maxpos
    end do

    call self%entropy()

  end subroutine struc_info_setup

!========================================================================================!
  subroutine struc_marginal_entropy(self,entropy)
    implicit none
    class(struc_info) :: self
    real(wp),intent(out),optional :: entropy
    real(wp) :: S,p
    integer :: i
    real(wp),parameter :: R = Rcal

    S = 0.0_wp

    if (allocated(self%maxprob)) then
      do i = 1,self%nmax
        p = self%maxprob(i)
        if(p > 0.0_wp)then 
        S = S-R*p*log(p)
        endif
      end do
    end if

    self%S = S
    if (present(entropy)) entropy = S

  end subroutine struc_marginal_entropy


!========================================================================================!
  function struc_discrete(self,var) result(dvar)
     implicit none
     class(struc_info) :: self
     real(wp),intent(in) :: var
     integer :: dvar
     real(wp) :: vartmp
     integer :: i,j,m
     dvar = 0
     vartmp = var
     if(var < self%maxintv(1,1))then
        vartmp = var + 2.0_wp*pi
     endif
     m = self%nmax
     if(var > self%maxintv(2,m))then
       vartmp = var - 2.0_wp*pi
     endif
     do i=1,m
       if(vartmp >= self%maxintv(1,i) .and. &
       &  vartmp <= self%maxintv(2,i))then
         dvar = i
         exit
       endif
     enddo
  end function struc_discrete

!========================================================================================!
  subroutine struc_prinfo(self,k,conv)
    implicit none
    class(struc_info) :: self
    integer,intent(in),optional :: k
    real(wp),intent(in),optional :: conv
    real(wp) :: convert
    integer :: i

    write(stdout,*)
    if(present(k))then
    write(stdout,'(a,i5,2x,a,i5,2x,a,f12.6)') '> discrete variable',k,'  maxima',self%nmax, &
    & '  marginal entropy',self%S
    else
    write(stdout,'(a,7x,a,i5,2x,a,f12.6)') '> discrete variable','  maxima',self%nmax, &
    & '  marginal entropy',self%S

    endif

    convert = 1.0_wp
    if(present(conv))then
       convert = conv
    endif
    if(self%nmax > 0)then
    write(stdout,'(2x,3a18,a16)') 'max','lower','upper','p'
    do i=1,self%nmax
       write(stdout,'(2x,3f18.6,f16.8)') self%maxpos(i)*convert,self%maxintv(1:2,i)*convert, &
       &                                 self%maxprob(i)   
    enddo
    endif
   
  end subroutine struc_prinfo

!========================================================================================!


  subroutine write_data_matrix_bin(datmat)
     implicit none
     integer,intent(in) :: datmat(:,:)
     integer :: n,m
     integer :: ich,i
     n = size(datmat,1)
     m = size(datmat,2)
     open(newunit=ich, file='ddata.bin', action='readwrite', form='unformatted', &
     & access='stream', status='replace')
     write(ich) n,m
     write(ich) datmat
     close(ich)
     open(newunit=ich, file='ddata.txt')
     write(ich,*) n,m
       write(ich,*) datmat
     close(ich)
  end subroutine write_data_matrix_bin


  subroutine read_data_matrix_bin(n,m,datmat)
     implicit none
     integer,allocatable,intent(out) :: datmat(:,:)
     integer,intent(out) :: n,m
     integer :: ich,i
     open(newunit=ich, file='ddata.bin', action='readwrite', form='unformatted', &
     & access='stream', status='old')
     read(ich) n,m
     allocate(datmat(n,m), source = 0)
     read(ich) datmat
     close(ich)
  end subroutine read_data_matrix_bin


!========================================================================================!
!========================================================================================!
end module discretize_module
