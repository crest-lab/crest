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

module mutual_information_module
  use crest_parameters
  use discretize_module
  implicit none
  private

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

subroutine mie_entropy( dataset, N, entropy)
!**********************************************************************
!* This subroutine calculates the Mutual Information Expansion (MIE)
!* entropy up to a given order N for a data set A={A_1, A_2, ..., A_M}
!* with M discrete variables.
!*
!**********************************************************************
   implicit none
   !> INPUT
   integer,intent(in) :: dataset(:,:)  !> discrete data
   integer,intent(in) :: N             !> max order of MIE
   !> OUTPUT
   real(wp),intent(out) :: entropy     !> MIE entropy
   !> LOCAL
   integer :: ndata,M
   integer :: maxorder  !> "true" maxorder because N <= M
   integer :: k
   real(wp) :: marginal_entropy, factor

   !> initialize
   entropy = 0.0_wp
   M = size(dataset,1)
   ndata = size(dataset,2)
   maxorder = min(N,M)

   do k=1,maxorder
     entropy_k = 0.0_wp
     factor = 1.0_wp

     ! call mie_factor(
     ! call n_order_entropy(

     entropy = entropy + factor*marginal_entropy
   enddo

end subroutine mie_entropy

!========================================================================================!
subroutine n_order_entropy


   !>-- loops over all n_order subsets (probably needs to be done recursively)
   
   !...
   ! call n_order_subset_entropy


end subroutine  n_order_entropy

!========================================================================================!
subroutine n_order_subset_entropy
  
end subroutine n_order_subset_entropy

!========================================================================================!
!========================================================================================!
end module mutual_information_module
