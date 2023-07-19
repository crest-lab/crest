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

  integer,parameter :: hugeint = huge(hugeint)
  real(wp),parameter :: hugeint_as_real = float(hugeint)

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
     marginal_entropy = 0.0_wp
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
subroutine n_order_subset_entropy(M,N,dataset,weights,mask,work,p,entropy)
   implicit none
   !> INPUT
   integer,intent(in) :: M,N
   integer,intent(in) :: dataset(M,N)
   logical,intent(in) :: mask(M)
   integer,intent(inout) :: work(N,2) !> a work arry to avoid allocating one in the subroutine
   real(wp),intent(inout) :: p(N)  !> p is also a work array
   real(wp),intent(in) :: weights(N)
   !> OUTPUT
   real(wp),intent(out) :: entropy
   !> LOCAL
   integer :: i,j,k,l,ref
   logical :: check
   real(wp) :: sumweights
   
   !> reset
   entropy = 0.0_wp 
   work(:,:) = 0
   p(:) = 0.0_wp
   k = 0
   sumweights = sum(weights(:))
   !> the first is always the first wrt the subset
   work(1,2) = 1
   work(1,1) = 1
   p(1) = p(1) + weights(1)
   k = k + 1
   !> loop over the data 
   do i=2,N
   !TODO
     do j=1,i-1 !> better only over the k mappings?
       do l=1,M
          if(mask(l))then
          
          endif
       enddo 
     enddo
   enddo

   !> calculate probabilities
   p = p / sumweights 

   !> calculate Shannon entropy
   do i=1,k
     entropy = entropy - Rcal * p(i) * log( p(i) )
   enddo

end subroutine n_order_subset_entropy

!========================================================================================!
!========================================================================================!
end module mutual_information_module
