! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
! Taken from github.com/grimme-lab/mctc-lib
!

!> Numerical storage size parameters for real and integer values
module mctc_env_accuracy
   implicit none
   public

   !> Single precision real numbers
   integer, parameter :: sp = selected_real_kind(6)

   !> Double precision real numbers
   integer, parameter :: dp = selected_real_kind(15)

   !> Wanted precision
   integer, parameter :: wp = dp

   !> Char length for integers
   integer, parameter :: i1 = selected_int_kind(2)

   !> Short length for integers
   integer, parameter :: i2 = selected_int_kind(4)

   !> Length of default integers
   integer, parameter :: i4 = selected_int_kind(9)

   !> Long length for integers
   integer, parameter :: i8 = selected_int_kind(18)


end module mctc_env_accuracy
