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

!> Central registry for error codes
module mctc_env_error
   implicit none
   private

   public :: mctc_stat, error_type
   public :: fatal_error


   !> Possible error codes
   type :: enum_stat

      !> Successful run
      integer :: success = 0

      !> Internal error:
      integer :: fatal = 1

   end type enum_stat

   !> Actual enumerator for return states
   type(enum_stat), parameter :: mctc_stat = enum_stat()


   !> Error message
   type :: error_type

      !> Error code
      integer :: stat

      !> Payload of the error
      character(len=:), allocatable :: message

   end type error_type


contains


!> A fatal error is encountered
subroutine fatal_error(error, message, stat)

   !> Instance of the error
   type(error_type), allocatable, intent(out) :: error

   !> A detailed message describing the error and (optionally) offering advice
   character(len=*), intent(in), optional :: message

   !> Overwrite of the error code
   integer, intent(in), optional :: stat

   allocate(error)

   if (present(stat)) then
      error%stat = stat
   else
      error%stat = mctc_stat%fatal
   end if

   if (present(message)) then
      error%message = message
   else
      error%message = "Fatal error"
   end if

end subroutine fatal_error


end module mctc_env_error
