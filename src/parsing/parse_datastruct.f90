!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022-2023 Philipp Pracht
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

module parse_datastruct
  use crest_parameters
  use parse_keyvalue
  use parse_block
  use iomod,only:lowercase
  implicit none

  public :: root_object
  type :: root_object
    !> filename (if any)
    character(len=:),allocatable :: filename

    !> number & list of root key-value objects
    integer :: nkv = 0
    type(keyvalue),allocatable :: kv_list(:)

    !> number & list of blocks
    integer :: nblk = 0
    type(datablock),allocatable :: blk_list(:)

  contains
    procedure :: new => root_deallocate
    procedure :: deallocate => root_deallocate
    procedure :: addkv => root_addkv
    procedure :: addblk => root_addblk
    procedure :: lowercase_keys => root_lowercase_keys
    procedure :: print => root_print
  end type root_object

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!
  subroutine root_deallocate(self)
    implicit none
    class(root_object) :: self
    self%nkv = 0
    if (allocated(self%kv_list)) deallocate (self%kv_list)
    self%nblk = 0
    if (allocated(self%blk_list)) deallocate (self%blk_list)
    return
  end subroutine root_deallocate

  subroutine root_addkv(self,kv)
    implicit none
    class(root_object) :: self
    type(keyvalue) :: kv
    type(keyvalue),allocatable :: newlist(:)
    integer :: i,j
    i = self%nkv
    j = i+1
    allocate (newlist(j))
    newlist(1:i) = self%kv_list(1:i)
    newlist(j) = kv
    call move_alloc(newlist,self%kv_list)
    self%nkv = j
  end subroutine root_addkv

  subroutine root_addblk(self,blk)
    implicit none
    class(root_object) :: self
    type(datablock) :: blk
    type(datablock),allocatable :: newlist(:)
    integer :: i,j
    i = self%nblk
    j = i+1
    allocate (newlist(j))
    newlist(1:i) = self%blk_list(1:i)
    newlist(j) = blk
    call move_alloc(newlist,self%blk_list)
    self%nblk = j
  end subroutine root_addblk

  subroutine root_lowercase_keys(self)
    implicit none
    class(root_object) :: self
    integer :: i,j
    do i = 1,self%nkv
      self%kv_list(i)%key = lowercase(self%kv_list(i)%key)
    end do
    do i = 1,self%nblk
      self%blk_list(i)%header = lowercase(self%blk_list(i)%header)
      do j = 1,self%blk_list(i)%nkv
        self%blk_list(i)%kv_list(j)%key = lowercase(self%blk_list(i)%kv_list(j)%key)
      end do
    end do
  end subroutine root_lowercase_keys

  subroutine root_print(self)
    class(root_object) :: self
    integer :: i
    if (allocated(self%filename)) then
      write (stdout,'(1x,a)') self%filename
    end if
    do i = 1,self%nkv
      if (i < self%nkv) then
        write (stdout,'(1x,a,a)') '├──',trim(self%kv_list(i)%print())
      else
        write (stdout,'(1x,a,a)') '└──',trim(self%kv_list(i)%print())
      end if
    end do
    do i = 1,self%nblk
      call self%blk_list(i)%print()
    end do
  end subroutine root_print

!========================================================================================!
end module parse_datastruct
