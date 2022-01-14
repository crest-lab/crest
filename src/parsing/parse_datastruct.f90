
module parse_datastruct
  use iso_fortran_env,only:wp => real64,sp => real32
  use parse_keyvalue
  use parse_block
  implicit none

  public :: root_object
  type :: root_object
    !> filename (if any)
    character(len=:),allocatable :: filename

    !> number & list of root key-value objects
    integer :: nkv = 0
    type(keyvalue),allocatable :: kv_list(:)

    !> number & list of blocks
    integer :: nblk
    type(datablock),allocatable :: blk_list(:)

  contains
    procedure :: new => root_deallocate
    procedure :: deallocate => root_deallocate
    procedure :: addkv => root_addkv
    procedure :: addblk => root_addblk
    procedure :: print => root_print
  end type root_object

contains

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
    j = i + 1
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
    j = i + 1
    allocate (newlist(j))
    newlist(1:i) = self%blk_list(1:i)
    newlist(j) = blk
    call move_alloc(newlist,self%blk_list)
    self%nblk = j
  end subroutine root_addblk

subroutine root_print(self)
   class(root_object) :: self
   integer :: i
   if(allocated(self%filename))then
   write(*,'(1x,a)')self%filename
   endif
   do i=1,self%nkv
     if(i < self%nkv)then
       write(*,'(1x,a,a)') '├──',trim(self%kv_list(i)%print() )
     else
       write(*,'(1x,a,a)') '└──',trim(self%kv_list(i)%print())
     endif
   enddo
   do i=1,self%nblk
     call self%blk_list(i)%print()
   enddo
end subroutine root_print


end module parse_datastruct
