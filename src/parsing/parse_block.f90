
module parse_block
  use iso_fortran_env,only:wp => real64,sp => real32
  use filemod
  use parse_keyvalue
  implicit none

  public :: datablock
  type :: datablock
    !> the data object's name
    character(len=:),allocatable :: header
    !> number & list of key value pairs
    integer :: nkv = 0
    type(keyvalue),allocatable :: kv_list(:)
    contains
     procedure :: addkv => blk_addkv
     procedure :: print => blk_print
     procedure :: deallocate => blk_deallocate
  end type datablock

  public :: parseblock
  type :: parseblock
    integer :: len = 0
    character(len=:),allocatable :: header
    character(len=:),allocatable :: content(:)
  contains
    procedure :: deallocate => deallocate_block
    procedure :: print => print_block
  end type parseblock

  public :: parse_infile_block
  public :: isheader

contains
!========================================================================================!
  subroutine read_datablock(file,i,blk)
    implicit none
    type(filetype),intent(inout)  :: file
    type(datablock),intent(inout) :: blk
    integer,intent(in) :: i

    type(parseblock) :: rawblk
    type(keyvalue) :: kvdum
    integer :: j,k,io

    call blk%deallocate()

    call parse_infile_block(file,i,rawblk)
 !   call rawblk%print()

    blk%header = rawblk%header
    call clearheader(blk%header)

    do j=1,rawblk%len
       call get_keyvalue(kvdum,rawblk%content(j),io)
       if(io == 0)then
         call blk%addkv(kvdum)
       endif
    enddo

  end subroutine read_datablock
!========================================================================================!
!> for given input file parse the next block
  subroutine parse_infile_block(file,i,rawblk)
    implicit none
    type(filetype),intent(inout)      :: file
    type(parseblock),intent(inout) :: rawblk
    integer,intent(in) :: i
    logical :: saveblock
    integer :: j,k,l
    character(len=:),allocatable :: src

    call rawblk%deallocate()

    src = repeat(' ',file%lwidth)

    !saveblock = .false.
    !iloop: do i = 1,file%nlines
    !  if (i < file%current_line) cycle
    !  if (.not. saveblock) then
        if (isheader(file%line(i))) then
          saveblock = .true.
          rawblk%header = file%line(i)
    !      cycle
        end if
    !  end if
      !> get blocklength
      k = i+1
      l = 0
      jloop: do j = k,file%nlines
        if (isheader(file%line(j))) then
          file%current_line = j
          exit jloop
        end if
        if (len_trim(file%line(j)) > 0) then
          l = l + 1
        end if
        if(j==file%nlines) file%current_line = j
      end do jloop
      !if (l < 1) exit iloop
      if(l < 1) return 
      !> get block
      rawblk%len = l
      allocate (rawblk%content(l),source=src)
      l = 0
      jloop2: do j = k,file%nlines
        if (isheader(file%line(j))) then
          saveblock = .false.
          return
        end if
        if (len_trim(file%line(j)) > 0) then
          l = l + 1
          rawblk%content(l) = file%line(j)
        end if
      end do jloop2
    !end do iloop

    return
  end subroutine parse_infile_block

!========================================================================================!
!> deallocate block data
  subroutine deallocate_block(self)
    implicit none
    class(parseblock) :: self

    self%len = 0
    if (allocated(self%header)) deallocate (self%header)
    if (allocated(self%content)) deallocate (self%content)
    return
  end subroutine deallocate_block

!========================================================================================!
!> deallocate block data
  subroutine print_block(self)
    implicit none
    class(parseblock) :: self
    integer :: i

    write (*,*)
    if (allocated(self%header)) then
      write (*,*) trim(self%header)
      if (allocated(self%content)) then
        do i = 1,self%len
          write (*,*) self%content(i)
        end do
      end if
    end if
    return
  end subroutine print_block

!=======================================================================================!
!> check if string is a block header
  function isheader(str)
    implicit none
    logical :: isheader
    character(len=*) :: str
    character(len=:),allocatable :: atmp
    integer :: l
    isheader = .false.
    atmp = adjustl(trim(str))
    l = len_trim(atmp)
    if(l < 1) return
    if ((atmp(1:1) == '[') .and. (atmp(l:l) == ']')) then
      isheader = .true.
    end if
    return
  end function isheader

  subroutine clearheader(hdr)
    implicit none
    character(len=*) :: hdr
    integer :: i,k,l
    character(len=:),allocatable :: atmp,btmp
    character(len=1) :: s
    atmp = adjustl(hdr)
    atmp = trim(atmp)
    !>remove whitespaces
    l = len_trim(atmp)
    btmp = ''
    do i = 1,l
      s = atmp(i:i)
      if (s == ' ') cycle
      btmp = btmp//s
    end do
    atmp = btmp
    !> cut off "[" and "]"
    k = len_trim(atmp) - 1
    atmp = atmp(2:k)
     hdr = trim(atmp)
    return
  end subroutine clearheader

!========================================================================================!
  subroutine blk_addkv(self,kv)
    implicit none
    class(datablock) :: self
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
  end subroutine blk_addkv

!========================================================================================!
  subroutine blk_deallocate(self)
    implicit none
    class(datablock) :: self
    self%header = ''
    self%nkv = 0
    if(allocated(self%kv_list)) deallocate(self%kv_list)
  end subroutine blk_deallocate

!========================================================================================!
subroutine blk_print(self)
   class(datablock) :: self
   integer :: i
   write(*,'(1x,"object:",1x,a)')self%header
   do i=1,self%nkv
     if(i < self%nkv)then
       write(*,'(1x,a,a)') '├──',trim(self%kv_list(i)%print() )
     else
       write(*,'(1x,a,a)') '└──',trim(self%kv_list(i)%print())
     endif
   enddo
end subroutine blk_print
!========================================================================================!
end module parse_block
