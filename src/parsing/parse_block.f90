
module parse_block
  use crest_parameters
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
    procedure :: print2 => blk_print2
    procedure :: fmt_header => blk_fmt_header
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

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

!> the following routines are type-bound procedures

!========================================================================================!
  subroutine blk_addkv(self,kv)
    implicit none
    class(datablock) :: self
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
  end subroutine blk_addkv

!========================================================================================!
  subroutine blk_deallocate(self)
    implicit none
    class(datablock) :: self
    self%header = ''
    self%nkv = 0
    if (allocated(self%kv_list)) deallocate (self%kv_list)
  end subroutine blk_deallocate

!========================================================================================!
  subroutine blk_fmt_header(self)
    implicit none
    class(datablock) :: self
    character(len=:),allocatable :: atmp
    integer :: l
    if (.not.allocated(self%header)) return

    self%header = adjustl(self%header)
    l = len_trim(self%header)
    if (self%header(1:1) .eq. '['.and. &
    &  self%header(l:l) .eq. ']') then
      atmp = self%header(2:l-1)
      call move_alloc(atmp,self%header)
    end if
  end subroutine blk_fmt_header

!========================================================================================!
  subroutine blk_print(self)
    class(datablock) :: self
    integer :: i
    write (stdout,'(1x,"object:",1x,a)') self%header
    do i = 1,self%nkv
      if (i < self%nkv) then
        write (stdout,'(1x,a,a)') '├──',trim(self%kv_list(i)%print())
      else
        write (stdout,'(1x,a,a)') '└──',trim(self%kv_list(i)%print())
      end if
    end do
  end subroutine blk_print

    subroutine blk_print2(self)
    class(datablock) :: self
    integer :: i
    if(index(self%header,'.').ne.0)then
      write (stdout,'("*",1x,a)') '[['//self%header//']]'
    else
      write (stdout,'("*",1x,a)') '['//self%header//']'
    endif
    do i = 1,self%nkv
      write (stdout,'("*",1x,a)') trim(self%kv_list(i)%print2())
    end do
  end subroutine blk_print2

!========================================================================================!

!> the following routines are only used in the fallback parsing routines
!> the standard implementation uses the toml-f library instead (see parse_toml.F90)

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
    blk%header = rawblk%header
    call clearheader(blk%header)

    do j = 1,rawblk%len
      call get_keyvalue(kvdum,rawblk%content(j),io)
      if (io == 0) then
        call blk%addkv(kvdum)
      end if
    end do

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
        l = l+1
      end if
      if (j == file%nlines) file%current_line = j
    end do jloop
    !if (l < 1) exit iloop
    if (l < 1) return
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
        l = l+1
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
!> print block data
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

  subroutine print_block2(self)
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
  end subroutine print_block2



!=======================================================================================!
!> check if string is a toml block header
  function isheader(str)
    implicit none
    logical :: isheader
    character(len=*) :: str
    character(len=:),allocatable :: atmp
    integer :: l
    isheader = .false.
    atmp = adjustl(trim(str))
    l = len_trim(atmp)
    if (l < 1) return
    if ((atmp(1:1) == '[').and.(atmp(l:l) == ']')) then
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
    k = len_trim(atmp)-1
    atmp = atmp(2:k)
    hdr = trim(atmp)
    return
  end subroutine clearheader

!========================================================================================!
end module parse_block
