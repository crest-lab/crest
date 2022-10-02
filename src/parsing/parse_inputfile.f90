
module parse_inputfile
  use iso_fortran_env,only:wp => real64,sp => real32
  use parse_keyvalue
  use parse_block
  use parse_datastruct
  use filemod
  implicit none
  private

  public :: parse_test
  public :: parse_input

contains
!========================================================================================!
  subroutine parse_test(fname)
    implicit none
    character(len=*) :: fname !> name of the input file
    type(root_object) :: dict
    call parse_input(fname,dict)
    call dict%print()
    return
  end subroutine parse_test
!========================================================================================!
  subroutine parse_input(fname,dict)
    implicit none

    character(len=*) :: fname !> name of the input file
    type(root_object),intent(out) :: dict
    type(filetype) :: file
    integer :: i,j,k,io
    logical :: get_root_kv
    type(keyvalue) :: kvdum
    !type(parseblock) :: rawblk
    type(datablock) :: blkdum 

    call dict%new()
!>--- open file to read and remove comments
    call file%open(trim(fname))
    dict%filename = trim(file%filename)
    call remove_comments(file)
    
    get_root_kv = .true.
!>--- the loop where the input file is read
    do i = 1,file%nlines
      if (file%current_line > i) cycle

      !> first, get key-value pairs of the root dict
      if(get_root_kv)then
        call get_keyvalue(kvdum,file%line(i),io)
        if(io == 0)then
          call dict%addkv(kvdum) !> add to dict
        endif
      endif

      if(isheader(file%line(i)))then
        get_root_kv = .false.
        call read_datablock(file,i,blkdum)
        !call blkdum%print()
        call dict%addblk(blkdum) !> add to dict
      endif
    end do

    call file%close()
    !file%filename = trim(fname)//'.tmp'
    !call file%flushclose()

    return
  end subroutine parse_input
!========================================================================================!
  subroutine remove_comments(file)
    implicit none
    type(filetype) :: file
    character(len=1),parameter :: com = '#'
    integer :: i
    do i = 1,file%nlines
      call clearcomment(file%f(i),com)
    end do
  end subroutine remove_comments
!========================================================================================!
end module parse_inputfile
