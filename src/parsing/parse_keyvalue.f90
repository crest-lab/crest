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

module parse_keyvalue
  use crest_parameters
  implicit none
  private

!> the key-value data structure used by the crest input parser
  public :: keyvalue
  type :: keyvalue
    character(len=:),allocatable :: key      !> the key as string
    character(len=:),allocatable :: rawvalue !> the value as string
    integer :: id = 0     !> type of the value (see enumerator below)
    real(wp) :: value_f   !> id=1, a float
    integer  :: value_i   !> id=2, an integer
    logical  :: value_b   !> id=3, a boolean
    character(len=:),allocatable :: value_c !> id=4, a string
    integer ::  na !> number of array elements
    character(len=:),allocatable :: value_rawa(:)
    real(wp),allocatable :: value_fa(:) !> id=5, an array of floats
    integer,allocatable  :: value_ia(:) !> id=6, an array of integers
    logical,allocatable  :: value_ba(:) !> id=7, an array of booleans
    character(len=:),allocatable :: value_ca(:) !> id=8, an array of strings/multiline string
  contains
    procedure :: print => print_kv
    procedure :: print2 => print_kv2
    procedure :: deallocate => deallocate_kv
    procedure :: set_valuestring => kv_set_valuestring
    procedure :: add_raw_array_string => kv_add_raw_array_string
  end type keyvalue

  character(len=*),parameter :: kv_indicator = '='  !> used for fallback parsing

!> enumerator definition (and parameter) for value type
  type,private :: enum_valuetypes
    integer :: raw = 0
    integer :: unknown = 0
    integer :: float = 1
    integer :: int = 2
    integer :: bool = 3
    integer :: string = 4
    integer :: float_array = 5
    integer :: int_array = 6
    integer :: bool_array = 7
    integer :: string_array = 8
    integer :: multiline_string = 8
    integer :: raw_array = 9
  end type enum_valuetypes
  type(enum_valuetypes),public,parameter :: valuetypes = enum_valuetypes()

  public :: get_keyvalue

!========================================================================================!
!========================================================================================!
contains !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine deallocate_kv(self)
    implicit none
    class(keyvalue) :: self
    if (allocated(self%key)) deallocate (self%key)
    if (allocated(self%rawvalue)) deallocate (self%rawvalue)
    self%id = 0
    if (allocated(self%value_c)) deallocate (self%value_c)
    self%na = 0
    if (allocated(self%value_rawa)) deallocate (self%value_rawa)
    if (allocated(self%value_fa)) deallocate (self%value_fa)
    if (allocated(self%value_ia)) deallocate (self%value_ia)
    if (allocated(self%value_ba)) deallocate (self%value_ba)
    if (allocated(self%value_ca)) deallocate (self%value_ca)
    return
  end subroutine

!========================================================================================!

!> The following routines are used only in the fallback implementation!
!> The standard implementation works via the toml-f library

!========================================================================================!
  subroutine get_keyvalue(kv,str,io)
    implicit none
    class(keyvalue) :: kv
    character(len=*) :: str
    integer,intent(out) :: io
    integer :: i
    call kv%deallocate()
    io = 0
    !>--- if the sting is no key-value pair, return with error
    if (.not.iskeyvalue(str)) then
      io = -1
      return
    end if

    !>--- get the key and its corresponding (raw) value
    call get_key(kv,str)
    call get_rawvalue(kv,str)

    !>--- check types of raw value and read
    !> multiline comment (must be read outside this subroutine)
    if (value_ismultiline(kv%rawvalue)) then
      io = 1
      return

      !> inline array
    elseif (value_isarray(kv%rawvalue)) then
      call get_rawarray(kv)
      call get_typearray(kv)

      !> inline string
    else if (value_isstring(kv%rawvalue)) then
      call get_cvalue(kv)

      !> anything else (boolean, float, or integer)
    else
      call get_valueautotype(kv)
    end if

    return
  end subroutine get_keyvalue
!========================================================================================!
  function iskeyvalue(str)
    implicit none
    logical :: iskeyvalue
    character(len=*) :: str
    character(len=:),allocatable :: atmp
    integer :: k
    iskeyvalue = .false.
    atmp = adjustl(str)
    if (len_trim(atmp) < 1) return
    if (atmp(1:1) == '[') return
    k = index(atmp,kv_indicator)
    if (k .ne. 0) then
      iskeyvalue = .true.
    end if
    return
  end function iskeyvalue
!========================================================================================!
  subroutine get_key(kv,str)
    implicit none
    class(keyvalue) :: kv
    character(len=*) :: str
    integer :: k,l
    character(len=:),allocatable :: atmp
    atmp = adjustl(str)
    k = index(atmp,kv_indicator)
    l = k-1
    if (atmp(1:1) == '"'.or.atmp(1:1) == "'") then
      l = k-2
      kv%key = trim(atmp(2:l))
    else
      kv%key = trim(atmp(1:l))
    end if
    return
  end subroutine get_key
!========================================================================================!
  subroutine get_rawvalue(kv,str)
    implicit none
    class(keyvalue) :: kv
    character(len=*) :: str
    integer :: k,l
    character(len=:),allocatable :: atmp
    atmp = adjustl(str)
    k = index(atmp,kv_indicator)
    k = k+len(kv_indicator)
    l = len_trim(str)
    kv%rawvalue = trim(adjustl(atmp(k:l)))
    return
  end subroutine get_rawvalue
!========================================================================================!
  function value_ismultiline(str)
    implicit none
    logical :: value_ismultiline
    character(len=*) :: str
    character(len=:),allocatable :: atmp
    value_ismultiline = .false.
    atmp = adjustl(str)
    if (len_trim(atmp) < 3) return
    if (atmp(1:3) == "'''".or.atmp(1:3) == '"""') then
      value_ismultiline = .true.
    end if
  end function value_ismultiline
!========================================================================================!
  function value_isstring(str)
    implicit none
    logical :: value_isstring
    character(len=*) :: str
    character(len=:),allocatable :: atmp
    integer :: l
    value_isstring = .false.
    atmp = adjustl(trim(str))
    l = len_trim(atmp)
    if (.not.value_ismultiline(str)) then
      if ((atmp(1:1) == "'".and.atmp(l:l) == "'").or. &
      &  (atmp(1:1) == '"'.and.atmp(l:l) == '"')) then
        value_isstring = .true.
      end if
    end if
  end function value_isstring
!========================================================================================!
  subroutine get_cvalue(kv)
    implicit none
    class(keyvalue) :: kv
    integer :: k,l
    character(len=:),allocatable :: atmp
    atmp = kv%rawvalue
    l = len_trim(atmp)-1
    kv%value_c = atmp(2:l)
    kv%id = valuetypes%string
    return
  end subroutine get_cvalue
!========================================================================================!
  function value_isarray(str)
    implicit none
    logical :: value_isarray
    character(len=*) :: str
    character(len=:),allocatable :: atmp
    integer :: l
    value_isarray = .false.
    atmp = adjustl(str)
    l = len_trim(atmp)
    if (atmp(1:1) == "[".and.atmp(l:l) == "]") then
      value_isarray = .true.
    end if
  end function value_isarray
!========================================================================================!
  subroutine get_rawarray(kv)
    implicit none
    class(keyvalue) :: kv
    integer :: k,l,opend,closd
    integer :: n,b
    character(len=:),allocatable :: atmp,btmp
    character(len=1) :: s
    atmp = kv%rawvalue
    l = len_trim(atmp)
    opend = 0
    n = 0
    do k = 1,l
      s = atmp(k:k)
      if (s == '[') opend = opend+1
      if (s == ']') opend = opend-1
      if (opend == 1.and.s == ',') then
        n = n+1
      end if
    end do
    if (allocated(kv%value_rawa)) deallocate (kv%value_rawa)
    atmp = repeat(' ',l)
    kv%na = n+1
    allocate (kv%value_rawa(kv%na),source=atmp)
    atmp = kv%rawvalue
    n = 1
    opend = 0
    btmp = ''
    do k = 1,l
      s = atmp(k:k)
      if (s == ']') then
        opend = opend-1
      end if
      if (opend == 1.and.s == ',') then
        kv%value_rawa(n) = btmp
        n = n+1
        btmp = ''
        cycle
      end if
      if (opend >= 1) then
        btmp = btmp//s
      end if
      if (s == '[') then
        opend = opend+1
      end if
      if (k == l) then
        kv%value_rawa(n) = btmp
      end if
    end do
    deallocate (atmp,btmp)
    return
  end subroutine get_rawarray
  subroutine get_typearray(kv)
    implicit none
    class(keyvalue) :: kv
    integer :: k,l
    integer :: n,b,io
    real(wp) :: num
    character(len=:),allocatable :: atmp,btmp
    character(len=1) :: s
    logical :: allt

    !> check if all elements of the raw array fit a specific type
    kv%id = 9 !> unspecified array type
    !> int
    allt = .true.
    do k = 1,kv%na
      atmp = trim(adjustl(kv%value_rawa(k)))
      allt = allt.and.is_int(atmp)
    end do
    if (allt) then
      kv%id = 6
      allocate (kv%value_ia(kv%na))
      do k = 1,kv%na
        atmp = trim(adjustl(kv%value_rawa(k)))
        read (atmp,*,iostat=io) num
        kv%value_ia(k) = nint(num)
      end do
      return
    end if

    !> float
    allt = .true.
    do k = 1,kv%na
      atmp = trim(adjustl(kv%value_rawa(k)))
      allt = allt.and.(is_int(atmp).or.is_float(atmp))
    end do
    if (allt) then
      kv%id = 5
      allocate (kv%value_fa(kv%na))
      do k = 1,kv%na
        atmp = trim(adjustl(kv%value_rawa(k)))
        read (atmp,*,iostat=io) num
        kv%value_fa(k) = num
      end do
      return
    end if

    !> bool
    allt = .true.
    do k = 1,kv%na
      atmp = trim(adjustl(kv%value_rawa(k)))
      allt = allt.and.is_bool(atmp)
    end do
    if (allt) then
      kv%id = 7
      allocate (kv%value_ba(kv%na))
      do k = 1,kv%na
        atmp = trim(adjustl(kv%value_rawa(k)))
        if (atmp == 'true') kv%value_ba(k) = .true.
        if (atmp == 'false') kv%value_ba(k) = .false.
      end do
      return
    end if

  end subroutine get_typearray
!========================================================================================!
  subroutine get_valueautotype(kv)
    implicit none
    class(keyvalue) :: kv
    integer :: k,l
    integer :: io
    real(wp) :: num
    logical :: isnumber
    character(len=1),parameter :: floatchar(4) = [',','.','e','E']
    character(len=:),allocatable :: atmp,btmp
    atmp = trim(adjustl(kv%rawvalue))
    !> bool
    if (is_bool(atmp)) then
      kv%id = valuetypes%bool
      if (atmp == 'true') kv%value_b = .true.
      if (atmp == 'false') kv%value_b = .false.
      return
    end if

    if (is_float(atmp)) then
      !>--- real
      read (atmp,*,iostat=io) num
      kv%id = valuetypes%float
      kv%value_f = num
      kv%value_i = nint(num) !> backed up

    else if (is_int(atmp)) then
      !>--- integer
      read (atmp,*,iostat=io) num
      kv%id = valuetypes%int
      kv%value_i = nint(num)
      kv%value_f = real(num,wp) !> backed up

    end if

    return
  end subroutine get_valueautotype
  function is_bool(str)
    character(len=*) :: str
    logical :: is_bool
    character(len=:),allocatable :: atmp
    is_bool = .false.
    atmp = trim(adjustl(str))
    if (atmp == 'true'.or.atmp == 'false') then
      is_bool = .true.
    end if
    deallocate (atmp)
  end function is_bool
  function is_int(str)
    character(len=*) :: str
    logical :: is_int
    character(len=:),allocatable :: atmp
    integer :: io
    real(wp) :: num
    logical :: isnumber
    character(len=1),parameter :: floatchar(4) = [',','.','e','E']
    is_int = .false.
    atmp = trim(adjustl(str))
    read (atmp,*,iostat=io) num
    isnumber = io == 0
    if (isnumber) then
      is_int = .true.
      if (any(index(atmp,floatchar(:)) .ne. 0)) is_int = .false.
    end if
    deallocate (atmp)
  end function is_int
  function is_float(str)
    character(len=*) :: str
    logical :: is_float
    character(len=:),allocatable :: atmp
    integer :: io
    real(wp) :: num
    logical :: isnumber
    character(len=1),parameter :: floatchar(4) = [',','.','e','E']
    is_float = .false.
    atmp = trim(adjustl(str))
    read (atmp,*,iostat=io) num
    isnumber = io == 0
    if (isnumber) then
      if (any(index(atmp,floatchar(:)) .ne. 0)) is_float = .true.
    end if
    deallocate (atmp)
  end function is_float

!========================================================================================!

  subroutine kv_set_valuestring(self)
    implicit none
    class(keyvalue) :: self
    character(len=20) :: atmp
    character(len=:),allocatable :: btmp
    integer :: i
    if (allocated(self%rawvalue)) deallocate (self%rawvalue)
    select case (self%id)
    case (valuetypes%float) !> float
      write (atmp,'(f20.10)') self%value_f
      atmp = adjustl(atmp)
    case (valuetypes%int) !> integer
      write (atmp,'(i0)') self%value_i
      atmp = adjustl(atmp)
    case (valuetypes%bool) !> boolean
      if (self%value_b) then
        atmp = 'true'
      else
        atmp = 'false'
      end if
    case (valuetypes%string) !> string
      btmp = self%value_c

    case (valuetypes%float_array) !> float array
      btmp = '['
      do i = 1,self%na-1
        write (atmp,*) self%value_fa(i)
        btmp = trim(btmp)//trim(atmp)//','
      end do
      write (atmp,*) self%value_fa(self%na)
      btmp = trim(btmp)//trim(atmp)//']'

    case (valuetypes%int_array) !> int array
      btmp = '['
      do i = 1,self%na-1
        write (atmp,'(i0)') self%value_ia(i)
        btmp = trim(btmp)//trim(atmp)//','
      end do
      write (atmp,'(i0)') self%value_ia(self%na)
      btmp = trim(btmp)//trim(atmp)//']'

    case (valuetypes%bool_array) !> bool array
      btmp = '['
      do i = 1,self%na-1
        if (self%value_ba(i)) then
          atmp = 'true'
        else
          atmp = 'false'
        end if
        btmp = trim(btmp)//trim(atmp)//','
      end do
      if (self%value_ba(self%na)) then
        atmp = 'true'
      else
        atmp = 'false'
      end if
      btmp = trim(btmp)//trim(atmp)//']'

    case (valuetypes%string_array) !> multiline comment
      btmp = '['
      do i = 1,self%na-1
        btmp = trim(btmp)//'"'//trim(self%value_rawa(i))//'",'
      end do
      btmp = trim(btmp)//'"'//trim(self%value_rawa(self%Na))//'"]'

    case (valuetypes%raw_array) !> unspecified array
      btmp = '['
      do i = 1,self%na-1
        btmp = trim(btmp)//trim(self%value_rawa(i))//','
      end do
      btmp = trim(btmp)//trim(self%value_rawa(self%Na))//']'

    case default
      atmp = ''
    end select
    if (allocated(btmp)) then
      call move_alloc(btmp,self%rawvalue)
    else
      self%rawvalue = trim(atmp)
    end if
  end subroutine kv_set_valuestring

!========================================================================================!
  subroutine kv_add_raw_array_string(self,str)
    implicit none
    class(keyvalue) :: self
    character(len=*),intent(in) :: str
    character(len=:),allocatable :: newlist(:)
    integer :: k,i,j,maxwid
    i = self%na
    maxwid = 0
    if (i > 0) then
      do j = 1,i
        k = len_trim(self%value_rawa(j))
        if (k > maxwid) maxwid = k
      end do
    end if
    if (len_trim(str) > maxwid) then
      maxwid = len_trim(str)
    end if
    j = i+1
    allocate (newlist(j),source=repeat(' ',maxwid))
    newlist(1:i) = self%value_rawa(1:i)
    newlist(j) = str
    call move_alloc(newlist,self%value_rawa)
    self%na = j
  end subroutine kv_add_raw_array_string

!========================================================================================!
  function print_kv(self)
    implicit none
    character(len=:),allocatable :: print_kv
    class(keyvalue) :: self
    character(len=25) :: typ
    character(len=1),parameter :: ws = ' '
    character(len=100) :: atmp
    character(len=25) :: abbrval
    character(len=20) :: abbrkey
    character(len=*),parameter :: outfmt = &
    & '(1x,"key:",1x,a20,"value:"1x,a25,"type:",1x,a)'

    atmp = ''
    print_kv = atmp

    select case (self%id)
    case (valuetypes%float) !> float
      write (abbrval,'(e20.10)') self%value_f
      abbrval = adjustl(abbrval)
      typ = 'float'
    case (valuetypes%int) !> integer
      write (abbrval,'(i20)') self%value_i
      abbrval = adjustl(abbrval)
      typ = 'integer'
    case (valuetypes%bool) !> boolean
      if (self%value_b) then
        abbrval = 'true'
      else
        abbrval = 'false'
      end if
      typ = 'bool'
    case (valuetypes%string) !> string
      call truncate20(self%value_c,abbrval)
      if (index(abbrval,NEW_LINE('C')) .ne. 0) then
        abbrval = '<multiline string>'
      end if
      typ = 'string'
    case (valuetypes%float_array) !> float array
      call truncate20(self%rawvalue,abbrval)
      write (typ,'("array,len=",i0)') self%na
      typ = 'array(float)'
    case (valuetypes%int_array) !> int array
      call truncate20(self%rawvalue,abbrval)
      write (typ,'("array,len=",i0)') self%na
      typ = 'array(int)'
    case (valuetypes%bool_array) !> bool array
      call truncate20(self%rawvalue,abbrval)
      write (typ,'("array,len=",i0)') self%na
      typ = 'array(bool)'
    case (valuetypes%string_array) !> multiline comment
      call truncate20(self%rawvalue,abbrval)
      write (typ,'("array,len=",i0)') self%na
      typ = 'array(string)'
    case (valuetypes%raw_array) !> unspecified array
      call truncate20(self%rawvalue,abbrval)
      write (typ,'("array,len=",i0)') self%na
      typ = 'array(unspecified)'
    case default
      call truncate20(self%rawvalue,abbrval)
      typ = 'unspecified'
    end select
    call truncate15(self%key,abbrkey)
    write (atmp,outfmt) abbrkey,abbrval,trim(typ)
    print_kv = atmp
  end function print_kv

  subroutine truncate20(ins,outs)
    character(len=*) :: ins
    character(len=25) :: outs
    integer :: l
    !write(outs,'(a)') trim(ins)
    if (len_trim(ins) > 20) then
      outs(1:20) = ins(1:20)
      outs(21:25) = '...  '
    else
      l = len_trim(ins)
      outs = ins(1:l)
    end if
  end subroutine truncate20
  subroutine truncate15(ins,outs)
    character(len=*) :: ins
    character(len=20) :: outs
    integer :: l
    !write(outs,'(a)') trim(ins)
    if (len_trim(ins) > 15) then
      outs(1:15) = ins(1:20)
      outs(16:20) = '...  '
    else
      l = len_trim(ins)
      outs = ins(1:l)
    end if
  end subroutine truncate15

  function print_kv2(self) result(outstr)
    implicit none
    character(len=:),allocatable :: outstr
    class(keyvalue) :: self
    select case (self%id)
    case (valuetypes%string)
      outstr = self%key//' = "'//self%rawvalue//'"'
    case default
      outstr = self%key//' = '//self%rawvalue
    end select
  end function print_kv2

!========================================================================================!
end module parse_keyvalue
