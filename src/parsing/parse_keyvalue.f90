
module parse_keyvalue
  use iso_fortran_env,only:wp => real64,sp => real32
  implicit none

  public :: keyvalue
  type :: keyvalue
    character(len=:),allocatable :: key
    character(len=:),allocatable :: rawvalue
    integer :: id = 0
    real(wp) :: value_f   !> id=1
    integer  :: value_i   !> id=2
    logical  :: value_b   !> id=3
    character(len=:),allocatable :: value_c !> id=4
    integer ::  na
    character(len=:),allocatable :: value_rawa(:)
    real(wp),allocatable :: value_fa(:) !> id=5
    integer,allocatable  :: value_ia(:) !> id=6
    logical,allocatable  :: value_ba(:) !> id=7
    character(len=:),allocatable :: value_cml(:) !> id=8
  contains
    procedure :: print => print_kv
    procedure :: deallocate => deallocate_kv
  end type keyvalue

  character(len=*),parameter :: kv_indicator = '='

  public :: get_keyvalue

contains
  subroutine deallocate_kv(self)
    implicit none
    class(keyvalue) :: self
    if (allocated(self%key)) deallocate (self%key)
    if (allocated(self%rawvalue)) deallocate (self%rawvalue)
    self%id = 0
    if (allocated(self%value_c)) deallocate (self%value_c)
    self%na = 0
    if (allocated(self%value_c)) deallocate (self%value_rawa)
    if (allocated(self%value_c)) deallocate (self%value_fa)
    if (allocated(self%value_c)) deallocate (self%value_ia)
    if (allocated(self%value_c)) deallocate (self%value_ba)
    if (allocated(self%value_c)) deallocate (self%value_cml)
    return
  end subroutine
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
    if (.not. iskeyvalue(str)) then
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
    if(len_trim(atmp) < 1) return
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
    l = k - 1
    if (atmp(1:1) == '"' .or. atmp(1:1) == "'") then
      l = k - 2
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
    k = k + len(kv_indicator)
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
    if(len_trim(atmp) < 3)return
    if (atmp(1:3) == "'''" .or. atmp(1:3) == '"""') then
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
    if (.not. value_ismultiline(str)) then
      if ((atmp(1:1) == "'" .and. atmp(l:l) == "'") .or. &
      &  (atmp(1:1) == '"' .and. atmp(l:l) == '"')) then
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
    l = len_trim(atmp) - 1
    kv%value_c = atmp(2:l)
    kv%id = 4
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
    if (atmp(1:1) == "[" .and. atmp(l:l) == "]") then
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
      if (s == '[') opend = opend + 1
      if (s == ']') opend = opend - 1
      if (opend == 1 .and. s == ',') then
        n = n + 1
      end if
    end do
    if (allocated(kv%value_rawa)) deallocate (kv%value_rawa)
    atmp = repeat(' ',l)
    kv%na = n + 1
    allocate (kv%value_rawa(kv%na),source=atmp)
    atmp = kv%rawvalue
    n = 1
    opend = 0
    btmp = ''
    do k = 1,l
      s = atmp(k:k)
      if (s == ']') then
        opend = opend - 1
      end if
      if (opend == 1 .and. s == ',') then
        kv%value_rawa(n) = btmp
        n = n + 1
        btmp = ''
        cycle
      end if
      if (opend >= 1) then
        btmp = btmp//s
      end if
      if (s == '[') then
        opend = opend + 1
      end if
      if (k == l) then
        kv%value_rawa(n) = btmp
      end if
    end do
    deallocate(atmp,btmp)
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
    
    !> int
    allt = .true.
    do k=1,kv%na
      atmp = trim(adjustl(kv%value_rawa(k)))
      allt = allt .and. is_int(atmp)
    enddo
    if(allt)then
      kv%id = 6
      allocate(kv%value_ia(kv%na))
      do k=1,kv%na
      atmp = trim(adjustl(kv%value_rawa(k)))
      read (atmp,*,iostat=io) num
      kv%value_ia(k) = nint(num)
      enddo
      return
    endif

    !> float
    allt = .true.
    do k=1,kv%na
      atmp = trim(adjustl(kv%value_rawa(k)))
      allt = allt .and. (is_int(atmp) .or. is_float(atmp))
    enddo
    if(allt)then
      kv%id = 5
      allocate(kv%value_fa(kv%na))
      do k=1,kv%na
      atmp = trim(adjustl(kv%value_rawa(k)))
      read (atmp,*,iostat=io) num
      kv%value_fa(k) = num
      enddo
      return
    endif

    !> bool
    allt = .true.
    do k=1,kv%na
      atmp = trim(adjustl(kv%value_rawa(k)))
      allt = allt .and. is_bool(atmp)   
    enddo 
    if(allt)then
      kv%id = 7
      allocate(kv%value_ba(kv%na))
      do k=1,kv%na
      atmp = trim(adjustl(kv%value_rawa(k)))
      if (atmp == 'true') kv%value_ba(k) = .true.
      if (atmp == 'false') kv%value_ba(k) = .false.
      enddo
      return 
    endif

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
      kv%id = 3
      if (atmp == 'true') kv%value_b = .true.
      if (atmp == 'false') kv%value_b = .false.
      return
    end if
    !> real
    if (is_float(atmp)) then
      read (atmp,*,iostat=io) num
      kv%id = 1
      kv%value_f = num

      !> integer
    else if (is_int(atmp)) then
      read (atmp,*,iostat=io) num
      kv%id = 2
      kv%value_i = nint(num)
    end if

    return
  end subroutine get_valueautotype
  function is_bool(str)
    character(len=*) :: str
    logical :: is_bool
    character(len=:),allocatable :: atmp
    is_bool = .false.
    atmp = trim(adjustl(str))
    if (atmp == 'true' .or. atmp == 'false') then
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
    case (1) !> float
      write (abbrval,'(f20.10)') self%value_f
      abbrval = adjustl(abbrval)
      typ = 'float'
    case (2) !> integer
      write (abbrval,'(i20)') self%value_i
      abbrval = adjustl(abbrval)
      typ = 'integer'
    case (3) !> boolean
      if (self%value_b) then
        abbrval = 'true'
      else
        abbrval = 'false'
      end if
      typ = 'bool'
    case (4) !> string
      call truncate20(self%value_c,abbrval)
      typ = 'string'
    case (5) !> float array
      call truncate20(self%rawvalue,abbrval)
      write (typ,'("array,len=",i0)') self%na
      typ = 'array(float)'
    case (6) !> int array
      call truncate20(self%rawvalue,abbrval)
      write (typ,'("array,len=",i0)') self%na
      typ = 'array(int)'
    case (7) !> bool array
      call truncate20(self%rawvalue,abbrval)
      write (typ,'("array,len=",i0)') self%na
      typ = 'array(bool)'
    case (8) !> multiline comment
      abbrval = "'''...'''"
      typ = 'multiline comment'
    case default
      call truncate20(self%rawvalue,abbrval)
      typ = 'unspecified'
    end select
    call truncate15(self%key,abbrkey)
    write (atmp,outfmt) abbrkey,abbrval,trim(typ)
    print_kv = atmp
  contains
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
  end function print_kv
!========================================================================================!
end module parse_keyvalue
