
module parse_calcdata
  use iso_fortran_env,only:wp => real64,sp => real64
  use calc_type,only:calcdata,calculation_settings
  use constraints,only:constraint

  use parse_block,only:datablock
  use parse_keyvalue,only:keyvalue
  use parse_datastruct,only:root_object

  implicit none
  private

  interface parse_setting
    module procedure :: parse_setting_auto
    module procedure :: parse_setting_float
    module procedure :: parse_setting_int
    module procedure :: parse_setting_c
    module procedure :: parse_setting_bool
  end interface parse_setting

  interface parse_calc
    module procedure :: parse_calc_auto
    module procedure :: parse_calc_float
    module procedure :: parse_calc_int
    module procedure :: parse_calc_c
    module procedure :: parse_calc_bool
  end interface parse_calc

  public :: parse_calculation_data
contains
!========================================================================================!
  subroutine parse_calculation_data(calc,dict,included)
    implicit none
    type(calcdata) :: calc
    type(root_object) :: dict
    type(datablock) :: blk
    type(calculation_settings) :: newjob
    type(constraint) :: newcstr
    integer :: i,j,k,l
    logical,intent(out) :: included

    included = .false.
    call calc%reset()

    do i = 1,dict%nblk
      call blk%deallocate()
      blk = dict%blk_list(i)
      if (blk%header == 'calculation') then
        included = .true.
        call parse_calcdat(blk,calc)
        !write(*,*) 'read [calculation]'
      else if (blk%header == '[calculation.level]') then
        call parse_leveldata(blk,newjob)
        call calc%add(newjob)
        included = .true.
        !write(*,*) 'read [calculation.level]'
      else if (blk%header == '[calculation.constraints]') then

        included = .true.
        !write(*,*) 'read [calculation.constraints]'
      end if
    end do
    if (included) then
      call calc%init()
    end if
    return
  end subroutine parse_calculation_data

!========================================================================================!
!> The following routines are used to
!> read information into the "calculation_settings" object
!>---------------------------------------------------------
  subroutine parse_leveldata(blk,job)
    implicit none
    type(datablock),intent(in) :: blk
    type(calculation_settings),intent(out) :: job
    integer :: i
    call job%deallocate()
    if (blk%header .ne. '[calculation.level]') return
    do i = 1,blk%nkv
      call parse_setting(job,blk%kv_list(i))
    end do
    return
  end subroutine parse_leveldata
  subroutine parse_setting_auto(job,kv)
    implicit none
    type(calculation_settings) :: job
    type(keyvalue) :: kv
    select case (kv%id)
    case (1) !> float
      call parse_setting(job,kv%key,kv%value_f)
    case (2) !> int
      call parse_setting(job,kv%key,kv%value_i)
    case (3) !> bool
      call parse_setting(job,kv%key,kv%value_b)
    case (4) !> string
      call parse_setting(job,kv%key,kv%value_c)
    end select
  end subroutine parse_setting_auto
  subroutine parse_setting_float(job,key,val)
    implicit none
    type(calculation_settings) :: job
    character(len=*) :: key
    real(wp) :: val
    !> calculation_settings actually has no
    !> float arguments to be read
    !> but maybe some were mistyped as such
    select case (key)
    case ('uhf')
      job%uhf = nint(val)
    case ('chrg','charge')
      job%chrg = nint(val)
    end select
    return
  end subroutine parse_setting_float
  subroutine parse_setting_int(job,key,val)
    implicit none
    type(calculation_settings) :: job
    character(len=*) :: key
    integer :: val
    select case (key)
    case ('uhf')
      job%uhf = val
    case ('chrg','charge')
      job%chrg = val
    case ('id')
      job%id = val
    end select
    return
  end subroutine parse_setting_int
  subroutine parse_setting_c(job,key,val)
    implicit none
    type(calculation_settings) :: job
    character(len=*) :: key
    character(len=*) :: val
    select case (key)
    case ('bin','binary')
      job%binary = val
    case ('flags')
      job%other = val
      !> don't.
      !case ('sys','syscall','systemcall')
      !  job%systemcall = val
    case ('calcspace','dir')
      job%calcspace = val
    case ('method')
      select case (val)
      case ('gfn-xtb','gfn','xtb')
        job%id = 10
      case default
        job%id = 0
      end select
    end select
    return
  end subroutine parse_setting_c
  subroutine parse_setting_bool(job,key,val)
    implicit none
    type(calculation_settings) :: job
    character(len=*) :: key
    logical :: val
    select case (key)
    case ('rdwbo')
      job%rdwbo = val
    case ('rddip','rddipole')
      job%rddip = val
    case ('dipgrad')
      job%rddipgrad = val
    end select
    return
  end subroutine parse_setting_bool
!========================================================================================!
!> The following routines are used to
!> read information into the "calcdata" object
!>---------------------------------------------------------
  subroutine parse_calcdat(blk,calc)
    implicit none
    type(datablock),intent(in) :: blk
    type(calcdata),intent(inout) :: calc
    integer :: i
    if (blk%header .ne. 'calculation') return
    do i = 1,blk%nkv
      call parse_calc(calc,blk%kv_list(i))
    end do
    return
  end subroutine parse_calcdat
  subroutine parse_calc_auto(calc,kv)
    implicit none
    type(calcdata) :: calc
    type(keyvalue) :: kv
    select case (kv%id)
    case (1) !> float
      call parse_calc(calc,kv%key,kv%value_f)
    case (2) !> int
      call parse_calc(calc,kv%key,kv%value_i)
    case (3) !> bool
      call parse_calc(calc,kv%key,kv%value_b)
    case (4) !> string
      call parse_calc(calc,kv%key,kv%value_c)
    end select
  end subroutine parse_calc_auto
  subroutine parse_calc_float(calc,key,val)
    implicit none
    type(calcdata) :: calc
    character(len=*) :: key
    real(wp) :: val
    select case (key)
    case default
      return
    end select
    return
  end subroutine parse_calc_float
  subroutine parse_calc_int(calc,key,val)
    implicit none
    type(calcdata) :: calc
    character(len=*) :: key
    integer :: val
    select case (key)
    case ('id','type')
      calc%id = val
    case default
      return
    end select
    return
  end subroutine parse_calc_int
  subroutine parse_calc_c(calc,key,val)
    implicit none
    type(calcdata) :: calc
    character(len=*) :: key
    character(len=*) :: val
    select case (key)
    case('elog')
      calc%elog = val
      calc%pr_energies = .true.
    case default
      return
    end select
    return
  end subroutine parse_calc_c
  subroutine parse_calc_bool(calc,key,val)
    implicit none
    type(calcdata) :: calc
    character(len=*) :: key
    logical :: val
    select case (key)
    case ('eprint')
      calc%pr_energies = val
    case default
      return
    end select
    return
  end subroutine parse_calc_bool

!========================================================================================!
end module parse_calcdata
