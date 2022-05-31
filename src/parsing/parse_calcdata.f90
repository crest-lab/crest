
module parse_calcdata
  use iso_fortran_env,only:wp => real64,sp => real64
  use calc_type,only:calcdata,calculation_settings
  use constraints
  use dynamics_module
  use metadynamics_module

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

  interface parse_md
    module procedure :: parse_md_auto
    module procedure :: parse_md_float
    module procedure :: parse_md_int
    module procedure :: parse_md_c
    module procedure :: parse_md_bool
  end interface parse_md

  interface parse_mtd
    module procedure :: parse_metadyn_auto
    module procedure :: parse_mtd_float
    module procedure :: parse_mtd_int
    module procedure :: parse_mtd_c
    module procedure :: parse_mtd_bool
  end interface parse_mtd

  public :: parse_calculation_data
  public :: parse_dynamics_data
contains
!========================================================================================!
  subroutine parse_calculation_data(calc,dict,included)
    implicit none
    type(calcdata) :: calc
    type(root_object) :: dict
    type(datablock) :: blk
    type(calculation_settings) :: newjob,newjob2
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
      else if (blk%header == '[calculation.mecp]') then
        !>-- setup
        if(allocated(calc%calcs))deallocate(calc%calcs)
        calc%ncalculations = 0
        calc%id = -1 
        call parse_leveldata(blk,newjob)
        !>-- S0 setup
        call parse_leveldata(blk,newjob)
        newjob%uhf = 0
        newjob%calcspace = 's0'
        call calc%add(newjob)
        !>-- S1 setup
        newjob%uhf = 2
        newjob%calcspace = 's1'
        call calc%add(newjob)
        included = .true.
      else if (blk%header == '[calculation.constraints]') then
        call parse_constraintdat(blk,calc)
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
    if ((blk%header .ne. '[calculation.level]') .and. &
    & (blk%header .ne. '[calculation.mecp]')) then
       return
    endif
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
    case( 'type' )
      select case( val )
      case( 'mecp' )
        calc%id = -1  
      case default
        calc%id = 1
      end select
    case ('elog')
      calc%elog = val
      calc%pr_energies = .true.
    case ('hess_update','hupdate')
      select case(val)
        case( 'bfgs' )
          calc%iupdat = 0
        case( 'powell' )
          calc%iupdat = 1
        case( 'sr1' )
          calc%iupdat = 2
        case( 'bofill' )
          calc%iupdat = 3
        case( 'schlegel' )
          calc%iupdat = 4
        case default
          calc%iupdat = 0
      end select
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
!> The following routines are used to
!> read information into the "constraint" object
!> and add it to a calculation data object
!>---------------------------------------------------------
  subroutine parse_constraintdat(blk,calc)
    implicit none
    type(datablock),intent(in) :: blk
    type(calcdata),intent(inout) :: calc
    logical :: success
    type(constraint) :: constr
    integer :: i
    if (blk%header .ne. '[calculation.constraints]') return
    do i = 1,blk%nkv
      call parse_constraint_auto(constr,blk%kv_list(i),success)
      if (success) then
        call calc%add(constr)
      end if
    end do
    return
  end subroutine parse_constraintdat
  subroutine parse_constraint_auto(constr,kv,success)
    implicit none
    type(keyvalue) :: kv
    type(constraint) :: constr
    logical,intent(out) :: success
    real(wp) :: dum1,dum2,dum3
    success = .false.
    select case (kv%key)
    case ('bond','bonds')
      select case (kv%id)
      case (4) !> string
        select case (kv%value_c)
        case ('all','allauto')
          call constr%dummyconstraint(11)
          success = .true.
        end select
      case (9) !> unspecified array
        call constr%analyzedummy(11,kv%na,kv%value_rawa)
        success = .true.
      end select
    case ('sphere')
      dum1 = kv%value_fa(3)  !> sphere radius
      dum2 = kv%value_fa(1)  !> prefactor
      dum3 = kv%value_fa(2)  !> exponent
      call constr%sphereconstraint(0,dum1,dum2,dum3,.false.)
      success = .true.
    case ('sphere_logfermi')
      dum1 = kv%value_fa(3)  !> sphere radius
      dum2 = kv%value_fa(1)  !> fermi temperature
      dum3 = kv%value_fa(2)  !> exponent factor
      call constr%sphereconstraint(0,dum1,dum2,dum3,.true.)
      success = .true.
    case ('gapdiff')
      dum1 = kv%value_fa(1)
      dum2 = kv%value_fa(2)
      call constr%gapdiffconstraint(dum1,dum2)
      success = .true.
    case ('gapdiff2', 'mecp')
      success = .true.
      if(kv%id==3)then
        if(kv%value_b)then
          dum1 = 10.0_wp
          dum2 = 0.005_wp
          dum3 = 0.20_wp
        else
          success = .false. 
        endif 
      else
      dum1 = kv%value_fa(1)
      dum2 = kv%value_fa(2)
      dum3 = kv%value_fa(3)
      endif
      call constr%gapdiffconstraint2(dum1,dum2,dum3)
    case default
      return
    end select

    return
  end subroutine parse_constraint_auto

!========================================================================================!
  subroutine parse_dynamics_data(mddat,dict,included)
    implicit none
    type(mddata) :: mddat
    type(root_object) :: dict
    type(datablock) :: blk
    type(calculation_settings) :: newjob
    type(constraint) :: newcstr
    integer :: i,j,k,l
    logical,intent(out) :: included

    included = .false.
    !call mddat%deallocate()

    do i = 1,dict%nblk
      call blk%deallocate()
      blk = dict%blk_list(i)
      if (blk%header == 'dynamics') then
        included = .true.
        call parse_mddat(blk,mddat)
      else if (blk%header == '[dynamics.meta]') then
        call parse_metadyn(blk,mddat)
        !call calc%add(newjob)
        included = .true.
      end if
    end do
    if (included) then
      mddat%requested = .true.
    end if
    return
  end subroutine parse_dynamics_data
  subroutine parse_mddat(blk,mddat)
    implicit none
    type(datablock),intent(in) :: blk
    type(mddata),intent(inout) :: mddat
    integer :: i
    if (blk%header .ne. 'dynamics') return
    do i = 1,blk%nkv
      call parse_md(mddat,blk%kv_list(i))
    end do
    return
  end subroutine parse_mddat
  subroutine parse_md_auto(mddat,kv)
    implicit none
    type(mddata) :: mddat
    type(keyvalue) :: kv
    select case (kv%id)
    case (1) !> float
      call parse_md(mddat,kv%key,kv%value_f)
    case (2) !> int
      call parse_md(mddat,kv%key,kv%value_i)
    case (3) !> bool
      call parse_md(mddat,kv%key,kv%value_b)
    case (4) !> string
      call parse_md(mddat,kv%key,kv%value_c)
    end select
  end subroutine parse_md_auto
  subroutine parse_md_float(mddat,key,val)
    implicit none
    type(mddata) :: mddat
    character(len=*) :: key
    real(wp) :: val
    select case (key)
    case ('length')
      mddat%length_ps = val
    case ('dump')
      mddat%dumpstep = val
    case ('hmass')
      mddat%md_hmass = val
    case ('tstep')
      mddat%tstep = val
    case default
      return
    end select
    return
  end subroutine parse_md_float
  subroutine parse_md_int(mddat,key,val)
    implicit none
    type(mddata) :: mddat
    character(len=*) :: key
    integer :: val
    real(wp) :: fval
    select case (key)
    case ('length','dump','hmass','tstep')
      fval = float(val)
      call parse_md(mddat,key,fval)
    case default
      return
    end select
    return
  end subroutine parse_md_int
  subroutine parse_md_c(mddat,key,val)
    implicit none
    type(mddata) :: mddat
    character(len=*) :: key
    character(len=*) :: val
    select case (key)
    case default
      return
    end select
    return
  end subroutine parse_md_c
  subroutine parse_md_bool(mddat,key,val)
    implicit none
    type(mddata) :: mddat
    character(len=*) :: key
    logical :: val
    select case (key)
    case default
      return
    end select
    return
  end subroutine parse_md_bool
!========================================================================================!
!> The following routines are used to
!> read information into the "metadynamics" object
!> and add it to a mol.dynamics data object
!>---------------------------------------------------------
  subroutine parse_metadyn(blk,mddat)
    implicit none
    type(datablock),intent(in) :: blk
    type(mddata),intent(inout) :: mddat
    logical :: success
    type(mtdpot) :: mtd
    integer :: i,k
    call mtd%deallocate()
    if (blk%header .ne. '[dynamics.meta]') return
    do i = 1,blk%nkv
      call parse_metadyn_auto(mtd,blk%kv_list(i),success)
    end do
    call mddat%add(mtd)
    return
  end subroutine parse_metadyn
  subroutine parse_metadyn_auto(mtd,kv,success)
    implicit none
    type(keyvalue) :: kv
    type(mtdpot) :: mtd
    logical,intent(out) :: success
    success = .false.
    select case (kv%id)
    case (1) !> float
      call parse_mtd(mtd,kv%key,kv%value_f)
      success = .true.
    case (2) !> int
      call parse_mtd(mtd,kv%key,kv%value_i)
      success = .true.
    case (3) !> bool
      call parse_mtd(mtd,kv%key,kv%value_b)
      success = .true.
    case (4) !> string
      call parse_mtd(mtd,kv%key,kv%value_c)
      success = .true.
    end select
  end subroutine parse_metadyn_auto
  subroutine parse_mtd_float(mtd,key,val)
    implicit none
    type(mtdpot) :: mtd
    character(len=*) :: key
    real(wp) :: val
    select case (key)
    case ('alpha')
      mtd%alpha = val
    case ('kpush')
      mtd%kpush = val
    case ('dump','dump_fs')
      mtd%cvdump_fs = val
    case ('dump_ps')
      mtd%cvdump_fs = val * 1000.0_wp
    case default
      return
    end select
    return
  end subroutine parse_mtd_float
  subroutine parse_mtd_int(mtd,key,val)
    implicit none
    type(mtdpot) :: mtd
    character(len=*) :: key
    integer :: val
    real(wp) :: fval
    select case (key)
    case ('type')
      mtd%mtdtype = val
    case ('dump','dump_fs','dump_ps')
      fval = float(val)
      call parse_mtd(mtd,key,fval)
    case default
      return
    end select
    return
  end subroutine parse_mtd_int
  subroutine parse_mtd_c(mtd,key,val)
    implicit none
    type(mtdpot) :: mtd
    character(len=*) :: key
    character(len=*) :: val
    select case (key)
    case ('type')
      select case (val)
      case ('rmsd','RMSD')
        mtd%mtdtype = cv_rmsd
      case default
        mtd%mtdtype = 0
      end select
    case default
      return
    end select
    return
  end subroutine parse_mtd_c
  subroutine parse_mtd_bool(mtd,key,val)
    implicit none
    type(mtdpot) :: mtd
    character(len=*) :: key
    logical :: val
    select case (key)
    case default
      return
    end select
    return
  end subroutine parse_mtd_bool

!========================================================================================!
end module parse_calcdata
