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

!> NOTE: This is work in progress, not all input conventions have been set yet
!========================================================================================!
!> Routines contained here are for parsing calculation and simulation settings that will 
!> enter separate setup objects.
!> This concerns mainly [calculation]/[calculation. ...] and [dynamics] blocks

module parse_calcdata
  use crest_parameters
  use crest_calculator,only:calcdata,calculation_settings,jobtype,constraint,scantype
  use dynamics_module
  use gradreader_module,only:gradtype,conv2gradfmt
  use tblite_api,only:xtblvl

  use parse_block,only:datablock
  use parse_keyvalue,only:keyvalue
  use parse_datastruct,only:root_object

  implicit none
  private

!>-- routines for parsing a calculation_settings object
  interface parse_setting  
    module procedure :: parse_setting_auto
    module procedure :: parse_setting_float
    module procedure :: parse_setting_int
    module procedure :: parse_setting_c
    module procedure :: parse_setting_bool
  end interface parse_setting

!>-- routines for parsing a calcdata object
  interface parse_calc
    module procedure :: parse_calc_auto
    module procedure :: parse_calc_float
    module procedure :: parse_calc_int
    module procedure :: parse_calc_c
    module procedure :: parse_calc_bool
  end interface parse_calc


!>-- routines for parsing a mddata object
  interface parse_md
    module procedure :: parse_md_auto
    module procedure :: parse_md_float
    module procedure :: parse_md_int
    module procedure :: parse_md_c
    module procedure :: parse_md_bool
  end interface parse_md

!>-- routines for parsing a mtdpot object
  interface parse_mtd
    module procedure :: parse_metadyn_auto
    module procedure :: parse_mtd_float
    module procedure :: parse_mtd_int
    module procedure :: parse_mtd_c
    module procedure :: parse_mtd_bool
  end interface parse_mtd

  public :: parse_calculation_data
  public :: parse_dynamics_data

!========================================================================================!
!========================================================================================!
contains !> MODULE PROCEDURES START HERE
!========================================================================================!
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

      else if (blk%header == 'calculation.level') then
        call parse_leveldata(blk,newjob)
        call newjob%autocomplete(calc%ncalculations+1)
        call calc%add(newjob)
        included = .true.

      else if (blk%header == 'calculation.mecp') then
        !>-- setup
        if (allocated(calc%calcs)) deallocate (calc%calcs)
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

      else if (blk%header == 'calculation.constraints') then
        call parse_constraintdat(blk,calc)
        included = .true.

      else if (blk%header == 'calculation.scans') then
        call parse_scandat(blk,calc)
        included = .true.

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
    if ((blk%header .ne. 'calculation.level').and. &
    & (blk%header .ne. 'calculation.mecp')) then
      return
    end if
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
    case (6,7) !> int/float array
      call parse_setting_array(job,kv,kv%key)
    end select
  end subroutine parse_setting_auto
  subroutine parse_setting_float(job,key,val)
    implicit none
    type(calculation_settings) :: job
    character(len=*) :: key
    real(wp) :: val
    select case (key)
    case ('uhf')
      job%uhf = nint(val)
    case ('chrg','charge')
      job%chrg = nint(val)
    case ('etemp')
      job%etemp = val
    case ('accuracy')
      job%accuracy = val
    case ('weight')
      job%weight = val
    case ('pressure')
      job%extpressure = val
    case ('proberad')
      job%proberad = val
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
    case ('maxscc')
      job%maxscc = val
    case ('tblite_level','tblite_hamiltonian')
      job%tblitelvl = val
    case ('lebedev')
      job%ngrid = val
    case('vdwSet')
      job%vdwset = val
    end select
    return
  end subroutine parse_setting_int
  subroutine parse_setting_c(job,key,val)
    implicit none
    type(calculation_settings) :: job
    character(len=*) :: key
    character(len=*) :: val
    select case (key)

    case ('method')
      select case (val)
      case ('gfn-xtb','gfn','xtb')
        job%id = jobtype%xtbsys
      case ('generic')
        job%id = jobtype%generic
      case ('orca')
        job%id = jobtype%orca
      case ('turbomole','tm')
        job%id = jobtype%turbomole
      case ('terachem')
        job%id = jobtype%terachem
      case ('tblite')
        job%id = jobtype%tblite
      case ('gfn0','gfn0-xtb')
        job%id = jobtype%gfn0
      case ('gfn0*','gfn0*-xtb')
        job%id = jobtype%gfn0occ
      case ('gfnff','gff','gfn-ff')
        job%id = jobtype%gfnff
      case ('xhcff')
        job%id = jobtype%xhcff
      case ('none')
        job%id = jobtype%unknown
      case default
        job%id = jobtype%unknown
      end select

    case ('bin','binary','script')
      job%binary = val

    case ('flags')
      job%other = val

      !> don't.
      !case ('sys','syscall','systemcall')
      !  job%systemcall = val

    case ('calcspace','dir')
      job%calcspace = val

    case ('gradfile')
      job%gradfile = val

    case ('gradtype')
      select case (val)
      case ('engrad','xtb','orca')
        job%gradtype = gradtype%engrad
      case ('turbomole','tm')
        job%gradtype = gradtype%turbomole
      case ('generic')
        job%gradtype = gradtype%unknown
      case default
        job%gradtype = gradtype%unknown
      end select

    case ('gradkey')
      job%gradkey = val

    case ('gradmt')
      job%gradfmt = conv2gradfmt(val)

    case ('efile')
      job%efile = val

    case ('tblite_level','tblite_hamiltonian')
      select case (val)
      case ('gfn2','gfn2-xtb')
        job%tblitelvl = xtblvl%gfn2
      case ('gfn1','gfn1-xtb')
        job%tblitelvl = xtblvl%gfn1
      case ('ipea1')
        job%tblitelvl = xtblvl%ipea1
      case default
        job%tblitelvl = xtblvl%unknown
      end select

    case ('gbsa','alpb','cpcm')
      job%solvmodel = key
      job%solvent = val

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
    case ('refresh')
      job%apiclean = val
    case ('print')
      if (val) job%prch = 999  !> the actual ID will be overwritten automatically
    end select
    return
  end subroutine parse_setting_bool
  subroutine parse_setting_array(job,kv,key)
    implicit none
    type(calculation_settings) :: job
    type(keyvalue) :: kv
    character(len=*) :: key
    select case (key)
    case ('config')
      call job%addconfig(kv%value_ia)
    end select
    return
  end subroutine parse_setting_array
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
    case ('maxcycle')
      calc%maxcycle = val
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
    case ('type')
      select case (val)
      case ('mecp')
        calc%id = -1
      case default
        calc%id = 1
      end select
    case ('elog')
      calc%elog = val
      calc%pr_energies = .true.
    case ('hess_update','hupdate')
      select case (val)
      case ('bfgs')
        calc%iupdat = 0
      case ('powell')
        calc%iupdat = 1
      case ('sr1')
        calc%iupdat = 2
      case ('bofill')
        calc%iupdat = 3
      case ('schlegel')
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
    if (blk%header .ne. 'calculation.constraints') return
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
    integer :: atm1,atm2,atm3,atm4
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
      case(5)  !> regular array
        call constr%rdbondconstraint(kv%na,kv%value_fa)
        success = .true.
      case (9) !> unspecified array
        call constr%analyzedummy(11,kv%na,kv%value_rawa)
        success = .true.
      case default
        success = .false.
      end select
    case ('dihedral')
      read (kv%value_rawa(1),*) atm1
      read (kv%value_rawa(2),*) atm2
      read (kv%value_rawa(3),*) atm3
      read (kv%value_rawa(4),*) atm4
      read (kv%value_rawa(5),*) dum1
      if (kv%na > 5) then
        read (kv%value_rawa(6),*) dum2
        call constr%dihedralconstraint(atm1,atm2,atm3,atm4,dum1,dum2)
      else
        call constr%dihedralconstraint(atm1,atm2,atm3,atm4,dum1)
      end if
      success = .true.
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
    case ('gapdiff2','mecp')
      success = .true.
      if (kv%id == 3) then
        if (kv%value_b) then
          dum1 = 10.0_wp
          dum2 = 0.005_wp
          dum3 = 0.20_wp
        else
          success = .false.
        end if
      else
        dum1 = kv%value_fa(1)
        dum2 = kv%value_fa(2)
        dum3 = kv%value_fa(3)
      end if
      call constr%gapdiffconstraint2(dum1,dum2,dum3)
    case default
      return
    end select

    return
  end subroutine parse_constraint_auto

!========================================================================================!
!> The following routines are used to
!> read information into the "scan" object
!> and add it to a calculation data object
!>---------------------------------------------------------
  subroutine parse_scandat(blk,calc)
    implicit none
    type(datablock),intent(in) :: blk
    type(calcdata),intent(inout) :: calc
    logical :: success
    type(scantype) :: scn
    integer :: i
    if (blk%header .ne. 'calculation.scans') return
    do i = 1,blk%nkv
      call parse_scan_auto(scn,blk%kv_list(i),success)
      if (success) then
        call calc%add(scn)
      end if
    end do
    return
  end subroutine parse_scandat
  subroutine parse_scan_auto(scn,kv,success)
    implicit none
    type(keyvalue) :: kv
    type(scantype) :: scn
    logical,intent(out) :: success
    real(wp) :: dum1,dum2,dum3
    integer :: atm1,atm2,atm3,atm4
    integer :: nsteps
    success = .false.
    call scn%deallocate()
    select case (kv%key)
    case ('bond','distance')
      scn%type = 1
      scn%n = 2
      allocate (scn%atms(2))
      read (kv%value_rawa(1),*) atm1
      scn%atms(1) = atm1
      read (kv%value_rawa(2),*) atm2
      scn%atms(2) = atm2
      read (kv%value_rawa(3),*) dum1
      scn%minval = dum1
      read (kv%value_rawa(4),*) dum2
      scn%maxval = dum2
      if (kv%na > 4) then
        read (kv%value_rawa(5),*) nsteps
        scn%steps = nsteps
      end if
      success = .true.
    case ('dihedral')
      scn%type = 3
      scn%n = 2
      allocate (scn%atms(4))
      read (kv%value_rawa(1),*) atm1
      scn%atms(1) = atm1
      read (kv%value_rawa(2),*) atm2
      scn%atms(2) = atm2
      read (kv%value_rawa(3),*) atm3
      scn%atms(3) = atm3
      read (kv%value_rawa(4),*) atm4
      scn%atms(4) = atm4
      if (kv%na > 4) then
        read (kv%value_rawa(5),*) nsteps
        scn%steps = nsteps
      end if
      if (kv%na > 6) then
        read (kv%value_rawa(6),*) dum1
        scn%minval = dum1
        read (kv%value_rawa(7),*) dum2
        scn%maxval = dum2
      end if
      success = .true.
    case default
      return
    end select

    return
  end subroutine parse_scan_auto

!========================================================================================!
!> The following routines are used to
!> read information into the "mddata" object
!>---------------------------------------------------------
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

    do i = 1,dict%nblk
      call blk%deallocate()
      blk = dict%blk_list(i)
      if (blk%header == 'dynamics') then
        included = .true.
        call parse_mddat(blk,mddat)
      else if (blk%header == 'dynamics.meta') then
        call parse_metadyn(blk,mddat)
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
    case ('length','length_ps')
      mddat%length_ps = val
    case ('dump')
      mddat%dumpstep = val
    case ('hmass')
      mddat%md_hmass = val
    case ('tstep')
      mddat%tstep = val
    case ('t','temp','temperature' )
      mddat%tsoll = val
      mddat%thermostat =.true.
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
    case ('length','length_ps','dump','hmass','tstep')
      fval = float(val)
      call parse_md(mddat,key,fval)
    case ('shake')
      if( val <= 0 )then
         mddat%shake = .false.
      else
         mddat%shake = .true.
         mddat%shk%shake_mode = min(val,2)
      endif
    case ('printstep')
       mddat%printstep = val
    case ('t','temp','temperature' )
      mddat%tsoll = float(val)
      mddat%thermostat =.true.
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
    case ('shake')
       mddat%shake = val
       if(val) mddat%shk%shake_mode=1
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
    if (blk%header .ne. 'dynamics.meta') return
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
      mtd%cvdump_fs = val*1000.0_wp
    case ('ramp')
      mtd%ramp = val
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
      case ('rmsd')
        mtd%mtdtype = cv_rmsd
      case default
        mtd%mtdtype = 0
      end select
    case ('biasfile')
       mtd%mtdtype = cv_rmsd_static
       mtd%biasfile = val
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
