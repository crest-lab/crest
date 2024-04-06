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
  use crest_data
  use crest_calculator,only:calcdata,calculation_settings,jobtype,constraint,scantype
  use dynamics_module
  use gradreader_module,only:gradtype,conv2gradfmt
  use tblite_api,only:xtblvl
  use strucrd,only:get_atlist,coord
  use axis_module

  use parse_block,only:datablock
  use parse_keyvalue,only:keyvalue,valuetypes
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

  subroutine parse_calculation_data(env,calc,dict,included)
    implicit none
    type(systemdata) :: env
    type(calcdata) :: calc
    type(root_object) :: dict
    type(datablock) :: blk
    type(calculation_settings) :: newjob,newjob2
    type(constraint) :: newcstr
    integer :: i,j,k,l
    logical,intent(out) :: included
    type(coord) :: moltmp

    included = .false.
    call calc%reset()
    call env%ref%to(moltmp)
    call axis(moltmp%nat,moltmp%at,moltmp%xyz) 

    do i = 1,dict%nblk
      call blk%deallocate()
      blk = dict%blk_list(i)
      if (blk%header == 'calculation') then
        included = .true.
        call parse_calcdat(env,blk,calc)

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

      else if (index(blk%header,'calculation.constraint').ne.0) then
        call parse_constraintdat(env,moltmp,blk,calc)
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

  subroutine parse_leveldata(blk,job)
!**********************************************************
!* The following routines are used to
!* read information into the "calculation_settings" object
!**********************************************************
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
    !> first, go through settings with fixed type expactations
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
    !> then, all others by key, automatic
    select case (kv%key)
    case default
      continue
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
    case ('uhf','multiplicity')
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
    case ('vdwset')
      job%vdwset = val
    end select
    return
  end subroutine parse_setting_int
  subroutine parse_setting_c(job,key,val)
    implicit none
    type(calculation_settings) :: job
    character(len=*) :: key
    character(len=*) :: val
    logical :: ex
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
      case ('gfn2','gfn2-xtb')
        job%id = jobtype%tblite
        job%tblitelvl = 2
      case ('gfn1','gfn1-xtb')
        job%id = jobtype%tblite
        job%tblitelvl = 1
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
      case ('lj','lennard-jones')
        job%id = jobtype%lj
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

    case ('orca_cmd')
      job%id = jobtype%orca
      job%ORCA%cmd = val
      job%binary = val
    case ('orca_template')
      job%id = jobtype%orca
      call job%ORCA%read(val)

    case ('gbsa','alpb','cpcm')
      job%solvmodel = key
      job%solvent = val

    case ('refine','refinement')
      select case (val)
      case ('sp','singlepoint')
        job%refine_lvl = refine%singlepoint
      case ('add','correction')
        job%refine_lvl = refine%correction
      case ('opt','optimization')
        job%refine_lvl = refine%geoopt
      case default
        job%refine_lvl = refine%non
      end select

    case('restartfile','topo','reftopo')
      inquire(file=val,exist=ex)
      if(ex)then
        job%restart = .true.
        job%restartfile = val
      else
        write(stderr,'(a,a,a)') 'specified restart file ',val,' does not exist'
        error stop
      endif
    case('refgeo','refxyz')
      inquire(file=val,exist=ex)
      if(ex)then
        job%refgeo = val
      else
        write(stderr,'(a,a,a)') 'specified reference geometry file ',val,' does not exist'
        error stop
      endif
    case('parametrisation')
      inquire(file=val,exist=ex)
      if(ex)then
        job%parametrisation  = val
      else
        write(stderr,'(a,a,a)') 'specified parametrisation file ',val,' does not exist'
        error stop
      endif

    case ('print')
      select case (val)
      case ('true','yes')
        job%pr = .true.
      case ('false','no')
        job%pr = .false.
      case ('append','cont','continuous')
        job%pr = .true.
        job%prappend = .true.
      end select
      if (job%pr) job%prch = 999  !> the actual ID will be generated automatically

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
    case ('rdgrad')
      job%rdgrad = val
    case ('refresh')
      job%apiclean = val
    case ('print')
      job%pr = val
      if (val) job%prch = 999  !> the actual ID will be generated automatically
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

  subroutine parse_calcdat(env,blk,calc)
!***********************************************
!* The following routines are used to
!* read information into the "calcdata" object
!***********************************************
    implicit none
    type(systemdata),intent(inout) :: env
    type(datablock),intent(in) :: blk
    type(calcdata),intent(inout) :: calc
    integer :: i
    if (blk%header .ne. 'calculation') return
    do i = 1,blk%nkv
      call parse_calc(env,calc,blk%kv_list(i))
    end do
    return
  end subroutine parse_calcdat
  subroutine parse_calc_auto(env,calc,kv)
    implicit none
    type(systemdata),intent(inout) :: env
    type(calcdata) :: calc
    type(keyvalue) :: kv
    select case (kv%id)
    case (1) !> float
      call parse_calc(env,calc,kv%key,kv%value_f)
    case (2) !> int
      call parse_calc(env,calc,kv%key,kv%value_i)
    case (3) !> bool
      call parse_calc(env,calc,kv%key,kv%value_b)
    case (4) !> string
      call parse_calc(env,calc,kv%key,kv%value_c)
    end select
    !> other, with multiple or raw type
    select case (kv%key)
    case ('optlev','ancopt_level')
      env%optlev = optlevnum(kv%rawvalue)
    end select
  end subroutine parse_calc_auto
  subroutine parse_calc_float(env,calc,key,val)
    implicit none
    type(systemdata),intent(inout) :: env
    type(calcdata) :: calc
    character(len=*) :: key
    real(wp) :: val
    select case (key)
    case ('converge_e','ethr_opt')
      calc%ethr_opt = val  !> optimization ΔE convergenve threshold (Ha)

    case ('converge_g','gthr_opt','rmsforce')
      calc%gthr_opt = val !> optimization RMS convergence threshold (Ha/a0)

    case ('maxerise')
      calc%maxerise = val !> optimization max E rise (Ha)

    case ('displ_opt','maxdispl')
      calc%maxdispl_opt = val !> optimization step size/scaling

    case ('hguess')
      calc%hguess = val  !> guess for the initial hessian

    case default
      return
    end select
    return
  end subroutine parse_calc_float
  subroutine parse_calc_int(env,calc,key,val)
    implicit none
    type(systemdata),intent(inout) :: env
    type(calcdata) :: calc
    character(len=*) :: key
    integer :: val
    select case (key)
    case ('id','type')
      calc%id = val

    case ('maxcycle')
      calc%maxcycle = val  !> optimization max cycles

    case default
      return
    end select
    return
  end subroutine parse_calc_int
  subroutine parse_calc_c(env,calc,key,val)
    implicit none
    type(systemdata),intent(inout) :: env
    type(calcdata) :: calc
    character(len=*) :: key
    character(len=*) :: val
    logical,allocatable :: atlist(:)
    select case (key)
    case ('type')
      select case (val)
      case ('mecp')
        calc%id = -1
      case default
        calc%id = 0
      end select
    case ('elog')
      calc%elog = val
      calc%pr_energies = .true.

    case ('hess_update','hupdate')
      select case (val) !> Hessian updates in geom. Opt.
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

    case ('opt','opt_engine','opt_algo')
      select case (val)
      case ('ancopt','rfo-anc')
        calc%opt_engine = 0
      case ('lbfgs','l-bfgs')
        calc%opt_engine = 1
      case ('rfo','rfo-cart')
        calc%opt_engine = 2
      case ('gd','gradient descent')
        calc%opt_engine = -1
      end select

    case ('freeze')
      call get_atlist(env%ref%nat,atlist,val,env%ref%at)
      calc%nfreeze = count(atlist)
      call move_alloc(atlist,calc%freezelist)

    case default
      return
    end select
    return
  end subroutine parse_calc_c
  subroutine parse_calc_bool(env,calc,key,val)
    implicit none
    type(systemdata),intent(inout) :: env
    type(calcdata) :: calc
    character(len=*) :: key
    logical :: val
    select case (key)
    case ('eprint')
      calc%pr_energies = val

    case ('exact_rf')
      calc%exact_rf = val

    case default
      return
    end select
    return
  end subroutine parse_calc_bool

!========================================================================================!

  subroutine parse_constraintdat(env,mol,blk,calc)
!*************************************************
!* The following routines are used to
!* read information into the "constraint" object
!* and add it to a calculation data object
!*************************************************
    implicit none
    type(systemdata),intent(inout) :: env
    type(coord),intent(inout) :: mol
    type(datablock),intent(in) :: blk
    type(calcdata),intent(inout) :: calc
    logical :: success
    type(constraint) :: constr
    integer :: i
    !type(coord) :: mol 
    logical,allocatable :: atlist(:)
    if (blk%header .ne. 'calculation.constraints' .and.  &
    & blk%header .ne. 'calculation.constraint') return
    success = .false.
    call constr%deallocate()
    do i = 1,blk%nkv
      call parse_constraint_auto(env,constr,blk%kv_list(i),success)

      select case (blk%kv_list(i)%key)
      case ('freeze')
        call get_atlist(env%ref%nat,atlist,blk%kv_list(i)%rawvalue,env%ref%at)
        calc%nfreeze = count(atlist)
        call move_alloc(atlist,calc%freezelist)
      end select
    end do
    if (success) then
      call constr%complete(mol)
      call calc%add(constr)
    end if
    return
  end subroutine parse_constraintdat
  subroutine parse_constraint_auto(env,constr,kv,success)
    implicit none
    type(systemdata) :: env
    type(keyvalue) :: kv
    type(constraint) :: constr
    logical,intent(inout) :: success
    real(wp) :: dum1,dum2,dum3,dum4
    real(wp) :: rabc(3)
    integer :: atm1,atm2,atm3,atm4,n,k,j
    logical,allocatable :: atlist(:)
    !success = .false.
    select case (kv%key)
    case ('type') !> the type of constraint
      select case (kv%value_c)
      case ('bond','bonds'); constr%type = 1
      case ('angle'); constr%type = 2
      case ('dihedral'); constr%type = 3
      case ('wall'); constr%type = 4
      case ('wall_logfermi','ellipsoid'); constr%type = 5
      case ('box'); constr%type = 6
      case ('bondrange'); constr%type = 8
      case ('gapdiff'); constr%type = -1
      case ('gapdiff2','mecp'); constr%type = -2
      end select
      if (constr%type /= 0) success = .true.

    case ('fc','k','forceconstant','params') !> force constants or parameters
      if (allocated(constr%fc)) deallocate (constr%fc)
      select case (kv%id)
      case (valuetypes%int)
        allocate (constr%fc(1),source=0.0_wp)
        constr%fc(1) = kv%value_i
      case (valuetypes%float)
        allocate (constr%fc(1),source=0.0_wp)
        constr%fc(1) = kv%value_f
      case (valuetypes%int_array)
        n = size(kv%value_ia,1)
        allocate (constr%fc(n),source=0.0_wp)
        constr%fc(:) = real(kv%value_ia(:))
      case (valuetypes%float_array)
        n = size(kv%value_fa,1)
        allocate (constr%fc(n),source=0.0_wp)
        constr%fc(:) = real(kv%value_fa(:))
      end select

    case ('atoms')
      if (allocated(constr%atms)) deallocate (constr%atms)
      if (kv%id == valuetypes%int_array) then
        n = size(kv%value_ia,1)
        allocate (constr%atms(n),source=0)
        constr%atms = kv%value_ia
      else
        allocate (atlist(env%ref%nat),source=.false.)
        call get_atlist(env%ref%nat,atlist,kv%rawvalue,env%ref%at)
        n = count(atlist)
        allocate (constr%atms(n),source=0)
        k = 0
        do j = 1,env%ref%nat
          if (atlist(j)) then
            k = k+1
            constr%atms(k) = j
          end if
        end do
        deallocate (atlist) 
      end if
      constr%n = n

    case ('ref','val') !> constrained value
      if (allocated(constr%ref)) deallocate (constr%ref)
      select case (kv%id)
      case (valuetypes%int)
        allocate (constr%ref(1),source=0.0_wp)
        constr%ref(1) = kv%value_i
      case (valuetypes%float)
        allocate (constr%ref(1),source=0.0_wp)
        constr%ref(1) = kv%value_f
      case (valuetypes%int_array)
        n = size(kv%value_ia,1)
        allocate (constr%ref(n),source=0.0_wp)
        constr%ref(:) = real(kv%value_ia(:))
      case (valuetypes%float_array)
        n = size(kv%value_fa,1)
        allocate (constr%ref(n),source=0.0_wp)
        constr%ref(:) = real(kv%value_fa(:))
      end select

    case ('wscal') !> scaling factor if the wall potential is automatically set up
      if(kv%id == valuetypes%int)then
        constr%wscal = max(0.0_wp, real(kv%value_i))
      elseif(kv%id == valuetypes%float)then
        constr%wscal = max(0.0_wp, kv%value_f)
      endif

!>--- the following are for specifiying keywords in a single line
!>--- I don't know it was wise to code them like this because it's hacky,
!>--- but i'll leave them so I don't get confused.
    case ('bond','bonds')
      select case (kv%id)
      case (4) !> string
        select case (kv%value_c)
        case ('all','allauto')
          call constr%dummyconstraint(11)
          success = .true.
        end select
      case (5)  !> regular array
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

    case ('ellipsoid','ellipsoid_logfermi')
      rabc(1:3) = kv%value_fa(1:3)
      if (index(kv%key,'logfermi') .ne. 0) then
        dum1 = 300.0_wp
        dum2 = 6.0_wp
        if (kv%na > 3) dum1 = kv%value_fa(4)
        if (kv%na > 4) dum2 = kv%value_fa(5)
        call constr%ellipsoid(0,atlist,rabc,dum1,dum2,.true.)
      else
        dum1 = 1.0_wp
        dum2 = 30.0_wp
        if (kv%na > 3) dum1 = kv%value_fa(4)
        if (kv%na > 4) dum2 = kv%value_fa(5)
        call constr%ellipsoid(0,atlist,rabc,dum1,dum2,.false.)
      end if
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

    case ('bondrange')
      atm1 = nint(kv%value_fa(1))
      atm2 = nint(kv%value_fa(2))
      dum1 = kv%value_fa(3)*aatoau
      dum1 = max(0.0_wp,dum1) !> can't be negative
      select case (kv%na)
      case (3)
        dum2 = huge(dum2)/3.0_wp !> some huge value
        call constr%bondrangeconstraint(atm1,atm2,dum1,dum2)
      case (4)
        dum2 = kv%value_fa(4)*aatoau
        call constr%bondrangeconstraint(atm1,atm2,dum1,dum2)
      case (5)
        dum3 = kv%value_fa(5)
        call constr%bondrangeconstraint(atm1,atm2,dum1,dum2,beta=dum3)
      case (6)
        dum4 = kv%value_fa(6)
        call constr%bondrangeconstraint(atm1,atm2,dum1,dum2,beta=dum3,T=dum4)
      case default
        error stop '**ERROR** wrong number of arguments in bondrange constraint'
      end select
      success = .true.
!>--------------
!>--------------
!>--------------
    case default
      return
    end select

    return
  end subroutine parse_constraint_auto

!========================================================================================!

  subroutine parse_scandat(blk,calc)
!*******************************************
!* The following routines are used to
!* read information into the "scan" object
!* and add it to a calculation data object
!*******************************************
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
      if (kv%id == valuetypes%float_array) then
        scn%atms(1) = nint(kv%value_fa(1))
        scn%atms(2) = nint(kv%value_fa(2))
        scn%minval = kv%value_fa(3)
        scn%maxval = kv%value_fa(4)
        if (kv%na > 4) then
          scn%steps = nint(kv%value_fa(5))
        end if
        success = .true.
      else if (kv%id == valuetypes%int_array) then
        scn%atms(1) = kv%value_ia(1)
        scn%atms(2) = kv%value_ia(2)
        scn%minval = real(kv%value_ia(3))
        scn%maxval = real(kv%value_ia(4))
        if (kv%na > 4) then
          scn%steps = kv%value_ia(5)
        end if
        success = .true.
      end if

    case ('dihedral')
      scn%type = 3
      scn%n = 2
      write (*,*) kv%value_rawa(:)
      write (*,*) kv%value_ia(:)
      write (*,*) kv%value_fa(:)
      allocate (scn%atms(4))

      if (kv%id == valuetypes%float_array) then
        scn%atms(1) = nint(kv%value_fa(1))
        scn%atms(2) = nint(kv%value_fa(2))
        scn%atms(3) = nint(kv%value_fa(3))
        scn%atms(4) = nint(kv%value_fa(4))
        if (kv%na > 4) then
          scn%steps = nint(kv%value_fa(5))
        end if
        if (kv%na > 6) then
          scn%minval = kv%value_fa(6)
          scn%maxval = kv%value_fa(7)
        end if
        success = .true.
      else if (kv%id == valuetypes%int_array) then
        scn%atms(1) = kv%value_ia(1)
        scn%atms(2) = kv%value_ia(2)
        scn%atms(3) = kv%value_ia(3)
        scn%atms(4) = kv%value_ia(4)
        if (kv%na > 4) then
          scn%steps = kv%value_ia(5)
        end if
        if (kv%na > 6) then
          scn%minval = real(kv%value_ia(6))
          scn%maxval = real(kv%value_ia(7))
        end if
        success = .true.
      end if

    case default
      return
    end select

    return
  end subroutine parse_scan_auto

!========================================================================================!

  subroutine parse_dynamics_data(env,mddat,dict,included)
!*********************************************
!* The following routines are used to
!* read information into the "mddata" object
!*********************************************
    implicit none
    type(systemdata) :: env
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
        call parse_mddat(env,blk,mddat)
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
  subroutine parse_mddat(env,blk,mddat)
    implicit none
    type(systemdata),intent(inout) :: env
    type(datablock),intent(in) :: blk
    type(mddata),intent(inout) :: mddat
    logical,allocatable :: atlist(:)
    integer :: i,j,nat
    if (blk%header .ne. 'dynamics') return
    nat = env%ref%nat
    allocate (atlist(nat),source=.false.)

    do i = 1,blk%nkv
      call parse_md(mddat,blk%kv_list(i))

      select case (blk%kv_list(i)%key)
      case ('includermsd','atlist+')
        call get_atlist(nat,atlist,blk%kv_list(i)%rawvalue,env%ref%at)
        if (.not.allocated(env%includeRMSD)) allocate (env%includeRMSD(nat),source=1)
        do j = 1,nat
          if (atlist(j)) env%includeRMSD(j) = 1
        end do

      case ('excludermsd','atlist-')
        call get_atlist(nat,atlist,blk%kv_list(i)%rawvalue,env%ref%at)
        if (.not.allocated(env%includeRMSD)) allocate (env%includeRMSD(nat),source=1)
        do j = 1,nat
          if (atlist(j)) env%includeRMSD(j) = 0
        end do

      end select

    end do
    deallocate (atlist)
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
    case default
      select case (kv%key)
      case ('active','active_levels')
        mddat%active_potentials = kv%value_ia
      end select
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
    case ('t','temp','temperature')
      mddat%tsoll = val
      mddat%thermostat = .true.
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
      if (val <= 0) then
        mddat%shake = .false.
      else
        mddat%shake = .true.
        mddat%shk%shake_mode = min(val,2)
      end if
    case ('printstep')
      mddat%printstep = val
    case ('t','temp','temperature')
      mddat%tsoll = float(val)
      mddat%thermostat = .true.
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
      if (val) mddat%shk%shake_mode = 1
    case default
      return
    end select
    return
  end subroutine parse_md_bool

!========================================================================================!

  subroutine parse_metadyn(blk,mddat)
!**************************************************
!* The following routines are used to
!* read information into the "metadynamics" object
!* and add it to a mol.dynamics data object
!***************************************************
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
!========================================================================================!
end module parse_calcdata
