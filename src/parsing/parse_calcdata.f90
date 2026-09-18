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
  use bh_module
  use gradreader_module,only:gradtype,conv2gradfmt
  use tblite_api,only:xtblvl,have_gxtb
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
  end interface parse_setting

!>-- routines for parsing a calcdata object
  interface parse_calc
    module procedure :: parse_calc_auto
  end interface parse_calc

!>-- routines for parsing a mddata object
  interface parse_md
    module procedure :: parse_md_auto
  end interface parse_md

!>-- routines for parsing a mtdpot object
  interface parse_mtd
    module procedure :: parse_metadyn_auto
  end interface parse_mtd

  public :: parse_calculation_data
  public :: parse_dynamics_data
  public :: parse_basinhopping_data

  character(len=*),parameter,private :: fmturk = '("unrecognized KEYWORD in ",a," : ",a)'
  character(len=*),parameter,private :: fmtura = '("unrecognized ARGUMENT : ",a)'

  external creststop

!========================================================================================!
!========================================================================================!
contains !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine parse_calculation_data(env,calc,dict,included,istat)
    implicit none
    type(systemdata) :: env
    type(calcdata) :: calc
    type(root_object) :: dict
    type(datablock) :: blk
    type(calculation_settings) :: newjob,newjob2
    type(constraint) :: newcstr
    integer :: i,j,k,l
    logical,intent(out)   :: included
    integer,intent(inout) :: istat
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
        call parse_calcdat(env,blk,calc,istat)

      else if (blk%header == 'calculation.level') then
        call parse_leveldata(env,blk,newjob,istat)
        call newjob%autocomplete(calc%ncalculations+1)
        call calc%add(newjob)
        included = .true.

      else if (blk%header == 'calculation.mecp') then
        !>-- setup
        if (allocated(calc%calcs)) deallocate (calc%calcs)
        calc%ncalculations = 0
        calc%id = -1
        call parse_leveldata(env,blk,newjob,istat)
        !>-- S0 setup
        call parse_leveldata(env,blk,newjob,istat)
        newjob%uhf = 0
        newjob%calcspace = 's0'
        call calc%add(newjob)
        !>-- S1 setup
        newjob%uhf = 2
        newjob%calcspace = 's1'
        call calc%add(newjob)
        included = .true.

      else if (index(blk%header,'calculation.constraint') .ne. 0) then
        call parse_constraintdat(env,moltmp,blk,calc,istat)
        included = .true.

      else if (blk%header == 'calculation.scans') then
        call parse_scandat(blk,calc,istat)
        included = .true.

      end if
    end do
    if (included) then
      call calc%init()
    end if
    return
  end subroutine parse_calculation_data

!========================================================================================!

  subroutine parse_leveldata(env,blk,job,istat)
!**********************************************************
!* The following routines are used to
!* read information into the "calculation_settings" object
!**********************************************************
    implicit none
    type(systemdata),intent(inout) :: env
    type(datablock),intent(in) :: blk
    type(calculation_settings),intent(out) :: job
    integer,intent(inout) :: istat
    logical :: rd
    integer :: i
    call job%deallocate()
    if ((blk%header .ne. 'calculation.level').and. &
    & (blk%header .ne. 'calculation.mecp')) then
      return
    end if
    do i = 1,blk%nkv
      call parse_setting_auto(env,job,blk%kv_list(i),rd)
      if (.not.rd) then
        istat = istat+1
        write (stdout,fmturk) '[['//blk%header//']]-block',blk%kv_list(i)%key
      end if
    end do

!>-- method-specific defaults, applied only after the entire block was read
!>   so that they do not depend on the order of the key/value pairs
    call parse_level_defaults(job)
    return
  end subroutine parse_leveldata

!========================================================================================!

  subroutine parse_level_defaults(job)
!************************************************************
!* Apply defaults that depend on the selected method but must
!* not overwrite anything the user specified explicitly.
!* Currently only the electronic temperature: g-xTB works with
!* integer occupations, i.e. Tel = 0 K, whereas the GFNn-xTB
!* Hamiltonians default to 300 K.
!************************************************************
    implicit none
    type(calculation_settings),intent(inout) :: job

    if (job%tblitelvl == xtblvl%gxtb.and..not.job%etemp_user_set) then
      job%etemp = 0.0_wp
    end if

  end subroutine parse_level_defaults

!========================================================================================!

  subroutine parse_setting_auto(env,job,kv,rd)
    implicit none
    type(systemdata),intent(inout) :: env
    type(calculation_settings) :: job
    type(keyvalue) :: kv
    logical,intent(out) :: rd
    logical :: ex
    integer :: n
    rd = .true.
    select case (kv%key)

!>--- floats
    case ('etemp')
      job%etemp = kv%value_f
      job%etemp_user_set = .true.
    case ('accuracy')
      job%accuracy = kv%value_f
    case ('weight')
      job%weight = kv%value_f
    case ('pressure')
      job%extpressure = kv%value_f
    case ('proberad','pvol_proberad')
      job%proberad = kv%value_f
    case ('radscal','pvol_radscal')
      job%pvradscal = kv%value_f
    case ('efield')
      n = size(kv%value_fa,1)
      if (n .ne. 3) then
        !>--- keyword was recognized, but invalid argument supplied
        write (stdout,fmtura) trim(kv%rawvalue)
        call creststop(status_config)
      end if
      allocate (job%efield(3),source=0.0_wp)
      job%efield(:) = kv%value_fa(:)

!>--- integers
    case ('uhf')
      job%uhf = kv%value_i
      call job%sync_multiplicity()
    case ('multiplicity','mult')
      call job%set_multiplicity(kv%value_i)
    case ('chrg','charge')
      job%chrg = kv%value_i
    case ('id')
      job%id = kv%value_i
    case ('maxscc')
      job%maxscc = kv%value_i
    case ('lebedev','pvol_ngrid')
      job%ngrid = kv%value_i
    case ('vdwset','pvol_vdwset')
      job%vdwset = kv%value_i
    case ('mlip_port','mlip_base_port')
      job%MPAR%BASE_PORT = kv%value_i
    case ('mlip_timeout','mlip_timeout_sec')
      job%MPAR%TIMEOUT_SEC = kv%value_i
    case ('config')
      call job%addconfig(kv%value_ia)

!>--- strings
    case ('method')
      select case (kv%value_c)
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
        job%tblitelvl = xtblvl%gfn2
      case ('gfn1','gfn1-xtb')
        job%id = jobtype%tblite
        job%tblitelvl = xtblvl%gfn1
      case ('ceh')
        job%id = jobtype%tblite
        job%tblitelvl = xtblvl%ceh
        job%rdgrad = .false.
        job%rdqat = .true.
        job%rddip = .true.
      case ('gfn0','gfn0-xtb')
        job%id = jobtype%gfn0
      case ('gfn0*','gfn0*-xtb')
        job%id = jobtype%gfn0occ
      case ('gfnff','gff','gfn-ff')
        job%id = jobtype%gfnff
      case ('pvol','libpvol','pv')
        job%id = jobtype%libpvol
      case ('gxtb','g-xtb','gxtb-xtb')
        if (have_gxtb) then
          job%id = jobtype%tblite
          job%tblitelvl = xtblvl%gxtb
        else
          job%id = jobtype%xtbsys
          job%other = '--gxtb'
        end if
      case ('gxtb_dev')
        if (have_gxtb) then
          call gxtb_dev_warning()
        else
          job%id = jobtype%xtbsys
          job%other = '--gxtb'
        end if
      case ('none')
        job%id = jobtype%unknown
      case ('lj','lennard-jones')
        job%id = jobtype%lj
      case ('modh')
        job%id = jobtype%approxg
      case ('rmsdbias','penalty')
        job%id = jobtype%penalty
        nullify (job%penalty%biaslist)
      case ('mlip','fmlip_relay')
        job%id = jobtype%mlip
      case ('solvation','solvator')
        job%id = jobtype%solvation
        if (.not.allocated(job%solv)) allocate (job%solv)
      case ('eeq','eeqbc','electrostatic')
        job%id = jobtype%electrostatic
        if (.not.allocated(job%eeq)) allocate (job%eeq)
        if (kv%value_c /= 'electrostatic') job%eeq%charge_model = kv%value_c
      case default
        job%id = jobtype%unknown
        !>--- keyword was recognized, but invalid argument supplied
        write (stdout,fmtura) kv%value_c
        call creststop(status_config)
      end select

    case ('bin','binary','script')
      job%binary = kv%value_c

    case ('flags')
      job%other = kv%value_c

      !> don't.
      !case ('sys','syscall','systemcall')
      !  job%systemcall = val

    case ('calcspace','dir')
      job%calcspace = kv%value_c

    case ('threads','ncores')
      job%threads = kv%value_i

    case ('gradfile')
      job%gradfile = kv%value_c

    case ('gradtype')
      select case (kv%value_c)
      case ('engrad','xtb','orca')
        job%gradtype = gradtype%engrad
      case ('turbomole','tm')
        job%gradtype = gradtype%turbomole
      case ('generic')
        job%gradtype = gradtype%unknown
      case default
        job%gradtype = gradtype%unknown
        !>--- keyword was recognized, but invalid argument supplied
        write (stdout,fmtura) kv%value_c
        call creststop(status_config)

      end select

    case ('gradkey')
      job%gradkey = kv%value_c

    case ('gradmt')
      job%gradfmt = conv2gradfmt(kv%value_c)

    case ('numgrad')
      job%numgrad = kv%value_b
    case ('gradstep')
      job%gradstep = kv%value_f

    case ('efile')
      job%efile = kv%value_c

    case ('tblite_level','tblite_hamiltonian')
      select case (kv%value_c)
      case ('gfn2','gfn2-xtb')
        job%tblitelvl = xtblvl%gfn2
      case ('gfn1','gfn1-xtb')
        job%tblitelvl = xtblvl%gfn1
      case ('ipea1')
        job%tblitelvl = xtblvl%ipea1
      case ('ceh')
        job%tblitelvl = xtblvl%ceh
        job%rdgrad = .false.
      case ('eeq','d4eeq')
        job%tblitelvl = xtblvl%eeq
        job%rdgrad = .false.
      case ('gxtb','g-xtb')
        job%tblitelvl = xtblvl%gxtb
      case default
        job%tblitelvl = xtblvl%unknown
        !>--- keyword was recognized, but invalid argument supplied
        write (stdout,fmtura) kv%value_c
        call creststop(status_config)
      end select

    case ('tblite_param')
      job%tbliteparam = kv%value_c
      job%tblitelvl = xtblvl%param

    case ('mlip_backend')
      job%MPAR%backend = kv%value_c

    case ('mlip_modelpath')
      job%MPAR%modelpath = kv%value_c

    case ('mlip_modelsize','mlip_modelname')
      job%MPAR%modelsize = kv%value_c

    case ('mlip_device')
      job%MPAR%device = kv%value_c

    case ('mlip_uma_model','mlip_umamodel')
      job%MPAR%umamodel = kv%value_c

    case ('mlip_uma_task','mlip_umatask')
      job%MPAR%umatask = kv%value_c

    case ('orca_cmd')
      job%id = jobtype%orca
      job%ORCA%cmd = kv%value_c
      job%binary = kv%value_c
    case ('orca_template')
      job%id = jobtype%orca
      if (job%ORCA%srckind == 2) then
        write (stdout,'(a)') '**ERROR** orca_template and orca_input are mutually exclusive!'
        call creststop(status_config)
      end if
      call job%ORCA%read(kv%value_c)
    case ('orca_input')
      job%id = jobtype%orca
      if (job%ORCA%srckind == 1) then
        write (stdout,'(a)') '**ERROR** orca_template and orca_input are mutually exclusive!'
        call creststop(status_config)
      end if
      call job%ORCA%build(kv%value_c)
    case ('orca_memory','orca_maxcore')
      job%id = jobtype%orca
      job%ORCA%maxcore = kv%value_i

    case ('gbsa','alpb','cpcm','cosmo','pcm')
      job%solvmodel = kv%key
      job%solvent = kv%value_c

!>--- composite solvation calculator (method = "solvation") settings
    case ('solvent')
      if (.not.allocated(job%solv)) allocate (job%solv)
      job%solv%solvent = kv%value_c
    case ('solv_model')
      if (.not.allocated(job%solv)) allocate (job%solv)
      job%solv%smodel = kv%value_c
    case ('solv_charges')
      if (.not.allocated(job%solv)) allocate (job%solv)
      job%solv%charge_model = kv%value_c
    case ('solv_hbond')
      if (.not.allocated(job%solv)) allocate (job%solv)
      job%solv%do_hbond = kv%value_b

!>--- charge-equilibration electrostatics calculator settings
    case ('eeq_model','charge_model')
      if (.not.allocated(job%eeq)) allocate (job%eeq)
      job%eeq%charge_model = kv%value_c

    case ('refine','refinement')
      select case (kv%value_c)
      case ('sp','singlepoint')
        job%refine_lvl = refine%singlepoint
      case ('add','correction')
        job%refine_lvl = refine%correction
      case ('opt','optimization')
        job%refine_lvl = refine%geoopt
      case ('freq','frequencies','deltag')
        job%refine_lvl = refine%deltaG
      case default
        job%refine_lvl = refine%non
        !>--- keyword was recognized, but invalid argument supplied
        write (stdout,fmtura) kv%value_c
        call creststop(status_config)

      end select

    case ('restartfile','topo','reftopo')
      inquire (file=kv%value_c,exist=ex)
      if (ex) then
        job%restart = .true.
        job%restartfile = kv%value_c
      else
        write (stderr,'(a,a,a)') 'specified restart file ',kv%value_c,' does not exist'
        call creststop(status_input)
      end if

    case ('refgeo','refxyz')
      inquire (file=kv%value_c,exist=ex)
      if (ex) then
        job%refgeo = kv%value_c
      else
        write (stderr,'(a,a,a)') 'specified reference geometry file ',kv%value_c,' does not exist'
        call creststop(status_input)
      end if

    case ('biasfile')
      inquire (file=kv%value_c,exist=ex)
      if (ex) then
        job%penalty%biasfile = kv%value_c
      else
        write (stderr,'(a,a,a)') 'specified bias file ',kv%value_c,' does not exist'
        call creststop(status_input)
      end if

    case ('penalty_kpush')
      job%penalty%kpush = kv%value_f

    case ('penalty_alpha')
      job%penalty%alpha = kv%value_f

    case ('parametrisation')
      inquire (file=kv%value_c,exist=ex)
      if (ex) then
        job%parametrisation = kv%value_c
      else
        write (stderr,'(a,a,a)') 'specified parametrisation file ',kv%value_c,' does not exist'
        call creststop(status_input)
      end if

    case ('refchrg','refcharges')
      inquire (file=kv%value_c,exist=ex)
      if (ex) then
        job%refcharges = kv%value_c
      else
        write (stderr,'(a,a,a)') 'specified reference charge file ',kv%value_c,' does not exist'
        call creststop(status_config)
      end if

    case ('fragment','fragments')
      !>--- user-defined fragmentation, e.g. fragments = ["1-432","433-458"]
      !>    each entry is an atom-range string defining one fragment.
      !>    A single string (one fragment) is also accepted.
      if (allocated(kv%value_ca)) then
        job%gff_fragments = kv%value_ca
      else if (allocated(kv%value_c)) then
        allocate (character(len=len(kv%value_c)) :: job%gff_fragments(1))
        job%gff_fragments(1) = kv%value_c
      else
        write (stderr,'(a)') 'could not read fragment definition'
        call creststop(status_config)
      end if

    case ('print')
      select case (kv%id)
      case (2)
        select case (kv%value_c)
        case ('true','yes')
          job%pr = .true.
        case ('false','no')
          job%pr = .false.
        case ('append','cont','continuous')
          job%pr = .true.
          job%prappend = .true.
        end select
      case (3)
        job%pr = kv%value_b
      end select
      if (job%pr) job%prch = 999  !> the actual ID will be generated automatically

    case ('getsasa')
      call get_atlist(env%ref%nat,job%getsasa,kv%value_c,env%ref%at)

    case ('pvol_model','pgrad','pvgrad')
      select case (kv%id)
      case (valuetypes%int)
        job%pvmodel = kv%value_i
      case (valuetypes%string)
        select case (kv%value_c)
        case ('xhcff')
          job%pvmodel = 0
        case ('analytic')
          job%pvmodel = 1
        case default
          write (stdout,fmtura) kv%value_c
          call creststop(status_config)
        end select
      case default
        write (stdout,fmtura) trim(kv%rawvalue)
        call creststop(status_config)
      end select

!>--- booleans
    case ('rdwbo')
      job%rdwbo = kv%value_b
    case ('rddip','rddipole')
      job%rddip = kv%value_b
    case ('rdqat','rdchrg')
      job%rdqat = kv%value_b
    case ('dumpq','dumpchrg')
      job%rdqat = kv%value_b
      job%dumpq = kv%value_b
    case ('dipgrad')
      job%rddipgrad = kv%value_b
    case ('rdgrad')
      job%rdgrad = kv%value_b
    case ('refresh')
      job%apiclean = kv%value_b
    case ('lmo','lmocent')
      job%getlmocent = kv%value_b
    case ('ceh_guess')
      job%ceh_guess = kv%value_b
    case ('spin_polarized')
      job%spin_polarized = kv%value_b

    case default
      !>--- keyword not correctly read/found
      rd = .false.
      continue
    end select
  end subroutine parse_setting_auto

!========================================================================================!

  subroutine parse_calcdat(env,blk,calc,istat)
!***********************************************
!* The following routines are used to
!* read information into the "calcdata" object
!***********************************************
    implicit none
    type(systemdata),intent(inout) :: env
    type(datablock),intent(in) :: blk
    type(calcdata),intent(inout) :: calc
    integer,intent(inout) :: istat
    integer :: i
    logical :: rd
    if (blk%header .ne. 'calculation') return
    do i = 1,blk%nkv
      call parse_calc_auto(env,calc,blk%kv_list(i),rd)
      if (.not.rd) then
        istat = istat+1
        write (stdout,fmturk) '[calculation]-block',blk%kv_list(i)%key
      end if
    end do
    return
  end subroutine parse_calcdat
  subroutine parse_calc_auto(env,calc,kv,rd)
    implicit none
    type(systemdata),intent(inout) :: env
    type(calcdata) :: calc
    type(keyvalue) :: kv
    logical,intent(out) :: rd
    logical,allocatable :: atlist(:)
    rd = .true.
    select case (kv%key)
    case ('optlev','ancopt_level')
      env%optlev = optlevnum(kv%rawvalue)

!>--- floats
    case ('converge_e','ethr_opt')
      calc%ethr_opt = kv%value_f  !> optimization ΔE convergenve threshold (Ha)

    case ('converge_g','gthr_opt','rmsforce')
      calc%gthr_opt = kv%value_f !> optimization RMS convergence threshold (Ha/a0)

    case ('maxerise')
      calc%maxerise = kv%value_f !> optimization max E rise (Ha)

    case ('displ_opt','maxdispl')
      calc%maxdispl_opt = kv%value_f !> optimization step size/scaling

    case ('hguess')
      calc%hguess = kv%value_f  !> guess for the initial hessian

    case ('opt_lval')
      calc%L = kv%value_f !> Parameters for smooth function for stepsize control within optimizer

    case ('opt_k')
      calc%k = kv%value_f

    case ('opt_shift')
      calc%shift = kv%value_f

    case ('scaling')
      calc%scaling = kv%value_f

    case ('doh_stepsize')
      calc%doh_stepsize = kv%value_f

    case ('chess_id_guess')
      calc%chess_id_guess = kv%value_f

!>--- integers
    case ('maxcycle')
      calc%maxcycle = kv%value_i  !> optimization max cycles

    case ('chess_steps')
      calc%hu_steps = kv%value_i

!>--- strings
    case ('id','type')
      !> (OLD setting) calculation type
      select case (kv%id)
      case (2)
        calc%id = kv%value_i
      case (4)
        select case (kv%value_c)
        case ('mecp')
          calc%id = -1
        case default
          calc%id = 0
        end select
      end select

    case ('elog')
      calc%elog = kv%value_c
      calc%pr_energies = .true.

    case ('hess_update','hupdate')
      select case (kv%value_c) !> Hessian updates in geom. Opt.
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
        !>--- keyword was recognized, but invalid argument supplied
        write (stdout,fmtura) kv%value_c
        call creststop(status_config)
      end select

    case ('opt','opt_engine','opt_algo')
      select case (kv%value_c)
      case ('ancopt','rfo-anc')
        calc%opt_engine = 0
      case ('lbfgs','l-bfgs')
        calc%opt_engine = 1
      case ('rfo','rfo-cart')
        calc%opt_engine = 2
      case ('newton','nr')
        calc%opt_engine = 3
      case ('gd','gradient descent')
        calc%opt_engine = -1
      case default
        !>--- keyword was recognized, but invalid argument supplied
        write (stdout,fmtura) kv%value_c
        call creststop(status_config)
      end select

    case ('hr_init','hr_initialization') !> here we set how the matrix for hessian reconstruction is initialized
      select case (kv%value_c)
      case ('identity')
        calc%initialize_hr_type = 0
      case ('gfnff','gfn-ff')
        calc%initialize_hr_type = 1
      case ('gfn0')
        calc%initialize_hr_type = 2
      case ('gfn1')
        calc%initialize_hr_type = 3
      case ('gfn2')
        calc%initialize_hr_type = 4
      case ('modhess')
        calc%initialize_hr_type = 5
      case default
        !>--- keyword was recognized, but invalid argument supplied
        write (stdout,fmtura) kv%value_c
        call creststop(status_config)
      end select

    case ('gs_hess_type')
      select case (kv%value_c)
      case ('gfnff','gfn-ff')
        calc%gs_hess_type = 1
      case ('gfn0')
        calc%gs_hess_type = 2
      case ('gfn1')
        calc%gs_hess_type = 3
      case ('gfn2')
        calc%gs_hess_type = 4
      case ('modhess')
        calc%gs_hess_type = 5
      case default
        !>--- keyword was recognized, but
        write (stdout,fmtura) kv%value_c
        call creststop(status_config)
      end select

    case ('modhess_type','mh_type') !> here we set how the matrix for hessian reconstruction is initialized
      select case (kv%value_c)      !>maybe need to add another keywort for crosstesting hr and geopt -> No,
      case ('lindh95')
        calc%mh_type = 0
      case ('lindh')
        calc%mh_type = 1
      case ('lindh07')
        calc%mh_type = 2
      case ('swart')
        calc%mh_type = 3
      case default
        !>--- keyword was recognized, but invalid argument supplied
        write (stdout,fmtura) kv%value_c
        call creststop(status_config)
      end select

    case ('hess_init','hess_initialization') !> here we set how the hessian for optimization
      select case (kv%value_c)
      case ('identity')
        calc%hess_init = 0
      case ('gfnff','gfn-ff')
        calc%hess_init = 1
      case ('gfn0')
        calc%hess_init = 2
      case ('gfn1')
        calc%hess_init = 3
      case ('gfn2')
        calc%hess_init = 4
      case ('modhess')
        calc%hess_init = 5
      case default
        !>--- keyword was recognized, but invalid argument supplied
        write (stdout,fmtura) kv%value_c
        call creststop(status_config)
      end select

    case ('hr_hess_update','hr_hu_update')
      select case (kv%value_c) !> Hessian updates in hessian reconstruction
      case ('bfgs')
        calc%hr_hu_type = 0
      case ('powell')
        calc%hr_hu_type = 1
      case ('sr1')
        calc%hr_hu_type = 2
      case ('bofill')
        calc%hr_hu_type = 3
      case ('schlegel')
        calc%hr_hu_type = 4
      case default
        !>--- keyword was recognized, but invalid argument supplied
        write (stdout,fmtura) kv%value_c
        call creststop(status_config)
      end select

    case ('freeze')
      call get_atlist(env%ref%nat,atlist,kv%value_c,env%ref%at)
      calc%nfreeze = count(atlist)
      call move_alloc(atlist,calc%freezelist)

!>--- booleans
    case ('eprint')
      calc%pr_energies = kv%value_b

    case ('exact_rf')
      calc%exact_rf = kv%value_b

    case ('chess')
      calc%do_HR = kv%value_b

    case ('full_chess') !> Do Hessian Reconstruct with all optimization steps
      calc%full_HR = kv%value_b

    case ('deform_opt_hess')
      calc%deform_opt_hess = kv%value_b

    case ("g_sampling") !> Do sampling on free energy surface as approximated by lindh95 hessian
      calc%g_sampling = kv%value_b

    case default
      rd = .false.
    end select
  end subroutine parse_calc_auto

!========================================================================================!

  subroutine parse_constraintdat(env,mol,blk,calc,istat)
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
    integer,intent(inout) :: istat
    logical :: success
    type(constraint) :: constr
    integer :: i
    logical :: rd
    if (blk%header .ne. 'calculation.constraints'.and.  &
    & blk%header .ne. 'calculation.constraint') return
    success = .false.
    call constr%deallocate()
    do i = 1,blk%nkv
      call parse_constraint_auto(env,calc,constr,blk%kv_list(i),success,rd)
      if (.not.rd) then
        istat = istat+1
        write (stdout,fmturk) '[['//blk%header//']]-block',blk%kv_list(i)%key
      end if
    end do
    if (success) then
      call constr%complete(mol)
      call calc%add(constr)
    end if
    return
  end subroutine parse_constraintdat
  subroutine parse_constraint_auto(env,calc,constr,kv,success,rd)
    implicit none
    type(systemdata) :: env
    type(keyvalue) :: kv
    type(constraint) :: constr
    logical,intent(inout) :: success
    type(calcdata),intent(inout) :: calc
    real(wp) :: dum1,dum2,dum3,dum4
    real(wp) :: rabc(3)
    integer :: atm1,atm2,atm3,atm4,n,k,j
    logical,allocatable :: atlist(:)
    logical,intent(out) :: rd
    rd = .true.
    select case (kv%key)

    case ('freeze')
      call get_atlist(env%ref%nat,atlist,kv%rawvalue,env%ref%at)
      calc%nfreeze = count(atlist)
      call move_alloc(atlist,calc%freezelist)

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
      case default
        !>--- keyword was recognized, but invalid argument supplied
        write (stdout,fmtura) kv%value_c
        call creststop(status_config)
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
      if (kv%id == valuetypes%int) then
        constr%wscal = max(0.0_wp,real(kv%value_i))
      elseif (kv%id == valuetypes%float) then
        constr%wscal = max(0.0_wp,kv%value_f)
      end if

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
        write (stdout,'(a)') '**ERROR** wrong number of arguments in bondrange constraint'
        call creststop(status_config)
      end select
      success = .true.
!>--------------
!>--------------
!>--------------
    case default
      rd = .false.
      return
    end select

    return
  end subroutine parse_constraint_auto

!========================================================================================!

  subroutine parse_scandat(blk,calc,istat)
!*******************************************
!* The following routines are used to
!* read information into the "scan" object
!* and add it to a calculation data object
!*******************************************
    implicit none
    type(datablock),intent(in) :: blk
    type(calcdata),intent(inout) :: calc
    integer,intent(inout) :: istat
    logical :: success
    type(scantype) :: scn
    integer :: i
    logical :: rd
    if (blk%header .ne. 'calculation.scans') return
    do i = 1,blk%nkv
      call parse_scan_auto(scn,blk%kv_list(i),success,rd)
      if (.not.rd) then
        istat = istat+1
        write (stdout,fmturk) '[['//blk%header//']]-block',blk%kv_list(i)%key
      end if
      if (success) then
        call calc%add(scn)
      end if
    end do
    return
  end subroutine parse_scandat
  subroutine parse_scan_auto(scn,kv,success,rd)
    implicit none
    type(keyvalue) :: kv
    type(scantype) :: scn
    logical,intent(out) :: success,rd
    real(wp) :: dum1,dum2,dum3
    integer :: atm1,atm2,atm3,atm4
    integer :: nsteps
    success = .false.
    rd = .true.
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
      rd = .false.
      return
    end select

    return
  end subroutine parse_scan_auto

!========================================================================================!

  subroutine parse_dynamics_data(env,mddat,dict,included,istat)
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
    integer,intent(inout) :: istat

    included = .false.

    do i = 1,dict%nblk
      call blk%deallocate()
      blk = dict%blk_list(i)
      if (blk%header == 'dynamics') then
        included = .true.
        call parse_mddat(env,blk,mddat,istat)
      else if (blk%header == 'dynamics.meta') then
        call parse_metadyn(env,blk,mddat,istat)
        included = .true.
      end if
    end do
    if (included) then
      mddat%requested = .true.
    end if
    return
  end subroutine parse_dynamics_data
  subroutine parse_mddat(env,blk,mddat,istat)
    implicit none
    type(systemdata),intent(inout) :: env
    type(datablock),intent(in) :: blk
    type(mddata),intent(inout) :: mddat
    integer,intent(inout) :: istat
    logical,allocatable :: atlist(:)
    integer :: i,j,nat
    logical :: rd
    if (blk%header .ne. 'dynamics') return
    nat = env%ref%nat
    allocate (atlist(nat),source=.false.)

    do i = 1,blk%nkv
      call parse_md(env,mddat,blk%kv_list(i),rd)
      if (.not.rd) then
        istat = istat+1
        write (stdout,fmturk) '['//blk%header//']-block',blk%kv_list(i)%key
      end if
    end do
    deallocate (atlist)
    return
  end subroutine parse_mddat
  subroutine parse_md_auto(env,mddat,kv,rd)
    implicit none
    type(systemdata),intent(inout) :: env
    type(mddata) :: mddat
    type(keyvalue) :: kv
    logical,intent(out) :: rd
    logical,allocatable :: atlist(:)
    integer :: nat,j
    logical :: ex
    rd = .true.

    select case (kv%key)
    case ('active','active_levels')
      mddat%active_potentials = kv%value_ia

    case ('includermsd','atlist+')
      nat = env%ref%nat
      call get_atlist(nat,atlist,kv%rawvalue,env%ref%at)
      if (.not.allocated(env%includeRMSD)) allocate (env%includeRMSD(nat),source=0)
      do j = 1,nat
        if (atlist(j)) env%includeRMSD(j) = 1
      end do

    case ('excludermsd','atlist-')
      nat = env%ref%nat
      call get_atlist(nat,atlist,kv%rawvalue,env%ref%at)
      if (.not.allocated(env%includeRMSD)) allocate (env%includeRMSD(nat),source=1)
      do j = 1,nat
        if (atlist(j)) env%includeRMSD(j) = 0
      end do

    case ('length','length_ps')
      mddat%length_ps = kv%value_f
      env%user_mdtime = .true.
    case ('dump')
      mddat%dumpstep = kv%value_f
    case ('hmass')
      mddat%md_hmass = kv%value_f
    case ('tstep')
      mddat%tstep = kv%value_f
    case ('t','temp','temperature')
      mddat%tsoll = kv%value_f
      mddat%thermostat = .true.

    case ('thermostat')
      select case (kv%value_c)
      case ('off','nve')
        mddat%thermotype = 'none'
      case ('berendsen','langevin','bbk')
        mddat%thermotype = trim(kv%value_c)
      case ('bussi','bussi-donaido-parinello','bussi-parinello','csvr')
        mddat%thermotype = 'bussi'
      case default
        write (stdout,fmtura) kv%value_c
        call creststop(status_config)
      end select
      mddat%thermostat = .true.

    case ('shake')
      select case (kv%id)
      case (valuetypes%int)
        if (kv%value_i <= 0) then
          mddat%shake = .false.
        else
          mddat%shake = .true.
          mddat%shk%shake_mode = min(kv%value_i,2)
        end if
      case (valuetypes%bool)
        mddat%shake = kv%value_b
        if (kv%value_b) mddat%shk%shake_mode = 1
      end select
    case ('printstep')
      mddat%printstep = kv%value_i
    case ('blocklength','blockl')
      mddat%blockl = kv%value_i

    case ('restart')
      inquire (file=trim(kv%value_c),exist=ex)
      if (ex) then
        mddat%restart = .true.
        mddat%restartfile = trim(kv%value_c)
      end if

    case default
      rd = .false.
      return
    end select

  end subroutine parse_md_auto

!========================================================================================!

  subroutine parse_metadyn(env,blk,mddat,istat)
!**************************************************
!* The following routines are used to
!* read information into the "metadynamics" object
!* and add it to a mol.dynamics data object
!***************************************************
    implicit none
    type(systemdata),intent(inout) :: env
    type(datablock),intent(in) :: blk
    type(mddata),intent(inout) :: mddat
    integer,intent(inout) :: istat
    logical :: success
    type(mtdpot) :: mtd
    integer :: i,k
    logical :: rd
    call mtd%deallocate()
    success = .false.
    if (blk%header .ne. 'dynamics.meta') return
    do i = 1,blk%nkv
      call parse_metadyn_auto(env,mtd,blk%kv_list(i),success,rd)
      if (.not.rd) then
        istat = istat+1
        write (stdout,fmturk) '[['//blk%header//']]-block',blk%kv_list(i)%key
      end if
    end do
    if (success) call mddat%add(mtd)
    return
  end subroutine parse_metadyn
  subroutine parse_metadyn_auto(env,mtd,kv,success,rd)
    implicit none
    type(systemdata),intent(inout) :: env
    type(keyvalue) :: kv
    type(mtdpot) :: mtd
    logical,intent(inout) :: success
    logical,intent(out) :: rd
    integer :: j,nat
    logical,allocatable :: atlist(:)
    rd = .true.

    select case (kv%key)
!>--- floats
    case ('alpha')
      mtd%alpha = kv%value_f
    case ('kpush')
      mtd%kpush = kv%value_f
    case ('dump','dump_fs')
      mtd%cvdump_fs = kv%value_f
    case ('dump_ps')
      mtd%cvdump_fs = kv%value_f*1000.0_wp
    case ('ramp')
      mtd%ramp = kv%value_f

!>--- strings
    case ('type')
      success = .true.
      select case (kv%id)
      case (valuetypes%int)
        mtd%mtdtype = kv%value_i
      case (valuetypes%string)
        select case (kv%value_c)
        case ('rmsd')
          mtd%mtdtype = cv_rmsd
        case default
          mtd%mtdtype = 0
        end select
      end select
    case ('biasfile')
      mtd%mtdtype = cv_rmsd_static
      mtd%biasfile = kv%value_c

    case ('includermsd','atlist+')
      nat = env%ref%nat
      call get_atlist(nat,atlist,kv%rawvalue,env%ref%at)
      if (.not.allocated(mtd%atinclude)) allocate (mtd%atinclude(nat),source=.false.)
      do j = 1,nat
        if (atlist(j)) mtd%atinclude(j) = .true.
      end do

    case ('excludermsd','atlist-')
      nat = env%ref%nat
      call get_atlist(nat,atlist,kv%rawvalue,env%ref%at)
      if (.not.allocated(mtd%atinclude)) allocate (mtd%atinclude(nat),source=.true.)
      do j = 1,nat
        if (atlist(j)) mtd%atinclude(j) = .false.
      end do

    case default
      rd = .false.
      return
    end select

  end subroutine parse_metadyn_auto

!========================================================================================!

  subroutine parse_basinhopping_data(env,bh,dict,included,istat)
!**********************************************
!* The following routines are used to
!* read information into the "bh_class" object
!**********************************************
    implicit none
    type(systemdata) :: env
    type(bh_class) :: bh
    type(root_object) :: dict
    type(datablock) :: blk
    type(calculation_settings) :: newjob
    type(constraint) :: newcstr
    integer :: i,j,k,l
    logical,intent(out) :: included
    integer,intent(inout) :: istat

    included = .false.

    do i = 1,dict%nblk
      call blk%deallocate()
      blk = dict%blk_list(i)
      if (blk%header == 'basinhopping') then
        included = .true.
        call parse_bh_class(env,blk,bh,istat)
      end if
    end do
    return
  end subroutine parse_basinhopping_data
  subroutine parse_bh_class(env,blk,bh,istat)
    implicit none
    type(systemdata),intent(inout) :: env
    type(datablock),intent(in) :: blk
    type(bh_class),intent(inout) :: bh
    integer,intent(inout) :: istat
    integer :: i,j,nat
    logical :: rd
    if (blk%header .ne. 'basinhopping') return

    do i = 1,blk%nkv
      call parse_bh_auto(env,bh,blk%kv_list(i),rd)
      if (.not.rd) then
        istat = istat+1
        write (stdout,fmturk) '['//blk%header//']-block',blk%kv_list(i)%key
      end if
    end do
    return
  end subroutine parse_bh_class
  subroutine parse_bh_auto(env,bh,kv,rd)
    implicit none
    type(systemdata),intent(inout) :: env
    type(bh_class) :: bh
    type(keyvalue) :: kv
    logical,intent(out) :: rd
    logical,allocatable :: atlist(:)
    integer :: n,j
    logical :: ex
    rd = .true.

    select case (kv%key)
    case ('maxiter') !> these are NOT the BH steps!
      bh%maxiter = max(1,kv%value_i)

    case ('maxsave')
      bh%maxsave = kv%value_i

    case ('seed')
      if (.not.allocated(bh%seed)) allocate (bh%seed)
      bh%seed = kv%value_i

    case ('step','stepsize')
      select case (kv%id)
      case (valuetypes%int)
        bh%stepsize(1) = real(kv%value_i)
      case (valuetypes%float)
        bh%stepsize(1) = kv%value_f
      case (valuetypes%float_array)
        n = min(size(kv%value_fa,1),3)
        bh%stepsize(1:n) = kv%value_fa(1:n)
      case default
        !>--- keyword was recognized, but invalid argument supplied
        write (stdout,fmtura) kv%rawvalue
        call creststop(status_config)
      end select

    case ('steps','maxsteps')  !> these are the BH steps
      bh%maxsteps = kv%value_i

    case ('steptype')
      select case (kv%value_c)
      case ('cartesian')
        bh%steptype = 0
      case ('internal')
        bh%steptype = 1
      case ('dihedral')
        bh%steptype = 2
      case ('intermol')
        bh%steptype = 3
      case default
        write (stdout,fmtura) trim(kv%value_c)
        call creststop(status_config)
      end select

    case ('temp','T')
      bh%temp = kv%value_f

    case ('parallel')
      bh%parallel = kv%value_b

    case default
      rd = .false.
      return
    end select

  end subroutine parse_bh_auto

!========================================================================================!
!========================================================================================!
end module parse_calcdata
