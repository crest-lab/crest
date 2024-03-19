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
!> Routines contained here are for parsing 'top level' settings that will
!> enter the env systemdata object

module parse_maindata
  use crest_parameters
  !> modules for data storage in crest
  use crest_data
  use crest_restartlog
  use strucrd,only:coord
  !> modules used for parsing the root_object
  !>
  use parse_keyvalue,only:keyvalue,valuetypes
  use parse_block,only:datablock
  use parse_datastruct,only:root_object
  !> Declarations
  implicit none
  public

!========================================================================================!
!========================================================================================!
contains   !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine parse_main_auto(env,kv)
    implicit none
    type(systemdata) :: env
    type(keyvalue) :: kv
    select case (kv%id)
    case (valuetypes%float) !> float
      call parse_main_float(env,kv%key,kv%value_f)
    case (valuetypes%int) !> int
      call parse_main_int(env,kv%key,kv%value_i)
    case (valuetypes%bool) !> bool
      call parse_main_bool(env,kv%key,kv%value_b)
    case (valuetypes%string) !> string
      call parse_main_c(env,kv%key,kv%value_c)
    end select
!> other, with multiple or raw type
    select case (kv%key)
    case ('optlev','ancopt_level')
      env%optlev = optlevnum(kv%rawvalue)
    end select
  end subroutine parse_main_auto
  subroutine parse_main_float(env,key,val)
    implicit none
    type(systemdata) :: env
    character(len=*) :: key
    real(wp) :: val
    select case (key)
    case ('wscal')
      env%potscal = val
      env%wallsetup = .true.
    case ('wpad')
      env%potpad = val
      env%wallsetup = .true.

    end select
    return
  end subroutine parse_main_float
  subroutine parse_main_int(env,key,val)
    implicit none
    type(systemdata) :: env
    character(len=*) :: key
    integer :: val
    select case (key)
    case ('threads','parallel')
      env%Threads = val
      env%autothreads = .true.
      env%threadssetmanual = .true.
    end select
    return
  end subroutine parse_main_int
  subroutine parse_main_c(env,key,val)
    implicit none
    type(systemdata) :: env
    character(len=*) :: key
    character(len=*) :: val
    type(coord) :: mol
    select case (key)
    case ('bin','binary')
      env%ProgName = val
    case ('runtype')
      select case (val)
      case ('none')
        env%crestver = crest_none
      case ('playground','test')
        env%preopt = .false.
        env%crestver = crest_test
      case ('singlepoint','sp')
        env%preopt = .false.
        env%crestver = crest_sp
      case ('numgrad')
        env%preopt = .false.
        env%crestver = crest_sp
        env%testnumgrad = .true.
      case ('ancopt','optimize','ohess')
        env%preopt = .false.
        env%crestver = crest_optimize
        env%optlev = 0.0_wp
        if(val.eq.'ohess')then
          env%crest_ohess=.true.
        endif   
      case ('ancopt_ensemble','optimize_ensemble','mdopt')
        env%preopt = .false.
        env%crestver = crest_mdopt2
        env%optlev = 0.0d0
      case ('screen_ensemble','screen')
        env%preopt = .false.
        env%crestver = crest_screen
      case ('ensemble_singlepoints','ensemblesp','mdsp')
        env%preopt = .false.
        env%crestver = crest_ensemblesp

      case ('md','mtd','metadynamics','dynamics')
        env%preopt = .false.
        env%crestver = crest_moldyn
      case ('scan')
        env%preopt = .false.
        env%crestver = crest_scanning
      case ('search_1')
        env%preopt = .false.
        env%crestver = crest_s1
        env%runver = crest_s1
      case ('mecp','mecp_search')
        env%preopt = .false.
        env%crestver = crest_mecp
        env%runver = crest_mecp
      case ('imtd-gc')
        env%preopt = .false.
        env%crestver = crest_imtd
        env%runver = 1
      case ('nci-mtd','nci')
        env%NCI = .true.
        env%runver = 4
        env%autozsort = .false.
        env%performCross = .false.
        env%rotamermds = .false.
      case ('entropy','imtd-stmd')
        env%crestver = crest_imtd  !> the entropy mode acts as subtype of the crest_imtd algo
        env%properties = abs(p_CREentropy)
        env%autozsort = .false.     !> turn off zsort (since we are not going to GC anyways)
        env%performCross = .false.  !> turn off GC
        env%entropic = .true.       !> indicator for this runtype
        env%Maxrestart = 1          !> turn off MTD iterations (just do one)
        env%rotamermds = .false.    !> turn off normMDs
        env%entropymd = .true.      !> special static MTDs
        env%runver = 111            !> version  for selection of MTD bias settings
        env%doNMR = .true.          !> we need equivalencies
        env%emtd%bhess = .false.    !> currently there is no BHESS version, TODO!
        call env%addjob(env%properties)
      case ('numhess','numerical hessian')
        env%preopt = .false.
        env%crestver = crest_numhessian
        env%runver = crest_numhessian
      case ('rigidconf')
        env%preopt = .false.
        env%crestver = crest_rigcon
        env%runver = crest_rigcon
      case default
        env%crestver = crest_imtd
      end select
    case ('ensemble_input','ensemble','input_ensemble')
      env%ensemblename = val
      env%inputcoords = val
    case ('input','structure')
      env%inputcoords = val
      call mol%open(val)
      call env%ref%load(mol)

    case ('constraints','xtbconstraints','xtbinput') !> equivalent to --cinp
      env%constraints = val
    case ('rigidconf_file')
      env%rigidconf_userfile = val

    case ('watlist','wat')
      env%potatlist = val
      env%wallsetup = .true.
    end select
    return
  end subroutine parse_main_c
  subroutine parse_main_bool(env,key,val)
    implicit none
    type(systemdata) :: env
    character(len=*) :: key
    logical :: val
    select case (key)
    case ('preopt')
      env%preopt = val
    case ('noopt')
      env%preopt = .not.val
    case ('topo')
      env%checktopo = val
    case ('notopo')
      env%checktopo = .not.val
    case ('restart')
      if (val) then
        call read_restart(env)
      end if
    case ('multilevelopt')
      env%multilevelopt = val 
    case ('refine_presort')
      env%refine_presort = val
    end select
    return
  end subroutine parse_main_bool
!========================================================================================!

  subroutine parse_main_blk(env,blk)
!**************************************
!* Some shorter blocks are not defined
!* in separate source files. They can
!* be found below.
!**************************************
    implicit none
    type(systemdata) :: env
    type(datablock) :: blk
    select case (blk%header)
    case ('cregen')
      call parse_cregen(env,blk)
    case ('confsolv')
      call parse_confsolv(env,blk)
    case ('thermo')
      call parse_thermo(env,blk)
    end select
  end subroutine parse_main_blk

!========================================================================================!
  subroutine parse_cregen(env,blk)
!****************************************
!* parse settings for the CREGEN routine
!****************************************
    implicit none
    type(systemdata) :: env
    type(datablock) :: blk
    type(keyvalue) :: kv
    integer :: i
!>--- parse the arguments
    do i = 1,blk%nkv
      kv = blk%kv_list(i)
      select case (kv%key)
      case ('ewin')
        env%ewin = kv%value_f
      case ('ethr')
        env%ethr = kv%value_f
      case ('rthr')
        env%rthr = kv%value_f
      case ('bthr')
        env%bthr2 = kv%value_f
      case ('eqv','nmr')
        env%doNMR = kv%value_b
      end select
    end do
  end subroutine parse_cregen

!========================================================================================!
  subroutine parse_confsolv(env,blk)
    use ConfSolv_module
    implicit none
    type(systemdata) :: env
    type(datablock) :: blk
    type(keyvalue) :: kv
    integer :: i
!>--- add ConfSolv as refinement level to give a ΔΔGsoln
    call env%addrefine(refine%ConfSolv)
    !env%refine_presort = .true.
    env%ewin = 100.0_wp

!>--- parse the arguments
    do i = 1,blk%nkv
      kv = blk%kv_list(i)
      select case (kv%key)
      case ('pid')
        if (.not.allocated(cs_pid)) allocate (cs_pid)
        cs_pid = kv%value_i
      case ('bin')
        cs_bin = trim(kv%value_c)
      case ('port')
        if (.not.allocated(cs_port)) allocate (cs_port)
        cs_port = kv%value_i
      case ('solvent')
        !> to define a single solvent like: solvent = ['water','O']
        if (kv%na == 2) then
          cs_solvent = trim(kv%value_rawa(1))
          cs_smiles = trim(kv%value_rawa(2))
        else if(index(kv%value_c,'.csv').ne.0)then
          cs_solvfile = kv%value_c
        else
          cs_solvent = kv%value_c
        end if
      case('solvent_csv','solvfile')
        cs_solvfile = kv%value_c
      case ('solvent_name')
        cs_solvent = kv%value_c
      case ('solvent_smiles')
        cs_smiles = kv%value_c
      case ('model_path','param','checkpoint')
        cs_param = kv%value_c
      end select
    end do
  end subroutine parse_confsolv

!========================================================================================!

  subroutine parse_thermo(env,blk)
!****************************************
!* parse settings for the Thermo routine
!****************************************
    implicit none
    type(systemdata) :: env
    type(datablock) :: blk
    type(keyvalue) :: kv
    integer :: i
!>--- parse the arguments
    do i = 1,blk%nkv
      kv = blk%kv_list(i)
      select case (kv%key)
      case ('ithr','freq_ithr','freq_invert')
        env%thermo%ithr = kv%value_f
      case ('fscal','freq_scal')
        env%thermo%fscal = kv%value_f
      case ('sthr','freq_interpol')
        env%thermo%sthr = kv%value_f
      case ('trange')
        if (kv%na >= 2) then
          env%thermo%trange(1) = minval(kv%value_fa(1:2),1)
          env%thermo%trange(2) = maxval(kv%value_fa(1:2),1) 
        endif  
        if(kv%na >= 3)then
          env%thermo%trange(3) = kv%value_fa(3)
        endif 
      case ('tstep')
        env%thermo%trange(3) = kv%value_f

      case ('input','coords')
        env%thermo%coords = kv%value_c
        if(allocated(env%thermo%vibfile)) env%properties = p_thermo
      case ('freq_input','vibs','hessian')
        env%thermo%vibfile = kv%value_c 
        if(allocated(env%thermo%coords)) env%properties = p_thermo
      end select
    end do
  end subroutine parse_thermo

!========================================================================================!
!========================================================================================!
end module parse_maindata
