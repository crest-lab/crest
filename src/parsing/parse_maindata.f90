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
  end subroutine parse_main_auto
  subroutine parse_main_float(env,key,val)
    implicit none
    type(systemdata) :: env
    character(len=*) :: key
    real(wp) :: val
    select case (key)

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
      case ('ancopt','optimize')
        env%preopt = .false.
        env%crestver = crest_optimize
      case ('ancopt_ensemble','optimize_ensemble')
        env%preopt = .false.
        env%crestver = crest_mdopt2
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
    case ('ensemble_input','ensemble')
      env%ensemblename = val
      env%inputcoords = val
    case ('input')
      env%inputcoords = val
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
    end select
    return
  end subroutine parse_main_bool
!========================================================================================!
  subroutine parse_main_blk(env,blk)
    implicit none
    type(systemdata) :: env
    type(datablock) :: blk
    select case (blk%header)
    case ('cregen')

    end select
  end subroutine parse_main_blk
!========================================================================================!
end module parse_maindata
