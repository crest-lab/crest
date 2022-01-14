!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022 Philipp Pracht
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
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> INPUT FILE PARSER FOR CREST
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
!> A supplement to the parseflags routine
!> in confparse.f90. This routine reads
!> an input file (very loosly based on the
!> TOML format)
!>
!> Input/Output:
!>  env   -  crest's systemdata object, which
!>           contains basically all information
!>           for the calculation
!>  iname -  name of the input file
!>-----------------------------------------------
subroutine parseinputfile(env,fname)
  use iso_fortran_env,only:wp=>real64
  !> modules for data storage in crest
  use crest_data
  use calc_type,only:calcdata
  use dynamics_module,only:mddata
 
  !> modules used for parsing the root_object
  !>   
  use parse_keyvalue,only:keyvalue
  use parse_block,only:datablock
  use parse_datastruct,only:root_object
  use parse_inputfile, only: parse_input
  use parse_calcdata, only: parse_calculation_data
  !> Declarations
  implicit none
  type(systemdata),intent(inout) :: env
  character(len=*),intent(in) :: fname
  type(root_object) :: dict
  type(keyvalue)  :: kv
  type(datablock) :: blk
  type(calcdata) :: newcalc
  logical :: ex,l1,l2
  integer :: i,j,k,l

!>--- check for the input files' existence
  inquire (file=fname,exist=ex)
  if (.not. ex) then
    return
  else
    write (*,*) 'reading ',trim(fname)
  end if

!>--- read the file into the object 'dict'
  call parse_input(fname,dict)
  call dict%print()

!>--- parse all root-level key-value pairs
  do i=1,dict%nkv
    kv = dict%kv_list(i)
    call parse_main_auto(env,kv)
  enddo

!>--- parse all objects
  do i=1,dict%nblk
    blk = dict%blk_list(i)
    call parse_main_blk(env,blk)
  enddo
!>--- check objects for a calculation setup
  call parse_calculation_data(newcalc,dict,l1)
  if(l1)then
    env%calc = newcalc
  endif


  call dict%deallocate()
  return
contains
!========================================================================================!
  subroutine parse_main_auto(env,kv)
    implicit none
    type(systemdata) :: env
    type(keyvalue) :: kv
    select case (kv%id)
    case (1) !> float
      call parse_main_float(env,kv%key,kv%value_f)
    case (2) !> int
      call parse_main_int(env,kv%key,kv%value_i)
    case (3) !> bool
      call parse_main_bool(env,kv%key,kv%value_b)
    case (4) !> string
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
     case('threads','parallel')
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
      case default
        env%crestver = crest_imtd
      end select
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
    case( 'noopt' )
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
end subroutine parseinputfile
