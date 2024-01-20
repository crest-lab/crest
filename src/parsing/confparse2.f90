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
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> INPUT FILE PARSER FOR CREST
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!

subroutine parseinputfile(env,fname)
!****************************************************
!* A supplement to the parseflags routine
!* in confparse.f90. This routine reads
!* an input file (TOML format via the toml-f library)
!*
!* Input/Output:
!*  env   -  crest's systemdata object, which
!*           contains basically all information
!*           for the calculation
!*  fname -  name of the input file
!****************************************************
  use crest_parameters
  !> modules for data storage in crest
  use crest_data
  use crest_calculator,only:calcdata
  use dynamics_module,only:mddata

  !> modules used for parsing the root_object
  use parse_keyvalue,only:keyvalue
  use parse_block,only:datablock
  use parse_datastruct,only:root_object
  use parse_maindata
  use parse_inputfile,only:parse_input
  use parse_calcdata,only:parse_calculation_data, &
  &                         parse_dynamics_data
  !> Declarations
  implicit none
  type(systemdata),intent(inout) :: env
  character(len=*),intent(in) :: fname
  type(root_object) :: dict
  type(keyvalue)  :: kv
  type(datablock) :: blk
  type(calcdata) :: newcalc
  type(mddata) :: mddat
  logical :: ex,l1,l2
  integer :: i,j,k,l

!>--- check for the input files' existence
  inquire (file=fname,exist=ex)
  if (.not.ex) then
    return
  else
    write (stdout,*) 'reading ',trim(fname)
    env%legacy = .false.
  end if

!>--- read the file into the object 'dict'
  call parse_input(fname,dict)
  call dict%print()

!>--- parse all root-level key-value pairs
  do i = 1,dict%nkv
    kv = dict%kv_list(i)
    call parse_main_auto(env,kv)
  end do
!>------------------------------------------------------
!> After this point I assume an input structure was
!> provided. Either in parse_main_auto or confparse.f90
!>------------------------------------------------------

!>--- parse all objects that write to env or global data
  do i = 1,dict%nblk
    blk = dict%blk_list(i)
    call parse_main_blk(env,blk)
  end do

!>--- check objects for a calculation setup
!     i.e., all [calculation] and [[calculation.*]] blocks
  call parse_calculation_data(env,newcalc,dict,l1)
  if (l1) then
    env%calc = newcalc
    call env_calcdat_specialcases(env)
  end if

!>--- check for molecular dynamics setup
!     i.e., all [dynamics] and [[dynamics.*]] blocks
  call parse_dynamics_data(env,mddat,dict,l1)
  if (l1) then
    env%mddat = mddat
  end if

!>--- check for lwONIOM setup (will be read at end of confparse)
  do i = 1,dict%nblk
    if (dict%blk_list(i)%header == 'lwoniom') then
      env%ONIOM_toml = trim(fname)
    end if
  end do

  call dict%deallocate()
  return
end subroutine parseinputfile

!========================================================================================!
!========================================================================================!
subroutine internal_constraint_repair(env,bondconst)
!**************************************************
!* subroutine internal_constraint_setup
!* add global settings from crests legacy keywords
!* to the new calculator setup
!**************************************************
  use crest_parameters
  use crest_data
  use parse_xtbinput
  implicit none
  type(systemdata) :: env
  logical,intent(in) :: bondconst
  integer :: i,j,k,l,n
  integer :: nat
  logical,allocatable :: atms(:)
  logical :: ex

!>--- don't need it for the legacy runtype
  if (env%legacy) return

!>--- check if we want to read the bondlengths file
  inquire (file='bondlengths',exist=ex)
  if (bondconst.and.ex) then
    call parse_xtbinputfile(env,'bondlengths')
  end if

!>--- wall potential setup for calculator version
  if (env%NCI.or.env%wallsetup) then
    write (stdout,'("> ",a)') 'Generating logfermi wall potential for the system.'
    if (allocated(env%potatlist)) then
      write (stdout,'("> ",a,f8.3,2x,a,f8.3,a,a)') 'wscal=',env%potscal, &
      & 'wpad=',env%potpad,'  atoms: ',trim(env%potatlist)
    else
      write (stdout,'("> ",a,f8.3,2x,a,f8.3,a,a)') 'wscal=',env%potscal, &
      & 'wpad=',env%potpad,'  atoms: ','all'
    end if
    call wallpot(env)
  end if

  return
end subroutine internal_constraint_repair

!========================================================================================!
!========================================================================================!
subroutine env_calcdat_specialcases(env)
!****************************************************
!* Some special treatments are sometimes necessary
!* depending on a choosen calculation level,
!* for example, a MD timestep adjustment for GFN-FF.
!* This routine takes care of that,
!* but some stuff (in particular MD settings) may be
!* overwritten later in the respective block.
!***************************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  implicit none
  type(systemdata) :: env
  integer :: i,j,k,l
  integer :: refine_lvl

  !> special case for GFN-FF calculations
  if (any(env%calc%calcs(:)%id == jobtype%gfnff)) then
    env%mdstep = 1.5d0
    env%hmass = 5.0d0
    env%cts%cbonds_md = .true.
    env%checkiso = .true.
  end if

  !> check if any refinement is to be done between opt.
  if (any(env%calc%calcs(:)%refine_lvl > 0)) then
    do i = 1,env%calc%ncalculations
      refine_lvl = env%calc%calcs(i)%refine_lvl
      if (refine_lvl <= 0) cycle
      if(allocated(env%refine_queue))then
        if (any(env%refine_queue(:) == refine_lvl)) cycle
      endif
      call env%addrefine(refine_lvl)
    end do
  end if

end subroutine env_calcdat_specialcases
