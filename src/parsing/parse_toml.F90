!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Philipp Pracht
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

module parse_toml
  use crest_parameters
  implicit none
  private

  public :: parse_tomlf
  public :: parse_toml_input_fallback

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

!> Handle calls to the toml-f library with preprocessor statements
#ifndef WITH_TOMLF

  subroutine parse_tomlf(fname,dict)
    use parse_datastruct,only:root_object
    implicit none
    character(len=*) :: fname !> name of the input file
    type(root_object),intent(inout) :: dict
    write (stdout,*) 'Error: Compiled without https://github.com/toml-f/toml-f support!'
    write (stdout,*) 'Use -DWITH_TOMLF=true in the setup to enable this function'
    error stop
  end subroutine parse_tomlf

#else /* WITH_TOMLF */

  subroutine parse_tomlf(fname,dict)
    use tomlf
    use parse_datastruct,only:root_object
    implicit none
    !> INPUT
    character(len=*) :: fname !> name of the input file
    !> OUTPUT
    type(root_object),intent(inout) :: dict
    !> LOCAL
    logical,parameter :: color = .true.
    type(toml_table),allocatable :: table
    type(toml_error),allocatable :: error
    type(toml_context) :: context
    type(toml_key),allocatable :: list(:)
    integer :: ikey,nkey
    type(toml_table),pointer :: child
    integer :: depth

!>--- initialize CREST's keyword dictionary object
    call dict%new()

!>--- read the toml file into table
    call toml_load(table,fname,error=error,context=context, &
    & config=toml_parser_config(color=color))
    call handle_tomlf_error(error)

!>--- iterate over keys (recursively)
    depth = 0
    call parse_tomlf_subtable(table,dict,depth)

!>--- do some formatting
    call dict%lowercase_keys()
    return
  end subroutine parse_tomlf
!========================================================================================!
!> handler for toml-f errors
  subroutine handle_tomlf_error(error)
    use tomlf,only:toml_error
    type(toml_error),intent(in),optional :: error
    if (present(error)) then
      write (stderr,'(a)') error%message
      stop 1
    end if
  end subroutine handle_tomlf_error
!=======================================================================================!
!>--- recursive handler for subtables
  recursive subroutine parse_tomlf_subtable(table,dict,depth,parentname)
    use tomlf
    use parse_datastruct,only:root_object
    use parse_keyvalue
    use parse_block
    implicit none
    !> INPUT
    type(toml_table) :: table
    type(root_object),intent(inout) :: dict
    integer :: depth,current_depth
    character(len=*),optional,intent(in) :: parentname
    !> LOCAL
    type(toml_error),allocatable :: error
    type(toml_context) :: context
    type(toml_key),allocatable :: list(:)
    integer :: ikey,nkey
    type(toml_table),pointer :: child
    type(toml_array),pointer :: array
    integer :: depth2,io,typ
    type(keyvalue) :: kvdum
    type(datablock) :: blkdum
    character(len=:),allocatable :: tablename

    current_depth = depth
!>--- prepare
    if (depth .ne. 0) then
      if (present(parentname)) then
        tablename = trim(parentname)//'.'//trim(table%key)
      else
        tablename = trim(table%key)
      end if
      if (.not.allocated(dict%blk_list(current_depth)%header)) then
        dict%blk_list(current_depth)%header = tablename
      else
        tablename = dict%blk_list(current_depth)%header
      end if
    end if

    call table%get_keys(list)
    nkey = size(list,1)
    do ikey = 1,nkey
      call get_value(table,list(ikey),child)
      call get_value(table,list(ikey),array)

!>--- "normal" key-value pairs
      if (.not.associated(child).and..not.associated(array)) then
        typ = parse_tomlf_checktype_normal(table,list(ikey))
        call parse_tomlf_kv_normal(table,list(ikey),typ,kvdum)
        if (current_depth == 0) then
          call dict%addkv(kvdum) !> add to root dict if global key value
        else
          call dict%blk_list(current_depth)%addkv(kvdum) !> add key-value to block
        end if

!>--- arrays
      else if (associated(array)) then
        typ = parse_tomlf_checktype_array(array)
        select case (typ)
        case (-1,-2) !> array of subtables or nested array
          call parse_tomlf_subtable_or_nested_array(array,dict,depth,tablename)

        case default !> typical value array
          call parse_tomlf_kv_array(array,typ,kvdum)
          if (current_depth == 0) then
            call dict%addkv(kvdum) !> add to root dict if global key value
          else
            call dict%blk_list(current_depth)%addkv(kvdum) !> add key-value to block
          end if
        end select

!>--- subtables
      else if (associated(child)) then
        call dict%addblk(blkdum) !> add new empty block to dict
        depth = depth+1
        call parse_tomlf_subtable(child,dict,depth,tablename)

      end if
    end do

  end subroutine parse_tomlf_subtable
!=======================================================================================!
!>-- routines to  convert toml-f key-value pairs into the crest datastructure
  function parse_tomlf_checktype_normal(table,key) result(typ)
    use tomlf
    use parse_keyvalue,only:valuetypes
    implicit none
    integer :: typ
    type(toml_table) :: table
    type(toml_key)   :: key
    integer :: io
    logical :: bool_dum
    integer :: int_dum
    real(wp) :: float_dum
    character(len=:),allocatable :: str_dum
    typ = valuetypes%unknown
    call get_value(table,key,bool_dum,stat=io)
    if (io == 0) then
      typ = valuetypes%bool
      return
    end if
    call get_value(table,key,int_dum,stat=io)
    if (io == 0) then
      typ = valuetypes%int
      return
    end if
    call get_value(table,key,float_dum,stat=io)
    if (io == 0) then
      typ = valuetypes%float
      return
    end if
    call get_value(table,key,str_dum,stat=io)
    if (io == 0) then
      typ = valuetypes%string
      return
    end if
  end function parse_tomlf_checktype_normal
  subroutine parse_tomlf_kv_normal(table,key,typ,kv)
    use tomlf
    use parse_keyvalue
    implicit none
    type(toml_table) :: table
    type(toml_key)   :: key
    integer :: typ
    type(keyvalue) :: kv

    call kv%deallocate()

    kv%id = typ
    kv%key = key%key
    select case (typ)
    case (valuetypes%bool)
      call get_value(table,key,kv%value_b)
    case (valuetypes%int)
      call get_value(table,key,kv%value_i)
      kv%value_f = real(kv%value_i)
    case (valuetypes%float)
      call get_value(table,key,kv%value_f)
      kv%value_i = nint(kv%value_f)
    case (valuetypes%string)
      call get_value(table,key,kv%value_c)
    end select

    call kv%set_valuestring()
  end subroutine parse_tomlf_kv_normal
!=======================================================================================!
!>-- routines to convert toml-f key-value pairs with an array type into the crest datastructure
  function parse_tomlf_checktype_array(tarray) result(typ)
    use tomlf
    use parse_keyvalue,only:valuetypes
    implicit none
    integer :: typ
    type(toml_array) :: tarray
    integer :: io,t_len,i
    logical,allocatable :: bool_dum(:)
    integer,allocatable :: int_dum(:)
    real(wp),allocatable :: float_dum(:)
    character(len=:),allocatable :: str_dum(:)

    type(toml_array),pointer ::  nested_array
    type(toml_table),pointer :: child

    typ = valuetypes%unknown
    t_len = len(tarray)
    !> Get a pointer to the 1st element to check type
    call get_value(tarray,1,nested_array)
    call get_value(tarray,1,child)
    if (.not.associated(nested_array).and..not.associated(child)) then

      call get_value(tarray,bool_dum,stat=io)
      if (io == 0) then
        typ = valuetypes%bool_array
        deallocate (bool_dum)
        return
      end if
      call get_value(tarray,int_dum,stat=io)
      if (io == 0) then
        typ = valuetypes%int_array
        deallocate (int_dum)
        return
      end if
      call get_value(tarray,float_dum,stat=io)
      if (io == 0) then
        typ = valuetypes%float_array
        deallocate (float_dum)
        return
      end if
      typ = valuetypes%raw_array

    else if (associated(nested_array)) then
      typ = -1
    else if (associated(child)) then
      typ = -2
    end if

  end function parse_tomlf_checktype_array
  subroutine parse_tomlf_kv_array(tarray,typ,kv)
    use tomlf
    use parse_keyvalue
    implicit none
    type(toml_array) :: tarray
    integer :: typ
    type(keyvalue) :: kv

    call kv%deallocate()

    kv%id = typ
    kv%key = tarray%key
    select case (typ)
    case (valuetypes%bool_array)
      call get_value(tarray,kv%value_ba)
      kv%na = size(kv%value_ba,1)

    case (valuetypes%int_array)
      call get_value(tarray,kv%value_ia)
      kv%na = size(kv%value_ia,1)

    case (valuetypes%float_array)
      call get_value(tarray,kv%value_fa)
      kv%na = size(kv%value_fa,1)

    case default
      !> mixed type arrays and arrays made of strings
      call parse_tomlf_rawarray_helper(tarray,typ,kv)

    end select

    call kv%set_valuestring()
  end subroutine parse_tomlf_kv_array
  subroutine parse_tomlf_rawarray_helper(tarray,typ,kv)
    use tomlf
    use parse_keyvalue
    use iomod,only:to_str
    implicit none
    type(toml_array) :: tarray
    integer :: typ,io
    type(keyvalue),intent(inout) :: kv
    character(len=:),allocatable :: str_dum
    integer :: int_dum
    real(wp) :: float_dum
    logical :: bool_dum
    integer :: t_len,i
    logical :: all_strings

    all_strings = .true.
!>--- iterate pointers
    t_len = len(tarray)
    do i = 1,t_len
      call get_value(tarray,i,bool_dum,stat=io)
      if (io == 0) then
        all_strings = .false.
        call kv%add_raw_array_string(to_str(bool_dum))
        cycle
      end if
      call get_value(tarray,i,int_dum,stat=io)
      if (io == 0) then
        all_strings = .false.
        call kv%add_raw_array_string(to_str(int_dum))
        cycle
      end if
      call get_value(tarray,i,float_dum,stat=io)
      if (io == 0) then
        all_strings = .false.
        call kv%add_raw_array_string(to_str(float_dum))
        cycle
      end if
      call get_value(tarray,i,str_dum,stat=io)
      if (io == 0) then
        call kv%add_raw_array_string(str_dum)
        cycle
      end if
    end do
    !if (kv%na .gt. 0) then
    !  do i = 1,kv%na
    !    write (*,*) kv%value_rawa(i)
    !  end do
    !end if
    if (all_strings) then
      kv%id = valuetypes%string_array
      kv%value_ca = kv%value_rawa
    end if
  end subroutine parse_tomlf_rawarray_helper
!========================================================================================!
!> special case that need recursive handling: an array of subtables or nested arrays
!> nested arrays will be added as duplicate key-value objects to the dict or latest blk
  recursive subroutine parse_tomlf_subtable_or_nested_array(tarray,dict,depth,parentname)
    use tomlf
    use parse_datastruct,only:root_object
    use parse_keyvalue
    use parse_block
    implicit none
    !> INPUT
    type(toml_array) :: tarray
    type(root_object),intent(inout) :: dict
    integer :: depth
    character(len=*),optional,intent(in) :: parentname
    !> LOCAL
    integer :: ikey,nkey
    type(toml_table),pointer :: child
    type(toml_array),pointer :: nested_array
    integer :: depth2,io,typ,t_len,i,current_depth
    type(keyvalue) :: kvdum
    type(datablock) :: blkdum
    character(len=:),allocatable :: tablename

!>--- prepare
    current_depth = depth
    if (depth .ne. 0) then
      if (present(parentname)) then
        tablename = trim(parentname)//'.'//trim(tarray%key)
      else
        tablename = trim(tarray%key)
      end if
    end if

    t_len = len(tarray)
!>--- iterate pointers
    do i = 1,t_len
      call get_value(tarray,i,child)
      call get_value(tarray,i,nested_array)
      if (associated(child)) then !> subtables
        blkdum%header = tablename
        call dict%addblk(blkdum)

        depth = depth+1
        call parse_tomlf_subtable(child,dict,depth,parentname)

      else if (associated(nested_array)) then !> nested arrays
        typ = parse_tomlf_checktype_array(nested_array)
        select case (typ)
        case (-1,-2)
          call parse_tomlf_subtable_or_nested_array(nested_array,dict,depth)
        case default
          call parse_tomlf_kv_array(nested_array,typ,kvdum)
          kvdum%key = tarray%key
          if (current_depth == 0) then
            call dict%addkv(kvdum) !> add to root dict if global key value
          else
            call dict%blk_list(current_depth)%addkv(kvdum) !> add key-value to block
          end if
        end select
      end if

      nullify (nested_array)
      nullify (child)
    end do

  end subroutine parse_tomlf_subtable_or_nested_array

!=======================================================================================!
#endif /* WITH_TOMLF */

!========================================================================================!
!========================================================================================!
!> This is the fallback reader for toml files >> IF TOML-F IS NOT USED <<.
!> It can handle basic toml files but is not as refined as the dedicated library!
!> For example, multiline comments are not read correctly
  subroutine parse_toml_input_fallback(fname,dict)
    use parse_keyvalue
    use parse_block
    use parse_datastruct
    use filemod
    implicit none

    character(len=*) :: fname !> name of the input file
    type(root_object),intent(out) :: dict
    type(filetype) :: file
    integer :: i,j,k,io
    logical :: get_root_kv
    type(keyvalue) :: kvdum
    type(datablock) :: blkdum

    call dict%new()
!>--- open file to read and remove comments
    call file%open(trim(fname))
    dict%filename = trim(file%filename)
    call remove_comments(file)

    get_root_kv = .true.
!>--- the loop where the input file is read
    do i = 1,file%nlines
      if (file%current_line > i) cycle

      !> first, get key-value pairs of the root dict
      if (get_root_kv) then
        call get_keyvalue(kvdum,file%line(i),io)
        if (io == 0) then
          call dict%addkv(kvdum) !> add to dict
        end if
      end if

      if (isheader(file%line(i))) then
        get_root_kv = .false.
        call read_datablock(file,i,blkdum)
        call blkdum%fmt_header()
        call dict%addblk(blkdum) !> add to dict
      end if
    end do

    call file%close()

    return
  end subroutine parse_toml_input_fallback
!========================================================================================!
  subroutine remove_comments(file)
    use filemod
    implicit none
    type(filetype) :: file
    character(len=1),parameter :: com = '#'
    integer :: i
    do i = 1,file%nlines
      call clearcomment(file%f(i),com)
    end do
  end subroutine remove_comments
!========================================================================================!
end module parse_toml
