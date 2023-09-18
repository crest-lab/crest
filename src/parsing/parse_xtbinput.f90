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

!> This is the fallback reader for xtb input files.
!> Defining constrains in this format is a bit easier than with toml
!> Furthermore, it should provide some higher degree of compatibility
!> between CREST and xTB.
!> The xtb-style keywords in CREST are limited to geometrical constraints
!> These are files that can be read with the --cinp option

module parse_xtbinput
  use crest_parameters
  use crest_data
  use parse_datastruct
  use parse_keyvalue
  use parse_block
  use parse_datastruct
  use filemod
  use iomod
  use strucrd
  use wall_setup 
  use constraints,only:constraint
  implicit none
  private
  !logical,parameter,private :: debug = .true.
  logical,parameter,private :: debug = .false.

  public :: parse_xtbinputfile
  interface parse_xtbinputfile
    module procedure :: parse_xtb_inputfile
    module procedure :: parse_xtb_input_fallback
  end interface parse_xtbinputfile

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine parse_xtb_inputfile(env,fname)
!*********************************************
!* Routine for parsing the input file fname
!* and storing information in env
!*********************************************
    implicit none
    type(systemdata) :: env
    character(len=*) :: fname

    type(root_object) :: dict
    type(datablock) :: blk
    logical :: ex
    character(len=:),allocatable :: hdr
    integer :: i,j,k,l

    inquire (file=fname,exist=ex)
    if (.not.ex) return

    call parse_xtb_input_fallback(fname,dict)
    !call dict%print()

    write (stdout,'(a,a,a)') 'Parsing xtb-type input file ',trim(fname), &
    & ' to set up calculators ...'
    !> iterate through the blocks and save the necessary information
    do i = 1,dict%nblk
      blk = dict%blk_list(i)
      hdr = trim(blk%header)
      select case (hdr)
      case ('constrain')
        call get_xtb_constraint_block(env,blk)
      case ('wall')
        call get_xtb_wall_block(env,blk)
      case ('fix')
        call get_xtb_fix_block(env,blk)
      case default
        write (stdout,'(a,a,a)') 'xtb-style input block: "$',trim(hdr),'" not defined for CREST'
      end select
    end do

    if (debug) stop
  end subroutine parse_xtb_inputfile

!========================================================================================!

  subroutine get_xtb_constraint_block(env,blk)
!********************************************************************
!* This is the fallback reader for xtb input files to set up a dict
!********************************************************************
    implicit none
    type(systemdata),intent(inout) :: env
    type(filetype) :: file
    integer :: i,j,k,io
    type(keyvalue) :: kv
    type(datablock),intent(in) :: blk
    real(wp) :: force_constant,dist,angl
    real(wp) :: rdum
    type(coord) :: mol
    type(coord) :: molref
    logical :: useref
    logical,allocatable :: pairwise(:)
    logical,allocatable :: atlist(:)
    integer :: i1,i2,i3,i4
    type(constraint) :: cons

    useref = .false.
!>--- a default xtb force constant (in Eh), must be read first, if present
    force_constant = 0.05
    do i = 1,blk%nkv
      kv = blk%kv_list(i)
      select case (kv%key)

      case ('force constant')
        read (kv%rawvalue,*,iostat=io) rdum
        if (io == 0) then
          if (debug) write (stdout,'(a,a,a)') 'read force constant: ',to_str(rdum),' Eh'
          force_constant = rdum
        end if

      case ('reference')
      !> a reference geometry (must be the same molecule as the input)
        call molref%open(kv%rawvalue)
        if(any(mol%at(:) .ne. molref%at(:)))then
         write (stdout,'(a,/,a)') '**ERROR** while reading xtb-style input:',&
         & '  Geometry provided as "reference=" appears not to be the same molecule as CREST input!'
         error stop
        endif 
        useref = .true.
        if(debug) write(stdout,'(a,a,a)') '> Using reference geometry "',kv%rawvalue,'"'

      end select
    end do

!>--- get reference input geometry
    call env%ref%to(mol)

!>--- then the common constraints: distance, angle, dihedral
    do i = 1,blk%nkv
      kv = blk%kv_list(i)
      select case (kv%key)
      case ('force constant','reference')
        !> already read above

      case ('distance','bond')
        if (kv%na .eq. 3) then
          read (kv%value_rawa(1),*,iostat=io) i1
          if (io == 0) read (kv%value_rawa(2),*,iostat=io) i2
          if (io == 0) then
            if (trim(kv%value_rawa(3)) .eq. 'auto') then
              if(useref)then
                dist = molref%dist(i1,i2)
              else 
                dist = mol%dist(i1,i2)
              endif
            else
              read (kv%value_rawa(3),*) dist
              dist = dist*aatoau
            end if
            call cons%deallocate()
            call cons%bondconstraint(i1,i2,dist,force_constant)
            if (debug) call cons%print(stdout)
            call env%calc%add(cons)
          end if
        end if

      case ('angle')
        if (kv%na .eq. 4) then
          read (kv%value_rawa(1),*,iostat=io) i1
          if (io == 0) read (kv%value_rawa(2),*,iostat=io) i2
          if (io == 0) read (kv%value_rawa(3),*,iostat=io) i3
          if (io == 0) then
            if (trim(kv%value_rawa(4)) .eq. 'auto') then
              if(useref)then
                angl = molref%angle(i1,i2,i3)*radtodeg
              else
                angl = mol%angle(i1,i2,i3)*radtodeg
              endif
            else
              read (kv%value_rawa(4),*) angl
            end if
            call cons%deallocate()
            call cons%angleconstraint(i1,i2,i3,angl,force_constant)
            if (debug) call cons%print(stdout)
            call env%calc%add(cons)
          end if
        end if

      case ('dihedral')
        if (kv%na .eq. 5) then
          read (kv%value_rawa(1),*,iostat=io) i1
          if (io == 0) read (kv%value_rawa(2),*,iostat=io) i2
          if (io == 0) read (kv%value_rawa(3),*,iostat=io) i3
          if (io == 0) read (kv%value_rawa(4),*,iostat=io) i4
          if (io == 0) then
            if (trim(kv%value_rawa(5)) .eq. 'auto') then
              if(useref)then
                angl = molref%dihedral(i1,i2,i3,i4)*radtodeg
              else
                angl = mol%dihedral(i1,i2,i3,i4)*radtodeg
              endif
            else
              read (kv%value_rawa(5),*) angl
            end if
            call cons%deallocate()
            call cons%dihedralconstraint(i1,i2,i3,i4,angl,force_constant)
            if (debug) call cons%print(stdout)
            call env%calc%add(cons)
          end if
        end if

      case ('atoms')
        if (.not.allocated(pairwise)) allocate (pairwise(mol%nat),source=.false.)
        call get_atlist(mol%nat, atlist, kv%rawvalue, mol%at)
        do j=1,mol%nat
          if(atlist(j)) pairwise(j) = .true.
        enddo
    
      case ('elements')
        if (.not.allocated(pairwise)) allocate (pairwise(mol%nat),source=.false.)
        if (kv%id == valuetypes%raw_array) then
          do j = 1,kv%na
            i1 = e2i(kv%value_rawa(j))
            do k=1,mol%nat
             if (i1 == mol%at(k)) pairwise(k) = .true.
            enddo
          enddo
        else
          i1 = e2i(kv%rawvalue)
          do j=1,mol%nat
            if (i1 == mol%at(j)) pairwise(j) = .true.
          enddo
        end if


      case default
        write (stdout,'(a,a,a)') 'xtb-style input key: "',kv%key,'" not defined for CREST'

      end select
    end do

!>--- if the pairwise section was allocated, set the distance constraints up here
    if (allocated(pairwise)) then
      do i = 1,mol%nat
        do j = 1,i-1
          if (pairwise(i).and.pairwise(j)) then
            if(useref)then
              dist = molref%dist(j,i)
            else 
              dist = mol%dist(j,i)
            endif
            call cons%deallocate()
            call cons%bondconstraint(j,i,dist,force_constant)
            if (debug) call cons%print(stdout)
            call env%calc%add(cons)
          end if
        end do
      end do
      deallocate (pairwise)
    end if

  end subroutine get_xtb_constraint_block

  subroutine get_xtb_wall_block(env,blk)
!**************************************
!* This is a reader for the $wall block
!***************************************
    implicit none
    type(systemdata),intent(inout) :: env
    type(filetype) :: file
    integer :: i,j,k,io
    type(keyvalue) :: kv
    type(datablock),intent(in) :: blk
    real(wp) :: force_constant,dist,angl
    real(wp) :: T,alpha,beta
    real(wp) :: rdum,rabc(3),r1,r2,r3
    type(coord) :: mol
    logical,allocatable :: pairwise(:)
    logical,allocatable :: atlist(:)
    integer :: i1,i2,i3,i4
    integer :: pot
    type(constraint) :: cons

!>--- asome defaults
    force_constant = 1.0_wp
    alpha = 30
    beta  = 6.0_wp
    T     = 300.0_wp
    pot = 1 !> 1= polynomial, 2= logfermi

!>--- get reference input geometry
    call env%ref%to(mol)

!>--- get the parameters first
    do i = 1,blk%nkv
      kv = blk%kv_list(i)
      select case (kv%key)
      case ('force constant')
        !> already read above
        read (kv%rawvalue,*,iostat=io) rdum
        if (io == 0) force_constant = rdum

      case('potential')
        if(trim(kv%rawvalue).eq.'logfermi')then
          pot = 2
        else
          pot = 1 
        endif

      case('alpha')
         read(kv%rawvalue,*,iostat=io) i1
         if(io == 0) alpha = i1

      case('beta')
         read(kv%rawvalue,*,iostat=io) i1
         if(io == 0) beta = i1

      case('temp')
         read(kv%rawvalue,*,iostat=io) i1
         if(io == 0) T = i1

      end select
   enddo

!>--- create the potentials
   do i = 1,blk%nkv
      kv = blk%kv_list(i)
      select case (kv%key)
      case('force constant','potential','alpha','beta','temp')
      !> created in higher prio loop already 

      case ('sphere')
      !> the sphere constraint si technically identical to the ellipsoid one, but
      !> with equal axis lengths in all 3 directions
         if (kv%na > 0)then
            if(trim(kv%value_rawa(1)) .eq. 'auto')then
               !> determine sphere
               call wallpot_core(mol,rabc, potscal=env%potscal)
               rdum = maxval(rabc(:))
               rabc(:) = rdum
            else            
              read(kv%rawvalue,*,iostat=io) rdum
              if(io == 0) rabc(:) = rdum
            endif
            call get_atlist(mol%nat, atlist, kv%rawvalue, mol%at)
            call cons%deallocate()
            select case( pot ) 
            case ( 1 ) !> polynomial
              call cons%ellipsoid( mol%nat, atlist, rabc, force_constant, alpha, .false.)
            case ( 2 ) !> logfermi
              call cons%ellipsoid( mol%nat, atlist, rabc, T, beta, .true.) 
            end select
            if (debug) call cons%print(stdout)
            call env%calc%add(cons)
         endif

      case('ellipsoid') 
         if (kv%na > 0)then
            if(trim(kv%value_rawa(1)) .eq. 'auto')then
               !> determine ellipsoid
              call wallpot_core(mol,rabc, potscal=env%potscal)
            else
              read (kv%value_rawa(1),*,iostat=io) r1
              if (io == 0) read (kv%value_rawa(2),*,iostat=io) r2
              if (io == 0) read (kv%value_rawa(3),*,iostat=io) r3
              if (io == 0)then
                 rabc(1) = r1
                 rabc(2) = r2
                 rabc(3) = r3
              endif 
            endif
            call get_atlist(mol%nat, atlist, kv%rawvalue, mol%at)
            call cons%deallocate() 
            select case( pot )
            case ( 1 ) !> polynomial
              call cons%ellipsoid( mol%nat, atlist, rabc, force_constant, alpha, .false.)
            case ( 2 ) !> logfermi
              call cons%ellipsoid( mol%nat, atlist, rabc, T, beta, .true.) 
            end select
            if (debug) call cons%print(stdout)
            call env%calc%add(cons)
         endif

      case default
        write (stdout,'(a,a,a)') 'xtb-style input key: "',kv%key,'" not defined for CREST'

      end select
    end do

  end subroutine get_xtb_wall_block

  subroutine get_xtb_fix_block(env,blk)
!**************************************
!* This is a reader for the $fix block
!***************************************
    implicit none
    type(systemdata),intent(inout) :: env
    type(filetype) :: file
    integer :: i,j,k,io
    type(keyvalue) :: kv
    type(datablock),intent(in) :: blk
    real(wp) :: force_constant,dist,angl
    real(wp) :: T,alpha,beta
    real(wp) :: rdum,rabc(3),r1,r2,r3
    type(coord) :: mol
    logical,allocatable :: pairwise(:)
    logical,allocatable :: atlist(:)
    integer :: i1,i2,i3,i4
    integer :: pot


!>--- get reference input geometry
    call env%ref%to(mol)

!>--- get the parameters first
    do i = 1,blk%nkv
      kv = blk%kv_list(i)
      select case (kv%key)

      case ('atoms')
      !> define frozen atoms via indices
        if (.not.allocated(pairwise)) allocate (pairwise(mol%nat),source=.false.)
        call get_atlist(mol%nat, atlist, kv%rawvalue, mol%at)
        do j=1,mol%nat
          if(atlist(j)) pairwise(j) = .true.
        enddo

      case ('elements')
      !> define frozen atoms via elements
        if (.not.allocated(pairwise)) allocate (pairwise(mol%nat),source=.false.)
        if (kv%id == valuetypes%raw_array) then
          do j = 1,kv%na
            i1 = e2i(kv%value_rawa(j))
            do k=1,mol%nat
             if (i1 == mol%at(k)) pairwise(k) = .true.
            enddo
          enddo
        else
          i1 = e2i(kv%rawvalue)
          do j=1,mol%nat
            if (i1 == mol%at(j)) pairwise(j) = .true.
          enddo
        end if

      case default
        write (stdout,'(a,a,a)') 'xtb-style input key: "',kv%key,'" not defined for CREST'

      end select
    end do

    if(allocated(pairwise))then
       i1 = count(pairwise)
       env%calc%nfreeze = i1
       if(debug)then
         write(stdout,'("> ",a)') 'Frozen atoms:'
         do i=1,mol%nat
           if(pairwise(i))write(stdout,'(1x,i0)',advance='no') i
         enddo
         write(stdout,*) 
       endif
       call move_alloc(pairwise,env%calc%freezelist)
    endif

  end subroutine get_xtb_fix_block

!========================================================================================!

  subroutine parse_xtb_input_fallback(fname,dict)
!********************************************************************
!* This is the fallback reader for xtb input files to set up a dict
!********************************************************************
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

!>--- all valid key-values must be in $-blocks, no root-level ones
    get_root_kv = .false.
!>--- the loop where the input file is read
    do i = 1,file%nlines
      if (file%current_line > i) cycle

      !> key-value pairs of the root dict (ignored for xtb)
      if (get_root_kv) then
        call get_keyvalue(kvdum,file%line(i),io)
        if (io == 0) then
          call dict%addkv(kvdum) !> add to dict
        end if
      end if

      !> the $-blocks
      if (isxtbheader(file%line(i))) then
        get_root_kv = .false.
        call read_xtbdatablock(file,i,blkdum)
        call dict%addblk(blkdum) !> add to dict
      end if
    end do

    call file%close()

    return
  end subroutine parse_xtb_input_fallback
!========================================================================================!
  subroutine remove_comments(file)
    use filemod
    implicit none
    type(filetype) :: file
    character(len=1),parameter :: com = '#'
    integer :: i
    do i = 1,file%nlines
      call clearcomment(file%f(i),com)
      call clearcomment(file%f(i),"$end")
    end do
  end subroutine remove_comments

  function isxtbheader(str)
    implicit none
    logical :: isxtbheader
    character(len=*) :: str
    character(len=:),allocatable :: atmp
    integer :: l
    isxtbheader = .false.
    atmp = adjustl(trim(str))
    l = len_trim(atmp)
    if (l < 1) return
    if ((atmp(1:1) == '$')) then
      isxtbheader = .true.
    end if
    return
  end function isxtbheader

!========================================================================================!

  subroutine read_xtbdatablock(file,i,blk)
    implicit none
    type(filetype),intent(inout)  :: file
    type(datablock),intent(inout) :: blk
    integer,intent(in) :: i

    character(len=:),allocatable :: rawline
    type(keyvalue) :: kvdum
    integer :: j,k,io

    call blk%deallocate()

    blk%header = file%line(i)
    call clearxtbheader(blk%header)

    do j = i+1,file%nlines
      rawline = file%line(j)
      if (isxtbheader(rawline)) exit
      call get_xtb_keyvalue(kvdum,rawline,io)
      if (io == 0) then
        call blk%addkv(kvdum)
      end if
    end do

  end subroutine read_xtbdatablock

!========================================================================================!

  subroutine get_xtb_keyvalue(kv,str,io)
    implicit none
    class(keyvalue) :: kv
    character(len=*) :: str
    integer,intent(out) :: io
    character(len=:),allocatable :: tmpstr
    character(len=:),allocatable :: ktmp
    character(len=:),allocatable :: vtmp
    integer :: i,j,k,na,plast
    integer :: l(3)
    call kv%deallocate()
    io = 0
    tmpstr = adjustl(lowercase(str))

    !> key-value conditions
    l(1) = index(tmpstr,'=')
    l(2) = index(tmpstr,':')
    l(3) = index(tmpstr,' ')

    k = 0
    if (l(1) .ne. 0) then
      k = l(1)
    else if (l(2) .ne. 0) then
      k = l(2)
    else if (l(3) .ne. 0) then
      k = l(3)
    end if

    if (k .eq. 0) then
      io = -1
      return
    end if

    ktmp = trim(adjustl(tmpstr(:k-1)))
    vtmp = trim(adjustl(tmpstr(k+1:)))
    kv%key = ktmp !> the key as string
    kv%rawvalue = vtmp !> value as unformatted string

    !> comma denotes an array of strings
    k = index(vtmp,',')
    if (k .ne. 0) then
      kv%id = valuetypes%raw_array
      j = len_trim(vtmp)
      na = 1
      !> count elements
      do i = 1,j
        if (vtmp(i:i) .eq. ',') na = na+1
      end do
      !> allocate
      kv%na = na
      allocate (kv%value_rawa(na),source=repeat(' ',j))
      plast = 1
      na = 1
      do i = 1,j
        if (na == kv%na) then !> for the last argument
          kv%value_rawa(na) = trim(adjustl(vtmp(plast:)))
        end if
        if (vtmp(i:i) .eq. ',') then
          kv%value_rawa(na) = trim(adjustl(vtmp(plast:i-1)))
          plast = i+1
          na = na+1
        end if
      end do

    end if
  end subroutine get_xtb_keyvalue

!========================================================================================!
!> for given input file parse the next block
  subroutine parse_xtbinfile_block(file,i,rawblk)
    implicit none
    type(filetype),intent(inout)      :: file
    type(parseblock),intent(inout) :: rawblk
    integer,intent(in) :: i
    logical :: saveblock
    integer :: j,k,l
    character(len=:),allocatable :: src

    call rawblk%deallocate()

    src = repeat(' ',file%lwidth)

    if (isxtbheader(file%line(i))) then
      saveblock = .true.
      rawblk%header = file%line(i)
      !      cycle
    end if
    !> get blocklength
    k = i+1
    l = 0
    jloop: do j = k,file%nlines
      if (isxtbheader(file%line(j))) then
        file%current_line = j
        exit jloop
      end if
      if (len_trim(file%line(j)) > 0) then
        l = l+1
      end if
      if (j == file%nlines) file%current_line = j
    end do jloop
    !if (l < 1) exit iloop
    if (l < 1) return
    !> get block
    rawblk%len = l
    allocate (rawblk%content(l),source=src)
    l = 0
    jloop2: do j = k,file%nlines
      if (isheader(file%line(j))) then
        saveblock = .false.
        return
      end if
      if (len_trim(file%line(j)) > 0) then
        l = l+1
        rawblk%content(l) = file%line(j)
      end if
    end do jloop2
    !end do iloop

    return
  end subroutine parse_xtbinfile_block

!=======================================================================================!

  subroutine clearxtbheader(hdr)
    implicit none
    character(len=*) :: hdr
    integer :: i,k,l
    character(len=:),allocatable :: atmp,btmp
    character(len=1) :: s
    atmp = adjustl(hdr)
    atmp = trim(atmp)
    !>remove whitespaces
    l = len_trim(atmp)
    k = index(hdr,'$')
    if (k > 0) then
      atmp(k:k) = ' '
      atmp = adjustl(atmp)
    end if
    hdr = trim(atmp)
    return
  end subroutine clearxtbheader

!========================================================================================!
end module parse_xtbinput
