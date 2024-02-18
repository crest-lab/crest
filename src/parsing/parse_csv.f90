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

!> This module implements a simple reader for csv files
!> NOT suitable for large csv tables

module parse_csv
  use filemod
  use iomod
  implicit none
  private
  !logical,parameter,private :: debug = .true.
  logical,parameter,private :: debug = .false.

  public :: parse_csv_file_column
  interface parse_csv_file_column
    module procedure :: parse_csv_file_columnname
    module procedure :: parse_csv_file_columnnumber
  end interface parse_csv_file_column

  public :: parse_csv_column_int
  interface parse_csv_column_int
    module procedure :: parse_csv_columnname_int
    module procedure :: parse_csv_columnnumber_int
  end interface parse_csv_column_int

  public :: parse_csv_column_real
  interface parse_csv_column_real
    module procedure :: parse_csv_columnname_real
    module procedure :: parse_csv_columnnumber_real
  end interface parse_csv_column_real

  public :: parse_csv_allcolumns
  interface parse_csv_allcolumns
    module procedure :: parse_csv_allcolumns_real
  end interface parse_csv_allcolumns

  public :: parse_csv_file_row
  interface parse_csv_file_row
    module procedure :: parse_csv_file_rownumber
  end interface parse_csv_file_row

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine parse_csv_file_columnname(fname,header,column)
!******************************************
!* Routine for parsing the csv file fname
!* and get a column as array of strings
!******************************************
    implicit none
    character(len=*),intent(in) :: fname
    character(len=*),intent(in) :: header
    character(len=:),intent(out),allocatable :: column(:)
    logical :: ex
    character(len=:),allocatable :: hdr
    integer :: i,j,k,l
    integer :: nrow,ncol,getcol
    type(filetype) :: file

    inquire (file=fname,exist=ex)
    if (.not.ex) return
    call file%open(fname)

    call csv_params(file,nrow,ncol)
    if (debug) write (*,*) 'nrow',nrow
    if (debug) write (*,*) 'ncol',ncol
    l = file%lwidth
    allocate (column(nrow),source=repeat(' ',l))

    !> first line should contain the header names
    hdr = file%line(1)
    if (debug) write (*,*) trim(hdr)
    getcol = csv_get_column_number(hdr,header,ncol)
    if (debug) write (*,*) 'trying to get column',getcol
    if (getcol > 0) then
      do i = 1,nrow
        column(i) = get_column_element(file%line(i),getcol)
        if (debug) write (*,*) trim(column(i))
      end do
    end if

    call file%close()
    !if (debug) stop
  end subroutine parse_csv_file_columnname

  subroutine parse_csv_file_columnnumber(fname,getcol,column)
!*********************************************
!* Routine for parsing the csv file fname
!* and get a column as array of strings
!*********************************************
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: getcol
    character(len=:),intent(out),allocatable :: column(:)
    logical :: ex
    character(len=:),allocatable :: hdr
    integer :: i,j,k,l,nrow,ncol
    type(filetype) :: file

    inquire (file=fname,exist=ex)
    if (.not.ex) return
    call file%open(fname)

    call csv_params(file,nrow,ncol)
    if (debug) write (*,*) 'nrow',nrow
    if (debug) write (*,*) 'ncol',ncol
    l = file%lwidth
    allocate (column(nrow),source=repeat(' ',l))

    if (debug) write (*,*) 'trying to get column',getcol
    if (getcol > 0.and.getcol <= ncol) then
      do i = 1,nrow
        !if (debug) write(*,*) file%line(i)
        column(i) = get_column_element(file%line(i),getcol)
        if (debug) write (*,*) trim(column(i))
      end do
    end if
    call file%close()
    !if (debug) stop
  end subroutine parse_csv_file_columnnumber


  subroutine parse_csv_file_rownumber(fname,getrow,row)
!*********************************************
!* Routine for parsing the csv file fname
!* and get a column as array of strings
!*********************************************
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: getrow
    character(len=:),intent(out),allocatable :: row(:)
    logical :: ex
    character(len=:),allocatable :: hdr
    integer :: i,j,k,l,nrow,ncol
    type(filetype) :: file

    inquire (file=fname,exist=ex)
    if (.not.ex) return
    call file%open(fname)

    call csv_params(file,nrow,ncol)
    if (debug) write (*,*) 'nrow',nrow
    if (debug) write (*,*) 'ncol',ncol
    l = file%lwidth
    allocate (row(nrow),source=repeat(' ',l))

    if (debug) write (*,*) 'trying to get row elements',getrow
    if (getrow > 0.and.getrow <= nrow) then
      write(*,*) ncol,l
      do i = 1,ncol
        write(*,*) trim(get_column_element(file%line(getrow),i) )
        !row(i) = get_column_element(file%line(getrow),i)
        if (debug) write (*,*) trim(row(i))
      end do
    end if
    call file%close()
    !if (debug) stop
 end subroutine parse_csv_file_rownumber

!========================================================================================!

  subroutine parse_csv_columnname_real(fname,header,column)
!*********************************************
!* Routine for parsing the csv file fname
!* and get a column as array of reals
!*********************************************
    implicit none
    character(len=*),intent(in) :: fname
    character(len=*),intent(in) ::header
    real(wp),intent(out),allocatable :: column(:)
    character(len=:),allocatable :: strcolumn(:)
    logical :: ex
    character(len=:),allocatable :: hdr
    integer :: i,j,k,l,io
    real(wp) :: dum

    if(debug) write(*,*) 'parse column as real'
    call parse_csv_file_column(fname,header,strcolumn)
    l = size(strcolumn,1)
    k = l-1
    allocate (column(k),source=0.0_wp)
    do i = 2,l
      k = i-1
      read (strcolumn(i),*,iostat=io) dum
      if (io == 0) column(k) = dum
      if (debug) write (*,*) dum
    end do
    deallocate (strcolumn)
    if (debug) stop
  end subroutine parse_csv_columnname_real

  subroutine parse_csv_columnnumber_real(fname,getcol,column)
!*********************************************
!* Routine for parsing the csv file fname
!* and get a column as array of reals
!*********************************************
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: getcol
    real(wp),intent(out),allocatable :: column(:)
    character(len=:),allocatable :: strcolumn(:)
    logical :: ex
    integer :: i,j,k,l,io
    real(wp) :: dum

    if(debug) write(*,*) 'parse column as real'
    call parse_csv_file_column(fname,getcol,strcolumn)
    l = size(strcolumn,1)
    k = l-1
    allocate (column(k),source=0.0_wp)
    do i = 2,l
      k = i-1
      read (strcolumn(i),*,iostat=io) dum
      if (io == 0) column(k) = dum
      if (debug) write (*,*) dum
    end do
    deallocate (strcolumn)
    if (debug) stop
  end subroutine parse_csv_columnnumber_real


  subroutine parse_csv_allcolumns_real(fname,columns,cols,rows)
!*********************************************
!* Routine for parsing the csv file fname
!* and get a matrix of all columns/rows
!*********************************************
    implicit none
    character(len=*),intent(in) :: fname
    real(wp),intent(out),allocatable :: columns(:,:)
    integer,intent(out),optional :: cols,rows
    character(len=:),allocatable :: strcolumn
    type(filetype) :: file
    logical :: ex
    integer :: getcol,i,j,k,l,io,nrow,ncol
    real(wp) :: dum

    inquire (file=fname,exist=ex)
    if (.not.ex) return
    call file%open(fname)

    call csv_params(file,nrow,ncol)
    if (debug) write (*,*) 'nrow',nrow
    if (debug) write (*,*) 'ncol',ncol
    l = file%lwidth
    strcolumn=repeat(' ',l)
    allocate (columns(ncol,nrow-1), source=0.0_wp)

    if (debug) write (*,*) 'trying to get column',getcol
    do i = 2,nrow
      k = i-1 !> to skip header
      do getcol=1,ncol
        strcolumn = get_column_element(file%line(i),getcol)
        if (debug) write (*,*) trim(strcolumn)
        read (strcolumn,*,iostat=io) dum
        if (io == 0) columns(getcol,k) = dum
        if (debug) write (*,*) dum
           
      end do
    enddo
    deallocate(strcolumn)

    if(present(cols))then
       cols = ncol
    endif
    if(present(rows))then
       rows = nrow-1  !> not counting the header lines
    endif

    call file%close()
  end subroutine parse_csv_allcolumns_real

!========================================================================================!

  subroutine parse_csv_columnname_int(fname,header,column)
!*********************************************
!* Routine for parsing the csv file fname
!* and get a column as array of integers
!*********************************************
    implicit none
    character(len=*),intent(in) :: fname
    character(len=*),intent(in) ::header
    integer,intent(out),allocatable :: column(:)
    character(len=:),allocatable :: strcolumn(:)
    logical :: ex
    character(len=:),allocatable :: hdr
    integer :: i,j,k,l,io

    call parse_csv_file_column(fname,header,strcolumn)
    l = size(strcolumn,1)
    k = l-1
    allocate (column(k),source=0)
    do i = 2,l
      k = i-1
      read (strcolumn(i),*,iostat=io) j
      if(io==0) column(k) = j
      if (debug) write (*,*) j
    end do
    deallocate (strcolumn)
    if (debug) stop
  end subroutine parse_csv_columnname_int

  subroutine parse_csv_columnnumber_int(fname,getcol,column)
!*********************************************
!* Routine for parsing the csv file fname
!* and get a column as array of integers
!*********************************************
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: getcol
    integer,intent(out),allocatable :: column(:)
    character(len=:),allocatable :: strcolumn(:)
    logical :: ex
    integer :: i,j,k,l,io

    call parse_csv_file_column(fname,getcol,strcolumn)
    l = size(strcolumn,1)
    k = l-1
    allocate (column(k),source=0)
    do i = 2,l
      k = i-1
      read (strcolumn(i),*,iostat=io) j
      if(io==0) column(k) = j
      if (debug) write (*,*) j
    end do
    deallocate (strcolumn)
    !if (debug) stop
  end subroutine parse_csv_columnnumber_int

!========================================================================================!

  function count_columns(str) result(columns)
!***************************************************
!* count columns of a row based on number of commas
!***************************************************
    implicit none
    character(len=*) :: str
    character(len=:),allocatable :: atmp
    integer :: columns
    integer :: l,commas,i
    atmp = adjustl(trim(str))
    l = len_trim(atmp)
    commas = 0
    columns = 1
    if (l < 1) return
    do i = 1,l
      if ((atmp(i:i) == ',')) then
        commas = commas+1
      end if
    end do
    columns = columns+commas
    return
  end function count_columns

!========================================================================================!

  function get_column_element(str,nr) result(col)
!***************************************************
!* count columns of a row based on number of commas
!***************************************************
    implicit none
    character(len=*) :: str
    integer,intent(in) :: nr
    character(len=:),allocatable :: atmp
    character(len=:),allocatable :: col
    integer :: l,k,commas
    atmp = adjustl(trim(str))
    l = len_trim(atmp)
    col = ''
    k = 0
    commas = 0
    if (l < 1) return
    do k = 1,l
      if ((atmp(k:k) == ',')) then
        if (commas+1 == nr) exit
        commas = commas+1
        col = ''
      end if
      col = col//atmp(k:k)
      !if(debug) write(*,*) col
    end do
    !> if we don't have enough columnd
    if (nr-commas > 1) col = ''
    !> some formatting
    if (col(1:1) == ',') col(1:1) = ' '
    col = trim(adjustl(col))
    return
  end function get_column_element

!========================================================================================!

  subroutine csv_params(file,rows,columns)
!*****************************************************************
!* analyse a csv file and get the max numbers of rows and columns
!*****************************************************************
    implicit none
    type(filetype) :: file
    integer,intent(out) :: rows,columns
    integer :: i,c
    rows = 0
    columns = 0
    rows = file%nlines
    do i = 1,rows
      c = count_columns(file%line(i))
      if (c > columns) columns = c
    end do
  end subroutine csv_params

!========================================================================================!

  function csv_get_column_number(line,header,maxelement) result(colnr)
!***************************************************
!* get the column number matching to a given header
!***************************************************
    implicit none
    integer :: colnr
    character(len=*),intent(in) :: line
    character(len=*),intent(in) :: header
    integer,intent(in) :: maxelement
    character(len=:),allocatable :: element
    integer :: i,j,k
    logical :: casesensitive
    character(len=26),parameter :: abc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    colnr = 0
    !> check if we have any uppercase characters in the header.
    !> if so, we treat the search case-sensitively
    casesensitive = .false.
    do i = 1,26
      if (index(header,abc(i:i)) .ne. 0) casesensitive = .true.
    end do
    !> iterate through column elements
    do i = 1,maxelement
      if (casesensitive) then
        element = get_column_element(line,i)
        element = lowercase(element)
        if (element .eq. lowercase(header)) then
          colnr = i
          exit
        end if
      else
        element = get_column_element(line,i)
        if (element .eq. header) then
          colnr = i
          exit
        end if
      end if
    end do
  end function csv_get_column_number

!========================================================================================!
!========================================================================================!
end module parse_csv
