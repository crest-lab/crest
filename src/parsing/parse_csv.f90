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

!> This module implements a reader for csv files

module parse_csv
  use filemod
  use iomod
  implicit none
  private
  !logical,parameter,private :: debug = .true.
  logical,parameter,private :: debug = .false.

  interface parse_csv_file_column
    module procedure :: parse_csv_file_columnname
    module procedure :: parse_csv_file_columnnumber
  end interface parse_csv_file_column

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine parse_csv_file_columnname(fname,header,column)
!*********************************************
!* Routine for parsing the csv file fname
!* and get a column as array of strings
!*********************************************
    implicit none
    character(len=*),intent(in) :: fname
    character(len=*),intent(in) :: header
    character(len=:),intent(out),allocatable :: column(:)
    logical :: ex
    character(len=:),allocatable :: hdr
    integer :: i,j,k,l
    integer :: nrow,ncol
    type(filetype) :: file

    inquire (file=fname,exist=ex)
    if (.not.ex) return
    call file%open(fname)

    call csv_params(file,nrow,ncol)
    

    if (debug) stop
  end subroutine parse_csv_file_columnname

  subroutine parse_csv_file_columnnumber(fname)
!*********************************************
!* Routine for parsing the csv file fname
!* and get a column as array of strings
!*********************************************
    implicit none
    character(len=*) :: fname

    logical :: ex
    character(len=:),allocatable :: hdr
    integer :: i,j,k,l
    type(filetype) :: file


    inquire (file=fname,exist=ex)
    if (.not.ex) return


    if (debug) stop
  end subroutine parse_csv_file_columnumber

!========================================================================================!

  function count_columns(str) result(commas)
!***************************************************
!* count columns of a row based on number of commas
!***************************************************
    implicit none
    character(len=*) :: str
    character(len=:),allocatable :: atmp
    integer :: l,commas
    atmp = adjustl(trim(str))
    l = len_trim(atmp)
    commas = 0
    if (l < 1) return
    if ((atmp(1:1) == ',')) then
      commas = commas + 1 
    end if
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
    integer :: l,k 
    atmp = adjustl(trim(str))
    l = len_trim(atmp)
    col = ''
    k = 0
    if (l < 1) return
    if ((atmp(1:1) == ',')) then
      commas = commas + 1
    end if
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
    cloumns = 0
    rows = file%nlines
    do i=1,rows
      c = count_columns( file%line(i) )
      if( c > columns ) columns = c  
    enddo 
  end subroutine csv_params 

!========================================================================================!
!========================================================================================!
end module parse_csv
