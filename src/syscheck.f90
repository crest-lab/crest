!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Philipp Pracht
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

!===========================================================================================!
! Quick and dirty module for checking for a correct setup of a script-like fortran program
!===========================================================================================!
module syscheck
   use iso_fortran_env, wp => real64

   implicit none

   public :: checkprog
   public :: checkprog_secret

   private

!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------

contains
!-------------------------------------------------------------------------------------------------------
! PRIVATE PROCEDURES
!-------------------------------------------------------------------------------------------------------
subroutine getpath(fname,path)
    implicit none
    character(len=*) :: fname
    character(len=*) :: path
    character(len=:),allocatable :: checkcall
    character(len=:),allocatable :: pipe
    integer :: rcode,ich,io


    pipe=' >/dev/null 2>/dev/null'

    checkcall='command -v '//trim(fname)//pipe
    call execute_command_line(checkcall,.true.,rcode)

    if(rcode.eq.0)then
      checkcall='command -v '//trim(fname)//' > pathout.tmp 2>/dev/null' 
      !call system(checkcall)
      call execute_command_line(trim(checkcall), exitstat=io)
      open(newunit=ich,file='pathout.tmp')
      read(ich,'(a)')path
      close(ich,status='delete')
    endif
    return
end subroutine getpath

!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------
! PUBLIC subroutines from the module
!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------
subroutine checkprog(fname,r)
    implicit none
    character(len=*) :: fname
    character(len=:),allocatable :: checkcall
    character(len=:),allocatable :: pipe
    integer :: rcode,r
    character(len=512) :: path
    

    pipe=' >/dev/null 2>/dev/null'

    checkcall='command -v '//trim(fname)//pipe
    call execute_command_line(checkcall,.true.,rcode)     

    write(*,'(4x,a,a,a)')'binary: "',trim(fname),'"'
    if(rcode.ne.0)then
      write(*,'(4x,a)')'status: not found'
      r=r+1
    else
      write(*,'(4x,a)')'status: present'
      call getpath(fname,path)
      write(*,'(4x,a,a)')'path  : ',trim(path)
    endif

    return
end subroutine checkprog
!-------------------------------------------------------------------------------------------------------
! the same as above but only provide a printout (at all) if the program is not present!
subroutine checkprog_secret(fname) 
    implicit none
    character(len=*) :: fname
    character(len=:),allocatable :: checkcall
    character(len=:),allocatable :: pipe
    integer :: rcode

    pipe=' >/dev/null 2>/dev/null'

    checkcall='command -v '//trim(fname)//pipe
    call execute_command_line(checkcall,.true.,rcode)

    if(rcode.ne.0)then
      write(*,'(4x,a,a,a)')'binary: "',trim(fname),'"'
      write(*,'(4x,a)')'status: not found'
    endif

    return
end subroutine checkprog_secret
!-------------------------------------------------------------------------------------------------------
end module syscheck
