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

module gradreader_module

  use iso_fortran_env,only:wp => real64
  implicit none
  !--- private module variables and parameters
  private

  !> Type enumerator
  type :: enum_filetype
    integer :: unknown = 0
    integer :: engrad = 1
    integer :: orca = 1
    integer :: turbomole = 2
  end type enum_filetype
  type(enum_filetype),parameter :: gradtype = enum_filetype()

  public :: gradtype
  public :: conv2gradfmt
  public :: rd_efile
  public :: rd_grad_engrad
  public :: rd_grad_generic
  public :: rd_grad_tm
  public :: write_engrad

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine write_engrad(fname,energy,grad)
!***************************************************
!* subroutine write_engrad
!* write energy and gradient in the .engrad format
!***************************************************
    implicit none
    character(len=*),intent(in) :: fname
    real(wp),intent(in) :: energy
    real(wp),intent(in) :: grad(:,:)
    integer :: c,ich,n,i,j
    real(wp) :: dum

    n = size(grad,2)
    open (newunit=ich,file=fname)
    write(ich,'(a)') '#'
    write(ich,'(a)') '# Atoms'
    write(ich,'(a)') '#'
    write(ich,'(5x,i0)') n
    write(ich,'(a)') '#'
    write(ich,'(a)') '# Energy ( Eh )'
    write(ich,'(a)') '#'
    write(ich,'(f25.15)') energy
    write(ich,'(a)') '#'
    write(ich,'(a)') '# Gradient ( Eh/a0 )'
    write(ich,'(a)') '#'
    do i=1,n
      do j=1,3
        write(ich,'(f25.15)') grad(j,i)
      enddo
    enddo
    close (ich)
  end subroutine write_engrad

!========================================================================================!

  subroutine rd_efile(fname,energy,iostatus)
!***************************************************
!* subroutine rd_efile
!* read the energy from a file containing only that
!***************************************************
    implicit none
    character(len=*),intent(in) :: fname
    real(wp),intent(out) :: energy
    integer,intent(out) :: iostatus
    integer :: c,iunit,n,i,j
    character(len=128) :: atmp
    real(wp) :: dum
    energy = 0.0_wp
    iostatus = 0
    open (newunit=iunit,file=fname)
    read (iunit,*,iostat=iostatus) energy
    close (iunit)
  end subroutine rd_efile

!========================================================================================!
  subroutine rd_grad_engrad(iunit,nat,energy,grad,iostatus)
!**************************************************
!* subrotuine rd_grad_engrad
!* read *.engrad file (used e.g. by xtb and orca)
!* all comments (#) are ignored
!* first number read should be number of atoms
!* followed by 3N lines á 1 float for the gradient
!**************************************************
    implicit none
    integer,intent(in) :: iunit
    integer,intent(in) :: nat
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: grad(3,nat)
    integer,intent(out) :: iostatus
    integer :: c,io,n,i,j
    character(len=128) :: atmp
    real(wp) :: dum

    iostatus = 0
    energy = 0.0_wp
    grad(:,:) = 0.0_wp

    c = 0
    do
      read (iunit,'(a)',iostat=io) atmp
      if (io < 0) exit
      atmp = adjustl(atmp)
      if (atmp(1:1) == '#') cycle
      if (c == 0) then
        read (atmp,*) n
        if (n /= nat) then
          iostatus = 2
          exit
        end if
        c = c+1
      else if (c == 1) then
        read (atmp,*) energy
        c = c+1
        cycle
      else if (c == 2) then
        backspace (iunit)
        do i = 1,n
          do j = 1,3
            read (iunit,*,iostat=io) dum
            if (io < 0) then
              iostatus = 3
              exit
            end if
            grad(j,i) = dum
          end do
        end do
        c = c+1
      else if (c >= 3) then
        exit
      end if
    end do

    return
  end subroutine rd_grad_engrad

!========================================================================================!
  subroutine rd_grad_tm(iunit,nat,energy,grad,iostatus)
!**************************************************
!* subrotuine rd_grad_tm
!* read a turbomole-type "gradient" file
!**************************************************
    implicit none
    integer,intent(in) :: iunit
    integer,intent(in) :: nat
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: grad(3,nat)
    integer,intent(out) :: iostatus
    integer :: c,io,n,i,j
    character(len=128) :: atmp
    character(len=20) :: btmp(8)
    real(wp) :: dum
    logical :: readblock

    iostatus = 0
    energy = 0.0_wp
    grad(:,:) = 0.0_wp

    c = 0
    readblock = .false.
    do
      read (iunit,'(a)',iostat=io) atmp
      if (io < 0) exit !> EOF exit
      atmp = adjustl(atmp)
      if (atmp(1:4) == '$end') readblock = .false.
      if( readblock ) then      
        if(index(atmp,'cycle').ne.0)then
          read(atmp,*) btmp(1:2),j,btmp(3:6),energy,btmp(7:8),dum
        elseif(c < nat)then !> skip coords
          c = c + 1  
        else !> read grad
          call rd_grad_n3(iunit,nat,grad,iostatus)
          exit
        endif
      endif
      if (atmp(1:5) == '$grad') readblock = .true.
    end do

    return
  end subroutine rd_grad_tm

!========================================================================================!
  subroutine rd_grad_generic(iunit,nat,grad,gradkey,gradfmt,iostatus)
!**************************************************
!* subrotuine rd_grad_generic
!* read unspecified gradient file
!* routine can look for a (case sensitive) keyword 
!* after which the gradient is read
!* NOTE: routine does not provide the energy
!************************************************** 
    implicit none
    integer,intent(in) :: iunit
    integer,intent(in) :: nat
    real(wp),intent(out) :: grad(3,nat)
    character(len=*),intent(in) :: gradkey
    integer,intent(in) :: gradfmt
    integer,intent(out) :: iostatus
    integer :: c,io,n,i,j
    character(len=128) :: atmp
    real(wp) :: dum

    iostatus = 0
    grad(:,:) = 0.0_wp

    c = 0
    do
      read (iunit,'(a)',iostat=io) atmp
      if (io < 0) exit
      atmp = adjustl(atmp)
      if (atmp(1:1) == '#') cycle
      if (index(atmp,gradkey) .ne. 0) then
        c = c+1
        if (gradfmt == 0) then
          call rd_grad_3n(iunit,nat,grad,iostatus)
        else
          call rd_grad_n3(iunit,nat,grad,iostatus)
        end if
      end if
      if (c >= 1) then
        exit
      end if
    end do

  end subroutine rd_grad_generic

!========================================================================================!
  subroutine rd_grad_3n(iunit,nat,grad,iostatus)
!***************************
!* subrotuine rd_grad_3n
!* read 3N lines into grad
!***************************
    implicit none
    integer,intent(in) :: iunit
    integer,intent(in) :: nat
    real(wp),intent(out) :: grad(3,nat)
    integer,intent(out) :: iostatus
    integer :: c,io,n,i,j
    character(len=128) :: atmp
    real(wp) :: dum

    iostatus = 0
    grad(:,:) = 0.0_wp

    c = 0
    do i = 1,n
      do j = 1,3
        read (iunit,*,iostat=io) dum
        if (io < 0) then
          iostatus = 3
          exit
        end if
        grad(j,i) = dum
      end do
    end do

    return
  end subroutine rd_grad_3n
  subroutine rd_grad_n3(iunit,nat,grad,iostatus)
!***************************************
!* subrotuine rd_grad_n3
!* read gradient as N lines á 3 entries
!***************************************
    implicit none
    integer,intent(in) :: iunit
    integer,intent(in) :: nat
    real(wp),intent(out) :: grad(3,nat)
    integer,intent(out) :: iostatus
    integer :: c,io,n,i,j
    character(len=128) :: atmp
    real(wp) :: dum(3)

    iostatus = 0
    grad(:,:) = 0.0_wp

    c = 0
    do i = 1,n
      read (iunit,*,iostat=io) dum(1:3)
      if (io < 0) then
        iostatus = 3
        exit
      end if
      grad(:,i) = dum(:)
    end do

    return
  end subroutine rd_grad_n3

!========================================================================================!
  function conv2gradfmt(str) result(ifmt)
!*************************************************
!* utility function, determine whether to
!* read 3N lines into grad or N lines á 3 entries
!*************************************************
    implicit none
    character(len=*),intent(in) :: str
    integer :: ifmt
    ifmt = 0
    select case (str)
    case ('3n','3N')
      ifmt = 0
    case ('N3','n3','n*3','N*3','n;3','N;3')
      ifmt = 1
    case default
      ifmt = 0
    end select
  end function conv2gradfmt


!========================================================================================!
!========================================================================================!
end module gradreader_module

