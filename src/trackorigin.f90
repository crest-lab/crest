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

!=====================================================================================!
! Add the string "origin" as a comment in every energy-line of an ensemble file
!=====================================================================================!
subroutine addorigin(filename,origin)
  use iso_fortran_env,only:wp => real64
  use iomod
  use strucrd
  implicit none
  character(len=*),intent(in) :: filename
  character(len=*),intent(in) :: origin
  integer :: i
  integer :: nat,nall
  integer,allocatable :: at(:)
  real(wp),allocatable :: xyz(:,:,:)
  real(wp),allocatable :: er(:)
  character(len=128),allocatable :: comments(:)

  call rdensembleparam(filename,nat,nall)
  if (nall < 1) then
    write (*,*) 'Warning: file ',trim(filename),' is empty!'
    return
  end if
  allocate (at(nat),comments(nall),xyz(3,nat,nall),er(nall))
  call rdensemble(filename,nat,nall,at,xyz,comments)

  do i = 1,nall
    er(i) = grepenergy(comments(i))
    comments(i) = '   !'//trim(origin)
  end do
  call wrensemble('dummy',nat,nall,at,xyz,er,comments)

  call remove(filename)
  call rename('dummy',filename)

  return
end subroutine addorigin

!--------------------------------------------------------------------------------------
! Add the string "origin" as a comment in every energy-line of an ensemble file
!--------------------------------------------------------------------------------------
subroutine addorigin2(filename,origin,nall)
  use iomod
  implicit none
  integer,intent(in) :: nall
  character(len=*),intent(in) :: filename
  character(len=*),intent(in) :: origin(nall)

  integer :: i,j,k
  integer :: io,from,to
  character(len=512) :: atmp

  open (newunit=from,file=filename)
  open (newunit=to,file='dummy')

  do k = 1,nall
    read (from,'(a)',iostat=io) atmp
    if (io < 0) exit
    read (atmp,*) j
    write (to,'(a)') trim(atmp)
    read (from,'(a)') atmp
    write (to,'(a,3x,a,a)') trim(atmp),'!',trim(origin(k))
    do i = 1,j
      read (from,'(a)',iostat=io) atmp
      if (io < 0) exit
      write (to,'(a)') trim(atmp)
    end do
  end do

  close (from)
  close (to)
  call remove(filename)
  call rename('dummy',filename)

end subroutine addorigin2

!====================================================================================!
! Add a timestamp indicator  as a comment in every energy-line of an ensemble file
! The whole length(=nall) of the trajectory is split in "nsplit" parts.
!====================================================================================!
subroutine addtimestamp(filename,nsplit)
  use iso_fortran_env,wp => real64
  use iomod
  use strucrd,only:rdensembleparam,rdensemble
  implicit none
  character(len=*),intent(in) :: filename
  integer,intent(in) :: nsplit

  integer :: i
  integer :: nat,nall
  integer :: tag3
  real(wp) :: tag1,tag2
  character(len=40),allocatable :: stamp(:)

  call rdensembleparam(filename,nat,nall)
  allocate (stamp(nall))
  do i = 1,nall
    tag1 = float(i)/float(nall)
    tag2 = tag1*float(nsplit)
    tag3 = nint(tag2)
    tag3 = max(1,tag3) !no 0 timestamp
    write (stamp(i),'(a,i0)') 't',tag3
  end do
  call addorigin2(filename,stamp,nall)
  deallocate (stamp)
  return
end subroutine addtimestamp

!--------------------------------------------------------------------------------------
! same as above, but only for a single xyz file
!--------------------------------------------------------------------------------------
subroutine addoriginXYZ(filename,origin)
  use iomod
  implicit none
  character(len=*),intent(in) :: filename
  character(len=*) :: origin
  integer :: i,j
  integer :: io,from,to
  character(len=512) :: atmp

  open (newunit=from,file=filename)
  open (newunit=to,file='dummy')

  read (from,'(a)',iostat=io) atmp
  read (atmp,*) j
  write (to,'(a)') trim(atmp)
  read (from,'(a)') atmp
  write (to,'(a,3x,a,a)') trim(atmp),'!',trim(origin)
  do i = 1,j
    read (from,'(a)',iostat=io) atmp
    if (io < 0) exit
    write (to,'(a)') trim(atmp)
  end do

  close (from)
  close (to)
  call remove(filename)
  call rename('dummy',filename)

end subroutine addoriginXYZ

!==================================================================================!
! Set origin as comment for several trajectories in subdirectories
!==================================================================================!
subroutine set_trj_origins(base,origin)
  use iomod
  implicit none

  character(len=*) :: base
  character(len=*) :: origin
  character(len=256) :: str,dir
  integer :: i
  logical :: ex
  i = 1
  do
    write (dir,'(a,i0)') trim(base),i
    ex = directory_exist(trim(dir))
    if (.not.ex) then
      exit
    else
      write (str,'(a,a,''xtb.trj'')') trim(dir),'/'
      write (dir,'(a,i0)') trim(origin),i
      call addorigin(trim(str),trim(dir))
      i = i+1
    end if
  end do
  return
end subroutine set_trj_origins

!==================================================================================!
! Set origin as comment for several trajectories in subdirectories
! each trajectory length is cut in "nsplit" parts
!==================================================================================!
subroutine set_trj_timestamps(base,nsplit)
  use iomod
  implicit none

  character(len=*) :: base
  integer,intent(in) :: nsplit
  character(len=256) :: str,dir
  integer :: i
  logical :: ex
  i = 1
  do
    write (dir,'(a,i0)') trim(base),i
    ex = directory_exist(trim(dir))
    if (.not.ex) then
      exit
    else
      write (str,'(a,a,''xtb.trj'')') trim(dir),'/'
      call addtimestamp(trim(str),nsplit)
      i = i+1
    end if
  end do
  return
end subroutine set_trj_timestamps

!==================================================================================!
! get the string "origin" as a substring from line, if seperated by "!"
!==================================================================================!
subroutine getorigin(line,origin)
  implicit none
  character(len=*),intent(in) :: line
  character(len=*),intent(inout) :: origin
  integer :: i,j
  character(len=1) :: digit

  origin = ''

  do i = 1,len(trim(line))
    digit = line(i:i)
    if (digit == '!') then
      j = i+1
      origin = trim(line(j:))
      exit
    end if
  end do
  return
end subroutine getorigin

!==========================================================================!
! subrotuine rdensemble_origin
! is a wrappen for the rdensemble routine from module strucrd
! that extracts the energy and the origin tag from the comment lines
!==========================================================================!
subroutine rdensemble_origin(fname,nat,nall,at,xyz,energy,origin)
  use iso_fortran_env,wp => real64
  use strucrd,only:rdensemble,grepenergy
  implicit none
  character(len=*) :: fname
  integer,intent(inout)  :: nat
  integer,intent(inout)  :: nall
  integer,intent(inout)  :: at(nat)
  real(wp),intent(inout) :: xyz(3,nat,nall)
  real(wp),intent(inout) :: energy(nall)
  character(len=*),intent(inout) :: origin(nall)

  character(len=256),allocatable :: comments(:)
  character(len=256) :: a
  integer :: i

  allocate (comments(nall))
  call rdensemble(fname,nat,nall,at,xyz,comments)

  do i = 1,nall
    a = comments(i)
    energy(i) = grepenergy(a)
    call getorigin(trim(a),origin(i))
  end do

  deallocate (comments)
  return
end subroutine rdensemble_origin

!=========================================================================!
!  get timestamp tag (as an integer) from a string tag
!=========================================================================!
subroutine origin2time(nall,orig,tims)
  implicit none
  integer :: nall
  character(len=*) :: orig(nall)
  integer :: tims(nall)
  character(len=40) :: atmp
  integer :: i
  do i = 1,nall
    atmp = trim(orig(i))
    if (atmp(1:1) == 't') then !should work only with timetag flag
      atmp = atmp(2:) !cut the "t"
      read (atmp,*) tims(i)
    else
      tims(i) = 1
    end if
  end do
  return
end subroutine origin2time
