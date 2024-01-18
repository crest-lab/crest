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

!> scratch directory handling (for simplicity just use shell commands)

!=========================================================================================!
subroutine scrdir(env)
  use crest_parameters
  use crest_data
  use iomod
  type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
  integer :: ich,io

  if (len_trim(env%scratchdir) .lt. 1) then
    call command('mktemp -d > tmpconf 2>/dev/null',io)
    open (newunit=ich,file='tmpconf')
    read (ich,'(a)',iostat=io) env%scratchdir
    if (io < 0) then   ! if mktemp failed and tmpconf is empty
      env%scratchdir = ''
      error stop 'Failed to create scratch directory!'
      return
    end if
    close (ich,status='delete')
  end if

  write (stdout,'(a,a)') 'Home directory          : ',trim(env%homedir)
  write (stdout,'(a,a)') 'Using scratch directory : ',trim((env%scratchdir))
  io = makedir(trim(env%scratchdir))

  call copy('.CHRG',trim(env%scratchdir)//'/'//'.CHRG')
  call copy('.UHF',trim(env%scratchdir)//'/'//'.UHF')

  write (stdout,'(a)',advance='no') 'Copying data to scratch directory ...'
  flush (stdout)
  call command('scp -r ./* '//trim(env%scratchdir)//'/')
  write (stdout,*) 'done.'

  call chdir(trim(env%scratchdir))

end subroutine scrdir

subroutine scrend(env)
  use crest_parameters
  use crest_data
  use iomod
  type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
  character(len=1024) :: crefi,crefi2
  logical :: ex

  if (len_trim(env%scratchdir) .lt. 1) then
    return
  end if

  write (stdout,'(/,a)',advance='no') 'Retrieving data from scratch directory ...'
  flush (stdout)
  call command('scp -r '//trim(env%scratchdir)//'/* '//trim(env%homedir)//'/')
  write (stdout,*) 'done.'

  if (.not.env%keepScratch) then
    write (stdout,'(a)',advance='no') 'Removing scratch directory ...'
    flush (stdout)
    call rmrf(env%scratchdir)
    write (stdout,*) 'done.'
  end if

  return
end subroutine scrend
