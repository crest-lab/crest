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
subroutine wsigint !> Ctrl+C
  use crest_parameters,only:stderr,stdout
  use crest_restartlog,only:dump_restart
  use ConfSolv_module
  integer :: myunit,io
  write (*,*)
  write (stderr,'(" recieved SIGINT, trying to terminate CREST...")')
  !call dump_restart()
  call cs_shutdown(io)
  call exit(1)
  error stop
end subroutine wsigint

subroutine wsigquit !> Ctrl+D or Ctrl+\
  use crest_parameters,only:stderr,stdout
  use crest_restartlog,only:dump_restart
  use ConfSolv_module
  integer :: myunit,io
  write (*,*)
  write (stderr,'(" recieved SIGQUIT, trying to terminate CREST...")')
  !call dump_restart()
  call cs_shutdown(io)
  call exit(1)
  error stop
end subroutine wsigquit

subroutine wsigterm !> Recieved by the "kill" pid command
  use crest_parameters,only:stderr,stdout
  use crest_restartlog,only:dump_restart
  use ConfSolv_module
  integer :: io
  write (stdout,*)
  write (stderr,'(" recieved SIGTERM, trying to terminate CREST...")')
  !call dump_restart()
  call cs_shutdown(io)
  call exit(1)
  error stop
end subroutine wsigterm

subroutine wsigkill
  use crest_parameters,only:stderr,stdout
  use crest_restartlog,only:dump_restart
  use ConfSolv_module
  integer :: io
  !call dump_restart()
  call cs_shutdown(io)
  error stop 'CREST recieved SIGKILL.'
end subroutine wsigkill

subroutine initsignal()
  external :: wSIGINT
  external :: wSIGTERM
  external :: wSIGKILL
  external :: wSIGQUIT

  call signal(2,wSIGINT)
  call signal(3,wSIGQUIT)
  call signal(9,wSIGKILL)
  call signal(15,wSIGTERM)
  call signal(69,wSIGINT)
end subroutine initsignal

