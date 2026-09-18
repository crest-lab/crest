!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023-2024 Philipp Pracht
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

subroutine creststop(io)
  use crest_parameters
  use crest_data
  implicit none
  integer,intent(in) :: io

  call graceful_shutdowns()

  select case (io)
  case (status_normal)
    write (stdout,*) 'CREST terminated normally.'
  case default
    write (stdout,*) 'CREST terminated abnormally.'
  case (status_error)
    write (stdout,*) 'CREST terminated with errors.'
  case (status_ioerr)
    write (stdout,*) 'CREST terminated with I/O errors.'
  case (status_args)
    write (stdout,*) 'CREST terminated due to invalid parameters.'
  case (status_input)
    write (stdout,*) 'CREST terminated due to failed input file read.'
  case (status_config)
    write (stdout,*) 'CREST terminated due to invalid configuration.'
  case (status_failed)
    write (stdout,*) 'CREST terminated with failures.'
  case (status_safety)
    write (stdout,*) 'Safety termination of CREST.'
  end select
  call exit(io)

end subroutine creststop

!================================================================================!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!================================================================================!
!&>
#if defined(__INTEL_LLVM_COMPILER)
  subroutine wsigint() bind(C,name="crest_wsigint") !> Ctrl+C
#else
  subroutine wsigint() !> Ctrl+C
#endif
    use crest_parameters,only:stderr,stdout
    write (*,*)
    write (stderr,'(" recieved SIGINT, trying to terminate CREST...")')
    !call dump_restart()
    call exit(130)
    error stop
  end subroutine wsigint

#if defined(__INTEL_LLVM_COMPILER)
  subroutine wsigquit() bind(C,name="crest_wsigquit") !> Ctrl+D
#else
  subroutine wsigquit() !> Ctrl+D
#endif
    use crest_parameters,only:stderr,stdout
    write (*,*)
    write (stderr,'(" recieved SIGQUIT, trying to terminate CREST...")')
    !call dump_restart()
    call exit(131)
    error stop
  end subroutine wsigquit

#if defined(__INTEL_LLVM_COMPILER)
  subroutine wsigterm() bind(C,name="crest_wsigterm") !> Recieved by the "kill" pid command
#else
  subroutine wsigterm() !> Recieved by the "kill" pid command
#endif
    use crest_parameters,only:stderr,stdout
    write (stdout,*)
    write (stderr,'(" recieved SIGTERM, trying to terminate CREST...")')
    !call dump_restart()
    call exit(143)
    error stop
  end subroutine wsigterm

#if defined(__INTEL_LLVM_COMPILER)
  subroutine wsigkill() bind(C,name="crest_wsigkill")
#else
  subroutine wsigkill()
#endif
    use crest_parameters,only:stderr,stdout
    !call dump_restart()
    call exit(137)
    error stop 'CREST recieved SIGKILL.'
  end subroutine wsigkill

  subroutine initsignal()
#if defined(__INTEL_LLVM_COMPILER)
   ! ifx: libifport's SIGNAL intrinsic crashes with ifx procedure thunks.
   ! Register handlers via ISO_C_BINDING → crest_install_signal() in signal.c.
   use iso_c_binding,only:c_int,c_funloc,c_funptr
   implicit none
   interface
     subroutine crest_install_signal(signum,handler) &
       bind(C,name='crest_install_signal')
       import :: c_int,c_funptr
       integer(c_int),value :: signum
       type(c_funptr),value :: handler
     end subroutine
     subroutine wsigint() bind(C,name='crest_wsigint')
     end subroutine
     subroutine wsigquit() bind(C,name='crest_wsigquit')
     end subroutine
     subroutine wsigterm() bind(C,name='crest_wsigterm')
     end subroutine
   end interface
   call crest_install_signal(2_c_int,c_funloc(wsigint))
   call crest_install_signal(3_c_int,c_funloc(wsigquit))
   call crest_install_signal(15_c_int,c_funloc(wsigterm))
   ! SIGKILL (9) cannot be caught; signal 69 is invalid — omit both.
#else
   external :: wSIGINT
   external :: wSIGTERM
   external :: wSIGKILL
   external :: wSIGQUIT
   call signal(2,wSIGINT)
   call signal(3,wSIGQUIT)
   call signal(9,wSIGKILL)
   call signal(15,wSIGTERM)
   call signal(69,wSIGINT)
#endif
  end subroutine initsignal

!=============================================================!
  subroutine graceful_shutdowns()
    use mlip_sc
    implicit none
    call mlips_shutdown()
  end subroutine graceful_shutdowns
!&<

