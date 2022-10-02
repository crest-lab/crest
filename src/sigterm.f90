subroutine wsigint
  use iso_fortran_env,only:error_unit
  write (error_unit,'(" recieved SIGINT, terminating...")')
  call exit(1)
  error stop
end subroutine wsigint

subroutine wsigterm
  use iso_fortran_env,only:error_unit
  write (error_unit,'(" recieved SIGTERM, terminating...")')
  call exit(1)
  error stop
end subroutine wsigterm

subroutine initsignal()
  external :: wSIGINT
  external :: wSIGTERM
!  external :: exit

!  call signal(2,exit)
  call signal(2,wSIGINT)
  call signal(15,wSIGTERM)
  call signal(69,wSIGINT)
end subroutine initsignal

integer(4) function h_abort_sigint(signum)
  integer(4) :: signum
  !DIRS$ ATTRIBUTES DEFAULT :: h_abort_sigint
  write (*,*) 'SIGINT recieved, terminating program ...'
  h_abort_sigint = 1
  error stop
end function h_abort_sigint

integer(4) function h_abort_sigterm(signum)
  integer(4) :: signum
  !DIRS$ ATTRIBUTES DEFAULT :: h_abort_sigterm
  write (*,*) 'SIGTERM recieved, terminating program ...'
  h_abort_sigterm = 1
  error stop
end function h_abort_sigterm

!subroutine initsignal2(io)
!        use ifport, only : ifsignal=>signal,SIGINT,SIGTERM
!  external :: h_abort_sigint
!  external :: h_abort_sigterm
!  integer(4) :: h_abort_sigint
!  integer(4) :: h_abort_sigterm
!  integer,intent(out) :: io
!  integer(4) :: io4
!  io = 0
!
!  io4 = ifsignal(SIGINT,h_abort_sigint,-1)
!  if (io4 == 1) io = 1
!  io4 = ifsingal(SIGTERM,h_abort_sigterm,-1)
!  if (io4 == 1) io = 1
!
!  return
!end subroutine initsignal2
