subroutine wsigint
   use iso_fortran_env, only : error_unit
   write(error_unit,'(" recieved SIGINT, terminating...")')
   error stop
end subroutine wsigint

subroutine wsigterm
   use iso_fortran_env, only : error_unit
   write(error_unit,'(" recieved SIGTERM, terminating...")')
   error stop
end subroutine wsigterm

subroutine initsignal()
    external :: wSIGINT
    external :: wSIGTERM

    call signal(2,wSIGINT)
    call signal(15,wSIGTERM)
end subroutine initsignal
