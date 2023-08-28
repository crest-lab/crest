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

!> A data storage module for hosting ConfSolv as a http server with requests
!> We need to know the server's PORT to create the request
!> CREST can either try to host the server on its own, or the user
!> must provide the PORT.
!> All of this assumes that the ConfSolv submodule was loaded

!========================================================================================!
module ConfSolv_module
  use crest_parameters
  use crest_data
  use strucrd
  use iomod
  implicit none
  public

  !> ConfSolv helper script PID
  integer,allocatable :: cs_pid
  !> ConfSolv helper script name
  character(len=:),allocatable :: cs_bin
  !> ConfSolv port server port
  integer,allocatable :: cs_port
  !> ConfSolv teardown instruction
  logical :: cs_teardown = .false.

  !> ConfSolv parameter location
  character(len=:),allocatable :: cs_param
  !> ConfSolv solvent   
  character(len=:),allocatable :: cs_solvent

!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!

  subroutine cs_deallocate()
     if(allocated(cs_pid))    deallocate(cs_pid)
     if(allocated(cs_bin))    deallocate(cs_bin)
     if(allocated(cs_port))   deallocate(cs_port)
     if(allocated(cs_param))  deallocate(cs_param)
     if(allocated(cs_solvent))deallocate(cs_solvent)
  end subroutine cs_deallocate

!=================================================!

  function cs_running() result(running)
    implicit none
    logical :: running
    integer :: io
    character(len=:),allocatable :: job
    running = .false.

    io = 0
    if(.not.allocated(cs_bin) .or. &
    &  .not.allocated(cs_pid) .or. &
    &  .not.allocated(cs_port))then
      running = .false.
      return
    endif

    job = trim(cs_bin)//' --test '//to_str(cs_port)

    call command(job, io)

  end function cs_running

!=================================================!
  subroutine cs_shutdown(io)
    implicit none
    integer,intent(out) :: io

    io = 0
    if(cs_teardown .and. allocated(cs_pid) .and. allocated(cs_port)) then
      write(stdout,'(a,i0)') 'Shutting down http://localhost/',cs_port
      call kill(cs_pid,9,io)
      deallocate (cs_pid)
    end if
    
  end subroutine cs_shutdown


!=================================================!

  subroutine cs_deploy()
     implicit none
     character(len=:),allocatable :: job
     integer :: io

     job = 'nohup '//trim(cs_bin)//' -l '//'> confsolv.out 2>/dev/null'//' &'

     call command(job, io)

     !> TODO read port and pid

  end subroutine cs_deploy

!=================================================!

  subroutine cs_write_config(ensname) 
     implicit none
     character(len=*) :: ensname 
     integer :: i,j,k,l,ich,io


     open(newunit=ich, file='config.json')
     write(ich,'(a)') '{'
     call wrjson(ich, 'port',to_str(cs_port))
     call wrjson(ich, 'pid', to_str(cs_pid))
     if(allocated(cs_param))then
     call wrjson(ich, 'trained_model_dir',cs_param)
     endif
     !call wrjson(ich, 'n_jobs', to_str(cs_n_jobs))
     !call wrjson(ich, 'no_ionic_solvents', to_str(cs_ionic_solvents))
     !call wrjson(ich, 'coords_path', cs_coords_path))
     !call wrjson(ich, 'energies_path', cs_energies_path)
     !call wrjson(ich, 'split_path', cs_split_path)
     !call wrjson(ich, 'seed', cs_seed)
     call wrjson(ich, 'ensemble', trim(ensname),.true.)
     write(ich,'(a)') '}'
     close(ich)
   contains
     subroutine wrjson(ch,key,val,fin)
        integer :: ch
        character(len=*) :: key
        character(len=*) :: val
        logical,intent(in),optional :: fin
        if(present(fin))then
          write(ch,'(a,a,a,a,a)') '"',trim(key),'": "',trim(val),'"'
        else
          write(ch,'(a,a,a,a,a)') '"',trim(key),'": "',trim(val),'",'
        endif
     end subroutine wrjson
  end subroutine cs_write_config


!========================================================================================!

end module ConfSolv_module

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!

subroutine confsolv_request(ensname,nall,esolv)
!***********************************************
!* Interface to ConfSolv via requests
!*
!* Input/Output:
!*  ensname - ensemble file
!*  nall    - number of structures in ensemble
!*  esolv   - δΔGsolv for each structure
!***********************************************
  use crest_parameters
  use crest_data
  use strucrd
  use iomod
  use ConfSolv_module
  implicit none
  !> INPUT
  character(len=*),intent(in) :: ensname
  integer,intent(in) :: nall
  !> OUTPUT
  real(wp),intent(out) :: esolv(nall)
  !> LOCAL
  integer :: i,j,k,l,io,ich
  logical :: pr,wr

  esolv(:) = 0.0_wp

!>--- some printout
  if(allocated(cs_pid) .and. allocated(cs_port))then
    write(stdout,'(a,i0,a,i0)') ' ConfSolv server (PID: ',cs_pid,') running at '//&
    & 'http://localhost:',cs_port
  else
   !error stop 'No running ConfSolv server found!'
  endif 
  if(allocated(cs_param))then
    write(stdout,'(a,a)') ' pyTorch checkpoint file located at ',cs_param
  endif
  if(allocated(cs_solvent))then
    write(stdout,'(a,a)') ' Requested δΔGsolv for ',cs_solvent 
  endif
  write(stdout,'(a,a)') ' Processing ensemble file ',trim(ensname)
 

!>---- doing the request
  write(stdout,'(a,a)') ' Writing config file and sending request ...'
  call cs_write_config(ensname)

  call command('echo "ConfSolv integration is a TODO"', io)
  write(*,*) io 

!>--- read δΔGsolv

   !TODO


  return
end subroutine confsolv_request



!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!
