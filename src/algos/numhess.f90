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

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> Implementation of a numerical Hessian calculation
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
!> Input/Output:
!>  env  -  crest's systemdata object
!>  tim  -  timer object
!>-----------------------------------------------
subroutine crest_numhess(env,tim)
  use iso_fortran_env,only: wp =>real64,stdout=>output_unit
  use crest_data
  use strucrd 
  use calc_type
  use calc_module
  use optimize_module
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,nat3
  logical :: pr,wr
!========================================================================================!
  type(calcdata) :: calc

  real(wp) :: energy
  real(wp),allocatable :: hess(:,:)
!========================================================================================!
  call tim%start(14,'numerical Hessian') 
!========================================================================================!
  !call system('figlet numhess')
  write(stdout,*)
  write(stdout,*) "                       _                   "
  write(stdout,*) " _ __  _   _ _ __ ___ | |__   ___  ___ ___ "
  write(stdout,*) "| '_ \| | | | '_ ` _ \| '_ \ / _ \/ __/ __|"
  write(stdout,*) "| | | | |_| | | | | | | | | |  __/\__ \__ \"
  write(stdout,*) "|_| |_|\__,_|_| |_| |_|_| |_|\___||___/___/"
  write(stdout,*)
                                           
  call env%ref%to(mol)
  write(stdout,*)
  write(stdout,*) 'Input structure:'
  call mol%append(stdout)
  write(stdout,*) 
!========================================================================================!

  calc = env%calc
  nat3 = mol%nat * 3 
  allocate(hess(nat3, nat3), source=0.0_wp) 

  write(stdout,*)
  write(stdout,'(1x,a)',advance='no') 'Calculating numerical Hessian ... '
  flush(stdout)

  call numhess(mol%nat,mol%at,mol%xyz,calc,hess,io)

  if(io.ne.0)then
    write(stdout,*) 'FAILED!'
  else 
    write(stdout,*) 'done.'
    write(stdout,'(1x,a)',advance='no') 'Will be written to file "numhess" ...'
    flush(stdout)
    open(newunit=ich,file='numhess')
    do i=1,nat3
      k = 0
      do j=1,nat3
        k = k + 1
        if( k .le. 4 )then
        write(ich,'(f16.8)',advance='no') hess(i,j)
        else
        write(ich,'(f16.8)') hess(i,j)
        k = 0 
        endif
      enddo
      if(k .ne. 0)then 
      write(ich,*)
      endif
    enddo
    close(ich)

    write(stdout,*) 'done.'
    write(stdout,*) 'Note: This is just the plain non-mass-weighted Hessian!'
    write(stdout,*) '      The Hessian is printed as nat*3 blocks of nat*3 entries.'
  endif
 
  deallocate(hess)
!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_numhess
