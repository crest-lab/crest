!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022 - 2023  Philipp Pracht, Gereon Feldmann
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
  use crest_parameters
  use crest_data
  use strucrd
  use calc_type
  use calc_module
  use optimize_module
  use hessian_tools
  use gradreader_module
  use xtb_sc
  implicit none

  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew

  integer :: i,j,io,ich,nat3,n_freqs,i0
  logical :: pr,wr

  character :: fname
!========================================================================================!
  type(calcdata) :: calc

  real(wp) :: energy
  real(wp),allocatable :: hess(:,:,:),freq(:,:),grad(:),grad1(:,:),grad2(:,:),heff(:,:)
!========================================================================================!
  call tim%start(14,'numerical Hessian')
!========================================================================================!
  !call system('figlet numhess')
  write (stdout,*)
  write (stdout,*) "                       _                   "
  write (stdout,*) " _ __  _   _ _ __ ___ | |__   ___  ___ ___ "
  write (stdout,*) "| '_ \| | | | '_ ` _ \| '_ \ / _ \/ __/ __|"
  write (stdout,*) "| | | | |_| | | | | | | | | |  __/\__ \__ \"
  write (stdout,*) "|_| |_|\__,_|_| |_| |_|_| |_|\___||___/___/"
  write (stdout,*)

  call env%ref%to(mol)
  write (stdout,*)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)
!========================================================================================!

  calc = env%calc

  nat3 = mol%nat*3

  if (calc%id .eq. -1) then
    n_freqs = calc%ncalculations+1
  else
    n_freqs = calc%ncalculations
  end if

  allocate (hess(nat3,nat3,calc%ncalculations),source=0.0_wp)
  allocate (freq(nat3,n_freqs),source=0.0_wp)

  if (calc%ncalculations .eq. 1) then
    write (stdout,'(1x,a)',advance='no') 'Calculating numerical Hessian ...'
  else
    write (stdout,'(1x,a)',advance='no') 'Calculating numerical Hessians ...'
  end if

  flush (stdout)

!========================================================================================!

  !>-- Computes numerical Hessians
  call numhess2(mol%nat,mol%at,mol%xyz,calc,hess,io) 

  write (stdout,*) 'done.'
  write (stdout,*)

  !> if calc type is set to -1: Computes the effective Hessian of the 
  !> first two given calculation levels
  if (calc%id .eq. -1) then

    if (calc%ncalculations .gt. 1) then

      allocate (heff(nat3,nat3),source=0.0_wp)
      allocate (grad1(3,mol%nat),grad2(3,mol%nat),source=0.0_wp)

      grad1 = calc%grdtmp(:,:,1)
      grad2 = calc%grdtmp(:,:,2)

      !>-- Computes effective Hessian
      call effective_hessian(mol%nat,nat3,grad1,grad2,hess(:,:,1),hess(:,:,2),heff)

      !>-- Printout
      call print_hessian(heff,nat3,'','effhess') 

      !>-- projection of translation and rotational modes + mass-weighting of Hessian
      call prj_mw_hess(mol%nat,mol%at,nat3,mol%xyz,heff)

      !>-- Comp. of Frequencies
      call frequencies(mol%nat,mol%at,mol%xyz,nat3,calc,heff,freq(:,n_freqs),io)

      !>-- Printout of vibspectrum
      call print_vib_spectrum(mol%nat,mol%at,nat3,mol%xyz,freq(:,n_freqs),'','vibspectrum')

      !>-- Printout of g98 format file
      call print_g98_fake(mol%nat,mol%at,nat3,mol%xyz,freq(:,n_freqs),heff,'','g98.out')

      deallocate (heff)
      deallocate (grad1,grad2)

    else

      write (stdout,*) 'At least two calculation level must be'
      write (stdout,*) 'given for the calculation of the effective Hessian.'
      write (stdout,*)

    end if

  end if

  !> Prints hessian of numerical hessian and does a prj and mass weighting as well for the
  !> computation of normal modes and frequencies
  do i = 1,calc%ncalculations

    if (io .ne. 0) then
      write (stdout,*) 'FAILED!'

    else

      !>-- Prints Hessian
      call print_hessian(hess(:,:,i),nat3,calc%calcs(i)%calcspace,'numhess')

      !>-- Projects and mass-weights the Hessian
      call prj_mw_hess(mol%nat,mol%at,nat3,mol%xyz,hess(:,:,i))

      !>-- Computes the Frequencies
      call frequencies(mol%nat,mol%at,mol%xyz,nat3,calc,hess(:,:,i),freq(:,i),io)

      if (io .ne. 0) then
        write (stdout,*) 'FAILED!'
      else

        !>-- Prints vibspectrum with artifical intensities
        call print_vib_spectrum(mol%nat,mol%at,nat3,mol%xyz,freq(:,i), &
        &    calc%calcs(i)%calcspace,'vibspectrum')

        !>-- Prints g98.out format file
        call print_g98_fake(mol%nat,mol%at,nat3,mol%xyz,freq(:,i),hess(:,:,i), &
        &    calc%calcs(i)%calcspace,'g98.out') 
      end if

    end if

  end do

  write (stdout,*) 'Note: The Hessian is printed as nat*3 blocks of nat*3 entries.'
  write (stdout,*)

  deallocate (hess)
!========================================================================================!
  call tim%stop(14)

  return

end subroutine crest_numhess
