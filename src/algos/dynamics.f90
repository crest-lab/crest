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


subroutine crest_moleculardynamics(env,tim)
  !use iso_fortran_env,only:wp => real64,stdout => output_unit
  use crest_parameters
  use crest_data
  use strucrd
  use calc_type
  use calc_module
  use dynamics_module
  use shake_module
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich
  logical :: pr,wr
!========================================================================================!
  type(calcdata) :: calc
  type(mddata) :: mddat
  type(shakedata) :: shk

  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)

  character(len=80) :: atmp
  character(len=*),parameter :: trjf='crest_dynamics.trj'
!========================================================================================!
  call tim%start(14,'molecular dynamics')
  call env%ref%to(mol)
  write (stdout,*)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)
!========================================================================================!

   !>--- parallelization settings
   if(env%autothreads)then
      call ompautoset(env%threads,8,env%omp,env%MAXRUN,env%threads) 
   endif

  pr = .true.
  mddat = env%mddat
  calc = env%calc
  !>--- check if we have any MD & calculation settings allocated
  if (.not. mddat%requested) then
    write (stdout,*) 'MD requested, but no MD settings present.'
    return
  else if (calc%ncalculations < 0) then
    write (stdout,*) 'MD requested, but no calculation settings present.'
    return
  end if

  !>--- init SHAKE?
  if (mddat%shake) then
    calc%calcs(1)%rdwbo = .true.
    allocate (grad(3,mol%nat),source=0.0_wp)
    call engrad(mol,calc,energy,grad,io)
    deallocate (grad)
    calc%calcs(1)%rdwbo = .false.

    shk%shake_mode = 2
    call move_alloc(calc%calcs(1)%wbo,shk%wbo)

    mddat%shk = shk
    call init_shake(mol%nat,mol%at,mol%xyz,mddat%shk,pr)
    mddat%nshake = mddat%shk%ncons
  end if

  !>--- complete real-time settings to steps
  call mdautoset(mddat,io)
  mddat%trajectoryfile = trjf

  !>--- run the MD
  call dynamics(mol,mddat,calc,pr,io)

  if (io == 0) then
    write (stdout,*) 'MD run completed successfully'
    write (stdout,*) 'Trajectory written to ',trjf
  else
    write (stdout,*) 'MD run terminated with error'
  end if
!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_moleculardynamics
