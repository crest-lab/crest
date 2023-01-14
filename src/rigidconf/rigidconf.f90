!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Philipp Pracht, Christopher Zurek, Christoph Bannwarth
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
!
! Routines were adapted from the xtb code (github.com/grimme-lab/xtb)
! under the Open-source software LGPL-3.0 Licencse.
!================================================================================!

!========================================================================================!
!========================================================================================!
subroutine crest_rigidconf(env,tim)
!************************************************************
!* Standalone runtype for conformer generation
!* based on chemoinformatic principles.
!*
!* Input:
!*    env  - CREST's systemdata
!*    tim  - CREST's timer object
!*
!************************************************************
  use crest_parameters
  use crest_data
  use strucrd
  implicit none
  !> INPUT/OUTPUT
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  !> LOCAL 
  type(coord) :: start_mol


!========================================================================================!
  call this_header()
  call tim%start(14,'rule-based conf.')

!>--- some calculation info should be printed out at this point. TODO

!========================================================================================!
!>--- get structure from reference
  call env%ref%to(start_mol)
  write (stdout,*)
  call smallhead('Input structure:')
  call start_mol % append(stdout)
  write (stdout,*)   

!========================================================================================!
!>--- pass the structure to the desired algorithm
  select case( env%rigidconf_algo ) 
  !case ( 1 ) !> "genetic crossing"-type algo

  !  

  case default !> straight-forward generation ("tree"-type algo)

    call rigidconf_tree(env,start_mol)

  end select

!========================================================================================!
  call tim%stop(14)
  return
!========================================================================================!
  contains
!========================================================================================!
  subroutine this_header
     implicit none
      write(stdout,'(/)')
      write(stdout,'(7x,"┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓")')
      write(stdout,'(7x,"┃          R I G I D C O N F           ┃ ")')
      write(stdout,'(7x,"┃     (name is work-in-progress)       ┃ ")')
      !write(stdout,'(7x,"┃               R i C o                ┃")')
      !write(stdout,'(7x,"┃              ConfAcc                 ┃")')
      write(stdout,'(7x,"┃    rule-based conformer generator    ┃")')
      write(stdout,'(7x,"┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛")')
      write(stdout,'(12x,"C.Zurek, C.Bannwarth, P.Pracht")')
      write(stdout,*)
  end subroutine this_header
end subroutine crest_rigidconf
!========================================================================================!
!========================================================================================!
