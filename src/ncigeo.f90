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

!========================================================================================!
! calculate settings for an elipsoid wall potential
!========================================================================================!
subroutine wallpot(env)
  use crest_parameters
  use crest_data
  use strucrd,only:coord,get_atlist
  use axis_module
  use crest_calculator
  use wall_setup
  implicit none
  type(systemdata) :: env
  real(wp) :: rabc(3)
  logical :: pr = .false.
  type(coord) :: mol
  type(constraint) :: constr
  logical,allocatable :: atms(:)
  real(wp),parameter :: pi43 = pi*(4.0_wp/3.0_wp)
  real(wp),parameter :: third = 1.0_wp/3.0_wp
  real(wp),parameter :: kdefault = 298.15_wp  !> xtb version doesn't use k
  real(wp),parameter :: betadefault = 50.0_wp !> polynomial default in xtb
!=========================================================!

!>--- read in coord
  call env%ref%to(mol)

!>--- calculate the surrounding ellipsoid
  call boxpot_core(mol,rabc,potscal=env%potscal,potpad=env%potpad)

!>--- write CMA transformed coord file and into memory
  env%ref%xyz = mol%xyz
  call mol%write('coord')

  if (.not.env%legacy) then
!>--- add constraint in calculator framwork
    allocate (atms(mol%nat),source=.false.)
    if(allocated(env%potatlist))then
      call get_atlist(mol%nat,atms,env%potatlist,mol%at) !> selected atoms
    else
      atms(:) = .true. !> all atoms
    endif
    call constr%ellipsoid(mol%nat,atms,rabc,kdefault,betadefault,.true.)
    deallocate (atms)
    call constr%print(stdout)
    call env%calc%add(constr)
  else
!>--- constraint in legacy framework

    allocate (env%cts%pots(10))
    env%cts%pots = ''
    env%cts%NCI = .true.
    write (env%cts%pots(1),'("$wall")')
    write (env%cts%pots(2),'(2x,"potential=logfermi")')
    write (env%cts%pots(3),'(2x,"ellipsoid:",1x,3(g0,",",1x),"all")') rabc
  end if
  return
end subroutine wallpot
