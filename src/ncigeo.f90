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
  use strucrd,only:rdnat,rdcoord,wrc0,coord
  use axis_module
  use crest_calculator
  use wall_setup
  implicit none
  type(systemdata) :: env

  integer :: nat
  integer,allocatable  :: at(:)
  real(wp),allocatable :: xyz(:,:)
  real(wp) :: eaxr(3),rabc(3)
  real(wp) :: rmax
  real(wp) :: sola,vtot
  real(wp) :: r,roff,boxr
  real(wp) :: pshape
  real(wp) :: volsum
  real(wp) :: natfac,erffac,erfscal
  logical :: pr = .false.

  type(coord) :: mol
  type(constraint) :: constr
  logical,allocatable :: atms(:)

  real(wp),parameter :: pi43 = pi * (4.0_wp / 3.0_wp)
  real(wp),parameter :: third = 1.0_wp / 3.0_wp
  real(wp),parameter :: kdefault = 1.0_wp  !> xtb version doesn't use k
  real(wp),parameter :: alphadefault = 30.0_wp !> polynomial default in xtb
!=========================================================!
  allocate (env%cts%pots(10))
  env%cts%pots = ''
  env%cts%NCI = .true.

!>--- read in coord
  call env%ref%to(mol)

!>--- calculate the surrounding ellipsoid
  call wallpot_core(mol,rabc, potscal=env%potscal) 

!>--- write CMA transformed coord file
  call wrc0('coord',env%nat,mol%at,mol%xyz)
  call wall_dummypot(rabc,mol%xyz,mol%at,mol%nat)

!>--- constraint in legacy framework
  write (env%cts%pots(1),'("$wall")')
  write (env%cts%pots(2),'(2x,"potential=polynomial")')
  write (env%cts%pots(3),'(2x,"ellipsoid:",1x,3(g0,",",1x),"all")') rabc

!>--- add constraint in calculator framwork
  if(.not.env%legacy)then
    allocate(atms(mol%nat), source=.true.) !> all atoms
    call constr%ellipsoid( mol%nat, atms,rabc, kdefault,alphadefault,.false.)
    deallocate(atms)
    !call constr%print(stdout)
    call env%calc%add(constr)
  endif
  return
end subroutine wallpot
