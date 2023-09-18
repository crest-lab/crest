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
!!============================================================================
!! only calculate the ellipsoide wall potential
!!============================================================================
!subroutine wallpot_calc(nat,at,xyz,rabc)
!  use iso_fortran_env,wp => real64
!  use axis_module
!  implicit none
!  integer,intent(in)     :: nat
!  integer,intent(inout)  :: at(nat)
!  real(wp),intent(inout) :: xyz(3,nat)
!  real(wp),intent(out)   :: rabc(3)     !potential ellipsoide-axis
!  real(wp) :: eaxr(3)
!  real(wp) :: rmax
!  real(wp) :: sola,vtot
!  real(wp) :: r,roff,boxr
!  real(wp) :: pshape
!  real(wp) :: volsum
!  real(wp) :: natfac,erffac,erfscal
!  logical :: pr = .false.
!  real(wp),parameter :: pi43 = 3.1415926540_wp * (4.0_wp / 3.0_wp)
!  real(wp),parameter :: pi = 3.1415926540_wp
!  real(wp),parameter :: third = 1.0_wp / 3.0_wp
!  rabc = 0.0d0
!  pshape = 1.0d0
!  eaxr = 0.0d0
!
!!--- CMA trafo
!  call axis(pr,nat,at,xyz,eaxr)
!  sola = sqrt(1.0d0 + (eaxr(1) - eaxr(3)) / ((eaxr(1) + eaxr(2) + eaxr(3)) / 3.0d0))
!  call getmaxdist(nat,xyz,at,rmax)
!  vtot = volsum(nat,at) !--- volume as sum of spherical atoms (crude approximation)
!!--- calculate ellipsoid
!  roff = sola * vtot / 1000.0d0
!  boxr = ((sola * vtot) / pi43)**third + roff + rmax * 0.5
!  r = (boxr**3 / (eaxr(1) * eaxr(2) * eaxr(3)))**third  ! volume of ellipsoid = volume of sphere
!  rabc = eaxr**pshape / sum((eaxr(1:3))**pshape)
!!--- scale pot size by number of atoms
!  natfac = 0.08_wp * nat - 0.08_wp * 50.0_wp
!  erffac = erf(natfac) * 0.25_wp
!  erfscal = 1.0_wp - erffac
!
!  rabc = eaxr * r * erfscal * 1.5_wp
!  return
!end subroutine wallpot_calc
!
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!c
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!subroutine getmaxdist(n,xyz,at,dist)
!  use iso_fortran_env,wp => real64
!  use miscdata, only: rcov 
!  implicit none
!  integer :: n
!  integer :: at(n)
!  integer :: i
!  real(wp) :: xyz(3,n)
!  real(wp) :: dist
!  real(wp) :: dum
!
!  dist = 0.0d0
!
!  do i = 1,n
!    dum = sqrt(xyz(1,i)**2 + xyz(2,i)**2 + xyz(3,i)**2)
!    dum = dum + rcov(at(i))
!    if (dum .gt. dist) dist = dum
!  end do
!
!end subroutine getmaxdist
!
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!c
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!subroutine dummypot(ellips,xyz,at,nat)
!  use iso_fortran_env,wp => real64
!  use crest_data
!  use strucrd,only:i2e
!  implicit none
!
!  integer :: i
!  integer :: nat,ich
!  integer :: at(nat)
!
!  real(wp) :: xyz(3,nat)
!  real(wp) :: r,dum
!  real(wp) :: ellips(3)
!
!  real(wp) :: point(3)
!
!!---- ellipsoide has to statisfy x²/a² + y²/b² + z²/c² = 1
!
!  call init_random_seed()
!  open (file='wall.coord',newunit=ich)
!  write (ich,'(a)') "$coord"
!
!  do i = 1,nat
!    write (ich,'(3F24.12,5x,a2)') xyz(1:3,i),i2e(at(i))
!  end do
!
!  point = 0.0d0
!  do i = 1,100
!    call random_number(r)
!    point(1) = r * ellips(1)
!    do
!      call random_number(r)
!      point(2) = r * ellips(2)
!      dum = (point(1)**2 / ellips(1)**2)  &
!      &    + (point(2)**2 / ellips(2)**2)
!      if (dum .lt. 1.0d0) exit
!    end do
!    dum = 1.0d0 - (point(1)**2 / ellips(1)**2) - (point(2)**2 / ellips(2)**2)
!    dum = dum * ellips(3)**2
!
!    point(3) = sqrt(dum)
!
!    write (ich,'(3F24.12,5x,a2)') point(1:3),'he'
!    write (ich,'(3F24.12,5x,a2)') - 1 * point(1:3),'he'
!    write (ich,'(3F24.12,5x,a2)') - 1 * point(1),point(2),point(3),'he'
!    write (ich,'(3F24.12,5x,a2)') point(1),-1 * point(2),point(3),'he'
!    write (ich,'(3F24.12,5x,a2)') point(1),point(2),-1 * point(3),'he'
!    write (ich,'(3F24.12,5x,a2)') - 1 * point(1),-1 * point(2),point(3),'he'
!    write (ich,'(3F24.12,5x,a2)') - 1 * point(1),point(2),-1 * point(3),'he'
!    write (ich,'(3F24.12,5x,a2)') point(1),-1 * point(2),-1 * point(3),'he'
!  end do
!
!  write (ich,'(a)') "$end"
!  close (ich)
!end subroutine dummypot
!
!SUBROUTINE init_random_seed()
!  INTEGER :: i,n,clock
!  INTEGER,DIMENSION(:),ALLOCATABLE :: seed
!
!  CALL RANDOM_SEED(size=n)
!  ALLOCATE (seed(n))
!
!  CALL SYSTEM_CLOCK(COUNT=clock)
!
!  seed = clock + 37 * (/(i - 1,i=1,n)/)
!  CALL RANDOM_SEED(PUT=seed)
!
!  DEALLOCATE (seed)
!END SUBROUTINE
!
!!---- set up  a box around the molecule. the return value is the box volume
!function getbox(n,xyz,box)
!  use iso_fortran_env,wp => real64
!  implicit none
!  real(wp) :: getbox
!  integer,intent(in)  :: n
!  real(wp),intent(in) :: xyz(3,n)
!  real(wp),intent(out) :: box(3,3)
!  integer :: i
!  box = 0.0_wp
!  getbox = 1.0_wp
!  do i = 1,3 !i are the X,Y,Z axis
!    box(i,1) = maxval(xyz(i,:))
!    box(i,2) = minval(xyz(i,:))
!    box(i,3) = abs(box(i,1) - box(i,2)) !side length
!    getbox = getbox * box(i,3)          !to volume
!  end do
!  return
!end function getbox
!!---- set up  a box around the molecule. the return value is the box volume
!function getbox2(n,xyz,at,box)
!  use iso_fortran_env,wp => real64
!  use miscdata, only: rcov
!  implicit none
!  real(wp) :: getbox2
!  integer,intent(in)  :: n
!  real(wp),intent(in) :: xyz(3,n)
!  real(wp),intent(out) :: box(3,3)
!  integer,intent(in) :: at(n)
!  integer :: i
!  real(wp),allocatable :: rat(:)
!  real(wp) :: rcovmax
!  do i = 1,n
!    rat(i) = rcov(at(i))
!  end do
!  rcovmax = maxval(rat)
!  box = 0.0_wp
!  getbox2 = 1.0_wp
!  do i = 1,3 !i are the X,Y,Z axis
!    box(i,1) = maxval(xyz(i,:))
!    box(i,2) = minval(xyz(i,:))
!    box(i,3) = abs(box(i,1) - box(i,2)) !side length
!    if (box(i,3) .lt. 1d-6) then
!      box(i,1) = rcovmax
!      box(i,2) = -rcovmax
!      box(i,3) = rcovmax * 2.0_wp
!    end if
!    getbox2 = getbox2 * box(i,3)          !to volume
!  end do
!  deallocate (rat)
!  return
!end function getbox2
!
!!---- get the volume simply as a sum of spherical atom volumes (crude approximation)
!function volsum(n,at)
!  use iso_fortran_env,wp => real64
!  use miscdata, only: vdw_d3
!  implicit none
!  real(wp) :: volsum
!  integer,intent(in)  :: n
!  integer,intent(in) :: at(n)
!  integer :: i
!  real(wp) :: r
!  !> using D3 vdW radii
!  volsum = 0.0_wp
!  do i = 1,n
!    r = vdw_d3(at(i))
!    volsum = volsum + (4.0_wp / 3.0_wp) * (3.14159265359_wp * r**3)
!  end do
!  return
!end function volsum
!
