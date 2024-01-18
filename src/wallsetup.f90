!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2023 Philipp Pracht
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

module wall_setup
  use crest_parameters
  use miscdata
  use strucrd
  use axis_module
  implicit none
  private

  public :: wallpot_core
  public :: wall_dummypot
  public :: boxpot_core

  real(wp),parameter,private :: pi43 = pi*(4.0_wp/3.0_wp)
  real(wp),parameter,private :: third = 1.0_wp/3.0_wp

  logical,parameter,private :: debug = .true.

!========================================================================================!
!========================================================================================!
contains !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine wallpot_core(mol_in,rabc,potscal,atlist)
!********************************************
!* calculate a surrounding ellipsoid for mol
!* the ellipsoid axes are returened as rabc
!********************************************
    implicit none
    !> INPUT
    type(coord) :: mol_in
    real(wp),intent(in),optional :: potscal
    logical,intent(in),optional :: atlist(mol_in%nat)
    !> OUTPUT
    real(wp),intent(out) :: rabc(3)
    !> LOCAL
    type(coord) :: mol
    integer :: nat
    real(wp) :: eaxr(3)
    real(wp) :: rmax
    real(wp) :: sola,vtot
    real(wp) :: r,roff,boxr
    real(wp) :: pshape,pscal
    real(wp) :: natfac,erffac,erfscal
    logical,parameter :: pr = .false.

    pshape = 1.0d0
    eaxr = 0.0d0
    if (present(potscal)) then
      pscal = potscal
    else
      pscal = 1.0_wp
    end if

!>-- obtain structure
    if (present(atlist)) then
      mol = mol_in%cutout(atlist)
    else
      mol = mol_in
    end if

!>--- CMA trafo
    nat = mol%nat
    call axis(pr,mol%nat,mol%at,mol%xyz,eaxr)
    sola = sqrt(1.0d0+(abs(eaxr(1)-eaxr(3)))/(abs(eaxr(1)+eaxr(2)+eaxr(3))/3.0d0))
    call getmaxdist(mol%nat,mol%xyz,mol%at,rmax)

!>--- volume as sum of speherical atoms (crude approximation)
    vtot = volsum(mol%nat,mol%at)

!>--- calculate ellipsoid
    roff = sola*vtot/1000.0d0
    boxr = ((sola*vtot)/pi43)**third+roff+rmax*0.5_wp
    r = (boxr**3/(eaxr(1)*eaxr(2)*eaxr(3)))**third  ! volume of ellipsoid = volume of sphere
    rabc = eaxr**pshape/sum((eaxr(1:3))**pshape)

!> scale pot size by number of atoms
!> pure empirics
    natfac = 0.08_wp*nat-0.08_wp*50.0_wp
    erffac = erf(natfac)*0.25_wp
    erfscal = 1.0_wp-erffac
!> erfscal is ~ 1.25 for systems <<50 atoms
!> erfscal is ~ 0.75 for systems >>50 atoms
    rabc = eaxr*r*pscal*erfscal*1.5_wp

    return
  end subroutine wallpot_core

!========================================================================================!

  subroutine boxpot_core(mol,rabc,potscal,potpad)
!***********************************************
!* simplified routine to set up wall potentials
!* draws a box around the molecule and adds 3AA
!***********************************************
    implicit none
    type(coord) :: mol
    real(wp),intent(out) :: rabc(3)
    real(wp),intent(in),optional :: potscal
    real(wp),intent(in),optional :: potpad
    real(wp) :: box(3,3),pscal,ppad
    real(wp) :: boxvol

    if (present(potscal)) then
      pscal = potscal
    else
      pscal = 1.0_wp
    end if
    if (present(potpad)) then
      ppad = potpad
    else
      ppad = 0.0_wp
    end if

!>--- calculate the box
    boxvol = getbox2(mol%nat,mol%xyz,mol%at,box)
    rabc(1:3) = box(1:3,3)/2.0_wp  !> radius, half the box length

!>--- add 3AA and scale further by a user-defined factor, if necessary
    rabc = (rabc+3.0_wp*aatoau+ppad)*pscal

    if (debug) call wall_dummypot(rabc,mol%xyz,mol%at,mol%nat)

  end subroutine boxpot_core

!========================================================================================!

  subroutine getmaxdist(n,xyz,at,dist)
    implicit none
    integer :: n
    integer :: at(n)
    integer :: i
    real(wp) :: xyz(3,n)
    real(wp) :: dist
    real(wp) :: dum

    dist = 0.0d0

    do i = 1,n
      dum = sqrt(xyz(1,i)**2+xyz(2,i)**2+xyz(3,i)**2)
      dum = dum+rcov(at(i))
      if (dum .gt. dist) dist = dum
    end do

  end subroutine getmaxdist

!========================================================================================!

  subroutine wall_dummypot(ellips,xyz,at,nat)
    implicit none

    integer :: i,j,num_points
    integer :: nat,ich
    integer :: at(nat)

    real(wp) :: xyz(3,nat)
    real(wp) :: r,dum
    real(wp) :: ellips(3),a,b,c
    real(wp) :: theta,phi,vol,surf,p
    real(wp) :: point(3)

!>--- ellipsoid has to statisfy x²/a² + y²/b² + z²/c² = 1
!>--- the volume is 4/3*pi*a*b*c
    call init_random_seed()
    open (file='wall.coord',newunit=ich)
    write (ich,'(a)') "$coord"

    do i = 1,nat
      write (ich,'(3F24.12,5x,a2)') xyz(1:3,i),i2e(at(i))
    end do

    point = 0.0d0
    ! Generate points on the surface of the ellipsoid
    a = ellips(1)
    b = ellips(2)
    c = ellips(3)
    num_points = 20
    vol = (4.0d0/3.0d0)*pi*a*b*c
    !write(*,*) vol,sqrt(vol),sqrt(sqrt(vol))
    !p = 1.6075d0
    !surf = 4.0d0*pi*( ((a**p)*(b**p) + (a**p)*(c**p) + (b**p)*(c**p))/3.0d0)**(1.0d0/p)
    !write(*,*) surf,sqrt(surf),sqrt(sqrt(surf))
    num_points = ceiling(sqrt(sqrt(vol)))
    num_points = max(num_points,8)
    do i = 1,num_points
      ! Generate evenly spaced spherical coordinates
      theta = real((i-1))*pi/real((num_points))
      do j = 1,num_points
        ! Generate evenly spaced spherical coordinates
        phi = 2.0d0*pi*real((j-1))/real((num_points))

        ! Calculate Cartesian coordinates using parametric equations
        point(1) = a*sin(theta)*cos(phi)
        point(2) = b*sin(theta)*sin(phi)
        point(3) = c*cos(theta)
        write (ich,'(3F24.12,5x,a2)') point(1:3),'he'
      end do
    end do

    write (ich,'(a)') "$end"
    close (ich)
  end subroutine wall_dummypot

  SUBROUTINE init_random_seed()
    INTEGER :: i,n,clock
    INTEGER,DIMENSION(:),ALLOCATABLE :: seed

    CALL RANDOM_SEED(size=n)
    ALLOCATE (seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    seed = clock+37*(/(i-1,i=1,n)/)
    CALL RANDOM_SEED(PUT=seed)

    DEALLOCATE (seed)
  END SUBROUTINE

!========================================================================================!

  function getbox(n,xyz,box)
!************************************
!* set up  a box around the molecule.
!* the return value is the box volume
!*************************************
    implicit none
    real(wp) :: getbox
    integer,intent(in)  :: n
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(out) :: box(3,3)
    integer :: i
    box = 0.0_wp
    getbox = 1.0_wp
    do i = 1,3 !i are the X,Y,Z axis
      box(i,1) = maxval(xyz(i,:))
      box(i,2) = minval(xyz(i,:))
      box(i,3) = abs(box(i,1)-box(i,2)) !side length
      getbox = getbox*box(i,3)          !to volume
    end do
    return
  end function getbox

!========================================================================================!

  function getbox2(n,xyz,at,box)
!************************************
!* set up  a box around the molecule.
!* the return value is the box volume
!*************************************
    implicit none
    real(wp) :: getbox2
    integer,intent(in)  :: n
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(out) :: box(3,3)
    integer,intent(in) :: at(n)
    integer :: i
    real(wp),allocatable :: rat(:)
    real(wp) :: rcovmax
    allocate (rat(n))
    do i = 1,n
      rat(i) = rcov(at(i))
    end do
    rcovmax = maxval(rat)*2.0_wp
    box = 0.0_wp
    getbox2 = 1.0_wp
    do i = 1,3 !i are the X,Y,Z axis
      box(i,1) = maxval(xyz(i,:))
      box(i,2) = minval(xyz(i,:))
      box(i,3) = abs(box(i,1)-box(i,2)) !side length
      if (box(i,3) .lt. rcovmax) then
        box(i,1) = rcovmax*2.0_wp
        box(i,2) = -rcovmax*2.0_wp
        box(i,3) = rcovmax*4.0_wp
      end if
      getbox2 = getbox2*box(i,3)          !to volume
    end do
    deallocate (rat)
    return
  end function getbox2

!========================================================================================!

  function volsum(n,at)
!***********************************************
!* get the volume simply as a sum of
!* spherical atom volumes (crude approximation)
!***********************************************
    implicit none
    real(wp) :: volsum
    integer,intent(in)  :: n
    integer,intent(in) :: at(n)
    integer :: i
    real(wp) :: r
    !> using D3 vdW radii
    volsum = 0.0_wp
    do i = 1,n
      r = vdw_d3(at(i))
      volsum = volsum+(4.0_wp/3.0_wp)*(3.14159265359_wp*r**3)
    end do
    return
  end function volsum

!========================================================================================!
!========================================================================================!
end module wall_setup
