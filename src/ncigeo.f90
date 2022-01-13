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
  use iso_fortran_env,wp => real64
  use crest_data
  use strucrd,only:rdnat,rdcoord,wrc0
  use axis_module
  implicit none
  !type(options) :: opt
  type(systemdata) :: env

  integer :: i,j,k,l
  integer :: nat
  integer,allocatable  :: at(:)
  real(wp),allocatable :: xyz(:,:)
  real(wp) :: eaxr(3),rabc(3)
  real(wp) :: rmax
  real(wp) :: sola,atot,vtot
  real(wp) :: r,roff,boxr
  real(wp) :: pshape
  real(wp) :: dumbox(3,3),getbox
  real(wp) :: volsum
  real(wp) :: natfac,erffac,erfscal
  logical :: pr = .false.

  real(wp),parameter :: pi43 = 3.1415926540_wp * (4.0_wp / 3.0_wp)
  real(wp),parameter :: pi = 3.1415926540_wp
  real(wp),parameter :: third = 1.0_wp / 3.0_wp

  allocate (env%cts%pots(10))
  env%cts%pots = ''
  env%cts%NCI = env%NCI

  pshape = 1.0d0

  nat = env%nat
  eaxr = 0.0d0

!--- read in coord
  allocate (xyz(3,nat),at(nat))
  call rdcoord('coord',nat,at,xyz)

!--- CMA trafo
  !call axistrf(nat,nat,at,xyz)
  !call axis(nat,at,xyz) 
  call axis(pr,nat,at,xyz,eaxr)
  sola = sqrt(1.0d0 + (abs(eaxr(1) - eaxr(3))) / (abs(eaxr(1) + eaxr(2) + eaxr(3)) / 3.0d0))
  call getmaxdist(nat,xyz,at,rmax)
  vtot = volsum(nat,xyz,at) !--- volume as sum of speherical atoms (crude approximation)

!--- calculate ellipsoid

  !env%potscal
  roff = sola * vtot / 1000.0d0
  !write(*,*) 'roff',roff
  boxr = ((sola * vtot) / pi43)**third + roff + rmax * 0.5_wp
  !write(*,*) 'boxr',boxr
  r = (boxr**3 / (eaxr(1) * eaxr(2) * eaxr(3)))**third  ! volume of ellipsoid = volume of sphere
  !write(*,*) 'r',r
  rabc = eaxr**pshape / sum((eaxr(1:3))**pshape)
  !write(*,*) 'rabc',rabc

  !> scale pot size by number of atoms
  !> pure empirics
  natfac = 0.08_wp * nat - 0.08_wp * 50.0_wp
  erffac = erf(natfac) * 0.25_wp
  erfscal = 1.0_wp - erffac
  !> erfscal is ~ 1.25 for systems <<50 atoms
  !> erfscal is ~ 0.75 for systems >>50 atoms
  rabc = eaxr * r * env%potscal * erfscal * 1.5_wp
  !write(*,*) 'rabc',rabc

!--- write CMA transformed coord file
  call wrc0('coord',env%nat,at,xyz)
  call dummypot(env,rabc,xyz,at,env%nat)

  deallocate (at,xyz)

  write (env%cts%pots(1),'("$wall")')
  write (env%cts%pots(2),'(2x,"potential=polynomial")')
  write (env%cts%pots(3),'(2x,"ellipsoid:",1x,3(g0,",",1x),"all")') rabc

  return
end subroutine wallpot
!============================================================================
! only calculate the ellipsoide wall potential
!============================================================================
subroutine wallpot_calc(nat,at,xyz,rabc)
  use iso_fortran_env,wp => real64
  use axis_module
  implicit none
  integer,intent(in)     :: nat
  integer,intent(inout)  :: at(nat)
  real(wp),intent(inout) :: xyz(3,nat)
  real(wp),intent(out)   :: rabc(3)     !potential ellipsoide-axis
  integer :: i,j,k,l
  real(wp) :: eaxr(3)
  real(wp) :: rmax
  real(wp) :: sola,atot,vtot
  real(wp) :: r,roff,boxr
  real(wp) :: pshape
  real(wp) :: dumbox(3,3),getbox
  real(wp) :: volsum
  real(wp) :: natfac,erffac,erfscal
  logical :: pr = .false.
  real(wp),parameter :: pi43 = 3.1415926540_wp * (4.0_wp / 3.0_wp)
  real(wp),parameter :: pi = 3.1415926540_wp
  real(wp),parameter :: third = 1.0_wp / 3.0_wp
  rabc = 0.0d0
  pshape = 1.0d0
  eaxr = 0.0d0

!--- CMA trafo
  !call axistrf(nat,nat,at,xyz)
  call axis(pr,nat,at,xyz,eaxr)
  sola = sqrt(1.0d0 + (eaxr(1) - eaxr(3)) / ((eaxr(1) + eaxr(2) + eaxr(3)) / 3.0d0))
  call getmaxdist(nat,xyz,at,rmax)
  vtot = volsum(nat,xyz,at) !--- volume as sum of speherical atoms (crude approximation)
!--- calculate ellipsoid
  roff = sola * vtot / 1000.0d0
  boxr = ((sola * vtot) / pi43)**third + roff + rmax * 0.5
  r = (boxr**3 / (eaxr(1) * eaxr(2) * eaxr(3)))**third  ! volume of ellipsoid = volume of sphere
  rabc = eaxr**pshape / sum((eaxr(1:3))**pshape)
!--- scale pot size by number of atoms
  natfac = 0.08_wp * nat - 0.08_wp * 50.0_wp
  erffac = erf(natfac) * 0.25_wp
  erfscal = 1.0_wp - erffac

  rabc = eaxr * r * erfscal * 1.5_wp
  return
end subroutine wallpot_calc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getmaxdist(n,xyz,at,dist)
  use iso_fortran_env,wp => real64
  implicit none
  integer :: n
  integer :: at(n)
  integer :: i,j,k,l
  real(wp) :: xyz(3,n)
  real(wp) :: dist
  real(wp) :: dum
  real(wp),allocatable :: rcov(:)

  allocate (rcov(94))
  call setrcov(rcov)
  dist = 0.0d0

  do i = 1,n
    dum = sqrt(xyz(1,i)**2 + xyz(2,i)**2 + xyz(3,i)**2)
    dum = dum + rcov(at(i))
    if (dum .gt. dist) dist = dum
  end do

  deallocate (rcov)

end subroutine getmaxdist

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine dummypot(env,ellips,xyz,at,nat)
  use iso_fortran_env,wp => real64
  use crest_data
  use strucrd,only:i2e
  implicit none
  !type(options) :: opt
  type(systemdata) :: env

  integer :: i,j,k,l
  integer :: imax
  integer :: nat,ich
  integer :: at(nat)

  real(wp) :: xyz(3,nat)
  real(wp) :: r,dum
  real(wp) :: ellips(3)

  real(wp) :: point(3)

!---- ellipsoide has to statisfy x²/a² + y²/b² + z²/c² = 1

  call init_random_seed()
  open (file='wall.coord',newunit=ich)
  write (ich,'(a)') "$coord"

  do i = 1,nat
    write (ich,'(3F24.12,5x,a2)') xyz(1:3,i),i2e(at(i))
  end do

  point = 0.0d0
  do i = 1,100
    call random_number(r)
    point(1) = r * ellips(1)
    do
      call random_number(r)
      point(2) = r * ellips(2)
      dum = (point(1)**2 / ellips(1)**2)  &
      &    + (point(2)**2 / ellips(2)**2)
      if (dum .lt. 1.0d0) exit
    end do
    dum = 1.0d0 - (point(1)**2 / ellips(1)**2) - (point(2)**2 / ellips(2)**2)
    dum = dum * ellips(3)**2

    point(3) = sqrt(dum)

    write (ich,'(3F24.12,5x,a2)') point(1:3),'he'
    write (ich,'(3F24.12,5x,a2)') - 1 * point(1:3),'he'
    write (ich,'(3F24.12,5x,a2)') - 1 * point(1),point(2),point(3),'he'
    write (ich,'(3F24.12,5x,a2)') point(1),-1 * point(2),point(3),'he'
    write (ich,'(3F24.12,5x,a2)') point(1),point(2),-1 * point(3),'he'
    write (ich,'(3F24.12,5x,a2)') - 1 * point(1),-1 * point(2),point(3),'he'
    write (ich,'(3F24.12,5x,a2)') - 1 * point(1),point(2),-1 * point(3),'he'
    write (ich,'(3F24.12,5x,a2)') point(1),-1 * point(2),-1 * point(3),'he'
  end do

  write (ich,'(a)') "$end"
  close (ich)
end subroutine dummypot

SUBROUTINE init_random_seed()
  INTEGER :: i,n,clock
  INTEGER,DIMENSION(:),ALLOCATABLE :: seed

  CALL RANDOM_SEED(size=n)
  ALLOCATE (seed(n))

  CALL SYSTEM_CLOCK(COUNT=clock)

  seed = clock + 37 * (/(i - 1,i=1,n)/)
  CALL RANDOM_SEED(PUT=seed)

  DEALLOCATE (seed)
END SUBROUTINE
!---- statically merge two molecules A and B and paste into C
subroutine quickdock(Aname,Bname,Cname)
  use iso_fortran_env,wp => real64
  use strucrd,only:rdnat,rdcoord,wrc0
  use crest_data,only:bohr
  use axis_module
  implicit none

  character(len=*) :: Aname
  character(len=*) :: Bname
  character(len=*) :: Cname

  integer :: i,j,k,l
  integer :: ich,ich2
  integer :: nA,nB

  real(wp),allocatable :: A(:,:)
  real(wp),allocatable :: B(:,:)
  integer,allocatable  :: atA(:)
  integer,allocatable  :: atB(:)
  real(wp),allocatable :: boxA(:,:)
  real(wp),allocatable :: boxB(:,:)
  real(wp) :: volA,volB

  integer :: nC
  integer,allocatable :: atC(:)
  real(wp),allocatable :: C(:,:)

  real(wp),allocatable :: dum(:)

  real(wp) :: getbox
  real(wp) :: X(3)
  real(wp) :: shift

!---- write
  write (*,'(2x,a,a,1x,a,a)') '-dock ',trim(Aname),trim(Bname), &
  & ' : merging structures to create input coord'
!---- some init stuff
  allocate (boxA(3,3),boxB(3,3),dum(3))

!---- Read in the files
  call rdnat(Aname,nA)
  call rdnat(Bname,nB)
  allocate (A(3,nA),atA(nA),B(3,nB),atB(nB))
  call rdcoord(Aname,nA,atA,A)
  call rdcoord(Bname,nB,atB,B)
!---- Transform to CMA
  call axis(nA,atA,A,A,dum)
  call axis(nB,atB,B,B,dum)

  volA = getbox(nA,A,boxA)
  volB = getbox(nB,B,boxB)

  !write(*,*) 'Volume of A:',volA
  !write(*,*) 'Volume of B:',volB

  call wrc0('A.new',nA,atA,A)
  call wrc0('B.new',nB,atB,B)
!---- shift molecules
  X = (/0.0_wp,0.0_wp,1.0_wp/) !translation along z
  shift = 2.5_wp / bohr
  if (volA .ge. volB) then  !shift the larger molecule above the smaller one
    call shiftXvec(boxA,boxB,X,shift)
    call shiftmol(A,nA,B,nB,X)
  else
    call shiftXvec(boxB,boxA,X,shift)
    call shiftmol(B,nB,A,nA,X)
  end if
!--- produce merged geometry C
  nC = nA + nB
  allocate (atC(nC),C(3,nC))
  do i = 1,nA
    atC(i) = atA(i)
    C(1:3,i) = A(1:3,i)
  end do
  do j = 1,nB
    k = j + nA !offset
    atC(k) = atB(j)
    C(1:3,k) = B(1:3,j)
  end do
  call wrc0(Cname,nC,atC,C)

  deallocate (C,atC)
  deallocate (atB,B,atA,A)
  deallocate (dum,boxB,boxA)
  return
end subroutine quickdock

!---- set up  a box around the molecule. the return value is the box volume
function getbox(n,xyz,box)
  use iso_fortran_env,wp => real64
  implicit none
  real(wp) :: getbox
  integer,intent(in)  :: n
  real(wp),intent(in) :: xyz(3,n)
  real(wp),intent(out) :: box(3,3)
  integer :: i,j
  box = 0.0_wp
  getbox = 1.0_wp
  do i = 1,3 !i are the X,Y,Z axis
    box(i,1) = maxval(xyz(i,:))
    box(i,2) = minval(xyz(i,:))
    box(i,3) = abs(box(i,1) - box(i,2)) !side length
    getbox = getbox * box(i,3)          !to volume
  end do
  return
end function getbox
!---- set up  a box around the molecule. the return value is the box volume
function getbox2(n,xyz,at,box)
  use iso_fortran_env,wp => real64
  implicit none
  real(wp) :: getbox2
  integer,intent(in)  :: n
  real(wp),intent(in) :: xyz(3,n)
  real(wp),intent(out) :: box(3,3)
  integer,intent(in) :: at(n)
  integer :: i,j
  real(wp),allocatable :: rcov(:),rat(:)
  real(wp) :: rcovmax
  allocate (rcov(94),rat(n))
  call setrcov(rcov)
  do i = 1,n
    rat(i) = rcov(at(i))
  end do
  rcovmax = maxval(rat)
  box = 0.0_wp
  getbox2 = 1.0_wp
  do i = 1,3 !i are the X,Y,Z axis
    box(i,1) = maxval(xyz(i,:))
    box(i,2) = minval(xyz(i,:))
    box(i,3) = abs(box(i,1) - box(i,2)) !side length
    if (box(i,3) .lt. 1d-6) then
      box(i,1) = rcovmax
      box(i,2) = -rcovmax
      box(i,3) = rcovmax * 2.0_wp
    end if
    getbox2 = getbox2 * box(i,3)          !to volume
  end do
  deallocate (rat,rcov)
  return
end function getbox2

!---- get the volume simply as a sum of spherical atom volumes (crude approximation)
function volsum(n,xyz,at)
  use iso_fortran_env,wp => real64
  implicit none
  real(wp) :: volsum
  integer,intent(in)  :: n
  real(wp),intent(in) :: xyz(3,n)
  integer,intent(in) :: at(n)
  integer :: i,j
  real(wp),allocatable :: rcov(:),rat(:)
  real(wp) :: r
  allocate (rcov(94))
  !--- D3 radii in Bohr
  rcov = (/ &
   &  2.18230009,1.73469996,3.49559999,3.09820008,3.21600008, &
   &  2.91030002,2.62249994,2.48169994,2.29959989,2.13739991, &
   &  3.70819998,3.48390007,4.01060009,3.79169989,3.50169992, &
   &  3.31069994,3.10459995,2.91479993,4.24109983,4.10349989, &
   &  3.89030004,3.76419997,3.72110009,3.44140005,3.54620004, &
   &  3.44210005,3.43269992,3.34619999,3.30080009,3.23090005, &
   &  3.95790005,3.86190009,3.66249990,3.52679992,3.36619997, &
   &  3.20959997,4.61759996,4.47639990,4.21960020,4.05970001, &
   &  3.85960007,3.75430012,3.56900001,3.46230006,3.39750004, &
   &  3.35249996,3.33080006,3.46199989,4.26230001,4.18739986, &
   &  4.01499987,3.89010000,3.73799992,3.58890009,5.05670023, &
   &  5.18139982,4.62610006,4.62010002,4.57019997,4.52710009, &
   &  4.48960018,4.45149994,4.42339993,4.12430000,4.24270010, &
   &  4.15409994,4.27939987,4.24499989,4.22079992,4.19859982, &
   &  4.01300001,4.24499989,4.09800005,3.98550010,3.89549994, &
   &  3.74900007,3.44560003,3.35249996,3.25640011,3.35990000, &
   &  4.31269979,4.27640009,4.11749983,4.00540018,3.86439991, &
   &  3.72160006,5.07959986,4.92939997,4.70429993,4.42519999, &
   &  4.45940018,4.39569998,4.35389996,4.43410015/)
  volsum = 0.0_wp
  do i = 1,n
    r = rcov(at(i))
    volsum = volsum + (4.0_wp / 3.0_wp) * (3.14159265359_wp * r**3)
  end do
  deallocate (rcov)
  return
end function volsum

!---- move molecule B above molecule A by a the vector X
subroutine shiftmol(A,nA,B,nB,X)
  use iso_fortran_env,wp => real64
  implicit none
  integer,intent(in)     :: nA
  real(wp),intent(inout) :: A(3,nA)
  integer,intent(in)     :: nB
  real(wp),intent(inout) :: B(3,nB)
  real(wp),intent(in)    :: X(3)
  integer :: i
  do i = 1,nB
    B(1:3,i) = B(1:3,i) + X(1:3)
  end do
  return
end subroutine shiftmol
!--- transform unit vector X to include coodrinates of mol boxes
subroutine shiftXvec(boxA,boxB,X,c)
  use iso_fortran_env,wp => real64
  implicit none
  real(wp),intent(inout)    :: X(3)
  real(wp),intent(in)       :: c
  real(wp) :: boxA(3,3)
  real(wp) :: boxB(3,3)
  X(1) = X(1) * (boxA(1,1) + abs(boxB(1,2)) + c)
  X(2) = X(2) * (boxA(2,1) + abs(boxB(2,2)) + c)
  X(3) = X(3) * (boxA(3,1) + abs(boxB(3,2)) + c)
  return
end subroutine shiftXvec
