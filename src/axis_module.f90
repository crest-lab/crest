!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2017 Stefan Grimme
! Copyright (C) 2021 Philipp Pracht
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

module axis_module
  use iso_fortran_env,only:wp => real64
  use atmasses,only:ams
  implicit none

  real(wp),parameter :: icm2MHz = 2.9979245d+4     !> cm⁻¹ to MHz
  real(wp),parameter :: MHz2icm = 1.0_wp / icm2MHz   !> MHz to cm⁻¹
  !************************************************************************
  !*   Aamu2icm   = conversion factor from Angstrom-amu to cm⁻¹
  !*              = (planck's constant*n*10**16)/(8*pi*pi*c)
  !*              = 6.62618*10**(-27)[erg-sec]*6.02205*10**23*10**16/
  !*                (8*(3.1415926535)**2*2.997925*10**10[cm/sec])
  !************************************************************************
  real(wp),parameter :: Aamu2icm = 16.8576522_wp

  !> 2π/3
  real(wp), parameter :: twothirdpi = 8.0_wp * atan(1.0_wp) / 3.0_wp

  public :: axis
  interface axis
    module procedure axis_0
    !> ARGS: nat,at,coord,rot,avmom,evec
    !> calculate rotational constants (rot) in MHz,
    !> the av. momentum (avmom in a.u.) and the trafo matrix (evec)
    module procedure axis_1
    !> ARGS: nat,at,coord,rot,avmom
    !> as axis_0, but omit evec
    module procedure axis_2
    !> ARGS: pr,nat,at,coord,eax
    !> calculate ellipsoid axes (eax) from rot constats, somehow, idk
    module procedure axis_3
    !> ARGS: nat,at,coord,coordout,rot
    !> as axis_0, but output only rot and write
    !> transformed, i.e., CMA shifted and rot-aligned coordinates
    !> to the output coordout
    module procedure axis_4
    !> ARGS: nat,at,coord
    !> as axis_3, but overwirte coord and
    !> doesn't output anything else
  end interface axis
  !> argument types:
  !> nat            -> integer
  !> at             -> integer,dimentsion(nat)
  !> coord,coordout -> real(wp),dimension(3,nat)
  !> rot, eax       -> real(wp),dimension(3)
  !> evec           -> real(wp),dimension(3,3)
  !> pr             -> logical

  public :: axistrf

  public :: cma
  interface cma
    module procedure CMAxyz
    module procedure CMAv
  end interface cma

contains
!========================================================================================!
!> subroutine axis_0
!> This is the original axis routine for calculating the
!> rotational constants of a molecule
!>
!> Input:    nat - number of atoms
!>            at - atom types
!>         coord - atomic coordinates in ANGSTROEM
!>
!> Output:   rot - rotational constants in MHz
!>         avmom - average momentum in a.u. (10⁻⁴⁷kg m²)
!>          evec - rot. matrix
!>
!========================================================================================!
  subroutine axis_0(nat,at,coord,rot,avmom,evec)
    implicit none
    integer :: nat
    integer :: at(nat)
    real(wp) :: coord(3,nat)
    real(wp) :: rot(3),avmom,evec(3,3)
    real(wp) :: a(3,3)
    real(wp) :: t(6),xyzmom(3),eig(3)
    !real(wp) :: x(nat),y(nat),z(nat)
    real(wp),allocatable :: x(:),y(:),z(:)
    real(wp) :: atmass
    integer :: i,j
    !************************************************************************
    !*     const1 =  10**40/(n*a*a)
    !*               n = avergadro's number
    !*               a = cm in an angstrom
    !*               10**40 is to allow units to be 10**(-40)gram-cm**2
    !*
    !************************************************************************
    real(wp),parameter :: const1 = 1.66053_wp
    !> first we move the molecule to the CMA
    !> this depends on the isotopic masses, and the cartesian geometry.
    !>
    allocate (x(nat),y(nat),z(nat),source=0.0_wp)
    call CMA(nat,at,coord,x,y,z)

    !************************************************************************
    !*    matrix for moments of inertia is of form
    !*
    !*   │ y²+z²                      │
    !*   │ -y*x       z²+x²           │ -i =0
    !*   │ -z*x       -z*y      x²+y² │
    !*
    !************************************************************************
    t = 0.0_wp
    do i = 1,6
      t(i) = dble(i) * 1.0d-10
    end do
    do i = 1,nat
      atmass = ams(at(i))
      t(1) = t(1) + atmass * (y(i)**2 + z(i)**2)
      t(2) = t(2) - atmass * x(i) * y(i)
      t(3) = t(3) + atmass * (z(i)**2 + x(i)**2)
      t(4) = t(4) - atmass * z(i) * x(i)
      t(5) = t(5) - atmass * y(i) * z(i)
      t(6) = t(6) + atmass * (x(i)**2 + y(i)**2)
      a(1,1) = t(1)
      a(2,1) = t(2)
      a(1,2) = t(2)
      a(2,2) = t(3)
      a(3,1) = t(4)
      a(1,3) = t(4)
      a(3,2) = t(5)
      a(2,3) = t(5)
      a(3,3) = t(6) 
    end do
    deallocate (z,y,x)

    evec = 0.0_wp
    eig = 0.0_wp
    !>--- old, numerical: rsp
    ! call rsp(t,3,3,eig,evec)
    !>--- new, analytical: eigvec3x3
    call eigvec3x3(a, eig, evec)

    do i = 1,3
      do j = 1,3
        if (abs(evec(i,j)) .lt. 1d-9) evec(i,j) = 0.0_wp
      end do
      if (eig(i) .lt. 3.d-4) then
        eig(i) = 0.d0
        rot(i) = 0.d0
      else
        rot(i) = icm2MHz * Aamu2icm / eig(i)
      end if
      xyzmom(i) = eig(i) * const1
    end do
    avmom = 1.d-47 * (xyzmom(1) + xyzmom(2) + xyzmom(3)) / 3.0_wp

    return
  end subroutine axis_0

!========================================================================================!
!> subroutine axis_1
!> format of axis_0 consistent with the original axis routine
!>-------------------------------------------
  subroutine axis_1(nat,at,coord,rot,avmom)
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: coord(3,nat)
    real(wp),intent(out) :: rot(3),avmom
    real(wp) :: dum(3,3)
    call axis_0(nat,at,coord,rot,avmom,dum)
    return
  end subroutine axis_1

!========================================================================================!
!> subroutine axis_2
!> format of axis consistent with the axis2 routine from the crest code.
!> I'm actually not sure whats going on there
!>---------------------------------------
  subroutine axis_2(pr,nat,at,coord,eax)
    implicit none
    logical,intent(in) :: pr
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(inout) :: coord(3,nat)
    real(wp),intent(out) :: eax(3)
    real(wp) :: rot(3),eig(3)
    real(wp) :: eps
    real(wp) :: sumw
    integer :: i
    real(wp),allocatable :: xyztmp(:,:)

    !call axis_0(nat,at,coord,rot,avmom,dum)
    allocate (xyztmp(3,nat))
    call axis_3(nat,at,coord,xyztmp,rot)
    coord = xyztmp
    deallocate (xyztmp)
    !> recover eig(3) from rot(3)
    !> this is needed because axis_0 outputs rot in MHz
    do i = 1,3
      if (rot(i) < 1.d-5) then
        eig(i) = 0.0_wp
      else
        eig(i) = icm2MHz * Aamu2icm / rot(i)
      end if
    end do

    eps = 1.d-9
    eig = 1.0_wp / (eig + eps)**0.25_wp !> ??? no idea
    sumw = sum(eig)
    eax = eig / sumw
    if (pr) then
      if (eax(2) .lt. 0.1_wp .and. eax(3) .lt. 0.1_wp) then
        eax(2:3) = 0.4_wp
        write (*,'(7x,''adjusting axis to finite length'')')
      end if
    end if
    sumw = sum(eax)
    eax = eax / sumw
    if (pr) write (*,'(7x,''unit ellipsoid axis a,b,c     :'',3f8.3)') eax

    return
  end subroutine axis_2

!========================================================================================!
!> subroutine axis_3
!> axis routine that orients the molecule along the
!> calculated principle axes and shifts it to CMA.
!> new geometry is written to coordout.
!>---------------------------------------------
  subroutine axis_3(nat,at,coord,coordout,rot)
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: coord(3,nat)
    real(wp),intent(out) :: coordout(3,nat)
    real(wp),intent(out) :: rot(3)
    real(wp) :: evec(3,3),avmom
    real(wp) :: xsum
    real(wp),allocatable :: coordtmp(:,:)
    real(wp),allocatable :: x(:),y(:),z(:)
    integer :: i,j,k

    !> call axis routine
    call axis_0(nat,at,coord,rot,avmom,evec)

    !> shift to CMA
    allocate (coordtmp(3,nat),x(nat),y(nat),z(nat))
    call CMA(nat,at,coord,x,y,z)
    coordtmp(1,:) = x
    coordtmp(2,:) = y
    coordtmp(3,:) = z

    !> do the trafo (chirality is preserved)
    xsum = calcxsum(evec)
    if (xsum .lt. 0.0_wp) then
      do j = 1,3
        evec(j,1) = -evec(j,1)
      end do
    end if

    do i = 1,nat
      do j = 1,3
        xsum = 0.0_wp
        do k = 1,3
          xsum = xsum + coordtmp(k,i) * evec(k,j)
        end do
        coordout(j,i) = xsum
      end do
    end do

    deallocate (z,y,x,coordtmp)

    return
  end subroutine axis_3

!========================================================================================!
!> subroutine axis_4
!> axis routine that orients the molecule along the
!> calculated principle axes and shifts it to CMA.
!> new geometry OVERWRITES input.
!>--------------------------------
  subroutine axis_4(nat,at,coord)
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(inout) :: coord(3,nat)
    real(wp) :: rot(3)
    real(wp),allocatable :: coordtmp(:,:)

    allocate (coordtmp(3,nat))
    !> call axis routine
    call axis_3(nat,at,coord,coordtmp,rot)
    coord = coordtmp
    deallocate (coordtmp)

    return
  end subroutine axis_4

!========================================================================================!
!> subroutine axistrf
!> as axis_3, but only the first nat0 atoms are taking for the
!> trafo and CMA shift. input coords are overwritten
!>--------------------------------------
  subroutine axistrf(nat,nat0,at,coord)
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: nat0
    integer,intent(in) :: at(nat)
    real(wp),intent(inout) :: coord(3,nat)
    real(wp) :: rot(3),evec(3,3)
    real(wp) :: xsum,avmom
    real(wp),allocatable :: coordtmp(:,:)
    integer,allocatable :: attmp(:)
    integer :: i,j,k

    !> call axis routine, only with the initial nat0 atoms
    allocate (attmp(nat0))
    allocate (coordtmp(3,nat0))
    attmp(1:nat0) = at(1:nat0)
    coordtmp(3,1:nat0) = coord(3,1:nat0)
    call axis_0(nat0,attmp,coordtmp,rot,avmom,evec)
    deallocate (coordtmp,attmp)

    !> shift to CMA of first nat0 atoms
    allocate (coordtmp(3,nat))
    coordtmp = coord
    call CMAtrf(nat,nat0,at,coordtmp)

    !> do the trafo
    xsum = calcxsum(evec)
    if (xsum .lt. 0.0_wp) then
      do j = 1,3
        evec(j,1) = -evec(j,1)
      end do
    end if
    do i = 1,nat
      do j = 1,3
        xsum = 0.0_wp
        do k = 1,3
          xsum = xsum + coordtmp(k,i) * evec(k,j)
        end do
        coord(j,i) = xsum
      end do
    end do
    deallocate (coordtmp)
    return
  end subroutine axistrf

!========================================================================================!
!> function calcxsum
!> calculates the determinant of evec(3,3)
!>---------------------------------
  real(wp) function calcxsum(evec)
    real(wp),intent(in) :: evec(3,3)
    calcxsum = evec(1,1) * (evec(2,2) * evec(3,3) - evec(3,2) * evec(2,3)) + &
  & evec(1,2) * (evec(2,3) * evec(3,1) - evec(2,1) * evec(3,3)) + &
  & evec(1,3) * (evec(2,1) * evec(3,2) - evec(2,2) * evec(3,1))
    return
  end function calcxsum

!========================================================================================!
!> subroutine CMA
!> calculate CMA-shifted coordinates x y z
!>--------------------------------------
  subroutine CMAxyz(nat,at,coord,x,y,z)
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: coord(3,nat)
    real(wp),intent(out) :: x(nat),y(nat),z(nat)
    integer :: i
    real(wp) :: sumw,sumwx,sumwy,sumwz,atmass
    sumw = 1.d-20
    sumwx = 0.d0
    sumwy = 0.d0
    sumwz = 0.d0
    do i = 1,nat
      atmass = ams(at(i))
      sumw = sumw + atmass
      sumwx = sumwx + atmass * coord(1,i)
      sumwy = sumwy + atmass * coord(2,i)
      sumwz = sumwz + atmass * coord(3,i)
    end do
    sumwx = sumwx / sumw
    sumwy = sumwy / sumw
    sumwz = sumwz / sumw
    do i = 1,nat
      x(i) = coord(1,i) - sumwx
      y(i) = coord(2,i) - sumwy
      z(i) = coord(3,i) - sumwz
    end do
    return
  end subroutine CMAxyz

!========================================================================================!
!> subroutine CMAtrf
!> calculate a shift to the first nat0 atoms' CMA
!>---------------------------------------
  subroutine CMAtrf(nat,nat0,at,coord)
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: nat0
    integer,intent(in) :: at(nat)
    real(wp),intent(inout) :: coord(3,nat)
    integer :: i
    real(wp) :: sumw,sumwx,sumwy,sumwz,atmass
    sumw = 1.d-20
    sumwx = 0.d0
    sumwy = 0.d0
    sumwz = 0.d0
    do i = 1,nat0  !> only these atoms!
      atmass = ams(at(i))
      sumw = sumw + atmass
      sumwx = sumwx + atmass * coord(1,i)
      sumwy = sumwy + atmass * coord(2,i)
      sumwz = sumwz + atmass * coord(3,i)
    end do
    sumwx = sumwx / sumw
    sumwy = sumwy / sumw
    sumwz = sumwz / sumw
    do i = 1,nat  !> overwrite input!
      coord(1,i) = coord(1,i) - sumwx
      coord(2,i) = coord(2,i) - sumwy
      coord(3,i) = coord(3,i) - sumwz
    end do
    return
  end subroutine CMAtrf

!========================================================================================!
!> subroutine CMAv
!> calculate a CMA coordinats and save them to vec
!>----------------------------------
  subroutine CMAv(nat,at,coord,vec)
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: coord(3,nat)
    real(wp),intent(out) :: vec(3)
    integer :: i
    real(wp) :: sumw,sumwx,sumwy,sumwz,atmass
    sumw = 1.d-20
    sumwx = 0.d0
    sumwy = 0.d0
    sumwz = 0.d0
    do i = 1,nat
      atmass = ams(at(i))
      sumw = sumw + atmass
      sumwx = sumwx + atmass * coord(1,i)
      sumwy = sumwy + atmass * coord(2,i)
      sumwz = sumwz + atmass * coord(3,i)
    end do
    vec(1) = sumwx / sumw
    vec(2) = sumwy / sumw
    vec(3) = sumwz / sumw
    return
  end subroutine CMAv

!========================================================================================!
!****************************************************************************************!
! The following routines were translated from older F70 code in drsp.f
!****************************************************************************************!
!========================================================================================!

!*******************************************************************
!*
!*   eispack diagonalization routines: to find the eigenvalues and
!*           eigenvectors (if desired) of a real symmetric packed matrix.
!* on input-      n  is the order of the matrix  a,
!*                a  contains the lower triangle of the real symmetric
!*                   packed matrix stored row-wise,
!*             matz  is an integer variable set equal to zero if only
!*                   eigenvalues are desired,  otherwise it is set to
!*                   any non-zero integer for both eigenvalues and
!*                   eigenvectors.
!* on output-     w  contains the eigenvalues in ascending order,
!*                z  contains the eigenvectors if matz is not zero,
!*
!*******************************************************************
  subroutine rsp(a,n,matz,w,z)
    implicit none
    integer :: n,matz
    real(wp) :: a(n * (n + 1) / 2)
    real(wp) :: w(n),z(n,n)
    real(wp) :: fv1(2 * n),fv2(2 * n)
    real(wp) :: eps,eta
    integer :: nv,nm,ierr
    integer :: i

    if (n .eq. 1) then
      z(1,1) = 1.0d0
      w(1) = a(1)
      return
    end if

    call epseta(eps,eta)

    nv = (n * (n + 1)) / 2
    nm = n
    call tred3(n,a,w,fv1,fv2)
    if (matz .ne. 0) then ! go to 10
      !>--- find eigenvalues only
      z = 0.0_wp
      do i = 1,n
        z(i,i) = 1.0_wp
      end do
      call tql2(nm,n,w,fv1,z,ierr,eps)
      if (ierr .ne. 0) return
      call trbak3(nm,n,a,n,z)
    else
      call tqlrat(n,w,fv2,ierr,eps)
    end if
    return
  end subroutine rsp

!*******************************************************************
!*  compute and return eta, the smallest representable number,
!*  and eps is the smallest number for which 1+eps.ne.1.
!*******************************************************************
  subroutine epseta(eps,eta)
    implicit none
    real(wp) :: eps,eta
    !> epsilon(X) returns the smallest number E
    !> of the same kind as X such that 1 + E > 1.
    !> I.e., it does exactly what this routine did
    intrinsic :: epsilon
    !eta = 1.0_wp
    !do
    !  if ((eta / 2.0_wp) .eq. 0.0_wp) exit
    !  if (eta .lt. 1.d-38) exit
    !  eta = eta / 2.0_wp
    !end do
    eta = epsilon(eta - 1.0_wp)
    !eps = 1.d0
    !do
    !  if ((1.0_wp + (eps / 2.0_wp)) .eq. 1.0_wp) exit
    !  if (eps .lt. 1.d-17) exit
    !  eps = eps / 2.0_wp
    !end do
    eps = epsilon(eps)
    return
  end subroutine epseta

!*******************************************************************
!*     this subroutine finds the eigenvalues and eigenvectors
!*     of a symmetric tridiagonal matrix by the ql method.
!*     the eigenvectors of a full symmetric matrix can also
!*     be found if  tred2  has been used to reduce this
!*     full matrix to tridiagonal form.
!*
!*     on input-
!*
!*        nm must be set to the row dimension of two-dimensional
!*          array parameters as declared in the calling program
!*          dimension statement,
!*
!*        n is the order of the matrix,
!*
!*        d contains the diagonal elements of the input matrix,
!*
!*        e contains the subdiagonal elements of the input matrix
!*          in its last n-1 positions.  e(1) is arbitrary,
!*
!*        z contains the transformation matrix produced in the
!*          reduction by  tred2, if performed.  if the eigenvectors
!*          of the tridiagonal matrix are desired, z must contain
!*          the identity matrix.
!*
!*      on output-
!*
!*        d contains the eigenvalues in ascending order.  if an
!*          error exit is made, the eigenvalues are correct but
!*          unordered for indices 1,2,...,ierr-1,
!*
!*        e has been destroyed,
!*
!*        z contains orthonormal eigenvectors of the symmetric
!*          tridiagonal (or full) matrix.  if an error exit is made,
!*          z contains the eigenvectors associated with the stored
!*          eigenvalues,
!*
!*        ierr is set to
!*          zero       for normal return,
!*          j          if the j-th eigenvalue has not been
!*                     determined after 30 iterations.
!*******************************************************************
  subroutine tql2(nm,n,d,e,z,ierr,eps)
    implicit none
    integer :: nm,n
    real(wp) :: d(*)
    real(wp) :: e(*)
    real(wp) :: z(nm,*)
    integer :: ierr
    real(wp) :: eps
    real(wp) :: f,b,h,g,p,r,c,s
    integer :: i,j,l,ii,m,l1,mml,k

    ierr = 0
    if (n .eq. 1) return

    do i = 2,n
      e(i - 1) = e(i)
    end do
    f = 0.0_wp
    b = 0.0_wp
    e(n) = 0.0_wp

    LLOOP: do l = 1,n
      j = 0
      h = eps * (abs(d(l)) + abs(e(l)))
      if (b .lt. h) b = h
!>    ********** look for small sub-diagonal element **********
      do m = l,n
        if (abs(e(m)) .le. b) exit !go to 30
!>    ********** e(n) is always zero, so there is no exit
!>               through the bottom of the loop **********
      end do

      JLOOP: do
        if (m .eq. l) exit JLOOP !go to 100
        if (j .eq. 30) then
          ierr = l
          return
        end if
        j = j + 1
!>    ********** form shift **********
        l1 = l + 1
        g = d(l)
        p = (d(l1) - g) / (2.0d0 * e(l))
        r = sqrt(p * p + 1.0d0)
        d(l) = e(l) / (p + sign(r,p))
        h = g - d(l)

        do i = l1,n
          d(i) = d(i) - h
        end do

        f = f + h
!>    ********** ql transformation **********
        p = d(m)
        c = 1.0d0
        s = 0.0d0
        mml = m - l
!>    ********** for i=m-1 step -1 until l do -- **********
        do ii = 1,mml
          i = m - ii
          g = c * e(i)
          h = c * p
          if (abs(p) .lt. abs(e(i))) then
            c = p / e(i)
            r = sqrt(c * c + 1.0d0)
            e(i + 1) = s * e(i) * r
            s = 1.0d0 / r
            c = c * s
          else
            c = e(i) / p
            r = sqrt(c * c + 1.0d0)
            e(i + 1) = s * p * r
            s = c / r
            c = 1.0d0 / r
          end if
          p = c * d(i) - s * g
          d(i + 1) = h + s * (c * g + s * d(i))
!>    ********** form vector **********
          do k = 1,n
            h = z(k,i + 1)
            z(k,i + 1) = s * z(k,i) + c * h
            z(k,i) = c * z(k,i) - s * h
          end do
        end do

        e(l) = s * p
        d(l) = c * p
        if (abs(e(l)) .gt. b) cycle JLOOP ! go to 40
        exit JLOOP
      end do JLOOP
      d(l) = d(l) + f
    end do LLOOP
!>    ********** order eigenvalues and eigenvectors **********
    do ii = 2,n
      i = ii - 1
      k = i
      p = d(i)

      do j = ii,n
        if (d(j) .ge. p) exit !go to 120
        k = j
        p = d(j)
      end do

      if (k .eq. i) return
      d(k) = d(i)
      d(i) = p

      do j = 1,n
        p = z(j,i)
        z(j,i) = z(j,k)
        z(j,k) = p
      end do
    end do

    return
  end subroutine tql2

!*******************************************************************
!*     this subroutine finds the eigenvalues of a symmetric
!*     tridiagonal matrix by the rational ql method.
!*
!*     on input-
!*
!*        n is the order of the matrix,
!*
!*        d contains the diagonal elements of the input matrix,
!*
!*        e2 contains the squares of the subdiagonal elements of the
!*          input matrix in its last n-1 positions.  e2(1) is arbitrary.
!*
!*      on output-
!*
!*        d contains the eigenvalues in ascending order.  if an
!*          error exit is made, the eigenvalues are correct and
!*          ordered for indices 1,2,...ierr-1, but may not be
!*          the smallest eigenvalues,
!*
!*        e2 has been destroyed,
!*
!*        ierr is set to
!*          zero       for normal return,
!*          j          if the j-th eigenvalue has not been
!*                     determined after 30 iterations.
!*
!*******************************************************************
  subroutine tqlrat(n,d,e2,ierr,eps)
    implicit none
    integer :: n
    real(wp) :: d(*)
    real(wp) :: e2(*)
    integer :: ierr
    real(wp) :: eps
    integer :: i,j,l,l1,m,ii,mml
    real(wp) :: f,b,h,c,r,s,g,p

    ierr = 0
    if (n .eq. 1) return

    do i = 2,n
      e2(i - 1) = e2(i)
    end do

    f = 0.0_wp
    b = 0.0_wp
    c = 0.0_wp
    e2(n) = 0.0_wp

    LLOOP: do l = 1,n
      j = 0
      h = eps * (abs(d(l)) + sqrt(e2(l)))
      if (b <= h) then !go to 20
        b = h
        c = b * b
      end if
!>    ********** look for small squared sub-diagonal element **********
      do m = l,n
        if (e2(m) .le. c) exit !go to 40
!>    ********** e2(n) is always zero, so there is no exit
!>               through the bottom of the loop **********
      end do

      JLOOP: do
        if (m .eq. l) exit JLOOP !go to 80
        if (j .eq. 30) then
          ierr = l !go to 130
          return
        end if
        j = j + 1
!>    ********** form shift **********
        l1 = l + 1
        s = sqrt(e2(l))
        g = d(l)
        p = (d(l1) - g) / (2.0_wp * s)
        r = sqrt(p * p + 1.0_wp)
        d(l) = s / (p + sign(r,p))
        h = g - d(l)

        do i = l1,n
          d(i) = d(i) - h
        end do

        f = f + h
!>    ********** rational ql transformation **********
        g = d(m)
        if (g .eq. 0.0_wp) g = b
        h = g
        s = 0.0_wp
        mml = m - l
!>    ********** for i=m-1 step -1 until l do -- **********
        do ii = 1,mml
          i = m - ii
          p = g * h
          r = p + e2(i)
          e2(i + 1) = s * r
          s = e2(i) / r
          d(i + 1) = h + s * (h + d(i))
          g = d(i) - e2(i) / g
          if (g .eq. 0.0_wp) g = b
          h = g * p / r
        end do
        e2(l) = s * g
        d(l) = h
!>    ********** guard against underflow in convergence test **********
        if (h .eq. 0.0_wp) exit JLOOP !go to 80
        if (abs(e2(l)) .le. abs(c / h)) exit JLOOP !go to 80
        e2(l) = h * e2(l)
        if (e2(l) .ne. 0.0_wp) cycle JLOOP !go to 50
        exit JLOOP
      end do JLOOP
      p = d(l) + f
!>    ********** order eigenvalues **********
      if (l /= 1) then !if (l .eq. 1) go to 100
!>    ********** for i=l step -1 until 2 do -- **********
        do ii = 2,l
          i = l + 2 - ii
          if (p .ge. d(i - 1)) then !go to 110
            d(i) = p
            cycle LLOOP
          end if
          d(i) = d(i - 1)
        end do
      end if
      i = 1
      d(i) = p
    end do LLOOP

    return
  end subroutine tqlrat

!*******************************************************************
!*     this subroutine forms the eigenvectors of a real symmetric
!*     matrix by back transforming those of the corresponding
!*     symmetric tridiagonal matrix determined by  tred3.
!*
!*     on input-
!*
!*        nm must be set to the row dimension of two-dimensional
!*          array parameters as declared in the calling program
!*          dimension statement,
!*
!*        n is the order of the matrix,
!*
!*        nv must be set to the dimension of the array parameter a
!*          as declared in the calling program dimension statement,
!*
!*        a contains information about the orthogonal transformations
!*          used in the reduction by  tred3  in its first
!*          n*(n+1)/2 positions,
!*
!*        m is the number of eigenvectors to be back transformed,
!*
!*        z contains the eigenvectors to be back transformed
!*          in its first m columns.
!*
!*     on output-
!*
!*        z contains the transformed eigenvectors
!*          in its first m columns.
!*
!*     note that trbak3 preserves vector euclidean norms.
!*
!*******************************************************************
  subroutine trbak3(nm,n,a,m,z)
    implicit none
    integer :: nm,n,m
    real(wp) :: a(*),z(nm,*)
    integer :: i,iz,ik,j,k,l
    real(wp) :: h,s

    if (m .eq. 0) return
    if (n .eq. 1) return

    do i = 2,n
      l = i - 1
      iz = (i * l) / 2
      ik = iz + i
      h = a(ik)
      if (h .eq. 0.0d0) exit
      do j = 1,m
        s = 0.0d0
        ik = iz
        do k = 1,l
          ik = ik + 1
          s = s + a(ik) * z(k,j)
        end do
!>    ********** double division avoids possible underflow **********
        s = (s / h) / h
        ik = iz
        do k = 1,l
          ik = ik + 1
          z(k,j) = z(k,j) - s * a(ik)
        end do
      end do
    end do

    return
  end subroutine trbak3

!*******************************************************************
!*     this subroutine reduces a real symmetric matrix, stored as
!*     a one-dimensional array, to a symmetric tridiagonal matrix
!*     using orthogonal similarity transformations.
!*
!*     on input-
!*
!*        n is the order of the matrix,
!*
!*        nv must be set to the dimension of the array parameter a
!*          as declared in the calling program dimension statement,
!*
!*        a contains the lower triangle of the real symmetric
!*          input matrix, stored row-wise as a one-dimensional
!*          array, in its first n*(n+1)/2 positions.
!*
!*     on output-
!*
!*        a contains information about the orthogonal
!*          transformations used in the reduction,
!*
!*        d contains the diagonal elements of the tridiagonal matrix,
!*
!*        e contains the subdiagonal elements of the tridiagonal
!*          matrix in its last n-1 positions.  e(1) is set to zero,
!*
!*        e2 contains the squares of the corresponding elements of e.
!*          e2 may coincide with e if the squares are not needed.
!*
!*******************************************************************
  subroutine tred3(n,a,d,e,e2)
    implicit none
    integer :: n
    real(wp) :: a(*),d(*),e(*),e2(*)
    integer :: i,j,k,l,iz,jk,ii
    real(wp) :: scale,h,hh,f,g
    !implicit double precision(a - h,o - z)
    !dimension a(*),d(*),e(*),e2(*)

    do ii = 1,n
      i = n + 1 - ii
      l = i - 1
      iz = (i * l) / 2
      h = 0.0d0
      scale = 0.0d0
      do k = 1,l
        iz = iz + 1
        d(k) = a(iz)
        scale = scale + abs(d(k))
      end do

      if (scale .ne. 0.d0) then !go to 20
        do k = 1,l
          d(k) = d(k) / scale
          h = h + d(k) * d(k)
        end do

        e2(i) = scale * scale * h
        f = d(l)
        g = -sign(sqrt(h),f)
        e(i) = scale * g
        h = h - f * g
        d(l) = f - g
        a(iz) = scale * d(l)
        if (l .eq. 1) go to 90
        f = 0.0d0

        do j = 1,l
          g = 0.0d0
          jk = (j * (j - 1)) / 2
!>    ********** form element of a*u **********
          k = 0
          do
            k = k + 1
            jk = jk + 1
            g = g + a(jk) * d(k)
            if (k .lt. j) cycle !go to 40
            exit
          end do
          if (k /= l) then
            do
              jk = jk + k
              k = k + 1
              g = g + a(jk) * d(k)
              if (k .lt. l) cycle !go to 50
              exit
            end do
!>    ********** form element of p **********
          end if
          e(j) = g / h
          f = f + e(j) * d(j)
        end do

        hh = f / (h + h)
        jk = 0
!>    ********** form reduced a **********
        do j = 1,l
          f = d(j)
          g = e(j) - hh * f
          e(j) = g

          do k = 1,j
            jk = jk + 1
            a(jk) = a(jk) - f * e(k) - g * d(k)
          end do
        end do

      else
        e(i) = 0.0d0
        e2(i) = 0.0d0
      end if

90    d(i) = a(iz + 1)
      a(iz + 1) = scale * sqrt(h)
    end do

    return
  end subroutine tred3


!========================================================================================!
!****************************************************************************************!
!> Analytical 3x3 eigenvalue solver to replace rsp routine
!> Taken from github.com/awvwgk/diag3x3
!****************************************************************************************!
!========================================================================================!
!> Calculates eigenvalues based on the trigonometric solution of A = pB + qI
pure subroutine eigval3x3(a, w)

   !> The symmetric input matrix
   real(wp), intent(in) :: a(3, 3)

   !> Contains eigenvalues on exit
   real(wp), intent(out) :: w(3)

   real(wp) :: q, p, r

   r = a(1, 2) * a(1, 2) + a(1, 3) * a(1, 3) + a(2, 3) * a(2, 3)
   q = (a(1, 1) + a(2, 2) + a(3, 3)) / 3.0_wp
   w(1) = a(1, 1) - q
   w(2) = a(2, 2) - q
   w(3) = a(3, 3) - q
   p = sqrt((w(1) * w(1) + w(2) * w(2) + w(3) * w(3) + 2*r) / 6.0_wp)
   r = (w(1) * (w(2) * w(3) - a(2, 3) * a(2, 3)) &
      & - a(1, 2) * (a(1, 2) * w(3) - a(2, 3) * a(1, 3)) &
      & + a(1, 3) * (a(1, 2) * a(2, 3) - w(2) * a(1, 3))) / (p*p*p) * 0.5_wp

   if (r <= -1.0_wp) then
      r = 0.5_wp * twothirdpi
   else if (r >= 1.0_wp) then
      r = 0.0_wp
   else
      r = acos(r) / 3.0_wp
   end if

   w(3) = q + 2 * p * cos(r)
   w(1) = q + 2 * p * cos(r + twothirdpi)
   w(2) = 3 * q - w(1) - w(3)

end subroutine eigval3x3


!> Calculates eigenvector using an analytical method based on vector cross
!  products.
pure subroutine eigvec3x3(a, w, q)

   !> The symmetric input matrix, destroyed while solving
   real(wp), intent(inout) :: a(3,3)

   !> Contains eigenvalues on exit
   real(wp), intent(out) :: w(3)

   !> Contains eigenvectors on exit
   real(wp), intent(out) :: q(3,3)

   !> Numerical precision
   real(wp), parameter :: eps = epsilon(1.0_wp)

   !> Local variables
   real(wp) :: norm, n1, n2, n3, precon
   integer :: i

   w(1) = max(abs(a(1, 1)), abs(a(1, 2)))
   w(2) = max(abs(a(1, 3)), abs(a(2, 2)))
   w(3) = max(abs(a(2, 3)), abs(a(3, 3)))
   precon = max(w(1), max(w(2), w(3)))

   ! null matrix
   if (precon < eps) then
      w(1) = 0.0_wp
      w(2) = 0.0_wp
      w(3) = 0.0_wp
      q(1, 1) = 1.0_wp
      q(2, 2) = 1.0_wp
      q(3, 3) = 1.0_wp
      q(1, 2) = 0.0_wp
      q(1, 3) = 0.0_wp
      q(2, 3) = 0.0_wp
      q(2, 1) = 0.0_wp
      q(3, 1) = 0.0_wp
      q(3, 2) = 0.0_wp
      return
   end if

   norm = 1.0_wp / precon

   a(1, 1) = a(1, 1) * norm
   a(1, 2) = a(1, 2) * norm
   a(2, 2) = a(2, 2) * norm
   a(1, 3) = a(1, 3) * norm
   a(2, 3) = a(2, 3) * norm
   a(3, 3) = a(3, 3) * norm

   ! Calculate eigenvalues
   call eigval3x3(a, w)

   ! Compute first eigenvector
   a(1, 1) = a(1, 1) - w(1)
   a(2, 2) = a(2, 2) - w(1)
   a(3, 3) = a(3, 3) - w(1)

   q(1, 1) = a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2)
   q(2, 1) = a(1, 3) * a(1, 2) - a(1, 1) * a(2, 3)
   q(3, 1) = a(1, 1) * a(2, 2) - a(1, 2) * a(1, 2)
   q(1, 2) = a(1, 2) * a(3, 3) - a(1, 3) * a(2, 3)
   q(2, 2) = a(1, 3) * a(1, 3) - a(1, 1) * a(3, 3)
   q(3, 2) = a(1, 1) * a(2, 3) - a(1, 2) * a(1, 3)
   q(1, 3) = a(2, 2) * a(3, 3) - a(2, 3) * a(2, 3)
   q(2, 3) = a(2, 3) * a(1, 3) - a(1, 2) * a(3, 3)
   q(3, 3) = a(1, 2) * a(2, 3) - a(2, 2) * a(1, 3)
   n1 = q(1, 1) * q(1, 1) + q(2, 1) * q(2, 1) + q(3, 1) * q(3, 1)
   n2 = q(1, 2) * q(1, 2) + q(2, 2) * q(2, 2) + q(3, 2) * q(3, 2)
   n3 = q(1, 3) * q(1, 3) + q(2, 3) * q(2, 3) + q(3, 3) * q(3, 3)

   norm = n1
   i = 1
   if (n2 > norm) then
      i = 2
      norm = n1
   end if
   if (n3 > norm) then
      i = 3
   end if

   if (i == 1) then
      norm = sqrt(1.0_wp / n1)
      q(1, 1) = q(1, 1) * norm
      q(2, 1) = q(2, 1) * norm
      q(3, 1) = q(3, 1) * norm
   else if (i == 2) then
      norm = sqrt(1.0_wp / n2)
      q(1, 1) = q(1, 2) * norm
      q(2, 1) = q(2, 2) * norm
      q(3, 1) = q(3, 2) * norm
   else
      norm = sqrt(1.0_wp / n3)
      q(1, 1) = q(1, 3) * norm
      q(2, 1) = q(2, 3) * norm
      q(3, 1) = q(3, 3) * norm
   end if

   ! Robustly compute a right-hand orthonormal set (ev1, u, v)
   if (abs(q(1, 1)) > abs(q(2, 1))) then
      norm = sqrt(1.0_wp / (q(1, 1) * q(1, 1) + q(3, 1) * q(3, 1)))
      q(1, 2) = -q(3, 1) * norm
      q(2, 2) = 0.0_wp
      q(3, 2) = +q(1, 1) * norm
   else
      norm = sqrt(1.0_wp / (q(2, 1) * q(2, 1) + q(3, 1) * q(3, 1)))
      q(1, 2) = 0.0_wp
      q(2, 2) = +q(3, 1) * norm
      q(3, 2) = -q(2, 1) * norm
   end if
   q(1, 3) = q(2, 1) * q(3, 2) - q(3, 1) * q(2, 2)
   q(2, 3) = q(3, 1) * q(1, 2) - q(1, 1) * q(3, 2)
   q(3, 3) = q(1, 1) * q(2, 2) - q(2, 1) * q(1, 2)

   ! Reset A
   a(1, 1) = a(1, 1) + w(1)
   a(2, 2) = a(2, 2) + w(1)
   a(3, 3) = a(3, 3) + w(1)

   ! A*U
   n1 = a(1, 1) * q(1, 2) + a(1, 2) * q(2, 2) + a(1, 3) * q(3, 2)
   n2 = a(1, 2) * q(1, 2) + a(2, 2) * q(2, 2) + a(2, 3) * q(3, 2)
   n3 = a(1, 3) * q(1, 2) + a(2, 3) * q(2, 2) + a(3, 3) * q(3, 2)

   ! A*V, note out of order computation
   a(3, 3) = a(1, 3) * q(1, 3) + a(2, 3) * q(2, 3) + a(3, 3) * q(3, 3)
   a(1, 3) = a(1, 1) * q(1, 3) + a(1, 2) * q(2, 3) + a(1, 3) * q(3, 3)
   a(2, 3) = a(1, 2) * q(1, 3) + a(2, 2) * q(2, 3) + a(2, 3) * q(3, 3)

   ! UT*(A*U) - l2*E
   n1 = q(1, 2) * n1 + q(2, 2) * n2 + q(3, 2) * n3 - w(2)
   ! UT*(A*V)
   n2 = q(1, 2) * a(1, 3) + q(2, 2) * a(2, 3) + q(3, 2) * a(3, 3)
   ! VT*(A*V) - l2*E
   n3 = q(1, 3) * a(1, 3) + q(2, 3) * a(2, 3) + q(3, 3) * a(3, 3) - w(2)

   if (abs(n1) >= abs(n3)) then
      norm = max(abs(n1), abs(n2))
      if (norm > eps) then
         if (abs(n1) >= abs(n2)) then
            n2 = n2 / n1
            n1 = sqrt(1.0_wp / (1.0_wp + n2 * n2))
            n2 = n2 * n1
         else
            n1 = n1 / n2
            n2 = sqrt(1.0_wp / (1.0_wp + n1 * n1))
            n1 = n1 * n2
         end if
         q(1, 2) = n2 * q(1, 2) - n1 * q(1, 3)
         q(2, 2) = n2 * q(2, 2) - n1 * q(2, 3)
         q(3, 2) = n2 * q(3, 2) - n1 * q(3, 3)
      end if
   else
      norm = max(abs(n3), abs(n2))
      if (norm > eps) then
         if (abs(n3) >= abs(n2)) then
            n2 = n2 / n3
            n3 = sqrt(1.0_wp / (1.0_wp + n2 * n2))
            n2 = n2 * n3
         else
            n3 = n3 / n2
            n2 = sqrt(1.0_wp / (1.0_wp + n3 * n3))
            n3 = n3 * n2
         end if
         q(1, 2) = n3 * q(1, 2) - n2 * q(1, 3)
         q(2, 2) = n3 * q(2, 2) - n2 * q(2, 3)
         q(3, 2) = n3 * q(3, 2) - n2 * q(3, 3)
      end if
   end if

   ! Calculate third eigenvector from cross product
   q(1, 3) = q(2, 1) * q(3, 2) - q(3, 1) * q(2, 2)
   q(2, 3) = q(3, 1) * q(1, 2) - q(1, 1) * q(3, 2)
   q(3, 3) = q(1, 1) * q(2, 2) - q(2, 1) * q(1, 2)

   w(1) = w(1) * precon
   w(2) = w(2) * precon
   w(3) = w(3) * precon

end subroutine eigvec3x3

!========================================================================================!
!========================================================================================!
end module axis_module
