!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2022 Philipp Pracht
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

module geo
  use iso_fortran_env,wp => real64
  implicit none
  public

  real(wp),parameter,private :: pi = acos(0.0_wp)*2.0_wp
  real(wp),parameter,private :: pi2 = pi*2.0_wp
 
!================================================================================!
!================================================================================!
contains  !> MODULE PROCEDURES START HERE
!================================================================================!
!================================================================================!
!> length of a vector
    pure function vec_len(v) result(l)
        implicit none
        real(wp), intent(in) :: v(3)
        real(wp) :: l
        l = sqrt(dot_product(v,v))
    end function
!================================================================================!
!> euclidean distance between two points
    pure function distance(p,q) result(l)
        implicit none
        real(wp), intent(in) :: p(3),q(3)
        real(wp) :: l,pq(3)
        pq(:) = p(:) - q(:)
        l = vec_len(pq)
    end function distance
!================================================================================!
!> normalize the vector k
  subroutine unitv(k)
    implicit none
    real(wp),intent(inout) :: k(3)
    k = k / vec_len(k)
    return
  end subroutine unitv
!================================================================================!
!> calculate angle theta between vectors u and v
  function tangle(u,v) result(angle)
    implicit none
    real(wp) :: angle
    real(wp),intent(in) :: u(3)
    real(wp),intent(in) :: v(3)
    real(wp) :: ulen
    real(wp) :: vlen
    real(wp) :: uv

    angle = 0.0_wp
    ulen = vec_len(u)
    vlen = vec_len(v)
    !> uv = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
    uv = dotp(u,v,3)
    !> acos returns a value between 0 and pi
    angle = acos(uv / (ulen * vlen))
    return
  end function tangle
!================================================================================!
!> calculate dot product between two vectors u and v
  function dotp(u,v,n) result(dotproduct)
    implicit none
    real(wp) :: dotproduct
    integer :: n     !dimension
    real(wp),intent(in) :: u(n) !first vector
    real(wp),intent(in) :: v(n) !second vector
    integer :: i
    dotproduct = 0.0_wp
    do i = 1,n
      dotproduct = dotproduct + u(i) * v(i)
    end do
    return
  end function dotp
!================================================================================!
!> calculate the cross product between two vectors u and v
  subroutine crosp(u,v,uv)
    implicit none
    real(wp),intent(in) :: u(3) !first vector
    real(wp),intent(in) :: v(3) !second vector
    real(wp),intent(out) :: uv(3)
    uv = 0.0_wp
    uv(1) = u(2) * v(3) - u(3) * v(2)
    uv(2) = u(3) * v(1) - u(1) * v(3)
    uv(3) = u(1) * v(2) - u(2) * v(1)
    return
  end subroutine crosp
!================================================================================!
!> calculate the normal between two vectors u and v spanning a plane
!> a third point (the origin) can optionally be provided
  function normal(u,v,o) result(nvec)
    implicit none
    real(wp) :: nvec(3)
    real(wp),intent(in) :: u(3) !> first vector
    real(wp),intent(in) :: v(3) !> second vector
    real(wp),intent(in),optional :: o(3) !> optional, origin
    real(wp) :: utmp(3), vtmp(3)
    nvec = 0.0_wp
    if(present(o))then
      utmp(1:3) = u(1:3) - o(1:3)
      vtmp(1:3) = v(1:3) - o(1:3)
    else
      utmp = u
      vtmp = v
    endif
    call crosp(utmp,vtmp,nvec)
    call unitv(nvec)
    return
  end function normal
!================================================================================!
!> calculate the dihedral angle between three vectors u, v, and w
!> spanning two half planes. 
!> For four given points A-B-C-D defining the dihedral angle between the planes
!> A-B-C and B-C-D, the vectors are u:A->B, v:B->C, and w:C->D
!> The result will be in radians between −π (-180°) and π (180°).
  function dihedral(u,v,w) result(dihed)
    implicit none
    real(wp) :: dihed
    real(wp),intent(in) :: u(3) !> first vector
    real(wp),intent(in) :: v(3) !> second vector
    real(wp),intent(in) :: w(3) !> third vector
    real(wp) :: uv(3), vw(3), uvw(3), tmp1, tmp2
    dihed = 0.0_wp
    call crosp(u, v, uv)
    call crosp(v, w, vw)
    call crosp(uv, vw, uvw) 
    tmp1 = dot_product(v, uvw) + 1.0d-10  !> small constant to avoid singularities
    tmp2 = vec_len(v) * dot_product(uv,vw)
    dihed = atan2(tmp1,tmp2) 
    return
  end function dihedral

!================================================================================!
!> shift an angle x into the interval −π (-180°) <= x < π (180°)
  function angleshift(x) result(xnew)
     implicit none  
      real(wp) :: x
      real(wp) :: xnew
      real(wp) :: f
      f = floor(0.5_wp + ( x / pi2 ))
      xnew = 2.0_wp*((x/pi2) - f)*pi
  end function angleshift


!================================================================================!
!================================================================================!
!> rotate around z axis
  subroutine rotRz(u,theta)
    implicit none
    real(wp) :: u(3) !first vector
    real(wp) :: v(3) !second vector
    real(wp) :: theta
    real(wp) :: Rz(3,3)

    Rz = 0.0_wp
    Rz = reshape(&
      &  [  cos(theta),  sin(theta),   0.0_wp,   &
      &    -sin(theta),  cos(theta),   0.0_wp,   &
      &         0.0_wp,      0.0_wp,   1.0_wp ], &
      &  [3,3])
    
    v = 0.0_wp
    v = matmul(Rz,u)

    u = v
    return
  end subroutine rotRz
!================================================================================!
  subroutine rotRz180(u)
    implicit none
    real(wp) :: u(3) !first vector
    real(wp) :: v(3) !second vector
    real(wp) :: Rz(3,3)

    Rz = 0.0_wp
    Rz(1,1) = -1.0_wp !cos(theta)
    Rz(2,1) = 0.0_wp  !sin(theta)
    Rz(1,2) = 0.0_wp  !-sin(theta)
    Rz(2,2) = -1.0_wp !cos(theta)
    Rz(3,3) = 1.0_wp

    v = 0.0_wp
    v = matmul(Rz,u)

    u = v
    return
  end subroutine rotRz180
!================================================================================!
!> rotate around y axis
  subroutine rotRy(u,theta)
    implicit none
    real(wp) :: u(3) !first vector
    real(wp) :: v(3) !second vector
    real(wp) :: theta
    real(wp) :: Ry(3,3)

    Ry = reshape(&
      &  [  cos(theta), 0.0_wp, sin(theta),   &
      &         0.0_wp, 1.0_wp,     0.0_wp,   &
      &    -sin(theta), 0.0_wp, cos(theta) ], &
      &  [3,3])

    v = 0.0_wp
    v = matmul(Ry,u)
    u = v
    return
  end subroutine rotRy
!================================================================================!
  subroutine rotRy180(u)
    implicit none
    real(wp) :: u(3) !first vector
    real(wp) :: v(3) !second vector
    real(wp) :: Ry(3,3)

    Ry = 0.0_wp
    Ry(1,1) = -1.0_wp !cos(theta)
    Ry(3,1) = 0.0_wp  !-sin(theta)
    Ry(3,1) = 0.0_wp  !sin(theta)
    Ry(3,3) = -1.0_wp !cos(theta)
    Ry(2,2) = 1.0_wp

    v = 0.0_wp
    v = matmul(Ry,u)
    u = v
    return
  end subroutine rotRy180
!================================================================================!
!> rotate around x axis
  subroutine rotRx(u,theta)
    implicit none
    real(wp) :: u(3) !first vector
    real(wp) :: v(3) !second vector
    real(wp) :: theta
    real(wp) :: Rx(3,3)

    Rx = reshape(&
      &  [ 1.0_wp,      0.0_wp,      0.0_wp,   &
      &    0.0_wp,  cos(theta),  sin(theta),   &
      &    0.0_wp, -sin(theta),  cos(theta) ], &
      &  [3,3])

    v = 0.0_wp
    v = matmul(Rx,u)
    u = v
    return
  end subroutine rotRx
!================================================================================!
!> rotate around x axis 180deg
  subroutine rotRx180(u)
    implicit none
    real(wp) :: u(3) !first vector
    real(wp) :: v(3) !second vector
    real(wp) :: Rx(3,3)

    Rx = 0.0_wp
    Rx(2,2) = -1.0_wp !cos(theta)
    Rx(2,3) = 0.0_wp  !-sin(theta)
    Rx(2,3) = 0.0_wp  !sin(theta)
    Rx(3,3) = -1.0_wp !cos(theta)
    Rx(1,1) = 1.0_wp

    v = 0.0_wp
    v = matmul(Rx,u)
    u = v
    return
  end subroutine rotRx180

!================================================================================!
!> rodriguez rotation of v, arount vector k  with angle theta
  subroutine rodrot(v,k,theta)
    implicit none
    real(wp) :: v(3)    !> point to be rotated
    real(wp) :: k(3)    !> the rotation axis !UNIT VECTOR!
    real(wp) :: theta   !> angle of rotation around k
    real(wp) :: vnew(3)
    real(wp) :: kxv(3)
    real(wp) :: kv
    integer :: i

    call unitv(k) !> make k into a unit vector (just in case)

    !> v' = v*cos(theta) + (k x v)sin(theta) + k(k*v)(1-cos(theta)
    vnew = 0
    call crosp(k,v,kxv)
    kv = dotp(k,v,3)
    do i = 1,3
      vnew(i) = v(i) * cos(theta)
      vnew(i) = vnew(i) + kxv(i) * sin(theta)
      vnew(i) = vnew(i) + k(i) * kv * (1 - cos(theta))
    end do
    v = vnew
    return
  end subroutine rodrot
!================================================================================!
!================================================================================!
end module geo
