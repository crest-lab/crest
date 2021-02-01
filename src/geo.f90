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
!rodriguez rotation of v, arount vector k  with angle theta
subroutine rodrot(v,k,theta)
     use iso_fortran_env , wp => real64
     implicit none
     real(wp) :: v(3)    !point to be rotated
     real(wp) :: k(3)    !the rotation axis !UNIT VECTOR!
     real(wp) :: theta   !angle of rotation around k
     real(wp) :: vnew(3)
     real(wp) :: kxv(3)
     real(wp) :: kv,dotp
     integer :: i

     call unitv(k) !make k into a unit vector (just in case)

     !v' = v*cos(theta) + (k x v)sin(theta) + k(k*v)(1-cos(theta)
     vnew = 0
     call crosp(k,v,kxv)
     kv=dotp(k,v,3)
     do i=1,3
     vnew(i) = v(i) * cos(theta)
     vnew(i) = vnew(i) + kxv(i)*sin(theta)
     vnew(i) = vnew(i) + k(i)*kv*(1-cos(theta))
     enddo
     v=vnew
     return
end subroutine rodrot

!normalize the vector k
subroutine unitv(k)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: k(3)
     real(wp) :: klen
     klen = sqrt( k(1)**2 + k(2)**2 + k(3)**2)
     k = k/klen
     return
end subroutine unitv

!calculate angle theta between vectors u and v
function tangle(u,v)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: tangle
     real(wp) :: u(3)
     real(wp) :: v(3)
     real(wp) :: ulen
     real(wp) :: vlen
     real(wp) :: uv,dotp

     tangle =0.0_wp
     ulen =sqrt(u(1)**2 + u(2)**2 + u(3)**2)
     vlen =sqrt(v(1)**2 + v(2)**2 + v(3)**2)
     !uv = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
     uv = dotp(u,v,3)
     !acos returns a value between 0 and pi
     tangle = acos( uv / ( ulen*vlen) )
     return
end function tangle

!calculate dot product between two vectors u and v
function dotp(u,v,n)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: dotp
     integer :: n     !dimension
     real(wp) :: u(n) !first vector
     real(wp) :: v(n) !second vector
     integer :: i
     dotp=0.0_wp
     do i=1,n
     dotp = dotp + u(i)*v(i)
     enddo
     return
end function dotp

!calculate the cross product between two vectors u and v
subroutine crosp(u,v,uv)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: u(3) !first vector
     real(wp) :: v(3) !second vector
     real(wp) :: uv(3)
     uv = 0
     uv(1) = u(2)*v(3) - u(3)*v(2)
     uv(2) = u(3)*v(1) - u(1)*v(3)
     uv(3) = u(1)*v(2) - u(2)*v(1)
     return
end subroutine crosp

!rotate around z axis
subroutine rotRz(u,theta)
 use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: u(3) !first vector
     real(wp) :: v(3) !second vector
     real(wp) :: theta
     real(wp) :: Rz(3,3)

     Rz = 0.0_wp
     Rz(1,1) = cos(theta)
     Rz(2,1) = sin(theta)
     Rz(1,2) = -sin(theta)
     Rz(2,2) = cos(theta)
     Rz(3,3) = 1.0_wp

     v = 0.0_wp

     v(1) = Rz(1,1)*u(1) + Rz(1,2)*u(2) + Rz(1,3)*u(3)
     v(2) = Rz(2,1)*u(1) + Rz(2,2)*u(2) + Rz(2,3)*u(3)
     v(3) = u(3)

     u = v
     return
end subroutine rotRz
subroutine rotRz180(u)
 use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: u(3) !first vector
     real(wp) :: v(3) !second vector
     real(wp) :: theta
     real(wp) :: Rz(3,3)

     Rz = 0.0_wp
     Rz(1,1) = -1.0_wp !cos(theta)
     Rz(2,1) = 0.0_wp  !sin(theta)
     Rz(1,2) = 0.0_wp  !-sin(theta)
     Rz(2,2) = -1.0_wp !cos(theta)
     Rz(3,3) = 1.0_wp

     v = 0.0_wp

     v(1) = Rz(1,1)*u(1) + Rz(1,2)*u(2) + Rz(1,3)*u(3)
     v(2) = Rz(2,1)*u(1) + Rz(2,2)*u(2) + Rz(2,3)*u(3)
     v(3) = u(3)

     u = v
     return
end subroutine rotRz180


!rotate around y axis
subroutine rotRy(u,theta)
 use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: u(3) !first vector
     real(wp) :: v(3) !second vector
     real(wp) :: theta
     real(wp) :: Ry(3,3)

     Ry = 0.0_wp
     Ry(1,1) = cos(theta)
     Ry(3,1) = -sin(theta)
     Ry(3,1) = sin(theta)
     Ry(3,3) = cos(theta)
     Ry(2,2) = 1.0_wp

     v = 0.0_wp

     v(1) = Ry(1,1)*u(1) + Ry(1,2)*u(2) + Ry(1,3)*u(3)
     v(2) = u(2)
     v(3) = Ry(3,1)*u(1) + Ry(3,2)*u(2) + Ry(3,3)*u(3) 

     u = v
     return
end subroutine rotRy
subroutine rotRy180(u)
 use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: u(3) !first vector
     real(wp) :: v(3) !second vector
     real(wp) :: theta
     real(wp) :: Ry(3,3)

     Ry = 0.0_wp
     Ry(1,1) = -1.0_wp !cos(theta)
     Ry(3,1) = 0.0_wp  !-sin(theta)
     Ry(3,1) = 0.0_wp  !sin(theta)
     Ry(3,3) = -1.0_wp !cos(theta)
     Ry(2,2) = 1.0_wp

     v = 0.0_wp

     v(1) = Ry(1,1)*u(1) + Ry(1,2)*u(2) + Ry(1,3)*u(3)
     v(2) = u(2)
     v(3) = Ry(3,1)*u(1) + Ry(3,2)*u(2) + Ry(3,3)*u(3)

     u = v
     return
end subroutine rotRy180


!rotate around x axis
subroutine rotRx(u,theta)
 use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: u(3) !first vector
     real(wp) :: v(3) !second vector
     real(wp) :: theta
     real(wp) :: Rx(3,3)

     Rx = 0.0_wp
     Rx(2,2) = cos(theta)
     Rx(2,3) = -sin(theta)
     Rx(2,3) = sin(theta)
     Rx(3,3) = cos(theta)
     Rx(1,1) = 1.0_wp

     v = 0.0_wp

    !v(1) = Rx(1,1)*u(1) + Rx(1,2)*u(2) + Rx(1,3)*u(3)
     v(1) = u(1)
     v(2) = Rx(2,1)*u(1) + Rx(2,2)*u(2) + Rx(2,3)*u(3)
     v(3) = Rx(3,1)*u(1) + Rx(3,2)*u(2) + Rx(3,3)*u(3)

     u = v
     return
end subroutine rotRx
!rotate around x axis 180deg
subroutine rotRx180(u)
 use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: u(3) !first vector
     real(wp) :: v(3) !second vector
     real(wp) :: theta
     real(wp) :: Rx(3,3)

     Rx = 0.0_wp
     Rx(2,2) = -1.0_wp !cos(theta)
     Rx(2,3) = 0.0_wp  !-sin(theta)
     Rx(2,3) = 0.0_wp  !sin(theta)
     Rx(3,3) = -1.0_wp !cos(theta)
     Rx(1,1) = 1.0_wp

     v = 0.0_wp

    !v(1) = Rx(1,1)*u(1) + Rx(1,2)*u(2) + Rx(1,3)*u(3)
     v(1) = u(1)
     v(2) = Rx(2,1)*u(1) + Rx(2,2)*u(2) + Rx(2,3)*u(3)
     v(3) = Rx(3,1)*u(1) + Rx(3,2)*u(2) + Rx(3,3)*u(3)

     u = v
     return
end subroutine rotRx180



!-------------------------------------------------------------------------------------------------------------
! (x,y,z) --> (r,φ,θ), where (0,0,0) is the position of the zatom "self"
! r = sqrt(x²+y²+z²) 
! φ = arctan(y/x)
! θ = arccos(z/r)
! on input self%zmat should contain cartesian coordinates
!subroutine cart2pol(self)
!   implicit none
!   class(zatom) :: self
!   real(wp),allocatable :: dum(:,:)
!   integer :: i
!   allocate(dum(3,self%nei))
!   do i=1,self%nei
!   dum(1:3,i)=self%zmat(1:3,i)-self%cart(1:3)
!   self%zmat(1,i)=sqrt(dum(1,i)**2 + dum(2,i)**2 + dum(3,i)**2)
!   self%zmat(2,i)= atan(dum(2,i)/dum(1,i))
!   self%zmat(3,i)= atan(dum(3,i)/self%zmat(1,i))
!   enddo
!   deallocate(dum)
!   return
!end subroutine cart2pol
!!-------------------------------------------------------------------------------------------------------------
!! (r,φ,θ) --> (x,y,z), where (0,0,0) is the position of the zatom "self"
!! x = r sin(θ) cos(φ)
!! y = r sin(θ) sin(φ)
!! z = r cos(θ)
!! on input self%zmat should contain polar coordinates
!subroutine pol2cart(self)
!   implicit none
!   class(zatom) :: self
!   real(wp),allocatable :: dum(:,:)
!   integer :: i
!   allocate(dum(3,self%nei))
!   do i=1,self%nei
!   dum(1:3,i)=self%zmat(1:3,i)
!   self%zmat(1,i)=dum(1,i)* cos(dum(2,i)) * sin(dum(3,i))
!   self%zmat(2,i)=dum(1,i)* sin(dum(2,i)) * sin(dum(3,i))
!   self%zmat(3,i)=dum(1,i)* cos(dum(3,i))
!   self%zmat(1:3,i)=self%zmat(1:3,i)+self%cart(1:3)
!   enddo
!   deallocate(dum)
!   return
!end subroutine pol2cart

