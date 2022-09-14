!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 Christoph Plett
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

! Adapted from
! Ján Busa, Jozef Dzurina, Edik Hayryan, Shura Hayryan, Chin-Kun Hu, Ján Plavka, Imrich Pokorný,
! Jaroslav Skrivánek, Ming-Chya Wu
! Comput. Phys. Commun. 165(2005)59

subroutine get_volume(zmol, rad)
   use iso_fortran_env, wp => real64
   use zdata
   implicit none
   type(Zmolecule), intent(inout) :: zmol
   real(wp), intent(in)  :: rad(zmol%nat)
   real(wp), allocatable :: xyz_rad(:, :)
   integer, allocatable  :: neigh_list(:)
   integer, allocatable  :: neigh_index(:)
   integer, allocatable  :: neigh_type(:)
   real(wp)              :: va_part(2)
   integer               :: i

   allocate (xyz_rad(zmol%nat, 4), neigh_list(zmol%nat), neigh_index(zmol%nat))
   allocate (neigh_type(zmol%nat**2))

   zmol%vtot = 0d0
   zmol%atot = 0d0

!--- Copying Input
   do i = 1, zmol%nat
      xyz_rad(i, 1:3) = zmol%xyz(1:3, i)
      xyz_rad(i, 4) = rad(i)
   end do

!--- Checking neighbors (different to usual CREST neighbors to account for more atoms)
   call create_neigh(zmol%nat, xyz_rad, neigh_list, &
   &                               neigh_index, neigh_type)

!--- Compute V and A
   do i = 1, zmol%nat
      call calcVA(i, xyz_rad, neigh_list, neigh_index, &
      &                        neigh_type, zmol%nat, va_part)
      zmol%vtot = zmol%vtot + va_part(1)
      zmol%atot = zmol%atot + va_part(2)
   end do

   deallocate (xyz_rad, neigh_type, neigh_index)
   deallocate (neigh_list)

   return
end subroutine get_volume

subroutine create_neigh(nat, xyz_rad, neigh_list, neigh_index, neigh_type)
   use iso_fortran_env, wp => real64
   implicit none

   integer, intent(in)  :: nat
   real(wp), intent(in) :: xyz_rad(nat, 4)
   integer              :: neigh_list(nat), neigh_index(nat), neigh_type(nat**2)

   integer              :: neigh_tmp(nat), dum(nat)
   integer              :: i, j
   real(wp)             :: x, y, z, d, ri, r

   neigh_index = 0
   neigh_index(1) = 1
   neigh_list = 0
   neigh_tmp = 0
   do i = 1, nat
      !--- Check, if there are neighbors and which
      neigh_list(i) = 0
      x = xyz_rad(i, 1)
      y = xyz_rad(i, 2)
      z = xyz_rad(i, 3)
      r = xyz_rad(i, 4)
      do j = 1, nat
         if (j .NE. i) then
            if (dabs(x - xyz_rad(j, 1)) .lt. r + xyz_rad(j, 4)) then
               d = dsqrt((x - xyz_rad(j, 1))**2 + (y - xyz_rad(j, 2))**2 + (z - xyz_rad(j, 3))**2)
               ri = xyz_rad(j, 4)
               if (d .lt. r + ri) then
                  if (d + r .LE. ri) then
                     neigh_list(i) = -1
                     exit
                  elseif (d + ri .gt. r) then
                     neigh_list(i) = neigh_list(i) + 1
                     neigh_tmp(neigh_list(i)) = j
                  end if
               end if
            end if
         end if
      end do
      dum = neigh_list !Somhow the first entry in neigh_list is overwritten in the following do cycle

      !--- No neighbors
      if (neigh_list(i) .LE. 0) then
         neigh_index(i + 1) = neigh_index(i)
      !--- Neighbors
      else
         if(i < nat) then
            neigh_index(i + 1) = neigh_index(i) + neigh_list(i)
         end if
         do j = 1, neigh_list(i)
            neigh_type(neigh_index(i) + j - 1) = neigh_tmp(j)
         end do
      end if
   end do
   neigh_list = dum

   return
end subroutine create_neigh

subroutine calcVA(num, xyz_rad, neigh_list, neigh_index, neigh_type, nat, va_part)
   use iso_fortran_env, wp => real64
   implicit none

   integer, intent(in)   :: num, nat
   real(wp), intent(in)  :: xyz_rad(nat, 4)
   integer, intent(in)   :: neigh_list(nat), neigh_index(nat), neigh_type(nat**2)
   real(wp), intent(out) :: va_part(2)
   real(wp)              :: circles(nat, 4), single_sphere(nat, 4)
   real(wp)              :: int_parts(nat**2, 3), av_part(2)
   real(wp)              :: rad

   integer               :: neigh_tmp(nat), nint_parts, npos
   integer               :: i, j
   real(wp), parameter   :: pi = 3.1415926540d0

   !--- No neighbors
   if (neigh_list(num) .eq. 0) then
      va_part(1) = 4d0*pi*xyz_rad(num, 4)**3/3.d0
      va_part(2) = 4d0*pi*xyz_rad(num, 4)**2
   !--- Subset
   elseif (neigh_list(num) .lt. 0) then
      va_part(1) = 0d0
      va_part(2) = 0d0
   !--- Neighbors exist
   else

      neigh_tmp(1) = num

      do i = 1, (neigh_list(num))
         neigh_tmp(i + 1) = neigh_type(neigh_index(num) + i - 1)
      end do
      do i = 1, neigh_list(num) + 1
         do j = 1, 4
            single_sphere(i, j) = xyz_rad(neigh_tmp(i), j)
         end do
      end do

      va_part(1) = 0d0
      va_part(2) = 0d0

      call generate_integration_parts(single_sphere, circles, nat, neigh_list(num), int_parts, nint_parts)

      npos = 0
      do i = 1, (neigh_list(num))
         if (circles(i, 4) .gt. 0) then
            npos = npos + 1
         end if
      end do

      rad = single_sphere(1, 4)
      !--- Selective integration as overlap was found
      if (npos .gt. 0) then
         call integrate(circles, int_parts, nat, nint_parts, rad, single_sphere(1, 3), av_part)
         va_part(1) = va_part(1) + av_part(1)
         va_part(2) = va_part(2) + av_part(2)
      !--- Complete integration
      else
         call integrate(circles, int_parts, nat, nint_parts, rad, single_sphere(1, 3), av_part)
         va_part(1) = va_part(1) + av_part(1) + 4d0*pi*single_sphere(1, 4)**3/3d0
         va_part(2) = va_part(2) + av_part(2) + 4d0*pi*single_sphere(1, 4)**2
      end if
   end if

   return
end subroutine calcVA

subroutine generate_integration_parts(single_sphere, circles, nat, num_neigh, int_parts, num_parts)

   use iso_fortran_env, wp => real64
   implicit none
   real(wp), intent(in)  :: single_sphere(nat, 4)
   real(wp), intent(out) :: circles(nat, 4)
   integer, intent(in)   :: nat, num_neigh
   real(wp), intent(out) :: int_parts(nat**2, 3)
   integer, intent(out)  :: num_parts

   integer               :: nna
   real(wp)              :: int_partsnew(nat**2, 3), rad, x, y, a, b, c, d
   integer               :: i, j, k
   real(wp), parameter   :: pi = 3.1415926540d0

   num_parts = 0

   !--- Create circles first
   rad = single_sphere(1, 4)
   do i = 1, (num_neigh)
      x = single_sphere(1, 1) - single_sphere(i + 1, 1)
      y = single_sphere(1, 2) - single_sphere(i + 1, 2)
      a = 8d0*rad**2*x
      b = 8d0*rad**2*y
      c = x**2 + y**2 + (single_sphere(1, 3) + rad - single_sphere(i + 1, 3))**2 - single_sphere(i + 1, 4)**2
      d = 4d0*rad**2*(x**2 + y**2 + (single_sphere(1, 3) - rad - single_sphere(i + 1, 3))**2 - single_sphere(i + 1, 4)**2)
      circles(i, 1) = -a/(2d0*c)
      circles(i, 2) = -b/(2d0*c)
      circles(i, 3) = dsqrt((a**2 + b**2 - 4d0*c*d)/(4d0*c**2))
      if (c .gt. 0) then
         circles(i, 4) = -1
      else
         circles(i, 4) = 1
      end if
   end do

   !--- And than integration parts
   !--- Only one circle
   if (num_neigh .eq. 1) then
      num_parts = 1
      int_parts(1, 1) = 1
      int_parts(1, 2) = 0d0
      int_parts(1, 3) = 2d0*pi*circles(1, 4)
   !--- More circles
   else
      do i = 1, (num_neigh)
         call make_parts(i, circles, nat, num_neigh, nna, int_partsnew)
         if (nna .gt. 0) then
            do j = 1, nna
               do k = 1, 3
                  int_parts(num_parts + j, k) = int_partsnew(j, k)
               end do
            end do
            num_parts = num_parts + nna
         end if
      end do
   end if
   return
end subroutine generate_integration_parts

!--- Create parts that are integratet later
subroutine make_parts(num, circles, nat, num_neigh, no_arc, int_partsnew)
   use iso_fortran_env, wp => real64
   implicit none

   integer, intent(in)   :: num
   real(wp), intent(in)  :: circles(nat, 4)
   integer, intent(in)   :: nat, num_neigh
   real(wp), intent(out) :: int_partsnew(nat**2, 3)
   integer, intent(out)  :: no_arc

   integer               :: i, j, m, counter_angles, counter_no_angles, counter
   real(wp)              :: c11, c12, c13, c21, c22, c23, dist, int11, int12, int21, int22, dum
   real(wp)              :: angles(nat**2), anglesnew(nat**2)
   logical               :: minmax
   real(wp), parameter   :: pi = 3.1415926540d0

   no_arc = 0
   counter_angles = 0

   c11 = circles(num, 1)
   c12 = circles(num, 2)
   c13 = circles(num, 3)
   do i = 1, (num_neigh)
      if (i .NE. num) then
         c21 = circles(i, 1)
         c22 = circles(i, 2)
         c23 = circles(i, 3)
         dist = dsqrt((c11 - c21)**2 + (c12 - c22)**2)
         if ((dist .lt. c23 + c13) .and. (dabs(c23 - c13) .lt. dist)) then
            !--- Two intersections
            call intersection(num, i, nat, circles, int11, int12, int21, int22)
            angles(counter_angles + 1) = int11
            angles(counter_angles + 2) = int12
            counter_angles = counter_angles + 2
         end if
      end if
   end do
   if (counter_angles .eq. 0) then
      counter_no_angles = 0
      do i = 1, (num_neigh)
         if (i .NE. num) then

            !--- Check overlapping circles
            dist = dsqrt((circles(num, 1) + circles(num, 3) - circles(i, 1))**2 + &
            &        (circles(num, 2) - circles(i, 2))**2)
            if (dist .lt. circles(i, 3)) then
               if (circles(i, 4) .gt. 0) then
                  counter = 1
               else
                  counter = 0
               end if
            elseif (dist .gt. circles(i, 3)) then
               if (circles(i, 4) .gt. 0) then
                  counter = 0
               else
                  counter = 1
               end if
            else
               dist = dsqrt((circles(num, 1) - circles(i, 1))**2 + (circles(num, 2) - circles(i, 2))**2)
               if (dist .lt. circles(i, 3)) then
                  if (circles(i, 4) .gt. 0) then
                     counter = 1
                  else
                     counter = 0
                  end if
               else
                  if (circles(i, 4) .gt. 0) then
                     counter = 0
                  else
                     counter = 1
                  end if
               end if
            end if

            counter_no_angles = counter_no_angles + counter
         end if
      end do
      if (counter_no_angles .eq. (num_neigh - 1)) then
         no_arc = 1
         int_partsnew(1, 1) = num
         int_partsnew(1, 2) = 0d0
         int_partsnew(1, 3) = 2d0*pi*circles(num, 4)
      end if
   else
      if (circles(num, 4) .gt. 0) then
         minmax = .true.
      else
         minmax = .false.
      end if
      !--- Sort angles
      do i = 1, (counter_angles - 1)
         counter = i
         dum = angles(i)
         do j = i + 1, counter_angles
            if (minmax) then
               if (dum .gt. angles(j)) then
                  counter = j
                  dum = angles(j)
               end if
            else
               if (dum .lt. angles(j)) then
                  counter = j
                  dum = angles(j)
               end if
            end if
         end do
         if (counter .NE. i) then
            angles(counter) = angles(i)
            angles(i) = dum
         end if
      end do

      !--- Remove equals
      m = 1
      anglesnew(1) = angles(1)
      do i = 2, counter_angles
         if (dabs(angles(i) - angles(i - 1)) .gt. 1d-12) then
            m = m + 1
            anglesnew(m) = angles(i)
         end if
      end do
      counter_angles = m
      do i = 1, m
         angles(i) = anglesnew(i)
      end do
      do i = 1, (counter_angles - 1)
         counter_no_angles = 0
         do j = 1, (num_neigh)
            if (j .NE. num) then
               c21 = c11 + c13*dcos((angles(i) + angles(i + 1))/2d0)
               c22 = c12 + c13*dsin((angles(i) + angles(i + 1))/2d0)
               !--- Check, if point is inside circle
               dist = dsqrt((c21 - circles(j, 1))**2 + (c22 - circles(j, 2))**2)
               if (dist .lt. circles(j, 3)) then
                  if (circles(j, 4) .gt. 0) then
                     counter_no_angles = counter_no_angles + 1
                  end if
               else
                  if (circles(j, 4) .LE. 0) then
                     counter_no_angles = counter_no_angles + 1
                  end if
               end if

            end if
         end do
         if (counter_no_angles .eq. (num_neigh - 1)) then
            no_arc = no_arc + 1
            int_partsnew(no_arc, 1) = num
            int_partsnew(no_arc, 2) = angles(i)
            int_partsnew(no_arc, 3) = angles(i + 1) - angles(i)
         end if
      end do
      counter_no_angles = 0
      do i = 1, (num_neigh)
         if (i .NE. num) then
            c21 = c11 + c13*dcos((angles(1) + 2d0*pi + angles(counter_angles))/2d0)
            c22 = c12 + c13*dsin((angles(1) + 2d0*pi + angles(counter_angles))/2d0)
            !--- Check, if point is inside circle
            dist = dsqrt((c21 - circles(i, 1))**2 + (c22 - circles(i, 2))**2)
            if (dist .lt. circles(i, 3)) then
               if (circles(i, 4) .gt. 0) then
                  counter_no_angles = counter_no_angles + 1
               end if
            else
               if (circles(i, 4) .LE. 0) then
                  counter_no_angles = counter_no_angles + 1
               end if
            end if

         end if
      end do
      if (counter_no_angles .eq. (num_neigh - 1)) then
         no_arc = no_arc + 1
         int_partsnew(no_arc, 1) = num
         int_partsnew(no_arc, 2) = angles(counter_angles)
         int_partsnew(no_arc, 3) = angles(1) + circles(num, 4)*2d0*pi - angles(counter_angles)
      end if
   end if

   return
end subroutine make_parts

!--- Computation of angles of two intersection points
subroutine intersection(point1, point2, nat, circles, int11, int12, int21, int22)

   use iso_fortran_env, wp => real64
   implicit none
   integer, intent(in)   :: point1, point2, nat
   real(wp), intent(in)  :: circles(nat, 4)
   real(wp), intent(out) :: int11, int12, int21, int22

   real(wp)              :: c11, c12, c13, c21, c22, c23, f1, f2, f3, f4
   real(wp), parameter   :: pi = 3.1415926540d0

   c11 = circles(point1, 1)
   c12 = circles(point1, 2)
   c13 = circles(point1, 3)
   c21 = circles(point2, 1)
   c22 = circles(point2, 2)
   c23 = circles(point2, 3)
   if (dabs(c21 - c11) .lt. 1d-12) then
      f1 = ((c13**2 - c23**2)/(c22 - c12) - (c22 - c12))/2d0
      f2 = dsqrt(c23**2 - f1**2)
      if (f1 .eq. 0) then
         int21 = 0d0
         int22 = pi
      elseif (f1 .gt. 0) then
         int21 = datan(dabs(f1/f2))
         int22 = pi - int21
      else
         int21 = pi + datan(dabs(f1/f2))
         int22 = 3d0*pi - int21
      end if
      f1 = f1 + c22 - c12
      if (f1 .eq. 0) then
         int11 = 0d0
         int12 = pi
      elseif (f1 .gt. 0) then
         int11 = datan(dabs(f1/f2))
         int12 = pi - int11
      else
         int11 = pi + datan(dabs(f1/f2))
         int12 = 3d0*pi - int11
      end if
   else
      f3 = ((c13**2 - c23**2 - (c22 - c12)**2)/(c21 - c11) - (c21 - c11))/2d0
      f4 = (c12 - c22)/(c21 - c11)
      f1 = (-f3*f4 + dsqrt((f4**2 + 1d0)*c23**2 - f3**2))/(f4**2 + 1d0)
      f2 = f3 + f4*f1
      if (f2 .eq. 0) then
         if (f1 .gt. 0) then
            int21 = pi/2d0
         else
            int21 = -pi/2d0
         end if
      elseif (f2 .gt. 0) then
         int21 = datan(f1/f2)
      else
         int21 = pi + datan(f1/f2)
      end if
      f1 = f1 + c22 - c12
      f2 = f2 + c21 - c11
      if (f2 .eq. 0) then
         if (f1 .gt. 0) then
            int11 = pi/2d0
         else
            int11 = -pi/2d0
         end if
      elseif (f2 .gt. 0) then
         int11 = datan(f1/f2)
      else
         int11 = pi + datan(f1/f2)
      end if
      f1 = (-f3*f4 - dsqrt((f4**2 + 1d0)*c23**2 - f3**2))/(f4**2 + 1d0)
      f2 = f3 + f4*f1
      if (f2 .eq. 0) then
         if (f1 .gt. 0) then
            int22 = pi/2d0
         else
            int22 = -pi/2d0
         end if
      elseif (f2 .gt. 0) then
         int22 = datan(f1/f2)
      else
         int22 = pi + datan(f1/f2)
      end if
      f1 = f1 + c22 - c12
      f2 = f2 + c21 - c11
      if (f2 .eq. 0) then
         if (f1 .gt. 0) then
            int12 = pi/2d0
         else
            int12 = -pi/2d0
         end if
      elseif (f2 .gt. 0) then
         int12 = datan(f1/f2)
      else
         int12 = pi + datan(f1/f2)
      end if
   end if
   if (int11 .lt. 0) int11 = int11 + 2d0*pi
   if (int12 .lt. 0) int12 = int12 + 2d0*pi
   if (int21 .lt. 0) int21 = int21 + 2d0*pi
   if (int22 .lt. 0) int22 = int22 + 2d0*pi

   return
end subroutine intersection

subroutine integrate(circles, int_parts, nat, nint_parts, rad, z1, av_part)
   use iso_fortran_env, wp => real64
   implicit none

   real(wp), intent(in)    :: circles(nat, 4), int_parts(nat**2, 3)
   integer, intent(in)     :: nat, nint_parts
   real(wp), intent(in)    :: rad, z1
   real(wp), intent(inout) :: av_part(2)

   integer                 :: i
   real(wp)                :: x, y, z, pre_V, xz, yz, pre_A, f
   real(wp)                :: pa, pd, p1, pb, pc, p2, v1, v2, v3, vJ1, vJ2, vJ3
   real(wp)                :: d_v, d_a, part1, part2, part3, part4
   real(wp), parameter     :: pi = 3.1415926540d0

   av_part(1) = 0d0
   av_part(2) = 0d0

   do i = 1, nint_parts
      x = circles(nint(int_parts(i, 1)), 1) !> int_parts is type real(wp) and therefore
      y = circles(nint(int_parts(i, 1)), 2) !> should not be used as an array index?
      z = circles(nint(int_parts(i, 1)), 3) !> added nint()
      xz = x*z
      yz = y*z
      pre_V = (4d0*rad**2 + x**2 + y**2 + z**2)/2d0
      pre_A = dsqrt(pre_V**2 - xz**2 - yz**2)
      f = z**2 - pre_V
      if (dabs(dabs(int_parts(i, 3)) - 2d0*pi) .lt. 1d-12) then
         v1 = 2d0*pi/pre_A
         v2 = 2d0*pi*pre_V/(pre_A**3)
         v3 = pi*(2d0*pre_V**2 + xz**2 + yz**2)/(pre_A**5)
         vJ1 = pi + f/2d0*v1
         vJ2 = (v1 + f*v2)/4d0
         vJ3 = (v2 + f*v3)/8d0
         d_v = (128d0*vJ3*rad**7 + 8d0*vJ2*rad**5 + &
         &  2d0*vJ1*rad**3)/3d0 - 8d0*rad**4*vJ2*(z1 + rad)
         d_a = 2d0*vJ1*rad**2
         if (int_parts(i, 3) .lt. 0) then
            d_v = -d_v
            d_a = -d_a
         end if
         av_part(1) = av_part(1) + d_v
         av_part(2) = av_part(2) + d_a
      else
         if (int_parts(i, 3) .lt. 0) then
            p2 = int_parts(i, 2) + int_parts(i, 3)
            p1 = int_parts(i, 2)
         else
            p1 = int_parts(i, 2) + int_parts(i, 3)
            p2 = int_parts(i, 2)
         end if
         v1 = 2d0*(pi/2d0 - datan((pre_V*dcos((p1 - p2)/2d0) + &
         &   xz*dcos((p2 + p1)/2d0) + yz*dsin((p2 + p1)/2d0))/ &
         &   (pre_A*dsin((p1 - p2)/2d0))))/pre_A
         pa = dsin(p1)
         pb = dcos(p1)
         pc = dsin(p2)
         pd = dcos(p2)
         part1 = (-xz*pa + yz*pb)/(pre_V + xz*pb + yz*pa)**1
         part2 = (-xz*pc + yz*pd)/(pre_V + xz*pd + yz*pc)**1
         part3 = (-xz*pa + yz*pb)/(pre_V + xz*pb + yz*pa)**2
         part4 = (-xz*pc + yz*pd)/(pre_V + xz*pd + yz*pc)**2
         v2 = (part1 - part2 + pre_V*v1)/(pre_A**2)
         v3 = (part3 - part4 + (part1 - part2)/pre_V + (2d0*pre_V**2 + xz**2 + yz**2)*v2/pre_V)/(2d0*pre_A**2)
         vJ1 = ((p1 - p2) + f*v1)/2d0
         vJ2 = (v1 + f*v2)/4d0
         vJ3 = (v2 + f*v3)/8d0
         d_v = (128d0*vJ3*rad**7 + 8d0*vJ2*rad**5 + &
         &  2d0*vJ1*rad**3)/3d0 - 8d0*rad**4*vJ2*(z1 + rad)
         d_a = 2d0*vJ1*rad**2
         if (int_parts(i, 3) .lt. 0) then
            d_v = -d_v
            d_a = -d_a
         end if
         av_part(1) = av_part(1) + d_v
         av_part(2) = av_part(2) + d_a
      end if
   end do

   return
end subroutine integrate
