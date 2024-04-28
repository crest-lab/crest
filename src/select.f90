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

!=================================================================!
! subroutine mrec
! molcount: number of total fragments (increased during search)
! xyz: overall Cart. coordinates
! nat: overall number of atoms
! at: atomic number array
! molvec: assignment vector of atom to fragment
!=================================================================!
subroutine mrec(molcount,xyz,nat,at,molvec)
  use crest_parameters
  use crest_cn_module
  implicit none
  real(wp) :: xyz(3,nat),cn(nat),bond(nat,nat)
  integer  :: nat,molvec(nat),i,molcount,at(nat)
  logical  :: taken(nat)
  molvec = 0
  molcount = 1
  taken = .false.
  cn = 0.0d0
  bond = 0.0d0
  call calc_ncoord(nat,at,xyz,cn,bond)
  do i = 1,nat
    if (.not.taken(i)) then
      molvec(i) = molcount
      taken(i) = .true.
      call neighbours(i,xyz,at,taken,nat,cn,bond,molvec,molcount)
      molcount = molcount+1
    end if
  end do
  molcount = molcount-1
end subroutine mrec

!=================================================================!
! subroutine mreclm
! a variant of the mrec routine with less allocate statements
! should be faster than the old version if several thousand
! calls are done.
!=================================================================!
subroutine mreclm(molcount,nat,at,xyz,molvec,bond,rcov,cn)
  use iso_fortran_env,only:wp => real64
  use crest_cn_module
  implicit none
  integer :: molcount,nat
  integer :: at(nat),molvec(nat)
  real(wp) :: xyz(3,nat)
  real(wp) :: cn(nat),bond(nat,nat)
  real(wp) :: bref(nat,nat)
  real(wp) :: rcov(*)
  integer :: i
  logical :: taken(nat)
  molvec = 0
  molcount = 1
  taken = .false.
  cn = 0.0d0
  bond = 0.0d0
  call calc_ncoord(nat,at,xyz,rcov,cn,400.0_wp,bond)
  bref = bond
  do i = 1,nat
    if (.not.taken(i)) then
      molvec(i) = molcount
      taken(i) = .true.
      call neighbours(i,xyz,at,taken,nat,cn,bond,molvec,molcount)
      molcount = molcount+1
    end if
  end do
  molcount = molcount-1
  bond = bref
end subroutine mreclm

!==================================================================================!
recursive subroutine neighbours(i,xyz,iat,taken,nat,cn,bond,molvec,molcnt)
  use crest_parameters
  implicit none
  real(wp) :: xyz(3,nat),tr,xi(3),cn(nat),bond(nat,nat)
  integer  :: i,nat,molcnt,molvec(nat),j,iat(nat),icn,k
  logical  :: taken(nat)
  tr = 2.0d0

  xi(1:3) = xyz(1:3,i)
  icn = nint(cn(i))
  do k = 1,icn
    j = maxloc(bond(:,i),1)
    bond(j,i) = 0.0d0
    if (i .eq. j) cycle
    if (.not.taken(j)) then
      molvec(j) = molcnt
      taken(j) = .true.
      call neighbours(j,xyz,iat,taken,nat,cn,bond,molvec,molcnt)
    end if
  end do
end subroutine neighbours

!========================================================================================!

subroutine neighdist(natoms,at,xyz,nb,dist)
  use crest_parameters
  use miscdata,only:rad => rcov
  implicit none
  integer :: at(natoms),natoms,nb(200,natoms)
  real(wp) :: xyz(3,natoms),dist(natoms,natoms)
  logical :: da
  integer :: iat,i,j,k,nn
  real(wp) :: dx,dy,dz,r,rco,r2,f

  nb = 0
  nn = min(natoms,2)-1

  do i = 1,natoms
    f = 1.0d0
    k = 0
    do while (k .lt. 1.and.f .lt. 1.5d0)
100   do iat = 1,natoms
        da = .false.
        do j = 1,k
          if (nb(j,i) .eq. iat) da = .true.
        end do
        dx = xyz(1,iat)-xyz(1,i)
        dy = xyz(2,iat)-xyz(2,i)
        dz = xyz(3,iat)-xyz(3,i)
        r2 = dx*dx+dy*dy+dz*dz
        r = sqrt(r2)
        dist(iat,i) = r
        if (iat .ne. i.and.(.not.da)) then
          rco = rad(at(i))+rad(at(iat))
!>-- critical step
          if (r .lt. f*rco.and.k .lt. 199) then
            k = k+1
            nb(k,i) = iat
          end if
        end if
      end do
      !if(k.lt.1.and.f.lt.1.5)then
      f = f*1.1d0
      !   goto 100
      !endif
    end do
    nb(200,i) = k
  end do

end subroutine neighdist

