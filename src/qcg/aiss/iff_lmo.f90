!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022-2026 Christoph Plett, Philipp Pracht
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
!
! Ported from the xtb IFF submodule (github.com/grimme-lab/xtb,
! src/iff/iff_lmo.f90), LGPL-3.0-or-later.
!================================================================================!

!> Off-center LMO point charges for the anisotropic IFF electrostatics
module aiss_iff_lmo
  use crest_parameters,only:wp,stdout
  use aiss_param
  implicit none
  private

  public :: setlmo

  integer,parameter :: maxlmo = 50000
  !> per LMO charge: type (1 atom, 2 sigma, 4 LP/pi) and anchor atoms
  integer :: lmoatom(4,maxlmo,2) = 0

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine setlmo(imol,n,nl,at,xyz,alp0,q,cn,c,lmo,dipol_read,verbose)
    !*****************************************************************
    !* Replace the LMO centers c(1:3,:) and types lmo(:) of fragment
    !* imol by off-center point charges c(1:4,:) for LPs, sigma holes
    !* and pi bonds. On output lmo(i) is the atom a charge belongs to
    !* and nl the number of charges.
    !*****************************************************************
    integer,intent(in)     :: imol,n
    real(wp),intent(in)    :: xyz(3,n),alp0(n),q(n),cn(n)
    integer,intent(in)     :: at(n)
    real(wp),intent(inout) :: c(4,10*n)
    integer,intent(inout)  :: lmo(10*n),nl
    real(wp),intent(in)    :: dipol_read
    logical,intent(in)     :: verbose

    integer  :: i,j,ii,ia1,ind(n),ltmp(10*n),nlp,npi,nsig,ia2,ia3
    real(wp) :: ex,ey,ez,f,norm(n),dip(3),di,ff,rr1
    real(wp) :: alp,r1,r2,qq,c1,c2,ex1,ey1,ez1,rr,rr2
    real(wp) :: r0_free,ctmp(4,10*n),prxs(3,3),eig(3),aa,vtmp(3),vnew(3)
    real(wp) :: s1,s2
    logical  :: sigbond,pibond

    lmoatom(1,1:maxlmo,imol) = 0
    ltmp = lmo
    ctmp = c
    c = 0

    !> a single atom has only LPs, which would break its symmetry
    if (n .eq. 1) then
      nl = 0
      return
    end if

    npi = 0
    nlp = 0
    nsig = 0
    i = 0
    do ii = 1,nl
      ex = ctmp(1,ii)
      ey = ctmp(2,ii)
      ez = ctmp(3,ii)
      do ia1 = 1,n
        norm(ia1) = sqrt((xyz(1,ia1)-ex)**2+(xyz(2,ia1)-ey)**2+(xyz(3,ia1)-ez)**2)
        ind(ia1) = ia1
      end do
      call qsort(norm,1,n,ind)
      ia3 = 0
      ia1 = ind(1)
      ia2 = ind(1)
      if (n .gt. 1) ia2 = ind(2)
      if (n .gt. 2) ia3 = ind(3)

!>-- LP: minus charge on the atom-LP axis, plus charge at the nucleus
      if (ltmp(ii) .eq. 2.and.lpatom(at(ia1)) .eq. 1) then
        i = i+1
        if (i .gt. maxlmo) error stop 'too many LMO'
        f = -par_pos1_lp
        c(1,i) = xyz(1,ia1)+f*(xyz(1,ia1)-ex)
        c(2,i) = xyz(2,ia1)+f*(xyz(2,ia1)-ey)
        c(3,i) = xyz(3,ia1)+f*(xyz(3,ia1)-ez)
        c(4,i) = -0.010*par_chrg_lp*alp0(ia1)/r0_atom(at(ia1))
        if (at(ia1) .gt. 10) c(4,i) = -0.010*par_chrg_lp*alp0(ia1)/r0_atom(at(ia1))/1.8
        lmo(i) = ia1
        lmoatom(1,i,imol) = 4
        lmoatom(2,i,imol) = ia1
        !> lmoatom(3:4) = 0 marks an LP
        lmoatom(3,i,imol) = 0
        lmoatom(4,i,imol) = 0
        i = i+1
        if (i .gt. maxlmo) error stop 'too many LMO'
        f = -par_pos2_lp
        c(1,i) = xyz(1,ia1)+f*(xyz(1,ia1)-ex)
        c(2,i) = xyz(2,ia1)+f*(xyz(2,ia1)-ey)
        c(3,i) = xyz(3,ia1)+f*(xyz(3,ia1)-ez)
        c(4,i) = -c(4,i-1)
        lmo(i) = ia1
        lmoatom(1,i,imol) = 1
        lmoatom(2,i,imol) = ia1
        nlp = nlp+1
      end if

      sigbond = (sigatom(at(ia1)) .eq. 1.or.sigatom(at(ia2)) .eq. 1) &
      &         .and.(ltmp(ii) .eq. 1).and.(n .ne. 1).and.at(ia1) .ne. 1 &
      &         .and.at(ia2) .ne. 1
      pibond = ltmp(ii) .ge. 3
      if (sigbond.or.pibond) then
        i = i+1
        if (i .gt. maxlmo) error stop 'too many LMO'
        alp = (alp0(ia1)*alp0(ia2))**0.5
        r0_free = (r0_atom(at(ia1))*r0_atom(at(ia2)))**0.5
        c(1:3,i) = ctmp(1:3,ii)
        lmo(i) = ia1
        r1 = sqrt((ex-xyz(1,ia1))**2+(ey-xyz(2,ia1))**2+(ez-xyz(3,ia1))**2)
        r2 = sqrt((ex-xyz(1,ia2))**2+(ey-xyz(2,ia2))**2+(ez-xyz(3,ia2))**2)
        c1 = r1/(r1+r2)
        c2 = r2/(r1+r2)
        if (sigbond) then
!>-- sigma hole: minus charge at the LMO center, two plus charges along the bond
          rr = sqrt((xyz(1,ia1)-xyz(1,ia2))**2+(xyz(2,ia1)-xyz(2,ia2))**2 &
          &        +(xyz(3,ia1)-xyz(3,ia2))**2)
          ex1 = xyz(1,ia2)-xyz(1,ia1)
          ey1 = xyz(2,ia2)-xyz(2,ia1)
          ez1 = xyz(3,ia2)-xyz(3,ia1)
          qq = 0.01*par_chrg_sig*alp/r0_free
          c(4,i) = -qq
          lmoatom(1,i,imol) = 2
          lmoatom(2,i,imol) = ia1
          lmoatom(3,i,imol) = ia2
          i = i+1
          if (i .gt. maxlmo) error stop 'too many LMO'
          c(1,i) = ex-par_pos_sig*ex1/rr
          c(2,i) = ey-par_pos_sig*ey1/rr
          c(3,i) = ez-par_pos_sig*ez1/rr
          c(4,i) = qq*c1
          lmo(i) = ia1
          lmoatom(1,i,imol) = 2
          lmoatom(2,i,imol) = ia1
          lmoatom(3,i,imol) = ia2
          i = i+1
          if (i .gt. maxlmo) error stop 'too many LMO'
          c(1,i) = ex+par_pos_sig*ex1/rr
          c(2,i) = ey+par_pos_sig*ey1/rr
          c(3,i) = ez+par_pos_sig*ez1/rr
          c(4,i) = qq*c2
          lmo(i) = ia2
          lmoatom(1,i,imol) = 2
          lmoatom(2,i,imol) = ia2
          lmoatom(3,i,imol) = ia1
          nsig = nsig+1
        else
!>-- pi: minus charge at the LMO center, plus charges on both atoms
          ff = 1.0
          if (ltmp(ii) .eq. 4) ff = 3.0
          qq = ff*0.01*par_chrg_pi*alp/r0_free
          c(4,i) = -qq
          lmoatom(1,i,imol) = 4
          lmoatom(2,i,imol) = ia1
          lmoatom(3,i,imol) = ia2
          lmoatom(4,i,imol) = ia3
          i = i+1
          if (i .gt. maxlmo) error stop 'too many LMO'
          c(1:3,i) = xyz(1:3,ia1)
          c(4,i) = qq*c1
          lmo(i) = ia1
          lmoatom(1,i,imol) = 1
          lmoatom(2,i,imol) = ia1
          i = i+1
          if (i .gt. maxlmo) error stop 'too many LMO'
          c(1:3,i) = xyz(1:3,ia2)
          c(4,i) = qq*c2
          lmo(i) = ia2
          lmoatom(1,i,imol) = 1
          lmoatom(2,i,imol) = ia2
          npi = npi+1
        end if
      end if
    end do

!>-- linear molecules: project off-center charges onto the molecular axis
    call inertia_axes(n,at,xyz,eig,prxs)
    if (eig(1) .lt. 3.0e-4_wp) then
      do j = 1,i
        if (lmoatom(1,j,imol) .le. 1) cycle
        ia1 = lmoatom(2,j,imol)
        vtmp(1:3) = c(1:3,j)-xyz(1:3,ia1)
        s2 = sqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
        s1 = dot_product(vtmp,prxs(1:3,1))
        vnew(1:3) = s1*prxs(1:3,1)
        c(1:3,j) = xyz(1:3,ia1)+vnew(1:3)

        !> rescale by 1+sin^2 of the projection angle (cos = 0 for a
        !> charge perpendicular to the axis)
        s1 = sqrt(vnew(1)**2+vnew(2)**2+vnew(3)**2)
        if (s1*s2 > 1.0e-10_wp) then
          aa = (vtmp(1)*vnew(1)+vtmp(2)*vnew(2)+vtmp(3)*vnew(3))/(s1*s2)
        else
          aa = 0.0_wp
        end if
        aa = 2.0_wp-aa**2
        ia2 = lmoatom(2,j+1,imol)
        if (ia1 .eq. ia2) then
          if (lmoatom(3,j,imol) .eq. 0) then
            c(4,j) = c(4,j)*aa
            c(4,j+1) = c(4,j+1)*aa
            cycle
          end if
          c(4,j) = c(4,j)*aa
          c(4,j+1) = c(4,j+1)*aa
          c(4,j+2) = c(4,j+2)*aa
        end if
      end do
    end if

    nl = i
    if (verbose) then
      write (stdout,'(''# of LP/pi/sigma LMOs    :'',3i5)') nlp,npi,nsig
      write (stdout,'(''# of off-center charges  :'', i5)') nl
    end if

    dip = 0
    do i = 1,nl
      dip(1:3) = dip(1:3)+c(1:3,i)*c(4,i)
    end do
    do i = 1,n
      dip(1:3) = dip(1:3)+xyz(1:3,i)*q(i)
    end do
    di = sqrt(dip(1)**2+dip(2)**2+dip(3)**2)
    if (verbose) write (stdout,'(''dipole moment (read/calc):'',2F8.3)') dipol_read,di
  end subroutine setlmo

!========================================================================================!

  subroutine inertia_axes(n,at,xyz,eig,evec)
    !*****************************************************************
    !* Principal moments (amu*Angstrom^2, ascending, below 3e-4
    !* zeroed) and axes, as xtb's axisvec for the linearity check
    !*****************************************************************
    integer,intent(in)   :: n,at(n)
    real(wp),intent(in)  :: xyz(3,n)
    real(wp),intent(out) :: eig(3),evec(3,3)
    real(wp) :: cma(3),sumw,eigw(3),coord(3,n)
    integer :: i
    !> axisvec works in Angstrom and does not divide by the total mass
    coord = xyz*autoang
    call axis_docking(n,at,coord,cma,sumw,eigw,evec)
    eig = eigw*sumw
    do i = 1,3
      if (eig(i) .lt. 3.0e-4_wp) eig(i) = 0.0_wp
    end do
  end subroutine inertia_axes

!========================================================================================!
end module aiss_iff_lmo
