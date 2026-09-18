!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2026 Philipp Pracht
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

module crest_surface
!********************************************************************************
!* Solvent-accessible-surface-area (SASA) component.
!*
!* Thin wrapper around tblite's surface integrator providing per-atom SASA and
!* the geometry-only nonpolar (cavity/dispersion) energy E = sum_i gamma_i S_i
!* together with its Cartesian gradient.  This is the surface-tension part of
!* the xTB ALPB/GBSA CDS term; the charge-dependent hydrogen-bond contribution
!* is left to the solvation composite, where atomic charges are available.
!********************************************************************************
  use crest_parameters
  use strucrd,only:coord
#ifdef WITH_TBLITE
  use mctc_io,only:structure_type,new
  use tblite_solvation_surface,only:surface_integrator,new_surface_integrator
  use tblite_solvation_data,only:get_vdw_rad_cosmo
#endif
  implicit none
  private

  public :: surface_engrad

#ifdef WITH_TBLITE
  real(wp),parameter :: probe_default = 1.0_wp*aatoau  !> probe radius (Bohr)
  integer,parameter  :: nang_default = 230             !> Lebedev grid points
#endif

contains
!========================================================================================!

  subroutine surface_engrad(mol,tension,energy,gradient,iostatus,sasa,dsdr,rad,probe)
    !***********************************************************************
    !* Nonpolar SASA energy and gradient.
    !*
    !*  mol      : molecular structure (Bohr)
    !*  tension  : per-atom surface tensions gamma_i
    !*  energy   : nonpolar energy sum_i gamma_i S_i (out)
    !*  gradient : 3,nat Cartesian gradient (out, overwritten)
    !*  sasa     : optional per-atom solvent-accessible surface area (out)
    !*  dsdr     : optional dS_i/dR_a (3,nat,nat); for charge-dependent terms
    !*             (e.g. hydrogen bonding) assembled by the caller
    !*  rad      : optional per-species vdW radii (default: COSMO radii)
    !*  probe    : optional solvent probe radius (default: 1.0 Angstrom)
    !***********************************************************************
    type(coord),intent(in)        :: mol
    real(wp),intent(in)           :: tension(:)
    real(wp),intent(out)          :: energy
    real(wp),intent(out)          :: gradient(:,:)
    integer,intent(out)           :: iostatus
    real(wp),intent(out),optional :: sasa(:)
    real(wp),intent(out),optional :: dsdr(:,:,:)
    real(wp),intent(in),optional  :: rad(:)
    real(wp),intent(in),optional  :: probe
#ifdef WITH_TBLITE
    type(structure_type)     :: struc
    type(surface_integrator) :: integ
    real(wp),allocatable :: vdw(:),surface(:),ds(:,:,:)
    real(wp) :: prb
    integer :: nat,iat,jat

    iostatus = 0
    nat = mol%nat
    energy = 0.0_wp
    gradient(:,:) = 0.0_wp

    call new(struc,mol%at,mol%xyz)

    ! ── per-species radii + probe + surface integrator ───────────────────────
    allocate (vdw(size(struc%num)))
    if (present(rad)) then
      vdw(:) = rad(:)
    else
      vdw(:) = get_vdw_rad_cosmo(struc%num)
    end if
    prb = probe_default
    if (present(probe)) prb = probe
    call new_surface_integrator(integ,struc%id,vdw,prb,nang_default)

    ! ── SASA and its Cartesian derivative ────────────────────────────────────
    allocate (surface(nat),ds(3,nat,nat))
    call integ%get_surface(struc,surface,ds)

    ! ── energy + gradient: E = sum_i gamma_i S_i ─────────────────────────────
    energy = sum(surface*tension)
    do iat = 1,nat
      do jat = 1,nat
        gradient(:,iat) = gradient(:,iat)+tension(jat)*ds(:,iat,jat)
      end do
    end do
    if (present(sasa)) sasa(:) = surface(:)
    if (present(dsdr)) dsdr(:,:,:) = ds(:,:,:)
#else
    iostatus = 1
    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
    if (present(sasa)) sasa(:) = 0.0_wp
#endif
  end subroutine surface_engrad

!========================================================================================!
end module crest_surface
