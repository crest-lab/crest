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

module crest_solvation
!********************************************************************************
!* Composite implicit-solvation contribution.
!*
!* Stitches the building-block components into a single, method-independent
!* solvation free energy and gradient that can be added on top of any parent
!* potential:
!*   - charges (+ dq/dR) from the EEQ/EEQ-BC electrostatic component
!*   - polar (electrostatic) free energy from the ddX continuum engine
!*   - nonpolar surface-tension term from the SASA component
!*   - optional charge-dependent hydrogen-bond term  h_i S_i q_i^2
!*
!* All settings and the cached (GFN2/ALPB) parameters are bundled in the
!* solvation_data object, which is stored on the calculation_settings and
!* passed together with the geometry to solvation_core, mirroring the
!* tblite/gfn0/gfnff calculators.
!********************************************************************************
  use crest_parameters
  use strucrd,only:coord
#ifdef WITH_DDX
  use crest_electrostatic,only:electrostatic_core
  use crest_ddx_pc,only:ddx_pc_engrad
  use crest_surface,only:surface_engrad
  use mctc_env,only:error_type
  use mctc_io,only:structure_type,new
  use tblite_solvation_cds,only:cds_input
  use tblite_solvation_data_cds,only:get_cds_param
  use tblite_solvation,only:get_solvent_data,solvent_data
#endif
  implicit none
  private

  public :: solvation_data
  public :: solvation_setup
  public :: solvation_core

!> Bundled solvation settings + cached parameters (stored on calculation_settings)
  type :: solvation_data
    !> user configuration
    character(len=:),allocatable :: charge_model  !> 'eeq' | 'eeqbc'
    character(len=:),allocatable :: smodel        !> 'cosmo' | 'cpcm' | 'pcm'
    character(len=:),allocatable :: solvent       !> e.g. 'water'
    logical  :: do_hbond = .true.
    !> derived/cached by solvation_setup
    logical  :: loaded = .false.
    real(wp) :: eps = 0.0_wp
    real(wp) :: probe = 0.0_wp
    real(wp),allocatable :: tension(:)  !> per atom
    real(wp),allocatable :: hbond(:)    !> per atom (scaled)
    real(wp),allocatable :: rad(:)      !> per species (D3)
  end type solvation_data

contains
!========================================================================================!

  subroutine solvation_setup(mol,solv,iostatus)
    !***********************************************************************
    !* Resolve the solvent dielectric constant and load the GFN2/ALPB CDS
    !* nonpolar parameters for the given molecule into the solvation_data
    !* object.  Idempotent: skips once loaded.
    !***********************************************************************
    type(coord),intent(in)            :: mol
    type(solvation_data),intent(inout) :: solv
    integer,intent(out)               :: iostatus
#ifdef WITH_DDX
    type(solvent_data) :: sdat
    real(wp),allocatable :: tension(:),hbond(:),rad(:)
    real(wp) :: probe
    character(len=:),allocatable :: solv_name

    iostatus = 0
    if (solv%loaded) return
    if (.not.allocated(solv%charge_model)) solv%charge_model = 'eeqbc'
    if (.not.allocated(solv%smodel)) solv%smodel = 'cpcm'
    if (.not.allocated(solv%solvent)) then
      iostatus = 1; return
    end if

    !> dielectric constant from the solvent database (tblite calls water 'water')
    solv_name = solv%solvent
    if (solv_name == 'h2o') solv_name = 'water'
    sdat = get_solvent_data(solv_name)
    if (sdat%eps <= 0.0_wp) then
      iostatus = 1; return
    end if
    solv%eps = sdat%eps

    !> GFN2/ALPB CDS nonpolar parameters
    allocate (tension(mol%nat),hbond(mol%nat))
    call solvation_cds_params(mol,solv_name,tension,hbond,rad,probe,iostatus)
    if (iostatus /= 0) return
    call move_alloc(tension,solv%tension)
    call move_alloc(hbond,solv%hbond)
    call move_alloc(rad,solv%rad)
    solv%probe = probe
    solv%loaded = .true.
#else
    iostatus = 1
#endif
  end subroutine solvation_setup

!========================================================================================!

  subroutine solvation_core(mol,chrg,solv,energy,gradient,iostatus)
    !***********************************************************************
    !* Total implicit-solvation energy and gradient for a prepared
    !* solvation_data object (call solvation_setup first).
    !*
    !*  mol      : molecular structure (Bohr)
    !*  chrg     : total molecular charge
    !*  solv     : bundled settings + cached parameters
    !*  energy   : total solvation free energy (out)
    !*  gradient : 3,nat Cartesian gradient (out, overwritten)
    !***********************************************************************
    type(coord),intent(in)            :: mol
    integer,intent(in)                :: chrg
    type(solvation_data),intent(in)   :: solv
    real(wp),intent(out)              :: energy
    real(wp),intent(out)              :: gradient(:,:)
    integer,intent(out)              :: iostatus
#ifdef WITH_DDX
    real(wp),allocatable :: qat(:),dqdr(:,:,:),gtmp(:,:),sasa(:),dsdr(:,:,:)
    real(wp) :: edum,etmp
    integer :: nat,iat,jat

    iostatus = 0
    nat = mol%nat
    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
    allocate (qat(nat),dqdr(3,nat,nat),gtmp(3,nat),sasa(nat),dsdr(3,nat,nat))

    ! ── charges + dq/dR from the electrostatic component ─────────────────────
    call electrostatic_core(mol,chrg,solv%charge_model,edum,gtmp,qat,iostatus,dqdr=dqdr)
    if (iostatus /= 0) return

    ! ── polar (ddX) free energy + gradient (incl. dq/dR chain term) ──────────
    call ddx_pc_engrad(mol,qat,solv%smodel,solv%eps,etmp,gradient,iostatus,dqdr=dqdr)
    if (iostatus /= 0) return
    energy = etmp

    ! ── nonpolar surface-tension term ────────────────────────────────────────
    call surface_engrad(mol,solv%tension,etmp,gtmp,iostatus,sasa=sasa,dsdr=dsdr, &
    &                   rad=solv%rad,probe=solv%probe)
    if (iostatus /= 0) return
    energy = energy+etmp
    gradient(:,:) = gradient+gtmp

    ! ── charge-dependent hydrogen-bond term:  E = sum_i h_i S_i q_i^2 ─────────
    if (solv%do_hbond) then
      energy = energy+sum(solv%hbond*sasa*qat**2)
      do iat = 1,nat
        do jat = 1,nat
          gradient(:,iat) = gradient(:,iat) &
          &  +solv%hbond(jat)*qat(jat)**2*dsdr(:,iat,jat) &            !> dS/dR part
          &  +solv%hbond(jat)*sasa(jat)*2.0_wp*qat(jat)*dqdr(:,iat,jat) !> dq/dR part
        end do
      end do
    end if
#else
    iostatus = 1
    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
#endif
  end subroutine solvation_core

!========================================================================================!

  subroutine solvation_cds_params(mol,solvent,tension,hbond,rad,probe,iostatus)
    !***********************************************************************
    !* Fetch the GFN2/ALPB CDS nonpolar parameters for a given solvent from
    !* tblite (no hard-coded numbers): per-atom surface tensions, per-atom
    !* scaled hydrogen-bond strengths, per-species (D3) radii and the probe.
    !*
    !* We might want to switch this parameter getter out for something else
    !* that is not hardcoded to the GFN2 calculator
    !***********************************************************************
    type(coord),intent(in)        :: mol
    character(len=*),intent(in)   :: solvent
    real(wp),intent(out)          :: tension(:)
    real(wp),intent(out)          :: hbond(:)
    real(wp),allocatable,intent(out) :: rad(:)
    real(wp),intent(out)          :: probe
    integer,intent(out)           :: iostatus
#ifdef WITH_DDX
    type(structure_type) :: struc
    type(cds_input)      :: inp
    type(error_type),allocatable :: error
    real(wp),allocatable :: hbspec(:)
    integer :: iat

    iostatus = 0
    call new(struc,mol%at,mol%xyz)
    inp%alpb = .true.
    inp%solvent = solvent
    call get_cds_param(inp,struc,'gfn2',error)
    if (allocated(error).or..not.allocated(inp%tension)) then
      iostatus = 1; return
    end if
    probe = inp%probe
    rad = inp%rad
    !> hydrogen-bond strength is scaled by the atomic surface (cf. tblite new_cds)
    allocate (hbspec(size(inp%hbond)))
    hbspec = inp%hbond/(4.0_wp*pi*(inp%rad+inp%probe)**2)
    do iat = 1,mol%nat
      tension(iat) = inp%tension(struc%id(iat))
      hbond(iat) = hbspec(struc%id(iat))
    end do
#else
    iostatus = 1
    tension(:) = 0.0_wp
    hbond(:) = 0.0_wp
    probe = 0.0_wp
#endif
  end subroutine solvation_cds_params

!========================================================================================!
end module crest_solvation
