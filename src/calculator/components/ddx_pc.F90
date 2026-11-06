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

module crest_ddx_pc
!********************************************************************************
!* Standalone ddX point-charge solvation engine.
!*
!* Drives the low-level ddX API (COSMO/CPCM/PCM) with an *externally* supplied
!* set of atomic point charges, returning the electrostatic solvation free
!* energy and its Cartesian gradient.  This is the method-independent analogue
!* of tblite's self-consistent reaction field: instead of variational xTB
!* charges it accepts fixed charges (e.g. EEQ/EEQ-BC) plus their geometry
!* derivative dq/dR, which adds the chain-rule term the SCRF case does not need.
!********************************************************************************
  use crest_parameters
  use strucrd,only:coord
#ifdef WITH_DDX
  use ddx,only:ddx_type,ddx_state_type,ddx_error_type,ddinit,allocate_state, &
  &            setup,fill_guess,fill_guess_adjoint,solve,solve_adjoint, &
  &            solvation_force_terms,check_error
  use ddx_core,only:ddx_electrostatics_type
  use ddx_multipolar_solutes,only:multipole_electrostatics,multipole_psi, &
  &                               multipole_force_terms
  use tblite_solvation_data,only:get_vdw_rad_cosmo
  use tblite_mesh_lebedev,only:grid_size
#endif
  implicit none
  private

  public :: ddx_pc_engrad

contains
!========================================================================================!

  subroutine ddx_pc_engrad(mol,q,smodel,eps,energy,gradient,iostatus,dqdr)
    !***********************************************************************
    !* Electrostatic solvation energy + gradient for fixed point charges.
    !*
    !*  mol      : molecular structure (Bohr)
    !*  q(nat)   : atomic point charges
    !*  smodel   : 'cosmo' | 'cpcm' | 'pcm'
    !*  eps      : solvent dielectric constant
    !*  energy   : solvation free energy (out)
    !*  gradient : 3,nat Cartesian gradient (out, overwritten)
    !*  dqdr     : optional dq_i/dR_a (3,nat,nat); adds the charge-response
    !*             term.  If absent the charges are treated as geometry fixed.
    !***********************************************************************
    type(coord),intent(in)        :: mol
    real(wp),intent(in)           :: q(:)
    character(len=*),intent(in)   :: smodel
    real(wp),intent(in)           :: eps
    real(wp),intent(out)          :: energy
    real(wp),intent(out)          :: gradient(:,:)
    integer,intent(out)           :: iostatus
    real(wp),intent(in),optional  :: dqdr(:,:,:)
#ifdef WITH_DDX
    type(ddx_type)                :: ddx
    type(ddx_state_type)          :: state
    type(ddx_error_type)          :: error
    type(ddx_electrostatics_type) :: elec
    real(wp),allocatable :: rvdw(:),multipoles(:,:),force(:,:),jmat(:,:),dedq(:)
    real(wp) :: feps,shift,sqrt4pi
    integer  :: model,nat,iat,jat,izp
    real(wp),parameter :: conv = 1.0e-9_wp

    iostatus = 0
    nat = mol%nat
    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
    sqrt4pi = sqrt(4.0_wp*pi)

    ! ── model mapping + dielectric scaling factor (cf. tblite) ───────────────
    select case (trim(smodel))
    case ('cosmo')
      model = 1; shift = -1.0_wp; feps = (eps-1.0_wp)/(eps+0.5_wp)
    case ('cpcm')
      model = 1; shift = -1.0_wp; feps = (eps-1.0_wp)/eps
    case ('pcm')
      model = 2; shift = 0.0_wp;  feps = 1.0_wp
    case default
      iostatus = 1; return
    end select

    ! ── van-der-Waals (COSMO) cavity radii ───────────────────────────────────
    allocate (rvdw(nat))
    do iat = 1,nat
      izp = mol%at(iat)
      rvdw(iat) = get_vdw_rad_cosmo(izp)
    end do

    ! ── ddX model + state setup ──────────────────────────────────────────────
    call ddinit(model,nat,mol%xyz,rvdw,eps,ddx,error,force=1,ngrid=grid_size(8), &
    &           lmax=1,eta=0.1_wp,shift=shift)
    if (error%flag /= 0) then; iostatus = 1; return; end if
    call allocate_state(ddx%params,ddx%constants,state,error)
    if (error%flag /= 0) then; iostatus = 1; return; end if

    ! ── build RHS from the (normalized) monopole distribution ────────────────
    allocate (multipoles(1,nat))
    multipoles(1,:) = q(:)/sqrt4pi
    elec%do_phi = .true.; elec%do_e = .true.; elec%do_g = .true.
    call multipole_electrostatics(ddx%params,ddx%constants,ddx%workspace, &
    &                             multipoles,0,elec,error)
    call multipole_psi(ddx%params,multipoles,0,state%psi)
    call setup(ddx%params,ddx%constants,ddx%workspace,state,elec,state%psi,error)
    if (error%flag /= 0) then; iostatus = 1; return; end if

    ! ── solve primal + adjoint linear systems ────────────────────────────────
    call fill_guess(ddx%params,ddx%constants,ddx%workspace,state,conv,error)
    call fill_guess_adjoint(ddx%params,ddx%constants,ddx%workspace,state,conv,error)
    call solve(ddx%params,ddx%constants,ddx%workspace,state,conv,error)
    call solve_adjoint(ddx%params,ddx%constants,ddx%workspace,state,conv,error)
    if (error%flag /= 0) then; iostatus = 1; return; end if

    ! ── energy: feps * 1/2 <xs|psi> ──────────────────────────────────────────
    energy = feps*0.5_wp*sum(state%xs*state%psi)

    ! ── explicit (fixed-charge) gradient ─────────────────────────────────────
    allocate (force(3,nat),source=0.0_wp)
    call solvation_force_terms(ddx%params,ddx%constants,ddx%workspace,state, &
    &                          elec,force,error)
    call multipole_force_terms(ddx%params,ddx%constants,ddx%workspace,state, &
    &                          0,multipoles,force,error)
    if (error%flag /= 0) then; iostatus = 1; return; end if
    gradient(:,:) = feps*force(:,:)

    ! ── charge-response term: sum_i (dE/dq_i) * dq_i/dR ──────────────────────
    if (present(dqdr)) then
      allocate (jmat(ddx%constants%ncav,nat),source=0.0_wp)
      call get_coulomb_matrix(mol%xyz,ddx%constants%ccav,jmat)
      ! dE/dq_i = 1/2 feps ( -[J^T zeta]_i + sqrt(4pi) xs_0i )
      allocate (dedq(nat))
      dedq = 0.5_wp*feps*(-matmul(transpose(jmat),state%zeta)+sqrt4pi*state%xs(1,:))
      do iat = 1,nat
        do jat = 1,nat
          gradient(:,iat) = gradient(:,iat)+dedq(jat)*dqdr(:,iat,jat)
        end do
      end do
    end if
#else
    iostatus = 1
    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
#endif
  end subroutine ddx_pc_engrad

#ifdef WITH_DDX
!========================================================================================!

  subroutine get_coulomb_matrix(xyz,ccav,jmat)
    !*********************************************************
    !* 1/r Coulomb matrix between atoms and cavity points.
    !*********************************************************
    real(wp),intent(in)    :: xyz(:,:)
    real(wp),intent(in)    :: ccav(:,:)
    real(wp),intent(inout) :: jmat(:,:)
    integer :: ic,jat
    real(wp) :: vec(3)
    do ic = 1,size(ccav,2)
      do jat = 1,size(xyz,2)
        vec(:) = ccav(:,ic)-xyz(:,jat)
        jmat(ic,jat) = 1.0_wp/sqrt(sum(vec**2))
      end do
    end do
  end subroutine get_coulomb_matrix
#endif

!========================================================================================!
end module crest_ddx_pc
