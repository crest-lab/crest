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

module crest_electrostatic
!********************************************************************************
!* Charge-equilibration electrostatic component.
!*
!* Thin wrapper around the multicharge library exposing the electrostatic
!* energy, its Cartesian gradient, the atomic partial charges and their
!* geometry derivative dq/dR for the EEQ (2019) and EEQ-BC (2025) models.
!* Used standalone and as the charge source for the ddX solvation engine.
!********************************************************************************
  use crest_parameters
  use strucrd,only:coord
#ifdef WITH_TBLITE
  use mctc_env,only:error_type
  use mctc_io,only:structure_type,new
  use mctc_cutoff,only:get_lattice_points
  use multicharge_model,only:mchrg_model_type
  use multicharge_param,only:new_eeq2019_model,new_eeqbc2025_model
#endif
  implicit none
  private

  public :: electrostatic_data
  public :: electrostatic_core

!> Bundled electrostatic settings (stored on calculation_settings)
  type :: electrostatic_data
    character(len=:),allocatable :: charge_model  !> 'eeq' | 'eeqbc'
  end type electrostatic_data

contains
!========================================================================================!

  subroutine electrostatic_core(mol,chrg,model,energy,gradient,qat,iostatus,dqdr)
    !***********************************************************************
    !* Charge-equilibration energy, gradient and atomic charges.
    !*
    !*  mol      : molecular structure (Bohr)
    !*  chrg     : total molecular charge
    !*  model    : 'eeq' (2019) | 'eeqbc' (2025)
    !*  energy   : electrostatic energy (out)
    !*  gradient : 3,nat Cartesian gradient (out, overwritten)
    !*  qat      : atomic partial charges (out)
    !*  dqdr     : optional dq_i/dR_a (3,nat,nat) charge derivative (out)
    !***********************************************************************
    type(coord),intent(in)        :: mol
    integer,intent(in)            :: chrg
    character(len=*),intent(in)   :: model
    real(wp),intent(out)          :: energy
    real(wp),intent(out)          :: gradient(:,:)
    real(wp),intent(out)          :: qat(:)
    integer,intent(out)           :: iostatus
    real(wp),intent(out),optional :: dqdr(:,:,:)
#ifdef WITH_TBLITE
    type(structure_type) :: struc
    class(mchrg_model_type),allocatable :: mchrg
    type(error_type),allocatable :: error
    real(wp),allocatable :: cn(:),dcndr(:,:,:),dcndL(:,:,:)
    real(wp),allocatable :: qloc(:),dqlocdr(:,:,:),dqlocdL(:,:,:)
    real(wp),allocatable :: trans(:,:),ener(:),sigma(:,:),dqdL(:,:,:)
    integer :: nat

    iostatus = 0
    nat = mol%nat
    energy = 0.0_wp
    gradient(:,:) = 0.0_wp

    call new(struc,mol%at,mol%xyz,charge=real(chrg,wp))

    select case (trim(model))
    case ('eeq')
      call new_eeq2019_model(struc,mchrg,error)
    case ('eeqbc')
      call new_eeqbc2025_model(struc,mchrg,error)
    case default
      iostatus = 1; return
    end select
    if (allocated(error)) then; iostatus = 1; return; end if

    ! ── coordination numbers + local charges (and their derivatives) ─────────
    allocate (cn(nat),qloc(nat))
    allocate (dcndr(3,nat,nat),dcndL(3,3,nat))
    allocate (dqlocdr(3,nat,nat),dqlocdL(3,3,nat))
    call get_lattice_points(struc%periodic,struc%lattice,mchrg%ncoord%cutoff,trans)
    call mchrg%ncoord%get_coordination_number(struc,trans,cn,dcndr,dcndL)
    call mchrg%local_charge(struc,trans,qloc,dqlocdr,dqlocdL)

    ! ── solve the equilibration: energy + gradient + charges (+ dq/dR) ───────
    allocate (ener(nat),source=0.0_wp)
    allocate (sigma(3,3),source=0.0_wp)
    !> branch on dqdr presence: forwarding a non-present optional to solve's
    !> contiguous dummy segfaults under ifx.  dqdr is only evaluated by solve
    !> when dqdL is present as well, so it has to be passed too.
    if (present(dqdr)) then
      allocate (dqdL(3,3,nat),source=0.0_wp)
      call mchrg%solve(struc,error,cn,qloc,dcndr,dcndL,dqlocdr,dqlocdL, &
      &                energy=ener,gradient=gradient,sigma=sigma,qvec=qat, &
      &                dqdr=dqdr,dqdL=dqdL)
    else
      call mchrg%solve(struc,error,cn,qloc,dcndr,dcndL,dqlocdr,dqlocdL, &
      &                energy=ener,gradient=gradient,sigma=sigma,qvec=qat)
    end if
    if (allocated(error)) then; iostatus = 1; return; end if
    energy = sum(ener)
#else
    iostatus = 1
    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
    qat(:) = 0.0_wp
    if (present(dqdr)) dqdr(:,:,:) = 0.0_wp
#endif
  end subroutine electrostatic_core

!========================================================================================!
end module crest_electrostatic
