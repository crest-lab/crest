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

!========================================================================================!
!========================================================================================!
module ttconf_settings_mod
!************************************************************
!* Settings container for the TTConf-light conformer search.
!*
!* This is a deliberately *minimal* module (it only depends on
!* crest_parameters) so that a ttconf_settings instance can be
!* stored directly on CREST's systemdata object without creating
!* a circular module dependency (ttconf_light_mod itself pulls in
!* crest_data via its worker routines).
!************************************************************
  use crest_parameters,only:wp
  implicit none
  private

  public :: ttconf_settings

  type :: ttconf_settings
    !> ── TT-cross sweep parameters ──────────────────────────────────────────
    integer  :: rank      = 3          !> TT rank r
    integer  :: sweeps    = 8          !> number of sweeps s
    integer  :: ngrid     = 6          !> dihedral grid points (360/ngrid degrees)
    integer  :: ninit     = 3          !> number of random initial tail seeds
    real(wp) :: ewin      = 6.0_wp     !> conformer energy window (kcal/mol)
    real(wp) :: kt        = 6.0_wp     !> maxvol energy->weight temperature (kcal/mol)
    logical  :: use_sweep = .true.     !> .false. -> brute-force oracle
    logical  :: sp_only   = .false.    !> .true. -> singlepoints only (no geometry opt)
    logical  :: use_cache = .true.     !> .false. -> disable the energy cache
    logical  :: excl_rings = .true.    !> .true. -> ring bonds are not TT variables
    logical  :: ring_sample = .false.  !> .true. -> sample ring templates as TT sites
    character(len=20) :: ring_method = 'mtd' !> ring-conformation generator to use
    !> ('mtd' = GFN-FF metadynamics; selectable, see ttconf_ringmtd_mod)
    integer  :: seed      = -1         !> RNG seed; <0 -> non-deterministic (default)
    character(len=20) :: preset = 'normal'
    character(len=:),allocatable :: userfile  !> optional dihedral spec (rigidconf fmt)
    !> user-forced TT-variable bonds (atom pairs); when set they REPLACE the
    !> automatic selection: (3, nbonds) = (atomA, atomB, npoints); npoints<=0 -> default grid
    integer,allocatable :: userbonds(:,:)
  contains
    procedure :: setpreset => ttconf_settings_setpreset
  end type ttconf_settings

!========================================================================================!
contains
!========================================================================================!

  subroutine ttconf_settings_setpreset(self,name,ok)
!************************************************************
!* Apply one of the paper's named presets to (r,s).
!*   fast     : r = s = 2
!*   normal   : r = 3, s = 8   (default)
!*   accurate : r = s = 6
!* "ok" (optional) reports whether the name was recognized.
!************************************************************
    implicit none
    class(ttconf_settings),intent(inout) :: self
    character(len=*),intent(in) :: name
    logical,intent(out),optional :: ok
    if (present(ok)) ok = .true.
    select case (trim(adjustl(name)))
    case ('fast')
      self%rank = 2; self%sweeps = 2; self%preset = 'fast'
    case ('normal','default','')
      self%rank = 3; self%sweeps = 8; self%preset = 'normal'
    case ('accurate')
      self%rank = 6; self%sweeps = 6; self%preset = 'accurate'
    case default
      if (present(ok)) ok = .false.
    end select
  end subroutine ttconf_settings_setpreset

!========================================================================================!
end module ttconf_settings_mod
!========================================================================================!
!========================================================================================!
