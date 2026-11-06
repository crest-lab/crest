!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Philipp Pracht, Christopher Zurek, Christoph Bannwarth
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
! Routines were adapted from the xtb code (github.com/grimme-lab/xtb)
! under the Open-source software LGPL-3.0 Licencse.
!================================================================================!

!========================================================================================!
!========================================================================================!
subroutine crest_rigidconf(env,tim)
!**********************************************************************
!* DEPRECATED runtype.
!*
!* "rigidconf" was a rule-based rigid-rotor conformer generator that
!* enumerated the dihedral grid and rebuilt every combination from the
!* z-matrix. That is exactly what the TTConf-light BRUTE-FORCE ORACLE
!* (crest_ttconf with use_sweep = .false.) does -- only better
!* integrated (shared classification, ring sites, topology screen,
!* multilevel optimization + CREGEN). To avoid maintaining two copies
!* of the same enumerate-and-rebuild path, rigidconf is now a thin
!* shim that forces the brute-force oracle and forwards to crest_ttconf.
!*
!* Input:
!*    env  - CREST's systemdata
!*    tim  - CREST's timer object
!**********************************************************************
  use crest_parameters
  use crest_data
  implicit none
  !> INPUT/OUTPUT
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim

!========================================================================================!
  call this_header()

! ── force the TTConf brute-force oracle (full grid enumeration) ─────────────
  env%ttconf%use_sweep = .false.

! ── carry over the optional user dihedral spec (rigidconf file format) ───────
  if (allocated(env%rigidconf_userfile)) env%ttconf%userfile = env%rigidconf_userfile

! ── hand off to the TTConf-light runtype ────────────────────────────────────
  call crest_ttconf(env,tim)

  return
!========================================================================================!
contains
!========================================================================================!
  subroutine this_header
    implicit none
    write (stdout,'(/)')
    write (stdout,'(7x,"┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓")')
    write (stdout,'(7x,"┃          R I G I D C O N F           ┃ ")')
    write (stdout,'(7x,"┃   (deprecated -> TTConf brute force) ┃ ")')
    write (stdout,'(7x,"┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛")')
    write (stdout,'(/,7x,a)') '** NOTE ** the rigidconf runtype is deprecated.'
    write (stdout,'(7x,a)')   'It now redirects to the TTConf-light brute-force oracle'
    write (stdout,'(7x,a,/)') '(equivalent: runtype = "ttconf" with bruteforce = true).'
  end subroutine this_header
end subroutine crest_rigidconf
!========================================================================================!
!========================================================================================!
