!=============================================================================!
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
!=============================================================================!

module molbuilder_classify
  use molbuilder_classify_type
  use molbuilder_classify_func
  implicit none
  private

  !> RE-EXEPORTS
  public :: coord_classify  !> the extended coord type
  public :: mol_ring        !> a single molecular-graph ring
  public :: setup_classify  !> setup a coord_classify from coord
  public :: atinfo_classify !> add atinfo string to a coord_classify
  public :: functional_group_classify !> try to determine some functional groups

!==============================================================================!
contains  !> MODULE PROCEDURES START HERE
!==============================================================================!

end module molbuilder_classify
