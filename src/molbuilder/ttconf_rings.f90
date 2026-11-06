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
module ttconf_rings_mod
!****************************************************************************************
!* Ring-conformation sampling for TTConf-light.
!*
!* A flexible ring is treated as just another tensor-train SITE: its grid is a
!* discrete set of pre-computed ring conformations ("templates", e.g. chair/boat/
!* twist puckers) rather than an evenly-spaced torsion angle. The TT machinery
!* (sweep / maxvol / cache / topology screen / optimization) is completely agnostic
!* to this; it only needs, per site, a grid size site_ngrid(site) and a way to turn a
!* chosen grid index into geometry.
!*
!* This module realizes that with three pieces:
!*   1. ttconf_ring_site / ttconf_ringset  - the per-ring template store, indexed by
!*      the TT site number it occupies.
!*   2. ttconf_ring_templates              - produces a ring's template z-matrix
!*      blocks from the registered conformer source (ttconf_ringsample_mod);
!*      with none registered it returns the input ring unchanged.
!*   3. ttconf_build_ringsites             - appends one TT site per flexible ring,
!*      growing the site_ngrid/site_step/site_bond site arrays, and
!*      ttconf_apply_active / %apply       - overlays the selected template onto a
!*      candidate z-matrix (called from ttconf_gen_check, right between
!*      construct_new_zmat and the Cartesian rebuild).
!*
!* The "active" ring set is held at module scope and consulted by the geometry
!* builder, which keeps the torsion code path (and its many call signatures)
!* untouched: with no ring set registered, %apply is a no-op.
!****************************************************************************************
  use crest_parameters,only:wp,stdout
  use strucrd,only:i2e
  use molbuilder_classify,only:coord_classify,mol_ring
  use ttconf_ringsample_mod,only:ttconf_ring_templates,ttconf_ring_flexible
  implicit none
  private

  !> one flexible ring, mapped onto a tensor-train site
  type,public :: ttconf_ring_site
    integer :: site = 0              !> the TT site index this ring occupies
    integer :: nrows = 0            !> number of z-matrix rows the ring controls
    integer,allocatable :: rows(:) !> the controlled z-matrix rows (atom indices)
    integer :: ntempl = 0          !> number of template conformations (= site_ngrid(site))
    real(wp),allocatable :: ztempl(:,:,:) !> (3, nrows, ntempl) z-matrix internals
  end type ttconf_ring_site

  !> the full collection of ring sites for one molecule
  type,public :: ttconf_ringset
    integer :: n = 0
    type(ttconf_ring_site),allocatable :: r(:)
  contains
    procedure :: apply => ttconf_ringset_apply
    procedure :: clear => ttconf_ringset_clear
    procedure :: add   => ttconf_ringset_add
  end type ttconf_ringset

  !> module-scoped "active" ring set consulted by the geometry builder.
  !> Empty (n=0) by default -> %apply is a no-op (torsion-only path).
  type(ttconf_ringset),save :: active_ringset

  public :: ttconf_set_ringset    !> register the active ring set
  public :: ttconf_clear_ringset  !> deregister it (back to torsion-only)
  public :: ttconf_apply_active   !> overlay the active ring templates onto a z-matrix
  public :: ttconf_build_ringsites!> append ring sites to the TT site arrays
  public :: ttconf_print_ringsites!> short summary of the registered ring sites

!========================================================================================!
contains
!========================================================================================!

  subroutine ttconf_ringset_clear(self)
    implicit none
    class(ttconf_ringset),intent(inout) :: self
    if (allocated(self%r)) deallocate (self%r)
    self%n = 0
  end subroutine ttconf_ringset_clear

  subroutine ttconf_ringset_add(self,rs)
    !************************************************************
    !* Append one ring site to the set (deep copy).
    !************************************************************
    implicit none
    class(ttconf_ringset),intent(inout) :: self
    type(ttconf_ring_site),intent(in) :: rs
    type(ttconf_ring_site),allocatable :: tmp(:)
    if (.not.allocated(self%r)) then
      allocate (self%r(1))
      self%r(1) = rs
    else
      allocate (tmp(self%n+1))
      tmp(1:self%n) = self%r(1:self%n)
      tmp(self%n+1) = rs
      call move_alloc(tmp,self%r)
    end if
    self%n = self%n+1
  end subroutine ttconf_ringset_add

  subroutine ttconf_ringset_apply(self,combi,zmat,zrow2site)
    !**********************************************************************
    !* Overlay the selected ring template onto a candidate z-matrix. For
    !* each ring site, combi(site) selects the template whose stored
    !* internal coordinates replace the ring atoms' z-matrix rows.
    !* (GMETRY2 then rebuilds the Cartesians, unchanged.)
    !*
    !* "zrow2site" (zmat row -> torsion TT-variable index, 0 if none) protects
    !* OTHER TT variables: a controlled row may also be an independent
    !* torsion site (e.g. an exocyclic substituent that rides along with
    !* the pucker but carries its own rotatable bond). For such a row the
    !* dihedral (component 3) was just set from the torsion grid index in
    !* construct_new_zmat -- we keep it, and only overlay the template's
    !* bond length / angle (components 1:2) so the atom follows the pucker
    !* without the ring clobbering its torsion.
    !**********************************************************************
    implicit none
    class(ttconf_ringset),intent(in) :: self
    integer,intent(in) :: combi(:)
    real(wp),intent(inout) :: zmat(:,:)
    integer,intent(in),optional :: zrow2site(:)
    integer :: i,k,t,a,hi
    do i = 1,self%n
      associate (rs => self%r(i))
        if (rs%site < 1.or.rs%site > size(combi)) cycle
        t = combi(rs%site)
        if (t < 1.or.t > rs%ntempl) cycle
        do k = 1,rs%nrows
          a = rs%rows(k)
! ── keep an independent torsion's dihedral; ride along in dist/angle ─────────
          hi = 3
          if (present(zrow2site)) then
            if (a >= 1.and.a <= size(zrow2site)) then
              if (zrow2site(a) /= 0) hi = 2
            end if
          end if
          zmat(1:hi,a) = rs%ztempl(1:hi,k,t)
        end do
      end associate
    end do
  end subroutine ttconf_ringset_apply

!========================================================================================!

  subroutine ttconf_set_ringset(rs)
    implicit none
    type(ttconf_ringset),intent(in) :: rs
    call ttconf_clear_ringset()
    if (rs%n > 0) then
      allocate (active_ringset%r(rs%n))
      active_ringset%r(1:rs%n) = rs%r(1:rs%n)
    end if
    active_ringset%n = rs%n
  end subroutine ttconf_set_ringset

  subroutine ttconf_clear_ringset()
    implicit none
    call active_ringset%clear()
  end subroutine ttconf_clear_ringset

  subroutine ttconf_apply_active(combi,zmat,zrow2site)
    !************************************************************
    !* Geometry seam: overlay the module-scoped active ring set
    !* onto "zmat". A no-op when no ring set is registered.
    !* "zrow2site" (optional, zmat row -> torsion TT index) shields
    !* independent torsion sites from the ring overlay.
    !************************************************************
    implicit none
    integer,intent(in) :: combi(:)
    real(wp),intent(inout) :: zmat(:,:)
    integer,intent(in),optional :: zrow2site(:)
    if (active_ringset%n < 1) return
    call active_ringset%apply(combi,zmat,zrow2site)
  end subroutine ttconf_apply_active

!========================================================================================!

  subroutine ttconf_build_ringsites(molc,nsite,site_ngrid,site_step,site_bond,ringset,nadded)
    !**************************************************************************
    !* Append one TT site per flexible ring, growing the parallel site arrays
    !* (site_ngrid/site_step/site_bond) and nsite. Each new site's grid size is the ring's
    !* template count; its site_step is 0 and its site_bond column is (0,0) (the
    !* geometry comes from the ring template overlay, not a torsion step). The
    !* returned "ringset" holds the template store for the geometry seam.
    !*
    !* In/out:
    !*   nsite           - number of TT sites (grows by the number of ring sites)
    !*   site_ngrid,site_step  - per-site grid size / torsion step (grown)
    !*   site_bond        - (2,nsite) bond atoms per site (grown; (0,0) for rings)
    !* Out:
    !*   ringset,nadded - the ring template store and the count of rings added
    !**************************************************************************
    implicit none
    type(coord_classify),intent(in) :: molc
    integer,intent(inout) :: nsite
    integer,allocatable,intent(inout) :: site_ngrid(:)
    real(wp),allocatable,intent(inout) :: site_step(:)
    integer,allocatable,intent(inout) :: site_bond(:,:)
    type(ttconf_ringset),intent(out) :: ringset
    integer,intent(out) :: nadded
    integer :: ir,nrows,ntempl
    integer,allocatable :: rows(:)
    real(wp),allocatable :: ztempl(:,:,:)
    type(ttconf_ring_site) :: rs

    call ringset%clear()
    nadded = 0
    if (molc%nrings < 1) return

    do ir = 1,molc%nrings
! ── flexibility filter: skip rigid (small/planar/aromatic) rings ─────────────
      if (.not.ttconf_ring_flexible(molc,molc%ringlist(ir))) cycle

      call ttconf_ring_templates(molc,molc%ringlist(ir),rows,nrows,ztempl,ntempl)
      if (ntempl < 1.or.nrows < 1) cycle

! ── append a new TT site for this ring ───────────────────────────────────────
      nsite = nsite+1
      site_ngrid = [site_ngrid,ntempl]                       !> grid size = #templates
      site_step = [site_step,0.0_wp]                            !> unused for a ring site
      site_bond = reshape([site_bond,0,0],[2,nsite])         !> (0,0): not a torsion bond

! ── stash the ring's template store, keyed on its site index ─────────────────
      rs%site = nsite
      rs%nrows = nrows
      if (allocated(rs%rows)) deallocate (rs%rows)
      rs%rows = rows
      rs%ntempl = ntempl
      if (allocated(rs%ztempl)) deallocate (rs%ztempl)
      rs%ztempl = ztempl
      call ringset%add(rs)
      nadded = nadded+1

      deallocate (rows,ztempl)
    end do
  end subroutine ttconf_build_ringsites

!========================================================================================!

  subroutine ttconf_print_ringsites(molc,ringset)
    !************************************************************
    !* Short summary of the registered ring sites.
    !************************************************************
    implicit none
    type(coord_classify),intent(in) :: molc
    type(ttconf_ringset),intent(in) :: ringset
    integer :: i,k,a
    character(len=:),allocatable :: line
    character(len=16) :: tok

    if (ringset%n < 1) return
    write (stdout,'(/,1x,a,i0,a)') 'Ring sampling: ',ringset%n, &
    &  ' ring site(s) registered as TT variables'
    do i = 1,ringset%n
!> list only the heavy controlled atoms (the ring skeleton); the controlled
!> rows also include the riding hydrogens, which would clutter the summary
      line = ''
      do k = 1,ringset%r(i)%nrows
        a = ringset%r(i)%rows(k)
        if (molc%at(a) == 1) cycle
        write (tok,'(a,i0)') trim(i2e(molc%at(a))),a
        line = trim(line)//' '//trim(tok)
      end do
      write (stdout,'(3x,a,i0,a,i0,a,a)') 'site ',ringset%r(i)%site, &
      &  '  (',ringset%r(i)%ntempl,' templates):',trim(line)
    end do
  end subroutine ttconf_print_ringsites

!========================================================================================!
!========================================================================================!
end module ttconf_rings_mod
!========================================================================================!
!========================================================================================!
