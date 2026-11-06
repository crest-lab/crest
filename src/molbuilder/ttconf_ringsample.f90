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
module ttconf_ringsample_mod
!****************************************************************************************
!* Ring-conformation SAMPLING for TTConf-light.
!*
!* This module answers two questions for a perceived ring:
!*   1. Is the ring flexible enough to be worth treating as a TT site?
!*      -> ttconf_ring_flexible (size + planarity filter)
!*   2. What template conformations should that site offer, and as which
!*      z-matrix internal coordinates?
!*      -> ttconf_ring_templates
!*
!* CENTRAL DESIGN POINT -- "re-measure, do not copy".
!* A ring atom's z-matrix row (dist,angle,dihedral) is meaningful only relative
!* to ITS OWN reference atoms zmap(a,1:3) in THIS molecule. A ring conformation
!* obtained elsewhere (a mirror image; a separately sampled cut-out) therefore
!* cannot be transplanted by copying internal coordinates -- only the geometry
!* transfers. So the provider produces each template as Cartesian coordinates
!* and then RE-MEASURES the internals (dist/angle/dihedral) of every ring atom
!* against the molecule's existing zmap references. Overlaying those internals
!* (ttconf_rings_mod) and rebuilding via GMETRY2 then reproduces the template's
!* ring shape exactly, while every substituent rides along on its own unchanged
!* internals (axial/equatorial follows the pucker automatically).
!*
!* CONFORMER SOURCE. An external conformer source is plugged in via a
!* registerable hook (ttconf_set_ring_sampler). This module owns the
!* level-agnostic geometry work (the cut-out, the fold-back, the re-measure)
!* while the actual GFN-FF MD/MTD + optimization + clustering (which needs env
!* and the high-level algos routines) is supplied by the algos layer and
!* registered before sampling. With no hook registered, the provider returns
!* the identity: the single template is the input ring geometry, passed through
!* the same re-measure/overlay/rebuild path, reproducing the input exactly.
!* (A mirror-through-the-mean-plane stand-in does NOT work: a reflection is an
!* IMPROPER operation, so overwriting only the ring rows while substituents ride
!* along on their fixed local internals breaks them; only a PROPER pucker change,
!* i.e. a real sampled conformer, transplants cleanly.) Keeping the source behind
!* the hook leaves molbuilder free of algos dependencies.
!****************************************************************************************
  use crest_parameters,only:wp,stdout,aatoau
  use strucrd,only:coord
  use molbuilder_classify,only:coord_classify,mol_ring
  implicit none
  private

  !> ring smaller than this is treated as rigid (never a TT site)
  integer,parameter :: min_ring_size = 4
  !> rings larger than this are always flexible
  integer,parameter :: large_ring_size = 6
  !> mean-plane RMSD (Angstrom) above which a small ring counts as non-planar
  real(wp),parameter :: planar_rmsd_thr = 0.1_wp
  !> two templates closer than this (max |dz| over all internals) are duplicates
  real(wp),parameter :: templ_dedup_thr = 1.0e-4_wp

  !> registerable conformer source: given an isolated cut-out molecule,
  !> return nconf candidate Cartesian geometries (3, frag%nat, nconf).
  abstract interface
    subroutine ttconf_ring_sampler_i(frag,confs,nconf)
      import :: wp,coord
      type(coord),intent(in) :: frag
      real(wp),allocatable,intent(out) :: confs(:,:,:)
      integer,intent(out) :: nconf
    end subroutine ttconf_ring_sampler_i
  end interface
  procedure(ttconf_ring_sampler_i),pointer,save :: ring_sampler => null()

  public :: ttconf_ring_flexible
  public :: ttconf_set_ring_sampler   !> register the conformer source
  public :: ttconf_clear_ring_sampler !> deregister it (back to identity)
  public :: ttconf_ring_templates

!========================================================================================!
contains
!========================================================================================!

  subroutine ttconf_set_ring_sampler(proc)
    !************************************************************
    !* Register the conformer source (supplied by algos).
    !************************************************************
    implicit none
    procedure(ttconf_ring_sampler_i) :: proc
    ring_sampler => proc
  end subroutine ttconf_set_ring_sampler

  subroutine ttconf_clear_ring_sampler()
    implicit none
    ring_sampler => null()
  end subroutine ttconf_clear_ring_sampler

!========================================================================================!

  logical function ttconf_ring_flexible(molc,ring) result(flex)
    !**************************************************************************
    !* Decide whether a ring is flexible enough to sample as a TT site.
    !*   size <  min_ring_size      -> rigid  (e.g. 3-membered)
    !*   size >  large_ring_size    -> flexible (macrocyclic-ish)
    !*   min..large                 -> flexible only if non-planar, judged by
    !*                                 the ring atoms' mean-plane RMSD in the
    !*                                 input geometry (aromatic/sp2 rings are
    !*                                 flat and fall below the threshold).
    !**************************************************************************
    implicit none
    type(coord_classify),intent(in) :: molc
    type(mol_ring),intent(in) :: ring
    real(wp) :: rmsd

    flex = .false.
    if (ring%size < min_ring_size) return
    if (ring%size > large_ring_size) then
      flex = .true.
      return
    end if
    rmsd = ring_plane_rmsd(molc,ring)
    flex = (rmsd > planar_rmsd_thr)
  end function ttconf_ring_flexible

!========================================================================================!

  subroutine ttconf_ring_templates(molc,ring,rows,nrows,ztempl,ntempl)
    !**************************************************************************
    !* Supply a ring's template conformations as the z-matrix internal
    !* coordinates of the atoms the ring controls.
    !*
    !* Out:
    !*   rows(nrows)            - the z-matrix rows the ring controls
    !*   ztempl(3,nrows,ntempl) - re-measured (dist,angle,dihedral) per template
    !*   ntempl                 - number of templates (= the site grid size)
    !*
    !* With no conformer source registered, ntempl = 1 (the identity, re-measured
    !* against molc%zmap). The dedup guard below collapses coincident templates,
    !* so a registered sampler can simply hand back more candidate geometries.
    !**************************************************************************
    implicit none
    type(coord_classify),intent(in) :: molc
    type(mol_ring),intent(in) :: ring
    integer,allocatable,intent(out) :: rows(:)
    integer,intent(out) :: nrows
    real(wp),allocatable,intent(out) :: ztempl(:,:,:)
    integer,intent(out) :: ntempl

    real(wp),allocatable :: conf(:,:,:)   !> (3,nat,nconf) candidate geometries
    real(wp),allocatable :: zt(:,:)       !> internals of one template
    integer :: nconf,t,maxconf

! ── candidate geometries + the z-matrix rows the ring controls ──────────────
    call get_ring_candidates(molc,ring,conf,nconf,rows,nrows)

    maxconf = max(nconf,1)
    allocate (ztempl(3,nrows,maxconf),source=0.0_wp)
    allocate (zt(3,nrows))

    ntempl = 0
    do t = 1,nconf
! ── re-measure each ring atom's internals on this candidate geometry ─────────
      call measure_ring_internals(molc,rows,nrows,conf(:,:,t),zt)
      if (ntempl > 0) then
        if (template_is_duplicate(ztempl(:,:,ntempl),zt,nrows)) cycle
      end if
      ntempl = ntempl+1
      ztempl(1:3,1:nrows,ntempl) = zt(1:3,1:nrows)
    end do
    if (ntempl < 1) ntempl = 1               !> always keep at least the identity

    deallocate (zt)
    if (allocated(conf)) deallocate (conf)
  end subroutine ttconf_ring_templates

!========================================================================================!

  subroutine get_ring_candidates(molc,ring,conf,nconf,rows,nrows)
    !**************************************************************************
    !* Conformer source for the provider, and the z-matrix rows the ring
    !* controls. Each returned conf(:,:,t) is a full-molecule-sized geometry.
    !*
    !* No hook registered: identity. rows = the ring atoms; the lone candidate
    !* is the input geometry, reproducing the ring exactly.
    !*
    !* Hook registered: build the isolated cut-out, let the registered sampler
    !* return cut-out conformers, then FOLD each one back -- the kept
    !* (cut-out) atoms take their sampled coordinates, the rest keep the input
    !* geometry. The controlled rows are EVERY cut-out atom whose three zmap
    !* references are themselves cut-out atoms: that is the ring plus its first-
    !* shell substituents (e.g. the ring hydrogens), all of which must follow a
    !* proper pucker change. Atoms reaching outside the cut-out keep their input
    !* internals and ride along. Because a controlled atom and all its references
    !* are kept, they share the cut-out's frame, so the (frame-invariant)
    !* internal-coordinate re-measure is self-consistent.
    !**************************************************************************
    implicit none
    type(coord_classify),intent(in) :: molc
    type(mol_ring),intent(in) :: ring
    real(wp),allocatable,intent(out) :: conf(:,:,:)
    integer,intent(out) :: nconf
    integer,allocatable,intent(out) :: rows(:)
    integer,intent(out) :: nrows
    type(coord) :: frag
    integer,allocatable :: o2c(:)
    real(wp),allocatable :: fconf(:,:,:)
    integer :: nat,t,a,k

    nat = molc%nat

    if (.not.associated(ring_sampler)) then
! ── no sampler: identity, ring atoms only ───────────────────────────────────
      nrows = ring%size
      allocate (rows(nrows))
      rows(:) = ring%atoms(:)
      allocate (conf(3,nat,1),source=0.0_wp)
      conf(1:3,1:nat,1) = molc%xyz(1:3,1:nat)
      nconf = 1
      return
    end if

! ── sampler registered: cut out, decide controlled rows, sample, fold back ──
    call build_ring_cutout(molc,ring,frag,o2c)

! ── controlled rows = kept atoms whose every zmap reference is also kept ────
    nrows = 0
    do a = 1,nat
      if (o2c(a) <= 0) cycle
      if (all_refs_kept(molc,o2c,a)) nrows = nrows+1
    end do
    allocate (rows(nrows))
    k = 0
    do a = 1,nat
      if (o2c(a) <= 0) cycle
      if (all_refs_kept(molc,o2c,a)) then
        k = k+1
        rows(k) = a
      end if
    end do

    call ring_sampler(frag,fconf,nconf)
    if (nconf < 1.or..not.allocated(fconf)) then
      allocate (conf(3,nat,1),source=0.0_wp)
      conf(1:3,1:nat,1) = molc%xyz(1:3,1:nat)
      nconf = 1
      return
    end if

    allocate (conf(3,nat,nconf),source=0.0_wp)
    do t = 1,nconf
      conf(1:3,1:nat,t) = molc%xyz(1:3,1:nat)         !> unsampled atoms stay put
      do a = 1,nat
        if (o2c(a) > 0) conf(1:3,a,t) = fconf(1:3,o2c(a),t)
      end do
    end do
    if (allocated(fconf)) deallocate (fconf)
  end subroutine get_ring_candidates

  logical function all_refs_kept(molc,o2c,a) result(kept)
    !************************************************************
    !* True if every (nonzero) zmap reference of atom a is a
    !* kept cut-out atom, so a's internals can be re-measured
    !* consistently from a sampled cut-out geometry.
    !************************************************************
    implicit none
    type(coord_classify),intent(in) :: molc
    integer,intent(in) :: o2c(:),a
    integer :: j,r
    kept = .true.
    do j = 1,3
      r = molc%zmap(a,j)
      if (r > 0) then
        if (o2c(r) <= 0) then
          kept = .false.
          return
        end if
      end if
    end do
  end function all_refs_kept

!========================================================================================!

  subroutine build_ring_cutout(molc,ring,frag,o2c)
    !**************************************************************************
    !* Build an isolated molecule ("frag") around a ring so it can be sampled
    !* on its own. The kept atoms are the ring atoms, every atom bonded to a
    !* ring atom (first shell), and every zmap reference of a ring atom (so the
    !* internal-coordinate re-measure always has its frame).
    !* Bonds from a kept atom to a NON-kept atom are saturated with a hydrogen
    !* placed along the original bond direction at a C-H-like distance.
    !*
    !* The kept atoms keep their input Cartesian coordinates; o2c(i) maps an
    !* original atom index i to its index in "frag" (0 if i is not kept), which
    !* the conformer re-measure uses to fold sampled geometries back in.
    !**************************************************************************
    implicit none
    type(coord_classify),intent(in) :: molc
    type(mol_ring),intent(in) :: ring
    type(coord),intent(out) :: frag
    integer,allocatable,intent(out) :: o2c(:)
    logical,allocatable :: keep(:)
    real(wp),allocatable :: fxyz(:,:)
    integer,allocatable :: fat(:)
    real(wp) :: u(3),d
    integer :: nat,i,a,j,r,nf,ncap
    real(wp),parameter :: rch = 1.09_wp*aatoau   !> cap X-H distance (Bohr)

    nat = molc%nat
    allocate (keep(nat),source=.false.)

! ── mark the kept atoms: ring + first shell + ring-atom zmap references ──────
    do i = 1,ring%size
      a = ring%atoms(i)
      keep(a) = .true.
      do j = 1,3
        r = molc%zmap(a,j)
        if (r > 0) keep(r) = .true.
      end do
      do j = 1,nat
        if (molc%bond(j,a) > 0) keep(j) = .true.
      end do
    end do

! ── count cap atoms: kept-atom bonds that leave the kept set ─────────────────
    ncap = 0
    do a = 1,nat
      if (.not.keep(a)) cycle
      do j = 1,nat
        if (molc%bond(j,a) > 0.and..not.keep(j)) ncap = ncap+1
      end do
    end do

    nf = count(keep)
    allocate (fxyz(3,nf+ncap),source=0.0_wp)
    allocate (fat(nf+ncap),source=0)
    allocate (o2c(nat),source=0)

! ── copy the kept atoms (keep their input coordinates) ──────────────────────
    nf = 0
    do a = 1,nat
      if (.not.keep(a)) cycle
      nf = nf+1
      o2c(a) = nf
      fat(nf) = molc%at(a)
      fxyz(1:3,nf) = molc%xyz(1:3,a)
    end do

! ── add a capping H for every bond crossing the cut ─────────────────────────
    do a = 1,nat
      if (.not.keep(a)) cycle
      do j = 1,nat
        if (molc%bond(j,a) > 0.and..not.keep(j)) then
          u = molc%xyz(1:3,j)-molc%xyz(1:3,a)
          d = sqrt(u(1)**2+u(2)**2+u(3)**2)
          if (d < 1.0e-6_wp) cycle
          nf = nf+1
          fat(nf) = 1
          fxyz(1:3,nf) = molc%xyz(1:3,a)+rch*u/d
        end if
      end do
    end do

    frag%nat = nf
    call move_alloc(fat,frag%at)
    call move_alloc(fxyz,frag%xyz)
    deallocate (keep)
  end subroutine build_ring_cutout

!========================================================================================!

  subroutine measure_ring_internals(molc,rows,nrows,xyz,zt)
    !**************************************************************************
    !* Re-measure the z-matrix internals (dist,angle,dihedral) of each ring
    !* atom against ITS molc%zmap references, but on the supplied Cartesian
    !* geometry "xyz". Missing references (zmap==0, e.g. ring atoms among the
    !* first z-matrix rows) keep the molecule's original internal value.
    !**************************************************************************
    implicit none
    type(coord_classify),intent(in) :: molc
    integer,intent(in) :: rows(:),nrows
    real(wp),intent(in) :: xyz(:,:)
    real(wp),intent(out) :: zt(:,:)
    type(coord) :: frag
    integer :: k,a,b,c,d

    frag = molc%as_coord()
    frag%xyz = xyz

    do k = 1,nrows
      a = rows(k)
      b = molc%zmap(a,1); c = molc%zmap(a,2); d = molc%zmap(a,3)
      if (b > 0) then
        zt(1,k) = frag%dist(a,b)
      else
        zt(1,k) = molc%zmat(1,a)
      end if
      if (c > 0) then
        zt(2,k) = frag%angle(a,b,c)
      else
        zt(2,k) = molc%zmat(2,a)
      end if
      if (d > 0) then
        zt(3,k) = frag%dihedral(a,b,c,d)
      else
        zt(3,k) = molc%zmat(3,a)
      end if
    end do
  end subroutine measure_ring_internals

!========================================================================================!

  logical function template_is_duplicate(za,zb,nrows) result(dup)
    !************************************************************
    !* True if two template internal-coordinate blocks coincide
    !* to within templ_dedup_thr (collapses degenerate puckers).
    !************************************************************
    implicit none
    real(wp),intent(in) :: za(:,:),zb(:,:)
    integer,intent(in) :: nrows
    dup = (maxval(abs(za(1:3,1:nrows)-zb(1:3,1:nrows))) < templ_dedup_thr)
  end function template_is_duplicate

!========================================================================================!
!> GEOMETRY HELPERS
!========================================================================================!

  subroutine ring_mean_plane(molc,ring,c0,nrm)
    !************************************************************
    !* Centroid c0 and unit normal nrm of a ring's best-fit
    !* (mean) plane. nrm is the eigenvector of the coordinate
    !* covariance with the smallest eigenvalue.
    !************************************************************
    implicit none
    type(coord_classify),intent(in) :: molc
    type(mol_ring),intent(in) :: ring
    real(wp),intent(out) :: c0(3),nrm(3)
    real(wp) :: cov(3,3),p(3)
    integer :: i,a,n

    n = ring%size
    c0 = 0.0_wp
    do i = 1,n
      c0 = c0+molc%xyz(1:3,ring%atoms(i))
    end do
    c0 = c0/real(n,wp)

    cov = 0.0_wp
    do i = 1,n
      a = ring%atoms(i)
      p = molc%xyz(1:3,a)-c0
      cov(1,1) = cov(1,1)+p(1)*p(1); cov(1,2) = cov(1,2)+p(1)*p(2)
      cov(1,3) = cov(1,3)+p(1)*p(3); cov(2,2) = cov(2,2)+p(2)*p(2)
      cov(2,3) = cov(2,3)+p(2)*p(3); cov(3,3) = cov(3,3)+p(3)*p(3)
    end do
    cov(2,1) = cov(1,2); cov(3,1) = cov(1,3); cov(3,2) = cov(2,3)

    call smallest_eigvec3(cov,nrm)
  end subroutine ring_mean_plane

  real(wp) function ring_plane_rmsd(molc,ring) result(rmsd)
    !************************************************************
    !* RMSD of the ring atoms from their mean plane (Angstrom).
    !* Equals sqrt(lambda_min/n) of the coordinate covariance.
    !************************************************************
    implicit none
    type(coord_classify),intent(in) :: molc
    type(mol_ring),intent(in) :: ring
    real(wp) :: c0(3),nrm(3),s
    integer :: i,a

    call ring_mean_plane(molc,ring,c0,nrm)
    rmsd = 0.0_wp
    do i = 1,ring%size
      a = ring%atoms(i)
      s = (molc%xyz(1,a)-c0(1))*nrm(1)+(molc%xyz(2,a)-c0(2))*nrm(2) &
      &  +(molc%xyz(3,a)-c0(3))*nrm(3)
      rmsd = rmsd+s*s
    end do
    rmsd = sqrt(rmsd/real(ring%size,wp))
  end function ring_plane_rmsd

  subroutine smallest_eigvec3(a,vec)
    !************************************************************
    !* Unit eigenvector of the smallest eigenvalue of a real
    !* symmetric 3x3 matrix. Eigenvalues come from the analytic
    !* (trigonometric) solution of the characteristic cubic; the
    !* eigenvector is the null space of (A - lambda I), found as
    !* the largest cross product of its rows.
    !************************************************************
    implicit none
    real(wp),intent(in) :: a(3,3)
    real(wp),intent(out) :: vec(3)
    real(wp) :: p1,q,p2,p,r,phi,e1,e3,lam
    real(wp) :: b(3,3),r1(3),r2(3),r3(3),c12(3),c13(3),c23(3)
    real(wp) :: n12,n13,n23
    real(wp),parameter :: pi = acos(-1.0_wp)

    p1 = a(1,2)**2+a(1,3)**2+a(2,3)**2
    if (p1 <= 1.0e-14_wp) then
! ── already diagonal: smallest diagonal entry, axis-aligned vector ──────────
      lam = min(a(1,1),a(2,2),a(3,3))
    else
      q = (a(1,1)+a(2,2)+a(3,3))/3.0_wp
      p2 = (a(1,1)-q)**2+(a(2,2)-q)**2+(a(3,3)-q)**2+2.0_wp*p1
      p = sqrt(p2/6.0_wp)
      b = a
      b(1,1) = b(1,1)-q; b(2,2) = b(2,2)-q; b(3,3) = b(3,3)-q
      b = b/p
      r = det3(b)/2.0_wp
      if (r <= -1.0_wp) then
        phi = pi/3.0_wp
      else if (r >= 1.0_wp) then
        phi = 0.0_wp
      else
        phi = acos(r)/3.0_wp
      end if
      e1 = q+2.0_wp*p*cos(phi)                 !> largest eigenvalue
      e3 = q+2.0_wp*p*cos(phi+2.0_wp*pi/3.0_wp) !> smallest eigenvalue
      lam = e3
    end if

! ── null space of (A - lam I): pick the best-conditioned row cross product ──
    b = a
    b(1,1) = b(1,1)-lam; b(2,2) = b(2,2)-lam; b(3,3) = b(3,3)-lam
    r1 = b(1,1:3); r2 = b(2,1:3); r3 = b(3,1:3)
    c12 = cross(r1,r2); c13 = cross(r1,r3); c23 = cross(r2,r3)
    n12 = dot(c12,c12); n13 = dot(c13,c13); n23 = dot(c23,c23)
    if (n12 >= n13.and.n12 >= n23) then
      vec = c12
    else if (n13 >= n23) then
      vec = c13
    else
      vec = c23
    end if
    p = sqrt(dot(vec,vec))
    if (p > 1.0e-12_wp) then
      vec = vec/p
    else
      vec = [0.0_wp,0.0_wp,1.0_wp]
    end if
  end subroutine smallest_eigvec3

  pure real(wp) function det3(a) result(d)
    implicit none
    real(wp),intent(in) :: a(3,3)
    d = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
    & -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
    & +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
  end function det3

  pure function cross(u,v) result(w)
    implicit none
    real(wp),intent(in) :: u(3),v(3)
    real(wp) :: w(3)
    w(1) = u(2)*v(3)-u(3)*v(2)
    w(2) = u(3)*v(1)-u(1)*v(3)
    w(3) = u(1)*v(2)-u(2)*v(1)
  end function cross

  pure real(wp) function dot(u,v) result(s)
    implicit none
    real(wp),intent(in) :: u(3),v(3)
    s = u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
  end function dot

!========================================================================================!
!========================================================================================!
end module ttconf_ringsample_mod
!========================================================================================!
!========================================================================================!
