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
module ttconf_light_mod
!************************************************************
!* Helper routines for the TTConf-light conformer search.
!* Kept in a module so that the coord_classify object (which
!* carries allocatable components) is passed with an explicit
!* interface.
!************************************************************
  use crest_parameters
  use molbuilder_classify_type,only:coord_classify
  use ttconf_rings_mod,only:ttconf_apply_active
  implicit none
  private

  public :: ttconf_identify_dihedrals
  public :: ttconf_force_userbonds
  public :: ttconf_print_dihedrals
  public :: ttconf_bruteforce_generate
  public :: ttconf_sweep
  public :: ttconf_print_settings

  !> energy book-keeping for infeasible (topology-broken) candidates
  real(wp),parameter :: ttl_einf = 1.0e6_wp   !> "infinity" energy
  real(wp),parameter :: ttl_feas = 1.0e5_wp   !> feasibility cutoff (E < ttl_feas)
  !> covalent-CN adjacency threshold
  real(wp),parameter :: ttl_athr = 0.5_wp

  !> ── energy cache ──────────────────────────────────────────────────────────
  !> Open-addressing hash map (int64 packed grid combination -> optimized
  !> energy). The TT-cross sweeps revisit the same grid combinations many
  !* times (overlapping head/tail sets, forward+backward passes); caching the
  !> optimized energy lets us skip the (dominant) re-optimization on revisits.
  !> Infeasible (topology-broken) combinations are cached as ttl_einf, so their
  !> structure generation + screening is skipped too.
  type :: ttconf_ecache
    integer :: n = 0               !> number of stored entries
    integer :: cap = 0             !> capacity (power of two)
    logical :: active = .false.    !> caching enabled (no key overflow yet)
    integer(int64),allocatable :: keys(:)
    real(wp),allocatable :: ener(:)
    logical,allocatable :: used(:)
    integer :: hits = 0            !> number of reused evaluations
    integer :: opt = 0             !> number of genuine optimizations cached
  contains
    procedure :: init    => ecache_init
    procedure :: lookup  => ecache_lookup
    procedure :: insert  => ecache_insert
    procedure :: destroy => ecache_destroy
  end type ttconf_ecache

!========================================================================================!
contains
!========================================================================================!

  subroutine ttconf_print_settings(env)
!************************************************************
!* Print a settings block for the TTConf-light run, in the
!* style of CREST's other runtype headers. Includes the TT
!* sampling parameters AND the basic global-optimization /
!* energy settings, which always matter for conformer search.
!************************************************************
    use crest_data,only:systemdata,optlevflag
    implicit none
    type(systemdata),intent(inout) :: env
    character(len=52) :: hl
    character(len=*),parameter :: cs = '(2x,"┃",1x,a26," :  ",a,t66,"┃")'
    character(len=*),parameter :: ci = '(2x,"┃",1x,a26," :  ",i0,t66,"┃")'
    real(wp) :: degstep
    integer :: n

    associate (t => env%ttconf)
      write (stdout,*)
      n = max(0,(60-21)/2)
      hl = repeat(' ',n)//'TTConf-light settings'
      write (stdout,'(2x,"┏",60("━"),"┓")')
      write (stdout,'(2x,"┃",a,t66,"┃")') hl
      write (stdout,'(2x,"┣",60("━"),"┫")')

      degstep = 360.0_wp/real(t%ngrid,wp)
      write (stdout,'(2x,"┃",1x,a,t66,"┃")') '> Conformer sampling'
      if (t%use_sweep) then
!>--- sweep mode: rank/sweeps/seeds/temperature (and the preset) all apply
        write (stdout,cs) 'sampling mode','TT-cross sweep'
        write (stdout,cs) 'preset',trim(t%preset)
        write (stdout,ci) 'TT rank (r)',t%rank
        write (stdout,ci) 'number of sweeps (s)',t%sweeps
        write (stdout,ci) 'initial tail seeds',t%ninit
        write (stdout,'(2x,"┃",1x,a26," :  ",f6.2," kcal/mol",t66,"┃")') 'maxvol temperature',t%kt
        write (stdout,cs) 'energy cache',trim(merge('on ','off',t%use_cache))
        write (stdout,cs) 'in-ring bonds', &
        &  trim(merge('excluded','TT-var  ',t%excl_rings))
        if (t%seed >= 0) then
          write (stdout,ci) 'RNG seed (fixed)',t%seed
        else
          write (stdout,cs) 'RNG seed','random'
        end if
      else
!>--- brute-force mode: only the grid and the energy window matter; the
!>--- sweep parameters (rank/sweeps/seeds/temperature) are not used
        write (stdout,cs) 'sampling mode','brute-force (oracle)'
      end if
      write (stdout,'(2x,"┃",1x,a26," :  ",i0,"  (",f5.1," deg steps)",t66,"┃")') &
      &  'dihedral grid points',t%ngrid,degstep
      write (stdout,'(2x,"┃",1x,a26," :  ",f6.2," kcal/mol",t66,"┃")') 'energy window (ewin)',t%ewin

    end associate

    write (stdout,'(2x,"┠",a,"┨")') repeat("─",60)
    write (stdout,'(2x,"┃",1x,a,t66,"┃")') '> Geometry optimization settings'
    if (env%ttconf%sp_only) then
      write (stdout,cs) 'grid evaluation','singlepoint only (no opt)'
    end if
    write (stdout,cs) 'optimization level',trim(optlevflag(env%optlev))
    if (.not.env%ttconf%sp_only.and.associated(env%calc)) then
      block
        use optimize_utils,only:get_optthr
        real(wp) :: ethr,gthr
        integer :: iolev
        iolev = nint(env%optlev)
        call get_optthr(env%ref%nat,iolev,env%calc,ethr,gthr)
        write (stdout,'(2x,"┃",1x,a26," :  ",es11.4," Eh",t66,"┃")') 'energy convergence',ethr
        write (stdout,'(2x,"┃",1x,a26," :  ",es11.4," Eh/a0",t66,"┃")') 'gradient convergence',gthr
      end block
    end if
    write (stdout,ci) 'OMP threads',env%threads

    write (stdout,'(2x,"┗",a,"┛")') repeat("━",60)
    write (stdout,*)
    if (associated(env%calc)) call env%calc%info(stdout)
  end subroutine ttconf_print_settings

!========================================================================================!
  subroutine ttconf_identify_dihedrals(molc,wbo,ttbond,ttwbo,ttverdict,nttvar,exclrings)
!************************************************************
!* Identify which of the candidate rotatable bonds (one per
!* unique single bond carrying a Z-matrix dihedral) actually
!* span the conformational space, following TTConf §3.3:
!*  - exclude STIFF bonds   (WBO >= 1.1: double/aromatic/amide)
!*  - exclude IN-RING bonds: both endpoints lie on a common
!*    ring (rotating them would distort/break the ring); these
!*    belong to a separate in-ring class handled elsewhere.
!*    Controlled by "exclrings" (optional, default .true.).
!*  - exclude ROTAMER bonds: a tetrahedral end atom carrying
!*    THREE neighbours of equal canonical priority (molc%prio),
!*    i.e. methyl/CF3/tert-butyl, whose rotation only permutes
!*    equivalent atoms. (A 2-fold case such as a para-symmetric
!*    phenyl is kept, matching the paper.)
!* Hydroxyl-type terminal rotors (a single substituent) are
!* kept, as in the paper.
!*
!* Output:
!*   ttbond(2,ndieder) - the two atoms of each candidate bond
!*   ttwbo(ndieder)    - WBO of each candidate bond
!*   ttverdict(ndieder)- 0=TT-variable, 1=stiff, 2=rotamer, 3=in-ring
!*   nttvar            - number of bonds with verdict 0
!*
!************************************************************
    implicit none
    !> INPUT
    type(coord_classify),intent(in) :: molc
    real(wp),intent(in) :: wbo(molc%nat,molc%nat)
    logical,intent(in),optional :: exclrings
    !> OUTPUT
    integer,allocatable,intent(out)  :: ttbond(:,:)
    real(wp),allocatable,intent(out) :: ttwbo(:)
    integer,allocatable,intent(out)  :: ttverdict(:)
    integer,intent(out) :: nttvar
    !> LOCAL
    integer :: nd,i,k,b,c
    logical :: dropring
    real(wp),parameter :: stiffthr = 1.1_wp   !> WBO threshold for "stiff" bonds

    dropring = .true.
    if (present(exclrings)) dropring = exclrings

    nd = molc%ndieder
    nttvar = 0
    if (nd < 1) then
      allocate (ttbond(2,0),ttwbo(0),ttverdict(0))
      return
    end if

    allocate (ttbond(2,nd),source=0)
    allocate (ttwbo(nd),source=0.0_wp)
    allocate (ttverdict(nd),source=0)

! ── find a representative Z-matrix row for each unique bond index ────────────
    do i = 1,molc%nat
      k = molc%ztod(i)
      if (k < 1) cycle
      if (ttbond(1,k) /= 0) cycle           !> bond k already characterized
      b = molc%zmap(i,1)                     !> dihedral i: atoms (i)-(b)-(c)-(d)
      c = molc%zmap(i,2)                     !> rotatable bond is b-c
      ttbond(1,k) = b
      ttbond(2,k) = c
      ttwbo(k) = wbo(b,c)
    end do

! ── classify each candidate bond ────────────────────────────────────────────
    do k = 1,nd
      b = ttbond(1,k)
      c = ttbond(2,k)
      if (b == 0.or.c == 0) then
        ttverdict(k) = 1                     !> degenerate, treat as stiff/skip
        cycle
      end if
      if (molc%term(b).or.molc%term(c)) then
        ttverdict(k) = 1                     !> improper: a terminal atom (one graph
        cycle                                !> partner, e.g. F-C) can't be a rotatable bond
      end if
      if (ttwbo(k) >= stiffthr) then
        ttverdict(k) = 1                     !> stiff (double/aromatic/amide)
      else if (dropring.and.ttconf_bond_in_ring(molc,b,c)) then
        ttverdict(k) = 3                     !> in-ring bond (separate ring class)
      else if (rotor_side_symmetric(molc,b,c).or. &
      &        rotor_side_symmetric(molc,c,b)) then
        ttverdict(k) = 2                     !> pure rotamer (methyl/CF3/...)
      else
        ttverdict(k) = 0                     !> genuine TT variable
        nttvar = nttvar+1
      end if
    end do

    return
  end subroutine ttconf_identify_dihedrals

!========================================================================================!
  subroutine ttconf_force_userbonds(molc,userbonds,ttbond,ttverdict,nttvar,bondngrid)
!**********************************************************************
!* REPLACE the automatic TT-variable selection with a user-specified
!* set of atom pairs. When the user lists bonds in [ttconf], only those
!* bonds (if they are indeed rotatable z-matrix bonds) are treated as
!* TT variables -- every bond the classifier had picked is dropped.
!*
!* For each user pair (A,B) the matching candidate z-matrix bond is
!* located in "ttbond"; if found it is set to verdict 0 (TT variable),
!* regardless of its automatic class. A pair that is not a rotatable
!* z-matrix bond is warned about and skipped ("if they are indeed in
!* the zmat"). The automatically selected bonds that are NOT re-listed
!* are demoted to verdict 4 (user-excluded).
!*
!* In/out:
!*   ttverdict - per candidate bond class (rewritten: 0 for user bonds,
!*               4 for the dropped automatic ones)
!*   nttvar    - number of TT variables (= valid user bonds)
!*   bondngrid - per candidate bond grid-point override (0 = default);
!*               set from the optional third entry of a user bond.
!*
!* Input:
!*   userbonds(3,:) - (atomA, atomB, npoints) as parsed from [ttconf] bonds
!**********************************************************************
    use strucrd,only:i2e
    implicit none
    type(coord_classify),intent(in) :: molc
    integer,intent(in)    :: userbonds(:,:)
    integer,intent(in)    :: ttbond(:,:)
    integer,intent(inout) :: ttverdict(:)
    integer,intent(inout) :: nttvar
    integer,intent(inout) :: bondngrid(:)
    integer :: u,k,a,b,np,kk,nd,origv
    character(len=10) :: ta,tb
    character(len=*),parameter :: cls(0:4) = &
    &  [character(len=9) :: 'TT-var','stiff','rotamer','in-ring','user-excl']

    nd = size(ttverdict)
    if (size(userbonds,2) < 1) return
    write (stdout,'(/,1x,a)') &
    &  'User-specified TT-variable bonds (replacing the automatic selection):'

! ── drop the automatic selection (verdict 0 -> 4, "user-excluded") ──────────
    do k = 1,nd
      if (ttverdict(k) == 0) ttverdict(k) = 4
    end do
    nttvar = 0

    do u = 1,size(userbonds,2)
      a = userbonds(1,u)
      b = userbonds(2,u)
      np = userbonds(3,u)
      write (ta,'(a,i0)') trim(i2e(label(molc,a))),a
      write (tb,'(a,i0)') trim(i2e(label(molc,b))),b

! ── locate the candidate bond matching the (A,B) pair (either order) ────────
      kk = 0
      do k = 1,nd
        if ((ttbond(1,k) == a.and.ttbond(2,k) == b).or. &
        &   (ttbond(1,k) == b.and.ttbond(2,k) == a)) then
          kk = k; exit
        end if
      end do

      if (kk == 0) then
        write (stdout,'(3x,a,a,1x,a)') trim(ta),'-'//trim(tb), &
        &  '** not a rotatable z-matrix bond -> ignored'
        cycle
      end if

! ── select it as TT variable (origv = 4 means it was auto-selected) ─────────
      origv = ttverdict(kk)
      ttverdict(kk) = 0
      nttvar = nttvar+1
      if (np > 0) bondngrid(kk) = np
      if (origv == 4) then
        write (stdout,'(3x,a,a,1x,a)') trim(ta),'-'//trim(tb), &
        &  'selected (was auto-detected)'
      else
        write (stdout,'(3x,a,a,1x,a,a,a)') trim(ta),'-'//trim(tb), &
        &  'selected (overriding class ',trim(cls(origv)),')'
      end if
    end do
    return
  contains
    integer function label(m,idx) result(z)
      type(coord_classify),intent(in) :: m
      integer,intent(in) :: idx
      z = 0
      if (idx >= 1.and.idx <= m%nat) z = m%at(idx)
    end function label
  end subroutine ttconf_force_userbonds

!========================================================================================!
! ── is the bond b-c a ring bond? (both atoms share a common ring) ───────────
  logical function ttconf_bond_in_ring(molc,b,c) result(onring)
!************************************************************
!* True iff a single ring of the molecule's ring library
!* contains BOTH b and c. Requiring a *common* ring keeps
!* inter-ring linkers (e.g. the biphenyl pivot, whose two
!* carbons sit in different rings) classified as genuine
!* rotatable bonds rather than in-ring bonds.
!************************************************************
    implicit none
    type(coord_classify),intent(in) :: molc
    integer,intent(in) :: b,c
    integer :: r
    onring = .false.
    if (molc%nrings < 1.or..not.allocated(molc%ringlist)) return
    do r = 1,molc%nrings
      if (any(molc%ringlist(r)%atoms == b).and. &
      &   any(molc%ringlist(r)%atoms == c)) then
        onring = .true.
        return
      end if
    end do
  end function ttconf_bond_in_ring

!========================================================================================!
! ── is rotation about atom-partner symmetric on the "atom" side? ────────────
  logical function rotor_side_symmetric(molc,atom,partner) result(sym)
    implicit none
    type(coord_classify),intent(in) :: molc
    integer,intent(in) :: atom,partner
    integer :: j,cnt,p0
    sym = .false.
    cnt = 0
    p0 = -1
    do j = 1,molc%nat
      if (j == atom.or.j == partner) cycle
      if (molc%bond(j,atom) <= 0) cycle
      cnt = cnt+1
      if (p0 < 0) then
        p0 = molc%prio(j)
      else if (molc%prio(j) /= p0) then
        return                             !> substituents differ -> not symmetric
      end if
    end do
!>-- a tetrahedral atom with exactly three equal-priority substituents
!>-- (methyl/CF3/tert-butyl) only yields rotamers upon rotation.
!>-- Fewer (e.g. a single O-H substituent, or a 2-fold phenyl) is kept.
    if (cnt == 3) sym = .true.
  end function rotor_side_symmetric

!========================================================================================!
! ── per-site label for the sweep printout (torsion bond vs. ring site) ──────
  subroutine ttconf_site_label(molc,site_bond,k,ng,skind,dlabel)
    use strucrd,only:i2e
    implicit none
    type(coord_classify),intent(in) :: molc
    integer,intent(in) :: site_bond(:,:),k,ng
    character(len=*),intent(out) :: skind,dlabel
    if (site_bond(1,k) > 0.and.site_bond(2,k) > 0) then
      skind = 'dihedral'
      write (dlabel,'(a,i0,a,a,i0)') trim(i2e(molc%at(site_bond(1,k)))), &
      &  site_bond(1,k),'-',trim(i2e(molc%at(site_bond(2,k)))),site_bond(2,k)
    else
      skind = 'ring    '
!>--- show the template count explicitly so it is not confused with the
!>--- block "candidates" count (= head rank x templates x tail rank)
      if (ng == 1) then
        dlabel = '(1 template)'
      else
        write (dlabel,'(a,i0,a)') '(',ng,' templates)'
      end if
    end if
  end subroutine ttconf_site_label

!========================================================================================!
  subroutine ttconf_print_dihedrals(molc,ttbond,ttwbo,ttverdict,nttvar)
!************************************************************
!* Pretty-print the identified tensor-train dihedral
!* variables and the discarded (stiff/rotamer) bonds.
!************************************************************
    use strucrd,only:i2e
    implicit none
    type(coord_classify),intent(in) :: molc
    integer,intent(in)  :: ttbond(:,:)
    real(wp),intent(in) :: ttwbo(:)
    integer,intent(in)  :: ttverdict(:)
    integer,intent(in)  :: nttvar
    integer :: k,nd,b,c
    character(len=12) :: verdict
    character(len=8)  :: ba,ca
    !> shared column layout (right-justified), spanning the z-matrix width (63)
    character(len=*),parameter :: hfmt = '(1x,a6,2x,a10,2x,a10,2x,a12,2x,a)'
    character(len=*),parameter :: dfmt = '(1x,i6,2x,a10,2x,a10,2x,f12.4,2x,a)'
    character(len=*),parameter :: gfmt = '(1x,i6,2x,a10,2x,a10,2x,a12,2x,a)'

    nd = size(ttverdict)
    write (stdout,*)
    call smallhead('Tensor-train dihedral identification')
    write (stdout,'(1x,a,i0,a)') 'Found ',nd,' candidate rotatable bond(s):'
    write (stdout,*)
    write (stdout,hfmt) 'bond','atom 1','atom 2','WBO','classification'
    write (stdout,'(1x,a)') repeat('-',61)
    do k = 1,nd
      b = ttbond(1,k)
      c = ttbond(2,k)
      select case (ttverdict(k))
      case (0); verdict = 'TT-variable'
      case (1); verdict = 'stiff'
      case (2); verdict = 'rotamer'
      case (3); verdict = 'in-ring'
      case (4); verdict = 'user-excl'
      case default; verdict = 'skip'
      end select
      if (b > 0.and.c > 0) then
        write (ba,'(a,i0)') trim(i2e(molc%at(b))),b
        write (ca,'(a,i0)') trim(i2e(molc%at(c))),c
        write (stdout,dfmt) k,trim(ba),trim(ca),ttwbo(k),trim(verdict)
      else
        write (stdout,gfmt) k,'n/a','n/a','-',trim(verdict)
      end if
    end do
    write (stdout,'(1x,a)') repeat('-',61)
    write (stdout,'(1x,a,i0,a)') '=> ',nttvar, &
    &  ' tensor-train variable(s) span the conformational space.'
    write (stdout,*)
    return
  end subroutine ttconf_print_dihedrals

!========================================================================================!
  subroutine ttconf_bruteforce_generate(mol,molc,zrow2site,site_ngrid,site_step,nsite, &
  &                                      rawfile,ngen,ncombi)
!************************************************************
!* Brute-force oracle: enumerate ALL grid
!* combinations of the TT-variable dihedrals, build each
!* structure from the Z-matrix, and keep only those whose
!* covalent connectivity matches the reference (no QM).
!* Survivors are written (xyz, Angstrom) to "rawfile".
!*
!* Input:
!*   mol      - reference geometry (Bohr)
!*   molc     - classified reference (holds zmat, zmap)
!*   zrow2site  - zmat row -> TT-variable index (0 if none)
!*   site_ngrid  - grid points per TT variable
!*   site_step    - step size (rad) per TT variable
!*   nsite     - number of TT variables
!*   rawfile  - output ensemble file name
!*
!* Output:
!*   ngen     - number of structures with intact topology
!*   ncombi   - total number of enumerated combinations
!*
!************************************************************
    use strucrd,only:coord
    use adjacency,only:wbo2adjacency
    implicit none
    !> INPUT
    type(coord),intent(in) :: mol
    type(coord_classify),intent(in) :: molc
    integer,intent(in) :: nsite
    integer,intent(in) :: zrow2site(molc%nat)
    integer,intent(in) :: site_ngrid(nsite)
    real(wp),intent(in) :: site_step(nsite)
    character(len=*),intent(in) :: rawfile
    !> OUTPUT
    integer,intent(out) :: ngen,ncombi
    !> LOCAL
    integer :: nat,i,ich
    integer(int8),allocatable :: combi(:)
    real(wp),allocatable :: zmat_new(:,:)
    type(coord) :: newmol
    real(wp),allocatable :: cn(:),Bmat(:,:)
    integer,allocatable :: Aref(:,:),Anew(:,:)
    logical :: sane,more
    real(wp),parameter :: athr = 0.5_wp        !> adjacency threshold for CN bonds
    integer,parameter  :: maxcombi = 200000    !> brute-force safety cap

    nat = molc%nat
    ngen = 0
    ncombi = 1
    do i = 1,nsite
      ncombi = ncombi*site_ngrid(i)
    end do
    write (stdout,'(1x,i0,a,i0,a)') nsite,' dihedral(s) on a grid of ',ncombi, &
    &  ' total combinations.'
    if (ncombi > maxcombi) then
      write (stdout,'(1x,a,i0,a)') '**WARNING** ',ncombi, &
      &  ' combinations exceed the brute-force limit; skipping generation.'
      return
    end if

! ── reference covalent connectivity (CN-based, no QM) ───────────────────────
    call mol%cn_to_bond(cn,Bmat,'cov')
    call wbo2adjacency(nat,Bmat,Aref,athr)

    allocate (zmat_new(3,nat),source=0.0_wp)
    allocate (combi(nsite),source=1_int8)

    open (newunit=ich,file=rawfile)
    do
! ── build the structure for the current grid combination ────────────────────
      call construct_new_zmat(nat,molc%zmat,combi,nsite,site_ngrid,site_step,zrow2site,zmat_new)
      call ttconf_apply_active(int(combi),zmat_new,zrow2site) !> ring overlay (no-op if none)
      call reconstruct_zmat_to_mol(nat,mol%at,zmat_new, &
      &     molc%zmap(:,1),molc%zmap(:,2),molc%zmap(:,3),newmol)

! ── keep only structures whose connectivity matches the reference ───────────
      call newmol%cn_to_bond(cn,Bmat,'cov')
      call wbo2adjacency(nat,Bmat,Anew,athr)
      sane = all(Anew == Aref)
      if (sane) then
        call newmol%append(ich)
        ngen = ngen+1
      end if

! ── advance the mixed-radix odometer; stop after the last combination ───────
      call combi_next(combi,site_ngrid,nsite,more)
      if (.not.more) exit
    end do
    close (ich)

    if (allocated(zmat_new)) deallocate (zmat_new)
    if (allocated(combi)) deallocate (combi)
    if (allocated(cn)) deallocate (cn)
    if (allocated(Bmat)) deallocate (Bmat)
    if (allocated(Aref)) deallocate (Aref)
    if (allocated(Anew)) deallocate (Anew)
    return
  contains
! ── increment a mixed-radix counter; more=.false. once it wraps around ──────
    subroutine combi_next(combi,site_ngrid,n,more)
      implicit none
      integer,intent(in) :: n
      integer(int8),intent(inout) :: combi(n)
      integer,intent(in) :: site_ngrid(n)
      logical,intent(out) :: more
      integer :: p
      more = .true.
      p = 1
      do
        combi(p) = combi(p)+1_int8
        if (combi(p) <= int(site_ngrid(p),int8)) return
        combi(p) = 1_int8
        p = p+1
        if (p > n) then
          more = .false.
          return
        end if
      end do
    end subroutine combi_next
  end subroutine ttconf_bruteforce_generate

!========================================================================================!
!========================================================================================!
!>                            the TT-cross sweep                                         !
!========================================================================================!
!========================================================================================!

  subroutine ttconf_gen_check(mol,molc,zrow2site,site_ngrid,site_step,nsite,combi,Aref, &
  &                           sane,newmol)
!************************************************************
!* Build the structure for one grid-index combination and
!* check that its covalent connectivity matches the
!* reference (CN-based, no QM).
!************************************************************
    use strucrd,only:coord
    use adjacency,only:wbo2adjacency
    implicit none
    type(coord),intent(in) :: mol
    type(coord_classify),intent(in) :: molc
    integer,intent(in) :: nsite
    integer,intent(in) :: zrow2site(molc%nat),site_ngrid(nsite),combi(nsite)
    real(wp),intent(in) :: site_step(nsite)
    integer,intent(in) :: Aref(molc%nat,molc%nat)
    logical,intent(out) :: sane
    type(coord),intent(out) :: newmol
    integer :: nat
    integer(int8) :: c8(nsite)
    real(wp),allocatable :: zmat_new(:,:),cn(:),Bmat(:,:)
    integer,allocatable :: Anew(:,:)

    nat = molc%nat
    allocate (zmat_new(3,nat),source=0.0_wp)
    c8(:) = int(combi(:),int8)

    call construct_new_zmat(nat,molc%zmat,c8,nsite,site_ngrid,site_step,zrow2site,zmat_new)

! ── ring-site geometry seam: overlay any selected ring templates ────────────
!  No-op unless ring sites are registered (ttconf_set_ringset). Torsion sites
!  are untouched; ring sites overwrite their atoms' z-matrix rows here.
    call ttconf_apply_active(combi,zmat_new,zrow2site)

    call reconstruct_zmat_to_mol(nat,mol%at,zmat_new, &
    &     molc%zmap(:,1),molc%zmap(:,2),molc%zmap(:,3),newmol)

    call newmol%cn_to_bond(cn,Bmat,'cov')

    call wbo2adjacency(nat,Bmat,Anew,ttl_athr)
    sane = all(Anew == Aref)

    deallocate (zmat_new)
  end subroutine ttconf_gen_check

!========================================================================================!
  subroutine ttconf_eval_batch(env,mol,molc,zrow2site,site_ngrid,site_step,nsite, &
  &                            combis,ncand,Aref,energies,rawunit,nadd,cache)
!************************************************************
!* Evaluate a batch of candidate dihedral combinations:
!* generate + topology-screen, optimize the survivors in
!* parallel (crest_oloop, at env%calc), return the optimized
!* energy of every candidate (ttl_einf if infeasible) and
!* append the optimized survivors to "rawunit".
!*
!* The energy cache short-circuits combinations already seen
!* in an earlier batch/sweep: a cached candidate skips both
!* structure generation/screening and optimization, and is
!* NOT re-appended to "rawunit" (it was written on first
!* evaluation). Newly evaluated combinations (feasible or
!* infeasible) are inserted into the cache.
!************************************************************
    use crest_data,only:systemdata
    use strucrd,only:coord
    use parallel_interface,only:crest_oloop,crest_sploop
    implicit none
    type(systemdata),intent(inout) :: env
    type(coord),intent(in) :: mol
    type(coord_classify),intent(in) :: molc
    integer,intent(in) :: nsite,ncand,rawunit
    integer,intent(in) :: zrow2site(molc%nat),site_ngrid(nsite)
    real(wp),intent(in) :: site_step(nsite)
    integer,intent(in) :: combis(nsite,ncand)
    integer,intent(in) :: Aref(molc%nat,molc%nat)
    real(wp),intent(out) :: energies(ncand)
    integer,intent(out) :: nadd
    type(ttconf_ecache),intent(inout) :: cache
    integer :: nat,j,i,p,dup,nsane
    logical :: sane,hit
    type(coord) :: newmol,tmp
    type(coord),allocatable :: strucs(:)
    real(wp),allocatable :: xyzc(:,:,:),eread(:)
    integer,allocatable :: map(:),dupof(:)
    integer(int64),allocatable :: keys(:)
    logical,allocatable :: haskey(:)
    real(wp) :: ecached

    nat = molc%nat
    energies(:) = ttl_einf
    nadd = 0
    allocate (xyzc(3,nat,ncand),source=0.0_wp)
    allocate (map(ncand),source=0)
    allocate (dupof(ncand),source=0)
    allocate (keys(ncand),source=0_int64)
    allocate (haskey(ncand),source=.false.)

! ── cache lookup, then generate + topology-screen the misses ────────────────
    nsane = 0
    do j = 1,ncand
      call ttconf_packkey(combis(:,j),nsite,site_ngrid,keys(j),haskey(j))
      !haskey(j) = .false.
      if (haskey(j)) then
        hit = cache%lookup(keys(j),ecached)
        if (hit) then
          energies(j) = ecached                  !> reuse (feasible or ttl_einf)
          cycle
        end if
      end if
      call ttconf_gen_check(mol,molc,zrow2site,site_ngrid,site_step,nsite, &
      &                     combis(:,j),Aref,sane,newmol)
      if (.not.sane) then
        energies(j) = ttl_einf                    !> infeasible / topology broken
        if (haskey(j)) call cache%insert(keys(j),ttl_einf,.false.)
        cycle
      end if
! ── within-batch duplicate? (maxvol can emit repeated index sets) ───────────
      dup = 0
      if (haskey(j)) then
        do p = 1,nsane
          if (haskey(map(p)).and.keys(map(p)) == keys(j)) then
            dup = p; exit
          end if
        end do
      end if
      if (dup > 0) then
        dupof(j) = dup                            !> reuse queued slot's energy
      else
        nsane = nsane+1
        map(nsane) = j
        xyzc(:,:,nsane) = newmol%xyz
      end if
    end do

! ── optimize the unique survivors in parallel and harvest energies ─────────
    if (nsane > 0) then
      allocate (eread(nsane),source=0.0_wp)
! ── marshal the unique survivors into a coord list (canonical API) ──────────
      allocate (strucs(nsane))
      do i = 1,nsane
        strucs(i)%nat = nat
        strucs(i)%at = mol%at
        strucs(i)%xyz = xyzc(:,:,i)
      end do
      if (env%ttconf%sp_only) then
        call crest_sploop(env,nsane,strucs,eread,silent=.true.)
      else
        call crest_oloop(env,nsane,strucs,.false.,silent=.true.,eread=eread)
      end if
      do i = 1,nsane
        j = map(i)
        energies(j) = eread(i)
        if (haskey(j)) call cache%insert(keys(j),eread(i),.true.)
        !> harvest the (possibly optimized) geometry back
        xyzc(:,:,i) = strucs(i)%xyz
        call tmp%deallocate()
        tmp%nat = nat
        allocate (tmp%at(nat)); tmp%at = mol%at
        allocate (tmp%xyz(3,nat)); tmp%xyz = xyzc(:,:,i)
        tmp%energy = eread(i)
        call tmp%append(rawunit)
        nadd = nadd+1
      end do
      deallocate (strucs)
! ── propagate energies to within-batch duplicates ──────────────────────────
      do j = 1,ncand
        if (dupof(j) > 0) energies(j) = eread(dupof(j))
      end do
      deallocate (eread)
    end if
    deallocate (xyzc,map,dupof,keys,haskey)
  end subroutine ttconf_eval_batch

!========================================================================================!
  subroutine ttconf_init_tails(mol,molc,zrow2site,site_ngrid,site_step,nsite,ngrid,rmax, &
  &                            ntarget,Aref,Rfull,ninit)
!************************************************************
!* Seed the tail index sets with topologically reasonable
!* structures: tail #1 is the reference (all grid index 1),
!* the rest are random valid combinations (no energy eval).
!* Fills up to ntarget seeds (capped at the array bound rmax).
!************************************************************
    use strucrd,only:coord
    implicit none
    type(coord),intent(in) :: mol
    type(coord_classify),intent(in) :: molc
    integer,intent(in) :: nsite,ngrid,rmax,ntarget
    integer,intent(in) :: zrow2site(molc%nat),site_ngrid(nsite)
    real(wp),intent(in) :: site_step(nsite)
    integer,intent(in) :: Aref(molc%nat,molc%nat)
    integer,intent(out) :: Rfull(nsite,rmax)
    integer,intent(out) :: ninit
    integer :: i,attempts,cc(nsite),ngoal
    logical :: sane
    type(coord) :: newmol
    real(wp) :: x

    ngoal = max(1,min(ntarget,rmax))
    Rfull(:,:) = 1
    Rfull(:,1) = 1            !> reference structure: always valid
    ninit = 1
    attempts = 0
    do while (ninit < ngoal.and.attempts < 200*rmax)

      attempts = attempts+1
      do i = 1,nsite
        call random_number(x)
        cc(i) = min(site_ngrid(i),1+int(x*real(site_ngrid(i),wp)))   !> per-site grid size
      end do

      call ttconf_gen_check(mol,molc,zrow2site,site_ngrid,site_step,nsite,cc,Aref,sane,newmol)
      if (sane) then
        ninit = ninit+1
        Rfull(:,ninit) = cc(:)
      end if

    end do
  end subroutine ttconf_init_tails

!========================================================================================!
  subroutine ttconf_seed_rng(seed)
!************************************************************
!* Seed the intrinsic PRNG for reproducible tail seeding.
!* A non-negative "seed" fills the whole seed array with a
!* deterministic, well-spread pattern (so runs from the same
!* input reproduce exactly). A negative "seed" defers to the
!* processor's default (non-deterministic) initialization.
!************************************************************
    implicit none
    integer,intent(in) :: seed
    integer :: n,i
    integer,allocatable :: sd(:)
    if (seed < 0) then
      call random_seed()                      !> processor default (random per run)
      return
    end if
    call random_seed(size=n)
    allocate (sd(n))
    do i = 1,n
      sd(i) = seed+37*(i-1)+1                  !> spread the value across the state
    end do
    call random_seed(put=sd)
    deallocate (sd)
  end subroutine ttconf_seed_rng

!========================================================================================!
  subroutine ttconf_maxvol(Amat,m,n,isel,nsel)
!*******************************************************************************
!* Row-selecting maximum-volume submatrix (Goreinov-Tyrtyshnikov maxvol).
!* Given A (m x n, m >= n), choose n row indices "isel" such that the n x n
!* submatrix  Ahat = A(isel,:)  has (locally) maximal absolute determinant,
!* i.e. maximal "volume". This is the ONLY genuine tensor-train linear algebra
!* in TTConf-light: in the TT-cross sweep it picks, from the r*ngrid x r block
!* of Boltzmann weights at a site, the r rows (head/tail multi-indices) that
!* best represent the block. Maximal volume is the right criterion because the
!* coefficient matrix  C = A * Ahat^-1  then satisfies |C_ij| <= 1 for all i,j,
!* which bounds the cross-interpolation error and keeps the skeleton
!* well-conditioned.
!*
!* Algorithm (two stages):
!*   1. Seed:   n linearly independent rows via rank-revealing QR with column
!*              pivoting applied to A^T (LAPACK dgeqp3). The first n pivots of
!*              A^T are n well-separated rows of A -> a non-degenerate start.
!*   2. Refine: iterate the "dominance" swap. With Z = A * Ahat^-1 (m x n), the
!*              selected rows reproduce the identity and every other Z_ij is an
!*              interpolation coefficient. If max|Z_ij| = |Z(ii,jj)| > 1+tol,
!*              swapping selected row jj for row ii multiplies |det Ahat| by
!*              exactly |Z(ii,jj)| > 1, strictly growing the volume. Iterate
!*              until max|Z| <= 1+tol (maxvol/dominance condition) or maxiter.
!*   Ahat^-1 is formed by LU factorization (dgetrf) followed by inversion
!*   (dgetri); both reuse the same workspace.
!*
!* For m <= n there is nothing to select: all m rows are returned.
!*
!* In:   Amat(m,n) - the (non-negative) weight block; m >= n is the useful case
!* Out:  isel(1:nsel) - selected row indices; nsel = min(m,n)
!*******************************************************************************
    implicit none
    integer,intent(in) :: m,n
    real(wp),intent(in) :: Amat(m,n)
    integer,intent(out) :: isel(*)
    integer,intent(out) :: nsel
    real(wp),allocatable :: At(:,:),Asub(:,:),Z(:,:),tau(:),work(:)
    integer,allocatable :: jpvt(:),ipiv(:)
    integer :: i,j,iter,ii,jj,lwork,info
    real(wp) :: zmax,wq(1)
    real(wp),parameter :: tol = 1.0e-2_wp   !> dominance tolerance (swap if |Z|>1+tol)
    integer,parameter  :: maxiter = 100
    external :: dgeqp3,dgetrf,dgetri

    nsel = min(m,n)
    if (m <= n) then
      do i = 1,m                               !> nothing to choose: take all rows
        isel(i) = i
      end do
      return
    end if

! ── stage 1: seed with n independent rows via column-pivoted QR of A^T ───────
    allocate (At(n,m)); At = transpose(Amat)   !> pivot the COLUMNS of A^T = rows of A
    allocate (jpvt(m),source=0)                !> 0 = column free to be pivoted (dgeqp3)
    allocate (tau(n))                          !> Householder scalars (unused afterwards)
    lwork = -1
    call dgeqp3(n,m,At,n,jpvt,tau,wq,lwork,info)   !> workspace query: optimal lwork in wq(1)
    lwork = max(int(wq(1)),4*m+64)             !> generous floor; also reused by dgetri below
    allocate (work(lwork))
    call dgeqp3(n,m,At,n,jpvt,tau,work,lwork,info) !> actual rank-revealing QR factorization
    do i = 1,n
      isel(i) = jpvt(i)                        !> first n pivot columns = first n chosen rows
    end do
    deallocate (At,jpvt,tau)

! ── stage 2: maxvol dominance iterations ────────────────────────────────────
    allocate (Asub(n,n),Z(m,n),ipiv(n))
    do iter = 1,maxiter
      do j = 1,n
        Asub(:,j) = Amat(isel(1:n),j)          !> gather Ahat = A(isel,:)  (n x n)
      end do
      call dgetrf(n,n,Asub,n,ipiv,info)        !> LU factorization of Ahat
      if (info /= 0) exit                      !> singular -> keep current rows
      call dgetri(n,Asub,n,ipiv,work,lwork,info)   !> invert: Asub <- Ahat^-1
      if (info /= 0) exit
      Z = matmul(Amat,Asub)                    !> Z = A * Ahat^-1  (m x n coefficients)
! ── locate the most "dominant" off-basis entry |Z(ii,jj)| ───────────────────
      zmax = 1.0_wp+tol                        !> only swaps beyond the tolerance count
      ii = 0; jj = 0
      do j = 1,n
        do i = 1,m
          if (abs(Z(i,j)) > zmax) then
            zmax = abs(Z(i,j)); ii = i; jj = j
          end if
        end do
      end do
      if (ii == 0) exit                        !> max|Z| <= 1+tol -> volume is dominant
      isel(jj) = ii                            !> swap row jj for row ii (grows |det Ahat|)
    end do
    deallocate (Asub,Z,ipiv,work)
  end subroutine ttconf_maxvol

!========================================================================================!
  subroutine ttconf_sweep(env,mol,molc,zrow2site,site_ngrid,site_step,nsite,ngrid, &
  &                       site_bond,rank,nsweeps,ninit,kt,rawfile,ntot,usecache,seed)
!************************************************************
!* The TT-cross sweep. Maintains head (left) and tail (right)
!* index sets per site; sweeps left<->right, at each site
!* evaluating the r x ngrid x r candidate block, then using
!* maxvol on the energy-weighted matrix to select the new
!* head (forward) or tail (backward) index sets. All optimized
!* survivors are streamed to "rawfile" for final CREGEN.
!*
!* Input:
!*   env,mol,molc,zrow2site,site_ngrid,site_step,nsite,ngrid
!*   site_bond  - (2,nsite) the two atoms of each TT-variable bond
!*   rank     - maximal TT rank r
!*   nsweeps  - number of (alternating) sweeps s
!*   ninit    - number of random initial tail seeds (capped at rank)
!*   kt       - maxvol energy->weight temperature (Hartree)
!*   rawfile  - ensemble file for the optimized survivors
!*   usecache - .true. to enable the energy cache (reuse repeated
!*              grid evaluations); .false. re-optimizes every candidate
!*   seed     - PRNG seed for the random tail seeding; <0 leaves the
!*              processor default (non-deterministic run-to-run)
!* Output:
!*   ntot     - number of structures written to rawfile
!************************************************************
    use crest_data,only:systemdata
    use strucrd,only:coord,i2e
    use adjacency,only:wbo2adjacency
    implicit none
    type(systemdata),intent(inout) :: env
    type(coord),intent(in) :: mol
    type(coord_classify),intent(in) :: molc
    integer,intent(in) :: nsite,ngrid,rank,nsweeps,ninit
    integer,intent(in) :: zrow2site(molc%nat),site_ngrid(nsite)
    integer,intent(in) :: site_bond(2,nsite)
    real(wp),intent(in) :: site_step(nsite),kt
    character(len=*),intent(in) :: rawfile
    logical,intent(in) :: usecache
    integer,intent(in) :: seed
    integer,intent(out) :: ntot
    integer :: nat,r,k,sweepi,hi,g,tj,jc,ncand,rL,rR,nsel,s,ridx,nr,nc
    integer :: rawunit,nadd,nfound,ng,ngmax
    character(len=24) :: dlabel
    character(len=8)  :: skind
    character(len=70) :: hdr
    character(len=70) :: sline
    logical :: fwd
    integer,allocatable :: Lfull(:,:),Rfull(:,:),Lnew(:,:),Rnew(:,:)
    integer,allocatable :: rhoL(:),rhoR(:)
    integer,allocatable :: combis(:,:),rrow(:),rcol(:),isel(:)
    real(wp),allocatable :: energies(:),Amat(:,:),cn(:),Bmat(:,:)
    integer,allocatable :: Aref(:,:)
    real(wp) :: emin
    type(ttconf_ecache) :: cache

    nat = molc%nat
    r = rank
!> largest per-site grid (ngrid for torsions, #templates for ring sites);
!> all working buffers are sized to this so any site's block fits
    ngmax = max(ngrid,maxval(site_ngrid))

! ── reference connectivity (computed once) ─────────────────────────────────
    call mol%cn_to_bond(cn,Bmat,'cov')
    call wbo2adjacency(nat,Bmat,Aref,ttl_athr)

    allocate (Lfull(nsite,r),Rfull(nsite,r),Lnew(nsite,r),Rnew(nsite,r),source=1)
    allocate (rhoL(0:nsite),source=1)
    allocate (rhoR(nsite+1),source=1)
    allocate (combis(nsite,r*ngmax*r),source=1)
    allocate (rrow(r*ngmax*r),rcol(r*ngmax*r),source=0)
    allocate (energies(r*ngmax*r),source=ttl_einf)
    allocate (isel(r*ngmax),source=0)

    open (newunit=rawunit,file=rawfile)
    ntot = 0

! ── initialize the energy cache (reuse repeated grid evaluations) ───────────
!  When disabled, the cache is left inactive: lookups always miss and inserts
!  are no-ops, so every candidate is re-optimized (within-batch dedup still
!  applies). ttconf_eval_batch needs a valid object either way.
    if (usecache) call cache%init()

! ── seed the PRNG (deterministic if seed >= 0, random otherwise) ────────────
    call ttconf_seed_rng(seed)

! ── initialize tail index sets ─────────────────────────────────────────────
    call ttconf_init_tails(mol,molc,zrow2site,site_ngrid,site_step,nsite,ngrid,r,ninit,Aref, &
    &                      Rfull,nfound)
    rhoR(:) = nfound
    rhoR(nsite+1) = 1
    rhoL(0) = 1

! ── sweeps ─────────────────────────────────────────────────────────────────
    do sweepi = 1,nsweeps
      fwd = (mod(sweepi,2) == 1)
      hdr = ''
      write (hdr,'(a,i0,a,i0,a)') '   sweep ',sweepi,'/',nsweeps, &
      &  trim(merge('  (forward) ','  (backward)',fwd))
!>--- one box per sweep: a double-ruled header cell with the sweep line,
!>--- a single-line separator, then the site lines build up below it
      write (stdout,'(/,1x,"╔",72("═"),"╗")')
      write (stdout,'(1x,"║ ",a," ║")') hdr
      write (stdout,'(1x,"╟",72("─"),"╢")')
      if (fwd) then
!>====== forward (left -> right): update head index sets ======================
        do k = 1,nsite
          rL = rhoL(k-1); rR = rhoR(k+1)
          ng = site_ngrid(k)                          !> this site's grid size
          jc = 0
          do tj = 1,rR
            do g = 1,ng
              do hi = 1,rL
                jc = jc+1
                if (k > 1) combis(1:k-1,jc) = Lfull(1:k-1,hi)
                combis(k,jc) = g
                if (k < nsite) combis(k+1:nsite,jc) = Rfull(k+1:nsite,tj)
                rrow(jc) = hi+(g-1)*rL
                rcol(jc) = tj
              end do
            end do
          end do

          ncand = jc
          call ttconf_site_label(molc,site_bond,k,ng,skind,dlabel)
          sline = ''
          write (sline,'(a,i0,a,i0,a,a,a,a,a,i0,a)') 'site ',k,'/',nsite, &
          &  '   ',trim(skind),'  ',trim(dlabel),'   ',ncand,' candidates'
          write (stdout,'(1x,"║ ",a," ║")') sline

          call ttconf_eval_batch(env,mol,molc,zrow2site,site_ngrid,site_step,nsite, &
          &     combis(:,1:ncand),ncand,Aref,energies(1:ncand),rawunit,nadd,cache)

          ntot = ntot+nadd
          nr = rL*ng; nc = rR
          allocate (Amat(nr,nc),source=0.0_wp)
          emin = ttl_einf
          do jc = 1,ncand
            if (energies(jc) < ttl_feas) emin = min(emin,energies(jc))
          end do
          if (emin >= ttl_feas) emin = 0.0_wp
          do jc = 1,ncand
            Amat(rrow(jc),rcol(jc)) = ttconf_weight(energies(jc),emin,kt)
          end do
          call ttconf_maxvol(Amat,nr,nc,isel,nsel)
          do s = 1,nsel
            ridx = isel(s)
            hi = mod(ridx-1,rL)+1
            g = (ridx-1)/rL+1
            if (k > 1) Lnew(1:k-1,s) = Lfull(1:k-1,hi)
            Lnew(k,s) = g
          end do
          rhoL(k) = nsel
          Lfull(1:k,1:nsel) = Lnew(1:k,1:nsel)
          deallocate (Amat)
        end do
      else
!>====== backward (right -> left): update tail index sets =====================
        do k = nsite,1,-1
          rL = rhoL(k-1); rR = rhoR(k+1)
          ng = site_ngrid(k)                          !> this site's grid size
          jc = 0
          do hi = 1,rL
            do g = 1,ng
              do tj = 1,rR
                jc = jc+1
                if (k > 1) combis(1:k-1,jc) = Lfull(1:k-1,hi)
                combis(k,jc) = g
                if (k < nsite) combis(k+1:nsite,jc) = Rfull(k+1:nsite,tj)
                rrow(jc) = g+(tj-1)*ng
                rcol(jc) = hi
              end do
            end do
          end do

          ncand = jc
          call ttconf_site_label(molc,site_bond,k,ng,skind,dlabel)
          sline = ''
          write (sline,'(a,i0,a,i0,a,a,a,a,a,i0,a)') 'site ',k,'/',nsite, &
          &  '   ',trim(skind),'  ',trim(dlabel),'   ',ncand,' candidates'
          write (stdout,'(1x,"║ ",a," ║")') sline
          call ttconf_eval_batch(env,mol,molc,zrow2site,site_ngrid,site_step,nsite, &
          &     combis(:,1:ncand),ncand,Aref,energies(1:ncand),rawunit,nadd,cache)
          ntot = ntot+nadd
          nr = ng*rR; nc = rL
          allocate (Amat(nr,nc),source=0.0_wp)
          emin = ttl_einf
          do jc = 1,ncand
            if (energies(jc) < ttl_feas) emin = min(emin,energies(jc))
          end do
          if (emin >= ttl_feas) emin = 0.0_wp
          do jc = 1,ncand
            Amat(rrow(jc),rcol(jc)) = ttconf_weight(energies(jc),emin,kt)
          end do
          call ttconf_maxvol(Amat,nr,nc,isel,nsel)
          do s = 1,nsel
            ridx = isel(s)
            g = mod(ridx-1,ng)+1
            tj = (ridx-1)/ng+1
            Rnew(k,s) = g
            if (k < nsite) Rnew(k+1:nsite,s) = Rfull(k+1:nsite,tj)
          end do
          rhoR(k) = nsel
          Rfull(k:nsite,1:nsel) = Rnew(k:nsite,1:nsel)
          deallocate (Amat)
        end do
      end if
!>--- close the per-sweep box
      write (stdout,'(1x,"╚",72("═"),"╝")')
    end do

    close (rawunit)

! ── energy-cache statistics ────────────────────────────────────────────────
    write (stdout,'(/,1x,a)') '> Energy cache'
    if (usecache) then
      write (stdout,'(3x,a,i0)') 'optimizations performed   : ',cache%opt
      write (stdout,'(3x,a,i0)') 'evaluations reused (cache) : ',cache%hits
    else
      write (stdout,'(3x,a)') 'disabled (every candidate re-optimized)'
    end if
    call cache%destroy()

    deallocate (Lfull,Rfull,Lnew,Rnew,rhoL,rhoR,combis,rrow,rcol,isel,energies)
    if (allocated(Aref)) deallocate (Aref)
  end subroutine ttconf_sweep

!========================================================================================!
! ── energy -> maxvol weight (low energy = large modulus) ────────────────────
  real(wp) function ttconf_weight(e,emin,kt) result(w)
    implicit none
    real(wp),intent(in) :: e,emin,kt
    if (e >= ttl_feas) then
      w = 0.0_wp                                !> infeasible / topology broken
    else
      w = exp(-(e-emin)/kt)
    end if
  end function ttconf_weight

!========================================================================================!
!>                         energy-cache implementation                                   !
!========================================================================================!
  subroutine ttconf_packkey(combi,nsite,site_ngrid,key,ok)
!************************************************************
!* Pack a grid-index combination into a single int64 key
!* via mixed-radix encoding (radix = site_ngrid per site). If
!* the product of grid sizes would overflow int64, ok is
!* returned .false. and the combination is left uncached.
!************************************************************
    implicit none
    integer,intent(in) :: nsite
    integer,intent(in) :: combi(nsite),site_ngrid(nsite)
    integer(int64),intent(out) :: key
    logical,intent(out) :: ok
    integer :: i
    integer(int64) :: base,dig
    integer(int64),parameter :: maxk = huge(0_int64)

    key = 0_int64
    ok = .true.
    do i = 1,nsite
      base = int(site_ngrid(i),int64)
      dig = int(combi(i)-1,int64)
      if (base <= 0_int64) then
        ok = .false.; return
      end if
! ── overflow guard for key = key*base + dig ──────────────────────────────────
      if (key > (maxk-dig)/base) then
        ok = .false.; return
      end if
      key = key*base+dig
    end do
  end subroutine ttconf_packkey

!========================================================================================!
  integer(int64) function ecache_hash(key) result(h)
!************************************************************************
!* Turn the packed grid-combination "key" into a well-scrambled 64-bit
!* hash. The goal is only that nearby keys (1,1,1,1) and (1,1,1,2) land
!* in very different table slots, so the open-addressing cache spreads
!* its entries evenly and stays fast. It is NOT cryptographic.
!*
!* Step by step (a "xorshift-multiply" mixer):
!*
!*  (0) mask = huge(int64) = 0x7FFF...FF, i.e. all 63 low bits set, sign
!*      bit clear. iand(x,mask) keeps a number non-negative, so the
!*      result is always a valid (>=0) array-index seed.
!*
!*  (1) h = iand(key,mask)
!*      Start from the key, sign bit cleared.
!*
!*  (2) h = ieor( h, ishft(h,-29) )
!*      ishft(h,-29) shifts h right by 29 bits; ieor is XOR. XOR-ing a
!*      value with a shifted copy of itself folds the HIGH bits down
!*      into the LOW bits, so information from the whole word starts to
!*      affect the low end (which is what the table index reads).
!*
!*  (3) h = iand( h*2654435761, mask )
!*      Multiply by a large odd constant (Knuth's 2^32/phi, the golden-
!*      ratio multiplier). Multiplication smears each input bit across
!*      many output bits; an odd multiplier is invertible, so no two
!*      keys collapse to the same value here. The overflow past 64 bits
!*      simply wraps (that wrap-around IS the mixing), and the mask then
!*      drops the sign bit again to stay non-negative.
!*
!*  (4) h = ieor( h, ishft(h,-17) )
!*      A second, shorter xorshift fold to finish scrambling the bits
!*      the multiply pushed upward back down into the low bits.
!*
!*  (5) h = iand(h,mask)
!*      Final clamp to a non-negative 63-bit value.
!*
!* The caller maps this to a slot with iand(h, cap-1) (cap is a power of
!* two), which keeps only the low bits -- hence all the effort above to
!* make the low bits depend on the entire key.
!************************************************************************
    implicit none
    integer(int64),intent(in) :: key
    integer(int64),parameter :: mask = huge(0_int64)   !> 0x7FFF...FF
    h = iand(key,mask)
    h = ieor(h,ishft(h,-29))
    h = iand(h*2654435761_int64,mask)
    h = ieor(h,ishft(h,-17))
    h = iand(h,mask)
  end function ecache_hash

!========================================================================================!
  subroutine ecache_init(self)
!************************************************************
!* Allocate the hash arrays and mark the cache active.
!************************************************************
    implicit none
    class(ttconf_ecache),intent(inout) :: self
    call self%destroy()
    self%cap = 1024
    allocate (self%keys(0:self%cap-1),source=0_int64)
    allocate (self%ener(0:self%cap-1),source=0.0_wp)
    allocate (self%used(0:self%cap-1),source=.false.)
    self%n = 0
    self%hits = 0
    self%opt = 0
    self%active = .true.
  end subroutine ecache_init

!========================================================================================!
  subroutine ecache_destroy(self)
    implicit none
    class(ttconf_ecache),intent(inout) :: self
    if (allocated(self%keys)) deallocate (self%keys)
    if (allocated(self%ener)) deallocate (self%ener)
    if (allocated(self%used)) deallocate (self%used)
    self%n = 0
    self%cap = 0
    self%active = .false.
  end subroutine ecache_destroy

!========================================================================================!
  logical function ecache_lookup(self,key,e) result(found)
!************************************************************
!* Linear-probe lookup. On a hit, return the stored energy
!* and bump the hit counter.
!************************************************************
    implicit none
    class(ttconf_ecache),intent(inout) :: self
    integer(int64),intent(in) :: key
    real(wp),intent(out) :: e
    integer :: idx
    found = .false.
    e = 0.0_wp
    if (.not.self%active) return
    idx = int(iand(ecache_hash(key),int(self%cap-1,int64)))
    do
      if (.not.self%used(idx)) return            !> empty slot -> not present
      if (self%keys(idx) == key) then
        e = self%ener(idx)
        found = .true.
        self%hits = self%hits+1
        return
      end if
      idx = iand(idx+1,self%cap-1)
    end do
  end function ecache_lookup

!========================================================================================!
  subroutine ecache_insert(self,key,e,isopt)
!************************************************************
!* Insert (or update) an entry, growing the table when the
!* load factor exceeds 0.5. "isopt" flags a genuine
!* optimization (vs. an infeasible/topology-broken combo),
!* counted for the final cache statistics.
!************************************************************
    implicit none
    class(ttconf_ecache),intent(inout) :: self
    integer(int64),intent(in) :: key
    real(wp),intent(in) :: e
    logical,intent(in) :: isopt
    integer :: idx
    if (.not.self%active) return
    if (2*(self%n+1) > self%cap) call ecache_grow(self)
    idx = int(iand(ecache_hash(key),int(self%cap-1,int64)))
    do
      if (.not.self%used(idx)) then
        self%used(idx) = .true.
        self%keys(idx) = key
        self%ener(idx) = e
        self%n = self%n+1
        if (isopt) self%opt = self%opt+1
        return
      else if (self%keys(idx) == key) then
        self%ener(idx) = e                        !> update in place
        return
      end if
      idx = iand(idx+1,self%cap-1)
    end do
  end subroutine ecache_insert

!========================================================================================!
! ── double the table capacity and rehash all live entries ───────────────────
  subroutine ecache_grow(self)
    implicit none
    type(ttconf_ecache),intent(inout) :: self
    integer(int64),allocatable :: oldkeys(:)
    real(wp),allocatable :: oldener(:)
    logical,allocatable :: oldused(:)
    integer :: oldcap,i,idx
    call move_alloc(self%keys,oldkeys)
    call move_alloc(self%ener,oldener)
    call move_alloc(self%used,oldused)
    oldcap = self%cap
    self%cap = self%cap*2
    allocate (self%keys(0:self%cap-1),source=0_int64)
    allocate (self%ener(0:self%cap-1),source=0.0_wp)
    allocate (self%used(0:self%cap-1),source=.false.)
    do i = 0,oldcap-1
      if (.not.oldused(i)) cycle
      idx = int(iand(ecache_hash(oldkeys(i)),int(self%cap-1,int64)))
      do
        if (.not.self%used(idx)) then
          self%used(idx) = .true.
          self%keys(idx) = oldkeys(i)
          self%ener(idx) = oldener(i)
          exit
        end if
        idx = iand(idx+1,self%cap-1)
      end do
    end do
    deallocate (oldkeys,oldener,oldused)
  end subroutine ecache_grow

!========================================================================================!
end module ttconf_light_mod
!========================================================================================!
