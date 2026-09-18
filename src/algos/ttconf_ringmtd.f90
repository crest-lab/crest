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
module ttconf_ringmtd_mod
!****************************************************************************************
!* Conformer source for TTConf-light ring sampling (algos layer).
!*
!* The level-agnostic geometry work (cut-out, fold-back, re-measure) lives in the
!* molbuilder module ttconf_ringsample_mod. The actual conformer GENERATION needs
!* env and the high-level search/optimization routines, which only exist this far
!* up the build -- so it is implemented here and REGISTERED into the molbuilder
!* module via its procedure-pointer hook.
!*
!* ttconf_ringmtd_provider runs, for one isolated ring cut-out, a short GFN-FF
!* MD/MTD + optimization + CREGEN in a private per-ring subdirectory (so all the
!* scratch files stay isolated) and returns a handful of representative ring
!* conformations. GFN-FF is used regardless of the run's target level.
!****************************************************************************************
  use crest_parameters,only:wp,stdout,sep,aatoau,autoaa
  use crest_data
  use strucrd,only:coord
  use ttconf_ringsample_mod,only:ttconf_set_ring_sampler,ttconf_clear_ring_sampler
  implicit none
  private

  !> env handle and per-run ring counter for the registered provider
  type(systemdata),pointer :: s_env => null()
  integer,save :: s_ringcount = 0

  public :: ttconf_enable_ring_sampling
  public :: ttconf_disable_ring_sampling

!========================================================================================!
contains
!========================================================================================!

  subroutine ttconf_enable_ring_sampling(env)
    !**************************************************************************
    !* Stash the env handle and register ONE ring-conformation provider as the
    !* sampler hook, selected by env%ttconf%ring_method. This is the single
    !* place that ties a concrete generator to a run; adding a new generator
    !* means writing a provider with the ttconf_ring_sampler_i signature (see
    !* ttconf_ringtemplate_provider below for a fill-in template) and adding a
    !* case here.
    !**************************************************************************
    implicit none
    type(systemdata),intent(inout),target :: env
    s_env => env
    s_ringcount = 0
    select case (trim(adjustl(env%ttconf%ring_method)))

    case ('mtd','gfnff','metadyn','')
      write (stdout,'(1x,a)') 'Ring sampling method: GFN-FF metadynamics'
      call ttconf_set_ring_sampler(ttconf_ringmtd_provider)

    case ('template','library','lib')
      write (stdout,'(1x,a)') 'Ring sampling method: template library (identity)'
      call ttconf_set_ring_sampler(ttconf_ringtemplate_provider)

    case default
      write (stdout,'(1x,a)') 'Unknown ring sampling method "'// &
      &  trim(adjustl(env%ttconf%ring_method))//'"; using GFN-FF metadynamics'
      call ttconf_set_ring_sampler(ttconf_ringmtd_provider)

    end select
  end subroutine ttconf_enable_ring_sampling

  subroutine ttconf_disable_ring_sampling()
    implicit none
    call ttconf_clear_ring_sampler()
    s_env => null()
    s_ringcount = 0
  end subroutine ttconf_disable_ring_sampling

!========================================================================================!

  subroutine ttconf_ringmtd_provider(frag,confs,nconf)
    !**************************************************************************
    !* Registered ring source: sample one ring cut-out and return its
    !* representative conformations as confs(3, frag%nat, nconf).
    !*
    !* Pipeline (all inside a private per-ring subdirectory so the scratch
    !* files stay isolated, GFN-FF regardless of the run's target level):
    !*   1. build a fresh GFN-FF calculator,
    !*   2. run a parallel block of nmtd RMSD-metadynamics on the cut-out, with
    !*      the bias height coupled to system size (kpush = k*Nat) and the MTD
    !*      length coupled to the cut-out size -- as the main iMTD-GC runtype,
    !*   3. optimize every snapshot (GFN-FF),
    !*   4. drop snapshots whose connectivity broke,
    !*   5. energy-window + energy-gap dedup, capped at a size-scaled maxtempl,
    !* returning the survivors (input frame). On any failure nconf=0 is
    !* returned and the caller falls back to the identity template.
    !**************************************************************************
    use crest_calculator,only:calcdata
    use dynamics_module,only:mddata,mtdpot,mdautoset,type_mtd,cv_rmsd
    use strucrd,only:coord,rdensemble
    use parallel_interface,only:crest_oloop
    use iomod,only:makedir
    implicit none
    type(coord),intent(in) :: frag
    real(wp),allocatable,intent(out) :: confs(:,:,:)
    integer,intent(out) :: nconf

    !> sampling controls
    integer,parameter  :: nmtd = 4          !> parallel MTDs per ring
    real(wp),parameter :: mtd_temp = 400.0_wp   !> MTD temperature (K)
    real(wp),parameter :: mtd_dump_fs = 500.0_wp   !> snapshot dump period (fs)
    real(wp),parameter :: ewin_kcal = 100.0_wp   !> conformer energy window
    real(wp),parameter :: egap_kcal = 0.25_wp    !> dedup energy gap
    real(wp),parameter :: autokcal = 627.50947_wp
    !> the nmtd bias settings: kpush = k*Nat couples the bias height to system
    !> size; the alpha widths are spread -- mirrors the iMTD-GC bias block
    real(wp),parameter :: kbase(nmtd) = [0.00125_wp,0.00250_wp,0.00500_wp,0.00500_wp]
    real(wp),parameter :: abias(nmtd) = [0.80_wp,1.00_wp,0.60_wp,1.30_wp]

    type(calcdata) :: gff,calcbak
    type(mddata) :: mddat,mddatbak
    type(mddata),allocatable :: mddats(:)
    type(mtdpot) :: pot
    type(coord) :: cmol
    type(coord),allocatable :: structures(:)
    character(len=64) :: dname
    character(len=128) :: ftmp
    character(len=512) :: thispath
    character(len=*),parameter :: trjf = 'crest_dynamics.trj.xyz'
    real(wp),allocatable :: xyz(:,:,:),eread(:)
    real(wp),allocatable :: esel(:)
    integer,allocatable :: at(:),order(:),keep(:)
    integer :: nat,nall,io,i,j,idx,nsel,nkeep,maxtempl
    real(wp) :: ebase,de,length_ps
    logical :: dup

    s_ringcount = s_ringcount+1
    nconf = 0
    nat = frag%nat
    write (stdout,'(3x,a,i0,a,i0,a,i0,a)') 'ring ',s_ringcount,': ',nmtd, &
    &  ' GFN-FF metadynamics runs on a ',nat,'-atom cut-out ...'
    flush (stdout)

! ── enter a private per-ring subdirectory ───────────────────────────────────
    write (dname,'(a,i0)') 'ttring_',s_ringcount
    call getcwd(thispath)
    io = makedir(trim(dname))
    call chdir(trim(dname))

! ── GFN-FF calculator (independent of the run's target level) ───────────────
    call gff%create('gfnff',chrg=s_env%chrg,uhf=s_env%uhf)
    gff%optlev = int(s_env%optlev)

! ── base MTD settings; length coupled to the cut-out size ───────────────────
    length_ps = min(40.0_wp,max(5.0_wp,0.5_wp*real(nat,wp)))
    maxtempl = max(2,min(10,nat/3))                 !> templates scale with size
    mddatbak = s_env%mddat
    call env_to_mddat(s_env)
    mddat = s_env%mddat
    mddat%requested = .true.
    mddat%shake = .false.            !> no WBO/SHAKE setup on the fragment
    mddat%length_ps = length_ps
    mddat%tsoll = mtd_temp
    mddat%dumpstep = mtd_dump_fs
    mddat%sdump = 0                  !> recomputed by mdautoset
    mddat%simtype = type_mtd
    call mdautoset(mddat,io)

! ── a parallel block of nmtd different RMSD-metadynamics (bias ~ k*Nat) ──────
    io = makedir('MDFILES')
    allocate (mddats(nmtd),source=mddat)
    do i = 1,nmtd
      mddats(i)%md_index = i
      write (ftmp,'(a,a,a,i0,a)') 'MDFILES',sep,'ring_',i,'.trj'
      mddats(i)%trajectoryfile = trim(ftmp)
      write (ftmp,'(a,a,a,i0,a)') 'MDFILES',sep,'ring_',i,'.mdrestart'
      mddats(i)%restartfile = trim(ftmp)
      mddats(i)%npot = 1
      if (allocated(mddats(i)%mtd)) deallocate (mddats(i)%mtd)
      if (allocated(mddats(i)%cvtype)) deallocate (mddats(i)%cvtype)
      pot%kpush = kbase(i)*real(nat,wp)
      pot%alpha = abias(i)
      pot%cvdump_fs = mtd_dump_fs
      pot%mtdtype = cv_rmsd
      allocate (mddats(i)%mtd(1),source=pot)
      allocate (mddats(i)%cvtype(1),source=cv_rmsd)
    end do

! ── run them in parallel; env briefly scoped to the GFN-FF calculator ───────
!  (crest_search_multimd reads only env%calc and env%mddat%requested, both of
!   which env_to_mddat/the copy below set, and merges all trajectories into
!   crest_dynamics.trj.xyz)
    cmol = frag
    call calcbak%copy(s_env%calc)
    call s_env%calc%copy(gff)
    call cmol%write('ring_cutout.xyz')
    call crest_search_multimd(s_env,cmol,mddats,nmtd)
    call s_env%calc%copy(calcbak)
    s_env%mddat = mddatbak
    deallocate (mddats)

! ── read snapshots as coord objects (already in Bohr) and optimize (GFN-FF) ──
    call rdensemble(trjf,nall,structures)
    if (nall < 1) goto 99
    allocate (at(nat),source=frag%at(1:nat))
    allocate (xyz(3,nat,nall),source=0.0_wp)
    allocate (eread(nall),source=0.0_wp)
    do i = 1,nall
      xyz(1:3,1:nat,i) = structures(i)%xyz(1:3,1:nat)   !> Bohr
      eread(i) = structures(i)%energy
    end do
    write(stdout,'(/,"> Optimizing ",i0," ring snapshots")') nall
    call crest_oloop(s_env,nat,nall,at,xyz,eread,.false.,gff)

! ── topology filter against the cut-out reference connectivity ──────────────
    allocate (keep(nall),source=0)
    nkeep = 0
    do i = 1,nall
      if (same_topology(nat,frag%at,frag%xyz,xyz(:,:,i))) then
        nkeep = nkeep+1
        keep(nkeep) = i
      end if
    end do
    if (nkeep < 1) goto 99

! ── energy-window + energy-gap dedup, capped at maxtempl ────────────────────
    allocate (order(nkeep))
    call esort_idx(nkeep,eread(keep(1:nkeep)),order)   !> ascending energy
    allocate (confs(3,nat,maxtempl),source=0.0_wp)
    allocate (esel(maxtempl),source=0.0_wp)
    nsel = 0
    ebase = eread(keep(order(1)))
    do i = 1,nkeep
      idx = keep(order(i))
      de = (eread(idx)-ebase)*autokcal
      if (de > ewin_kcal) exit                          !> sorted: rest are higher
      dup = .false.
      do j = 1,nsel
        if (abs(eread(idx)-esel(j))*autokcal < egap_kcal) then
          dup = .true.
          exit
        end if
      end do
      if (dup) cycle
      nsel = nsel+1
      esel(nsel) = eread(idx)
      confs(1:3,1:nat,nsel) = xyz(1:3,1:nat,idx)
      if (nsel >= maxtempl) exit
    end do
    nconf = nsel
    write (stdout,'(3x,a,i0,a,i0,a)') '      -> ',nall, &
    &  ' snapshots reduced to ',nconf,' ring template(s)'

99  continue
! ── always return to the original directory; tidy scratch ───────────────────
    call chdir(trim(thispath))
    !block
    !  use strucrd
    !  call wrensemble('test.xyz',nat,nconf,at,confs*autoaa)
    !end block
    call rmrf(trim(dname))
    if (nconf < 1.and.allocated(confs)) deallocate (confs)
  end subroutine ttconf_ringmtd_provider

!========================================================================================!

  subroutine ttconf_ringtemplate_provider(frag,confs,nconf)
    !**************************************************************************
    !* Reference (identity) ring-conformation provider, and a copy-paste
    !* skeleton for writing an alternative one.
    !*
    !* Any provider plugged into the ttconf_ringsample hook must obey ONE
    !* contract:
    !*   in : frag   - the isolated ring cut-out (coord, coordinates in BOHR)
    !*   out: confs  - (3, frag%nat, nconf) candidate geometries, SAME atom
    !*                 count and order as frag, also in BOHR
    !*   out: nconf  - number of candidates; nconf = 0 makes the caller fall
    !*                 back to the identity template (safe no-op)
    !* The molbuilder side then re-measures internals and transplants them, so
    !* a provider only has to PRODUCE GEOMETRIES -- it never touches z-matrices.
    !* Run context (settings, charge, threads, ...) is available through the
    !* module-scoped s_env pointer, exactly as ttconf_ringmtd_provider uses it.
    !*
    !* Extension ideas: a canonical pucker library keyed on the ring size, a
    !* user-supplied ensemble file, or a different force field / MLIP sampler.
    !* Register it from ttconf_enable_ring_sampling via a new ring_method case.
    !*
    !* The body returns the identity (one template = the input cut-out).
    !**************************************************************************
    implicit none
    type(coord),intent(in) :: frag
    real(wp),allocatable,intent(out) :: confs(:,:,:)
    integer,intent(out) :: nconf

    allocate (confs(3,frag%nat,1),source=0.0_wp)
    confs(1:3,1:frag%nat,1) = frag%xyz(1:3,1:frag%nat)
    nconf = 1
  end subroutine ttconf_ringtemplate_provider

!========================================================================================!

  logical function same_topology(nat,at,xyzref,xyzcand) result(same)
    !************************************************************
    !* Compare a candidate's covalent connectivity to a
    !* reference geometry's (CN-based, no QM).
    !************************************************************
    use strucrd,only:coord
    use adjacency,only:wbo2adjacency
    implicit none
    integer,intent(in) :: nat,at(nat)
    real(wp),intent(in) :: xyzref(3,nat),xyzcand(3,nat)
    type(coord) :: m
    real(wp),allocatable :: cn(:),bmat(:,:)
    integer,allocatable :: aref(:,:),acand(:,:)
    real(wp),parameter :: athr = 0.02_wp

    m%nat = nat
    m%at = at
    m%xyz = xyzref
    call m%cn_to_bond(cn,bmat,'cov')
    call wbo2adjacency(nat,bmat,aref,athr)
    m%xyz = xyzcand
    call m%cn_to_bond(cn,bmat,'cov')
    call wbo2adjacency(nat,bmat,acand,athr)
    same = all(aref == acand)
  end function same_topology

!========================================================================================!

  subroutine esort_idx(n,e,idx)
    !************************************************************
    !* Return idx so that e(idx(1:n)) is ascending (small n).
    !************************************************************
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: e(n)
    integer,intent(out) :: idx(n)
    integer :: i,j,k
    do i = 1,n
      idx(i) = i
    end do
    do i = 2,n
      k = idx(i)
      j = i-1
      do while (j >= 1)
        if (e(idx(j)) <= e(k)) exit
        idx(j+1) = idx(j)
        j = j-1
      end do
      idx(j+1) = k
    end do
  end subroutine esort_idx

!========================================================================================!
!========================================================================================!
end module ttconf_ringmtd_mod
!========================================================================================!
!========================================================================================!
