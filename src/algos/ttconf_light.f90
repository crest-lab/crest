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
subroutine crest_ttconf(env,tim)
!************************************************************
!* Standalone runtype for the "TTConf-light" conformer
!* search (tensor-train conformer sampling).
!*
!* This is the LIGHT reimplementation of TTConf
!* (Zurek et al., J. Chem. Theory Comput. 2025, 21, 1459):
!* it realizes the TT-cross *sweep* as a sampling heuristic,
!* NOT a full tensor-train algebra library.
!*
!* Input:
!*    env  - CREST's systemdata
!*    tim  - CREST's timer object
!*
!************************************************************
  use crest_parameters
  use crest_data
  use strucrd
  implicit none
  !> INPUT/OUTPUT
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  !> LOCAL
  type(coord) :: mol

!========================================================================================!
  call this_header()
  call tim%start(14,'TTConf-light conformer search')

!========================================================================================!
!>--- get reference structure
  call env%ref%to(mol)
  write (stdout,*)
  call smallhead('Input structure:')
  call mol%append(stdout)
  write (stdout,*)

!========================================================================================!
!>--- pass the structure to the setup/sweep machinery
  call ttconf_light_core(env,mol)

!========================================================================================!
  call tim%stop(14)
  return
!========================================================================================!
contains
!========================================================================================!
  subroutine this_header
    implicit none
    write (stdout,'(/)')
    write (stdout,'(3x,"╔",57("═"),"╗")')
    !write (stdout,'(3x,"║",57x,"║")')
    write (stdout,'(3x,"║  ",a,"  ║")') "d888888P d888888P  a88888b.                   .8888b "
    write (stdout,'(3x,"║  ",a,"  ║")') "   88       88    d8'   `88                   88   "" "
    write (stdout,'(3x,"║  ",a,"  ║")') "   88       88    88        .d8888b. 88d888b. 88aaa  "
    write (stdout,'(3x,"║  ",a,"  ║")') "   88       88    88        88'  `88 88'  `88 88     "
    write (stdout,'(3x,"║  ",a,"  ║")') "   88       88    Y8.   .88 88.  .88 88    88 88     "
    write (stdout,'(3x,"║  ",a,"  ║")') "   dP       dP     Y88888P' `88888P' dP    dP dP     "
    write (stdout,'(3x,"║",57x,"║")')
    write (stdout,'(3x,"║",5x,a,4x,"║")') 'Tensor-train conformer search (reimplementation)'
    !write (stdout,'(3x,"║",57x,"║")')
    write (stdout,'(3x,"╚",57("═"),"╝")')
    write (stdout,'(/,4x,a)') 'A lightweight reimplementation of the TTConf algorithm by'
    write (stdout,'(4x,a)') 'Zurek et al., J. Chem. Theory Comput. 2025, 21, 1459-1475.'
    write (stdout,'(4x,a)') 'https://doi.org/10.1021/acs.jctc.4c01275'
  end subroutine this_header
end subroutine crest_ttconf

!========================================================================================!
!========================================================================================!
subroutine ttconf_light_core(env,mol)
!************************************************************
!* Core driver of the TTConf-light algorithm.
!*
!* Setup: build BAT/internal coordinates via the modern
!* coord_classify engine, obtain real-valued bond orders
!* from a one-off GFN0 singlepoint (independent of env%calc),
!* and identify the relevant rotatable dihedral angles that
!* become the tensor-train variables.
!*
!* Sampling: lay the TT variables on a torsional grid and
!* generate candidate structures by one of two routes:
!*   - the TT-cross sweep (default), a sub-exponential heuristic
!*     that visits only a low-rank cross of the grid, or
!*   - the brute-force oracle, which enumerates ALL grid
!*     combinations (exhaustive reference; exponential cost).
!* Surviving structures are topology-screened, optimized, and
!* CREGEN-filtered into the final conformer ensemble.
!*
!* Input:
!*    env  - CREST's systemdata
!*    mol  - reference geometry (in Bohr by convention)
!*
!************************************************************
  use crest_parameters
  use crest_data
  use strucrd
  use crest_calculator,only:calcdata,engrad
  use molbuilder_classify_type,only:coord_classify,setup_classify
  use ttconf_light_mod,only:ttconf_identify_dihedrals,ttconf_print_dihedrals, &
  &                         ttconf_bruteforce_generate,ttconf_sweep, &
  &                         ttconf_print_settings,ttconf_force_userbonds
  use ttconf_rings_mod,only:ttconf_ringset,ttconf_build_ringsites, &
  &                         ttconf_set_ringset,ttconf_clear_ringset, &
  &                         ttconf_print_ringsites
  use ttconf_ringmtd_mod,only:ttconf_enable_ring_sampling,ttconf_disable_ring_sampling
  use ttconf_ringsample_mod,only:ttconf_ring_flexible
  use iomod
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env
  type(coord),intent(in) :: mol
  !> LOCAL
  type(coord_classify) :: molc
  type(calcdata) :: wbocalc
  real(wp),allocatable :: wbo(:,:)
  real(wp),allocatable :: grad(:,:)
  real(wp) :: energy
  integer  :: io
  !> tensor-train variable bookkeeping
  integer,allocatable :: ttbond(:,:)   !> (2,ndieder): the two atoms of each candidate bond
  real(wp),allocatable :: ttwbo(:)     !> WBO of each candidate bond
  integer,allocatable :: ttverdict(:)  !> 0=TT-variable, 1=stiff, 2=rotamer
  integer :: nttvar
  type(coord) :: molref
  !> TT-variable grid + candidate-structure generation
  integer,allocatable :: bond2tt(:)    !> candidate bond index -> TT-variable index (or 0)
  integer,allocatable :: site_bond(:,:)  !> (2,nsite): the two atoms of each TT-variable bond
  integer,allocatable :: zrow2site(:)    !> zmat row -> TT-variable index (or 0)
  integer,allocatable :: site_ngrid(:)    !> grid points per TT variable
  real(wp),allocatable :: site_step(:)     !> step size (rad) per TT variable
  integer,allocatable :: bondngrid(:)     !> per-candidate-bond grid override (0=default)
  type(ttconf_ringset) :: ringset      !> ring-site template store (optional)
  integer :: nrsite                    !> number of ring sites appended
  integer :: nsite,i,k,j,ngen,ncombi
  logical :: multilevel(6)
  character(len=*),parameter :: rawfile = 'crest_ttconf_raw.xyz'
  !> TT-cross sweep settings (read from env%ttconf, see ttconf_settings_mod)
  integer  :: ngrid,ttrank,ttsweeps,ttninit
  logical  :: use_sweep,use_cache
  integer  :: ttseed                   !> PRNG seed (<0 -> non-deterministic)
  real(wp) :: ttkt                     !> maxvol temperature in Hartree

!========================================================================================!
!>--- (0) unpack settings and print the settings block
  ngrid = env%ttconf%ngrid
  ttrank = env%ttconf%rank
  ttsweeps = env%ttconf%sweeps
  ttninit = env%ttconf%ninit
  use_sweep = env%ttconf%use_sweep
  use_cache = env%ttconf%use_cache
  ttseed = env%ttconf%seed
  ttkt = env%ttconf%kt*kcaltoau          !> kcal/mol -> Hartree
  env%ewin = env%ttconf%ewin                 !> hand the window to CREGEN
  call ttconf_print_settings(env)

!========================================================================================!
!>--- (1) real-valued bond orders from a dedicated GFN0 singlepoint.
!>    This is deliberately decoupled from env%calc, so the production
!>    level (e.g. GFN-FF) stays free; GFN-FF would not yield clean WBOs.
  call smallhead('Bond orders (GFN0-xTB)')
  call molref%copy(mol)                   !> mutable deep copy (engrad takes a target)
  call wbocalc%create('gfn0',chrg=env%chrg,uhf=env%uhf)
  wbocalc%calcs(1)%rdwbo = .true.
  allocate (grad(3,molref%nat),source=0.0_wp)
  call engrad(molref,wbocalc,energy,grad,io)
  if (io /= 0.or..not.allocated(wbocalc%calcs(1)%wbo)) then
    write (stdout,'(1x,a)') '**ERROR** could not obtain WBOs from GFN0 singlepoint'
    return
  end if
  call move_alloc(wbocalc%calcs(1)%wbo,wbo)
  write (stdout,'(1x,a,f16.8,a)') 'GFN0 singlepoint energy:',energy,' Eh'
  write (stdout,'(1x,a)') 'Real-valued Wiberg bond orders obtained.'

!========================================================================================!
!>--- (2) build the molecular graph + per-atom info from the WBO-seeded topology
  call setup_classify(mol,molc,wbo)
  call molc%print_rings(stdout)

! ── covalency check: TTConf is a single (covalent) molecule feature ──────────
  if (molc%nfrag > 1) then
    write (stdout,'(/,1x,a,i0,a)') '**WARNING** the input contains ',molc%nfrag, &
    &  ' covalent fragments. TTConf is designed for a single covalent molecule;'
    write (stdout,'(1x,a)') '            intermolecular degrees of freedom are '// &
    &  'NOT sampled and results may be meaningless.'
  end if

!========================================================================================!
!>--- (3) build the Z-matrix and its dihedral mapping.
!>    NOTE: plain (non-natural) ordering is used on purpose: get_zmat(.true.)
!>    creates extra improper/pyramidal-H ztod groups (via prune+hpyrad) that
!>    are not proper rotatable bonds. Natural ordering would give a better TT
!>    rank but must keep this dihedral mapping consistent.
  call molc%get_zmat(.false.)
  call molc%print_zmat(stdout)

!========================================================================================!
!>--- (4) identify the relevant (rotatable) dihedral angles = TT variables
  call ttconf_identify_dihedrals(molc,wbo,ttbond,ttwbo,ttverdict,nttvar, &
  &                              exclrings=env%ttconf%excl_rings)
!>    user-specified bonds (if any) REPLACE the automatic selection; the
!>    bondngrid override carries any per-bond grid-point counts to step (5)
  allocate (bondngrid(molc%ndieder),source=0)
  if (allocated(env%ttconf%userbonds)) then
    call ttconf_force_userbonds(molc,env%ttconf%userbonds,ttbond,ttverdict, &
    &                           nttvar,bondngrid)
  end if
  call ttconf_print_dihedrals(molc,ttbond,ttwbo,ttverdict,nttvar)
!>    (ring sampling can still contribute sites even with no rotatable torsions)
  if (nttvar < 1.and..not. (env%ttconf%ring_sample.and.molc%nrings > 0)) then
    write (stdout,'(1x,a)') 'No rotatable dihedral angles found. Nothing to do.'
    if (allocated(wbo)) deallocate (wbo)
    if (allocated(grad)) deallocate (grad)
    return
  end if

!========================================================================================!
!>--- (5) build the TT-variable grid: renumber verdict-0 bonds 1..nsite and map
!>    every Z-matrix row to its TT variable (0 if it belongs to none).
  nsite = nttvar
  allocate (bond2tt(molc%ndieder),source=0)
  j = 0
  do k = 1,molc%ndieder
    if (ttverdict(k) == 0) then
      j = j+1
      bond2tt(k) = j
    end if
  end do
  allocate (zrow2site(molc%nat),source=0)
  do i = 1,molc%nat
    k = molc%ztod(i)
    if (k >= 1) zrow2site(i) = bond2tt(k)
  end do
  allocate (site_bond(2,nsite),source=0)
  do k = 1,molc%ndieder
    if (ttverdict(k) == 0) site_bond(:,bond2tt(k)) = ttbond(:,k)
  end do
  allocate (site_ngrid(nsite),source=ngrid)
  allocate (site_step(nsite),source=(2.0_wp*pi/real(ngrid,wp)))
! ── apply per-bond grid overrides from user-specified bonds (if any) ────────
  do k = 1,molc%ndieder
    if (ttverdict(k) == 0.and.bondngrid(k) > 0) then
      site_ngrid(bond2tt(k)) = bondngrid(k)
      site_step(bond2tt(k)) = 2.0_wp*pi/real(bondngrid(k),wp)
    end if
  end do

!========================================================================================!
!>--- (5b) register ring sites as extra TT variables (optional, rings/puckers).
!>    The TT machinery (sweep/maxvol/cache/screen/opt) is agnostic to what a
!>    "site" is: it only needs nsite, the per-site grid size site_ngrid(site), and
!>    the combination->geometry builder. A ring is simply a site whose grid is a
!>    discrete set of pre-computed ring templates (puckers) instead of an
!>    evenly-spaced torsion angle. ttconf_build_ringsites appends one site per
!>    flexible ring (growing site_ngrid/site_step/site_bond and nsite); the
!>    template conformations come from ttconf_rings_mod::ttconf_ring_templates,
!>    and the matching overlay (ttconf_apply_active) lives in ttconf_gen_check.
  call ttconf_clear_ringset()
  if (env%ttconf%ring_sample) then
    call ttconf_enable_ring_sampling(env)   !> register the GFN-FF cut-out sampler
    call ttconf_build_ringsites(molc,nsite,site_ngrid,site_step,site_bond,ringset,nrsite)
    call ttconf_disable_ring_sampling()
    if (nrsite > 0) then
      call ttconf_print_ringsites(molc,ringset)
      call ttconf_set_ringset(ringset)   !> activate the geometry overlay
    end if
  else if (molc%nrings > 0) then
! ── flexible ring present but ring sampling off: hint at the CLI flag ────────
    do k = 1,molc%nrings
      if (ttconf_ring_flexible(molc,molc%ringlist(k))) then
        write (stdout,'(/,1x,a)') 'Note: a flexible ring was detected but ring '// &
        &  'sampling is off. Pass -ttrings to sample ring conformations.'
        exit
      end if
    end do
  end if

!========================================================================================!
!>--- (6) generate the candidate structures (topology-screened) into an ensemble
!>    file, via either the TT-cross sweep or the brute-force oracle.
  write (stdout,*)
  ngen = 0
  if (nsite < 1) then
    write (stdout,'(1x,a)') 'No TT variables (torsions or rings). Nothing to do.'
  else if (use_sweep) then
    call smallhead('Tensor-train cross sweep')
    write (stdout,'(1x,a,i0,a,i0,a)') 'TT rank r = ',ttrank,', sweeps s = ', &
    &  ttsweeps,'.'
    call ttconf_sweep(env,mol,molc,zrow2site,site_ngrid,site_step,nsite,ngrid, &
    &                 site_bond,ttrank,ttsweeps,ttninit,ttkt,rawfile,ngen,use_cache,ttseed)
    write (stdout,'(1x,a,i0,a)') 'TT-cross sweep '// &
    &  trim(merge('evaluated','optimized',env%ttconf%sp_only))//' ',ngen, &
    &  ' structures (with intact topology).'
  else
    call smallhead('Brute-force conformer generation (oracle)')
    call ttconf_bruteforce_generate(mol,molc,zrow2site,site_ngrid,site_step,nsite, &
    &                               rawfile,ngen,ncombi)
    write (stdout,'(1x,a,i0,a,i0,a)') 'Generated ',ngen,' of ',ncombi, &
    &  ' grid structures with intact topology.'
  end if

!========================================================================================!
!>--- (7) optimize + CREGEN-filter the generated ensemble
  if (ngen > 0) then
    write (stdout,*)
    if (env%ttconf%sp_only) then
!>--- singlepoint mode: evaluate the raw ensemble at env%calc WITHOUT any
!>    geometry optimization, then CREGEN-sort/deduplicate by those energies.
      call smallhead('Singlepoints and sorting')
      block
        use parallel_interface,only:crest_sploop
        use utilities,only:checkname_xyz
        use strucrd,only:coord
        integer :: spnat,spnall,i
        integer,allocatable :: spat(:)
        real(wp),allocatable :: spxyz(:,:,:),speread(:)
        type(coord),allocatable :: spmols(:)
        character(len=128) :: inpnam,outnam
        call rdensembleparam(rawfile,spnat,spnall)
        if (spnall > 0) then
          allocate (spat(spnat),spxyz(3,spnat,spnall),speread(spnall))
          call rdensemble(rawfile,spnat,spnall,spat,spxyz,speread)
          spxyz = spxyz/bohr                   !> crest_sploop expects Bohr
!>--- marshal into a coord list (canonical crest_sploop API)
          allocate (spmols(spnall))
          do i = 1,spnall
            spmols(i)%nat = spnat
            spmols(i)%at = spat
            spmols(i)%xyz = spxyz(:,:,i)
          end do
          call crest_sploop(env,spnall,spmols,speread)
          deallocate (spmols)
          spxyz = spxyz/angstrom               !> ensemble file must be Angstrom
!>--- seed the crest_rotamers_* sequence, then sort it like the opt path does
!>    (sort_and_check sets env%nat and runs CREGEN in conformer-search mode)
          call checkname_xyz(crefile,inpnam,outnam)
          call wrensemble(trim(inpnam),spnat,spnall,spat,spxyz,speread)
          deallocate (spat,spxyz,speread)
          call sort_and_check(env,trim(inpnam))
        end if
      end block
    else
      call smallhead('Optimization and sorting')
      call optlev_to_multilev(env%optlev,multilevel)
      call crest_multilevel_oloop(env,rawfile,multilevel,0)
    end if
!>--- standard CREST cleanup: latest crest_rotamers_*.xyz -> crest_rotamers.xyz,
!>    discard the numbered intermediates (and other scratch files)
    call V2terminating()
    call remove(rawfile)
    write (stdout,'(/,1x,a,1x,a)') 'Final conformer ensemble on file', &
    &  '<'//conformerfile//'>'
  else
    write (stdout,'(1x,a)') '**WARNING** no structures with intact topology generated.'
  end if

  write (stdout,*)
   write(stdout,'(1x,a)') 'Refer to this as the "CREST reimplementation of TTConf"!'
  write (stdout,'(1x,a)') 'Please make sure to cite the original TTConf method appropirately:'
  call drawbox(stdout,'',procedual=0,width=62,charset=4,ltab=1) 
  call drawbox(stdout,'C.Zurek, R.A.Malleav, A.C.Paul, N.van Staalduinen, et al.', width=62, procedual=1, charset=4, padl=1,ltab=1)
  call drawbox(stdout,'J. Chem. Theory Comput. 2025, 21, 3, 1459-1475.', width=62, procedual=1, charset=4, padl=1,ltab=1) 
  call drawbox(stdout,'https://doi.org/10.1021/acs.jctc.4c01275', width=62, procedual=1, charset=4, padl=1,ltab=1)  
  call drawbox(stdout,'',procedual=2,width=62,charset=4,ltab=1)


!========================================================================================!
  call ttconf_clear_ringset()            !> deregister the ring-template overlay
  if (allocated(wbo)) deallocate (wbo)
  if (allocated(grad)) deallocate (grad)
  if (allocated(bond2tt)) deallocate (bond2tt)
  if (allocated(site_bond)) deallocate (site_bond)
  if (allocated(zrow2site)) deallocate (zrow2site)
  if (allocated(site_ngrid)) deallocate (site_ngrid)
  if (allocated(site_step)) deallocate (site_step)
  if (allocated(bondngrid)) deallocate (bondngrid)
  return
end subroutine ttconf_light_core
!========================================================================================!
!========================================================================================!
