!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Philipp Pracht
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

!================================================================================!
! THIS MODULE CONTAINS THE MAIN AND GLOBAL DATA STORAGE OF CREST
!================================================================================!

module crest_data
  use iso_fortran_env,wp => real64,dp => int64
  use crest_calculator,only:calcdata
  use dynamics_module,only:mddata
  use strucrd,only:coord
  use crest_type_timer,only:timer
  use lwoniom_module, only: lwoniom_input
  implicit none

  public :: systemdata
  public :: timer   !> RE-EXPORT from crest_type_timer
  public :: protobj
  public :: constra
  public :: optlevflag,optlevnum,optlevmap_alt
  public :: optlev_to_multilev

  !> basename for the CRE files
  character(len=14),parameter,public :: crefile = 'crest_rotamers'
  !> basename for the conformer file
  character(len=20),parameter,public :: conformerfile = 'crest_conformers.xyz'
  !> basename for the conformer file
  character(len=16),parameter,public :: conformerfilebase = 'crest_conformers'
  !> basename for the conformer file
  character(len=19),parameter,public :: clusterfile = 'crest_clustered.xyz'
  !> basename for the cluster file
  character(len=*),parameter,public :: ensemblefile = 'crest_ensemble.xyz'
  !> basename for the log file
  character(len=*),parameter,public :: ensembleelog = 'ensemble_energies.log'
  !> basename for the mecp file
  character(len=*),parameter,public :: mecpensemble = 'crest_mecp_search.xyz'

  logical,public :: gui = .false.

!============================================================================!
! runtype variables
!============================================================================!
!&<
!>---- runtypes (variables for crest%crestver)
  integer,parameter,public :: crest_none     = 0
  integer,parameter,public :: crest_mfmdgc   = 1
  integer,parameter,public :: crest_imtd     = 2
  integer,parameter,public :: crest_imtd2    = 22
  integer,parameter,public :: crest_mdopt    = 3
  integer,parameter,public :: crest_mdopt2   = 33
  integer,parameter,public :: crest_screen   = 4
  integer,parameter,public :: crest_nano     = 7
  integer,parameter,public :: crest_compr    = 77
  integer,parameter,public :: crest_msreac   = 9
  integer,parameter,public :: crest_pka      = 14
  integer,parameter,public :: crest_solv     = 15
!>> runtypes with IDs between use non-legacy routines  <<!
  integer,parameter,public :: crest_sp         = 264
  integer,parameter,public :: crest_optimize   = 265
  integer,parameter,public :: crest_moldyn     = 266
  integer,parameter,public :: crest_s1         = 267
  integer,parameter,public :: crest_mecp       = 268
  integer,parameter,public :: crest_numhessian = 269
  integer,parameter,public :: crest_scanning   = 270
  integer,parameter,public :: crest_rigcon     = 271
  integer,parameter,public :: crest_trialopt   = 272
  integer,parameter,public :: crest_ensemblesp = 273
!>> <<!
  integer,parameter,public :: crest_test       = 456

!>---- tools (variables for env%properties)
  integer,parameter,public :: p_none         =  0
  integer,parameter,public :: p_cregen       = -1
  integer,parameter,public :: p_compare      = -2
  integer,parameter,public :: p_protonate    = -3
  integer,parameter,public :: p_deprotonate  = -4
  integer,parameter,public :: p_tautomerize  = -5
  integer,parameter,public :: p_zsort        = -6
  integer,parameter,public :: p_tautomerize2 = -555
  integer,parameter,public :: p_isomerize    = -92
  integer,parameter,public :: p_reactorset   = -312
  integer,parameter,public :: p_CREentropy   = -45
  integer,parameter,public :: p_rrhoaverage  = -4450
  integer,parameter,public :: p_cluster      = -70
  integer,parameter,public :: p_propcalc     = -666
  integer,parameter,public :: p_acidbase     = -788
  integer,parameter,public :: p_ligand       = -355
  integer,parameter,public :: p_gesc1        = -9224
  integer,parameter,public :: p_gesc2        = -9225
  integer,parameter,public :: p_thermo       = -3654
  integer,parameter,public :: p_useonly      = -227
  integer,parameter,public :: p_qcg          = 37

!>--- refinement levels (typically after multilevel opt.)
  type ,private:: refine_type
    integer :: non         = 0
    integer :: singlepoint = 1
    integer :: correction  = 2
    integer :: geoopt      = 3
    integer :: ConfSolv    = 4
  end type refine_type
  type(refine_type), parameter,public :: refine = refine_type()

!&>

!========================================================================================!
!========================================================================================!

  private

!========================================================================================!
!========================================================================================!

  type :: constra
    integer :: ndim
    logical :: used = .false.
    logical :: NCI = .false.
    logical :: ggrid = .false.
    character(len=:),allocatable :: gbsagrid
    character(len=128),allocatable :: sett(:)
    character(len=128),allocatable :: buff(:)
    character(len=128),allocatable :: pots(:)
    !--- data of bondlength constraints
    logical :: cbonds_global = .false.
    logical :: cbonds_md = .false.
    integer :: n_cbonds = 0
    character(len=128),allocatable :: cbonds(:)
    !--- data of nano-reactor
    logical :: ureactor = .false.
    integer :: nrctrl = 0
    character(len=128),allocatable :: rctrl(:)
    !--- dispersion scaling (>xtb 6.4.0)
    logical :: dispscal_md = .false.
    logical :: dispscal_global = .false.
    real(wp) :: dscal = 1.0_wp
    !--- RMSD corrected GFN/DFT hypersurface:  gESC method
    character(len=:),allocatable :: rmsdpotfile
    logical :: usermsdpot = .false.
    logical :: gesc_heavy = .false.
  contains
    procedure :: allocate => allocate_constraints
    procedure :: deallocate => deallocate_constraints
  end type constra

!========================================================================================!

  type :: protobj
    integer :: nfrag
    integer :: newchrg
    integer :: iter
    real(wp) :: popthr
    real(wp) :: ewin
    integer :: swchrg        !switch element charge
    integer :: swat          !switch element element
    logical :: swelem        !switch element to add to lmo lp pair?
    logical :: allowFrag
    logical :: threshsort    !use ewin threshold
    logical :: protdeprot    !currently unused!
    logical :: deprotprot    !(tautomerize) do first deprotonation and then protonation

    logical :: strictPDT = .false.  ! strict mode (i.e. bond constraints) for (de)protonation,tautomerization
    logical :: fixPDT = .false.  ! extension to the strict mode, fix heavy atom positions
    logical :: ABcorrection = .false.

    integer,allocatable :: atmap(:)

    !--- ligandtool
    integer :: centeratom = 0
    integer :: ligand = 0
    logical :: isatom = .false.
    character(len=:),allocatable :: infile
    character(len=:),allocatable :: newligand

    !--- pka
    integer :: h_acidic = 0  !which h atom to remove in pka script
    integer :: pka_mode = 0  !what to do in the pka calc.
    character(len=:),allocatable :: pka_baseinp  !if a base file is read in instead
    character(len=:),allocatable :: pka_acidensemble  !
    character(len=:),allocatable :: pka_baseensemble  !
    logical :: rdCFER = .false.
    character(len=:),allocatable :: cferfile

    integer :: divers = 1    !number of structures read from given ensemble for extended taut. mode
    logical :: alldivers = .false.  !use all structures of given ensemble for extended taut mode
    logical,allocatable :: blacklist(:) !a blacklist of atoms to disallow deprotonation from
  end type protobj

!========================================================================================!

!>--- ENTROPY mode setting object
  type :: entropyMTD
    integer :: nMDs           ! number of static MTDs
    integer :: nBias          ! number of Bias structures
    real(wp) :: nbiasgrow = 1.2d0
    integer :: iter = 4       ! number of iterations
    integer :: iterlast = 1   ! document how many iterations were done
    integer :: nconflast = 0  ! number of conformers in last iteraton
    real(wp) :: lenfac        ! length factor of each static MTD
    real(wp) :: temperature
    real(wp) :: kpush = 0.0005d0  ! kpush per atom !SG
    real(wp),allocatable :: klist(:)
    integer :: nklist = 1
    integer :: maxfallback = 10
    real(wp) :: alpha = 1.0d0     ! exponent
    real(wp) :: confthr = -0.01d0  ! stop iterations if we get less than this fraction of new conformers
    real(wp) :: sconvthr = -0.01d0 ! stop iterations based on estimated entroy change
    integer  :: rmax = 5          ! max ring size to exclude from bias
    integer :: katoms             ! number of atoms in bias

    real(wp) :: mtdramp = 0.03_wp
    real(wp) :: sapprox
    real(wp) :: sapproxlast
    character(len=:),allocatable :: atomlist
    logical,allocatable :: atomlist2(:)
    real(wp) :: trange(3)
    integer :: nt = 0
    real(wp),allocatable ::  cpoft(:)  !Cp(T)
    real(wp),allocatable ::  soft(:)   !S(T)
    real(wp),allocatable :: hoft(:)    !H(T)-H(0)
    !-- reference structure for BHESS
    logical :: bhess = .true.
    character(len=:),allocatable :: fromfile
    integer :: nat
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
  end type entropyMTD

!========================================================================================!
!>--- thermodynamics evaluation data
  type :: thermodata
    real(wp) :: ithr = -50.0_wp  !> imaginary mode inversion (in xtb -20.0)
    real(wp) :: fscal = 1.0_wp   !> frequency scaling
    real(wp) :: sthr = 25.0_wp   !> rot/vib interpol threshold (in xtb 50.0)
    real(wp) :: trange(3)
    integer  :: ntemps
    character(len=:),allocatable :: coords !> name of structure file for standalone use
    character(len=:),allocatable :: vibfile !> name of file containing frequencies or hessian
    real(wp),allocatable :: temps(:)
    real(wp) :: ptot = 0.9d0       !> population sthreshold
    integer  :: pcap = 50000       !> limit number of structures
    logical :: avbhess = .false.   !> use bhess in the msRRHO average calc. for all structures (expensive!)
    logical :: constrhess = .false. !> apply constraints in rrhoav?
    logical :: printpop   = .false. !> print a file with populations at different T
  contains
    procedure :: get_temps => thermo_get_temps
    procedure :: read_temps => thermo_read_temps
  end type thermodata

!========================================================================================!
!>--- storage of the reference (input structure)
  type :: refdata
    integer :: nat
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    integer :: ichrg = 0
    integer :: uhf = 0
    integer :: ntopo
    integer,allocatable :: topo(:)
    real(wp),allocatable :: charges(:)
    real(wp),allocatable :: wbo(:,:)
  contains
    procedure :: rdcharges => read_charges
    procedure :: to => ref_to_mol
    procedure :: load => ref_load_mol
  end type refdata

!========================================================================================!
!========================================================================================!
!>--- GENERAL data
  type :: systemdata
    integer :: crestver          !> Runtype-variable
    integer :: runver            !> additional runtype-variable
    integer :: properties        !> additional stuff before or after the confsearch
    integer :: properties2       !> backup variable
    integer,allocatable :: pqueue(:) !> property job queue
    integer :: npq = 0           !> number of jobs in priority queue

    real(wp) :: level            !> increase/decrease factor of #modes in V1
    real(wp) :: thresholds(8)    !> CREGEN thresholds
    real(wp) :: ewin             !> Energy window (6 kcal default)
    real(wp) :: ethr             !> dE between conformers (0.05 kcal default)
    real(wp) :: ethrpurge        !> dE between conformers (purgemode, 0.2 kcal default)
    real(wp) :: couthr           !> CREGEN coulomb energy threshold
    real(wp) :: rthr
    real(wp) :: bthr
    real(wp) :: bthr2            !> alternative BTHR (relative value BTHR threshold)
    real(wp) :: bthrmax          !> max BTHR due to anisotropy
    real(wp) :: bthrshift        !> erf-shift of BTHR
    real(wp) :: athr
    real(wp) :: pthr
    real(wp) :: pthrsum
    real(wp) :: tboltz
    logical  :: cgf(6)           !> collection of CREGEN options

    real(wp) :: mdtemps(10)      !> different temperatures for the QMDFF-MDs in V1
    real(wp) :: mdtime           !> MD length (V1&2)
    real(wp) :: elowest          !> Energy of the lowest conformer
    real(wp) :: eprivious        !> Energy of the priviously lowest conformer
    real(wp) :: gcmax            !> Max. number of structures for GC
    real(wp) :: gcmaxparent = 1.2d7 !> Max. number of parent pairs for GC

    integer :: icount            !> Mode counter in V1
    integer :: mdmode            !> MD mode in V1 xtb-MD or QMDFF
    integer :: nmodes            !> number of NMs to follow (V1)
    integer :: temps             !> number of MDs in V1, and NORMMDs in V2
    integer :: snapshots         !> number of snapshots from MD (V1)
    integer :: Maxrestart        !> max number of restarts in V1 and V3
    integer :: nreset = 0        !> tracker of algo-iterations
    integer :: nrotammds         !> Number of additional normal MDs in V2
    integer :: maxcompare        !> maximal number of (lowest) conformers to compare in -compare
    integer :: tsplit            !> number of time splitting for entropy S(t=inf) extrapolation

    integer  :: nat              !> number of atoms
    integer  :: chrg             !> molecular charge
    integer  :: uhf              !>  nα - nβ electrons
    integer  :: MAXRUN           !> number of parallel xtb jobs
    integer  :: omp              !> OMP/MKL_NUM_THREADS
    integer  :: Threads          !> Total number of threads (=omp*MAXRUN)
    integer  :: rednat           !> reduced Nat, if atom list is given
    real(wp) :: optlev
    real(wp) :: forceconst       !> forceconstant (mainly for gff iMTD-GC)

    real(wp) :: dummypercent

    !>--- various names and flags
    character(len=128) :: ensemblename   !> ensemble input name for SCREEN,MDOPT and CREGEN
    character(len=128) :: ensemblename2  !> another ensemble input name
    character(len=128) :: fixfile
    character(len=512) :: constraints    !> name of the constraint file
    character(len=20)  :: solvent        !> the solvent
    character(len=:),allocatable :: solv !> the entrie gbsa flag including solvent
    character(len=20)  :: gfnver         !> GFN version
    character(len=20)  :: gfnver2        !> GFN version (multilevel)
    character(len=20)  :: lmover         !> GFN version for LMO computation in xtb_lmo subroutine
    character(len=512) :: ProgName       !> name of the xtb executable (+ path)
    character(len=512) :: ProgIFF        !> name of xtbiff for QCG-mode
    character(len=512) :: homedir        !> original directory from which calculation was started
    character(len=512) :: scratchdir     !> path to the scratch directory
    character(len=:),allocatable :: cmd
    character(len=:),allocatable :: inputcoords
    character(len=:),allocatable :: wbofile
    character(len=:),allocatable :: atlist
    character(len=:),allocatable :: chargesfilename

    !>--- METADYN data
    real(wp) :: hmass
    real(wp) :: mdtemp
    real(wp) :: nmdtemp
    real(wp) :: mdstep
    real(wp) :: mdlenfac
    real(wp) :: tmtd
    real(wp) :: flexi
    integer  :: shake
    integer  :: mddumpxyz
    integer  :: mdskip
    integer  :: mddump
    integer  :: maxopt
    real(wp) :: hlowopt
    real(wp) :: microopt
    real(wp) :: s6opt

    integer :: nmetadyn
    real(wp),allocatable :: metadfac(:)
    real(wp),allocatable :: metadexp(:)
    integer,allocatable  :: metadlist(:)
    real(wp) :: mtd_kscal = 1.0_wp !> globally scale kpush for all metadynamics

    character(len=:),allocatable :: mtdstaticfile
    integer :: nstatic

    integer,allocatable  :: includeRMSD(:)
    logical,allocatable  :: excludeTOPO(:)

    !>--- property data objects
    type(protobj) :: ptb

    !>--- saved constraints
    type(constra) :: cts

    !>--- NCI mode data
    real(wp) :: potscal = 1.0_wp
    real(wp) :: potpad  = 0.0_wp
    character(len=:),allocatable :: potatlist

    !>--- Nanoreactor data
    real(wp) :: rdens     !reactor density
    real(wp) :: tempfermi = 6000.0d0 !logfermi temperature

    !>--- Entropy static MTDs (umbrella sampling) object
    type(entropyMTD) :: eMTD
    real(wp) :: XH3 = 0
    real(wp) :: kappa = 1.5_wp   !> vM-kernel discretization

    !>--- thermo data
    type(thermodata) :: thermo

    !>--- reference structure data (the input structure)
    type(refdata) :: ref

    !>--- QCG data
    integer :: qcg_runtype = 0      !> Default is grow, 1= ensemble & opt, 2= e_solv, 3= g_solv
    integer :: nsolv = 0            !> Number of solventmolecules
    integer :: nqcgclust = 0        !> Number of cluster to be taken
    integer :: max_solv = 0         !> Maximal number of solvents added, if none is given
    integer :: ensemble_method = -1 !> Default -1 for qcgmtd, 0= crest, 1= standard MD, 2= MTD
    character(len=:), allocatable :: directed_file !name of the directed list
    character(len=64), allocatable :: directed_list(:,:) !How many solvents at which atom to add
    integer, allocatable :: directed_number(:) !Numbers of solvents added per defined atom
    character(len=20) :: ensemble_opt         !> Method for ensemble optimization in qcg mode
    character(len=20) :: freqver              !> Method for frequency computation in qcg mode
    real(wp)          :: freq_scal            !> Frequency scaling factor
    character(len=:),allocatable :: solu_file,solv_file !> solute  and solvent input file
    character(len=5) :: docking_qcg_flag = '--qcg'

    !>--- clustering data
    integer  :: maxcluster = 0  !> maximum number of clusters to be generated
    integer  :: nclust = 0      !> fixed number of clusters (unly used if !=0)
    integer  :: pccap = 100     !> maximum number of principal components used for clustering
    real(wp) :: pcthr = 0.85d0
    real(wp) :: pcmin = 0.05d0
    real(wp) :: csthr = 0.80d0
    character(len=:),allocatable :: pcmeasure
    integer :: clustlev = 0    !> clustering level

    !>--- additional structure generation settings
    logical :: doOHflip = .true.
    integer :: maxflip = 1000

    !>--- external RMSD bias to optimizations
    character(len=:),allocatable :: biasfile
    real(wp) :: rthr2 = 0.3_wp    !>  Discard all structures with a bias smaller than this
    real(wp) :: kshift = 3.0_wp   !>  Shift of the k_i (in kcal/mol)
    integer  :: kshiftnum = 4     !>  try 5 different kshift (if not specified otherwise)
    real(wp) :: gescoptlev = 2.0_wp

    !>--- DFT driver arguments [DEPRECATED]
    character(len=:),allocatable :: dftrcfile !i.e. ~/.dftrc
    logical :: hardcutDFT = .false.  !use hard cut-off criteria for aoforce consideration
    real(wp) :: harcutpthr = 0.5_wp  !take no less than 50% population of ensemble
    integer  :: hardcutnst = 5       !take no more than 5 structures

    !================================================!
    !>--- Calculation settings for newer implementations (version >= 3.0)
    type(calcdata) :: calc
    type(mddata)   :: mddat
    !>--- rigidconf data
    integer :: rigidconf_algo = 0
    integer :: rigidconf_toposource = 0
    character(len=:),allocatable :: rigidconf_userfile
    !>--- refinement queue
    integer,allocatable :: refine_queue(:)
    !>--- lwONIOM input
    character(len=:),allocatable    :: ONIOM_toml
    type(lwoniom_input),allocatable :: ONIOM_input
    !================================================!

    !>--- msreact mode settings
    logical :: msnoiso =.false. ! print only dissociated structures in msreact
    logical :: msiso =.false. ! only print non-dissociated structures in msreact
    logical :: msmolbar =.false. ! sort out duplicates by molbar
    logical :: msinchi =.false. ! sort out duplicates by inchi
    logical :: mslargeprint=.false. ! dont remove temporary files
    logical :: msattrh=.true. ! add attractive potential for H-atoms
    integer :: msnbonds = 3 ! distance of bonds up to nonds are stretched
    integer :: msnshifts = 0 ! number of random shifts applied to whole mol
    integer :: msnshifts2 = 0 ! number of random shifts applied to whole mol
    integer :: msnfrag = 0 !number of fragments that are printed in msreact mode
    character(len=80) :: msinput = '' !name of input file for msreact

    !>--- general logical data
    logical :: allrot = .true.       !> use all rotational constants for check instead of mean?
    logical :: altopt = .false.
    logical :: autothreads           !> automatically determine threads
    logical :: autozsort             !> do the ZSORT in the beginning ?
    logical :: allowrestart = .true. !> allow restart in crest algos?
    logical :: better                !> found a better conformer and restart in V1
    logical :: cff                   !> CFF used in QCG-energy calculation
    logical :: cluster = .false.     !> perform a clustering analysis
    logical :: checktopo = .true.    !> perform topolgy check in CREGEN
    logical :: checkiso = .false.    !> perform E/Z isomerization check in CREGEN
    logical :: chargesfile = .false. !> use a given charges file for gfnff
    logical :: compareens            !> try to correlate 2 given Ensemble files
    logical :: confgo                !> perform only the CREGEN routine ?
    logical :: constrain_solu        !> constrain the solute
    logical :: crest_ohess = .false. !> append numerical Hessian after optimization
    logical :: doNMR                 !> determine NMR equivalencies in CREGEN ?
    logical :: dryrun = .false.      !> dryrun to print settings
    logical :: ENSO                  !> some options for usage of CREST within ENSO
    logical :: ens_const = .false.   !> constrain solute also in Ensemble generation
    logical :: entropic = .false.    !> entropy mode
    logical :: entropymd = .false.   !> entropy mode static mtds
    logical :: esort = .false.       !> legacy option in old cregen
    logical :: ext                   !> external
    logical :: extLFER = .false.     !> read in external LFER parameters
    logical :: FINAL_GFN2_OPT = .false.
    logical :: fullcre = .false.     !> calculate exact rotamer degeneracies
    logical :: gbsa                  !> use gbsa
    logical :: gcmultiopt            !> 2 level optimization for GC in V2
    logical :: heavyrmsd = .false.   !> use only heavy atoms for RMSD in CREGEN?
    logical :: inplaceMode = .true.  !> in-place mode: optimization dirs are created "on-the-fly"
    logical :: iterativeV2           !> iterative version of V2 (= V3)
    logical :: iru                   !> re-use previously found conformers as bias in iterative approach
    logical :: keepModef             !> keep MODEF* dirs in V1 ?
    logical :: keepScratch = .false. !> keep scratch directory or delete it?
    logical :: legacy = .false.       !> switch between the original system call routines of crest and newer, e.g. tblite implementations
    logical :: metadynset            !> is the number of MTDs already set (V2) ?
    logical :: methautocorr          !> try to automatically include Methyl equivalencies in CREGEN ?
    logical :: multilevelopt =.true. !> perform the multileveloptimization
    logical :: newcregen = .false.   !> use the CREGEN rewrite
    logical :: NCI                   !> NCI special usage
    logical :: niceprint             !> make a nice progress-bar printout
    logical :: noconst = .false.     !> no constrain of solute during QCG Growth
    logical :: onlyZsort             !> do only the ZSORT routine ?
    logical :: optpurge = .false.    !> MDOPT purge application
    logical :: outputsdf = .false.   !> write output ensemble as sdf?
    logical :: pcaexclude = .false.  !> exclude user set atoms from PCA?
    logical :: pclean                !> cleanup option for property mode
    logical :: performCross          !> perform the GC in V1/V2 ?
    logical :: performMD             !> perform the MD in V1 ?
    logical :: performModef          !> perform the MF in V1 ?
    logical :: performMTD            !> perform the MTD in V2 ?
    logical :: preactormtd           !> prepare reactor mtd?
    logical :: preactorpot           !> prepare reactor logfermi?
    logical :: preopt                !> do a GFNn-xTB pre-optimization of the geometry?
    logical :: presp = .false.       !> do a Sp calculation before starting a job?
    logical :: printscoords          !> write scoord.* files in CREGEN ?
    logical :: QCG                   !> QCG special usage
    logical :: qcg_flag = .false.    !> QCG-parsing logical, only true, if qcg exclusive flags used
    logical :: qcg_restart = .false. !> QCG, only true, if results from previous run are found
    logical :: nopreopt = .false.    !> Switch off preoptimization for QCG
    logical :: quick                 !> quick-run option (mainly for testing)
    logical :: readbias = .false.    !> read MTD parameters from file
    logical :: reftopo = .true.      !> use a reference topology from the given input structure
    logical :: relax = .false.       !> was the --relax function used for protonation site search?
    logical :: restartopt            !> restart in the second step of the multilevel opt (V2) ?
    logical :: reweight = .false.    !> reweight structures on the fly after optimizations (i.e. do SPs)?
    logical :: riso = .false.        !> take only isomers in reactor mode
    logical :: rotamermds            !> do additional MDs after second  multilevel OPT step in V2 ?
    logical :: refine_presort = .false.  !> run CREGEN at the beginning of crest_refine?
    logical :: refine_esort   = .false.  !> if CREGEN is run after crest_refine, only sort energy?
    logical :: sameRandomNumber = .false. !> QCG related, choose same random number for iff
    logical :: scallen               !> scale the automatically determined MD length by some factor?
    logical :: scratch               !> use scratch directory
    logical :: setgcmax = .false.    !> adjust the maxmimum number of structures taken into account for GC?
    logical :: sdfformat = .false.   !> was the SDF format used as input file?
    logical :: slow                  !> slowmode (counterpart to quick mode)
    logical :: solv_md = .false.     !> switches on QCG-ensemblerun instead of CFF
    logical :: staticmtd = .false.   !> do a static MTD instead of normal MDs
    logical :: subRMSD               !> include only the selected substructure into the CREGEN RMSD
    logical :: superquick            !> very crude quick-run option
    logical :: threadssetmanual      !> are #CPUs set with the '-T' flag ?
    logical :: trackorigin           !> track the origin of a conformation?
    logical :: testnumgrad = .false. !> test numerical gradient in singlepoint
    logical :: use_xtbiff = .false.  !> use xtbiff for QCG?
    logical :: user_enslvl = .false. !> true if user set qcg enslvl
    logical :: user_temp = .false.   !> true if user set the MD temp
    logical :: user_mdtime = .false. !> true if mdtime set by user
    logical :: user_mdstep = .false. !> true if mdstep is set by user
    logical :: user_nclust = .false. !> true if number of cluster is set by user (only QCG)
    logical :: user_dumxyz = .false. !> true if dumpxyz is set by user
    logical :: user_wscal = .false.  !> true if wscal is set by user
    logical :: useqmdff              !> use QMDFF in V2?
    logical :: water = .false.       !> true if water is used as solvent (only QCG)
    logical :: wallsetup = .false.   !> set up a wall potential?
    logical :: wbotopo = .false.     !> set up topo with WBOs
  contains
    procedure :: allocate => allocate_metadyn
    procedure :: deallocate => deallocate_metadyn
    procedure :: addjob => add_to_pqueue
    procedure :: checkhy => pqueue_hybrid
    procedure :: rmhy => pqueue_removehybrid
    procedure :: addrefine => add_to_refinequeue
    procedure :: wrtCHRG => wrtCHRG
  end type systemdata

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine allocate_metadyn(self,n)
    implicit none
    class(systemdata) :: self
    integer,intent(in)  :: n
    self%nmetadyn = n
    if (.not.allocated(self%metadfac)) then
      allocate (self%metadfac(n),source=0.0_wp)
    end if
    if (.not.allocated(self%metadexp)) then
      allocate (self%metadexp(n),source=0.0_wp)
    end if
    if (.not.allocated(self%metadlist)) then
      allocate (self%metadlist(n),source=0)
    end if
    return
  end subroutine allocate_metadyn
!========================================================================================!
  subroutine deallocate_metadyn(self)
    implicit none
    class(systemdata) :: self
    if (allocated(self%metadfac)) deallocate (self%metadfac)
    if (allocated(self%metadexp)) deallocate (self%metadexp)
    if (allocated(self%metadlist)) deallocate (self%metadlist)
  end subroutine deallocate_metadyn
!========================================================================================!
  subroutine allocate_constraints(self,n)
    implicit none
    class(constra) :: self
    integer,intent(in)  :: n
    self%ndim = n
    allocate (self%sett(n))
    allocate (self%buff(n))
    self%sett = ''
    self%buff = ''
  end subroutine allocate_constraints
!========================================================================================!
  subroutine deallocate_constraints(self)
    implicit none
    class(constra) :: self
    if (allocated(self%sett)) deallocate (self%sett)
    if (allocated(self%buff)) deallocate (self%buff)
  end subroutine deallocate_constraints

!========================================================================================!
!========================================================================================!

  subroutine add_to_pqueue(self,pjob)
    implicit none
    class(systemdata) :: self
    integer :: pjob
    integer :: idum
    integer,allocatable :: qdum(:)
    idum = self%npq
    self%npq = self%npq+1 !add a job to the queue
    if (.not.allocated(self%pqueue)) then
      allocate (self%pqueue(1))
      self%pqueue(1) = pjob
    else
      allocate (qdum(self%npq))
      qdum(1:idum) = self%pqueue(1:idum)
      qdum(self%npq) = pjob
      call move_alloc(qdum,self%pqueue)
    end if
    return
  end subroutine add_to_pqueue

  subroutine pqueue_hybrid(self)
!> check the queue for requested hybrid reoptimization (e.g. '-gfn2@gff')
!> and sort it accordingly
    implicit none
    class(systemdata) :: self
    integer :: i,pjob,n,k
    integer,allocatable :: qdum(:)
    n = self%npq
    if (self%npq .gt. 1) then
      if (any(self%pqueue .ge. 50.and.self%pqueue .lt. 60)) then
        allocate (qdum(n),source=0)
        do i = 1,n
          pjob = self%pqueue(i)
          if (pjob .ge. 50.and.pjob .lt. 60) then
            qdum(1) = pjob
          end if
        end do
        k = 2
        do i = 1,n
          pjob = self%pqueue(i)
          if (pjob == 0) cycle
          if (pjob .lt. 50.or.pjob .ge. 60) then
            qdum(k) = pjob
            k = k+1
          end if
        end do
        call move_alloc(qdum,self%pqueue)
      end if
    end if
    return
  end subroutine pqueue_hybrid

  subroutine pqueue_removehybrid(self)
    implicit none
    class(systemdata) :: self
    integer :: i,pjob,n
    n = self%npq
    if (self%npq .gt. 1) then
      do i = 1,n
        pjob = self%pqueue(i)
        if (pjob .ge. 50.and.pjob .lt. 60) then
          self%pqueue(i) = 0
          self%npq = self%npq-1
        end if
      end do
    end if
    return
  end subroutine pqueue_removehybrid


  subroutine add_to_refinequeue(self,refinetype)
    implicit none
    class(systemdata) :: self
    integer :: refinetype
    integer :: idum
    integer,allocatable :: qdum(:)
    if( refinetype <= 0 ) return
    if (.not.allocated(self%refine_queue)) then
      allocate (self%refine_queue(1))
      self%refine_queue(1) = refinetype
    else
      idum = size(self%refine_queue,1)
      allocate (qdum(idum+1))
      qdum(1:idum) = self%refine_queue(1:idum)
      qdum(idum+1) = refinetype
      call move_alloc(qdum,self%pqueue)
    end if
    return
  end subroutine add_to_refinequeue

!========================================================================================!
!========================================================================================!
!> write a .CHRG (and .UHF) file in the specified dir, but only if it is needed
  subroutine wrtCHRG(self,dir)
    implicit none
    class(systemdata) :: self
    character(len=*) :: dir
    character(len=:),allocatable :: path
    integer :: ich,k,i
    k = len_trim(dir)
    if (self%chrg .ne. 0) then
      if (k > 0) then
        path = trim(dir)//'/'//'.CHRG'
      else
        path = '.CHRG'
      end if
      open (newunit=ich,file=path)
      write (ich,*) self%chrg
      if (self%QCG) then
        write (ich,*) self%chrg
        write (ich,*) '0'
      end if
      close (ich)
    end if
    if (self%UHF .ne. 0) then
      if (k > 0) then
        path = trim(dir)//'/'//'.UHF'
      else
        path = '.UHF'
      end if
      open (newunit=ich,file=path)
      write (ich,*) self%uhf
      if (self%QCG) then
        write (ich,*) self%uhf
        write (ich,*) '0'
      end if
      close (ich)
    end if
    if (self%chargesfile.and.allocated(self%ref%charges)) then
      if (k > 0) then
        path = trim(dir)//'/'//'charges'
      else
        path = 'charges'
      end if
      open (newunit=ich,file=path)
      do i = 1,self%ref%nat
        write (ich,'(1x,f16.8)') self%ref%charges(i)
      end do
      close (ich)
    end if
    return
  end subroutine wrtCHRG

!========================================================================================!
!> read atomic charges from a file (one line per atom)
  subroutine read_charges(self,chargesfilename,totchrg)
    implicit none
    class(refdata) :: self
    character(len=*) :: chargesfilename
    integer :: ich,io,i
    real(wp) :: dum,tot
    integer :: totchrg
    if (allocated(self%charges)) deallocate (self%charges)
    if (self%nat > 0) then
      allocate (self%charges(self%nat),source=0.0_wp)
      open (newunit=ich,file=chargesfilename)
      do i = 1,self%nat
        read (ich,*,iostat=io) dum
        if (io == 0) then
          self%charges(i) = dum
        end if
      end do
      close (io)
    end if
    tot = 0.0_wp
    do i = 1,self%nat
      tot = tot+self%charges(i)
    end do
    totchrg = nint(tot)
    return
  end subroutine read_charges

  subroutine ref_to_mol(self,mol)
    implicit none
    class(refdata) :: self
    type(coord) :: mol
    mol%nat = self%nat
    if(allocated(self%at)) mol%at = self%at
    if(allocated(self%xyz)) mol%xyz = self%xyz
    mol%chrg = self%ichrg
    mol%uhf = self%uhf
    return
  end subroutine ref_to_mol

  subroutine ref_load_mol(self,mol)
    implicit none
    class(refdata) :: self
    type(coord) :: mol
    self%nat    = mol%nat 
    self%at     = mol%at  
    self%xyz    = mol%xyz 
    self%ichrg  = mol%chrg
    self%uhf    = mol%uhf 
    return
  end subroutine ref_load_mol

!========================================================================================!
!========================================================================================!

  function optlevflag(optlev) result(flag)
    implicit none
    real(wp),intent(in) :: optlev
    character(len=:),allocatable :: flag
    flag = ''
    if (optlev <= 3.0d0) flag = 'extreme'
    if (optlev <= 2.0d0) flag = 'very tight'
    if (optlev <= 1.0d0) flag = 'tight'
    if (optlev <= 0.0d0) flag = 'normal'
    if (optlev <= -1.0d0) flag = 'loose'
    if (optlev <= -2.0d0) flag = 'very loose'
    if (optlev <= -3.0d0) flag = 'crude'
    if (optlev <= -4.0d0) flag = 'lax'
    return
  end function optlevflag

  function optlevnum(flag) result(optlev)
    implicit none
    real(wp) :: optlev
    character(len=*):: flag
    optlev = 0.0_wp
    if (index(flag,'crude') .ne. 0) optlev = -3.0d0
    if (index(flag,'loose') .ne. 0) optlev = -1.0d0
    if (index(flag,'vloose') .ne. 0) optlev = -2.0d0
    if (index(flag,'sloppy') .ne. 0) optlev = -2.0d0
    if (index(flag,'normal') .ne. 0) optlev = 0.0d0
    if (index(flag,'tight') .ne. 0) optlev = 1.0d0
    if (index(flag,'verytight') .ne. 0) optlev = 2.0d0
    if (index(flag,'vtight') .ne. 0) optlev = 2.0d0
    if (index(flag,'extreme') .ne. 0) optlev = 3.0d0 
    if (index(flag,'3') .ne. 0) optlev = 3.0d0
    if (index(flag,'2') .ne. 0) optlev = 2.0d0
    if (index(flag,'1') .ne. 0) optlev = 1.0d0
    if (index(flag,'0') .ne. 0) optlev = 0.0d0
    if (index(flag,'-3') .ne. 0) optlev = -3.0d0
    if (index(flag,'-2') .ne. 0) optlev = -2.0d0
    if (index(flag,'-1') .ne. 0) optlev = -1.0d0
    return
  end function optlevnum

  function optlevmap_alt(optin) result(optout)
!> a mapping from env%optlev to integer
!> that is used in the multilevel_oloop selection
    implicit none
    integer :: optout
    real(wp),intent(in) :: optin
    integer :: k
    integer,parameter :: omap(7) = [1,2,3,4,5,6,6]
    k = nint(optin)+4
    optout = omap(k)
  end function optlevmap_alt

  subroutine optlev_to_multilev(optlev,multilev)
    implicit none
    real(wp),intent(in) :: optlev
    logical,intent(out) :: multilev(6)
    integer :: j
    if (optlev <= 3.0d0)then !> "extreme" thresholds
      multilev(:) = .false.
      multilev(6) = .true.
      multilev(4) = .true.
      multilev(1) = .true.
    endif
    j = optlevmap_alt(optlev)
    j = max(j-1, 1)  !> j is reduced by one
    if (optlev <= 2.0d0)then  !> "normal" to "vtight"
     multilev(:) = .false.
     multilev(1) = .true.
     multilev(j) = .true.
    endif
    if (optlev <= -1.0d0)then !> "loose" to "crude"
     multilev(:) = .false.
     multilev(j) = .true.
    endif
  end subroutine optlev_to_multilev

!========================================================================================!
!========================================================================================!

  subroutine thermo_get_temps(self)
    implicit none
    class(thermodata) :: self
    integer :: i,nt
    real(wp) :: dum1
    if (allocated(self%temps)) then
      deallocate (self%temps)
    end if
    dum1 = (self%trange(2)-self%trange(1))
    nt = nint(dum1/self%trange(3))+1
    if (nt < 1) nt = 1
    self%ntemps = nt
    allocate (self%temps(nt))
    dum1 = self%trange(1)
    do i = 1,nt
      self%temps(i) = max(dum1,1.0_wp) !> never do 0 K
      dum1 = dum1+self%trange(3)
    end do
    return
  end subroutine thermo_get_temps

  subroutine thermo_read_temps(self,fname)
    implicit none
    class(thermodata) :: self
    character(len=*) :: fname
    integer :: i,nt,io,ich
    real(wp) :: dum1
    character(80) :: atmp
    write (*,*) 'reading from ',trim(fname)
    if (allocated(self%temps)) then
      deallocate (self%temps)
    end if
    open (newunit=ich,file=fname)
    nt = 0
    do
      read (ich,'(a)',iostat=io) atmp
      if (io < 0) exit  !EOF
      read (atmp,*,iostat=io) dum1
      if (io == 0) then
        nt = nt+1
      end if
    end do
    self%ntemps = nt
    write (*,*) nt,' temperatures'
    allocate (self%temps(nt))
    rewind ich
    i = 0
    do
      read (ich,'(a)',iostat=io) atmp
      if (io < 0) exit  !EOF
      read (atmp,*,iostat=io) dum1
      if (io == 0) then
        i = i+1
        self%temps(i) = dum1
      end if
    end do
    close (ich)
    write (*,*) self%ntemps
    write (*,*) self%temps
    return
  end subroutine thermo_read_temps

!========================================================================================!
!========================================================================================!
end module crest_data
