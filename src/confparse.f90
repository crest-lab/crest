!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2022 Philipp Pracht
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
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> ARGUMENT PARSER FOR CREST
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
!> The parseflags routine does at one expects from
!> its name, but also some other things as setting
!> several defaults for the program.
!> This is the program initialization.
!>
!> Input/Output:
!>  env  -  crest's systemdata object, which
!>          contains basically all information
!>          for the calculation
!>  arg  -  an array of the command line args
!>          that were read in the beginning
!>  nra  -  number of command line args in "arg"
!>-----------------------------------------------
subroutine parseflags(env,arg,nra)
  use iso_fortran_env,wp => real64
  use crest_data
  use crest_calculator
  use iomod
  use utilities
  use strucrd
  use dynamics_module
  use optimize_module
  use parse_inputfile
  use crest_restartlog
  use parse_hybrid

  implicit none
  type(systemdata),intent(inout) :: env
  integer,intent(in) :: nra
  real(wp),allocatable :: xx(:),floats(:)
  character(len=256),allocatable :: strings(:)
  character(len=*) :: arg(nra)
  character(len=1024) :: cmd
  character(len=512) :: atmp,btmp
  character(len=:),allocatable :: ctmp,dtmp
  integer :: i,j,k,l,io,ich,idum
  real(wp) :: rdum
  integer :: ctype
  logical :: ex,bondconst,presetok
  logical :: user_finopt
  character(len=:),allocatable :: argument
  logical,allocatable :: processedarg(:)
  logical,allocatable :: atlist(:)
  character(len=:),allocatable :: arg1,arg2,arg3
  character(len=:),allocatable :: hybrid_quality,hybrid_workhorse
  character(len=4) :: hybrid_mode

  allocate (xx(10),floats(3),strings(3))
  ctmp = ''
  dtmp = ''

  allocate (processedarg(nra),source=.false.)
!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> Set the defaults
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
  do i = 1,nra
    if (any((/'--GUI','--gui'/) == trim(arg(i)))) gui = .true.
    if ('-niceprint' == trim(arg(i))) env%niceprint = .true.
    if (any((/character(9)::'-version','--version'/) == trim(arg(i)))) then
      call confscript_head(.true.)
      stop
    end if
  end do

!=========================================================================================!
!>--- print the program header and command line input
  call get_command(cmd)
  if (.not.gui) then
    call confscript_head(.false.)

    write (stdout,'(/,1x,a)') 'Command line input:'
    write (stdout,'(1x,a,a,/)') '$ ',trim(cmd)
  end if
  env%cmd = trim(cmd)

!=========================================================================================!
!>--- check if help is requested or citations shall be diplayed
  do i = 1,nra
    if (processedarg(i)) cycle
    if (any((/character(6)::'-h','-H','--h','--H','--help'/) == trim(arg(i)))) then
      if (nra > i) then
        ctmp = trim(arg(i+1))
        if (ctmp(1:1) .ne. '-') then
          call confscript_morehelp(ctmp)
        end if
      end if
      call confscript_help()
    end if
    if (any((/character(10)::'-cite','--cite','--citation'/) == trim(arg(i)))) then
      call crestcite()
    end if
    if (index(arg(i),'-newversion') .ne. 0) then !> as in CREST version >= 3.0
      env%legacy = .false.
      !processedarg(i) = .true.
    end if
    if (index(arg(i),'-legacy') .ne. 0) then  !> as in CREST version <3.0
      env%legacy = .true.
      !processedarg(i) = .true.
    end if
    if (index(arg(i),'-dry') .ne. 0) then   !> "dry" run to print settings
      env%dryrun = .true.
      processedarg(i) = .true.
    end if
  end do

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!>    D E F A U L T   S E T T I N G S
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
!&<
!>--- parallelization stuff
  env%Threads = 1                !> total number of threads
  env%MAXRUN = 1                 !> number of parallel xtb jobs
  env%omp = 1                    !> # of OMP_NUM_THREADS and MKL_NUMTHREADS to be used
  env%autothreads = .true.       !> automatically determine optimal parameters omp and MAXRUN
  env%threadssetmanual = .false. !> did the user set the #threads manually?

  env%scratch = .false.          !> use scratch directory?
  call getcwd(env%homedir)       !> original directory
  env%scratchdir = ''            !> directory that shall be used for scratch

!>--- xtb settings
  env%ProgName = 'xtb'           !> the name of the xtb executable used per default
  env%ProgIFF = 'xtbiff'         !> name of the IFF that is per default xtbiff, only for  QCG
  env%optlev = 2.0d0             !> optimization level for the GFN-xTB optimizations in ALL steps
  env%gbsa = .false.             !> use GBSA (or ALPB)
  env%solv = ''                  !> if gbsa is used, the entrie flag will be written into here
  env%chrg = 0                   !> molecular chrg
  env%uhf = 0                    !> nα-nβ electrons
  env%gfnver = '--gfn2'          !> selct the GFN verison as complete flag(!)
  env%gfnver2 = ''               !> a second level, used for multilevel post-optimization
  env%ensemble_opt = '--gff'     !> qcg specific method for ensemble search and optimization

!--- cregen settings
  env%confgo = .false.           !> perform confg (cregen) subroutine only
  env%methautocorr = .false.     !> for the "-metac" flag
  env%printscoords = .false.     !> write scoord.* files?
  env%doNMR = .false.            !> option for the very last confg call
  env%elowest = 0.0d0            !> energy of the lowest conformer
  env%ENSO = .false.             !> some options for CREST usage within ENSO
  env%subRMSD = .false.          !> use only the RMSD of the selected atoms (e.g. if given bei atomlist+/-) in CREGEN

  env%newcregen = .true.         !> use the re-written version of CREGEN?
  env%checktopo = .true.       !> check topology (neighborlists) in CREGEN
  env%checkiso = .false.       !> check for E/Z C=C isomerizations

!>--- important: DEFAULT THRESHOLDS
  env%ewin = 6.0d0    !> EWIN - energy window in kcal/mol (for confg and confcross)
  env%rthr = 0.125d0  !> RTHR - RMSD thr in Angstroem
  env%ethr = 0.05d0   !> ETHR - E threshold in kcal
  env%ethrpurge = 0.20d0   !> ETHRPURGE - E threshold in kcal (purgemode)
  env%couthr = 0.1d0    !> COUTHR - CREGEN coulomb energy threshold
  env%thresholds(4) = 15.0d0   !> BTHR - rot.const. thr., A B and C in MHz
  env%bthr2 = 0.01d0   !> BTHR2 - rot.const relative deviation
  env%bthrmax = 0.025d0  !> max value of BTHR2 for anisotropic rot. const.
  env%bthrshift = 0.5d0    !> BHTR2 anisotropy error-function shift
  env%athr = 0.04d0   !> ATHR -to det. int. rotation. equal atoms for NMR, CRITICAL!
  env%pthr = 0.05d0   !> PTHR - population thr (I don't really know why we have this)
  env%pthrsum = 0.85d0   !> PTHRSUM - sum of populations threshold
  env%tboltz = 298.15d0 !> T - Temperature for Boltzmann weights
  env%esort = .false.  !> legacy option for energy sorting only in cregen

  !---set default logical options for confg
  env%cgf(1) = .false.         !> DEBUG
  env%cgf(2) = .true.          !> NEWFILE
  env%cgf(3) = .false.         !> ANAL
  env%cgf(4) = .false.         !> HEAVY - compare only heavy atoms + OH in RMSD
  env%cgf(5) = .true.          !> RMSDCHK
  env%cgf(6) = .false.         !> write confg output to file <tmp> instead of <confg.out>

!>--- general runtype settings (shared for V1 & V2 and other functionallities)
  env%autozsort = .false.      !> zsort at the beginning
  env%onlyZsort = .false.      !> perform only the zsort routine
  env%performCross = .true.    !> do the GC in V1 and V2
  env%slow = .false.
  env%setgcmax = .false.       !> adjust max. number of structures for GC?
  env%quick = .false.          !> use loose options for a quick conformation search
  env%superquick = .false.     !> very crude variant of quick-mode
  env%niceprint = .false.      !> progressbar printout for some of the steps
  env%multilevelopt = .true.   !> perform multilevel optimization
  env%trackorigin = .true.     !> for v2 track generation step by default
  env%compareens = .false.     !> compare two given ensembles
  env%maxcompare = 100         !> maximum number of (lowest) structures to compare when using "-compare"
  env%QCG = .false.          !> special QCG usage

!>--- The following settings are mainly for v.1 (MF-MD-GC)
  env%level = 1             !> full number of modes
  env%performMD = .true.    !do the MD in V1
  env%keepModef = .true.    !keep the MODEF* directories ate the end?
!>---some md defaults
  env%mdmode = 0            !1=qmdff, 0=normal md
  env%mdtime = -1.0d0       !dummy argument, the actual MD length is set depending on the number of atoms and number of meta-MDs
  env%snapshots = 100       !number of snapshots to be taken from the MD
  env%temps = -1            !dummy argument for normMDs
  env%shake = 2             !SHAKE on
  env%mdtemps = 0
  env%nrotammds = -1    !number of normMDs (-1 means that a default will be set automatically)
  env%mdtemps(1:3) = (/500,400,300/) !default temperatures for the 3 QMDFF MDs

!>--- Mixed MD settings required by iMTD-GC (V2)
  env%hmass = 2.0d0       ! if hmass=0, hmass is taken from the .xtbrc
  env%mdtemp = 300.d0     ! Due to Guiding Force the Temperature is not so important
  env%nmdtemp = 400.0d0   ! base temperature for additional normal MDs
  env%mdstep = 5.0d0      ! 4 fs
  env%shake = 2           ! shake 1 makes it more stable but requires mdstep 2.0
  env%mddumpxyz = 100     ! if not set by the user mddumpxyz is adjusted in subroutine
  env%mdskip = 1          ! 1 --> no structure is skipped
  env%mddump = 1000       ! Vbias dump in fs, e.g. every 2 ps new Vbias
  env%scallen = .false.   ! scale md length?
  env%mdlenfac = 1.0d0     ! md length scaling factor
  env%hlowopt = 0.005d0   ! ancopt
  env%microopt = 20d0     !  "
  env%s6opt = 20.0d0      !  "
  env%Maxrestart = 5       !maximum number of restarts

!>--- Settings for MTD-GC (V2)
  env%restartopt = .false.  !> jump to second iteration of the Multilevel optimization (V2 only)
  env%rotamermds = .true.   !> do some additional mds for the lowermost conformers in V2 (after first step of multilevel optimization)
  env%gcmultiopt = .true.   !> optimize in two steps after GC (loose/vtight) in V2 ? !SG
  env%performMTD = .true.   !> do the MTD in V2
  env%metadynset = .false.  !> is the metadyn prepared? (V2)
  env%useqmdff = .false.    !> use qmdff for the MDs?
  env%iru = .false.         !> re-use previously found conformers as bias in iterative approach
!>--- array to determine if RMSD are included
  env%keepModef = .false.   !> delete intermediate Directories
  env%nmetadyn = 0          !> number of METADYNs (dummy argument at this point; set later)

  env%forceconst = 0.02_wp !> force constant (mainly for GFN-FF iMTD-GC)
  bondconst = .false.      !> constrain all bonds
  user_finopt = .false.    !> -fin_opt_gfn2 given explicitly (QCG)
  ctype = 0                !> bond constraint type

  env%NCI = .false.      !> use specialized NCI mode?
  env%potscal = 1.0d0    !> scale automatically set wall potential by this factor

!>--- get CHRG and UHF if the respective files are present
  inquire (file='.CHRG',exist=ex)
  if (any(index(arg,'-chrg') .ne. 0)) ex = .false.
  if (ex) then
    write(stdout,*) '**ERROR** CREST will not read .CHRG files following version 3.0.2'
    write(stdout,*) 'Please use --chrg or the input file specifications.'
    error stop
  end if
  inquire (file='.UHF',exist=ex)
  if (any(index(arg,'-uhf') .ne. 0)) ex = .false.
  if (any(index(arg,'-mult') .ne. 0)) ex = .false.
  if (ex) then
    write(stdout,*) '**ERROR** CREST will not read .UHF files following version 3.0.2'
    write(stdout,*) 'Please use --uhf or the input file specifications.'
    error stop
  end if

!>--- options for constrained conformer sampling
  env%fixfile = 'none selected'

!>--- options for possible property calculations, mainly protonation/deprotonation/taut. tool
  env%protb%ewin = 30.0_wp  !> 30 kcal for protonation
  env%protb%swat = 0
  env%protb%swchrg = 0
  env%protb%iter = 1             !> number of iteration cycles for tautomerization
  env%protb%swelem = .false.     !> replace H⁺ in protonation routine by something else?
  env%protb%allowFrag = .false.  !> allow dissociated Structures?
  env%protb%threshsort = .false. !> use ewin threshold window
  env%protb%protdeprot = .false. !> (tautomerize) do first protonation and then deprotonation
  env%protb%deprotprot = .false. !> (tautomerize) do first deprotonation and then protonation
  env%protb%strictPDT = .false.  !> strict mode (i.e. bond constraints) for (de)protonation,tautomerization
  env%protb%ffopt = .true. !> FF pre-optimization (CREST>3.0.2)
  env%pclean = .false.       !> cleanup option for property mode

!>--- options for principal component analysis (PCA) and clustering
  env%pcmeasure = 'dihedral'

!>--- Standard THERMO options
  env%thermo%trange(1) = 278.15d0  !> T start
  env%thermo%trange(2) = 380.0d0   !> T stop (approx.)
  env%thermo%trange(3) = 10.0d0    !> T step
  env%thermo%ptot = 0.9d0          !> for hessians take x% conformers
  env%thermo%pcap = 50000          !> limit number of structures
  env%thermo%sthr = 25.0d0         !> rotor cutoff
  env%thermo%fscal = 1.0d0         !> frequency scaling factor
  env%thermo%emodel = 'grimme'     !> Svib treatment

!>--- other things
  env%crest_ohess = .false.

!>--- options for QCG
  env%cff = .true.
  env%nqcgclust = 0
  env%freq_scal = 0.75
  env%freqver = '--gfn2'
  env%max_solv = 150
  env%solv_file = ''
  env%solu_file = ''

!>--- options for msreact
  env%msei = .true.
  env%mscid = .false.
  env%msiso = .false.  ! msiso and msnoiso are mutually exclusive !!!
  env%msnoiso = .false.
  env%msmolbar = .false.
  env%mslargeprint=.false. ! dont remove temporary files
  env%msattrh=.true. ! add attractive potential for H-atoms
  env%msnbonds = 3 ! distance of bonds up to nonds are stretched
  env%msnshifts = 0 ! number of random shifts applied to whole mol
  env%msnshifts2 = 0 ! number of random shifts applied to whole mol
  env%msnfrag = 0 !number of fragments that are printed in msreact mode

!&>

!=========================================================================================!
!=========================================================================================!
!> MAIN RUNTYPE SELECTION VIA CMD
!=========================================================================================!
!=========================================================================================!
!>--- get the CREST version/runtype
  env%crestver = crest_none !> no runtype selected — must be set explicitly
  env%runver = 1            !> default
  env%properties = p_none   !> additional calculations/options before or after confsearch
  env%properties2 = p_none  !> backup for env%properties
  env%iterativeV2 = .true.  !> iterative crest V2 version
  env%preopt = .true.
!>--- check for (TOML) input file
  call find_input_file(arg,nra,idum)
  if (idum .ne. 0) then
    call parseinputfile(env,trim(arg(idum)))
    processedarg(idum) = .true.
  end if

!>--- first arg loop
  do i = 1,nra
    if (processedarg(i)) cycle
    argument = trim(arg(i))
    arg1 = ''
    if (i+1 .le. nra) arg1 = trim(arg(i+1))
    arg2 = ''
    if (i+2 .le. nra) arg2 = trim(arg(i+2))
    arg3 = ''
    if (i+3 .le. nra) arg3 = trim(arg(i+3))
    if (argument(1:2) == '--') then
      argument = argument(2:)
    end if
    if (argument .ne. '') then
      select case (argument) !> RUNTYPES

      case ('-v1') !> confscript version 1 (MF-MD-GC)
        processedarg(i) = .true.
        env%crestver = crest_mfmdgc
        write (stdout,'(2x,a,'' : MF-MD-GC'')') trim(arg(i))
        env%mdtime = 40.0d0       !> simulation length of the MD, 40ps total (2*20ps)(default for QMDFF would be 500)
        env%temps = 1             !> number of default MD cycles
        env%Maxrestart = 15
        env%performModef = .true. !> do the MF in V1
        env%trackorigin = .false. !> for v1 there is not much insight from this
        call parseflags_deprecated(argument)
        call creststop(status_safety)
        exit

      case ('-v2') !> confscript version 2 (MTD-GC)
        processedarg(i) = .true.
        env%crestver = crest_imtd
        write (stdout,'(2x,a,'' : MTD-GC'')') trim(arg(i))
        env%iterativeV2 = .false.  !> iterative crest V2 version
        env%Maxrestart = 1       !> for non-iterative MTD-GC only
        exit

      case ('-v3','-v2i','-imtdgc') !> confscript version 2 but iterativ (iMTD-GC)
        processedarg(i) = .true.
        env%crestver = crest_imtd
        env%iterativeV2 = .true.
        write (stdout,'(2x,a,'' : iMTD-GC'')') trim(arg(i))
        exit

      case ('-entropy') !> sMTD-iMTD+entropy extrapol
        processedarg(i) = .false.  !> THIS IS IMPORTANT, WE READ FURTHER ENTROPY SETTINGS BELOW
        env%crestver = crest_imtd
        env%iterativeV2 = .true.
        write (stdout,'(2x,a,'' : iMTD-sMTD + entropy extrapolation'')') trim(arg(i))
        exit

      case ('-v4') !> sMTD-iMTD (same as entropy mode)
        processedarg(i) = .true.
        env%crestver = crest_imtd2
        env%iterativeV2 = .true.
        env%entropymd = .true.
        env%rotamermds = .false.
        env%performCross = .false.
        env%emtd%maxfallback = 1
        write (stdout,'(2x,a,'' : iMTD-sMTD'')') trim(arg(i))
        exit

      case ('-mdopt','-purge') !> MDOPT
        processedarg(i) = .true.
        env%crestver = crest_mdopt
        atmp = ''
        env%preopt = .false.
        env%ensemblename = 'none selected'
        if (nra .ge. (i+1)) atmp = adjustl(arg(i+1))
        if ((atmp(1:1) /= '-').and.(len_trim(atmp) .ge. 1)) then
          processedarg(i+1) = .true.
          env%ensemblename = trim(atmp)
          call xyz2coord(env%ensemblename,'coord') !> write coord from lowest structure
          env%inputcoords = env%ensemblename !> just for a printout
        end if
        exit

      case ('-screen')  !> SCREEN
        processedarg(i) = .true.
        env%crestver = crest_screen
        atmp = ''
        env%ensemblename = 'none selected'
        if (nra .ge. (i+1)) atmp = adjustl(arg(i+1))
        if ((atmp(1:1) /= '-').and.(len_trim(atmp) .ge. 1)) then
          processedarg(i+1) = .true.
          env%ensemblename = trim(atmp)
          call xyz2coord(env%ensemblename,'coord') !write coord from lowest structure
          env%inputcoords = env%ensemblename !just for a printout
        end if
        exit

      case ('-mdsp','-ensemblesp') !> Singlepoints along ensemble
        processedarg(i) = .true.
        env%crestver = crest_ensemblesp
        atmp = ''
        env%preopt = .false.
        env%ensemblename = 'none selected'
        if (nra .ge. (i+1)) atmp = adjustl(arg(i+1))
        if ((atmp(1:1) /= '-').and.(len_trim(atmp) .ge. 1)) then
          processedarg(i+1) = .true.
          env%ensemblename = trim(atmp)
          call xyz2coord(env%ensemblename,'coord') !> write coord from lowest structure
          env%inputcoords = env%ensemblename !> just for a printout
        end if
        exit

      case ('-mdhess','-ensemblehess') !> Hessians + thermochemistry along ensemble
        processedarg(i) = .true.
        env%crestver = crest_ensemblehess
        atmp = ''
        env%preopt = .false.
        env%ensemblename = 'none selected'
        if (nra .ge. (i+1)) atmp = adjustl(arg(i+1))
        if ((atmp(1:1) /= '-').and.(len_trim(atmp) .ge. 1)) then
          processedarg(i+1) = .true.
          env%ensemblename = trim(atmp)
          call xyz2coord(env%ensemblename,'coord') !> write coord from lowest structure
          env%inputcoords = env%ensemblename !> just for a printout
        end if
        exit

      case ('-pka','-pKa')  !> pKa calculation script
        processedarg(i) = .true.
        env%crestver = crest_pka
        env%runver = 33
        !env%relax=.true.
        env%performCross = .false.    !skip the genetic crossing
        env%trackorigin = .false.
        env%Maxrestart = 1
        env%protb%ewin = 15.0d0
        env%gbsa = .true.
        env%solv = '--alpb h2o'
        env%protb%h_acidic = 0
        call pka_argparse(arg(i+1),env%protb%h_acidic)
        if (env%protb%h_acidic == -2) env%protb%pka_baseinp = arg1

      case ('-compare')   !> flag for comparing two ensembles, analysis tool
        processedarg(i) = .true.
        env%compareens = .true.
        env%crestver = crest_sorting
        env%sortmode = 'compare'
        env%preopt = .false.
        env%ensemblename = 'none selected'
        env%ensemblename2 = 'none selected'
        if (nra .ge. (i+2)) then
          atmp = adjustl(arg(i+1))
          btmp = adjustl(arg(i+2))
        else
          write (stdout,'(a,a)') trim(arg(i)),' requires two arguments:'
          write (stdout,'(2x,a,a)') trim(arg(i)),' [ensemble1] [ensemble2]'
          error stop
        end if
        if ((atmp(1:1) /= '-').and.(len_trim(atmp) .ge. 1).and. &
        &  (btmp(1:1) /= '-').and.(len_trim(btmp) .ge. 1)) then
          env%ensemblename = trim(atmp)
          env%ensemblename2 = trim(btmp)
          processedarg(i+1) = .true.
          processedarg(i+2) = .true.
        end if
        write (stdout,'(1x,a,1x,a,1x,a)') trim(arg(i)),trim(env%ensemblename),trim(env%ensemblename2)
        exit

      case ('-protonate') !> protonation tool
        processedarg(i) = .true.
        env%properties = p_protonate
        env%crestver = crest_protonate
        env%legacy = .false.
        write (stdout,'(2x,a,'' : automated protonation script'')') trim(arg(i))
        exit

      case ('-deprotonate') !> deprotonation tool
        processedarg(i) = .true.
        env%properties = p_deprotonate
        env%crestver = crest_deprotonate
        env%legacy = .false. 
        write (stdout,'(2x,a,'' : automated deprotonation script'')') trim(arg(i))
        exit

      case ('-tautomerize') !> tautomerization tool
        processedarg(i) = .true.
        env%properties = p_tautomerize
        env%crestver = crest_tautomerize
        env%legacy = .false. 
        write (stdout,'(2x,a,'' : automated tautomerization script'')') trim(arg(i))
        exit

      case ('-isomerize','-stereomers') !> isomerization tool
        processedarg(i) = .true.
        env%properties = p_isomerize
        write (stdout,'(2x,a,'' : automated stereoisomerization script'')') trim(arg(i))
        write (stdout,'(2x,''Note: Use of GFN-FF required for stereoisomer generation.'')')
        exit

      case ('-forall','-for') !> property mode with ensemble as input (deprecated)
        call parseflags_deprecated(argument)
        exit

      case ('-rrhoav')  !> Hessians along given ensemble and average
        processedarg(i) = .true.
        env%properties = p_rrhoaverage
        atmp = ''
        env%ensemblename = 'none selected'
        if (nra .ge. (i+1)) atmp = adjustl(arg(i+1))
        if ((atmp(1:1) /= '-').and.(len_trim(atmp) .ge. 1)) then
          env%ensemblename = trim(atmp)
          processedarg(i+1) = .true.
        end if
        inquire (file=env%ensemblename,exist=ex)
        if (.not.ex) then
          write (stdout,'(1x,a,a,a)') 'invalid ensemble file <',trim(env%ensemblename),'>. exit.'
          error stop
        end if
        exit

      case ('-reactor')  !> xtb nanoreactor workarounds
        processedarg(i) = .true.
        env%preopt = .false.
        env%crestver = crest_nano
        exit

      case ('-solvtool','-qcg')
        processedarg(i) = .true.
        !> Set solute file if present
        if (i == 2) then
          env%solu_file = trim(arg(i-1))
          inquire (file=env%solu_file,exist=ex)
          if (ex) then
            processedarg(i-1) = .true.
          end if
        end if
        !> Set solvent file if prensent
        !> If it is another argument, it doesent matter as solvent file is checke in solvtool
        if (nra >= i+1) then
          env%solv_file = arg1
          inquire (file=env%solv_file,exist=ex)
          if (ex) then
            processedarg(i+1) = .true.
          end if
        end if
        !> Set QCG defaults
        env%preopt = .false.
        env%crestver = crest_solv
        env%QCG = .true.
        env%runver = 3
        env%performCross = .false.
        env%optlev = 0 !> If QCG is invoked, optlevel default is normal
        env%properties = p_qcg
        env%ewin = 3.0d0
        env%doOHflip = .false. !> Switch off OH-flip
        if (env%iterativeV2) env%iterativeV2 = .false.
        exit
        !env%legacy = .true. !> force legacy routines for now

      case ('-compress')
        processedarg(i) = .true.
        env%crestver = crest_compr
        env%runver = 77
        env%mdstep = 2.5d0
        env%mddump = 2000
        env%autozsort = .false.
        exit

      case ('-msreact')
        processedarg(i) = .true.
        env%crestver = crest_msreac
        env%preopt = .false.
        env%presp = .true.
        env%ewin = 200.0d0 !> 200 kcal for msreact

      case ('-splitfile')
        processedarg(i) = .true.
        ctmp = arg1
        processedarg(i+1) = .true.
        k = huge(j)
        l = 1
        if (nra >= i+2) then
          read (arg(i+2),*,iostat=io) j
          if (io == 0) then
            k = j
            processedarg(i+2) = .true.
          end if
        end if
        if (nra >= i+3) then
          read (arg(i+3),*,iostat=io) j
          if (io == 0) then
            l = j
            processedarg(i+3) = .true.
          end if
        end if
        call splitfile(ctmp,k,l)
        call creststop(status_normal)

      case ('-printaniso')
        processedarg(i) = .true.
        ctmp = arg1
        inquire (file=ctmp,exist=ex)
        if (ex) then
          processedarg(i+1) = .true.
          call printaniso(ctmp,0.01_wp,0.025_wp,0.5_wp)
        end if
        stop

      case ('-rotalign')
        processedarg(i) = .true.
        ctmp = arg1
        inquire (file=ctmp,exist=ex)
        if (ex) then
          processedarg(i+1) = .true.
          call rotalign_tool(ctmp)
        end if
        stop

      case ('-printboltz')
        processedarg(i) = .true.
        if (nra >= i+2) then
          ctmp = arg1
          dtmp = arg2
          call prbweight(ctmp,dtmp)
          processedarg(i+1) = .true.
          processedarg(i+2) = .true.
        else
          ctmp = arg1
          call prbweight(ctmp,'')
          processedarg(i+1) = .true.
        end if

      case ('-wbotopo','-usewbo')  !> try to use a WBO file in topology analysis
        processedarg(i) = .true.
        ctmp = arg1
        if (ctmp(1:1) .ne. '-'.and.(nra >= i+1)) then
          env%wbofile = trim(ctmp)
          processedarg(i+1) = .true.
        else
          env%wbofile = 'wbo'
        end if
        env%wbotopo = .true.

      case ('-testtopo')
        processedarg(i) = .true.
        ctmp = arg1
        inquire (file=ctmp,exist=ex)
        if (i+2 .le. nra) then
          dtmp = arg2
          if (dtmp(1:1) == '-') then
            dtmp = 'default'
          else
            processedarg(i+2) = .true.
          end if
        end if
        if (ex) then
          processedarg(i+1) = .true.
          call testtopo(ctmp,env,dtmp)
        end if

      case ('-resortensemble')
        processedarg(i) = .true.
        ctmp = arg1
        inquire (file=ctmp,exist=ex)
        if (ex) then
          processedarg(i+1) = .true.
          call resort_ensemble(ctmp)
        end if
        stop

      case ('-thermo','-thermotool')
        processedarg(i) = .true.
        env%properties = p_thermo
        ctmp = trim(arg(1))  ! first argument to read the structure
        if (ctmp(1:1) .ne. '-') then
          env%inputcoords = trim(ctmp)
          env%thermo%coords = trim(ctmp)
        end if
        ctmp = arg1 ! second argument to read the vibspectrum
        if (ctmp(1:1) .ne. '-') then
          processedarg(i+1) = .true.
          env%thermo%vibfile = trim(ctmp)
        end if

      case ('-rmsd','-rmsdheavy','-hrmsd')
        processedarg(i) = .true.
        if ((argument == '-rmsdheavy').or.(argument == '-hrmsd')) then
          env%sortmode = 'hrmsd'
        else
          env%sortmode = 'rmsd'
        end if
        ctmp = arg1
        dtmp = arg2
        env%preopt = .false.
        env%crestver = crest_sorting
        inquire (file=ctmp,exist=ex)
        if (ex) then
          env%inputcoords = ctmp
          env%ensemblename = ctmp
          processedarg(i+1) = .true.
        end if
        inquire (file=dtmp,exist=ex)
        if (ex) then
          env%ensemblename2 = dtmp
          processedarg(i+2) = .true.
        end if

      case ('-irmsd','-irmsd_noinv')
        processedarg(i) = .true.
        ctmp = arg1
        dtmp = arg2
        env%preopt = .false.
        env%crestver = crest_sorting
        env%sortmode = 'irmsd'
        inquire (file=ctmp,exist=ex)
        if (ex) then
          env%inputcoords = ctmp
          env%ensemblename = ctmp
          processedarg(i+1) = .true.
        end if
        inquire (file=dtmp,exist=ex)
        if (ex) then
          env%ensemblename2 = dtmp
          processedarg(i+2) = .true.
        end if
        if (index(argument,'_noinv') .ne. 0) then
          env%iinversion = 2
        end if

      case ('-hungarian','-hungarianheavy','-hhungarian','-lsap','-hlsap','-lsapheavy')
        processedarg(i:i+2) = .true.
        ctmp = arg1
        dtmp = arg2
        if ((argument == '-hungarianheavy').or.(argument == '-hhungarian').or. &
           &(argument == '-lsapheavy').or.(argument == '-hlsap')) then
          call quick_hungarian_match(ctmp,dtmp,.true.)
        else
          call quick_hungarian_match(ctmp,dtmp,.false.)
        end if
        stop

      case ('-symmetries')
        processedarg(i) = .true.
        ctmp = arg1
        inquire (file=ctmp,exist=ex)
        if (ex) then
          processedarg(i+1) = .true.
          call ensemble_analsym(trim(ctmp),.true.)
        end if
        call exit(0)

      case ('-exlig','-exligand','-exchligand')
        processedarg(i) = .true.
        env%properties = p_ligand
        env%protb%infile = trim(arg(1))
        ctmp = arg1
        processedarg(i+1) = .true.
        env%protb%newligand = trim(ctmp)
        read (arg(i+2),*,iostat=io) j
        if (io == 0) then
          env%protb%centeratom = j
          processedarg(i+2) = .true.
        end if
        read (arg(i+3),*,iostat=io) j
        if (io == 0) then
          env%protb%ligand = j
          processedarg(i+3) = .true.
        end if
        exit

      case ("-acidbase","-ab",'-abprep','-pkaprep','-gdissprep')  !-- acid base correction
        processedarg(i) = .true.
        !> crest --ab <acid.xyz> <base.xyz> --chrg <acidchrg>
        env%properties = p_acidbase
        if (index(arg(i),'prep') .ne. 0) then
          call pka_argparse2(env,arg(i+1),arg(i+2),env%protb%pka_mode)
        else
          ctmp = arg1
          inquire (file=ctmp,exist=ex)
          if (ex) then
            processedarg(i+1) = .true.
            env%protb%pka_acidensemble = trim(ctmp)
            write (stdout,'(1x,a,a)') 'File used for the acid: ',trim(ctmp)
          end if
          ctmp = arg2
          inquire (file=ctmp,exist=ex)
          if (ex) then
            processedarg(i+2) = .true.
            env%protb%pka_baseensemble = trim(ctmp)
            write (stdout,'(1x,a,a)') 'File used for the base: ',trim(ctmp)
          end if
        end if
        env%solv = '--alpb h2o'
        env%gfnver = '--gfn2'

      case ('-redoextrapol')
        processedarg(i) = .true.
        ctmp = arg1
        processedarg(i+1) = .true.
        read (arg(i+2),*,iostat=io) j
        if (io == 0) then
          processedarg(i+2) = .true.
          call redo_extrapol(ctmp,j)
        else
          call redo_extrapol(ctmp,0)
        end if
        stop

      case ('-sp') !> singlepoint calculation (uses new calculator routines)
        processedarg(i) = .true.
        env%crestver = crest_sp
        env%preopt = .false.
        env%legacy = .false.
        write (stdout,'(2x,a,t15,a)') argument//':','Singlepoint energy calculation runtype'
        exit

      case ('-opt','-optimize','-ancopt','-ohess') !> ANCOPT structure optimization (uses new calculator routines)
        processedarg(i) = .true.
        env%preopt = .false.
        env%crestver = crest_optimize
        env%legacy = .false.
        if (argument .eq. '-ohess') then
          env%crest_ohess = .true.
          write (stdout,'(2x,a,t15,a)') argument//':','Geometry optimization + frequency calculation runtype'
        else
          write (stdout,'(2x,a,t15,a)') argument//':','Geometry optimization runtype'
        end if
        !if (i+1 .le. nra) then
        env%optlev = optlevnum(arg(i+1),iostat=io)
        if (io == 0) processedarg(i+1) = .true.
        !end if

        exit

      case ('-hess','-numhess') !> Numerical hessian
        processedarg(i) = .true.
        env%preopt = .false.
        env%crestver = crest_numhessian
        env%legacy = .false.
        write (stdout,'(2x,a,t15,a)') argument//':','Frequency calculation runtype'
        exit

      case ('-trialopt')  !> test optimization with topocheck
        processedarg(i) = .true.
        env%preopt = .false.
        env%crestver = crest_trialopt
        write (stdout,'(2x,a,t15,a)') argument//':','Trial geometry optimization'
        exit

      case ('-ttconf')  !> tensor-train (light) conformer search
        processedarg(i) = .true.
        env%preopt = .true.
        env%crestver = crest_ttc
        env%runver = crest_ttc
        env%legacy = .false.
        !>--- optional argument: a preset (fast/normal/accurate) OR "bruteforce"
        if (arg1 .ne. ''.and.arg1(1:1) .ne. '-') then
          select case (arg1)
          case ('bruteforce','brute','oracle')
            env%ttconf%use_sweep = .false.
            processedarg(i+1) = .true.
          case default
            call env%ttconf%setpreset(arg1,presetok)
            if (presetok) processedarg(i+1) = .true.
          end select
        end if
        if (env%ttconf%use_sweep) then
          write (stdout,'(2x,a,t15,a)') argument//':', &
            & 'Tensor-train conformer search ['//trim(env%ttconf%preset)//']'
        else
          write (stdout,'(2x,a,t15,a)') argument//':', &
            & 'Tensor-train conformer search [brute-force]'
        end if
        exit

      case ('-dynamics','-dyn') !> molecular dynamics (uses new calculator routines)
        processedarg(i) = .true.
        env%preopt = .false.
        env%crestver = crest_moldyn
        env%legacy = .false.
        write (stdout,'(2x,a,t15,a)') argument//':','Molecular dynamics simulation'
        exit

      case ('-sort','-cregen')
        processedarg(i) = .true.
        env%preopt = .false.
        env%crestver = crest_sorting
        env%autozsort = .false.
        if (argument == '-cregen') then
          env%sortmode = 'cregen'
          env%confgo = .true.
        end if
        ctmp = arg1
        inquire (file=ctmp,exist=ex)
        if (ex) then
          processedarg(i+1) = .true.
          env%inputcoords = ctmp
          env%ensemblename = ctmp
        end if
        if (argument == '-sort'.and.nra >= i+2) then
          ctmp = arg2
          if (ctmp(1:1) .ne. '-') then
            processedarg(i+2) = .true.
            env%sortmode = trim(ctmp)
          end if
        end if

      case ('-bh','-GMIN')
        processedarg(i) = .true.
        env%crestver = crest_bh
        write (stdout,'(2x,a,t15,a)') argument//':','Basin-hopping global optimization'
        exit

      case ('-SANDBOX')
        processedarg(i) = .true.
        !>--- readl vs readl_old test suite
        !>-----
        stop
      case ('-PLAYGROUND','-TEST')
        processedarg(i) = .true.
        env%preopt = .false.
        env%crestver = crest_test
        exit
      case default
        continue
      end select !> RUNTYPES
    end if
  end do

!=========================================================================================!
!=========================================================================================!
!=========================================================================================!
!>--- options for the xtb Nano-reactor
  if (any((/crest_nano,crest_compr/) == env%crestver)) then
    env%rdens = -1.0_wp       !> reference density (i.e., =-1 means it is not set here)
    env%preactormtd = .false.  !> prepare reactor mtd?
    env%preactorpot = .false.  !> prepare reactor logfermi?
  end if

!>--- turn off autozsort for additional applications preceeding conf.searches
  if (env%properties .lt. 0) then
    env%autozsort = .false.
  end if

!>--- options for topology related applications.
  if (.not.allocated(env%wbofile)) then
    env%wbotopo = .false.
    env%wbofile = ''
  end if
!>    E.g. stereoisomer sampling or rotamer enhancement for entropy calc.
  if (env%properties == p_isomerize) then
    env%gfnver = '-gff'    !stereoisomer generation works only with GFN-FF!
  end if
  env%tsplit = 5   !timeframe splitting in conformational search for Entropy extrapolation

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!>    I N P U T   C O O R D I N A T E S
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
!> If env%inputcoords is initialized anywhere prior to
!> this point, it will be taken as the input.
!> Otherwise, the first cmd argument after "crest" will
!> taken for the input coordinates
  if (allocated(env%inputcoords)) then
    call inputcoords(env,env%inputcoords)
  else
    call inputcoords(env,trim(arg(1)))
    processedarg(1) = .true.
  end if
!> For sorting runtypes, fall back to the input file if no ensemble was set explicitly
  if (env%crestver == crest_sorting.and.len_trim(env%ensemblename) == 0) then
    env%ensemblename = env%inputcoords
  end if

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> parse the input flags
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
  do i = 1,nra
    if (processedarg(i)) cycle
    argument = trim(arg(i))
    arg1 = ''
    if (i+1 .le. nra) arg1 = trim(arg(i+1))
    arg2 = ''
    if (i+2 .le. nra) arg2 = trim(arg(i+2))
    arg3 = ''
    if (i+3 .le. nra) arg3 = trim(arg(i+3))
    if (argument(1:2) == '--') then
      argument = argument(2:)
    end if
    if (argument .ne. '') then

!========================================================================================!
!-------- switch between legacy (systemcall) and new code (API) implementations
!========================================================================================!
      select case (argument)
      case ('-legacy')      !> switch to old xtb-call version where possible
        processedarg(i) = .true.
        env%legacy = .true.
      case ('-newversion')  !> switch to newer implementations (CREST >3.0)
        processedarg(i) = .true.
        env%legacy = .false.
      end select
!========================================================================================!
!------- flags exclusively for V2 (MTD-GC) V2i/V3 (iMTD-GC) V4 (iMTD-sMTD)
!========================================================================================!
      if (any((/crest_imtd,crest_imtd2,11/) == env%crestver)) then
        select case (argument) !> V2

        case ('-quick')                           !> performing quick conformational search
          processedarg(i) = .true.
          env%quick = .true.
          env%runver = 2
          env%ewin = 5.0d0
          if (env%optlev > 1.0d0) env%optlev = 1.0d0    !> optlev tight for quick run

        case ('-mdskip')                          !> set skipping structures in -mdopt
          processedarg(i) = .true.
          call readl(arg1,xx,j)
          if (j > 0) then
            env%mdskip = nint(xx(1))
            processedarg(i+1) = .true.
          else
            call parseflags_missing(argument)
          end if

        case ('-nomtd')                           !> Don't do the MTD in V2
          processedarg(i) = .true.
          env%performMTD = .false.

        case ('-restartopt')                      !> go to step 2 of multilevel optimization immideatly
          processedarg(i) = .true.
          env%restartopt = .true.
          env%autozsort = .false.

        case ('-norotmd')                         !> don't do the regular mds after step 2 in multilevel optimization of V2
          processedarg(i) = .true.
          env%rotamermds = .false.

        case ('-rotmd')
          processedarg(i) = .true.
          env%rotamermds = .true.

        case ('-gcmopt')                          !> GC multilevel optimization activate in V2
          processedarg(i) = .true.
          env%gcmultiopt = .true.

        case ('-gcsopt')                          !> GC single level optimization in V2
          processedarg(i) = .true.
          env%gcmultiopt = .false.

        case ('-nogcmopt')                        !> GC single level optimization in V2
          processedarg(i) = .true.
          env%gcmultiopt = .false.

        case ('-qmdff')                           !> use QMDFF for the MDs in V2?
          processedarg(i) = .true.
          env%useqmdff = .true.
          call parseflags_deprecated(argument)

        case ('-nci')                             !> NCI special mode
          processedarg(i) = .true.
          write (stdout,'(2x,a,1x,a)') trim(arg(i)),' : Special NCI mode for non-covalently bound complexes or clusters.'
          env%NCI = .true.
          env%runver = 4
          env%autozsort = .false.
          env%performCross = .false.
          env%rotamermds = .false.

        case ('-squick','-superquick')            !> extremely crude quick mode
          processedarg(i) = .true.
          write (stdout,'(2x,a,1x,a)') trim(arg(i)),' : very crude quick-mode (no NORMMD, no GC, crude opt.)'
          env%rotamermds = .false.      !> no NORMMD
          env%performCross = .false.    !> no GC
          env%quick = .true.            !> MTD settings from the quick-mode
          env%superquick = .true.       !> use user-set opt level in Multilevel opt.
          env%runver = 5
          if (env%optlev > 0.0d0) env%optlev = 0.0d0    !> user-set opt level
          env%ewin = 5.0d0              !> smaller energy window

        case ('-mquick','-megaquick')   !> extremely crude quick mode pt.2
          processedarg(i) = .true.
          write (stdout,'(2x,a,1x,a)') trim(arg(i)),' : very crude quick-mode (no NORMMD, no GC, crude opt.)'
          env%rotamermds = .false.      !> no NORMMD
          env%performCross = .false.    !> no GC
          env%quick = .true.            !> MTD settings from the quick-mode
          env%superquick = .true.       !> use user-set opt level in Multilevel opt.
          env%Maxrestart = 1            !> only one MTD iteration
          env%runver = 6
          if (env%optlev > 0.0d0) env%optlev = 0.0d0  !> user-set opt level
          env%ewin = 2.5d0              !> smaller energy window

        case ('-extensive')   !> counterpart to quick mode
          processedarg(i) = .true.
          env%slow = .true.
          env%quick = .false.
          env%superquick = .false.
          env%optlev = 0.0d0
          env%ewin = 8.0d0
          env%runver = 8

        case ('-static','-staticmtd')
          processedarg(i) = .true.
          env%staticmtd = .true.

        case default
          continue
        end select !> V2
        !--- iterative version of V2
        if (env%iterativeV2) then
          select case (argument) !> V2i
          case ('-mrest')                  !> set max number of restarts
            processedarg(i) = .true.
            call readl(arg1,xx,j)
            if (j > 0) then
              env%Maxrestart = nint(xx(1))
              processedarg(i+1) = .true.
            else
              call parseflags_missing(argument)
            end if

          case ('-iru')                    !> re-use previously found conformers as bias in iterative approach
            processedarg(i) = .true.
            env%iru = .true.

          case ('-singlerun')              !> QCG special mode
            processedarg(i) = .true.
            write (stdout,'(2x,a,1x,a)') trim(arg(i)),' : run mode with only a single MTD and no iterations (for testing)'
            env%runver = 45
            env%Maxrestart = 1
            env%rotamermds = .false.

          case default
            continue
          end select !> V2i
          !-----
        end if
      end if
!========================================================================================!
!------- Settings for MDOPT and SCREEN
!========================================================================================!
      if (env%crestver == crest_mdopt.or.env%crestver == crest_screen) then
        select case (argument) !> SCREEN
        case ('-purge')        !> Purge special application
          processedarg(i) = .true.
          env%optpurge = .true.

        case ('-ethrpurge','-ethrp')
          processedarg(i) = .true.
          read (arg(i+1),*,iostat=io) rdum
          if (io == 0) then
            env%ethrpurge = rdum
            processedarg(i+1) = .true.
          end if

        case default
          continue
        end select !> SCREEN
      end if
!========================================================================================!
!------- Settings for NANO-REACTOR
!========================================================================================!
      if (env%crestver == crest_nano) then
        select case (argument) !> RCTR
        case ('-genpot')
          processedarg(i) = .true.
          if (i+1 .le. nra) then
            atmp = arg1
            if (atmp(1:1) .ne. '-') then
              call readl(arg1,xx,j)
              env%rdens = xx(1)
              processedarg(i+1) = .true.
            end if
          end if
          env%properties = p_reactorset
          env%preactorpot = .true.

        case ('-genmtd')
          processedarg(i) = .true.
          env%properties = p_reactorset
          env%mdtime = 20.0d0
          if (i+1 .le. nra) then
            atmp = arg1
            if (atmp(1:1) .ne. '-') then
              call readl(arg1,xx,j)
              env%mdtime = xx(1)
              processedarg(i+1) = .true.
            end if
          end if
          env%nmetadyn = 1
          if (.not.allocated(env%metadfac)) then
            allocate (env%metadfac(1))
            allocate (env%metadexp(1))
            allocate (env%metadlist(1))
          end if
          env%metadlist(1) = nint(env%mdtime)
          env%metadexp(1) = 1.00_wp
          env%metadfac(1) = 0.04_wp
          env%preactormtd = .true.

        case ('-fragopt')
          processedarg(i) = .true.
          env%restartopt = .true.

        case ('-iso')
          processedarg(i) = .true.
          env%riso = .true.

        case default
          continue
        end select !> RCTR
      end if
!========================================================================================!
!------- Flags for QCG
!========================================================================================!
      if (env%QCG) then
        select case (argument) !> QCG
        case ('-keeptmp')
          processedarg(i) = .true.
          env%keepModef = .true.

        case ('-nomtd')                           !> Don't do the MTD in V2
          processedarg(i) = .true.
          env%performMTD = .false.

        case ('-fixsolute')                       !> Fix the solute after CMA trafo
          processedarg(i) = .true.
          env%constrain_solu = .true.

        case ('-nofix')                           !> No fixing of the solute after CMA trafo
          processedarg(i) = .true.
          env%noconst = .true.

        case ('-restartopt')                      !> go to step 2 of multilevel optimization immideatly
          processedarg(i) = .true.
          env%restartopt = .true.
          env%autozsort = .false.

        case ('-norotmd')                         !> don't do the regular mds after step 2 in multilevel optimization of V2
          processedarg(i) = .true.
          env%rotamermds = .false.

        end select !> QCG
      end if

!========================================================================================!
!------- Flags for msreact
!========================================================================================!
      if (env%crestver == crest_msreac) then
        select case (argument) !> msreact
        case ('-msei')
          processedarg(i) = .true.
          env%msei = .true.

        case ('-mscid')
          processedarg(i) = .true.
          env%mscid = .true.
          env%msei = .false.

        case ('-msnoiso') !> filter out non fragmentated structures in msreact
          processedarg(i) = .true.
          env%msnoiso = .true.

        case ('-msiso') !> filter out fragmentated structures in msreact
          processedarg(i) = .true.
          env%msiso = .true.

        case ('-msnbonds') ! give number of bonds up to which bias potential is added between atoms default 3
          processedarg(i) = .true.
          call readl(arg1,xx,j)
          if (j > 0) then
            env%msnbonds = xx(1)
            processedarg(i+1) = .true.
          else
            call parseflags_missing(argument)
          end if

        case ('-msnshifts') ! give number of times atoms are randomly shifted before optimization
          processedarg(i) = .true.
          call readl(arg1,xx,j)
          if (j > 0) then
            env%msnshifts = xx(1)
            processedarg(i+1) = .true.
          else
            call parseflags_missing(argument)
          end if

        case ('-msnshifts2') ! give number of times atoms are randomly shifted before applying the constrained optimization default 0
          processedarg(i) = .true.
          call readl(arg1,xx,j)
          if (j > 0) then
            env%msnshifts2 = xx(1)
            processedarg(i+1) = .true.
          else
            call parseflags_missing(argument)
          end if

        case ('-msnfrag') ! give number of structures that should be generated
          processedarg(i) = .true.
          call readl(arg1,xx,j)
          if (j > 0) then
            env%msnfrag = xx(1)
            processedarg(i+1) = .true.
          else
            call parseflags_missing(argument)
          end if

        case ('-msmolbar') !> filter out structures with same molbar code in msreact
          processedarg(i) = .true.
          env%msmolbar = .true.

        case ('-msinchi') !> filter out structures with same inchi code in msreact
          processedarg(i) = .true.
          env%msinchi = .true.

        case ('-msnoattrh') !> add attractive potential for H-atoms
          processedarg(i) = .true.
          env%msattrh = .false.

        case ('-mslargeprint') !> additional printouts and keep MSDIR
          processedarg(i) = .true.
          env%mslargeprint = .true.

        case ('-msinput') ! msreact input file
          processedarg(i) = .true.
          ctmp = arg1
          if (ctmp(1:1) .ne. '-') then
            env%msinput = trim(ctmp)
            processedarg(i+1) = .true.
          else

          end if
        end select !> msreact
      end if
!========================================================================================!
!------- other general flags
!========================================================================================!
      select case (argument) !> ARGPARSER1
      case ('-dry')             !> "dry" run to print settings
        processedarg(i) = .true.
        env%dryrun = .true.

      case ('-ttnocache')       !> disable the TTConf-light energy cache
        processedarg(i) = .true.
        env%ttconf%use_cache = .false.
        write (stdout,'(2x,a,1x,a)') trim(arg(i)), &
          & ' : TTConf-light energy cache disabled.'

      case ('-ttsp')            !> TTConf-light: singlepoints only (no opt)
        processedarg(i) = .true.
        env%ttconf%sp_only = .true.
        write (stdout,'(2x,a,1x,a)') trim(arg(i)), &
          & ' : TTConf-light evaluates singlepoints only (no geometry optimization).'

      case ('-ttringbonds')     !> allow in-ring bonds as TTConf TT variables
        processedarg(i) = .true.
        env%ttconf%excl_rings = .false.
        write (stdout,'(2x,a,1x,a)') trim(arg(i)), &
          & ' : TTConf-light treats in-ring bonds as TT variables.'

      case ('-ttrings')         !> sample ring templates as TTConf TT sites
        processedarg(i) = .true.
        env%ttconf%ring_sample = .true.
        write (stdout,'(2x,a,1x,a)') trim(arg(i)), &
          & ' : TTConf-light samples ring templates as TT sites.'

      case ('-ttseed')          !> fix the TTConf-light RNG seed (reproducible)
        processedarg(i) = .true.
        ctmp = arg1
        if (ctmp(1:1) .ne. '-'.and.len_trim(ctmp) .ge. 1) then
          read (ctmp,*,iostat=io) idum
          if (io .eq. 0) then
            env%ttconf%seed = idum
            processedarg(i+1) = .true.
          end if
        end if
        write (stdout,'(2x,a,1x,a,i0,a)') trim(arg(i)), &
          & ' : TTConf-light RNG seed = ',env%ttconf%seed,'.'

      case ('-ttrank')          !> TTConf-light TT rank r
        processedarg(i) = .true.
        ctmp = arg1
        if (ctmp(1:1) .ne. '-'.and.len_trim(ctmp) .ge. 1) then
          read (ctmp,*,iostat=io) idum
          if (io .eq. 0) then
            env%ttconf%rank = idum
            processedarg(i+1) = .true.
          end if
        end if
        write (stdout,'(2x,a,1x,a,i0,a)') trim(arg(i)), &
          & ' : TTConf-light TT rank = ',env%ttconf%rank,'.'

      case ('-ttsweeps')        !> TTConf-light number of sweeps s
        processedarg(i) = .true.
        ctmp = arg1
        if (ctmp(1:1) .ne. '-'.and.len_trim(ctmp) .ge. 1) then
          read (ctmp,*,iostat=io) idum
          if (io .eq. 0) then
            env%ttconf%sweeps = idum
            processedarg(i+1) = .true.
          end if
        end if
        write (stdout,'(2x,a,1x,a,i0,a)') trim(arg(i)), &
          & ' : TTConf-light sweeps = ',env%ttconf%sweeps,'.'

      case ('-ttgrid')          !> TTConf-light dihedral grid points
        processedarg(i) = .true.
        ctmp = arg1
        if (ctmp(1:1) .ne. '-'.and.len_trim(ctmp) .ge. 1) then
          read (ctmp,*,iostat=io) idum
          if (io .eq. 0) then
            env%ttconf%ngrid = idum
            processedarg(i+1) = .true.
          end if
        end if
        write (stdout,'(2x,a,1x,a,i0,a)') trim(arg(i)), &
          & ' : TTConf-light dihedral grid points = ',env%ttconf%ngrid,'.'

      case ('-ttninit')         !> TTConf-light number of random tail seeds
        processedarg(i) = .true.
        ctmp = arg1
        if (ctmp(1:1) .ne. '-'.and.len_trim(ctmp) .ge. 1) then
          read (ctmp,*,iostat=io) idum
          if (io .eq. 0) then
            env%ttconf%ninit = idum
            processedarg(i+1) = .true.
          end if
        end if
        write (stdout,'(2x,a,1x,a,i0,a)') trim(arg(i)), &
          & ' : TTConf-light initial tail seeds = ',env%ttconf%ninit,'.'

      case ('-ttewin')          !> TTConf-light conformer energy window (kcal/mol)
        processedarg(i) = .true.
        ctmp = arg1
        if (ctmp(1:1) .ne. '-'.and.len_trim(ctmp) .ge. 1) then
          read (ctmp,*,iostat=io) rdum
          if (io .eq. 0) then
            env%ttconf%ewin = rdum
            processedarg(i+1) = .true.
          end if
        end if
        write (stdout,'(2x,a,1x,a,f0.2,a)') trim(arg(i)), &
          & ' : TTConf-light energy window = ',env%ttconf%ewin,' kcal/mol.'

      case ('-ttkt')            !> TTConf-light maxvol weight temperature (kcal/mol)
        processedarg(i) = .true.
        ctmp = arg1
        if (ctmp(1:1) .ne. '-'.and.len_trim(ctmp) .ge. 1) then
          read (ctmp,*,iostat=io) rdum
          if (io .eq. 0) then
            env%ttconf%kt = rdum
            processedarg(i+1) = .true.
          end if
        end if
        write (stdout,'(2x,a,1x,a,f0.2,a)') trim(arg(i)), &
          & ' : TTConf-light maxvol temperature = ',env%ttconf%kt,' kcal/mol.'

      case ('-ttringmethod')    !> TTConf-light ring-conformation generator
        processedarg(i) = .true.
        ctmp = arg1
        if (ctmp(1:1) .ne. '-'.and.len_trim(ctmp) .ge. 1) then
          env%ttconf%ring_method = trim(ctmp)
          processedarg(i+1) = .true.
        end if
        write (stdout,'(2x,a,1x,a,a,a)') trim(arg(i)), &
          & ' : TTConf-light ring method = ',trim(env%ttconf%ring_method),'.'

      case ('-nozs')
        processedarg(i) = .true.
        env%autozsort = .false.   !> turn off automatic zsort (default)

      case ('-zs')
        processedarg(i) = .true.
        env%autozsort = .true.    !> turn on automatic zsort

      case ('-nocross')
        processedarg(i) = .true.
        env%performCross = .false.    !> skip the genetic crossing
        write (stdout,'(2x,a,1x,a)') trim(arg(i)),' : skipping GC part.'

      case ('-cross')
        processedarg(i) = .true.
        env%performCross = .true.     !> do the genetic crossing
        env%autozsort = .true.

      case ('-keepdir','-keeptmp')     !> Do not delete temporary directories at the end
        processedarg(i) = .true.
        env%keepModef = .true.

      case ('-opt','-optlev')             !> settings for optimization level of GFN-xTB
        processedarg(i) = .true.
        if (i+1 .le. nra) then
          env%optlev = optlevnum(arg(i+1),iostat=io)
          if (io == 0) processedarg(i+1) = .true.
        end if
        write (stdout,'(2x,a,1x,a)') trim(arg(i)),optlevflag(env%optlev)

      case ('-gfn','-gfn1','-gfn2','-gfn0','-gff','-gfnff')
        processedarg(i) = .true.
        ctmp = argument
        if (argument == '-gfn') then
          processedarg(i+1) = .true.
          dtmp = arg1
          ctmp = ctmp//dtmp
        end if
        if (env%properties == p_isomerize) then
          ctmp = 'stereoisomers'
        end if
        select case (ctmp) !> GFN
        case ('-gfn1')
          env%gfnver = '--gfn1'
          write (stdout,'(2x,a,'' : Use of GFN1-xTB requested.'')') env%gfnver
        case ('-gfn2')
          env%gfnver = '--gfn2'
          write (stdout,'(2x,a,'' : Use of GFN2-xTB requested.'')') env%gfnver
        case ('-gfn0')
          env%gfnver = '--gfn0'
          write (stdout,'(2x,a,'' : Use of GFN0-xTB requested.'')') env%gfnver
        case ('-gff','-gfnff')
          env%gfnver = '--gff'
          write (stdout,'(2x,a,'' : Use of GFN-FF requested.'')') '--gfnff'
          env%mdstep = 1.5d0
          env%hmass = 5.0d0
          ctype = 5 !> bond constraint activated
          if (any((/crest_imtd,crest_imtd2/) == env%crestver)) then
            bondconst = .true.
          end if
          env%cts%cbonds_md = .true.
          env%checkiso = .true.
        case ('stereoisomers')
          env%gfnver = '--gff'
        case default
          env%gfnver = '--gfn2'
        end select !> GFN

      case ('-gxtb')
        processedarg(i) = .true.
        env%gfnver = '--gxtb'
        write (stdout,'(2x,a,'' : Use of g-xTB requested.'')') env%gfnver
        call gxtb_syscall_warning()

      case ('-gxtb_dev')
        processedarg(i) = .true.
        env%gfnver = '--gxtb'
        write (stdout,'(2x,a)') 'Note: --gxtb_dev is deprecated, redirecting to --gxtb.'

      case ('-mlip') !> quick setup of an fmlip-relay ML potential: --mlip <uma|maceoff>
        processedarg(i) = .true.
        env%legacy = .false. !> new calculators only!
        select case (trim(arg1))
        case ('uma')
          env%gfnver = 'uma'
          processedarg(i+1) = .true.
          write (stdout,'(2x,a)') '--mlip uma : FairChem UMA (omol task) via fmlip-relay requested.'
        case ('maceoff','mace-off','mace_off')
          env%gfnver = 'maceoff'
          processedarg(i+1) = .true.
          write (stdout,'(2x,a)') '--mlip maceoff : MACE-OFF23 (medium) via fmlip-relay requested.'
        case default
          write (stdout,'(/,2x,a)') '┌─[ --mlip: missing or unknown backend ]'
          write (stdout,'(2x,a)')   '│ The --mlip flag selects an fmlip-relay ML potential'
          write (stdout,'(2x,a)')   '│ and expects one backend keyword:'
          write (stdout,'(2x,a)')   '│   uma      FairChem UMA foundation model (omol task)'
          write (stdout,'(2x,a)')   '│   maceoff  MACE-OFF23 organic force field (medium)'
          write (stdout,'(2x,a)')   '│ For finer control (model size, task, device, custom'
          write (stdout,'(2x,a)')   '│ checkpoints) use the TOML [[calculation.level]] block.'
          write (stdout,'(2x,a,/)') '└────────────────────────────────────────────────'
          call creststop(status_config)
        end select

      case ('-orca') !> set up a single ORCA driver level: --orca <template> <exe>
        processedarg(i) = .true.
        env%legacy = .false. !> new calculators only!
        ! ── both the template and the executable are mandatory ──────────
        if (arg1 == '' .or. arg1(1:1) == '-' .or. &
        &   arg2 == '' .or. arg2(1:1) == '-') then
          write (stdout,'(/,2x,a)') '┌─[ --orca: missing argument(s) ]'
          write (stdout,'(2x,a)')   '│ The --orca flag sets up an ORCA driver and expects'
          write (stdout,'(2x,a)')   '│ exactly two arguments, in this order:'
          write (stdout,'(2x,a)')   '│   1) an ORCA input template'
          write (stdout,'(2x,a)')   '│   2) the path to the ORCA executable'
          write (stdout,'(2x,a)')   '│ The template must define ORCA EnGrad jobs ONLY'
          write (stdout,'(2x,a)')   '│ (no opt/freq/md) - CREST drives the SP/opt loop.'
          write (stdout,'(2x,a)')   '│ Example:'
          write (stdout,'(2x,a)')   '│   crest struc.xyz --orca ORCA.in /path/to/orca'
          write (stdout,'(2x,a,/)') '└────────────────────────────────────────────────'
          call creststop(status_config)
        end if
        env%orca_template = trim(arg1)
        env%orca_cmd = trim(arg2)
        processedarg(i+1) = .true.
        processedarg(i+2) = .true.
        env%gfnver = '--orca'
        write (stdout,'(2x,a)') '--orca : Use of ORCA driver requested.'
        write (stdout,'(4x,a,1x,a)') 'template  :',trim(env%orca_template)
        write (stdout,'(4x,a,1x,a)') 'executable:',trim(env%orca_cmd)

      case ('-refine','-rsp','-ropt') !> add one refinement step (via cmd only one is possible)
        processedarg(i) = .true.
        env%legacy = .false. !> new calculators only!
        if (i+1 .le. nra) then
          env%gfnver2 = arg1
          write (stdout,'(2x,a,1x,a,a)') argument,trim(env%gfnver2), &
          & ' : adding refinement step (singlepoint on optimized structures)'
          processedarg(i+1) = .true.
        end if

      case ('-rerank') !> post-search SP reranking at a higher level (pqueue job 52)
        processedarg(i) = .true.
        env%legacy = .false.
        if (i+1 .le. nra) then
          env%rerank_lvl = lowercase(trim(arg1))
          processedarg(i+1) = .true.
          call env%addjob(p_prop_multilevel+2)
          write (stdout,'(2x,a,1x,a,a)') argument,trim(env%rerank_lvl), &
            & ' : post-search SP re-ranking of conformer ensemble'
        end if

      case ('-reopt') !> post-search geometry re-optimization (standalone, pqueue job 53)
        processedarg(i) = .true.
        env%legacy = .false.
        if (i+1 .le. nra) then
          env%reopt_lvl = lowercase(trim(arg1))
          processedarg(i+1) = .true.
          call env%addjob(p_prop_multilevel+3)
          write (stdout,'(2x,a,1x,a,a)') argument,trim(env%reopt_lvl), &
            & ' : post-search re-optimization of conformer ensemble'
        end if

      case ('-finalhess')
        processedarg(i) = .true.
        env%legacy = .false.
        call env%addjob(p_prop_finalhess)
        write (stdout,'(2x,a,a)') argument, &
          & ' : post-search Hessian + free-energy re-ranking of conformer ensemble'

      case default !> catch composite method arguments: A@B, A//B, A/sp/B, A/opt/B
        if (argument(1:1) == '-') then
          call parse_hybrid_argument(argument(2:),hybrid_quality,hybrid_workhorse, &
            &                        hybrid_mode,io)
          if (io == 0) then
            processedarg(i) = .true.
            env%legacy = .false.
            call setup_hybrid_calc(env,trim(hybrid_workhorse), &
              &                    trim(hybrid_quality),trim(hybrid_mode))
          end if
        end if

      case ('-charges') !> read charges from file for GFN-FF calcs.
        processedarg(i) = .true.
        ctmp = arg1
        if ((len_trim(ctmp) < 1).or.(ctmp(1:1) == '-')) then
          ctmp = 'charges'
        else
          processedarg(i+1) = .true.
        end if
        inquire (file=ctmp,exist=ex)
        if (ex) then
          env%chargesfilename = ctmp
          env%chargesfile = .true.
          write (stdout,'(2x,a,a,a)') '-charges: file <',trim(ctmp),'> used for atomic charges'
          call env%ref%rdcharges(env%chargesfilename,idum)
          if (idum .ne. env%chrg) then
            write (stdout,'(12x,a,i0)') 'with total summed up molecular charge: ',idum
            env%chrg = idum
            env%ref%ichrg = idum
          end if
        end if

      case ('-efield')  !> electric field in V/Ang, only compatibe with tblite
        processedarg(i) = .true.
        if (.not.allocated(env%ref%efield)) allocate (env%ref%efield(3),source=0.0_wp)
        if (i+3 .le. nra) then
          ctmp = arg1
          read (ctmp,*,iostat=io) env%ref%efield(1)
          if (io == 0) processedarg(i+1) = .true.
          ctmp = arg2
          read (ctmp,*,iostat=io) env%ref%efield(2)
          if (io == 0) processedarg(i+2) = .true.
          ctmp = arg3
          read (ctmp,*,iostat=io) env%ref%efield(3)
          if (io == 0) processedarg(i+3) = .true.
          write (stdout,'("  --efield: ",3(1x,es10.3)," V/Å")') env%ref%efield(1:3)
        else
          call parseflags_missing(argument)
        end if

      case ('-ceh_guess')
        processedarg(i) = .true.
        env%ceh_guess = .true.

      case ('-spin-polarized','-spinpol')
        processedarg(i) = .true.
        env%spin_polarized = .true.

      case ('-dscal','-dispscal','-dscal_global','-dispscal_global')
        processedarg(i) = .true.
        env%cts%dispscal_md = .true.
        if (index(argument,'_global') .ne. 0) then
          env%cts%dispscal_global = .true.
        end if
        if (nra .ge. i+1) then
          ctmp = arg1
          read (ctmp,*,iostat=io) rdum
          if (io .eq. 0) then
            env%cts%dscal = rdum
            processedarg(i+1) = .true.
          end if
        end if
        if (.not.env%legacy) call parseflags_deprecated(argument)

      case ('-mtd_kscal','-mtdkscal')
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          env%mtd_kscal = xx(1)
          processedarg(i+1) = .true.
        else
          call parseflags_missing(argument)
        end if

      case ('-norestart')
        processedarg(i) = .true.
        env%allowrestart = .false.

      case ('-readbias')
        processedarg(i) = .true.
        env%readbias = .true.

      case ('-useonly')
        processedarg(i) = .true.
        env%properties = p_useonly
        env%autozsort = .false.
        env%dummypercent = 1.0_wp
        if (nra .ge. i+1) then
          atmp = adjustl(arg(i+1))
          if (atmp(1:1) .ne. '-') then
            read (atmp,*) env%dummypercent
            processedarg(i+1) = .true.
          end if
        end if

      case ('-gbsa','-g','-alpb')   !> use GBSA or ALPB implicit solvation
        processedarg(i) = .true.
        env%gbsa = .true.
        atmp = adjustl(arg(i+1))
        if (atmp(1:1) .ne. '-'.and.atmp(1:1) .ne. ' ') then
          env%solvent = arg(i+1)
          processedarg(i+1) = .true.
          if (trim(argument) == '-alpb') then
            env%solv = '--alpb '//trim(env%solvent)
          else
            env%solv = '--gbsa '//trim(env%solvent)
          end if
        else
          call parseflags_missing(argument)
        end if
        write (stdout,'(2x,a,a)') trim(env%solv),' : implicit solvation'

      case ('-chrg')                  !> create a .CHRG file
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          processedarg(i+1) = .true.
          open (newunit=ich,file='.CHRG')
          env%chrg = nint(xx(1))
          env%ref%ichrg = env%chrg
          write (ich,'(i0)') nint(xx(1))
          close (ich)
          write (stdout,'(2x,a,1x,a)') trim(arg(i)),arg1
        else
          call parseflags_missing(argument)
        end if

      case ('-uhf')                    !> create a .UHF file
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          processedarg(i+1) = .true.
          open (newunit=ich,file='.uhf')
          env%uhf = nint(xx(1))
          env%ref%uhf = env%uhf
          write (ich,'(i0)') nint(xx(1))
          close (ich)
          write (stdout,'(2x,a,1x,a)') trim(arg(i)),arg1
        else
          call parseflags_missing(argument)
        end if

      case ('-mult','-multiplicity')   !> spin multiplicity (2S+1); stored as uhf = mult-1
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          processedarg(i+1) = .true.
          env%uhf = nint(xx(1))-1
          env%ref%uhf = env%uhf
          open (newunit=ich,file='.uhf')
          write (ich,'(i0)') env%uhf
          close (ich)
          write (stdout,'(2x,a,1x,a)') trim(arg(i)),arg1
        else
          call parseflags_missing(argument)
        end if

      case ('-grad')
        processedarg(i) = .true.
        env%gradsp = .true.

      case ('-nograd')
        processedarg(i) = .true.
        env%gradsp = .false.

      case ('-len','-mdlen','-mdtime') !> set md length in ps
        processedarg(i) = .true.
        atmp = arg(i+1)
        call to_lower(atmp)
        j = index(atmp,'x')
        env%user_mdtime = .true.
        if (j .ne. 0) then             !> scaling of the md length
          btmp = atmp(j+1:)
          env%scallen = .true.
          call readl(btmp,xx,j)
          if (j > 0) then
            env%mdlenfac = xx(1)
            processedarg(i+1) = .true.
          else
            call parseflags_missing(argument)
          end if
        else                           !> direct setting of the md length
          call readl(arg1,xx,j)
          if (j > 0) then
            processedarg(i+1) = .true.
            env%mdtime = xx(1)
            write (stdout,'(2x,a,1x,a,1x,a)') trim(arg(i)),arg1, &
            &    '(MD length in ps)'
          else
            call parseflags_missing(argument)
          end if
        end if

      case ('-shake')                           !> set shake
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          env%shake = nint(xx(1))
          processedarg(i+1) = .true.
        else
          call parseflags_missing(argument)
        end if

      case ('-tstep')                           !> set MD timestep in fs
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          env%mdstep = xx(1)
          env%user_mdstep = .true.
          processedarg(i+1) = .true.
        else
          call parseflags_missing(argument)
        end if

      case ('-vbdump')                          !> Vbias dump in ps
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          xx(2) = xx(1)*1000
          env%mddump = nint(xx(2))
          processedarg(i+1) = .true.
        else
          call parseflags_missing(argument)
        end if

      case ('-tnmd')                            !> temperature for additional normal MDs
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          env%nmdtemp = xx(1)
          processedarg(i+1) = .true.
        else
          call parseflags_missing(argument)
        end if

      case ('-mdskip')                          !> set skipping structures in -mdopt
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          env%mdskip = nint(xx(1))
          processedarg(i+1) = .true.
        else
          call parseflags_missing(argument)
        end if

      case ('-mddump')                          !> set dumpstep for writing structures from MD
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          env%user_dumxyz = .true.
          env%mddumpxyz = nint(xx(1))
          processedarg(i+1) = .true.
        else
          call parseflags_missing(argument)
        end if

      case ('-mdscal','-lenscal')       !> scale md length
        processedarg(i) = .true.
        env%scallen = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          env%mdlenfac = xx(1)
          processedarg(i+1) = .true.
        else
          call parseflags_missing(argument)
        end if

      case ('-mdtemp')                          !> set MTD temperature (V2 version)
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          env%mdtemp = xx(1)
          env%user_temp = .true.
          processedarg(i+1) = .true.
        else
          call parseflags_missing(argument)
        end if

      case ('-nmtd') !> set number of MTDs
        processedarg(i) = .true.
        env%runver = 787878
        call readl(arg1,xx,j)
        if (j > 0) then
          processedarg(i+1) = .true.
          env%nmetadyn = nint(xx(1))
        else
          call parseflags_missing(argument)
        end if

      case ('-gcmax','-setgcmax')       !> set maximum number of structures for GC
        processedarg(i) = .true.
        env%setgcmax = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          env%gcmax = xx(1)
          processedarg(i+1) = .true.
        else
          call parseflags_missing(argument)
        end if

      case ('-xnam')                    !> select a name for the xTB executeable
        processedarg(i) = .true.
        env%ProgName = arg1
        write (stdout,'(2x,''-xnam :'')',advance='no')
        write (stdout,'(1x,''xtb executable was set to: "'',a,''"'')') trim(env%ProgName)
        processedarg(i+1) = .true.

      case ('-niceprint')               !> progres bar printout
        processedarg(i) = .true.
        env%niceprint = .true.

      case ('-origin')                  !> track the origin (i.e. the generation step) of each conformer
        processedarg(i) = .true.
        env%trackorigin = .true.
        write (stdout,'(2x,a,1x,a)') trim(arg(i)),': tracking conformer origins.'

      case ('-constrain')               !> provide a list of atoms to write a .xcontrol.sample
        processedarg(i) = .true.
        ctmp = arg1
        call quick_constrain_file('coord',env%nat,env%ref%at,ctmp)
        processedarg(i+1) = .true.

      case ('-freeze')               !> freeze atoms during optimization/MD
        processedarg(i) = .true.
        if (i+1 .le. nra) then
          ctmp = arg1
          call get_atlist(env%ref%nat,atlist,ctmp,env%ref%at)
          env%calc%nfreeze = count(atlist)
          call move_alloc(atlist,env%calc%freezelist)
          processedarg(i+1) = .true.
          write (stdout,'(2x,a,1x,i0,1x,a)') '-freeze :',env%calc%nfreeze,'atoms frozen'
        else
          call parseflags_missing(argument)
        end if

      case ('-nocbonds')
        processedarg(i) = .true.
        bondconst = .false.
        env%cts%cbonds_global = .false.
        env%cts%cbonds_md = .false.
        inquire (file='bondlengths',exist=ex)
        if (ex) call remove('bondlengths')

      case ('-cbonds','-cbonds_md','-cbonds_ez')  !> constrain all bonds
        processedarg(i) = .true.
        ctmp = arg1
        if (ctmp(1:1) .ne. '-') then
          read (ctmp,*,iostat=io) rdum
          if (io .eq. 0) then
            env%forceconst = rdum
            processedarg(i+1) = .true.
          end if
        end if
        ctype = 1
        bondconst = .true.
        env%cts%cbonds_global = .true.
        if (index(argument,'_md') .ne. 0) then !> if the bond constraint shall be present only in the MDs/MTDs
          env%cts%cbonds_md = .true.
          env%cts%cbonds_global = .false.
        end if
        if (index(argument,'_ez') .ne. 0) then !> if the only E/Z shall be constrained
          ctype = 5
        end if

      case ('-cmetal','-cmetal_md')            !> constrain transition metal coordination sites
        processedarg(i) = .true.
        ctmp = arg1
        if (ctmp(1:1) .ne. '-') then
          read (ctmp,*,iostat=io) rdum
          if (io .eq. 0) then
            env%forceconst = rdum
            processedarg(i+1) = .true.
          end if
        end if
        ctype = 2
        bondconst = .true.
        env%cts%cbonds_global = .true.
        if (index(argument,'_md') .ne. 0) then
          env%cts%cbonds_md = .true.
          env%cts%cbonds_global = .false.
        end if

      case ('-cheavy','-fixheavy','-cheavy_md')  !> constrain all heavy atom bonds
        processedarg(i) = .true.
        ctmp = arg1
        if (ctmp(1:1) .ne. '-') then
          read (ctmp,*,iostat=io) rdum
          if (io .eq. 0) then
            env%forceconst = rdum
            processedarg(i+1) = .true.
          end if
        end if
        ctype = 3
        bondconst = .true.
        env%cts%cbonds_global = .true.
        if (index(argument,'_md') .ne. 0) then
          env%cts%cbonds_md = .true.
          env%cts%cbonds_global = .false.
        end if

      case ('-clight','-fixhyd','-clight_md')  !> constraint all X-H bonds
        processedarg(i) = .true.
        ctmp = arg1
        if (ctmp(1:1) .ne. '-') then
          read (ctmp,*,iostat=io) rdum
          if (io .eq. 0) then
            env%forceconst = rdum
            processedarg(i+1) = .true.
          end if
        end if
        ctype = 4
        bondconst = .true.
        env%cts%cbonds_global = .true.
        if (index(argument,'_md') .ne. 0) then
          env%cts%cbonds_md = .true.
          env%cts%cbonds_global = .false.
        end if

      case ('-cfile','-cinp','-C','-c')     !> specify the constrain file
        processedarg(i) = .true.
        ctmp = arg1
        if (ctmp(1:1) .ne. '-') then
          env%constraints = trim(ctmp)
          processedarg(i+1) = .true.
          write (stdout,'(2x,a,1x,a)') argument//' :',trim(ctmp)
        end if

      case ('-fc','-forceconstant')
        processedarg(i) = .true.
        ctmp = arg1
        if (i+1 .le. nra) then
          call readl(arg1,xx,j)
          if (j > 0) then
            env%forceconst = xx(1)
            processedarg(i+1) = .true.
          else
            call parseflags_missing(argument)
          end if
        end if
        write (stdout,'(2x,a,f6.4,a)') '-fc ',env%forceconst,': selected force constant in Eh'

      case ('-nomlo','-no-multilevel')   !> turn off multilevel optimization
        processedarg(i) = .true.
        env%multilevelopt = .false.

      case ('-normmd')  !> set number of normMDs
        processedarg(i) = .true.
        env%rotamermds = .true.
        if (i+1 .le. nra) then
          call readl(arg1,xx,j)
          if (j > 0) then
            env%nrotammds = nint(xx(1))  !> how many lowest conformers?
            processedarg(i+1) = .true.
          else
            call parseflags_missing(argument)
          end if
        end if
        if (i+2 .le. nra) then
          call readl(arg2,xx,j)
          if (j > 0) then
            env%temps = nint(xx(1))      !> how many different temperatures
            processedarg(i+2) = .true.
          end if
        end if

      case ('-rmsdpot','-gesc')
        processedarg(i) = .true.
        ctmp = arg1
        inquire (file=ctmp,exist=ex)
        if (ex) then
          env%cts%usermsdpot = .true.
          call getcwd(atmp)
          env%cts%rmsdpotfile = trim(atmp)//'/'//ctmp
          write (stdout,'(2x,a,a,a,a)') argument,': using <',ctmp,'> as bias'
          processedarg(i+1) = .true.
        else
          write (stdout,'(a,a)') argument,': Warning! File could not be found!'
        end if
        call parseflags_deprecated(argument)

      case ('-mergebias','-mergebias+','-gesc+')
        processedarg(i) = .true.
        env%properties = -9224
        if (index(argument,'+') > 0) env%properties = p_gesc2
        ctmp = arg1
        inquire (file=ctmp,exist=ex)
        if (ex) then
          env%biasfile = ctmp
          processedarg(i+1) = .true.
        end if
        env%autozsort = .false.
        call parseflags_deprecated(argument)

      case ('-gescopt')
        processedarg(i) = .true.
        env%gescoptlev = optlevnum(arg(i+1),iostat=io)
        if (io == 0) processedarg(i+1) = .true.

      case ('-gescheavy','-heavygesc','-gesc_heavy')
        processedarg(i) = .true.
        env%cts%gesc_heavy = .true.

      case ('-rthr2') !> bias rmsd threshold
        processedarg(i) = .true.
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0) then
          env%rthr2 = rdum
          processedarg(i+1) = .true.
        else
          call parseflags_missing(argument)
        end if
        call parseflags_deprecated(argument)

      case ('-kshift')
        processedarg(i) = .true.
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0) then
          env%kshift = rdum
          processedarg(i+1) = .true.
        else
          call parseflags_missing(argument)
        end if
        env%kshiftnum = 1
        call parseflags_deprecated(argument)

      case ('-hflip')
        processedarg(i) = .true.
        env%doOHflip = .true.

      case ('-noflip')
        processedarg(i) = .true.
        env%doOHflip = .false.

      case ('-maxflip')
        processedarg(i) = .true.
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0.and.(index(arg(i+1),'-') .eq. 0)) then
          env%maxflip = nint(rdum)
          processedarg(i+1) = .true.
        else
          call parseflags_missing(argument)
        end if

      case ('-osdf')
        processedarg(i) = .true.
        env%outputsdf = .true.
        write (stdout,'(2x,a," :",1x,a)') trim(arg(i)), &
        & "output ensemble requested in sdf format"

      case ('-wscal')                           !> scale size of wall potential
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          env%potscal = xx(1)
          env%user_wscal = .true.
          processedarg(i+1) = .true.
        else
          call parseflags_missing(argument)
        end if

      case ('-wpad')                            !> scale size of wall potential
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          env%potpad = xx(1)
          processedarg(i+1) = .true.
        else
          call parseflags_missing(argument)
        end if

      case ('-watoms','-wat')
        processedarg(i) = .true.
        ctmp = arg(i+1)
        if (ctmp(1:1) .ne. '-') then
          env%potatlist = trim(ctmp)
          write (stdout,*) env%potatlist
          processedarg(i+1) = .true.
        end if

      case ('-wall')
        processedarg(i) = .true.
        env%wallsetup = .true.
        write (stdout,'(2x,a,1x,a)') '--wall:','requesting setup of wall potential'

      case ('-wallxl','-wall-xl')
        processedarg(i) = .true.
        env%wallsetup = .true.
        env%potscal = 1.5_wp
        write (stdout,'(2x,a,1x,a)') '--wall-xl:','requesting setup of wall potential (x1.5 size)'

      case ('-wallxxl','-wall-xxl')
        processedarg(i) = .true.
        env%wallsetup = .true.
        env%potscal = 2.0_wp
        write (stdout,'(2x,a,1x,a)') '--wall-xxl:','requesting setup of wall potential (x2.0 size)'

      case ('-alkylize')
        processedarg(i) = .true.
        write (stdout,'(2x,a,1x)',advance='no') '--alkylize'
        env%alkylize = .true.
        if (nra >= i+1) then
          ctmp = arg1
          select case (ctmp)
          case ('full','sample')
            env%alkylizeskip = .false.
            write (stdout,'(a,1x)',advance='no') ctmp
            processedarg(i+1) = .true.
          case ('simple')
            write (stdout,'(a,1x)',advance='no') ctmp
            processedarg(i+1) = .true.
          end select
        end if
        write (stdout,'(a)') ': automatic alkyl group dispatch'

!========================================================================================!
!------ flags for parallelization / disk space
!========================================================================================!
      case ('-T','-P','-parallel')  !> set total number of OMP threads, this replaces -P and -O entirely
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          processedarg(i+1) = .true.
        else
          call parseflags_missing(argument)
        end if
        if (index(arg(i+1),'-') .ne. 0) xx = 0d0
        env%Threads = nint(xx(1))
        env%autothreads = .true.
        env%threadssetmanual = .true.
        write (stdout,'(2x,a,1x,i0,1x,a)') trim(arg(i)),nint(xx(1)), &
        &     '(CPUs/Threads selected)'

      case ('-omp-nested')  !> allow nested OpenMP threading (ON by default)
        processedarg(i) = .true.
        env%omp_allow_nested = .true.
        write (stdout,'(2x,a,1x,a)') trim(arg(i)),': nested OpenMP threading enabled'
      case ('-no-omp-nested','-omp-nonested')  !> forbid nested OpenMP threading
        processedarg(i) = .true.
        env%omp_allow_nested = .false.
        write (stdout,'(2x,a,1x,a)') trim(arg(i)),': nested OpenMP threading disabled'

      case ('-inplace')     !> activate "in-place" mode for optimizations (ON by default)
        processedarg(i) = .true.
        env%inplaceMode = .true.
      case ('-enso')             !> compare two given ensembles
        processedarg(i) = .true.
        env%ENSO = .true.

      case ('-compare')          !> compare two given ensembles
        processedarg(i) = .true.
        env%compareens = .true.

      case ('-maxcomp')          !> maximum number of lowest structures to compare with "-compare"
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          processedarg(i+1) = .true.
          env%maxcompare = nint(xx(1))
        else
          call parseflags_missing(argument)
        end if

      case ('-ewin')             !> set energy threshold in kcal/mol
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          processedarg(i+1) = .true.
          env%ewin = abs(xx(1))
          if (any((/p_protonate,p_deprotonate,p_tautomerize/) == env%properties)) then
            env%protb%ewin = abs(xx(1))
          end if
          write (stdout,'(2x,a,1x,a)') trim(arg(i)),arg1
        else
          call parseflags_missing(argument)
        end if

      case ('-rthr')             !> set RMSD thr
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          processedarg(i+1) = .true.
          env%rthr = xx(1)
          write (stdout,'(2x,a,1x,a)') trim(arg(i)),arg1
        else
          call parseflags_missing(argument)
        end if

      case ('-ethr')             !> set E thr
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          processedarg(i+1) = .true.
          env%ethr = xx(1)
          write (stdout,'(2x,a,1x,a)') trim(arg(i)),arg1
        else
          call parseflags_missing(argument)
        end if

      case ('-bthr')             !> set rot const thr
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          processedarg(i+1) = .true.
          env%thresholds(4) = xx(1)  !> legacy
          env%bthr2 = xx(1)
          write (stdout,'(2x,a,1x,a)') trim(arg(i)),arg1
        else
          call parseflags_missing(argument)
        end if

      case ('-allrot')           !> use all rotational constants for comparison, instead of mean
        processedarg(i) = .true.
        env%allrot = .true.

      case ('-athr')             !> set int. rotation. equal atoms for NMR thr
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          processedarg(i+1) = .true.
          env%athr = xx(1)
          write (stdout,'(2x,a,1x,a)') trim(arg(i)),arg1
        else
          call parseflags_missing(argument)
        end if

      case ('-pthr')             !> set population thr
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          processedarg(i+1) = .true.
          rdum = min(1.0_wp,xx(1)) !--> pthr <= 1
          rdum = max(0.0_wp,rdum)  !--> pthr >= 0
          env%pthr = rdum
          write (stdout,'(2x,a,1x,f6.4)') trim(arg(i)),rdum !arg1
        else
          call parseflags_missing(argument)
        end if

      case ('-eqv')
        processedarg(i) = .true.
        env%doNMR = .true. !> option for the very last confg call

      case ('-zsort')
        processedarg(i) = .true.
        env%onlyZsort = .true.                                 !perform only the zsort subroutine
        env%autozsort = .true.                                 ! CB: needs to be set to run zsort
        write (stdout,'(2x,a,1x,a)') trim(arg(i)),' : only using the ZSORT subroutine.'

      case ('-metac')                                        !automatic complete of mag. and chem. methyl equivalencies
        processedarg(i) = .true.
        env%methautocorr = .true.

      case ('-esort')          !> cregen legacy option
        processedarg(i) = .true.
        env%esort = .true.

      case ('-debug')
        processedarg(i) = .true.
        env%cgf(1) = .true.    !> debug option for confg

      case ('-nowr')
        processedarg(i) = .true.
        env%cgf(2) = .false.   !> newfile option for confg

      case ('-eqan')
        processedarg(i) = .true.
        env%cgf(3) = .true.    !> equivalence analysis on (for NMR)

      case ('-noeqan')
        processedarg(i) = .true.
        env%cgf(3) = .false.   !> equivalence analysis off (for nmr)

      case ('-rot')
        processedarg(i) = .true.
        env%cgf(5) = .false.   !> just rotamer check

      case ('-nmr')            !> NMR mode for confscript
        processedarg(i) = .true.
        env%doNMR = .true.
        env%optlev = 2.0d0

      case ('-fullcre')
        processedarg(i) = .true.
        env%doNMR = .true.
        env%fullcre = .true.

      case ('-heavy')
        processedarg(i) = .true.
        env%cgf(4) = .true.   !> perform just the heavy atom RMSD
        env%heavyrmsd = .true.

      case ('-temp')
        processedarg(i) = .true.
        ctmp = arg1
        if (index(ctmp,'-') .eq. 0) then
          call readl(arg1,xx,j)
          env%tboltz = xx(1)
          processedarg(i+1) = .true.
        end if

      case ('-prsc')                !> write scoord files
        processedarg(i) = .true.
        env%printscoords = .true.

      case ('-noprsc')              !> don't write scoord files
        processedarg(i) = .true.
        env%printscoords = .false.

      case ('-subrmsd')             !> use only the RMSD for atoms that are included in the MTD
        processedarg(i) = .true.
        env%subRMSD = .true.

      case ('-noopt')               !> skip the pre-optimization with GFNn-xTB before the confsearch
        processedarg(i) = .true.
        env%preopt = .false.

      case ('-topo','-topocheck')
        processedarg(i) = .true.
        env%checktopo = .true.

      case ('-notopo','-notopocheck','-noreftopo')
        processedarg(i) = .true.
        env%checktopo = .false.
        ctmp = arg1
        if (ctmp(1:1) .ne. '-') then
          processedarg(i+1) = .true.
          call parse_topo_excl(env,ctmp)
          if (allocated(env%excludeTOPO)) then
            env%checktopo = .true.
          end if
        end if
        env%reftopo = .false.

      case ('-ezcheck','-checkez')
        processedarg(i) = .true.
        env%checkiso = .true.

      case ('-noezcheck','-nocheckez')
        processedarg(i) = .true.
        env%checkiso = .false.

      case ('-inversion')
        processedarg(i) = .true.
        ctmp = lowercase(arg1)
        select case (ctmp)
        case ('auto')
          env%iinversion = 0
        case ('on')
          env%iinversion = 1
        case ('off')
          env%iinversion = 2
        case default
          write (stdout,'(a,a,a,a)') 'invalid argument for ',argument,': ',trim(ctmp)
          stop
        end select
        processedarg(i+1) = .true.
!========================================================================================!
!-------- PROPERTY CALCULATION related flags
!========================================================================================!
      case ('-protonate')             !> protonation tool
        processedarg(i) = .true.
        env%properties = p_protonate
        env%autozsort = .false.
        env%protb%threshsort = .true.
        ctmp = arg1
        if (ctmp(1:1) .ne. '-') then
          processedarg(i+1) = .true.
          read (ctmp,*,iostat=io) idum
          if (io .eq. 0) env%protb%amount = idum
        end if

      case ('-swel')                  !> switch out H+ to something else in protonation script
        processedarg(i) = .true.
        if (env%properties .eq. -3) then
          call swparse(arg(i+1),env%protb)
          processedarg(i+1) = .true.
        end if

      case ('-deprotonate')           !> deprotonation tool
        processedarg(i) = .true.
        env%properties = p_deprotonate
        env%autozsort = .false.
        env%protb%threshsort = .true.
        ctmp = arg1
        if (ctmp(1:1) .ne. '-') then
          processedarg(i+1) = .true.
          read (ctmp,*,iostat=io) idum
          if (io .eq. 0) env%protb%amount = idum
        end if

      case ('-tautomerize')           !> tautomerization tool
        processedarg(i) = .true.
        env%properties = p_tautomerize
        env%autozsort = .false.
        env%protb%threshsort = .true.

      case ('-tautomerize2','-exttautomerize')
        processedarg(i) = .true.
        if (env%properties == p_propcalc) then
          env%properties = p_tautomerize2
        else
          call env%addjob(abs(p_tautomerize2))
        end if
        env%autozsort = .false.
        env%protb%threshsort = .true.
        env%runver = 33
        env%relax = .true.
        env%performCross = .false.  !> skip the genetic crossing
        env%trackorigin = .false.
        env%Maxrestart = 1

      case ('-relax')
        processedarg(i) = .true.
        env%runver = 33
        env%relax = .true.
        env%performCross = .false.  !> skip the genetic crossing
        env%trackorigin = .false.
        env%Maxrestart = 1

      case ('-trev','-tdp')
        processedarg(i) = .true.
        env%protb%deprotprot = .true. !> switch to deprotonation-first mode in tautomerization

      case ('-iter')                !> number of Protonation/Deprotonation cycles in Tautomerization
        processedarg(i) = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          processedarg(i+1) = .true.
          env%protb%iter = nint(xx(1))
        end if

      case ('-texcl','-blacklist')
        processedarg(i) = .true.
        ctmp = arg1
        processedarg(i+1) = .true.

      case ('-strict')
        processedarg(i) = .true.
        env%protb%strictPDT = .true.

      case ('-verystrict','-vstrict')
        processedarg(i) = .true.
        env%protb%strictPDT = .false.
        env%protb%fixPDT = .true.

      case ('-fstrict')
        processedarg(i) = .true.
        env%protb%strictPDT = .true.
        env%protb%fixPDT = .true.

      case ('-corr','-abcorr')
        processedarg(i) = .true.
        env%protb%strictPDT = .true.
        env%protb%fixPDT = .true.
        env%protb%ABcorrection = .true.

      case ('-pkaensemble')
        processedarg(i) = .true.
        env%preopt = .false.
        env%presp = .false.
        call pka_argparse2(env,arg(i+1),arg(i+2),env%protb%pka_mode)
        processedarg(i+1) = .true.
        processedarg(i+2) = .true.

      case ('-pkaparam')
        processedarg(i) = .true.
        env%protb%rdcfer = .true.
        if (i+1 .le. nra) then
          ctmp = arg1
          if (ctmp(1:1) .ne. '-') then
            processedarg(i+1) = .true.
            env%protb%cferfile = ctmp
          end if
        end if

!========================================================================================!
!--------- ENTROPY related settings
!========================================================================================!
      case ('-entropy','-entropic')  !> new, specialized calculation of molecular entropies
        processedarg(i) = .true.
        write (stdout,'(2x,a,'' : enhanced ensemble entropy calculation'')') argument
        if (env%properties == p_propcalc) then
          !>--- for standalone use
          env%properties = p_CREentropy

        elseif (env%confgo.and.env%properties == -1) then
          !>--- as extension for CREGEN
          env%entropic = .true.
          env%fullcre = .true.

        else if (env%crestver == crest_imtd) then
          !>--- works as an extensiton to the conformational search
          env%properties = abs(p_CREentropy)
          env%autozsort = .false.     !> turn off zsort (since we are not going to GC anyways)
          env%performCross = .false.  !> turn off GC
          env%entropic = .true.       !> indicator for this runtype
          env%Maxrestart = 1          !> turn off MTD iterations (just do one)
          env%rotamermds = .false.    !> turn off normMDs
          env%entropymd = .true.      !> special static MTDs
          call read_bhess_ref(env,'coord')
        end if
        env%runver = 111             !> version  for selection of MTD bias settings
        env%doNMR = .true.           !> we need equivalencies
        if (i+1 .le. nra) then
          ctmp = arg1    !> second argument can be the temperature
          if (index(ctmp,'-') .eq. 0) then
            processedarg(i+1) = .true.
            call readl(arg1,xx,j)
            env%tboltz = xx(1)
          end if
        end if
        call env%addjob(env%properties)

      case ('-scthr','-entropy_cthr')
        processedarg(i) = .true.
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0) then
          processedarg(i+1) = .true.
          env%emtd%confthr = rdum
        else
          call parseflags_missing(argument)
        end if

      case ('-ssthr','-entropy_sthr')
        processedarg(i) = .true.
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0) then
          processedarg(i+1) = .true.
          env%emtd%sconvthr = rdum
        else
          call parseflags_missing(argument)
        end if

      case ('-rrhoav')             !> see above in the first specification of -rrhoav
        processedarg(i) = .true.
        env%properties = p_rrhoaverage
        call read_bhess_ref(env,'coord')

      case ('-avbhess')
        processedarg(i) = .true.
        env%thermo%avbhess = .true. !> use bhess in rrhoav for all structures (expensive)

      case ('-avchess')
        processedarg(i) = .true.
        env%thermo%constrhess = .true.   !> apply constraints during rrhoav routine

      case ('-printpop')
        processedarg(i) = .true.
        env%thermo%printpop = .true. !> print a file with free energy pop. at different T

      case ('-noref')              !> dont use a bhess reference
        processedarg(i) = .true.
        env%emtd%bhess = .false.

      case ('-ref')
        processedarg(i) = .true.
        env%emtd%bhess = .true.
        inquire (file=arg1,exist=ex)
        if (ex) then
          processedarg(i+1) = .true.
          call read_bhess_ref(env,arg1)
        end if

      case ('-pcap')
        processedarg(i) = .true.
        read (arg(i+1),*,iostat=io) j
        if (io == 0.and.(index(arg(i+1),'-') .eq. 0)) then
          processedarg(i+1) = .true.
          env%thermo%pcap = j
        end if

      case ('-ptot')
        processedarg(i) = .true.
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0.and.(index(arg(i+1),'-') .eq. 0)) then
          processedarg(i+1) = .true.
          if (rdum > 1.0d0) rdum = 1.0d0
          env%thermo%ptot = rdum
        end if

      case ('-ithr')
        processedarg(i) = .true.
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0) then
          processedarg(i+1) = .true.
          if (rdum > 0.0d0) rdum = 0.0
          env%thermo%ithr = rdum
        end if

      case ('-rotorcut','-sthr')
        processedarg(i) = .true.
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0.and.(index(arg(i+1),'-') .eq. 0)) then
          processedarg(i+1) = .true.
          if (rdum < 0.0d0) rdum = 0.0d0
          env%thermo%sthr = rdum
        end if

      case ('-fscal')
        processedarg(i) = .true.
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0.and.(index(arg(i+1),'-') .eq. 0)) then
          processedarg(i+1) = .true.
          env%thermo%fscal = rdum
        end if

      case ('-trange')    !> provide a range of temperatures (min max step) for entropy evaluation
        processedarg(i) = .true.
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0.and.(index(arg(i+1),'-') .eq. 0)) then
          processedarg(i+1) = .true.
          env%thermo%trange(1) = rdum  !> T start
        end if
        read (arg(i+2),*,iostat=io) rdum
        if (io == 0.and.(index(arg(i+2),'-') .eq. 0)) then
          processedarg(i+2) = .true.
          env%thermo%trange(2) = rdum  !> T stop (approx.)
        end if
        read (arg(i+3),*,iostat=io) rdum
        if (io == 0.and.(index(arg(i+3),'-') .eq. 0)) then
          processedarg(i+3) = .true.
          env%thermo%trange(3) = rdum  !> T step
        end if

      case ('-tread')   !> read a file with temperatures (one per line) for entropy evaluation
        processedarg(i) = .true.
        ctmp = arg1
        inquire (file=ctmp,exist=ex)
        if (ex) then
          processedarg(i+1) = .true.
          call env%thermo%read_temps(ctmp)
        end if
!=========================================================================================!
!-------- QCG-Related flags
!=========================================================================================!
      case ('-nopreopt')
        processedarg(i) = .true.
        env%nopreopt = .true.
        env%qcg_flag = .true.

      case ('-xtbiff')
        processedarg(i) = .true.
        env%use_xtbiff = .true.

      case ('-grow')
        processedarg(i) = .true.
        env%qcg_runtype = 0
        env%qcg_flag = .true.

      case ('-ensemble')
        processedarg(i) = .true.
        env%qcg_runtype = 1
        env%qcg_flag = .true.

      case ('-esolv')
        processedarg(i) = .true.
        env%qcg_runtype = 2
        env%qcg_flag = .true.

      case ('-gsolv')
        processedarg(i) = .true.
        env%qcg_runtype = 3
        env%qcg_flag = .true.

      case ('-nsolv')
        processedarg(i) = .true.
        env%qcg_flag = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          processedarg(i+1) = .true.
          env%nsolv = NINT(xx(1))
        else
          call parseflags_missing(argument)
        end if

      case ('-maxsolv')
        processedarg(i) = .true.
        env%qcg_flag = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          processedarg(i+1) = .true.
          env%max_solv = NINT(xx(1))
        else
          call parseflags_missing(argument)
        end if

      case ('-normdock')
        processedarg(i) = .true.
        env%docking_qcg_flag = ''

      case ('-fin_opt_gfn2')
        processedarg(i) = .true.
        env%final_gfn2_opt = .true.
        user_finopt = .true.

      case ('-no_fin_opt_gfn2')
        processedarg(i) = .true.
        env%final_gfn2_opt = .false.

      case ('-directed') !> specify the directed list
        processedarg(i) = .true.
        env%qcg_flag = .true.
        ctmp = arg1
        if (ctmp(1:1) .ne. '-') then
          processedarg(i+1) = .true.
          env%directed_file = trim(ctmp)
          write (stdout,'(2x,a,1x,a)') trim(argument)//' :',trim(ctmp)
        end if

      case ('-nclus')
        processedarg(i) = .true.
        env%qcg_flag = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          processedarg(i+1) = .true.
          env%nqcgclust = NINT(xx(1))
          env%user_nclust = .true.
        else
          call parseflags_missing(argument)
        end if

      case ('-freqscal')
        processedarg(i) = .true.
        env%qcg_flag = .true.
        call readl(arg1,xx,j)
        if (j > 0) then
          processedarg(i+1) = .true.
          env%freq_scal = (xx(1))
        else
          call parseflags_missing(argument)
        end if

      case ('-qcgmtd')
        processedarg(i) = .true.
        env%ensemble_method = -1
        env%qcg_flag = .true.

      case ('-ncimtd')
        processedarg(i) = .true.
        env%ensemble_method = 0
        env%qcg_flag = .true.

      case ('-md')
        processedarg(i) = .true.
        env%ensemble_method = 1
        env%qcg_flag = .true.
        if (.not.env%user_enslvl) then
          env%ensemble_opt = '--gfn2'
        end if

      case ('-mtd')
        processedarg(i) = .true.
        env%ensemble_method = 2
        env%qcg_flag = .true.
        if (.not.env%user_enslvl) then
          env%ensemble_opt = '--gfn2'
        end if

      case ('-samerand')
        processedarg(i) = .true.
        env%sameRandomNumber = .true.
        env%qcg_flag = .true.

      case ('-nocff')
        processedarg(i) = .true.
        env%cff = .false.
        env%qcg_flag = .true.

      case ('-enslvl')
        processedarg(i) = .true.
        ctmp = arg(i+1)
        processedarg(i+1) = .true.
        env%user_enslvl = .true.
        env%qcg_flag = .true.
        if (arg(i+1) == '-gfn') then
          processedarg(i+2) = .true.
          dtmp = arg2
          ctmp = trim(ctmp)//dtmp
        end if
        select case (ctmp)
        case ('gfn1')
          env%ensemble_opt = '--gfn1'
          write (stdout,'(2x, a)') 'Use of GFN1-xTB for ensemble search requested.'
        case ('gfn2')
          env%ensemble_opt = '--gfn2'
          write (stdout,'(2x, a)') 'Use of GFN2-xTB for ensemble search requested.'
        case ('gfn0')
          env%ensemble_opt = '--gfn0'
          write (stdout,'(2x, a)') 'Use of GFN0-xTB for ensemble search requested.'
        case ('gff','gfnff')
          env%ensemble_opt = '--gff'
          write (stdout,'(2x, a)') 'Use of GFN-FF for ensemble search requested.'
        end select

      case ('-freqlvl')
        processedarg(i) = .true.
        ctmp = arg(i+1)
        processedarg(i+1) = .true.
        env%qcg_flag = .true.
        if (arg(i+1) == '-gfn') then
          processedarg(i+2) = .true.
          dtmp = arg2
          ctmp = trim(ctmp)//dtmp
        end if
        select case (ctmp)
        case ('gfn1')
          env%freqver = '--gfn1'
          write (stdout,'(2x, a)') 'Use of GFN1-xTB for frequency computation requested.'
        case ('gfn2')
          env%freqver = '--gfn2'
          write (stdout,'(2x, a)') 'Use of GFN2-xTB for frequency computation requested.'
        case ('gfn0')
          env%freqver = '--gfn0'
          write (stdout,'(2x, a)') 'Use of GFN0-xTB for frequency computation requested.'
        case ('gff','gfnff')
          env%freqver = '--gff'
          write (stdout,'(2x, a)') 'Use of GFN-FF for frequency computation requested.'
        end select
!========================================================================================!
!-------- PRINCIPAL COMPONENT analysis and CLUSTERING flags
!========================================================================================!
      case ('-cluster')
        processedarg(i) = .true.
        write (stdout,'(2x,a,'' : ensemble clustering'')') trim(arg(i))
        if (env%properties == p_propcalc) then
          !>--- for standalone use
          env%properties = p_cluster
        elseif (env%crestver == crest_sorting) then
          !>--- standalone clustering via --sort --cluster
          env%sortmode = 'cluster'
        else if (any((/crest_imtd,crest_imtd2/) == env%crestver)) then
          !>--- works as an extensiton to the conformational search
          env%properties = abs(p_cluster)
        elseif (env%QCG) then
          env%properties = abs(p_cluster)
        end if
        env%doNMR = .true.     !> we need equivalencies
        call env%addjob(env%properties)
        if (i+1 .le. nra) then !second argument a distinct number of clusters
          read (arg(i+1),*,iostat=io) j
          if (io == 0.and.(index(arg(i+1),'-') .eq. 0)) then
            processedarg(i+1) = .true.
            env%nclust = j
          else
            env%nclust = 0
            if ((index(arg(i+1),'-') .eq. 0)) then
              processedarg(i+1) = .true.
              ctmp = arg1
              select case (ctmp)
              case ('loose')
                env%clustlev = -1
                write (stdout,'(2x,a,'' loose : using loose clustering setting'')') trim(arg(i))
              case ('normal')
                env%clustlev = 0
                write (stdout,'(2x,a,'' normal : using normal clustering setting'')') trim(arg(i))
              case ('tight')
                env%clustlev = 1
                write (stdout,'(2x,a,'' tight : using tight clustering setting'')') trim(arg(i))
              case ('vtight','verytight')
                env%clustlev = 2
                write (stdout,'(2x,a,'' vtight : using very tight clustering setting'')') trim(arg(i))
              case ('incremental','incr')
                env%clustlev = 10
                write (stdout,'(2x,a,'' incremental : using incremental clustering settings'')') trim(arg(i))
              case ('tightincremental','tightincr')
                env%clustlev = 11
                write (stdout,'(2x,a,'' tightincremental : using incremental clustering settings'')') &
                &    trim(arg(i))
              case ('vtightincremental','vtightincr')
                env%clustlev = 12
                write (stdout,'(2x,a,'' vtightincremental : using incremental clustering settings'')') &
                 & trim(arg(i))
              end select
            end if
          end if
        end if

      case ('-pccap')
        processedarg(i) = .true.
        if (i+1 .le. nra) then !> second argument is the max. number of PCs
          read (arg(i+1),*,iostat=io) j
          if (io == 0.and.(index(arg(i+1),'-') .eq. 0)) then
            processedarg(i+1) = .true.
            env%pccap = j
          end if
        end if

      case ('-nopcmin')
        processedarg(i) = .true.
        env%pcmin = 0.0d0

      case ('-pctype','-pctyp')
        processedarg(i) = .true.
        if (i+1 .le. nra) then
          ctmp = arg1
          if (ctmp(1:1) .ne. '-') then
            processedarg(i+1) = .true.
            env%pcmeasure = ctmp
          end if
        end if

      case ('-pcaex','-pcaexclude')
        processedarg(i) = .true.
        if (i+1 .le. nra) then
          ctmp = arg1
          if (ctmp(1:1) .ne. '-') then
            processedarg(i+1) = .true.
            env%atlist = ctmp
            env%pcaexclude = .true.
          end if
        end if
!========================================================================================!
!---------- PROPERTY MODE
!========================================================================================!
      case ('-prop')
        processedarg(i) = .true.
!>----------------------------------------------------------------
!> NOTE: These flags are outdated and using them is discouraged!
!>----------------------------------------------------------------
        if ((env%properties == p_none.or.    &
        &  env%properties == p_propcalc)) then         !property selection
          ctmp = arg1
          processedarg(i+1) = .true.
          PROPARG:select case(ctmp)
          case ('hess')                  !hessian calculation to free energies for all conformers
          env%properties2 = p_prop_hess
          case ('ohess')                 !optimization+hessian calculation
          env%properties2 = p_prop_ohess
          case ('autoir','autoIR')       !automated IR averaging for populated (-pthr) conformers
          env%properties2 = p_prop_autoir
          case ('reopt')                 !reoptimize only conformers at vtight level
          env%properties2 = p_prop_reopt
          case ('singlepoint','sp')      !singlepoint calculation and ensemble sorting
          env%properties2 = p_prop_rerank
          env%pclean = .true.
          case ('dipole')                !singlepoint calculation and dipole grepping
          env%properties2 = p_prop_dipole
          env%pclean = .true.
          case default
          env%properties2 = 0
          end select PROPARG
          if (env%properties2 .ne. 0) then
            call env%addjob(env%properties2)
          end if
        end if
        call parseflags_deprecated(argument)

      case ('-dftrc')                            !provide dft-rc file (including path)
        processedarg(i) = .true.
        atmp = ''
        if (nra .ge. (i+1)) atmp = adjustl(arg(i+1))
        if ((atmp(1:1) /= '-').and.(len_trim(atmp) .ge. 1)) then
          processedarg(i+1) = .true.
          env%dftrcfile = trim(atmp)
        end if
        call parseflags_deprecated(argument)

      case ('-hardcut')                          !cut DFT populations hard
        processedarg(i) = .true.
        env%hardcutDFT = .true.
        call parseflags_deprecated(argument)

      case ('-pclean')                           !cleanup option for property mode, i.e., remove PROP/
        processedarg(i) = .true.
        env%pclean = .true.
!========================================================================================!
      case ('-scratch')
        processedarg(i) = .true.
        !use a scratch directory to perform the calculation in
        env%scratch = .true.
        atmp = ''
        if (nra .ge. (i+1)) atmp = adjustl(arg(i+1))
        if ((atmp(1:1) /= '-').and.(len_trim(atmp) .ge. 1)) then
          processedarg(i+1) = .true.
          env%scratchdir = trim(atmp)
        end if

      case ('-keepscratch')
        processedarg(i) = .true.
        env%keepScratch = .true.
      end select !> ARGPARSER1
!========================================================================================!
    end if
  end do
!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> END OF ARGPARSER LOOP
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
  deallocate (strings,floats,xx)

  call parseflags_cli_summary(nra,arg,processedarg)

!>----- additional checks and settings
  if (env%crestver .eq. crest_solv) bondconst = .false.

  if (env%qcg_flag.and.env%crestver .ne. crest_solv) then
    error stop 'At least one flag is only usable for QCG runtype. Exit.'
  end if

  if (env%autozsort.and.env%crestver .eq. crest_solv) then
    error stop 'Z sorting of the input is unavailable for -qcg runtyp.'
  end if

!>--- avoid 0 potscal
  if (env%potscal < 1.0d-5) env%potscal = 1.0_wp

!>--- automatic wall potential for the LEGACY version
  if ((env%NCI.or.env%wallsetup).and.env%legacy) then
    call wallpot(env)
    if (env%wallsetup) then
      write (stdout,'(2x,a)') 'Automatically generated ellipsoide potential:'
    else
      write (stdout,'(2x,a)') 'Automatically generated ellipsoide potential for NCI mode:'
    end if
    call write_cts_NCI_pr(stdout,env%cts)
    write (stdout,*)
  end if

!>--- automatic bond constraint setup
  if (bondconst) then
    select case (ctype)
    case (1)
      call autoBondConstraint('coord',env%forceconst,env%wbofile)
    case (2)
      call autoMetalConstraint('coord',env%forceconst,env%wbofile)
    case (3)
      call autoHeavyConstraint('coord',env%forceconst)
    case (4)
      call autoHydrogenConstraint('coord',env%forceconst)
    case (5)
      call autoBondConstraint_withEZ('coord',env%forceconst,env%wbofile)
    end select
  end if

!>--- additional parsing of $setblock, .constrains and .confscriptrc file
  call parseRC2(env,bondconst)

!>--- internal constraint check-up
  call internal_constraint_repair(env,bondconst)

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> settings after user input parsing
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
  if ((any((/crest_imtd,crest_imtd2,crest_pka,crest_compr,11/) == env%crestver)).and.  &
  &  .not.env%confgo) then
    call defaultGF(env) !set Guiding Force default if none was read
  end if

  !-- increase gbsa grid for GFNn-xTB calculations (not for the FF)
  if ((env%gfnver .ne. '--gff').and.(env%gbsa)) then
    env%cts%ggrid = .true.
    env%cts%gbsagrid = 'tight'
  end if

  if (env%gfnver == '--gff') then
    env%hmass = 5.0d0
    if (bondconst) then
      env%autozsort = .false.
    end if
  end if

  !>--- defaults for QCG gfnff ensemble search
  if (env%crestver == crest_solv) then
    if (env%ensemble_opt .EQ. '--gff') then
      env%hmass = 5.0d0
      ctype = 5 !bond constraint
      bondconst = .true.
      env%cts%cbonds_md = .true.
      env%checkiso = .true.
      env%lmover = '--gfn2'
    end if
    if ((env%gfnver .EQ. '--gff').OR.(env%gfnver .EQ. '--gfn0')) then
      env%lmover = '--gfn2'
    else
      env%lmover = env%gfnver
    end if
  end if
  !> the extra GFN2 re-optimization is redundant only if the ensemble level is GFN2 already
  if (env%ensemble_opt == '--gfn2'.and.env%final_gfn2_opt) then
    if (env%crestver == crest_solv.and.user_finopt) write (stdout,'(1x,a)') &
      & 'Note: -fin_opt_gfn2 is ignored because the ensemble level is already GFN2-xTB.'
    env%final_gfn2_opt = .false.
  end if

  if (env%useqmdff) then
    env%autozsort = .false.
  end if

  if (.not.env%preopt.and.env%crestver .ne. crest_trialopt) then
    if (allocated(env%ref%topo)) deallocate (env%ref%topo)
  end if

!>-- turn off niceprint if we are not writing to terminal
  if (env%niceprint) then
    env%niceprint = myisatty(output_unit)
  end if

!>-- driver for optimization along trajectory, additional settings
  if (.not.any((/crest_mfmdgc,crest_imtd,crest_imtd2,crest_compr/) == env%crestver) &
      & .OR.(env%qcg_runtype .GT. 0.and.env%ensemble_method .EQ. 0)) then
    env%autozsort = .false.
    env%trackorigin = .false.
    env%confgo = .false.
  end if

!>-- final settings for property mode (-prop)
  if (env%properties .eq. p_none) then
    env%properties = env%properties2
  end if
  if (env%properties .eq. p_propcalc) then
    env%autozsort = .false.
  end if

!>-- some more zsort checks
  if (.not.env%onlyZsort.and.env%autozsort) then
    call zsortwarning2(env) !turn autozsort off when a .constrains file is present.
  end if
  if (env%autozsort) then
    if (allocated(env%ref%topo)) deallocate (env%ref%topo)
  end if
  if (env%sdfformat) then
    env%autozsort = .false.
  end if

!>--- 2023/08/19 moved zsort to a standalone property tool
  if (env%autozsort) then
    env%properties = p_zsort
  end if

!>--- convert ProgName to absolute path (to make legacy routines more stable)
  ctmp = absolute_filepath(trim(env%ProgName))
  env%ProgName = ctmp

!>--- for legacy runtypes, check if xtb is present
  if (env%legacy.or.env%QCG) then
    call checkprog_silent(env%ProgName,.true.,iostat=io)
    if (io /= 0) error stop
    write (stdout,'(/,a,a)') 'Selected path to xtb binary: ',trim(env%Progname)
  end if

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> FALLBACK setup of new calculator (important CREST >3.0 handling things)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
  if (.not.env%legacy.and.env%calc%ncalculations == 0) then
    write (stdout,'(/,a)',advance='no') '> Setting up backup calculator ...'
    flush (stdout)
    call env2calc_setup(env)
    write (stdout,*) 'done.'
  end if

!>--- pass on other settings (from cli) to new calculator
  if (.not.env%legacy) then
    call env2calc_modify(env)
  end if

!>--- important printouts
  if (.not.env%legacy) then

    if (env%crestver .ne. crest_sorting) then
      call env%calc%info(stdout)
    end if

    call print_frozen(env)
  end if

  return
end subroutine parseflags

subroutine parseflags_missing(arg)
  use crest_parameters
  implicit none
  character(len=*),intent(in) :: arg
  write (stdout,'(a)') trim(arg)//' requires a valid argument'
end subroutine parseflags_missing

subroutine parseflags_deprecated(arg)
  !**********************************************
  !* Print a deprecation error and stop.
  !**********************************************
  use crest_parameters
  use crest_data
  implicit none
  character(len=*),intent(in) :: arg
  write (stdout,'(/,a)') repeat('!',60)
  write (stdout,'(a)') '  DEPRECATED FLAG: '//trim(arg)
  write (stdout,'(a)') '  This flag has been removed. Please update your input.'
  write (stdout,'(a,/)') repeat('!',60)
  call creststop(status_safety)
end subroutine parseflags_deprecated

subroutine parseflags_cli_summary(nra,args,processedarg)
  use crest_parameters
  use crest_data
  implicit none
  integer,intent(in) :: nra
  character(len=*),intent(in) :: args(nra)
  logical,intent(in) :: processedarg(nra)
  integer :: ii,nprocessed
  nprocessed = count(processedarg,1)
  if (nprocessed == nra) then
    write (stdout,'(/,a)') '> All CLI arguments successfully processed.'
  else
    write (stdout,'(70("-"))')
    write (stdout,'(a,/)') '** WARNING ** Some CLI arguments were not correctly processed:'
    do ii = 1,nra
      if (processedarg(ii)) cycle
      write (stdout,'(1x,a,i4,a,t20,a)') 'Argument',ii,': ',trim(args(ii))
    end do

    write (stdout,'(/,a)') 'Please check your command line input for sanity.'
    write (stdout,'(70("-"))')
    write (stdout,*)
    call creststop(status_safety)
  end if

end subroutine parseflags_cli_summary

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!=========================================================================================!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!=========================================================================================!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!

subroutine parseRC2(env,bondconst)
!*****************************************
!* subroutine parseRC2
!* Parse the "confscriptrc" and set-block
!* in the input file
!* NOTE: this routine is very old, chances
!* are that some things are obsolete now
!******************************************
  use iso_fortran_env,wp => real64
  use crest_data
  use iomod
  use utilities
  use parse_xtbinput
  implicit none

  type(systemdata),intent(inout) :: env

  integer :: i,j,k
  character(len=512),allocatable :: cfiles(:)
  character(len=256) :: atmp,btmp
  character(len=512) :: dg,argument
  integer,allocatable :: atlist(:)
  logical :: ex,ex1,ex2
  logical :: create,atomlistused
  logical :: bondconst

!>--- check for constraint file
  ex1 = .false.
  inquire (file=env%constraints,exist=ex1)

!>--- do we have a user-set constraint to all bonds?
  if (bondconst) then
    inquire (file='bondlengths',exist=ex2)
    if (ex2) then
      call rd_cbonds('bondlengths',env)
      if (.not.env%cts%cbonds_md) then
        env%cts%cbonds_global = .true.
      end if
    end if
  end if

  if (ex1) then
    write (stdout,'(/,1x,a,a,a)') '<',trim(env%constraints),'> file present.'
    env%cts%used = .true.
  else
    env%cts%used = .false.
    return
  end if

!>--- read the data
  call read_constrainbuffer(env%constraints,env%cts)
  call sort_constraints(env%cts)
  write (stdout,*) 'content of the constraining file (sorted):'
  if (env%cts%ndim .gt. 20) then
    write (stdout,'(1x,a)') '<skipped due to length of constraining file>'
  else
    do i = 1,env%cts%ndim
      if (trim(env%cts%sett(i)) .ne. '') then
        write (stdout,'(''>'',1x,a)') trim(env%cts%sett(i))
      end if
    end do
  end if
  if (.not.env%legacy) then
    call parse_xtbinputfile(env,env%constraints)
  end if

!>--- some settings
  create = .false.
  atomlistused = .false.
  allocate (atlist(env%nat))

!>--- parse for special arguments that are used by CREST also
  do i = 1,env%cts%ndim
    btmp = env%cts%sett(i)
    if (trim(btmp) .eq. '') cycle
    atmp = btmp
    call to_lower(atmp)  !convert to lower-case for case-insensitivity
    if (index(atmp,'atomlist+') .ne. 0) then
      create = .true.
      atomlistused = .true.
      dg = atmp
      call split_set_args(dg,argument)
      call parse_atlist_new(trim(argument),env%rednat,env%nat,env%ref%at,atlist)
      write (stdout,'(2x,a)') trim(adjustl(btmp))
      write (stdout,'(5x,a,i0)') '# of atoms considered for RMSDs:',env%rednat
      env%includeRMSD = atlist !includeRMSD contains only the atoms that are included in RMSD
    end if
    if (index(atmp,'atomlist-') .ne. 0) then
      create = .true.
      atomlistused = .true.
      dg = atmp
      call split_set_args(dg,argument)
      call parse_atlist_new(trim(argument),j,env%nat,env%ref%at,atlist)
      env%rednat = env%nat-j
      write (stdout,'(2x,a)') trim(adjustl(btmp))
      write (stdout,'(3x,a,i0)') '# of atoms considered for RMSDs:',env%rednat
      env%includeRMSD = atlist !includeRMSD contains the atoms that are NOT included in RMSD
      do k = 1,env%nat
        if (env%includeRMSD(k) .lt. 1) then   !therefore the values have to be "inverted"
          env%includeRMSD(k) = 1
        else
          env%includeRMSD(k) = 0
        end if
      end do
    end if
    if ((index(atmp,'$metadyn') .ne. 0)) then
      do j = i+1,env%cts%ndim
        btmp = env%cts%sett(j)
        if (index(btmp,'$') .ne. 0) exit    !--- exit $metadyn-block
        if (index(btmp,'atoms:') .ne. 0) then
          create = .true.
          atomlistused = .true.
          dg = btmp
          call split_set_args(dg,argument)
          call parse_atlist_new(trim(argument),env%rednat,env%nat,env%ref%at,atlist)
          write (stdout,'(2x,a)') trim(adjustl(btmp))
          write (stdout,'(5x,a,i0)') '# of atoms considered for RMSDs:',env%rednat
          env%includeRMSD = atlist !includeRMSD contains only the atoms that are included in RMSD
        end if
      end do
    end if
    if (index(btmp,'reference=') .ne. 0) then
      call rdarg(btmp,'reference=',env%fixfile)
      write (stdout,'(1x,a,1x,a)') 'fix file:',trim(env%fixfile)
    end if
    if ((index(atmp,'$wall') .ne. 0)) then
      if (env%NCI) then
        env%cts%sett(i) = ''
        env%cts%pots = ''
        write (env%cts%pots(1),'(a)') '$wall'
        k = 2
        do j = i+1,env%cts%ndim
          btmp = env%cts%sett(j)
          if (index(btmp,'$') .ne. 0) exit    !--- exit $wall-block
          env%cts%sett(j) = ''
          write (env%cts%pots(k),'(a)') trim(btmp)
          k = k+1
        end do

        write (stdout,'(/,2x,a)') 'Automatically generated ellipsoide potential overwritten by:'
        call write_cts_NCI(6,env%cts)
        write (stdout,*)

      end if
    end if
  end do

  if (.not.atomlistused) then
    atlist = 1
    env%includeRMSD = atlist
  end if

  deallocate (atlist)

  return
end subroutine parseRC2

!========================================================================================!

subroutine inputcoords(env,arg)
!***********************************************************************
!* subroutine inputcoords
!* Convert given input coordinate file into a "coord" file (TM format)
!* If the input is in SDF format, document the info to convert the
!* final ensemble back into this format
!***********************************************************************
  use crest_parameters
  use crest_data
  use strucrd
  use axis_module
  use zdata
  use iomod
  implicit none
  !> Input
  type(systemdata) :: env
  character(len=*) :: arg
  !> Local variables
  logical :: ex,ex2
  character(len=:),allocatable :: inputfile
  character(len=:),allocatable :: arg2
  type(coord) :: mol
  type(zmolecule) :: zmol
  integer :: i,idiff
  integer,allocatable :: tmpinclude(:)
  real(wp),allocatable :: savelat(:,:)   !> preserve input lattice across coord round-trip

!>--- Redirect for QCG input reading
  if (env%QCG) then
    !Input coordinates are processed in solvtool.f90 file during solvtool subroutine
    return
  end if
!>---

  inquire (file=arg,exist=ex)
  inquire (file='coord',exist=ex2)
  if (.not.ex.and..not.ex2) then
    if (env%dryrun) then
      write (stdout,*) 'No (valid) input file, but ignoring for dry run.'
      return
    else
      error stop 'No (valid) input file! exit.'
    end if
  end if
  if (ex2) then
!>-- save coord as reference
    call copy('coord','coord.original')
  end if
  if (ex.and.arg(1:1) .ne. '-') then
    call mol%open(arg)
    if (allocated(mol%lat)) savelat = mol%lat
    call mol%write('coord')
    call mol%write('crest_input_copy.xyz')
    call mol%deallocate()
    inputfile = trim(arg)
  else
    inputfile = 'coord'
  end if
  if (.not.allocated(env%inputcoords)) env%inputcoords = inputfile

!>-- if the input was a SDF file, special handling
  env%sdfformat = .false.
  call checkcoordtype(inputfile,i)
  if (any((/31,32/) == i)) then
    env%sdfformat = .true.
    env%outputsdf = .true.
  end if

!>-- after this point there should always be an coord file present
  if (.not.allocated(env%inputcoords)) env%inputcoords = 'coord'
  call mol%open('coord')
  if (allocated(savelat)) call move_alloc(savelat,mol%lat)
!>-- shift to CMA and/or align according to rot.const. We have to be careful about this.
  if (allocated(mol%lat)) then
    !> periodic systems must not be rotated/translated
    continue
  else if (any((/crest_sp,crest_optimize,crest_numhessian,crest_trialopt, &
  &         crest_ensemblesp,crest_ensemblehess/) == env%crestver)) then
    !> some runtypes should only do a CMA translation, but no rotation
    call CMAtrf(mol%nat,mol%nat,mol%at,mol%xyz)
  else if (env%crestver == crest_solv) then
    !> runtypes like qcg must not modify input coordinates!
    continue
  else
    !> all other can align with rot. axis
    call axis(mol%nat,mol%at,mol%xyz)
  end if
!>-- overwrite coord
  call mol%write('coord')

!>-- get the number of atoms and the reduced number of atoms if some of
!>-- them are excluded from the RMSD calc in V2. Initially they are the same
  env%nat = mol%nat
  env%rednat = env%nat
!>-- reference geo save
  env%ref%nat = mol%nat
  env%ref%at = mol%at
  env%ref%xyz = mol%xyz
  if (allocated(mol%lat)) then
    env%ref%lat = mol%lat
!>-- the genetic crossing (GC) step operates on bare internal coordinates and is
!>-- not defined under periodic boundary conditions -> disable it for PBC input
    if (env%performCross) then
      env%performCross = .false.
      write (stdout,'(a)') '> Periodic system detected: genetic crossing (GC) step disabled.'
    end if
  end if
  env%ref%ichrg = env%chrg
  env%ref%uhf = env%uhf
!>-- topology save
  if (any((/crest_mfmdgc,crest_imtd,crest_imtd2,crest_trialopt/) == env%crestver)) then
    if (.not.env%autozsort) then
      env%ref%ntopo = mol%nat*(mol%nat+1)/2
      allocate (env%ref%topo(env%ref%ntopo))
      call quicktopo(mol%nat,mol%at,mol%xyz,env%ref%ntopo,env%ref%topo)
    end if
  end if
  call mol%deallocate()

!>-- for protonation/deprotonation applications get ref. number of fragments
!>-- also get some other structure based info
  call simpletopo_file('coord',zmol,.false.,.false.,'')
  env%protb%nfrag = zmol%nfrag
  call zmol%deallocate()

!>--- Repair logic of includeRMSD array (especially for something like QCG)
  if (.not.allocated(env%includeRMSD)) then
    allocate (env%includeRMSD(env%ref%nat))
    env%includeRMSD(:) = 1
  else
    !> assuming if the current includeRMSD is smaller than the
    !> current system we have *appended* some atoms
    idiff = size(env%includeRMSD,1)
    if (idiff < env%ref%nat) then
      allocate (tmpinclude(env%ref%nat),source=1)
      do i = 1,idiff
        tmpinclude(i) = env%includeRMSD(i)
      end do
      !deallocate(env%includeRMSD)
      call move_alloc(tmpinclude,env%includeRMSD)
    end if
  end if

  return
end subroutine inputcoords

!========================================================================================!
!========================================================================================!

!========================================================================================!
!========================================================================================!
