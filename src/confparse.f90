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
  use lwoniom_module
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
  logical :: ex,bondconst
  character(len=:),allocatable :: argument

  allocate (xx(10),floats(3),strings(3))
  ctmp = ''
  dtmp = ''
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

    write (*,'(/,1x,a)') 'Command line input:'
    write (*,'(1x,a,a,/)') '$ ',trim(cmd)
  end if
  env%cmd = trim(cmd)

!=========================================================================================!
!>--- check if help is requested or citations shall be diplayed
  do i = 1,nra
    if (any((/character(6)::'-h','-H','--h','--H','--help'/) == trim(arg(i)))) then
      if(nra > i)then
        ctmp=trim(arg(i+1))
        if(ctmp(1:1).ne.'-')then
          call confscript_morehelp(ctmp)
        endif
      endif 
      call confscript_help()
    end if
    if (any((/character(10)::'-cite','--cite','--citation'/) == trim(arg(i)))) then
      call crestcite()
    end if
    if (index(arg(i),'-newversion') .ne. 0) then !> as in CREST version >= 3.0
      env%legacy = .false.
    end if
    if (index(arg(i),'-legacy') .ne. 0) then  !> as in CREST version <3.0
      env%legacy = .true.
    end if
    if (index(arg(i),'-dry') .ne. 0) then   !> "dry" run to print settings
      env%dryrun = .true.
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
  env%maxcompare = 10          !> maximum number of (lowest) conformers to compare when using "-compare"
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
  ctype = 0                !> bond constraint type

  env%NCI = .false.      !> use specialized NCI mode?
  env%potscal = 1.0d0    !> scale automatically set wall potential by this factor

!>--- get CHRG and UHF if the respective files are present
  inquire (file='.CHRG',exist=ex)
  if (any(index(arg,'-chrg') .ne. 0)) ex = .false.
  if (ex) then
    call rdshort('.CHRG',env%chrg)
    if (env%chrg .ne. 0) then
      write (*,'(2x,a,i0)') 'molecular charge read from .CHRG:  ',env%chrg
    end if
  end if
  inquire (file='.UHF',exist=ex)
  if (any(index(arg,'-uhf') .ne. 0)) ex = .false.
  if (ex) then
    call rdshort('.UHF',env%uhf)
    if (env%uhf .ne. 0) then
      write (*,'(2x,a,i0)') 'nα-nβ electrons read from .UHF:  ',env%uhf
    end if
  end if

!>--- options for constrained conformer sampling
  env%fixfile = 'none selected'

!>--- options for possible property calculations, mainly protonation/deprotonation/taut. tool
  env%ptb%popthr = 0.01_wp  !> = 1% population
  env%ptb%ewin = 30.0_wp  !> 30 kcal for protonation
  env%ptb%swat = 0
  env%ptb%swchrg = 0
  env%ptb%iter = 1           !> number of iteration cycles for tautomerization
  env%ptb%swelem = .false. !> replace H⁺ in protonation routine by something else?
  env%ptb%allowFrag = .false. !> allow dissociated Structures?
  env%ptb%threshsort = .false. !> use ewin threshold window
  env%ptb%protdeprot = .false. !> (tautomerize) do first protonation and then deprotonation
  env%ptb%deprotprot = .false. !> (tautomerize) do first deprotonation and then protonation
  env%ptb%strictPDT = .false.  !> strict mode (i.e. bond constraints) for (de)protonation,tautomerization
  env%pclean = .false.       !> cleanup option for property mode

!>--- options for principal component analysis (PCA) and clustering
  env%pcmeasure = 'dihedral'

!>--- thermo options
  env%thermo%trange(1) = 278.15d0  !> T start
  env%thermo%trange(2) = 380.0d0   !> T stop (approx.)
  env%thermo%trange(3) = 10.0d0    !> T step
  env%thermo%ptot = 0.9d0   !> for hessians take x% conformers
  env%thermo%pcap = 50000   !> limit number of structures
  env%thermo%sthr = 25.0d0  !> rotor cutoff
  env%thermo%fscal = 1.0d0   !> frequency scaling factor

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
  env%crestver = crest_imtd !> confscript version (v.1 = MF-MD-GC, v.2 = MTD)
  env%runver = 1            !> default
  env%properties = p_none   !> additional calculations/options before or after confsearch
  env%properties2 = p_none  !> backup for env%properties
  env%iterativeV2 = .true.  !> iterative crest V2 version
  env%preopt = .true.
!>--- check for input file
  do i = 1,nra
    argument = trim(arg(i))
    if (argument == '--input'.or.argument == '-i') then
      call parseinputfile(env,trim(arg(i+1)))
      exit
    end if
    if (i == 1.and.index(argument,'.toml') .ne. 0) then
      call parseinputfile(env,trim(arg(1)))
      exit
    end if
  end do
!>--- first arg loop
  do i = 1,nra
    argument = trim(arg(i))
    if (argument(1:2) == '--') then
      argument = argument(2:)
    end if
    if (argument .ne. '') then
      select case (argument) !> RUNTYPES

      case ('-v1') !> confscript version 1 (MF-MD-GC)
        env%crestver = crest_mfmdgc
        write (*,'(2x,a,'' : MF-MD-GC'')') trim(arg(i))
        env%mdtime = 40.0d0       !> simulation length of the MD, 40ps total (2*20ps)(default for QMDFF would be 500)
        env%temps = 1             !> number of default MD cycles
        env%Maxrestart = 15
        env%performModef = .true. !> do the MF in V1
        env%trackorigin = .false. !> for v1 there is not much insight from this
        exit

      case ('-v2') !> confscript version 2 (MTD-GC)
        env%crestver = crest_imtd
        write (*,'(2x,a,'' : MTD-GC'')') trim(arg(i))
        env%iterativeV2 = .false.  !> iterative crest V2 version
        env%Maxrestart = 1       !> for non-iterative MTD-GC only
        exit

      case ('-v3','-v2i') !> confscript version 2 but iterativ (iMTD-GC)
        env%crestver = crest_imtd
        env%iterativeV2 = .true.
        write (*,'(2x,a,'' : iMTD-GC'')') trim(arg(i))
        exit

      case ('-v4') !> sMTD-iMTD (same as entropy mode)
        env%crestver = crest_imtd2
        env%iterativeV2 = .true.
        env%entropymd = .true.
        env%rotamermds = .false.
        env%performCross = .false.
        env%emtd%maxfallback = 1
        write (*,'(2x,a,'' : iMTD-sMTD'')') trim(arg(i))
        exit

      case ('-mdopt','-purge') !> MDOPT
        env%crestver = crest_mdopt
        atmp = ''
        env%preopt = .false.
        env%ensemblename = 'none selected'
        if (nra .ge. (i+1)) atmp = adjustl(arg(i+1))
        if ((atmp(1:1) /= '-').and.(len_trim(atmp) .ge. 1)) then
          env%ensemblename = trim(atmp)
        end if
        call xyz2coord(env%ensemblename,'coord') !> write coord from lowest structure
        env%inputcoords = env%ensemblename !> just for a printout
        exit

      case ('-screen')  !> SCREEN
        env%crestver = crest_screen
        atmp = ''
        env%ensemblename = 'none selected'
        if (nra .ge. (i+1)) atmp = adjustl(arg(i+1))
        if ((atmp(1:1) /= '-').and.(len_trim(atmp) .ge. 1)) then
          env%ensemblename = trim(atmp)
        end if
        call xyz2coord(env%ensemblename,'coord') !write coord from lowest structure
        env%inputcoords = env%ensemblename !just for a printout
        exit

      case ('-mdsp','-ensemblesp') !> Singlepoints along ensemble
        env%crestver = crest_ensemblesp
        atmp = ''
        env%preopt = .false.
        env%ensemblename = 'none selected'
        if (nra .ge. (i+1)) atmp = adjustl(arg(i+1))
        if ((atmp(1:1) /= '-').and.(len_trim(atmp) .ge. 1)) then
          env%ensemblename = trim(atmp)
        end if
        call xyz2coord(env%ensemblename,'coord') !> write coord from lowest structure
        env%inputcoords = env%ensemblename !> just for a printout
        exit


      case ('-pka','-pKa')  !> pKa calculation script
        env%crestver = crest_pka
        env%runver = 33
        !env%relax=.true.
        env%performCross = .false.    !skip the genetic crossing
        env%trackorigin = .false.
        env%Maxrestart = 1
        env%ptb%ewin = 15.0d0
        env%gbsa = .true.
        env%solv = '--alpb h2o'
        env%ptb%h_acidic = 0
        call pka_argparse(arg(i+1),env%ptb%h_acidic)
        if (env%ptb%h_acidic == -2) env%ptb%pka_baseinp = trim(arg(i+1))

      case ('-compare')   !> flag for comparing two ensembles, analysis tool
        env%compareens = .true.
        env%crestver = 5
        env%properties = p_compare
        env%ensemblename = 'none selected'
        env%ensemblename2 = 'none selected'
        if (nra .ge. (i+2)) then
          atmp = adjustl(arg(i+1))
          btmp = adjustl(arg(i+2))
        else
          write (*,'(a,a)') trim(arg(i)),' requires two arguments:'
          write (*,'(2x,a,a)') trim(arg(i)),' [ensemble1] [ensemble2]'
          error stop
        end if
        if ((atmp(1:1) /= '-').and.(len_trim(atmp) .ge. 1).and. &
        &  (btmp(1:1) /= '-').and.(len_trim(btmp) .ge. 1)) then
          env%ensemblename = trim(atmp)
          env%ensemblename2 = trim(btmp)
        end if
        write (*,'(1x,a,1x,a,1x,a)') trim(arg(i)),trim(env%ensemblename),trim(env%ensemblename2)
        exit

      case ('-protonate') !> protonation tool
        env%properties = p_protonate
        write (*,'(2x,a,'' : automated protonation script'')') trim(arg(i))
        exit

      case ('-deprotonate') !> deprotonation tool
        env%properties = p_deprotonate
        write (*,'(2x,a,'' : automated deprotonation script'')') trim(arg(i))
        exit

      case ('-tautomerize') !> tautomerization tool
        env%properties = p_tautomerize
        write (*,'(2x,a,'' : automated tautomerization script'')') trim(arg(i))
        exit

      case ('-isomerize','-stereomers') !> isomerization tool
        env%properties = p_isomerize
        write (*,'(2x,a,'' : automated stereoisomerization script'')') trim(arg(i))
        write (*,'(2x,''Note: Use of GFN-FF required for stereoisomer generation.'')')
        exit

      case ('-forall','-for') !> property mode with ensemble as input
        env%properties = p_propcalc
        atmp = ''
        env%ensemblename = 'none selected'
        if (nra .ge. (i+1)) atmp = adjustl(arg(i+1))
        if ((atmp(1:1) /= '-').and.(len_trim(atmp) .ge. 1)) then
          env%ensemblename = trim(atmp)
        end if
        inquire (file=env%ensemblename,exist=ex)
        if (.not.ex) then
          write (*,'(1x,a,a,a)') 'invalid ensemble file <',trim(env%ensemblename),'>. exit.'
          error stop
        end if
        call xyz2coord(env%ensemblename,'coord') !> write coord from lowest structure
        env%inputcoords = env%ensemblename !> just for a printout
        if (argument == '-forall') then
          env%ptb%alldivers = .true.
        end if
        exit

      case ('-rrhoav')  !> Hessians along given ensemble and average
        env%properties = p_rrhoaverage
        atmp = ''
        env%ensemblename = 'none selected'
        if (nra .ge. (i+1)) atmp = adjustl(arg(i+1))
        if ((atmp(1:1) /= '-').and.(len_trim(atmp) .ge. 1)) then
          env%ensemblename = trim(atmp)
        end if
        inquire (file=env%ensemblename,exist=ex)
        if (.not.ex) then
          write (*,'(1x,a,a,a)') 'invalid ensemble file <',trim(env%ensemblename),'>. exit.'
          error stop
        end if
        exit

      case ('-reactor')  !> xtb nanoreactor workarounds
        env%preopt = .false.
        env%crestver = crest_nano
        exit

      case ('-solvtool','-qcg')
        !> Set solute file if present
        if(i == 2) env%solu_file = trim(arg(i-1))
        !> Set solvent file if prensent
        !> If it is another argument, it doesent matter as solvent file is checke in solvtool
        if (nra >= i+1) env%solv_file = trim(arg(i+1))
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

      case ('-compress')
        env%crestver = crest_compr
        env%runver = 77
        env%mdstep = 2.5d0
        env%mddump = 2000
        env%autozsort = .false.
        exit

    case ('-msreact')
        env%crestver = crest_msreac
        env%preopt = .false.
        env%presp = .true.
        env%ewin = 200.0d0 !> 200 kcal for msreact

      case ('-splitfile')
        ctmp = trim(arg(i+1))
        k = huge(j)
        l = 1
        if (nra >= i+2) then
          read (arg(i+2),*,iostat=io) j
          if (io == 0) then
            k = j
          end if
        end if
        if (nra >= i+3) then
          read (arg(i+3),*,iostat=io) j
          if (io == 0) then
            l = j
          end if
        end if
        call splitfile(ctmp,k,l)
        stop

      case ('-printaniso')
        ctmp = trim(arg(i+1))
        inquire (file=ctmp,exist=ex)
        if (ex) then
          call printaniso(ctmp,0.01_wp,0.025_wp,0.5_wp)
        end if
        stop

      case ('-printboltz')
        if (nra >= i+2) then
          ctmp = trim(arg(i+1))
          dtmp = trim(arg(i+2))
          call prbweight(ctmp,dtmp)
        else
          ctmp = trim(arg(i+1))
          call prbweight(ctmp,'')
        end if

      case ('-wbotopo','-usewbo')  !> try to use a WBO file in topology analysis
        ctmp = trim(arg(i+1))
        if (ctmp(1:1) .ne. '-'.and.(nra >= i+1)) then
          env%wbofile = trim(ctmp)
        else
          env%wbofile = 'wbo'
        end if
        env%wbotopo = .true.

      case ('-testtopo')
        ctmp = trim(arg(i+1))
        inquire (file=ctmp,exist=ex)
        if (i+2 .le. nra) then
          dtmp = trim(arg(i+2))
          if (dtmp(1:1) == '-') then
            dtmp = 'default'
          end if
        end if
        if (ex) then
          call testtopo(ctmp,env,dtmp)
        end if

      case ('-resortensemble')
        ctmp = trim(arg(i+1))
        inquire (file=ctmp,exist=ex)
        if (ex) then
          call resort_ensemble(ctmp)
        end if
        stop

      case ('-thermo','-thermotool')
        env%properties = p_thermo
        ctmp = trim(arg(1))  ! first argument to read the structure
        if (ctmp(1:1) .ne. '-') then
          env%inputcoords = trim(ctmp)
          env%thermo%coords = trim(ctmp)
        end if
        ctmp = trim(arg(i+1)) ! second argument to read the vibspectrum
        if (ctmp(1:1) .ne. '-') then
          env%thermo%vibfile = trim(ctmp)
        end if

      case ('-rmsd','-rmsdheavy','-hrmsd')
        ctmp = trim(arg(i+1))
        dtmp = trim(arg(i+2))
        if ((argument == '-rmsdheavy').or.(argument=='-hrmsd')) then
          call quick_rmsd_tool(ctmp,dtmp,.true.)
        else
          call quick_rmsd_tool(ctmp,dtmp,.false.)
        end if
        stop

      case ('-symmetries')
        ctmp = trim(arg(i+1))
        inquire (file=ctmp,exist=ex)
        if (ex) then
          call ensemble_analsym(trim(ctmp),.true.)
        end if
        stop

      case ('-exlig','-exligand','-exchligand')
        env%properties = p_ligand
        env%ptb%infile = trim(arg(1))
        ctmp = trim(arg(i+1))
        env%ptb%newligand = trim(ctmp)
        read (arg(i+2),*,iostat=io) j
        if (io == 0) then
          env%ptb%centeratom = j
        end if
        read (arg(i+3),*,iostat=io) j
        if (io == 0) then
          env%ptb%ligand = j
        end if
        exit

      case ("-acidbase","-ab",'-abprep','-pkaprep','-gdissprep')  !-- acid base correction
        !> crest --ab <acid.xyz> <base.xyz> --chrg <acidchrg>
        env%properties = p_acidbase
        if (index(arg(i),'prep') .ne. 0) then
          call pka_argparse2(env,arg(i+1),arg(i+2),env%ptb%pka_mode)
        else
          ctmp = trim(arg(i+1))
          inquire (file=ctmp,exist=ex)
          if (ex) then
            env%ptb%pka_acidensemble = trim(ctmp)
            write (*,'(1x,a,a)') 'File used for the acid: ',trim(ctmp)
          end if
          ctmp = trim(arg(i+2))
          inquire (file=ctmp,exist=ex)
          if (ex) then
            env%ptb%pka_baseensemble = trim(ctmp)
            write (*,'(1x,a,a)') 'File used for the base: ',trim(ctmp)
          end if
        end if
        env%solv = '--alpb h2o'
        env%gfnver = '--gfn2'

      case ('-redoextrapol')
        ctmp = trim(arg(i+1))
        read (arg(i+2),*,iostat=io) j
        if (io == 0) then
          call redo_extrapol(ctmp,j)
        else
          call redo_extrapol(ctmp,0)
        end if
        stop

      case ('-sp') !> singlepoint calculation (uses new calculator routines)
        env%crestver = crest_sp
        env%preopt = .false.
        env%legacy = .false.
        exit

      case ('-opt','-optimize','-ancopt','-ohess') !> ANCOPT structure optimization (uses new calculator routines)
        env%preopt = .false.
        env%crestver = crest_optimize
        env%legacy = .false.
        if(argument.eq.'-ohess') env%crest_ohess=.true.
        exit

      case ('-trialopt')  !> test optimization with topocheck
        env%preopt = .false.
        env%crestver = crest_trialopt
        exit

      case ('-dynamics','-dyn') !> molecular dynamics (uses new calculator routines)
        env%preopt = .false.
        env%crestver = crest_moldyn
        env%legacy = .false.
        exit

      case ('-SANDBOX')
        !>--- IMPLEMENT HERE WHATEVER YOU LIKE, FOR TESTING
        !>-----
        stop
      case ('-PLAYGROUND','-TEST')
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
  end if
!========================================================================================!
!> after this point there should always be a "coord" file present
!========================================================================================!
  allocate (env%includeRMSD(env%nat))
  env%includeRMSD = 1

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> parse the input flags
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
  do i = 1,nra
    argument = trim(arg(i))
    if (argument(1:2) == '--') then
      argument = argument(2:)
    end if
    if (argument .ne. '') then
!========================================================================================!
!-------- switch between legacy (systemcall) and new code (API) implementations
!========================================================================================!
      select case (argument)
      case ('-legacy')      !> switch to old xtb-call version where possible
        env%legacy = .true.
      case ('-newversion')  !> switch to newer implementations (CREST >3.0)
        env%legacy = .false.
      end select
!========================================================================================!
!------- flags exclusively for V2 (MTD-GC) V2i/V3 (iMTD-GC) V4 (iMTD-sMTD)
!========================================================================================!
      if (any((/crest_imtd,crest_imtd2,11/) == env%crestver)) then
        select case (argument) !> V2
        case ('-mdtemp')                          !> set MTD temperature (V2 version)
          call readl(arg(i+1),xx,j)
          env%mdtemp = xx(1)
          env%user_temp = .true.
        case ('-quick')                           !> performing quick conformational search
          env%quick = .true.
          env%runver = 2
          env%ewin = 5.0d0
          if (env%optlev > 1.0d0) env%optlev = 1.0d0    !> optlev tight for quick run
        case ('-shake')                           !> set shake
          call readl(arg(i+1),xx,j)
          env%shake = nint(xx(1))
        case ('-tstep')                           !> set MD timestep in fs
          call readl(arg(i+1),xx,j)
          env%mdstep = xx(1)
          env%user_mdstep = .true.
        case ('-vbdump')                          !> Vbias dump in ps
          call readl(arg(i+1),xx,j)
          xx(2) = xx(1)*1000
          env%mddump = nint(xx(2))
        case ('-mdskip')                          !> set skipping structures in -mdopt
          call readl(arg(i+1),xx,j)
          env%mdskip = nint(xx(1))
        case ('-mddump')                          !> set dumpstep for writing structures from MD
          call readl(arg(i+1),xx,j)
          env%mddumpxyz = nint(xx(1))
        case ('-nomtd')                           !> Don't do the MTD in V2
          env%performMTD = .false.
        case ('-restartopt')                      !> go to step 2 of multilevel optimization immideatly
          env%restartopt = .true.
          env%autozsort = .false.
        case ('-norotmd')                         !> don't do the regular mds after step 2 in multilevel optimization of V2
          env%rotamermds = .false.
        case ('-rotmd')
          env%rotamermds = .true.
        case ('-tnmd')                            !> temperature for additional normal MDs
          call readl(arg(i+1),xx,j)
          env%nmdtemp = xx(1)
        case ('-gcmopt')                          !> GC multilevel optimization activate in V2
          env%gcmultiopt = .true.
        case ('-gcsopt')                          !> GC single level optimization in V2
          env%gcmultiopt = .false.
        case ('-nogcmopt')                        !> GC single level optimization in V2
          env%gcmultiopt = .false.
        case ('-qmdff')                           !> use QMDFF for the MDs in V2?
          env%useqmdff = .true.
        case ('-nci')                             !> NCI special mode
          write (*,'(2x,a,1x,a)') trim(arg(i)),' : Special NCI mode for non-covalently bound complexes or clusters.'
          env%NCI = .true.
          env%runver = 4
          env%autozsort = .false.
          env%performCross = .false.
          env%rotamermds = .false.
        case ('-squick','-superquick')            !> extremely crude quick mode
          write (*,'(2x,a,1x,a)') trim(arg(i)),' : very crude quick-mode (no NORMMD, no GC, crude opt.)'
          env%rotamermds = .false.      !> no NORMMD
          env%performCross = .false.    !> no GC
          env%quick = .true.            !> MTD settings from the quick-mode
          env%superquick = .true.       !> use user-set opt level in Multilevel opt.
          env%runver = 5
          if (env%optlev > 0.0d0) env%optlev = 0.0d0    !> user-set opt level
          env%ewin = 5.0d0              !> smaller energy window

        case ('-mquick','-megaquick')   !> extremely crude quick mode pt.2
          write (*,'(2x,a,1x,a)') trim(arg(i)),' : very crude quick-mode (no NORMMD, no GC, crude opt.)'
          env%rotamermds = .false.      !> no NORMMD
          env%performCross = .false.    !> no GC
          env%quick = .true.            !> MTD settings from the quick-mode
          env%superquick = .true.       !> use user-set opt level in Multilevel opt.
          env%Maxrestart = 1            !> only one MTD iteration
          env%runver = 6
          if (env%optlev > 0.0d0) env%optlev = 0.0d0  !> user-set opt level
          env%ewin = 2.5d0              !> smaller energy window

        case ('-extensive')   !> counterpart to quick mode
          env%slow = .true.
          env%quick = .false.
          env%superquick = .false.
          env%optlev = 0.0d0
          env%ewin = 8.0d0
          env%runver = 8
        case ('-static','-staticmtd')
          env%staticmtd = .true.
        case default
          continue
        end select !> V2
        !--- iterative version of V2
        if (env%iterativeV2) then
          select case (argument) !> V2i
          case ('-mrest')                  !> set max number of restarts
            call readl(arg(i+1),xx,j)
            env%Maxrestart = nint(xx(1))
          case ('-iru')                    !> re-use previously found conformers as bias in iterative approach
            env%iru = .true.
          case ('-keepdir','-keeptmp')     !> Do not delete METADYN and NORMMD directories
            env%keepModef = .true.
          case ('-singlerun')              !> QCG special mode
            write (*,'(2x,a,1x,a)') trim(arg(i)),' : run mode with only a single MTD and no iterations (for testing)'
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
          env%optpurge = .true.
        case ('-ethrpurge','-ethrp')
          read (arg(i+1),*,iostat=io) rdum
          if (io == 0) env%ethrpurge = rdum
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
          if (i+1 .le. nra) then
            atmp = trim(arg(i+1))
            if (atmp(1:1) .ne. '-') then
              call readl(arg(i+1),xx,j)
              env%rdens = xx(1)
            end if
          end if
          env%properties = p_reactorset
          env%preactorpot = .true.
        case ('-genmtd')
          env%properties = p_reactorset
          env%mdtime = 20.0d0
          if (i+1 .le. nra) then
            atmp = trim(arg(i+1))
            if (atmp(1:1) .ne. '-') then
              call readl(arg(i+1),xx,j)
              env%mdtime = xx(1)
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
          env%restartopt = .true.
        case ('-iso')
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
        case ('-keepdir','-keeptmp')
          env%keepModef = .true.
        case ('-tstep')                           !> set MD timestep in fs
          call readl(arg(i+1),xx,j)
          env%mdstep = xx(1)
          env%user_mdstep = .true.
        case ('-vbdump')                          !> Vbias dump in ps
          call readl(arg(i+1),xx,j)
          xx(2) = xx(1)*1000
          env%mddump = nint(xx(2))
        case ('-mdskip')                          !> set skipping structures in -mdopt
          call readl(arg(i+1),xx,j)
          env%mdskip = nint(xx(1))
        case ('-mddump')                          !> set dumpstep for writing structures out of the md
          env%user_dumxyz = .true.
          call readl(arg(i+1),xx,j)
          env%mddumpxyz = nint(xx(1))
        case ('-nomtd')                           !> Don't do the MTD in V2
          env%performMTD = .false.
        case ('-wscal')                           !> scale size of wall potential
          call readl(arg(i+1),xx,j)
          env%potscal = xx(1)
          env%user_wscal = .true.
        case ('-fixsolute')                       !> Fix the solute after CMA trafo
          env%constrain_solu = .true.
        case ('-nofix')                           !> No fixing of the solute after CMA trafo
          env%noconst = .true.
        case ('-restartopt')                      !> go to step 2 of multilevel optimization immideatly
          env%restartopt = .true.
          env%autozsort = .false.
        case ('-norotmd')                         !> don't do the regular mds after step 2 in multilevel optimization of V2
          env%rotamermds = .false.
        case ('-mdtemp')                          !> set MTD temperature (V2 version)
          call readl(arg(i+1),xx,j)
          env%mdtemp = xx(1)
          env%user_temp = .true.
        case ('-tnmd')                            !> temperature for additional normal MDs
          call readl(arg(i+1),xx,j)
          env%nmdtemp = xx(1)
        case ('-shake')                           !> set shake
          call readl(arg(i+1),xx,j)
          env%shake = nint(xx(1))
        end select !> QCG
      end if

!========================================================================================!
!------- Flags for msreact
!========================================================================================!
      if (env%crestver == crest_msreac) then
        select case (argument) !> msreact
        case('-msnoiso') !> filter out non fragmentated structures in msreact
          env%msnoiso=.true.
        case('-msiso') !> filter out fragmentated structures in msreact
          env%msiso=.true.
        case('-msnbonds') ! give number of bonds up to which bias potential is added between atoms default 3
          call readl(arg(i + 1),xx,j)
          env%msnbonds = xx(1)
        case('-msnshifts') ! give number of times atoms are randomly shifted before optimization
          call readl(arg(i + 1),xx,j)
          env%msnshifts = xx(1)
        case('-msnshifts2') ! give number of times atoms are randomly shifted before applying the constrained optimization default 0 
          call readl(arg(i + 1),xx,j)
          env%msnshifts2 = xx(1)
        case('-msnfrag') ! give number of structures that should be generated
          call readl(arg(i + 1),xx,j)
          env%msnfrag = xx(1)
        case('-msmolbar') !> filter out structures with same molbar code in msreact
          env%msmolbar=.true.
        case('-msinchi') !> filter out structures with same inchi code in msreact
            env%msinchi=.true.
        case('-msnoattrh') !> add attractive potential for H-atoms
          env%msattrh=.false.
        case('-mslargeprint') !> additional printouts and keep MSDIR
          env%mslargeprint=.true.
        case('-msinput') ! give number of times atoms are randomly shifted before applying the constrained optimization default 0 
            ctmp = trim(arg(i+1))
            if (ctmp(1:1) .ne. '-') then
              env%msinput = trim(ctmp)
            end if
        end select !> msreact
      end if
!========================================================================================!
!------- other general flags
!========================================================================================!
      select case (argument) !> ARGPARSER1
      case ('-dry')             !> "dry" run to print settings
        env%dryrun = .true.
      case ('-nozs')
        env%autozsort = .false.   !> turn off automatic zsort (default)
      case ('-zs')
        env%autozsort = .true.    !> turn on automatic zsort
      case ('-nocross')
        env%performCross = .false.    !> skip the genetic crossing
        write (*,'(2x,a,1x,a)') trim(arg(i)),' : skipping GC part.'
      case ('-cross')
        env%performCross = .true.     !> do the genetic crossing
        env%autozsort = .true.
      case ('-opt','-optlev')             !> settings for optimization level of GFN-xTB
        env%optlev = optlevnum(arg(i+1))
        write (*,'(2x,a,1x,a)') trim(arg(i)),optlevflag(env%optlev)
      case ('-gfn','-gfn1','-gfn2','-gfn0','-gff','-gfnff')
        ctmp = argument
        if (argument == '-gfn') then
          dtmp = trim(arg(i+1))
          ctmp = ctmp//dtmp
        end if
        if (env%properties == p_isomerize) then
          ctmp = 'stereoisomers'
        end if
        select case (ctmp) !> GFN
        case ('-gfn1')
          env%gfnver = '--gfn1'
          write (*,'(2x,a,'' : Use of GFN1-xTB requested.'')') env%gfnver
        case ('-gfn2')
          env%gfnver = '--gfn2'
          write (*,'(2x,a,'' : Use of GFN2-xTB requested.'')')  env%gfnver
        case ('-gfn0')
          env%gfnver = '--gfn0'
          write (*,'(2x,a,'' : Use of GFN0-xTB requested.'')') env%gfnver
        case ('-gff','-gfnff')
          env%gfnver = '--gff'
          write (*,'(2x,a,'' : Use of GFN-FF requested.'')') '--gfnff'
          env%mdstep = 1.5d0
          env%hmass = 5.0d0
          ctype = 5 !> bond constraint activated
          if (any((/crest_imtd,crest_imtd2/) == env%crestver)) then
            bondconst = .true.
          endif
          env%cts%cbonds_md = .true.
          env%checkiso = .true.
        case ('stereoisomers')
          env%gfnver = '--gff'
        case default
          env%gfnver = '--gfn2'
        end select !> GFN
      case ('-gfn2@gfn0','-gfn2@gfn1','-gfn2@gff','-gfn2@ff','-gfn2@gfnff')
        if (.not.env%legacy) then !TODO
          write (*,'("> ",a,1x,a)') argument,'option not yet available with new calculator'
          error stop
        end if
        select case (argument) !> GFN2ON
        case ('-gfn2@gfn0')
          env%gfnver = '--gfn0'
        case ('-gfn2@gfn1')
          env%gfnver = '--gfn1'
        case ('-gfn2@gff','-gfn2@ff','-gfn2@gfnff')
          env%gfnver = '--gff'
          env%mdstep = 2.0d0
        case default
          env%gfnver = '--gfn2'
        end select !> GFN2ON
        env%gfnver2 = '--gfn2'
        call env%addjob(51)
        call env%checkhy()
        env%reweight = .false.
      case ('-gfn2//gfnff')
        if (.not.env%legacy) then !TODO
          write (*,'("> ",a,1x,a)') argument,'option not yet available with new calculator'
          error stop
        end if
        env%gfnver = '--gff'
        env%mdstep = 2.0d0
        env%gfnver2 = '--gfn2'
        env%reweight = .true.
        env%mdstep = 2.0d0
        env%hmass = 4.0d0
        ctype = 1 !> bond constraint
        bondconst = .true.
        env%cts%cbonds_md = .true.
        env%checkiso = .true.
        if (index(arg(i+1),'opt') .ne. 0) then
          env%altopt = .true.
          write (*,'(2x,a,a)') argument,' : GFN-FF MDs + GFN2 opt.'
        else
          write (*,'(2x,a,a)') argument,' : energy reweighting'
        end if

      case('-refine','-rsp','-ropt') !> add one refinement step (via cmd only one is possible)
        env%legacy = .false. !> new calculators only!
        if(nra >= i+1)then
          env%gfnver2 = trim(arg(i+1))
          write (*,'(2x,a,1x,a,a)') argument,trim(env%gfnver2), &
          & ' : adding refinement step (singlepoint on optimized structures)'
        endif
        
      case ('-charges') !> read charges from file for GFN-FF calcs.
        ctmp = trim(arg(i+1))
        if ((len_trim(ctmp) < 1).or.(ctmp(1:1) == '-')) then
          ctmp = 'charges'
        end if
        inquire (file=ctmp,exist=ex)
        if (ex) then
          env%chargesfilename = ctmp
          env%chargesfile = .true.
          write (*,'(2x,a,a,a)') '-charges: file <',trim(ctmp),'> used for atomic charges'
          call env%ref%rdcharges(env%chargesfilename,idum)
          if (idum .ne. env%chrg) then
            write (*,'(12x,a,i0)') 'with total summed up molecular charge: ',idum
            env%chrg = idum
            env%ref%ichrg = idum
          end if
        end if
      case ('-dscal','-dispscal','-dscal_global','-dispscal_global')
        env%cts%dispscal_md = .true.
        if (index(argument,'_global') .ne. 0) then
          env%cts%dispscal_global = .true.
        end if
        if (nra .ge. i+1) then
          ctmp = trim(arg(i+1))
          read (ctmp,*,iostat=io) rdum
          if (io .eq. 0) env%cts%dscal = rdum
        end if
      case ('-mtd_kscal','-mtdkscal')
         call readl(arg(i+1),xx,j)
         env%mtd_kscal = xx(1)
      case ('-norestart')
        env%allowrestart = .false.
      case ('-readbias')
        env%readbias = .true.
      case ('-useonly')
        env%properties = p_useonly
        env%autozsort = .false.
        env%dummypercent = 1.0_wp
        if (nra .ge. i+1) then
          atmp = adjustl(arg(i+1))
          if (atmp(1:1) .ne. '-') then
            read (atmp,*) env%dummypercent
          end if
        end if

      case ('-gbsa','-g','-alpb')   !> use GBSA or ALPB implicit solvation
        env%gbsa = .true.
        atmp = adjustl(arg(i+1))
        if (atmp(1:1) .ne. '-'.and.atmp(1:1) .ne. ' ') then
          env%solvent = arg(i+1)
          if (trim(argument) == '-alpb') then
            env%solv = '--alpb '//trim(env%solvent)
          else
            env%solv = '--gbsa '//trim(env%solvent)
          end if
        end if
        write (*,'(2x,a,a)') trim(env%solv),' : implicit solvation'

      case ('-chrg')                  !> create a .CHRG file
        call readl(arg(i+1),xx,j)
        open (newunit=ich,file='.CHRG')
        env%chrg = nint(xx(1))
        env%ref%ichrg = env%chrg
        write (ich,'(i0)') nint(xx(1))
        close (ich)
        write (*,'(2x,a,1x,a)') trim(arg(i)),trim(arg(i+1))

      case ('-uhf')                    !> create a .UHF file
        call readl(arg(i+1),xx,j)
        open (newunit=ich,file='.UHF')
        env%uhf = nint(xx(1))
        env%ref%uhf = env%uhf
        write (ich,'(i0)') nint(xx(1))
        close (ich)
        write (*,'(2x,a,1x,a)') trim(arg(i)),trim(arg(i+1))

      case ('-len','-mdlen','-mdtime') !> set md length in ps
        atmp = arg(i+1)
        call to_lower(atmp)
        j = index(atmp,'x')
        env%user_mdtime = .true.
        if (j .ne. 0) then             !> scaling of the md length
          btmp = atmp(j+1:)
          env%scallen = .true.
          call readl(btmp,xx,j)
          env%mdlenfac = xx(1)
        else                           !> direct setting of the md length
          call readl(arg(i+1),xx,j)
          env%mdtime = xx(1)
          write (*,'(2x,a,1x,a,1x,a)') trim(arg(i)),trim(arg(i+1)), &
          &    '(MD length in ps)'
        end if
      case ('-mdscal','-lenscal')       !> scale md length
        env%scallen = .true.
        call readl(arg(i+1),xx,j)
        env%mdlenfac = xx(1)
      case ('-nmtd') !> set number of MTDs
        env%runver = 787878
        call readl(arg(i+1),xx,j)
        env%nmetadyn = nint(xx(1))
      case ('-gcmax','-setgcmax')       !> set maximum number of structures for GC
        env%setgcmax = .true.
        call readl(arg(i+1),xx,j)
        env%gcmax = xx(1)
      case ('-xnam')                    !> select a name for the xTB executeable
        env%ProgName = trim(arg(i+1))
        write (*,'(2x,''-xnam :'')')
        write (*,'(5x,''xtb executable was set to: "'',a,''"'')') trim(env%ProgName)
      case ('-niceprint')               !> progres bar printout
        env%niceprint = .true.
      case ('-origin')                  !> track the origin (i.e. the generation step) of each conformer
        env%trackorigin = .true.
        write (*,'(2x,a,1x,a)') trim(arg(i)),': tracking conformer origins.'
      case ('-constrain')               !> provide a list of atoms to write a .xcontrol.sample
        ctmp = trim(arg(i+1))
        call quick_constrain_file('coord',env%nat,env%ref%at,ctmp)
      case ('-nocbonds')
        bondconst = .false.
        env%cts%cbonds_global = .false.
        env%cts%cbonds_md = .false.
        inquire (file='bondlengths',exist=ex)
        if (ex) call remove('bondlengths')
      case ('-cbonds','-cbonds_md','-cbonds_ez')  !> constrain all bonds
        ctmp = trim(arg(i+1))
        if (ctmp(1:1) .ne. '-') then
          read (ctmp,*,iostat=io) rdum
          if (io .eq. 0) env%forceconst = rdum
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
        ctmp = trim(arg(i+1))
        if (ctmp(1:1) .ne. '-') then
          read (ctmp,*,iostat=io) rdum
          if (io .eq. 0) env%forceconst = rdum
        end if
        ctype = 2
        bondconst = .true.
        env%cts%cbonds_global = .true.
        if (index(argument,'_md') .ne. 0) then
          env%cts%cbonds_md = .true.
          env%cts%cbonds_global = .false.
        end if
      case ('-cheavy','-fixheavy','-cheavy_md')  !> constrain all heavy atom bonds
        ctmp = trim(arg(i+1))
        if (ctmp(1:1) .ne. '-') then
          read (ctmp,*,iostat=io) rdum
          if (io .eq. 0) env%forceconst = rdum
        end if
        ctype = 3
        bondconst = .true.
        env%cts%cbonds_global = .true.
        if (index(argument,'_md') .ne. 0) then
          env%cts%cbonds_md = .true.
          env%cts%cbonds_global = .false.
        end if
      case ('-clight','-fixhyd','-clight_md')  !> constraint all X-H bonds
        ctmp = trim(arg(i+1))
        if (ctmp(1:1) .ne. '-') then
          read (ctmp,*,iostat=io) rdum
          if (io .eq. 0) env%forceconst = rdum
        end if
        ctype = 4
        bondconst = .true.
        env%cts%cbonds_global = .true.
        if (index(argument,'_md') .ne. 0) then
          env%cts%cbonds_md = .true.
          env%cts%cbonds_global = .false.
        end if
      case ('-cfile','-cinp')                   !> specify the constrain file
        ctmp = trim(arg(i+1))
        if (ctmp(1:1) .ne. '-') then
          env%constraints = trim(ctmp)
          write (*,'(2x,a,1x,a)') '--cinp :',trim(ctmp)
        end if
      case ('-fc','-forceconstant')
        ctmp = trim(arg(i+1))
        if (i+1 >= nra) then
          call readl(arg(i+1),xx,j)
          env%forceconst = xx(1)
        end if
        write (*,'(2x,a,f6.4,a)') '-fc ',env%forceconst,': selected force constant in Eh'
      case ('-nomlo','-no-multilevel')   !> turn off multilevel optimization
        env%multilevelopt = .false.
      case ('-normmd')  !> set number of normMDs
        env%rotamermds = .true.
        if (i+1 .le. nra) then
          call readl(arg(i+1),xx,j)
          env%nrotammds = nint(xx(1))  !> how many lowest conformers?
        end if
        if (i+2 .le. nra) then
          call readl(arg(i+2),xx,j)
          env%temps = nint(xx(1))      !> how many different temperatures
        end if
      case ('-rmsdpot','-gesc')
        ctmp = trim(arg(i+1))
        inquire (file=ctmp,exist=ex)
        if (ex) then
          env%cts%usermsdpot = .true.
          call getcwd(atmp)
          env%cts%rmsdpotfile = trim(atmp)//'/'//ctmp
          write (*,'(2x,a,a,a,a)') argument,': using <',ctmp,'> as bias'
        else
          write (*,'(a,a)') argument,': Warning! File could not be found!'
        end if
      case ('-mergebias','-mergebias+','-gesc+')
        env%properties = -9224
        if (index(argument,'+') > 0) env%properties = p_gesc2
        ctmp = trim(arg(i+1))
        inquire (file=ctmp,exist=ex)
        if (ex) then
          env%biasfile = ctmp
        end if
        env%autozsort = .false.
      case ('-gescopt')
        env%gescoptlev = optlevnum(arg(i+1))
      case ('-gescheavy','-heavygesc','-gesc_heavy')
        env%cts%gesc_heavy = .true.
      case ('-rthr2') !> bias rmsd threshold
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0) env%rthr2 = rdum
      case ('-kshift')
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0) env%kshift = rdum
        env%kshiftnum = 1
      case ('-hflip')
        env%doOHflip = .true.
      case ('-noflip')
        env%doOHflip = .false.
      case ('-maxflip')
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0.and.(index(arg(i+1),'-') .eq. 0)) then
          env%maxflip = nint(rdum)
        end if
      case ('-osdf')
        env%outputsdf = .true.
        write (*,'(2x,a," :",1x,a)') trim(arg(i)), &
        & "output ensemble requested in sdf format"

      case ('-wscal')                           !> scale size of wall potential
        call readl(arg(i+1),xx,j)
        env%potscal = xx(1)
      case ('-wpad')                            !> scale size of wall potential
        call readl(arg(i+1),xx,j)
        env%potpad = xx(1)
      case ('-watoms','-wat')
        ctmp = arg(i+1)
        if(ctmp(1:1) .ne. '-')then
           env%potatlist = trim(ctmp)
           write(*,*) env%potatlist 
        endif
      case ('-wall')
        env%wallsetup = .true.
        write (*,'(2x,a,1x,a)') '--wall:','requesting setup of wall potential'
      case ('-wallxl','-wall-xl')
        env%wallsetup = .true.
        env%potscal = 1.5_wp
        write (*,'(2x,a,1x,a)') '--wall-xl:','requesting setup of wall potential (x1.5 size)'
      case ('-wallxxl','-wall-xxl')
        env%wallsetup = .true.
        env%potscal = 2.0_wp
        write (*,'(2x,a,1x,a)') '--wall-xxl:','requesting setup of wall potential (x2.0 size)'
!========================================================================================!
!------ flags for parallelization / disk space
!========================================================================================!
      case ('-T','-P','-parallel')  !> set total number of OMP threads, this replaces -P and -O entirely
        call readl(arg(i+1),xx,j)
        if (index(arg(i+1),'-') .ne. 0) xx = 0d0
        env%Threads = nint(xx(1))
        env%autothreads = .true.
        env%threadssetmanual = .true.
        write (*,'(2x,a,1x,i0,1x,a)') trim(arg(i)),nint(xx(1)), &
        &     '(CPUs/Threads selected)'
      case ('-inplace')     !> activate "in-place" mode for optimizations (ON by default)
        env%inplaceMode = .true.
!========================================================================================!
!------- CREGEN related flags
!========================================================================================!
      case ('-cregen','-oldcregen')  !> CREGEN standalone use
        env%confgo = .true.
        env%properties = p_cregen
        env%autozsort = .false.
        atmp = ''
        env%ensemblename = 'none selected'
        if (nra .ge. (i+1)) atmp = adjustl(arg(i+1))
        if ((atmp(1:1) /= '-').and.(len_trim(atmp) .ge. 1)) then
          env%ensemblename = trim(atmp)
        end if
        if (index(env%ensemblename,'none selected') .ne. 0) then
          write (*,'(2x,a,1x,a)') trim(arg(i)),': CREGEN standalone usage.'
        else
          write (*,'(2x,a,1x,a,a,a)') trim(arg(i)),': CREGEN standalone usage. Sorting file <', &
          & trim(env%ensemblename),'>'
        end if
        if (trim(arg(i)) .eq. '-oldcregen') then
          write (*,'(3x,a)') 'Using the old version of the CREGEN subroutine.'
          env%newcregen = .false.
        end if
      case ('-oldcr')
        write (*,'(3x,a)') 'Using the old version of the CREGEN subroutine.'
        env%newcregen = .false.
        env%ethr = 0.1d0 !> ETHR old value
      case ('-enso')             !> compare two given ensembles
        env%ENSO = .true.
      case ('-compare')          !> compare two given ensembles
        env%compareens = .true.
      case ('-maxcomp')          !> maximum number of lowest conformers to compare with "-compare"
        call readl(arg(i+1),xx,j)
        env%maxcompare = nint(xx(1))
      case ('-ewin')             !> set energy threshold in kcal/mol
        call readl(arg(i+1),xx,j)
        env%ewin = abs(xx(1))
        if (any((/p_protonate,p_deprotonate,p_tautomerize/) == env%properties)) then
          env%ptb%ewin = abs(xx(1))
        end if
        write (*,'(2x,a,1x,a)') trim(arg(i)),trim(arg(i+1))
      case ('-rthr')             !> set RMSD thr
        call readl(arg(i+1),xx,j)
        env%rthr = xx(1)
        write (*,'(2x,a,1x,a)') trim(arg(i)),trim(arg(i+1))
      case ('-ethr')             !> set E thr
        call readl(arg(i+1),xx,j)
        env%ethr = xx(1)
        write (*,'(2x,a,1x,a)') trim(arg(i)),trim(arg(i+1))
      case ('-bthr')             !> set rot const thr
        call readl(arg(i+1),xx,j)
        env%thresholds(4) = xx(1)  !> legacy
        env%bthr2 = xx(1)
        write (*,'(2x,a,1x,a)') trim(arg(i)),trim(arg(i+1))
      case ('-allrot')           !> use all rotational constants for comparison, instead of mean
        env%allrot = .true.
      case ('-athr')             !> set int. rotation. equal atoms for NMR thr
        call readl(arg(i+1),xx,j)
        env%athr = xx(1)
        write (*,'(2x,a,1x,a)') trim(arg(i)),trim(arg(i+1))
      case ('-pthr')             !> set population thr
        call readl(arg(i+1),xx,j)
        rdum = min(1.0_wp,xx(1)) !--> pthr <= 1
        rdum = max(0.0_wp,rdum)  !--> pthr >= 0
        env%pthr = rdum
        write (*,'(2x,a,1x,f6.4)') trim(arg(i)),rdum !trim(arg(i+1))
      case ('-eqv')
        env%doNMR = .true. !> option for the very last confg call
      case ('-zsort')
        env%onlyZsort = .true.                                 !perform only the zsort subroutine
        env%autozsort = .true.                                 ! CB: needs to be set to run zsort
        write (*,'(2x,a,1x,a)') trim(arg(i)),' : only using the ZSORT subroutine.'
      case ('-metac')                                        !automatic complete of mag. and chem. methyl equivalencies
        env%methautocorr = .true.
      case ('-esort')          !> cregen legacy option
        env%esort = .true.
      case ('-debug')
        env%cgf(1) = .true.    !> debug option for confg
      case ('-nowr')
        env%cgf(2) = .false.   !> newfile option for confg
      case ('-eqan')
        env%cgf(3) = .true.    !> equivalence analysis on (for NMR)
      case ('-noeqan')
        env%cgf(3) = .false.   !> equivalence analysis off (for nmr)
      case ('-rot')
        env%cgf(5) = .false.   !> just rotamer check
      case ('-nmr')            !> NMR mode for confscript
        env%doNMR = .true.
        env%optlev = 2.0d0
      case ('-fullcre')
        env%doNMR = .true.
        env%fullcre = .true.
      case ('-heavy')
        env%cgf(4) = .true.   !> perform just the heavy atom RMSD
        env%heavyrmsd = .true.
      case ('-temp')
        ctmp = trim(arg(i+1))
        if (index(ctmp,'-') .eq. 0) then
          call readl(arg(i+1),xx,j)
          env%tboltz = xx(1)
        end if
      case ('-prsc')                !> write scoord files
        env%printscoords = .true.
      case ('-noprsc')              !> don't write scoord files
        env%printscoords = .false.
      case ('-subrmsd')             !> use only the RMSD for atoms that are included in the MTD
        env%subRMSD = .true.
      case ('-noopt')               !> skip the pre-optimization with GFNn-xTB before the confsearch
        env%preopt = .false.
      case ('-topo','-topocheck')
        env%checktopo = .true.
      case ('-notopo','-notopocheck')
        env%checktopo = .false.
        ctmp = trim(arg(i+1))
        if (ctmp(1:1) .ne. '-') then
          call parse_topo_excl(env,ctmp)
          if (allocated(env%excludeTOPO)) then
            env%checktopo = .true.
          end if
        end if
      case ('-noreftopo')
        env%reftopo = .false.
      case ('-ezcheck','-checkez')
        env%checkiso = .true.
      case ('-noezcheck','-nocheckez')
        env%checkiso = .false.
!========================================================================================!
!-------- PROPERTY CALCULATION related flags
!========================================================================================!
      case ('-protonate')             !> protonation tool
        env%properties = p_protonate
        env%autozsort = .false.
        env%ptb%threshsort = .true.
      case ('-swel')                  !> switch out H+ to something else in protonation script
        if (env%properties .eq. -3) then
          call swparse(arg(i+1),env%ptb)
        end if
      case ('-deprotonate')           !> deprotonation tool
        env%properties = p_deprotonate
        env%autozsort = .false.
        env%ptb%threshsort = .true.
      case ('-tautomerize')           !> tautomerization tool
        env%properties = p_tautomerize
        env%autozsort = .false.
        env%ptb%threshsort = .true.
      case ('-tautomerize2','-exttautomerize')
        if (env%properties == p_propcalc) then
          env%properties = p_tautomerize2
        else
          call env%addjob(abs(p_tautomerize2))
        end if
        env%autozsort = .false.
        env%ptb%threshsort = .true.
        env%runver = 33
        env%relax = .true.
        env%performCross = .false.  !> skip the genetic crossing
        env%trackorigin = .false.
        env%Maxrestart = 1
      case ('-relax')
        env%runver = 33
        env%relax = .true.
        env%performCross = .false.  !> skip the genetic crossing
        env%trackorigin = .false.
        env%Maxrestart = 1
      case ('-trev','-tdp')
        env%ptb%deprotprot = .true. !> switch to deprotonation-first mode in tautomerization
      case ('-iter')                !> number of Protonation/Deprotonation cycles in Tautomerization
        call readl(arg(i+1),xx,j)
        env%ptb%iter = nint(xx(1))
      case ('-texcl','-blacklist')
        ctmp = trim(arg(i+1))
      case ('-strict')
        env%ptb%strictPDT = .true.
      case ('-verystrict','-vstrict')
        env%ptb%strictPDT = .false.
        env%ptb%fixPDT = .true.
      case ('-fstrict')
        env%ptb%strictPDT = .true.
        env%ptb%fixPDT = .true.
      case ('-corr','-abcorr')
        env%ptb%strictPDT = .true.
        env%ptb%fixPDT = .true.
        env%ptb%ABcorrection = .true.
      case ('-pkaensemble')
        env%preopt = .false.
        env%presp = .false.
        call pka_argparse2(env,arg(i+1),arg(i+2),env%ptb%pka_mode)
      case ('-pkaparam')
        env%ptb%rdcfer = .true.
        if (i+1 .le. nra) then
          ctmp = trim(arg(i+1))
          if (ctmp(1:1) .ne. '-') then
            env%ptb%cferfile = ctmp
          end if
        end if
!========================================================================================!
!--------- ENTROPY related settings
!========================================================================================!
      case ('-entropy','-entropic')  !> new, specialized calculation of molecular entropies
        write (*,'(2x,a,'' : enhanced ensemble entropy calculation'')') trim(arg(i))
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
          ctmp = trim(arg(i+1))    !> second argument can be the temperature
          if (index(ctmp,'-') .eq. 0) then
            call readl(arg(i+1),xx,j)
            env%tboltz = xx(1)
          end if
        end if
        call env%addjob(env%properties)

      case ('-scthr','-entropy_cthr')
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0) env%emtd%confthr = rdum
      case ('-ssthr','-entropy_sthr')
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0) env%emtd%sconvthr = rdum
      case ('-rrhoav')             !> see above in the first specification of -rrhoav
        env%properties = p_rrhoaverage
        call read_bhess_ref(env,'coord')
      case ('-avbhess')
        env%thermo%avbhess = .true. !> use bhess in rrhoav for all structures (expensive)
      case ('-avchess')
        env%thermo%constrhess = .true.   !> apply constraints during rrhoav routine
      case ('-printpop')
        env%thermo%printpop = .true. !> print a file with free energy pop. at different T
      case ('-noref')              !> dont use a bhess reference
        env%emtd%bhess = .false.
      case ('-ref')
        env%emtd%bhess = .true.
        inquire (file=trim(arg(i+1)),exist=ex)
        if (ex) then
          call read_bhess_ref(env,trim(arg(i+1)))
        end if
      case ('-pcap')
        read (arg(i+1),*,iostat=io) j
        if (io == 0.and.(index(arg(i+1),'-') .eq. 0)) then
          env%thermo%pcap = j
        end if
      case ('-ptot')
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0.and.(index(arg(i+1),'-') .eq. 0)) then
          if (rdum > 1.0d0) rdum = 1.0d0
          env%thermo%ptot = rdum
        end if
      case ('-ithr')
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0) then
          if (rdum > 0.0d0) rdum = 0.0
          env%thermo%ithr = rdum
        end if
      case ('-rotorcut','-sthr')
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0.and.(index(arg(i+1),'-') .eq. 0)) then
          if (rdum < 0.0d0) rdum = 0.0d0
          env%thermo%sthr = rdum
        end if
      case ('-fscal')
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0.and.(index(arg(i+1),'-') .eq. 0)) then
          env%thermo%fscal = rdum
        end if
      case ('-trange')    !> provide a range of temperatures (min max step) for entropy evaluation
        read (arg(i+1),*,iostat=io) rdum
        if (io == 0.and.(index(arg(i+1),'-') .eq. 0)) then
          env%thermo%trange(1) = rdum  !> T start
        end if
        read (arg(i+2),*,iostat=io) rdum
        if (io == 0.and.(index(arg(i+2),'-') .eq. 0)) then
          env%thermo%trange(2) = rdum  !> T stop (approx.)
        end if
        read (arg(i+3),*,iostat=io) rdum
        if (io == 0.and.(index(arg(i+3),'-') .eq. 0)) then
          env%thermo%trange(3) = rdum  !> T step
        end if
      case ('-tread')   !> read a file with temperatures (one per line) for entropy evaluation
        ctmp = trim(arg(i+1))
        inquire (file=ctmp,exist=ex)
        if (ex) then
          call env%thermo%read_temps(ctmp)
        end if
!=========================================================================================!
!-------- QCG-Related flags
!=========================================================================================!
      case ('-nopreopt')
        env%nopreopt = .true.
        env%qcg_flag = .true.
      case ('-xtbiff')
        env%use_xtbiff = .true.
      case ('-grow')
        env%qcg_runtype = 0
        env%qcg_flag = .true.
      case ('-ensemble')
        env%qcg_runtype = 1
        env%qcg_flag = .true.
      case ('-esolv')
        env%qcg_runtype = 2
        env%qcg_flag = .true.
      case ('-gsolv')
        env%qcg_runtype = 3
        env%qcg_flag = .true.
      case ('-nsolv')
        env%qcg_flag = .true.
        call readl(arg(i+1),xx,j)
        env%nsolv = NINT(xx(1))
      case ('-maxsolv')
        env%qcg_flag = .true.
        call readl(arg(i+1),xx,j)
        env%max_solv = NINT(xx(1))
      case ('-normdock')
        env%docking_qcg_flag = ''
      case ('-fin_opt_gfn2')
        env%final_gfn2_opt = .true.
      case ('-directed') !> specify the directed list
        env%qcg_flag = .true.
        ctmp = trim(arg(i + 1))
        if (ctmp(1:1) .ne. '-') then
          env%directed_file = trim(ctmp)
          write (*,'(2x,a,1x,a)') trim(argument)//' :',trim(ctmp)
        end if
      case ('-nclus')
        env%qcg_flag = .true.
        call readl(arg(i+1),xx,j)
        env%nqcgclust = NINT(xx(1))
        env%user_nclust = .true.
      case ('-freqscal')
        env%qcg_flag = .true.
        call readl(arg(i+1),xx,j)
        env%freq_scal = (xx(1))
      case ('-qcgmtd')
        env%ensemble_method = -1
        env%qcg_flag = .true.
      case ('-ncimtd')
        env%ensemble_method = 0
        env%qcg_flag = .true.
      case ('-md')
        env%ensemble_method = 1
        env%qcg_flag = .true.
        if (.not.env%user_enslvl) then
          env%ensemble_opt = '--gfn2'
        end if
      case ('-mtd')
        env%ensemble_method = 2
        env%qcg_flag = .true.
        if (.not.env%user_enslvl) then
          env%ensemble_opt = '--gfn2'
        end if
      case ('-samerand')
        env%sameRandomNumber = .true.
        env%qcg_flag = .true.
      case ('-nocff')
        env%cff = .false.
        env%qcg_flag = .true.
      case ('-enslvl')
        ctmp = arg(i+1)
        env%user_enslvl = .true.
        env%qcg_flag = .true.
        if (arg(i+1) == 'gfn') then
          dtmp = trim(arg(i+2))
          ctmp = trim(ctmp)//dtmp
        end if
        select case (ctmp)
        case ('gfn1')
          env%ensemble_opt = '--gfn1'
          write (*,'(2x, a)') 'Use of GFN1-xTB for ensemble search requested.'
        case ('gfn2')
          env%ensemble_opt = '--gfn2'
          write (*,'(2x, a)') 'Use of GFN2-xTB for ensemble search requested.'
        case ('gfn0')
          env%ensemble_opt = '--gfn0'
          write (*,'(2x, a)') 'Use of GFN0-xTB for ensemble search requested.'
        case ('gff','gfnff')
          env%ensemble_opt = '--gff'
          write (*,'(2x, a)') 'Use of GFN-FF for ensemble search requested.'
        end select

      case ('-freqlvl')
        ctmp = arg(i+1)
        env%qcg_flag = .true.
        if (arg(i+1) == 'gfn') then
          dtmp = trim(arg(i+2))
          ctmp = trim(ctmp)//dtmp
        end if
        select case (ctmp)
        case ('gfn1')
          env%freqver = '--gfn1'
          write (*,'(2x, a)') 'Use of GFN1-xTB for frequency computation requested.'
        case ('gfn2')
          env%freqver = '--gfn2'
          write (*,'(2x, a)') 'Use of GFN2-xTB for frequency computation requested.'
        case ('gfn0')
          env%freqver = '--gfn0'
          write (*,'(2x, a)') 'Use of GFN0-xTB for frequency computation requested.'
        case ('gff','gfnff')
          env%freqver = '--gff'
          write (*,'(2x, a)') 'Use of GFN-FF for frequency computation requested.'
        end select
!========================================================================================!
!-------- PRINCIPAL COMPONENT analysis and CLUSTERING flags
!========================================================================================!
      case ('-cluster')
        write (*,'(2x,a,'' : ensemble clustering'')') trim(arg(i))
        if (env%properties == p_propcalc) then
          !>--- for standalone use
          env%properties = p_cluster
        elseif (env%confgo.and.env%properties == p_cregen) then
          !>--- as extension for CREGEN
          env%cluster = .true.
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
            env%nclust = j
          else
            env%nclust = 0
            if ((index(arg(i+1),'-') .eq. 0)) then
              ctmp = trim(arg(i+1))
              select case (ctmp)
              case ('loose')
                env%clustlev = -1
                write (*,'(2x,a,'' loose : using loose clustering setting'')') trim(arg(i))
              case ('normal')
                env%clustlev = 0
                write (*,'(2x,a,'' normal : using normal clustering setting'')') trim(arg(i))
              case ('tight')
                env%clustlev = 1
                write (*,'(2x,a,'' tight : using tight clustering setting'')') trim(arg(i))
              case ('vtight','verytight')
                env%clustlev = 2
                write (*,'(2x,a,'' vtight : using very tight clustering setting'')') trim(arg(i))
              case ('incremental','incr')
                env%clustlev = 10
                write (*,'(2x,a,'' incremental : using incremental clustering settings'')') trim(arg(i))
              case ('tightincremental','tightincr')
                env%clustlev = 11
                write (*,'(2x,a,'' tightincremental : using incremental clustering settings'')') &
                &    trim(arg(i))
              case ('vtightincremental','vtightincr')
                env%clustlev = 12
                write (*,'(2x,a,'' vtightincremental : using incremental clustering settings'')') &
                 & trim(arg(i))
              end select
            end if
          end if
        end if
      case ('-pccap')
        if (i+1 .le. nra) then !> second argument is the max. number of PCs
          read (arg(i+1),*,iostat=io) j
          if (io == 0.and.(index(arg(i+1),'-') .eq. 0)) then
            env%pccap = j
          end if
        end if
      case ('-nopcmin')
        env%pcmin = 0.0d0
      case ('-pctype','-pctyp')
        if (i+1 .le. nra) then
          ctmp = trim(arg(i+1))
          if (ctmp(1:1) .ne. '-') then
            env%pcmeasure = ctmp
          end if
        end if
      case ('-pcaex','-pcaexclude')
        if (i+1 .le. nra) then
          ctmp = trim(arg(i+1))
          if (ctmp(1:1) .ne. '-') then
            env%atlist = ctmp
            env%pcaexclude = .true.
          end if
        end if
!========================================================================================!
!---------- PROPERTY MODE
!========================================================================================!
      case ('-prop')
!>----------------------------------------------------------------
!> NOTE: These flags are outdated and using them is discouraged!
!>----------------------------------------------------------------
        if ((env%properties == p_none.or.    &
        &  env%properties == p_propcalc)) then         !property selection
          ctmp = trim(arg(i+1))
          PROPARG:select case(ctmp)
          case ('hess')                  !hessian calculation to free energies for all conformers
          env%properties2 = 1
          case ('ohess')                 !optimization+hessian calculation
          env%properties2 = 10
          case ('autoir','autoIR')       !automated IR averaging for populated (-pthr) conformers
          env%properties2 = 2
          case ('b973c')                 !B97-3c optimization (xtb driver for ancopt)
          env%properties2 = 3
          case ('b973cIR')               !B97-3c optimization + IR spectra average
          env%properties2 = 4
          case ('dft')                   !DFT (custom) job, read from dftrc
          env%properties2 = 100
          case ('dftOPT')                !DFT (custom) optimization (xtb driver for ancopt)
          env%properties2 = 5
          case ('dftIR')                 !DFT (custom) optimization + IR spectra average
          env%properties2 = 6
          case ('dftSP')                 !DFT (custom) singlepoint
          env%properties2 = 7
          env%harcutpthr = 0.75
          case ('dftFREQ')               !DFT (custom) optimization + frequencies
          env%properties2 = 8
          case ('reopt')                 !reoptimize only conformers at vtight level
          env%properties2 = 20
          case ('TEST')                  !testSTUFF
          env%properties2 = -9999
          case ('singlepoint','sp')      !singlepoint calculation and ensemble sorting
          env%properties2 = 999
          env%pclean = .true.
          case ('dipole')      !singlepoint calculation and dipole grepping
          env%properties2 = 998
          env%pclean = .true.
          case default
          env%properties2 = 0
          end select PROPARG
          if (env%properties2 .ne. 0) then
            call env%addjob(env%properties2)
          end if
        end if
      case ('-dftrc')                            !provide dft-rc file (including path)
        atmp = ''
        if (nra .ge. (i+1)) atmp = adjustl(arg(i+1))
        if ((atmp(1:1) /= '-').and.(len_trim(atmp) .ge. 1)) then
          env%dftrcfile = trim(atmp)
        end if
      case ('-hardcut')                          !cut DFT populations hard
        env%hardcutDFT = .true.
      case ('-pclean')                           !cleanup option for property mode, i.e., remove PROP/
        env%pclean = .true.
!========================================================================================!
      case ('-scratch')
        !use a scratch directory to perform the calculation in
        env%scratch = .true.
        atmp = ''
        if (nra .ge. (i+1)) atmp = adjustl(arg(i+1))
        if ((atmp(1:1) /= '-').and.(len_trim(atmp) .ge. 1)) then
          env%scratchdir = trim(atmp)
        end if
      case ('-keepscratch')
        env%keepScratch = .true.
      case default
        continue
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

!>----- additional checks and settings
  if (env%crestver .eq. crest_solv) bondconst = .false.

  if (env%qcg_flag.and.env%crestver .ne. crest_solv) then
    error stop 'At least one flag is only usable for QCG runtype. Exit.'
  end if

  if (env%autozsort.and.env%crestver .eq. crest_solv) then
    error stop 'Z sorting of the input is unavailable for -qcg runtyp.'
  end if

!>--- automatic wall potential for the LEGACY version
  if (env%NCI.or.env%wallsetup .and. env%legacy) then
    call wallpot(env)
    if (env%wallsetup) then
      write (*,'(2x,a)') 'Automatically generated ellipsoide potential:'
    else
      write (*,'(2x,a)') 'Automatically generated ellipsoide potential for NCI mode:'
    end if
    call write_cts_NCI_pr(stdout,env%cts)
    write (*,*)
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
  if (env%ensemble_opt == '--gfn2' .or. env%gfnver == '--gfn2') &
          & env%final_gfn2_opt = .false. !Prevent additional opt.

  if (env%useqmdff) then
    env%autozsort = .false.
  end if

  if (.not.env%preopt .and. env%crestver.ne.crest_trialopt) then
    if (allocated(env%ref%topo)) deallocate (env%ref%topo)
  end if

!>-- turn off niceprint if we are not writing to terminal
  if(env%niceprint)then
    env%niceprint = myisatty(output_unit)
  endif

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
  if(env%autozsort)then
    env%properties = p_zsort
  endif 

!>--- for legacy runtypes, check if xtb is present
  if(env%legacy)then
    call checkprog_silent(env%ProgName,.true.,iostat=io)
    if(io /= 0 ) error stop
  endif

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> FALLBACK setup of new calculator (important CREST >3.0 handling things)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
  if (.not.env%legacy.and.env%calc%ncalculations == 0) then
    write (stdout,'(/,a)',advance='no') '> Setting up backup calculator ...'
    flush (stdout)
    call env2calc_setup(env)
    write(stdout,*) 'done.'
    call env%calc%info(stdout)
  end if
!>--- pass on opt-level to new calculator
  if(.not.env%legacy)then
     env%calc%optlev = nint(env%optlev)
  endif

!>--- ONIOM setup from toml file
  if (allocated(env%ONIOM_toml))then
    allocate(env%calc%ONIOM)
    call ONIOM_read_toml(env%ONIOM_toml,env%nat,env%ref%at,env%ref%xyz,env%calc%ONIOM)     
    call env%calc%ONIOMexpand()
  endif

!>--- important printouts
  if( .not.env%legacy)then
    call print_frozen(env)
  endif

  return
end subroutine parseflags

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
    write (*,'(/,1x,a,a,a)') '<',trim(env%constraints),'> file present.'
    env%cts%used = .true.
  else
    env%cts%used = .false.
    return
  end if 

!>--- read the data
  call read_constrainbuffer(env%constraints,env%cts)
  call sort_constraints(env%cts)
  write (*,*) 'content of the constraining file (sorted):'
  if (env%cts%ndim .gt. 20) then
    write (*,'(1x,a)') '<skipped due to length of constraining file>'
  else
    do i = 1,env%cts%ndim
      if (trim(env%cts%sett(i)) .ne. '') then
        write (*,'(''>'',1x,a)') trim(env%cts%sett(i))
      end if
    end do
  end if
  if(.not.env%legacy)then
    call parse_xtbinputfile(env,env%constraints)
  endif

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
      write (*,'(2x,a)') trim(adjustl(btmp))
      write (*,'(5x,a,i0)') '# of atoms considered for RMSDs:',env%rednat
      env%includeRMSD = atlist !includeRMSD contains only the atoms that are included in RMSD
    end if
    if (index(atmp,'atomlist-') .ne. 0) then
      create = .true.
      atomlistused = .true.
      dg = atmp
      call split_set_args(dg,argument)
      call parse_atlist_new(trim(argument),j,env%nat,env%ref%at,atlist)
      env%rednat = env%nat-j
      write (*,'(2x,a)') trim(adjustl(btmp))
      write (*,'(3x,a,i0)') '# of atoms considered for RMSDs:',env%rednat
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
          write (*,'(2x,a)') trim(adjustl(btmp))
          write (*,'(5x,a,i0)') '# of atoms considered for RMSDs:',env%rednat
          env%includeRMSD = atlist !includeRMSD contains only the atoms that are included in RMSD
        end if
      end do
    end if
    if (index(btmp,'reference=') .ne. 0) then
      call rdarg(btmp,'reference=',env%fixfile)
      write (*,'(1x,a,1x,a)') 'fix file:',trim(env%fixfile)
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

        write (*,'(/,2x,a)') 'Automatically generated ellipsoide potential overwritten by:'
        call write_cts_NCI(6,env%cts)
        write (*,*)

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
  integer :: i

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
    call mol%write('coord')
    call mol%write('crest_input_copy.xyz')
    call mol%deallocate()
    inputfile = trim(arg)
  else
    inputfile = 'coord'
  end if
  if(.not.allocated(env%inputcoords)) env%inputcoords = inputfile

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
!>-- shift to CMA and align according to rot.const.
  if (env%crestver /= crest_solv .and. env%crestver /= crest_sp &
  & .and. env%crestver /= crest_optimize) call axis(mol%nat,mol%at,mol%xyz)
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
  env%ptb%nfrag = zmol%nfrag
  call zmol%deallocate()

  return
end subroutine inputcoords

!========================================================================================!
!========================================================================================!
