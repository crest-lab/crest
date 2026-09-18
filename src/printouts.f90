!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2024 Philipp Pracht
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

!=========================================================================================!
!=========================================================================================!
!> CREST PRINTOUT ROUTINES
!=========================================================================================!
!=========================================================================================!
subroutine confscript_head(vers)
!*******************************
!* Print program header section
!*******************************
  implicit none
  logical,intent(in) :: vers
  logical :: niceprint
  include 'crest_metadata.fh' !> this file should be created by meson or CMake

  niceprint = .true.
  call box3(version,date,commit,author)
  write (*,*)

  if (vers) then
    write (*,*) "crest ",trim(version)
    stop
  end if

  write (*,'(3x,''Cite work conducted with this code as'')')
  write (*,'(/,3x,''• P.Pracht, F.Bohle, S.Grimme, PCCP, 2020, 22, 7169-7192.'')')
  write (*,'(  3x,''• S.Grimme, JCTC, 2019, 15, 2847-2862.'')')
  write (*,'(  3x,''• P.Pracht, S.Grimme, C.Bannwarth, F.Bohle, S.Ehlert,'')')
  write (*,'(  3x,''  G.Feldmann, J.Gorges, M.Müller, T.Neudecker, C.Plett,'')')
  write (*,'(  3x,''  S.Spicher, P.Steinbach, P.Wesołowski, F.Zeller,'')')
  write (*,'(  3x,''  J. Chem. Phys., 2024, 160, 114110.'')')
  write (*,'(/,3x,''for works involving QCG cite'')')
  write (*,'(/,3x,''• S.Spicher, C.Plett, P.Pracht, A.Hansen, S.Grimme,'')')
  write (*,'(  3x,''  JCTC, 2022, 18 (5), 3174-3189.'')')
  write (*,'(  3x,''• C.Plett, S. Grimme,'')')
  write (*,'(  3x,''  Angew. Chem. Int. Ed. 2023, 62, e202214477.'')')
  write (*,'(/,3x,''for works involving MECP screening cite'')')
  write (*,'(/,3x,''• P.Pracht, C.Bannwarth, JCTC, 2022, 18 (10), 6370-6385.'')')
  write (*,'(/,3x,''for works involving the TTConf reimplementaiton cite'')')
  write (*,'(/,3x,''• C.Zurek, et al., JCTC, 2025, 21, 1459-1475.'')')
  write (*,*)

  write (*,'(3x,a)') 'Original code'
  write (*,'(4x,a)') ' P.Pracht, S.Grimme, Universität Bonn, MCTC'
  write (*,'(3x,a)') 'with help from (alphabetical order):'
  write (*,'(4x,a)') ' C.Bannwarth, F.Bohle, S.Ehlert, G.Feldmann, J.Gorges,'
  write (*,'(4x,a)') ' S.Grimme, C.Plett, P.Pracht, L.Rindt, S.Spicher,'
  write (*,'(4x,a)') ' P.Steinbach, P.Wesolowski, F.Zeller'
  write (*,*)

  write (*,'(3x,a)') 'Online documentation is available at'
  write (*,'(3x,a)') 'https://crest-lab.github.io/crest-docs/'
  write (*,*)

  call disclaimer()
end subroutine confscript_head

subroutine box3(version,date,commit,author)
!***************
!* Print banner
!***************
  implicit none
  character(len=*) :: version
  character(len=*) :: date
  character(len=*) :: commit
  character(len=*) :: author
  character(len=200) :: logo(13)
  character(len=200) :: info(2)
  integer,parameter :: pad_left = 7
  integer :: i,lcount
  write (*,*)
  !write (logo(1),'(''╔════════════════════════════════════════════╗'')')
  !write (logo(2),'(''║            ___ ___ ___ ___ _____           ║'')')
  !write (logo(3),'(''║           / __| _ \ __/ __|_   _|          ║'')')
  !write (logo(4),'(''║          | (__|   / _|\__ \ | |            ║'')')
  !write (logo(5),'(''║           \___|_|_\___|___/ |_|            ║'')')
  !write (logo(6),'(''║                                            ║'')')
  !write (logo(7),'(''║  Conformer-Rotamer Ensemble Sampling Tool  ║'')')
  !write (logo(8),'(''║          based on the xTB methods          ║'')')
  !write (logo(9),'(''║                                            ║'')')
  !write (logo(10),'("╚════════════════════════════════════════════╝")')

  write (logo(1),'(''╔════════════════════════════════════════════════╗'')')
  write (logo(2),'(''║                                                ║'')')
  write (logo(3),'(''║     ██████╗██████╗ ███████╗███████╗████████╗   ║'')')
  write (logo(4),'(''║    ██╔════╝██╔══██╗██╔════╝██╔════╝╚══██╔══╝   ║'')')
  write (logo(5),'(''║    ██║     ██████╔╝█████╗  ███████╗   ██║      ║'')')
  write (logo(6),'(''║    ██║     ██╔══██╗██╔══╝  ╚════██║   ██║      ║'')')
  write (logo(7),'(''║    ╚██████╗██║  ██║███████╗███████║   ██║      ║'')')
  write (logo(8),'(''║     ╚═════╝╚═╝  ╚═╝╚══════╝╚══════╝   ╚═╝      ║'')')
  write (logo(9),'(''║                                                ║'')')
  write (logo(10),'(''║    Conformer-Rotamer Ensemble Sampling Tool    ║'')')
  write (logo(11),'(''║            based on the xTB methods            ║'')')
  write (logo(12),'(''║                                                ║'')')
  write (logo(13),'(''╚════════════════════════════════════════════════╝'')')

  do i = 1,13
    write (*,'(a,a)') repeat(" ",pad_left),trim(logo(i))
  end do
  write (*,'(a,'' Version '',a,'', '',a)') repeat(" ",pad_left),trim(version),trim(date)
  if (author(1:2) .eq. "'@") then
    write (*,'(a," commit (",a,") compiled by ",a)') repeat(" ",pad_left),commit,"'usr"//author(2:)
  else
    write (*,'(a," commit (",a,") compiled by ",a)') repeat(" ",pad_left),commit,author
  end if
end subroutine box3

subroutine disclaimer

  write (*,'(3x,a)') 'This program is distributed in the hope that it will be useful,'
  write (*,'(3x,a)') 'but WITHOUT ANY WARRANTY; without even the implied warranty of'
  write (*,'(3x,a)') 'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the'
  write (*,'(3x,a)') 'GNU Lesser General Public License (LGPL) for more details.'

end subroutine disclaimer

!=========================================================================================!

subroutine help_section(title)
  !*************************************
  !* Print a colored section header.   *
  !*************************************
  use iomod,only:colorify
  use crest_parameters,only:stdout
  implicit none
  character(len=*),intent(in) :: title
  integer :: n
  n = len_trim(title)
  write (stdout,'(/,1x,a)') colorify(trim(title),'yellow')
  write (stdout,'(1x,a)') colorify(repeat('─',n),'yellow')
end subroutine help_section

subroutine help_opt(flag,fw,desc)
  !*********************************************************************
  !* Print one colored flag + description, padding the flag column to  *
  !* fw characters wide so descriptions align.                         *
  !*********************************************************************
  use iomod,only:colorify
  use crest_parameters,only:stdout
  implicit none
  character(len=*),intent(in) :: flag,desc
  integer,intent(in) :: fw
  integer :: fl,pad
  fl = len_trim(flag)
  pad = max(fw-fl,1)
  write (stdout,'(a,a,a,a,a)') '   ',colorify(trim(flag),'green'), &
    & repeat(' ',pad),' : ',trim(desc)
end subroutine help_opt

!=========================================================================================!

subroutine confscript_help()
  use iomod,only:colorify
  use crest_parameters,only:stdout
  implicit none

  write (stdout,'(/,1x,a)') colorify(repeat('─',76),'gold')
  write (stdout,'(1x,a,a)') colorify('Usage:','yellow'),'  crest [INPUT] [OPTIONS]'
  write (stdout,'(1x,a)') colorify(repeat('─',76),'gold')
  write (stdout,*)
  write (stdout,'(1x,a)') 'The '//colorify('[INPUT]','blue')//' argument CAN be a coordinate file in the'
  write (stdout,'(1x,a)') 'TM (coord, Bohr) or Xmol (*.xyz, Ang.) format.'
  write (stdout,'(1x,a)') 'If no such file is present as the first argument, crest will'
  write (stdout,'(1x,a)') 'automatically search for a file called "coord" in the TM format.'
  write (stdout,*)
  write (stdout,'(1x,a)') colorify('Versions >3.0 allow specifying detailed input instructions via','green')
  write (stdout,'(1x,a)') colorify('input files in the TOML format.','green')
  write (stdout,'(1x,a)') colorify('*.toml files can be ','green')//colorify(' [INPUT]','blue')// &
    colorify(' or specified via "--input <file>"','green')
  write (stdout,*)
  call confscript_morehelp2()
  stop '   [-h] displayed. exit.'
end subroutine confscript_help

subroutine confscript_morehelp(flag)
  use iomod,only:colorify
  use crest_parameters,only:stdout
  implicit none
  character(len=*),intent(in) :: flag
  integer :: fw

  write (stdout,'(/,1x,a)') colorify(repeat('─',76),'gold')
  write (stdout,*)
  select case (flag)

    ! ── General / technical ──────────────────────────────────────────────
  case default
    fw = 16
    call help_section('Run modes:')
    call help_opt('-sp',fw,'Single-point energy calculation')
    call help_opt('-opt/-optimize',fw,'Geometry optimization')
    call help_opt('-hess/-numhess',fw,'Numerical Hessian / vibrational frequencies')
    call help_opt('-md/-dynamics',fw,'Molecular dynamics simulation')
    call help_opt('-v3/-imtdgc',fw,'iMTD-GC conformational search  (see '//colorify('--help conf','gold')//')')
    call help_opt('-v4/-entropy',fw,'Entropy/free-energy sampling  (see '//colorify('--help conf','gold')//')')
    call help_opt('-mdopt',fw,'Ensemble optimization (no sorting)')
    call help_opt('-screen',fw,'Ensemble screening')
    call help_opt('-protonate',fw,'Protonation site search')
    call help_opt('-deprotonate',fw,'Deprotonation site search')
    call help_opt('-tautomerize',fw,'Tautomer generation')
    call help_opt('-qcg',fw,'Quantum Cluster Growth workflows  (see '//colorify('--help qcg','gold')//')')
    call help_opt('-msreact',fw,'MS fragment generator  (see '//colorify('--help msreact','gold')//')')
    !call help_opt('-bh/-GMIN',fw,'Basin-hopping global optimization')
    call help_opt('-sort/-cregen',fw,'Ensemble sorting, comparison, clustering  (see '//colorify('--help compare','gold')//')')
    call help_opt('-ttcon',fw,'TTConf conformational sampling (see '//colorify('--help conf','gold')//')')
    write (stdout,*)
    fw = 22
    call help_section('Method selection:')
    call help_opt('-gfn2',fw,'Use GFN2-xTB  [default]')
    call help_opt('-gfn1',fw,'Use GFN1-xTB')
    call help_opt('-gfn0',fw,'Use GFN0-xTB')
    call help_opt('-gff/-gfnff',fw,'Use GFN-FF  (bond constraints applied automatically)')
    call help_opt('-gxtb',fw,'Use g-xTB  (requires special build)')
    call help_opt('-orca <tmpl> <exe>',fw, &
    & 'Use an ORCA driver: <tmpl> is an ORCA input template (EnGrad jobs only), <exe> the ORCA executable')
    call help_opt('-mlip <backend>',fw, &
    & 'Use an fmlip-relay ML potential: <backend> is uma (omol task) or maceoff (MACE-OFF23 medium)')
    call help_opt('-A//B  -A/sp/B',fw,'sampling at B; inline SP re-ranking at A  (e.g. -gfn2//gfnff)')
    call help_opt('-A/opt/B',fw,'sampling at B; inline geometry refinement at A')
    call help_opt('-A@B',fw,'sampling at B; post-search re-optimization of ensemble at A')
    write (stdout,'(9x,a)') 'A, B in: gfn0, gfn1, gfn2, gxtb, gfnff'
    call help_opt('-rerank <method>',fw,'Post-search SP re-ranking of conformer ensemble at <method>')
    call help_opt('-reopt <method>',fw,'Post-search geometry optimization of conformer ensemble at <method>')
    call help_opt('-finalhess',fw,'Post-search Hessians on final ensemble, re-sort by Gibbs free energy')
    call help_opt('-refine <method>',fw,'Post-process conformers at a higher level (single step)')
    call help_opt('-optlev <level>',fw,'Optimization convergence level for ALL semiempirical calculations')
    write (stdout,'(9x,a)') '<level> = crude, vloose, loose, normal, tight, vtight, extreme'
    call help_opt('-dscal [<factor>]',fw,'Scale dispersion energy in MD/MTD simulations')
    write (stdout,*)
    fw = 22
    call help_section('Molecular system:')
    call help_opt('-T <int>',fw,'Number of CPU threads (or read from OMP_NUM_THREADS)')
    call help_opt('-chrg <int>',fw,"Molecular charge")
    call help_opt('-uhf <int>',fw,'Unpaired electrons (N_alpha - N_beta)')
    call help_opt('-mult <int>',fw,'Spin multiplicity (2S+1 = uhf+1)')
    call help_opt('-spin-polarized',fw,'Spin-polarized treatment (open-shell, tblite only)')
    call help_opt('-g/-gbsa <solvent>',fw,'GBSA implicit solvation')
    call help_opt('-alpb <solvent>',fw,'ALPB implicit solvation')
    call help_opt('-efield <x> <y> <z>',fw,'External electric field in V/Ang along x, y, z')
    call help_opt('-charges [<file>]',fw,'Read atomic partial charges from file  [default: "charges"]')
    write (stdout,*)
    fw = 22
    call help_section('Technical:')
    call help_opt('--input <file>',fw,'Specify TOML input file with detailed settings')
    call help_opt('-xnam <bin>',fw,'Path to the xtb executable (when using xtb as backend)')
    call help_opt('-noopt',fw,'Skip pre-optimization of the input structure')
    call help_opt('-niceprint',fw,'Show progress bar during optimizations')
    call help_opt('-dry',fw,'Parse args, print resolved settings, then exit')
    call help_opt('-legacy',fw,'Force CREST < 3.0 behavior')
    write (stdout,*)
    fw = 22
    call help_section('Constraints (applied to ALL calculations):')
    call help_opt('-cinp <file>',fw,'Read constraints file  (xtb format; formerly ".constrains")')
    call help_opt('-freeze <atoms>',fw,'Freeze atoms in optimizations/MD  (indices, ranges, or element symbols)')
    call help_opt('-cbonds [<fc>]',fw,'Constrain all bonds globally  (set up from topology)')
    call help_opt('-cbonds_md [<fc>]',fw,'Constrain all bonds during MDs/MTDs only')
    call help_opt('-nocbonds',fw,'Disable automatic bond constraints')
    call help_opt('-fc <float>',fw,'Global force constant for bond constraints')
    write (stdout,*)

    ! ── Ensemble comparison / CREGEN ─────────────────────────────────────
  case ('compare','cregen')
    fw = 22
    call help_section('Ensemble sorting (CREGEN):')
    call help_opt('-cregen [file]',fw,'Run CREGEN standalone to sort an ensemble file')
    call help_opt('-sort <file>',fw,'Ensemble sorting (default: iRMSD-based conformer assignment)')
    call help_opt('-sort <file> cregen',fw,'CREGEN sorting via --sort')
    write (stdout,*)
    call help_section('CREGEN thresholds:')
    call help_opt('-ewin <real>',fw,'Energy window in kcal/mol  [default: 6.0]')
    call help_opt('-rthr <real>',fw,'RMSD threshold in Ang  [default: 0.125]')
    call help_opt('-ethr <real>',fw,'Energy threshold in kcal/mol  [default: 0.05]')
    call help_opt('-bthr <real>',fw,'Rotational constant threshold  [default: 0.01 = 1%]')
    call help_opt('-pthr <real>',fw,'Boltzmann population threshold (0-1)  [default: 0.05]')
    call help_opt('-temp <real>',fw,'Boltzmann temperature in K  [default: 298.15]')
    write (stdout,*)
    call help_section('CREGEN algorithm options:')
    call help_opt('-topo/-notopo',fw,'Enable/disable topology change check')
    call help_opt('-ezcheck',fw,'Enable E/Z double-bond isomer check')
    call help_opt('-heavy',fw,'Use heavy-atom-only RMSD')
    call help_opt('-allrot',fw,'Use all three rotational constants (A, B, C)')
    call help_opt('-eqv/-nmr',fw,'NMR nuclear equivalence analysis (requires rotamers)')
    write (stdout,*)
    call help_section('CREGEN output:')
    call help_opt('-prsc',fw,'Write scoord.* file for each conformer')
    call help_opt('-nowr',fw,"Skip writing the sorted ensemble file")
    call help_opt('-osdf',fw,'Also write output ensemble in SDF format')
    write (stdout,*)
    call help_section('Structural comparison tools:')
    call help_opt('-rmsd <f1> <f2>',fw,'RMSD between two structures')
    call help_opt('-rmsdheavy <f1> <f2>',fw,'Heavy-atom RMSD between two structures')
    call help_opt('-irmsd <f1> <f2>',fw,'iRMSD (iterative RMSD with atom permutation)')
    call help_opt('-compare <f1> <f2>',fw,'Compare two ensembles for structural overlap')
    write (stdout,'(11x,a)') colorify('-maxcomp <int>','green')//' : max structures per ensemble  [default: 100]'
    write (stdout,*)
    call help_section('PCA / k-Means clustering (CCEGEN):')
    call help_opt('-sort <file> --cluster',fw,'Standalone PCA clustering of an ensemble')
    call help_opt('-cluster [<N>]',fw,'PCA + k-Means clustering, optionally fix the cluster count')
    call help_opt('-cluster <level>',fw,'Clustering level: loose, normal, tight, vtight')
    call help_opt('-cluster incremental',fw,'Test cluster sizes incrementally (large ensembles)')
    call help_opt('-pctype <type>',fw,'PCA descriptor: dihedral (default), zmat, cartesian, cma')
    call help_opt('-pccap <int>',fw,'Cap the number of principal components  [default: 100]')
    call help_opt('-nopcmin',fw,'Include all PCs (disable minimum eigenvalue cutoff)')
    call help_opt('-pcaexclude <atoms>',fw,'Exclude specified atoms from the PCA')
    write (stdout,*)

    ! ── Conformer search / sampling ──────────────────────────────────────
  case ('conf','sampling')
    fw = 20
    call help_section('Conformer search algorithms:')
    call help_opt('-v3/-imtdgc',fw,'iMTD-GC (iterative MTD-GC)  [default]')
    call help_opt('-v4',fw,'iMTD-sMTD (entropy-focused search)')
    call help_opt('-entropy',fw,'Same as -v4, specialized for conformational entropy')
    call help_opt('-ttconf [<preset>]',fw,'TTConf conformational sampling, where <preset> is one of [fast,normal,accurate] or [bruteforce]')

    write (stdout,*)
    call help_section('TTConf-light parameters:')
    call help_opt('-ttrank <int>',fw,'TT rank r  [default: 3]')
    call help_opt('-ttsweeps <int>',fw,'Number of TT-cross sweeps s  [default: 8]')
    call help_opt('-ttgrid <int>',fw,'Dihedral grid points (360/n degree steps)  [default: 6]')
    call help_opt('-ttninit <int>',fw,'Number of random initial tail seeds  [default: 3]')
    call help_opt('-ttewin <float>',fw,'Conformer energy window in kcal/mol  [default: 6.0]')
    call help_opt('-ttkt <float>',fw,'Maxvol weight temperature in kcal/mol  [default: 6.0]')
    call help_opt('-ttseed <int>',fw,'Fix the RNG seed (reproducible run)')
    call help_opt('-ttsp',fw,'Evaluate singlepoints only (no geometry optimization)')
    call help_opt('-ttnocache',fw,'Disable the TTConf energy cache')
    call help_opt('-ttringbonds',fw,'Treat in-ring bonds as TT variables')
    call help_opt('-ttrings',fw,'Sample ring conformations as TT sites')
    call help_opt('-ttringmethod <str>',fw,'Ring-conformation generator  [default: mtd]')

    write (stdout,*)
    call help_section('MD / MTD parameters:')
    call help_opt('-len/-mdlen <t>[x]',fw,'MD/MTD length in ps; append "x" for a scaling factor')
    call help_opt('-tstep <float>',fw,'MD timestep in fs  [default: 5 fs]')
    call help_opt('-shake <int>',fw,'SHAKE mode: 0=off, 1=X-H only, 2=all bonds  [default: 1]')
    call help_opt('-mdtemp <float>',fw,'Temperature for MTD runs in K')
    call help_opt('-tnmd <float>',fw,'Temperature for extra normal MDs in K')
    call help_opt('-mddump <int>',fw,'Trajectory dump interval in fs  [default: 100]')
    call help_opt('-vbdump <real>',fw,'Vbias dump frequency in ps  [default: 1.0]')
    call help_opt('-nmtd <int>',fw,'Number of MTD simulations per cycle')
    write (stdout,*)
    call help_section('Search control:')
    call help_opt('-cross/-nocross',fw,'Enable/disable genetic structure crossing  [cross=default]')
    call help_opt('-gcmax <int>',fw,'Max structures fed into genetic crossing')
    call help_opt('-nozs',fw,'Disable z-matrix sorting')
    call help_opt('-normmd [<n> [<T>]]',fw,'Run additional unbiased MDs on lowest conformers')
    call help_opt('-quick/-squick/-mquick',fw,'Progressively reduced search settings')
    call help_opt('-origin',fw,'Track conformer origin step  [default]')
    call help_opt('-keepdir',fw,'Keep temporary working directories')
    call help_opt('-NCI',fw,'NCI cluster mode (flat-bottom wall + specialised MTD settings)')
    call help_opt('-wscal <float>',fw,'Scale wall potential sphere radius')
    call help_opt('-hflip/-noflip',fw,'OH proton flip after MTD  [default: OFF]')
    call help_opt('-maxflip <int>',fw,'Max OH flip attempts  [default: 1000]')
    write (stdout,*)

    ! ── Thermochemistry / entropy ─────────────────────────────────────────
  case ('thermo','entropy')
    fw = 28
    call help_section('Thermostatistical options:')
    call help_opt('-trange <Tmin> <Tmax> <Tstep>',fw,'Temperature range in K for entropy output')
    write (stdout,'(9x,a)') '[default: 280-380 K in 10 K steps]'
    call help_opt('-tread <file>',fw,'Read temperatures (one per line) from file')
    call help_opt('-fscal <float>',fw,'Frequency scaling factor  [default: 1.0]')
    call help_opt('-sthr/-rotorcut <float>',fw,'Rotor cutoff in cm^-1 (free-rotor interpolation)  [default: 25.0]')
    call help_opt('-ithr <float>',fw,'Imaginary mode inversion cutoff  [default: -50.0 cm^-1]')
    call help_opt('-ptot <float>',fw,'Cumulative population threshold for msRRHO  [default: 0.9]')
    call help_opt('-pcap <int>',fw,'Max structures used in property calculations')
    call help_opt('-printpop',fw,'Print Boltzmann populations at every temperature')
    call help_opt('-avbhess',fw,'Use Boltzmann-averaged Hessian in rrhoav')
    write (stdout,*)

    ! ── QCG ──────────────────────────────────────────────────────────────
  case ('qcg')
    fw = 20
    call help_section('Quantum Cluster Growth (QCG)')
    write (stdout,'(1x,a)') 'General usage:  crest <solute> -qcg <solvent> [options]'
    write (stdout,'(1x,a)') 'Options (in addition to general / iMTD-GC options):'
    write (stdout,*)
    call help_section('Cluster growth:')
    call help_opt('-grow',fw,'Cluster generation run type')
    call help_opt('-nsolv <int>',fw,'Number of solvent molecules to add')
    call help_opt('-fixsolute',fw,'Fix the solute during growth  (auto for water)')
    call help_opt('-nofix',fw,'Do not fix the solute  (override for water)')
    call help_opt('-nopreopt',fw,'Skip pre-optimization')
    call help_opt('-xtbiff',fw,'Use xTB-IFF standalone for solvent docking')
    call help_opt('-normdock',fw,'More extensive docking during growth')
    call help_opt('-maxsolv',fw,'Convergence limit if -nsolv not given  [default: 150]')
    call help_opt('-wscal <float>',fw,'Scaling factor for outer wall potential')
    call help_opt('-samerand',fw,'Use same random seed for every xtbiff run')
    call help_opt('-directed <file>',fw,'Directed solvation at positions in <file>')
    call help_opt('-fin_opt_gfn2',fw,'Final GFN2-xTB optimization for grow and ensemble  [default]')
    call help_opt('-no_fin_opt_gfn2',fw,'Skip it, energies at the -enslvl level instead')
    write (stdout,*)
    call help_section('Ensemble generation:')
    call help_opt('-ensemble',fw,'Ensemble generation run type')
    call help_opt('-qcgmtd',fw,'NCI-MTD CREST ensemble generation  [default]')
    call help_opt('-ncimtd',fw,'NCI-MTD CREST ensemble generation  (alias)')
    call help_opt('-mtd',fw,'MTD for QCG ensemble generation')
    call help_opt('-md',fw,'Normal MD for QCG ensemble search')
    call help_opt('-enslvl [method]',fw,'Method for ensemble search (all GFN methods supported)')
    call help_opt('-clustering',fw,'Clustering for ensemble search  (qcgmtd/ncimtd only)')
    write (stdout,*)
    call help_section('Solvation free energy:')
    call help_opt('-esolv',fw,'Solvation energy  (reference cluster generation)')
    call help_opt('-gsolv',fw,'Solvation free energy  (reference cluster generation)')
    call help_opt('-nclus',fw,'Clusters for reference generation  [default: 4]')
    call help_opt('-nocff',fw,'Switch off the CFF algorithm')
    call help_opt('-freqscal',fw,'Frequency scale factor  (output only)')
    call help_opt('-freqlvl [method]',fw,'Method for frequency computation')
    call help_opt('-keepdir',fw,'Keep temporary directories')
    write (stdout,*)

    ! ── MSReact ──────────────────────────────────────────────────────────
  case ('msreact')
    fw = 22
    call help_section('Mass spectral fragment generator (msreact)')
    write (stdout,'(1x,a)') 'General usage:  crest <input> -msreact [options]'
    write (stdout,*)
    call help_opt('-msnoattrh',fw,'Deactivate H–LMO attractive potential')
    call help_opt('-msnshifts <int>',fw,'n optimizations with randomly shifted atoms  [default: 0]')
    call help_opt('-msnshifts2 <int>',fw,'Same but with bond-repulsive potential  [default: 0]')
    call help_opt('-msnbonds <int>',fw,'Max bond distance for repulsive potential  [default: 3]')
    call help_opt('-msmolbar',fw,'Deduplicate by molbar codes  (requires "molbar")')
    call help_opt('-msinchi',fw,'Deduplicate by InChI codes  (requires "obabel")')
    call help_opt('-msnfrag <int>',fw,'Number of fragments to print  (random selection)')
    call help_opt('-msiso',fw,'Print only non-dissociated structures (isomers)')
    call help_opt('-msnoiso',fw,'Print only dissociated structures')
    call help_opt('-mslargeprint',fw,'Keep all temporary files and MSDIR')
    call help_opt('-chrg <int>',fw,"Molecular charge")
    call help_opt('-ewin <real>',fw,'Energy window for fragment sorting in kcal/mol  [default: 200.0]')
    call help_opt('-msinput <file>',fw,'Read special settings from input file')
    write (stdout,*)
    fw = 22
    call help_section('msreact input file keywords:')
    call help_opt('fragdist <real>',fw,'Inter-fragment distance increase  [default: 0.0 Ang]')
    call help_opt('atomshift <real>',fw,'Random atom displacement  [default: 0.75 Ang]')
    call help_opt('distthr_attr <real>',fw,'H–LMO attraction distance cutoff  [default: 4.0 Ang]')
    call help_opt('fc_rep <real>',fw,'Repulsive potential force constant  [default: 0.5]')
    call help_opt('fc_attr <real>',fw,'H–LMO attractive force constant  [default: -0.5]')
    call help_opt('etemp <real>',fw,'Electronic temperature in xTB optimizations')
    write (stdout,*)

    ! ── Standalone tools ─────────────────────────────────────────────────
  case ('other')
    fw = 26
    call help_section('Single-structure calculations:')
    call help_opt('-sp',fw,'Single-point energy')
    call help_opt('-opt/-optimize',fw,'Geometry optimization')
    call help_opt('-hess/-numhess',fw,'Numerical Hessian / vibrational frequencies')
    call help_opt('-dynamics/-dyn',fw,'Stand-alone MD run')
    call help_opt('-thermo <file>',fw,'Thermochemistry from existing Hessian data')
    write (stdout,'(9x,a)') '(also requires "vibspectrum" in TM format)'
    write (stdout,*)
    call help_section('Ensemble tools:')
    call help_opt('-mdopt <file>',fw,'Optimize every structure in an ensemble (XYZ)')
    call help_opt('-mdsp <file>',fw,'Singlepoints for every structure in an ensemble (XYZ)')
    call help_opt('-mdhess <file>',fw,'Hessians + thermochemistry for every structure in an ensemble (XYZ)')
    call help_opt('-screen <file>',fw,'Multi-level energy screening of an ensemble')
    call help_opt('-entropy [<T>]',fw,'Conformational entropy from ensemble')
    call help_opt('-symmetries',fw,'Symmetry analysis of all structures in an ensemble')
    call help_opt('-printboltz',fw,'Print Boltzmann population weights')
    call help_opt('-splitfile <f> [i] [j]',fw,'Split ensemble into per-structure directories (SPLIT/)')
    write (stdout,'(9x,a)') 'For sorting, clustering, RMSD, and ensemble comparison tools'
    write (stdout,'(9x,a)') 'see '//colorify('--help compare','gold')
    write (stdout,*)
    call help_section('Protonation / tautomerization:')
    call help_opt('-protonate',fw,"Find a molecule's protomers  (LMO π/LP-center approach)")
    call help_opt('-deprotonate',fw,"Find a molecule's deprotomers")
    call help_opt('-tautomerize',fw,'Find prototropic tautomers  (protonation + deprotonation)')
    write (stdout,'(9x,a)') colorify('-trev','green')//'         : deprotonate first, then protonate  (reverse order)'
    write (stdout,'(9x,a)') colorify('-iter <int>','green')//'    : number of prot/deprot cycles  [default: 2]'
    write (stdout,*)
    call help_section('Miscellaneous:')
    call help_opt('-zsort',fw,'Z-matrix sorting of the input coord file')
    call help_opt('-testtopo <file>',fw,'Topology / bond connectivity analysis')
    call help_opt('-constrain <atoms>',fw,'Write example constraint file ".xcontrol.sample"')
    write (stdout,*)

    ! ── TOML input files ─────────────────────────────────────────────────
  case ('toml')
    call help_section('TOML input files')
    write (stdout,'(1x,a)') 'CREST (v3+) accepts a TOML file as a flexible alternative to CLI flags.'
    write (stdout,'(1x,a)') 'Pass it as the first argument or explicitly with --input:'
    write (stdout,*)
    write (stdout,'(3x,a)') colorify('crest structure.xyz --input settings.toml','green')
    write (stdout,'(3x,a)') colorify('crest settings.toml','green')//'  (structure path given inside the file)'
    write (stdout,*)
    call help_section('Minimal example:')
    write (stdout,'(3x,a)') colorify('input','yellow')//'   = "struc.xyz"'
    write (stdout,'(3x,a)') colorify('runtype','yellow')//' = "iMTD-GC"'
    write (stdout,'(3x,a)') colorify('threads','yellow')//' = 4'
    write (stdout,*)
    write (stdout,'(3x,a)') colorify('[calculation]','yellow')
    write (stdout,'(5x,a)') colorify('[[calculation.level]]','yellow')
    write (stdout,'(7x,a)') 'method = "gfn2"'
    write (stdout,'(7x,a)') 'chrg   = 0'
    write (stdout,'(7x,a)') 'gbsa   = "h2o"'
    write (stdout,*)
    call help_section('Key root-level settings:')
    fw = 20
    call help_opt('input / structure',fw,'Input coordinate file')
    call help_opt('runtype',fw,'Workflow to run  (e.g. "iMTD-GC", "optimize", "md", "singlepoint")')
    call help_opt('threads',fw,'Number of CPU threads')
    call help_opt('preopt',fw,'Pre-optimize input structure  (true/false)')
    call help_opt('constraints',fw,'Path to an xtb-format constraint file')
    write (stdout,*)
    call help_section('Main blocks:')
    ! ── padding = 30 - visible_len, so ' — ' aligns at column 30 ──
    write (stdout,'(3x,a,a)') colorify('[calculation]','yellow'), &
      & repeat(' ',17)//' — method, charge, solvent, …'
    write (stdout,'(3x,a,a)') colorify('  [[calculation.level]]','yellow'), &
      & repeat(' ',7)//' — one or more calculation levels'
    write (stdout,'(3x,a,a)') colorify('  [[calculation.constraint]]','yellow'), &
      & repeat(' ',2)//' — geometric constraints'
    write (stdout,'(3x,a,a)') colorify('[dynamics]','yellow'), &
      & repeat(' ',20)//' — MD length, timestep, temperature, …'
    write (stdout,'(3x,a,a)') colorify('  [[dynamics.meta]]','yellow'), &
      & repeat(' ',11)//' — metadynamics bias settings'
    write (stdout,'(3x,a,a)') colorify('[cregen]','yellow'), &
      & repeat(' ',22)//' — ensemble sorting thresholds'
    write (stdout,'(3x,a,a)') colorify('[thermo]','yellow'), &
      & repeat(' ',22)//' — thermochemistry settings'
    write (stdout,*)
    write (stdout,'(1x,a)') 'Full TOML keyword reference:'
    write (stdout,'(3x,a)') colorify('https://crest-lab.github.io/crest-docs/','blue')
    write (stdout,*)

  end select
  call confscript_morehelp2()
  stop '   [-h] displayed. exit.'
end subroutine confscript_morehelp

subroutine confscript_morehelp2
  use iomod,only:colorify
  use crest_parameters,only:stdout
  implicit none
  write (stdout,'(/,1x,a)') 'For detailed help on option groups, use:'
  write (stdout,'(3x,a)') colorify('--help general','gold')//'    '// &
    & colorify('--help compare','gold')//'    '//colorify('--help conf','gold')
  write (stdout,'(3x,a)') colorify('--help thermo','gold')//'     '// &
    & colorify('--help qcg','gold')//'        '//colorify('--help msreact','gold')
  write (stdout,'(3x,a)') colorify('--help other','gold')//'      '// &
    & colorify('--help toml','gold')
  write (stdout,*)
  write (stdout,'(1x,a,a)') 'View literature references with ',colorify('--cite','green')
  write (stdout,'(1x,a)') 'For detailed documentation refer to:'
  write (stdout,'(3x,a)') colorify('https://crest-lab.github.io/crest-docs/','blue')
  write (stdout,*)
end subroutine confscript_morehelp2

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine crestcite
  write (*,*)
  write (*,'(4x,''MAIN REFERENCES:'')')
  write (*,'(/5x,''• P. Pracht, F. Bohle, S. Grimme,'')')
  write (*,'( 5x,''  PCCP, 2020, 22, 7169-7192.'')')
  write (*,'(/5x,''• S.Grimme, JCTC, 2019, 15, 2847-2862.'')')
  write (*,'(/5x,''• P.Pracht, S.Grimme, Chem. Sci., 2021, 12, 6551-6568.'')')
  write (*,'(/5x,''• S.Spicher, C.Plett, P.Pracht, A.Hansen, S.Grimme,'')')
  write (*,'( 5x,''  JCTC, 2022, 18 (5), 3174-3189.'')')
  write (*,'(/5x,''• P.Pracht, C.Bannwarth, JCTC, 2022, 18 (10), 6370-6385.'')')
  write (*,'(/3x,''• P.Pracht, S.Grimme, C.Bannwarth, F.Bohle, S.Ehlert,'')')
  write (*,'( 3x,''  G.Feldmann, J.Gorges, M.Müller, T.Neudecker, C.Plett,'')')
  write (*,'( 3x,''  S.Spicher, P.Steinbach, P.Wesołowski, F.Zeller,'')')
  write (*,'( 3x,''  J. Chem. Phys., 2024, 160, 114110.'')')

  write (*,'(/,/)')
  write (*,'(4x,''GFNn-xTB references:'')')
  write (*,'(5x,''GFN1-xTB'')')
  write (*,'(5x,''• S.Grimme, C.Bannwarth, P.Shushkov, JCTC, 2017,'')')
  write (*,'(5x,''  13, 1989-2009. DOI: 10.1021/acs.jctc.7b00118'')')
  write (*,'(5x,''GFN2-xTB'')')
  write (*,'(5x,''• C.Bannwarth, S.Ehlert and S.Grimme., JCTC, 2019,'')')
  write (*,'(5x,''  15, 1652-1671. DOI: 10.1021/acs.jctc.8b01176'')')
  write (*,'(5x,''GFN0-xTB'')')
  write (*,'(5x,''• P.Pracht, E.Caldeweyher, S.Ehlert, S.Grimme, 2019,'')')
  write (*,'(5x,''  ChemRxiv preprint, DOI: 10.26434/chemrxiv.8326202.v1'')')
  write (*,'(5x,''GFN-FF'')')
  write (*,'(5x,''• S.Spicher, S.Grimme, Angew. Chem. Int. Ed., 2020,'')')
  write (*,'(5x,''  132, 2-11, DOI: 10.1002/ange.202004239'')')

  write (*,'(/,/)')
  write (*,'(4x,''related references:'')')
  write (*,'(5x, ''• S.Grimme, C.Bannwarth, S.Dohm, A.Hansen,'')')
  write (*,'(5x, ''  J.Pisarek, P.Pracht, J.Seibert, F.Neese,'')')
  write (*,'(5x, ''  Angew. Chem. Int. Ed., 2017, 56, 14763-14769'')')
  write (*,'(/5x,''• P.Pracht, C.A.Bauer, S.Grimme, JCC, 2017,'')')
  write (*,'(5x, ''  38, 2618–2631, DOI: 10.1002/jcc.24922'')')
  write (*,'(/5x,''• P.Pracht, R.Wilcken, A.Udvarhelyi, S.Rodde, S.Grimme,'')')
  write (*,'(5x, ''  JCAMD, 2018, 32, 1139-1149.'')')
  write (*,'(/5x,''• P.Pracht, S.Grimme, JPCA, 2021, 125, 5681-5692'')')
  write (*,'(/5x,''• J.Gorges, S.Grimme, A.Hansen, P.Pracht,'')')
  write (*,'(5x, ''  PCCP, 2022,24, 12249-12259.'')')
  write (*,'(/5x,''• C.Zurek, R.A.Malleav, A.C.Paul, N.van Staalduinen, et al.,'')')
  write (*,'(5x, ''  J. Chem. Theory Comput. 2025, 21, 3, 1459-1475.'')')

  write (*,'(/,/)')
  write (*,'(3x,''Please cite work conducted with this code appropriately.'')')
  stop '   [--cite] displayed. exit.'
end subroutine crestcite

!========================================================================================!

subroutine crestcrest
  write (*,'(7x,''|                                            |'')')
  write (*,'(7x,''|              ###############               |'')')
  write (*,'(7x,''|              #//////|@@@@@@#               |'')')
  write (*,'(7x,''|              #\\\\\\|@@@@@@#               |'')')
  write (*,'(7x,''|              #//////|@@@@@@#               |'')')
  write (*,'(7x,''|              #@@@@@@|\\\\\\#               |'')')
  write (*,'(7x,''|              #@@@@@@|//////#               |'')')
  write (*,'(7x,''|               #@@@@@|\\\\\#                |'')')
  write (*,'(7x,''|                ###@@|//###                 |'')')
  write (*,'(7x,''|                   #####                    |'')')
  write (*,'(7x,''|                                            |'')')
  write (*,'(7x,''|                 C R E S T                  |'')')
end subroutine crestcrest

!========================================================================================!

subroutine prchd
  write (*,*)
  write (*,'(7x,''========================================'')')
  write (*,'(7x,''|             C R E G E N              |'')')
  write (*,'(7x,''|     conformer/rotamer generation     |'')')
  write (*,'(7x,''| & NMR symmetry/equivalence analysis  |'')')
  write (*,'(7x,''|      SG, Universitaet Bonn, MCTC     |'')')
  write (*,'(7x,''|     Fri Aug 11 13:16:16 CEST 2017    |'')')
  write (*,'(7x,''========================================'')')
  write (*,*)
end

!========================================================================================!

subroutine header_stereo
  write (*,*)
  write (*,'(5x,''========================================'')')
  write (*,'(5x,''|   automated stereoisomer generator   |'')')
  write (*,'(5x,''|              P.Pracht                |'')')
  write (*,'(5x,''|       Universitaet Bonn, MCTC        |'')')
  write (*,'(5x,''|    Wed 9. Oct 11:20:34 CEST 2019     |'')')
  write (*,'(5x,''========================================'')')
  write (*,*)

  write (*,*) 'NOTE: This is a work-in-progress project!'
end subroutine header_stereo

!========================================================================================!

subroutine prreactorhd
  write (*,*)
  write (*,'(7x,''========================================'')')
  write (*,'(7x,''|          GFNn-xTB NANOREACTOR        |'')')
  write (*,'(7x,''|      SG, Universitaet Bonn, MCTC     |'')')
  write (*,'(7x,''========================================'')')
  write (*,'(/,7x,''JCTC, 2019, 15, 2847-2862.'')')
  write (*,*)
end

!========================================================================================!

subroutine zsortwarning2(env)
  use crest_data
  implicit none
  type(systemdata) :: env
  logical :: ex
  inquire (file=env%constraints,exist=ex)
  if (ex.and.env%autozsort) then
    write (*,*) '==========================================='
    write (*,*) 'WARNING:'
    write (*,*) 'The input coordinate file would be sorted'
    write (*,*) 'by zsort and a constraining file is'
    write (*,*) 'present. To avoid constrainment of the'
    write (*,*) 'wrong atoms zsort will be turned off.'
    write (*,*) 'This also might influence the results.'
    write (*,*) '==========================================='
    write (*,*)
    env%autozsort = .false.
  end if
end subroutine zsortwarning2

!========================================================================================!

subroutine msreact_head()
  implicit none
  write (*,*)
  write (*,'(2x,''========================================'')')
  write (*,'(2x,''|                                      |'')')
  write (*,'(2x,''|               MSREACT                |'')')
  write (*,'(2x,''| automated MS fragment generator      |'')')
  write (*,'(2x,''|                                      |'')')
  write (*,'(2x,''|       University of Bonn, MCTC       |'')')
  write (*,'(2x,''========================================'')')
  write (*,'(2x,'' S. Grimme, P. Pracht, J. Gorges.'')')
  write (*,*)
  write (*,'(3x,''Cite work conducted with this code as'')')
  write (*,'(/,3x,''Philipp Pracht, Stefan Grimme, Christoph Bannwarth, Fabian Bohle, Sebastian Ehlert, Gereon Feldmann,'')')
  write (*,'(3x,''Johannes Gorges, Marcel Müller, Tim Neudecker, Christoph Plett, Sebastian Spicher, Pit Steinbach,'')')
  write (*,'(3x,''Patryk A. Wesolowski, and Felix Zeller J. Chem. Phys., 2024, submitted.'')')
  write (*,*)
end subroutine msreact_head

!========================================================================================!

subroutine smallhead(str)
!**********************************************
!> convert a string in a small header printout
!**********************************************
  use crest_parameters,only:stdout
  implicit none
  character(len=*) :: str
  integer :: strlen
  character(len=:),allocatable :: str2
  strlen = len_trim(str)
  str2 = repeat('-',strlen)
  write (stdout,'(1x,a)') trim(str2)
  write (stdout,'(1x,a)') trim(str)
  write (stdout,'(1x,a)') trim(str2)
  return
end subroutine smallhead
subroutine smallheadline(line)
  implicit none
  character(len=*) :: line
  integer :: lw
  lw = len_trim(line)
  write (*,'(/,1x,a)') repeat('=',lw)
  write (*,'(1x,a)') trim(line)
  write (*,'(1x,a)') repeat('=',lw)
  return
end subroutine smallheadline
subroutine underline(str)
  implicit none
  character(len=*) :: str
  integer :: strlen
  character(len=:),allocatable :: str2
  strlen = len_trim(str)
  str2 = repeat('-',strlen)
  write (*,'(1x,a)') trim(str)
  write (*,'(1x,a)') trim(str2)
  return
end subroutine underline

!========================================================================================!

function str_center_align(str,ilen) result(res)
  character(len=*),intent(in) :: str
  integer,intent(in) :: ilen
  character(len=:),allocatable :: res
  integer :: i,slen,str_len,pad_left,pad_right

  str_len = len_trim(str)
  if (str_len >= ilen) then
    slen = str_len+2
  else
    slen = ilen
  end if
  write (*,*) slen,ilen,str_len
  pad_left = (slen-str_len)/2
  pad_right = slen-pad_left-str_len

  res = ""
  write (*,*) pad_left
  do i = 1,pad_left
    res = res//" "
  end do
!res = repeat(" ",pad_left)
  res = res//trim(str)
end function str_center_align

!========================================================================================!

subroutine mtdwarning(lenv)
  implicit none
  real*8 :: lenv
!> a warning when the MTD length exceeds 200ps
  write (*,*)
  write (*,'(a,f5.1,a)') "! WARNING: the estimated MTD time exceeds ",lenv," ps."
  write (*,'(a       )') "! Because the estimate is uncertain, the program restricts"
  write (*,'(a,f5.1,a)') "! this to ",lenv," ps and continues. The user may"
  write (*,'(a       )') "! re-run crest with manual setting by '-mdlen <time>' and"
  write (*,'(a       )') "! check the results carefully."
  write (*,*)

end subroutine mtdwarning

!========================================================================================!

subroutine printiter
  implicit none
  write (*,*)
  write (*,'(80("*"))')
  write (*,'("**",20x,"N E W   I T E R A T I O N  C Y C L E",20x,"**")')
  write (*,'(80("*"))')
end subroutine printiter
subroutine printiter2(i)
  implicit none
  integer :: i
  write (*,*)
  write (*,'(80("*"))')
  write (*,'("**",21x,"I T E R A T I O N    C Y C L E    ",i3,18x,"**")') i
  write (*,'(80("*"))')
end subroutine printiter2
subroutine printiter3(text,i)
  implicit none
  character(len=*),intent(in) :: text
  integer,intent(in) :: i
  character(len=128) :: atmp
  write (atmp,'(a,1x,i0)') trim(text),i
  !call largehead(trim(atmp))
  call construct_boxed_headline(trim(atmp),80,.true.)
end subroutine printiter3

!========================================================================================!

subroutine largehead(str)
!**********************************************
!* convert a string in a large header printout
!**********************************************
  implicit none
  character(len=*) :: str
  call construct_large_headline('*',str)
  return
end subroutine largehead
subroutine largehead2(str)
  implicit none
  character(len=*) :: str
  call construct_large_headline('+',str)
  return
end subroutine largehead2
subroutine construct_large_headline(symb,str)
  implicit none
  character(len=1) :: symb
  character(len=*) :: str
  integer :: strlen,strlen2
  integer :: k
  integer :: i,j
  character(len=128) :: str2
  character(len=128) :: str3
  strlen = len_trim(str)
  str2 = repeat(symb,80)
  strlen2 = len_trim(str2)
  if (strlen .ge. strlen2) then
    str2 = ''
    do i = 1,strlen+6
      str2 = trim(str2)//symb
    end do
    strlen = strlen+6
  end if
  k = strlen2-strlen
  j = k/2-2
  str3 = symb//symb
  str3(j+1:) = trim(str)
  str3(strlen2-1:strlen2) = symb//symb
  write (*,*)
  write (*,'(a)') trim(str2)
  write (*,'(a)') trim(str3)
  write (*,'(a)') trim(str2)
  return
end subroutine construct_large_headline
subroutine construct_boxed_headline(str,width,bold)
  use crest_parameters,only:stdout
  implicit none
  character(len=*),intent(in) :: str
  integer,intent(in) :: width
  logical,intent(in) :: bold
  integer :: strlen,strlen2
  integer :: k,i,j,jj
  integer :: wid
  wid = max(width,len_trim(str)+4)-2
  wid = width-2
  strlen = len_trim(str)
  if (strlen > wid) wid = strlen
  write (stdout,*)
  if (bold) then
    write (stdout,'(a)') "┏"//repeat("━",wid)//"┓"
  else
    write (stdout,'(a)') "┌"//repeat("─",wid)//"┐"
  end if
  strlen2 = wid+2
  k = strlen2-strlen
  j = k/2
  jj = k-j-2
  if (bold) then
    write (stdout,'(a)') "┃"//repeat(" ",j)//trim(str)//repeat(" ",jj)//"┃"
    write (stdout,'(a)') "┗"//repeat("━",wid)//"┛"
  else
    write (stdout,'(a)') "│"//repeat(" ",j)//trim(str)//repeat(" ",jj)//"│"
    write (stdout,'(a)') "└"//repeat("─",wid)//"┘"
  end if
  return
end subroutine construct_boxed_headline

!========================================================================================!

subroutine print_crest_metadata()
!********************************
!* print metadata from include
!********************************
  include 'crest_metadata.fh'
  integer :: l
  write (*,'(2x,a,t22,":   ",a)') 'CREST version    ',version
  write (*,'(2x,a,t22,":   ",a)') 'timestamp        ',date
  write (*,'(2x,a,t22,":   ",a)') 'commit           ',commit
  l = len_trim(author)
  if (author(1:2) .eq. "'@") then
    write (*,'(2x,a,t22,":   ",a)') 'compiled by      ',"usr"//author(2:l-1)
  else
    write (*,'(2x,a,t22,":   ",a)') 'compiled by      ',author(2:l-1)
  end if
  write (*,'(2x,a,t22,":   ",a)') 'Fortran compiler ',fcompiler
  write (*,'(2x,a,t22,":   ",a)') 'C compiler       ',ccompiler
  write (*,'(2x,a,t22,":   ",a)') 'build system     ',bsystem
  write (*,'(2x,a,t22,":   ",a)') '-DWITH_TOMLF     ',tomlfvar
  write (*,'(2x,a,t22,":   ",a)') '-DWITH_GFN0      ',gfn0var
  write (*,'(2x,a,t22,":   ",a)') '-DWITH_GFNFF     ',gfnffvar
  write (*,'(2x,a,t22,":   ",a)') '-DWITH_TBLITE    ',tblitevar
  write (*,'(2x,a,t22,":   ",a)') '-DWITH_DDX       ',ddxvar
  write (*,'(2x,a,t22,":   ",a)') '-DWITH_LIBPVOL   ',libpvolvar
  write (*,'(2x,a,t22,":   ",a)') '-DWITH_LWONIOM   ',lwoniomvar
  write (*,'(2x,a,t22,":   ",a)') '-DWITH_FMLIP_RELAY',fmliprelayvar
end subroutine print_crest_metadata

subroutine cat_mod(ch,pre,fname,post)
  implicit none
  integer :: ch
  character(len=*) :: pre
  character(len=*) :: fname
  character(len=*) :: post
  character(len=256) :: adum
  integer :: ich,io
  open (newunit=ich,file=fname)
  do
    read (ich,'(a)',iostat=io) adum
    if (io < 0) exit
    write (ch,'(a,a,a)') pre,trim(adum),post
  end do
  close (ich)
  return
end subroutine cat_mod

subroutine checkbinary(env)
  use crest_data
  use iomod,only:checkprog
  implicit none
  type(systemdata) :: env
  integer :: r
  r = 0
  call checkprog(trim(env%ProgName),r)
  if (r .ne. 0) then
    write (*,'(4x,a)') 'Warning! The xtb binary was not found and hence CREST might crash'
  end if
  if (env%crestver .eq. crest_solv) then
    call checkprog(trim('xtbiff'),r)
    if (r .ne. 0) then
      write (*,'(4x,a)') 'Warning! The xtbiff binary was not found and hence the qcg mode in CREST will probably crash'
    end if
  end if
  return
end subroutine checkbinary

!========================================================================================!
!========================================================================================!

subroutine wrGUIpercent(current,maxv,interval)
!********************************************
!* printout percent calculation for GUI mode
!********************************************
  use iso_fortran_env,wp => real64
  implicit none
  integer :: current,maxv,interval
  real(wp) :: perc,inc,trc
  integer :: i,interval2

  if (current == maxv) then
    write (*,'(1x,f6.2,a)') 100.0d0,' percent done. finished loop.'
    return
  end if

  interval2 = interval
  inc = float(maxv)/float(interval)
  if (inc .gt. 1.0d0) interval2 = maxv
  trc = 0.0d0
  do i = 1,interval2
    trc = floor(float(i)*inc)
    if (current == nint(trc)) then
      perc = (float(current)/float(maxv))*100.0d0
      write (*,'(1x,f6.2,a)') perc,' percent done'
      exit
    end if
  end do

  return
end subroutine wrGUIpercent

!=======================================================================================!

subroutine print_frozen(env)
  use crest_parameters
  use crest_data
  implicit none
  type(systemdata) :: env
  integer :: i
  if (env%calc%nfreeze > 0.and.allocated(env%calc%freezelist)) then
    write (stdout,'(/,a)') repeat('-',50)
    write (stdout,'(a)') ' FROZEN ATOMS:'
    do i = 1,env%ref%nat
      if (env%calc%freezelist(i)) then
        write (stdout,'(1x,i0)',advance='no') i
      end if
    end do
    write (stdout,'(/,a)') repeat('-',50)
  end if
end subroutine print_frozen

!========================================================================================!
!========================================================================================!

subroutine progbar(percent,bar)
  use crest_parameters
  implicit none
  real(wp),intent(in) :: percent
  character(len=52),intent(inout) :: bar
  integer :: i
  integer :: done,notdone

  bar = '['

  done = nint(percent/2)
  notdone = 50-done

  do i = 1,done
    bar = trim(bar)//'#'
  end do

  do i = 1,notdone
    bar = trim(bar)//'-'
  end do

  bar = trim(bar)//']'

end subroutine progbar

subroutine printprogbar(percent)
  use crest_parameters
  implicit none
  real(wp),intent(in) :: percent
  character(len=52) :: bar

  if (percent > 0.0_wp) then
    call progbar(percent,bar)
  else
    call progbar(0.0_wp,bar)
  end if
  write (0,FMT="(A1,A52,2x,F6.2,A)",ADVANCE="NO") achar(13), &
  & bar,percent,'% finished.'

  flush (0)
end subroutine printprogbar
!========================================================================================!
!========================================================================================!

subroutine gxtb_dev_warning
  use crest_parameters
  use crest_data,only:status_safety
  write (stdout,*)
  write (stdout,'(a)') "Note: '--gxtb_dev' is deprecated. Use '--gxtb'."
  write (stdout,'(a)') "g-xTB via the tblite API is available with this build."
  write (stdout,*)
  call creststop(status_safety)
end subroutine gxtb_dev_warning

!========================================================================================!

subroutine gxtb_syscall_warning
!*************************************************************
!* Warn the user that this build runs g-xTB by invoking the
!* xtb binary via system calls (slow, temporary workaround).
!* Only prints in the non-WITH_GXTB call path (have_gxtb=F),
!* i.e. when no native tblite g-xTB build is available.
!* Points the user to the static xtb binary of the g-xTB
!* development release and to the (much faster) static CREST
!* binary with a tblite g-xTB build available upon request.
!*************************************************************
  use crest_parameters,only:stdout
  use tblite_api,only:have_gxtb
  implicit none
  integer,parameter :: w = 70  !> inner box width (ASCII content)
  character(len=*),parameter :: hz = '─'
  character(len=*),parameter :: vt = '│'
  if (have_gxtb) return
  write (stdout,*)
  write (stdout,'(1x,a)') '╭'//repeat(hz,w)//'╮'
  call boxline('Note: g-xTB is currently invoked via the xtb binary (system')
  call boxline('calls). This is a temporary implementation and may be slow.')
  call boxline('')
  call boxline('* A static xtb binary providing g-xTB is available from the')
  call boxline('  g-xTB development release at:')
  call boxline('    https://github.com/grimme-lab/g-xtb')
  call boxline('')
  call boxline('* A static CREST binary with a native tblite development build')
  call boxline('  of g-xTB (much faster than the bespoke xtb system call) is')
  call boxline('  available upon request.')
  write (stdout,'(1x,a)') '╰'//repeat(hz,w)//'╯'
  write (stdout,*)
contains
  subroutine boxline(text)
    character(len=*),intent(in) :: text
    character(len=w) :: buf
    buf = '  '//text
    write (stdout,'(1x,a)') vt//buf//vt
  end subroutine boxline
end subroutine gxtb_syscall_warning

!========================================================================================!
!========================================================================================!

subroutine gxtb_solvation_guard(env)
!**************************************************************************
!* g-xTB does not (yet) support dedicated implicit solvation. This guard
!* is meant to be called right after input parsing (CLI and/or TOML),
!* once the dedicated env%calc object is populated.
!*
!* It scans all requested calculation levels (both the native tblite
!* g-xTB level and the xtb-subprocess driver) as well as the legacy CLI
!* selection. If a g-xTB level is combined with a request for implicit
!* solvation, CREST is terminated safely with an explanatory message.
!*
!* State of the g-xTB implementation: 05/2026 (development version) — no
!* solvation support, neither via the native tblite interface nor via the
!* xtb subprocess driver.
!*
!*   env  -  CREST system data (read-only here; holds env%calc)
!**************************************************************************
  use crest_parameters
  use crest_data
  use crest_calculator,only:jobtype
  use tblite_api,only:xtblvl,have_gxtb
  use iomod,only:drawbox
  implicit none
  type(systemdata),intent(in) :: env
  integer :: i
  logical :: usesgxtb,wantssolv
  integer,parameter :: bw = 72   !> box width

  usesgxtb = .false.
  wantssolv = .false.

  ! ── legacy CLI selection (env-level) ────────────────────────────────
  if (index(env%gfnver,'gxtb') .ne. 0) usesgxtb = .true.
  if (env%gbsa) wantssolv = .true.

  ! ── per-calculation levels collected in env%calc ────────────────────
  if (associated(env%calc)) then
    do i = 1,env%calc%ncalculations
      associate (job => env%calc%calcs(i))
        !> native tblite g-xTB level
        if (job%id == jobtype%tblite.and.job%tblitelvl == xtblvl%gxtb) then
          usesgxtb = .true.
        end if
        !> g-xTB run as an xtb subprocess (--gxtb flag or gxtb binary)
        if (allocated(job%other)) then
          if (index(job%other,'gxtb') .ne. 0) usesgxtb = .true.
        end if
        if (allocated(job%binary)) then
          if (index(job%binary,'gxtb') .ne. 0) usesgxtb = .true.
        end if
        !> any implicit solvation request on this level
        if (allocated(job%solvmodel).or.allocated(job%solvent)) then
          wantssolv = .true.
        end if
      end associate
    end do
  end if

  !> the unsupported combination is g-xTB AND implicit solvation
  if (.not. (usesgxtb.and.wantssolv)) return

  ! ── unsupported combination → print box and terminate safely ────────
  write (stdout,*)
  call drawbox(stdout,'',width=bw,charset=6,procedual=0)
  call drawbox(stdout,'g-xTB: IMPLICIT SOLVATION NOT SUPPORTED',width=bw, &
  &            charset=6,procedual=1)
  call drawbox(stdout,'',width=bw,charset=6,procedual=3)
  call gsg_line('')
  call gsg_line('You requested g-xTB together with an implicit solvation model.')
  call gsg_line('')
  call gsg_line('The current g-xTB implementation (development version, 05/2026)')
  call gsg_line('does not support dedicated implicit solvation yet, neither via')
  call gsg_line('the native tblite interface nor via the xtb subprocess driver.')
  call gsg_line('')
  call gsg_line('Please run g-xTB in the gas phase (drop --alpb / --gbsa or the')
  call gsg_line('corresponding TOML solvation keywords), or pick a different')
  call gsg_line('level of theory for solvated calculations.')
  call drawbox(stdout,'',width=bw,charset=6,procedual=2)
  write (stdout,*)

  call creststop(status_safety)

contains
  subroutine gsg_line(txt)
    !> emit a single left-aligned content line of the warning box
    use iomod,only:drawbox
    character(len=*),intent(in) :: txt
    call drawbox(stdout,txt,width=bw,charset=6,procedual=1,padl=2)
  end subroutine gsg_line
end subroutine gxtb_solvation_guard

!========================================================================================!
!========================================================================================!

subroutine crest_no_runtype_selected()
  !*****************************************************
  !* Print an error when no runtype has been selected, *
  !* list the available main runtypes, and stop.       *
  !*****************************************************
  use crest_parameters,only:stdout
  use crest_data,only:status_safety
  use iomod,only:colorify
  implicit none
  write (stdout,*)
  write (stdout,'(1x,a)') repeat('=',60)
  write (stdout,'(1x,a)') 'No runtype was selected.'
  write (stdout,'(1x,a)') 'Please choose one of the main runtypes listed below.'
  write (stdout,'(1x,a)') repeat('=',60)
  write (stdout,*)
  write (stdout,'(3x,a)') 'Main runtypes:'
  write (stdout,*)
  write (stdout,'(5x,a,t30,a)') '--sp','Single-point energy calculation'
  write (stdout,'(5x,a,t30,a)') '--opt','Structure optimization'
  write (stdout,'(5x,a,t30,a)') '--md','Molecular dynamics simulation'
  write (stdout,'(5x,a,t30,a)') '--imtdgc/--v3','iMTD-GC conformational search '//colorify('(PREVIOUS DEFAULT)','gold')
  write (stdout,'(5x,a,t30,a)') '--entropy','Entropy/free-energy sampling'
  write (stdout,'(5x,a,t30,a)') '--mdopt','Ensemble optimization (no sorting)'
  write (stdout,'(5x,a,t30,a)') '--screen','Ensemble screening'
  write (stdout,'(5x,a,t30,a)') '--protonate','Protonation site search'
  write (stdout,'(5x,a,t30,a)') '--deprotonate','Deprotonation site search'
  write (stdout,'(5x,a,t30,a)') '--tautomerize','Tautomer generation'
  write (stdout,'(5x,a,t30,a)') '--qcg','QCG workflows'
  write (stdout,'(5x,a,t30,a)') '--msreact','MSREACT workflows'
  !write (stdout,'(5x,a,t30,a)') '--bh','Basin-hopping global optimization'
  write (stdout,'(5x,a,t30,a)') '--sort','Ensemble sorting (CREGEN)'
  write (stdout,'(5x,a,t30,a)') '--ttconf','TTConf conformational sampling'

  write (stdout,*)
  write (stdout,'(3x,a)') 'For TOML input files use:  crest --input <file.toml>'
  write (stdout,'(3x,a)') 'For the full option list:  crest --help'
  write (stdout,*)
  call creststop(status_safety)
end subroutine crest_no_runtype_selected

!========================================================================================!

subroutine crest_output_summary(env)
!**************************************************
!* Print a table of output files written by CREST *
!* and a short description of their contents.     *
!* Which files are listed depends on the runtype. *
!**************************************************
  use crest_parameters,only:stdout
  use crest_data
  implicit none
  type(systemdata),intent(in) :: env
  character(len=72),parameter :: hbar = repeat('-',80)
  logical :: lexists

  select case (env%crestver)
  case (crest_imtd,crest_imtd2,crest_ttc,crest_screen,crest_mdopt, &
      & crest_sorting,crest_optimize,crest_trialopt,crest_rigcon, &
      & crest_moldyn,crest_bh,crest_bhpt,crest_protonate,crest_deprotonate, &
      & crest_tautomerize,crest_numhessian,crest_ensemblehess)
    write (stdout,*)
    write (stdout,*)
    write (stdout,'(2x,a)') 'Important files written by CREST during this run'
    write (stdout,'(1x,a)') hbar
  end select

  select case (env%crestver)
    ! ── iMTD-GC / sMTD-iMTD conformer search ─────
  case (crest_imtd,crest_imtd2)
    call wfe('crest_conformers.xyz','unique conformers (1 per rotamer group), energy-sorted')
    call wfe('crest_rotamers.xyz','all structures including rotamers, with Boltzmann weights')
    call wfe('crest_best.xyz','lowest-energy conformer')
    call wfe('cregen.full','full CREGEN output (written when terminal output is abbreviated)')
    call wfe('crest.restart','restart/checkpoint file for the iMTD-GC algorithm')

    ! ── TTConf-light (tensor-train) conformer search ──
  case (crest_ttc)
    call wfe('crest_conformers.xyz','unique conformers (1 per rotamer group), energy-sorted')
    call wfe('crest_rotamers.xyz','all structures including rotamers, with Boltzmann weights')
    call wfe('crest_best.xyz','lowest-energy conformer')
    call wfe('cregen.full','full CREGEN output (written when terminal output is abbreviated)')

    ! ── ensemble screening / MDOPT ────────────────
  case (crest_screen)
    call wfe('crest_ensemble.xyz','ensemble output')
    call wfe('crest_conformers.xyz','unique conformers (1 per rotamer group), energy-sorted')
    call wfe('crest_rotamers.xyz','all structures including rotamers')
    call wfe('crest_best.xyz','lowest-energy structure')

  case (crest_mdopt)
    call wfe('crest_ensemble.xyz','ensemble output with optimized geometries')
    if (allocated(env%refine_queue)) &
    & call wfe('crest_ensemble.xyz.refine','refined ensemble (after the requested --refine step)')

    ! ── standalone CREGEN sorting ─────────────────
  case (crest_sorting)
    call wfe('crest_conformers.xyz','unique conformers (1 per rotamer group), energy-sorted')
    call wfe('crest_rotamers.xyz','all structures including rotamers')
    call wfe('crest_best.xyz','lowest-energy structure')
    call wfe('cregen.full','full CREGEN output (written when terminal output is abbreviated)')

    ! ── geometry optimization ─────────────────────
  case (crest_optimize,crest_trialopt,crest_rigcon)
    call wfe('crestopt.xyz','final optimized geometry')
    call wfe('crestopt.log.xyz','step-by-step optimization trajectory')
    if (env%crest_ohess) then
      call wfe('numhess','Hessian matrix in Turbomole format (2nd derivatives, a.u.)')
      call wfe('vibspectrum','vibrational frequencies (cm⁻¹) and IR intensities (km/mol)')
      call wfe('g98.out','frequencies and normal modes in Gaussian output format')
    end if

    ! ── numerical Hessian ─────────────────────────
  case (crest_numhessian)
    call wfe('numhess','Hessian matrix in Turbomole format (2nd derivatives, a.u.)')
    call wfe('vibspectrum','vibrational frequencies (cm⁻¹) and IR intensities (km/mol)')
    call wfe('g98.out','frequencies and normal modes in Gaussian 98 output format')
    call wfe('dipgrad','Cartesian dipole gradient ∂μ/∂x (a.u.), if available')

    ! ── ensemble Hessians ─────────────────────────
  case (crest_ensemblehess)
    call wfe('crest_ensemble.xyz','input ensemble with Gibbs free energies as comments')
    call wfe('crest.energies','plain list of total Gibbs free energies (Eh)')

    ! ── molecular dynamics ────────────────────────
  case (crest_moldyn)
    call wfe('crest_dynamics.trj.xyz','MD trajectory (snapshots)')

    ! ── basin-hopping ─────────────────────────────
  case (crest_bh,crest_bhpt)
    call wfe('crest_conformers.xyz','unique conformers found by basin-hopping, energy-sorted')
    call wfe('crest_rotamers.xyz','all accepted structures including rotamers')
    call wfe('crest_best.xyz','lowest-energy structure')

    ! ── protonation / deprotonation / tautomers ───
  case (crest_protonate)
    call wfe('protonated.xyz','unique site candidates, energy-sorted')
    call wfe('crest_best.xyz','lowest-energy structure')
  case (crest_deprotonate)
    call wfe('deprotonated.xyz','unique site candidates, energy-sorted')
    call wfe('crest_best.xyz','lowest-energy structure')
  case (crest_tautomerize)
    call wfe('tautomers.xyz','unique site candidates, energy-sorted')
    call wfe('crest_best.xyz','lowest-energy structure')

  case default
    ! no dedicated file summary for this runtype

  end select

  ! ── if PCA clustering was performed, report the cluster file ─────────────────
  select case (env%crestver)
  case (crest_imtd,crest_imtd2,crest_sorting)
    inquire (file=clusterfile,exist=lexists)
    if (lexists) call wfe(clusterfile,'representative structures from PCA/k-means clustering')
  end select

  select case (env%crestver)
  case (crest_imtd,crest_imtd2,crest_ttc,crest_screen,crest_mdopt, &
      & crest_sorting,crest_optimize,crest_trialopt,crest_rigcon, &
      & crest_moldyn,crest_bh,crest_bhpt,crest_protonate,crest_deprotonate, &
      & crest_tautomerize,crest_numhessian,crest_ensemblehess)

    write (stdout,'(1x,a)') hbar
  end select

contains
  subroutine wfe(fname,descr)
    character(len=*),intent(in) :: fname,descr
    write (stdout,'(1x,a,t30,''│ '',a)') fname,trim(descr)
  end subroutine wfe
end subroutine crest_output_summary
