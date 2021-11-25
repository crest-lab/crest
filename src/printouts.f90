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

!===================================================================================================!
! Print Header
!===================================================================================================!
subroutine confscript_head
      implicit none
      character(len=40),parameter:: date='Mon 16. Aug 09:59:32 CEST 2021'
      character(len=10),parameter:: version='2.11.1'
      logical :: niceprint
      
      niceprint=.true.
      write(*,*)
      write(*,'(7x,''=============================================='')')
      write(*,'(7x,''|                                            |'')')
      write(*,'(7x,''|                 C R E S T                  |'')')
      write(*,'(7x,''|                                            |'')')
      write(*,'(7x,''|  Conformer-Rotamer Ensemble Sampling Tool  |'')')
      write(*,'(7x,''|          based on the GFN methods          |'')')
      write(*,'(7x,''|             P.Pracht, S.Grimme             |'')')
      write(*,'(7x,''|          Universitaet Bonn, MCTC           |'')')
      write(*,'(7x,''=============================================='')')
      write(*,'(7x,''Version '',a,'', '',a)')trim(version),trim(date)
      write(*,'(2x,''Using the xTB program. Compatible with xTB version 6.4.0'')')
      write(*,*)
      write(*,'(3x,''Cite work conducted with this code as'')')
      write(*,'(/,3x,''P. Pracht, F. Bohle, S. Grimme, PCCP, 2020, 22, 7169-7192.'')')
      write(*,'(/,3x,''and  S. Grimme, JCTC, 2019, 15, 2847-2862.'')')
      write(*,*)

      write(*,'(3x,a)')'with help from:'
      write(*,'(3x,a)')'C.Bannwarth, F.Bohle, S.Ehlert, S.Grimme,'
      write(*,'(3x,a)')'C. Plett, P.Pracht, S. Spicher'
      write(*,*)

      call disclaimer()

end subroutine confscript_head

subroutine disclaimer

      write(*,'(3x,a)')'This program is distributed in the hope that it will be useful,'
      write(*,'(3x,a)')'but WITHOUT ANY WARRANTY; without even the implied warranty of'
      write(*,'(3x,a)')'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.'

end subroutine disclaimer

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Confscript help printout
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine confscript_help()
      use crest_data
      implicit none

      write(*,*)
      write(*,'(1x, ''usage        :'')')
      write(*,'(1x, ''crest [input] [options]'')')
      write(*,*)
      write(*,'(1x, ''The FIRST argument CAN be a coordinate file in the'')')
      write(*,'(1x, ''TM (coord, Bohr) or Xmol (*.xyz, Ang.) format.'')')
      write(*,'(1x, ''If no such file is present as the first argument crest will'')')
      write(*,'(1x, ''automatically search for a file called "coord" in the TM format.'')')
      write(*,'(/,1x,''General and technical options:'')')
      write(*,'(5x,''-v1                : use the MF-MD-GC workflow. (OUTDATED)'')')
      write(*,'(5x,''-v2                : use the MTD-GC workflow. (OUTDATED)'')')
      write(*,'(5x,''-v3 (or -v2i)      : use the iMTD-GC workflow. [default]'')')
      write(*,'(5x,''-v4                : use the iMTD-sMTD workflow.'')')
      write(*,'(5x,''-entropy           : the same workflow as with "-v4", specialized'')')
      write(*,'(5x,''                     for the calculation of conformational entropy.'')')
      write(*,'(5x,''-xnam <"bin">      : specify name of the xtb binary'')')
      write(*,'(5x,''                     that should be used.'')')
      write(*,'(5x,''-niceprint         : progress bar printout for optimizations'')')
      write(*,'(5x,''-dry               : perform a "dry run". Only prints the settings'')')
      write(*,'(5x,''                     that would be applied with the CMD input'')')
      write(*,'(5x,''                     and stops the run before any calculations'')')
      write(*,'(5x,''-T <int>           : Set total number of CPUs(threads)'')')
      write(*,'(5x,''                     to be used. Parallel settings are then'')')
      write(*,'(5x,''                     determined automatically for each step.'')')
      write(*,'(5x,''                     If not set by "-T", this number is read'')')
      write(*,'(5x,''                     from the OMP_NUM_THREADS global variable.'')')
  write(*,*)
  write(*,*)

      write(*,'(1x,''Calculation options:'')')
      write(*,'(5x,''-g <string>        : use GBSA implicit solvent'')')
      write(*,'(5x,''                     for solvent <string>'')')
      write(*,'(5x,''-alpb <string>     : use ALPB implicit solvent'')')
      write(*,'(5x,''                     for solvent <string>'')')
      write(*,'(5x,''-chrg <int>        : set the molecules´ charge'')')
      write(*,'(5x,''-uhf <int>         : set <int>=Nα-Nβ electrons'')')
      write(*,'(5x,''-charges <file>    : copy a existing atomic charges file for'')')
      write(*,'(5x,''                     all optimizations (can be used by GFN-FF)'')')
      write(*,'(5x,''-nozs              : do not perform z-mat sorting [default]'')')
      write(*,'(5x,''-opt <"lev">       : set opt. level for ALL GFN-xTB'')')
      write(*,'(5x,''                     optimizations.'')')
      write(*,'(5x,''                     <lev>=vloose,loose,normal,tight,vtight'')')
      write(*,'(5x,''                     [default: vtight]'')')
      write(*,'(5x,''-gfn1              : use GFN1-xTB'')')
      write(*,'(5x,''-gfn2              : use GFN2-xTB [default]'')')
      !write(*,'(5x,''-gfn0              : use GFN0-xTB'')')
      write(*,'(5x,''-gff, -gfnff       : use GFN-FF (requires xtb 6.3 or newer)'')')
      write(*,'(5x,''                     (for GFN-FF searches bond constraints are applied automatically)'')')
      write(*,'(5x,''-gfn2//gfnff       : GFN2-xTB//GFN-FF composite mode)'')')
      write(*,'(3x,''Adding additional constraints to the calculations:'')')
      write(*,'(3x,''The user is able to include additional constraints to ALL'')')
      write(*,'(3x,''xtb calculations that are conducted by CREST.'')')
      write(*,'(5x,''-cinp <file>       : read in a file containing the constraints.'')')
      write(*,'(5x,''                     constraints have to be in the same format as in xtb.'')')
      write(*,'(5x,''                     (this was done previously via the ".constrains" file)'')')
      write(*,'(5x,''-cbonds            : define automatic bond constraints (set up from topology)'')')
      write(*,'(5x,''-nocbonds          : turn of -cbonds (for GFN-FF, mainly. see above)'')')
      write(*,'(5x,''-fc <float>        : define force constant for defined constraints (-cbonds)'')')
  write(*,*)
  write(*,*)

      write(*,'(1x,''Options for ensemble comparisons:'')')
      write(*,'(5x,''-cregen [file]     : use ONLY the CREGEN subroutine'')')
      write(*,'(5x,''                     to sort a given ensemble file.'')')
      write(*,'(5x,''-ewin <real>       : set energy window in kcal/mol,'')')
      write(*,'(5x,''                     [default: 6.0 kcal/mol]'')')
      write(*,'(5x,''-rthr <real>       : set RMSD threshold in Ang,'')')
      write(*,'(5x,''                     [default: 0.125 Ang]'')')
      write(*,'(5x,''-ethr <real>       : set E threshold in kcal/mol,'')')
      write(*,'(5x,''                     [default: 0.05 kcal/mol]'')')
      write(*,'(5x,''-bthr <real>       : set Rot. const. threshold ,'')')
      write(*,'(5x,''                     [default: 0.01 (= 1%)]'')')
      write(*,'(5x,''-pthr <real>       : Boltzmann population threshold'')')
      write(*,'(5x,''                     [default: 0.05 (= 5%)]'')')
      write(*,'(5x,''-temp <real>       : set temperature in cregen, [default: 298.15 K]'')')
      write(*,'(5x,''-prsc              : create a scoord.* file for each conformer'')')
      write(*,'(5x,''-nowr              : don´t write new ensemble files'')')
      write(*,'(5x,''-eqv,-nmr,-entropy : compare nuclear equivalences (requires rotamers)'')')
      write(*,'(5x,''-cluster <int>     : PCA and k-Means clustering of sorted ensemble.'')')
      write(*,'(5x,''                     Works as extenstion to the CREGEN sorting.'')')
      write(*,'(5x,''                     <int> is the number of clusters to be formed.'')')
      write(*,'(5x,''-notopo            : turn off any topology checks in CREGEN.'')')
  write(*,*)
  write(*,*)
      write(*,'(1x,''Options for the iMTD-GC workflows:'')')
      write(*,'(5x,''-cross             : do the GC part [default]'')')
      write(*,'(5x,''-nocross           : don´t do the GC part'')')
      write(*,'(5x,''-shake <int>       : set SHAKE mode for MD'')')
      write(*,'(5x,''                     (0=off,1=H-only,2=all bonds) [default: 2]'')')
      write(*,'(5x,''-tstep <int>       : set MD time step in fs'')')
      write(*,'(5x,''                     [default: 5 fs]'')')
      write(*,'(5x,''-mdlen/-len <real> : set MD length (all MTDs) in ps.'')')
      write(*,'(5x,''                     Also possible are multiplicative factors'')')
      write(*,'(5x,''                     for the default MD length with "x<real>"'')')
      write(*,'(5x,''-mddump <int>      : xyz dumpstep to Trajectory in fs'')')
      write(*,'(5x,''                     [default: 100 fs]'')')
      write(*,'(5x,''-vbdump <real>     : set Vbias dump frequency in ps'')')
      write(*,'(5x,''                     [default: 1.0 ps]'')')
      write(*,'(5x,''-tnmd <real>       : set temperature for additional normal MDs'')')
      write(*,'(5x,''                     [default: 400 K]'')')
      write(*,'(5x,''-norotmd           : don´t do the regular MDs after the '')')
      write(*,'(5x,''                     second multilevel optimization step'')')
      write(*,'(5x,''-hflip/-noflip     : turn on/off a small enhancement routine to'')')
      write(*,'(5x,''                     rotate OH groups after MTD. [default: OFF]'')')
      write(*,'(5x,''-maxflip <int>     : max. number of new structures by the above'')')
      write(*,'(5x,''                     enhancement routine. [default: 1000]'')')
      write(*,'(5x,''-quick             : perform a search with reduced settings'')')
      write(*,'(5x,''                     for a crude ensemble.'')')
      write(*,'(5x,''-squick            : perform a even further reduced search'')')
      write(*,'(5x,''-mquick            : perform a search with maximum reduced settings'')')
      write(*,'(5x,''                     (do not reduce the settings more than that)'')')
      write(*,'(5x,''-origin            : track the step of generation for'')')
      write(*,'(5x,''                     each conformer/rotamer. [default]'')')
      write(*,'(5x,''-keepdir           : keep sub-directories of the conformer'')')
      write(*,'(5x,''                     generation step.'')')
      write(*,'(5x,''-nci               : generate an ellipsoide potential around the'')')
      write(*,'(5x,''                     input structure and add it to the MTD simulation.'')')
      write(*,'(5x,''                     This can be used to find aggregates of NCI complexes.'')')
      write(*,'(5x,''-pscal <real>      : scale the ellipsoide potential axes by factor <real>.'')')
  write(*,*)
  write(*,*)
      write(*,'(1x,''Thermostatistical options (used in entropy mode):'')')
      write(*,'(5x,''-trange <lower> <upper> <step>   : entropies are calculated for different temperatures.'')')
      write(*,'(5x,''                                   these are calculated in a temperature range from'')')
      write(*,'(5x,''                                   <lower> to <upper> with <step> in between.'')')
      write(*,'(5x,''                                   [default: 280K-380K in 10K steps]'')')
      write(*,'(5x,''-fscal <float>     : frequency scaling factor. [default: 1.0]'')')
      write(*,'(5x,''-sthr <float>      : vibrational/rotational entropy interpolation threshold (τ)'')')
      write(*,'(5x,''                     [default: 25.0 cm^-1]'')')
      write(*,'(5x,''-ithr <float>      : imaginary mode inversion cutoff [default: -50.0 cm^-1]'')')
      write(*,'(5x,''-ptot <float>      : sum of population for structures considered in msRRHO average.'')')
      write(*,'(5x,''                     [default: 0.9 (=90%)]'')')
    

      write(*,*)
      write(*,*)
      
      write(*,'(1x,''Quantum Cluster Growth (QCG)'')')
      write(*,'(1x,''General usage :'')')
      write(*,'(5x,''<solute> -qcg <solvent> [options]'')')
      write(*,'(1x,''options:'')')
      write(*,'(5x,''-nopreopt          : do not perform preoptimization (only for qcg).'')')
      write(*,'(5x,''-grow              : cluster generation'')')
      write(*,'(5x,''-nsolv <INT>       : number of solvent molecules to add'')')
      write(*,'(5x,''-samerand          : use same random number for every xtbiff run'')')
      write(*,'(5x,''-ensemble          : ensemble generation'')')
      write(*,'(5x,''-ncimtd            : NCI-MTD CREST ensemble generation (Default)'')')
      write(*,'(5x,''-mtd               : MTD for QCG ensemble generation'')')
      write(*,'(5x,''-md                : normal MD for QCG ensemble search'')')
      write(*,'(5x,''-enslvl [method]   : define a method for ensemble search. All gfn methods are supported'')')
      write(*,'(5x,''-esolv             : reference cluster generation and comp. of solvation energy'')')
      write(*,'(5x,''-gsolv             : reference cluster generation and comp. of solvation free energy'')')
      write(*,'(5x,''-nclus             : defines how many clusters are taken for reference cluster generation'')')
      write(*,'(5x,''                   : default 4'')')
      write(*,'(5x,''-nocff             : switches off the CFF algorithm'')')
      write(*,'(5x,''-freqscal          : defines frequency scale factor. Only for outprint'')')
      write(*,'(5x,''-freqlvl [method]  : define a method for frequency computation. All gfn versions are supported'')')


  write(*,*)
  write(*,*)
      write(*,'(1x,''Other tools for standalone use:'')')
      write(*,'(5x,''-zsort             : use only the zsort subroutine'')')
      write(*,'(5x,''                     to sort the z-matrix of the input'')')
      write(*,'(5x,''                     coord file.'')')
      write(*,'(5x,''-mdopt <file>      : optimize along trajectory or'')')
      write(*,'(5x,''                     ensemble file in the XYZ format.'')')
      write(*,'(5x,''                     Each point on the file is optimized.'')')
      write(*,'(5x,''-screen <file>     : optimize along ensemble file'')')
      write(*,'(5x,''                     in the XYZ format. A multilevel'')')
      write(*,'(5x,''                     optimization is performed with continiously'')')
      write(*,'(5x,''                     increasing thresholds. After each step'')')
      write(*,'(5x,''                     the ensemble file is sorted.'')')
      write(*,'(5x,''-protonate         : find a molecules protomes by using a'')')
      write(*,'(5x,''                     LMO π- or LP-center approach.'')')
      write(*,'(5x,''-deprotonate       : find a molecules deprotomers.'')')
      write(*,'(5x,''-tautomerize       : combine the protonation and deprotonation'')')
      write(*,'(5x,''                     to find prototropic tautomers.'')')
      write(*,'(6x,''↳ -trev           : do first the deprotonation and then the'')')
      write(*,'(8x,''                  protonation in the -tautomerize mode, i.e.,'')')
      write(*,'(8x,''                  reverse of the default procedure.'')')
      write(*,'(6x,''↳ -iter <int>     : set number of protonation/deprotonation cycles'')')
      write(*,'(8x,''                  in the tautomerization script. [default: 2]'')')
      write(*,'(5x,''-compare <f1> <f2> : compare two ensembles <f1> and <f2>.'')')
      write(*,'(5x,''                     Both ensembles must have the same'')')
      write(*,'(5x,''                     order of atoms of the molecule and'')')
      write(*,'(5x,''                     should contain rotamers.'')')
      write(*,'(6x,''↳ -maxcomp <int>  : Selcect the lowest <int> conformers'')')
      write(*,'(8x,''                  out of each ensemble to be compared'')')
      write(*,'(8x,''                  with "-compare". [default: 10]'')')
      write(*,'(5x,''-testtopo <file>   : Analyze some stuctural info (topology) for a given file.'')')
      write(*,'(5x,''-constrain <atoms> : write example file ".xcontrol.sample" for constraints'')')
      write(*,'(5x,''                     in crest. (see -cinp option above)'')')
      write(*,'(5x,''-thermo <file>     : Calculate thermo data for given structure. Also requires vibrational'')')
      write(*,'(5x,''                     frequencies in the TM format, saved as file called "vibspectrum"'')')
      write(*,'(5x,''-rmsd,-rmsdheavy <file1> <file2>  : Calculate RMSD or heavy atom RMSD between two structures.'')')
      write(*,'(5x,''                                    Input coords are automatically transformed to Angstroem.'')')      
      write(*,'(5x,''-splitfile <file> [from] [to]     : Split an ensemble from <file> into seperate directories'')')
      write(*,'(5x,''                                    for each structure. [from] and [to] can be used to select'')')      
      write(*,'(5x,''                                    specific structures from the file.'')')      
      write(*,'(5x,''                                    The new directories are collected in the SPLIT directory.'')')      

  write(*,*)
  write(*,*)
      write(*,*) 'View literature references with [--cite]'
      write(*,*) 'Also refer to:'
      write(*,*) 'https://xtb-docs.readthedocs.io/en/latest/crest.html'
      write(*,*)

      stop '   [-h] displayed. exit.'
end subroutine confscript_help

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


subroutine crestcite
     write(*,*)
     write(*,'(4x,''MAIN REFERENCES:'')')
     write(*,'(/,7x,''P. Pracht, F. Bohle, S. Grimme,'')')
     write(*,'(  7x,''PCCP, 2020, 22, 7169-7192.'')')
     write(*,'(/,7x,''S. Grimme, JCTC, 2019, 15, 2847-2862.'')')
     write(*,'(/7x,''P. Pracht, S. Grimme, Chem. Sci., 2021,'')')
     write(*,'(7x,''DOI: 10.1039/d1sc00621e'')')
     write(*,'(/,/)')
     write(*,'(4x,''GFNn-xTB references:'')')
     write(*,'(6x,''GFN1-xTB'')')
     write(*,'(7x,''S. Grimme, C. Bannwarth, P. Shushkov, JCTC, 2017,'')')
     write(*,'(7x,''13, 1989-2009. DOI: 10.1021/acs.jctc.7b00118'')')
     write(*,'(6x,''GFN2-xTB'')')
     write(*,'(7x,''C. Bannwarth, S. Ehlert and S. Grimme., JCTC, 2019,'')')
     write(*,'(7x,''15, 1652-1671. DOI: 10.1021/acs.jctc.8b01176'')')
     write(*,'(6x,''GFN0-xTB'')')
     write(*,'(7x,''P. Pracht, E. Caldeweyher, S. Ehlert, S. Grimme, 2019,'')')
     write(*,'(7x,''ChemRxiv preprint, DOI: 10.26434/chemrxiv.8326202.v1'')')
     write(*,'(6x,''GFN-FF'')')
     write(*,'(7x,''S. Spicher, S. Grimme, Angew. Chem. Int. Ed., 2020,'')')
     write(*,'(7x,''132, 2-11, DOI: 10.1002/ange.202004239'')')
    
     write(*,'(/,/)')
     write(*,'(4x,''related references:'')')
     write(*,'(7x, ''S. Grimme, C. Bannwarth, S. Dohm, A. Hansen,'')')
     write(*,'(7x, ''J. Pisarek, P. Pracht, J. Seibert, F. Neese,'')')
     write(*,'(7x, ''Angew. Chem. Int. Ed., 2017, 56, 14763-14769'')')
     write(*,'(/7x,''P. Pracht, C.A. Bauer, S. Grimme, JCC, 2017,'')')
     write(*,'(7x, ''38, 2618–2631, DOI: 10.1002/jcc.24922'')')
     write(*,'(/7x,''P. Pracht, R. Wilcken, A. Udvarhelyi, S. Rodde, S. Grimme,'')')
     write(*,'(7x, ''JCAMD, 2018, 32, 1139-1149, DOI: 10.1007/s10822-018-0145-7'')')
 

     write(*,'(/,/)')
     write(*,'(3x,''Please cite work conducted with this code appropriately.'')')
     stop '   [--cite] displayed. exit.'
end subroutine crestcite




!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine crestcrest
     write(*,'(7x,''|                                            |'')')    
     write(*,'(7x,''|              ###############               |'')')
     write(*,'(7x,''|              #//////|@@@@@@#               |'')')
     write(*,'(7x,''|              #\\\\\\|@@@@@@#               |'')')
     write(*,'(7x,''|              #//////|@@@@@@#               |'')')
     write(*,'(7x,''|              #@@@@@@|\\\\\\#               |'')')
     write(*,'(7x,''|              #@@@@@@|//////#               |'')')
     write(*,'(7x,''|               #@@@@@|\\\\\#                |'')') 
     write(*,'(7x,''|                ###@@|//###                 |'')') 
     write(*,'(7x,''|                   #####                    |'')') 
     write(*,'(7x,''|                                            |'')')
     write(*,'(7x,''|                 C R E S T                  |'')')
end subroutine crestcrest


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine prchd
      write(*,*)
      write(*,'(7x,''========================================'')')
      write(*,'(7x,''|             C R E G E N              |'')')
      write(*,'(7x,''|     conformer/rotamer generation     |'')')
      write(*,'(7x,''| & NMR symmetry/equivalence analysis  |'')')
      write(*,'(7x,''|      SG, Universitaet Bonn, MCTC     |'')')
      write(*,'(7x,''|     Fri Aug 11 13:16:16 CEST 2017    |'')')
      write(*,'(7x,''========================================'')')
      write(*,*)
end
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


subroutine header_stereo
      write(*,*)
      write(*,'(5x,''========================================'')')
      write(*,'(5x,''|   automated stereoisomer generator   |'')')
      write(*,'(5x,''|              P.Pracht                |'')')
      write(*,'(5x,''|       Universitaet Bonn, MCTC        |'')')
      write(*,'(5x,''|    Wed 9. Oct 11:20:34 CEST 2019     |'')')
      write(*,'(5x,''========================================'')')
      write(*,*)

      write(*,*) 'NOTE: This is a work-in-progress project!'
end subroutine header_stereo


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


subroutine prreactorhd
      write(*,*)
      write(*,'(7x,''========================================'')')
      write(*,'(7x,''|          GFNn-xTB NANOREACTOR        |'')')
      write(*,'(7x,''|      SG, Universitaet Bonn, MCTC     |'')')
      write(*,'(7x,''========================================'')')
      write(*,'(/,7x,''JCTC, 2019, 15, 2847-2862.'')')
      write(*,*)
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine zsortwarning2(env)
      use crest_data
      implicit none
      type(systemdata) :: env
      logical :: ex
      inquire(file=env%constraints,exist=ex)
      if(ex.and.env%autozsort)then
      write(*,*)'==========================================='
      write(*,*) 'WARNING:'
      write(*,*) 'The input coordinate file would be sorted'
      write(*,*) 'by zsort and a constraining file is'
      write(*,*) 'present. To avoid constrainment of the'
      write(*,*) 'wrong atoms zsort will be turned off.'
      write(*,*) 'This also might influence the results.'
      write(*,*)'==========================================='
      write(*,*)
      env%autozsort=.false.
      endif
end subroutine zsortwarning2

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine qcg_head()
  implicit none
  write(*,*)
  write(*,'(2x,''========================================'')')
  write(*,'(2x,''|           ----------------           |'')')
  write(*,'(2x,''|                 Q C G                |'')')
  write(*,'(2x,''|           ----------------           |'')')
  write(*,'(2x,''|        Quantum Cluster Growth        |'')')
  write(*,'(2x,''|       University of Bonn, MCTC       |'')')
  write(*,'(2x,''========================================'')')
  write(*,'(2x,'' S. Grimme, S. Spicher, unpublished.'')')
  write(*,*)
end subroutine qcg_head

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!convert a string in a small header printout
subroutine smallhead(str)
      implicit none
      character(len=*) :: str
      integer :: strlen
      integer :: i,j
      character(len=:),allocatable :: str2
      strlen=len_trim(str)
      str2=repeat('-',strlen) 
      write(*,'(1x,a)') trim(str2)
      write(*,'(1x,a)') trim(str)
      write(*,'(1x,a)') trim(str2)
      return
end subroutine smallhead
subroutine smallheadline(line)
      implicit none
      character(len=*) :: line
      integer :: lw
      lw=len_trim(line)
      write(*,'(/,1x,a)') repeat('=',lw)
      write(*,'(1x,a)') trim(line)
      write(*,'(1x,a)') repeat('=',lw)
      return
end subroutine smallheadline
subroutine underline(str)
      implicit none
      character(len=*) :: str
      integer :: strlen
      integer :: i,j
      character(len=:),allocatable :: str2
      strlen=len_trim(str)
      str2=repeat('-',strlen) 
      write(*,'(1x,a)') trim(str)
      write(*,'(1x,a)') trim(str2)
      return
end subroutine underline

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! a warning when the MTD length exceeds 200ps
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

subroutine mtdwarning(lenv)
        implicit none
        real*8 :: lenv

       write(*,*)
       write(*,'(a,f5.1,a)')"! WARNING: the estimated MTD time exceeds ",lenv," ps."
       write(*,'(a       )')"! Because the estimate is uncertain, the program restricts"
       write(*,'(a,f5.1,a)')"! this to ",lenv," ps and continues. The user may"
       write(*,'(a       )')"! re-run crest with manual setting by '-mdlen <time>' and"
       write(*,'(a       )')"! check the results carefully." 
       write(*,*)

end subroutine mtdwarning

!--------------------------------------------------------------------------------------------
! iteration cycler printout
!--------------------------------------------------------------------------------------------
subroutine printiter
      implicit none
      write(*,*)
      write(*,'(''*******************************************************************************************'')')
      write(*,'(''**                        N E W    I T E R A T I O N    C Y C L E                        **'')')
      write(*,'(''*******************************************************************************************'')')
end subroutine printiter
!--------------------------------------------------------------------------------------------
subroutine printiter2(i)
      implicit none
      integer :: i
      write(*,*)
      write(*,'(90("*"))')
      write(*,'("**",26x,"I T E R A T I O N    C Y C L E    ",i3,23x,"**")')i
      write(*,'(90("*"))')
end subroutine printiter2

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!convert a string in a large header printout
subroutine largehead(str)
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
      strlen=len_trim(str)
      str2=repeat(symb,90)
      strlen2=len_trim(str2)
      if(strlen.ge.strlen2)then
        str2=''
        do i=1,strlen+6
          str2=trim(str2)//symb
        enddo
        strlen=strlen+6
      endif
      k=strlen2-strlen
      j= k/2 - 2
      str3=symb//symb
      str3(j+1:)=trim(str)
      str3(strlen2-1:strlen2)=symb//symb
      write(*,*)
      write(*,'(a)') trim(str2)
      write(*,'(a)') trim(str3)
      write(*,'(a)') trim(str2)
      return
end subroutine construct_large_headline

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Confscript dry-run printout
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine crest_dry(env)
      use iso_fortran_env, wp => real64
      use crest_data
      use iomod
      implicit none

      type(systemdata),intent(inout) :: env
      !type(options),intent(inout)    :: opt
      character(len=1024) :: dumstr
      character(len=:),allocatable :: dum
      logical :: cregenpr = .false.
      logical :: mdsetpr  = .false.
      logical :: constpr  = .false.
      logical :: jobpr    = .true.
      logical :: xtbpr    = .true.
      logical :: techpr   = .true.

      call largehead('D R Y    R U N')
      write(*,'(1x,a)') 'Dry run was requested.'
      write(*,'(1x,a)') 'Running CREST with the chosen cmd arguments will result in the following settings:'

      write(*,'(/,1x,a,a)') 'Input file : ',trim(env%inputcoords)

      write(*,'(/,1x,a)') 'Job type :'
      if(env%onlyZsort)then
         write(*,'(2x,a)',advance='no')'1.'
         write(*,'(2x,a)')'Standalone use of ZSORT routine.'
         xtbpr=.false.
         techpr=.false.
         jobpr=.false.
      else if( env%properties.lt.0 )then
         jobpr=.false.
         write(*,'(2x,a)',advance='no')'1.'
         select case( env%properties )
         case( -1 )
             write(*,'(2x,a)')'Standalone use of CREGEN sorting routine.'
             cregenpr=.true.
             xtbpr=.false.
         case( -2 )
             write(*,'(2x,a)')'Comparison of two conformer-rotamer ensembles'
             cregenpr=.true.
             xtbpr=.false.
         case( -3 )
           write(*,'(2x,a)')'Automated protonation'
         case( -4 )
           write(*,'(2x,a)')'Automated deprotonation'
         case( -5 )
           write(*,'(2x,a)')'Automated tautomerization'
         case( -666 )
            write(*,'(2x,a)')'"Property" calculation (-prop) for a given ensemble.'
         case default
            write(*,'(2x,a)') '<undefined>'
         end select
      else
         cregenpr=.true.
         write(*,'(2x,a)',advance='no')'1.'
         select case( env%crestver )
           case( 1 )
              write(*,'(2x,a)')'Conformational search via the MF-MD-GC algo'
           case( 2 )
              if(env%properties==45)then
              write(*,'(2x,a)')'Conformational search via the iMTD-sMTD algo'
              else if(env%iterativeV2)then
              write(*,'(2x,a)')'Conformational search via the iMTD-GC algo'
              else
              write(*,'(2x,a)')'Conformational search via the MTD-GC algo'
              endif
              mdsetpr=.true.
           case( crest_imtd2 )  
              write(*,'(2x,a)')'Conformational search cia the iMTD-sMTD algo (-v4)'
              mdsetpr=.true.
           case( 3 )
              write(*,'(2x,a)')'Reoptimization of all structures in a given ensemble (-mdopt)'
           case( 4 )
              write(*,'(2x,a)')'Reoptimization and sorting of all structures in a given ensemble (-screen)'
           case( 7 )
              write(*,'(2x,a)')'GFNn-xTB nano reactor (-reactor)'  
           case( crest_pka )
              write(*,'(2x,a)')'GFN2-xTB/ALPB(H2O) pKa calculation (-pka)'  
           case default
              write(*,'(2x,a)') '<undefined>'
         end select
         if( env%properties.gt.0 )then
             select case( env%properties )
             case( 45 )
               write(*,'(2x,a,2x,a)')'2.','Calculation of molecular ensemble entropy (-entropy)'
             case default    
               write(*,'(2x,a,2x,a)')'2.','Additional "property" calculation requested (-prop)'
             end select
         endif
      endif


      if(jobpr)then
      write(*,'(/,1x,a)') 'Job settings'
      write(*,'(2x,a,l6)') 'sort Z-matrix        : ',env%autozsort
      if(env%crestver .eq. 2)then
      select case( env%runver )
       case( 2 )
      write(*,'(2x,a,a)')  'MTD-GC modified mode : ','"-quick"'
       case( 4 )
      write(*,'(2x,a,a)')  'MTD-GC modified mode : ','"-nci"'
       case( 5 )
      write(*,'(2x,a,a)')  'MTD-GC modified mode : ','"-squick"'
       case( 6 )
      write(*,'(2x,a,a)')  'MTD-GC modified mode : ','"-mquick"'
       case( 111  )
      write(*,'(2x,a,a)')  'iMTD-sMTD mode       : ','"-entropy"'
       case default
        continue
      end select
      endif
      if( env%properties.gt.0 )then
      select case( env%properties )
      case( 1 )
        dum='"hess"'
      case( 10 )
        dum='"ohess"'
      case( 2 )
        dum='"autoIR (GFN)"'
      case( 20 )
        dum='"reopt"'
      case(3:7,100 )
        dum='DFT'  
      case( 45 )
        dum='none'  
      case default
        dum='<undefined>'
      end select
      if(dum.ne.'none')then
      write(*,'(2x,a,a)')  'PROP mode (-prop)    : ',dum
      endif
      endif
      endif

      if(cregenpr)then
      write(*,'(/,1x,a)') 'CRE settings'
      write(*,'(2x,a,f10.4)')'energy window         (-ewin) :',env%ewin
      write(*,'(2x,a,f10.4)')'RMSD threshold        (-rthr) :',env%rthr          !RTHR - RMSD thr in Angstroem
      write(*,'(2x,a,f10.4)')'energy threshold      (-ethr) :',env%ethr          !ETHR - E threshold in kcal
      write(*,'(2x,a,f10.2)')'rot. const. threshold (-bthr) :',env%bthr2         !BTHR - rot const  thr    
      write(*,'(2x,a,f10.2)')'T (for boltz. weight) (-temp) :',env%tboltz  
      endif


      if(mdsetpr)then
      write(*,'(/,1x,a)') 'General MD/MTD settings'
      if(env%mdtime.gt.0.0d0)then
      write(*,'(2x,a,f10.1)')'simulation length [ps]    (-len) :',env%mdtime
      else
      write(*,'(2x,a,a)')    'simulation length [ps]    (-len) : ','<system dependent>'
      endif
      write(*,'(2x,a,f10.1)')'time step [fs]          (-tstep) :',env%mdstep
      write(*,'(2x,a,i10)')  'shake mode              (-shake) :',env%shake
      write(*,'(2x,a,f10.2)')'MTD temperature [K]    (-mdtemp) :',env%mdtemp
      write(*,'(2x,a,i10)')  'trj dump step  [fs]    (-mddump) :',env%mddumpxyz
      write(*,'(2x,a,f10.1)')'MTD Vbias dump [ps]    (-vbdump) :',real(env%mddump)/1000.0d0
      endif

      if(env%cts%used)then
      write(*,'(/,1x,a)') 'Constrainment info'
      write(*,'(2x,a,l7)')  'applying constraints?  : ',env%cts%used
      write(*,'(2x,a,a)')   'constraining file      : ',trim(env%constraints)
      write(*,'(2x,a)')     'file content :'
      call cat_mod(6,'  > ',env%constraints,'')
      endif

      if(xtbpr)then
      write(*,'(/,1x,a)') 'XTB settings'
      write(*,'(2x,a,a)')  'binary name        (-xnam) : ',trim(env%ProgName)
      call checkbinary(env)
      write(*,'(2x,a,a)')  'GFN method         (-gfn)  : ',trim(env%gfnver)
      write(*,'(2x,a,i0)') '(final) opt level  (-opt)  : ',nint(env%optlev)
      if(env%gbsa)then
        if(index(env%solv,'--alpb').ne.0)then
      write(*,'(2x,a,a)')  'Implicit solvation (-alpb) : ',trim(env%solvent)
        else
      write(*,'(2x,a,a)')  'Implicit solvation (-gbsa) : ',trim(env%solvent)
        endif
      endif
      if(env%chrg.ne.0.0d0)then
      write(*,'(2x,a,i0)') 'Molecular charge   (-chrg) : ',env%chrg
      endif
      if(env%uhf.ne.0)then
      write(*,'(2x,a,i0)') 'UHF (nα-nβ elec.)  (-uhf)  : ',env%uhf
      endif
      endif

      if(techpr)then
      call getcwd(dumstr)
      write(*,'(/,1x,a)') 'Technical settings'
      write(*,'(2x,a,a)')  'working directory : ',trim(dumstr)
      write(*,'(2x,a,i0)') 'CPUs (threads)     (-T) : ',env%threads

      endif 
      
      write(*,'(/)')
      stop 'normal dry run termination.'
end subroutine crest_dry

subroutine cat_mod(ch,pre,fname,post)
       implicit none
       integer :: ch
       character(len=*) :: pre
       character(len=*) :: fname
       character(len=*) :: post
       character(len=256) :: adum
       integer :: ich,i,j,k,l,io
           
       open(newunit=ich,file=fname)
       do
         read(ich,'(a)',iostat=io) adum
         if(io<0)exit
         write(ch,'(a,a,a)')pre,trim(adum),post
       enddo
       close(ich)
       return
end subroutine cat_mod

subroutine checkbinary(env)
        use crest_data
        use syscheck
        implicit none
        type(systemdata) :: env
        !type(options) :: opt
        integer :: r

        r = 0
        call checkprog(trim(env%ProgName),r)
  
        if(r .ne. 0) then
        write(*,'(4x,a)') 'Warning! The xtb binary was not found and hence CREST might crash'
        endif

        if(env%crestver .eq. crest_solv) then
            call checkprog(trim('xtbiff'),r)
            if(r .ne. 0)then
                write(*,'(4x,a)') 'Warning! The xtbiff binary was not found and hence the qcg mode in CREST will probably crash'
            end if
        end if
        return
end subroutine checkbinary

!==============================================================================!
!  QCG-printouts
!==============================================================================!

!____________________________________________________________________________

subroutine print_qcg_grow()
  implicit none
  write(*,*)
  write(*,'(2x,''========================================='')')
  write(*,'(2x,''|   quantum cluster growth: GROW        |'')')
  write(*,'(2x,''========================================='')')
  write(*,*)
end subroutine print_qcg_grow

!____________________________________________________________________________

subroutine pr_qcg_fastgrow()
  implicit none
  write(*,*)
  write(*,'(2x,''========================================='')')
  write(*,'(2x,''|   quantum cluster growth: FASTGROW    |'')')
  write(*,'(2x,''========================================='')')
  write(*,*)
end subroutine pr_qcg_fastgrow

!____________________________________________________________________________

subroutine print_qcg_ensemble()
  implicit none
  write(*,*)
  write(*,'(2x,''========================================='')')
  write(*,'(2x,''|   quantum cluster growth: ENSEMBLE    |'')')
  write(*,'(2x,''========================================='')')
  write(*,*)
end subroutine print_qcg_ensemble

!____________________________________________________________________________

subroutine print_qcg_opt()
  implicit none
  write(*,*)
  write(*,'(2x,''========================================='')')
  write(*,'(2x,''|   quantum cluster growth: OPT         |'')')
  write(*,'(2x,''========================================='')')
  write(*,*)
  write(*,'(2x,''Very tight post optimization of lowest cluster'')')
end subroutine print_qcg_opt

!____________________________________________________________________________

subroutine pr_qcg_fill()
  implicit none
  write(*,*)
  write(*,'(2x,''========================================='')')
  write(*,'(2x,''|   quantum cluster growth: CFF         |'')')
  write(*,'(2x,''========================================='')')
  write(*,*)
  write(*,'(2x,''CUT-FREEZE-FILL Algorithm to generate reference solvent cluster'')')
end subroutine pr_qcg_fill

!____________________________________________________________________________

subroutine pr_qcg_freq()
  implicit none
  write(*,*)
  write(*,'(2x,''========================================='')')
  write(*,'(2x,''|          Frequency evaluation         |'')')
  write(*,'(2x,''========================================='')')
  write(*,*)
end subroutine pr_qcg_freq

!____________________________________________________________________________

subroutine pr_eval_solute()
  implicit none
  write(*,*)
  write(*,'(2x,''________________________________________________________________________'')')
  write(*,*)
  write(*,'(2x,''__________________     Solute Cluster Generation   _____________________'')')
  write(*,*)
  write(*,'(2x,''________________________________________________________________________'')')
  write(*,*)
end subroutine pr_eval_solute

!____________________________________________________________________________

subroutine pr_eval_solvent()
  implicit none
  write(*,*)
  write(*,*)
  write(*,'(2x,''________________________________________________________________________'')')
  write(*,*)
  write(*,'(2x,''_________________     Solvent Cluster Generation   _____________________'')')
  write(*,*)
  write(*,'(2x,''________________________________________________________________________'')')
  write(*,*)
end subroutine pr_eval_solvent

!____________________________________________________________________________

subroutine pr_eval_eval()
  implicit none
  write(*,*)
  write(*,*)
  write(*,'(2x,''________________________________________________________________________'')')
  write(*,*)
  write(*,'(2x,''_________________________     Evaluation    ____________________________'')')
  write(*,*)
  write(*,'(2x,''________________________________________________________________________'')')
  write(*,*)
  write(*,*)
end subroutine pr_eval_eval

!____________________________________________________________________________

subroutine pr_freq_energy()
  implicit none
  write(*,'(2x,"#       H(T)       SVIB      SROT       STRA      G(T)")')
  write(*,'(2x,"     [kcal/mol]    [      cal/mol/K        ]    [kcal/mol]")')
  write(*,'(2x,"--------------------------------------------------------")')
end subroutine pr_freq_energy  

!____________________________________________________________________________

subroutine pr_eval_1(G,H)
  use iso_fortran_env, only : wp => real64
  implicit none
  real(wp),intent(in)       :: G,H
  write(*,'(2x,"-----------------------------------------------------")')
  write(*,'(2x,"Gsolv and Hsolv ref. state: [1 M gas/solution] ")')
  write(*,'(2x,"G_solv (incl.RRHO)      =",F8.2," kcal/mol")') G
  write(*,'(2x,"H_solv (incl.RRHO)      =",F8.2," kcal/mol")') H
  write(*,'(2x,"-----------------------------------------------------")')
  write(*,*)
end subroutine pr_eval_1  


!____________________________________________________________________________

subroutine pr_eval_2(srange,G,scal)
  use iso_fortran_env, only : wp => real64
  implicit none
! Dummy  
  integer,intent(in)        :: srange
  real(wp),intent(in)       :: G(srange)
  real(wp),intent(in)       :: scal(srange)
! Stack
  integer                   :: i
  write(*,'(2x,"-----------------------------------------------------")')
  write(*,'(2x,"Solvation free energies with scaled translational")')
  write(*,'(2x,"and rotational degrees of freedom: Gsolv (scaling)")')
  do i = 1,srange
    write(*,'(10x,">>",2x,f8.2," (",f4.2,")",4x,"<<")') G(i),scal(i)
  end do
  write(*,'(2x,"-----------------------------------------------------")')
end subroutine pr_eval_2

!____________________________________________________________________________

subroutine pr_eval_3(srange,freqscal,scal,G)
  use iso_fortran_env, only : wp => real64
  implicit none
! Dummy  
  integer,intent(in)        :: srange
  integer,intent(in)        :: freqscal
  real(wp),intent(in)       :: scal
  real(wp),intent(in)       :: G(srange)
  write(*,*)
  write(*,'(2x,"==================================================")')
  write(*,'(2x,"|  Gsolv with SCALED RRHO contributions: ",f4.2,4x"|")') scal
  write(*,'(2x,"|  [1 bar gas/ 1 M solution]                     |")')
  write(*,'(2x,"|                                                |")')
  write(*,'(2x,"|  G_solv (incl.RRHO)+dV(T)=",F8.2," kcal/mol    |")') G(freqscal)
  write(*,'(2x,"==================================================")')
  write(*,*)
end subroutine pr_eval_3

!____________________________________________________________________________

subroutine pr_fill_energy()
  implicit none
  write(*,'(x,'' Size'',2x,''Cluster '',2x,''E /Eh '',7x,''De/kcal'',3x,&
          &''Detot/kcal'',2x,''Opt'',4x)')
end subroutine pr_fill_energy 

!____________________________________________________________________________

subroutine pr_ensemble_energy()
  implicit none
  write(*,*)
  write(*,'(x,'' Cluster'',3x,''E /Eh '',7x,&
           &''Density'',2x,''Efix'',7x,''R   av/act.'',1x,&
           &''Surface'',3x,''Opt'',4x)')
end subroutine pr_ensemble_energy

!____________________________________________________________________________

subroutine pr_qcg_esolv()
  implicit none
  write(*,*)
  write(*,'(2x,''========================================='')')
  write(*,'(2x,''|   quantum cluster growth: ESOLV       |'')')
  write(*,'(2x,''|                                       |'')')
end subroutine pr_qcg_esolv

!____________________________________________________________________________

subroutine pr_grow_energy()
  implicit none
  write(*,'(x,'' Size'',2x,''E /Eh '',6x,''De/kcal'',3x,''Detot/kcal'',2x,&
           &''Density'',3x,''Efix'',9x,''R   av/act.'',1x,&
           &''Surface'',3x,''Opt'',4x)')
end subroutine pr_grow_energy


!==============================================================================!
!  printout percent calculation for GUI mode
!==============================================================================!
subroutine wrGUIpercent(current,maxv,interval)
       use iso_fortran_env, wp => real64
       implicit none
       integer :: current,maxv,interval
       real(wp) :: perc,inc,trc
       integer :: i,interval2

       if(current==maxv)then
         write(*,'(1x,f6.2,a)') 100.0d0,' percent done. finished loop.'
         return
       endif

       interval2 = interval
       inc=float(maxv)/float(interval)
       if(inc.gt.1.0d0) interval2=maxv
       trc=0.0d0
       do i=1,interval2
         trc = floor(float(i) * inc)
         if(current == nint(trc))then
            perc = (float(current) / float(maxv)) *100.0d0
            write(*,'(1x,f6.2,a)') perc,' percent done'
            exit
         endif    
       enddo

       return      
end subroutine wrGUIpercent

