!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021-2025 Christoph Plett, Philipp Pracht
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

module qcg_printouts
  use crest_parameters,only:stdout,wp

  use crest_data
  use iomod
  implicit none
  public

contains

  subroutine qcg_head()
    implicit none
    write (stdout,*)
    write (stdout,'(10x,''========================================'')')
    write (stdout,'(10x,''|           ----------------           |'')')
    write (stdout,'(10x,''|                 Q C G                |'')')
    write (stdout,'(10x,''|           ----------------           |'')')
    write (stdout,'(10x,''|        Quantum Cluster Growth        |'')')
    write (stdout,'(10x,''|       University of Bonn, MCTC       |'')')
    write (stdout,'(10x,''========================================'')')
    write (stdout,'(10x,''    S. Grimme, S. Spicher, C. Plett.'')')
    write (stdout,*)
    write (stdout,'(10x,''Cite work conducted with this code as:'')')
    write (stdout,'(/,9x,''S. Spicher, C. Plett, P. Pracht, A. Hansen,'')')
    write (stdout,'(9x,''S. Grimme, JCTC, 2022, 18, 3174-3189.'')')
    write (stdout,*)
  end subroutine qcg_head

  subroutine write_qcg_setup(env)
    implicit none
    type(systemdata) :: env

    write (stdout,*)
    write (stdout,'(2x,''========================================='')')
    write (stdout,'(2x,''|   quantum cluster growth: INPUT       |'')')
    write (stdout,'(2x,''========================================='')')
    write (stdout,*)
    select case (env%qcg_runtype)
    case (0)
      write (stdout,'(2x,''QCG: Only Cluster Generation'')')
    case (1)
      write (stdout,'(2x,''QCG: Cluster + Ensemble Generation'')')
      if (env%ensemble_method .eq. 0) write (stdout,'(2x,''Ensemble generated via CREST'')')
      if (env%ensemble_method .eq. 1) write (stdout,'(2x,''Ensemble generated via MD Simulation'')')
      if (env%ensemble_method .eq. 2) write (stdout,'(2x,''Ensemble generated via MetaDynamic'')')
    case (2)
      write (stdout,'(2x,''QCG: Calculation of delta E_solv'')')
      if (env%ensemble_method .eq. 0) write (stdout,'(2x,''Ensemble generated via CREST'')')
      if (env%ensemble_method .eq. 1) write (stdout,'(2x,''Ensemble generated via MD Simulation'')')
      if (env%ensemble_method .eq. 2) write (stdout,'(2x,''Ensemble generated via MetaDynamic'')')
    case (3)
      write (stdout,'(2x,''QCG: Calculation of delta G_solv'')')
      if (env%ensemble_method .eq. 0) write (stdout,'(2x,''Ensemble generated via CREST'')')
      if (env%ensemble_method .eq. 1) write (stdout,'(2x,''Ensemble generated via MD Simulation'')')
      if (env%ensemble_method .eq. 2) write (stdout,'(2x,''Ensemble generated via MetaDynamic'')')
    case default
      continue
    end select
    write (stdout,*)
    write (stdout,'(2x,''input parameters     '')')
    write (stdout,'(2x,''solute                 : '',a)') trim(env%solu_file)
    write (stdout,'(2x,''charge                 : '',i0)') env%chrg
    write (stdout,'(2x,''uhf                    : '',i0)') env%uhf
    write (stdout,'(2x,''solvent                : '',a)') trim(env%solv_file)
    if (env%nsolv .ne. 0) then
      write (stdout,'(2x,''# of solvents to add   : '',i0)') env%nsolv
    else if (env%nsolv .eq. 0) then
      write (stdout,'(2x,''# of solvents to add   : until convergence, but maximal'',1x,i4)') env%max_solv
    end if
    if (env%nqcgclust .ne. 0) then
      write (stdout,'(2x,''# of cluster generated : '',i0)') env%nqcgclust
    else
      write (stdout,'(2x,''Cluster generated that are above 10 % populated '')')
    end if

    write (stdout,'(2x,''# of CPUs used         : '',i0)') env%Threads
    if (env%solvent .eq. '') then
      write (stdout,'(2x,''No gbsa/alpb model''  )')
    else
      write (stdout,'(2x,''Solvation model        : '',a)') env%solvent
    end if
    write (stdout,'(2x,''xtb opt level          : '',a)') trim(optlevflag(env%optlev))
    write (stdout,'(2x,''System temperature [K] : '',F5.1)') env%tboltz
    write (stdout,'(2x,''RRHO scaling factor    : '',F4.2)') env%freq_scal
    write (stdout,*)

  end subroutine write_qcg_setup

!========================================================================================!
!========================================================================================!
!>  QCG-printouts
!==============================================================================!
!========================================================================================!

  subroutine print_qcg_grow()
    implicit none
    write (stdout,*)
    write (stdout,'(2x,''========================================='')')
    write (stdout,'(2x,''|   quantum cluster growth: GROW        |'')')
    write (stdout,'(2x,''========================================='')')
    write (stdout,*)
  end subroutine print_qcg_grow
  subroutine pr_qcg_fastgrow()
    implicit none
    write (stdout,*)
    write (stdout,'(2x,''========================================='')')
    write (stdout,'(2x,''|   quantum cluster growth: FASTGROW    |'')')
    write (stdout,'(2x,''========================================='')')
    write (stdout,*)
  end subroutine pr_qcg_fastgrow
  subroutine print_qcg_ensemble()
    implicit none
    write (stdout,*)
    write (stdout,'(2x,''========================================='')')
    write (stdout,'(2x,''|   quantum cluster growth: ENSEMBLE    |'')')
    write (stdout,'(2x,''========================================='')')
    write (stdout,*)
  end subroutine print_qcg_ensemble
  subroutine print_qcg_opt()
    implicit none
    write (stdout,*)
    write (stdout,'(2x,''========================================='')')
    write (stdout,'(2x,''|   quantum cluster growth: OPT         |'')')
    write (stdout,'(2x,''========================================='')')
    write (stdout,*)
    write (stdout,'(2x,''Very tight post optimization of lowest cluster'')')
  end subroutine print_qcg_opt
  subroutine pr_qcg_fill()
    implicit none
    write (stdout,*)
    write (stdout,'(2x,''========================================='')')
    write (stdout,'(2x,''|   quantum cluster growth: CFF         |'')')
    write (stdout,'(2x,''========================================='')')
    write (stdout,*)
    write (stdout,'(2x,''CUT-FREEZE-FILL Algorithm to generate reference solvent cluster'')')
  end subroutine pr_qcg_fill
  subroutine pr_qcg_freq()
    implicit none
    write (stdout,*)
    write (stdout,'(2x,''========================================='')')
    write (stdout,'(2x,''|          Frequency evaluation         |'')')
    write (stdout,'(2x,''========================================='')')
    write (stdout,*)
  end subroutine pr_qcg_freq
  subroutine pr_eval_solute()
    implicit none
    write (stdout,*)
    write (stdout,'(2x,''________________________________________________________________________'')')
    write (stdout,*)
    write (stdout,'(2x,''__________________     Solute Cluster Generation   _____________________'')')
    write (stdout,*)
    write (stdout,'(2x,''________________________________________________________________________'')')
    write (stdout,*)
  end subroutine pr_eval_solute
  subroutine pr_eval_solvent()
    implicit none
    write (stdout,*)
    write (stdout,*)
    write (stdout,'(2x,''________________________________________________________________________'')')
    write (stdout,*)
    write (stdout,'(2x,''_________________     Solvent Cluster Generation   _____________________'')')
    write (stdout,*)
    write (stdout,'(2x,''________________________________________________________________________'')')
    write (stdout,*)
  end subroutine pr_eval_solvent
  subroutine pr_eval_eval()
    implicit none
    write (stdout,*)
    write (stdout,*)
    write (stdout,'(2x,''________________________________________________________________________'')')
    write (stdout,*)
    write (stdout,'(2x,''_________________________     Evaluation    ____________________________'')')
    write (stdout,*)
    write (stdout,'(2x,''________________________________________________________________________'')')
    write (stdout,*)
    write (stdout,*)
  end subroutine pr_eval_eval
  subroutine pr_freq_energy()
    implicit none
    write (stdout,'(2x,"#       H(T)       SVIB      SROT       STRA      G(T)")')
    write (stdout,'(2x,"     [kcal/mol]    [      cal/mol/K        ]    [kcal/mol]")')
    write (stdout,'(2x,"--------------------------------------------------------")')
  end subroutine pr_freq_energy
  subroutine pr_eval_1(G,H)
    implicit none
    real(wp),intent(in)       :: G,H
    write (stdout,'(2x,"-----------------------------------------------------")')
    write (stdout,'(2x,"Gsolv and Hsolv ref. state: [1 M gas/solution] ")')
    write (stdout,'(2x,"G_solv (incl.RRHO)      =",F12.2," kcal/mol")') G
    write (stdout,'(2x,"H_solv (incl.RRHO)      =",F12.2," kcal/mol")') H
    write (stdout,'(2x,"-----------------------------------------------------")')
    write (stdout,*)
  end subroutine pr_eval_1
  subroutine pr_eval_2(srange,G,scal)
    implicit none
! Dummy
    integer,intent(in)        :: srange
    real(wp),intent(in)       :: G(srange)
    real(wp),intent(in)       :: scal(srange)
! Stack
    integer                   :: i
    write (stdout,'(2x,"-----------------------------------------------------")')
    write (stdout,'(2x,"Solvation free energies with scaled translational")')
    write (stdout,'(2x,"and rotational degrees of freedom: Gsolv (scaling)")')
    do i = 1,srange
      write (stdout,'(10x,">>",f10.2," (",f4.2,")",4x,"<<")') G(i),scal(i)
    end do
    write (stdout,'(2x,"-----------------------------------------------------")')
  end subroutine pr_eval_2
  subroutine pr_eval_3(srange,freqscal,scal,G)
    implicit none
! Dummy
    integer,intent(in)        :: srange
    integer,intent(in)        :: freqscal
    real(wp),intent(in)       :: scal
    real(wp),intent(in)       :: G(srange)
    write (stdout,*)
    write (stdout,'(2x,"==================================================")')
    write (stdout,'(2x,"|  Gsolv with SCALED RRHO contributions: ",f4.2,4x"|")') scal
    write (stdout,'(2x,"|  [1 bar gas/ 1 M solution]                     |")')
    write (stdout,'(2x,"|                                                |")')
    write (stdout,'(2x,"|  G_solv (incl.RRHO)+dV(T)=",F10.2," kcal/mol  |")') G(freqscal)
    write (stdout,'(2x,"==================================================")')
    write (stdout,*)
  end subroutine pr_eval_3
  subroutine pr_fill_energy()
    implicit none
    write (stdout,'(x,'' Size'',2x,''Cluster '',2x,''E /Eh '',7x,''De/kcal'',3x,&
            &''Detot/kcal'',2x,''Opt'',4x)')
  end subroutine pr_fill_energy
  subroutine pr_ensemble_energy()
    implicit none
    write (stdout,*)
    write (stdout,'(x,'' Cluster'',3x,''E /Eh '',7x,&
             &''Density'',2x,''Efix'',7x,''R   av/act.'',1x,&
             &''Surface'',3x,''Opt'',4x)')
  end subroutine pr_ensemble_energy
  subroutine pr_qcg_esolv()
    implicit none
    write (stdout,*)
    write (stdout,'(2x,''========================================='')')
    write (stdout,'(2x,''|   quantum cluster growth: ESOLV       |'')')
    write (stdout,'(2x,''|                                       |'')')
  end subroutine pr_qcg_esolv
  subroutine pr_grow_energy()
    implicit none
    write (stdout,'(x,'' Size'',7x,''E'',8x,''De'',7x,''Detot'',6x,&
             &''Density'',5x,''Eatom'',4x,''av. R'', 1x,'' Rlast'',3x,&
             &''Volume'',4x,''Opt'')')
    write (stdout,'(12x,''[Eh]'',4x,''[kcal]'',5x,''[kcal]'',5x,&
             &''[u/Å^3]'',5x,''[kcal]'',3x,''[bohr]'', 1x,''[bohr]'',1x,&
             &''[bohr^3]'')')

  end subroutine pr_grow_energy

  subroutine pr_freq_file(ich)
    implicit none
    integer :: ich
    write (ich,'(2x,"#       H(T)       SVIB      SROT       STRA      G(T)")')
    write (ich,'(2x,"     [kcal/mol]    [      cal/mol/K        ]    [kcal/mol]")')
    write (ich,'(2x,"--------------------------------------------------------")')
  end subroutine pr_freq_file

!========================================================================================!

  subroutine xtbiff_print_deprecated()
    external creststop
    write (stdout,*)
    write (stdout,*) 'WARNING WARNING WARNING'
    write (stdout,*) ' The use of xtbiff in QCG is deprecated and is disabled'
    write (stdout,*) ' following CREST 3.1.0, in favor of the aISS algorithm.'
    write (stdout,*) ' This requires a current version of the xtb program.'
    write (stdout,*)
    call creststop(status_safety)
  end subroutine xtbiff_print_deprecated

!========================================================================================!
!========================================================================================!
end module qcg_printouts
