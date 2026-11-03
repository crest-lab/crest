!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2024 Philipp Pracht
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

!> Routines for parsing and setting up hybrid two-level method combinations
!> expressed as CLI arguments of the form A@B, A//B, A/sp/B, A/opt/B.

module parse_hybrid
  use crest_parameters
  use crest_data
  use crest_calculator,only:calcdata,calculation_settings,jobtype
  use iomod,only:lowercase
  implicit none
  private

  public :: parse_hybrid_argument
  public :: is_valid_method
  public :: setup_hybrid_calc
  public :: level_add_cli_solvation

  external creststop

contains

!========================================================================================!

  subroutine parse_hybrid_argument(str,quality,workhorse,mode,iostat)
  !***************************************************************
  !* Parse a composite method argument of the form:
  !*   A@B     → quality=A, workhorse=B, mode='at'  (post-search opt)
  !*   A//B    → quality=A, workhorse=B, mode='sp'  (inline SP refine)
  !*   A/sp/B  → quality=A, workhorse=B, mode='sp'
  !*   A/opt/B → quality=A, workhorse=B, mode='opt' (inline geo-opt refine)
  !* Returns iostat=0 on success, non-zero if not recognised.
  !*
  !* Input:
  !*   str  - argument string (without leading dash)
  !* Output:
  !*   quality   - the higher-level method token (left of separator)
  !*   workhorse - the fast workhorse method token (right of separator)
  !*   mode      - 'at', 'sp', or 'opt'
  !*   iostat    - 0 on success
  !***************************************************************
    implicit none
    character(len=*),intent(in) :: str
    character(len=:),intent(out),allocatable :: quality,workhorse
    character(len=4),intent(out) :: mode
    integer,intent(out) :: iostat
    integer :: k
    character(len=:),allocatable :: s

    iostat = 1
    quality = ''
    workhorse = ''
    mode = ''
    s = lowercase(trim(str))

  ! ── try each separator in order of specificity ────────────────
    k = index(s,'/opt/')
    if (k > 0) then
      quality   = s(1:k-1)
      workhorse = s(k+5:)
      mode = 'opt'
    end if
    if (mode == '') then
      k = index(s,'/sp/')
      if (k > 0) then
        quality   = s(1:k-1)
        workhorse = s(k+4:)
        mode = 'sp'
      end if
    end if
    if (mode == '') then
      k = index(s,'//')
      if (k > 0) then
        quality   = s(1:k-1)
        workhorse = s(k+2:)
        mode = 'sp'
      end if
    end if
    if (mode == '') then
      k = index(s,'@')
      if (k > 0) then
        quality   = s(1:k-1)
        workhorse = s(k+1:)
        mode = 'at'
      end if
    end if

    if (mode == '') return

    if (.not.is_valid_method(quality)) return
    if (.not.is_valid_method(workhorse)) return
    iostat = 0
  end subroutine parse_hybrid_argument

!========================================================================================!

  logical function is_valid_method(token)
  !*************************************************
  !* Returns .true. if token names a supported method.
  !*************************************************
    implicit none
    character(len=*),intent(in) :: token
    select case (lowercase(trim(token)))
    case ('gfn0','gfn1','gfn2','gxtb','gfnff','gff')
      is_valid_method = .true.
    case default
      is_valid_method = .false.
    end select
  end function is_valid_method

!========================================================================================!

  subroutine setup_hybrid_calc(env,calc,workhorse_str,quality_str,mode)
  !*******************************************************************
  !* Set up a two-level calcdata from a hybrid method pair.
  !* Called after all CLI flags are parsed, so charge, spin and the
  !* implicit solvation from env apply to both levels.
  !*  workhorse → primary calc (refine_lvl=0, runs during iMTD-GC)
  !*  quality   → secondary calc:
  !*               mode='sp'  → refine_lvl=singlepoint (inline SP)
  !*               mode='opt' → refine_lvl=geoopt      (inline opt)
  !*               mode='at'  → refine_lvl=post_opt    (post-search + pqueue 51)
  !*
  !* Input:
  !*   workhorse_str - method string for the fast workhorse (e.g. 'gfnff')
  !*   quality_str   - method string for the quality level  (e.g. 'gfn2')
  !*   mode          - 'sp', 'opt', or 'at'
  !* Output:
  !*   calc          - the two levels are appended to it
  !*******************************************************************
    implicit none
    type(systemdata),intent(inout) :: env
    type(calcdata),intent(inout) :: calc
    character(len=*),intent(in) :: workhorse_str,quality_str,mode
    type(calculation_settings) :: cal_work,cal_qual
    integer :: rlvl

    select case (trim(mode))
    case ('sp');  rlvl = refine%singlepoint
    case ('opt'); rlvl = refine%geoopt
    case ('at');  rlvl = refine%post_opt
    case default; rlvl = refine%singlepoint
    end select

  ! ── workhorse (active during normal search) ────────────────────
    call cal_work%create(workhorse_str)
    cal_work%chrg = env%chrg
    cal_work%uhf  = env%uhf
    call cal_work%sync_multiplicity()
    cal_work%spin_polarized = env%spin_polarized
    call hybrid_solvation(cal_work,workhorse_str)
    cal_work%refine_lvl = 0
    call cal_work%autocomplete(1)
    call calc%add(cal_work)

  ! ── quality level ──────────────────────────────────────────────
    call cal_qual%create(quality_str)
    cal_qual%chrg = env%chrg
    cal_qual%uhf  = env%uhf
    call cal_qual%sync_multiplicity()
    cal_qual%spin_polarized = env%spin_polarized
    call hybrid_solvation(cal_qual,quality_str)
    cal_qual%refine_lvl = rlvl
    call cal_qual%autocomplete(2)
    call calc%add(cal_qual)

  ! ── register refine_queue / pqueue ────────────────────────────
    if (trim(mode) == 'at') then
      call env%addjob(51)
      call env%checkhy()
      write(stdout,'(2x,a,"@",a,a)') trim(quality_str),trim(workhorse_str), &
        & ' : post-search re-optimization of conformer ensemble'
    else
  !   Mirror what env2calc does for --refine: populate refine_queue now
  !   so the inline refine path is active even without a TOML input file.
      call env%addrefine(rlvl)
      if (trim(mode) == 'opt') then
        write(stdout,'(2x,a,"/opt/",a,a)') trim(quality_str),trim(workhorse_str), &
          & ' : inline geometry refinement'
      else
        write(stdout,'(2x,a,"//",a,a)') trim(quality_str),trim(workhorse_str), &
          & ' : inline singlepoint re-ranking'
      end if
    end if

  contains

    subroutine hybrid_solvation(cal,method)
      type(calculation_settings),intent(inout) :: cal
      character(len=*),intent(in) :: method
      logical :: ok
      call level_add_cli_solvation(env,cal,ok)
      if (.not.ok) then
        write (stdout,'(a,a,a)') '**ERROR** implicit solvation (',trim(env%solv), &
          & ') is not available for '//trim(method)//' in the composite method'
        call creststop(status_config)
      end if
    end subroutine hybrid_solvation
  end subroutine setup_hybrid_calc

!========================================================================================!

  subroutine level_add_cli_solvation(env,cal,ok)
  !*******************************************************
  !* Copy the CLI implicit solvation (--alpb/--gbsa) into
  !* one calculation level. ok is .false. if the level's
  !* method has no implicit solvation in the calculator
  !* layer (e.g. g-xTB via the xtb binary, ORCA, MLIPs),
  !* so the caller can stop instead of silently running
  !* in the gas phase. Without --alpb/--gbsa: no-op.
  !*******************************************************
    implicit none
    type(systemdata),intent(in) :: env
    type(calculation_settings),intent(inout) :: cal
    logical,intent(out) :: ok
    ok = .true.
    if (.not.env%gbsa) return
    select case (cal%id)
    case (jobtype%tblite,jobtype%gfn0,jobtype%gfn0occ,jobtype%gfnff)
      continue
    case default
      ok = .false.
      return
    end select
    if (index(env%solv,'gbsa') .ne. 0) then
      cal%solvmodel = 'gbsa'
    else
      cal%solvmodel = 'alpb'
    end if
    cal%solvent = trim(env%solvent)
  end subroutine level_add_cli_solvation

!========================================================================================!
end module parse_hybrid
