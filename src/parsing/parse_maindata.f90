!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022-2023 Philipp Pracht
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

!> NOTE: This is work in progress, not all input conventions have been set yet
!========================================================================================!
!> Routines contained here are for parsing 'top level' settings that will
!> enter the env systemdata object

module parse_maindata
  use crest_parameters
  !> modules for data storage in crest
  use crest_data
  use strucrd,only:coord
  use molecule_parameters,only:extxyz_units_global
  !> modules used for parsing the root_object
  use parse_keyvalue,only:keyvalue,valuetypes
  use parse_block,only:datablock
  use parse_datastruct,only:root_object
  !> Declarations
  implicit none
  public

  character(len=*),parameter,private :: fmturk = '("unrecognized KEYWORD in ",a," : ",a)'
  character(len=*),parameter,private :: fmtura = '("unrecognized ARGUMENT : ",a)'

  external creststop

!========================================================================================!
!========================================================================================!
contains   !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine parse_main_auto(env,kv,istat)
    implicit none
    type(systemdata) :: env
    type(keyvalue) :: kv
    integer,intent(inout) :: istat
    integer :: istat_ref
    logical :: rd
    istat_ref = istat
    rd = .false.
    select case (kv%id)
    case (valuetypes%float)  !>--- float
      call parse_main_float(env,kv%key,kv%value_f,rd)
    case (valuetypes%int)    !>--- int
      call parse_main_int(env,kv%key,kv%value_i,rd)
    case (valuetypes%bool)   !>--- bool
      call parse_main_bool(env,kv%key,kv%value_b,rd)
    case (valuetypes%string) !>--- string
      call parse_main_c(env,kv%key,kv%value_c,rd)
    end select
!>--- other, with multiple or raw type
    if (.not.rd) then
      select case (kv%key)
      case ('optlev','ancopt_level')
        env%optlev = optlevnum(kv%rawvalue)

      case ('split')
        if (kv%id .ne. valuetypes%int_array.or. &
          & kv%na < 3) then
          write (stdout,'(a)') '**ERROR** "split" must be a list of at least 3 integers'
          call creststop(status_config)
        end if
        call env%addsplitqueue(kv%value_ia)

      case default
        istat = istat+1
      end select
    end if
!>--- if none of the options was recognizeda and istat increased as a consequence, print that
    if (istat > istat_ref) then
      write (stdout,fmturk) 'main section',kv%key
    end if

  end subroutine parse_main_auto
  subroutine parse_main_float(env,key,val,rd)
    implicit none
    type(systemdata) :: env
    character(len=*) :: key
    real(wp) :: val
    logical,intent(out) :: rd
    rd = .true.
    select case (key)
    case ('wscal')
      env%potscal = val
      env%wallsetup = .true.
    case ('wpad')
      env%potpad = val
      env%wallsetup = .true.
    case default
      rd = .false.
    end select
    return
  end subroutine parse_main_float
  subroutine parse_main_int(env,key,val,rd)
    implicit none
    type(systemdata) :: env
    character(len=*) :: key
    logical,intent(out) :: rd
    integer :: val
    rd = .true.
    select case (key)
    case ('threads','parallel')
      env%Threads = val
      env%autothreads = .true.
      env%threadssetmanual = .true.
    case default
      rd = .false.
    end select
    return
  end subroutine parse_main_int
  subroutine parse_main_c(env,key,val,rd)
    implicit none
    type(systemdata) :: env
    character(len=*) :: key
    character(len=*) :: val
    logical,intent(out) :: rd
    type(coord) :: mol
    rd = .true.
    select case (key)
    case ('bin','binary')
      env%ProgName = val
    case ('runtype')
      select case (val)
      case ('none')
        env%crestver = crest_none
      case ('playground','test')
        env%preopt = .false.
        env%crestver = crest_test
      case ('singlepoint','sp')
        env%preopt = .false.
        env%crestver = crest_sp
      case ('numgrad')
        env%preopt = .false.
        env%crestver = crest_sp
        env%testnumgrad = .true.
      case ('ancopt','optimize','ohess')
        env%preopt = .false.
        env%crestver = crest_optimize
        env%optlev = 0.0_wp
        if (val .eq. 'ohess') then
          env%crest_ohess = .true.
        end if
      case ('ancopt_ensemble','optimize_ensemble','mdopt')
        env%preopt = .false.
        env%crestver = crest_mdopt2
        env%optlev = 0.0d0
      case ('screen_ensemble','screen')
        env%preopt = .false.
        env%crestver = crest_screen
      case ('ensemble_singlepoints','ensemblesp','mdsp')
        env%preopt = .false.
        env%crestver = crest_ensemblesp

      case ('md','mtd','metadynamics','dynamics')
        env%preopt = .false.
        env%crestver = crest_moldyn
      case ('scan')
        env%preopt = .false.
        env%crestver = crest_scanning
      case ('search_1')
        env%preopt = .true.
        env%crestver = crest_s1
        env%runver = crest_s1
      case ('mecp','mecp_search')
        env%preopt = .false.
        env%crestver = crest_mecp
        env%runver = crest_mecp
      case ('imtd-gc','mtd_search')
        env%preopt = .true.
        env%crestver = crest_imtd
        env%runver = 1
      case ('mtd_search_quick')
        env%preopt = .true.
        env%crestver = crest_imtd
        env%quick = .true.
        env%runver = 2
        env%ewin = 5.0d0
        env%optlev = 1.0d0    !> optlev tight for quick run
      case ('mtd_search_mquick')
        env%preopt = .true.
        env%crestver = crest_imtd
        env%rotamermds = .false.      !> no NORMMD
        env%performCross = .false.    !> no GC
        env%quick = .true.            !> MTD settings from the quick-mode
        env%superquick = .true.       !> use user-set opt level in Multilevel opt.
        env%Maxrestart = 1            !> only one MTD iteration
        env%runver = 6
        env%optlev = 0.0d0  !> user-set opt level
        env%ewin = 2.5d0              !> smaller energy window
      case ('mtd_search_squick')
        env%preopt = .true.
        env%crestver = crest_imtd
        env%rotamermds = .false.      !> no NORMMD
        env%performCross = .false.    !> no GC
        env%quick = .true.            !> MTD settings from the quick-mode
        env%superquick = .true.       !> use user-set opt level in Multilevel opt.
        env%runver = 5
        env%optlev = 0.0d0            !> user-set opt level
        env%ewin = 5.0d0              !> smaller energy window
      case ('nci-mtd','nci','nci_search')
        env%NCI = .true.
        env%runver = 4
        env%autozsort = .false.
        env%performCross = .false.
        env%rotamermds = .false.
        env%crestver = crest_imtd
      case ('bh','gmin')
        env%crestver = crest_bh
      case ('entropy','imtd-smtd','entropy_search')
        env%crestver = crest_imtd  !> the entropy mode acts as subtype of the crest_imtd algo
        env%properties = abs(p_CREentropy)
        env%autozsort = .false.     !> turn off zsort (since we are not going to GC anyways)
        env%performCross = .false.  !> turn off GC
        env%entropic = .true.       !> indicator for this runtype
        env%Maxrestart = 1          !> turn off MTD iterations (just do one)
        env%rotamermds = .false.    !> turn off normMDs
        env%entropymd = .true.      !> special static MTDs
        env%runver = 111            !> version  for selection of MTD bias settings
        env%doNMR = .true.          !> we need equivalencies
        env%emtd%bhess = .false.    !> currently there is no BHESS version, TODO!
        call env%addjob(env%properties)
      case ('numhess','numerical hessian')
        env%preopt = .false.
        env%crestver = crest_numhessian
        env%runver = crest_numhessian
      case ('rigidconf')
        env%preopt = .true.
        env%crestver = crest_rigcon
        env%runver = crest_rigcon

      case ('ttconf')
        env%preopt = .true.
        env%crestver = crest_ttc
        env%runver = crest_ttc

      case ('protonate')
        env%properties = p_protonate
        env%crestver = crest_protonate

      case ('deprotonate')
        env%properties = p_deprotonate
        env%crestver = crest_deprotonate

      case ('tautomerize')
        env%properties = p_tautomerize
        env%crestver = crest_tautomerize

      case ('thermo')
        env%properties = p_thermo
        env%crestver = crest_none
        env%preopt = .false.

      case ('cregen','sort')
        env%preopt    = .false.
        env%crestver  = crest_sorting
        env%autozsort = .false.
        if (val .eq. 'cregen') then
          env%sortmode = 'cregen'
          env%confgo   = .true.
        end if

      case default
        !>--- keyword was recognized, but invalid argument supplied
        write (stdout,fmtura) val
        call creststop(status_config)

      end select
    case ('ensemble_input','ensemble','input_ensemble')
      env%ensemblename = val
      env%inputcoords = val
    case ('input','structure','coord','coords')
      env%inputcoords = val
      call mol%open(val)
      call env%ref%load(mol)

    case ('constraints','xtbconstraints','xtbinput') !> equivalent to --cinp
      env%constraints = val
    case ('rigidconf_file')
      env%rigidconf_userfile = val

    case ('watlist','wat')
      env%potatlist = val
      env%wallsetup = .true.
    case ('extxyz_units')
      extxyz_units_global = trim(val)
    case default
      rd = .false.
    end select
    return
  end subroutine parse_main_c
  subroutine parse_main_bool(env,key,val,rd)
    implicit none
    type(systemdata) :: env
    character(len=*) :: key
    logical,intent(out) :: rd
    logical :: val
    rd = .true.
    select case (key)
    case ('preopt')
      env%preopt = val
    case ('noopt')
      env%preopt = .not.val
    case ('topo')
      env%checktopo = val
    case ('notopo')
      env%checktopo = .not.val
    case ('multilevelopt')
      env%multilevelopt = val
    case ('refine_presort')
      env%refine_presort = val

    case ('omp_nested')
      env%omp_allow_nested = val
    case default
      rd = .false.
    end select
    return
  end subroutine parse_main_bool
!========================================================================================!

  subroutine parse_main_blk(env,blk,istat)
!**************************************
!* Some shorter blocks are not defined
!* in separate source files. They can
!* be found below.
!**************************************
    implicit none
    type(systemdata) :: env
    type(datablock) :: blk
    integer,intent(inout) :: istat
    select case (blk%header)
    case ('cregen')
      call parse_cregen(env,blk,istat)
    case ('thermo')
      call parse_thermo(env,blk,istat)
    case ('protonation')
      call parse_protonation(env,blk,istat)
    case ('ttconf')
      call parse_ttconf(env,blk,istat)
    end select
  end subroutine parse_main_blk

!========================================================================================!
  subroutine parse_ttconf(env,blk,istat)
!*******************************************************
!* parse detailed settings for the TTConf-light runtype
!* ([ttconf] block). A "preset" is applied first, then
!* individual keys may override it.
!*******************************************************
    implicit none
    type(systemdata) :: env
    type(datablock) :: blk
    type(keyvalue) :: kv
    integer,intent(inout) :: istat
    integer :: i
    logical :: ok
!>--- apply a preset first (if present), so explicit keys override it
    do i = 1,blk%nkv
      if (blk%kv_list(i)%key == 'preset') then
        call env%ttconf%setpreset(blk%kv_list(i)%value_c,ok)
        if (.not.ok) write (stdout,'(1x,a)') &
        &  '**WARNING** unknown [ttconf] preset: '//trim(blk%kv_list(i)%value_c)
      end if
    end do
!>--- parse the remaining keys
    do i = 1,blk%nkv
      kv = blk%kv_list(i)
      select case (kv%key)
      case ('preset')                  !> already handled above
        continue
      case ('rank','r')
        env%ttconf%rank = kv%value_i
      case ('sweeps','s')
        env%ttconf%sweeps = kv%value_i
      case ('grid','ngrid')
        env%ttconf%ngrid = kv%value_i
      case ('ninit')
        env%ttconf%ninit = kv%value_i
      case ('ewin')
        env%ttconf%ewin = kv%value_f
      case ('kt','temperature')
        env%ttconf%kt = kv%value_f
      case ('bruteforce','oracle')
        env%ttconf%use_sweep = .not.kv%value_b
      case ('sp','singlepoint','sponly') !> singlepoints only (no geometry opt)
        env%ttconf%sp_only = kv%value_b
      case ('cache','ecache')
        env%ttconf%use_cache = kv%value_b
      case ('ringbonds')              !> treat in-ring bonds as TT variables?
        env%ttconf%excl_rings = .not.kv%value_b
      case ('ringsample')             !> sample ring templates as TT sites?
        env%ttconf%ring_sample = kv%value_b
      case ('ringmethod','ringsampler') !> which ring-conformation generator
        env%ttconf%ring_method = trim(kv%value_c)
      case ('seed')
        env%ttconf%seed = kv%value_i
      case ('bonds','userbonds')      !> force atom pairs to be TT variables
        call parse_ttconf_bonds(env,kv)
      case default
        istat = istat+1
        write (stdout,fmturk) '[ttconf]-block',kv%key
      end select
    end do
  end subroutine parse_ttconf

!========================================================================================!
  subroutine parse_ttconf_bonds(env,kv)
!*******************************************************
!* Parse the [ttconf] "bonds" key: atom pairs (optionally
!* with a per-bond grid-point count) that the user wants
!* treated as TT variables, e.g.
!*    bonds = [[1,2], [3,4,12]]
!* Each entry is [A,B] or [A,B,npoints]. They are stored as
!* env%ttconf%userbonds(3,:) = (A,B,npoints), npoints = 0
!* meaning "use the default grid". A single pair may also be
!* given flat, e.g. bonds = [1,2].
!*
!* NOTE: the toml reader expands the nested form [[..],[..]]
!* into one int-array key per sub-array, so this routine may
!* be called repeatedly for the same key -- bonds therefore
!* ACCUMULATE (each call appends) rather than overwrite.
!*******************************************************
    implicit none
    type(systemdata),intent(inout) :: env
    type(keyvalue),intent(in) :: kv
    integer :: k,a,b,np
    character(len=:),allocatable :: s

    if (kv%id == valuetypes%int_array) then
!>--- a flat pair: bonds = [1,2] or [1,2,12] (also each [[..]] sub-array)
      if (kv%na < 2) then
        write (stdout,'(1x,a)') '**WARNING** [ttconf] bonds entry needs >=2 atoms; ignored'
        return
      end if
      np = 0
      if (kv%na >= 3) np = abs(kv%value_ia(3))
      call append_userbond(env,kv%value_ia(1),kv%value_ia(2),np)
    else if (kv%id == valuetypes%raw_array) then
!>--- a raw list of sub-arrays: bonds = [[1,2],[3,4,12]]
      do k = 1,kv%na
        s = trim(adjustl(kv%value_rawa(k)))
        call read_bracketed_triplet(s,a,b,np)
        if (a <= 0.or.b <= 0) then
          write (stdout,'(1x,a,a)') '**WARNING** could not parse [ttconf] bond entry: ', &
          &  trim(kv%value_rawa(k))
          cycle
        end if
        call append_userbond(env,a,b,np)
      end do
    else
      write (stdout,'(1x,a)') '**WARNING** [ttconf] bonds must be a list of atom pairs; ignored'
    end if
  end subroutine parse_ttconf_bonds

!========================================================================================!
  subroutine append_userbond(env,a,b,np)
!*******************************************************
!* Append one (atomA, atomB, npoints) bond to the growing
!* env%ttconf%userbonds(3,:) list.
!*******************************************************
    implicit none
    type(systemdata),intent(inout) :: env
    integer,intent(in) :: a,b,np
    integer,allocatable :: tmp(:,:)
    integer :: n
    if (.not.allocated(env%ttconf%userbonds)) then
      allocate (env%ttconf%userbonds(3,1))
      env%ttconf%userbonds(:,1) = [a,b,np]
      return
    end if
    n = size(env%ttconf%userbonds,2)
    allocate (tmp(3,n+1))
    tmp(:,1:n) = env%ttconf%userbonds(:,1:n)
    tmp(:,n+1) = [a,b,np]
    call move_alloc(tmp,env%ttconf%userbonds)
  end subroutine append_userbond

!========================================================================================!
  subroutine read_bracketed_triplet(str,a,b,np)
!*******************************************************
!* Read up to three integers from a "[A,B]" / "[A,B,N]"
!* bracketed substring. Missing values come back as 0.
!*******************************************************
    implicit none
    character(len=*),intent(in) :: str
    integer,intent(out) :: a,b,np
    character(len=:),allocatable :: clean
    integer :: i,io
    a = 0; b = 0; np = 0
    clean = ''
    do i = 1,len_trim(str)
      select case (str(i:i))
      case ('[',']')
        cycle                          !> drop brackets
      case (',')
        clean = clean//' '             !> commas -> blanks for list-directed read
      case default
        clean = clean//str(i:i)
      end select
    end do
    read (clean,*,iostat=io) a,b,np    !> try three
    if (io /= 0) then
      np = 0
      read (clean,*,iostat=io) a,b     !> fall back to two
      if (io /= 0) then
        a = 0; b = 0
      end if
    end if
    np = abs(np)
  end subroutine read_bracketed_triplet

!========================================================================================!
  subroutine parse_cregen(env,blk,istat)
!****************************************
!* parse settings for the CREGEN routine
!****************************************
    implicit none
    type(systemdata) :: env
    type(datablock) :: blk
    type(keyvalue) :: kv
    integer,intent(inout) :: istat
    integer :: i
!>--- parse the arguments
    do i = 1,blk%nkv
      kv = blk%kv_list(i)
      select case (kv%key)
      case ('ewin')
        env%ewin = kv%value_f
      case ('ethr')
        env%ethr = kv%value_f
      case ('rthr')
        env%rthr = kv%value_f
      case ('bthr')
        env%bthr2 = kv%value_f
      case ('eqv','nmr')
        env%doNMR = kv%value_b
      case default
        !>--- unrecognized keyword
        istat = istat+1
        write (stdout,fmturk) '[cregen]-block',kv%key
      end select
    end do
  end subroutine parse_cregen

!========================================================================================!

  subroutine parse_thermo(env,blk,istat)
!****************************************
!* parse settings for the Thermo routine
!****************************************
    implicit none
    type(systemdata) :: env
    type(datablock) :: blk
    type(keyvalue) :: kv
    integer,intent(inout) :: istat
    integer :: i
!>--- parse the arguments
    do i = 1,blk%nkv
      kv = blk%kv_list(i)
      select case (kv%key)
      case ('ithr','freq_ithr','freq_invert')
        env%thermo%ithr = kv%value_f
      case ('fscal','freq_scal')
        env%thermo%fscal = kv%value_f
      case ('sthr','freq_interpol')
        env%thermo%sthr = kv%value_f
      case ('trange')
        if (kv%na >= 2) then
          env%thermo%trange(1) = minval(kv%value_fa(1:2),1)
          env%thermo%trange(2) = maxval(kv%value_fa(1:2),1)
        end if
        if (kv%na >= 3) then
          env%thermo%trange(3) = kv%value_fa(3)
        end if
      case ('tstep')
        env%thermo%trange(3) = kv%value_f

      case ('input','coords')
        env%thermo%coords = kv%value_c
        if (allocated(env%thermo%vibfile)) env%properties = p_thermo
      case ('freq_input','vibs','hessian')
        env%thermo%vibfile = kv%value_c
        if (allocated(env%thermo%coords)) env%properties = p_thermo

      case ('entropy_model','svib_model')
        select case (kv%value_c)
        case ('grimme')
          env%thermo%emodel = kv%value_c
        case ('truhlar') 
          env%thermo%emodel = kv%value_c 
          env%thermo%sthr = 100.0_wp
        case default
          write (stdout,fmtura) trim(kv%rawvalue)
          call creststop(status_config)
        end select

      case default
        !>--- unrecognized keyword
        istat = istat+1
        write (stdout,fmturk) '[thermo]-block',kv%key
      end select
    end do
  end subroutine parse_thermo

!========================================================================================!
  subroutine parse_protonation(env,blk,istat)
!******************************************
!* parse settings for protonation settings
!******************************************
    implicit none
    type(systemdata) :: env
    type(datablock) :: blk
    type(keyvalue) :: kv
    integer,intent(inout) :: istat
    integer :: i
    external :: swparse
!>--- parse the arguments
    do i = 1,blk%nkv
      kv = blk%kv_list(i)
      select case (kv%key)
      case ('ewin')
        env%protb%ewin = kv%value_f
      case ('swel','ion')
        call swparse(kv%value_c,env%protb)
      case ('ffopt')
        env%protb%ffopt = kv%value_b
      case ('freezeopt')
        env%protb%hnewopt = kv%value_b
      case ('finalopt')
        env%protb%finalopt = kv%value_b

      case ('activelmo')
        env%protb%active_lmo(1:) = kv%value_ba(1:)
      case ('pi')
        env%protb%active_lmo(1) = kv%value_b
      case ('lp')
        env%protb%active_lmo(2) = kv%value_b
      case ('delpi','delocpi')
        env%protb%active_lmo(3) = kv%value_b
      case default
        !>--- unrecognized keyword
        istat = istat+1
        write (stdout,fmturk) '[protonation]-block',kv%key
      end select
    end do
  end subroutine parse_protonation

!========================================================================================!
!========================================================================================!
end module parse_maindata
