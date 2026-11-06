!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022-2026 Christoph Plett, Philipp Pracht
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
!
! Ported from the xtb docking submodule (github.com/grimme-lab/xtb,
! src/prog/dock.f90, src/docking/set_module.f90), LGPL-3.0-or-later.
!================================================================================!

!> Standalone aISS docking runtype (replacement for "xtb dock A B")
module aiss_standalone
  use crest_parameters,only:wp,stdout
  use crest_data
  use strucrd,only:coord,get_atlist
  use crest_calculator,only:calcdata
  use aiss_param
  use aiss_dock
  use parse_keyvalue,only:keyvalue
  use parse_block,only:datablock
  implicit none
  private

  public :: parse_dock_cli
  public :: parse_dock_toml
  public :: crest_aiss_dock

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine parse_dock_cli(env,nra,arg,processedarg)
    !*****************************************************************
    !* If -dock/--dock is present, take over all docking options
    !* (xtb dock syntax) before the general parser sees them:
    !*   crest A --dock B [opts]   or   crest --dock A B [opts]
    !* Presets (--fast, --qcg) are applied before explicit values.
    !*****************************************************************
    type(systemdata),intent(inout) :: env
    integer,intent(in) :: nra
    character(len=*),intent(in) :: arg(nra)
    logical,intent(inout) :: processedarg(nra)
    type(aiss_input) :: dk
    character(len=:),allocatable :: flag,val
    integer :: i,j,idock,npos,io
    logical :: ex,preset_done
    real(wp) :: rdum
    integer :: idum

    idock = 0
    do i = 1,nra
      if (trim(arg(i)) == '-dock'.or.trim(arg(i)) == '--dock') idock = i
    end do
    if (idock == 0) return
    dk = env%dock
    env%crestver = crest_dock
    processedarg(idock) = .true.
    write (stdout,'(2x,a,1x,a)') trim(arg(idock)),': standalone aISS docking'

!>--- positional structure files after --dock
    npos = 0
    do j = idock+1,min(idock+2,nra)
      if (arg(j)(1:1) == '-') exit
      inquire (file=trim(arg(j)),exist=ex)
      if (.not.ex) exit
      npos = npos+1
      processedarg(j) = .true.
      if (npos == 1) dk%fileB = trim(arg(j))
      if (npos == 2) then
        dk%fileA = dk%fileB
        dk%fileB = trim(arg(j))
      end if
    end do
    !> crest --dock A B: A also serves as the general input structure
    if (allocated(dk%fileA)) env%inputcoords = dk%fileA

!>--- presets first (xtb: --fast and --qcg share one switch, first wins)
    preset_done = .false.
    do i = 1,nra
      if (processedarg(i)) cycle
      flag = nodash(arg(i))
      if (flag /= 'fast'.and.flag /= 'qcg') cycle
      processedarg(i) = .true.
      if (preset_done) cycle
      if (flag == 'fast') call dk%cfg%set_fast()
      if (flag == 'qcg') call dk%cfg%set_qcg()
      preset_done = .true.
    end do

    do i = 1,nra
      if (processedarg(i)) cycle
      if (arg(i)(1:1) /= '-') cycle
      flag = nodash(arg(i))
      val = ''
      if (i < nra) val = trim(arg(i+1))
      select case (flag)
      !>-- switches
      case ('pocket','nopocket','onlypocket','stack','nostack','onlystack', &
      &     'angular','noangular','onlyangular')
        call set_grids(dk%cfg,flag)
      case ('cs')
        dk%cfg%cssym = .true.
      case ('atm')
        dk%cfg%fulle = .true.
      case ('ensemble','ens')
        dk%cfg%docking_ens = .true.
      case ('samerand')
        dk%cfg%samerand = .true.
      case ('verbose')
        dk%cfg%debug = .true.
      case ('attractive')
        dk%directed_type = p_atom_att
      case ('repulsive')
        dk%directed_type = p_atom_pot
      !>-- options with a value
      case ('stepr','stepa','etemp','nfinal','maxgen','maxparent','nstack','nfrag1', &
      &     'chrg','chrg1','chrg2','uhf','uhf1','uhf2','i','input','optlvl','gfn', &
      &     'aiss_lmo','lmo','o','opt','fix','directed')
        if (i == nra) then
          write (stdout,'(a,a)') '**ERROR** missing value for ',trim(arg(i))
          call creststop(status_input)
        end if
        processedarg(i+1) = .true.
        select case (flag)
        case ('stepr')
          read (val,*,iostat=io) rdum
          if (io == 0) dk%cfg%stepr = max(1.0_wp,rdum)
        case ('stepa')
          read (val,*,iostat=io) rdum
          if (io == 0) dk%cfg%stepa = max(1.0_wp,rdum)
        case ('etemp')
          read (val,*,iostat=io) rdum
          if (io == 0) dk%cfg%etemp = rdum
        case ('nfinal')
          read (val,*,iostat=io) idum
          if (io == 0) dk%cfg%n_opt = idum
          dk%nfinal_set = .true.
        case ('maxgen')
          read (val,*,iostat=io) idum
          if (io == 0) dk%cfg%maxgen = idum
        case ('maxparent')
          read (val,*,iostat=io) idum
          if (io == 0) dk%cfg%maxparent = idum
        case ('nstack')
          read (val,*,iostat=io) idum
          if (io == 0) dk%cfg%mxcma = idum
        case ('nfrag1')
          read (val,*,iostat=io) idum
          if (io == 0) dk%cfg%nfrag1 = idum
        case ('chrg','chrg1','chrg2','uhf','uhf1','uhf2')
          read (val,*,iostat=io) idum
          if (io == 0) call set_chrguhf(dk,flag,idum)
        case ('i','input')
          dk%xcontrol = val
        case ('optlvl','gfn')
          dk%cfg%optlvl = level_name(val)
        case ('aiss_lmo','lmo')
          dk%cfg%lmolvl = level_name(val)
        case ('o','opt')
          dk%cfg%optlev = nint(optlevnum(val,io))
        case ('fix')
          dk%fix = val
        case ('directed')
          dk%directed = val
        end select
      case default
        cycle
      end select
      processedarg(i) = .true.
    end do
    env%dock = dk
  end subroutine parse_dock_cli

!========================================================================================!

  subroutine parse_dock_toml(env,blk,istat)
    !*****************************************************************
    !* [dock] block of the TOML input; a "preset" (fast/qcg) is
    !* applied first so that explicit keys override it
    !*****************************************************************
    type(systemdata),intent(inout) :: env
    type(datablock),intent(in) :: blk
    integer,intent(inout) :: istat
    type(keyvalue) :: kv
    integer :: i

    do i = 1,blk%nkv
      if (blk%kv_list(i)%key /= 'preset') cycle
      select case (blk%kv_list(i)%value_c)
      case ('fast')
        call env%dock%cfg%set_fast()
      case ('qcg')
        call env%dock%cfg%set_qcg()
      case default
        write (stdout,'(1x,a)') '**WARNING** unknown [dock] preset: '//trim(blk%kv_list(i)%value_c)
      end select
    end do

    do i = 1,blk%nkv
      kv = blk%kv_list(i)
      select case (kv%key)
      case ('preset')
      case ('molecule_a','file_a','host')
        env%dock%fileA = kv%value_c
      case ('molecule_b','file_b','guest','partner')
        env%dock%fileB = kv%value_c
      case ('xcontrol','input')
        env%dock%xcontrol = kv%value_c
      case ('level','method','optlvl')
        env%dock%cfg%optlvl = level_name(kv%value_c)
      case ('lmo','lmolvl','lmo_level')
        env%dock%cfg%lmolvl = level_name(kv%value_c)
      case ('optlevel')
        if (kv%id == 2) then
          env%dock%cfg%optlev = kv%value_i
        else
          env%dock%cfg%optlev = nint(optlevnum(kv%value_c))
        end if
      case ('maxgen')
        env%dock%cfg%maxgen = kv%value_i
      case ('maxparent')
        env%dock%cfg%maxparent = kv%value_i
      case ('nfinal','n_opt')
        env%dock%cfg%n_opt = kv%value_i
        env%dock%nfinal_set = .true.
      case ('nstack')
        env%dock%cfg%mxcma = kv%value_i
      case ('nfrag1')
        env%dock%cfg%nfrag1 = kv%value_i
      case ('stepr')
        env%dock%cfg%stepr = max(1.0_wp,toreal(kv))
      case ('stepa')
        env%dock%cfg%stepa = max(1.0_wp,toreal(kv))
      case ('etemp')
        env%dock%cfg%etemp = toreal(kv)
      case ('pocket')
        env%dock%cfg%pocket_grid = kv%value_b
      case ('stack')
        env%dock%cfg%stack_grid = kv%value_b
      case ('angular')
        env%dock%cfg%angular_grid = kv%value_b
      case ('cs','cssym')
        env%dock%cfg%cssym = kv%value_b
      case ('atm')
        env%dock%cfg%fulle = kv%value_b
      case ('ensemble')
        env%dock%cfg%docking_ens = kv%value_b
      case ('samerand')
        env%dock%cfg%samerand = kv%value_b
      case ('verbose')
        env%dock%cfg%debug = kv%value_b
      case ('chrg','uhf')
        call pair_int(kv,kv%key)
      case ('fix')
        env%dock%fix = kv%value_c
      case ('directed')
        env%dock%directed = kv%value_c
      case ('directed_type')
        if (kv%value_c == 'repulsive') then
          env%dock%directed_type = p_atom_pot
        else
          env%dock%directed_type = p_atom_att
        end if
      case ('directed_scal','directed_fc')
        env%dock%directed_fc = toreal(kv)
      case ('solvent')
        env%dock%cfg%solvent = kv%value_c
        if (.not.allocated(env%dock%cfg%solvmodel)) env%dock%cfg%solvmodel = 'alpb'
      case ('solvmodel','solvation')
        env%dock%cfg%solvmodel = kv%value_c
      case default
        istat = istat+1
        write (stdout,'("unrecognized KEYWORD in ",a," : ",a)') '[dock]-block',kv%key
      end select
    end do

  contains

    function toreal(k) result(r)
      type(keyvalue),intent(in) :: k
      real(wp) :: r
      if (k%id == 2) then
        r = real(k%value_i,wp)
      else
        r = k%value_f
      end if
    end function toreal

    subroutine pair_int(k,key)
      !> "chrg = [a,b]" for A and B, a scalar sets molecule A
      type(keyvalue),intent(in) :: k
      character(len=*),intent(in) :: key
      if (allocated(k%value_ia)) then
        if (size(k%value_ia) >= 1) call set_chrguhf(env%dock,key//'1',k%value_ia(1))
        if (size(k%value_ia) >= 2) call set_chrguhf(env%dock,key//'2',k%value_ia(2))
      else
        call set_chrguhf(env%dock,key//'1',k%value_i)
      end if
    end subroutine pair_int

  end subroutine parse_dock_toml

!========================================================================================!

  subroutine crest_aiss_dock(env,tim)
    !*****************************************************************
    !* Runtype driver: read both fragments and the xtb-style detailed
    !* input, dock B onto A and report the best structure (best.xyz)
    !* and all optimized ones (final_structures.xyz)
    !*****************************************************************
    type(systemdata),intent(inout) :: env
    type(timer),intent(inout)      :: tim
    type(aiss_input) :: dk
    type(coord) :: molA,molB,comb,best
    type(aiss_wall),allocatable :: walls(:)
    type(calcdata) :: ucons
    logical,allocatable :: frozen(:),dirsel(:),atl(:)
    integer,allocatable :: directed(:)
    real(wp) :: eint
    integer :: i,io,T,Tn
    logical :: ex

    call tim%start(14,'aISS docking')
    dk = env%dock
    if (.not.allocated(dk%fileA).and.allocated(env%inputcoords)) dk%fileA = env%inputcoords
    if (.not.allocated(dk%fileA).or..not.allocated(dk%fileB)) then
      write (stdout,'(a)') '**ERROR** docking needs two structures: crest A --dock B'
      call creststop(status_input)
    end if
    do i = 1,2
      if (i == 1) inquire (file=dk%fileA,exist=ex)
      if (i == 2) inquire (file=dk%fileB,exist=ex)
      if (.not.ex) then
        write (stdout,'(a)') '**ERROR** structure file not found for docking'
        call creststop(status_input)
      end if
    end do
    call molA%open(dk%fileA)
    call molB%open(dk%fileB)

!>--- settings not given explicitly fall back to the general CREST ones
    if (.not.allocated(dk%cfg%optlvl)) dk%cfg%optlvl = aiss_level_from_flag(env%gfnver)
    if (.not.allocated(dk%cfg%lmolvl)) dk%cfg%lmolvl = trim(env%aiss_lmolvl)
    if (env%sameRandomNumber) dk%cfg%samerand = .true.
    if (env%gbsa.and..not.allocated(dk%cfg%solvent)) then
      if (allocated(env%solv)) then
        if (len_trim(env%solv) > 0) then
          dk%cfg%solvmodel = 'alpb'
          if (index(env%solv,'gbsa') > 0) dk%cfg%solvmodel = 'gbsa'
          dk%cfg%solvent = trim(env%solvent)
        end if
      end if
    end if

    comb%nat = molA%nat+molB%nat
    comb%at = [molA%at,molB%at]
    allocate (comb%xyz(3,comb%nat))
    comb%xyz(:,1:molA%nat) = molA%xyz
    comb%xyz(:,molA%nat+1:) = molB%xyz
    allocate (frozen(comb%nat),dirsel(comb%nat),source=.false.)
    allocate (walls(0))

    if (allocated(dk%xcontrol)) then
      call read_xcontrol(dk%xcontrol,dk,comb,walls,frozen,dirsel,ucons)
    end if
    if (allocated(dk%fix)) then
      call get_atlist(comb%nat,atl,dk%fix,comb%at)
      frozen = frozen.or.atl
    end if
    if (allocated(dk%directed)) then
      call get_atlist(comb%nat,atl,dk%directed,comb%at)
      dirsel = dirsel.or.atl
    end if
    directed = pack([(i,i=1,comb%nat)],dirsel)

    call print_dock_settings(dk,molA,molB,size(walls),count(frozen),size(directed), &
    &                        ucons%nconstraints)

    call new_ompautoset(env,'max',1,T,Tn)
    call aiss_dock_run(molA,molB,dk%cfg,best,eint,io,walls=walls,frozen=frozen, &
    &                  directed=directed,directed_type=dk%directed_type, &
    &                  directed_fc=dk%directed_fc,ucons=ucons,iunit=stdout)
    call tim%stop(14)
    if (io /= 0) then
      write (stdout,'(/,a)') '**ERROR** aISS docking failed'
      call creststop(status_failed)
    end if
    write (stdout,'(/,1x,a)') 'Lowest structure written to <best.xyz>, all optimized'
    write (stdout,'(1x,a)') 'structures to <final_structures.xyz> (energies in Eh).'
  end subroutine crest_aiss_dock

!========================================================================================!

  subroutine read_xcontrol(fname,dk,comb,walls,frozen,dirsel,ucons)
    !*****************************************************************
    !* xtb-style detailed input of xtb dock:
    !*  $wall      polynomial|logfermi ellipsoid/sphere walls
    !*             (alpha=, temp=, beta=, atom lists)
    !*  $fix       atoms: frozen atoms
    !*  $directed  atoms:/elements:, attractive|repulsive,
    !*             scaling factor=
    !*  $dock      the options of the command line (key=value/switch)
    !*  $constrain passed on to the CREST constraint parser
    !* Atom lists refer to the complex (A first, then B).
    !*****************************************************************
    use parse_xtbinput,only:parse_xtbinputfile,get_xtb_constraint_block
    use parse_datastruct,only:root_object
    character(len=*),intent(in)    :: fname
    type(aiss_input),intent(inout) :: dk
    type(coord),intent(inout)      :: comb
    type(aiss_wall),allocatable,intent(inout) :: walls(:)
    logical,intent(inout)          :: frozen(:),dirsel(:)
    type(calcdata),intent(inout)   :: ucons
    type(root_object),allocatable,target :: dict
    character(len=512) :: line
    character(len=:),allocatable :: blk,key,val
    logical,allocatable :: atl(:)
    type(aiss_wall) :: w
    real(wp) :: r(3),temp,beta
    integer :: ich,io,ie,alpha,i,nr,nw0
    logical :: ex,logfermi

    inquire (file=fname,exist=ex)
    if (.not.ex) then
      write (stdout,'(a,a)') '**ERROR** detailed input not found: ',fname
      call creststop(status_input)
    end if

    !> $wall settings are global in xtb (all walls, any order)
    alpha = 30
    logfermi = .false.
    temp = 300.0_wp
    beta = 6.0_wp
    nw0 = size(walls)
    blk = ''
    open (newunit=ich,file=fname)
    do
      read (ich,'(a)',iostat=io) line
      if (io /= 0) exit
      ie = index(line,'#')
      if (ie > 0) line = line(:ie-1)
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '$') then
        blk = lower(trim(line(2:)))
        cycle
      end if
      ie = scan(line,'=:')
      if (ie > 0) then
        key = lower(trim(adjustl(line(:ie-1))))
        val = trim(adjustl(line(ie+1:)))
      else
        key = lower(trim(line))
        val = ''
      end if
      select case (blk)
      case ('wall')
        select case (key)
        case ('alpha')
          read (val,*,iostat=io) alpha
        case ('potential')
          logfermi = (index(lower(val),'logfermi') > 0)
        case ('temp')
          read (val,*,iostat=io) temp
        case ('beta')
          read (val,*,iostat=io) beta
        case ('sandwich')
          write (stdout,'(a)') '**WARNING** sandwich walls are not supported by the docking, ignored'
        case ('ellipsoid','sphere')
          nr = 3
          if (key == 'sphere') nr = 1
          call split_radii(val,nr,r,io)
          if (io /= 0) then
            write (stdout,'(a,a)') '**WARNING** could not read wall: ',trim(line)
            cycle
          end if
          if (nr == 1) r(2:3) = r(1)
          w%radius = r
          if (allocated(w%list)) deallocate (w%list)
          call get_atlist(comb%nat,atl,after_commas(val,nr),comb%at)
          if (.not.all(atl)) w%list = pack([(i,i=1,comb%nat)],atl)
          walls = [walls,w]
        end select
      case ('fix')
        if (key == 'atoms'.or.key == 'elements') then
          call get_atlist(comb%nat,atl,val,comb%at)
          frozen = frozen.or.atl
        end if
      case ('directed')
        select case (key)
        case ('atoms','elements')
          call get_atlist(comb%nat,atl,val,comb%at)
          dirsel = dirsel.or.atl
        case ('attractive')
          dk%directed_type = p_atom_att
        case ('repulsive')
          dk%directed_type = p_atom_pot
        case ('scaling factor')
          read (val,*,iostat=io) dk%directed_fc
        end select
      case ('dock')
        call set_dock_key(dk,key,val)
      end select
    end do
    close (ich)
    walls(nw0+1:)%alpha = alpha
    walls(nw0+1:)%logfermi = logfermi
    walls(nw0+1:)%temp = temp
    walls(nw0+1:)%beta = beta

!>--- $constrain blocks through the CREST xtb-input parser
    allocate (dict)
    call parse_xtbinputfile(fname,dict)
    do i = 1,dict%nblk
      if (trim(dict%blk_list(i)%header) == 'constrain') then
        call get_xtb_constraint_block(ucons,comb,dict%blk_list(i))
      end if
    end do
  end subroutine read_xcontrol

!========================================================================================!

  subroutine set_dock_key(dk,key,val)
    !> $dock block entries (same names as the command-line options)
    type(aiss_input),intent(inout) :: dk
    character(len=*),intent(in) :: key,val
    integer :: io,idum
    real(wp) :: rdum
    select case (key)
    case ('fast')
      call dk%cfg%set_fast()
    case ('qcg')
      call dk%cfg%set_qcg()
    case ('pocket','nopocket','onlypocket','stack','nostack','onlystack', &
    &     'angular','noangular','onlyangular')
      call set_grids(dk%cfg,key)
    case ('cs')
      dk%cfg%cssym = .true.
    case ('atm')
      dk%cfg%fulle = .true.
    case ('ensemble')
      dk%cfg%docking_ens = .true.
    case ('stepr','stepa')
      read (val,*,iostat=io) rdum
      if (io /= 0) return
      if (key == 'stepr') dk%cfg%stepr = max(1.0_wp,rdum)
      if (key == 'stepa') dk%cfg%stepa = max(1.0_wp,rdum)
    case ('nfinal','maxgen','maxparent','nstack')
      read (val,*,iostat=io) idum
      if (io /= 0) return
      select case (key)
      case ('nfinal')
        dk%cfg%n_opt = idum
      case ('maxgen')
        dk%cfg%maxgen = idum
      case ('maxparent')
        dk%cfg%maxparent = idum
      case ('nstack')
        dk%cfg%mxcma = idum
      end select
    end select
  end subroutine set_dock_key

!========================================================================================!

  subroutine set_grids(cfg,flag)
    !> xtb dock search-grid switches, applied in the given order
    type(aiss_config),intent(inout) :: cfg
    character(len=*),intent(in) :: flag
    select case (flag)
    case ('pocket')
      cfg%pocket_grid = .true.
    case ('nopocket')
      cfg%pocket_grid = .false.
    case ('stack')
      cfg%stack_grid = .true.
    case ('nostack')
      cfg%stack_grid = .false.
    case ('angular')
      cfg%angular_grid = .true.
    case ('noangular')
      cfg%angular_grid = .false.
    case ('onlypocket')
      cfg%pocket_grid = .true.
      cfg%stack_grid = .false.
      cfg%angular_grid = .false.
    case ('onlystack')
      cfg%pocket_grid = .false.
      cfg%stack_grid = .true.
      cfg%angular_grid = .false.
    case ('onlyangular')
      cfg%pocket_grid = .false.
      cfg%stack_grid = .false.
      cfg%angular_grid = .true.
    end select
  end subroutine set_grids

!========================================================================================!

  subroutine set_chrguhf(dk,key,ival)
    type(aiss_input),intent(inout) :: dk
    character(len=*),intent(in) :: key
    integer,intent(in) :: ival
    select case (key)
    case ('chrg','chrg1')
      dk%cfg%chrg(1) = ival
      dk%chrg_set(1) = .true.
    case ('chrg2')
      dk%cfg%chrg(2) = ival
      dk%chrg_set(2) = .true.
    case ('uhf','uhf1')
      dk%cfg%uhf(1) = ival
      dk%uhf_set(1) = .true.
    case ('uhf2')
      dk%cfg%uhf(2) = ival
      dk%uhf_set(2) = .true.
    end select
  end subroutine set_chrguhf


!========================================================================================!

  subroutine print_dock_settings(dk,molA,molB,nwalls,nfrozen,ndirected,ncons)
    type(aiss_input),intent(in) :: dk
    type(coord),intent(in) :: molA,molB
    integer,intent(in) :: nwalls,nfrozen,ndirected,ncons
    character(len=*),parameter :: f1 = '(2x,a,t30,a)'
    character(len=*),parameter :: f2 = '(2x,a,t30,i0,1x,a,i0,a,i0,a)'
    write (stdout,'(/,1x,a)') 'aISS docking settings'
    write (stdout,f2) 'molecule A:',molA%nat,'atoms (chrg ',dk%cfg%chrg(1),', uhf ',dk%cfg%uhf(1),')'
    write (stdout,f1) '  file',dk%fileA
    write (stdout,f2) 'molecule B:',molB%nat,'atoms (chrg ',dk%cfg%chrg(2),', uhf ',dk%cfg%uhf(2),')'
    write (stdout,f1) '  file',dk%fileB
    write (stdout,f1) 'final optimizations:',dk%cfg%optlvl
    write (stdout,f1) 'LMO/IFF precomputation:',dk%cfg%lmolvl
    if (allocated(dk%cfg%solvent)) write (stdout,f1) 'implicit solvation:', &
    &  trim(dk%cfg%solvmodel)//' '//trim(dk%cfg%solvent)
    if (allocated(dk%xcontrol)) write (stdout,f1) 'detailed input:',dk%xcontrol
    if (nwalls > 0) write (stdout,'(2x,a,t30,i0)') 'wall potentials:',nwalls
    if (nfrozen > 0) write (stdout,'(2x,a,t30,i0)') 'frozen atoms:',nfrozen
    if (ndirected > 0) write (stdout,'(2x,a,t30,i0)') 'directed docking atoms:',ndirected
    if (ncons > 0) write (stdout,'(2x,a,t30,i0)') 'user constraints:',ncons
  end subroutine print_dock_settings

!========================================================================================!

  function nodash(a) result(f)
    !> option name without leading dashes, lower case
    character(len=*),intent(in) :: a
    character(len=:),allocatable :: f
    integer :: k
    k = verify(a,'-')
    if (k == 0) then
      f = ''
    else
      f = lower(trim(a(k:)))
    end if
  end function nodash

  function level_name(a) result(lvl)
    !> 'gfn2','--gfn2','2','ff','gff',... -> gfn2/gfn1/gfn0/gfnff
    character(len=*),intent(in) :: a
    character(len=:),allocatable :: lvl,b
    b = nodash(a)
    select case (b)
    case ('ff','gff','gfnff')
      lvl = 'gfnff'
    case ('0','gfn0')
      lvl = 'gfn0'
    case ('1','gfn1')
      lvl = 'gfn1'
    case default
      lvl = 'gfn2'
    end select
  end function level_name

  function lower(a) result(b)
    character(len=*),intent(in) :: a
    character(len=len(a)) :: b
    integer :: k
    b = a
    do k = 1,len(a)
      if (a(k:k) >= 'A'.and.a(k:k) <= 'Z') b(k:k) = achar(iachar(a(k:k))+32)
    end do
  end function lower

  subroutine split_radii(val,nr,r,io)
    character(len=*),intent(in) :: val
    integer,intent(in) :: nr
    real(wp),intent(out) :: r(3)
    integer,intent(out) :: io
    r = 0.0_wp
    read (val,*,iostat=io) r(1:nr)
  end subroutine split_radii

  function after_commas(val,n) result(tail)
    !> the text after the n-th comma (the atom list of a wall)
    character(len=*),intent(in) :: val
    integer,intent(in) :: n
    character(len=:),allocatable :: tail
    integer :: k,ic
    tail = val
    do k = 1,n
      ic = index(tail,',')
      if (ic == 0) then
        tail = 'all'
        return
      end if
      tail = tail(ic+1:)
    end do
    tail = trim(adjustl(tail))
    if (len_trim(tail) == 0) tail = 'all'
  end function after_commas

!========================================================================================!
end module aiss_standalone
