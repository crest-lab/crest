!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022 Philipp Pracht
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
!> New implementation of protonation routines
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
subroutine crest_new_protonate(env,tim)
!********************************************************************
!* Standalone runtype for protonating a molecule
!*
!* Input/Output:
!*  env  -  crest's systemdata object
!*  tim  -  timer object
!********************************************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
  use parallel_interface
  use cregen_interface
  use iomod
  use utilities,only:binomial
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,T,Tn,np
  logical :: pr,wr
  character(len=80) :: atmp
!========================================================================================!
  type(calcdata) :: calc
  real(wp) :: accuracy,etemp

  real(wp) :: energy,dip
  real(wp),allocatable :: grad(:,:)
  type(calcdata),allocatable :: tmpcalc
  type(calcdata),allocatable :: tmpcalc_ff
  type(calculation_settings) :: tmpset
  real(wp),allocatable :: protxyz(:,:)
  integer :: natp,pstep,npnew
  integer,allocatable :: atp(:)
  real(wp),allocatable :: xyzp(:,:,:)
  real(wp),allocatable :: ep(:)
  logical,allocatable :: atlist(:)
  character(len=*),parameter :: basename = 'protonate_'
!========================================================================================!
  write (stdout,*)
  !call system('figlet singlepoint')
  write (stdout,*) "                 _                    _        "
  write (stdout,*) " _ __  _ __ ___ | |_ ___  _ __   __ _| |_ ___  "
  write (stdout,*) "| '_ \| '__/ _ \| __/ _ \| '_ \ / _` | __/ _ \ "
  write (stdout,*) "| |_) | | | (_) | || (_) | | | | (_| | ||  __/ "
  write (stdout,*) "| .__/|_|  \___/ \__\___/|_| |_|\__,_|\__\___| "
  write (stdout,*) "|_|                                            "
  write (stdout,*) "-----------------------------------------------"
  write (stdout,*) "  automated protonation site screening script  "
  write (stdout,*) "       revised version (c) P.Pracht 2024       "
  write (stdout,*) 'Cite as:'
  write (stdout,*) '  P.Pracht, C.A.Bauer, S.Grimme'
  write (stdout,*) '  JCC, 2017, 38, 2618–2631.'
  write (stdout,*)

!========================================================================================!
  call new_ompautoset(env,'max',0,T,Tn)
  call ompprint_intern()
!========================================================================================!
  call env%ref%to(mol)
  write (stdout,*)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)
!========================================================================================!
!>--- The first step is to perfom a LMO calculation to identify suitable protonation sites
!>--- in this version of the program this step is always done with GFN0-xTB.

  call tim%start(14,'LMO center calculation')
  write (stdout,'(a)') repeat('-',80)
  write (stdout,'(a)')
  write (stdout,'(a)',advance='no') '> Setting up GFN0-xTB for LMO center calculation ... '
  flush (stdout)

  allocate (grad(3,mol%nat),source=0.0_wp)
!>--- GFN0 job adapted from global settings
  allocate (tmpcalc)
  call env2calc(env,tmpcalc,mol)
  tmpcalc%calcs(1)%id = jobtype%gfn0
  tmpcalc%calcs(1)%rdwbo = .true.
  tmpcalc%calcs(1)%getlmocent = .true.
  tmpcalc%calcs(1)%chrg = env%chrg
  tmpcalc%ncalculations = 1
  write (stdout,'(a)') 'done.'
!>--- and then start it
  write (stdout,'(a)',advance='yes') '> Performing singlepoint calculation ... '
  flush (stdout)
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<!
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<!
  call engrad(mol,tmpcalc,energy,grad,io)
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<!
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<!
  call tim%stop(14)
  write (stdout,'(a)') '> done.'
  write (atmp,'(a)') '> Total wall time for calculations'
  call tim%write_timing(stdout,14,trim(atmp),.true.)
  write (stdout,'(a)') repeat('-',80)
  write (stdout,*)
  if (io /= 0) then
    write (stdout,*)
    write (stdout,*) 'WARNING: Calculation exited with error!'
  end if
!>--- check LMO center
  np = tmpcalc%calcs(1)%nprot
  if (np > 0) then
    write (stdout,'(a,i0,a)') '> ',np,' π- or LP-centers identified as protonation candidates:'
    call move_alloc(tmpcalc%calcs(1)%protxyz,protxyz)
    write (stdout,'(1x,a5,1x,a,5x,a)') 'LMO','type','center(xyz/Ang)'
    do i = 1,np
      select case (nint(protxyz(4,i)))
      case (1)
        write (stdout,'(1x,i5,1x,a,3F12.5)') i,'LP   ',protxyz(1:3,i)*autoaa
      case (2)
        write (stdout,'(1x,i5,1x,a,3F12.5)') i,'π    ',protxyz(1:3,i)*autoaa
      case (3)
        write (stdout,'(1x,i5,1x,a,3F12.5)') i,'del.π',protxyz(1:3,i)*autoaa
      case default
        write (stdout,'(1x,i5,1x,a,3F12.5)') i,'???  ',protxyz(1:3,i)*autoaa
      end select
    end do
  else
    write (stdout,*)
    write (stdout,*) 'WARNING: No suitable protonation sites found!'
    write (stdout,*) '         Confirm whether you expect π- or LP-centers for your molecule!'
    write (stdout,*)
    env%iostatus_meta = status_failed
    return
  end if
  deallocate (tmpcalc)
  deallocate (grad)

!========================================================================================!
!>--- If we reached this point, we have candidate positions for our protonation!
  pstep = 0
  write (stdout,'(a)',advance='yes') '> Generating candidate structures ... '
  flush (stdout)
  if (env%protb%amount > np) then
    env%protb%amount = np
    write (stdout,'(a,i0,a)') '> A maximum of ',np,' protonation sites can be used'
  end if
  natp = mol%nat+env%protb%amount
  npnew = binomial(np,env%protb%amount)
  write (stdout,'(a,i0,a)') '> Up to ',npnew,' new structures'
  allocate (atp(natp),source=0)
  allocate (xyzp(3,natp,npnew),ep(npnew),source=0.0_wp)

  call protonation_candidates(env,mol,natp,np,protxyz,atp,xyzp,npnew)
!>--- NOTE: after this the global charge (env%chrg) and all charges saved in the calc levels have been increased

  if (npnew < 1) then
    write (stdout,*)
    write (stdout,'(a)') '> WARNING: No remaining protonation sites after applying user defined conditions!'
    write (stdout,'(a)') '>          Modify the search criteria and check your input structure for sanity.'
    env%iostatus_meta = status_failed
    return
  end if

  write (atmp,'(a,i0,a)') basename,pstep,'.xyz'
  write (stdout,'(a,a,a,i0,a)') '> Write ',trim(atmp),' with ',npnew,' candidates ... '

  call wrensemble(basename//'0.xyz',natp,npnew,atp,xyzp(:,:,1:npnew)*autoaa,ep(1:npnew))

  write (stdout,'(a)') '> done.'
  write (stdout,*)

!========================================================================================!
!>--- Optimize candidates, optional FF pre-step
  if (env%protb%ffopt) then
    call smallhead('Protomer Ensemble FF Pre-Optimization')
    !>--- set up a temporary calculation object
    write (stdout,'(a)') '> LMO centers can be very close to atoms, leading to initial'
    write (stdout,'(a)') '  extremely high-energy candidates which may undergo unwanted'
    write (stdout,'(a)') '  chemical changes in optimizations. Classical force-fields'
    write (stdout,'(a)') '  with defined bond-topology can help circumvent this issue.'
    write (stdout,'(a)') '> Setting up force-field structure pre-optimization ...'
    allocate (tmpcalc_ff)
    tmpcalc_ff%optnewinit = .true.
    env%calc%optnewinit = .true.
    tmpcalc_ff%optnewinit = .true.
    call tmpset%create('gfnff')
    tmpset%chrg = env%chrg
    call tmpcalc_ff%add(tmpset)
    tmpcalc_ff%maxcycle = 10000
    call tmpcalc_ff%info(stdout)

    !>--- Run optimizations
    call tim%start(15,'Ensemble optimization (FF)')
    call print_opt_data(tmpcalc_ff,stdout)
    write (stdout,'(a,i0,a)') '> ',npnew,' structures to optimize ...'
    call crest_oloop(env,natp,npnew,atp,xyzp(:,:,1:npnew),ep(1:npnew),.false.,tmpcalc_ff)
    call tim%stop(15)

    deallocate (tmpcalc_ff)

    pstep = pstep+1
    write (atmp,'(a,i0,a)') basename,pstep,'.xyz'
    write (stdout,'(a,a,a)') '> Write ',trim(atmp),' with optimized structures ... '
    call wrensemble(trim(atmp),natp,npnew,atp,xyzp(:,:,1:npnew)*autoaa,ep(1:npnew))
    deallocate (xyzp,atp) !> clear this space to re-use it

    !>--- sorting
    write (stdout,'(a)') '> Sorting structures by energy to remove failed opts. ...'
    env%ewin = 1000.0_wp
    call newcregen(env,7,trim(atmp))
    call rename(trim(atmp)//'.sorted',trim(atmp))
    write (stdout,'(a)') '> Sorted file was renamed to '//trim(atmp)
    write (stdout,'(a)') '> WARNING: These are force-field energies which are '
    write (stdout,'(a)') '           NOT(!) accurate for bond formation and breaking!'
    write (stdout,*)

    !>--- re-read sorted ensemble
    call rdensemble(trim(atmp),natp,npnew,atp,xyzp)
    xyzp = xyzp*aatoau !> don't forget to restore BOHR
  end if

!========================================================================================!
!>--- H-position-only optimization (only makes sense after FF preoptimization)
  if (env%protb%hnewopt.and.env%protb%ffopt) then
    call smallhead('Protomer Ensemble Frozen-Atom Optimization')
    !>--- create temporary calculation from our intended level of theory
    write (stdout,'(a)') '> Setting up frozen structure optimization ...'
    allocate (tmpcalc)
    call env2calc(env,tmpcalc,mol)
    !>--- freeze all atoms, EXCEPT the new one
    allocate (atlist(natp),source=.true.)
    atlist(natp) = .false. !> the new one is always last
    do i = 1,natp
      if (atp(i) == 1) then
        atlist(i) = .false. !> additionally un-freeze all H's (this seems to be beneficial)
      end if
    end do
    tmpcalc%nfreeze = count(atlist)
    tmpcalc%optnewinit = .true.
    write (stdout,'(a,i0,a)') '> ',tmpcalc%nfreeze,' frozen atoms. All H non-frozen.'
    call move_alloc(atlist,tmpcalc%freezelist)
    call tmpcalc%info(stdout)

    !>--- run opt
    call tim%start(16,'Ensemble optimization (frozen)')
    write (stdout,'(a,i0,a)') '> ',npnew,' structures to optimize ...'
    call crest_oloop(env,natp,npnew,atp,xyzp(:,:,1:npnew),ep(1:npnew),.false.,tmpcalc)
    call tim%stop(16)

    pstep = pstep+1
    write (atmp,'(a,i0,a)') basename,pstep,'.xyz'
    write (stdout,'(a,a,a)') '> Write ',trim(atmp),' with optimized structures ... '
    call wrensemble(trim(atmp),natp,npnew,atp,xyzp(:,:,1:npnew)*autoaa,ep(1:npnew))
    deallocate (xyzp,atp) !> clear this space to re-use it
    deallocate (tmpcalc)

!    call tim%start(17,'Ensemble refinement')
!    call crest_refine(env,trim(atmp),trim(atmp))
!    call tim%stop(17)

    !>--- sorting
    write (stdout,'(a)') '> Sorting structures by energy to remove failed opts. ...'
    env%ewin = env%protb%ewin*5.0_wp  !> large energy threshold
    call newcregen(env,7,trim(atmp))
    call rename(trim(atmp)//'.sorted',trim(atmp))
    write (stdout,'(a)') '> Sorted file was renamed to '//trim(atmp)
    write (stdout,*)

    !>--- re-read sorted ensemble
    call rdensemble(trim(atmp),natp,npnew,atp,xyzp)
    xyzp = xyzp*aatoau !> don't forget to restore BOHR

  end if

!========================================================================================!
!>--- Optimize with global settings
  if (env%protb%finalopt.and.env%protb%ffopt) then
    call smallhead('Final Protomer Ensemble Optimization')
    allocate (tmpcalc)
    call env2calc(env,tmpcalc,mol)
    call tmpcalc%info(stdout)
    tmpcalc%optnewinit = .true.
    call tim%start(20,'Ensemble optimization')
    call print_opt_data(env%calc,stdout)
    write (stdout,'(a,i0,a)') '> ',npnew,' structures to optimize ...'
    call crest_oloop(env,natp,npnew,atp,xyzp(:,:,1:npnew),ep,.false.,tmpcalc)
    call tim%stop(20)

    pstep = pstep+1
    write (atmp,'(a,i0,a)') basename,pstep,'.xyz'
    write (stdout,'(a,a,a)') '> Write ',trim(atmp),' with optimized structures ... '
    call wrensemble(trim(atmp),natp,npnew,atp,xyzp(:,:,1:npnew)*autoaa,ep(1:npnew))
    deallocate (xyzp,atp)

    call tim%start(17,'Ensemble refinement')
    call crest_refine(env,trim(atmp),trim(atmp))
    call tim%stop(17)

!>--- sorting
    write (stdout,'(a)') '> Sorting structures by energy ...'
    env%ewin = env%protb%ewin*3.0_wp
    call newcregen(env,7,trim(atmp))
    call rename(trim(atmp)//'.sorted',trim(atmp))
  end if

!========================================================================================!
!> Remove doubly generated tautomers
  call protonation_prep_canonical(env,mol,trim(atmp))

!========================================================================================!
!>--- move final ensemble to protonated.xyz
  call rename(trim(atmp),'protonated.xyz')
  write (stdout,'(a)') '> Sorted file was renamed to protonated.xyz'

  write (stdout,'(a)') '> All other temporary protonate_*.xyz files are removed by default.'
  if (.not.env%keepmodef) then
    call rmrf('protonate_*.xyz')
  end if
!========================================================================================!
  return
end subroutine crest_new_protonate

!========================================================================================!
!========================================================================================!

subroutine protonation_candidates(env,mol,natp,np,protxyz,at,xyz,npnew)
!********************************************************
!* generate protonation/ionization candidate structures
!* The outputs are at and xyz, the latter being in Bohr
!********************************************************
  use crest_data
  use crest_parameters
  use strucrd,only:coord
  use utilities,only:binomial,get_combinations
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env
  type(coord),intent(in) :: mol
  integer,intent(in) :: natp
  integer,intent(in) :: np
  real(wp),intent(in) :: protxyz(4,np)
  !> OUTPUT
  integer,intent(out)  :: at(natp)
  real(wp),intent(out) :: xyz(3,natp,npnew)
  integer,intent(inout) :: npnew
  !> LOCAL
  integer :: i,j,jj,k,l,ii,c,nc,kk
  integer :: ati,ichrg,ctype
  integer,allocatable :: combi(:,:),tmp(:)

  if (natp .ne. mol%nat+env%protb%amount) then
    write (stdout,'(a)') 'WARNING: Inconsistent number of atoms in protonation routine'
    call creststop(status_args)
  end if

  if (env%protb%swelem) then
!>--- User-defined monoatomic ion
    ichrg = env%protb%swchrg*env%protb%amount
    ati = env%protb%swat
  else
!>--- DEFAULT: H⁺
    ichrg = 1*env%protb%amount
    ati = 1
  end if
  write (stdout,'(a,i0)') '> Increasing the molecular charge by ',ichrg
  call env%calc%increase_charge(ichrg)
  env%chrg = env%chrg+ichrg

!>--- Check if we have some other conditions
  if (any(.not.env%protb%active_lmo(:))) then
    write (stdout,'(a)',advance='no') '> User-defined: IGNORING '
    if (.not.env%protb%active_lmo(1)) write (stdout,'(a)',advance='no') 'π '
    if (.not.env%protb%active_lmo(2)) write (stdout,'(a)',advance='no') 'LP '
    if (.not.env%protb%active_lmo(3)) write (stdout,'(a)',advance='no') 'deloc.π '
    write (stdout,'(a)') 'LMOs ...'
  end if

!>--- check combinations
  k = env%protb%amount
  nc = binomial(np,k)
  allocate (tmp(k),combi(k,nc),source=0)
  c = 0
  call get_combinations(np,k,nc,c,combi,tmp,0)

!>--- Populate
  npnew = 0
  ii = 0
  COMBILOOP: do i = 1,nc
!>--- Enforce further constraints, conditions, etc.
    ADDLOOP1: do j = 1,k
      jj = combi(j,i)
      ctype = nint(protxyz(4,jj))
      if (.not.env%protb%active_lmo(ctype)) cycle COMBILOOP
    end do ADDLOOP1

!>--- passed checks in ADDLOOP1 means we can add this config
    ii = ii+1 !> counter of actually created structures
    do j = 1,mol%nat
      xyz(1:3,j,ii) = mol%xyz(1:3,j)
      at(j) = mol%at(j)
    end do
    kk = mol%nat
    ADDLOOP2: do j = 1,k
      jj = combi(j,i)
      kk = kk+1
      xyz(1:3,kk,ii) = protxyz(1:3,jj)
      at(kk) = ati
    end do ADDLOOP2
  end do COMBILOOP
  npnew = ii

  deallocate (combi,tmp)
end subroutine protonation_candidates

!========================================================================================!

subroutine protonation_prep_canonical(env,refmol,fname)
  use crest_parameters
  use crest_data
  use strucrd
  use crest_calculator
  use tblite_api,only:xtblvl
  use canonical_mod
  use iomod,only:remove
  use adjacency
  use cregen_interface
  implicit none
  type(systemdata) :: env
  type(coord),intent(in) :: refmol 
  character(len=*),intent(in) :: fname

  type(calcdata),allocatable :: tmpcalc
  type(calculation_settings) :: ceh
  type(coord),allocatable :: structures(:)

  real(wp),allocatable :: cn(:),Bmat(:,:)
  integer,allocatable :: frag(:),Amat(:,:)
  real(wp),allocatable :: grad(:,:)
  real(wp) :: energy
  integer :: nat,nall,i,j,k,io,ich,refnfrag

  type(canonical_sorter),allocatable :: canon(:)
  integer,allocatable :: group(:)
  character(len=*),parameter :: outfile = 'unique.xyz'

  call smallhead('Remove duplicates via canonical atom identity')

!>--- reference molecule (unprotonated)
!  call env%ref%to(refmol) !> DO NOT USE THIS, CREGEN MAY OVERWRITE
  call refmol%cn_to_bond(cn,Bmat,'cov')
  call wbo2adjacency(refmol%nat,Bmat,Amat,0.01_wp)
  allocate (frag(refmol%nat),source=0)
  call setup_fragments(refmol%nat,Amat,frag)
  refnfrag = maxval(frag(:),1)
  deallocate (frag,Amat,Bmat)

!>--- read ensemble
  call rdensemble(fname,nall,structures)
  nat = structures(1)%nat

!>--- run singlepoints to document WBOs and charges (doesn't need to be parallel)
  allocate (canon(nall))

  write (stdout,'(a,i0,a)') '> Setting up canonical atom order for ',nall,' structures via CN-based molecular graphs ...'
  call crest_oloop_pr_progress(env,nall,0)
  do i = 1,nall
    call canon(i)%init(structures(i),invtype='apsp+')
    call canon(i)%stereo(structures(i))
    !write(*,*)
    !call canon(i)%rankprint(structures(i))
    call canon(i)%shrink()

    call crest_oloop_pr_progress(env,nall,i)
  end do
  call crest_oloop_pr_progress(env,nall,-1)

!>--- grouping loop
  allocate (group(nall),source=0)
  do i = 1,nall
    k = maxval(group(:),1)
    if (group(i) == 0) then
      k = k+1
      group(i) = k
      !write(*,*)'structure ',i,' is the reference of ',k
    end if
    do j = i+1,nall
      if (group(j) == 0) then
        if (canon(i)%compare(canon(j))) then
          group(j) = k
          !write(*,*)'structure ',j,' is in group ',k
        end if
      end if
    end do
  end do

  write (stdout,'(a,i0,a)') '> ',maxval(group,1),' groups identified based on canonical atom identifier algorithm'

!>--- more rules (setting group(i) to zero will discard the structure)
  k = 0
  do i = 1,nall
    if (canon(i)%nfrag .ne. refnfrag) then
      group(i) = 0
      k = k+1
    end if
  end do
  if (k > 0) then
    write (stdout,'(a,i0,a,i0,a)') '> ',k,' structures discared due to fragment change (!=',refnfrag,')'
  end if

!>--- dump to new file
  open (newunit=ich,file=outfile)
  GROUPLOOP: do i = 1,maxval(group,1)
    STRUCTLOOP: do j = 1,nall
      if (group(j) == i) then !> since the structures are energy sorted, writing the first is taking the lowest in the group
        call structures(j)%append(ich)
        cycle GROUPLOOP
      end if
    end do STRUCTLOOP
  end do GROUPLOOP
  close (ich)
!>--- "sort" again for printout
  call newcregen(env,7,trim(outfile))
  call remove(outfile)
  call rename(outfile//'.sorted',fname)

  deallocate (group)
  deallocate (canon)
end subroutine protonation_prep_canonical

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> New implementation of DEprotonation routines
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
subroutine crest_new_deprotonate(env,tim)
!********************************************************************
!* Standalone runtype for deprotonating a molecule
!*
!* Input/Output:
!*  env  -  crest's systemdata object
!*  tim  -  timer object
!********************************************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
  use parallel_interface
  use cregen_interface
  use iomod
  use utilities,only:binomial
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,T,Tn,np
  logical :: pr,wr
  character(len=80) :: atmp
!========================================================================================!
  type(calcdata) :: calc
  real(wp) :: accuracy,etemp

  real(wp) :: energy,dip
  real(wp),allocatable :: grad(:,:)
  type(calcdata),allocatable :: tmpcalc
  type(calcdata),allocatable :: tmpcalc_ff
  type(calculation_settings) :: tmpset
  integer :: natp,pstep,npnew
  integer,allocatable :: atp(:)
  real(wp),allocatable :: xyzp(:,:,:)
  real(wp),allocatable :: ep(:)
  logical,allocatable :: atlist(:)
  character(len=*),parameter :: basename = 'deprotonate_'
!========================================================================================!
  write (stdout,*)
  !call system('figlet singlepoint')
  write (stdout,*) "     _                      _                    _        "
  write (stdout,*) "  __| | ___ _ __  _ __ ___ | |_ ___  _ __   __ _| |_ ___  "
  write (stdout,*) " / _` |/ _ \ '_ \| '__/ _ \| __/ _ \| '_ \ / _` | __/ _ \ "
  write (stdout,*) "| (_| |  __/ |_) | | | (_) | || (_) | | | | (_| | ||  __/ "
  write (stdout,*) " \__,_|\___| .__/|_|  \___/ \__\___/|_| |_|\__,_|\__\___| "
  write (stdout,*) "           |_|                                            "
  write (stdout,*) "--------------------------------------------------------- "
  write (stdout,*) "      automated deprotonation site screening script       "
  write (stdout,*) "            revised version (c) P.Pracht 2024             "
  write (stdout,*) 'Cite as:'
  write (stdout,*) '  P.Pracht, R.Wilcken, A.Udvarhelyi, S.Rodde, S.Grimme'
  write (stdout,*) '  JCAMD, 2018, 32, 1139-1149.'
  write (stdout,*)

!========================================================================================!
  call new_ompautoset(env,'max',0,T,Tn)
  call ompprint_intern()
!========================================================================================!
  call env%ref%to(mol)
  write (stdout,*)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)
!========================================================================================!
!>--- The first step is to sort the input structure so that all heavy atoms are first and
!>--- all hydrogen atoms follow after that. This is because some of the (parallel) routines
!>--- rely on identical atom order.

  write (stdout,'(a)') repeat('-',80)
  write (stdout,'(a)')
  write (stdout,'(a)') '> Re-sorting input structure ... '
  call deprotonation_sort(mol,np)
  call mol%append(stdout)
  write (stdout,*)
!>--- check
  if (np == 0) then
    write (stdout,*)
    write (stdout,*) 'WARNING: No deprotonation sites found!'
    write (stdout,*)
    return
  end if

!========================================================================================!
!>--- If we reached this point, we have candidate positions for our protonation!
  pstep = 0
  write (stdout,'(a)',advance='yes') '> Generating candidate structures ... '
  flush (stdout)
  if (env%protb%amount > np) then
    env%protb%amount = np
    write (stdout,'(a,i0,a)') '> A maximum of ',np,' deprotonation sites can be used'
  end if
  natp = mol%nat-env%protb%amount
  npnew = binomial(np,env%protb%amount)
  write (stdout,'(a,i0,a)') '> Up to ',npnew,' new structures'
  allocate (atp(natp),source=0)
  allocate (xyzp(3,natp,npnew),ep(npnew),source=0.0_wp)

  call deprotonation_candidates(env,mol,natp,np,atp,xyzp,npnew)
!>--- NOTE: after this the global charge (env%chrg) and all charges saved in the calc levels have been decreased

  if (npnew < 1) then
    write (stdout,*)
    write (stdout,'(a)') '> WARNING: No remaining deprotonation sites after applying user defined conditions!'
    write (stdout,'(a)') '>          Modify the search criteria and check your input structure for sanity.'
    return
  end if

  write (atmp,'(a,i0,a)') basename,pstep,'.xyz'
  write (stdout,'(a,a,a,i0,a)') '> Write ',trim(atmp),' with ',npnew,' candidates ... '

  call wrensemble(basename//'0.xyz',natp,npnew,atp,xyzp(:,:,1:npnew)*autoaa,ep(1:npnew))

  write (stdout,'(a)') '> done.'
  write (stdout,*)

!========================================================================================!
!>--- Optimize candidates, optional FF pre-step
  if (env%protb%ffopt) then
    call smallhead('Deprotomer Ensemble FF Pre-Optimization')
    !>--- set up a temporary calculation object
    write (stdout,'(a)') '> Removal of some hydrogen atoms may lead to initial'
    write (stdout,'(a)') '  extremely high-energy candidates which can undergo unwanted'
    write (stdout,'(a)') '  chemical changes in optimizations. Classical force-fields'
    write (stdout,'(a)') '  with defined bond-topology can help circumvent this issue.'
    write (stdout,'(a)') '> Setting up force-field structure pre-optimization ...'
    allocate (tmpcalc_ff)
    tmpcalc_ff%optnewinit = .true.
    env%calc%optnewinit = .true.
    tmpcalc_ff%optnewinit = .true.
    call tmpset%create('gfnff')
    tmpset%chrg = env%chrg
    call tmpcalc_ff%add(tmpset)
    tmpcalc_ff%maxcycle = 10000
    call tmpcalc_ff%info(stdout)

    !>--- Run optimizations
    call tim%start(15,'Ensemble optimization (FF)')
    call print_opt_data(tmpcalc_ff,stdout)
    write (stdout,'(a,i0,a)') '> ',npnew,' structures to optimize ...'
    call crest_oloop(env,natp,npnew,atp,xyzp(:,:,1:npnew),ep(1:npnew),.false.,tmpcalc_ff)
    call tim%stop(15)

    deallocate (tmpcalc_ff)

    pstep = pstep+1
    write (atmp,'(a,i0,a)') basename,pstep,'.xyz'
    write (stdout,'(a,a,a)') '> Write ',trim(atmp),' with optimized structures ... '
    call wrensemble(trim(atmp),natp,npnew,atp,xyzp(:,:,1:npnew)*autoaa,ep(1:npnew))
    deallocate (xyzp,atp) !> clear this space to re-use it

    !>--- sorting
    write (stdout,'(a)') '> Sorting structures by energy to remove failed opts. ...'
    env%ewin = 1000.0_wp
    call newcregen(env,7,trim(atmp))
    call rename(trim(atmp)//'.sorted',trim(atmp))
    write (stdout,'(a)') '> Sorted file was renamed to '//trim(atmp)
    write (stdout,'(a)') '> WARNING: These are force-field energies which are '
    write (stdout,'(a)') '           NOT(!) accurate for bond formation and breaking!'
    write (stdout,*)

    !>--- re-read sorted ensemble
    call rdensemble(trim(atmp),natp,npnew,atp,xyzp)
    xyzp = xyzp*aatoau !> don't forget to restore BOHR
  end if

!========================================================================================!
!>--- H-position-only optimization (only makes sense after FF preoptimization)
  if (env%protb%hnewopt.and.env%protb%ffopt) then
    call smallhead('Deprotomer Ensemble Frozen-Atom Optimization')
    !>--- create temporary calculation from our intended level of theory
    write (stdout,'(a)') '> Setting up frozen structure optimization ...'
    allocate (tmpcalc)
    call env2calc(env,tmpcalc,mol)
    !>--- freeze all atoms, EXCEPT the new one
    allocate (atlist(natp),source=.true.)
    atlist(natp) = .false. !> the new one is always last
    do i = 1,natp
      if (atp(i) == 1) then
        atlist(i) = .false. !> additionally un-freeze all H's (this seems to be beneficial)
      end if
    end do
    tmpcalc%nfreeze = count(atlist)
    tmpcalc%optnewinit = .true.
    write (stdout,'(a,i0,a)') '> ',tmpcalc%nfreeze,' frozen atoms. All H non-frozen.'
    call move_alloc(atlist,tmpcalc%freezelist)
    call tmpcalc%info(stdout)

    !>--- run opt
    call tim%start(16,'Ensemble optimization (frozen)')
    write (stdout,'(a,i0,a)') '> ',npnew,' structures to optimize ...'
    call crest_oloop(env,natp,npnew,atp,xyzp(:,:,1:npnew),ep(1:npnew),.false.,tmpcalc)
    call tim%stop(16)

    pstep = pstep+1
    write (atmp,'(a,i0,a)') basename,pstep,'.xyz'
    write (stdout,'(a,a,a)') '> Write ',trim(atmp),' with optimized structures ... '
    call wrensemble(trim(atmp),natp,npnew,atp,xyzp(:,:,1:npnew)*autoaa,ep(1:npnew))
    deallocate (xyzp,atp) !> clear this space to re-use it
    deallocate (tmpcalc)

!    call tim%start(17,'Ensemble refinement')
!    call crest_refine(env,trim(atmp),trim(atmp))
!    call tim%stop(17)

    !>--- sorting
    write (stdout,'(a)') '> Sorting structures by energy to remove failed opts. ...'
    env%ewin = env%protb%ewin*5.0_wp  !> large energy threshold
    call newcregen(env,7,trim(atmp))
    call rename(trim(atmp)//'.sorted',trim(atmp))
    write (stdout,'(a)') '> Sorted file was renamed to '//trim(atmp)
    write (stdout,*)

    !>--- re-read sorted ensemble
    call rdensemble(trim(atmp),natp,npnew,atp,xyzp)
    xyzp = xyzp*aatoau !> don't forget to restore BOHR

  end if

!========================================================================================!
!>--- Optimize with global settings
  if (env%protb%finalopt.and.env%protb%ffopt) then
    call smallhead('Final Deprotomer Ensemble Optimization')
    allocate (tmpcalc, source=env%calc)
    call tmpcalc%info(stdout)
    tmpcalc%optnewinit = .true.
    call tim%start(20,'Ensemble optimization')
    call print_opt_data(env%calc,stdout)
    write (stdout,'(a,i0,a)') '> ',npnew,' structures to optimize ...'
    call crest_oloop(env,natp,npnew,atp,xyzp(:,:,1:npnew),ep(1:npnew),.false.,tmpcalc)
    call tim%stop(20)

    pstep = pstep+1
    write (atmp,'(a,i0,a)') basename,pstep,'.xyz'
    write (stdout,'(a,a,a)') '> Write ',trim(atmp),' with optimized structures ... '
    call wrensemble(trim(atmp),natp,npnew,atp,xyzp(:,:,1:npnew)*autoaa,ep(1:npnew))
    deallocate (xyzp,atp)

    call tim%start(17,'Ensemble refinement')
    call crest_refine(env,trim(atmp),trim(atmp))
    call tim%stop(17)

!>--- sorting
    write (stdout,'(a)') '> Sorting structures by energy ...'
    env%ewin = env%protb%ewin*3.0_wp
    call newcregen(env,7,trim(atmp))
    call rename(trim(atmp)//'.sorted',trim(atmp))
  end if

!========================================================================================!
!> Remove doubly generated tautomers
  call protonation_prep_canonical(env,mol,trim(atmp))

!========================================================================================!
!>--- move final ensemble to deprotonated.xyz
  call rename(trim(atmp),'deprotonated.xyz')
  write (stdout,'(a)') '> Sorted file was renamed to deprotonated.xyz'

  write (stdout,'(a)') '> All other temporary deprotonate_*.xyz files are removed by default.'
  if (.not.env%keepmodef) then
    call rmrf('deprotonate_*.xyz')
  end if
!========================================================================================!
  return
end subroutine crest_new_deprotonate

!========================================================================================!
!========================================================================================!

subroutine deprotonation_sort(mol,nhyd)
!*******************************************************
!* Sort "mol" so that all heavy atoms are first
!* and all hydrogens follow. Also returns the number of
!* hydrogen atoms and writes a file with the original
!* atom order
!*******************************************************
  use crest_parameters
  use strucrd
  implicit none
  type(coord),intent(inout) :: mol
  integer,intent(out) :: nhyd
  integer :: i,j,k,l,io,ich
  integer,allocatable :: atnew(:)
  real(wp),allocatable :: xyznew(:,:)

  nhyd = 0
  open (newunit=ich,file='original.atomorder')
  write (ich,'(a15,a15)') '<atom was at>','<atom is now>'

  allocate (atnew(mol%nat),source=0)
  allocate (xyznew(3,mol%nat),source=0.0_wp)

  k = 0
  !>--- first heavy atoms
  do i = 1,mol%nat
    if (mol%at(i) .ne. 1) then
      k = k+1
      atnew(k) = mol%at(i)
      xyznew(1:3,k) = mol%xyz(1:3,i)
      write (ich,'(i15,i15)') i,k
    end if
  end do
  !>--- then hydrogens
  do i = 1,mol%nat
    if (mol%at(i) .eq. 1) then
      k = k+1
      nhyd = nhyd+1
      atnew(k) = mol%at(i)
      xyznew(1:3,k) = mol%xyz(1:3,i)
      write (ich,'(i15,i15)') i,k
    end if
  end do

  call move_alloc(atnew,mol%at)
  call move_alloc(xyznew,mol%xyz)

  close (ich)
end subroutine deprotonation_sort

!========================================================================================!

subroutine deprotonation_candidates(env,mol,natp,np,at,xyz,npnew)
!********************************************************
!* generate deprotonation/ionization candidate structures
!* The outputs are at and xyz, the latter being in Bohr
!********************************************************
  use crest_data
  use crest_parameters
  use strucrd,only:coord
  use utilities,only:binomial,get_combinations
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env
  type(coord),intent(in) :: mol
  integer,intent(in) :: natp  !> total number of atoms
  integer,intent(in) :: np    !> number of H⁺
  !> OUTPUT
  integer,intent(out)  :: at(natp)
  real(wp),intent(out) :: xyz(3,natp,npnew)
  integer,intent(inout) :: npnew
  !> LOCAL
  integer :: i,j,jj,k,l,ii,c,nc,kk
  integer :: ati,ichrg,ctype,nhvy
  integer,allocatable :: combi(:,:),tmp(:)

  if (natp .ne. mol%nat-env%protb%amount) then
    write (stdout,'(a)') 'WARNING: Inconsistent number of atoms in deprotonation routine'
    call creststop(status_args)
  end if

  if (env%protb%swelem) then
!>--- User-defined monoatomic ion
    ichrg = env%protb%swchrg*env%protb%amount
    ati = env%protb%swat
  else
!>--- DEFAULT: H⁺
    ichrg = 1*env%protb%amount
    ati = 1
  end if
  write (stdout,'(a,i0)') '> Decreasing the molecular charge by ',ichrg
  call env%calc%decrease_charge(ichrg)
  env%chrg = env%chrg-ichrg

!>--- check combinations
  k = env%protb%amount
  nc = binomial(np,k)
  allocate (tmp(k),combi(k,nc),source=0)
  c = 0
  call get_combinations(np,k,nc,c,combi,tmp,0)

!>--- Check other conditions?
  ! TODO

!>--- Populate heavy atoms
  nhvy = mol%nat-np
  at(1:nhvy) = mol%at(1:nhvy)

  npnew = 0
  ii = 0
  COMBILOOP: do i = 1,nc
    ii = ii+1
    kk = 0
    ADDLOOP1: do j = 1,mol%nat
      do l = 1,k
        jj = combi(l,i)+nhvy
        if (j == jj) cycle ADDLOOP1
      end do
      kk = kk+1
      xyz(1:3,kk,ii) = mol%xyz(1:3,j)
      at(kk) = mol%at(j)
    end do ADDLOOP1

  end do COMBILOOP
  npnew = ii

  deallocate (combi,tmp)
end subroutine deprotonation_candidates

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> New implementation of tautomer routines
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
subroutine crest_new_tautomerize(env,tim)
!********************************************************************
!* Standalone runtype for generating molecule tautomers
!*
!* Input/Output:
!*  env  -  crest's systemdata object
!*  tim  -  timer object
!********************************************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
  use parallel_interface
  use cregen_interface
  use iomod
  use utilities,only:binomial
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,T,Tn,np,npadd,npremove
  logical :: pr,wr
  character(len=80) :: atmp
!========================================================================================!
  type(calcdata) :: calc
  real(wp) :: accuracy,etemp

  real(wp) :: energy,dip
  real(wp),allocatable :: grad(:,:)
  type(calcdata),allocatable :: tmpcalc
  type(calcdata),allocatable :: tmpcalc_ff
  type(calculation_settings) :: tmpset
  integer :: natp,pstep,npnew
  integer :: tautiter,structiter
  integer,allocatable :: atp(:)
  real(wp),allocatable :: xyzp(:,:,:)
  real(wp),allocatable :: ep(:)
  logical,allocatable :: atlist(:)
  real(wp),allocatable :: protxyz(:,:)
  character(len=*),parameter :: basename = 'tautomerize_'
!========================================================================================!
  write (stdout,*)
  !call system('figlet singlepoint')
  write (stdout,*) " _              _                            _         "
  write (stdout,*) "| |_ __ _ _   _| |_ ___  _ __ ___   ___ _ __(_)_______ "
  write (stdout,*) "| __/ _` | | | | __/ _ \| '_ ` _ \ / _ \ '__| |_  / _ \"
  write (stdout,*) "| || (_| | |_| | || (_) | | | | | |  __/ |  | |/ /  __/"
  write (stdout,*) " \__\__,_|\__,_|\__\___/|_| |_| |_|\___|_|  |_/___\___|"
  write (stdout,*) "                                                       "
  write (stdout,*) "-------------------------------------------------------"
  write (stdout,*) "         automated tautomer screening script           "
  write (stdout,*) "          revised version (c) P.Pracht 2024            "
  write (stdout,*) 'Cite as:'
  write (stdout,*) '  P.Pracht, R.Wilcken, A.Udvarhelyi, S.Rodde, S.Grimme'
  write (stdout,*) '  JCAMD, 2018, 32, 1139-1149.'
  write (stdout,*)

!========================================================================================!
  call new_ompautoset(env,'max',0,T,Tn)
  call ompprint_intern()
!========================================================================================!
  call env%ref%to(mol)
  write (stdout,*)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)
!========================================================================================!
!>--- The first step is to sort the input structure so that all heavy atoms are first and
!>--- all hydrogen atoms follow after that. This is because some of the (parallel) routines
!>--- rely on identical atom order.

  write (stdout,'(a)') repeat('-',80)
  write (stdout,'(a)')
  write (stdout,'(a)') '> Re-sorting input structure ... '
  call deprotonation_sort(mol,npremove)
  call mol%append(stdout)
  write (stdout,*)
!>--- check
  if (npremove == 0) then
    write (stdout,*)
    write (stdout,*) 'WARNING: No deprotonation sites found, therefore no tautomer search possible!'
    write (stdout,*)
    return
  end if

!========================================================================================!
!========================================================================================!
!>--- Tautomerization can be run iteratively to do further permutations
!  TAUTLOOP : do tautiter=1,env%protb%iter
!========================================================================================!
!========================================================================================!
!>--- Then, we need to perfom a LMO calculation to identify suitable protonation sites
!>--- in this version of the program this step is always done with GFN0-xTB.

  call tim%start(14,'LMO center calculation')
  write (stdout,'(a)') repeat('-',80)
  write (stdout,'(a)')
  write (stdout,'(a)',advance='no') '> Setting up GFN0-xTB for LMO center calculation ... '
  flush (stdout)

  allocate (grad(3,mol%nat),source=0.0_wp)
!>--- GFN0 job adapted from global settings
  allocate (tmpcalc)
  call env2calc(env,tmpcalc,mol)
  tmpcalc%calcs(1)%id = jobtype%gfn0
  tmpcalc%calcs(1)%rdwbo = .true.
  tmpcalc%calcs(1)%getlmocent = .true.
  tmpcalc%calcs(1)%chrg = env%chrg
  tmpcalc%ncalculations = 1
  write (stdout,'(a)') 'done.'
!>--- and then start it
  write (stdout,'(a)',advance='yes') '> Performing singlepoint calculation ... '
  flush (stdout)
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<!
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<!
  call engrad(mol,tmpcalc,energy,grad,io)
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<!
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<!
  call tim%stop(14)
  write (stdout,'(a)') '> done.'
  write (atmp,'(a)') '> Total wall time for calculations'
  call tim%write_timing(stdout,14,trim(atmp),.true.)
  write (stdout,'(a)') repeat('-',80)
  write (stdout,*)
  if (io /= 0) then
    write (stdout,*)
    write (stdout,*) 'WARNING: Calculation exited with error!'
  end if
!>--- check LMO center
  npadd = tmpcalc%calcs(1)%nprot
  if (npadd > 0) then
    write (stdout,'(a,i0,a)') '> ',np,' π- or LP-centers identified as protonation candidates:'
    call move_alloc(tmpcalc%calcs(1)%protxyz,protxyz)
    write (stdout,'(1x,a5,1x,a,5x,a)') 'LMO','type','center(xyz/Ang)'
    do i = 1,npadd
      select case (nint(protxyz(4,i)))
      case (1)
        write (stdout,'(1x,i5,1x,a,3F12.5)') i,'LP   ',protxyz(1:3,i)*autoaa
      case (2)
        write (stdout,'(1x,i5,1x,a,3F12.5)') i,'π    ',protxyz(1:3,i)*autoaa
      case (3)
        write (stdout,'(1x,i5,1x,a,3F12.5)') i,'del.π',protxyz(1:3,i)*autoaa
      case default
        write (stdout,'(1x,i5,1x,a,3F12.5)') i,'???  ',protxyz(1:3,i)*autoaa
      end select
    end do
  else
    write (stdout,*)
    write (stdout,*) 'WARNING: No suitable protonation sites found!'
    write (stdout,*) '         Confirm whether you expect π- or LP-centers for your molecule!'
    write (stdout,*)
    env%iostatus_meta = status_failed
    return
  end if
  deallocate (tmpcalc)
  deallocate (grad)

!========================================================================================!
!>--- If we reached this point, we have candidate positions for our protonation!
  pstep = 0
  write (stdout,'(a)',advance='yes') '> Generating candidate structures ... '
  flush (stdout)
  !> no change in atom numbers for tautomers!
  natp = mol%nat
  !> however, the number of permutations can get gigantic!
  !npnew = binomial(npadd,env%protb%amount)*binomial(npremove,env%protb%amount)
  !> for now, refer to a single-tautomer-permutation (moving one H)
  npnew = npadd*npremove
  write (stdout,'(a,i0,a)') '> Up to ',npnew,' new structures'
  allocate (atp(natp),source=0)
  allocate (xyzp(3,natp,npnew),ep(npnew),source=0.0_wp)

  !> generate initial structures
  call tautomer_candidates(env,mol,natp,npadd,npremove,protxyz,atp,xyzp,npnew)

  if (npnew < 1) then
    write (stdout,*)
    write (stdout,'(a)') '> WARNING: No remaining tautomer candidates after applying user defined conditions!'
    write (stdout,'(a)') '>          Modify the search criteria and check your input structure for sanity.'
    return
  end if

  write (atmp,'(a,i0,a)') basename,pstep,'.xyz'
  write (stdout,'(a,a,a,i0,a)') '> Write ',trim(atmp),' with ',npnew,' candidates ... '

  call wrensemble(basename//'0.xyz',natp,npnew,atp,xyzp(:,:,1:npnew)*autoaa,ep(1:npnew))

  write (stdout,'(a)') '> done.'
  write (stdout,*)

!========================================================================================!
!>--- Optimize candidates, optional FF pre-step
  if (env%protb%ffopt) then
    call smallhead('Tautomer Ensemble FF Pre-Optimization')
    !>--- set up a temporary calculation object
    write (stdout,'(a)') '> Removal/addition of hydrogen atoms may lead to initial'
    write (stdout,'(a)') '  extremely high-energy candidates which can undergo unwanted'
    write (stdout,'(a)') '  chemical changes in optimizations. Classical force-fields'
    write (stdout,'(a)') '  with defined bond-topology can help circumvent this issue.'
    write (stdout,'(a)') '> Setting up force-field structure pre-optimization ...'
    allocate (tmpcalc_ff)
    tmpcalc_ff%optnewinit = .true.
    env%calc%optnewinit = .true.
    tmpcalc_ff%optnewinit = .true.
    call tmpset%create('gfnff')
    tmpset%chrg = env%chrg
    call tmpcalc_ff%add(tmpset)
    tmpcalc_ff%maxcycle = 10000
    tmpcalc_ff%anopt=.true.
    call tmpcalc_ff%info(stdout)

    !>--- Run optimizations
    call tim%start(15,'Ensemble optimization (FF)')
    call print_opt_data(tmpcalc_ff,stdout)
    write (stdout,'(a,i0,a)') '> ',npnew,' structures to optimize ...'
    call crest_oloop(env,natp,npnew,atp,xyzp(:,:,1:npnew),ep(1:npnew),.false.,tmpcalc_ff)
    call tim%stop(15)

    deallocate (tmpcalc_ff)

    pstep = pstep+1
    write (atmp,'(a,i0,a)') basename,pstep,'.xyz'
    write (stdout,'(a,a,a)') '> Write ',trim(atmp),' with optimized structures ... '
    call wrensemble(trim(atmp),natp,npnew,atp,xyzp(:,:,1:npnew)*autoaa,ep(1:npnew))
    deallocate (xyzp,atp) !> clear this space to re-use it

    !>--- sorting
    write (stdout,'(a)') '> Sorting structures by energy to remove failed opts. ...'
    env%ewin = 1000.0_wp
    call newcregen(env,7,trim(atmp))
    call rename(trim(atmp)//'.sorted',trim(atmp))
    write (stdout,'(a)') '> Sorted file was renamed to '//trim(atmp)
    write (stdout,'(a)') '> WARNING: These are force-field energies which are '
    write (stdout,'(a)') '           NOT(!) accurate for bond formation and breaking!'
    write (stdout,*)

    !> Remove doubly generated tautomers
    call protonation_prep_canonical(env,mol,trim(atmp))

    !>--- re-read sorted ensemble
    call rdensemble(trim(atmp),natp,npnew,atp,xyzp)
    xyzp = xyzp*aatoau !> don't forget to restore BOHR
  end if

!========================================================================================!
!>--- H-position-only optimization (only makes sense after FF preoptimization)
  if (env%protb%hnewopt.and.env%protb%ffopt) then
    call smallhead('Tautomer Ensemble Frozen-Atom Optimization')
    !>--- create temporary calculation from our intended level of theory
    write (stdout,'(a)') '> Setting up frozen structure optimization ...'
    allocate (tmpcalc)
    call env2calc(env,tmpcalc,mol)
    !>--- freeze all atoms, EXCEPT the new one
    allocate (atlist(natp),source=.true.)
    atlist(natp) = .false. !> the new one is always last
    do i = 1,natp
      if (atp(i) == 1) then
        atlist(i) = .false. !> additionally un-freeze all H's (this seems to be beneficial)
      end if
    end do
    tmpcalc%nfreeze = count(atlist)
    tmpcalc%optnewinit = .true.
    write (stdout,'(a,i0,a)') '> ',tmpcalc%nfreeze,' frozen atoms. All H non-frozen.'
    call move_alloc(atlist,tmpcalc%freezelist)
    call tmpcalc%info(stdout)

    !>--- run opt
    call tim%start(16,'Ensemble optimization (frozen)')
    write (stdout,'(a,i0,a)') '> ',npnew,' structures to optimize ...'
    call crest_oloop(env,natp,npnew,atp,xyzp(:,:,1:npnew),ep(1:npnew),.false.,tmpcalc)
    call tim%stop(16)

    pstep = pstep+1
    write (atmp,'(a,i0,a)') basename,pstep,'.xyz'
    write (stdout,'(a,a,a)') '> Write ',trim(atmp),' with optimized structures ... '
    call wrensemble(trim(atmp),natp,npnew,atp,xyzp(:,:,1:npnew)*autoaa,ep(1:npnew))
    deallocate (xyzp,atp) !> clear this space to re-use it
    deallocate (tmpcalc)

!    call tim%start(17,'Ensemble refinement')
!    call crest_refine(env,trim(atmp),trim(atmp))
!    call tim%stop(17)

    !>--- sorting
    write (stdout,'(a)') '> Sorting structures by energy to remove failed opts. ...'
    env%ewin = env%protb%ewin*3.0_wp  !> large energy threshold
    call newcregen(env,7,trim(atmp))
    call rename(trim(atmp)//'.sorted',trim(atmp))
    write (stdout,'(a)') '> Sorted file was renamed to '//trim(atmp)
    write (stdout,*)

    !>--- re-read sorted ensemble
    call rdensemble(trim(atmp),natp,npnew,atp,xyzp)
    xyzp = xyzp*aatoau !> don't forget to restore BOHR

  end if

!========================================================================================!
!>--- Optimize with global settings
  if (env%protb%finalopt.and.env%protb%ffopt) then
    call smallhead('Final Tautomer Ensemble Optimization')
    allocate (tmpcalc)
    call env2calc(env,tmpcalc,mol)
    call tmpcalc%info(stdout)
    tmpcalc%optnewinit = .true.
    call tim%start(20,'Ensemble optimization')
    call print_opt_data(env%calc,stdout)
    write (stdout,'(a,i0,a)') '> ',npnew,' structures to optimize ...'
    call crest_oloop(env,natp,npnew,atp,xyzp(:,:,1:npnew),ep,.false.,tmpcalc)
    call tim%stop(20)

    pstep = pstep+1
    write (atmp,'(a,i0,a)') basename,pstep,'.xyz'
    write (stdout,'(a,a,a)') '> Write ',trim(atmp),' with optimized structures ... '
    call wrensemble(trim(atmp),natp,npnew,atp,xyzp(:,:,1:npnew)*autoaa,ep(1:npnew))
    deallocate (xyzp,atp)

    call tim%start(17,'Ensemble refinement')
    call crest_refine(env,trim(atmp),trim(atmp))
    call tim%stop(17)

!>--- sorting
    write (stdout,'(a)') '> Sorting structures by energy ...'
    env%ewin = env%protb%ewin*2.0_wp
    call newcregen(env,7,trim(atmp))
    call rename(trim(atmp)//'.sorted',trim(atmp))
  end if

!========================================================================================!
!> Remove doubly generated tautomers
  call protonation_prep_canonical(env,mol,trim(atmp))

!========================================================================================!
!========================================================================================!
!  enddo TAUTLOOP
!========================================================================================!
!========================================================================================!
!>--- move final ensemble to deprotonated.xyz
  call rename(trim(atmp),'tautomers.xyz')
  write (stdout,'(a)') '> Sorted file was renamed to tautomers.xyz'

  write (stdout,'(a)') '> All other temporary tautomerize_*.xyz files are removed by default.'
  if (.not.env%keepmodef) then
    call rmrf('tautomerize_*.xyz')
  end if

  return
end subroutine crest_new_tautomerize

!=========================================================================================!

subroutine tautomer_candidates(env,mol,natp,npadd,npremove,protxyz,at,xyz,npnew)
!********************************************************
!* generate tautomer candidate structures
!* The outputs are at and xyz, the latter being in Bohr
!********************************************************
  use crest_data
  use crest_parameters
  use strucrd,only:coord
  use utilities,only:binomial,get_combinations
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env
  type(coord),intent(in) :: mol
  integer,intent(in) :: natp
  integer,intent(in) :: npadd,npremove
  real(wp),intent(in) :: protxyz(4,npadd)
  !> OUTPUT
  integer,intent(out)  :: at(natp)
  real(wp),intent(out) :: xyz(3,natp,npnew)
  integer,intent(inout) :: npnew
  !> LOCAL
  integer :: i,j,jj,k,l,ii,c,nc,kk,nhvy
  integer :: ati,ichrg,ctype
  integer,allocatable :: combiadd(:,:),tmpadd(:),combiremove(:,:),tmpremove(:)

  if (natp .ne. mol%nat) then
    write (stdout,'(a)') 'WARNING: Inconsistent number of atoms in protonation routine'
    call creststop(status_args)
  end if

  ati = 1  !> always refer to Hydrogen for tautomers

!>--- Check if we have some other conditions
  if (any(.not.env%protb%active_lmo(:))) then
    write (stdout,'(a)',advance='no') '> User-defined: IGNORING '
    if (.not.env%protb%active_lmo(1)) write (stdout,'(a)',advance='no') 'π '
    if (.not.env%protb%active_lmo(2)) write (stdout,'(a)',advance='no') 'LP '
    if (.not.env%protb%active_lmo(3)) write (stdout,'(a)',advance='no') 'deloc.π '
    write (stdout,'(a)') 'LMOs ...'
  end if

  nc = npadd*npremove
  nhvy = mol%nat-npremove

!>--- Populate
  npnew = 0
  ii = 0
  COMBILOOP: do i = 1,npremove
    l = i+nhvy
    ADDLOOP1: do jj = 1,npadd
      ctype = nint(protxyz(4,jj))
!>--- Enforce further constraints, conditions, etc.
      if (.not.env%protb%active_lmo(ctype)) cycle COMBILOOP

!>--- passed checks means we can add this config
      ii = ii+1 !> counter of actually created structures
      kk = 0
      do j = 1,nhvy  !> first all heavy atoms (xyz should be sorted!)
        kk=kk+1
        xyz(1:3,kk,ii) = mol%xyz(1:3,j)
        at(kk) = mol%at(j)
      end do
      ADDLOOP2: do j = nhvy+1,natp !> then all "old" H except the moving one
        if (j == l) cycle ADDLOOP2 !> cycle the moving H
        kk = kk+1
        xyz(1:3,kk,ii) = mol%xyz(1:3,j)
        at(kk) = ati
      end do ADDLOOP2
      !> finally, add the new H from the protonation center (it is always the last H)
      xyz(1:3,natp,ii) = protxyz(1:3,jj)
      at(natp) = ati
    end do ADDLOOP1
  end do COMBILOOP
  npnew = ii

end subroutine tautomer_candidates

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!

