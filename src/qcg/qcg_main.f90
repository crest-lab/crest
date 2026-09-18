!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Christoph Plett, Sebastian Spicher, Philipp Pracht
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

!> This file contains routines related to QCG and microsolvation

subroutine crest_solvtool(env,tim)
!***********************************************
!* Main driver for all QCG runtypes
!***********************************************
  use crest_parameters,only:wp,autokcal
  use qcg_printouts
  use crest_data
  use iomod
  use qcg_coord_type
  use strucrd
  implicit none

  type(systemdata):: env    ! MAIN STORAGE OS SYSTEM DATA
  type(timer):: tim
  !> Information about solvent, solute and cluster
  type(coord_qcg) :: solute,solvent,cluster,cluster_backup
  type(ensemble) :: full_ensemble,solvent_ensemble

  integer :: progress,io
  integer :: chrg_tmp,uhf_tmp
  character(len=512) :: thispath

!--- Molecule settings
  solute%nmol = 1
  solvent%nmol = 1
  cluster%nmol = 1

  progress = 0
  call getcwd(thispath)

  !>-----------------------------------
  call qcg_head()
  !>-----------------------------------

!> Check, if xtb is present
  call checkprog_silent(env%ProgName,.true.,iostat=io)
  if (io /= 0) error stop 'No xtb found'

!> Check, if xtbiff is present (if it is required)
  if (env%use_xtbiff) then
    call xtbiff_print_deprecated()
  else
    write (stdout,*)
    write (stdout,*) ' This program uses the the aISS algorithm as implemnted in xtb.'
    write (stdout,*) '      The aISS method requires xtb version 6.6.0 or newer.'
    write (stdout,*) '           Tested with xtb version 6.7.1 (902b313)'
    !write (stdout,*) '  xTB-IFF can still be used with the --xtbiff flag.'
    write (stdout,*)
  end if

!------------------------------------------------------------------------------
!   Setup
!------------------------------------------------------------------------------
  call write_qcg_setup(env) !Just an outprint of setup
  call read_qcg_input(env,solute,solvent) !Reading mol. data and determining r,V,A
  call qcg_setup(env,solute,solvent)
  call qcg_restart(env,progress,solute,solvent,cluster,full_ensemble,&
         & solvent_ensemble,cluster_backup)

!-----------------------------------------------------------------------------
!   Grow
!-----------------------------------------------------------------------------
  if (progress .le. env%qcg_runtype.and.progress .eq. 0) then
    cluster = solute
    call qcg_grow(env,solute,solvent,cluster,tim)
    if (.not.env%cff) then
      allocate (cluster_backup%at(cluster%nat))
      allocate (cluster_backup%xyz(3,cluster%nat))
      cluster_backup = cluster
    end if
    progress = progress+1
    call chdirdbug(thispath)
  end if

!>-- monomer references at the energy level of the ensemble, also needed on restart
  if (env%qcg_runtype .ge. 2.and.progress .le. env%qcg_runtype) then
    call qcg_update_reference(env,solute,.true.)
    call qcg_update_reference(env,solvent,.false.)
  end if

!------------------------------------------------------------------------------
!   Ensemble search
!------------------------------------------------------------------------------
  if (progress .le. env%qcg_runtype.and.progress .eq. 1) then
    call print_qcg_ensemble()
    call qcg_ensemble(env,solute,solvent,cluster,full_ensemble,tim,'ensemble')
    progress = progress+1
    call chdirdbug(thispath)
  end if

!------------------------------------------------------------------------------
!   Solvent cluster generation
!------------------------------------------------------------------------------
  if (progress .le. env%qcg_runtype.and.progress .eq. 2) then !esolv
    call pr_eval_solvent()
    if (env%cff) then !CFF
      call qcg_cff(env,solute,solvent,cluster,full_ensemble,&
             & solvent_ensemble,tim)
    else !Normal ensemble generation
      call print_qcg_ensemble()
      call cluster%deallocate
      allocate (cluster%at(cluster_backup%nat))
      allocate (cluster%xyz(3,cluster_backup%nat))
      cluster = cluster_backup
      deallocate (cluster_backup%at)
      deallocate (cluster_backup%xyz)
      env%solv_md = .true.
      !>-- pure solvent cluster, neutral closed-shell
      chrg_tmp = env%chrg
      uhf_tmp = env%uhf
      env%chrg = 0
      env%uhf = 0
      call qcg_ensemble(env,solute,solvent,cluster,solvent_ensemble,&
             & tim,'solvent_ensemble')
      env%chrg = chrg_tmp
      env%uhf = uhf_tmp
    end if
    call qcg_check_reference(env,solute)
    call pr_qcg_esolv()
    write (stdout,'(2x,"|",7x,F10.2," kcal/mol ",12x,"|")') &
           &   full_ensemble%g-solvent_ensemble%g-(solute%energy*autokcal)
    write (stdout,'(2x,''========================================='')')
    call chdirdbug(thispath)
    progress = progress+1
  end if

!------------------------------------------------------------------------------
!   Frequency computation and evaluation
!------------------------------------------------------------------------------
  if (progress .le. env%qcg_runtype.and.progress .eq. 3) then !gsolv
    call qcg_freq(env,tim,solute,solvent,full_ensemble,solvent_ensemble)
    call qcg_eval(env,solute,full_ensemble,solvent_ensemble)
    progress = progress+1
  end if

!------------------------------------------------------------------------------
!   Cleanup and deallocation
!------------------------------------------------------------------------------
  if (env%scratchdir .ne. 'qcg_tmp') call qcg_cleanup(env)
  if (.not.env%keepModef) call rmrf('qcg_tmp')
  call solute%deallocate
  call solvent%deallocate
  call cluster%deallocate
  call full_ensemble%deallocate
  call solvent_ensemble%deallocate
  return
end subroutine crest_solvtool

!==============================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!==============================================================================!

subroutine qcg_setup(env,solu,solv)
  use iso_fortran_env,wp => real64
  use crest_data
  use iomod
  use qcg_coord_type
  use strucrd
  use axis_module
  implicit none

  type(systemdata):: env
  type(coord_qcg) :: solv,solu

  integer :: io,f,r
  integer :: num_O,num_H,i
  integer :: chrg_tmp,uhf_tmp
  character(len=*),parameter :: outfmt = '(1x,1x,a,1x,f14.7,a,1x)'
  logical :: e_there,tmp,used_tmp,gbsa_tmp
  character(len=512) :: thispath,tmp_grow
  character(len=40)  :: solv_tmp
  character(len=80)  :: atmp
  character(len=20)  :: gfnver_tmp

  call getcwd(thispath)

  ! Remove scratch dir, if present
  inquire (file='./qcg_tmp/solute_properties/solute',exist=tmp)
  if (tmp) call rmrf('qcg_tmp') !User given scratch dir will be removed anyway after run

  ! Make scratch directories
  if (env%scratchdir .eq. '') then !check if scratch was not set
    env%scratchdir = 'qcg_tmp'
    io = makedir('qcg_tmp')
  end if
  if (env%fixfile /= 'none selected') then
    call copysub(env%fixfile,env%scratchdir)
  end if
  call chdirdbug(env%scratchdir)

  f = makedir('solute_properties')
  if (env%fixfile /= 'none selected') then
    call copysub(env%fixfile,env%scratchdir)
  end if
  r = makedir('solvent_properties')

  if (.not.env%nopreopt) then
    write (stdout,*)
    write (stdout,'(2x,''========================================='')')
    write (stdout,'(2x,''|            Preoptimization            |'')')
    write (stdout,'(2x,''========================================='')')
  end if

  solv_tmp = env%solv
  gbsa_tmp = env%gbsa
  env%solv = ''
  env%gbsa = .false.

!---- Properties solute
  call chdirdbug('solute_properties')
  call env%wrtCHRG('') !Write three lines in QCG mode, but xtb anyway only reads first one

!---- Geometry preoptimization solute
  if (env%final_gfn2_opt) then !If GFN2 final opt, solute also GFN2 optimized
    gfnver_tmp = env%gfnver
    env%gfnver = '--gfn2'
  end if

  if ((.not.env%nopreopt).and.(solu%nat /= 1)) then
    call xtb_opt_qcg(env,solu,.true.)
  end if

!--- Axistrf
  call axistrf(solu%nat,solu%nat,solu%at,solu%xyz)
  call solu%write('solute')

!---- SP-Computation solute
  call xtb_sp_qcg(env,'solute',e_there,solu%energy)
  solu%elevel = env%gfnver

  if (.not.e_there) then
    write (stdout,*) 'Total Energy of solute not found'
  else
    write (stdout,outfmt) 'Total Energy of solute: ',solu%energy,' Eh'
  end if

  if (env%final_gfn2_opt) then !If GFN2 final opt, solute also GFN2 optimized
    env%gfnver = gfnver_tmp
  end if

  call chdirdbug(thispath)

! No constraints for solvent possible
  used_tmp = env%cts%used
  env%cts%used = .false.
! The solvent is neutral closed-shell, charge and uhf belong to the solute
  chrg_tmp = env%chrg
  uhf_tmp = env%uhf
  env%chrg = 0
  env%uhf = 0

!---- Properties solvent
  call chdirdbug(env%scratchdir)
  call chdirdbug('solvent_properties')
  !No charges for solvent written. This is currently not possible

!---- Geometry preoptimization solvent
  if ((.not.env%nopreopt).and.(solv%nat /= 1)) then
    call xtb_opt_qcg(env,solv,.false.)
  end if
  call solv%write('solvent')

!---- SP-Computation solvent
  call xtb_sp_qcg(env,'solvent',e_there,solv%energy)
  solv%elevel = env%gfnver

  if (.not.e_there) then
    write (stdout,'(1x,a)') 'Total Energy of solvent not found'
  else
    write (stdout,outfmt) 'Total energy of solvent:',solv%energy,' Eh'
  end if

  call chdirdbug(thispath)

!---- Overwriting solute and solvent in original folder
  call solu%write('solute')
  call solv%write('solvent')

  num_O = 0
  num_H = 0
!--- Check, if water is solvent
  if (solv%nat .eq. 3) then
    do i = 1,solv%nat
      if (solv%at(i) .eq. 8) num_O = num_O+1
      if (solv%at(i) .eq. 1) num_H = num_H+1
    end do
  end if
  if (num_O .eq. 1.AND.num_H .eq. 2) then
    env%water = .true.
    if (.not.env%noconst) env%constrain_solu = .true.
  end if

  env%solv = solv_tmp
  env%gbsa = gbsa_tmp
  env%cts%used = used_tmp
  env%chrg = chrg_tmp
  env%uhf = uhf_tmp

end subroutine qcg_setup

!==============================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!==============================================================================!

subroutine qcg_update_reference(env,mol,is_solute)
  !*****************************************************************
  !* Re-evaluate the gas-phase monomer reference energy at the QCG
  !* energy level (see qcg_energy_level), if it was obtained at a
  !* different one. The growth needs mol%energy at env%gfnver, while
  !* Esolv/Gsolv and CFF need it at the level of the ensemble SPs.
  !*
  !* mol       : solute or solvent; only %energy and %elevel change,
  !*             the geometry is kept
  !* is_solute : .true. for the solute, .false. for the solvent
  !*****************************************************************
  use crest_parameters
  use crest_data
  use iomod
  use qcg_coord_type
  use qcg_utils
  implicit none

  type(systemdata) :: env
  type(coord_qcg),intent(inout) :: mol
  logical,intent(in) :: is_solute

  type(coord_qcg) :: ref
  integer :: io,chrg_tmp,uhf_tmp
  logical :: e_there,used_tmp,gbsa_tmp
  real(wp) :: eref
  character(len=512) :: thispath
  character(len=40)  :: solv_tmp
  character(len=20)  :: gfnver_tmp,elevel
  character(len=:),allocatable :: molname
  character(len=*),parameter :: outfmt = '(1x,1x,a,1x,f14.7,a,1x)'

  elevel = qcg_energy_level(env)
  if (trim(mol%elevel) == trim(elevel)) return

  if (is_solute) then
    molname = 'solute'
  else
    molname = 'solvent'
  end if
  call getcwd(thispath)

!>-- gas phase at the energy level, same settings as in qcg_setup
  solv_tmp = env%solv
  gbsa_tmp = env%gbsa
  used_tmp = env%cts%used
  gfnver_tmp = env%gfnver
  env%solv = ''
  env%gbsa = .false.
  env%gfnver = elevel
  chrg_tmp = env%chrg
  uhf_tmp = env%uhf
  if (.not.is_solute) then
    env%cts%used = .false.
    env%chrg = 0
    env%uhf = 0
  end if

!>-- separate directory, the files of qcg_setup are still needed (CFF)
  call chdirdbug(env%scratchdir)
  call chdirdbug(molname//'_properties')
  io = makedir('reference')
  call chdirdbug('reference')
  if (is_solute) call env%wrtCHRG('')

!>-- work on a copy, mol%xyz is used later for the monomer frequencies
  ref = mol
  if ((.not.env%nopreopt).and.(ref%nat /= 1)) then
    call xtb_opt_qcg(env,ref,is_solute)
  end if
  call ref%write(molname)
  call xtb_sp_qcg(env,molname,e_there,eref)

  env%solv = solv_tmp
  env%gbsa = gbsa_tmp
  env%cts%used = used_tmp
  env%gfnver = gfnver_tmp
  env%chrg = chrg_tmp
  env%uhf = uhf_tmp
  call chdirdbug(thispath)

  if (.not.e_there) then
    write (stdout,'(1x,a,a,a,a)') 'ERROR: reference energy of the ',molname, &
      & ' could not be computed at the level ',trim(elevel)
    call creststop(status_safety)
  end if

  write (stdout,*)
  write (stdout,'(2x,a,a,a,a,a,a)') 'Reference energy of the ',molname, &
    & ' re-evaluated: ',trim(mol%elevel),' -> ',trim(elevel)
  write (stdout,outfmt) 'Total Energy of '//molname//':',eref,' Eh'

  mol%energy = eref
  mol%elevel = elevel
  call ref%deallocate()

end subroutine qcg_update_reference

!==============================================================================!

subroutine qcg_check_reference(env,mol)
  !***********************************************
  !* Stop if the monomer reference energy was not
  !* evaluated at the QCG energy level.
  !***********************************************
  use crest_parameters
  use crest_data
  use qcg_coord_type
  use qcg_utils
  implicit none
  type(systemdata) :: env
  type(coord_qcg),intent(in) :: mol
  character(len=20) :: elevel

  elevel = qcg_energy_level(env)
  if (trim(mol%elevel) /= trim(elevel)) then
    write (stdout,'(1x,a,a,a,a)') 'ERROR: solute reference energy evaluated at ', &
      & trim(mol%elevel),' but the cluster energies at ',trim(elevel)
    call creststop(status_safety)
  end if
end subroutine qcg_check_reference

!==============================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!==============================================================================!

subroutine read_qcg_input(env,solu,solv)
  use crest_parameters
  use crest_data
  use iomod
  use qcg_coord_type
  use strucrd
  use atmasses
  use qcg_utils
  implicit none

  type(systemdata)               :: env
  type(coord_qcg),intent(inout) :: solu,solv
  logical                        :: pr
  real(wp),parameter            :: third = 1.0d0/3.0d0
  integer                        :: i
  real(wp)                       :: r_solu,r_solv

  pr = .true.

!--- Read in solu and solv coordinates and make solute and solvent file in WD
  call inputcoords_qcg(env,solu,solv)

!--- CMA-Trafo
  call cma_shifting(solu,solv)

!--- Setting solute charge and uhf to input
  solu%chrg = env%chrg
  solu%uhf = env%uhf

!--- Getting r, V, A
  write (stdout,*)
  write (stdout,*) 'Solute geometry'
  call get_sphere(.true.,solu,.true.) !r,V,A of solute
  write (stdout,*) 'Solvent geometry'
  call get_sphere(.true.,solv,.true.) !r,V,A of solvent

  r_solu = solu%vtot**third
  r_solv = solv%vtot**third
  write (stdout,*)
  write (stdout,'(2x,''radius of solute    : '',f8.2)') r_solu
  write (stdout,'(2x,''radius of solvent   : '',f8.2)') r_solv

!--- Determine masses (for later density computation)
  do i = 1,solu%nat
    solu%mass = solu%mass+ams(solu%at(i))
  end do
  do i = 1,solv%nat
    solv%mass = solv%mass+ams(solv%at(i))
  end do
  solu%mass = solu%mass*amutokg
  solv%mass = solv%mass*amutokg

!--- If directed docking is requested, it is read in here:
  if (allocated(env%directed_file)) then
    call read_directed_input(env)
  end if

end subroutine read_qcg_input

!==============================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!==============================================================================!

!> Read input for directed docking
subroutine read_directed_input(env)
  use crest_parameters
  use crest_data
  implicit none

  type(systemdata)           :: env

  integer                    :: nlines
  integer                    :: io,ich,i,i_check
  integer                    :: index
  character(len=512)         :: dum
  character(len=1),parameter :: delim_space = ' ',delim_tab = achar(9)

  open (newunit=ich,file=env%directed_file)
  !First check number of lines
  nlines = 0
  do
    read (ich,*,iostat=io)
    if (io /= 0) exit
    nlines = nlines+1
  end do
  !Allocate directed list
  !First entry is the atom number, Second how many solvents to add to this atom
  allocate (env%directed_list(nlines,2))
  allocate (env%directed_number(nlines),source=0)
  !Now read lines into directed_list
  rewind (ich)
  do i = 1,nlines
    read (ich,'(A)') dum
    !> Remove leading tab and spaces first
    dum = adjustl(dum) !Leading spaces are removed
    index = SCAN(trim(dum),delim_tab)
    if (index == 1) then !Leading tab -> remove it
      dum = dum(2:)
    end if
    index = SCAN(trim(dum),delim_space)
    if (index == 0) then !No space = check for tab
      index = SCAN(trim(dum),delim_tab)
    end if
    if (index == 0) then  !Second value is missing
      write (stdout,'(a,1x,i0)') "No second value found in directed list on line",i
      error stop
    end if
    env%directed_list(i,1) = dum(1:index-1)
    env%directed_list(i,2) = dum(index+1:)
    !Remove multiple spaces
    env%directed_list(i,2) = adjustl(env%directed_list(i,2))
    !Check, if spaces are still in second argument (e.g. a third number is giveb)
    index = SCAN(trim(env%directed_list(i,2)),delim_space)
    if (index == 0) index = SCAN(trim(dum),delim_tab)
    if (index /= 0) then
      write (stdout,'(a,1x,i0)') "Too many values at line",i
      error stop
    end if
    !> Make array with which solvent molecule at which atom to add
    read (env%directed_list(i,2),*,iostat=io) env%directed_number(i)
    env%directed_number(i) = sum(env%directed_number)
    if (io /= 0) then
      write (stdout,'(a,1x,i0)') "Second value is no number in line",i
      error stop
    end if
  end do
  close (ich)
  write (stdout,*) 'Performing directed docking'
  do i = 1,nlines
    write (stdout,'(a,1x,a,1x,a,1x,a)') 'Docking',trim(env%directed_list(i,2)),&
           & 'solvent molecules at',trim(env%directed_list(i,1))
  end do

end subroutine read_directed_input

!==============================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!==============================================================================!

subroutine qcg_grow(env,solu,solv,clus,tim)
  use crest_parameters
  use crest_data
  use qcg_printouts
  use iomod
  use qcg_coord_type
  use strucrd
  use qcg_utils
  implicit none

  type(systemdata)           :: env
  type(coord_qcg)            :: solu,solv,clus
  type(timer)                :: tim

  integer                    :: minE_pos,m
  integer                    :: iter = 1
  integer                    :: i,j,io,v
  integer                    :: max_cycle
  integer                    :: nat_backup
  logical                    :: e_there,high_e,success,neg_E
  real(wp)                   :: etmp(500)
  real(wp),allocatable       :: e_each_cycle(:)
  real(wp)                   :: dens,dum,efix
  real(wp)                   :: e_diff = 0.0_wp
  real(wp),allocatable       :: E_inter(:)
  real(wp)                   :: shr = 0.0_wp
  real(wp)                   :: shr_av = 0.0_wp
  real(wp)                   :: mean = 0.0_wp
  real(wp)                   :: mean_old = 0.0_wp
  real(wp)                   :: mean_diff = 0.0_wp
  character(len=*),parameter :: outfmt = '(1x,1x,a,1x,f14.7,a,1x)'
  character(len=512)         :: thispath,resultspath
  character(len=20)          :: gfnver_tmp
  integer                    :: ich99,ich15,ich88
  character(len=LEN(env%solv)) :: solv_tmp
  logical                    :: gbsa_tmp

  if (env%nsolv .gt. 0) then
    allocate (e_each_cycle(env%nsolv))
    allocate (E_inter(env%nsolv))
  else
    allocate (e_each_cycle(env%max_solv))
    allocate (E_inter(env%max_solv))
  end if

  call tim%start(5,'QCG Grow')

  call pr_eval_solute()
  call print_qcg_grow()
  call getcwd(thispath)
  io = makedir('grow')
  call chdirdbug('grow') !Results directory

!--- Output Files
  open (newunit=ich99,file='qcg_energy.dat')
  write (ich99,'(i0,2F20.8)') 0,solu%energy,solv%energy
  open (newunit=ich15,file='qcg_grow.xyz')  ! for molden movie
  open (newunit=ich88,file='qcg_conv.dat')  ! for convergence check
  write (ich88,'(''   #   Energy       Run. Aver.   Diff / au.'')')

  call getcwd(resultspath)
  call chdirdbug(thispath)

  if (env%water) then
    if (.not.env%user_wscal) then
      if (solu%nat .lt. 18) then
        env%potscal = 0.7_wp
      else
        env%potscal = 0.8_wp
      end if
      write (stdout,*)
      write (stdout,'(2x,''Water as solvent recognized,&
              & adjusting scaling factor for outer wall pot to '',F4.2)')&
             & env%potscal
      write (stdout,*)
    end if
  end if
  if (env%constrain_solu) write (stdout,'(2x,''Constraining solute during Growth '')')

  call get_ellipsoid(env,solu,solv,clus,.true.)
  call pr_grow_energy()

  call chdirdbug(env%scratchdir)
  v = makedir('tmp_grow')
  if (env%fixfile /= 'none selected') then
    call copysub(env%fixfile,'tmp_grow')
  end if
  call chdirdbug('tmp_grow')
  call solu%write('solute')
  call solv%write('solvent')
  call env%wrtCHRG('') !Write .CHRG file for docking

  call ellipsout('solute_cavity.coord',clus%nat,clus%at,clus%xyz,solu%ell_abc)
  solv%ell_abc = clus%ell_abc

  clus%chrg = solu%chrg
  clus%uhf = solu%uhf

  if (env%nsolv .gt. 0) then
    max_cycle = env%nsolv !User set number of solvents to add
  else
    max_cycle = env%max_solv !No solvent number set
  end if

!--------------------------------------------------------
! Start Loop
!--------------------------------------------------------
  GROW_LOOP: do iter = 1,max_cycle
    e_there = .false.
    success = .false.
    high_e = .false.
    neg_E = .false.
!---- Computation
    if (iter .gt. 1) then
      call get_ellipsoid(env,solu,solv,clus,.false.)
    end if

    call both_ellipsout('twopot_1.coord',clus%nat,clus%at,clus%xyz,&
           & clus%ell_abc,solu%ell_abc)

    CHECKWALL: do while (.not.success) !For restart with larger wall pot
      if (iter .eq. 1) then
        call xtb_dock(env,'solute','solvent',solu,solv)
        call check_dock(neg_E)

!-- If Interaction Energy is not negativ and existent, wall pot. too small and increase
        if (neg_E) then
          success = .true.
        else
          if (env%potscal .lt. 1.0_wp) then
            write (stdout,*) '  Wall Potential too small, increasing size by 5 %'
            solv%ell_abc = solv%ell_abc*1.05_wp
            env%potscal = env%potscal*1.05_wp
            if (env%potscal .gt. 1.0_wp) env%potscal = 1.0_wp
            write (stdout,'(''   New scaling factor '',F4.2)') env%potscal
          else
            success = .true.
          end if
        end if
      else
        call xtb_dock(env,'cluster.coord','solvent',solu,clus)
        call check_dock(neg_E)

        if (neg_E) then
          success = .true.
        else
          if (env%potscal .lt. 1.0_wp) then
            write (stdout,*) '  Wall Potential too small, increasing size by 5 %'
            clus%ell_abc = clus%ell_abc*1.05_wp
            env%potscal = env%potscal*1.05_wp
            if (env%potscal .gt. 1.0_wp) env%potscal = 1.0_wp
            write (stdout,'(''   New scaling factor '',F4.2)') env%potscal
          else
            success = .true.
          end if
        end if
      end if
    end do CHECKWALL

!--- Increase cluster size
    nat_backup = clus%nat
    call clus%deallocate
    clus%nat = nat_backup+solv%nat
    allocate (clus%at(clus%nat))
    allocate (clus%xyz(3,clus%nat))
    clus%nmol = clus%nmol+1

!--- Select xtb-IFF stucture to proceed
    call rdcoord('best.xyz',clus%nat,clus%at,clus%xyz,clus%energy)

    call remove('cluster.coord')
    call clus%write('cluster.coord')
    call both_ellipsout('twopot_2.coord',clus%nat,clus%at,clus%xyz,&
           & clus%ell_abc,solu%ell_abc)

    success = .false.

!--- Cluster restart, if interaction energy not negativ (wall pot. too small)
    gfnver_tmp = env%gfnver !> backup original level of theory
    do while (.not.success)
!--- Cluster optimization
      if (env%cts%used) then
        call write_reference(env,solu,clus) !new fixed file
      end if

!--- Interaction energy
      !gfnver_tmp = env%gfnver
      env%gfnver = env%lmover
      gbsa_tmp = env%gbsa
      solv_tmp = env%solv
      env%gbsa = .false.
      env%solv = ''
      call get_interaction_E(env,solu,solv,clus,iter,E_inter)
      env%gbsa = gbsa_tmp
      env%solv = solv_tmp
      if (E_inter(iter) .lt. 0) then
        success = .true.
      else
        if (env%potscal .lt. 1.0_wp) then
          write (stdout,*) '  Interaction Energy positiv, increasing outer wall pot by 5 %'
          clus%ell_abc = clus%ell_abc*1.05_wp
          env%potscal = env%potscal*1.05_wp
          if (env%potscal .gt. 1.0_wp) env%potscal = 1.0_wp
          write (stdout,'('' New scaling factor '',F4.2)') env%potscal
        else
          success = .true.
        end if
      end if
    end do
    env%gfnver = gfnver_tmp

!--- For output
    !Energy already read from xyz file
    e_each_cycle(iter) = clus%energy

!--- Calclulate fix energy + diff. energy
    efix = clus%energy/sqrt(real(clus%nat))
    dum = solu%energy
    if (iter .gt. 1) dum = e_each_cycle(iter-1)
    e_diff = e_diff+autokcal*(e_each_cycle(iter)-solv%energy-dum)
    call ellipsout('cluster_cavity.coord',clus%nat,clus%at,clus%xyz,clus%ell_abc)
    call both_ellipsout('twopot_cavity.coord',clus%nat,clus%at,clus%xyz,&
           & clus%ell_abc,solu%ell_abc)

!--- Density calculations
    call get_sphere(.false.,clus,.false.) !V, A of new cluster
    dens = 0.001*(solu%mass+iter*solv%mass)/(1.0d-30*clus%vtot*bohr**3)

!--- Movie file
    write (ich15,*) clus%nat
    write (ich15,'('' SCF done '',2F16.8)') autokcal*(e_each_cycle(iter)-solv%energy-dum)
    do j = 1,clus%nat
      write (ich15,'(a,1x,3F24.10)') i2e(clus%at(j)),clus%xyz(1:3,j)*bohr
    end do

!--- Output
    ! dist of new mol from solute for output
    call analyze_cluster(iter,clus%nat,solu%nat,solv%nat,clus%xyz,clus%at,shr_av,shr)

    write (stdout,'(x,i4,F13.6,1x,f7.2,3x,es9.2,5x,f6.3,3x,f8.3,3x,2f6.1,2x,f8.1,3x,a,x)') &
          & iter,e_each_cycle(iter),autokcal*(e_each_cycle(iter)-solv%energy-dum),&
          & e_diff,dens,efix,shr_av,shr,clus%vtot,trim(optlevflag(env%optlev))
    write (ich99,'(i4,F20.10,3x,f8.1)') iter,e_each_cycle(iter),clus%vtot

!--- Calculate moving average
    mean_old = mean
    do i = 0,iter-1
      mean = mean+E_inter(iter-i)
    end do
    mean = mean/iter
    mean_diff = mean-mean_old
    write (ich88,'(i5,1x,3F13.8)') iter,E_inter(iter)*autokcal,mean,mean_diff

!--- Check if converged when no nsolv was given
    if (env%nsolv .eq. 0) then
      if (abs(mean_diff) .lt. 1.0d-4.and.iter .gt. 5) then
        env%nsolv = iter
        exit
      end if
      if (iter .eq. env%max_solv) then
        write (stdout,'(1x,''No convergence could be reached upon adding'',1x,i4,1x,&
                & ''solvent molecules.'')') env%max_solv
        write (stdout,*) ' Proceeding.'
        env%nsolv = env%max_solv
        exit
      end if
    end if
!-----------------------------------------------
! End loop
!-----------------------------------------------
  end do GROW_LOOP

  if (env%nsolv .eq. 0) env%nsolv = iter !if no env%solv was given

  if (env%gfnver .ne. '--gfn2'.and.env%final_gfn2_opt) then
    gfnver_tmp = env%gfnver
    env%gfnver = '--gfn2'
    write (stdout,'(2x,''Final gfn2 optimization'')')
    call opt_cluster(env,solu,clus,'cluster.coord',.false.)
    call rdcoord('xtbopt.coord',clus%nat,clus%at,clus%xyz)
    call clus%write('cluster.coord')
    call grepval('xtb_sp.out','| TOTAL ENERGY',e_there,clus%energy)
    if (.not.e_there) then
      write (stdout,'(1x,a)') 'Total Energy of cluster not found.'
    else
      write (stdout,'(2x,''Total gfn2-energy of cluster/Eh:'',f20.6)') clus%energy
    end if
    env%gfnver = gfnver_tmp
  end if

  call wrxyz('cluster.xyz',clus%nat,clus%at,clus%xyz*bohr)

!--- One optimization without Wall Potential and with implicit model
  gfnver_tmp = env%gfnver
  if (env%final_gfn2_opt) env%gfnver = '--gfn2'
  call opt_cluster(env,solu,clus,'cluster.xyz',.true.)
  env%gfnver = gfnver_tmp
  call rename('xtbopt.xyz','cluster_optimized.xyz')
  call copysub('cluster_optimized.xyz',resultspath)

!--- output and files
  write (stdout,*)
  write (stdout,'(2x,''Growth finished after '',i0,'' solvents added'')') env%nsolv
  write (stdout,'(2x,''Results can be found in grow directory'')')
  write (stdout,'(2x,''Energy list in file <qcg_energy.dat>'')')
  write (stdout,'(2x,''Interaction energy in file <qcg_conv.dat>'')')
  write (stdout,'(2x,''Growing process in <qcg_grow.xyz>'')')
  write (stdout,'(2x,''Final geometry after grow in <cluster.coord> and <cluster.xyz>'')')
  write (stdout,'(2x,''Final geometry optimized without wall potential in <cluster_optimized.xyz>'')')
  write (stdout,'(2x,''Potentials and geometry written in <cluster_cavity.coord> and <twopot_cavity.coord>'')')

  close (ich99)
  close (ich88)
  close (ich15)

!--- Saving results and cleanup
  call copysub('cluster.coord',resultspath)
  call copysub('cluster.xyz',resultspath)
  call copysub('twopot_cavity.coord',resultspath)
  call copysub('cluster_cavity.coord',resultspath)
  call copysub('solute_cavity.coord',resultspath)
!  call rename('xcontrol','wall_potential')
  env%constrain_solu = .false.
  call write_wall(env,solu%nat,solu%ell_abc,clus%ell_abc,'wall_potential')
  call copysub('wall_potential',resultspath)

  call chdirdbug(thispath)
  call chdirdbug(env%scratchdir)
  if (.not.env%keepModef) call rmrf('tmp_grow')

  deallocate (e_each_cycle,E_inter)

  call tim%stop(5)

end subroutine qcg_grow

!==============================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!==============================================================================!
!
subroutine qcg_ensemble(env,solu,solv,clus,ens,tim,fname_results)
  use crest_parameters
  use crest_data
  use cregen_interface
  use crest_calculator
  use iomod
  use parse_xtbinput
  use qcg_coord_type
  use qcg_printouts
  use qcg_utils
  use strucrd
  use utilities
  use parallel_interface,only:crest_sploop !explicit interface needed (optional arguments)
  implicit none

  type(systemdata) :: env
  type(coord_qcg)  :: solu,solv,clus
  type(ensemble)   :: ens,dum
  type(timer)      :: tim

  integer :: i,j,k,io,f,r,ich,T,Tn,minpos
  character(len=512) :: thispath,resultspath,tmppath,tmppath2
  character(len=512) :: scratchdir_tmp
  character(len=512) :: jobcall
  character(len=256) :: inpnam,outnam
  character(len=80)  :: fname,pipe,to
  character(len=*)   :: fname_results
  character(len=64)  :: comment
  character(len=:),allocatable :: gfnver_tmp
  character(len=:),allocatable :: solv_tmp
  logical              :: gbsa_tmp,ex,mdfail,e_there,checkiso_tmp,cbonds_tmp
  real(wp),allocatable :: e_fix(:),e_clus(:)
  real(wp)             :: S,H,G,dens,shr,shr_av
  real(wp)             :: sasa
  real(wp)             :: newtemp,newmdtime,newmdstep,newhmass
  real(wp)             :: newmetadlist,newmetadexp,newmetadfac
  real(wp)             :: optlev_tmp
  real(wp)             :: e0
  real(wp),allocatable :: de(:)
  real(wp),allocatable :: p(:)
  integer              :: ich98,ich65,ich48
  logical              :: not_param = .false.
  type(timer)          :: tim_dum !Dummy timer to avoid double counting
  type(calcdata) :: calc_tmp
  logical,parameter :: debug = .true.

  if (.not.env%solv_md) then
    call tim%start(6,'QCG Solute-Ensemble')
  else
    call tim%start(7,'QCG Solvent-Ensemble')
  end if

  call tim_dum%init(20)

!--- Setting up directories
  call getcwd(thispath)
  f = makedir(fname_results)
  call chdirdbug(fname_results)
  call getcwd(resultspath)
  call chdirdbug(thispath)

!--- Setting defaults and backups
  env%cts%NCI = .true.  !Activating to have wall pot. written in coord file for xtb
  optlev_tmp = env%optlev
  env%optlev = 0.0d0
  gbsa_tmp = env%gbsa
  solv_tmp = env%solv
  env%gbsa = .false.
  env%solv = ''

!--- Setting up potential constraints
  allocate (env%cts%pots(10))
  env%cts%pots = ''
  write (env%cts%pots(1),'("$wall")')
  write (env%cts%pots(2),'(2x,"potential=polynomial")')
  write (env%cts%pots(3),'(2x,"ellipsoid:",1x,3(g0,",",1x),"all")') clus%ell_abc
  if (.not.env%solv_md) write (env%cts%pots(4),'(2x,"ellipsoid:",1x,3(g0,",",1x),"1-",i0)')&
         & solu%ell_abc,solu%nat

  if (env%cts%used) then
    call write_reference(env,solu,clus) !new fixed file
    call copysub(env%fixfile,env%scratchdir)
  end if

  call chdirdbug(env%scratchdir)
  scratchdir_tmp = env%scratchdir
  if (.not.env%solv_md) then
    io = makedir('tmp_MTD')
    call copysub('.CHRG','tmp_MTD')
    call copysub('.UHF','tmp_MTD')
    if (env%cts%used) call copysub(env%fixfile,'tmp_MTD')
    call chdirdbug('tmp_MTD')
  else
    io = makedir('tmp_solv_MTD')
    call chdirdbug('tmp_solv_MTD')
  end if
  call getcwd(tmppath2)
  call clus%write('crest_input')

  if (env%solv_md) then
    call wr_cluster_cut('crest_input',solu%nat,solv%nat,env%nsolv,&
           & 'solute_cut.coord','solvent_shell.coord')
    call remove('crest_input')
    call copy('solvent_shell.coord','crest_input')
    call clus%open('solvent_shell.coord')
  end if

  !For newcregen: If env%crestver .eq. crest_solv .and. .not. env%QCG then conffile .eq. .true.
  env%QCG = .false.
  call inputcoords(env,'crest_input')
  call defaultGF(env)         !Setting MTD parameter

!--- Special constraints for gff to safeguard stability
  if (env%ensemble_opt .eq. '--gff') then
    checkiso_tmp = env%checkiso
    env%checkiso = .true.
    cbonds_tmp = env%cts%cbonds_md
    env%cts%cbonds_md = .true.
    call autoBondConstraint_withEZ('coord',env%forceconst,env%wbofile)
    call rd_cbonds('bondlengths',env)
  end if

  gfnver_tmp = env%gfnver
  write (stdout,*) 'Method for ensemble search: ',env%ensemble_opt
!  if (env%ens_const) write(stdout,*) '  Solute fixed during ensemble generation'
  env%gfnver = env%ensemble_opt  !Setting method for ensemble search

  !----------------------------------------------------------------
  ! Case selection of normal Crest, MD or MTD
  !----------------------------------------------------------------

  !> Parse contraints (wall potentials etc.) into new calculator
  !> if we are using it.
  if (.not.env%legacy) then
    calc_tmp = env%calc
    call qcg_envcalc_reinit(env,clus,.true.,.true.)
  end if

  ENSEMBLEGEN:select case(env%ensemble_method)
  case (-1:0) !qcgmtd/Crest runtype

  !> Some custom Defaults for running the standard search
  !General settings:
  if (.not.env%user_mdstep) then
    if (env%ensemble_opt .EQ. '--gff') then
      env%mdstep = 1.5d0
    else
      env%mdstep = 5.0d0
    end if
  end if
  !Runtype specific settings:
  if (env%ensemble_method == 0) then
    if (.not.env%user_dumxyz) then
      env%mddumpxyz = 200
    end if
    if (.not.env%user_mdtime) then
      env%mdtime = 10.0
    end if
  else if (env%ensemble_method == -1) then
    if (.not.env%user_dumxyz) then
      env%mddumpxyz = 50
    end if
    if (.not.env%user_mdtime) then
      env%mdtime = 5.0
    end if
    env%nmdtemp = 100
    env%MaxRestart = 6
  end if

  env%iterativeV2 = .true.  !Safeguards more precise ensemble search
  write (stdout,*) 'Starting ensemble cluster generation by CREST routine'
  call confscript2i(env,tim_dum) !Calling ensemble search
  call copy('crest_rotamers.xyz','crest_rotamers_0.xyz')

  case (1:2) ! Single MD or MTD
  call xtb_md_ensemble_qcg(env,solu,solv,clus,resultspath)

  end select ENSEMBLEGEN

  env%QCG = .true.

!--- Optimization with gfn2 if necessary
  if (env%final_gfn2_opt.and.env%gfnver .ne. '--gfn2') then
    write (stdout,'(2x,a)') 'GFN2-xTB optimization'
    env%gfnver = '--gfn2'

    if (.not.env%legacy) then
      !> reinit calculator with GFN2
      call qcg_envcalc_reinit(env,clus,.true.,.true.)
      call checkname_xyz(crefile,inpnam,outnam)
      call crest_multilevel_wrap(env,trim(inpnam),0)
    else
      call rmrf('OPTIM')
      call multilevel_opt(env,99)
    end if
    write (stdout,*)
  end if

!--- Final optimization without wall potentials
  env%optlev = 1.0d0    !Higher precision for less scattering
  env%cts%NCI = .false.  !Dactivating the wall pot.
  env%cts%pots = ''
  deallocate (env%cts%pots)

  if (.not.env%legacy) then
    !> wall potential was turned off, add any other constraint back in
    call qcg_envcalc_reinit(env,clus,.true.,.true.)
    call checkname_xyz(crefile,inpnam,outnam)
    call crest_multilevel_wrap(env,trim(inpnam),0)
  else
    call rmrf('OPTIM')
    call multilevel_opt(env,99)
  end if

!--- Clustering to exclude similar structures if requested with -cluster
  if (env%properties == 70) then
    write (stdout,'(3x,''Clustering the remaining structures'')')
    call checkname_xyz(crefile,inpnam,outnam)
    call ccegen(env,.false.,inpnam)
    call move(trim(clusterfile),trim(outnam))
  end if

!--- Energy sorting and removal of dublicates
  env%gbsa = gbsa_tmp
  env%solv = solv_tmp
  call newcregen(env,0)
  call checkname_xyz(crefile,inpnam,outnam)
  call copy(inpnam,'ensemble.xyz')
  call ens%open('ensemble.xyz') !Read in ensemble
  call clus%deallocate()
  clus%nat = ens%nat
  allocate (clus%at(clus%nat))
  allocate (clus%xyz(3,clus%nat))
  write (stdout,'(1x,i0,a)') ens%nall,' structures remaining.'
  write (stdout,*)

!-------------------------------------------------------------
!      SP with Implicit solvation model and without wall potentials
!-------------------------------------------------------------
  if (env%legacy) then
    !> old, I/O-heavy version
    call ens_sp_with_io(env,ens,clus,resultspath)
  else
    !> use internal parallel loop, but remember to convert to Bohrs for that
    clus%at(:) = ens%at(:)
    clus%xyz(1:3,1:clus%nat) = ens%xyz(1:3,1:ens%nat,1)*aatoau
    call qcg_envcalc_reinit(env,clus,.true.,.true.)

    !> build the coord list (in Bohr) for the internal parallel loop
    block
      type(coord),allocatable :: structures(:)
      integer :: isp
      allocate (structures(ens%nall))
      do isp = 1,ens%nall
        structures(isp)%nat = ens%nat
        structures(isp)%at = ens%at
        structures(isp)%xyz = ens%xyz(1:3,1:ens%nat,isp)*aatoau
      end do
      call crest_sploop(env,ens%nall,structures,ens%er)
      deallocate (structures)
    end block
  end if

!-------------------------------------------------------------
!      Processing results
!-------------------------------------------------------------
  env%gfnver = gfnver_tmp
  allocate (e_fix(ens%nall),source=0.0_wp)
  allocate (e_clus(ens%nall),source=0.0_wp)

  call pr_ensemble_energy()

  open (newunit=ich98,file='cluster_energy.dat')
  write (ich98,'(3x,''#'',9x,''Energy [Eh]'',6x,''SASA'')')

!--- Fixation energy of optimization
  do i = 1,ens%nall
    if (env%legacy) then
      !> old I/O-heady version
      call chdirdbug('OPTIM')
      write (to,'("TMPCONF",i0)') i
      call chdirdbug(to)
      call grepval('xtb.out','         :: add. restraining',e_there,e_fix(i))
      call chdirdbug(tmppath2)
      call rdxmolselec('full_ensemble.xyz',i,clus%nat,clus%at,clus%xyz)
    else
      !> quicker version, simply load from 'ens'
      call ens%get_mol(i,clus)
    end if

    call get_sphere(.false.,clus,.false.)
    dens = 0.001*(solu%mass+env%nsolv*solv%mass)/(1.0d-30*clus%vtot*bohr**3)
    if (env%solv_md) then
      call analyze_cluster(env%nsolv-1,clus%nat,solv%nat,solv%nat,clus%xyz,clus%at,shr_av,shr)
    else
      call analyze_cluster(env%nsolv,clus%nat,solu%nat,solv%nat,clus%xyz,clus%at,shr_av,shr)
    end if
    write (ich98,'(i4,F20.10,3x,f8.1)') env%nsolv,ens%er(i),clus%atot
    write (stdout,'(x,i4,4x,F13.6,2x,f6.3,1x,f8.3,2x,2f6.1,3x,f8.1,3x,a)') &
          & i,ens%er(i),dens,e_fix(i),shr_av,shr,clus%atot,trim(optlevflag(env%optlev))
    e_fix(i) = e_fix(i)*autokcal/sqrt(real(clus%nat,wp))
  end do
  close (ich98)
  call copysub('cluster_energy.dat',resultspath)

!--- Checking Boltzmann weighting
  write (stdout,*)
  call remove('full_ensemble.xyz')
  call qcg_dump_sorted_ensemble(ens,ens%er,'full_ensemble.xyz')
  e_clus = ens%er*autokcal
  call sort_min(ens%nall,1,1,e_clus)
  ens%er = e_clus/autokcal !Overwrite ensemble energy with sorted one
  allocate (de(ens%nall),source=0.0d0)
  allocate (p(ens%nall),source=0.0d0)
  e0 = e_clus(1)
  de(1:ens%nall) = (e_clus(1:ens%nall)-e0)
  call qcg_boltz(env,ens%nall,de,p)
  k = 0
  if (.not.env%user_nclust) env%nqcgclust = 0 !Needed for solvent ensemble
  if (env%nqcgclust .eq. 0) then
    do i = 1,ens%nall !Count how many are above 10%
      if ((p(i)) .gt. 0.1) then
        k = k+1
      end if
    end do
    if ((k .eq. 0).or.(k .gt. 10)) then
      k = 10 !If too many structures are relevant, set it 10
    else if ((k .lt. 4).and.(ens%nall .ge. 4)) then
      k = 4 !If too less structures are relevant, set it 4
    else if (ens%nall .gt. 0) then
      k = ens%nall
    else
      error stop 'No structure left. Something went wrong.'
    end if
    write (stdout,'(2x,a,1x,i0)') 'Conformers taken:',k
    env%nqcgclust = k
  else
    if (env%nqcgclust .gt. ens%nall) then
      k = ens%nall !Input larger than remaining structures
      write (stdout,'(''Less than '',1x,i0,1x,''structures remain'')') env%nqcgclust
      write (stdout,'(''Only '',1x,i0,1x,''structures are taken'')') ens%nall
      if (env%cff) env%nqcgclust = ens%nall !Only for CFF, else a second qcg_ensemble run starts for solvent
    else
      write (stdout,'(''Taking '',1x,i0,1x,''structures'')') env%nqcgclust
      k = env%nqcgclust !user input
    end if
  end if

  open (newunit=ich65,file='final_ensemble.xyz')
  do i = 1,k
    open (newunit=ich48,file='full_population.dat')
    write (ich48,'(2x, ''cluster'',2x,''E_norm [Eh]'',2x, ''De [kcal]'', 4x, ''p'')')
    do j = 1,ens%nall
      if (j .lt. 10) then
        write (ich48,'(5x,i0,3x,f11.6,5x,f6.4,3x,f6.4)') j,e_clus(j)/autokcal,de(j),p(j)
      else
        write (ich48,'(5x,i0,2x,f11.6,5x,f6.4,3x,f6.4)') j,e_clus(j)/autokcal,de(j),p(j)
      end if
    end do
    close (ich48)

!--- Take k energetic least structures (written at beginning of file)
    call rdxmolselec('full_ensemble.xyz',i,clus%nat,clus%at,clus%xyz)
    call wrxyz(ich65,clus%nat,clus%at,clus%xyz*bohr,ens%er(i))
  end do
  close (ich65)

  call ens%deallocate()
  call ens%open('final_ensemble.xyz')
  ens%er = e_clus(1:k)/autokcal

!--- Getting G,S,H
  write (stdout,*)
  write (stdout,'(2x,70("-"))')
  write (stdout,'(2x,70("-"))')
  write (stdout,'(2x,''Boltz. averaged energy of final cluster:'')')
  call aver(.true.,env,ens%nall,e_clus(1:ens%nall),S,H,G,sasa,.false.)
  write (stdout,'(7x,''G /Eh     :'',f15.8)') G/autokcal
  write (stdout,'(7x,''T*S /kcal :'',f15.8)') S

  ens%g = G
  ens%s = S

  deallocate (e_fix)
  deallocate (e_clus)

!---Folder management
  call rename('cregen.out.tmp','thermo_data')
  call copysub('thermo_data',resultspath)
  call copysub('crest_best.xyz',resultspath)
  call copysub('cre_members.out',resultspath)
  call copysub('full_ensemble.xyz',resultspath)
  call copysub('final_ensemble.xyz',resultspath)
  call copysub('population.dat',resultspath)
  call copysub('full_population.dat',resultspath)

!---Deleting ensemble tmp
  call chdirdbug(thispath)
  call chdirdbug(env%scratchdir)
  if (.not.env%keepModef) call rmrf(tmppath2)
!----Outprint
  write (stdout,*)
  write (stdout,'(2x,"Ensemble generation finished.")')
  write (stdout,'(2x,"Results can be found in the [ensemble] directory:")')
  write (stdout,'(2x,"--> What?                      --> Where?")')
  write (stdout,'(2x,"Lowest energy conformer        crest_best.xyz")')
  write (stdout,'(2x,"List of full ensemble          full_ensemble.xyz")')
  write (stdout,'(2x,"List of used ensemble          final_ensemble.xyz")')
  write (stdout,'(2x,"Ensemble thermodyn data        thermo_data")')
  write (stdout,'(2x,"Population of selected         population.dat")')
  write (stdout,'(2x,"Population of full ensemble    full_population.dat")')

  !>--- restore settings
  env%gfnver = gfnver_tmp
  env%optlev = optlev_tmp
  if (env%ensemble_opt .eq. '--gff') then
    env%cts%cbonds_md = cbonds_tmp
    env%checkiso = checkiso_tmp
  end if

  call tim_dum%clear

  if (.not.env%solv_md) then
    call tim%stop(6)
  else
    call tim%stop(7)
  end if

end subroutine qcg_ensemble

!==============================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!==============================================================================!

subroutine qcg_cff(env,solu,solv,clus,ens,solv_ens,tim)
  use crest_parameters
  use crest_data
  use qcg_printouts
  use iomod
  use qcg_coord_type
  use strucrd
  use qcg_utils
  implicit none

  type(systemdata)           :: env
  type(timer)                :: tim
  type(coord_qcg)            :: solu,solv,clus
  type(ensemble),intent(inout) :: solv_ens
  type(ensemble),intent(in)    :: ens

  integer                    :: i,j,k,iter
  integer                    :: io,r
  integer                    :: n_ini,nopen
  integer                    :: ipos,dum
  integer                    :: v_ratio
  integer                    :: minE_pos,m,nat_tot
  integer                    :: nat_frag1 !number of atoms larger fragment (=solvent shell)
  integer                    :: conv(env%nqcgclust+1)
  integer                    :: solv_added,minpos
  character(len=512)         :: thispath,resultspath,tmppath,tmppath2
  character(len=64)          :: comment
  character(len=20)          :: to
  real(wp),allocatable      :: e_empty(:),inner_ell_abc(:,:)
  real(wp),allocatable      :: outer_ell_abc(:,:)
  real(wp),allocatable      :: e_cur(:,:)
  real(wp)                   :: e_cluster(env%nqcgclust)
  real(wp)                   :: S,H,G
  real(wp)                   :: sasa,tmp_optlev
  real(wp)                   :: etmp(500)
  real(wp)                   :: e_fix(env%nqcgclust),e_norm(env%nqcgclust)
  real(wp)                   :: dum_e,de
  real(wp)                   :: de_tot(env%nqcgclust)
  real(wp)                   :: shr = 0
  real(wp)                   :: shr_av = 0
  real(wp)                   :: dens,atotS
  logical                    :: ex,skip,e_there
  logical                    :: all_converged
  logical,allocatable        :: converged(:),nothing_added(:)

  character(len=20)          :: gfnver_tmp
  real(wp)                   :: optlev_tmp
  integer                    :: ich98,ich31
  integer                    :: chrg_tmp,uhf_tmp

  call tim%start(8,'QCG CFF')

  allocate (e_empty(env%nqcgclust))
  allocate (converged(env%nqcgclust))
  allocate (nothing_added(env%nqcgclust))
  allocate (outer_ell_abc(env%nqcgclust,3))
  allocate (inner_ell_abc(env%nqcgclust,3))

  v_ratio = nint(solu%vtot/solv%vtot)
  allocate (e_cur(env%nsolv+v_ratio,env%nqcgclust),source=0.0d0)

!--- Setting defaults (same as ensemble optimization to have comparable structures)
  optlev_tmp = env%optlev
  env%optlev = 1.0d0    !Increaseing percision for ensemble search to minimze scattering
  gfnver_tmp = env%gfnver
  env%gfnver = qcg_energy_level(env) !CFF always at the level of the ensemble SPs
!--- Only pure solvent clusters are computed here, neutral closed-shell
  chrg_tmp = env%chrg
  uhf_tmp = env%uhf
  env%chrg = 0
  env%uhf = 0
  nothing_added = .false.

  dum = 0
  converged = .false.
  all_converged = .false.
  nat_tot = clus%nat-solu%nat!*env%nqcgclust

  if (solu%vtot/solv%vtot .lt. 1.0d0) then
    skip = .true.
  else
    skip = .false.
  end if

!--- Folder management
  call getcwd(thispath)
  r = makedir('solvent_ensemble')
  call chdirdbug('solvent_ensemble')
  call getcwd(resultspath)
  call chdirdbug(thispath)
  call chdirdbug(env%scratchdir)
  call getcwd(tmppath)
  io = makedir('tmp_CFF')
  call chdirdbug('tmp_CFF')
  call getcwd(tmppath2)
  call chdirdbug(tmppath)
  call chdirdbug('solvent_properties')
  call copysub('solvent',tmppath2)
  call chdirdbug(tmppath2)

!--- SP of each cluster
!    (works with legacy and calculator version since xtb_sp_qcg switches automatically)
  call ens%write('ensemble.xyz')
  do i = 1,env%nqcgclust
    call ens%get_mol(i,clus)
    clus%nmol = clus%nat/solv%nat

    write (to,'("TMPCFF",i0)') i
    io = makedir(trim(to))
    call copysub('solvent',to)
    call chdirdbug(to)

    call clus%write('cluster.coord')
    call wr_cluster_cut('cluster.coord',solu%nat,solv%nat,env%nsolv,'solute_cut.coord','solvent_shell.coord')
    call xtb_sp_qcg(env,'solvent_shell.coord',ex,e_empty(i))
    if (env%legacy) then
      call grepval('xtb.out','| TOTAL ENERGY',ex,e_empty(i))
    end if
    call copy('solvent_shell.coord','solvent_cluster.coord')
    call copy('solvent_cluster.coord','filled_cluster.coord')
    call get_ellipsoid(env,solu,solv,clus,.false.) !solu, to have same cavity to fill solvent in
    outer_ell_abc(i,1:3) = clus%ell_abc(1:3)
    inner_ell_abc(i,1:3) = solu%ell_abc(1:3)
    call chdirdbug(tmppath2)
  end do

  if (skip) write (stdout,'(2x,''solute smaller than solvent, cff skipped'')')

  clus%nat = clus%nat-solu%nat
  n_ini = clus%nat

!--- If solvent molecules are added
  if (.not.skip) then
    call pr_qcg_fill()
    write (stdout,'(2x,''now adding solvents to fill cluster...'')')
    call pr_fill_energy()
    write (stdout,'(2x,70("-"))')
    nat_frag1 = env%nsolv*solv%nat

    iter = 0
!--- Main cycle for addition of solvent molecules
    convergence: do while (.not.all_converged)
      k = 0
      iter = iter+1
      !--- Setting array, with only numbers of dirs that are not converged
      do i = 1,env%nqcgclust
        if (.not.converged(i)) then
          k = k+1
          conv(k) = i
        else
          cycle
        end if
      end do
      conv(k+1:env%nqcgclust) = 0
      conv(env%nqcgclust+1) = k !How many jobs are open
      nopen = k

      call chdirdbug(tmppath2)

!--- Solvent addition to the cluster---------------------------------------------
      call ensemble_dock(env,outer_ell_abc,nat_frag1,'solvent_cluster.coord',&
              &'solvent',clus%nat,solv%nat,nopen,'TMPCFF',conv)
!--------------------------------------------------------------------------------

      nat_frag1 = nat_frag1+solv%nat

      !--- Increase cluster size
      deallocate (clus%at)
      deallocate (clus%xyz)
      clus%nat = clus%nat+solv%nat
      allocate (clus%at(clus%nat))
      allocate (clus%xyz(3,clus%nat))
      clus%nmol = clus%nmol+1

      do i = 1,env%nqcgclust
        if (.not.converged(i)) then
          write (to,'("TMPCFF",i0)') i
          call chdirdbug(to)
          call remove('xtbrestart')
          call remove('xcontrol')

          call rdcoord('best.xyz',clus%nat,clus%at,clus%xyz,e_cur(iter,i))
          call clus%write('solvent_cluster.coord')

          !--- Check if converged
          call fill_take(env,solv%nat,clus%nat,inner_ell_abc(i,1:3),ipos)
          if (ipos .eq. 0) then
            converged(i) = .true.
            write (stdout,'(2x,a,i0,a)') &
              & "no more solvents can be placed inside cavity of cluster: ",i, &
              & ", taking previous."
            if (iter .eq. 1) nothing_added(i) = .true.
          end if
          call chdirdbug(tmppath2)

        else
          cycle
        end if
      end do

!--- Check, if a structure was converged and iff was not necessary
      k = 0
      do i = 1,env%nqcgclust
        if (.not.converged(i)) then
          k = k+1
          conv(k) = i
        else
          cycle
        end if
      end do
      conv(k+1:env%nqcgclust) = 0
      conv(env%nqcgclust+1) = k !How many jobs are open
      nopen = k

!--- Parallel optimization-------------------------------------------------------------------
      ! for some reason pre-processing with constraint is coupled to the pr flag
      ! also, I don't think this call does anything useful...
      ! I implemented e_cur readout, this makes sense to me at least
      if (nopen .gt. 0) then
        call cff_opt(.false.,env,'solvent_cluster.coord',n_ini,nopen&
                &,'TMPCFF',conv,nothing_added,e_cur(iter,:))
      end if
!----------------------------------------------------------------------------------------------
      de_tot(:) = 0.0_wp
      do i = 1,env%nqcgclust
        if (.not.converged(i)) then
          write (to,'("TMPCFF",i0)') i
          call chdirdbug(to)
          dum_e = e_empty(i)
          if (iter .gt. 1) dum_e = e_cur(iter-1,i)
          de = autokcal*(e_cur(iter,i)-solv%energy-dum_e)
          de_tot(i) = de_tot(i)+de
          !---- Check if solvent added is repulsive
          if (de .gt. 0) then
            converged(i) = .true.
            write (stdout,'(2x,"adding solvent is repulsive for cluster: ",i0,a)') i, &
            & ", taking previous one instead."
            if (iter .eq. 1) nothing_added(i) = .true.
          else !Only if the addition was not repulsive
            call copy('solvent_cluster.coord','filled_cluster.coord')
            write (stdout,'(i4,5x,i3,1x,F13.6,3x,f7.2,5x,f7.2,4x,a)') &
               & iter+env%nsolv,i,e_cur(iter,i),de,de_tot(i),&
               & trim(optlevflag(env%optlev))
          end if
          call chdirdbug(tmppath2)
        end if
      end do

      !--- Check if everything is converged
      dum = 0
      do i = 1,env%nqcgclust
        if (converged(i)) then
          dum = dum+1
        end if
      end do

      if (dum .eq. env%nqcgclust) then
        all_converged = .true.
      else
        nat_tot = nat_tot+solv%nat
      end if

      write (stdout,'(2x,70("-"))')
      !--- Or if maximum solvent is added
      if (iter .eq. v_ratio) then
        write (stdout,'(2x,''volume filled'')')
        all_converged = .true.
        call copy('solvent_cluster.coord','filled_cluster.coord')
      end if

    end do convergence

  end if

  !Now in every TMPPath the final cluster file filled_cluster.coord is present

!---------------------------------------------------------------------
!     Final Optimization
!---------------------------------------------------------------------

  tmp_optlev = env%optlev
  if (env%optlev .lt. 1.0) env%optlev = 1.0d0 !higher accuracy

  nopen = env%nqcgclust !final optimization is done for all clusters
  if (.not.skip) then
    call cff_opt(.true.,env,'filled_cluster.coord',n_ini,nopen,&
           & 'TMPCFF',conv,nothing_added,e_cluster)
  else
    n_ini = 0 !If this is 0, no constraining will be done (optimization of total system)
    nothing_added = .true.
    call cff_opt(.true.,env,'filled_cluster.coord',n_ini,nopen,'TMPCFF',&
           & conv,nothing_added,e_cluster)
  end if
  env%optlev = tmp_optlev

  call pr_ensemble_energy()

  solv_ens%nall = env%nqcgclust
  solv_ens%nat = nat_tot

!--- Getting results--------------------------------------------------------------
  open (newunit=ich31,file='crest_rotamers_0.xyz')
  open (newunit=ich98,file='cluster_energy.dat')
  write (ich98,'(3x,''#'',11x,''Energy [Eh]'',6x,''SASA'')')

  do i = 1,env%nqcgclust
    write (to,'("TMPCFF",i0)') i
    call chdirdbug(to)
    call copy('xtbopt.coord','final_cluster.coord')

!--- Reading structure
    call clus%open('final_cluster.coord')

!--- Getting energy and calculating properties
    if (env%legacy) then
      call grepval('xtb_sp.out','| TOTAL ENERGY',e_there,e_cluster(i))
      call grepval('xtb_sp.out','         :: add. restraining',e_there,e_fix(i))
    end if
    e_fix(i) = e_fix(i)*autokcal/sqrt(real(clus%nat))
    call get_sphere(.false.,clus,.false.)
    if (clus%nat .gt. n_ini) then
      solv_added = (clus%nat-(n_ini))/solv%nat
    else
      solv_added = 0
    end if
    dens = 0.001*((clus%nat/solv%nat)*solv%mass)/(1.0d-30*clus%vtot*bohr**3)
    call analyze_cluster(solv_added,clus%nat,n_ini,solv%nat,clus%xyz,clus%at,shr_av,shr)
    e_norm(i) = e_cluster(i)*env%nsolv/(clus%nat/solv%nat)
    atotS = clus%atot*env%nsolv/(clus%nat/solv%nat)

!--- Writing outputfiles
    write (ich31,'(2x,i0)') clus%nat
    write (ich31,'(2x,a,f18.8,2x,a)') 'energy=',e_cluster(i)
    do j = 1,clus%nat
      write (ich31,'(1x,a2,1x,3f20.10)') i2e(clus%at(j),'nc'),clus%xyz(1:3,j)*bohr
    end do

    write (ich98,'(''No'',i4,F20.10,3x,f8.1)') i,e_norm(i),atotS

!--- Print to screen
    write (stdout,'(x,i4,4x,F13.6,2x,f6.3,1x,f8.3,2x,2f6.1,3x,f8.1,3x,a)') &
            & i,e_norm(i),dens,e_fix(i),shr_av,shr,atotS,trim(optlevflag(env%optlev))

    call chdirdbug(tmppath2)
  end do

  close (ich98)
  close (ich31)

  call solv_ens%deallocate()
  call solv_ens%open('crest_rotamers_0.xyz')

  solv_ens%er = e_cluster
  call copy('crest_rotamers_0.xyz','crest_ensemble.xyz')

!--- crest_best structure
  minpos = minloc(solv_ens%er,dim=1)
  write (to,'("TMPCFF",i0)') minpos
  call chdirdbug(to)
  call clus%open('final_cluster.coord')

  clus%xyz = clus%xyz*bohr
  call chdirdbug(tmppath2)
  write (comment,'(a,F20.8)') 'energy=',solv_ens%er(minpos)
  call wrxyz('crest_best.xyz',clus%nat,clus%at,clus%xyz,comment)

!--- Boltz. average-------------------------------------------------------------------------
  write (stdout,*)
  write (stdout,'(2x,70("-"))')
  write (stdout,'(2x,70("-"))')
  write (stdout,'(2x,''Boltz. averaged energy of final cluster:'')')
  e_cluster = solv_ens%er*autokcal
  e_norm = e_norm*autokcal
  call sort_min(env%nqcgclust,1,1,e_norm)
  call aver(.true.,env,solv_ens%nall,e_norm(1:env%nqcgclust),S,H,G,sasa,.false.)
  write (stdout,'(7x,''G /Eh     :'',f15.8)') G/autokcal
  write (stdout,'(7x,''T*S /kcal :'',f15.8)') S
  solv_ens%er = e_norm/autokcal !normalized energy needed for final evaluation

  solv_ens%g = G
  solv_ens%s = S

!--- Cleanup
  call copysub('crest_ensemble.xyz',resultspath)
  call copysub('cluster_energy.dat',resultspath)
  call copysub('crest_best.xyz',resultspath)
  call copysub('population.dat',resultspath)
  call chdirdbug(tmppath)
  if (.not.env%keepModef) call rmrf('tmp_CFF')
  call chdirdbug(thispath)

!--- Printouts
  write (stdout,*)
  write (stdout,'(2x,"Solvent cluster generation finished.")')
  write (stdout,'(2x,"Results can be found in [solvent_cluster] directory")')
  write (stdout,'(2x,"--> What?       --> Where?")')
  write (stdout,'(2x,"Structures      crest_ensemble.xyz")')
  write (stdout,'(2x,"Energies        cluster_energy.dat")')
  write (stdout,'(2x,"Population      population.dat")')

  env%gfnver = gfnver_tmp
  env%optlev = optlev_tmp
  env%chrg = chrg_tmp
  env%uhf = uhf_tmp

  deallocate (e_empty)
  deallocate (converged)
  deallocate (outer_ell_abc)
  deallocate (inner_ell_abc)

  call tim%stop(8)

end subroutine qcg_cff

!==============================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!==============================================================================!

subroutine qcg_freq(env,tim,solu,solv,solu_ens,solv_ens)
  use crest_parameters
  use crest_data
  use qcg_printouts
  use iomod
  use qcg_coord_type
  use strucrd
  use qcg_utils
  implicit none

  type(systemdata)           :: env
  type(timer)                :: tim
  type(coord_qcg)            :: solu,solv,clus
  type(ensemble)             :: solu_ens,solv_ens

  integer                    :: r,io,f,g,h
  integer                    :: i
  character(len=512)         :: thispath,resultspath,tmppath,tmppath2
  character(len=80)          :: to
  character(len=20)          :: gfnver_tmp
  real(wp)                   :: optlev_tmp
  real(wp)                   :: gt(3)
  real(wp)                   :: ht(3)
  real(wp)                   :: svib(3)
  real(wp)                   :: srot(3)
  real(wp)                   :: stra(3)
  integer                    :: ich65,ich56,ich33,ich81
  logical                    :: opt
  integer                    :: chrg_tmp,uhf_tmp
  type(coord_qcg) :: tmpmol

  call tim%start(9,'QCG Frequencies')

  call pr_qcg_freq()

!--- Setting defaults (same as ensemble optimization and cff to have comparable structures)
  optlev_tmp = env%optlev
  env%optlev = 1.0d0    !Increaseing percision for ensemble search to minimze scattering
  gfnver_tmp = env%gfnver
  env%gfnver = env%freqver  !Setting method

!--- Folder management
  call getcwd(thispath)
  r = makedir('frequencies')
  call chdirdbug('frequencies')
  call getcwd(resultspath)
  call chdirdbug(thispath)
  call chdirdbug(env%scratchdir)
  call getcwd(tmppath)
  io = makedir('tmp_freq')
  call copysub('.CHRG','tmp_freq')
  call copysub('.UHF','tmp_freq')
  call chdirdbug('tmp_freq')
  call getcwd(tmppath2)
  f = makedir('tmp_solu')
  call copysub('.CHRG','tmp_solu')
  call copysub('.UHF','tmp_solu')
  g = makedir('tmp_solv')
  h = makedir('tmp_gas1') !One solute molecule
  call copysub('.CHRG','tmp_gas1')
  call copysub('.UHF','tmp_gas1')

!----------------------------------------------------------------------------
!   frequencies for solute molecule
!----------------------------------------------------------------------------
  write (stdout,'(1x,a)') 'processing  SOLUTE MOLECULE'
  call chdirdbug('tmp_gas1')
  call solu%write('solute.coord')
  call chdirdbug(tmppath2)
  opt = .false.
  call ens_freq(env,'solute.coord',1,'tmp_gas',opt)
  call chdirdbug('tmp_gas1')
  call rdtherm('xtb_freq.out',ht(3),svib(3),srot(3),stra(3),gt(3))
  solu%gt = gt(3)
  solu%ht = ht(3)
  solu%svib = svib(3)
  solu%srot = srot(3)
  solu%stra = stra(3)

  call chdirdbug(tmppath2)

!----------------------------------------------------------------------------
!   frequencies for solute cluster
!----------------------------------------------------------------------------
  write (stdout,'(/,1x,a)') 'processing  SOLUTE CLUSTER'

!--- Folder setup for cluster
  call chdirdbug('tmp_solu')
  call solu_ens%write('solute_ensemble.xyz')

!--- All cluster are of the same size
  call clus%deallocate()
  clus%nat = solu_ens%nat
  allocate (clus%at(clus%nat))
  allocate (clus%xyz(3,clus%nat))
  clus%xyz = 0
  clus%nmol = env%nsolv+1 !clus%nat/clus%at

  do i = 1,solu_ens%nall
    !call rdxmolselec('solute_ensemble.xyz',i,clus%nat,clus%at,clus%xyz)
    call solu_ens%get_mol(i,clus)

!--- Solute cluster
    write (to,'("TMPFREQ",i0)') i
    io = makedir(trim(to))
    call copysub('.UHF',to)
    call copysub('.CHRG',to)
    call chdirdbug(to)
    call clus%write("cluster.xyz")

    call chdirdbug(tmppath2)

    !--- Solvent cluster (only if cff, than the solvent shell is taken, which was fixed all the time)
    if (env%cff) then
      call chdirdbug('tmp_solv')
      write (to,'("TMPFREQ",i0)') i
      io = makedir(trim(to))
      call chdirdbug(to)
      call clus%write('cluster.coord')
      call wr_cluster_cut('cluster.coord',solu%nat,solv%nat,env%nsolv,&
             & 'solute_cut.coord','solvent_cut.coord')

      call chdirdbug(tmppath2)
    end if
    call chdirdbug('tmp_solu')

  end do

!> Frequency calculation
  opt = .true.
  call ens_freq(env,'cluster.xyz',solu_ens%nall,'TMPFREQ',opt)
  call chdirdbug(tmppath2)

!----------------------------------------------------------------------------
!   frequencies for solvent cluster
!----------------------------------------------------------------------------
  write (stdout,'(/,1x,a)') 'processing  SOLVENT CLUSTER'
!--- pure solvent clusters are neutral closed-shell
  chrg_tmp = env%chrg
  uhf_tmp = env%uhf
  env%chrg = 0
  env%uhf = 0
  if (env%cff) then
    call chdirdbug('tmp_solv')
    call ens_freq(env,'solvent_cut.coord',solu_ens%nall,'TMPFREQ',opt)
    call chdirdbug(tmppath2)
  end if

  call clus%deallocate()

  !--- Frequencies solvent cluster (only, if not cff was used)
  if (.not.env%cff) then
    call chdirdbug('tmp_solv')
    call solv_ens%write('solvent_ensemble.xyz')

    do i = 1,solv_ens%nall
      call solv_ens%get_mol(i,tmpmol)
      write (to,'("TMPFREQ",i0)') i
      io = makedir(trim(to))
      call copysub('.UHF',to)
      call copysub('.CHRG',to)
      call chdirdbug(to)
      call tmpmol%write("solv_cluster.xyz")

      call chdirdbug(tmppath2)
      call chdirdbug('tmp_solv')
    end do
!> Frequency calculation
    call ens_freq(env,'solv_cluster.xyz',solv_ens%nall,'TMPFREQ',opt)
    call chdirdbug(tmppath2)
  end if
  env%chrg = chrg_tmp
  env%uhf = uhf_tmp

!----------------------------------------------------------------------------
!   Data read out
!----------------------------------------------------------------------------

!--- Solute in gas phase
  write (stdout,*)
  write (stdout,*) '  Solute Gas properties'
  call pr_freq_energy()
  open (newunit=ich56,file='solute.dat')
  call pr_freq_file(ich56)
  write (stdout,'(2x,5f10.2)') ht(3),svib(3),srot(3),stra(3),gt(3)
  write (ich56,'(2x,5f10.2)') ht(3),svib(3),srot(3),stra(3),gt(3)
  close (ich56)

!--- Solute cluster
  write (stdout,*)
  write (stdout,*) '  Solute cluster properties'
  open (newunit=ich33,file='solute_cluster.dat')

  call chdirdbug('tmp_solu')

  allocate (solu_ens%gt(solu_ens%nall))
  allocate (solu_ens%ht(solu_ens%nall))
  allocate (solu_ens%svib(solu_ens%nall))
  allocate (solu_ens%srot(solu_ens%nall))
  allocate (solu_ens%stra(solu_ens%nall))

  call pr_freq_energy()
  call pr_freq_file(ich33)

  do i = 1,solu_ens%nall
    write (to,'("TMPFREQ",i0)') i
    call chdirdbug(to)
    call rdtherm('xtb_freq.out',ht(1),svib(1),srot(1),stra(1),gt(1))
    write (stdout,'(2x,i0,2x,5f10.2)') i,ht(1),svib(1),srot(1),stra(1),gt(1)
    write (ich33,'(2x,i0,2x,5f10.2)') i,ht(1),svib(1),srot(1),stra(1),gt(1)
    solu_ens%gt(i) = gt(1)
    solu_ens%ht(i) = ht(1)
    solu_ens%svib(i) = svib(1)
    solu_ens%srot(i) = srot(1)
    solu_ens%stra(i) = stra(1)

    call chdirdbug(tmppath2)
    call chdirdbug('tmp_solu')
  end do
  close (ich33)

!--- Solvent cluster
  write (stdout,*)
  write (stdout,*) '  Solvent cluster properties'
  call chdirdbug(tmppath2)
  open (newunit=ich81,file='solvent_cluster.dat')

  call chdirdbug('tmp_solv')

  allocate (solv_ens%gt(solv_ens%nall))
  allocate (solv_ens%ht(solv_ens%nall))
  allocate (solv_ens%svib(solv_ens%nall))
  allocate (solv_ens%srot(solv_ens%nall))
  allocate (solv_ens%stra(solv_ens%nall))

  call pr_freq_energy()
  call pr_freq_file(ich81)

  do i = 1,solv_ens%nall
    write (to,'("TMPFREQ",i0)') i
    call chdirdbug(to)
    call rdtherm('xtb_freq.out',ht(2),svib(2),srot(2),stra(2),gt(2))
    write (stdout,'(2x,i0,2x,5f10.2)') i,ht(2),svib(2),srot(2),stra(2),gt(2)
    write (ich81,'(2x,i0,2x,5f10.2)') i,ht(2),svib(2),srot(2),stra(2),gt(2)
    solv_ens%gt(i) = gt(2)
    solv_ens%ht(i) = ht(2)
    solv_ens%svib(i) = svib(2)
    solv_ens%srot(i) = srot(2)
    solv_ens%stra(i) = stra(2)
    call chdirdbug(tmppath2)
    call chdirdbug('tmp_solv')
  end do
  close (ich81)

!--- Saving results
  call chdirdbug(tmppath2)
  call copysub('solute.dat',resultspath)
  call copysub('solute_cluster.dat',resultspath)
  call copysub('solvent_cluster.dat',resultspath)

!--- Deleting tmp directory
  call chdirdbug(tmppath)
  if (.not.env%keepModef) call rmrf(tmppath2)
  call chdirdbug(thispath)

  env%gfnver = gfnver_tmp
  env%optlev = optlev_tmp

  call tim%stop(9)

end subroutine qcg_freq

!==============================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!==============================================================================!

subroutine qcg_eval(env,solu,solu_ens,solv_ens)
  use crest_parameters
  use crest_data
  use qcg_printouts
  use iomod
  use qcg_coord_type
  use strucrd
  use qcg_utils
  implicit none

  type(systemdata)   :: env
  type(coord_qcg)    :: solu
  type(ensemble)     :: solu_ens,solv_ens

  character(len=512) :: thispath

  integer            :: i,j
  integer            :: srange
  integer            :: freqscal
  real(wp)           :: g1(solu_ens%nall)
  real(wp)           :: g2(solv_ens%nall)
  real(wp)           :: g3
  real(wp)           :: Gsolv(20)
  real(wp)           :: Hsolv
  real(wp)           :: G_solute(20)
  real(wp)           :: H_solute
  real(wp)           :: G_solvent(20)
  real(wp)           :: H_solvent
  real(wp)           :: G_mono(20)
  real(wp)           :: H_mono
  real(wp)           :: S(20)
  real(wp)           :: volw
  real(wp)           :: sasa
  real(wp)           :: dum,dum1,dum2
  real(wp)           :: e_solute(solu_ens%nall)
  real(wp)           :: e_solvent(solv_ens%nall)
  real(wp)           :: scal(20)
  integer            :: ich23

  call pr_eval_eval()
  call qcg_check_reference(env,solu)

  call getcwd(thispath)

  freqscal = nint(env%freq_scal/0.05)
  srange = 20
  do i = 1,srange
    scal(i) = 0.05*i
  end do

!--- Solute Cluster
  !H_solv
  do i = 1,solu_ens%nall
    e_solute(i) = solu_ens%er(i)*autokcal+solu_ens%ht(i)
  end do
  call aver(.false.,env,solu_ens%nall,e_solute,dum1,H_solute,dum2,sasa,.false.)
  !G_solv
  do i = 1,srange
    do j = 1,solu_ens%nall
      g1(j) = solu_ens%ht(j)-(env%tboltz*(solu_ens%svib(j)+scal(i)*(solu_ens%srot(j)+solu_ens%stra(j)))/1000)
      e_solute(j) = solu_ens%er(j)*autokcal+g1(j)
    end do
    call aver(.false.,env,solu_ens%nall,e_solute,S(i),dum,G_solute(i),sasa,.false.)
  end do

!--- Solvent Cluster
  !H_solv
  do i = 1,solv_ens%nall
    e_solvent(i) = solv_ens%er(i)*autokcal+solv_ens%ht(i)
  end do
  call aver(.false.,env,solv_ens%nall,e_solvent,dum1,H_solvent,dum2,sasa,.false.)

  !G_solv
  do i = 1,srange
    do j = 1,solv_ens%nall
      g2(j) = solv_ens%ht(j)- &
              & (env%tboltz*(solv_ens%svib(j)+scal(i)*(solv_ens%srot(j)+solv_ens%stra(j)))/1000)
      e_solvent(j) = solv_ens%er(j)*autokcal+g2(j)
    end do
    call aver(.false.,env,solv_ens%nall,e_solvent,S(i),dum,G_solvent(i),sasa,.false.)
  end do

!--- Solute gas phase
  H_mono = solu%energy*autokcal+solu%ht
  do i = 1,srange
    g3 = solu%ht-(env%tboltz*(solu%svib+scal(i)*(solu%srot+solu%stra))/1000)
    G_mono(i) = solu%energy*autokcal+g3
  end do

  Gsolv(1:20) = G_solute(1:20)-G_solvent(1:20)-G_mono(1:20)
  Hsolv = H_solute-H_solvent-H_mono

!--- Calculate Volume work and include
  volw = (env%tboltz*8.31451/1000./4.184)*log(24.47d0*env%tboltz/298.15)
  Gsolv(1:20) = Gsolv(1:20)-volw
  Hsolv = Hsolv-volw
  call pr_eval_1(Gsolv(20),Hsolv)
  call pr_eval_2(srange,Gsolv,scal)
  call pr_eval_3(srange,freqscal,env%freq_scal,Gsolv)

! Save Result
  open (newunit=ich23,file='frequencies/result.dat')
  write (ich23,'("Solvation Free Energy [kcal/mol] :")')
  write (ich23,'(f12.2)') Gsolv(freqscal)
  close (ich23)

end subroutine qcg_eval

!==============================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!==============================================================================!

subroutine qcg_restart(env,progress,solu,solv,clus,solu_ens,solv_ens,clus_backup)
  use crest_parameters
  use crest_data
  use iomod
  use qcg_coord_type
  use strucrd
  use qcg_utils
  implicit none

  type(systemdata)           :: env
  type(coord_qcg)            :: solu,solv,clus,clus_backup
  type(ensemble)             :: solu_ens,solv_ens
  integer                    :: progress

  integer                    :: i
  character(len=512)         :: thispath
  character(len=6)           :: counter
  character(len=7)           :: counter2
  character(len=8)           :: counter3
  logical                    :: grow,solu_ensemble,solv_ensemble
  logical                    :: solv_cff,solv_present,freq,tmp,ex
  real(wp),allocatable       :: xyz(:,:)

  grow = .false.
  solu_ensemble = .false.
  solv_ensemble = .false.
  solv_cff = .false.
  solv_present = .false.
  freq = .false.
  tmp = .false.

  inquire (file='./grow/cluster.coord',exist=grow)
  inquire (file='./ensemble/final_ensemble.xyz',exist=solu_ensemble)
  inquire (file='./solvent_ensemble/final_ensemble.xyz',exist=solv_ensemble)
  inquire (file='./solvent_ensemble/crest_ensemble.xyz',exist=solv_cff)
  inquire (file='./frequencies/result.dat',exist=freq)

  if (solv_cff.or.solv_ensemble) solv_present = .true.

  call getcwd(thispath)

!---------------------------------------------------------------------------------
!        Check, if everything needed is present
!---------------------------------------------------------------------------------

  if (freq.and.((.not.grow).or.(.not.solu_ensemble).or.(.not.solv_ensemble))) then
    progress = 0
    call rmrf('frequencies')
    freq = .false.
  end if

  if (solv_present.and.((.not.grow).or.(.not.solu_ensemble))) then
    progress = 0
    call rmrf('solvent_ensemble')
    solv_present = .false.
    solv_cff = .false.
    solv_ensemble = .false.
  end if

  if (solu_ensemble.and.(.not.grow)) then
    progress = 0
    call rmrf('ensemble')
    solu_ensemble = .false.
  end if

!-------------------------------------------------------------
!           Data read out
!-------------------------------------------------------------

!--- Grow process
  if (grow) then
    env%qcg_restart = .true.
    call chdirdbug('grow')
    call clus%open('cluster.coord')

    clus%nmol = (clus%nat-solu%nat)/solv%nat+1
    allocate (xyz(3,clus%nat))
    xyz = clus%xyz
    call get_ellipsoid(env,solu,solv,clus,.true.)
    clus%xyz = xyz !Needed, because get_ellipsoid performs axistransformation and not fitting potential
    deallocate (xyz)

    if (.not.env%cff) then
      allocate (clus_backup%at(clus%nat))
      allocate (clus_backup%xyz(3,clus%nat))
      clus_backup = clus
    end if

    if (clus%nmol-1 .ge. env%nsolv) then
      progress = 1
      env%nsolv = clus%nmol-1
      write (stdout,*)
      write (stdout,*)
      write (stdout,'(''Found cluster with '',i0,'' solvents'')') env%nsolv
      call chdirdbug(thispath)
    else
      error stop 'The found cluster is smaller than nsolv. Please restart the whole computaion by removing the grow directory'
      !Future implementation continue grow process
      call chdirdbug(thispath)
      if (solu_ensemble) call rmrf('ensemble')
      if (solv_ensemble) call rmrf('solvent_ensemble')
      if (freq) call rmrf('frequencies')
      solu_ensemble = .false.
      solv_ensemble = .false.
      freq = .false.
      progress = 0
    end if
  end if

!--- Solute Ensemble
  if (solu_ensemble) then
    call chdirdbug('ensemble')
    call solu_ens%open('final_ensemble.xyz')
    call rdensemble('final_ensemble.xyz',solu_ens%nat,solu_ens%nall,solu_ens%at,solu_ens%xyz,solu_ens%er)
    env%nqcgclust = solu_ens%nall
    write (stdout,'("  Ensemble of solute-cluster found.")')
    write (stdout,'("  Taking all ", i0, " structures")') env%nqcgclust
    call grepval('population.dat','Ensemble free energy [Eh]:',ex,solu_ens%G)
    solu_ens%G = solu_ens%G*autokcal
    write (stdout,*) 'Solute Ensmeble Free E [kcal/mol]',solu_ens%G
    call chdirdbug(thispath)
    progress = 2
  end if

!--- Solvent Ensemble
  if (solv_present) then
    call chdirdbug('solvent_ensemble')
    write (stdout,'("  Ensemble of solvent-cluster found.")')

    !--- Case CFF
    if (solv_cff) then
      call solv_ens%open('crest_ensemble.xyz')
      do i = 1,solv_ens%nall
        if (i .le. 9) then
          write (counter,'(''No   '',i1)') i
          call grepval('cluster_energy.dat',counter,ex,solv_ens%er(i))
        else if (i .le. 99) then
          write (counter2,'(''No   '',i2)') i
          call grepval('cluster_energy.dat',counter2,ex,solv_ens%er(i))
        else
          write (counter3,'(''No   '',i3)') i
          call grepval('cluster_energy.dat',counter3,ex,solv_ens%er(i))
        end if
        write (stdout,*) 'Energy of cluster',i,solv_ens%er(i)
      end do
    end if

    !--- Case MD/Crest run
    if (solv_ensemble) then
      call solv_ens%open('final_ensemble.xyz')
      call rdensemble('final_ensemble.xyz',solv_ens%nat,solv_ens%nall,solv_ens%at,solv_ens%xyz,solv_ens%er)
    end if
    call grepval('population.dat','Ensemble free energy [Eh]:',ex,solv_ens%G)
    solv_ens%G = solv_ens%G*autokcal
    write (stdout,*) 'solvent ensmeble free E [kcal/mol]',solv_ens%G
    call chdirdbug(thispath)
    progress = 3
  end if

!--- Frequencies
  if (freq) then
    write (stdout,*)
    write (stdout,*)
    write (stdout,*) '  Nothing to do'
    progress = 4
  end if

end subroutine qcg_restart

!==============================================================================!
!
subroutine qcg_cleanup(env)
  use crest_data
  implicit none
  type(systemdata)      :: env
  character(len=280)    :: thispath
  logical               :: tmp
  call getcwd(thispath)
  call chdirdbug(env%scratchdir)
  inquire (file='./solute_properties/solute',exist=tmp)
  if (tmp) then
    call rmrf('solute_properties')
    call rmrf('solvent_properties')
  end if
end subroutine qcg_cleanup

!==============================================================================!
