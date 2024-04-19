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
!===================================================================!
! This file contains routines related to QCG and microsolvation
!===================================================================!
!======================================================!
! main routine
!======================================================!
subroutine crest_solvtool(env, tim)
   use iso_fortran_env, wp => real64
   use crest_data
   use iomod
   use zdata
   use strucrd
   implicit none

   type(systemdata):: env    ! MAIN STORAGE OS SYSTEM DATA
   type(timer):: tim
   !> Information about solvent, solute and cluster
   type(zmolecule) :: solute, solvent, cluster, cluster_backup 
   type(ensemble) :: full_ensemble, solvent_ensemble

   integer :: progress,io
   character(len=512) :: thispath

   real(wp), parameter         :: eh = 627.509541d0

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
   if(io /= 0 ) error stop 'No xtb found'

!> Check, if xtbiff is present (if it is required)
   if (env%use_xtbiff) then
      call checkprog_silent(env%ProgIFF,.true.,iostat=io)
      if(io /= 0 ) error stop 'No xtbiff found'
   else
      write (*, *)
      write (*, *) '  The use of the aISS algorithm is requested (recommend).'
      write (*, *) '  This requires xtb version 6.6.0 or newer.'
      write (*, *) '  xTB-IFF can still be used with the --xtbiff flag.'
      write (*, *)
   end if

!------------------------------------------------------------------------------
!   Setup
!------------------------------------------------------------------------------

   call write_qcg_setup(env) !Just an outprint of setup
   call read_qcg_input(env, solute, solvent) !Reading mol. data and determining r,V,A
   call qcg_setup(env, solute, solvent)
   call qcg_restart(env, progress, solute, solvent, cluster, full_ensemble,&
          & solvent_ensemble, cluster_backup)

!-----------------------------------------------------------------------------
!   Grow
!-----------------------------------------------------------------------------
   if (progress .le. env%qcg_runtype .and. progress .eq. 0) then
      cluster = solute
      call qcg_grow(env, solute, solvent, cluster, tim)
      if (.not. env%cff) then
         allocate (cluster_backup%at(cluster%nat))
         allocate (cluster_backup%xyz(3, cluster%nat))
         cluster_backup = cluster
      end if
      progress = progress + 1
      call chdir(thispath)
   end if

!------------------------------------------------------------------------------
!   Ensemble search
!------------------------------------------------------------------------------
   if (progress .le. env%qcg_runtype .and. progress .eq. 1) then
      call print_qcg_ensemble()
      call qcg_ensemble(env, solute, solvent, cluster, full_ensemble, tim, 'ensemble')
      progress = progress + 1
      call chdir(thispath)
   end if

!------------------------------------------------------------------------------
!   Solvent cluster generation
!------------------------------------------------------------------------------
   if (progress .le. env%qcg_runtype .and. progress .eq. 2) then !esolv
      call pr_eval_solvent()
      if (env%cff) then !CFF
         call qcg_cff(env, solute, solvent, cluster, full_ensemble,&
                & solvent_ensemble, tim)
      else !Normal ensemble generation
         call print_qcg_ensemble()
         call cluster%deallocate
         allocate (cluster%at(cluster_backup%nat))
         allocate (cluster%xyz(3, cluster_backup%nat))
         cluster = cluster_backup
         deallocate (cluster_backup%at)
         deallocate (cluster_backup%xyz)
         env%solv_md = .true.
         call qcg_ensemble(env, solute, solvent, cluster, solvent_ensemble,&
                & tim, 'solvent_ensemble')
      end if
      call pr_qcg_esolv()
      write (*, '(2x,"|",9x,F8.2," kcal/mol ",12x,"|")') &
             &   full_ensemble%g - solvent_ensemble%g - (solute%energy*eh)
      write (*, '(2x,''========================================='')')
      call chdir(thispath)
      progress = progress + 1
   end if

!------------------------------------------------------------------------------
!   Frequency computation and evaluation
!------------------------------------------------------------------------------
   if (progress .le. env%qcg_runtype .and. progress .eq. 3) then !gsolv
      call qcg_freq(env, tim, solute, solvent, full_ensemble, solvent_ensemble)
      call qcg_eval(env, solute, full_ensemble, solvent_ensemble)

      progress = progress + 1
   end if

   !<----------------------------------
!  call tim%stop(2) !stop a timer

!------------------------------------------------------------------------------
!   Cleanup and deallocation
!------------------------------------------------------------------------------
   if (env%scratchdir .ne. 'qcg_tmp') call qcg_cleanup(env)
   if (.not. env%keepModef) call rmrf('qcg_tmp')
   call solute%deallocate
   call solvent%deallocate
   call cluster%deallocate
   call full_ensemble%deallocate
   call solvent_ensemble%deallocate
   return
end subroutine crest_solvtool

subroutine qcg_setup(env, solu, solv)
   use iso_fortran_env, wp => real64
   use crest_data
   use iomod
   use zdata
   use strucrd
   use axis_module
   implicit none

   type(systemdata):: env
   type(zmolecule) :: solv, solu

   integer :: io, f, r
   integer :: num_O, num_H, i
   character(len=*), parameter :: outfmt = '(1x,1x,a,1x,f14.7,a,1x)'
   logical :: e_there, tmp, used_tmp
   character(len=512) :: thispath, tmp_grow
   character(len=40)  :: solv_tmp
   character(len=80)  :: atmp
   character(len=20)  :: gfnver_tmp

   call getcwd(thispath)

   ! Remove scratch dir, if present
   inquire (file='./qcg_tmp/solute_properties/solute', exist=tmp)
   if (tmp) call rmrf('qcg_tmp') !User given scratch dir will be removed anyway after run

   ! Make scratch directories
   if (env%scratchdir .eq. '') then !check if scratch was not set
      env%scratchdir = 'qcg_tmp'
      io = makedir('qcg_tmp')
   end if
   if(env%fixfile /= 'none selected') then
      call copysub(env%fixfile, env%scratchdir)
   end if
   call chdir(env%scratchdir)

   f = makedir('solute_properties')
   if(env%fixfile /= 'none selected') then
      call copysub(env%fixfile, env%scratchdir)
   end if
   r = makedir('solvent_properties')

   if (.not. env%nopreopt) then
      write (*, *)
      write (*, '(2x,''========================================='')')
      write (*, '(2x,''|            Preoptimization            |'')')
      write (*, '(2x,''========================================='')')
   end if

   solv_tmp = env%solv
   env%solv = ''

!---- Properties solute
   call chdir('solute_properties')
   call env%wrtCHRG('') !Write three lines in QCG mode, but xtb anyway only reads first one

!---- Geometry preoptimization solute
   if (env%final_gfn2_opt) then !If GFN2 final opt, solute also GFN2 optimized
      gfnver_tmp = env%gfnver
      env%gfnver = '--gfn2'
   end if

   if ((.not. env%nopreopt) .and. (solu%nat /= 1)) then
      call xtb_opt_qcg(env, solu, .true.)
   end if

!--- Axistrf
   call axistrf(solu%nat, solu%nat, solu%at, solu%xyz)
   call wrc0('solute', solu%nat, solu%at, solu%xyz)

!---- LMO/SP-Computation solute
   if (env%use_xtbiff) then
      write (*, *) 'Generating LMOs for solute'
      call xtb_lmo(env, 'solute')
   else
      call xtb_sp_qcg(env, 'solute')
   end if

   if (env%final_gfn2_opt) then !If GFN2 final opt, solute also GFN2 optimized
      env%gfnver = gfnver_tmp
   end if

   call grepval('xtb.out', '| TOTAL ENERGY', e_there, solu%energy)
   if (.not. e_there) then
      write (*, *) 'Total Energy of solute not found'
   else
      write (*, outfmt) 'Total Energy of solute: ', solu%energy, ' Eh'
   end if

   if (env%use_xtbiff) then
      call rename('xtblmoinfo', 'solute.lmo')
   end if

   call chdir(thispath)

! No constraints for solvent possible
   used_tmp = env%cts%used
   env%cts%used = .false.

!---- Properties solvent
   call chdir(env%scratchdir)
   call chdir('solvent_properties')
   !No charges for solvent written. This is currently not possible

!---- Geometry preoptimization solvent
   if ((.not. env%nopreopt) .and. (solv%nat /= 1)) then
      call xtb_opt_qcg(env, solv, .false.)
   end if
   call wrc0('solvent', solv%nat, solv%at, solv%xyz)

!---- LMO-Computation solvent
   if (env%use_xtbiff) then
      write (*, *) 'Generating LMOs for solvent'
      call xtb_lmo(env, 'solvent')!,solv%chrg)
   else
      call xtb_sp_qcg(env, 'solvent')
   end if

   call grepval('xtb.out', '| TOTAL ENERGY', e_there, solv%energy)
   if (.not. e_there) then
      write (*, '(1x,a)') 'Total Energy of solvent not found'
   else
      write (*, outfmt) 'Total energy of solvent:', solv%energy, ' Eh'
   end if

   if (env%use_xtbiff) then
      call rename('xtblmoinfo', 'solvent.lmo')
   end if

   call chdir(thispath)

!---- Overwriting solute and solvent in original folder
   call wrc0('solute', solu%nat, solu%at, solu%xyz)
   call wrc0('solvent', solv%nat, solv%at, solv%xyz)

   num_O = 0
   num_H = 0
!--- Check, if water is solvent
   if (solv%nat .eq. 3) then
      do i = 1, solv%nat
         if (solv%at(i) .eq. 8) num_O = num_O + 1
         if (solv%at(i) .eq. 1) num_H = num_H + 1
      end do
   end if
   if (num_O .eq. 1 .AND. num_H .eq. 2) then
      env%water = .true.
      if (.not. env%noconst) env%constrain_solu = .true.
   end if

   env%solv = solv_tmp
   env%cts%used = used_tmp

end subroutine qcg_setup

subroutine read_qcg_input(env, solu, solv)
   use iso_fortran_env, wp => real64
   use crest_data
   use iomod
   use zdata
   use strucrd
   use atmasses
   implicit none

   type(systemdata)               :: env
   type(zmolecule), intent(inout) :: solu, solv
   logical                        :: pr
   real(wp), parameter            :: amutokg = 1.66053886E-27
   real(wp), parameter            :: third = 1.0d0/3.0d0
   integer                        :: i
   real(wp)                       :: r_solu, r_solv

   pr = .true.

!--- Read in solu and solv coordinates and make solute and solvent file in WD
   call inputcoords_qcg(env, solu, solv)

!--- CMA-Trafo
   call cma_shifting(solu, solv)

!--- Setting solute charge and uhf to input
   solu%chrg = env%chrg
   solu%uhf = env%uhf

!--- Getting r, V, A
   write (*, *)
   write (*, *) 'Solute geometry'
   call get_sphere(.true., solu, .true.) !r,V,A of solute
   write (*, *) 'Solvent geometry'
   call get_sphere(.true., solv, .true.) !r,V,A of solvent

   r_solu = solu%vtot**third
   r_solv = solv%vtot**third
   write (*, *)
   write (*, '(2x,''radius of solute    : '',f8.2)') r_solu
   write (*, '(2x,''radius of solvent   : '',f8.2)') r_solv

!--- Determine masses (for later density computation)
   do i = 1, solu%nat
      solu%mass = solu%mass + ams(solu%at(i))
   end do
   do i = 1, solv%nat
      solv%mass = solv%mass + ams(solv%at(i))
   end do
   solu%mass = solu%mass*amutokg
   solv%mass = solv%mass*amutokg

!--- If directed docking is requested, it is read in here:
   if(allocated(env%directed_file)) then
      if (env%use_xtbiff) error stop 'xTB-IFF does not support directed docking. &
              &Please use the aISS algorithm of xtb.'
      call read_directed_input(env)
   end if

end subroutine read_qcg_input

!> Read input for directed docking
subroutine read_directed_input(env)
   use iso_fortran_env, wp => real64
   use crest_data
   implicit none

   type(systemdata)           :: env

   integer                    :: nlines
   integer                    :: io, ich, i, i_check
   integer                    :: index
   character(len=512)         :: dum
   character(len=1), parameter :: delim_space = ' ', delim_tab = achar(9)

     open (newunit=ich, file=env%directed_file) 
        !First check number of lines
        nlines = 0
        do
           read(ich,*,iostat=io)
           if (io /= 0) exit
           nlines = nlines + 1
        end do
        !Allocate directed list 
        !First entry is the atom number, Second how many solvents to add to this atom
        allocate(env%directed_list(nlines,2))
        allocate(env%directed_number(nlines), source = 0)
        !Now read lines into directed_list
        rewind(ich)
        do i=1, nlines
           read(ich,'(A)') dum
           !> Remove leading tab and spaces first
           dum = adjustl(dum) !Leading spaces are removed
           index = SCAN(trim(dum), delim_tab) 
           if (index == 1) then !Leading tab -> remove it
              dum = dum(2:)
           end if
           index = SCAN(trim(dum), delim_space)
           if (index == 0) then !No space = check for tab
              index = SCAN(trim(dum), delim_tab)
           end if
           if (index == 0) then  !Second value is missing
              write(*,'(a,1x,i0)') "No second value found in directed list on line", i
              error stop
           end if
           env%directed_list(i, 1) = dum(1:index-1)
           env%directed_list(i, 2) = dum(index+1:)
           !Remove multiple spaces
           env%directed_list(i, 2) = adjustl(env%directed_list(i, 2))
           !Check, if spaces are still in second argument (e.g. a third number is giveb)
           index = SCAN(trim(env%directed_list(i, 2)), delim_space)
           if (index == 0) index = SCAN(trim(dum), delim_tab)
           if (index /= 0) then
              write(*,'(a,1x,i0)') "Too many values at line", i
              error stop
           end if
           !> Make array with which solvent molecule at which atom to add
           read(env%directed_list(i,2), *, iostat=io) env%directed_number(i)
           env%directed_number(i) = sum(env%directed_number)
           if (io/= 0) then
              write(*,'(a,1x,i0)') "Second value is no number in line", i
              error stop
           end if
        end do
     close(ich)
     write(*,*) 'Performing directed docking'
     do i=1, nlines
        write(*,'(a,1x,a,1x,a,1x,a)') 'Docking', trim(env%directed_list(i,2)),&
               & 'solvent molecules at', trim(env%directed_list(i,1))
     end do

end subroutine read_directed_input


subroutine qcg_grow(env, solu, solv, clus, tim)
   use crest_parameters
   use crest_data
   use iomod
   use zdata
   use strucrd
   implicit none

   type(systemdata)           :: env
   type(zmolecule)            :: solu, solv, clus
   type(timer)                :: tim

   integer                    :: minE_pos, m
   integer                    :: iter = 1
   integer                    :: i, j, io, v
   integer                    :: max_cycle
   logical                    :: e_there, high_e, success, neg_E
   real(wp)                   :: etmp(500)
   real(wp), allocatable       :: e_each_cycle(:)
   real(wp)                   :: dens, dum, efix
   real(wp)                   :: e_diff = 0.0_wp
   real(wp), parameter         :: eh = 627.509541d0
   real(wp), allocatable       :: E_inter(:)
   real(wp)                   :: shr = 0.0_wp
   real(wp)                   :: shr_av = 0.0_wp
   real(wp)                   :: mean = 0.0_wp
   real(wp)                   :: mean_old = 0.0_wp
   real(wp)                   :: mean_diff = 0.0_wp
   character(len=*), parameter :: outfmt = '(1x,1x,a,1x,f14.7,a,1x)'
   character(len=512)         :: thispath, resultspath
   character(len=20)          :: gfnver_tmp
   integer                    :: ich99, ich15, ich88
   character(len=LEN(env%solv)) :: solv_tmp
   logical                    :: gbsa_tmp

   interface
      subroutine both_ellipsout(fname, n, at, xyz, r1, r2)
         use iso_fortran_env, only: wp => real64
         use strucrd, only: i2e
         implicit none

         integer            :: i, j
         integer            :: n, at(n)
         real(wp)           :: dum(3)
         real(wp)           :: rx, ry, rz
         real(wp)           :: xyz(3, n), r1(3)
         real(wp), optional :: r2(3)
         real               :: x, y, z, f, rr
         character(len=*)   :: fname
         integer            :: ich11
      end subroutine both_ellipsout
   end interface

   if (env%nsolv .gt. 0) then
      allocate (e_each_cycle(env%nsolv))
      allocate (E_inter(env%nsolv))
   else
      allocate (e_each_cycle(env%max_solv))
      allocate (E_inter(env%max_solv))
   end if

   call tim%start(5, 'Grow')

   call pr_eval_solute()
   call print_qcg_grow()
   call getcwd(thispath)
   io = makedir('grow')
   call chdir('grow') !Results directory

!--- Output Files
   open (newunit=ich99, file='qcg_energy.dat')
   write (ich99, '(i0,2F20.8)') 0, solu%energy, solv%energy
   open (newunit=ich15, file='qcg_grow.xyz')  ! for molden movie
   open (newunit=ich88, file='qcg_conv.dat')  ! for convergence check
   write (ich88, '(''   #   Energy       Run. Aver.   Diff / au.'')')

   call getcwd(resultspath)
   call chdir(thispath)

   if (env%water) then
      if (.not. env%user_wscal) then
         if (solu%nat .lt. 18) then
            env%potscal = 0.7_wp
         else
            env%potscal = 0.8_wp
         end if
         write (*, *)
         write (*, '(2x,''Water as solvent recognized, &
                 & adjusting scaling factor for outer wall pot to '',F4.2)')&
                & env%potscal
         write (*, *)
      end if
   end if
   if (env%constrain_solu) write (*, '(2x,''Constraining solute during Growth '')')

   call get_ellipsoid(env, solu, solv, clus, .true.)
   call pr_grow_energy()

   call chdir(env%scratchdir)
   v = makedir('tmp_grow')
   if(env%fixfile /= 'none selected') then
      call copysub(env%fixfile, 'tmp_grow')
   end if
   if (env%use_xtbiff) then
     call copy('solute_properties/solute.lmo', 'tmp_grow/solute.lmo')
     call copy('solvent_properties/solvent.lmo', 'tmp_grow/solvent.lmo')
   end if
   call chdir('tmp_grow')
   call wrc0('solute', solu%nat, solu%at, solu%xyz)
   call wrc0('solvent', solv%nat, solv%at, solv%xyz)

   call ellipsout('solute_cavity.coord', clus%nat, clus%at, clus%xyz, solu%ell_abc)
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
   do iter = 1, max_cycle

      e_there = .false.
      success = .false.
      high_e = .false.
      neg_E = .false.
!---- LMO-Computation
      if (iter .gt. 1) then
         call get_ellipsoid(env, solu, solv, clus, .false.)
         if (env%use_xtbiff) then
            call xtb_lmo(env, 'xtbopt.coord')!,clus%chrg)
            call grepval('xtb.out', '| TOTAL ENERGY', e_there, clus%energy)
            if (.not. e_there) then
               write (*, '(1x,a)') 'Total Energy of cluster LMO computation not found'
            end if
            call rename('xtblmoinfo', 'cluster.lmo')
         end if
      end if

      call both_ellipsout('twopot_1.coord', clus%nat, clus%at, clus%xyz,&
             & clus%ell_abc, solu%ell_abc)

      do while (.not. success) !For restart with larger wall pot
         if (iter .eq. 1) then
            if (env%use_xtbiff) then
               call xtb_iff(env, 'solute.lmo', 'solvent.lmo', solu, solv)
               !solu for nat of core pot. solv for outer ellips
               call check_iff(neg_E)
            else
               call xtb_dock(env, 'solute', 'solvent', solu, solv)
               call check_dock(neg_E)
            end if

!-- If Interaction Energy is not negativ and existent, wall pot. too small and increase
            if (neg_E) then
               success = .true.
            else
               if (env%potscal .lt. 1.0_wp) then
                  write (*, *) '  Wall Potential too small, increasing size by 5 %'
                  solv%ell_abc = solv%ell_abc*1.05_wp
                  env%potscal = env%potscal*1.05_wp
                  if (env%potscal .gt. 1.0_wp) env%potscal = 1.0_wp
                  write (*, '(''   New scaling factor '',F4.2)') env%potscal
               else
                  success = .true.
               end if
            end if
         else
            if (env%use_xtbiff) then
               call xtb_iff(env, 'cluster.lmo', 'solvent.lmo', solu, clus)
               call check_iff(neg_E)
            else
               call xtb_dock(env, 'cluster.coord', 'solvent', solu, clus)
               call check_dock(neg_E)
            end if

            if (neg_E) then
               success = .true.
            else
               if (env%potscal .lt. 1.0_wp) then
                  write (*, *) '  Wall Potential too small, increasing size by 5 %'
                  clus%ell_abc = clus%ell_abc*1.05_wp
                  env%potscal = env%potscal*1.05_wp
                  if (env%potscal .gt. 1.0_wp) env%potscal = 1.0_wp
                  write (*, '(''   New scaling factor '',F4.2)') env%potscal
               else
                  success = .true.
               end if
            end if
         end if
      end do

!--- Increase cluster size
      call clus%deallocate
      clus%nat = clus%nat + solv%nat
      allocate (clus%at(clus%nat))
      allocate (clus%xyz(3, clus%nat))
      clus%nmol = clus%nmol + 1

!--- Select xtb-IFF stucture to proceed
      if (env%use_xtbiff) then
         call rdxtbiffE('xtbscreen.xyz', m, clus%nat, etmp) !Get energy of screening
         minE_pos = minloc(etmp(1:m), dim=1) !Get minimum of those
         !Read the struc into clus%xyz
         call rdxmolselec('xtbscreen.xyz', minE_pos, clus%nat, clus%at, clus%xyz) 
      else
         call rdcoord('best.xyz', clus%nat, clus%at, clus%xyz, clus%energy)
      end if

      call remove('cluster.coord')
      call wrc0('cluster.coord', clus%nat, clus%at, clus%xyz)
      call both_ellipsout('twopot_2.coord', clus%nat, clus%at, clus%xyz,&
             & clus%ell_abc, solu%ell_abc)

      success = .false.

!--- Cluster restart, if interaction energy not negativ (wall pot. too small)
      do while (.not. success)
!--- Cluster optimization
         if (env%cts%used) then
            call write_reference(env, solu, clus) !new fixed file
         end if

         if (env%use_xtbiff) then
            call opt_cluster(env, solu, clus, 'cluster.coord', .false.)
            call rdcoord('xtbopt.coord', clus%nat, clus%at, clus%xyz)
         end if

!--- Interaction energy
         gfnver_tmp = env%gfnver
         env%gfnver = env%lmover
         gbsa_tmp = env%gbsa
         solv_tmp = env%solv
         env%gbsa = .false.
         env%solv = ''
         call get_interaction_E(env, solu, solv, clus, iter, E_inter)
         env%gbsa = gbsa_tmp
         env%solv = solv_tmp
         if (E_inter(iter) .lt. 0) then
            success = .true.
         else
            if (env%potscal .lt. 1.0_wp) then
               write (*, *) '  Interaction Energy positiv, increasing outer wall pot by 5 %'
               clus%ell_abc = clus%ell_abc*1.05_wp
               env%potscal = env%potscal*1.05_wp
               if (env%potscal .gt. 1.0_wp) env%potscal = 1.0_wp
               write (*, '('' New scaling factor '',F4.2)') env%potscal
            else
               success = .true.
            end if
         end if
      end do
      env%gfnver = gfnver_tmp

!--- For output
      if (env%use_xtbiff) then
         call grepval('xtb.out', '| TOTAL ENERGY', e_there, clus%energy)
         call wrc0('optimized_cluster.coord', clus%nat, clus%at, clus%xyz)
         if (.not. e_there) then
            write (*, '(1x,a)') 'Total Energy of cluster not found.'
         end if
      else
         !Energy already read from xyz file
      end if
      e_each_cycle(iter) = clus%energy

!--- Calclulate fix energy + diff. energy
      efix = clus%energy/sqrt(float(clus%nat))
      dum = solu%energy
      if (iter .gt. 1) dum = e_each_cycle(iter - 1)
      e_diff = e_diff + eh*(e_each_cycle(iter) - solv%energy - dum)
      call ellipsout('cluster_cavity.coord', clus%nat, clus%at, clus%xyz, clus%ell_abc)
      call both_ellipsout('twopot_cavity.coord', clus%nat, clus%at, clus%xyz,&
             & clus%ell_abc, solu%ell_abc)

!--- Density calculations
      call get_sphere(.false., clus, .false.) !V, A of new cluster
      dens = 0.001*(solu%mass + iter*solv%mass)/(1.0d-30*clus%vtot*bohr**3)

!--- Movie file
      write (ich15, *) clus%nat
      write (ich15, '('' SCF done '',2F16.8)') eh*(e_each_cycle(iter) - solv%energy - dum)
      do j = 1, clus%nat
         write (ich15, '(a,1x,3F24.10)') i2e(clus%at(j)), clus%xyz(1:3, j)*bohr
      end do

!--- Output
     ! dist of new mol from solute for output
      call analyze_cluster(iter, clus%nat, solu%nat, solv%nat, clus%xyz, clus%at, shr_av, shr)

      write (*, '(x,i4,F13.6,1x,f7.2,3x,f8.2,6x,f6.3,3x,f8.3,3x,2f6.1,2x,f8.1,3x,a,x)') &
            & iter, e_each_cycle(iter), eh*(e_each_cycle(iter) - solv%energy - dum),&
            & e_diff, dens, efix, shr_av, shr, clus%vtot, trim(optlevflag(env%optlev))
      write (ich99, '(i4,F20.10,3x,f8.1)') iter, e_each_cycle(iter), clus%vtot

!--- Calculate moving average
      mean_old = mean
      do i = 0, iter - 1
         mean = mean + E_inter(iter - i)
      end do
      mean = mean/iter
      mean_diff = mean - mean_old
      write (ich88, '(i5,1x,3F13.8)') iter, E_inter(iter)*eh, mean, mean_diff

!--- Check if converged when no nsolv was given
      if (env%nsolv .eq. 0) then
         if (abs(mean_diff) .lt. 1.0d-4 .and. iter .gt. 5) then
            env%nsolv = iter
            exit
         end if
         if (iter .eq. env%max_solv) then
            write (*, '(1x,''No convergence could be reached upon adding'',1x,i4,1x,&
                    & ''solvent molecules.'')') env%max_solv
            write (*, *) ' Proceeding.'
            env%nsolv = env%max_solv
            exit
         end if
      end if
!-----------------------------------------------
! End loop
!-----------------------------------------------
   end do

   if (env%nsolv .eq. 0) env%nsolv = iter !if no env%solv was given

   if (env%gfnver .ne. '--gfn2' .and. env%final_gfn2_opt) then
      gfnver_tmp = env%gfnver
      env%gfnver = '--gfn2'
      write (*, '(2x,''Final gfn2 optimization'')')
      call opt_cluster(env, solu, clus, 'cluster.coord', .false.)
      call rdcoord('xtbopt.coord', clus%nat, clus%at, clus%xyz)
      call wrc0('cluster.coord', clus%nat, clus%at, clus%xyz)
      call grepval('xtb_sp.out', '| TOTAL ENERGY', e_there, clus%energy)
      if (.not. e_there) then
         write (*, '(1x,a)') 'Total Energy of cluster not found.'
      else
         write (*, '(2x,''Total gfn2-energy of cluster/Eh:'',f20.6)') clus%energy
      end if
      env%gfnver = gfnver_tmp
   end if

   call wrxyz('cluster.xyz', clus%nat, clus%at, clus%xyz*bohr)

!--- One optimization without Wall Potential and with implicit model
   gfnver_tmp = env%gfnver
   if (env%final_gfn2_opt) env%gfnver = '--gfn2'
   call opt_cluster(env, solu, clus, 'cluster.xyz', .true.)
   env%gfnver = gfnver_tmp
   call rename('xtbopt.xyz', 'cluster_optimized.xyz')
   call copysub('cluster_optimized.xyz', resultspath)

!--- output and files
   write (*, *)
   write (*, '(2x,''Growth finished after '',i0,'' solvents added'')') env%nsolv
   write (*, '(2x,''Results can be found in grow directory'')')
   write (*, '(2x,''Energy list in file <qcg_energy.dat>'')')
   write (*, '(2x,''Interaction energy in file <qcg_conv.dat>'')')
   write (*, '(2x,''Growing process in <qcg_grow.xyz>'')')
   write (*, '(2x,''Final geometry after grow in <cluster.coord> and <cluster.xyz>'')')
   write (*, '(2x,''Final geometry optimized without wall potential in <cluster_optimized.xyz>'')')
   write (*, '(2x,''Potentials and geometry written in <cluster_cavity.coord> and <twopot_cavity.coord>'')')

   close (ich99)
   close (ich88)
   close (ich15)

!--- Saving results and cleanup
   call copysub('cluster.coord', resultspath)
   call copysub('cluster.xyz', resultspath)
   call copysub('twopot_cavity.coord', resultspath)
   call copysub('cluster_cavity.coord', resultspath)
   call copysub('solute_cavity.coord', resultspath)
!  call rename('xcontrol','wall_potential')
   env%constrain_solu = .false.
   call write_wall(env, solu%nat, solu%ell_abc, clus%ell_abc, 'wall_potential')
   call copysub('wall_potential', resultspath)

   call chdir(thispath)
   call chdir(env%scratchdir)
   if (.not. env%keepModef) call rmrf('tmp_grow')

   deallocate (e_each_cycle, E_inter)

   call tim%stop(5)

end subroutine qcg_grow

subroutine qcg_ensemble(env, solu, solv, clus, ens, tim, fname_results)
   use crest_parameters
   use crest_data
   use iomod
   use zdata
   use strucrd
   use utilities
   use cregen_interface
   implicit none

   type(systemdata)           :: env
   type(zmolecule)            :: solu, solv, clus
   type(ensemble)             :: ens, dum
   type(timer)                :: tim

   integer                    :: i, j, k
   integer                    :: io, f, r, ich,T,Tn
   integer                    :: minpos
   character(len=512)         :: thispath, resultspath, tmppath, tmppath2
   character(len=512)         :: scratchdir_tmp
   character(len=512)         :: jobcall
   character(len=256)         :: inpnam, outnam
   character(len=80)          :: fname, pipe, to
   character(len=*)           :: fname_results
   character(len=64)          :: comment
   character(len=20)          :: gfnver_tmp
   character(len=LEN(env%solv)) :: solv_tmp
   logical                    :: gbsa_tmp
   logical                    :: ex, mdfail, e_there
   logical                    :: checkiso_tmp, cbonds_tmp
   real(wp), allocatable      :: e_fix(:), e_clus(:)
   real(wp), parameter         :: eh = 627.509541d0
   real(wp)                   :: S, H, G, dens, shr, shr_av
   real(wp)                   :: sasa
   real(wp)                   :: newtemp, newmdtime, newmdstep, newhmass
   real(wp)                   :: newmetadlist, newmetadexp, newmetadfac
   real(wp)                   :: optlev_tmp
   real(wp)                   :: e0
   real(wp), allocatable      :: de(:)
   real(wp), allocatable      :: p(:)
   integer                    :: ich98, ich65, ich48
   logical                    :: not_param = .false.
   type(timer)                :: tim_dum !Dummy timer to avoid double counting

   interface
      subroutine aver(pr, env, runs, e_tot, S, H, G, sasa, a_present, a_tot)
         use iso_fortran_env, only: wp => real64
         use crest_data

         implicit none
         type(systemdata), intent(in)     :: env
         integer, intent(in)             :: runs
         real(wp), intent(inout)         :: e_tot
         real(wp), intent(in), optional  :: a_tot
         real(wp), intent(out)           :: S
         real(wp), intent(out)           :: H
         real(wp), intent(out)           :: G
         real(wp), intent(out)           :: sasa
         logical, intent(in)             :: pr, a_present
         dimension e_tot(runs)
         dimension a_tot(runs)
      end subroutine aver
   end interface

   if (.not. env%solv_md) then
      call tim%start(6, 'Solute-Ensemble')
   else
      call tim%start(7, 'Solvent-Ensemble')
   end if

   call tim_dum%init(20)

!--- Setting up directories
   call getcwd(thispath)
   f = makedir(fname_results)
   call chdir(fname_results)
   call getcwd(resultspath)
   call chdir(thispath)

!--- Setting defaults
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
   write (env%cts%pots(1), '("$wall")')
   write (env%cts%pots(2), '(2x,"potential=polynomial")')
   write (env%cts%pots(3), '(2x,"ellipsoid:",1x,3(g0,",",1x),"all")') clus%ell_abc
   if (.not. env%solv_md) write (env%cts%pots(4), '(2x,"ellipsoid:",1x,3(g0,",",1x),"1-",i0)')&
          & solu%ell_abc, solu%nat

   if (env%cts%used) then
      call write_reference(env, solu, clus) !new fixed file
      call copysub(env%fixfile, env%scratchdir)
   end if

   call chdir(env%scratchdir)
   scratchdir_tmp = env%scratchdir
   if (.not. env%solv_md) then
      io = makedir('tmp_MTD')
      call copysub('.CHRG', 'tmp_MTD')
      call copysub('.UHF', 'tmp_MTD')
      if (env%cts%used) call copysub(env%fixfile, 'tmp_MTD')
      call chdir('tmp_MTD')
   else
      io = makedir('tmp_solv_MTD')
      call chdir('tmp_solv_MTD')
   end if
   call getcwd(tmppath2)
   call wrc0('crest_input', clus%nat, clus%at, clus%xyz)

   if (env%solv_md) then
      call wr_cluster_cut('crest_input', solu%nat, solv%nat, env%nsolv,&
             & 'solute_cut.coord', 'solvent_shell.coord')
      call remove('crest_input')
      call copy('solvent_shell.coord', 'crest_input')
      deallocate (clus%at)
      deallocate (clus%xyz)
      call rdnat('solvent_shell.coord', clus%nat)
      allocate (clus%at(clus%nat), clus%xyz(3, clus%nat))
      call rdcoord('solvent_shell.coord', clus%nat, clus%at, clus%xyz)
   end if

   !For newcregen: If env%crestver .eq. crest_solv .and. .not. env%QCG then conffile .eq. .true.
   env%QCG = .false. 
   call inputcoords(env, 'crest_input')
   call defaultGF(env)         !Setting MTD parameter

!--- Special constraints for gff to safeguard stability
   if (env%ensemble_opt .eq. '--gff') then
      checkiso_tmp = env%checkiso
      env%checkiso = .true.
      cbonds_tmp = env%cts%cbonds_md
      env%cts%cbonds_md = .true.
      call autoBondConstraint_withEZ('coord', env%forceconst, env%wbofile)
      call rd_cbonds('bondlengths', env)
   end if

   gfnver_tmp = env%gfnver
   write (*, *) '  Method for ensemble search:', env%ensemble_opt
!  if (env%ens_const) write(*,*) '  Solute fixed during ensemble generation'
   env%gfnver = env%ensemble_opt  !Setting method for ensemble search

   !----------------------------------------------------------------
   ! Case selection of normal Crest, MD or MTD
   !----------------------------------------------------------------

   select case (env%ensemble_method)
   case (-1:0) !qcgmtd/Crest runtype

      !Defaults
      !General settings:
      if (.not. env%user_mdstep) then
         if (env%ensemble_opt .EQ. '--gff') then
            env%mdstep = 1.5d0
         else
            env%mdstep = 5.0d0
         end if
      end if
      !Runtype specific settings:
      if(env%ensemble_method == 0) then
         if (.not. env%user_dumxyz) then
            env%mddumpxyz = 200
         end if
         if (.not. env%user_mdtime) then
            env%mdtime = 10.0
         end if
      else if(env%ensemble_method == -1) then
         if (.not. env%user_dumxyz) then
            env%mddumpxyz = 50
         end if
         if (.not. env%user_mdtime) then
            env%mdtime = 5.0
         end if
         env%nmdtemp = 100
         env%MaxRestart = 6
      endif

      env%iterativeV2 = .true.  !Safeguards more precise ensemble search
      write (*, *) 'Starting ensemble cluster generation by CREST routine'
      call confscript2i(env, tim_dum) !Calling ensemble search
      call copy('crest_rotamers.xyz', 'crest_rotamers_0.xyz')

   case (1:2) ! Single MD or MTD

      !---- Setting threads
      call new_ompautoset(env,'auto',1,T,Tn)

      !--- Setting new defaults for MD/MTD in qcg
      if (env%mdtemp .lt. 0.0d0) then
         newtemp = 400.00d0
      else if (.not. env%user_temp) then
         newtemp = 298.0
      else
         newtemp = env%mdtemp
      end if

      if (.not. env%user_mdtime) then
         newmdtime = 100.0 !100.0
      else
         newmdtime = env%mdtime
      end if

      if (.not. env%user_dumxyz) then
         env%mddumpxyz = 1000
      end if

      if (.not. env%user_mdstep) then
         if (env%ensemble_opt .ne. '--gff') then
            newmdstep = 4.0d0
         else
            newmdstep = 1.5d0
         end if
      else
         newmdstep = env%mdstep
      end if

      if (env%ensemble_opt .ne. '--gff') then
         newhmass = 4.0
      else
         newhmass = 5.0
      end if

      if (.not. allocated(env%metadfac)) then
         allocate (env%metadfac(1))
         allocate (env%metadexp(1))
         allocate (env%metadlist(1))
      end if
      newmetadfac = 0.02_wp
      newmetadexp = 0.1_wp
      newmetadlist = 10.0_wp

      fname = 'coord'
      pipe = ' > xtb.out 2>/dev/null'

      !--- Writing constraining file xcontrol
      !--- Providing xcontrol overwrites constraints in coord file

      open (newunit=ich, file='xcontrol')
      if (env%cts%NCI) then
         do i = 1, 10
            if (trim(env%cts%pots(i)) .ne. '') then
               write (ich, '(a)') trim(env%cts%pots(i))
            end if
         end do
      end if

      if (.not. env%solv_md) then
         write (ich, '(a)') '$constrain'
         write (ich, '(2x,a,i0)') 'atoms: 1-', solu%nat
         write (ich, '(2x,a)') 'force constant=0.5'
         write (ich, '(2x,a,a)') 'reference=ref.coord'
      end if

      write (ich, '(a)') '$md'
      write (ich, '(2x,a,f10.2)') 'hmass=', newhmass
      write (ich, '(2x,a,f10.2)') 'time=', newmdtime
      write (ich, '(2x,a,f10.2)') 'temp=', newtemp
      write (ich, '(2x,a,f10.2)') 'step=', newmdstep
      write (ich, '(2x,a,i0)') 'shake=', env%shake
      write (ich, '(2x,a,i0)') 'dump=', env%mddumpxyz
      write (ich, '(2x,a)') 'dumpxyz=500.0'

      if (env%ensemble_method .EQ. 2) then
         write (ich, '(a)') '$metadyn'
         write (ich, '(2x,a,i0,a,i0)') 'atoms: ', solu%nat + 1, '-', clus%nat
         write (ich, '(2x,a,f10.2)') 'save=', newmetadlist
         write (ich, '(2x,a,f10.2)') 'kpush=', newmetadfac
         write (ich, '(2x,a,f10.2)') 'alp=', newmetadexp
      end if

      if (env%cts%cbonds_md) call write_cts_CBONDS(ich, env%cts)

      close (ich)

!--- Writing jobcall
      write (jobcall, '(a,1x,a,1x,a,'' --md --input xcontrol '',a,1x,a,a)') &
            &     trim(env%ProgName), trim(fname), trim(env%gfnver), trim(env%solv), pipe
!--- slightly different jobcall for QMDFF usage
      if (env%useqmdff) then
         write (jobcall, '(a,1x,a,1x,a,'' --md --input xcontrol --qmdff'',a,1x,a,a)') &
            &     trim(env%ProgName), trim(fname), trim(env%gfnver), trim(env%solv), pipe
      end if

!--- MD
      if (env%ensemble_method .EQ. 1) then
         call normalMD(fname, env, 1, newtemp, newmdtime)
         write (*, *) 'Starting MD with the settings:'
         write (*, '(''     MD time /ps        :'',f8.1)') newmdtime
         write (*, '(''     MD Temperature /K  :'',f8.1)') newtemp
         write (*, '(''     dt /fs             :'',f8.1)') newmdstep
         write (tmppath, '(a,i0)') 'NORMMD1'

         r = makedir(tmppath)
         call copysub('xcontrol', tmppath)
         call chdir(tmppath)
         call copy('coord', 'ref.coord')
         call chdir(tmppath2)

         call command('cd '//trim(tmppath)//' && '//trim(jobcall), io)

         inquire (file=trim(tmppath)//'/'//'xtb.trj', exist=ex)
         if (.not. ex .or. io .ne. 0) then
            write (*, '(a,i0,a)') '*Warning: MD seemingly failed (no xtb.trj)*'
         else
            write (*, *) '*MD finished*'
         end if

         if (env%trackorigin) then
            call set_trj_origins('NORMMD', 'md')
         end if
         call chdir('NORMMD1')
      end if

!--- MTD

      if (env%ensemble_method .EQ. 2) then
         call MetaMD(env, 1, newmdtime, env%metadfac(1), env%metadexp(1), &
            &               env%metadlist(1))
         write (*, '(a,i4,a)') 'Starting Meta-MD with the settings:'
         write (*, '(''     MTD time /ps       :'',f8.1)') newmdtime
         write (*, '(''     dt /fs             :'',f8.1)') newmdstep
         write (*, '(''     MTD Temperature /K  :'',f8.1)') newtemp
         write (*, '(''     dumpstep(trj) /fs  :'',i8)') env%mddumpxyz
         write (*, '(''     Vbias factor k /Eh :'',f8.4)') newmetadfac
         write (*, '(''     Vbias exp α /bohr⁻²:'',f8.2)') newmetadexp

         write (tmppath, '(a,i0)') 'METADYN1'
         r = makedir(tmppath)
         call copysub('xcontrol', tmppath)
         call chdir(tmppath)
         call copy('coord', 'ref.coord')

         call chdir(tmppath2)

         call command('cd '//trim(tmppath)//' && '//trim(jobcall), io)

         inquire (file=trim(tmppath)//'/'//'xtb.trj', exist=ex)
         if (.not. ex .or. io .ne. 0) then
            write (*, '(a,i0,a)') '*Warning: Meta-MTD seemingly failed (no xtb.trj)*'
         else
            write (*, *) '*MTD finished*'
         end if

         if (env%trackorigin) then
            call set_trj_origins('METADYN', 'mtd')
         end if

         call chdir('METADYN1')

      end if

      call rename('xtb.trj', 'crest_rotamers_0.xyz')
      call copysub('crest_rotamers_0.xyz', tmppath2)
      call dum%open('crest_rotamers_0.xyz')

!--- M(T)D stability check
      call minigrep('xtb.out', 'M(T)D is unstable, emergency exit', mdfail)
      if (dum%nall .eq. 1) then
         call copysub('xtb.out', resultspath)
         write (*, *) 'ERROR : M(T)D results only in one structure'
         if (mdfail) then
            write (*, *) '        It was unstable'
         else
            write (*, *) '        The M(T)D time step might be too large or the M(T)D time too short.'
         end if
         call copysub('xtb.out', resultspath)
         error stop '         Please check the xtb.out file in the ensemble folder'
      end if
      if (mdfail) then
         write (*, *)
         write (*, *) '   WARNING: The M(T)D was unstable.'
         write (*, *) '            Please check the xtb.out file in the ensemble folder.'
         write (*, *)
         call copysub('xtb.out', resultspath)
      end if
      call dum%deallocate
      call chdir(tmppath2)
      call wrc0('coord', clus%nat, clus%at, clus%xyz)
      call inputcoords(env, 'coord') !Necessary

!--- Optimization
      call print_qcg_opt
      !if (env%gfnver .eq. '--gfn2') 
      call multilevel_opt(env, 99)

   end select

   env%QCG = .true.

!--- Optimization with gfn2 if necessary
      if (env%final_gfn2_opt) then
      gfnver_tmp = env%gfnver
!      if (env%gfnver .ne. '--gfn2') then
      write (*, '(2x,a)') 'GFN2-xTB optimization'
      env%gfnver = '--gfn2'
      call rmrf('OPTIM')
      call multilevel_opt(env, 99)
   end if

!--- Final optimization without potentials
   call rmrf('OPTIM')
   env%optlev = 1.0d0    !Higher precision for less scattering
   env%cts%NCI = .false.  !Dactivating the wall pot.
   env%cts%pots = ''
   deallocate (env%cts%pots)
   call multilevel_opt(env, 99)

  !Clustering to exclude similar structures if requested with -cluster
  if (env%properties == 70) then
    write(*,'(3x,''Clustering the remaining structures'')')
    call checkname_xyz(crefile,inpnam,outnam)
    call ccegen(env, .false. , inpnam)
    call move(trim(clusterfile),trim(outnam))
  end if

!--- Energy sorting and removal of dublicates
   env%gbsa = gbsa_tmp
   env%solv = solv_tmp
   call newcregen(env, 0)
   call checkname_xyz(crefile, inpnam, outnam)
   call copy(inpnam, 'ensemble.xyz')
   call ens%open('ensemble.xyz') !Read in ensemble
   call clus%deallocate()
   clus%nat = ens%nat
   allocate (clus%at(clus%nat))
   allocate (clus%xyz(3, clus%nat))

!-------------------------------------------------------------
!      SP with GBSA model and without wall potentials
!-------------------------------------------------------------

   !--- Write folder with xyz-coordinates
   do i = 1, ens%nall
      call rdxmolselec('ensemble.xyz', i, clus%nat, clus%at, clus%xyz)
      write (to, '("TMPSP",i0)') i
      io = makedir(trim(to))
      call copysub('.UHF', to)
      call copysub('.CHRG', to)
      call chdir(to)
      call wrxyz('cluster.xyz', clus%nat, clus%at, clus%xyz*bohr)
      call chdir(tmppath2)
   end do
   !--- SP
   write (*, *)
   call ens_sp(env, 'cluster.xyz', ens%nall, 'TMPSP')
   !--- Getting energy
   do i = 1, ens%nall
      call rdxmolselec('ensemble.xyz', i, clus%nat, clus%at, clus%xyz)
      write (to, '("TMPSP",i0)') i
      call chdir(to)
      call grepval('xtb_sp.out', '| TOTAL ENERGY', e_there, ens%er(i))
      call chdir(tmppath2)
   end do

   if (.not. e_there) then
      write (*, *)
      write (*, *) 'Energy not found. Error in xTB computations occured'
      call chdir(to)
      call minigrep('xtb_sp.out', 'solv_model_loadInternalParam', not_param)
      call chdir(tmppath2)
      if (not_param) then
         write (*, *) '  !!!WARNIG: CHOSEN SOLVENT NOT PARAMETERIZED &
         & FOR IMPLICIT SOLVATION MODEL!!!'
         write (*, '(''  CHECK IF '',A,'' IS AVAILABLE IN xTB'')') env%solv
         write (*, *) '  PLEASE RESTART THE ENSEMBLE GENERATION WITH AVAILABLE&
                &  PARAMETERIZATION IF YOU NEED ENERGIES'
         call copysub('crest_conformers.xyz', resultspath)
         write (*, *) '  The enesemble can be found in the <ensemble> directory&
                & as <crest_conformers.xyz>'
         error stop
      end if
   end if

   env%gfnver = gfnver_tmp
   call ens%write('full_ensemble.xyz')

!--- crest_best structure
   minpos = minloc(ens%er, dim=1)
   write (to, '("TMPSP",i0)') minpos
   call chdir(to)
   call rdxmol('cluster.xyz', clus%nat, clus%at, clus%xyz)
   call chdir(tmppath2)
   write (comment, '(F20.8)') ens%er(minpos)
   inquire (file='crest_best.xyz', exist=ex)
   if (ex) then
      call rmrf('crest_best.xyz') !remove crest_best from
   end if
   call wrxyz('crest_best.xyz', clus%nat, clus%at, clus%xyz, comment)

!-------------------------------------------------------------
!      Processing results
!-------------------------------------------------------------

   allocate (e_fix(ens%nall))
   allocate (e_clus(ens%nall))

   call pr_ensemble_energy()

   open (newunit=ich98, file='cluster_energy.dat')
   write (ich98, '(3x,''#'',9x,''Energy [Eh]'',6x,''SASA'')')

!--- Fixation energy of optimization
   do i = 1, ens%nall
      call chdir('OPTIM')
      write (to, '("TMPCONF",i0)') i
      call chdir(to)
      call grepval('xtb.out', '         :: add. restraining', e_there, e_fix(i))
      call chdir(tmppath2)

      call rdxmolselec('full_ensemble.xyz', i, clus%nat, clus%at, clus%xyz)
      call get_sphere(.false., clus, .false.)
      dens = 0.001*(solu%mass + env%nsolv*solv%mass)/(1.0d-30*clus%vtot*bohr**3)
      if (env%solv_md) then
         call analyze_cluster(env%nsolv - 1, clus%nat, solv%nat, solv%nat, clus%xyz, clus%at, shr_av, shr)
      else
         call analyze_cluster(env%nsolv, clus%nat, solu%nat, solv%nat, clus%xyz, clus%at, shr_av, shr)
      end if
      write (ich98, '(i4,F20.10,3x,f8.1)') env%nsolv, ens%er(i), clus%atot
      write (*, '(x,i4,4x,F13.6,2x,f6.3,1x,f8.3,2x,2f6.1,3x,f8.1,3x,a)') &
            & i, ens%er(i), dens, e_fix(i), shr_av, shr, clus%atot, trim(optlevflag(env%optlev))
      e_fix(i) = e_fix(i)*eh/sqrt(float(clus%nat))
   end do
   close (ich98)
   call copysub('cluster_energy.dat', resultspath)

!--- Checking Boltzmann weighting
   write (*, *)
   call remove('full_ensemble.xyz')
   call sort_ensemble(ens, ens%er, 'full_ensemble.xyz')
   e_clus = ens%er*eh
   call sort_min(ens%nall, 1, 1, e_clus)
   ens%er = e_clus/eh !Overwrite ensemble energy with sorted one
   allocate (de(ens%nall), source=0.0d0)
   allocate (p(ens%nall), source=0.0d0)
   e0 = e_clus(1)
   de(1:ens%nall) = (e_clus(1:ens%nall) - e0)
   call qcg_boltz(env, ens%nall, de, p)
   k = 0
   if (.not. env%user_nclust) env%nqcgclust = 0 !Needed for solvent ensemble
   if (env%nqcgclust .eq. 0) then
      do i = 1, ens%nall !Count how many are above 10%
         if ((p(i)) .gt. 0.1) then
            k = k + 1
         end if
      end do
      if ((k .eq. 0) .or. (k .gt. 10)) then
         k = 10 !If too many structures are relevant, set it 10
      else if ((k .lt. 4) .and. (ens%nall .ge. 4)) then
         k = 4 !If too less structures are relevant, set it 4
      else if (ens%nall .gt. 0) then 
         k=ens%nall
      else
         error stop 'No structure left. Something went wrong.' 
      end if
      write (*, '(2x,a,1x,i0)') 'Conformers taken:', k
      env%nqcgclust = k
   else
      if (env%nqcgclust .gt. ens%nall) then
         k = ens%nall !Input larger than remaining structures
         write (*, '(''Less than '',1x,i0,1x,''structures remain'')') env%nqcgclust
         write (*, '(''Only '',1x,i0,1x,''structures are taken'')') ens%nall
         if (env%cff) env%nqcgclust = ens%nall !Only for CFF, else a second qcg_ensemble run starts for solvent
      else
         write (*, '(''Taking '',1x,i0,1x,''structures'')') env%nqcgclust
         k = env%nqcgclust !user input
      end if
   end if

   open (newunit=ich65, file='final_ensemble.xyz')
   do i = 1, k
      open (newunit=ich48, file='full_population.dat')
      write (ich48, '(2x, ''cluster'',2x,''E_norm [Eh]'',2x, ''De [kcal]'', 4x, ''p'')')
      do j = 1, ens%nall
         if (j .lt. 10) then
            write (ich48, '(5x,i0,3x,f11.6,5x,f6.4,3x,f6.4)') j, e_clus(j)/eh, de(j), p(j)
         else
            write (ich48, '(5x,i0,2x,f11.6,5x,f6.4,3x,f6.4)') j, e_clus(j)/eh, de(j), p(j)
         end if
      end do
      close (ich48)

!--- Take k energetic least structures (written at beginning of file)
      call rdxmolselec('full_ensemble.xyz', i, clus%nat, clus%at, clus%xyz)
      call wrxyz(ich65, clus%nat, clus%at, clus%xyz*bohr, ens%er(i))
   end do
   close (ich65)

   call ens%deallocate()
   call ens%open('final_ensemble.xyz')
   ens%er = e_clus(1:k)/eh

!--- Getting G,S,H
   write (*, *)
   write (*, '(2x,''------------------------------------------------------------------------'')')
   write (*, '(2x,''------------------------------------------------------------------------'')')
   write (*, '(2x,''Boltz. averaged energy of final cluster:'')')
   call aver(.true., env, ens%nall, e_clus(1:ens%nall), S, H, G, sasa, .false.)
   write (*, '(7x,''G /Eh     :'',F14.8)') G/eh
   write (*, '(7x,''T*S /kcal :'',f8.3)') S

   ens%g = G
   ens%s = S

   deallocate (e_fix)
   deallocate (e_clus)

!---Folder management
   call rename('cregen.out.tmp', 'thermo_data')
   call copysub('thermo_data', resultspath)
   call copysub('crest_best.xyz', resultspath)
   call copysub('cre_members.out', resultspath)
   call copysub('full_ensemble.xyz', resultspath)
   call copysub('final_ensemble.xyz', resultspath)
   call copysub('population.dat', resultspath)
   call copysub('full_population.dat', resultspath)

!---Deleting ensemble tmp
   call chdir(thispath)
   call chdir(env%scratchdir)
   if (.not. env%keepModef) call rmrf(tmppath2)
!----Outprint
   write (*, *)
   write (*, '(2x,''Ensemble generation finished.'')')
   write (*, '(2x,''Results can be found in ensemble directory'')')
   write (*, '(2x,''Lowest energy conformer in file <crest_best.xyz>'')')
   write (*, '(2x,''List of full ensemble in file <full_ensemble.xyz>'')')
   write (*, '(2x,''List of used ensemble in file <final_ensemble.xyz>'')')
   write (*, '(2x,''Thermodynamical data in file <thermo_data>'')')
   write (*, '(2x,''Population of full ensemble in file <full_population.dat>'')')
   write (*, '(2x,''Population in file <population.dat>'')')

   env%gfnver = gfnver_tmp
   env%optlev = optlev_tmp
   if (env%ensemble_opt .eq. '--gff') then
      env%cts%cbonds_md = cbonds_tmp
      env%checkiso = checkiso_tmp
   end if

   call tim_dum%clear

   if (.not. env%solv_md) then
      call tim%stop(6)
   else
      call tim%stop(7)
   end if

end subroutine qcg_ensemble

subroutine qcg_cff(env, solu, solv, clus, ens, solv_ens, tim)
   use crest_parameters
   use crest_data
   use iomod
   use zdata
   use strucrd

   implicit none

   type(systemdata)           :: env
   type(timer)                :: tim
   type(zmolecule)            :: solu, solv, clus
   type(ensemble)             :: solv_ens
   type(ensemble), intent(in)  :: ens

   integer                    :: i, j, k, iter
   integer                    :: io, r
   integer                    :: nsolv, n_ini
   integer                    :: ipos, dum
   integer                    :: v_ratio
   integer                    :: minE_pos, m, nat_tot
   integer                    :: nat_frag1 !number of atoms larger fragment (=solvent shell)
   integer                    :: conv(env%nqcgclust + 1)
   integer                    :: solv_added, minpos
   character(len=512)         :: thispath, resultspath, tmppath, tmppath2
   character(len=64)          :: fname_lmo1, fname_lmo2, comment
   character(len=20)          :: to
   real(wp), allocatable      :: e_empty(:), inner_ell_abc(:, :)
   real(wp), allocatable      :: outer_ell_abc(:, :)
   real(wp), allocatable      :: e_cur(:, :)
   real(wp)                   :: e_cluster(env%nqcgclust)
   real(wp), parameter        :: eh = 627.509541d0
   real(wp)                   :: S, H, G
   real(wp)                   :: sasa, tmp_optlev
   real(wp)                   :: etmp(500)
   real(wp)                   :: e_fix(env%nqcgclust), e_norm(env%nqcgclust)
   real(wp)                   :: dum_e, de
   real(wp)                   :: de_tot(env%nqcgclust)
   real(wp)                   :: shr = 0
   real(wp)                   :: shr_av = 0
   real(wp)                   :: dens, atotS
   logical                    :: ex, skip, e_there
   logical                    :: all_converged
   logical, allocatable        :: converged(:), nothing_added(:)

   character(len=20)          :: gfnver_tmp
   real(wp)                   :: optlev_tmp
   integer                    :: ich98, ich31

   interface
      subroutine aver(pr, env, runs, e_tot, S, H, G, sasa, a_present, a_tot)
         use iso_fortran_env, only: wp => real64
         use crest_data

         implicit none
         type(systemdata), intent(in)     :: env
         integer, intent(in)             :: runs
         real(wp), intent(inout)         :: e_tot
         real(wp), intent(in), optional  :: a_tot
         real(wp), intent(out)           :: S
         real(wp), intent(out)           :: H
         real(wp), intent(out)           :: G
         real(wp), intent(out)           :: sasa
         logical, intent(in)             :: pr, a_present
         dimension e_tot(runs)
         dimension a_tot(runs)
      end subroutine aver
   end interface

   call tim%start(8, 'CFF')

   allocate (e_empty(env%nqcgclust))
   allocate (converged(env%nqcgclust))
   allocate (nothing_added(env%nqcgclust))
   allocate (outer_ell_abc(env%nqcgclust, 3))
   allocate (inner_ell_abc(env%nqcgclust, 3))

   v_ratio = nint(solu%vtot/solv%vtot)
   allocate (e_cur(env%nsolv + v_ratio, env%nqcgclust), source=0.0d0)

!--- Setting defaults (same as ensemble optimization to have comparable structures)
   optlev_tmp = env%optlev
   env%optlev = 1.0d0    !Increaseing percision for ensemble search to minimze scattering
   gfnver_tmp = env%gfnver
   if (env%final_gfn2_opt) then
      env%gfnver = '--gfn2'
   else
      env%gfnver = env%ensemble_opt !CFF always with ensemble method
   end if
   nothing_added = .false.

   dum = 0
   converged = .false.
   all_converged = .false.
   nat_tot = clus%nat - solu%nat!*env%nqcgclust

   if (solu%vtot/solv%vtot .lt. 1.0d0) then
      skip = .true.
   else
      skip = .false.
   end if

!--- Folder management
   call getcwd(thispath)
   r = makedir('solvent_ensemble')
   call chdir('solvent_ensemble')
   call getcwd(resultspath)
   call chdir(thispath)
   call chdir(env%scratchdir)
   call getcwd(tmppath)
   io = makedir('tmp_CFF')
   call chdir('tmp_CFF')
   call getcwd(tmppath2)
   call chdir(tmppath)
   call chdir('solvent_properties')
   if (env%use_xtbiff) then
      call copysub('solvent.lmo', tmppath2)
   else
      call copysub('solvent', tmppath2)
   end if
   call chdir(tmppath2)

!--- SP of each cluster
   call ens%write('ensemble.xyz')
   do i = 1, env%nqcgclust
      call rdxmolselec('ensemble.xyz', i, clus%nat, clus%at, clus%xyz)
      clus%nmol = clus%nat/solv%nat
      write (to, '("TMPCFF",i0)') i
      io = makedir(trim(to))
      if (env%use_xtbiff) then
         call copysub('solvent.lmo', to)
      else
         call copysub('solvent', to)
      end if
      call chdir(to)
      call wrc0('cluster.coord', clus%nat, clus%at, clus%xyz)
      call wr_cluster_cut('cluster.coord', solu%nat, solv%nat, env%nsolv, 'solute_cut.coord', 'solvent_shell.coord')
      call xtb_sp_qcg(env, 'solvent_shell.coord')
      call grepval('xtb.out', '| TOTAL ENERGY', ex, e_empty(i))
      call copy('solvent_shell.coord', 'solvent_cluster.coord')
      call copy('solvent_cluster.coord', 'filled_cluster.coord')
      call get_ellipsoid(env, solu, solv, clus, .false.) !solu, to have same cavity to fill solvent in
      outer_ell_abc(i, 1:3) = clus%ell_abc(1:3)
      inner_ell_abc(i, 1:3) = solu%ell_abc(1:3)
      call chdir(tmppath2)
   end do

   if (skip) write (*, '(2x,''solute smaller than solvent, cff skipped'')')

   clus%nat = clus%nat - solu%nat
   n_ini = clus%nat

!--- If solvent molecules are added
   if (.not. skip) then
      call pr_qcg_fill()
      write (*, '(2x,''now adding solvents to fill cluster...'')')
      call pr_fill_energy()
      write (*, '(2x,''------------------------------------------------------------------------'')')
      nat_frag1 = env%nsolv*solv%nat

      iter = 0
!--- Main cycle for addition of solvent molecules
      convergence: do while (.not. all_converged)
         k = 0
         iter = iter + 1
         !--- Setting array, with only numbers of dirs that are not converged
         do i = 1, env%nqcgclust
            if (.not. converged(i)) then
               k = k + 1
               conv(k) = i
               conv(env%nqcgclust + 1) = k !How many jobs are open
            else
               cycle
            end if
         end do
         conv(k + 1:env%nqcgclust) = 0

         if (env%use_xtbiff) then
!----------- LMO computation for solvent cluster---------------------------------------------------
            call ensemble_lmo(env, 'solvent_cluster.coord', solv, conv(env%nqcgclust + 1),&
                   & 'TMPCFF', conv)
!--------------------------------------------------------------------------------------------------

            do i = 1, env%nqcgclust
               if (.not. converged(i)) then
                  write (to, '("TMPCFF",i0)') i
                  call chdir(to)
                  call rename('xtblmoinfo', 'solvent_cluster.lmo')
                  call chdir(tmppath2)
               else
                  cycle
               end if
            end do
         end if

         call chdir(tmppath2)

         fname_lmo1 = 'solvent_cluster.lmo'
         fname_lmo2 = 'solvent.lmo'

!--- Solvent addition to the cluster---------------------------------------------
         if (env%use_xtbiff) then
            call ensemble_iff(env, outer_ell_abc, nat_frag1, fname_lmo1, fname_lmo2,&
                    &conv(env%nqcgclust + 1), 'TMPCFF', conv)
         else
            call ensemble_dock(env, outer_ell_abc, nat_frag1, 'solvent_cluster.coord',&
                    &'solvent', clus%nat, solv%nat, conv(env%nqcgclust + 1), 'TMPCFF', conv)
         end if
!--------------------------------------------------------------------------------

         nat_frag1 = nat_frag1 + solv%nat

         !--- Increase cluster size
         deallocate (clus%at)
         deallocate (clus%xyz)
         clus%nat = clus%nat + solv%nat
         allocate (clus%at(clus%nat))
         allocate (clus%xyz(3, clus%nat))
         clus%nmol = clus%nmol + 1

         do i = 1, env%nqcgclust
            if (.not. converged(i)) then
               write (to, '("TMPCFF",i0)') i
               call chdir(to)
               call remove('xtbrestart')
               call remove('xcontrol')

               if (env%use_xtbiff) then
                  !--- Select xtb-IFF stucture to proceed
                  call rdxtbiffE('xtbscreen.xyz', m, clus%nat, etmp) !Get energy of screening
                  minE_pos = minloc(etmp(1:m), dim=1)            !Get minimum of those
                  call rdxmolselec('xtbscreen.xyz', minE_pos, clus%nat, clus%at, clus%xyz) !Read the struc into clus%xyz
                  call wrc0('solvent_cluster.coord', clus%nat, clus%at, clus%xyz)
               else
                  call rdcoord('best.xyz', clus%nat, clus%at, clus%xyz, e_cur(iter, i))
                  call wrc0('solvent_cluster.coord', clus%nat, clus%at, clus%xyz)
               end if

               !--- Check if converged
               call fill_take(env, solv%nat, clus%nat, inner_ell_abc(i, 1:3), ipos)
               if (ipos .eq. 0) then
                  converged(i) = .true.
                  write (*, '(2x,''no more solvents can be placed inside cavity of cluster: '',i0)') i
                  write (*, '(2x,''previous cluster taken...'')')
                  if (iter .eq. 1) nothing_added(i) = .true.
               end if
               call chdir(tmppath2)

            else
               cycle
            end if
         end do

!--- Check, if a structure was converged and iff was not necessary
         k = 0
         do i = 1, env%nqcgclust
            if (.not. converged(i)) then
               k = k + 1
               conv(k) = i
               conv(env%nqcgclust + 1) = k !How many jobs are open
            else
               cycle
            end if
         end do
         conv(k + 1:env%nqcgclust) = 0

!      if(env%use_xtbiff) then
!--- Parallel optimization-------------------------------------------------------------------
         call cff_opt(.false., env, 'solvent_cluster.coord', n_ini, conv(env%nqcgclust + 1)&
                 &, 'TMPCFF', conv, nothing_added)
!----------------------------------------------------------------------------------------------
!      end if

         do i = 1, env%nqcgclust
            if (.not. converged(i)) then
               write (to, '("TMPCFF",i0)') i
               call chdir(to)
               if (env%use_xtbiff) then
                  call copy('xtbopt.coord', 'solvent_cluster.coord')
                  call grepval('xtb_sp.out', '| TOTAL ENERGY', e_there, e_cur(iter, i))
               end if
               dum_e = e_empty(i)
               if (iter - nsolv .gt. 1) dum_e = e_cur(iter - 1, i)
               de = eh*(e_cur(iter, i) - solv%energy - dum_e)
               de_tot(i) = de_tot(i) + de
               !---- Check if solvent added is repulsive
               if (de .gt. 0) then
                  converged(i) = .true.
                  write (*, '(2x,''adding solvent is repulsive for cluster: '',i0)') i
                  write (*, '(2x,''previous cluster taken...'')')
                  if (iter .eq. 1) nothing_added(i) = .true.
               else !Only if the addition was not repulsive
                  call copy('solvent_cluster.coord', 'filled_cluster.coord')
                  write (*, '(i4,5x,i3,1x,F13.6,3x,f7.2,5x,f7.2,4x,a)') &
                     & iter + env%nsolv, i, e_cur(iter, i), de, de_tot(i),&
                     & trim(optlevflag(env%optlev))
               end if
               call chdir(tmppath2)
            end if
         end do

         !--- Check if everything is converged
         dum = 0
         do i = 1, env%nqcgclust
            if (converged(1)) then
               dum = dum + 1
            end if
         end do

         if (dum .eq. env%nqcgclust) then
            all_converged = .true.
         else
            nat_tot = nat_tot + solv%nat
         end if

         write (*, '(2x,''------------------------------------------------------------------------'')')
         !--- Or if maximum solvent is added
         if (iter - nsolv .eq. v_ratio) then
            write (*, '(2x,''volume filled'')')
            all_converged = .true.
            call copy('solvent_cluster.coord', 'filled_cluster.coord')
         end if

      end do convergence

   end if

   !Now in every TMPPath the final cluster file filled_cluster.coord is present

!---------------------------------------------------------------------
!     Final Optimization
!---------------------------------------------------------------------

   tmp_optlev = env%optlev
   if (env%optlev .lt. 1.0) env%optlev = 1.0d0 !higher accuracy

   if (.not. skip) then
      call cff_opt(.true., env, 'filled_cluster.coord', n_ini, conv(env%nqcgclust + 1),&
             & 'TMPCFF', conv, nothing_added)
   else
      n_ini = 0 !If this is 0, no constraining will be done (optimization of total system)
      nothing_added = .true.
      call cff_opt(.true., env, 'filled_cluster.coord', n_ini, env%nqcgclust, 'TMPCFF',&
             & conv, nothing_added)
   end if
   env%optlev = tmp_optlev

   call pr_ensemble_energy()

   solv_ens%nall = env%nqcgclust
   solv_ens%nat = nat_tot

!--- Getting results--------------------------------------------------------------
   open (newunit=ich31, file='crest_rotamers_0.xyz')
   open (newunit=ich98, file='cluster_energy.dat')
   write (ich98, '(3x,''#'',11x,''Energy [Eh]'',6x,''SASA'')')

   do i = 1, env%nqcgclust
      write (to, '("TMPCFF",i0)') i
      call chdir(to)
      call copy('xtbopt.coord', 'final_cluster.coord')

!--- Reading structure
      call clus%deallocate()
      call rdnat('final_cluster.coord', clus%nat)
      allocate (clus%at(clus%nat), clus%xyz(3, clus%nat))
      call rdcoord('final_cluster.coord', clus%nat, clus%at, clus%xyz)

!--- Getting energy and calculating properties
      call grepval('xtb_sp.out', '| TOTAL ENERGY', e_there, e_cluster(i))
      call grepval('xtb_sp.out', '         :: add. restraining', e_there, e_fix(i))
      e_fix(i) = e_fix(i)*eh/sqrt(float(clus%nat))
      call get_sphere(.false., clus, .false.)
      if (clus%nat .gt. n_ini) then
         solv_added = (clus%nat - (n_ini))/solv%nat
      else
         solv_added = 0
      end if
      dens = 0.001*((clus%nat/solv%nat)*solv%mass)/(1.0d-30*clus%vtot*bohr**3)
      call analyze_cluster(solv_added, clus%nat, n_ini, solv%nat, clus%xyz, clus%at, shr_av, shr)
      e_norm(i) = e_cluster(i)*env%nsolv/(clus%nat/solv%nat)
      atotS = clus%atot*env%nsolv/(clus%nat/solv%nat)

!--- Writing outputfiles
      write (ich31, '(2x,i0)') clus%nat
      write (ich31, '(2x,f18.8,2x,a)') e_cluster(i)
      do j = 1, clus%nat
         write (ich31, '(1x,a2,1x,3f20.10)') i2e(clus%at(j), 'nc'), clus%xyz(1:3, j)*bohr
      end do

      write (ich98, '(''No'',i4,F20.10,3x,f8.1)') i, e_norm(i), atotS

!--- Print to screen
      write (*, '(x,i4,4x,F13.6,2x,f6.3,1x,f8.3,2x,2f6.1,3x,f8.1,3x,a)') &
              & i, e_norm(i), dens, e_fix(i), shr_av, shr, atotS, trim(optlevflag(env%optlev))

      call chdir(tmppath2)
   end do

   close (ich98)
   close (ich31)

   call solv_ens%deallocate()
   call solv_ens%open('crest_rotamers_0.xyz')

   solv_ens%er = e_cluster
   call copy('crest_rotamers_0.xyz', 'crest_ensemble.xyz')

!--- crest_best structure
   minpos = minloc(solv_ens%er, dim=1)
   write (to, '("TMPCFF",i0)') minpos
   call chdir(to)
   call clus%deallocate
   call rdnat('final_cluster.coord', clus%nat)
   allocate (clus%at(clus%nat), clus%xyz(3, clus%nat))
   call rdcoord('final_cluster.coord', clus%nat, clus%at, clus%xyz)
   clus%xyz = clus%xyz*bohr
   call chdir(tmppath2)
   write (comment, '(F20.8)') solv_ens%er(minpos)
   call wrxyz('crest_best.xyz', clus%nat, clus%at, clus%xyz, comment)

!--- Boltz. average-------------------------------------------------------------------------
   write (*, *)
   write (*, '(2x,''------------------------------------------------------------------------'')')
   write (*, '(2x,''------------------------------------------------------------------------'')')
   write (*, '(2x,''Boltz. averaged energy of final cluster:'')')
   e_cluster = solv_ens%er*eh
   e_norm = e_norm*eh
   call sort_min(env%nqcgclust, 1, 1, e_norm)
   call aver(.true., env, solv_ens%nall, e_norm(1:env%nqcgclust), S, H, G, sasa, .false.)
   write (*, '(7x,''G /Eh     :'',F14.8)') G/eh
   write (*, '(7x,''T*S /kcal :'',f8.3)') S
   solv_ens%er = e_norm/eh !normalized energy needed for final evaluation

   solv_ens%g = G
   solv_ens%s = S

!--- Cleanup
   call copysub('crest_ensemble.xyz', resultspath)
   call copysub('cluster_energy.dat', resultspath)
   call copysub('crest_best.xyz', resultspath)
   call copysub('population.dat', resultspath)
   call chdir(tmppath)
   if (.not. env%keepModef) call rmrf('tmp_CFF')
   call chdir(thispath)

!--- Printouts
   write (*, *)
   write (*, '(2x,''Solvent cluster generation finished.'')')
   write (*, '(2x,''Results can be found in solvent_cluster directory'')')
   write (*, '(2x,''Structures in file <crest_ensemble.xyz>'')')
   write (*, '(2x,''Energies in file <cluster_energy.dat>'')')
   write (*, '(2x,''Population in file <population.dat>'')')

   env%gfnver = gfnver_tmp
   env%optlev = optlev_tmp

   deallocate (e_empty)
   deallocate (converged)
   deallocate (outer_ell_abc)
   deallocate (inner_ell_abc)

   call tim%stop(8)

end subroutine qcg_cff

subroutine qcg_freq(env, tim, solu, solv, solu_ens, solv_ens)
   use crest_parameters
   use crest_data
   use iomod
   use zdata
   use strucrd

   implicit none

   type(systemdata)           :: env
   type(timer)                :: tim
   type(zmolecule)            :: solu, solv, clus
   type(ensemble)             :: solu_ens, solv_ens

   integer                    :: r, io, f, g, h
   integer                    :: i
   character(len=512)         :: thispath, resultspath, tmppath, tmppath2
   character(len=80)          :: to
   character(len=20)          :: gfnver_tmp
   real(wp)                   :: optlev_tmp
   real(wp)                   :: gt(3)
   real(wp)                   :: ht(3)
   real(wp)                   :: svib(3)
   real(wp)                   :: srot(3)
   real(wp)                   :: stra(3)
   integer                    :: ich65, ich56, ich33, ich81
   logical                    :: opt

   call tim%start(9, 'Frequencies')

   call pr_qcg_freq()

!--- Setting defaults (same as ensemble optimization and cff to have comparable structures)
   optlev_tmp = env%optlev
   env%optlev = 1.0d0    !Increaseing percision for ensemble search to minimze scattering
   gfnver_tmp = env%gfnver
   env%gfnver = env%freqver  !Setting method

!--- Folder management
   call getcwd(thispath)
   r = makedir('frequencies')
   call chdir('frequencies')
   call getcwd(resultspath)
   call chdir(thispath)
   call chdir(env%scratchdir)
   call getcwd(tmppath)
   io = makedir('tmp_freq')
   call copysub('.CHRG', 'tmp_freq')
   call copysub('.UHF', 'tmp_freq')
   call chdir('tmp_freq')
   call getcwd(tmppath2)
   f = makedir('tmp_solu')
   call copysub('.CHRG', 'tmp_solu')
   call copysub('.UHF', 'tmp_solu')
   g = makedir('tmp_solv')
   h = makedir('tmp_gas1') !One solute molecule
   call copysub('.CHRG', 'tmp_gas1')
   call copysub('.UHF', 'tmp_gas1')

!--- Frequencies solute molecule
   write (*, *) '  SOLUTE MOLECULE'
   call chdir('tmp_gas1')
   call wrc0('solute.coord', solu%nat, solu%at, solu%xyz)
   call chdir(tmppath2)
   opt = .false.
   call ens_freq(env, 'solute.coord', 1, 'tmp_gas', opt)
   call chdir('tmp_gas1')
   call rdtherm('xtb_freq.out', ht(3), svib(3), srot(3), stra(3), gt(3))
   solu%gt = gt(3)
   solu%ht = ht(3)
   solu%svib = svib(3)
   solu%srot = srot(3)
   solu%stra = stra(3)

   call chdir(tmppath2)

!--- Folder setup for cluster
   call chdir('tmp_solu')
   call solu_ens%write('solute_ensemble.xyz')

!--- All cluster are of the same size
   call clus%deallocate()
   clus%nat = solu_ens%nat
   allocate (clus%at(clus%nat))
   allocate (clus%xyz(3, clus%nat))
   clus%xyz = 0
   clus%nmol = env%nsolv + 1 !clus%nat/clus%at

   do i = 1, solu_ens%nall
      call rdxmolselec('solute_ensemble.xyz', i, clus%nat, clus%at, clus%xyz)

!--- Solute cluster
      write (to, '("TMPFREQ",i0)') i
      io = makedir(trim(to))
      call copysub('.UHF', to)
      call copysub('.CHRG', to)
      call chdir(to)
      open (newunit=ich65, file='cluster.xyz')
      call wrxyz(ich65, clus%nat, clus%at, clus%xyz*bohr)
      close (ich65)

      call chdir(tmppath2)

      !--- Solvent cluster (only if cff, than the solvent shell is taken, which was fixed all the time)
      if (env%cff) then
         call chdir('tmp_solv')
         write (to, '("TMPFREQ",i0)') i
         io = makedir(trim(to))
         call chdir(to)
         call wrc0('cluster.coord', clus%nat, clus%at, clus%xyz)
         call wr_cluster_cut('cluster.coord', solu%nat, solv%nat, env%nsolv,&
                & 'solute_cut.coord', 'solvent_cut.coord')

         call chdir(tmppath2)
      end if
      call chdir('tmp_solu')

   end do

   write (*, *) '  SOLUTE CLUSTER'

!> Frequency calculation
   opt = .true.
   call ens_freq(env, 'cluster.xyz', solu_ens%nall, 'TMPFREQ', opt)
   call chdir(tmppath2)

   write (*, *) '  SOLVENT CLUSTER'
   if (env%cff) then
      call chdir('tmp_solv')
      call ens_freq(env, 'solvent_cut.coord', solu_ens%nall, 'TMPFREQ', opt)
      call chdir(tmppath2)
   end if

   call clus%deallocate()

   !--- Frequencies solvent cluster (only, if not cff was used)
   if (.not. env%cff) then
      call chdir('tmp_solv')
      call solv_ens%write('solvent_ensemble.xyz')

      do i = 1, solv_ens%nall
         write (to, '("TMPFREQ",i0)') i
         io = makedir(trim(to))
         call copysub('.UHF', to)
         call copysub('.CHRG', to)
         call chdir(to)
         open (newunit=ich65, file='solv_cluster.xyz')
         call wrxyz(ich65, solv_ens%nat, solv_ens%at, solv_ens%xyz(:, :, i))
         close (ich65)
         call chdir(tmppath2)
         call chdir('tmp_solv')
      end do
!> Frequency calculation
      call ens_freq(env, 'solv_cluster.xyz', solv_ens%nall, 'TMPFREQ', opt)
      call chdir(tmppath2)
   end if

!----------------------------------------------------------------------------
!   Data read out
!----------------------------------------------------------------------------

!--- Solute in gas phase
   write (*, *)
   write (*, *) '  Solute Gas properties'
   call pr_freq_energy()
   open (newunit=ich56, file='solute.dat')
   call pr_freq_file(56)
   write (*, '(2x,5f10.2)') ht(3), svib(3), srot(3), stra(3), gt(3)
   write (ich56, '(2x,5f10.2)') ht(3), svib(3), srot(3), stra(3), gt(3)
   close (ich56)

!--- Solute cluster
   write (*, *)
   write (*, *) '  Solute cluster properties'
   open (newunit=ich33, file='solute_cluster.dat')

   call chdir('tmp_solu')

   allocate (solu_ens%gt(solu_ens%nall))
   allocate (solu_ens%ht(solu_ens%nall))
   allocate (solu_ens%svib(solu_ens%nall))
   allocate (solu_ens%srot(solu_ens%nall))
   allocate (solu_ens%stra(solu_ens%nall))

   call pr_freq_energy()
   call pr_freq_file(ich33)

   do i = 1, solu_ens%nall
      write (to, '("TMPFREQ",i0)') i
      call chdir(to)
      call rdtherm('xtb_freq.out', ht(1), svib(1), srot(1), stra(1), gt(1))
      write (*, '(2x,i0,2x,5f10.2)') i, ht(1), svib(1), srot(1), stra(1), gt(1)
      write (ich33, '(2x,i0,2x,5f10.2)') i, ht(1), svib(1), srot(1), stra(1), gt(1)
      solu_ens%gt(i) = gt(1)
      solu_ens%ht(i) = ht(1)
      solu_ens%svib(i) = svib(1)
      solu_ens%srot(i) = srot(1)
      solu_ens%stra(i) = stra(1)

      call chdir(tmppath2)
      call chdir('tmp_solu')
   end do
   close (ich33)

!--- Solvent cluster
   write (*, *)
   write (*, *) '  Solvent cluster properties'
   call chdir(tmppath2)
   open (newunit=ich81, file='solvent_cluster.dat')

   call chdir('tmp_solv')

   allocate (solv_ens%gt(solv_ens%nall))
   allocate (solv_ens%ht(solv_ens%nall))
   allocate (solv_ens%svib(solv_ens%nall))
   allocate (solv_ens%srot(solv_ens%nall))
   allocate (solv_ens%stra(solv_ens%nall))

   call pr_freq_energy()
   call pr_freq_file(ich81)

   do i = 1, solv_ens%nall
      write (to, '("TMPFREQ",i0)') i
      call chdir(to)
      call rdtherm('xtb_freq.out', ht(2), svib(2), srot(2), stra(2), gt(2))
      write (*, '(2x,i0,2x,5f10.2)') i, ht(2), svib(2), srot(2), stra(2), gt(2)
      write (ich81, '(2x,i0,2x,5f10.2)') i, ht(2), svib(2), srot(2), stra(2), gt(2)
      solv_ens%gt(i) = gt(2)
      solv_ens%ht(i) = ht(2)
      solv_ens%svib(i) = svib(2)
      solv_ens%srot(i) = srot(2)
      solv_ens%stra(i) = stra(2)
      call chdir(tmppath2)
      call chdir('tmp_solv')
   end do
   close (ich81)

!--- Saving results
   call chdir(tmppath2)
   call copysub('solute.dat', resultspath)
   call copysub('solute_cluster.dat', resultspath)
   call copysub('solvent_cluster.dat', resultspath)

!--- Deleting tmp directory
   call chdir(tmppath)
   if (.not. env%keepModef) call rmrf(tmppath2)
   call chdir(thispath)

   env%gfnver = gfnver_tmp
   env%optlev = optlev_tmp

   call tim%stop(9)

end subroutine qcg_freq

subroutine qcg_eval(env, solu, solu_ens, solv_ens)
   use iso_fortran_env, wp => real64
   use crest_data
   use iomod
   use zdata
   use strucrd

   implicit none

   type(systemdata)           :: env
   type(zmolecule)            :: solu
   type(ensemble)             :: solu_ens, solv_ens

   character(len=512)         :: thispath

   integer                    :: i, j
   integer                    :: srange
   integer                    :: freqscal
   real(wp)                   :: g1(solu_ens%nall)
   real(wp)                   :: g2(solv_ens%nall)
   real(wp)                   :: g3
   real(wp)                   :: Gsolv(20)
   real(wp)                   :: Hsolv
   real(wp)                   :: G_solute(20)
   real(wp)                   :: H_solute
   real(wp)                   :: G_solvent(20)
   real(wp)                   :: H_solvent
   real(wp)                   :: G_mono(20)
   real(wp)                   :: H_mono
   real(wp)                   :: S(20)
   real(wp)                   :: volw
   real(wp)                   :: sasa
   real(wp)                   :: dum, dum1, dum2
   real(wp)                   :: e_solute(solu_ens%nall)
   real(wp)                   :: e_solvent(solv_ens%nall)
   real(wp)                   :: scal(20)
   integer                    :: ich23
   real(wp), parameter         :: eh = 627.509541d0

   interface
      subroutine aver(pr, env, runs, e_tot, S, H, G, sasa, a_present, a_tot)
         use iso_fortran_env, only: wp => real64
         use crest_data

         implicit none
         type(systemdata), intent(in)     :: env
         integer, intent(in)             :: runs
         real(wp), intent(inout)         :: e_tot
         real(wp), intent(in), optional  :: a_tot
         real(wp), intent(out)           :: S
         real(wp), intent(out)           :: H
         real(wp), intent(out)           :: G
         real(wp), intent(out)           :: sasa
         logical, intent(in)             :: pr, a_present
         dimension e_tot(runs)
         dimension a_tot(runs)
      end subroutine aver
   end interface

   call pr_eval_eval()

   call getcwd(thispath)

   freqscal = nint(env%freq_scal/0.05)
   srange = 20
   do i = 1, srange
      scal(i) = 0.05*i
   end do

!--- Solute Cluster
   !H_solv
   do i = 1, solu_ens%nall
      e_solute(i) = solu_ens%er(i)*eh + solu_ens%ht(i)
   end do
   call aver(.false., env, solu_ens%nall, e_solute, dum1, H_solute, dum2, sasa, .false.)
   !G_solv
   do i = 1, srange
      do j = 1, solu_ens%nall
         g1(j) = solu_ens%ht(j) - (env%tboltz*(solu_ens%svib(j) + scal(i)*(solu_ens%srot(j) + solu_ens%stra(j)))/1000)
         e_solute(j) = solu_ens%er(j)*eh + g1(j)
      end do
      call aver(.false., env, solu_ens%nall, e_solute, S(i), dum, G_solute(i), sasa, .false.)
   end do

!--- Solvent Cluster
   !H_solv
   do i = 1, solv_ens%nall
      e_solvent(i) = solv_ens%er(i)*eh + solv_ens%ht(i)
   end do
   call aver(.false., env, solv_ens%nall, e_solvent, dum1, H_solvent, dum2, sasa, .false.)

   !G_solv
   do i = 1, srange
      do j = 1, solv_ens%nall
         g2(j) = solv_ens%ht(j) - &
                 & (env%tboltz*(solv_ens%svib(j) + scal(i)*(solv_ens%srot(j) + solv_ens%stra(j)))/1000)
         e_solvent(j) = solv_ens%er(j)*eh + g2(j)
      end do
      call aver(.false., env, solv_ens%nall, e_solvent, S(i), dum, G_solvent(i), sasa, .false.)
   end do

!--- Solute gas phase
   H_mono = solu%energy*eh + solu%ht
   do i = 1, srange
      g3 = solu%ht - (env%tboltz*(solu%svib + scal(i)*(solu%srot + solu%stra))/1000)
      G_mono(i) = solu%energy*eh + g3
   end do

   Gsolv(1:20) = G_solute(1:20) - G_solvent(1:20) - G_mono(1:20)
   Hsolv = H_solute - H_solvent - H_mono

!--- Calculate Volume work and include
   volw = (env%tboltz*8.31451/1000./4.184)*log(24.47d0*env%tboltz/298.15)
   Gsolv(1:20) = Gsolv(1:20) - volw
   Hsolv = Hsolv - volw
   call pr_eval_1(Gsolv(20), Hsolv)
   call pr_eval_2(srange, Gsolv, scal)
   call pr_eval_3(srange, freqscal, env%freq_scal, Gsolv)

! Save Result
   open (newunit=ich23, file='frequencies/result.dat')
   write (ich23, '("Solvation Free Energy [kcal/mol] :")')
   write (ich23, '(f8.2)') Gsolv(freqscal)
   close (ich23)

end subroutine qcg_eval

subroutine write_qcg_setup(env)
   use crest_data
   use iomod
   implicit none

   type(systemdata) :: env

   write (*, *)
   write (*, '(2x,''========================================='')')
   write (*, '(2x,''|   quantum cluster growth: INPUT       |'')')
   write (*, '(2x,''========================================='')')
   write (*, *)
   select case (env%qcg_runtype)
   case (0)
      write (*, '(2x,''QCG: Only Cluster Generation'')')
   case (1)
      write (*, '(2x,''QCG: Cluster + Ensemble Generation'')')
      if (env%ensemble_method .eq. 0) write (*, '(2x,''Ensemble generated via CREST'')')
      if (env%ensemble_method .eq. 1) write (*, '(2x,''Ensemble generated via MD Simulation'')')
      if (env%ensemble_method .eq. 2) write (*, '(2x,''Ensemble generated via MetaDynamic'')')
   case (2)
      write (*, '(2x,''QCG: Calculation of delta E_solv'')')
      if (env%ensemble_method .eq. 0) write (*, '(2x,''Ensemble generated via CREST'')')
      if (env%ensemble_method .eq. 1) write (*, '(2x,''Ensemble generated via MD Simulation'')')
      if (env%ensemble_method .eq. 2) write (*, '(2x,''Ensemble generated via MetaDynamic'')')
   case (3)
      write (*, '(2x,''QCG: Calculation of delta G_solv'')')
      if (env%ensemble_method .eq. 0) write (*, '(2x,''Ensemble generated via CREST'')')
      if (env%ensemble_method .eq. 1) write (*, '(2x,''Ensemble generated via MD Simulation'')')
      if (env%ensemble_method .eq. 2) write (*, '(2x,''Ensemble generated via MetaDynamic'')')
   case default
      continue
   end select
   write (*, *)
   write (*, '(2x,''input parameters     '')')
   write (*, '(2x,''solute                 : '',a)') trim(env%solu_file)
   write (*, '(2x,''charge                 : '',i0)') env%chrg
   write (*, '(2x,''uhf                    : '',i0)') env%uhf
   write (*, '(2x,''solvent                : '',a)') trim(env%solv_file)
   if (env%nsolv .ne. 0) then
      write (*, '(2x,''# of solvents to add   : '',i0)') env%nsolv
   else if (env%nsolv .eq. 0) then
      write (*, '(2x,''# of solvents to add   : until convergence, but maximal'',1x,i4)') env%max_solv
   end if
   if (env%nqcgclust .ne. 0) then
      write (*, '(2x,''# of cluster generated : '',i0)') env%nqcgclust
   else
      write (*, '(2x,''Cluster generated that are above 10 % populated '')')
   end if

   write (*, '(2x,''# of CPUs used         : '',i0)') env%Threads
   if (env%solvent .eq. '') then
      write (*, '(2x,''No gbsa/alpb model''  )')
   else
      write (*, '(2x,''Solvation model        : '',a)') env%solvent
   end if
   write (*, '(2x,''xtb opt level          : '',a)') trim(optlevflag(env%optlev))
   write (*, '(2x,''System temperature [K] : '',F5.1)') env%tboltz
   write (*, '(2x,''RRHO scaling factor    : '',F4.2)') env%freq_scal
   write (*, *)
   if (env%use_xtbiff) write (*, '(2x,''Use of xTB-IFF standalone requested'')')

end subroutine write_qcg_setup

subroutine get_sphere(pr, zmol, r_logical)
   use crest_parameters, only : wp
   use zdata

   implicit none
   type(zmolecule), intent(inout) :: zmol
   type(zmolecule) :: dum
   logical        :: pr
   logical        :: r_logical !Determines wether r is overwritten or not
   real(wp), parameter :: pi43 = 3.1415926540d0*4.0d0/3.0d0
   real(wp), parameter :: pi = 3.1415926540d0
   real(wp), parameter :: third = 1.0d0/3.0d0
   real(wp), parameter :: bohr = 0.52917726d0

   integer :: i
   real*8 :: rad(zmol%nat), xyz_tmp(3, zmol%nat)
   real(wp), allocatable :: rcov(:)

   allocate (rcov(94))
   rcov = (/ &
       &  2.18230009, 1.73469996, 3.49559999, 3.09820008, 3.21600008, &
       &  2.91030002, 2.62249994, 2.48169994, 2.29959989, 2.13739991, &
       &  3.70819998, 3.48390007, 4.01060009, 3.79169989, 3.50169992, &
       &  3.31069994, 3.10459995, 2.91479993, 4.24109983, 4.10349989, &
       &  3.89030004, 3.76419997, 3.72110009, 3.44140005, 3.54620004, &
       &  3.44210005, 3.43269992, 3.34619999, 3.30080009, 3.23090005, &
       &  3.95790005, 3.86190009, 3.66249990, 3.52679992, 3.36619997, &
       &  3.20959997, 4.61759996, 4.47639990, 4.21960020, 4.05970001, &
       &  3.85960007, 3.75430012, 3.56900001, 3.46230006, 3.39750004, &
       &  3.35249996, 3.33080006, 3.46199989, 4.26230001, 4.18739986, &
       &  4.01499987, 3.89010000, 3.73799992, 3.58890009, 5.05670023, &
       &  5.18139982, 4.62610006, 4.62010002, 4.57019997, 4.52710009, &
       &  4.48960018, 4.45149994, 4.42339993, 4.12430000, 4.24270010, &
       &  4.15409994, 4.27939987, 4.24499989, 4.22079992, 4.19859982, &
       &  4.01300001, 4.24499989, 4.09800005, 3.98550010, 3.89549994, &
       &  3.74900007, 3.44560003, 3.35249996, 3.25640011, 3.35990000, &
       &  4.31269979, 4.27640009, 4.11749983, 4.00540018, 3.86439991, &
       &  3.72160006, 5.07959986, 4.92939997, 4.70429993, 4.42519999, &
       &  4.45940018, 4.39569998, 4.35389996, 4.43410015/)

   do i = 1, zmol%nat
      rad(i) = bohr*rcov(zmol%at(i))*1.40 ! scale factor adjusted to rough
      xyz_tmp(1:3, i) = bohr*zmol%xyz(1:3, i)
   end do

   dum = zmol
   dum%xyz = xyz_tmp

   call get_volume(dum, rad)

   zmol%atot = dum%atot/bohr**2
   zmol%vtot = dum%vtot/bohr**3

   if (r_logical) then
      zmol%rtot = zmol%vtot*3.0/4.d0/pi
      zmol%rtot = zmol%rtot**(1.d0/3.d0)
   end if

   if (pr) then
      if (r_logical) then
         write (*, '(2x,''molecular radius (Bohr**1):'',F8.2)') zmol%rtot
      end if
      write (*, '(2x,''molecular area   (Bohr**2):'',F8.2)') zmol%atot
      write (*, '(2x,''molecular volume (Bohr**3):'',F8.2)') zmol%vtot
   end if

   deallocate (rcov)
end subroutine get_sphere

subroutine cma_shifting(solu, solv)
   use iso_fortran_env, wp => real64
   use crest_data
   use iomod
   use zdata
   use strucrd
   use axis_module, only: cma
   implicit none

   type(zmolecule)    :: solu, solv

   integer            :: i

   call cma(solu%nat, solu%at, solu%xyz, solu%cma)
   call cma(solv%nat, solv%at, solv%xyz, solv%cma)

   do i = 1, solu%nat
      solu%xyz(1:3, i) = solu%xyz(1:3, i) - solu%cma(1:3)
   end do
   do i = 1, solv%nat
      solv%xyz(1:3, i) = solv%xyz(1:3, i) - solv%cma(1:3)
   end do

end subroutine cma_shifting

subroutine get_ellipsoid(env, solu, solv, clus, pr1)
   use iso_fortran_env, wp => real64
   use crest_data
   use iomod
   use zdata
   use strucrd
   use axis_module
   implicit none

   type(systemdata)   :: env
   type(zmolecule)    :: solu, solv, clus
   type(zmolecule)    :: dummy_solu, dummy_solv
   real(wp)           :: rabc_solu(3), rabc_solv(3)
   real(wp)           :: aniso, sola
   real(wp)           :: rmax_solu, rmax_solv
   real(wp)           :: boxr, roff, r
   character(len=10) :: fname
   logical            :: ex, pr, pr1

   real(wp), parameter :: pi43 = 3.1415926540d0*4.0d0/3.0d0
   real(wp), parameter :: pi = 3.1415926540d0
   real(wp), parameter :: third = 1.0d0/3.0d0

   pr = .false. !Outprint deactivated

   fname = 'eaxis.qcg'
   inquire (file=fname, exist=ex)

   if (pr1) then !First time called
!--- Moving all coords to the origin (transformation)
      call axistrf(solu%nat, solu%nat, solu%at, solu%xyz)
!  call axistrf(solv%nat,solv%nat,solv%at,solv%xyz)  !Not done in original QCG code
      call axistrf(clus%nat, solu%nat, clus%at, clus%xyz)

!--- Overwrite solute and solvent coord in original file with transformed and optimized ones
      call wrc0('solute', solu%nat, solu%at, solu%xyz)
      call wrc0('solvent', solv%nat, solv%at, solv%xyz)

!--- Getting axis
      write (*, *) 'Solute:'
      call axis(pr1, solu%nat, solu%at, solu%xyz, solu%eax)
      write (*, *) 'Solvent:'
      call axis(pr1, solv%nat, solv%at, solv%xyz, solv%eax)
      write (*, *)
   end if

!--- Computing anisotropy factor of solute and solvent
   sola = sqrt(1.+(solu%eax(1) - solu%eax(3))/((solu%eax(1) + solu%eax(2) + solu%eax(3))/3.))
   aniso = sqrt(1.+(solv%eax(1) - solv%eax(3))/((solv%eax(1) + solv%eax(2) + solv%eax(3))/3.)) ! =1 for a spherical system

!--- Get maximum intramoleclar distance of solute and solvent
   call getmaxrad(solu%nat, solu%at, solu%xyz, rmax_solu)
   call getmaxrad(solv%nat, solv%at, solv%xyz, rmax_solv)

!--- Getting V and A of dummies
   dummy_solu = solu
   dummy_solv = solv !Why is dummy_solv%vtot different to solv%vtot
   call get_sphere(.false., dummy_solu, .false.)
   call get_sphere(.false., dummy_solv, .false.)

!--- Computation of outer Wall
   roff = sola*dummy_solu%vtot/1000
   boxr = ((0.5*aniso*clus%nmol*dummy_solv%vtot + dummy_solu%vtot)/pi43)**third + roff + rmax_solv*0.5 !0.5 both
   r = (boxr**3/(solu%eax(1)*solu%eax(2)*solu%eax(3)))**third       ! volume of ellipsoid = volume of sphere
   rabc_solv = solu%eax*r                              ! outer solvent wall

!--- Computation of inner wall
   roff = sola*dummy_solu%vtot/1000
   boxr = ((sola*dummy_solu%vtot)/pi43)**third + roff + rmax_solu*0.1 !0.1 before
   r = (boxr**3/(solu%eax(1)*solu%eax(2)*solu%eax(3)))**third       ! volume of ellipsoid = volume of sphere
   rabc_solu = solu%eax*r
   dummy_solu%ell_abc(1) = solu%eax(1)**2/sum((solu%eax(1:3))**2)
   dummy_solu%ell_abc(2) = solu%eax(2)**2/sum((solu%eax(1:3))**2)
   dummy_solu%ell_abc(3) = solu%eax(3)**2/sum((solu%eax(1:3))**2)
   rabc_solu = dummy_solu%ell_abc*r

   solu%aniso = sola
   solv%aniso = aniso
   solu%ell_abc = rabc_solu
   clus%ell_abc = rabc_solv*env%potscal

   if (pr1) then
      write (*, '(2x,''solvent anisotropy            :'',4f10.3)') aniso
      write (*, '(2x,''solute anisotropy             :'',4f10.3)') sola
      write (*, '(2x,''roff inner wall               :'',4f10.3)') roff
      write (*, '(2x,''solute max dist               :'',4f10.3)') rmax_solu
      write (*, '(2x,''solvent max dist              :'',4f10.3)') rmax_solv
      write (*, '(2x,''inner unit axis               :'',3f10.3)') dummy_solu%ell_abc(1:3)
      write (*, '(2x,''inner ellipsoid/Bohr          :'',3f10.3)') rabc_solu(1:3)
      write (*, '(2x,''scaling factor outer ellipsoid:'',3f10.3)') env%potscal
      write (*, '(2x,''outer ellipsoid/Bohr          :'',3f10.3)') clus%ell_abc(1:3)
      if (env%potscal .gt. 1.0_wp) write &
           &(*, '(2x,''!!!WARNING: A SCALING FACTOR LARGER 1.0 IS ONLY FOR MICROSOLVATION RECOMMENDED'')')
      write (*, *)
   end if

end subroutine get_ellipsoid

subroutine getmaxrad(n, at, xyz, r)
   use crest_parameters, only : wp
   implicit none
   real(wp) :: xyz(3, n), r
   integer :: n, at(n)

   real(wp) :: rx, ry, rz, rr
   integer :: i, j
   real(wp), allocatable :: rcov(:)

   allocate (rcov(94))
   rcov = (/ &
    &  2.18230009, 1.73469996, 3.49559999, 3.09820008, 3.21600008, &
    &  2.91030002, 2.62249994, 2.48169994, 2.29959989, 2.13739991, &
    &  3.70819998, 3.48390007, 4.01060009, 3.79169989, 3.50169992, &
    &  3.31069994, 3.10459995, 2.91479993, 4.24109983, 4.10349989, &
    &  3.89030004, 3.76419997, 3.72110009, 3.44140005, 3.54620004, &
    &  3.44210005, 3.43269992, 3.34619999, 3.30080009, 3.23090005, &
    &  3.95790005, 3.86190009, 3.66249990, 3.52679992, 3.36619997, &
    &  3.20959997, 4.61759996, 4.47639990, 4.21960020, 4.05970001, &
    &  3.85960007, 3.75430012, 3.56900001, 3.46230006, 3.39750004, &
    &  3.35249996, 3.33080006, 3.46199989, 4.26230001, 4.18739986, &
    &  4.01499987, 3.89010000, 3.73799992, 3.58890009, 5.05670023, &
    &  5.18139982, 4.62610006, 4.62010002, 4.57019997, 4.52710009, &
    &  4.48960018, 4.45149994, 4.42339993, 4.12430000, 4.24270010, &
    &  4.15409994, 4.27939987, 4.24499989, 4.22079992, 4.19859982, &
    &  4.01300001, 4.24499989, 4.09800005, 3.98550010, 3.89549994, &
    &  3.74900007, 3.44560003, 3.35249996, 3.25640011, 3.35990000, &
    &  4.31269979, 4.27640009, 4.11749983, 4.00540018, 3.86439991, &
    &  3.72160006, 5.07959986, 4.92939997, 4.70429993, 4.42519999, &
    &  4.45940018, 4.39569998, 4.35389996, 4.43410015/)

   r = 0
   do i = 1, n - 1
      do j = i + 1, n
         rx = xyz(1, i) - xyz(1, j)
         ry = xyz(2, i) - xyz(2, j)
         rz = xyz(3, i) - xyz(3, j)
         rr = sqrt(rx**2 + ry**2 + rz**2) + rcov(at(i)) + rcov(at(j))
         if (rr .gt. r) r = rr
      end do
   end do

   deallocate (rcov)

end subroutine getmaxrad

subroutine ellipsout(fname, n, at, xyz, r1)
   use iso_fortran_env, only: wp => real64
   use strucrd, only: i2e
   implicit none

   integer            :: i
   integer            :: n, at(n)
   real(wp)           :: xyz(3, n), r1(3)
   real(wp)           :: x, y, z, f, rr
   character(len=*)   :: fname
   integer            :: ich11

   open (newunit=ich11, file=fname)
   write (ich11, '(a)') '$coord'
   do i = 1, n
      write (ich11, '(3F24.14,6x,a)') xyz(1, i), xyz(2, i), xyz(3, i), i2e(at(i))
   end do
   do i = 1, 500
      call random_number(x)
      call random_number(f)
      if (f .gt. 0.5) x = -x
      call random_number(y)
      call random_number(f)
      if (f .gt. 0.5) y = -y
      call random_number(z)
      call random_number(f)
      if (f .gt. 0.5) z = -z
      rr = sqrt(x*x + y*y + z*z)
      x = x*r1(1)/rr
      y = y*r1(2)/rr
      z = z*r1(3)/rr
      write (ich11, '(3F24.14,6x,a2)') x, y, z, 'he'
   end do
   write (ich11, '(a)') '$end'
   close (ich11)

end subroutine ellipsout

subroutine both_ellipsout(fname, n, at, xyz, r1, r2)
   use iso_fortran_env, only: wp => real64
   use strucrd, only: i2e
   implicit none

   integer            :: i
   integer            :: n, at(n)
   real(wp)           :: xyz(3, n), r1(3)
   real(wp), optional :: r2(3)
   real(wp)           :: x, y, z, f, rr
   character(len=*)   :: fname
   integer            :: ich11

   open (newunit=ich11, file=fname)
   write (ich11, '(a)') '$coord'
   do i = 1, n
      write (ich11, '(3F24.14,6x,a)') xyz(1, i), xyz(2, i), xyz(3, i), i2e(at(i))
   end do
   do i = 1, 500
      call random_number(x)
      call random_number(f)
      if (f .gt. 0.5) x = -x
      call random_number(y)
      call random_number(f)
      if (f .gt. 0.5) y = -y
      call random_number(z)
      call random_number(f)
      if (f .gt. 0.5) z = -z
      rr = sqrt(x*x + y*y + z*z)
      x = x*r1(1)/rr
      y = y*r1(2)/rr
      z = z*r1(3)/rr
      write (ich11, '(3F24.14,6x,a2)') x, y, z, 'he'
   end do
   if (present(r2)) then
      do i = 1, 100
         call random_number(x)
         call random_number(f)
         if (f .gt. 0.5) x = -x
         call random_number(y)
         call random_number(f)
         if (f .gt. 0.5) y = -y
         call random_number(z)
         call random_number(f)
         if (f .gt. 0.5) z = -z
         rr = sqrt(x*x + y*y + z*z)
         x = x*r2(1)/rr
         y = y*r2(2)/rr
         z = z*r2(3)/rr
         write (ich11, '(3F24.14,6x,a2)') x, y, z, 'b'
      end do
   end if
   write (ich11, '(a)') '$end'
   close (ich11)

end subroutine both_ellipsout

subroutine get_interaction_E(env, solu, solv, clus, iter, E_inter)
   use iso_fortran_env, wp => real64
   use crest_data
   use iomod
   use zdata
   use strucrd
   implicit none

   type(systemdata)            :: env
   type(zmolecule), intent(in) :: solu, solv, clus
   real(wp)                    :: e_cluster, e_solute, e_solvent
   real(wp)                    :: E_inter(env%nsolv)           ! interaction energy
   integer                     :: iter
   logical                     :: e_there

   call remove('cluster.coord')

!--- Prepare input coordinate files
   call wrc0('cluster.coord', clus%nat, clus%at, clus%xyz)
   call wr_cluster_cut('cluster.coord', solu%nat, solv%nat, iter, 'solute_cut.coord', 'solvent_cut.coord')

!--- Perform single point calculations and recieve energies
   call xtb_sp_qcg(env, 'solute_cut.coord')
   call grepval('xtb.out', '| TOTAL ENERGY', e_there, e_solute)
   if (.not. e_there) write (*, *) 'Solute energy not found'
   call xtb_sp_qcg(env, 'solvent_cut.coord')
   call grepval('xtb.out', '| TOTAL ENERGY', e_there, e_solvent)
   if (.not. e_there) write (*, *) 'Solvent energy not found'
   call xtb_sp_qcg(env, 'cluster.coord')
   call grepval('xtb.out', '| TOTAL ENERGY', e_there, e_cluster)
   if (.not. e_there) write (*, *) 'Cluster energy not found'

   E_inter(iter) = e_cluster - e_solute - e_solvent

end subroutine get_interaction_E

subroutine analyze_cluster(nsolv, n, nS, nM, xyz, at, av, last)
   use iso_fortran_env, only: wp => real64
   use axis_module, only: cma
   implicit none
   real(wp) xyz(3, n)
   real(wp) av, last
   integer n, nS, nM, nsolv, at(n)
   real(wp) xyzM(3, nM)
   integer atm(nM)
   real(wp) xyzS(3, nS)
   integer atS(nS)
   real(wp) x1(3), x2(3), r
   integer i, is, ie

   if (nsolv .eq. 1) return
   xyzS(1:3, 1:nS) = xyz(1:3, 1:nS)
   atS(1:nS) = at(1:nS)
   call cma(nS, atS, xyzS, x1)

   av = 0
   do i = 1, nsolv
      is = nS + (i - 1)*nM + 1
      ie = is + nM - 1
      xyzM(1:3, 1:nM) = xyz(1:3, is:ie)
      atM(1:nM) = at(is:ie)
      call cma(nM, atM, xyzM, x2)
      r = sqrt((x1(1) - x2(1))**2 + (x1(2) - x2(2))**2 + (x1(3) - x2(3))**2)
      if (i .lt. nsolv) then
         av = av + r
      else
         last = r
      end if
   end do
   av = av/float(nsolv - 1)
end subroutine analyze_cluster

subroutine aver(pr, env, runs, e_tot, S, H, G, sasa, a_present, a_tot)
   use iso_fortran_env, only: wp => real64
   use crest_data

   implicit none
!---- Dummy
   type(systemdata), intent(in)     :: env
   integer, intent(in)             :: runs
   real(wp), intent(inout)         :: e_tot
   real(wp), intent(in), optional  :: a_tot
   real(wp), intent(out)           :: S
   real(wp), intent(out)           :: H
   real(wp), intent(out)           :: G
   real(wp), intent(out)           :: sasa
!---- Stack
   logical, intent(in)             :: pr, a_present
   integer                         :: j, jmin
   real(wp)                        :: A
   real(wp)                        :: e0
   real(wp), allocatable           :: de(:)
   real(wp), allocatable           :: p(:)
   real(wp)                        :: pmax
   real(wp)                        :: eav
   real(wp)                        :: area
   real(wp)                        :: beta
   real(wp)                        :: temp
   integer                         :: ich48
   real(wp), parameter              :: eh = 627.509541d0
   dimension e_tot(runs)
   dimension a_tot(runs)

   temp = env%tboltz
   allocate (de(runs), source=0.0d0)
   allocate (p(runs), source=0.0d0)

   beta = 1./(temp*8.314510/4.184/1000.+1.d-14)
   e0 = e_tot(1)
   de(1:runs) = (e_tot(1:runs) - e0)
   call qcg_boltz(env, runs, de, p)

   A = 0
   eav = 0
   pmax = 0
   area = 0
   do j = 1, runs
      A = A + p(j)*log(p(j) + 1.d-12)
      eav = eav + p(j)*e_tot(j)
      if (p(j) .gt. pmax) then
         pmax = p(j)
         jmin = j
      end if
      if (a_present) area = area + p(j)*a_tot(j)
   end do
   sasa = area
   S = (1./beta)*A
   H = eav
   G = eav + S
   if (pr) then
      open (newunit=ich48, file='population.dat')
      write (ich48, '(2x, ''cluster'',2x,''E_norm [Eh]'',2x, ''De [kcal]'', 4x, ''p'')')
      do j = 1, runs
         if (j .lt. 10) then
            write (ich48, '(5x,i0,3x,f11.6,5x,f6.4,3x,f6.4)') j, e_tot(j)/eh, de(j), p(j)
         else
            write (ich48, '(5x,i0,2x,f11.6,5x,f6.4,3x,f6.4)') j, e_tot(j)/eh, de(j), p(j)
         end if
      end do
      write (ich48, *)
      write (ich48, '(''Ensemble free energy [Eh]:'', f20.10)') G/eh
      close (ich48)
   end if

   deallocate (de, p)

end subroutine aver

subroutine qcg_boltz(env, n, e, p)
   use iso_fortran_env, only: wp => real64
   use crest_data
   implicit none
   type(systemdata), intent(in)    :: env
   integer, intent(in)              :: n
   real(wp), intent(in)             :: e(*)
   real(wp), intent(out)            :: p(*)
   integer                         :: i
   real(wp)                        :: temp
   real(wp)                        :: f, hsum, esum

   temp = env%tboltz
   f = 8.314*temp/4.184d+3
   esum = 0
   do i = 1, n
      esum = esum + exp(-e(i)/f)
   end do
   hsum = 0
   do i = 1, n
      p(i) = exp(-e(i)/f)/esum
   end do
end subroutine qcg_boltz

subroutine fill_take(env, n2, n12, rabc, ipos)
   use iso_fortran_env, only: wp => real64
   use crest_data
   use strucrd
   use axis_module, only: cma
   implicit none

   type(systemdata)      :: env
   integer, intent(in)   :: n2, n12
   real(wp), intent(in)   :: rabc(3)
   integer, intent(out)  :: ipos
   integer               :: i, m, n21
   integer               :: at2(n2), at12(n12)
   integer               :: counter
   real(wp)              :: xyz2(3, n2), xyz12(3, n12)
   real(wp)              :: etmp(100)
   real(wp)              :: eabc
   real(wp)              :: cma2(3)
   real(wp), allocatable  :: dist(:)

   eabc = 0
   counter = 0
   n21 = n12 - n2 + 1
   if (env%use_xtbiff) then
      call rdxtbiffE('xtbscreen.xyz', m, n12, etmp)
   else
      call rdxtbiffE('best.xyz', m, n12, etmp)
   end if

   allocate (dist(m), source=0.0d0)
   dist = 0.0d0

   do i = 1, m
      if (env%use_xtbiff) then
         call rdxmolselec('xtbscreen.xyz', i, n12, at12, xyz12)
      else
         call rdxmolselec('final_structures.xyz', i, n12, at12, xyz12)
      end if

      at2(1:n2) = at12(n21:n12)
      xyz2(1:3, 1:n2) = xyz12(1:3, n21:n12)
      call cma(n2, at2, xyz2, cma2)
      call calc_dist(cma2, rabc, dist(i), eabc)
      if (eabc .gt. 1.0d0) then
         dist(i) = 1.0d42
         counter = counter + 1
      end if
   end do

   ipos = minloc(dist(1:m), dim=1)

   if (counter .eq. m) ipos = 0

   deallocate (dist)
end subroutine fill_take

subroutine calc_dist(xyz, rabc, dist, eabc)
   use iso_fortran_env, only: wp => real64
   implicit none

   real(wp), intent(in)    :: xyz(3)
   real(wp), intent(in)    :: rabc(3)
   real(wp), intent(out)   :: dist
   real(wp), intent(out)   :: eabc
   real(wp)                :: center(3), rc(3)

   center = 0.d0
   rc = (xyz(1:3) - center)
   dist = norm2(rc)
   eabc = sum((xyz(1:3)**2)/(rabc(1:3)**2))
end subroutine calc_dist

subroutine sort_min(i, j, col, A)
   use iso_fortran_env, only: wp => real64
   implicit none
   integer, intent(in)   :: i, j, col
   real*8, intent(inout) :: A(i, j)
   real*8                :: buf(j)
   integer               :: nsize, irow, krow
! dimension A(i,j)
   nsize = i

   do irow = 1, nsize
      krow = minloc(A(irow:nsize, col), dim=1) + irow - 1
      buf(:) = A(irow, :)
      A(irow, :) = A(krow, :)
      A(krow, :) = buf(:)
   end do
end subroutine sort_min

subroutine sort_ensemble(ens, e_ens, fname)
   use iso_fortran_env, only: wp => real64
   use crest_data
   use strucrd
   implicit none
   type(ensemble)       :: ens
   real(wp)             :: e_ens(ens%nall), dum(ens%nall)
   character(len=*)     :: fname
   integer              :: ich
   integer              :: i, e_min

   dum = e_ens

   open (newunit=ich, file=fname)

   do i = 1, ens%nall
      e_min = minloc(dum, dim=1)
      call wrxyz(ich, ens%nat, ens%at, ens%xyz(:, :, e_min), e_ens(e_min))
      dum(e_min) = 0.0d0
   end do
   close (ich)

end subroutine sort_ensemble

subroutine rdtherm(fname, ht, svib, srot, stra, gt)
   use iso_fortran_env, only: wp => real64
   use crest_data
   use iomod

   implicit none
! Dummy
   real(wp), intent(out)  :: ht
   real(wp), intent(out)  :: gt
   real(wp), intent(out)  :: svib
   real(wp), intent(out)  :: srot
   real(wp), intent(out)  :: stra
! Stack
   integer                :: nn
   integer                :: io
   integer                :: counter
   integer                :: hg_line
   real(wp)               :: xx(20)
   logical                :: ende
   character(len=*)       :: fname
   character(len=128)     :: a
   real(wp), parameter     :: eh = 627.509541d0
   integer                :: ich

   ende = .false.
   counter = 0
   hg_line = 0

   open (newunit=ich, file=fname)
   do while (.not. ende)
      read (ich, '(a)', iostat=io) a
      if (io .lt. 0) then
         ende = .true.
         cycle
      end if
      if (index(a, 'G(T)/Eh ') .ne. 0) then
         hg_line = counter
      end if
      if (index(a, '  VIB  ') .ne. 0) then
         call readl(a, xx, nn)
         svib = xx(5)
         if (svib .eq. 0.0d0) then
            call readl(a, xx, nn)
            svib = xx(4)
         end if
      end if
      if (index(a, '  ROT  ') .ne. 0) then
         call readl(a, xx, nn)
         srot = xx(4)
      end if
      if (index(a, '  TR   ') .ne. 0) then
         call readl(a, xx, nn)
         stra = xx(4)
      end if
      if (counter .eq. hg_line + 2) then
         call readl(a, xx, nn)
         ht = xx(3)*eh
         gt = xx(5)*eh
      end if
      counter = counter + 1
   end do
   close (ich)
end subroutine rdtherm

subroutine pr_freq_file(ich)
   implicit none
   integer :: ich
   write (ich, '(2x,"#       H(T)       SVIB      SROT       STRA      G(T)")')
   write (ich, '(2x,"     [kcal/mol]    [      cal/mol/K        ]    [kcal/mol]")')
   write (ich, '(2x,"--------------------------------------------------------")')
end subroutine pr_freq_file

subroutine qcg_restart(env, progress, solu, solv, clus, solu_ens, solv_ens, clus_backup)
   use iso_fortran_env, wp => real64
   use crest_data
   use iomod
   use zdata
   use strucrd

   implicit none

   type(systemdata)           :: env
   type(zmolecule)            :: solu, solv, clus, clus_backup
   type(ensemble)             :: solu_ens, solv_ens
   integer                    :: progress

   integer                    :: i
   character(len=512)         :: thispath
   character(len=6)           :: counter
   character(len=7)           :: counter2
   character(len=8)           :: counter3
   logical                    :: grow, solu_ensemble, solv_ensemble
   logical                    :: solv_cff, solv_present, freq, tmp, ex
   real(wp), allocatable       :: xyz(:, :)
   real(wp), parameter         :: eh = 627.509541d0

   grow = .false.
   solu_ensemble = .false.
   solv_ensemble = .false.
   solv_cff = .false.
   solv_present = .false.
   freq = .false.
   tmp = .false.

   inquire (file='./grow/cluster.coord', exist=grow)
   inquire (file='./ensemble/final_ensemble.xyz', exist=solu_ensemble)
   inquire (file='./solvent_ensemble/final_ensemble.xyz', exist=solv_ensemble)
   inquire (file='./solvent_ensemble/crest_ensemble.xyz', exist=solv_cff)
   inquire (file='./frequencies/result.dat', exist=freq)

   if (solv_cff .or. solv_ensemble) solv_present = .true.

   call getcwd(thispath)

!---------------------------------------------------------------------------------
!        Check, if everything needed is present
!---------------------------------------------------------------------------------

   if (freq .and. ((.not. grow) .or. (.not. solu_ensemble) .or. (.not. solv_ensemble))) then
      progress = 0
      call rmrf('frequencies')
      freq = .false.
   end if

   if (solv_present .and. ((.not. grow) .or. (.not. solu_ensemble))) then
      progress = 0
      call rmrf('solvent_ensemble')
      solv_present = .false.
      solv_cff = .false.
      solv_ensemble = .false.
   end if

   if (solu_ensemble .and. (.not. grow)) then
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
      call chdir('grow')
      call rdnat('cluster.coord', clus%nat)
      allocate (clus%at(clus%nat), clus%xyz(3, clus%nat))
      call rdcoord('cluster.coord', clus%nat, clus%at, clus%xyz)
      clus%nmol = (clus%nat - solu%nat)/solv%nat + 1
      allocate (xyz(3, clus%nat))
      xyz = clus%xyz
      call get_ellipsoid(env, solu, solv, clus, .true.)
      clus%xyz = xyz !Needed, because get_ellipsoid performs axistransformation and not fitting potential
      deallocate (xyz)

      if (.not. env%cff) then
         allocate (clus_backup%at(clus%nat))
         allocate (clus_backup%xyz(3, clus%nat))
         clus_backup = clus
      end if

      if (clus%nmol - 1 .ge. env%nsolv) then
         progress = 1
         env%nsolv = clus%nmol - 1
         write (*, *)
         write (*, *)
         write (*, '(''Found cluster with '',i0,'' solvents'')') env%nsolv
         call chdir(thispath)
      else
         error stop 'The found cluster is smaller than nsolv. Please restart the whole computaion by removing the grow directory'
         !Future implementation continue grow process
         call chdir(thispath)
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
      call chdir('ensemble')
      call solu_ens%open('final_ensemble.xyz')
      call rdensemble('final_ensemble.xyz', solu_ens%nat, solu_ens%nall, solu_ens%at, solu_ens%xyz, solu_ens%er)
      env%nqcgclust = solu_ens%nall
      write (*, '("  Ensemble of solute-cluster found.")')
      write (*, '("  Taking all ", i0, " structures")') env%nqcgclust
      call grepval('population.dat', 'Ensemble free energy [Eh]:', ex, solu_ens%G)
      solu_ens%G = solu_ens%G*eh
      write (*, *) 'Solute Ensmeble Free E [kcal/mol]', solu_ens%G
      call chdir(thispath)
      progress = 2
   end if

!--- Solvent Ensemble
   if (solv_present) then
      call chdir('solvent_ensemble')
      write (*, '("  Ensemble of solvent-cluster found.")')

      !--- Case CFF
      if (solv_cff) then
         call solv_ens%open('crest_ensemble.xyz')
         do i = 1, solv_ens%nall
            if (i .le. 9) then
               write (counter, '(''No   '',i1)') i
               call grepval('cluster_energy.dat', counter, ex, solv_ens%er(i))
            else if (i .le. 99) then
               write (counter2, '(''No   '',i2)') i
               call grepval('cluster_energy.dat', counter2, ex, solv_ens%er(i))
            else
               write (counter3, '(''No   '',i3)') i
               call grepval('cluster_energy.dat', counter3, ex, solv_ens%er(i))
            end if
            write (*, *) 'Energy of cluster', i, solv_ens%er(i)
         end do
      end if

      !--- Case MD/Crest run
      if (solv_ensemble) then
         call solv_ens%open('final_ensemble.xyz')
         call rdensemble('final_ensemble.xyz', solv_ens%nat, solv_ens%nall, solv_ens%at, solv_ens%xyz, solv_ens%er)
      end if
      call grepval('population.dat', 'Ensemble free energy [Eh]:', ex, solv_ens%G)
      solv_ens%G = solv_ens%G*eh
      write (*, *) 'solvent ensmeble free E [kcal/mol]', solv_ens%G
      call chdir(thispath)
      progress = 3
   end if

!--- Frequencies
   if (freq) then
      write (*, *)
      write (*, *)
      write (*, *) '  Nothing to do'
      progress = 4
   end if

end subroutine qcg_restart

subroutine qcg_cleanup(env)
   use crest_data

   implicit none

   type(systemdata)      :: env
   character(len=280)    :: thispath
   logical               :: tmp

   call getcwd(thispath)
   call chdir(env%scratchdir)
   inquire (file='./solute_properties/solute', exist=tmp)
   if (tmp) then
      call rmrf('solute_properties')
      call rmrf('solvent_properties')
   end if

end subroutine qcg_cleanup

subroutine write_reference(env, solu, clus)
   use iso_fortran_env, wp => real64
   use crest_data
   use zdata, only: zmolecule
   use iomod
   use strucrd

   implicit none
   type(systemdata):: env    ! MAIN STORAGE OS SYSTEM DATA
   type(zmolecule)            :: solu, clus
   type(zmolecule)            :: ref_mol, ref_clus

   ref_mol = solu
   call rdcoord(env%solu_file, ref_mol%nat, ref_mol%at, ref_mol%xyz) !original solute coordinates
   call remove(env%fixfile)
   ref_clus = clus
   ref_clus%xyz(1:3, 1:solu%nat) = solu%xyz
   call wrc0(env%fixfile, ref_clus%nat, ref_clus%at, ref_clus%xyz)

end subroutine write_reference


!========================================================================================!
!> Convert given QCG coordinate files into (TM format)
!> Write "solute" and "solvent" coordinate files
!========================================================================================!
subroutine inputcoords_qcg(env, solute, solvent)
  use iso_fortran_env,only:wp => real64
  use crest_data
  use strucrd
  use zdata
  use iomod
  implicit none

  type(systemdata), intent(inout) :: env
  type(zmolecule), intent(out) :: solute, solvent

  logical :: ex11,ex21,solu,solv
  type(coord) :: mol
  type(zmolecule) :: zmol,zmol1
  integer :: i

!--------------------Checking for input-------------!

  !Solute
  inquire (file=env%solu_file,exist=ex11)
  inquire (file='solute',exist=solu)
  if (solu) call copy('solute','solute.old') !Backup solute file
  if ((.not. ex11) .and. (.not. solu)) then
    error stop 'No (valid) solute file! exit.'
  else if ((.not. ex11) .and. (solu)) then
    env%solu_file = 'solute'
  end if

  !Solvent
  inquire (file=env%solv_file,exist=ex21)
  inquire (file='solvent',exist=solv)
  if (solu) call copy('solvent','solvent.old') !Backup solvent file
  if ((.not. ex21) .and. (.not. solv)) then
    error stop 'No (valid) solvent file! exit.'
  else if ((.not. ex11) .and. (solu)) then
    env%solu_file = 'solvent'
  end if

!---------------Handling solute---------------------!
  call mol%open(env%solu_file)
  call mol%write('solute')
  solute%nat = mol%nat
  solute%at = mol%at
  solute%xyz = mol%xyz
  call mol%deallocate()

  !--- if the input was a SDF file, special handling
  env%sdfformat = .false.
  call checkcoordtype(env%solu_file,i)
  if (i == 31.or.i == 32) then
    !Add sdf stuff here, if somebody needs it
  end if

!---------------Handling solvent---------------------!

  call mol%open(env%solv_file)
  call mol%write('solvent')
  solvent%nat = mol%nat
  solvent%at = mol%at
  solvent%xyz = mol%xyz
  call mol%deallocate()

  !--- if the input was a SDF file, special handling
  env%sdfformat = .false.
  call checkcoordtype(env%solv_file,i)
  if (i == 31.or.i == 32) then
    !Add sdf stuff here, if somebody needs it
  end if

  return
end subroutine inputcoords_qcg
