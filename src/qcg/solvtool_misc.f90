!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 Christoph Plett, Sebastian Spicher, Philipp Pracht
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

!--------------------------------------------------------------------------------------------
! A quick single point xtb calculation without wbo
!--------------------------------------------------------------------------------------------
subroutine xtb_sp_qcg(env, fname)
   use iso_fortran_env, only: wp => real64
   use iomod
   use crest_data
   implicit none
   character(len=*) :: fname
   type(systemdata) :: env
   character(len=512) :: jobcall
   character(*), parameter :: pipe = ' > xtb.out 2> /dev/null'
   integer :: io
   call remove('gfnff_topo')
   call remove('energy')
   call remove('charges')
   call remove('xtbrestart')
!---- setting threads
   if (env%autothreads) then
      call ompautoset(env%threads, 7, env%omp, env%MAXRUN, 1) !set the global OMP/MKL variables for the xtb jobs
   end if
!---- jobcall
   write (jobcall, '(a,1x,a,1x,a,'' --sp '',a,1x,a)') &
   &     trim(env%ProgName), trim(fname), trim(env%gfnver), trim(env%solv), trim(pipe)

   call command(trim(jobcall), io)
!---- cleanup
   call remove('energy')
   call remove('charges')
   call remove('xtbrestart')
   call remove('xtbtopo.mol')
   call remove('gfnff_topo')
end subroutine xtb_sp_qcg

!--------------------------------------------------------------------------------------------
! A quick single xtb optimization gets zmol and overwrites it with optimized stuff
!--------------------------------------------------------------------------------------------
subroutine xtb_opt_qcg(env, zmol, constrain)
   use iso_fortran_env, only: wp => real64
   use iomod
   use crest_data
   use zdata
   use strucrd

   implicit none
   type(systemdata), intent(in) :: env
   type(zmolecule), intent(inout) :: zmol

   character(:), allocatable :: fname
   character(len=512) :: jobcall
   logical :: constrain
   logical :: const
   character(*), parameter :: pipe = ' > xtb_opt.out 2> /dev/null'
   integer :: io

!--- Write coordinated
   fname = 'coord'
   call wrc0(fname, zmol%nat, zmol%at, zmol%xyz) !write coord for xtbopt routine

!---- setting threads
   if (env%autothreads) then
      call ompautoset(env%threads, 7, env%omp, env%MAXRUN, 1) !set the global OMP/MKL variables for the xtb jobs
   end if
!---- jobcall & Handling constraints
   if(constrain .AND. env%cts%used) then
     call write_constraint(env, fname, 'xcontrol')
     call wrc0('coord.ref', zmol%nat, zmol%at, zmol%xyz) !write coord for xtbopt routine
     write (jobcall, '(a,1x,a,1x,a,'' --opt --input xcontrol '',a,1x,a)') &
     &     trim(env%ProgName), trim(fname), trim(env%gfnver), trim(env%solv), trim(pipe)
   else
     write (jobcall, '(a,1x,a,1x,a,'' --opt '',a,1x,a)') &
     &     trim(env%ProgName), trim(fname), trim(env%gfnver), trim(env%solv), trim(pipe)
   end if

   call command(trim(jobcall), io)
!---- cleanup
   call rdcoord('xtbopt.coord', zmol%nat, zmol%at, zmol%xyz)
   call remove('energy')
   call remove('charges')
   call remove('xtbrestart')
   call remove('xtbtopo.mol')
   call remove('gfnff_topo')
end subroutine xtb_opt_qcg

!___________________________________________________________________________________
!
! An xTB single point calculation and lmo generation on all available threads
!___________________________________________________________________________________

subroutine xtb_lmo(env, fname)!,chrg)
   use iso_fortran_env, only: wp => real64
   use iomod
   use crest_data
   use zdata
   implicit none
   type(systemdata) :: env
   character(len=*), intent(in)     :: fname
   character(len=80)               :: pipe
   character(len=512)              :: jobcall

   pipe = ' > xtb.out 2>/dev/null'

!---- setting threads
   if (env%autothreads) then
      call ompautoset(env%threads, 7, env%omp, env%MAXRUN, 1) !set the global OMP/MKL variables for the xtb jobs
   end if

!---- jobcall, special gbsa treatment not needed, as the entire flag is included in env%solv
   write (jobcall, '(a,1x,a,1x,a,'' --sp --lmo '',a)') &
   &     trim(env%ProgName), trim(fname), trim(env%lmover), trim(pipe)
   call command(trim(jobcall))

!--- cleanup
   call remove('wbo')
   call remove('charges')
   call remove('xtbrestart')
   call remove('xtbscreen.xyz')
   call remove('lmocent.coord')
   call remove('coordprot.0')
end subroutine xtb_lmo

!___________________________________________________________________________________
!
! An xTB-IFF calculation on all available threads
!___________________________________________________________________________________

subroutine xtb_iff(env, file_lmo1, file_lmo2, solu, clus)
   use iso_fortran_env, only: wp => real64
   use iomod
   use crest_data
   use zdata

   implicit none

   type(systemdata)                :: env
   type(zmolecule), intent(in)     :: solu, clus
   character(len=80)               :: pipe
   character(len=512)              :: jobcall
   character(len=*)                :: file_lmo1, file_lmo2

!--- Option setting
   pipe = ' > iff.out 2>/dev/null'

!--- Setting threads
!         if(env%autothreads)then
   call ompautoset(env%threads, 7, env%omp, env%MAXRUN, 1) !set the global OMP/MKL variables for the xtb jobs
!         endif

!--- Jobcall
   if (env%sameRandomNumber) then
      write (jobcall, '(a,1x,a,1x,a,'' -nfrag1 '',i3,'' -ellips '',3f9.3,'' -qcg -test '',a)') &
      &     trim(env%ProgIFF), trim(file_lmo1), trim(file_lmo2), solu%nat, clus%ell_abc, trim(pipe)
   else
      write (jobcall, '(a,1x,a,1x,a,'' -nfrag1 '',i3,'' -ellips '',3f9.3,'' -qcg '',a)') &
      &     trim(env%ProgIFF), trim(file_lmo1), trim(file_lmo2), solu%nat, clus%ell_abc, trim(pipe)
!            &     trim(env%ProgIFF),trim(solvent_file),trim(solute_file),solu%nat,clus%ell_abc,trim(pipe)
   end if
   call command(trim(jobcall))

!--- Cleanup
   call remove('xtbiff_bestsofar.xyz')
   call remove('xtbiff_genstart.xyz')
   call remove('xtbrestart')

end subroutine xtb_iff

!___________________________________________________________________________________
!
! An xTB docking on all available threads
!___________________________________________________________________________________

subroutine xtb_dock(env, fnameA, fnameB, solu, clus)
   use iso_fortran_env, only: wp => real64
   use iomod
   use crest_data
   use zdata

   implicit none

   type(systemdata)                :: env
   type(zmolecule), intent(in)     :: solu, clus
   character(len=*), intent(in)    :: fnameA, fnameB
   character(len=80)               :: pipe
   character(len=512)              :: jobcall
   integer                         :: i, ich

   call remove('xtb_dock.out')
   call remove('xcontrol')

   pipe = ' 2>/dev/null'

!---- writing wall pot in xcontrol
   call write_wall(env, solu%nat, solu%ell_abc, clus%ell_abc, 'xcontrol')

!---- Write directed stuff, if requested
   if (allocated(env%directed_file)) then
      do i=1, size(env%directed_number)
         if &
         & ((i==1 .and. env%directed_number(i) >= clus%nmol) .OR. &
         & (env%directed_number(i) >= clus%nmol .and. env%directed_number(i-1) < clus%nmol)) &
         & then
            open(newunit=ich, file='xcontrol', status='old', position='append', action='write')
                 write(ich,'("$directed")')
                 write(ich,'(a,1x,a)') 'atoms:', trim(env%directed_list(i,1))
                 write(ich,'("$end")')
         end if
      end do
   end if

!--- Setting threads
!   if(env%autothreads)then
   call ompautoset(env%threads, 7, env%omp, env%MAXRUN, 1) !set the global OMP/MKL variables for the xtb jobs
!   endif

   !--- Jobcall docking
   write (jobcall, '(a,1x,''dock'',1x,a,1x,a,1x,a,1x,f4.2,1x,''--nfrag1'',1x,i0,1x,a,1x&
           & ''--input xcontrol > xtb_dock.out'',a)') &
           &     trim(env%ProgName), trim(fnameA), trim(fnameB), trim(env%gfnver),&
           &     env%optlev, solu%nat, trim(env%docking_qcg_flag), trim(pipe)
   call command(trim(jobcall))

! cleanup
   call remove('wbo')
   call remove('charges')
   call remove('xtbrestart')

end subroutine xtb_dock

!___________________________________________________________________________________
!
! An xTB optimization on all available threads
!___________________________________________________________________________________

subroutine opt_cluster(env, solu, clus, fname, without_pot)
   use iso_fortran_env, only: wp => real64
   use iomod
   use crest_data
   use zdata

   implicit none

   type(systemdata)                :: env
   type(zmolecule), intent(in)     :: solu, clus
   character(len=*), intent(in)     :: fname
   logical, optional, intent(in)   :: without_pot
   character(len=80)               :: pipe
   character(len=512)              :: jobcall

   if (env%niceprint) then
      call printprogbar(0.0_wp)
   end if

   call remove('xtb.out')
   pipe = ' 2>/dev/null'

!---- writing wall pot in xcontrol
   if (.not. without_pot) then
      call write_wall(env, solu%nat, solu%ell_abc, clus%ell_abc, 'xcontrol')
   end if

!--- Setting threads
   if (env%autothreads) then
      call ompautoset(env%threads, 7, env%omp, env%MAXRUN, 1) !set the global OMP/MKL variables for the xtb jobs
   end if

!--- Jobcall optimization
   if (.not. without_pot) then
      write (jobcall, '(a,1x,a,1x,a,'' --opt '',f4.2,'' --input xcontrol > xtb_opt.out'',a)') &
      &     trim(env%ProgName), trim(fname), trim(env%gfnver), env%optlev, trim(pipe)
   else
      write (jobcall, '(a,1x,a,1x,a,'' --opt '',f4.2,1x,a,'' > xtb_opt.out'',a)') &
      &     trim(env%ProgName), trim(fname), trim(env%gfnver), env%optlev, trim(env%solv), trim(pipe)
   end if
   call command(trim(jobcall))

! cleanup
   call remove('wbo')
   call remove('charges')
   call remove('xtbrestart')

!--- Jobcall SP for gbsa model
   if (.not. without_pot) then
      write (jobcall, '(a,1x,a,1x,a,'' --sp '',a,'' > xtb_sp.out'',a)') &
      &    trim(env%ProgName), 'xtbopt.coord', trim(env%gfnver), trim(env%solv), trim(pipe)
   end if
   call command(trim(jobcall))

! cleanup
   call remove('wbo')
   call remove('charges')
   call remove('xtbrestart')

end subroutine opt_cluster

!___________________________________________________________________________________
!
! xTB LMO calculation performed in parallel
!___________________________________________________________________________________

subroutine ensemble_lmo(env, fname, self, NTMP, TMPdir, conv)
   use iso_fortran_env, only: wp => real64
   use iomod
   use crest_data
   use zdata

   implicit none
   type(systemdata)                :: env
   type(zmolecule), intent(in)      :: self
   character(len=*), intent(in)     :: fname      !file base name
   character(len=*), intent(in)     :: TMPdir     !directory name
   integer, intent(in)              :: NTMP       !number of structures to be optimized
   integer, intent(in)              :: conv(env%nqcgclust + 1)
   integer                         :: i, k
   integer                         :: vz
   character(len=20)               :: pipe
   character(len=512)              :: thispath, tmppath
   character(len=1024)             :: jobcall
   real(wp)                        :: percent

! setting the threads for correct parallelization
   if (env%autothreads) then
      call ompautoset(env%threads, 7, env%omp, env%MAXRUN, NTMP)
   end if

   pipe = '2>/dev/null'

   !create the system call (it is the same for every optimization)

   write (jobcall, '(a,1x,a,1x,a,'' --sp --lmo --chrg '',f4.1,1x,a,'' >xtb_lmo.out'')') &
   &     trim(env%ProgName), trim(fname), trim(env%lmover), self%chrg, trim(pipe)
   k = 0 !counting the finished jobs
!___________________________________________________________________________________

!$omp parallel &
!$omp shared( vz,jobcall,NTMP,percent,k,TMPdir,conv )
!$omp single
   do i = 1, NTMP
      vz = i
      !$omp task firstprivate( vz ) private( tmppath )
      write (tmppath, '(a,i0)') trim(TMPdir), conv(vz)
      call command('cd '//trim(tmppath)//' && '//trim(jobcall))
      !$omp critical
      k = k + 1
      percent = float(k)/float(NTMP)*100
      !$omp end critical
      !$omp end task
   end do
!$omp taskwait
!$omp end single
!$omp end parallel

!___________________________________________________________________________________

   call getcwd(thispath)
   do i = 1, NTMP
      write (tmppath, '(a,i0)') trim(TMPdir), i
      call chdir(trim(tmppath))
      call remove('xtbrestart')
      call chdir(trim(thispath))
   end do

end subroutine ensemble_lmo

!___________________________________________________________________________________
!
! xTB-IFF calculation performed in parallel
!___________________________________________________________________________________

subroutine ensemble_iff(env, outer_ell_abc, nfrag1, frag1_file, frag2_file, NTMP, TMPdir, conv)
   use iso_fortran_env, only: wp => real64
   use iomod
   use crest_data
   use zdata

   implicit none
   type(systemdata)                :: env

   character(len=*), intent(in)     :: TMPdir     !directory name
   integer, intent(in)              :: NTMP       !number of structures to be optimized
   integer, intent(in)              :: nfrag1     !#atoms of larger fragment
   integer, intent(in)              :: conv(env%nqcgclust + 1)
   real(wp), intent(in)             :: outer_ell_abc(env%nqcgclust, 3)

   integer                         :: i, k
   integer                         :: vz
   character(len=20)               :: pipe
   character(len=512)              :: tmppath
   character(len=1024)             :: jobcall
   character(len=64), intent(in)    :: frag1_file
   character(len=64), intent(in)    :: frag2_file
   character(len=64)               :: frag1
   character(len=64)               :: frag2
   real(wp)                        :: percent

! some options
   pipe = '2>/dev/null'
   frag1 = 'solvent_cluster.lmo'
   frag2 = 'solvent.lmo'

! setting the threads for correct parallelization
   if (env%autothreads) then
      call ompautoset(env%threads, 7, env%omp, env%MAXRUN, NTMP)
   end if

   k = 0 !counting the finished jobs

!___________________________________________________________________________________

!$omp parallel &
!$omp shared( vz,NTMP,percent,k,TMPdir,conv )
!$omp single
   do i = 1, NTMP
      vz = i
      !$omp task firstprivate( vz ) private( tmppath,jobcall )
! create the system call
      write (jobcall, '(a,1x,a,1x,a,'' -nfrag1 '',i3,'' -ellips '',3f9.3,'' -qcg '',a,'' >iff.out'')') &
&         trim(env%ProgIFF), trim(frag1_file), trim(frag2_file), nfrag1, outer_ell_abc(conv(vz), 1:3)*0.9, trim(pipe)
      write (tmppath, '(a,i0)') trim(TMPdir), conv(vz)
      call command('cd '//trim(tmppath)//' && '//trim(jobcall))
      !$omp critical
      k = k + 1
      percent = float(k)/float(NTMP)*100
      !$omp end critical
      !$omp end task
   end do
!$omp taskwait
!$omp end single
!$omp end parallel

!___________________________________________________________________________________

end subroutine ensemble_iff

!___________________________________________________________________________________
!
! xTB docking calculation performed in parallel
!___________________________________________________________________________________

subroutine ensemble_dock(env, outer_ell_abc, nfrag1, frag1_file, frag2_file, n_shell&
        &, n_solvent, NTMP, TMPdir, conv)
   use iso_fortran_env, only: wp => real64
   use iomod
   use crest_data
   use zdata

   implicit none
   type(systemdata)                :: env

   character(len=*), intent(in)     :: TMPdir     !directory name
   integer, intent(in)              :: NTMP       !number of structures to be optimized
   integer, intent(in)              :: nfrag1     !#atoms of larger fragment
   integer, intent(in)              :: conv(env%nqcgclust + 1)
   real(wp), intent(in)             :: outer_ell_abc(env%nqcgclust, 3)
   integer, intent(in)              :: n_shell, n_solvent

   integer                         :: i, k
   integer                         :: vz
   character(len=20)               :: pipe
   character(len=1024)             :: jobcall
   character(len=512)              :: thispath, tmppath
   character(len=*), intent(in)     :: frag1_file
   character(len=*), intent(in)     :: frag2_file
   character(len=64)               :: frag1
   character(len=64)               :: frag2
   real(wp)                        :: percent
   character(len=2)                :: flag
   integer                         :: ich31

! some options
   pipe = '2>/dev/null'
   frag1 = 'solvent_cluster.coord'
   frag2 = 'solvent'
   call getcwd(thispath)

! setting the threads for correct parallelization
   if (env%autothreads) then
      call ompautoset(env%threads, 7, env%omp, env%MAXRUN, NTMP)
   end if

   write (jobcall, '(a,1x,''dock'',1x,a,1x,a,1x,a,1x,f4.2,1x,''--nfrag1'',1x,i0,1x,&
           & ''--input xcontrol --fast > xtb_dock.out '',a)') &
           & trim(env%ProgName), trim(frag1_file), trim(frag2_file),&
           & trim(env%gfnver), env%optlev, nfrag1, trim(pipe)

   flag = '$'
   do i = 1, NTMP
      vz = i
      write (tmppath, '(a,i0)') trim(TMPdir), conv(i)
      call chdir(trim(tmppath))
      open (newunit=ich31, file='xcontrol')
      write (ich31, '(a,"fix")') trim(flag)
      write (ich31, '(3x,"atoms: 1-",i0)') n_shell !Initial number of atoms (starting solvent shell)
      write (ich31, '(a,"wall")') trim(flag)
      write (31, '(3x,"potential=polynomial")')
      write (ich31, '(3x,"ellipsoid:",1x,3(g0,",",1x),i0,"-",i0)') outer_ell_abc(conv(vz), :), &
              & n_shell + 1, n_shell + n_solvent !Initial number of atoms (starting solvent shell)
      close (ich31)
      call chdir(trim(thispath))
   end do

   k = 0 !counting the finished jobs

!___________________________________________________________________________________

!$omp parallel &
!$omp shared( vz,NTMP,percent,k,TMPdir,conv,n_shell,n_solvent,jobcall )
!$omp single
   do i = 1, NTMP
      vz = i
      !$omp task firstprivate( vz ) private( tmppath )
      write (tmppath, '(a,i0)') trim(TMPdir), conv(vz)
      call command('cd '//trim(tmppath)//' && '//trim(jobcall))
      !$omp critical
      k = k + 1
      percent = float(k)/float(NTMP)*100
      !$omp end critical
      !$omp end task
   end do

!$omp taskwait
!$omp end single
!$omp end parallel

!___________________________________________________________________________________
   call chdir(trim(thispath))

end subroutine ensemble_dock

!___________________________________________________________________________________
!
! xTB CFF optimization performed in parallel
!___________________________________________________________________________________

subroutine cff_opt(postopt, env, fname, n12, NTMP, TMPdir, conv, nothing_added)
   use iso_fortran_env, only: wp => real64
   use iomod
   use crest_data
   use strucrd
   implicit none

   type(systemdata)                :: env
   character(len=*), intent(in)     :: fname      !file base name
   character(len=*), intent(in)     :: TMPdir     !directory name
   integer, intent(inout)           :: NTMP       !number of structures to be optimized
   integer, intent(inout)           :: conv(env%nqcgclust + 1)
   logical, intent(in)              :: postopt
   logical, intent(in)              :: nothing_added(env%nqcgclust)
   integer                         :: i, k, n12
   integer                         :: vz
   integer                         :: ich31
   character(len=20)               :: pipe
   character(len=512)              :: thispath, tmppath
   character(len=1024)             :: jobcall
   character(len=2)                :: flag
   real(wp)                        :: percent

! setting the threads for correct parallelization
   if (env%autothreads) then
      call ompautoset(env%threads, 7, env%omp, env%MAXRUN, NTMP)
   end if

   if (postopt) then
      write (*, '(2x,''Starting optimizations + SP  of structures'')')
      write (*, '(2x,i0,'' jobs to do.'')') NTMP
   end if

! postopt eq true => post opt run, which has to be performed in every directory !!!
   if (postopt) then
      k = 0
      NTMP = env%nqcgclust
      do i = 1, env%nqcgclust
         k = k + 1
         conv(k) = i
         conv(env%nqcgclust + 1) = k
      end do
   end if
   pipe = '2>/dev/null'

   call getcwd(thispath)
   do i = 1, NTMP
      write (tmppath, '(a,i0)') trim(TMPdir), conv(i)
      call chdir(trim(tmppath))
      open (newunit=ich31, file='xcontrol')
      if (n12 .ne. 0) then
         flag = '$'
         write (ich31, '(a,"fix")') trim(flag)
         write (ich31, '(3x,"atoms: 1-",i0)') n12 !Initial number of atoms (starting solvent shell)
      end if
      close (ich31)
      if (postopt .and. nothing_added(i)) call remove('xcontrol')
      call chdir(trim(thispath))
   end do

!--- Jobcall
   write (jobcall, '(a,1x,a,1x,a,'' --input xcontrol --opt '',i0,1x,a,'' >xtb.out'')') &
   &    trim(env%ProgName), trim(fname), trim(env%gfnver), nint(env%optlev), trim(pipe)

   if (NTMP .lt. 1) then
      write (*, '(2x,"No structures to be optimized")')
      return
   end if

   k = 0 !counting the finished jobs
   if (postopt) call printprogbar(0.0_wp)
!___________________________________________________________________________________

!$omp parallel &
!$omp shared( vz,jobcall,NTMP,percent,k,TMPdir,conv )
!$omp single
   do i = 1, NTMP
      vz = i
      !$omp task firstprivate( vz ) private( tmppath )
      write (tmppath, '(a,i0)') trim(TMPdir), conv(vz)
      call command('cd '//trim(tmppath)//' && '//trim(jobcall))
      !$omp critical
      k = k + 1
      percent = float(k)/float(NTMP)*100
      if (postopt) then
         call printprogbar(percent)
      end if
      !$omp end critical
      !$omp end task
   end do
!$omp taskwait
!$omp end single
!$omp end parallel

!__________________________________________________________________________________

   do i = 1, NTMP
      write (tmppath, '(a,i0)') trim(TMPdir), conv(i)
      call chdir(trim(tmppath))
      call remove('xtbrestart')
      call chdir(trim(thispath))
   end do

   !create the system call for sp (needed for gbsa model)
   write (jobcall, '(a,1x,a,1x,a,'' --sp '',a,1x,a,'' >xtb_sp.out'')') &
   &    trim(env%ProgName), 'xtbopt.coord', trim(env%gfnver), trim(env%solv), trim(pipe)

   if (NTMP .lt. 1) then
      write (*, '(2x,"Nothing to do")')
      return
   end if

   k = 0 !counting the finished jobs
   if (postopt) call printprogbar(0.0_wp)
!___________________________________________________________________________________

!$omp parallel &
!$omp shared( vz,jobcall,NTMP,percent,k,TMPdir,conv )
!$omp single
   do i = 1, NTMP
      vz = i
      !$omp task firstprivate( vz ) private( tmppath )
      write (tmppath, '(a,i0)') trim(TMPdir), conv(vz)
      call command('cd '//trim(tmppath)//' && '//trim(jobcall))
      !$omp critical
      k = k + 1
      percent = float(k)/float(NTMP)*100
      if (postopt) then
         call printprogbar(percent)
      end if
      !$omp end critical
      !$omp end task
   end do
!$omp taskwait
!$omp end single
!$omp end parallel

!___________________________________________________________________________________

   do i = 1, NTMP
      write (tmppath, '(a,i0)') trim(TMPdir), conv(i)
      call chdir(trim(tmppath))
      call remove('xtbrestart')
      !call remove('xcontrol')
      call chdir(trim(thispath))
   end do

   if (postopt) then
      write (*, *) ''
      write (*, '(2x,"done.")')
   end if

end subroutine cff_opt

!___________________________________________________________________________________
!
! xTB SP performed in parallel
!___________________________________________________________________________________

subroutine ens_sp(env, fname, NTMP, TMPdir)
   use iso_fortran_env, only: wp => real64
   use iomod
   use crest_data
   use strucrd
   implicit none

   type(systemdata)                :: env
   character(len=*), intent(in)     :: fname      !file base name
   character(len=*), intent(in)     :: TMPdir     !directory name
   integer, intent(inout)           :: NTMP       !number of structures to be optimized

   integer                         :: i, k
   integer                         :: vz
   character(len=20)               :: pipe
   character(len=512)              :: thispath, tmppath
   character(len=1024)             :: jobcall
   real(wp)                        :: percent

! setting the threads for correct parallelization
   if (env%autothreads) then
      call ompautoset(env%threads, 7, env%omp, env%MAXRUN, NTMP)
   end if

   write (*, '(2x,''Single point computation with GBSA model'')')
   write (*, '(2x,i0,'' jobs to do.'')') NTMP

   pipe = '2>/dev/null'

   call getcwd(thispath)

   if (NTMP .lt. 1) then
      write (*, '(2x,"No structures to be optimized")')
      return
   end if

!--- Jobcall
   write (jobcall, '(a,1x,a,1x,a,'' --sp '',a,1x,a,'' > xtb_sp.out'')') &
   &    trim(env%ProgName), trim(fname), trim(env%gfnver), trim(env%solv), trim(pipe)

   k = 0 !counting the finished jobs
   call printprogbar(0.0_wp)
!___________________________________________________________________________________

!$omp parallel &
!$omp shared( vz,NTMP,percent,k,TMPdir,jobcall )
!$omp single
   do i = 1, NTMP
      vz = i
      !$omp task firstprivate( vz ) private( tmppath )
      call initsignal()
      !$omp critical
      write (tmppath, '(a,i0)') trim(TMPdir), vz
      !$omp end critical
      call command('cd '//trim(tmppath)//' && '//trim(jobcall))
      !$omp critical
      k = k + 1
      percent = float(k)/float(NTMP)*100
      call printprogbar(percent)
      !$omp end critical
      !$omp end task
   end do
!$omp taskwait
!$omp end single
!$omp end parallel

!__________________________________________________________________________________

   do i = 1, NTMP
      write (tmppath, '(a,i0)') trim(TMPdir), i
      call chdir(trim(tmppath))
      call remove('xtbrestart')
      call chdir(trim(thispath))
   end do
   write (*, *) ''
   write (*, '(2x,"done.")')

end subroutine ens_sp

!___________________________________________________________________________________
!
! xTB Freq compuatation performed in parallel
!___________________________________________________________________________________

subroutine ens_freq(env, fname, NTMP, TMPdir, opt)
   use iso_fortran_env, only: wp => real64
   use iomod
   use crest_data
   use strucrd
   implicit none

   type(systemdata)                :: env
   character(len=*), intent(in)     :: fname      !file base name
   character(len=*), intent(in)     :: TMPdir     !directory name
   integer, intent(inout)           :: NTMP       !number of structures to be optimized

   integer                         :: i, k
   integer                         :: vz
   character(len=20)               :: pipe
   character(len=512)              :: thispath, tmppath
   character(len=1024)             :: jobcall
   real(wp)                        :: percent
   logical                         :: opt

! setting the threads for correct parallelization
   if (env%autothreads) then
      call ompautoset(env%threads, 7, env%omp, env%MAXRUN, NTMP)
   end if

   write (*, '(2x,''Starting reoptimizations + Frequency computation of ensemble'')')
   write (*, '(2x,i0,'' jobs to do.'')') NTMP

   pipe = '2>/dev/null'

   call getcwd(thispath)

   if (NTMP .lt. 1) then
      write (*, '(2x,"No structures to be optimized")')
      return
   end if

   k = 0 !counting the finished jobs
   call printprogbar(0.0_wp)

!--- Jobcall
   if (.not. opt) then
      write (jobcall, '(a,1x,a,1x,a,'' --hess '',a,'' >xtb_freq.out'')') &
       &    trim(env%ProgName), trim(fname), trim(env%gfnver), trim(pipe)
   else
      write (jobcall, '(a,1x,a,1x,a,'' --ohess '',a,'' >xtb_freq.out'')') &
      &    trim(env%ProgName), trim(fname), trim(env%gfnver), trim(pipe)
   end if

!___________________________________________________________________________________

!$omp parallel &
!$omp shared( vz,NTMP,percent,k,TMPdir,jobcall )
!$omp single
   do i = 1, NTMP
      vz = i
      !$omp task firstprivate( vz ) private( tmppath )
      write (tmppath, '(a,i0)') trim(TMPdir), i
      call command('cd '//trim(tmppath)//' && '//trim(jobcall))
      !$omp critical
      k = k + 1
      percent = float(k)/float(NTMP)*100
      call printprogbar(percent)
      !$omp end critical
      !$omp end task
   end do
!$omp taskwait
!$omp end single
!$omp end parallel

!__________________________________________________________________________________

   do i = 1, NTMP
      write (tmppath, '(a,i0)') trim(TMPdir), i
      call chdir(trim(tmppath))
      call remove('xtbrestart')
      call chdir(trim(thispath))
   end do
   write (*, *) ''
   write (*, '(2x,"done.")')

end subroutine ens_freq

!============================================================!
! Read the Energies from a xtbiff output
!============================================================!

subroutine rdxtbiffE(fname, m, n, e)

   implicit none
   integer :: m, n
   character*(*) :: fname
   real*8 :: e(*)

   character*128 :: line
   real*8 :: xx(10)
   integer :: ich, i, j, nn

   open (newunit=ich, file=fname)

   j = 1
10 continue
   read (ich, '(a)', end=999) line
   read (ich, '(a)') line
   call readl(line, xx, nn)
   e(j) = xx(1)
   do i = 1, n
      read (ich, '(a)') line
   end do
   j = j + 1
   goto 10

999 close (ich)
   m = j - 1
end

!============================================================!
! subroutine wr_cluster_cut
! Cuts a cluster file and and writes the parts
!
! On Input: fname          - name of the coord file
!           n1             - number of atoms fragment1
!           n2             - number of atmos fragment2
!           iter           - number of solvent molecules
!           fname_solu_cut - name of outputfile fragment1
!           fname_solv_cut - name of outputfile fragment2
!
!============================================================!

subroutine wr_cluster_cut(fname_cluster, n1, n2, iter, fname_solu_cut, fname_solv_cut)
   use iso_fortran_env, only: wp => real64
   use strucrd

   implicit none
   integer, intent(in)         :: n1, n2, iter
   real(wp)                    :: xyz1(3, n1)
   real(wp)                    :: xyz2(3, n2*iter)
   integer                     :: at1(n1), at2(n2*iter)
   character(len=*), intent(in) :: fname_cluster, fname_solu_cut, fname_solv_cut
   character(len=256)         :: atmp
   character(len=2)           :: a2
   integer                     :: ich, i, k, stat, io, io2

   ich = 142
   open (unit=ich, file=fname_cluster, iostat=stat)
   read (ich, '(a)') atmp
   k = 1
   do i = 1, n1
      read (ich, '(a)', iostat=io) atmp
      if (io < 0) exit
      atmp = adjustl(atmp)
      call coordline(atmp, a2, xyz1(1:3, k), io2)
      at1(k) = e2i(a2)
      k = k + 1
   end do
   k = 1
   do i = 1, n2*iter
      read (ich, '(a)', iostat=io) atmp
      if (io < 0) exit
      atmp = adjustl(atmp)
      call coordline(atmp, a2, xyz2(1:3, k), io2)
      at2(k) = e2i(a2)
      k = k + 1
   end do

   call wrc0(fname_solu_cut, n1, at1, xyz1)
   call wrc0(fname_solv_cut, n2*iter, at2, xyz2)
   close (ich)

end subroutine wr_cluster_cut

subroutine check_iff(neg_E)
   use iso_fortran_env, only: wp => real64
   use crest_data

   implicit none
   integer              :: io, ich
   real(wp)             :: int_E
   character(len=50)    :: tmp
   logical, intent(out) :: neg_E

   logical :: ex
   character(len=*), parameter :: filename = 'xtbscreen.xyz'

   neg_E = .false.
   int_E = 0.0_wp

   inquire (file=filename, exist=ex)
   if (.not. ex) return

   open (newunit=ich, file=filename, status="old", iostat=io)
   if (io == 0) read (ich, '(a)', iostat=io)
   if (io == 0) read (ich, '(a)', iostat=io) tmp
   close (ich)
   if (io /= 0) return

   tmp = adjustl(tmp(11:))
   read (tmp, *, iostat=io) int_E
   neg_E = io == 0 .and. int_E < 0.0_wp

end subroutine check_iff

!----------------------------------------------------------------------------
! write a wall potential in a file used as xtb input

subroutine write_wall(env,n1,rabc1,rabc12,fname)
  use iso_fortran_env, only : wp => real64
  use crest_data

  implicit none

  type(systemdata)     :: env
  integer, intent(in)  :: n1
  real(wp),intent(in)  :: rabc1(3),rabc12(3)
  character (len=8)    :: flag
  character(len=*)    :: fname

  open(unit=31,file=fname)
  flag='$'
  write(31,'(a,"wall")') trim(flag)
  write(31,'(3x,"potential=polynomial")')
  write(31,'(3x,"ellipsoid:",1x,3(g0,",",1x),"all")') rabc12
  write(31,'(3x,"ellipsoid:",1x,3(g0,",",1x),"1-",i0)') rabc1,n1
  if(env%constrain_solu) then
    write(31,'("$fix")')
    write(31,'(3x,"atoms: 1-",i0)') n1
  end if
  call write_cts(31,env%cts)
  call write_cts_biasext(31,env%cts)
  if(env%cts%used) then !Only, if user set constrians is an $end written
     write(31,'(a)') '$end'
  end if

  close(31)

end subroutine write_wall

subroutine check_dock(neg_E)
   use iso_fortran_env, only: wp => real64
   use crest_data
   use iomod, only: minigrep, grepval

   implicit none
   real(wp)             :: int_E
   logical, intent(out) :: neg_E
   logical :: ex
   character(len=*), parameter :: filename = 'xtbscreen.xyz'

   neg_E = .false.
   int_E = 0.0_wp

   call minigrep('xtb_dock.out', '  Lowest Interaction Energy: ********** kcal/mol', ex)
   if (ex) return

   call grepval('xtb_dock.out', 'Lowest Interaction Energy:', ex, int_E)

   if (ex .and. int_E < 0.0_wp) neg_E = .true.

end subroutine check_dock

subroutine write_constraint(env,coord_name,fname)
  use iso_fortran_env, only : wp => real64
  use crest_data
  use iomod

  implicit none

  type(systemdata)     :: env
  character(len=*),intent(in)    :: fname, coord_name

  call copysub(coord_name, 'coord.ref')
  open(unit=31,file=fname)
  call write_cts(31,env%cts)
  call write_cts_biasext(31,env%cts)
  if(env%cts%used) then !Only, if user set constrians is an $end written
     write(31,'(a)') '$end'
  end if

  close(31)

end subroutine write_constraint
