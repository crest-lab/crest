!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 Sebastian Spicher, Philipp Pracht
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
  type(zmolecule) :: solute, solvent, cluster !Information about solvent, solute and cluster
  type(ensemble) :: full_ensemble, solvent_ensemble
 
  integer :: progress
  integer :: io
  character(len=512) :: thispath

!--- Molecule settings
  solute%nmol = 1
  solvent%nmol = 1
  cluster%nmol = 1


  progress = 0
  call getcwd(thispath)


  !>-----------------------------------
  call qcg_head()
  call tim%start(2,'QCG') !start a timer
  !>-----------------------------------

  call write_qcg_setup(env) !Just an outprint of setup
  call read_qcg_input(env,solute,solvent) !Reading mol. data and determining r,V,A

  if(progress .le. env%qcg_runtype) then !grow
    call qcg_setup(env,solute,solvent)
!    call chdir(env%scratchdir)
    cluster = solute
    call qcg_grow(env,solute,solvent,cluster)
    progress=progress+1
    call chdir(thispath)
  end if

  if(progress .le. env%qcg_runtype) then !ensemble
    call print_qcg_ensemble()
    call qcg_ensemble(env,solute,solvent,cluster,full_ensemble,tim)
    progress=progress+1
    call chdir(thispath)
  end if

  if(progress .le. env%qcg_runtype) then !esolv
    if(env%cff)then
      call qcg_cff(env,solute,solvent,cluster,full_ensemble)
    else
!      call qcg_ensemble()
    end if
    call chdir(thispath)
    progress=progress+1
  end if

  if(progress .le. env%qcg_runtype) then !gsolv

      progress=progress+1
  end if

  !<----------------------------------
  call tim%stop(2) !stop a timer

  return
end subroutine crest_solvtool


subroutine write_qcg_setup(env)
  use crest_data
  use iomod
  implicit none

  type(systemdata) :: env

  write(*,*)
  write(*,'(2x,''========================================='')')
  write(*,'(2x,''|   quantum cluster growth: INPUT       |'')')
  write(*,'(2x,''========================================='')')
  write(*,*)
  select case( env%qcg_runtype )
    case( 0 )
      write(*,'(2x,''QCG: Only Cluster Generation'')')
    case( 1 )
      write(*,'(2x,''QCG: Cluster + Ensemble Generation'')')
!      if(opt%performMD) write(*,'(2x,''Ensemble generated via MD Simulation'')')
!      if(opt%performMTD) write(*,'(2x,''Ensemble generated via MetaDynamic'')')
      !if(opt%performCREST) 
      write(*,'(2x,''Ensemble generated via CREST'')')
    case( 2 )
      write(*,'(2x,''QCG: Calculation of delta E_solv'')')
!      if(opt%performMD) write(*,'(2x,''Ensemble generated via MD Simulation'')')
!      if(opt%performMTD) write(*,'(2x,''Ensemble generated via MetaDynamic'')')
      !if(opt%performCREST) 
      write(*,'(2x,''Ensemble generated via CREST'')')
    case( 3 )
      write(*,'(2x,''QCG: Calculation of delta G_solv'')')
!      if(opt%performMD) write(*,'(2x,''Ensemble generated via MD Simulation'')')
!      if(opt%performMTD) write(*,'(2x,''Ensemble generated via MetaDynamic'')')
      !if(opt%performCREST) 
      write(*,'(2x,''Ensemble generated via CREST'')')
    case default
      continue
  end select
  write(*,*)
  write(*,'(2x,''input parameters     '')')
  write(*,'(2x,''solute                 : '',a)') trim(env%solu_file)
  write(*,'(2x,''charge                 : '',i0)') env%chrg
  write(*,'(2x,''uhf                    : '',i0)') env%uhf
  write(*,'(2x,''solvent                : '',a)') trim(env%solv_file)
  if(env%nsolv.ne.0)then
    write(*,'(2x,''# of solvents to add   : '',i0)') env%nsolv
  else if(env%nsolv.eq.0)then
    write(*,'(2x,''# of solvents to add   : until convergence'')')
  end if
!  write(*,'(2x,''# of cluster generated : '',i0)') sys%NRUN
  write(*,'(2x,''# of CPUs used         : '',i0)') env%Threads
  write(*,'(2x,''gbsa model             : '',a)') env%solvent
  write(*,'(2x,''xtb opt level          : '',i0)') nint(env%optlev)
!  write(*,'(2x,''System temperature [K] : '',F5.1)') sys%mdtemp
!  write(*,'(2x,''RRHO scaling factor    : '',F4.2)') sys%freq_scal

end subroutine write_qcg_setup
 


subroutine qcg_setup(env,solu,solv)
  use iso_fortran_env, wp => real64
  use crest_data
  use iomod
  use zdata
  use strucrd
  implicit none

  type(systemdata):: env 
  type(zmolecule) :: solv, solu 

  integer :: io, f, r, v
  character(len=*),parameter :: outfmt = '(1x,1x,a,1x,f14.7,a,1x)'
  logical :: e_there
  character(len=512) :: thispath, tmp_grow, printout

  call getcwd(thispath)

  if(env%scratchdir .eq. '') then !check if scratch was not set
      env%scratchdir='qcg_tmp'
!      env%scratch = .true.
      io = makedir('qcg_tmp')
  end if

!  call scrdir(env)

  call copysub('solute',trim(env%scratchdir))
  call copysub('solvent',env%scratchdir)
  call copysub('.CHRG',env%scratchdir)
  call copysub('.UHF',env%scratchdir)
  call chdir(env%scratchdir)


  f = makedir('solute_properties')
  r = makedir('solvent_properties')
  v = makedir('tmp_grow')


  call chdir('tmp_grow')
  call getcwd(tmp_grow)
  call chdir(thispath)
  call chdir(env%scratchdir)

  call copysub('solute','solute_properties')
  call copysub('.CHRG','solute_properties')
  call copysub('.UHF','solute_properties')
  call copysub('solvent','solvent_properties')
  call copysub('.CHRG','solvent_properties')
  call copysub('.UHF','solvent_properties')

  call copysub('.CHRG',tmp_grow)
  call copysub('.UHF',tmp_grow)

  call remove('solute')
  call remove('solvent')
  call remove('.CHRG')
  call remove('.UHF')

!---- Properties solute
  call chdir('solute_properties')
  env%ref%nat = solu%nat
  env%ref%at = solu%at
  env%ref%xyz = solu%xyz

!---- Geometry preoptimization solute
  if (.not. env%nopreopt) then
    call wrc0('coord',solu%nat,solu%at,solu%xyz) !write coord for xtbopt routine
    call xtbopt(env)
    call rdcoord('coord',solu%nat,solu%at,solu%xyz)
    call remove('coord')
  end if
  call wrc0('solute',solu%nat,solu%at,solu%xyz)

!---- LMO-Computation solute
  write(*,*) 'Generating LMOs for solute'
  call xtb_lmo(env,'solute',solu%chrg)
  call grepval('xtb.out','| TOTAL ENERGY',e_there,solu%energy)
  if (e_there .eq. .false.) then
      write(*,*) 'Total Energy of solute not found'
  else
      write(*,outfmt) 'Total Energy of solute: ', solu%energy, ' Eh'
  end if
  call rename('xtblmoinfo','solute.lmo')
  call copysub('solute.lmo',tmp_grow)

  call chdir(thispath)

!---- Properties solvent
  call chdir(env%scratchdir)
  call chdir('solvent_properties')
  env%ref%nat = solv%nat
  env%ref%at = solv%at
  env%ref%xyz = solv%xyz

!---- Geometry preoptimization solvent
  if(.not. env%nopreopt) then
    call wrc0('coord',solv%nat,solv%at,solv%xyz) !write coord for xtbopt routine
    call xtbopt(env)
    call rdcoord('coord',solv%nat,solv%at,solv%xyz)
    call remove('coord')
  end if
  call wrc0('solvent',solv%nat,solv%at,solv%xyz)

!---- LMO-Computation solvent
  write(*,*) 'Generating LMOs for solvent'
  call xtb_lmo(env,'solvent',solv%chrg)

  call grepval('xtb.out','| TOTAL ENERGY',e_there,solv%energy)
  if (e_there .eq. .false.) then
      write(*,'(1x,a)') 'Total Energy of solvent not found'
  else
      write(*,outfmt) 'Total energy of solvent:', solv%energy,' Eh'
  end if
  call rename('xtblmoinfo','solvent.lmo')
  call copysub('solvent.lmo',tmp_grow)
  call chdir(thispath)

!---- Overwriting solute and solvent in original folder  
  call wrc0('solute',solu%nat,solu%at,solu%xyz)
  call wrc0('solvent',solv%nat,solv%at,solv%xyz)

end subroutine qcg_setup


subroutine read_qcg_input(env,solu,solv)
  use iso_fortran_env, wp => real64
  use crest_data
  use iomod
  use zdata
  use strucrd
  use atmasses
  implicit none

  type(systemdata)               :: env 
  type(zmolecule), intent(inout) :: solu,solv
  logical                        :: pr 
  real(wp),parameter             :: amutokg= 1.66053886E-27
  real(wp),parameter             :: third  = 1.0d0/3.0d0
  integer                        :: i
  real(wp)                       :: r_solu, r_solv

  pr = .true.

!--- Read in nat, at, xyz
  call simpletopo_file('solute',solu,.false.,'')
  allocate(solu%xyz(3,solu%nat))
  call rdcoord('solute',solu%nat,solu%at,solu%xyz)
  call simpletopo_file('solvent',solv,.false.,'')
  allocate(solv%xyz(3,solv%nat))
  call rdcoord('solvent',solv%nat,solv%at,solv%xyz)

!--- CMA-Trafo
  call cma_shifting(env,solu,solv) !During qcg_grow also axistrf done, therefore unnecessary

!--- Setting solute charge and uhf to input
  solu%chrg=env%chrg
  solu%uhf=env%uhf

!--- Getting r, V, A 
  if(pr) then
      write(*,*) 'Solute geometry'
  end if
  call get_sphere(.true.,solu,.true.)  !r,V,A of solute
  if(pr) then
      write(*,*) 'Solvent geometry'
  end if
  call get_sphere(.true.,solv,.true.) !r,V,A of solvent

  r_solu=solu%vtot**third    
  r_solv=solv%vtot**third    
  write(*,'(2x,''radius of solute    : '',f8.2)') r_solu     
  write(*,'(2x,''radius of solvent   : '',f8.2)') r_solv 

!--- Determine masses (for later density computation)
  do i=1,solu%nat
     solu%mass=solu%mass+ams(solu%at(i))
  end do
  do i=1,solv%nat
     solv%mass=solv%mass+ams(solv%at(i))
  enddo
  solu%mass=solu%mass*amutokg
  solv%mass=solv%mass*amutokg

end subroutine read_qcg_input


      subroutine get_sphere(pr,zmol,r_logical)
      use iso_fortran_env, wp => real64
      use zdata

      implicit none
      type(zmolecule), intent(inout) :: zmol
      logical        :: pr
      logical        :: r_logical !Determines wether r is overwritten or not 
  real(wp),parameter :: pi43   = 3.1415926540d0*4.0d0/3.0d0
  real(wp),parameter :: pi     = 3.1415926540d0
  real(wp),parameter :: third  = 1.0d0/3.0d0
  real(wp),parameter :: bohr   = 0.52917726d0

      integer i,j,k
      real*8 rad(zmol%nat),xyz_tmp(3,zmol%nat)
      real(wp),allocatable :: rcov(:)

     allocate(rcov(94))
     rcov = (/ &
      &  2.18230009,  1.73469996,  3.49559999,  3.09820008,  3.21600008, &
      &  2.91030002,  2.62249994,  2.48169994,  2.29959989,  2.13739991, &
      &  3.70819998,  3.48390007,  4.01060009,  3.79169989,  3.50169992, &
      &  3.31069994,  3.10459995,  2.91479993,  4.24109983,  4.10349989, &
      &  3.89030004,  3.76419997,  3.72110009,  3.44140005,  3.54620004, &
      &  3.44210005,  3.43269992,  3.34619999,  3.30080009,  3.23090005, &
      &  3.95790005,  3.86190009,  3.66249990,  3.52679992,  3.36619997, &
      &  3.20959997,  4.61759996,  4.47639990,  4.21960020,  4.05970001, &
      &  3.85960007,  3.75430012,  3.56900001,  3.46230006,  3.39750004, &
      &  3.35249996,  3.33080006,  3.46199989,  4.26230001,  4.18739986, &
      &  4.01499987,  3.89010000,  3.73799992,  3.58890009,  5.05670023, &
      &  5.18139982,  4.62610006,  4.62010002,  4.57019997,  4.52710009, &
      &  4.48960018,  4.45149994,  4.42339993,  4.12430000,  4.24270010, &
      &  4.15409994,  4.27939987,  4.24499989,  4.22079992,  4.19859982, &
      &  4.01300001,  4.24499989,  4.09800005,  3.98550010,  3.89549994, &
      &  3.74900007,  3.44560003,  3.35249996,  3.25640011,  3.35990000, &
      &  4.31269979,  4.27640009,  4.11749983,  4.00540018,  3.86439991, &
      &  3.72160006,  5.07959986,  4.92939997,  4.70429993,  4.42519999, &
      &  4.45940018,  4.39569998,  4.35389996,  4.43410015/)

      do i=1,zmol%nat
         rad(i)=bohr*rcov(zmol%at(i))*1.40 ! scale factor adjusted to rough
                                       ! surfaces are systematically sm
                                       ! for stacked alkane chains and 
         xyz_tmp(1:3,i)=bohr*zmol%xyz(1:3,i)
      enddo

      call arvo(zmol%nat,xyz_tmp,rad,zmol%atot,zmol%vtot)  ! does the job, see arvo.f

      zmol%atot = zmol%atot / bohr**2
      zmol%vtot = zmol%vtot / bohr**3
    if(r_logical) then
      zmol%rtot = zmol%vtot*3.0 / 4.d0 / pi
      zmol%rtot = zmol%rtot**(1.d0/3.d0)
    end if

      if(pr)then
          if(r_logical) then
            write(*,'(2x,''molecular radius (Bohr**1):'',F8.2)') zmol%rtot
          end if
          write(*,'(2x,''molecular area   (Bohr**2):'',F8.2)') zmol%atot
          write(*,'(2x,''molecular volume (Bohr**3):'',F8.2)') zmol%vtot   
      endif

      deallocate(rcov)
      end subroutine get_sphere


subroutine cma_shifting(env,solu,solv)
  use iso_fortran_env, wp => real64
  use crest_data
  use iomod
  use zdata
  use strucrd
  implicit none

  type(systemdata)   :: env 
  type(zmolecule)    :: solu, solv

  integer            :: i

  call cma(solu%nat,solu%at,solu%xyz,solu%cma)
  call cma(solv%nat,solv%at,solv%xyz,solv%cma)

  do i=1,solu%nat
     solu%xyz(1:3,i)=solu%xyz(1:3,i)-solu%cma(1:3)
  end do
  do i=1,solv%nat
     solv%xyz(1:3,i)=solv%xyz(1:3,i)-solv%cma(1:3)
  end do

end subroutine cma_shifting



subroutine qcg_grow(env,solu,solv,clus)
  use iso_fortran_env, wp => real64
  use crest_data
  use iomod
  use zdata
  use strucrd
  implicit none

  type(systemdata)           :: env 
  type(zmolecule)            :: solu, solv, clus
  integer                    :: minE_pos, m
  real(wp)                   :: etmp(500)
  real(wp)                   :: e_each_cycle(env%nsolv)
  character(len=*),parameter :: outfmt = '(1x,1x,a,1x,f14.7,a,1x)'
  logical                    :: e_there
!  real(wp),parameter         :: bohr   = 0.52917726d0
  real(wp)                   :: dens,dum,efix
  real(wp)                   :: e_diff = 0
  integer                    :: iter=1
  real(wp),parameter         :: eh = 627.509541d0
  integer                    :: i,j,io
  real(wp)                   :: E_inter(env%nsolv) 
  real(wp)                   :: shr = 0
  real(wp)                   :: shr_av = 0
  real(wp)                   :: mean = 0
  real(wp)                   :: mean_old = 0
  real(wp)                   :: mean_diff = 0
  real(wp),parameter         :: au = 627.509541d0

  character(len=512)           :: thispath, resultspath


  call print_qcg_grow()
  call getcwd(thispath)
  io = makedir('grow')
  call chdir('grow') !Results directory

!--- Output Files
! open output file: qcg_energy.dat
  open(unit=99,file='qcg_energy.dat')
  write(99,'(i4,2F20.8)') env%nsolv,solu%energy,solv%energy
  open(unit=15,file='qcg_grow.xyz')  ! for molden movie
  open(unit=88,file='qcg_conv.dat')  ! for convergence check
  write(88,'(''   #   Energy       Run. Aver.   Diff / au.'')')

  call getcwd(resultspath)
  call chdir(thispath)

  call get_ellipsoid(env, solu, solv, clus,.true.)
  call pr_grow_energy()

  call chdir(env%scratchdir)
  call chdir('tmp_grow')

  call ellipsout('solute_cavity.coord',clus%nat,clus%at,clus%xyz,solu%ell_abc)
  solv%ell_abc = clus%ell_abc

  clus%chrg = solu%chrg
  clus%uhf = solu%uhf

!--------------------------------------------------------
! Start Loop
!--------------------------------------------------------

do iter=1, env%nsolv

!---- LMO-Computation
  if (iter .gt. 1) then
  call get_ellipsoid(env, solu, solv, clus,.false.)
  call xtb_lmo(env,'xtbopt.coord',clus%chrg)
  call grepval('xtb.out','| TOTAL ENERGY',e_there,clus%energy)
  if (e_there .eq. .false.) then
      write(*,'(1x,a)') 'Total Energy of cluster LMO computation not found'
  end if
  call rename('xtblmoinfo','cluster.lmo')
  end if

  call both_ellipsout('twopot_1.coord',clus%nat,clus%at,clus%xyz,clus%ell_abc,solu%ell_abc)


  if(iter .eq. 1) then
    call xtb_iff(env, 'solute.lmo', 'solvent.lmo', solu, solv) !solu for nat of core pot. solv for outer ellips
  else
    call xtb_iff(env, 'cluster.lmo', 'solvent.lmo', solu, clus)
  end if

!--- Increase cluster size
  clus%nat = clus%nat + solv%nat
  deallocate(clus%at)
  deallocate(clus%xyz)
  allocate(clus%at(clus%nat))
  allocate(clus%xyz(3,clus%nat))
  clus%nmol = clus%nmol + 1

!--- Select xtb-IFF stucture to proceed
  call rdxtbiffE('xtbscreen.xyz',m,clus%nat,etmp) !Get energy of screening
  minE_pos = minloc(etmp(1:m),dim=1)            !Get minimum of those
  call rdxmolselec('xtbscreen.xyz', minE_pos, clus%nat,clus%at,clus%xyz) !Read the struc into clus%xyz
  call remove ('cluster.coord')
  call wrc0 ('cluster.coord',clus%nat,clus%at,clus%xyz)

  call both_ellipsout('twopot_2.coord',clus%nat,clus%at,clus%xyz,clus%ell_abc,solu%ell_abc)

!--- Cluster optimization  
  call opt_cluster(env,solu,clus,'cluster.coord')

  call rdcoord('xtbopt.coord',clus%nat,clus%at,clus%xyz)
  call grepval('xtb.out','| TOTAL ENERGY',e_there,clus%energy)
  e_each_cycle(iter) = clus%energy

  if (e_there .eq. .false.) then
      write(*,'(1x,a)') 'Total Energy of cluster not found.'
!  else
!      write(*,outfmt) 'Total energy of new cluster:', clus%energy,' Eh'
  end if

!--- Calclulate fix energy + diff. energy
!     efix=clus%energy*eh/sqrt(float(clus%nat))
     efix=clus%energy/sqrt(float(clus%nat))
     dum=solu%energy
     if(iter.gt.1) dum=e_each_cycle(iter-1)
     e_diff = e_diff + eh*(e_each_cycle(iter)-solv%energy-dum)


  call ellipsout('cluster_cavity.coord',clus%nat,clus%at,clus%xyz,clus%ell_abc)
  call both_ellipsout('twopot_cavity.coord',clus%nat,clus%at,clus%xyz,clus%ell_abc,solu%ell_abc)

!--- Density calculations         
  call get_sphere(.false.,clus,.false.) !V, A of new cluster
  dens=0.001*(solu%mass+iter*solv%mass)/(1.0d-30*clus%vtot*bohr**3)

!--- Movie file
  write(15,*) clus%nat
  write(15,'('' SCF done '',2F16.8)') eh*(e_each_cycle(iter)-solv%energy-dum) 
  do j=1,clus%nat
      write(15,'(a2,3F24.10)')asym(clus%at(j)),clus%xyz(1:3,j)*bohr 
  enddo

!----------------------------------
! Interaction energy and Output
!----------------------------------

    call get_interaction_E (env,solu,solv,clus,iter,E_inter)
     call analyze_cluster(iter,clus%nat,solu%nat,solv%nat,clus%xyz,clus%at,shr_av,shr)   ! dist of new mol from solute for output

    write(*,'(x,i4,F13.6,1x,f7.2,3x,f7.2,5x,f6.3,3x,f8.3,3x,2f6.1,2x,f8.1,3x,a,x)') &
          & iter,e_each_cycle(iter),eh*(e_each_cycle(iter)-solv%energy-dum),e_diff,dens,efix,shr_av,shr,&
          & clus%vtot,trim(optlevflag(env%optlev))
     write(99,'(i4,F20.10,3x,f8.1)') iter,e_each_cycle(iter),clus%vtot


!--- Calculate mean energy difference between current and last cycle         
!     do i=0,iter-1
!         mean = mean + E_inter(iter-i)
!     end do
     if(iter .gt. 1) mean = mean*(iter-1) !Getting the sum of energies back for next step
     mean = mean + E_inter(iter)
     mean = mean / iter
     mean_diff = mean - mean_old
     mean_old = mean
     write(88,'(i5,1x,3F13.8)') iter,E_inter(iter)*au,mean,mean_diff

!--- Check if converged
     if(abs(mean_diff).lt.1.0d-4.and.iter.gt.10) exit
     if(iter.gt.500.and.env%nsolv.eq.0) then
         write(*,*) 'No convergence could be reached upon adding 500 solvent molecules.'
         write(*,*) 'Continue with further treatment'
     end if
enddo

  
  if (env%nsolv .eq. 0) env%nsolv = iter !if no env%solv was given

!--- output and files
  write(*,*)
  write(*,'(2x,''Growth finished after '',i0,'' solvents added'')') env%nsolv 
  write(*,'(2x,''Results can be found in grow directory'')')
  write(*,'(2x,''Energy list on file <qcg_energy.dat>'')')
  write(*,'(2x,''Interaction energy on file <qcg_conv.dat>'')')
  write(*,'(2x,''Growing process on <qcg_grow.xyz>'')')
  write(*,'(2x,''Final geometry after grow in <cluster.coord>'')')
  write(*,'(2x,''Potentials and geometry written in <cluster_cavity.coord> and <twopot_cavity.coord>'')')

  open(unit=27,file='.nsolv')
  write(27,'(i4)') iter
  close(27)
!  call touch('.qcgok')
  close(99)
  close(88)
  close(15)

  call copysub ('cluster.coord', resultspath)
  call copysub ('twopot_cavity.coord', resultspath)
  call copysub ('cluster_cavity.coord', resultspath)
  call copysub ('solute_cavity.coord', resultspath)
  call rename('xcontrol','wall_potential')
  call copysub ('wall_potential', resultspath)



  call chdir(thispath)
  call chdir(env%scratchdir)
  call rmrf('tmp_grow')

end subroutine qcg_grow


subroutine qcg_ensemble(env,solu,solv,clus,ens,tim)
  use iso_fortran_env, wp => real64
  use crest_data
  use iomod
  use zdata
  use strucrd
  implicit none

  type(systemdata)           :: env 
  type(zmolecule)            :: solu, solv, clus
  type(ensemble)             :: ens,dum
  type(timer)                :: tim
  integer                    :: i
  integer                    :: io,f,r,ich
  character(len=512)         :: thispath,resultspath,tmppath,tmppath2
  character(len=512)         :: scratchdir_tmp
  character(len=512)         :: jobcall
  character(len=256)         :: inpnam, outnam
  character(len=80)          :: fname,pipe
  character(len=20)          :: gfnver_tmp
  real(wp), allocatable      :: e_fix(:), e_clus(:)
  real(wp),parameter         :: eh  = 627.509541d0
  real(wp),parameter         :: au  = 627.509541d0
  real(wp)                   :: S,H,G,e_average,dens
  real(wp)                   :: sasa
  logical                    :: g_found, ex, mdfail
  real(wp)                   :: newtemp,newmass,newmdtime,newmdstep, newhmass
  real(wp)                   :: newmetadlist, newmetadexp, newmetadfac



  call getcwd(thispath)
  f=makedir('ensemble')
  call chdir('ensemble')
  call getcwd(resultspath)
  call chdir(thispath)

    !--- Setting defaults
       env%cts%NCI = .true.  !Activating to have wall pot. written in coord file for xtb
       env%optlev = 2.0d0    !Increaseing percision for ensemble search to minimze scattering

    !--- Setting up potential constraints
       allocate( env%cts%pots(10))
       env%cts%pots=''
       write(env%cts%pots(1),'("$wall")')
       write(env%cts%pots(2),'(2x,"potential=polynomial")')
       write(env%cts%pots(3),'(2x,"ellipsoid:",1x,3(g0,",",1x),"all")') clus%ell_abc  
       write(env%cts%pots(4),'(2x,"ellipsoid:",1x,3(g0,",",1x),"1-",i0)') solu%ell_abc, solu%nat
       call chdir(env%scratchdir)
       scratchdir_tmp = env%scratchdir
!       env%scratchdir=''
       io = makedir('tmp_MTD')
       call copysub('.CHRG',env%scratchdir)
       call copysub('.UHF',env%scratchdir)
       call chdir('tmp_MTD')
       call getcwd(tmppath2)
       call wrc0('crest_input',clus%nat,clus%at,clus%xyz)
       call inputcoords(env,'crest_input')
       call iV2defaultGF(env)         !Setting MTD parameter

       gfnver_tmp = env%gfnver
       write(*,*) 'Method for ensemble search:', env%ensemble_opt
       env%gfnver = env%ensemble_opt  !Setting method for ensemble search


!----------------------------------------------------------------
! Case selection of normal Crest, MD or MTD
!----------------------------------------------------------------

  select case(env%ensemble_method)
    case(0)
       env%iterativeV2 = .true.  !Safeguards more precise ensemble search

       write(*,*) 'Starting ensemble cluster generation by CREST routine'


       call confscript2i(env,tim) !Calling ensemble search


    !For later: Use of entropic routine from Crest
    !                env%properties= 45
    !                env%autozsort = .false.    !turn off zsort (since we are not going to GC anyways)
    !                env%performCross=.false.   !turn off GC
    !                env%entropic = .true.      !indicator for this runtype
    !                env%Maxrestart=1           !turn off MTD iterations (just do one)
    !                env%rotamermds=.false.     !turn off normMDs
    !                env%entropymd=.true.       !special static MTDs
    !                call read_bhess_ref(env,'coord')
    !            env%runver=111             !version  for selection of MTD bias settings
    !            env%doNMR=.true.           !we need equivalencies
    !             call env%addjob(env%properties)

case( 1:2 ) ! Single MD or MTD

    !    write(*,*) 'Performing normal MD'

!---- Setting threads
         if(env%autothreads)then
            call ompautoset(env%threads,7,env%omp,env%MAXRUN,1) !set the global OMP/MKL variables for the xtb jobs
         endif

!--- Setting new defaults for MD/MTD in qcg

         if(env%mdtemp.lt.0.0d0)then
           newtemp=400.00d0 
         else if(env%user_temp .eq. .false.) then
           !env%mdtemp = 149.00
           !newtemp=env%mdtemp
           newtemp=149.0
         else
            newtemp=env%mdtemp
         endif
         
         if(env%user_mdtime .eq. .false.)then
             !env%mdtime = 10.0
             newmdtime = 10.0
         else
             newmdtime = env%mdtime
         end if

         if(env%user_mdstep .eq. .false.) then
            if(env%ensemble_opt .ne. '--gff') then
                !env%mdstep = 4.0d0
                newmdstep = 4.0d0
            else
                newmdstep = 1.5d0
            end if
         else
            newmdstep = env%mdstep
         end if
         if(env%ensemble_opt .ne. '--gff') then
            newhmass = 4.0
         else
            newmass = 5.0
         end if


        if(.not.allocated(env%metadfac))then
           allocate(env%metadfac(1))
           allocate(env%metadexp(1))
           allocate(env%metadlist(1))
        endif
        !env%metadfac(1) = 0.02_wp
        !env%metadexp(1) = 0.1_wp
        !env%metadlist(1) = 10
        newmetadfac = 0.02_wp
        newmetadexp = 0.1_wp
        newmetadlist = 10.0_wp

         fname='coord'
         pipe=' > xtb.out 2>/dev/null'


!--- Writing constrining file xcontrol
!--- Providing xcontrol overwrites constraints in coord file

     open(newunit=ich,file='xcontrol')
     if(env%cts%NCI)then
      do i=1,10
         if(trim(env%cts%pots(i)).ne.'')then
            write(ich,'(a)') trim(env%cts%pots(i))
         endif
      enddo
     endif

     write(ich,'(a)')'$constrain'
     write(ich,'(2x,a,i0)')'atoms: 1-',solu%nat
     write(ich,'(2x,a)')'force constant=0.5'
     write(ich,'(2x,a,a)')'reference=ref.coord'

     write(ich,'(a)') '$md'
     write(ich,'(2x,a,f10.2)') 'hmass=',newhmass
     write(ich,'(2x,a,f10.2)') 'time=',newmdtime
     write(ich,'(2x,a,f10.2)') 'temp=',newtemp
     write(ich,'(2x,a,f10.2)') 'step=',newmdstep
     write(ich,'(2x,a,i0)') 'shake=',env%shake
     write(ich,'(2x,a,i0)') 'dump=',env%mddumpxyz
     write(ich,'(2x,a)') 'dumpxyz=500.0'

    if(env%ensemble_method .EQ. 2)then
       write(ich,'(a)')'$metadyn'
       write(ich,'(2x,a,i0,a,i0)') 'atoms: ',solu%nat+1,'-',clus%nat
       write(ich,'(2x,a,f10.2)') 'save=',newmetadlist
       write(ich,'(2x,a,f10.2)') 'kpush=',newmetadfac
       write(ich,'(2x,a,f10.2)') 'alp=',newmetadexp
    end if
    close(ich)


!--- Writing jobcall

         write(jobcall,'(a,1x,a,1x,a,'' --md --input xcontrol '',a,1x,a,a)') &
         &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),pipe
      !--- slightly different jobcall for QMDFF usage
         if(env%useqmdff)then
           write(jobcall,'(a,1x,a,1x,a,'' --md --input xcontrol --qmdff'',a,1x,a,a)') &
           &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),pipe
         endif

    !--- MD

    if(env%ensemble_method .EQ. 1) then


         call normalMD(fname,env,1,newtemp,newmdtime)

             write(*,*) 'Starting MD with the settings:'
             write(*,'(''     MD time /ps        :'',f8.1)')newmdtime
             write(*,'(''     MD Temperature /K  :'',f8.1)')newtemp
             write(*,'(''     dt /fs             :'',f8.1)')newmdstep
         write(tmppath,'(a,i0)')'NORMMD1'

     r = makedir(tmppath)
     call copysub('xcontrol',tmppath)
     call chdir(tmppath)
     call copy('coord','ref.coord')

     call chdir(tmppath2) 


         call execute_command_line('cd '//trim(tmppath)//' && '//trim(jobcall), exitstat=io)
         inquire(file=trim(tmppath)//'/'//'xtb.trj',exist=ex)
         if(.not.ex .or. io.ne.0)then
         write(*,'(a,i0,a)')'*Warning: MD seemingly failed (no xtb.trj)*'
         else 
         write(*,*)'*MD finished*'
         end if

        if(env%trackorigin)then
        call set_trj_origins('NORMMD','md')
        endif

        call chdir('NORMMD1')
    end if

     !--- MTD

    if(env%ensemble_method .EQ. 2)then

     call MetaMD(env,1,env%mdtime,env%metadfac(1),env%metadexp(1), &
       &               env%metadlist(1))


             write(*,'(a,i4,a)') 'Starting Meta-MD with the settings:'
             write(*,'(''     MTD time /ps       :'',f8.1)')newmdtime
             write(*,'(''     dt /fs             :'',f8.1)')newmdstep
             write(*,'(''     MTD Temperature /K  :'',f8.1)')newtemp
             write(*,'(''     dumpstep(trj) /fs  :'',i8)')env%mddumpxyz
!             write(*,'(''     dumpstep(Vbias)/ps :'',f8.1)')float(env%mddump)/1000d0
             write(*,'(''     Vbias factor k /Eh :'',f8.4)')newmetadfac
             write(*,'(''     Vbias exp α /bohr⁻²:'',f8.2)')newmetadexp

         write(tmppath,'(a,i0)')'METADYN1'
     r = makedir(tmppath)
     call copysub('xcontrol',tmppath)
     call chdir(tmppath)
     call copy('coord','ref.coord')

    call chdir(tmppath2) 

         call execute_command_line('cd '//trim(tmppath)//' && '//trim(jobcall), exitstat=io)
         inquire(file=trim(tmppath)//'/'//'xtb.trj',exist=ex)
         if(.not.ex .or. io.ne.0)then
            write(*,'(a,i0,a)')'*Warning: Meta-MTD seemingly failed (no xtb.trj)*'
         else
            write(*,*)'*MTD finished*'
         end if 


        if(env%trackorigin)then
        call set_trj_origins('METADYN','mtd')
        endif

        call chdir('METADYN1')

    end if

         call rename('xtb.trj','crest_rotamers_0.xyz')
         call copysub('crest_rotamers_0.xyz',tmppath2)
         call dum%open('crest_rotamers_0.xyz')
         call minigrep('xtb.out','MD is unstable, emergency exit',mdfail)
         if(dum%nall .eq. 1) then
             call copysub('xtb.out',resultspath)
             write(*,*) 'ERROR : M(T)D results only in one structure'
             if(mdfail)then
                write(*,*) '        It was unstable'
             else
                write(*,*)'         Probably the M(T)D time was too small'
             end if
             call copysub('xtb.out',resultspath)
             error stop '         Please check the xtb.out file in the ensemble folder'
          end if
         if(mdfail)then
             write(*,*)
             write(*,*) '   WARNING: The MD was unstable.'
             write(*,*) '            Please check the xtb.out file in the ensemble folder.'
             write(*,*)
             call copysub('xtb.out',resultspath)
         end if
         call dum%deallocate
         call chdir(tmppath2)
       call wrc0('coord',clus%nat,clus%at,clus%xyz)
       call inputcoords(env,'coord')


!--- Optimization
          call print_qcg_opt
          if(env%multilevelopt)then
           if(env%optlev >= 1.0d0)then   
           call multilevel_opt(env,6)
           else
           call multilevel_opt(env,99)    
           endif
          endif
  end select


    !---Entropy and enthalpy (Only ensemble, not for use for free solvation energy
!       call grepval('cregen.out.tmp','ensemble average energy (kcal)        :',g_found,e_average)
!       call grepval('cregen.out.tmp','ensemble entropy (cal/mol K)          :',g_found,S)
!       call grepval('cregen.out.tmp','ensemble free energy (kcal/mol)       :   ',g_found,G)

 !       write(*,*)'ensemble average energy (kcal)        :',e_average
  !      write(*,*)'ensemble entropy (cal/mol K)          :',S
   !     write(*,*)'ensemble free energy (kcal/mol)       :',G


    !---Folder management
        call rename('cregen.out.tmp','thermo_data')
        call copysub('thermo_data',resultspath)
        call copysub('crest_conformers.xyz',resultspath)
        call copysub('crest_best.xyz',resultspath)
        call copysub('cre_members.out',resultspath)
        call checkname_xyz(crefile,inpnam,outnam)
        call rename(inpnam,"crest_rotamers.xyz")
        call copysub('crest_rotamers.xyz',resultspath)


    !---Total free enthalpy
       call ens%open('crest_rotamers.xyz') !Read in ensemble
       allocate(e_fix(ens%nall))
       allocate(e_clus(ens%nall))
       do i=1, ens%nall
          call rdxmolselec('crest_rotamers.xyz',i,clus%nat,clus%at,clus%xyz)
          call get_sphere(.false.,clus,.false.)
          dens=0.001*(solu%mass+env%nsolv*solv%mass)/(1.0d-30*clus%vtot*bohr**3)
!          call analyze_cluster(env%nsolv,clus%nat,solu%nat,solv%nat,clus%xyz,clus%at,shr_av,shr) 
!      e_fix(i) = ens%er(i)*eh/sqrt(float(clus%nat))
       end do
       e_clus =  ens%er * au

!--- Getting G,S,H
  write(*,*)
  write(*,'(2x,''------------------------------------------------------------------------'')')
  write(*,'(2x,''------------------------------------------------------------------------'')')
  write(*,'(2x,''Boltz. averaged energy of final cluster:'')')
  call aver(.false.,env,ens%nall,e_clus(1:ens%nall),S,H,G,sasa,.false.)
  write(*,'(7x,''G /Eh     :'',F14.8)') G/au
  write(*,'(7x,''T*S /kcal :'',f8.3)') S

      call chdir(resultspath)
  open(unit=75,file='cluster.dat')
  write(75,*) '#  E_tot [kcal/mol]'
  do i = 1,ens%nall
     write(75,'(i3,2x,F20.10)') i,e_clus(i)
  end do
  close(75)



       deallocate(e_fix)
       deallocate(e_clus)


       
    !---Deleting ensemble tmp
        call chdir(thispath)
        call chdir(env%scratchdir)
        call rmrf('tmp_MTD')

    !---Outprint
        write(*,*)
        write(*,'(2x,''Ensemble generation finished.'')') 
        write(*,'(2x,''Results can be found in ensemble directory'')')
        write(*,'(2x,''Lowest energy conformer on file <crest_best.xyz>'')')
        write(*,'(2x,''List of conformers on file <crest_conformers.xyz>'')')
        write(*,'(2x,''List of rotamers on file <crest_rotamers.xyz>'')')
        write(*,'(2x,''Thermodynamical data on file <thermo_data>'')')

       env%gfnver = gfnver_tmp

  end subroutine qcg_ensemble


subroutine qcg_cff(env,solu,solv,clus,ens)
  use iso_fortran_env, wp => real64
  use crest_data
  use iomod
  use zdata
  use strucrd
  implicit none

  type(systemdata)           :: env 
  type(zmolecule)            :: solu, solv, clus
  type(ensemble)             :: ens
  integer                    :: i
  integer                    :: io,f,r,ich
  integer                    :: nsolv
  integer                    :: ipos, dum
  integer                    :: minE_pos, m
  integer                    :: open_process
!  integer                    :: conv(env%nqcgclust+1)
  character(len=512)         :: thispath,resultspath,tmppath,tmppath2
  character(len=512)         :: scratchdir_tmp
  character(len=512)         :: jobcall
  character(len=256)         :: inpnam, outnam
  character(len=80)          :: fname,pipe
  character(len=20)          :: to
  real(wp), allocatable      :: e_fix(:), e_clus(:), e_empty(:), outer_ell_abc(:,:), inner_ell_abc(:,:)
  real(wp),parameter         :: eh  = 627.509541d0
  real(wp),parameter         :: au  = 627.509541d0
  real(wp)                   :: S,H,G,e_average
  real(wp)                   :: sasa
  real(wp)                   :: etmp(500)
  logical                    :: found, ex,skip
  logical                    :: all_converged = .false.
  logical,allocatable        :: converged(:)
  real(wp)                   :: newtemp,newmass,newmdtime,newmdstep, newhmass
  real(wp)                   :: newmetadlist, newmetadexp, newmetadfac

  allocate(e_empty(env%nqcgclust))
  allocate(converged(env%nqcgclust))
  allocate(outer_ell_abc(env%nqcgclust,3))
  allocate(inner_ell_abc(env%nqcgclust,3))

  dum = 0
  converged = .false.
  open_process = env%nqcgclust

  call pr_eval_solvent()

  if(solu%vtot/solv%vtot.lt.1.0d0)then
    skip=.true.
  end if

!--- Folder management
  call getcwd(thispath)
  call chdir(env%scratchdir)
  call getcwd(tmppath)
  io = makedir('tmp_CFF')
  call chdir('tmp_CFF')
  call getcwd(tmppath2)
  call chdir(tmppath)
  call chdir('solvent_properties')
  call copysub('solvent.lmo',tmppath2)
  call chdir(tmppath2)

!--- SP of each cluster
  call ens%write('ensemble.xyz')
  do i=1, env%nqcgclust
      call rdxmolselec('ensemble.xyz',i,clus%nat,clus%at,clus%xyz)
!-------------------------------------------------------
!      call axistrf(clus%nat,clus%nat,clus%at,clus%xyz) 
!---------------------------------------------------------
      call get_ellipsoid(env, solu, solv, clus,.false.)
      outer_ell_abc(i,1:3) = clus%ell_abc(1:3)*0.9 !Scaling for cff, outer ellipse wall for each cluster
      inner_ell_abc(i,1:3) = solu%ell_abc(1:3)
      write(to,'("TMPCFF",i0)') i
      io = makedir(trim(to))
      call copysub('solvent.lmo',to)
      call chdir(to)
      call wrc0('cluster.coord',clus%nat,clus%at,clus%xyz)
      call wr_cluster_cut('cluster.coord',solu%nat,solv%nat,env%nsolv,'solute_cut.coord','solvent_shell.coord')
      call xtbsp3(env,'solvent_shell.coord')
      call grepval('xtb.out','| TOTAL ENERGY',ex,e_empty(i))
      call copy('solvent_shell.coord','solvent_cluster.coord')
      if(skip) then
        call copy('colvent_cluster.coord','filled_cluster.coord')
        write(*,'(2x,''solute smaller than solvent, cff skipped'')')
      end if
      call chdir(tmppath2)
  end do

!--- If solvent molecules are added
  if(.not. skip) then
    call pr_qcg_fill()
    write(*,'(2x,''now adding solvents to fill cluster...'')')
    call pr_fill_energy()
    write(*,*)

  !--- Main cycle for addition of solvent molecules
    convergence: do while (all_converged .eq. .false.)
    do i=1,env%nqcgclust
       if(converged(i) .eq. .false.) then
         write(to,'("TMPCFF",i0)') i
         call chdir(to)

!--- LMO computation for solvent cluster---------------------------------------------------
!         call xtb_lmo(env,'solvent_cluster.coord',0) !solvent shells are always uncharged
         call ensemble_lmo(env,'solvent_cluster.coord',solv,open_process,'TMPCFF',conv)
!con fehlt noch
!-------------------------------------------------------------------------------------------
         call rename('xtblmoinfo','solvent_cluster.lmo')
         call remove('xtbrestart')

      !--- read in solvent cluster in folder
!         deallocate(clus%at)
         call clus%deallocate()
!         deallocate(clus%xyz)
         call simpletopo_file('solvent_cluster.coord',clus,.false.,'')

      !--- set ellipsoid potential to the one in folder
         clus%ell_abc(1:3)=outer_ell_abc(i,1:3)

!--- Solvent addition to the cluster----------------------------------------------------------
         call xtb_iff(env, 'solvent_cluster.lmo', 'solvent.lmo', clus, clus)!number of atoms and ell_abc of clus
!----------------------------------------------------------------------------------------------
         call remove('xtbrestart')
         call remove('xcontrol')

      !--- Increase cluster size
         clus%nat = clus%nat + solv%nat
         deallocate(clus%at)
         deallocate(clus%xyz)
         allocate(clus%at(clus%nat))
         allocate(clus%xyz(3,clus%nat))
         clus%nmol = clus%nmol + 1

      !--- Select xtb-IFF stucture to proceed
         call rdxtbiffE('xtbscreen.xyz',m,clus%nat,etmp) !Get energy of screening
         minE_pos = minloc(etmp(1:m),dim=1)            !Get minimum of those
         call rdxmolselec('xtbscreen.xyz', minE_pos, clus%nat,clus%at,clus%xyz) !Read the struc into clus%xyz

      !--- Check if converged
         call fill_take(solv%nat,clus%nat,inner_ell_abc(i,1:3),ipos)
         if (ipos.eq.0) then
            converged(i)=.true.
            open_process = open_process - 1
            write(*,'(2x,''no more solvents can be placed inside caviy of cluster: '',i0)') i
            write(*,'(2x,''previous cluster taken...'')')
         else
            call wrc0('filled_cluster.coord',clus%nat,clus%at,clus%xyz)
         end if

         call chdir(tmppath2)

       else
           cycle
       end if
    end do

  !--- Check if everything is converged
    dum = 0
    do i=1, env%nqcgclust
       if (converged(1) .eq. .true.) then
          dum = dum + 1
       end if
    end do
    if(dum .eq. env%nqcgclust) all_converged = .true.

  end do convergence

  end if

  !Now in every TMPPath one coord file with filled_cluster.coord is present

  !OPTIMIZATION


  deallocate(e_empty)

end subroutine qcg_cff


subroutine get_ellipsoid(env,solu,solv,clus,pr1)
  use iso_fortran_env, wp => real64
  use crest_data
  use iomod
  use zdata
  use strucrd
  implicit none

  type(systemdata)   :: env 
  type(zmolecule)    :: solu, solv, clus
  type(zmolecule)    :: dummy_solu, dummy_solv
  real(wp)           :: eax_solu(3), eax_solv(3)
  real(wp)           :: rabc_solu(3), rabc_solv(3)
  real(wp)           :: aniso, sola
  real(wp)           :: rmax_solu, rmax_solv
  real(wp)           :: boxr, roff, r
  real(wp)           :: ell_solu(3), ell_solv(3)
  character (len=10) :: fname
  logical            :: ex,pr,pr1

  real(wp),parameter :: pi43   = 3.1415926540d0*4.0d0/3.0d0
  real(wp),parameter :: pi     = 3.1415926540d0
  real(wp),parameter :: third  = 1.0d0/3.0d0

  pr = .false. !Outprint deactivated

  fname = 'eaxis.qcg'
  inquire(file=fname,exist=ex)

!  if(.not.ex) then   !If no eaxis file was found

!--- Moving all coords to the origin (transformation)
  call axistrf(solu%nat,solu%nat,solu%at,solu%xyz) 
!  call axistrf(solv%nat,solv%nat,solv%at,solv%xyz)  !Not done in original QCG code
  call axistrf(clus%nat,solu%nat,clus%at,clus%xyz)

!--- Overwrite solute and solvent coord in original file with transformed and optimized ones
  call wrc0('solute',solu%nat,solu%at,solu%xyz)
  call wrc0('solvent',solv%nat,solv%at,solv%xyz)

!--- Getting axis
  call axis2(pr1,solu%nat,solu%at,solu%xyz,eax_solu)
  call axis2(pr1,solv%nat,solv%at,solv%xyz,eax_solv)

!--- Computing anisotropy factor of solute and solvent
  sola=sqrt(1.+(eax_solu(1)-eax_solu(3))/((eax_solu(1)+eax_solu(2)+eax_solu(3))/3.))
  aniso=sqrt(1.+(eax_solv(1)-eax_solv(3))/((eax_solv(1)+eax_solv(2)+eax_solv(3))/3.)) ! =1 for a spherical system

!--- Get maximum intramoleclar distance of solute and solvent
  call getmaxrad(solu%nat,solu%at,solu%xyz,rmax_solu)
  call getmaxrad(solv%nat,solv%at,solv%xyz,rmax_solv)

!--- Getting V and A of dummies
  dummy_solu = solu
  dummy_solv = solv !Why is dummy_solv%vtot different to solv%vtot
  call get_sphere(.false.,dummy_solu,.false.)
  call get_sphere(.false.,dummy_solv,.false.)

!--- Write everything into a file
!  open(unit=87,file=fname)
!  open(unit=87,file=fname)
!  write(87,'(3F20.14)') eax_solu
!  write(87,'(4F20.14)') dummy_solu%vtot,dummy_solv%vtot,aniso,sola,rmax_solu,rmax_solv
!  close(87)
!  else
!    if(pr) write(*,'(2x,''ellipsoid axis read from file: '',a)') trim(fname)
!        open(unit=87,file=fname)
!        read(87,*) eax_solu(1:3)
!        read(87,*) dummy_solu%vtot,dummy_solv%vtot,aniso,sola,rmax_solu,rmax_solv
!        close(87)
!      end if
!  end if

!--- Computation of outer Wall        
      roff = sola * dummy_solu%vtot / 1000
!      write(*,*) 'roff:',roff
      boxr=((0.5*aniso*clus%nmol*dummy_solv%vtot+dummy_solu%vtot)/pi43)**third+roff+rmax_solv*0.5 !0.5 both
!      write(*,*) 'boxr', boxr
      r=(boxr**3/(eax_solu(1)*eax_solu(2)*eax_solu(3)))**third       ! volume of ellipsoid = volume of sphere
!      write(*,*) 'r:', r
      rabc_solv= eax_solu * r                              ! outer solvent wall
  
!--- Computation of inner wall
      roff= sola * dummy_solu%vtot / 1000 
      boxr=((sola*dummy_solu%vtot)/pi43)**third+roff+rmax_solu*0.1 !0.1 before
      r=(boxr**3/(eax_solu(1)*eax_solu(2)*eax_solu(3)))**third       ! volume of ellipsoid = volume of sphere
      rabc_solu = eax_solu * r
      dummy_solu%ell_abc(1) = eax_solu(1)**2/sum((eax_solu(1:3))**2)
      dummy_solu%ell_abc(2) = eax_solu(2)**2/sum((eax_solu(1:3))**2)
      dummy_solu%ell_abc(3) = eax_solu(3)**2/sum((eax_solu(1:3))**2)
      rabc_solu = dummy_solu%ell_abc * r


      if(pr1) then
        write(*,'(2x,''solvent anisotropy  :'',4f10.3)') aniso
        write(*,'(2x,''solute anisotropy   :'',4f10.3)') sola
        write(*,'(2x,''roff inner wall     :'',4f10.3)') roff
        write(*,'(2x,''solute max dist     :'',4f10.3)') rmax_solu
        write(*,'(2x,''solvent max dist    :'',4f10.3)') rmax_solv
        write(*,'(2x,''inner unit axis     :'',3f10.3)') dummy_solu%ell_abc(1:3)
        write(*,'(2x,''inner ellipsoid/Bohr:'',3f10.3)') rabc_solu(1:3)
        write(*,'(2x,''outer ellipsoid/Bohr:'',3f10.3)') rabc_solv(1:3)
        write(*,*)
      endif


      solu%aniso = sola
      solv%aniso = aniso
      solu%ell_abc = rabc_solu
!      solv%ell_abc = rabc_solv
      clus%ell_abc = rabc_solv

  end subroutine get_ellipsoid

!___________________________________________________________________________________
!
! determine max interatomic distance + 2* vdW          
!___________________________________________________________________________________

subroutine getmaxrad(n,at,xyz,r)
  use iso_fortran_env, wp => real64
  implicit none
  real(wp) :: xyz(3,n),r
  integer :: n,at(n)

  real(wp) :: rx,ry,rz,rr
  integer :: i,j
  real(wp),allocatable :: rcov(:)

  allocate(rcov(94))
     rcov = (/ &
      &  2.18230009,  1.73469996,  3.49559999,  3.09820008,  3.21600008, &
      &  2.91030002,  2.62249994,  2.48169994,  2.29959989,  2.13739991, &
      &  3.70819998,  3.48390007,  4.01060009,  3.79169989,  3.50169992, &
      &  3.31069994,  3.10459995,  2.91479993,  4.24109983,  4.10349989, &
      &  3.89030004,  3.76419997,  3.72110009,  3.44140005,  3.54620004, &
      &  3.44210005,  3.43269992,  3.34619999,  3.30080009,  3.23090005, &
      &  3.95790005,  3.86190009,  3.66249990,  3.52679992,  3.36619997, &
      &  3.20959997,  4.61759996,  4.47639990,  4.21960020,  4.05970001, &
      &  3.85960007,  3.75430012,  3.56900001,  3.46230006,  3.39750004, &
      &  3.35249996,  3.33080006,  3.46199989,  4.26230001,  4.18739986, &
      &  4.01499987,  3.89010000,  3.73799992,  3.58890009,  5.05670023, &
      &  5.18139982,  4.62610006,  4.62010002,  4.57019997,  4.52710009, &
      &  4.48960018,  4.45149994,  4.42339993,  4.12430000,  4.24270010, &
      &  4.15409994,  4.27939987,  4.24499989,  4.22079992,  4.19859982, &
      &  4.01300001,  4.24499989,  4.09800005,  3.98550010,  3.89549994, &
      &  3.74900007,  3.44560003,  3.35249996,  3.25640011,  3.35990000, &
      &  4.31269979,  4.27640009,  4.11749983,  4.00540018,  3.86439991, &
      &  3.72160006,  5.07959986,  4.92939997,  4.70429993,  4.42519999, &
      &  4.45940018,  4.39569998,  4.35389996,  4.43410015/)

  r=0
  do i=1,n-1
     do j=i+1,n
     rx=xyz(1,i)-xyz(1,j)
     ry=xyz(2,i)-xyz(2,j)
     rz=xyz(3,i)-xyz(3,j)                                                                    
     rr=sqrt(rx**2+ry**2+rz**2)+rcov(at(i))+rcov(at(j))
     if(rr.gt.r) r=rr
     enddo
  enddo

  deallocate(rcov)

end subroutine getmaxrad


subroutine ellipsout(fname,n,at,xyz,r1)
  use iso_fortran_env, only : wp => real64
  implicit none

  integer            :: i,j
  integer            :: n,at(n)
  real(wp)           :: dum(3)
  real(wp)           :: rx,ry,rz 
  real(wp)           :: xyz(3,n),r1(3)
  real               :: x,y,z,f,rr
  character(len=*)   :: fname
  character(len=2)   :: asym
  

  open(unit=11,file=fname)
  write(11,'(a)')'$coord'
  do i=1,n
     write(11,'(3F24.14,6x,a2)') xyz(1,i),xyz(2,i),xyz(3,i),asym(at(i))
  enddo
  do i=1,500
     call random_number(x)
     call random_number(f)
     if(f.gt.0.5) x=-x
     call random_number(y)
     call random_number(f)
     if(f.gt.0.5) y=-y
     call random_number(z)
     call random_number(f)
     if(f.gt.0.5) z=-z
     rr=sqrt(x*x+y*y+z*z)
     x=x*r1(1)/rr
     y=y*r1(2)/rr
     z=z*r1(3)/rr
     write(11,'(3F24.14,6x,a2)') x,y,z,asym(2)
  enddo
  write(11,'(a)')'$end'
  close(11)

end subroutine ellipsout

subroutine both_ellipsout(fname,n,at,xyz,r1,r2)
  use iso_fortran_env, only : wp => real64
  implicit none

  integer            :: i,j
  integer            :: n,at(n)
  real(wp)           :: dum(3)
  real(wp)           :: rx,ry,rz 
  real(wp)           :: xyz(3,n),r1(3)
  real(wp), optional :: r2(3)
  real               :: x,y,z,f,rr
  character(len=*)   :: fname
  character(len=2)   :: asym
  

  open(unit=11,file=fname)
  write(11,'(a)')'$coord'
  do i=1,n
     write(11,'(3F24.14,6x,a2)') xyz(1,i),xyz(2,i),xyz(3,i),asym(at(i))
  enddo
  do i=1,500
     call random_number(x)
     call random_number(f)
     if(f.gt.0.5) x=-x
     call random_number(y)
     call random_number(f)
     if(f.gt.0.5) y=-y
     call random_number(z)
     call random_number(f)
     if(f.gt.0.5) z=-z
     rr=sqrt(x*x+y*y+z*z)
     x=x*r1(1)/rr
     y=y*r1(2)/rr
     z=z*r1(3)/rr
     write(11,'(3F24.14,6x,a2)') x,y,z,asym(2)
  enddo
  if (present(r2)) then
     do i=1,100
     call random_number(x)
     call random_number(f)
     if(f.gt.0.5) x=-x
     call random_number(y)
     call random_number(f)
     if(f.gt.0.5) y=-y
     call random_number(z)
     call random_number(f)
     if(f.gt.0.5) z=-z
     rr=sqrt(x*x+y*y+z*z)
     x=x*r2(1)/rr
     y=y*r2(2)/rr
     z=z*r2(3)/rr
     write(11,'(3F24.14,6x,a2)')x,y,z,asym(5)
     enddo
  end if
  write(11,'(a)')'$end'
  close(11)

end subroutine both_ellipsout

subroutine get_interaction_E(env,solu,solv,clus,iter,E_inter)
  use iso_fortran_env, wp => real64
  use crest_data
  use iomod
  use zdata
  use strucrd
  implicit none

  type(systemdata)            :: env 
  type(zmolecule), intent(in) :: solu, solv, clus
  type(zmolecule)             :: dummy
  character (len=10)          :: fname
  real(wp)                    :: e_cluster, e_solute, e_solvent
  real(wp)                    :: E_inter(env%nsolv)           ! interaction energy
  integer                     :: iter
  logical                     :: e_there

  call remove('cluster.coord')

!--- Prepare input coordinate files
  call wrc0('cluster.coord',clus%nat,clus%at,clus%xyz)
  call wr_cluster_cut('cluster.coord',solu%nat,solv%nat,iter,'solute_cut.coord','solvent_cut.coord')

!--- Perform single point calculations and recieve energies
  call xtbsp3(env,'cluster.coord')
  call grepval('xtb.out','| TOTAL ENERGY',e_there,e_cluster)
  if(.not.e_there) write(*,*) 'Cluster energy not found'
  call xtbsp3(env,'solute_cut.coord')
  call grepval('xtb.out','| TOTAL ENERGY',e_there,e_solute)
  if(.not.e_there) write(*,*) 'Solute energy not found'
  call xtbsp3(env,'solvent_cut.coord')
  call grepval('xtb.out','| TOTAL ENERGY',e_there,e_solvent)
  if(.not.e_there) write(*,*) 'Solvent energy not found'

  E_inter(iter) = e_cluster - e_solute - e_solvent  
 
end subroutine get_interaction_E


subroutine analyze_cluster(nsolv,n,nS,nM,xyz,at,av,last)
  use iso_fortran_env, only : wp => real64
  implicit none
  real(wp) xyz(3,n)
  real(wp) av,last
  integer n,nS,nM,nsolv,at(n)
  real(wp) xyzM(3,nM)
  integer atm(nM)
  real(wp) xyzS(3,nS)
  integer atS(nS)
  real(wp) x1(3),x2(3),r
  integer i,is,ie
  
  if(nsolv.eq.1) return
  xyzS(1:3,1:nS)=xyz(1:3,1:nS)
   atS(1:nS)=at(1:nS)
  call cma(nS,atS,xyzS,x1)   
  
  av=0
  do i=1,nsolv
     is=nS+(i-1)*nM+1
     ie=is+nM-1
     xyzM(1:3,1:nM)=xyz(1:3,is:ie)
      atM(1:nM)=at(is:ie)
     call cma(nM,atM,xyzM,x2)   
     r=sqrt((x1(1)-x2(1))**2+(x1(2)-x2(2))**2+(x1(3)-x2(3))**2)
     if(i.lt.nsolv)then
        av=av+r        
     else
        last=r
     endif
  enddo
  av=av/float(nsolv-1)
end subroutine analyze_cluster

subroutine aver(pr,env,runs,e_tot,S,H,G,sasa,a_present,a_tot)
      use iso_fortran_env, only : wp => real64
      use crest_data

      implicit none
!---- Dummy      
      type(systemdata),intent(in)     :: env
      integer, intent(in)             :: runs
      real(wp), intent(inout)         :: e_tot
      real(wp), intent(in), optional  :: a_tot
      real(wp), intent(out)           :: S
      real(wp), intent(out)           :: H
      real(wp), intent(out)           :: G
      real(wp), intent(out)           :: sasa
!---- Stack      
      logical, intent(in)             :: pr,a_present
      integer                         :: i,j,jmin
      real(wp)                        :: A
      real(wp)                        :: e0
      real(wp), allocatable           :: de(:)
      real(wp), allocatable           :: p(:)
      real(wp), allocatable           :: deg(:)
      real(wp)                        :: pmax
      real(wp)                        :: eav
      real(wp)                        :: mv     ! mean value
      real(wp)                        :: md     ! mean deviation
      real(wp)                        :: sd     ! standard deviation
      real(wp)                        :: area
      real(wp)                        :: beta
      real(wp)                        :: temp
      dimension e_tot(runs)
      dimension a_tot(runs)
  
      temp=env%mdtemp

      allocate(de(runs), source = 0.0d0)
      allocate(p(runs), source = 0.0d0)
!      allocate(deg(runs))
      deg=1 !degeneracies
  
      beta = 1./(temp*8.314510/4.184/1000.+1.d-14)
!      call sort_min(runs,1,1,e_tot)
!      e0=minval(e_tot(1:runs))
      e0=e_tot(1)
      de(1:runs)=(e_tot(1:runs)-e0)
      call qcg_boltz(env,runs,de,p)

!      call entropy_boltz(runs,env%mdtemp,de,deg,p)

      if(pr)then
        write(*,'(7x,''De /kcal  :'',200f9.4)') de(1:runs)
        write(*,'(7x,''p         :'',200f9.4)')  p(1:runs)
      end if
      A=0
      eav=0
      pmax=0
      area=0
      do j = 1,runs
         A=A+p(j)*log(p(j)+1.d-12)
         eav=eav+p(j)*e_tot(j)
         if(p(j).gt.pmax) then
           pmax=p(j)
           jmin=j
         end if
!         if(present (a_tot) ) area=area+p(j)*a_tot(j)
         if(a_present) area=area+p(j)*a_tot(j)
      end do
      write(*,*)'end it'
      sasa= area
      S=(1./beta)*A
      H=eav
      G=eav+S
  
      deallocate(de,p)
  
    end subroutine aver

    subroutine qcg_boltz(env,n,e,p)
      use iso_fortran_env , only : wp => real64
      use  crest_data
      implicit none      
      type(systemdata),intent(in)    :: env
      integer,intent(in)              :: n 
      real(wp),intent(in)             :: e(*)
      real(wp),intent(out)            :: p(*)
      integer                         :: i
      real(wp)                        :: temp
      real(wp)                        :: f,hsum,esum

      temp=env%mdtemp
      f=8.314*temp/4.184d+3
      esum=0
      do i=1,n
         esum=esum+exp(-e(i)/f)
      end do
      hsum=0
      do i=1,n
         p(i)=exp(-e(i)/f)/esum
      end do
    end subroutine qcg_boltz

subroutine fill_take(n2,n12,rabc,ipos)
  use iso_fortran_env, only : wp => real64
  use crest_data
  use strucrd

  implicit none
  
  integer, intent(in)   :: n2,n12
  real(wp),intent(in)   :: rabc(3)
  integer, intent(out)  :: ipos
  integer               :: i,j,m,n21
  integer               :: at2(n2),at12(n12)
  integer               :: counter
  real(wp)              :: xyz2(3,n2),xyz12(3,n12)
  real(wp)              :: etmp(100)
  real(wp)              :: eabc
  real(wp)              :: cma2(3)
  real(wp),allocatable  :: dist(:)
  
  dist=0
  eabc=0
  counter=0
  n21=n12-n2+1
  call rdxtbiffE('xtbscreen.xyz',m,n12,etmp)
  
  allocate(dist(m),source=0.0d0)
  
  do i=1,m
     call rdxmolselec('xtbscreen.xyz',i,n12,at12,xyz12)
     at2(1:n2)=at12(n21:n12)
     xyz2(1:3,1:n2)=xyz12(1:3,n21:n12)
     call cma(n2,at2,xyz2,cma2)
     call calc_dist(cma2,rabc,dist(i),eabc)
     if(eabc.gt.1.0d0)then
       dist(i) = 1.0d42
       counter = counter + 1
     end if
  end do
  
  ipos=minloc(dist(1:m),dim=1)
  
  if(counter.eq.m) ipos = 0
  
  deallocate(dist)
end subroutine fill_take

subroutine calc_dist(xyz,rabc,dist,eabc)
  use iso_fortran_env, only : wp => real64
  implicit none
  
  real(wp), intent(in)    :: xyz(3)
  real(wp), intent(in)    :: rabc(3)
  real(wp), intent(out)   :: dist
  real(wp), intent(out)   :: eabc
  real(wp)                :: center(3),rc(3)
  
  center=0.d0
  rc=(xyz(1:3)-center)
  dist = norm2(rc)
  eabc = sum((xyz(1:3)**2)/(rabc(1:3)**2))
end subroutine calc_dist


