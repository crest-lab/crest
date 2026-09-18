!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022-2023  Philipp Pracht
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

!> A collection of routines to set up OMP-parallel runs of MDs and optimizations.

!========================================================================================!
!========================================================================================!
!> Interfaces to handle optional arguments
!========================================================================================!
!========================================================================================!
module parallel_interface
!*******************************************************
!* module to load an interface to the parallel routines
!* mandatory to handle any optional input arguments
!*******************************************************
  implicit none
  interface
    subroutine crest_sploop(env,nall,structures,eread,silent)
      use crest_parameters,only:wp,stdout,sep
      use crest_calculator
      use omp_lib
      use crest_data
      use strucrd
      implicit none
      type(systemdata),intent(inout) :: env
      integer,intent(in) :: nall
      type(coord),intent(inout) :: structures(nall)
      real(wp),intent(inout),optional :: eread(nall)
      logical,intent(in),optional :: silent
    end subroutine crest_sploop
  end interface

  !> crest_oloop is generic: the coord-list form is canonical (PBC-capable),
  !> the flat (nat,nall,at,xyz) form is a legacy adapter onto it.
  interface crest_oloop
    subroutine crest_oloop_struc(env,nall,structures,dump,customcalc,eread,silent)
      use crest_parameters,only:wp,stdout,sep
      use crest_calculator
      use omp_lib
      use crest_data
      use strucrd
      use optimize_module
      use iomod,only:makedir,directory_exist,remove
      implicit none
      type(systemdata),target,intent(inout) :: env
      integer,intent(in) :: nall
      type(coord),intent(inout) :: structures(nall)
      logical,intent(in) :: dump
      type(calcdata),intent(in),target,optional :: customcalc
      real(wp),intent(inout),optional :: eread(nall)
      logical,intent(in),optional :: silent
    end subroutine crest_oloop_struc

    subroutine crest_oloop_xyz(env,nat,nall,at,xyz,eread,dump,customcalc)
      use crest_parameters,only:wp
      use crest_calculator
      use crest_data
      use strucrd
      implicit none
      type(systemdata),target,intent(inout) :: env
      integer,intent(in) :: nat,nall
      integer,intent(in) :: at(nat)
      real(wp),intent(inout) :: xyz(3,nat,nall)
      real(wp),intent(inout) :: eread(nall)
      logical,intent(in) :: dump
      type(calcdata),intent(in),target,optional :: customcalc
    end subroutine crest_oloop_xyz
  end interface crest_oloop

  interface
    subroutine crest_hessloop(env,nat,nall,at,xyz,eread,gt_out,stot_out)
      use crest_parameters,only:wp,stdout,sep
      use crest_calculator
      use omp_lib
      use crest_data
      use strucrd
      use thermochem_module
      use iomod,only:makedir,directory_exist,remove
      implicit none
      type(systemdata),intent(inout) :: env
      real(wp),intent(inout) :: xyz(3,nat,nall)
      integer,intent(in)  :: at(nat)
      real(wp),intent(inout) :: eread(nall)
      integer,intent(in) :: nat,nall
      real(wp),optional,intent(out) :: gt_out(:,:)    !> (nall, nt_full)
      real(wp),optional,intent(out) :: stot_out(:,:)  !> (nall, nt_full)
    end subroutine crest_hessloop
  end interface

end module parallel_interface

!========================================================================================!
!========================================================================================!
!> Routines for concurrent singlepoint evaluations
!========================================================================================!
!========================================================================================!
subroutine crest_sploop(env,nall,structures,eread,silent)
!****************************************************************
!* subroutine crest_sploop
!* Concurrent singlepoint evaluations for a list of structures
!* passed as an array of coord objects. Each coord carries its
!* own %lat/%chrg/%uhf, so periodic and heterogeneous systems
!* are handled. xyz must be in Bohr.
!* Energies are stored in structures(i)%energy; the optional
!* eread array, if present, additionally receives them.
!* silent - suppress the progress bar (optional, default .false.)
!****************************************************************
  use crest_parameters,only:wp,stdout,sep
  use crest_calculator
  use omp_lib
  use crest_data
  use strucrd
  use iomod,only:makedir,directory_exist,remove
  use term_ui,only:progress_init,progress_update,progress_finish
  implicit none
  type(systemdata),intent(inout) :: env
  integer,intent(in) :: nall
  type(coord),intent(inout) :: structures(nall)
  real(wp),intent(inout),optional :: eread(nall)
  logical,intent(in),optional :: silent

  type(coord),allocatable :: mols(:)
  integer :: i,j,io,c,k,z,zcopy
  logical :: ex,quiet
  type(calcdata),allocatable :: calculations(:)
  real(wp) :: energy
  real(wp),allocatable :: grad(:,:)
  integer :: thread_id,vz,job
  character(len=80) :: atmp
  real(wp) :: percent,runtime
  type(timer) :: profiler
  integer :: T,Tn
  logical :: nested

!>--- check if we have any calculation settings allocated
  if (env%calc%ncalculations < 1) then
    write (stdout,*) 'no calculations allocated'
    return
  end if

!>--- silent mode? (suppress progress bar + summary printout)
  quiet = .false.
  if (present(silent)) quiet = silent

!>--- prepare calculation objects for parallelization (one per thread)
  call new_ompautoset(env,'auto_nested',nall,T,Tn)
  nested = env%omp_allow_nested
  if (.not.quiet) call ompautoset_summary(env,'singlepoints',T,Tn)

!>--- prepare objects for parallelization
  allocate (calculations(T))
  allocate (mols(T))
  do i = 1,T
    call calculations(i)%copy(env%calc)
    do j = 1,env%calc%ncalculations
      !>--- directories and io preparation
      ex = directory_exist(env%calc%calcs(j)%calcspace)
      if (.not.ex) then
        io = makedir(trim(env%calc%calcs(j)%calcspace))
      end if
      write (atmp,'(a,"_",i0)') sep,i
      calculations(i)%calcs(j)%calcspace = env%calc%calcs(j)%calcspace//trim(atmp)
      if (allocated(calculations(i)%calcs(j)%calcfile)) deallocate (calculations(i)%calcs(j)%calcfile)
      if (allocated(calculations(i)%calcs(j)%systemcall)) deallocate (calculations(i)%calcs(j)%systemcall)
      call calculations(i)%calcs(j)%printid(i,j)
    end do
    calculations(i)%pr_energies = .false.
  end do

!>--- timer initialization
  call profiler%init(1)
  call profiler%start(1)

!>--- initialize progress bar
  if (.not.quiet) then
    call progress_init(env%ps,nall,width=50,prefix=" ↳ ", &
      &                suffix="",show_time=.true.,show_eta=.false.)
    call progress_update(env%ps,0,nall)
  end if

!>--- shared variables
  c = 0  !> counter of successful evaluations
  k = 0  !> counter of total evaluations (fail+success)
  z = 0  !> counter to process structures in order (1...nall)
!>--- pre-start server-based calculators before forking OMP threads
  call preinit_mlip_parallel(calculations,T)
!>--- loop over the structures
  !$omp parallel &
  !$omp shared(env,calculations,nall,structures,c,k,z,mols,nested,Tn)
  !$omp single
  do i = 1,nall

    call initsignal()
    vz = i
    !$omp task firstprivate( vz ) private(i,j,job,energy,grad,io,thread_id,zcopy)
    call initsignal()

    !>--- OpenMP nested region threads
    if (nested) call ompmklset(Tn)

    thread_id = OMP_GET_THREAD_NUM()
    job = thread_id+1
    !>--- deep-copy this structure into the thread-local working mol
    !$omp critical
    z = z+1
    zcopy = z
    call mols(job)%copy(structures(zcopy))
    !$omp end critical

    allocate (grad(3,mols(job)%nat),source=0.0_wp)

    !>-- energy+gradient call
    call engrad(mols(job),calculations(job),energy,grad,io)

    !$omp critical
    if (io == 0) then
      c = c+1
      structures(zcopy)%energy = energy
    else
      structures(zcopy)%energy = 0.0_wp
    end if
    k = k+1
    !>--- print progress
    if (.not.quiet) call progress_update(env%ps,k,nall)
    !$omp end critical

    deallocate (grad)
    !$omp end task
  end do
  !$omp taskwait
  !$omp end single
  !$omp end parallel

!>--- energies are stored in the structures; optionally also return them
  if (present(eread)) then
    do i = 1,nall
      eread(i) = structures(i)%energy
    end do
  end if

!>--- finalize progress printout
  if (.not.quiet) call progress_finish(env%ps)

!>--- stop timer
  call profiler%stop(1)

!>--- prepare some summary printout
  if (.not.quiet) then
    percent = float(c)/float(nall)*100.0_wp
    write (atmp,'(f5.1,a)') percent,'% success)'
    write (stdout,'(">",1x,i0,a,i0,a,a)') c,' of ',nall,' structures successfully evaluated (', &
    &     trim(adjustl(atmp))
    write (atmp,'(">",1x,a,i0,a)') 'Total runtime for ',nall,' singlepoint calculations:'
    call profiler%write_timing(stdout,1,trim(atmp),.true.)
    runtime = profiler%get(1)
    write (atmp,'(f16.3,a)') runtime/real(nall,wp),' sec'
    write (stdout,'(a,a,a)') '> Corresponding to approximately ',trim(adjustl(atmp)), &
    &                       ' per processed structure'
  end if

  call profiler%clear()
  deallocate (calculations)
  if (allocated(mols)) deallocate (mols)
  return
end subroutine crest_sploop

!========================================================================================!
!========================================================================================!
!> Routines for concurrent singlepoint evaluations
!========================================================================================!
!========================================================================================!
subroutine crest_hessloop(env,nat,nall,at,xyz,eread,gt_out,stot_out)
!***************************************************************
!* subroutine crest_hessloop
!* Concurrent numerical Hessian evaluations for an ensemble.
!* Input eread is overwritten with Gibbs free energies.
!* xyz must be in Bohrs.
!* eread contains only the gt@RT on output!
!* Optional gt_out/stot_out return G and S at all temperatures
!* from env%thermo; requires pre-allocated (nall,nt) arrays.
!*
!* Parallelization is enabled using numhess1 (OpenMP-compatible).
!*
!***************************************************************
  use crest_parameters,only:wp,stdout,sep
  use crest_calculator
  use omp_lib
  use crest_data
  use strucrd
  use optimize_module
  use thermochem_module
  use iomod,only:makedir,directory_exist,remove
  use term_ui,only:progress_init,progress_update,progress_finish
  implicit none
  type(systemdata),intent(inout) :: env
  real(wp),intent(inout) :: xyz(3,nat,nall)
  integer,intent(in)  :: at(nat)
  real(wp),intent(inout) :: eread(nall)
  integer,intent(in) :: nat,nall
  real(wp),optional,intent(out) :: gt_out(:,:)    !> (nall, nt_full)
  real(wp),optional,intent(out) :: stot_out(:,:)  !> (nall, nt_full)

  type(coord),allocatable :: mols(:)
  integer :: i,j,k,l,io,ich,ich2,c,z,job_id,zcopy,nat3
  logical :: pr,wr,ex
  type(calcdata),allocatable :: calculations(:)
  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:),grads(:,:,:)
  real(wp),allocatable :: freqs(:,:),hess(:,:,:)
  integer :: thread_id,vz,job
  character(len=80) :: atmp
  real(wp) :: percent,runtime

  integer :: nt,nrt
  real(wp),allocatable :: temps(:,:),et(:,:),ht(:,:),gt(:,:),stot(:,:)
  real(wp) :: ithr,sthr,fscal,rt
  character(len=:),allocatable :: emodel

  type(timer) :: profiler
  integer :: T,Tn  !> threads and threads per core
  logical :: nested
  real(wp),parameter :: big = 10e10

!>--- check if we have any calculation settings allocated
  if (env%calc%ncalculations < 1) then
    write (stdout,*) 'no calculations allocated'
    return
  end if

!>--- prepare calculation objects for parallelization (one per thread)
  call new_ompautoset(env,'auto_nested',nall,T,Tn)
  nested = env%omp_allow_nested
  call ompautoset_summary(env,'Hessians',T,Tn)

!>--- prepare objects for parallelization (one working copy per parallel job)
  allocate (calculations(T))!,source=env%calc)
  allocate (mols(T))
  nat3 = nat*3
  allocate (freqs(nat3,T),source=0.0_wp)
  allocate (hess(nat3,nat3,T),source=0.0_wp)
  do i = 1,T
    call calculations(i)%copy(env%calc)
    do j = 1,env%calc%ncalculations
      !calculations(i)%calcs(j) = env%calc%calcs(j)
      !>--- directories and io preparation
      ex = directory_exist(env%calc%calcs(j)%calcspace)
      if (.not.ex) then
        io = makedir(trim(env%calc%calcs(j)%calcspace))
      end if
      write (atmp,'(a,"_",i0)') sep,i
      calculations(i)%calcs(j)%calcspace = env%calc%calcs(j)%calcspace//trim(atmp)
      if (allocated(calculations(i)%calcs(j)%calcfile)) deallocate (calculations(i)%calcs(j)%calcfile)
      if (allocated(calculations(i)%calcs(j)%systemcall)) deallocate (calculations(i)%calcs(j)%systemcall)
      call calculations(i)%calcs(j)%printid(i,j)
    end do
    calculations(i)%pr_energies = .false.
    allocate (mols(i)%at(nat),mols(i)%xyz(3,nat))
  end do

!>--- thermo settings
  !> inversion threshold
  ithr = env%thermo%ithr
  !> frequency scaling factor
  fscal = env%thermo%fscal
  !> RR-HO interpolation (or cut-off)
  sthr = env%thermo%sthr
  !> Svib model
  emodel = env%thermo%emodel
  if (.not.allocated(env%thermo%temps)) then
    call env%thermo%get_temps()
  end if
  if (present(gt_out)) then
    nt = env%thermo%ntemps
    allocate (temps(nt,T),et(nt,T),ht(nt,T),gt(nt,T),stot(nt,T),source=0.0_wp)
    do i = 1,T
      temps(:,i) = env%thermo%temps(:)
    end do
    rt = env%thermo%get_close_rt(nrt)
  else
    nt = 1
    allocate (temps(nt,T),et(nt,T),ht(nt,T),gt(nt,T),stot(nt,T),source=0.0_wp)
    rt = env%thermo%get_close_rt(nrt)
    temps = rt
    nrt = 1
  end if

!>--- printout directions and timer initialization
  pr = .false. !> stdout printout
  wr = .false. !> write crestopt.log.xyz
  call profiler%init(1)
  call profiler%start(1)

!>--- initialize progress bar
  call progress_init(env%ps,nall,width=50,prefix=" ↳ ", &
    &                suffix="",show_time=.true.,show_eta=.false.)
  call progress_update(env%ps,0,nall)

!>--- shared variables
  allocate (grads(3,nat,T),source=0.0_wp)
  c = 0  !> counter of successfull optimizations
  k = 0  !> counter of total optimization (fail+success)
  z = 0  !> counter to perform optimization in right order (1...nall)
  eread(:) = 0.0_wp
  grads(:,:,:) = 0.0_wp
!>--- pre-start server-based calculators before forking OMP threads
  call preinit_mlip_parallel(calculations,T)
!>--- loop over ensemble
  !$omp parallel &
  !$omp shared(env,calculations,nat,nall,at,xyz,eread,grads,c,k,z,pr,wr,nrt) &
  !$omp shared(mols,nested,Tn,freqs,hess,temps,et,ht,gt,stot,nat3,ithr,fscal,sthr,nt,emodel)
  !$omp single
  do i = 1,nall

    call initsignal()
    vz = i
    !$omp task firstprivate( vz ) private(i,j,job,energy,io,thread_id,zcopy)
    call initsignal()

    !>--- OpenMP nested region threads
    if (nested) call ompmklset(Tn)

    thread_id = OMP_GET_THREAD_NUM()
    job = thread_id+1
    !>--- modify calculation spaces
    !$omp critical
    z = z+1
    zcopy = z
    mols(job)%nat = nat
    mols(job)%at(:) = at(:)
    mols(job)%xyz(:,:) = xyz(:,:,z)
    !$omp end critical

    !>-- engery+gradient call first, for setup
    call engrad(mols(job),calculations(job),energy,grads(:,:,job),io)
    !>-- then, numerical hessian
    call numhess1(mols(job)%nat,mols(job)%at,mols(job)%xyz, &
                  calculations(job),hess(:,:,job),io)
    !!$omp critical
    if (io .eq. 0) then
      call prj_mw_hess(mols(job)%nat,mols(job)%at,nat3,mols(job)%xyz,hess(:,:,job))
      !>-- Computes the Frequencies
      call frequencies(mols(job)%nat,mols(job)%at,mols(job)%xyz, &
                       nat3,hess(:,:,job),freqs(:,job),io)
    end if

    if (io .eq. 0) then
      call calcthermo(mols(job)%nat,mols(job)%at,mols(job)%xyz,   &
                      freqs(:,job),.false.,ithr,fscal,sthr,nt,    &
                      temps(:,job),et(:,job),ht(:,job),gt(:,job), &
                      stot(:,job),emodel=emodel)
    end if
    !!$omp end critical

    !$omp critical
    if (io == 0) then
      c = c+1
      eread(zcopy) = gt(nrt,job)
      if (present(gt_out))   gt_out(zcopy,:)   = gt(:,job)
      if (present(stot_out)) stot_out(zcopy,:) = stot(:,job)
    else
      eread(zcopy) = big
    end if
    k = k+1
    call progress_update(env%ps,k,nall)
    !$omp end critical
    !$omp end task
  end do
  !$omp taskwait
  !$omp end single
  !$omp end parallel

!>--- finalize progress printout
  call progress_finish(env%ps)

!>--- stop timer
  call profiler%stop(1)

!>--- prepare some summary printout
  percent = float(c)/float(nall)*100.0_wp
  write (atmp,'(f5.1,a)') percent,'% success)'
  write (stdout,'(">",1x,i0,a,i0,a,a)') c,' of ',nall,' structures successfully evaluated (', &
  &     trim(adjustl(atmp))
  write (atmp,'(">",1x,a,i0,a)') 'Total runtime for ',nall,' frequency calculations:'
  call profiler%write_timing(stdout,1,trim(atmp),.true.)
  runtime = profiler%get(1)
  write (atmp,'(f16.3,a)') runtime/real(nall,wp),' sec'
  write (stdout,'(a,a,a)') '> Corresponding to approximately ',trim(adjustl(atmp)), &
  &                       ' per processed structure'

  deallocate (grads)
  call profiler%clear()
  deallocate (calculations)
  if (allocated(mols)) deallocate (mols)
  if (allocated(freqs)) deallocate (freqs)
  if (allocated(hess)) deallocate (hess)
  return
end subroutine crest_hessloop

!========================================================================================!
!========================================================================================!
!> Routines for concurrent geometry optimization
!========================================================================================!
!========================================================================================!
subroutine crest_oloop_struc(env,nall,structures,dump,customcalc,eread,silent)
!*******************************************************************************
!* subroutine crest_oloop_struc
!* Concurrent geometry optimizations for a list of coord objects.
!* Optimized geometries and energies are written back into structures;
!* the optional eread array, if present, additionally receives the energies.
!* Each coord carries its own %lat/%chrg/%uhf, so periodic and heterogeneous
!* systems are handled.
!*
!* dump       - dump an ensemble file (NOT in the input order)
!* customcalc - customized (optional) calculation level data
!* silent     - suppress the progress bar and summary printout (optional,
!*              default .false.); used when the caller drives many small
!*              batches and prints its own progress (e.g. TTConf-light)
!*
!* IMPORTANT: structures xyz must be in Bohr(!)
!******************************************************************************
  use crest_parameters,only:wp,stdout,sep
  use crest_calculator
  use omp_lib
  use crest_data
  use strucrd
  use optimize_module
  use iomod,only:makedir,directory_exist,remove
  use term_ui,only:progress_init,progress_update,progress_finish
  implicit none
  type(systemdata),target,intent(inout) :: env
  integer,intent(in) :: nall
  type(coord),intent(inout) :: structures(nall)
  logical,intent(in) :: dump
  type(calcdata),intent(in),target,optional :: customcalc
  real(wp),intent(inout),optional :: eread(nall)
  logical,intent(in),optional :: silent
  logical :: quiet

  type(coord),allocatable :: mols(:)
  type(coord),allocatable :: molsnew(:)
  integer :: i,j,io,ich,ich2,c,k,z,zcopy
  logical :: pr,wr,ex
  type(calcdata),allocatable :: calculations(:)
  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)
  integer :: thread_id,vz,job
  character(len=80) :: atmp
  real(wp) :: percent,runtime
  type(calcdata),pointer :: mycalc
  type(timer) :: profiler
  integer :: T,Tn  !> threads and threads per core
  logical :: nested

!>--- check which calc to use
  if (present(customcalc)) then
    mycalc => customcalc
  else
    mycalc => env%calc
  end if

!>--- silent mode? (suppress progress bar + summary printout)
  quiet = .false.
  if (present(silent)) quiet = silent

!>--- check if we have any calculation settings allocated
  if (mycalc%ncalculations < 1) then
    write (stdout,*) 'no calculations allocated'
    return
  end if

!>--- prepare calculation objects for parallelization (one per thread)
  call new_ompautoset(env,'auto_nested',nall,T,Tn)
  nested = env%omp_allow_nested
  if (.not.quiet) call ompautoset_summary(env,'optimizations',T,Tn)

!>--- prepare objects for parallelization
  allocate (calculations(T))
  allocate (mols(T),molsnew(T))
  do i = 1,T
    call calculations(i)%copy(mycalc)
    do j = 1,mycalc%ncalculations
      !>--- directories and io preparation
      ex = directory_exist(mycalc%calcs(j)%calcspace)
      if (.not.ex) then
        io = makedir(trim(mycalc%calcs(j)%calcspace))
      end if
      if (calculations(i)%calcs(j)%id == jobtype%tblite) then
        calculations(i)%optnewinit = .true.
      end if
      write (atmp,'(a,"_",i0)') sep,i
      calculations(i)%calcs(j)%calcspace = mycalc%calcs(j)%calcspace//trim(atmp)
      if (allocated(calculations(i)%calcs(j)%calcfile)) deallocate (calculations(i)%calcs(j)%calcfile)
      if (allocated(calculations(i)%calcs(j)%systemcall)) deallocate (calculations(i)%calcs(j)%systemcall)
      call calculations(i)%calcs(j)%printid(i,j)
    end do
    calculations(i)%pr_energies = .false.
  end do

!>--- printout directions and timer initialization
  pr = .false. !> stdout printout
  wr = .false. !> write crestopt.log.xyz
  if (dump) then
    open (newunit=ich,file=ensemblefile)
    open (newunit=ich2,file=ensembleelog)
  end if
  call profiler%init(1)
  call profiler%start(1)

!>--- initialize progress bar
  if (.not.quiet) then
    call progress_init(env%ps,nall,width=50,prefix=" ↳ ", &
      &                suffix="",show_time=.true.,show_eta=.false.)
    call progress_update(env%ps,0,nall)
  end if

!>--- shared variables
  c = 0  !> counter of successfull optimizations
  k = 0  !> counter of total optimization (fail+success)
  z = 0  !> counter to perform optimization in right order (1...nall)
!>--- pre-start server-based calculators before forking OMP threads
  call preinit_mlip_parallel(calculations,T)
!>--- loop over ensemble
  !$omp parallel &
  !$omp shared(env,calculations,nall,structures,c,k,z,pr,wr,dump) &
  !$omp shared(ich,ich2,mols,molsnew,nested,Tn)
  !$omp single
  do i = 1,nall

    call initsignal()
    vz = i
    !$omp task firstprivate( vz ) private(j,job,energy,grad,io,atmp,gnorm,thread_id,zcopy)
    call initsignal()

    !>--- OpenMP nested region threads
    if (nested) call ompmklset(Tn)

    thread_id = OMP_GET_THREAD_NUM()
    job = thread_id+1
    !>--- deep-copy this structure into the thread-local working mol
    !$omp critical
    z = z+1
    zcopy = z
    call mols(job)%copy(structures(zcopy))
    !$omp end critical

    allocate (grad(3,mols(job)%nat),source=0.0_wp)

    !>-- geometry optimization
    call optimize_geometry(mols(job),molsnew(job),calculations(job),energy,grad,pr,wr,io)

    !$omp critical
    if (io == 0) then
      !>--- successful optimization (io==0)
      c = c+1
      structures(zcopy)%xyz = molsnew(job)%xyz
      structures(zcopy)%energy = energy
      if (dump) then
        gnorm = norm2(grad)
        write (atmp,'(1x,"energy=",f16.10,1x,"g norm=",f12.8)') energy,gnorm
        molsnew(job)%comment = trim(atmp)
        call molsnew(job)%append(ich)
        call calc_eprint(calculations(job),energy,calculations(job)%etmp,gnorm,ich2)
      end if
    else if (io == calculations(job)%maxcycle.and.calculations(job)%anopt) then
      !>--- allow partial optimization?
      c = c+1
      structures(zcopy)%xyz = molsnew(job)%xyz
      structures(zcopy)%energy = energy
    else
      structures(zcopy)%energy = 1.0_wp
    end if
    k = k+1
    !>--- print progress
    if (.not.quiet) call progress_update(env%ps,k,nall)
    !$omp end critical

    deallocate (grad)
    !$omp end task
  end do
  !$omp taskwait
  !$omp end single
  !$omp end parallel

!>--- finalize progress printout
  if (.not.quiet) call progress_finish(env%ps)

!>--- stop timer
  call profiler%stop(1)

!>--- energies are stored in the structures; optionally also return them
  if (present(eread)) then
    do i = 1,nall
      eread(i) = structures(i)%energy
    end do
  end if

!>--- prepare some summary printout
  if (.not.quiet) then
    percent = float(c)/float(nall)*100.0_wp
    write (atmp,'(f5.1,a)') percent,'% success)'
    write (stdout,'(">",1x,i0,a,i0,a,a)') c,' of ',nall,' structures successfully optimized (', &
    &     trim(adjustl(atmp))
    write (atmp,'(">",1x,a,i0,a)') 'Total runtime for ',nall,' optimizations:'
    call profiler%write_timing(stdout,1,trim(atmp),.true.)
    runtime = profiler%get(1)
    write (atmp,'(f16.3,a)') runtime/real(nall,wp),' sec'
    write (stdout,'(a,a,a)') '> Corresponding to approximately ',trim(adjustl(atmp)), &
    &                       ' per processed structure'
  end if

!>--- close files (if they are open)
  if (dump) then
    close (ich)
    close (ich2)
  end if

  call profiler%clear()
  deallocate (calculations)
  if (allocated(mols)) deallocate (mols)
  if (allocated(molsnew)) deallocate (molsnew)
  return
end subroutine crest_oloop_struc

!========================================================================================!
!> Flat-array adapter for crest_oloop. Marshals the (nat,nall,at,xyz) ensemble
!> into a coord list, runs the optimization loop, and copies the optimized
!> coordinates and energies back. Kept for the legacy (non-periodic) callers;
!> new code should use the coord-list crest_oloop directly.
!========================================================================================!
subroutine crest_oloop_xyz(env,nat,nall,at,xyz,eread,dump,customcalc)
  use crest_parameters,only:wp
  use crest_calculator
  use crest_data
  use strucrd
  use parallel_interface,only:crest_oloop
  implicit none
  type(systemdata),target,intent(inout) :: env
  integer,intent(in) :: nat,nall
  integer,intent(in) :: at(nat)
  real(wp),intent(inout) :: xyz(3,nat,nall)
  real(wp),intent(inout) :: eread(nall)
  logical,intent(in) :: dump
  type(calcdata),intent(in),target,optional :: customcalc

  type(coord),allocatable :: structures(:)
  integer :: i

  allocate (structures(nall))
  do i = 1,nall
    structures(i)%nat = nat
    structures(i)%at = at
    structures(i)%xyz = xyz(1:3,1:nat,i)
  end do

  if (present(customcalc)) then
    call crest_oloop(env,nall,structures,dump,customcalc,eread)
  else
    call crest_oloop(env,nall,structures,dump,eread=eread)
  end if

  do i = 1,nall
    xyz(1:3,1:nat,i) = structures(i)%xyz(1:3,1:nat)
  end do
  deallocate (structures)
  return
end subroutine crest_oloop_xyz

!========================================================================================!
!========================================================================================!
!> Routines for parallel MDs
!========================================================================================!
!========================================================================================!

subroutine crest_search_multimd(env,mol,mddats,nsim)
!*****************************************************
!* subroutine crest_search_multimd
!* this runs #nsim MDs on the same structure (mol)
!*****************************************************
  use crest_parameters,only:wp,stdout,sep
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use iomod,only:makedir,directory_exist,remove
  use omp_lib
  implicit none
  type(systemdata),intent(inout) :: env
  type(mddata) :: mddats(nsim)
  integer :: nsim
  type(coord) :: mol
  type(coord),allocatable :: moltmps(:)
  integer :: i,j,io,ich
  logical :: pr,ex,nested
  integer :: T,Tn
  real(wp) :: percent
  character(len=80) :: atmp
  character(len=*),parameter :: mdir = 'MDFILES'

  type(calcdata),allocatable :: calculations(:)
  integer :: vz,job,thread_id
  real(wp) :: etmp
  real(wp),allocatable :: grdtmp(:,:)
  type(timer) :: profiler
!===========================================================!
!>--- check if we have any MD & calculation settings allocated
  if (.not.env%mddat%requested) then
    write (stdout,*) 'MD requested, but no MD settings present.'
    return
  else if (env%calc%ncalculations < 1) then
    write (stdout,*) 'MD requested, but no calculation settings present.'
    return
  end if

!>--- prepare calculation containers for parallelization (one per thread)
  call new_ompautoset(env,'auto_nested',nsim,T,Tn)
  nested = env%omp_allow_nested
  call ompautoset_summary(env,'MTD/MD runs',T,Tn)

  allocate (calculations(T),source=env%calc)
  allocate (moltmps(T),source=mol)
  allocate (grdtmp(3,mol%nat),source=0.0_wp)
  do i = 1,T
    moltmps(i)%nat = mol%nat
    moltmps(i)%at = mol%at
    moltmps(i)%xyz = mol%xyz
    do j = 1,env%calc%ncalculations
      calculations(i)%calcs(j) = env%calc%calcs(j)
      !>--- directories and io preparation
      ex = directory_exist(env%calc%calcs(j)%calcspace)
      if (.not.ex) then
        io = makedir(trim(env%calc%calcs(j)%calcspace))
      end if
      write (atmp,'(a,"_",i0)') sep,i
      calculations(i)%calcs(j)%calcspace = env%calc%calcs(j)%calcspace//trim(atmp)
      if (allocated(calculations(i)%calcs(j)%calcfile)) deallocate (calculations(i)%calcs(j)%calcfile)
      if (allocated(calculations(i)%calcs(j)%systemcall)) deallocate (calculations(i)%calcs(j)%systemcall)
      call calculations(i)%calcs(j)%printid(i,j)
    end do
    calculations(i)%pr_energies = .false.
    !>--- initialize the calculations
    call engrad(moltmps(i),calculations(i),etmp,grdtmp,io)
  end do

  !>--- other settings
  pr = .false.
  call profiler%init(nsim)

!>--- pre-start server-based calculators before forking OMP threads
  call preinit_mlip_parallel(calculations,T)
  !>--- run the MDs
  !$omp parallel &
  !$omp shared(env,calculations,mddats,mol,pr,percent,ich, nsim, moltmps, nested,Tn) &
  !!$omp single
  !$omp private(vz,i,job,thread_id,io,ex)
  !$omp do
  do i = 1,nsim

    call initsignal()
    vz = i

    !>--- OpenMP nested region threads
    if (nested) call ompmklset(Tn)

    !!$omp task firstprivate( vz ) private( job,thread_id,io,ex )
    call initsignal()

    thread_id = OMP_GET_THREAD_NUM()
    job = thread_id+1
    !$omp critical
    moltmps(job)%nat = mol%nat
    moltmps(job)%at = mol%at
    moltmps(job)%xyz = mol%xyz
    !>--- carry the lattice (PBC) through to the per-thread MD copy
    if (allocated(mol%lat)) moltmps(job)%lat = mol%lat
    !$omp end critical
    !>--- startup printout (thread safe)
    call parallel_md_block_printout(mddats(vz),vz)

    !>--- the acutal MD call with timing
    call profiler%start(vz)
    call dynamics(moltmps(job),mddats(vz),calculations(job),pr,io)
    call profiler%stop(vz)

    !>--- finish printout (thread safe)
    call parallel_md_finish_printout(mddats(vz),vz,io,profiler)
    !!$omp end task
  end do
  !!$omp taskwait
  !$omp end parallel

  !>--- collect trajectories into one
  call collect(nsim,mddats)

  call profiler%clear()
  deallocate (calculations)
  if (allocated(moltmps)) deallocate (moltmps)
  return
contains
  subroutine collect(n,mddats)
    implicit none
    integer :: n
    type(mddata) :: mddats(n)
    logical :: ex
    integer :: i,io,ich,ich2
    character(len=:),allocatable :: atmp
    character(len=256) :: btmp
    open (newunit=ich,file='crest_dynamics.trj.xyz')
    do i = 1,n
      atmp = mddats(i)%trajectoryfile
      inquire (file=atmp,exist=ex)
      if (ex) then
        open (newunit=ich2,file=atmp)
        io = 0
        do while (io == 0)
          read (ich2,'(a)',iostat=io) btmp
          if (io == 0) then
            write (ich,'(a)') trim(btmp)
          end if
        end do
        close (ich2)
      end if
    end do
    close (ich)
    return
  end subroutine collect
end subroutine crest_search_multimd

!========================================================================================!
subroutine crest_search_multimd_init(env,mol,mddat,nsim)
!*******************************************************
!* subroutine crest_search_multimd_init
!* This routine will initialize a copy of env%mddat
!* and save it to the local mddat. If we are about to
!* run RMSD metadynamics, the required number of
!* simulations (#nsim) is returned
!*******************************************************
  use crest_parameters,only:wp,stdout
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use iomod,only:makedir,directory_exist,remove
  use omp_lib
  implicit none
  type(systemdata),intent(inout) :: env
  type(mddata) :: mddat
  type(coord) :: mol
  integer,intent(inout) :: nsim
  integer :: i,io
  logical :: pr
!=======================================================!
  type(calcdata),target :: calc
  type(shakedata) :: shk

  real(wp) :: energy
  real(wp),allocatable :: grad(:,:)
  character(len=*),parameter :: mdir = 'MDFILES'
!======================================================!

  !>--- check if we have any MD & calculation settings allocated
  mddat = env%mddat
  if (.not.mddat%requested) then
    write (stdout,*) 'MD requested, but no MD settings present.'
    return
  else if (env%calc%ncalculations < 1) then
    write (stdout,*) 'MD requested, but no calculation settings present.'
    return
  end if

  !>--- init SHAKE?
  if (mddat%shake) then
    if (allocated(env%ref%wbo)) then
      shk%wbo = env%ref%wbo
    else
      calc = env%calc
      calc%calcs(1)%rdwbo = .true.
      allocate (grad(3,mol%nat),source=0.0_wp)
      call engrad(mol,calc,energy,grad,io)
      deallocate (grad)
      calc%calcs(1)%rdwbo = .false.

      shk%shake_mode = env%mddat%shk%shake_mode
      call move_alloc(calc%calcs(1)%wbo,shk%wbo)
    end if

    if (calc%nfreeze > 0) then
      shk%freezeptr => calc%freezelist
    else
      nullify (shk%freezeptr)
    end if

    shk%shake_mode = env%shake
    mddat%shk = shk
    call init_shake(mol%nat,mol%at,mol%xyz,mddat%shk,pr)
    mddat%nshake = mddat%shk%ncons
  end if
  !>--- complete real-time settings to steps
  call mdautoset(mddat,io)

  !>--- (optional)  MTD initialization
  if (nsim < 0) then
    mddat%simtype = type_mtd  !>-- set runtype to MTD

    call defaultGF(env)
    write (stdout,*) 'list of applied metadynamics Vbias parameters:'
    do i = 1,env%nmetadyn
      write (stdout,'(''$metadyn '',f10.5,f8.3,i5)') env%metadfac(i),env%metadexp(i)
    end do
    write (stdout,*)

    !>--- how many simulations
    nsim = env%nmetadyn
  end if

  return
end subroutine crest_search_multimd_init

!========================================================================================!
subroutine crest_search_multimd_init2(env,mddats,nsim)
  use crest_parameters,only:wp,stdout,sep
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use iomod,only:makedir,directory_exist,remove
  use omp_lib
  implicit none
  type(systemdata),intent(inout) :: env
  type(mddata) :: mddats(nsim)
  integer :: nsim
  integer :: i,io,j
  logical :: ex
!========================================================!
  type(mtdpot),allocatable :: mtds(:)

  character(len=80) :: atmp
  character(len=*),parameter :: mdir = 'MDFILES'

  !>--- parallel MD setup
  ex = directory_exist(mdir)
  if (ex) then
    call rmrf(mdir)
  end if
  io = makedir(mdir)
  do i = 1,nsim
    mddats(i)%md_index = i
    write (atmp,'(a,i0,a)') 'crest_',i,'.trj'
    mddats(i)%trajectoryfile = mdir//sep//trim(atmp)
    write (atmp,'(a,i0,a)') 'crest_',i,'.mdrestart'
    mddats(i)%restartfile = mdir//sep//trim(atmp)
  end do

  allocate (mtds(nsim))
  do i = 1,nsim
    if (mddats(i)%simtype == type_mtd) then
      mtds(i)%kpush = env%metadfac(i)
      mtds(i)%alpha = env%metadexp(i)
      mtds(i)%cvdump_fs = float(env%mddump)
      mtds(i)%mtdtype = cv_rmsd

      mddats(i)%npot = 1
      allocate (mddats(i)%mtd(1),source=mtds(i))
      allocate (mddats(i)%cvtype(1),source=cv_rmsd)
      !> if necessary exclude atoms from RMSD bias
      if (sum(env%includeRMSD) /= env%ref%nat) then
        if (.not.allocated(mddats(i)%mtd(1)%atinclude)) &
        & allocate (mddats(i)%mtd(1)%atinclude(env%ref%nat),source=.true.)
        do j = 1,env%ref%nat
          if (env%includeRMSD(j) .ne. 1) mddats(i)%mtd(1)%atinclude(j) = .false.
        end do
      end if
    end if
  end do
  if (allocated(mtds)) deallocate (mtds)

  return
end subroutine crest_search_multimd_init2

!========================================================================================!
subroutine crest_search_multimd2(env,mols,mddats,nsim)
!*******************************************************************
!* subroutine crest_search_multimd2
!* this runs #nsim MDs on #nsim selected different structures (mols)
!*******************************************************************
  use crest_parameters,only:wp,stdout,sep
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  use iomod,only:makedir,directory_exist,remove
  use omp_lib
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env
  type(mddata) :: mddats(nsim)
  integer :: nsim
  type(coord) :: mols(nsim)
  type(coord),allocatable :: moltmps(:)
  integer :: i,j,io,ich
  logical :: pr,ex,nested
  integer :: T,Tn
  real(wp) :: percent
  character(len=80) :: atmp
  character(len=*),parameter :: mdir = 'MDFILES'

  type(calcdata),allocatable :: calculations(:)
  integer :: vz,job,thread_id
  type(timer) :: profiler
!===========================================================!
!>--- check if we have any MD & calculation settings allocated
  if (.not.env%mddat%requested) then
    write (stdout,*) 'MD requested, but no MD settings present.'
    return
  else if (env%calc%ncalculations < 1) then
    write (stdout,*) 'MD requested, but no calculation settings present.'
    return
  end if

!>--- prepare calculation objects for parallelization (one per thread)
  call new_ompautoset(env,'auto_nested',nsim,T,Tn)
  nested = env%omp_allow_nested
  call ompautoset_summary(env,'MTD/MD runs',T,Tn)

  allocate (calculations(T),source=env%calc)
  allocate (moltmps(T),source=mols(1))
  do i = 1,T
    do j = 1,env%calc%ncalculations
      calculations(i)%calcs(j) = env%calc%calcs(j)
      !>--- directories and io preparation
      ex = directory_exist(env%calc%calcs(j)%calcspace)
      if (.not.ex) then
        io = makedir(trim(env%calc%calcs(j)%calcspace))
      end if
      write (atmp,'(a,"_",i0)') sep,i
      calculations(i)%calcs(j)%calcspace = env%calc%calcs(j)%calcspace//trim(atmp)
      if (allocated(calculations(i)%calcs(j)%calcfile)) deallocate (calculations(i)%calcs(j)%calcfile)
      if (allocated(calculations(i)%calcs(j)%systemcall)) deallocate (calculations(i)%calcs(j)%systemcall)
      call calculations(i)%calcs(j)%printid(i,j)
    end do
    calculations(i)%pr_energies = .false.
  end do

!>--- other settings
  pr = .false.
  call profiler%init(nsim)

!>--- pre-start server-based calculators before forking OMP threads
  call preinit_mlip_parallel(calculations,T)
!>--- run the MDs
  !$omp parallel &
  !$omp shared(env,calculations,mddats,mols,pr,percent,ich, moltmps,profiler, nested,Tn)
  !$omp single
  do i = 1,nsim

    call initsignal()
    vz = i

    !>--- OpenMP nested region threads
    if (nested) call ompmklset(Tn)

    !$omp task firstprivate( vz ) private( job,thread_id,io,ex )
    call initsignal()

    thread_id = OMP_GET_THREAD_NUM()
    job = thread_id+1
    !$omp critical
    moltmps(job)%nat = mols(vz)%nat
    moltmps(job)%at = mols(vz)%at
    moltmps(job)%xyz = mols(vz)%xyz
    !>--- carry the per-structure lattice (PBC); moltmps was sourced from mols(1)
    if (allocated(mols(vz)%lat)) moltmps(job)%lat = mols(vz)%lat
    !$omp end critical
    !>--- startup printout (thread safe)
    call parallel_md_block_printout(mddats(vz),vz)

    !>--- the acutal MD call with timing
    call profiler%start(vz)
    call dynamics(moltmps(job),mddats(vz),calculations(job),pr,io)
    call profiler%stop(vz)

    !>--- finish printout (thread safe)
    call parallel_md_finish_printout(mddats(vz),vz,io,profiler)
    !$omp end task
  end do
  !$omp taskwait
  !$omp end single
  !$omp end parallel

!>--- collect trajectories into one
  call collect(nsim,mddats)

  call profiler%clear()
  deallocate (calculations)
  if (allocated(moltmps)) deallocate (moltmps)
  return
contains
  subroutine collect(n,mddats)
    implicit none
    integer :: n
    type(mddata) :: mddats(n)
    logical :: ex
    integer :: i,io,ich,ich2
    character(len=:),allocatable :: atmp
    character(len=256) :: btmp
    open (newunit=ich,file='crest_dynamics.trj.xyz')
    do i = 1,n
      atmp = mddats(i)%trajectoryfile
      inquire (file=atmp,exist=ex)
      if (ex) then
        open (newunit=ich2,file=atmp)
        io = 0
        do while (io == 0)
          read (ich2,'(a)',iostat=io) btmp
          if (io == 0) then
            write (ich,'(a)') trim(btmp)
          end if
        end do
        close (ich2)
      end if
    end do
    close (ich)
    return
  end subroutine collect
end subroutine crest_search_multimd2

!========================================================================================!
subroutine parallel_md_block_printout(MD,vz)
!***********************************************
!* subroutine parallel_md_block_printout
!* This will print information about the MD/MTD
!* simulation. The execution is omp threadsave
!***********************************************
  use crest_parameters,only:wp,stdout,sep
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  use iomod,only:to_str,drawbox
  implicit none
  type(mddata),intent(in) :: MD
  integer,intent(in) :: vz
  character(len=60) :: atmp
  integer,parameter :: bw = 54
  !$omp critical

  ! ── title ────────────────────────────────────────────────────
  if (MD%simtype == type_md) then
    write (atmp,'(a,1x,i3)') 'starting MD',vz
  else if (MD%simtype == type_mtd) then
    if (MD%cvtype(1) == cv_rmsd_static) then
      write (atmp,'(a,1x,i3)') 'starting static MTD',vz
    else
      write (atmp,'(a,1x,i4)') 'starting MTD',vz
    end if
  end if
  call drawbox(stdout,trim(atmp),width=bw,charset=4,ltab=1,procedual=0)
  call drawbox(stdout,trim(atmp),width=bw,charset=4,ltab=1,procedual=1)
  call drawbox(stdout,trim(atmp),width=bw,charset=4,ltab=1,procedual=3)

  ! ── simulation parameters ────────────────────────────────────
  write (stdout,'(1x,"│   MD simulation time   :",f8.1," ps",16x,"│")') MD%length_ps
  write (stdout,'(1x,"│   target T             :",f8.1," K",17x,"│")') MD%tsoll
  write (stdout,'(1x,"│   timestep dt          :",f8.1," fs",16x,"│")') MD%tstep
  write (stdout,'(1x,"│   dump interval(trj)   :",f8.1," fs",16x,"│")') MD%dumpstep
  if (MD%shake.and.MD%shk%shake_mode > 0) then
    if (MD%shk%shake_mode == 2) then
      write (stdout,'(1x,"│   SHAKE algorithm      :",a5," (all bonds)",10x,"│")') to_str(MD%shake)
    else
      write (stdout,'(1x,"│   SHAKE algorithm      :",a5," (H only)",13x,"│")') to_str(MD%shake)
    end if
  end if
  if (allocated(MD%active_potentials)) then
    write (stdout,'(1x,"│   active potentials    :",i4," potential(s)",10x,"│")') size(MD%active_potentials,1)
  end if
  if (MD%simtype == type_mtd) then
    if (MD%cvtype(1) == cv_rmsd) then
      write (stdout,'(1x,"│   dump interval(Vbias) :",f8.2," ps",16x,"│")') &
          & MD%mtd(1)%cvdump_fs/1000.0_wp
    end if
    write (stdout,'(1x,"│   Vbias prefactor (k)  :",f8.4," Eh",16x,"│")') MD%mtd(1)%kpush
    if (MD%cvtype(1) == cv_rmsd.or.MD%cvtype(1) == cv_rmsd_static) then
      write (stdout,'(1x,"│   Vbias exponent (α)   :",f8.4," bohr⁻²",12x,"│")') MD%mtd(1)%alpha
    else
      write (stdout,'(1x,"│   Vbias exponent (α)   :",f8.4,19x,"│")') MD%mtd(1)%alpha
    end if
    if (allocated(MD%mtd(1)%atinclude)) then
      write (stdout,'(1x,"│   # active atoms       :",i9," atoms",12x,"│")') count(MD%mtd(1)%atinclude,1)
    end if
  end if

  call drawbox(stdout,'',width=bw,charset=4,ltab=1,procedual=2)

  !$omp end critical

end subroutine parallel_md_block_printout

subroutine parallel_md_finish_printout(MD,vz,io,profiler)
!*******************************************
!* subroutine parallel_md_finish_printout
!* This will print information termination
!* info about the MD/MTD simulation
!*******************************************
  use crest_parameters,only:wp,stdout,sep
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  implicit none
  type(mddata),intent(in) :: MD
  integer,intent(in) :: vz,io
  type(timer),intent(inout) :: profiler
  character(len=40) :: atmp
  character(len=80) :: btmp

  !$omp critical

  if (MD%simtype == type_mtd) then
    if (MD%cvtype(1) == cv_rmsd_static) then
      write (atmp,'(a)') '*sMTD'
    else
      write (atmp,'(a)') '*MTD'
    end if
  else
    write (atmp,'(a)') '*MD'
  end if
  if (io == 0) then
    write (btmp,'(a,1x,i3,a)') trim(atmp),vz,' completed successfully'
  else
    write (btmp,'(a,1x,i3,a)') trim(atmp),vz,' terminated EARLY'
  end if
  call profiler%write_timing(stdout,vz,trim(btmp))

  !$omp end critical

end subroutine parallel_md_finish_printout
!========================================================================================!
!========================================================================================!
