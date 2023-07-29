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

!> a collection of routines to set up parallel runs of
!> MDs and optimizations.


!========================================================================================!

subroutine crest_oloop(env,nat,nall,at,xyz,eread,dump)
  use crest_parameters, only: wp,stdout
  use omp_lib
  use crest_data
  use strucrd
  use calc_type
  use calc_module
  use optimize_module
  use iomod,only:makedir,directory_exist,remove
  implicit none
  type(systemdata),intent(inout) :: env
  real(wp),intent(inout) :: xyz(3,nat,nall)
  integer,intent(in)  :: at(nat)
  real(wp),intent(inout) :: eread(nall)
  integer,intent(in) :: nat,nall
  logical,intent(in) :: dump

  type(coord),allocatable :: mols(:)
  type(coord),allocatable :: molsnew(:)
  integer :: i,j,k,l,io,ich,ich2,c,z,job_id
  logical :: pr,wr,ex
  type(calcdata),allocatable :: calculations(:)
  type(calcdata) :: calc
  integer :: T  !> number of parallel running instances
  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)
  integer :: thread_id,vz,job
  character(len=80) :: atmp
  real(wp) :: percent
  character(len=52) :: bar

  !>--- check if we have any calculation settings allocated
  calc = env%calc
  if (calc%ncalculations < 1) then
    write (stdout,*) 'no calculations allocated'
    return
  end if

  !>--- prepare objects for parallelization
  T = env%threads
  allocate (calculations( T ),source=env%calc)
  allocate (mols( T ), molsnew( T ))
  do i = 1,T
    do j = 1,calc%ncalculations
      calculations(i)%calcs(j) = env%calc%calcs(j)
      !>--- directories
      ex = directory_exist(env%calc%calcs(j)%calcspace)
      if (.not. ex) then
        io = makedir(trim(env%calc%calcs(j)%calcspace))
      end if
      write (atmp,'(a,"_",i0)') sep,i
      calculations(i)%calcs(j)%calcspace = env%calc%calcs(j)%calcspace//trim(atmp)
    end do
    calculations(i)%pr_energies = .false.
    allocate (mols(i)%at(nat), mols(i)%xyz(3,nat))
    allocate (molsnew(i)%at(nat), molsnew(i)%xyz(3,nat))
  end do

  !>--- printout directions
  pr = .false. !> stdout printout
  wr = .false. !> write crestopt.log
  if (dump)then
     open (newunit=ich,file=ensemblefile)
     open (newunit=ich2,file=ensembleelog)
  endif
  if (env%niceprint) then
      percent = 0.0_wp
      call progbar(percent,bar)
      call printprogbar(percent,bar)
  endif
  !>--- shared variables
  allocate (grad(3,nat),source=0.0_wp)
   
  c = 0
  k = 0
  z = 0
  eread(:) = 0.0_wp 
  !>--- loop over ensemble
  !$omp parallel &
  !$omp shared(env,calculations,nat,nall,at,xyz,c,k,z,pr,wr,dump,percent,bar,ich,ich2,mols,molsnew)
  !$omp single
  do i = 1,nall

    call initsignal()
    vz = i
    !$omp task firstprivate( vz ) private(j,job,calc,energy,grad,io,atmp,gnorm,thread_id)
    call initsignal()

    thread_id = OMP_GET_THREAD_NUM()
    job = thread_id + 1
    !>--- modify calculation spaces
    !$omp critical
    z = z+1
    mols(job)%nat = nat
    mols(job)%at(:) = at(:)
    mols(job)%xyz(:,:) = xyz(:,:,z)

    molsnew(job)%nat = nat
    molsnew(job)%at(:) = at(:)
    molsnew(job)%xyz(:,:) = xyz(:,:,z)
    !$omp end critical

    !>--- first energy&gradient calculation
    call engrad(mols(job),calculations(job),energy,grad,io)

    !>-- geopetry optimization
    call optimize_geometry(mols(job),molsnew(job),calculations(job),energy,grad,pr,wr,io)

    !$omp critical
    if (io == 0) then
      !>--- successful optimization (io==0)
      c = c + 1
      if (dump) then
        gnorm = norm2(grad)
        write (atmp,'(1x,"Etot=",f16.10,1x,"g norm=",f12.8)') energy,gnorm
        molsnew(job)%comment = trim(atmp)
        call molsnew(job)%append(ich)
        call calc_eprint(calculations(job),energy,calculations(job)%etmp,ich2)
      end if
    end if
    k = k + 1
    if (env%niceprint) then
      percent = float(k) / float(nall) * 100.0_wp
      call progbar(percent,bar)
      call printprogbar(percent,bar)
    else
      write (stdout,'(1x,i0)',advance='no') k
      flush (stdout)
    end if

    !$omp end critical
    !$omp end task
  end do
  !$omp taskwait
  !$omp end single
  !$omp end parallel

  if (.not. env%niceprint) then
    write (stdout,'(/,1x,a)') 'done.'
  else
    write (stdout,*)
  end if
 
  write(stdout,'(1x,i0,a,i0,a)')c,' of ',nall,' structures successfully optimized.' 

  if (dump)then
     close (ich)
     close(ich2)
  endif 

  deallocate (grad)
  deallocate (calculations)
  if(allocated(mols))deallocate(mols)
  if(allocated(molsnew))deallocate(molsnew) 
  return
end subroutine crest_oloop

!========================================================================================!
!========================================================================================!

subroutine crest_search_multimd(env,mol,mddats,nsim)
!*****************************************************
!* this runs nsim MDs on the same structure (mol)
!*****************************************************
  use crest_parameters,only:wp,stdout
  use crest_data
  use strucrd
  use calc_type
  use calc_module
  use dynamics_module
  use shake_module
  use iomod,only:makedir,directory_exist,remove
  use omp_lib
  implicit none
  type(systemdata),intent(inout) :: env
  type(mddata) :: mddats(nsim)
  integer :: nsim
  type(coord) :: mol
  type(coord),allocatable :: moltmps(:)
  integer :: i,j,io,ich
  logical :: pr,ex
  type(calcdata) :: calc
  type(mddata) :: mddat
  integer :: T
  real(wp) :: percent
  character(len=52) :: bar
  character(len=80) :: atmp
  character(len=*),parameter :: mdir = 'MDFILES' 
 
  type(calcdata),allocatable :: calculations(:)
  integer :: vz,job,thread_id
  real(wp) :: etmp
  real(wp),allocatable :: grdtmp(:,:)
!===========================================================!


  !>--- set threads
  if (env%autothreads) then
    call ompautoset(env%threads,7,env%omp,env%MAXRUN,nsim)
  end if

  !>--- check if we have any MD & calculation settings allocated
  mddat = env%mddat
  calc = env%calc
  if (.not. mddat%requested) then
    write (stdout,*) 'MD requested, but no MD settings present.'
    return
  else if (calc%ncalculations < 1) then
    write (stdout,*) 'MD requested, but no calculation settings present.'
    return
  end if

  !>--- prepare calculation objects for parallelization (one per thread)
  T = env%threads
  allocate (calculations( T ),source=env%calc)
  allocate (moltmps( T ), source=mol )
  allocate (grdtmp(3,mol%nat), source =0.0_wp) 
  do i = 1,T
    moltmps(i)%nat = mol%nat
    moltmps(i)%at = mol%at
    moltmps(i)%xyz = mol%xyz
    do j = 1,calc%ncalculations
      calculations(i)%calcs(j) = env%calc%calcs(j)
      !>--- directories
      ex = directory_exist(env%calc%calcs(j)%calcspace)
      if (.not. ex) then
        io = makedir(trim(env%calc%calcs(j)%calcspace))
      end if
      write (atmp,'(a,"_",i0)') sep,i
      calculations(i)%calcs(j)%calcspace = env%calc%calcs(j)%calcspace//trim(atmp)
    end do
    calculations(i)%pr_energies = .false.
    !>--- initialize the calculations
    !call calculations(i)%info(stdout)
    call engrad(moltmps(i),calculations(i),etmp,grdtmp,io)
  end do

  !>--- other settings
  pr = .false.

  !>--- run the MDs
  !$omp parallel &
  !$omp shared(env,calculations,mddats,mol,pr,percent,bar,ich, moltmps)
  !$omp single
  do i = 1,nsim

    call initsignal()
    vz = i
    !$omp task firstprivate( vz ) private( job,thread_id,io,ex )
    call initsignal()

    thread_id = OMP_GET_THREAD_NUM()
    job = thread_id + 1
    !$omp critical
    moltmps(job)%nat = mol%nat
    moltmps(job)%at = mol%at
    moltmps(job)%xyz = mol%xyz

    if(mddats(vz)%simtype == type_md)then
    write (stdout,'(a,i4,a)') 'Starting MD',vz,' with the settings:'
    else if(mddats(vz)%simtype == type_mtd)then
    write (stdout,'(a,i4,a)') 'Starting MTD',vz,' with the settings:'
    endif 
    write (stdout,'(''     MD time /ps        :'',f8.1)') mddats(vz)%length_ps
    write (stdout,'(''     dt /fs             :'',f8.1)') mddats(vz)%tstep
    write (stdout,'(''     dumpstep(trj) /fs  :'',f8.1)') mddats(vz)%dumpstep
    if(mddats(vz)%simtype == type_mtd)then
      if(mddats(vz)%cvtype(1) == cv_rmsd)then
        write (stdout,'(''     dumpstep(Vbias) /ps:'',f8.2)') &
        & mddats(vz)%mtd(1)%cvdump_fs/1000.0_wp
        write (stdout,'(''     Vbias factor k /Eh :'',f8.4)') &
        &  mddats(vz)%mtd(1)%kpush
        write (stdout,'(''     Vbias exp α /bohr⁻²:'',f8.4)') &
        &  mddats(vz)%mtd(1)%alpha
      endif
    endif
    !$omp end critical

    !>--- the acutal MD call
    call dynamics(moltmps(job),mddats(vz),calculations(job),pr,io)

    !$omp critical
    if (io == 0) then
      write (stdout,'(a,i0,a)') '*MD ',vz,' completed successfully'
    else
      write (stdout,'(a,i0,a)') '*MD ',vz,' terminated early'
    end if
    !deallocate(moltmp%at,moltmp%xyz)
    !$omp end critical
    !$omp end task
  end do
  !$omp taskwait
  !$omp end single
  !$omp end parallel

  !>--- collect trajectories into one
  call collect(nsim,mddats)

  deallocate (calculations)
  if(allocated(moltmps)) deallocate(moltmps)
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
    open(newunit=ich,file='crest_dynamics.trj')
    do i=1,n
      atmp = mddats(i)%trajectoryfile
      inquire(file=atmp,exist=ex)
      if(ex)then
        open(newunit=ich2,file=atmp)
        io=0
        do while( io == 0 )
          read(ich2,'(a)',iostat=io) btmp
          if(io == 0)then
          write(ich,'(a)') trim(btmp)
          endif
        enddo
        close(ich2) 
      endif 
    enddo
    close(ich)
    return
end subroutine collect
end subroutine crest_search_multimd
!========================================================================================!
subroutine crest_search_multimd_init(env,mol,mddat,nsim)
  use crest_parameters,only:wp,stdout
  use crest_data
  use strucrd
  use calc_type
  use calc_module
  use dynamics_module
  use shake_module
  use iomod,only:makedir,directory_exist,remove
  use omp_lib
  implicit none
  type(systemdata),intent(inout) :: env
  type(mddata) :: mddat
  type(coord) :: mol
  integer,intent(inout) :: nsim
  integer :: i,io
  logical :: pr
!========================================================================================!
  type(calcdata) :: calc
  type(shakedata) :: shk

  real(wp) :: energy
  real(wp),allocatable :: grad(:,:)
  character(len=*),parameter :: mdir = 'MDFILES' 
 
!========================================================================================!

  !>--- check if we have any MD & calculation settings allocated
  mddat = env%mddat
  calc = env%calc
  if (.not. mddat%requested) then
    write (stdout,*) 'MD requested, but no MD settings present.'
    return
  else if (calc%ncalculations < 1) then
    write (stdout,*) 'MD requested, but no calculation settings present.'
    return
  end if

  !>--- init SHAKE?
  if (mddat%shake) then
    calc%calcs(1)%rdwbo = .true.
    allocate (grad(3,mol%nat),source=0.0_wp)
    call engrad(mol,calc,energy,grad,io)
    deallocate (grad)
    calc%calcs(1)%rdwbo = .false.

    shk%shake_mode = env%mddat%shk%shake_mode
    call move_alloc(calc%calcs(1)%wbo,shk%wbo)

    mddat%shk = shk
    call init_shake(mol%nat,mol%at,mol%xyz,mddat%shk,pr)
    mddat%nshake = mddat%shk%ncons
  end if
  !>--- complete real-time settings to steps
  call mdautoset(mddat,io)

  !>--- (optional)  MTD initialization
  if( nsim < 0 )then
    mddat%simtype = type_mtd  !>-- set runtype to MTD

    call defaultGF(env)
    write(stdout,*)'list of applied metadynamics Vbias parameters:'
    do i=1,env%nmetadyn
       write(stdout,'(''$metadyn '',f10.5,f8.3,i5)') env%metadfac(i)/env%rednat,env%metadexp(i)
    enddo
    write(stdout,*)

    !>--- how many simulations
    nsim = env%nmetadyn 
  endif

  return
end subroutine crest_search_multimd_init
subroutine crest_search_multimd_init2(env,mddats,nsim)
  use crest_parameters, only: wp,stdout
  use crest_data
  use strucrd
  use calc_type
  use calc_module
  use dynamics_module
  use metadynamics_module
  use shake_module
  use iomod,only:makedir,directory_exist,remove
  use omp_lib
  implicit none
  type(systemdata),intent(inout) :: env
  type(mddata) :: mddats(nsim)
  integer :: nsim
  integer :: i,io
  logical :: ex
!========================================================================================!
  type(mtdpot),allocatable :: mtds(:)
  
  character(len=80) :: atmp
  character(len=*),parameter :: mdir = 'MDFILES' 
 
  !>--- parallel MD setup
  ex = directory_exist(mdir)
  if (ex) then
    call rmrf(mdir)
  endif
  io = makedir(mdir)
  do i = 1,nsim
    mddats(i)%md_index = i
    write (atmp,'(a,i0,a)') 'crest_',i,'.trj'
    mddats(i)%trajectoryfile = mdir//sep//trim(atmp)
    write (atmp,'(a,i0,a)') 'crest_',i,'.mdrestart'
    mddats(i)%restartfile = mdir//sep//trim(atmp)
  end do

  allocate(mtds(nsim))
  do i=1,nsim
   if(mddats(i)%simtype == type_mtd)then
     mtds(i)%kpush = env%metadfac(i)/env%rednat
     mtds(i)%alpha = env%metadexp(i)
     mtds(i)%cvdump_fs = float(env%mddump)
     mtds(i)%mtdtype = cv_rmsd 

     mddats(i)%npot = 1
     allocate(mddats(i)%mtd(1), source=mtds(i))
     allocate(mddats(i)%cvtype(1), source=cv_rmsd)
   endif 
  enddo
  if(allocated(mtds))deallocate(mtds) 

  return
end subroutine crest_search_multimd_init2
!========================================================================================!
subroutine crest_search_multimd2(env,mols,mddats,nsim)
!*******************************************************************
!* this runs nsim MDs on nsim selected different structures (mols)
!*******************************************************************
  use crest_parameters,only:wp,stdout
  use crest_data
  use strucrd
  use calc_type
  use calc_module
  use dynamics_module
  use shake_module
  use iomod,only:makedir,directory_exist,remove
  use omp_lib
  implicit none
  type(systemdata),intent(inout) :: env
  type(mddata) :: mddats(nsim)
  integer :: nsim
  type(coord) :: mols(nsim)
  type(coord),allocatable :: moltmps(:)
  integer :: i,j,io,ich
  logical :: pr,ex
  type(calcdata) :: calc
  type(mddata) :: mddat
  integer :: T
  real(wp) :: percent
  character(len=52) :: bar
  character(len=80) :: atmp
  character(len=*),parameter :: mdir = 'MDFILES' 
 
  type(calcdata),allocatable :: calculations(:)
  integer :: vz,job,thread_id
!===========================================================!
  !>--- set threads
  if (env%autothreads) then
    call ompautoset(env%threads,7,env%omp,env%MAXRUN,nsim)
  end if

  !>--- check if we have any MD & calculation settings allocated
  mddat = env%mddat
  calc = env%calc
  if (.not. mddat%requested) then
    write (stdout,*) 'MD requested, but no MD settings present.'
    return
  else if (calc%ncalculations < 1) then
    write (stdout,*) 'MD requested, but no calculation settings present.'
    return
  end if

  !>--- prepare calculation objects for parallelization (one per thread)
  T = env%threads
  allocate (calculations( T ),source=env%calc)
  allocate (moltmps( T ), source=mols(1) )
  do i = 1,T
    do j = 1,calc%ncalculations
      calculations(i)%calcs(j) = env%calc%calcs(j)
      !>--- directories
      ex = directory_exist(env%calc%calcs(j)%calcspace)
      if (.not. ex) then
        io = makedir(trim(env%calc%calcs(j)%calcspace))
      end if
      write (atmp,'(a,"_",i0)') sep,i
      calculations(i)%calcs(j)%calcspace = env%calc%calcs(j)%calcspace//trim(atmp)
    end do
    calculations(i)%pr_energies = .false.
  end do

  !>--- other settings
  pr = .false.

  !>--- run the MDs
  !$omp parallel &
  !$omp shared(env,calculations,mddats,mols,pr,percent,bar,ich, moltmps)
  !$omp single
  do i = 1,nsim

    call initsignal()
    vz = i
    !$omp task firstprivate( vz ) private( job,thread_id,io,ex )
    call initsignal()

    thread_id = OMP_GET_THREAD_NUM()
    job = thread_id + 1
    !$omp critical
    moltmps(job)%nat = mols(vz)%nat
    moltmps(job)%at  = mols(vz)%at
    moltmps(job)%xyz = mols(vz)%xyz

    if(mddats(vz)%simtype == type_md)then
    write (stdout,'(a,i4,a)') 'Starting MD',vz,' with the settings:'
    else if(mddats(vz)%simtype == type_mtd)then
    write (stdout,'(a,i4,a)') 'Starting MTD',vz,' with the settings:'
    endif 
    write (stdout,'(''     MD time /ps        :'',f8.1)') mddats(vz)%length_ps
    write (stdout,'(''     target T /K        :'',f8.1)') mddats(vz)%tsoll
    write (stdout,'(''     dt /fs             :'',f8.1)') mddats(vz)%tstep
    write (stdout,'(''     dumpstep(trj) /fs  :'',f8.1)') mddats(vz)%dumpstep
    if(mddats(vz)%simtype == type_mtd)then
      if(mddats(vz)%cvtype(1) == cv_rmsd)then
        write (stdout,'(''     dumpstep(Vbias) /ps:'',f8.2)') &
        & mddats(vz)%mtd(1)%cvdump_fs/1000.0_wp
        write (stdout,'(''     Vbias factor k /Eh :'',f8.4)') &
        &  mddats(vz)%mtd(1)%kpush
        write (stdout,'(''     Vbias exp α /bohr⁻²:'',f8.4)') &
        &  mddats(vz)%mtd(1)%alpha
      endif
    endif
    !$omp end critical

    !>--- the acutal MD call
    call dynamics(moltmps(job),mddats(vz),calculations(job),pr,io)

    !$omp critical
    if (io == 0) then
      write (stdout,'(a,i0,a)') '*MD ',vz,' completed successfully'
    else
      write (stdout,'(a,i0,a)') '*MD ',vz,' terminated with early'
    end if
    !deallocate(moltmp%at,moltmp%xyz)
    !$omp end critical
    !$omp end task
  end do
  !$omp taskwait
  !$omp end single
  !$omp end parallel

  !>--- collect trajectories into one
  call collect(nsim,mddats)

  deallocate (calculations)
  if(allocated(moltmps)) deallocate(moltmps)
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
    open(newunit=ich,file='crest_dynamics.trj')
    do i=1,n
      atmp = mddats(i)%trajectoryfile
      inquire(file=atmp,exist=ex)
      if(ex)then
        open(newunit=ich2,file=atmp)
        io=0
        do while( io == 0 )
          read(ich2,'(a)',iostat=io) btmp
          if(io == 0)then
          write(ich,'(a)') trim(btmp)
          endif
        enddo
        close(ich2) 
      endif 
    enddo
    close(ich)
    return
end subroutine collect
end subroutine crest_search_multimd2
!========================================================================================!
