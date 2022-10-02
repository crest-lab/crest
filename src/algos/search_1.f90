
subroutine crest_search_1(env,tim)
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use crest_data
  use strucrd
  use calc_type
  use calc_module
  use dynamics_module
  use shake_module
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol
!===========================================================!
  type(mddata) :: mddat

  type(mddata),allocatable :: mddats(:)
  integer :: nsim

  character(len=:),allocatable :: ensnam
  integer :: nat,nall
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  logical :: dump

!===========================================================!
!>--- printout header
  write (stdout,*)
  write (stdout,'(10x,"┍",49("━"),"┑")')
  write (stdout,'(10x,"│",13x,a,12x,"│")') "CREST SAMPLING ALGORITHM"
  write (stdout,'(10x,"┕",49("━"),"┙")')
  write (stdout,*)

!===========================================================!
!>--- setup
  call env%ref%to(mol)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)

!===========================================================!
!>--- Dynamics

  write(stdout,*)
  write(stdout,'(1x,a)') '------------------------------'
  write(stdout,'(1x,a)') 'Molecular Dynamics Simulations'
  write(stdout,'(1x,a)') '------------------------------'

  call crest_search_multimd_init(env,mol,mddat,nsim)
  allocate (mddats(nsim), source=mddat)
  call crest_search_multimd_init2(env,mddats,nsim)

  call tim%start(2,'MD simulations')
  call crest_search_multimd(env,mol,mddats,nsim)
  call tim%stop(2)
  !>--- a file called crest_dynamics.trj should have been written
  ensnam = 'crest_dynamics.trj'

!==========================================================!
!>--- Reoptimization of trajectories

  write(stdout,*)
  write(stdout,'(1x,a)') '---------------------'
  write(stdout,'(1x,a)') 'Ensemble Optimization'
  write(stdout,'(1x,a)') '---------------------'

  !>--- read ensemble
  call rdensembleparam(ensnam,nat,nall)
  if (nall .lt. 1) then
    write(stdout,*) 'empty ensemble file'
    return
  endif
  allocate (xyz(3,nat,nall),at(nat),eread(nall))
  call rdensemble(ensnam,nat,nall,at,xyz,eread)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: crest_oloop requires coordinates in Bohrs
  xyz = xyz / bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
 
  write(stdout,'(1x,a,i0,a,a,a)')'Optimizing all ',nall, &
  & ' structures from file "',trim(ensnam),'" ...'
  !>--- set threads
  if (env%autothreads) then
    call ompautoset(env%threads,7,env%omp,env%MAXRUN,nall)
  end if
  !>--- optimize
  call tim%start(3,'geom. optimization')
  dump = .true.
  call crest_oloop(env,nat,nall,at,xyz,eread,dump)
  call tim%stop(3)

 

!==========================================================!
  return
end subroutine crest_search_1

!========================================================================================!
subroutine crest_search_multimd(env,mol,mddats,nsim)
  use iso_fortran_env,only:wp => real64,stdout => output_unit
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
  type(coord) :: mol,moltmp
  integer :: i,j,io,ich
  logical :: pr,ex
!========================================================================================!
  type(calcdata) :: calc
  type(mddata) :: mddat

  real(wp) :: percent
  character(len=52) :: bar
  character(len=80) :: atmp
  character(len=*),parameter :: mdir = 'MDFILES' 
 
  type(calcdata),allocatable :: calculations(:)
 ! type(mddata),allocatable :: mddats(:)
  integer :: vz,job,thread_id
!========================================================================================!


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
  allocate (calculations(env%threads),source=env%calc)
  do i = 1,env%threads
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
  !$omp shared(env,calculations,mddats,mol,pr,percent,bar,ich)
  !$omp single
  do i = 1,nsim

    call initsignal()
    vz = i
    !$omp task firstprivate( vz ) private( job,moltmp,thread_id,io,ex )
    call initsignal()

    thread_id = OMP_GET_THREAD_NUM()
    job = thread_id + 1
    !$omp critical
    allocate(moltmp%at(mol%nat),moltmp%xyz(3,mol%nat))
    moltmp%nat = mol%nat
    moltmp%at = mol%at
    moltmp%xyz = mol%xyz

    write (stdout,'(a,i4,a)') 'Starting MD',vz,' with the settings:'
    write (stdout,'(''     MD time /ps        :'',f8.1)') mddats(vz)%length_ps
    write (stdout,'(''     dt /fs             :'',f8.1)') mddats(vz)%tstep
    write (stdout,'(''     dumpstep(trj) /fs  :'',f8.1)') mddats(vz)%dumpstep
    !$omp end critical

    !>--- the acutal MD call
    call dynamics(moltmp,mddats(vz),calculations(job),pr,io)

    !$omp critical
    if (io == 0) then
      write (stdout,'(a,i0,a)') '*MD ',vz,' completed successfully'
    else
      write (stdout,'(a,i0,a)') '*MD ',vz,' terminated with error'
    end if
    deallocate(moltmp%at,moltmp%xyz)
    !$omp end critical
    !$omp end task
  end do
  !$omp taskwait
  !$omp end single
  !$omp end parallel

  !>--- collect trajectories into one
  call collect(nsim,mddats)

  deallocate (calculations)
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
  use iso_fortran_env,only:wp => real64,stdout => output_unit
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
  integer,intent(out) :: nsim
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

    shk%shake_mode = 2
    call move_alloc(calc%calcs(1)%wbo,shk%wbo)

    mddat%shk = shk
    call init_shake(mol%nat,mol%at,mol%xyz,mddat%shk,pr)
    mddat%nshake = mddat%shk%ncons
  end if
  !>--- complete real-time settings to steps
  call mdautoset(mddat,io)

  !>--- MTD initial setup
  call defaultGF(env)
  write(stdout,*)'list of applied metadynamics Vbias parameters:'
  do i=1,env%nmetadyn
     write(stdout,'(''$metadyn '',f10.5,f8.3,i5)') env%metadfac(i)/env%rednat,env%metadexp(i)
  enddo
  write(stdout,*)

  !>--- how many simulations
  nsim = env%nmetadyn 

  return
end subroutine crest_search_multimd_init
subroutine crest_search_multimd_init2(env,mddats,nsim)
  use iso_fortran_env,only:wp => real64,stdout => output_unit
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
  !allocate (mddats(nsim))!,source=mddat)
  do i = 1,nsim
    mddats(i)%md_index = i
    write (atmp,'(a,i0,a)') 'crest_',i,'.trj'
    mddats(i)%trajectoryfile = mdir//sep//trim(atmp)
    write (atmp,'(a,i0,a)') 'crest_',i,'.mdrestart'
    mddats(i)%restartfile = mdir//sep//trim(atmp)
  end do

  allocate(mtds(nsim))
  do i=1,nsim
  mddats(i)%simtype = type_mtd 
  mtds(i)%kpush = env%metadfac(i)/env%rednat
  mtds(i)%alpha = env%metadexp(i)
  mtds(i)%cvdump_fs = float(env%mddump)
  mtds(i)%mtdtype = cv_rmsd 

  mddats(i)%npot = 1
  allocate(mddats(i)%mtd(1), source=mtds(i))
  allocate(mddats(i)%cvtype(1), source=cv_rmsd)
  enddo

  return
end subroutine crest_search_multimd_init2

