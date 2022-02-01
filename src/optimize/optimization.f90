
subroutine crest_optimization(env,tim)
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use crest_data
  use strucrd
  use calc_type
  use calc_module
  use optimize_module
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich
  logical :: pr,wr
!========================================================================================!
  type(calcdata) :: calc

  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)

  character(len=80) :: atmp
!========================================================================================!
  call tim%start(14,'geometry optimization')
  call env%ref%to(mol)
  write (stdout,*)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)
!========================================================================================!

  allocate (grad(3,mol%nat),source=0.0_wp)
  calc = env%calc
  !>--- check if we have any calculation settings allocated
  if (calc%ncalculations < 1) then
    write (stdout,*) 'no calculations allocated'
    return
  end if

  !>--- first energy&gradient calculation
  call engrad(mol,calc,energy,grad,io)

  !>-- geopetry optimization
  pr = .true. !> stdout printout
  wr = .true. !> write crestopt.log
  call optimize_geometry(mol,molnew,calc,energy,grad,pr,wr,io)

  if (io == 0) then
    write (stdout,*) 'geometry successfully optimized!'
    write (stdout,*)
    call molnew%append(stdout)
    write (stdout,*)
    write (stdout,*) 'optimized geometry written to crestopt.xyz'
    gnorm = norm2(grad)
    write (atmp,'(1x,"Etot=",f16.10,1x,"g norm=",f12.8)') energy,gnorm
    molnew%comment = trim(atmp)
    open (newunit=ich,file='crestopt.xyz')
    call molnew%append(ich)
    close (ich)
  end if

  deallocate (grad)
!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_optimization

!========================================================================================!
!> subroutine crest_ensemble_optimization
!> Read in an ensemble file and optimize all structures
!>
!>------------------------------------------------------
subroutine crest_ensemble_optimization(env,tim)
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use crest_data
  use strucrd
  use calc_type
  use calc_module
  use optimize_module
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,c
  logical :: pr,wr,ex
!========================================================================================!
  type(calcdata) :: calc

  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)

  character(len=:),allocatable :: ensnam
  integer :: nat,nall
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  character(len=80) :: atmp
  real(wp) :: percent
  character(len=52) :: bar
!========================================================================================!
  write (*,*)
!>--- check for the ensemble file
  inquire (file=env%ensemblename,exist=ex)
  if (ex) then
    ensnam = env%ensemblename
  else
    write (stdout,*) 'no ensemble file provided.'
    return
  end if

!>--- start the timer
  call tim%start(14,'test implementation')

!>---- read the input ensemble
  call rdensembleparam(ensnam,nat,nall)
  if (nall .lt. 1) return
  allocate (xyz(3,nat,nall),at(nat),eread(nall))
  call rdensemble(ensnam,nat,nall,at,xyz,eread)
!>--- Important: crest_oloop requires coordinates in Bohrs
  xyz = xyz/bohr 

!>--- set OMP parallelization
  if (env%autothreads) then
    !>--- usually, one thread per xtb job
    call ompautoset(env%threads,7,env%omp,env%MAXRUN,nall)
  end if

!========================================================================================!
  !>--- printout header
  write (stdout,*)
  write (stdout,'(10x,"┍",49("━"),"┑")')
  write (stdout,'(10x,"│",14x,a,14x,"│")') "ENSEMBLE OPTIMIZATION"
  write (stdout,'(10x,"┕",49("━"),"┙")')
  write (stdout,*)

  !>--- call the loop
  call crest_oloop(env,nat,nall,at,xyz,eread,.true.)

  deallocate (eread,at,xyz)
!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_ensemble_optimization

!========================================================================================!

subroutine crest_oloop(env,nat,nall,at,xyz,eread,dump)
  use iso_fortran_env,only:wp => real64,stdout => output_unit
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

  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,c,job_id
  logical :: pr,wr,ex
!========================================================================================!
  type(calcdata),allocatable :: calculations(:)
  type(calcdata) :: calc

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
  allocate(calculations(env%threads),source=env%calc)
  do i=1,env%threads
   do j = 1,calc%ncalculations
      calculations(i)%calcs(j) = env%calc%calcs(j)

      !>--- directories
      ex = directory_exist(env%calc%calcs(j)%calcspace)
      if (.not. ex) then
        io = makedir(trim(env%calc%calcs(j)%calcspace))
      end if
      write (atmp,'(a,"_",i0)') sep,i
      !write (atmp,'("_",i0)') i
      calculations(i)%calcs(j)%calcspace = env%calc%calcs(j)%calcspace//trim(atmp)
      !write(*,*) calculations(i)%calcs(j)%calcspace
    end do
    calculations(i)%pr_energies = .false.
  enddo

  !>--- shared variables
  allocate(mol%at(nat),mol%xyz(3,nat))
  allocate(molnew%at(nat),molnew%xyz(3,nat))
  allocate (grad(3,nat),source=0.0_wp)
  pr = .false. !> stdout printout
  wr = .false. !> write crestopt.log

  if (dump) open (newunit=ich,file=ensemblefile)
  c = 0
  k = 0
  !>--- loop over ensemble
  !$omp parallel &
  !$omp shared(env,calculations,nat,nall,at,xyz,c,k,pr,wr,dump,percent,bar,ich) 
  !$omp single
  do i = 1,nall

    call initsignal()
    vz=i
    !$omp task firstprivate( vz ) private(j,job,mol,molnew,calc,energy,grad,io,atmp,gnorm,thread_id) 
    call initsignal()

    thread_id = OMP_GET_THREAD_NUM()
    job = thread_id + 1
    !>--- modify calculation spaces
    !$omp critical
    mol%nat = nat
    molnew%nat = nat
    mol%at = at
    molnew%at = at
    mol%xyz = xyz(1:3,1:nat,vz)
    molnew%xyz = xyz(1:3,1:nat,vz)
    !$omp end critical

    !>--- first energy&gradient calculation
    call engrad(mol,calculations(job),energy,grad,io)

    !>-- geopetry optimization
    call optimize_geometry(mol,molnew,calculations(job),energy,grad,pr,wr,io)

    !$omp critical
    if (io == 0) then
      c = c + 1
      if (dump) then
        gnorm = norm2(grad)
        write (atmp,'(1x,"Etot=",f16.10,1x,"g norm=",f12.8)') energy,gnorm
        molnew%comment = trim(atmp)
        call molnew%append(ich)
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

  if (dump) close (ich)

  deallocate (grad)
  deallocate(molnew%xyz,molnew%at)
  deallocate(mol%xyz,mol%at)
  deallocate(calculations)
  return
end subroutine crest_oloop
