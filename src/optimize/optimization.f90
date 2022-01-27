
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
!========================================================================================!
  write(*,*)
!>--- check for the ensemble file
  inquire(file=env%ensemblename,exist=ex)
  if(ex)then
     ensnam = env%ensemblename
  else
    write(stdout,*) 'no ensemble file provided.'
    return
  endif

!>--- start the timer
  call tim%start(14,'test implementation')


!>---- read the input ensemble
  call rdensembleparam(ensnam,nat,nall)
  if (nall .lt. 1) return
  allocate (xyz(3,nat,nall),at(nat),eread(nall))
  call rdensemble(ensnam,nat,nall,at,xyz,eread)

!========================================================================================!

  allocate (grad(3,nat),source=0.0_wp)
  calc = env%calc
  !>--- check if we have any calculation settings allocated
  if (calc%ncalculations < 1) then
    write (stdout,*) 'no calculations allocated'
    return
  end if

  !>--- shared variables
  mol%nat = nat
  mol%at = at
  molnew%nat = nat
  molnew%at = at
  pr = .false. !> stdout printout
  wr = .false. !> write crestopt.log

  open (newunit=ich,file=ensemblefile)
  c = 0
  !>--- loop over ensemble
  do i = 1,nall

    mol%xyz = xyz(:,:,i)
    molnew%xyz = xyz(:,:,i)

    !>--- first energy&gradient calculation
    call engrad(mol,calc,energy,grad,io)

    !>-- geopetry optimization
    call optimize_geometry(mol,molnew,calc,energy,grad,pr,wr,io)

    !$omp critical
    if (io == 0) then
      c = c+1
      gnorm = norm2(grad)
      write (atmp,'(1x,"Etot=",f16.10,1x,"g norm=",f12.8)') energy,gnorm
      molnew%comment = trim(atmp)
      call molnew%append(ich)
    end if
    !$omp end critical

  end do

  close (ich)

  deallocate (grad)
  deallocate(eread,at,xyz)
!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_ensemble_optimization

