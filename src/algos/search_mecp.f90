
subroutine crest_search_mecp(env,tim)
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
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,m
  logical :: pr,wr
!===========================================================!
  type(calcdata) :: calc
  type(mddata) :: mddat
  type(shakedata) :: shk

  type(mddata),allocatable :: mddats(:)
  integer :: nsim

  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)
  character(len=:),allocatable :: ensnam
  integer :: nat,nall
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  logical :: dump,ex
  character(len=80) :: atmp

!===========================================================!
!>--- printout header
  write (stdout,*)
  write (stdout,'(10x,"┍",49("━"),"┑")')
  write (stdout,'(10x,"│",10x,a,10x,"│")') "CREST MECP SAMPLING ALGORITHM"
  write (stdout,'(10x,"┕",49("━"),"┙")')
  write (stdout,*)

!===========================================================!
!>--- setup
  call env%ref%to(mol)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)

!>--- check calculation setup
  ex = env%calc%ncalculations > 1
  if(.not.ex)then
    write(stdout,*) 'not enough calculation levels specified for MECP search.'
    error stop
  endif
  call print_gapcons(env%calc)

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
!>--- rename ensemble and sort
  call rename(ensemblefile,mecpensemble)
  call newcregen(env,12)
 

!==========================================================!
  return
contains
subroutine print_gapcons(calc)
  use calc_type
  use constraints
  implicit none

  type(calcdata) :: calc
  integer :: i,t
  logical :: ex 

  if(calc%nconstraints < 1)then
     write(stdout,*) 'no gap constraint provided'
  endif

  ex = .false.
  do i=1,calc%nconstraints
     t = calc%cons(i)%type
     select case( t )
     case( -1 )
      ex = .true.
      write (stdout,'(1x,a1x,"[",a,"]")') 'nonadiabatic gap constraint','σ*ΔE²/(ΔE+α)'
      write (stdout,'(" σ=",f8.5," α=",f8.5)') calc%cons(i)%fc(1:2)
     case( -2 )
      ex = .true.
      write (stdout,'(1x,a)') 'nonadiabatic gap constraint'
      write (stdout,'(1x,"Vgap = [",a,"] with")') &
      &     'σ*(exp(-β|ΔE|)+C) * ΔE²/(|ΔE|+α)'
      write (stdout,'(" σ = ",f8.5,/," α = ",f8.5,/," C = ",f8.5,/," β = ",f8.5)') &
      & calc%cons(i)%fc(1:3),27.2114_wp
     case default
       continue
     end select
  enddo

  if(.not.ex)then
    write(stdout,*) 'no gap constraint provided'
    error stop
  else
    write(stdout,*) 
  endif

end subroutine print_gapcons

end subroutine crest_search_mecp


