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

subroutine crest_optimization(env,tim)
!***********************************************
!* subroutine crest_optimization
!* This routine implements a standalone runtype
!* to perform geometry optimization for the
!* specified input file (read from env%ref).
!* An optional inline refinement step is run
!* after geometry optimization when refine_queue
!* is allocated (e.g. set by --A//B or --refine).
!***********************************************
  use crest_parameters,only:wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
  use iomod,only:colorify
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,T,Tn
  logical :: pr,wr
!========================================================================================!
  type(calcdata) :: calc

  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)

  character(len=80) :: atmp
  character(len=*),parameter :: partial = '∂E/∂'

! ══════════════════════════════════════════════════════════════════════════════
  write (stdout,*)
  write (stdout,*) " -------------------------------------------- "
  write (stdout,*) "  ####  #####  ##### # #    # # ###### ###### "
  write (stdout,*) " #    # #    #   #   # ##  ## #     #  #      "
  write (stdout,*) " #    # #    #   #   # # ## # #    #   #####  "
  write (stdout,*) " #    # #####    #   # #    # #   #    #      "
  write (stdout,*) " #    # #        #   # #    # #  #     #      "
  write (stdout,*) "  ####  #        #   # #    # # ###### ###### "
  write (stdout,*) " -------------------------------------------- "
  write (stdout,*)

!========================================================================================!
  call new_ompautoset(env,'max',0,T,Tn)
  call ompprint_intern()
  call tim%start(14,'Geometry optimization')
!========================================================================================!
  call env%ref%to(mol)
  write (stdout,*)
  call smallhead('Input structure:')
  call mol%append(stdout)
  write (stdout,*)

!========================================================================================!

  allocate (grad(3,mol%nat),source=0.0_wp)
  call calc%copy(env%calc)
!>--- check if we have any calculation settings allocated
  if (calc%ncalculations < 1) then
    write (stdout,*) 'no calculations allocated'
    return
  else
    call calc%info(stdout)
  end if
  write (stdout,'(a)') repeat('-',80)

!>-- geometry optimization
  pr = .true. !> stdout printout
  wr = .true. !> write crestopt.log.xyz
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
  call optimize_geometry(mol,molnew,calc,energy,grad,pr,wr,io)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
  if (io == 0) then
!>--- print out the results
    write (stdout,*) 'SUCCESS!'
    write (stdout,*) 'geometry successfully optimized!'
    write (stdout,*)
    write (stdout,'(a)') repeat('-',80)

    write (stdout,*)
    write (stdout,'(1x,a)') 'FINAL CALCULATION SUMMARY'
    write (stdout,'(1x,a)') '========================='
    call calculation_summary(calc,mol,energy,grad,molnew)

    write (stdout,*)
    gnorm = norm2(grad)

! ── optional inline refinement (e.g. --A//B or --refine) ─────────
    if (allocated(env%refine_queue)) then
      do i = 1,size(env%refine_queue,1)
        select case (env%refine_queue(i))
        case (refine%singlepoint)
          write (stdout,'(a)') '> Running SP refinement ...'
          calc%refine_stage = refine%singlepoint
          call engrad(molnew,calc,energy,grad,io)
          calc%refine_stage = 0
          gnorm = norm2(grad)
          if (io == 0) write (stdout,'(1x,a,f20.10,a)') 'Refined energy: ',energy,' Eh'
        case (refine%geoopt)
          write (stdout,'(a)') '> Re-optimizing at higher level ...'
          mol = molnew
          calc%refine_stage = refine%geoopt
          call optimize_geometry(mol,molnew,calc,energy,grad,pr,wr,io)
          calc%refine_stage = 0
          gnorm = norm2(grad)
        end select
      end do
    end if

    write (stdout,'(a)') '> Optimized geometry written to crestopt.xyz'
    write (atmp,'(1x,"Etot=",f16.10,1x,"g norm=",f12.8)') energy,gnorm
    molnew%comment = trim(atmp)

    open (newunit=ich,file='crestopt.xyz')
    call molnew%append(ich)
    close (ich)
  else
    write (stdout,*) 'geometry optimization FAILED!'
    env%iostatus_meta = status_failed
  end if

  write (stdout,*)
  write (stdout,'(a)') repeat('=',42)
  write (stdout,'(1x,a,f20.10,a)') 'TOTAL ENERGY   ',energy,' Eh'
  write (stdout,'(1x,a,f20.10,a)') 'GRADIENT NORM  ',norm2(grad),' Eh/a0'
  write (stdout,'(a)') repeat('=',42)

  if (io /= 0) then
    write (stdout,*) 'WARNING: geometry optimization FAILED!'
  end if

  deallocate (grad)
!========================================================================================!
  call tim%stop(14)

!========================================================================================!
!>--- append numerical hessian calculation
  if (io == 0.and.env%crest_ohess) then
    call env%ref%load(molnew)      !> load the optimized geometry
    call crest_numhess(env,tim) !> run the numerical hessian
  end if

!========================================================================================!

!========================================================================================!
!>--- append deform opt hessian calculation
  if (io == 0.and.calc%deform_opt_hess) then !.and. calc%do_HR )then
    call env%ref%load(molnew)      !> load the optimized geometry
    call deform_opt_hess(calc,molnew) !> run the hessian reconstruction
  end if

!========================================================================================!

  return
end subroutine crest_optimization

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!
subroutine crest_ensemble_optimization(env,tim)
!***********************************************
!* subroutine crest_ensemble_optimization
!* This routine implements a standalone runtype
!* to perform geometry optimizations along an
!* ensemble or trajectory file.
!***********************************************
  use crest_parameters,only:wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
  use parallel_interface
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,c,T,Tn
  logical :: pr,wr,ex
!========================================================================================!
  type(calcdata) :: calc

  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)

  character(len=:),allocatable :: ensnam
  integer :: nat,nall
  type(coord),allocatable :: structures(:)
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
    write (stdout,*) '**ERROR** no ensemble file provided.'
    env%iostatus_meta = status_input
    return
  end if

!>--- start the timer
  call tim%start(14,'Ensemble optimization')

!>---- read the input ensemble as a list of coord objects (in Bohr)
  call rdensemble(ensnam,nall,structures)
  if (nall .lt. 1) then
    write (stdout,*) '**ERROR** empty ensemble file.'
    env%iostatus_meta = status_input
    return
  end if

!>--- set OMP parallelization
  call new_ompautoset(env,'auto',nall,T,Tn)

!========================================================================================!
!>--- printout header
  write (stdout,*)
  write (stdout,'(10x,"┍",49("━"),"┑")')
  write (stdout,'(10x,"│",14x,a,14x,"│")') "ENSEMBLE OPTIMIZATION"
  write (stdout,'(10x,"┕",49("━"),"┙")')
  write (stdout,*)
  write (stdout,'(1x,a,i0,a,1x,a)') 'Optimizing all ',nall,' structures of file',trim(ensnam)
!>--- call the loop (optimized geometries and energies are written into structures)
  call crest_oloop(env,nall,structures,.false.)

!========================================================================================!
!>--- output
  write (stdout,'(/,a,a,a)') 'Rewriting ',ensemblefile,' in the correct order'// &
  & ' (failed optimizations are assigned an energy of +1.0)'
!>--- write the optimized ensemble in extxyz format
  do i = 1,nall
    structures(i)%wrextxyz = .true.
  end do
  call wrensemble(ensemblefile,nall,structures)
  deallocate (structures)
  write (stdout,'(/,a,a,a)') 'Optimized ensemble written to <',ensemblefile,'>'

!========================================================================================!
!>--- (optional) refinement step
  if (allocated(env%refine_queue)) then
    write (stdout,*)
    call crest_refine(env,ensemblefile,ensemblefile//'.refine')
    write (stdout,'(/,a,a,a)') 'Refined ensemble written to <',ensemblefile,'.refine>'
  end if

!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_ensemble_optimization

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!
subroutine crest_ensemble_screening(env,tim)
!****************************************************
!* subroutine crest_ensemble_screening
!* This routine implements a standalone runtype
!* to perform geometry optimizations along an
!* ensemble in a multilevel step and sort in between
!****************************************************
  use crest_parameters,only:wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
  use iomod
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,c,T,Tn
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
  logical :: multilevel(6)
!========================================================================================!
  write (*,*)
!>--- check for the ensemble file
  inquire (file=env%ensemblename,exist=ex)
  if (ex) then
    ensnam = env%ensemblename
  else
    write (stdout,*) '**ERROR** no ensemble file provided.'
    env%iostatus_meta = status_input
    return
  end if

!>--- start the timer
  call tim%start(14,'Ensemble screening')

!>---- read the input ensemble
  call rdensembleparam(ensnam,nat,nall)
  if (nall .lt. 1) then
    write (stdout,*) '**ERROR** empty ensemble file.'
    env%iostatus_meta = status_input
    return
  end if
!>--- set OMP parallelization
  call new_ompautoset(env,'auto',nall,T,Tn)

!========================================================================================!
!>--- printout header
  write (stdout,*)
  write (stdout,'(10x,"┍",48("━"),"┑")')
  write (stdout,'(10x,"│",15x,a,15x,"│")') "ENSEMBLE SCREENING"
  write (stdout,'(10x,"┕",48("━"),"┙")')
  write (stdout,*)
  write (stdout,'(1x,''Multilevel optimization and structure screening.'')')
  write (stdout,*)
  write (stdout,'(1x,a,a)') 'Input file: ','<'//trim(ensnam)//'>'
  write (stdout,'(1x,a,i0,a)') 'Containing ',nall,' structures.'

!>--- call the loop
  call rmrfw('crest_rotamers_')
  call optlev_to_multilev(3.0d0,multilevel)
  call crest_multilevel_oloop(env,ensnam,multilevel,0)
  if (env%iostatus_meta .ne. 0) return

!>--- printout
  call catdel('cregen.out.tmp')
  write (stdout,'(/,1x,a,1x,a)') 'Final ensemble on file','<'//trim(ensemblefile)//'>'

  call rename(conformerfile,trim(ensemblefile))

!>--- clean up
  call screen_cleanup

!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_ensemble_screening

