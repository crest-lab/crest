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
!* specified input file (read from env%ref)
!***********************************************
  use crest_parameters,only:wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
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
  character(len=*),parameter :: partial = '∂E/∂'
!========================================================================================!
  call tim%start(14,'Geometry optimization')
  call env%ref%to(mol)
  write (stdout,*)
  call smallhead('Input structure:')
  call mol%append(stdout)
  write (stdout,*)

!========================================================================================!

  allocate (grad(3,mol%nat),source=0.0_wp)
  calc = env%calc
  !>--- check if we have any calculation settings allocated
  if (calc%ncalculations < 1) then
    write (stdout,*) 'no calculations allocated'
    return
  else
    call calc%info( stdout )
  end if
  write(stdout,'(a)') repeat('-',80)

  !>--- first energy&gradient calculation
  call engrad(mol,calc,energy,grad,io)

  !>-- geopetry optimization
  pr = .true. !> stdout printout
  wr = .true. !> write crestopt.log
  call optimize_geometry(mol,molnew,calc,energy,grad,pr,wr,io)

  if (io == 0) then
    write (stdout,*) 'geometry successfully optimized!'
    write (stdout,*)
    write(stdout,'(a)') repeat('-',80)
    write (stdout,*)
    call smallhead( 'Output structure:') 
    call molnew%append(stdout)
    write (stdout,*)
    write (stdout,*) 'optimized geometry written to crestopt.xyz'
    gnorm = norm2(grad)
    write (atmp,'(1x,"Etot=",f16.10,1x,"g norm=",f12.8)') energy,gnorm
    molnew%comment = trim(atmp)
    open (newunit=ich,file='crestopt.xyz')
    call molnew%append(ich)
    close (ich)
  else
    write (stdout,*) 'geometry optimization FAILED!'
  end if

!========================================================================================!
!>--- print out the results
   if(any(calc%calcs(:)%rdwbo))then
   write(stdout,*)
   write(stdout,*) 'Connectivity information (bond order):'
   do k=1,calc%ncalculations
     if(calc%calcs(k)%rdwbo)then
       write(stdout,'("> ",a,i0)') 'Calculation level ',k
       write(stdout,'(a12,a12,a10)') 'Atom A','Atom B','BO(A-B)'
       do i=1,mol%nat
         do j=1,i-1
           if(calc%calcs(k)%wbo(i,j) > 0.0002_wp)then
            write(stdout,*) i,j,calc%calcs(k)%wbo(i,j)
           endif
         enddo
       enddo
     endif
   enddo
   endif
   write(stdout,*)

   write(stdout,'(a)') repeat('-',80)
   write(stdout,'(a)') '> Final molecular gradient ( Eh/a0 ):'
   write(stdout,'(13x,a,13x,a,13x,a)')partial//'x',partial//'y',partial//'z'
   do i = 1,mol%nat
      write (stdout,'(3f18.8)') grad(1:3,i)
   end do
   write(stdout,'(a,f18.8,a)') '> Gradient norm:',norm2(grad),' Eh/α'

   if(calc%ncalculations > 1)then
   write(stdout,*)
   write(stdout,'(a)') '> Individual energies and gradient norms:'
     do k=1,calc%ncalculations
       write(stdout,'(1x,a,i0,2f18.8)') 'calculation ',k,calc%etmp(k),norm2(calc%grdtmp(:,:,k))
     enddo
     if(calc%nconstraints > 0)then
       write(stdout,'(1x,a)') '(+ constraints contribution)'
     endif
   endif

   write(stdout,*)
   write(stdout,'(a)') repeat('=',40)
   write(stdout,'(1x,a,f20.10,a)') 'TOTAL ENERGY ',energy,' Eh'
   write(stdout,'(1x,a,f20.10,a)') 'GRADIENT NORM',norm2(grad),' Eh/α'
   write(stdout,'(a)') repeat('=',40)
   
   if(io /= 0)then
    write (stdout,*) 'WARNING: geometry optimization FAILED!'
   endif

  deallocate (grad)
!========================================================================================!
  call tim%stop(14)
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
  call tim%start(14,'Ensemble optimization')

!>---- read the input ensemble
  call rdensembleparam(ensnam,nat,nall)
  if (nall .lt. 1) return
  allocate (xyz(3,nat,nall),at(nat),eread(nall))
  call rdensemble(ensnam,nat,nall,at,xyz,eread)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: crest_oloop requires coordinates in Bohrs
  xyz = xyz / bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

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
  write (stdout,'(1x,a,i0,a,1x,a)') 'Optimizing all ',nall,' structures of file',trim(ensnam)
  !>--- call the loop
  call crest_oloop(env,nat,nall,at,xyz,eread,.true.)

  deallocate (eread,at,xyz)
!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_ensemble_optimization

