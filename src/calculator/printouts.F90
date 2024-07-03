! This file is part of crest.
!
! Copyright (C) 2024 Philipp Pracht
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

module crest_calculator_printout
!>--- types and readers
  use calc_type
  use crest_parameters
  use strucrd
  use gradreader_module
  implicit none
!=========================================================================================!
!>--- private module variables and parameters
  private
  character(len=*),parameter :: partial = '∂E/∂'

  public :: calculation_summary

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine calculation_summary(calc,mol,energy,grad,molnew,iounit,print)
    implicit none
    type(calcdata),intent(inout)    :: calc
    type(coord),intent(in)          :: mol
    real(wp),intent(in),optional    :: energy
    real(wp),intent(in),optional    :: grad(3,mol%nat)
    type(coord),intent(in),optional :: molnew
    integer,intent(in),optional     :: iounit
    logical,intent(in),optional     :: print

    integer :: i,j,k,l
    real(wp) :: dip
    real(wp),allocatable :: cn(:)
    logical :: pr,skiplen
    integer :: iunit

!>--- optional arguments to local settings
    if (present(iounit)) then
      iunit = iounit
    else
      iunit = stdout
    end if

    if (present(print)) then
      pr = print
    else
      pr = .true.
    end if

!>--- we might want to skip some printout for large systems...
    skiplen = mol%nat > 100

!>--- Bond orders
    if (any(calc%calcs(:)%rdwbo)) then
      write (iunit,*)
      write (iunit,*) 'Connectivity information (bond order):'
      do k = 1,calc%ncalculations
        if (calc%calcs(k)%rdwbo) then
          write (iunit,'("> ",a,i3,a:)') 'Calculation level ',k,'; '//calc%calcs(k)%shortflag
          write (iunit,'(a12,a12,a10)') 'Atom A','Atom B','BO(A-B)'
          do i = 1,mol%nat
            do j = 1,i-1
              if (calc%calcs(k)%wbo(i,j) > 0.01_wp) then
                write (iunit,*) i,j,calc%calcs(k)%wbo(i,j)
              end if
            end do
          end do
        end if
      end do
      write (iunit,*)
      write (iunit,'(a)') repeat('-',80)
    end if

!>--- atomic charges
    if(any(calc%calcs(:)%rdqat))then
      if(present(molnew))then
        call mol%get_cn(cn,'cov')
      else
        call mol%get_cn(cn,'cov')
      endif
      write (iunit,*)
      write (iunit,*) 'Atomic charges and covalent CN:'
      do k = 1,calc%ncalculations
        if (calc%calcs(k)%rdqat.and.allocated(calc%calcs(k)%qat)) then
          write (iunit,'("> ",a,i3,a)') 'Calculation level ',k,': '//calc%calcs(k)%shortflag
          write (iunit,'(3x,a6,1x,a6,a12,a12)') 'atom','type','charge','covCN'
          do j=1,mol%nat
            write (iunit,'(3x,i6,1x,a6,2f12.6)') j,i2e(mol%at(j),'nc'), &
            & calc%calcs(k)%qat(j), cn(j)
          enddo  
          write (iunit,*)
          if(calc%calcs(k)%dumpq)then
            call dumpq(k,mol%nat,calc%calcs(k)%qat)
          endif
        end if
      end do
      write (iunit,'(a)') repeat('-',80)
    endif

!>--- Dipole moments
    if (any(calc%calcs(:)%rddip)) then
      write (iunit,*)
      write (iunit,*) 'Molecular dipole moments (a.u.):'
      do k = 1,calc%ncalculations
        if (calc%calcs(k)%rddip) then
          dip = norm2(calc%calcs(k)%dipole)*autod
          write (iunit,'("> ",a,i3,a:)') 'Calculation level ',k,': '//calc%calcs(k)%shortflag
          write (iunit,'(2x,19x,3x,a10,a10,a10,a14)') 'x','y','z','tot (Debye)' 
          write (iunit,'(2x,a,3x,4f10.3)') 'total dipole moment',calc%calcs(k)%dipole(:),dip
          write (iunit,*)
        end if
      end do
      write (iunit,'(a)') repeat('-',80)
    end if

!>--- gradients
    if (all(calc%calcs(:)%rdgrad.eqv..false.)) then
      write (iunit,*)
      write (iunit,'(a)') '> No gradients calculated'
    else if (present(grad)) then
      write (iunit,*)
      write (iunit,'(a)') '> Final molecular gradient ( Eh/a0 ):'
      write (iunit,'(13x,a,13x,a,13x,a)') partial//'x',partial//'y',partial//'z'
      k = mol%nat-5
      do i = 1,mol%nat
        if (.not.skiplen) then
          write (iunit,'(3f18.8)') grad(1:3,i)
        else
          if (i <= 5.or.i >= k) then
            write (iunit,'(3f18.8)') grad(1:3,i)
          elseif (i == 6) then
            write (iunit,'(a18,a18,a18)') '<...>','<...>','<...>'
          end if
        end if
      end do
      write (iunit,'(a,f18.8,a)') '> Gradient norm:',norm2(grad),' Eh/a0'
      if (skiplen) then
        write (iunit,'(a)') "> (printout partially skipped due to system size)"
        if (present(energy)) then
          write (iunit,'(a)') "> (see crest.engrad for full gradient)"
        end if
      end if
    end if

!>--- individual energies and norms
    if (calc%ncalculations > 1) then
      write (iunit,*)
      write (iunit,'(a)') '> Individual energies and gradient norms:'
      write (iunit,'(24x,a18,a18,3x,a)') 'Energy/Eh','||Grad||/Eh/a0','Method'
      do k = 1,calc%ncalculations
        if (calc%calcs(k)%rdgrad) then
          write (iunit,'(2x,a,i3,a,2f18.8,4x,a)') 'Calculation level ',k,':', &
          & calc%etmp(k),norm2(calc%grdtmp(:,:,k)),calc%calcs(k)%shortflag
        else
          write (iunit,'(2x,a,i3,a,f18.8,a18,4x,a)') 'Calculation level ',k,':', &
          & calc%etmp(k),'-',calc%calcs(k)%shortflag
        end if
      end do
      if (calc%nconstraints > 0) then
        write (iunit,'(1x,a)') '(+ constraints contribution)'
      end if
    end if

!>---- write gradient
    if (present(grad).and.present(energy)) then
      write (iunit,*)
      write (iunit,'(a)',advance='no') '> Writing crest.engrad ...'
      flush (iunit)
      call write_engrad('crest.engrad',energy,grad)
      write (iunit,*) 'done.'
    end if

!>--- print output structure (if present)
    if (present(molnew)) then
      write (iunit,*)
      write (iunit,'(a)') repeat('-',80)
      write (iunit,*)
      write (iunit,'(1x,a)') '-----------------'
      write (iunit,'(1x,a)') 'Output structure:'
      write (iunit,'(1x,a)') '-----------------'
      write (iunit,'(i6)') molnew%nat
      write (iunit,*)
      k = molnew%nat-5
      do i = 1,molnew%nat
        if (.not.skiplen) then
          write (iunit,'(a3,3f18.8)') i2e(molnew%at(i),'nc'),molnew%xyz(1:3,i)*autoaa
        else
          if (i <= 5.or.i >= k) then
            write (iunit,'(a2,3f18.8)') i2e(molnew%at(i),'nc'),molnew%xyz(1:3,i)*autoaa
          elseif (i == 6) then
            write (iunit,'(a18,a18,a18)') '<...>','<...>','<...>'
          end if
        end if
      end do
      if (skiplen) then
        write (iunit,'(a)') "> (printout partially skipped due to system size)"
        write (iunit,'(a)') "> (see crestopt.xyz for full structure)"
      end if
    end if

  end subroutine calculation_summary

!========================================================================================!

  subroutine dumpq(id,nat,q)
!********************************
!* write atomic charges to file
!********************************
    implicit none
    integer,intent(in) :: id,nat
    real(wp),intent(in) :: q(nat)
    integer :: i,ich
    character(len=50) :: atmp

    write(atmp,'("charges.",i0)') id
    open(newunit=ich,file=trim(atmp))
    do i=1,nat
       write(ich,'(F20.10)') q(i)
    enddo
    close(ich)
  end subroutine dumpq 

!========================================================================================!
!========================================================================================!
end module crest_calculator_printout
