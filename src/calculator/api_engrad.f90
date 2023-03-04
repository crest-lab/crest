!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2023 Philipp Pracht
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

!====================================================!
! module api_engrad
! a collection of engrad calls for different APIs
! this builds the communication between CRESTs
! calculation settings and the api setups
!====================================================!

module api_engrad

  use iso_fortran_env,only:wp => real64,stdout=>output_unit
  use strucrd
  use calc_type
  use iomod,only:makedir,directory_exist,remove
  !> APIs
  use tblite_api
  use gfn0_api
!=========================================================================================!
  implicit none
  !--- private module variables and parameters
  private
  integer :: i,j,k,l,ich,och,io
  logical :: ex


  public :: tblite_engrad
  public :: gfn0_engrad, gfn0occ_engrad

!=========================================================================================!
!=========================================================================================!
contains    !>--- Module routines start here
!=========================================================================================!
!=========================================================================================!

  subroutine tblite_engrad(mol,calc,energy,grad,iostatus)
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus

    character(len=:),allocatable :: cpath
    logical :: loadnew
    iostatus = 0
    
    !>--- setup system call information
    !$omp critical
    call tblite_init(calc,loadnew)
    !> tblite printout handling
    inquire(unit=calc%ctx%unit, opened=ex)
    if((calc%ctx%unit .ne. stdout) .and. ex)then
      close(calc%ctx%unit)
    endif
    if (allocated(calc%calcspace)) then
      ex = directory_exist(calc%calcspace)
      if (.not. ex) then
        io = makedir(trim(calc%calcspace))
      end if
      cpath = calc%calcspace//sep//'tblite.out'
    else
      cpath = 'tblite.out'
    end if
    !write(*,*) 'here', calc%ctx%unit,trim(cpath),stdout
    open(newunit=calc%ctx%unit, file=cpath, status='replace')
    !write(*,*) 'doe',calc%ctx%unit
    deallocate (cpath)
    !> populate parameters and wavefunction
    if(loadnew)then
      call tblite_setup(mol,calc%chrg,calc%uhf,calc%tblitelvl,calc%etemp, &
      &    calc%ctx,calc%wfn,calc%tbcalc)
      call tblite_addsettings(calc%tbcalc,calc%maxscc,calc%rdwbo,calc%saveint)
    endif
    !$omp end critical
    !>--- do the engrad call
    call initsignal()
    call tblite_singlepoint(mol,calc%chrg,calc%uhf,calc%accuracy, &
    & calc%ctx,calc%wfn,calc%tbcalc,energy,grad,calc%tbres,iostatus)
    if(iostatus /= 0) return

    !>--- postprocessing, getting other data
    !$omp critical
    call tblite_wbos(calc,mol,iostatus)
    !$omp end critical

    return
  contains
    subroutine tblite_init(calc,loadnew)
      implicit none
      type(calculation_settings),intent(inout) :: calc
      logical,intent(out) :: loadnew
      loadnew = .false.
      if(.not.allocated(calc%wfn))then
      allocate(calc%wfn)
      loadnew = .true.
      endif
      if(.not.allocated(calc%tbcalc))then
      allocate(calc%tbcalc)
      loadnew=.true.
      endif
      if( .not.allocated(calc%ctx) )then
      allocate(calc%ctx)
      loadnew = .true.
      endif
      if( .not.allocated(calc%tbres) )then
      allocate(calc%tbres)
      loadnew=.true.
      endif
      if( calc%apiclean ) loadnew = .true.
    end subroutine tblite_init
    subroutine tblite_wbos(calc,mol,iostatus)      
      implicit none
      type(calculation_settings),intent(inout) :: calc
      type(coord),intent(in) :: mol
      integer,intent(out) :: iostatus
      iostatus = 0  
      if(.not.calc%rdwbo) return
      if(allocated(calc%wbo))deallocate(calc%wbo)
      allocate(calc%wbo( mol%nat, mol%nat), source=0.0_wp)
      call tblite_getwbos(calc%tbcalc,calc%wfn,calc%tbres,mol%nat,calc%wbo)
    end subroutine tblite_wbos
  end subroutine tblite_engrad


!========================================================================================!

  subroutine gfn0_engrad(mol,calc,g0calc,energy,grad,iostatus)
!> This is the GFN0 engrad call that uses the standard implementation
    implicit none
    !> INPUT
    type(coord) :: mol
    type(calculation_settings) :: calc
    type(gfn0_data),intent(inout),allocatable  :: g0calc
    !> OUTPUT  
    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus
    !> LOCAL
    type(gfn0_results) :: res
    character(len=:),allocatable :: cpath
    logical :: loadnew
    logical :: pr
    iostatus = 0
    pr = .false.
    !>--- setup system call information
    !$omp critical
    call gfn0_init(calc,g0calc,loadnew)
    !> printout handling
    inquire(unit=calc%prch, opened=ex)
    if((calc%prch .ne. stdout) .and. ex)then
      close(calc%prch)
    endif
    if (allocated(calc%calcspace)) then
      ex = directory_exist(calc%calcspace)
      if (.not. ex) then
        io = makedir(trim(calc%calcspace))
      end if
      cpath = calc%calcspace//sep//'gfn0.out'
    else
      cpath = 'gfn0.out'
    end if
    if(( calc%prch .ne. stdout)) then
    open(newunit=calc%prch, file=cpath)
    pr = .true.
    endif
    deallocate (cpath)
    !> populate parameters and wavefunction
    if(loadnew)then
      call gfn0_setup(mol,calc%chrg,calc%uhf,g0calc)
      call gfn0_init2(mol,calc,g0calc)
    endif
    call gfn0_init3(mol,calc,g0calc) 
    !$omp end critical
    !>--- do the engrad call
    call initsignal()
    call gfn0_sp(mol,calc%chrg,calc%uhf,g0calc,energy,grad,iostatus,res)
    if(iostatus /= 0) return
    if(pr)then
      call gfn0_print(calc%prch,g0calc,res)
    endif

    !>--- postprocessing, getting other data
    !$omp critical
    call gfn0_wbos(calc,mol,iostatus)
    !$omp end critical

    return
  contains
    subroutine gfn0_init(calc,g0calc,loadnew)
      implicit none
      type(calculation_settings),intent(inout) :: calc
      type(gfn0_data),intent(inout),allocatable  :: g0calc
      logical,intent(out) :: loadnew
      loadnew = .false.
      if(.not.allocated(g0calc))then
      allocate(g0calc)
      loadnew = .true.
      endif
      if( calc%apiclean ) loadnew = .true.
    end subroutine gfn0_init
    subroutine gfn0_init2(mol,calc,g0calc)
      implicit none
      type(coord),intent(in) :: mol
      type(calculation_settings),intent(inout) :: calc
      type(gfn0_data),intent(inout)  :: g0calc
      if(allocated(calc%solvent) .and. allocated(calc%solvmodel))then
      call gfn0_addsettings(mol,g0calc,calc%solvent,calc%solvmodel)
      endif
      call gfn0_addsettings(mol,g0calc,loadwbo=calc%rdwbo)
    end subroutine gfn0_init2
    subroutine gfn0_init3(mol,calc,g0calc)
      implicit none
      type(coord),intent(in) :: mol
      type(calculation_settings),intent(inout) :: calc
      type(gfn0_data),intent(inout),allocatable  :: g0calc
      integer :: nel,uhf
      nel = g0calc%wfn%nel
      uhf = calc%uhf
      call g0calc%wfn%refresh_occu(nel, uhf)
      call gfn0_addsettings(mol,g0calc,etemp=calc%etemp)
    end subroutine gfn0_init3
    subroutine gfn0_wbos(calc,mol,iostatus)
      implicit none
      type(calculation_settings),intent(inout) :: calc
      type(coord),intent(in) :: mol
      integer,intent(out) :: iostatus
      iostatus = 0
      if(.not.calc%rdwbo) return
      if(allocated(calc%wbo))deallocate(calc%wbo)
      allocate(calc%wbo( mol%nat, mol%nat), source=0.0_wp)
      call gfn0_getwbos(calc%g0calc,mol%nat,calc%wbo)
    end subroutine gfn0_wbos
  end subroutine gfn0_engrad

!========================================================================================!

  subroutine gfn0occ_engrad(mol,calc,g0calc,energy,grad,iostatus)
!> This is the GFN0 engrad call in which a config can be specified
    implicit none
    !> INPUT
    type(coord) :: mol
    type(calculation_settings) :: calc
    type(gfn0_data),intent(inout),allocatable  :: g0calc
    !> OUTPUT  
    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(:,:)
    integer,intent(out) :: iostatus
    !> LOCAL
    type(gfn0_results) :: res
    character(len=:),allocatable :: cpath
    logical :: loadnew,pr
    iostatus = 0
    pr = .false.
    !>--- setup system call information
    !$omp critical
    call gfn0_init(calc,g0calc,loadnew)
    !> printout handling
    inquire(unit=calc%prch, opened=ex)
    if((calc%prch .ne. stdout) .and. ex)then
      close(calc%prch)
    endif
    if (allocated(calc%calcspace)) then
      ex = directory_exist(calc%calcspace)
      if (.not. ex) then
        io = makedir(trim(calc%calcspace))
      end if
      cpath = calc%calcspace//sep//'gfn0.out'
    else
      cpath = 'gfn0.out'
    end if
    if( calc%prch .ne. stdout) then
    open(newunit=calc%prch, file=cpath)
    pr = .true.
    endif
    deallocate (cpath)
    !> populate parameters and wavefunction
    if(loadnew)then
      call gfn0_setup(mol,calc%chrg,calc%uhf,g0calc)
      call gfn0_init2(mol,calc,g0calc)
    endif
    call gfn0_init3(mol,calc,g0calc)
    !$omp end critical
    !>--- do the engrad call
    call initsignal()
    call gfn0_sp_occ(mol,calc%chrg,calc%uhf,calc%occ,g0calc, &
    &    energy,grad,iostatus,res)
    if(iostatus /= 0) return
    if(pr)then
      call gfn0_print(calc%prch,g0calc,res)
    endif

    !>--- postprocessing, getting other data
    !$omp critical
    call gfn0_wbos(calc,mol,iostatus)
    !$omp end critical

    return
  contains
    subroutine gfn0_init(calc,g0calc,loadnew)
      implicit none
      type(calculation_settings),intent(inout) :: calc
      logical,intent(out) :: loadnew
      type(gfn0_data),intent(inout),allocatable  :: g0calc
      integer :: nel,nao,nlev
      loadnew = .false.
      if(.not.allocated(g0calc))then
      allocate(g0calc)
      loadnew = .true.
      endif
      if( calc%apiclean ) loadnew = .true.
    end subroutine gfn0_init
    subroutine gfn0_init2(mol,calc,g0calc)
      implicit none
      type(coord),intent(in) :: mol
      type(calculation_settings),intent(inout) :: calc
      type(gfn0_data),intent(inout),allocatable  :: g0calc
      integer :: nel,nao,nlev
      if(allocated(calc%solvent) .and. allocated(calc%solvmodel))then
      call gfn0_addsettings(mol,g0calc,calc%solvent,calc%solvmodel)
      endif
      call gfn0_addsettings(mol,g0calc,etemp=calc%etemp,loadwbo=calc%rdwbo)
    end subroutine gfn0_init2
    subroutine gfn0_init3(mol,calc,g0calc)
      implicit none
      type(coord),intent(in) :: mol
      type(calculation_settings),intent(inout) :: calc
      type(gfn0_data),intent(inout),allocatable  :: g0calc
      integer :: nel,nao,nlev
      if(.not.allocated(calc%occ))then
        nel = g0calc%wfn%nel
        nao = g0calc%basis%nao 
        if(allocated(calc%occ))deallocate(calc%occ)
        allocate(calc%occ(nao), source=0.0_wp)
        call gfn0_gen_occ(nel,nao,calc%config,calc%occ)
      endif
    end subroutine gfn0_init3
    subroutine gfn0_wbos(calc,mol,iostatus)
      implicit none
      type(calculation_settings),intent(inout) :: calc
      type(coord),intent(in) :: mol
      integer,intent(out) :: iostatus
      iostatus = 0
      if(.not.calc%rdwbo) return
      if(allocated(calc%wbo))deallocate(calc%wbo)
      allocate(calc%wbo( mol%nat, mol%nat), source=0.0_wp)
      call gfn0_getwbos(calc%g0calc,mol%nat,calc%wbo)
    end subroutine gfn0_wbos
  end subroutine gfn0occ_engrad



!========================================================================================!
end module api_engrad
