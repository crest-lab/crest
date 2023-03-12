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

module api_helpers

  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use strucrd
  use calc_type
  use iomod,only:makedir,directory_exist,remove
  !> APIs
  use tblite_api
  use gfn0_api
!=========================================================================================!
  implicit none
  public

!=========================================================================================!
!=========================================================================================!
contains    !>--- Module routines start here
!=========================================================================================!
!=========================================================================================!

  subroutine tblite_init(calc,loadnew)
    implicit none
    type(calculation_settings),intent(inout) :: calc
    logical,intent(out) :: loadnew
    loadnew = .false.
    if (.not.allocated(calc%wfn)) then
      allocate (calc%wfn)
      loadnew = .true.
    end if
    if (.not.allocated(calc%tbcalc)) then
      allocate (calc%tbcalc)
      loadnew = .true.
    end if
    if (.not.allocated(calc%ctx)) then
      allocate (calc%ctx)
      loadnew = .true.
    end if
    if (.not.allocated(calc%tbres)) then
      allocate (calc%tbres)
      loadnew = .true.
    end if
    if (calc%apiclean) loadnew = .true.
  end subroutine tblite_init
  subroutine tblite_wbos(calc,mol,iostatus)
    implicit none
    type(calculation_settings),intent(inout) :: calc
    type(coord),intent(in) :: mol
    integer,intent(out) :: iostatus
    iostatus = 0
    if (.not.calc%rdwbo) return
    if (allocated(calc%wbo)) deallocate (calc%wbo)
    allocate (calc%wbo(mol%nat,mol%nat),source=0.0_wp)
    call tblite_getwbos(calc%tbcalc,calc%wfn,calc%tbres,mol%nat,calc%wbo)
  end subroutine tblite_wbos

!========================================================================================!

  subroutine gfn0_init(calc,g0calc,loadnew)
    implicit none
    type(calculation_settings),intent(inout) :: calc
    type(gfn0_data),intent(inout),allocatable  :: g0calc
    logical,intent(out) :: loadnew
    loadnew = .false.
    if (.not.allocated(g0calc)) then
      allocate (g0calc)
      loadnew = .true.
    end if
    if (calc%apiclean) loadnew = .true.
  end subroutine gfn0_init
  subroutine gfn0_init2(mol,calc,g0calc)
    implicit none
    type(coord),intent(in) :: mol
    type(calculation_settings),intent(inout) :: calc
    type(gfn0_data),intent(inout)  :: g0calc
    if (allocated(calc%solvent).and.allocated(calc%solvmodel)) then
      call gfn0_addsettings(mol,g0calc,calc%solvent,calc%solvmodel)
    end if
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
    call g0calc%wfn%refresh_occu(nel,uhf)
    call gfn0_addsettings(mol,g0calc,etemp=calc%etemp)
  end subroutine gfn0_init3
  subroutine gfn0_wbos(calc,mol,iostatus)
    implicit none
    type(calculation_settings),intent(inout) :: calc
    type(coord),intent(in) :: mol
    integer,intent(out) :: iostatus
    iostatus = 0
    if (.not.calc%rdwbo) return
    if (allocated(calc%wbo)) deallocate (calc%wbo)
    allocate (calc%wbo(mol%nat,mol%nat),source=0.0_wp)
    call gfn0_getwbos(calc%g0calc,mol%nat,calc%wbo)
  end subroutine gfn0_wbos

!========================================================================================!

  subroutine gfn0occ_init(calc,g0calc,loadnew)
    implicit none
    type(calculation_settings),intent(inout) :: calc
    logical,intent(out) :: loadnew
    type(gfn0_data),intent(inout),allocatable  :: g0calc
    integer :: nel,nao,nlev
    loadnew = .false.
    if (.not.allocated(g0calc)) then
      allocate (g0calc)
      loadnew = .true.
    end if
    if (calc%apiclean) loadnew = .true.
  end subroutine gfn0occ_init
  subroutine gfn0occ_init2(mol,calc,g0calc)
    implicit none
    type(coord),intent(in) :: mol
    type(calculation_settings),intent(inout) :: calc
    type(gfn0_data),intent(inout),allocatable  :: g0calc
    integer :: nel,nao,nlev
    if (allocated(calc%solvent).and.allocated(calc%solvmodel)) then
      call gfn0_addsettings(mol,g0calc,calc%solvent,calc%solvmodel)
    end if
    call gfn0_addsettings(mol,g0calc,etemp=calc%etemp,loadwbo=calc%rdwbo)
  end subroutine gfn0occ_init2
  subroutine gfn0occ_init3(mol,calc,g0calc)
    implicit none
    type(coord),intent(in) :: mol
    type(calculation_settings),intent(inout) :: calc
    type(gfn0_data),intent(inout),allocatable  :: g0calc
    integer :: nel,nao,nlev
    if (.not.allocated(calc%occ)) then
      nel = g0calc%wfn%nel
      nao = g0calc%basis%nao
      if (allocated(calc%occ)) deallocate (calc%occ)
      allocate (calc%occ(nao),source=0.0_wp)
      call gfn0_gen_occ(nel,nao,calc%config,calc%occ)
    end if
  end subroutine gfn0occ_init3

!========================================================================================!

  subroutine gfnff_init(calc,loadnew)
    implicit none
    type(calculation_settings),intent(inout) :: calc
    logical,intent(out) :: loadnew
    integer :: nel,nao,nlev
    loadnew = .false.
    if (.not.allocated(calc%ff_dat)) then
      allocate (calc%ff_dat)
      loadnew = .true.
    end if
    if (allocated(calc%solvent)) then
      if (.not.allocated(calc%ff_dat%solvent)) then
        allocate (calc%ff_dat%solvent,source=trim(calc%solvent))
        calc%ff_dat%solvent = calc%solvent
      end if
    end if
    if (calc%apiclean) loadnew = .true.
  end subroutine gfnff_init

!========================================================================================!
end module api_helpers
