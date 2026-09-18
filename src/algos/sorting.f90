!================================================================================!
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
!================================================================================!

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> Implementation for standalone sorting
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
!> Input/Output:
!>  env  -  crest's systemdata object
!>  tim  -  timer object
!>-----------------------------------------------
subroutine crest_sort(env,tim)
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use cregen_interface
  use iomod,only:catdel
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich
  logical :: pr,wr
  external :: CCEGEN
!========================================================================================!
  integer :: nall
  type(coord),allocatable :: structures(:)
  integer,allocatable :: groups(:)

!========================================================================================!
  select case (env%sortmode)
  case default
    write (stdout,'(a,a,a)',advance='no') '> Read ensemble ',trim(env%ensemblename),' ... '
    flush (stdout)
    call rdensemble(env%ensemblename,nall,structures)
    allocate (groups(nall),source=0)
    write (stdout,'(i0,a)') nall,' structures!'
  case ('irmsd','rmsd','hrmsd')
    write (stdout,'(a,a)',advance='no') '> Reading files ',trim(env%ensemblename)
    flush (stdout)
    write (stdout,'(a,a)') ' and ',trim(env%ensemblename2)
  case ('cluster')
    continue !> ccegen reads the ensemble internally
  end select
  write (stdout,*)

  env%confgo = .true.
!========================================================================================!
  call tim%start(11,'Sorting')

  select case (env%sortmode)

  case ('rmsd')
    call quick_rmsd_tool(trim(env%ensemblename),trim(env%ensemblename2),.false.)
    stop

  case ('hrmsd')
    call quick_rmsd_tool(trim(env%ensemblename),trim(env%ensemblename2),.true.)
    stop

  case ('irmsd')
    call irmsd_tool(trim(env%ensemblename),trim(env%ensemblename2),env%iinversion)
    stop

  case ('compare')
    call compare_ensembles(env)

  case ('isort')
!>--- Assigning structures to conformers based on RTHR,with canonical atom IDs
    call underline('Assigning conformers based on iRMSD and RTHR')
    call cregen_irmsd_sort(env,nall,structures,groups,allcanon=.true.,printlvl=2)

  case ('isort_noid')
!>--- Assigning structures to conformers based on RTHR, WITHOUT canonical atom IDs
    call underline('Assigning conformers based on iRMSD and RTHR')
    call cregen_irmsd_sort(env,nall,structures,groups,allcanon=.false.,printlvl=2)

  case ('all','allpair')
!>--- all unique pairs of the ensemble (only suitable for small ensembles)
    call underline('Running all unique pair RMSDs incl. atom permutation')
    call cregen_irmsd_all(nall,structures,printlvl=2,iinversion=env%iinversion)

  case ('cregen')
!>--- the original CREGEN procedure
    call newcregen(env,structurelist=structures)
    call catdel('cregen.out.tmp')

  case ('cluster')
!>--- PCA and k-means clustering
    call CCEGEN(env,.true.,env%ensemblename)

  case default
!>--- all unique pairs of the ensemble (only suitable for small ensembles)
    call cregen_irmsd_sort(env,nall,structures,groups,allcanon=.true.,printlvl=2)
  end select

!========================================================================================!
  call tim%stop(11)
  if (allocated(structures)) deallocate (structures)
  return
end subroutine crest_sort

!=========================================================================================!

subroutine irmsd_tool(fname1,fname2,iinversion)
!*******************************************************
!* irmsd_tool
!* Standalone implementation to compare two structures
!* with the iRMSD method.
!* This implementation should be called only on its own,
!* for ensemble-based processing see the CREGEN file
!*******************************************************
  use crest_parameters
  use strucrd
  use axis_module
  use irmsd_module
  use canonical_mod
  implicit none
  character(len=*),intent(in) :: fname1
  character(len=*),intent(in) :: fname2
  integer,intent(in) :: iinversion
  type(coord) :: mol,ref
  real(wp) :: rmsdval,tmpd(3),tmpdist
  integer :: i,ich
  type(rmsd_cache) :: rcache
  type(canonical_sorter) :: canmol
  type(canonical_sorter) :: canref
  logical :: mirror
  logical,parameter :: debug = .false.

  write (stdout,*) 'iRMSD algorithm'
  write (stdout,*) 'reference: ',trim(fname1)
  write (stdout,*) 'processed: ',trim(fname2)
  write (stdout,*)

  !> read the geometries
  call ref%open(trim(fname1))
  call mol%open(trim(fname2))

  !> move ref to CMA and align rotational axes
  call axis(ref%nat,ref%at,ref%xyz)

  !> allocate memory
  call rcache%allocate(ref%nat)

  !> canonical atom ranks
  call canref%init(ref,invtype='apsp+',heavy=.false.)
  !call canref%add_h_ranks(ref)
  rcache%stereocheck = .not. (canref%hasstereo(ref))
  call canref%shrink()
  write (stdout,*) 'false enantiomers possible?: ',rcache%stereocheck
  select case (iinversion)
  case (0)
    mirror = .true.
  case (1)
    mirror = .true.
    rcache%stereocheck = .true.
  case (2)
    mirror = .false.
    rcache%stereocheck = .false.
  end select
  write (stdout,*) 'allow inversion?:            ',mirror

  call canmol%init(mol,invtype='apsp+',heavy=.false.)
  !call canmol%add_h_ranks(mol)
  call canmol%shrink()

  !> check if we can work with the determined ranks
  if (checkranks(ref%nat,canref%rank,canmol%rank)) then
    write (stdout,*) 'using canonical atom identities as rank backend'
    rcache%rank(:,1) = canref%rank(:)
    rcache%rank(:,2) = canmol%rank(:)
    if (debug) then
      write (*,*) 'iRMSD ranks:'
      write (*,*) 'atom',' rank('//fname1//')',' rank('//fname2//')'
      do i = 1,ref%nat
        write (*,*) i,rcache%rank(i,1),rcache%rank(i,2)
      end do
      write (*,*)
    end if
  else
    !> if not, fall back to atom types
    write (stdout,*) 'using atom types as rank backend'
    call fallbackranks(ref,mol,ref%nat,rcache%rank)
  end if

  call min_rmsd(ref,mol,rcache=rcache,rmsdout=rmsdval,align=.true.)

  !> write the rotated and shifted coordinates to one file
  open (newunit=ich,file='irmsd.xyz')
  call ref%append(ich)
  call mol%append(ich)
  close (ich)
  write (stdout,*)
  write (stdout,*) 'aligned structures written to irmsd.xyz'
  write (stdout,*)

  do i = 1,mol%nat
    tmpd(:) = (mol%xyz(:,i)-ref%xyz(:,i))**2
    tmpdist = sqrt(sum(tmpd(:)))*autoaa
    if (tmpdist > 0.01_wp) then
      write (*,*) i,mol%at(i),tmpdist
    end if
  end do

  rmsdval = rmsdval*autoaa
  write (*,'(1x,a,f16.8)') 'Calculated iRMSD (Å):',rmsdval

  return
end subroutine irmsd_tool

