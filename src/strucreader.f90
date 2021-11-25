!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020 Philipp Pracht
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

!=====================================================================================================!
! STRUCRD is a module for reading and writing molecular structures.
! 
! The source is organized as follows:
!   0. Variable declarations
!   1. Routines for reading and writing ensemble files/trajectories in the XYZ format
!   2. Routines for reading single structures in various formats
!   3. Routines for writing structures in various formats
!   4. Utility routines mainly used only within the module
!
! Currently supported formats:
!   .xyz (Xmol) files and trajectories (read and write)
!   coord (turbomole) files (read and write)
!   .sdf/.mol files (V2000, read only)
!   .pdb files (in development)
!
!=====================================================================================================!
module strucrd
      use iso_fortran_env, only : wp => real64
      use iso_c_binding
      implicit none

    !=========================================================================================! 
    !--- private module variables and parameters
      private
       integer :: i,j,k,l,ich,och,io
       logical :: ex

     !--- some constants and name mappings
       real(wp),parameter :: bohr  = 0.52917726_wp
       real(wp),parameter :: autokcal = 627.509541_wp
       !-- filetypes as integers
       integer,parameter :: tmcoord = 1  
       integer,parameter :: xmol    = 2
       integer,parameter :: sdf     = 3  !currently unused
       integer,parameter :: sdfV2000 = 31
       integer,parameter :: sdfV3000 = 32
       integer,parameter :: pdbfile     = 4  !currently unused
       ! [...]

    !--- private utility subroutines 
      private :: upperCase,lowerCase
      private :: convertlable,fextension,sgrep
 
    !=========================================================================================!
    !--- public subroutines
      public :: i2e
      public :: asym
        interface asym
        module procedure i2e
        end interface asym
      public :: e2i  
      public :: grepenergy
      public :: checkcoordtype

      public :: rdnat       !-- procedure to read number of atoms Nat
      public :: rdcoord     !-- read an input file, determine format automatically
      public :: rdxmol      !-- read a file in the Xmol (.xyz) format specifically
      public :: rdxmolselec !-- read only a certain structure in Xmol file

      public :: wrc0  !-- write a TM coord file
        interface wrc0
        module procedure wrc0_file
        module procedure wrc0_channel
        end interface wrc0
      public :: wrcoord  
        interface wrcoord
        module procedure wrc0_file
        module procedure wrc0_channel
        end interface wrcoord

      public :: wrxyz !-- write a XYZ coord file
        interface wrxyz
        module procedure wrxyz_file
        module procedure wrxyz_file_mask
        module procedure wrxyz_channel_energy
        module procedure wrxyz_channel
        end interface wrxyz

      public :: xyz2coord
      public :: coord2xyz  


      public :: openensembledummy
      public :: rdensembleparam   !-- read Nat and Nall for a XYZ trajectory
      public :: rdensemble        !-- read a XYZ trajectory
       interface rdensemble
       module procedure rdensemble_conf1
       module procedure rdensemble_conf2
       module procedure rdensemble_conf3

       module procedure rdensemble_mixed2
       end interface rdensemble
       
     public :: wrensemble
       interface wrensemble
       module procedure wrensemble_conf_energy
       module procedure wrensemble_conf_energy_comment
       end interface wrensemble

      public :: pdbdata
      public :: coord
      public :: ensemble  
      public :: coordline

!=========================================================================================!
   !coord class. contains a single structure in the PDB format.
   !coordinates by definition are in Angstroem.
   type :: pdbdata

   !--- data
       integer :: nat = 0
       integer :: frag = 0
   !--- arrays
       integer,allocatable  :: athet(:) !ATOM (1) or HETATM (2)
       character(len=4),allocatable :: pdbat(:) !PDB atom specifier
       character(len=3),allocatable :: pdbas(:) !PDB amino acid specifier
       integer,allocatable :: pdbfrag(:) !PDB fragment specifier
       character(len=1),allocatable :: pdbgrp(:)  !PDB group specifier
       real(wp),allocatable :: pdbocc(:) !PDB occupancy
       real(wp),allocatable :: pdbtf(:)  !PDB temperature factor

       contains
           procedure :: deallocate => deallocate_pdb !clear memory space
           procedure :: allocate => allocate_pdb  

   end type pdbdata
!=========================================================================================!
   !coord class. contains a single structure
   !by convention coordinates are in Bohr for a single structure!
   type :: coord

   !--- data
       integer :: nat  = 0
       real(wp) :: energy = 0.0_wp
   !--- arrays    
       integer,allocatable  :: at(:)    !atom types as integer, dimension will be at(nat)
       real(wp),allocatable :: xyz(:,:) !coordinates, dimension will be xyz(3,nat)
       character(len=:),allocatable :: comment !a comment line

   !--- (optional) PDB data
       type(pdbdata) :: pdb   

       contains
           procedure :: deallocate => deallocate_coord !clear memory space
           procedure :: open => opencoord !read an coord file
           procedure :: write => writecoord !write
           procedure :: append => appendcoord !append
           procedure :: get => getcoord !allocate & fill with data

   end type coord
!=========================================================================================!
   !ensemble class. contains all structures of an ensemble
   !by convention coordinates are in Angström for an ensemble!
   type :: ensemble

   !--- data
       integer :: nat  = 0             !number of total atoms
       integer :: nall = 0             !number of structures
       integer,allocatable :: vnat(:)     !used instead of nat if not all structures have the same      number of atoms, in which case nat will be  =maxval(vnat,1)

       integer,allocatable  :: at(:)      !atom types as integer, dimension will be at(nat)
       real(wp),allocatable :: xyz(:,:,:) !coordinates, dimension will be xyz(3,nat,nall)
       real(wp),allocatable :: er(:)   !energy of each structure, dimension will be eread(nall)

       real(wp)            :: g         !gibbs free energy
       real(wp)            :: s         !entropy
       real(wp),allocatable :: gt(:)    !gibbs free energy of each member
       real(wp),allocatable :: ht(:)    !enthalpy of each member
       real(wp),allocatable :: svib(:)  !vibrational entropy of each member
       real(wp),allocatable :: srot(:)  !rotational entropy of each member
       real(wp),allocatable :: stra(:)  !translational entropy of each member



       contains
           procedure :: deallocate => deallocate_ensembletype !clear memory space
           procedure :: open => openensemble !read an ensemble file
           procedure :: write => write_ensemble !write to file

   end type ensemble
!=========================================================================================!


contains
!=====================================================================================================!
!=====================================================================================================!
!  1. ROUTINES FOR READING ENTIRE ENSEMBLES (OR TRAJECTORIES)
!=====================================================================================================!
!=====================================================================================================!

subroutine openensembledummy(fname) !DUMMY FUNCTION FOR IMPLEMENTATION TESTING
      implicit none
      character(len=*),intent(in) :: fname
      integer :: nat
      integer,allocatable :: at(:)
      real(wp),allocatable :: xyz(:,:,:)
      real(wp),allocatable :: eread(:)
      integer :: nall
      logical :: conform

      type(coord) :: struc
      type(ensemble) :: trj
      character(len=:),allocatable :: dummy
      character(len=6) :: sym
      real(wp) :: crd(3)

      !inquire(file=fname,exist=ex)
      !if(.not.ex)then
      !  error stop 'ensemble file does not exist.'
      !endif

      call rdensembleparam(fname,nat,nall)

      write(*,*)trim(fname),nat,nall
      write(*,*) fextension(fname)

      if(nat>0 .and. nall>0)then
          allocate(at(nat),xyz(3,nat,nall),eread(nall))
          call rdensemble(fname,nat,nall,at,xyz,eread)

          call wrensemble('new.dummy.xyz',nat,nall,at,xyz,eread)

          call wrcoord('dum.coord',nat,at,xyz(:,:,1)/bohr)
          deallocate(eread,xyz,at)
      else
        error stop 'format error while reading ensemble file.'
      endif

      !call struc%open(fname)

      !call wrc0(6,struc%nat,struc%at,struc%xyz)
      !write(*,*)
      !call wrxyz(6,struc%nat,struc%at,struc%xyz)
      return
end subroutine openensembledummy

!==================================================================!
! subroutine rdensembleparam
! read a ensemble file and get some information from
! it:
! On Input: fname - name of the file, should be in
!                   the Xmol (*.xyz) format.
!
! On Output: nat  - number of atoms
!                   (if different sized structures are present,
!                    nat is the largest)
!            nall - number of structures
!            conform - (optional) do all structures
!                      have the same number of atoms?
!=================================================================!
subroutine rdensembleparam(fname,nat,nall,conform)
      implicit none
      character(len=*),intent(in) :: fname
      integer,intent(out) :: nat
      integer,intent(out) :: nall
      logical,optional :: conform
      logical :: conformdum
      integer :: dum,iosum
      integer :: natref 
      real(wp) :: x,y,z
      character(len=10) :: str
      conformdum = .true.
      nat = 0
      nall = 0
      inquire(file=fname,exist=ex)
      if(.not.ex)return
      open(newunit=ich,file=fname)
      do
        read(ich,*,iostat=io) dum
          if( io < 0 ) exit
          if( io > 0 ) cycle
          if( nat == 0 ) natref = dum
        read(ich,*,iostat=io)
          if( io < 0 ) exit
        iosum=0
        do i=1,dum
          read(ich,*,iostat=io) str,x,y,z
          if(io < 0) exit
          iosum=iosum+io 
        enddo
        if(iosum > 0) cycle
        nat = max(dum,nat)
        if(dum.ne.natref) conformdum =.false.
        nall = nall + 1
      enddo
      close(ich)
      if(present(conform))conform=conformdum
      return
end subroutine rdensembleparam

!==================================================================!
! subroutine rdensemble_conf1
! read a conformer ensemble/a MD trajectory, i.e., 
! all structures have the same number and order of atoms.
! version 1 also reads the energy
!=================================================================!
subroutine rdensemble_conf1(fname,nat,nall,at,xyz,eread)
      implicit none
      character(len=*),intent(in) :: fname
      integer,intent(in) :: nat
      integer,intent(in) :: nall
      integer :: at(nat)
      real(wp) :: xyz(3,nat,nall)
      real(wp) :: eread(nall)
      integer :: i,j
      integer :: io,ich
      integer :: dum
      real(wp) :: dum2
      character(len=512) :: line
      character(len=6) :: sym
 
      eread = 0.0_wp
      
      open(newunit=ich,file=fname)
      do i=1,nall
        read(ich,*,iostat=io) dum
          if( io < 0 ) exit
          if( io > 0 ) cycle
        read(ich,'(a)',iostat=io) line
          if( io < 0 ) exit
          eread(i) = grepenergy(line)
        do j=1,dum
          read(ich,'(a)',iostat=io)line
          if(io < 0) exit
          call coordline(line,sym,xyz(1:3,j,i))
          at(j) = e2i(sym)
        enddo
      enddo
      close(ich)

      if(io<0)then
          error stop 'error while reading ensemble file.'
      endif

      return
end subroutine rdensemble_conf1      

!==================================================================!
! subroutine rdensemble_conf2
! read a conformer ensemble/a MD trajectory, i.e.,
! all structures have the same number and order of atoms.
! version 2 does not read the energy
!=================================================================!
subroutine rdensemble_conf2(fname,nat,nall,at,xyz)
      implicit none
      character(len=*),intent(in) :: fname
      integer,intent(in) :: nat
      integer,intent(in) :: nall
      integer :: at(nat)
      real(wp) :: xyz(3,nat,nall)
      integer :: i,j
      integer :: io,ich
      integer :: dum
      character(len=512) :: line
      character(len=6) :: sym
      io=0
      open(newunit=ich,file=fname)
      do i=1,nall
        read(ich,*,iostat=io) dum
          if( io < 0 ) exit
          if( io > 0 ) cycle
        read(ich,'(a)',iostat=io) line
          if( io < 0 ) exit
        do j=1,dum
          read(ich,'(a)',iostat=io)line
          if(io < 0) exit
          call coordline(line,sym,xyz(1:3,j,i))
          at(j) = e2i(sym)
        enddo
      enddo
      close(ich)

      if(io<0)then
          error stop 'error while reading ensemble file.'
      endif

      return
end subroutine rdensemble_conf2


!==================================================================!
! subroutine rdensemble_conf3
! read a conformer ensemble/a MD trajectory, i.e.,
! all structures have the same number and order of atoms.
! version 3 saves the comment line for each structure
!=================================================================!
subroutine rdensemble_conf3(fname,nat,nall,at,xyz,comments)
      implicit none
      character(len=*),intent(in) :: fname
      integer,intent(in) :: nat
      integer,intent(in) :: nall
      integer :: at(nat)
      real(wp) :: xyz(3,nat,nall)
      character(len=*) :: comments(nall)
      integer :: i,j
      integer :: io,ich
      integer :: dum
      character(len=512) :: line
      character(len=6) :: sym
      io=0
      open(newunit=ich,file=fname)
      do i=1,nall
        read(ich,*,iostat=io) dum
          if( io < 0 ) exit
          if( io > 0 ) cycle
        read(ich,'(a)',iostat=io) line
          if( io < 0 ) exit
          comments(i)=trim(line)
        do j=1,dum
          read(ich,'(a)',iostat=io)line
          if(io < 0) exit
          call coordline(line,sym,xyz(1:3,j,i))
          at(j) = e2i(sym)
        enddo
      enddo
      close(ich)

      if(io<0)then
          error stop 'error while reading ensemble file.'
      endif

      return
end subroutine rdensemble_conf3


!==================================================================!
! subroutine rdensemble_mixed2
! read an ensemble of mixed strcutres, i.e., all stuctures
! can have a diferent number and order of atoms.
! version 2 does not read energies
!=================================================================!
subroutine rdensemble_mixed2(fname,natmax,nall,nats,ats,xyz)
      implicit none
      character(len=*),intent(in) :: fname
      integer,intent(in) :: natmax
      integer,intent(in) :: nall
      integer  :: nats(nall)
      integer  :: ats(natmax,nall)
      real(wp) :: xyz(3,natmax,nall)
      integer :: i,j
      integer :: io,ich
      integer :: dum
      character(len=512) :: line
      character(len=6) :: sym
      open(newunit=ich,file=fname)
      do i=1,nall
        read(ich,*,iostat=io) dum
          if( io < 0 ) exit
          if( io > 0 ) cycle
          nats(i) = dum
        read(ich,'(a)',iostat=io) line
          if( io < 0 ) exit
        do j=1,dum
          read(ich,'(a)',iostat=io)line
          if(io < 0) exit
          call coordline(line,sym,xyz(1:3,j,i))
          ats(j,i) = e2i(sym)
        enddo
      enddo
      close(ich)

      if(io<0)then
          error stop 'error while reading ensemble file.'
      endif

      return
end subroutine rdensemble_mixed2

!=================================================================!
! subroutine wrensemble_conf
! write a ensemble file/a trajectory from memory. 
!=================================================================!
subroutine wrensemble_conf(fname,nat,nall,at,xyz)
      implicit none
      character(len=*),intent(in) :: fname
      integer,intent(in) :: nat
      integer,intent(in) :: nall
      integer :: at(nat)
      real(wp) :: xyz(3,nat,nall)
      real(wp) :: er(nall)
      integer :: i,j
      integer :: io,ich
      integer :: dum
      real(wp) :: dum2
      character(len=512) :: line
      character(len=6) :: sym

      open(newunit=ich,file=fname,status='replace')
      do i=1,nall
         call wrxyz(ich,nat,at,xyz(:,:,i))
      enddo
      close(ich)

      return
end subroutine wrensemble_conf

!=================================================================!
! subroutine wrensemble_conf_energy
! write a ensemble file/a trajectory from memory. 
!=================================================================!
subroutine wrensemble_conf_energy(fname,nat,nall,at,xyz,er)
      implicit none
      character(len=*),intent(in) :: fname
      integer,intent(in) :: nat
      integer,intent(in) :: nall
      integer :: at(nat)
      real(wp) :: xyz(3,nat,nall)
      real(wp) :: er(nall)
      integer :: i,j
      integer :: io,ich
      integer :: dum
      real(wp) :: dum2
      character(len=512) :: line
      character(len=6) :: sym

      open(newunit=ich,file=fname,status='replace')
      do i=1,nall
         call wrxyz(ich,nat,at,xyz(:,:,i),er(i))
      enddo
      close(ich)

      return
end subroutine wrensemble_conf_energy

!=================================================================!
! subroutine wrensemble_conf_energy_comment
! write a ensemble file/a trajectory from memory.
!=================================================================!
subroutine wrensemble_conf_energy_comment(fname,nat,nall,at,xyz,er,comments)
      implicit none
      character(len=*),intent(in) :: fname
      integer,intent(in) :: nat
      integer,intent(in) :: nall
      integer :: at(nat)
      real(wp) :: xyz(3,nat,nall)
      real(wp) :: er(nall)
      character(len=*) :: comments(nall)
      integer :: i,j
      integer :: io,ich
      integer :: dum
      real(wp) :: dum2
      character(len=512) :: line
      character(len=6) :: sym

      open(newunit=ich,file=fname,status='replace')
      do i=1,nall
         write(line,'(2x,f18.8,2x,a)') er(i),trim(comments(i))
         call wrxyz(ich,nat,at,xyz(:,:,i),trim(line))
      enddo
      close(ich)

      return
end subroutine wrensemble_conf_energy_comment


!==================================================================!
! subroutine write_ensemble
! wrapper to write an ensemble from the "ensemble" class
!==================================================================!
subroutine write_ensemble(self,fname)
      implicit none
      class(ensemble) :: self
      character(len=*),intent(in) :: fname
      call wrensemble_conf_energy(fname,self%nat,self%nall,self%at,self%xyz,self%er)
      return
end subroutine write_ensemble

!==================================================================!
! subroutine deallocate_ensembletype
! is used to clear memory for the ensemble type
!==================================================================!
subroutine deallocate_ensembletype(self)
      implicit none
      class(ensemble) :: self
      self%nat = 0
      self%nall = 0
      if(allocated(self%vnat))deallocate(self%vnat)
      if(allocated(self%at))deallocate(self%at)
      if(allocated(self%xyz))deallocate(self%xyz)
      if(allocated(self%er))deallocate(self%er)
      if(allocated(self%gt))deallocate(self%gt)
      if(allocated(self%ht))deallocate(self%ht)
      if(allocated(self%svib))deallocate(self%svib)
      if(allocated(self%srot))deallocate(self%srot)
      if(allocated(self%stra))deallocate(self%stra)

      return
end subroutine deallocate_ensembletype

!==================================================================!
! subroutine openensemble
! is the open procedure for the "ensemble" class.
! a ensemble (trajectory) fname is read into a new ensemble object
!==================================================================!
subroutine openensemble(self,fname)
      implicit none
      class(ensemble) :: self
      character(len=*),intent(in) :: fname
      integer :: nat
      integer,allocatable :: at(:)
      real(wp),allocatable :: xyz(:,:,:)
      real(wp),allocatable :: eread(:)
      integer :: nall
      logical :: conform

      inquire(file=fname,exist=ex)
      if(.not.ex)then
        error stop 'ensemble file does not exist.'
      endif

      call rdensembleparam(fname,nat,nall)

      if(nat>0 .and. nall>0)then
          call self%deallocate()
          allocate(at(nat),xyz(3,nat,nall),eread(nall))
          call rdensemble(fname,nat,nall,at,xyz,eread)

          self%nat=nat
          self%nall=nall
          call move_alloc(at,self%at)
          call move_alloc(xyz,self%xyz)
          call move_alloc(eread,self%er)
      else
        error stop 'format error while reading ensemble file.'
      endif

      return
end subroutine openensemble




!=====================================================================================================!
!=====================================================================================================!
!  2. ROUTINES FOR READING SINGLE STRUCTURES (COORDS)
!=====================================================================================================!
!=====================================================================================================!

!============================================================!
! subroutine checkcoordtype
! try to identify the filetype of the coord type.
! first based on file extension, if that fails by
! a keyword within the file.
!============================================================!
subroutine checkcoordtype(fname,typint)
     implicit none
     character(len=*) :: fname
     integer,intent(out) :: typint
     typint = 0
     !-- check file extension first
     select case( fextension(fname) )
      case( '.coord','.COORD' )
       typint=tmcoord 
      case( '.xyz','.XYZ','.trj','.TRJ','.sorted' )
       typint=xmol
      case( '.sd','.sdf','.SDF','.mol','.MOL' )
       typint=sdf
       if(sgrep(fname,'V2000'))then
       typint=sdfV2000
       endif
       if(sgrep(fname,'V3000'))then
       typint=sdfV3000
       endif
      case( '.pdb','.PDB' )
       typint=pdbfile
      case default
       typint=0
     end select
     if(typint.ne.0)return !-- file extension was recognized
     !-- grep for keywords otherwise
     if(sgrep(fname,'$coord'))then
        typint=tmcoord
     else !--no match found
        typint=0
     endif
     return
end subroutine checkcoordtype


!============================================================!
! subroutine rdnat
! read number of atoms "nat" form file
!
! On Input: fname  - name of the coord file
!           ftype  - (OPTIONAL) format of the input coord file
!                    if ftype is not present, it is determined
! On Output: nat   - number of atoms
!============================================================!
subroutine rdnat(fname,nat,ftype)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(out) :: nat
    integer,optional :: ftype
    integer :: ftypedum
    integer :: ich,i,j,io
    logical :: ex
    character(len=256) :: atmp
    nat=0
    inquire(file=fname,exist=ex)
    if(.not.ex)then
      error stop 'file does not exist.'
    endif
    if(present(ftype))then
     ftypedum=ftype
    else
     call checkcoordtype(fname,ftypedum)
    endif 
    open(newunit=ich,file=fname) 
    select case( ftypedum )
    !--- *.xyz files
      case( xmol )
        read(ich,*,iostat=io) nat
    !--- TM coord file    
      case( tmcoord )
        do
          read(ich,'(a)',iostat=io) atmp
          if(io < 0) exit
          atmp = adjustl(atmp)
          if(index(atmp,"$coord").eq.1)exit
        enddo
        do
          read(ich,'(a)',iostat=io) atmp
          if(io < 0) exit
          atmp = adjustl(atmp) 
          if(atmp(1:1) == '$') exit
          nat=nat+1
        enddo
    !--- sdf V2000 (or *.mol) file        
      case( sdfV2000 )
        do i=1,3 !-- first three comment lines
         read(ich,'(a)',iostat=io) atmp
         if(io < 0) exit
        enddo 
        read(ich,'(a)',iostat=io) atmp
        if(index(atmp,'V2000').ne.0)then
         read(atmp,'(i3)')nat !- first argument is nat
        endif 
    !--- sdf V3000 file
      case( sdfV3000 )
        do
         read(ich,'(a)',iostat=io) atmp
         if(io < 0) exit
         if((index(atmp,'V30').ne.0) .and. &
         &  (index(atmp,'COUNTS').ne.0)) then
           j=index(atmp,'COUNTS') + 6
           k=len_trim(atmp)
           atmp=atmp(j:k)       
           atmp=adjustl(atmp)
           read(atmp,*) nat
         endif
        enddo
    !--- pdb file    
      case( pdbfile )
         !write(*,*) 'PDB file format not supported yet.'
         nat=0 
         do
          read(ich,'(a)',iostat=io) atmp
          if(io < 0) exit
          if((index(atmp,'ATOM').eq.1) .or. &
          &  (index(atmp,'HETATM').eq.1)) then
           nat=nat+1
          endif
         enddo
      case default
        continue
    end select
    close(ich)
    return
end subroutine rdnat

!============================================================!
! subroutine rdcoord
! read in a structure. The format is determined automatically
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!
! On Output: at   - atom number as integer
!            xyz  - coordinates (always in Bohr)
!            energy - (OPTIONAL) if present, try to get energy
!                      mainly from xyz files
!============================================================!
subroutine rdcoord(fname,nat,at,xyz,energy)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nat
    integer,intent(inout)  :: at(nat)
    real(wp),intent(inout) :: xyz(3,nat)
    real(wp),optional :: energy
    character(len=256) :: atmp
    integer :: ftype
    type(pdbdata) :: pdbdummy

    !--- determine the file type
    call checkcoordtype(fname,ftype)

    select case( ftype )
     case( tmcoord )  !-- TM coord file, is already in Bohr
       call rdtmcoord(fname,nat,at,xyz)
     case( xmol )     !-- XYZ file, is Angström, needs conversion
       if(present(energy))then
        call rdxmol(fname,nat,at,xyz,atmp)   
        energy = grepenergy(atmp)
       else
        call rdxmol(fname,nat,at,xyz)
       endif
       xyz=xyz/bohr
     case( sdfV2000 )      !-- SDF/MOL V2000 file, also Angström
       call rdsdf(fname,nat,at,xyz)    
       xyz=xyz/bohr
     case( sdfV3000 )     !-- SDF V3000 file, Angström
       call rdsdfV3000(fname,nat,at,xyz)
       xyz=xyz/bohr
     case( pdbfile )  !-- PDB file, Angström
       !error stop 'PDB file format not supported yet.'
       call rdPDB(fname,nat,at,xyz,pdbdummy)
       xyz=xyz/bohr
       call pdbdummy%deallocate()  
     case default
      continue
    end select

    return
end subroutine rdcoord

!============================================================!
! subroutine rdtmcoord
! read a struncture in the TM coord style.
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!
! On Output: at   - atom number as integer
!            xyz  - coordinates (always in Bohr)
!============================================================!
subroutine rdtmcoord(fname,nat,at,xyz)
     implicit none
     character(len=*),intent(in) :: fname
     integer,intent(in) :: nat
     integer,intent(inout)  :: at(nat)
     real(wp),intent(inout) :: xyz(3,nat)
     character(len=6) :: sym
     integer :: ich,io,i
     character(len=256) :: atmp  
     open(newunit=ich,file=fname)
     do                                     
       read(ich,'(a)',iostat=io) atmp
       if(io < 0) exit
       atmp = adjustl(atmp)
       if(index(atmp,"$coord").eq.1)exit
     enddo
     do i=1,nat
       read(ich,'(a)',iostat=io) atmp
       if(io < 0) exit
       atmp = adjustl(atmp) 
       if(atmp(1:1) == '$') exit
       call coordline(atmp,sym,xyz(1:3,i))
       at(i) = e2i(sym)
     enddo
     close(ich)
     return
end subroutine rdtmcoord

!============================================================!
! subroutine rdxmol
! read a struncture in the *.xyz (Xmol) style.
! The commentary (second) line is ignored
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!
! On Output: at   - atom number as integer
!            xyz  - coordinates (in Angström)
!            comment - (OPTIONAL) commentary line of the file
!============================================================!
subroutine rdxmol(fname,nat,at,xyz,comment)
     implicit none
     character(len=*),intent(in) :: fname
     integer,intent(in) :: nat
     integer,intent(inout)  :: at(nat)
     real(wp),intent(inout) :: xyz(3,nat)
     character(len=*),optional :: comment
     character(len=6) :: sym
     integer :: ich,io,i
     integer :: dum
     character(len=256) :: atmp  
     open(newunit=ich,file=fname)
     read(ich,*,iostat=io) dum
     if( nat .ne. dum)then
         error stop 'error while reading input coordinates'
     endif
     read(ich,'(a)') atmp !--commentary line
     if(present(comment)) comment=trim(adjustl(atmp))
     do i=1,nat
       read(ich,'(a)',iostat=io) atmp
       if(io < 0) exit
       atmp = adjustl(atmp) 
       call coordline(atmp,sym,xyz(1:3,i))
       at(i) = e2i(sym)
     enddo
     close(ich)
     return
end subroutine rdxmol

!============================================================!
! subroutine rdsdf
! read a struncture in the .sdf/.mol V2000 style.
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!
! On Output: at   - atom number as integer
!            xyz  - coordinates (in Angström)
!            comment - (OPTIONAL) commentary line of the file
!============================================================!
subroutine rdsdf(fname,nat,at,xyz,comment)
     implicit none
     character(len=*),intent(in) :: fname
     integer,intent(in) :: nat
     integer,intent(inout)  :: at(nat)
     real(wp),intent(inout) :: xyz(3,nat)
     character(len=*),optional :: comment
     character(len=6) :: sym
     integer :: ich,io,i
     integer :: dum
     character(len=256) :: atmp  
     open(newunit=ich,file=fname)
     read(ich,'(a)',iostat=io) atmp
     read(ich,'(a)',iostat=io) atmp
     read(ich,'(a)',iostat=io) atmp
     if(present(comment)) comment=trim(adjustl(atmp))
     read(ich,'(i3)',iostat=io) dum
     if( nat .ne. dum)then
         error stop 'error while reading input coordinates'
     endif
     do i=1,nat
       read(ich,'(a)',iostat=io) atmp
       if(io < 0) exit
       atmp = adjustl(atmp) 
       call coordline(atmp,sym,xyz(1:3,i))
       at(i) = e2i(sym)
     enddo
     close(ich)
     return                                        
end subroutine rdsdf

!============================================================!
! subroutine rdsdfV3000
! read a struncture in the .sdf/.mol V3000 style.
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!
! On Output: at   - atom number as integer
!            xyz  - coordinates (in Angström)
!            comment - (OPTIONAL) commentary line of the file
!============================================================!
subroutine rdsdfV3000(fname,nat,at,xyz,comment)
     implicit none
     character(len=*),intent(in) :: fname
     integer,intent(in) :: nat
     integer,intent(inout)  :: at(nat)
     real(wp),intent(inout) :: xyz(3,nat)
     character(len=*),optional :: comment
     character(len=6) :: sym
     integer :: ich,io,i,j,k,l
     integer :: dum
     character(len=256) :: atmp  
     character(len=32) :: btmp
     open(newunit=ich,file=fname)
     read(ich,'(a)',iostat=io) atmp
     read(ich,'(a)',iostat=io) atmp
     read(ich,'(a)',iostat=io) atmp
     if(present(comment)) comment=trim(adjustl(atmp))
     do
         read(ich,'(a)',iostat=io) atmp
         if(io < 0) exit
         if((index(atmp,'V30').ne.0) .and. &
         &  (index(atmp,'COUNTS').ne.0)) then
           j=index(atmp,'COUNTS') + 6
           k=len_trim(atmp)
           atmp=atmp(j:k)
           atmp=adjustl(atmp)
           read(atmp,*) dum
         endif
         if((index(atmp,'V30').ne.0) .and. &
         &  (index(atmp,'ATOM').ne.0)) then
         exit
         endif
     enddo
     if( nat .ne. dum)then
         error stop 'error while reading input coordinates'
     endif
     do i=1,nat
       read(ich,'(a)',iostat=io) atmp
       if(io < 0) exit
       write(btmp,'(i0)') i
       l = len_trim(btmp) + 1
       j = index(atmp,'V30') + 3
       k = len_trim(atmp)
       atmp = atmp(j:k)
       atmp = adjustl(atmp) 
       atmp = atmp(l:k)
       call coordline(atmp,sym,xyz(1:3,i))
       at(i) = e2i(sym)
     enddo
     close(ich)
     return                                        
end subroutine rdsdfV3000

!============================================================!
! subroutine rdPDB
! read a struncture in the .PDB style.
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!
! On Output: at   - atom number as integer
!            xyz  - coordinates (in Angström)
!            pdb  - pdbdata object
!============================================================!
subroutine rdPDB(fname,nat,at,xyz,pdb)
     implicit none
     character(len=*),intent(in) :: fname
     integer,intent(in) :: nat
     integer,intent(inout)  :: at(nat)
     real(wp),intent(inout) :: xyz(3,nat)
     type(pdbdata) :: pdb
     character(len=2) :: sym
     integer :: ich,io,i,j,k,l
     integer :: dum
     character(len=256) :: atmp
     character(len=32) :: btmp
     character(len=6) :: dum1
     character(len=1) :: dum2,dum3,pdbgp
     character(len=3) :: pdbas
     character(len=2) :: dum4
     character(len=4) :: pdbat
     real(wp) :: r1,r2
     call pdb%allocate(nat)
     open(newunit=ich,file=fname)
     k=0
     do
       read(ich,'(a)',iostat=io) atmp
       if(io < 0) exit
       if((index(atmp,'ATOM').eq.1) .or. &
       &  (index(atmp,'HETATM').eq.1)) then
        k=k+1
        read(atmp,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)') &
        &  dum1,i,pdbat,dum2,pdbas,pdbgp,j,dum3,xyz(1:3,k),r1,r2,sym,dum4     
        at(k)=e2i(sym)
        pdb%pdbat(k) = pdbat
        pdb%pdbas(k) = pdbas
        pdb%pdbgrp(k) = pdbgp
        pdb%pdbfrag(k) = j
        pdb%pdbocc(k) = r1
        pdb%pdbtf(k) = r2
        endif
       enddo
     close(ich)
     return
end subroutine rdPDB     

!============================================================!
! subroutine rdxmolselec
! Read a file with multiple structures in the *.xyz (Xmol) style.
! Picks one structure.
! The commentary (second) line is ignored
!
! On Input: fname  - name of the coord file
!           m      - position of the desired structure
!           nat    - number of atoms
!
! On Output: at   - atom number as integer
!            xyz  - coordinates (in Angström)
!============================================================!

subroutine rdxmolselec(fname,m,nat,at,xyz,comment)
     implicit none
     character(len=*),intent(in) :: fname
     integer,intent(in) :: nat, m
     integer,intent(inout)  :: at(nat)
     real(wp),intent(inout) :: xyz(3,nat)
     character(len=*),optional :: comment
     character(len=6) :: sym
     integer :: ich,io,i,j
     integer :: dum
     character(len=256) :: atmp  

     open(newunit=ich,file=fname)

     do j=1,m
       read(ich,*,iostat=io) dum
       if( nat .ne. dum)then
           error stop 'error while reading input coordinates'
       endif
       read(ich,'(a)') atmp !--commentary line
       if(present(comment)) comment=trim(adjustl(atmp))
       do i=1,nat
         read(ich,'(a)',iostat=io) atmp
         if(io < 0) exit
         atmp = adjustl(atmp) 
         call coordline(atmp,sym,xyz(1:3,i))
         at(i) = e2i(sym)
       enddo
     end do
     close(ich)
     xyz = xyz/bohr
     return
end subroutine rdxmolselec

!==================================================================!
! subroutine deallocate_coord
! is used to clear memory for the coord type
!==================================================================!
subroutine deallocate_coord(self)
      implicit none
      class(coord) :: self
      self%nat = 0
      if(allocated(self%at))deallocate(self%at)
      if(allocated(self%xyz))deallocate(self%xyz)
      call self%pdb%deallocate()
      return
end subroutine deallocate_coord

!==================================================================!
! subroutine deallocate_pdb
! is used to clear memory for the pdbdata type
!==================================================================!
subroutine deallocate_pdb(self)
      implicit none
      class(pdbdata) :: self
      self%nat = 0
      self%frag = 0
      if(allocated(self%athet))  deallocate(self%athet) 
      if(allocated(self%pdbat))  deallocate(self%pdbat) 
      if(allocated(self%pdbas))  deallocate(self%pdbas) 
      if(allocated(self%pdbfrag))deallocate(self%pdbfrag)
      if(allocated(self%pdbgrp)) deallocate(self%pdbgrp) 
      if(allocated(self%pdbocc)) deallocate(self%pdbocc) 
      if(allocated(self%pdbtf))  deallocate(self%pdbtf)  
      return
end subroutine deallocate_pdb

!==================================================================!
! subroutine allocate_pdb
! is used to clear memory for the pdbdata type
!==================================================================!
subroutine allocate_pdb(self,nat)
      implicit none
      class(pdbdata) :: self
      integer :: nat
      call deallocate_pdb(self)
      self%nat = nat
      allocate(self%athet(nat)) 
      allocate(self%pdbat(nat)) 
      allocate(self%pdbas(nat)) 
      allocate(self%pdbfrag(nat))
      allocate(self%pdbgrp(nat)) 
      allocate(self%pdbocc(nat)) 
      allocate(self%pdbtf(nat))  
      return
end subroutine allocate_pdb

!==================================================================!
! subroutine opencoord
! is the open procedure for the "coord" class.
!==================================================================!
subroutine opencoord(self,fname)
      implicit none
      class(coord) :: self
      character(len=*),intent(in) :: fname
      integer :: nat,i
      integer,allocatable :: at(:)
      real(wp),allocatable :: xyz(:,:)
      integer :: nall
      logical :: conform
      integer :: ftype

      inquire(file=fname,exist=ex)
      if(.not.ex)then
        error stop 'coord file does not exist.'
      endif

      call self%deallocate()

      call checkcoordtype(fname,ftype)
      call rdnat(fname,nat)

      if(nat>0)then
          allocate(at(nat),xyz(3,nat))
          if(ftype==pdbfile)then
          call rdPDB(fname,nat,at,xyz,self%pdb)
          xyz=xyz/bohr
          else
          call rdcoord(fname,nat,at,xyz)
          endif

          self%nat=nat
          call move_alloc(at,self%at)
          call move_alloc(xyz,self%xyz)
      else
        error stop 'format error while reading coord file.'
      endif

      return
end subroutine opencoord
!==================================================================!
! subroutine getcoord
! allocate "coord" class and fill with data
!==================================================================!
subroutine getcoord(self,convfac,nat,at,xyz)
      implicit none
      class(coord) :: self
      real(wp),intent(in) :: convfac          
      integer,intent(in)  :: nat
      integer,intent(in)  :: at(nat)
      real(wp),intent(in) :: xyz(3,nat)
      call self%deallocate()
      allocate(self%at(nat))
      allocate(self%xyz(3,nat))
      self%nat=nat
      self%at=at
      self%xyz=xyz/convfac
      return
end subroutine getcoord



!=====================================================================================================!
!=====================================================================================================!
!  3. ROUTINES FOR WRITING STRUCTURES AND CONVERTING THEM
!=====================================================================================================!
!=====================================================================================================!

!============================================================!
! subroutine wrc0_file
! this is the typical quick write routine for TM coord files
! version for writing directly to a new file
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!           at   - atom number as integer
!           xyz  - coordinates (in Bohr)
!
! On Output: file written to "fname"
!============================================================!
subroutine wrc0_file(fname,nat,at,xyz)
     implicit none
     character(len=*) :: fname
     integer :: nat
     integer :: at(nat)
     real(wp) ::  xyz(3,nat)
     open(newunit=ich,file=fname,status='replace')
     write(ich,'(''$coord'')')
     do j=1,nat
        write(ich,'(3F24.12,5x,a2)') xyz(1:3,j),i2e(at(j),'lc')
     enddo
     write(ich,'(''$end'')')
     close(ich)
     return
end subroutine wrc0_file

!============================================================!
! subroutine wrc0_channel
! this is the typical quick write routine for TM coord files
! version for writing to a output channel
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!           at   - atom number as integer
!           xyz  - coordinates (in Bohr)
!
! On Output: file written to "fname"
!============================================================!
subroutine wrc0_channel(ch,nat,at,xyz)
     implicit none
     !character(len=*) :: fname
     integer :: ch
     integer :: nat
     integer :: at(nat)
     real(wp) ::  xyz(3,nat)
     !open(newunit=ich,file=fname)
     write(ch,'(''$coord'')')
     do j=1,nat
        write(ch,'(3F24.12,5x,a2)') xyz(1:3,j),i2e(at(j),'lc')
     enddo
     write(ch,'(''$end'')')
     !close(ch)
     return
end subroutine wrc0_channel

!============================================================!
! subroutine wrxyz_file
! this is the typical quick write routine for TM coord files
! version for writing directly to a new file
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!           at   - atom number as integer
!           xyz  - coordinates (in Angström)
!           comment - (OPTIONAL) comment line
!
! On Output: file written to "fname"
!============================================================!
subroutine wrxyz_file(fname,nat,at,xyz,comment)
     implicit none
     character(len=*) :: fname
     integer :: nat
     integer :: at(nat)
     real(wp) ::  xyz(3,nat)
     character(len=*),optional :: comment
     open(newunit=ich,file=fname,status='replace')
     write(ich,'(2x,i0)')nat
     if(present(comment))then
         write(ich,'(a)') trim(comment)
     else
         write(ich,*)
     endif
     do j=1,nat
      write(ich,'(1x,a2,1x,3f20.10)')i2e(at(j),'nc'),xyz(1:3,j)
     enddo
     close(ich)
     return
end subroutine wrxyz_file


!============================================================!
! subroutine wrxyz_file_mask
! this is the typical quick write routine for TM coord files
! version for writing directly to a new file
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!           at   - atom number as integer
!           xyz  - coordinates (in Angström)
!           mask - a mask to determine to write which atoms
!           comment - (OPTIONAL) comment line
!
! On Output: file written to "fname"
!============================================================!
subroutine wrxyz_file_mask(fname,nat,at,xyz,mask,comment)
     implicit none
     character(len=*) :: fname
     integer :: nat
     integer :: at(nat)
     real(wp) ::  xyz(3,nat)
     logical :: mask(nat)
     integer :: maskednat
     character(len=*),optional :: comment
     open(newunit=ich,file=fname,status='replace')
     maskednat = count(mask(:))
     write(ich,'(2x,i0)')maskednat
     if(present(comment))then
         write(ich,'(a)') trim(comment)
     else
         write(ich,*)
     endif
     do j=1,nat
      if(mask(j))then
       write(ich,'(1x,a2,1x,3f20.10)')i2e(at(j),'nc'),xyz(1:3,j)
      endif
     enddo
     close(ich)
     return
end subroutine wrxyz_file_mask


!============================================================!
! subroutine wrxyz_channel
! this is the typical quick write routine for xyz files
! version for writing to a output channel
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!           at   - atom number as integer
!           xyz  - coordinates (in Angström)
!           comment - (OPTIONAL) the comment line
!
! On Output: file written to "fname"
!============================================================!
subroutine wrxyz_channel(ch,nat,at,xyz,comment)
     implicit none
     !character(len=*) :: fname
     integer :: ch
     integer :: nat
     integer :: at(nat)
     real(wp) ::  xyz(3,nat)
     character(len=*),optional :: comment
     !open(newunit=ich,file=fname)
     write(ch,'(2x,i0)') nat
     if(present(comment))then
       write(ch,'(a)') trim(comment)
     else
       write(ch,*)
     endif
     do j=1,nat
      write(ch,'(1x,a2,1x,3f20.10)')i2e(at(j),'nc'),xyz(1:3,j)
     enddo
     return
end subroutine wrxyz_channel

!============================================================!
! subroutine wrxyz_channel
! this is the typical quick write routine for xyz files
! version for writing to a output channel
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!           at   - atom number as integer
!           xyz  - coordinates (in Angström)
!           er   - energy
!
! On Output: file written to "fname"
!============================================================!
subroutine wrxyz_channel_energy(ch,nat,at,xyz,er)
     implicit none
     !character(len=*) :: fname
     integer :: ch
     integer :: nat
     integer :: at(nat)
     real(wp) ::  xyz(3,nat)
     real(wp) :: er
     !open(newunit=ich,file=fname)
     write(ch,'(2x,i0)') nat
     write(ch,'(2x,f18.8)') er
     do j=1,nat
      write(ch,'(1x,a2,1x,3f20.10)')i2e(at(j),'nc'),xyz(1:3,j)
     enddo
     return
end subroutine wrxyz_channel_energy

!============================================================!
! subroutine xyz2coord
! simple conversion of a xyz to a coord file.
!
! On Input: iname  - name of the xyz file
!           oname  - name of the coord file
!
! On Output: file written to "oname"
!============================================================!
subroutine xyz2coord(iname,oname)
     implicit none
     character(len=*) :: iname
     character(len=*) :: oname
     type(coord) :: struc
     call struc%open(iname)
     call wrc0(oname,struc%nat,struc%at,struc%xyz)
     call struc%deallocate()
     return
end subroutine xyz2coord

!============================================================!
! subroutine coord2xyz
! simple conversion of a coord to a xyz file.
!
! On Input: iname  - name of the coord file
!           oname  - name of the xyz file
!
! On Output: file written to "oname"
!============================================================!
subroutine coord2xyz(iname,oname)
     implicit none
     character(len=*) :: iname
     character(len=*) :: oname
     type(coord) :: struc
     call struc%open(iname)
     struc%xyz=struc%xyz*bohr !to Angström
     call wrxyz(oname,struc%nat,struc%at,struc%xyz)
     call struc%deallocate()
     return
end subroutine coord2xyz

!==================================================================!
! subroutine writecoord
! is the write procedure for the "coord" class.
!==================================================================!
subroutine writecoord(self,fname)
      implicit none
      class(coord) :: self
      character(len=*),intent(in) :: fname
      if(.not.allocated(self%xyz))then
          write(*,*) 'Cannot write ',trim(fname),'. not allocated'
      endif
      if(index(fname,'.xyz').ne.0)then
        self%xyz=self%xyz*bohr !to Angström
        call wrxyz(fname,self%nat,self%at,self%xyz)
        self%xyz = self%xyz/bohr !back
      else
       call wrc0(fname,self%nat,self%at,self%xyz)
      endif
      return
end subroutine writecoord

!==================================================================!
! subroutine appendcoord
! is the write procedure for the "coord" class.
! coords will be written out in XYZ format!
!==================================================================!
subroutine appendcoord(self,io)
      implicit none
      class(coord) :: self
      integer :: io
        self%xyz=self%xyz*bohr !to Angström
        if(allocated(self%comment))then
            call wrxyz(io,self%nat,self%at,self%xyz,trim(self%comment))
        else
            call wrxyz(io,self%nat,self%at,self%xyz)
        endif
        self%xyz = self%xyz/bohr !back
      return
end subroutine appendcoord

!=====================================================================================================!
!=====================================================================================================!
!  4. GENERAL UTILITY ROUTINES
!=====================================================================================================!
!=====================================================================================================!

!============================================================!
! read a line of coordinates and determine by itself
! if the format is x,y,z,at or at,x,y,z
!============================================================!
subroutine coordline(line,sym,xyz)
      implicit none
      character(len=*) :: line
      character(len=*) :: sym
      real(wp) :: xyz(3)
      character(len=len_trim(line)) :: dum
      integer :: io

      read(line,*,iostat=io) xyz(1:3),sym
      if(io.ne.0)then
        read(line,*,iostat=io) sym,xyz(1:3)
        if(io.ne.0)then
          error stop 'error while reading coord line'
        endif
      endif    

      return
end subroutine coordline

!============================================================!
! convert a string into uppercase
!============================================================!
function upperCase(s)
! wandelt Kleinbuchstaben in Großbuchstaben um
    implicit none
    character(len=*), intent(in) :: s
    character(len=:),allocatable :: sout
    character(len=:),allocatable :: upperCase
    integer :: ic, i
    character(26), Parameter :: high = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
    sout = s
    do i = 1, LEN_TRIM(s)
        ic = INDEX(low, s(i:i))
        if (ic > 0) sout(i:i) = high(ic:ic)
    end do
    call move_alloc(sout,upperCase)
end function upperCase

!============================================================!
! convert a string into lowercase
!============================================================!
function lowerCase(s)
! wandelt Kleinbuchstaben in Großbuchstaben um
    implicit none
    character(len=*), intent(in) :: s
    character(len=:),allocatable :: sout
    character(len=:),allocatable :: lowerCase
    integer :: ic, i
    character(26), Parameter :: high = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
    sout = s
    do i = 1, LEN_TRIM(s)
        ic = INDEX(high, s(i:i))
        if (ic > 0) sout(i:i) = low(ic:ic)
    end do
    call move_alloc(sout,lowerCase)
end function lowerCase

!============================================================!
! split element lable if some isotope indicator was given
! and convert to uppercase
!============================================================!
function convertlable(s)
    implicit none
    character(len=*), intent(in) :: s
    character(len=:),allocatable :: sout
    character(len=:),allocatable :: convertlable
    integer :: ic, i
    character(14),parameter :: lab = '0123456789*_+-'
    character(26),parameter :: high = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26),parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
    sout = s
    do i=1,len_trim(s)
        ic = index(lab, s(i:i))
        if(ic > 0) sout(i:i) = ' '
        ic = index(low, s(i:i))
        if(ic > 0) sout(i:i) = high(ic:ic)
    end do
    sout = trim(adjustl(sout))
    call move_alloc(sout,convertlable)
end function convertlable

!============================================================!
! e2i is used to map the element (as a string) to integer
!============================================================!
integer function e2i(cin)
    implicit none
    character(len=*), intent(in) :: cin
    character(len=:),allocatable :: c
    integer :: iout
    c = trim(convertlable(cin))
    select case( c )
     case( 'H' ); iout = 1
     case( 'D' ); iout = 1
     case( 'T' ); iout = 1
     case( 'HE'); iout = 2
     case( 'LI'); iout = 3
     case( 'BE'); iout = 4
     case( 'B' ); iout = 5
     case( 'C' ); iout = 6
     case( 'N' ); iout = 7
     case( 'O' ); iout = 8
     case( 'F' ); iout = 9
     case( 'NE'); iout = 10  
     case( 'NA'); iout = 11
     case( 'MG'); iout = 12
     case( 'AL'); iout = 13
     case( 'SI'); iout = 14
     case( 'P' ); iout = 15
     case( 'S' ); iout = 16
     case( 'CL'); iout = 17
     case( 'AR'); iout = 18
     case( 'K' ); iout = 19
     case( 'CA'); iout = 20
     case( 'SC'); iout = 21
     case( 'TI'); iout = 22
     case( 'V' ); iout = 23
     case( 'CR'); iout = 24
     case( 'MN'); iout = 25
     case( 'FE'); iout = 26
     case( 'CO'); iout = 27
     case( 'NI'); iout = 28
     case( 'CU'); iout = 29
     case( 'ZN'); iout = 30
     case( 'GA'); iout = 31
     case( 'GE'); iout = 32
     case( 'AS'); iout = 33
     case( 'SE'); iout = 34
     case( 'BR'); iout = 35
     case( 'KR'); iout = 36
     case( 'RB'); iout = 37
     case( 'SR'); iout = 38
     case( 'Y' ); iout = 39
     case( 'ZR'); iout = 40
     case( 'NB'); iout = 41
     case( 'MO'); iout = 42
     case( 'TC'); iout = 43
     case( 'RU'); iout = 44
     case( 'RH'); iout = 45
     case( 'PD'); iout = 46
     case( 'AG'); iout = 47
     case( 'CD'); iout = 48
     case( 'IN'); iout = 49
     case( 'SN'); iout = 50
     case( 'SB'); iout = 51
     case( 'TE'); iout = 52
     case( 'I' ); iout = 53
     case( 'XE'); iout = 54
     case( 'CS'); iout = 55
     case( 'BA'); iout = 56
     case( 'LA'); iout = 57
     case( 'CE'); iout = 58
     case( 'PR'); iout = 59
     case( 'ND'); iout = 60
     case( 'PM'); iout = 61
     case( 'SM'); iout = 62
     case( 'EU'); iout = 63
     case( 'GD'); iout = 64
     case( 'TB'); iout = 65
     case( 'DY'); iout = 66
     case( 'HO'); iout = 67
     case( 'ER'); iout = 68
     case( 'TM'); iout = 69
     case( 'YB'); iout = 70
     case( 'LU'); iout = 71
     case( 'HF'); iout = 72
     case( 'TA'); iout = 73
     case( 'W' ); iout = 74
     case( 'RE'); iout = 75
     case( 'OS'); iout = 76
     case( 'IR'); iout = 77
     case( 'PT'); iout = 78
     case( 'AU'); iout = 79
     case( 'HG'); iout = 80
     case( 'TL'); iout = 81
     case( 'PB'); iout = 82
     case( 'BI'); iout = 83
     case( 'PO'); iout = 84
     case( 'AT'); iout = 85
     case( 'RN'); iout = 86
     case default; iout=0
    end select     
    e2i = iout
end function e2i

!============================================================!
! i2e is used to map the element (as a integer) to a string
!============================================================!
character(len=2) function i2e(iin,oformat)
    implicit none
    integer,intent(in) :: iin
    character(len=:),allocatable :: c
    character(len=*),optional :: oformat
    select case( iin )
     case( 1 ); c = 'H'   
     case( 2 ); c = 'HE'  
     case( 3 ); c = 'LI'  
     case( 4 ); c = 'BE'  
     case( 5 ); c = 'B'   
     case( 6 ); c = 'C'   
     case( 7 ); c = 'N'   
     case( 8 ); c = 'O'   
     case( 9 ); c = 'F'   
     case( 10); c = 'NE'  
     case( 11); c = 'NA'  
     case( 12); c = 'MG'  
     case( 13); c = 'AL'  
     case( 14); c = 'SI'  
     case( 15); c = 'P'   
     case( 16); c = 'S'   
     case( 17); c = 'CL'  
     case( 18); c = 'AR'  
     case( 19); c = 'K'   
     case( 20); c = 'CA'  
     case( 21); c = 'SC'  
     case( 22); c = 'TI'  
     case( 23); c = 'V'   
     case( 24); c = 'CR'  
     case( 25); c = 'MN'  
     case( 26); c = 'FE'  
     case( 27); c = 'CO'  
     case( 28); c = 'NI'  
     case( 29); c = 'CU'  
     case( 30); c = 'ZN'  
     case( 31); c = 'GA'  
     case( 32); c = 'GE'  
     case( 33); c = 'AS'  
     case( 34); c = 'SE'  
     case( 35); c = 'BR'  
     case( 36); c = 'KR'  
     case( 37); c = 'RB'  
     case( 38); c = 'SR'  
     case( 39); c = 'Y'  
     case( 40); c = 'ZR'  
     case( 41); c = 'NB'  
     case( 42); c = 'MO'  
     case( 43); c = 'TC'  
     case( 44); c = 'RU'  
     case( 45); c = 'RH'  
     case( 46); c = 'PD'  
     case( 47); c = 'AG'  
     case( 48); c = 'CD'  
     case( 49); c = 'IN'  
     case( 50); c = 'SN'  
     case( 51); c = 'SB'  
     case( 52); c = 'TE'  
     case( 53); c = 'I'  
     case( 54); c = 'XE'  
     case( 55); c = 'CS'  
     case( 56); c = 'BA'  
     case( 57); c = 'LA'  
     case( 58); c = 'CE'  
     case( 59); c = 'PR'  
     case( 60); c = 'ND'  
     case( 61); c = 'PM'  
     case( 62); c = 'SM'  
     case( 63); c = 'EU'  
     case( 64); c = 'GD'  
     case( 65); c = 'TB'  
     case( 66); c = 'DY'  
     case( 67); c = 'HO'  
     case( 68); c = 'ER'  
     case( 69); c = 'TM'  
     case( 70); c = 'YB'  
     case( 71); c = 'LU'  
     case( 72); c = 'HF'  
     case( 73); c = 'TA'  
     case( 74); c = 'W'  
     case( 75); c = 'RE'  
     case( 76); c = 'OS'  
     case( 77); c = 'IR'  
     case( 78); c = 'PT'  
     case( 79); c = 'AU'  
     case( 80); c = 'HG'  
     case( 81); c = 'TL'  
     case( 82); c = 'PB'  
     case( 83); c = 'BI'  
     case( 84); c = 'PO'  
     case( 85); c = 'AT'  
     case( 86); c = 'RN'  
     case default; c = 'XX'
    end select     
    i2e = trim(c)
    if(present(oformat))then
      select case( oformat )
       case( 'lc','lowercase' )
         i2e = lowerCase(trim(c))
       case( 'nc','nicecase' )
         if(len_trim(c).gt.1)then
           c(2:2) = lowerCase(c(2:2))
           i2e = trim(c)
         endif
       case default
         continue
      end select
    endif    
end function i2e


!============================================================!
! get the file extension
!============================================================!
function fextension(s)
    implicit none
    character(len=*), intent(in) :: s !filename
    character(len=:),allocatable :: sout
    character(len=:),allocatable :: fextension !output
    integer :: ic, i
    sout = trim(adjustl(s))
    i=len_trim(sout)
    ic = index(sout,'.',.true.)
    if(ic.ne.0)then
      fextension = sout(ic:i) 
    else
      fextension='none'
    endif
    return
end function fextension

!============================================================!
! grep for a keyword within the file
!============================================================!
function sgrep(fname,key)
    implicit none
    character(len=*), intent(in) :: fname
    character(len=*), intent(in) :: key
    logical :: sgrep
    character(len=256) :: atmp
    integer :: ic, i, io
    sgrep=.false.
    open(newunit=ic,file=fname)
    do
      read(ic,'(a)',iostat=io) atmp
      if(io < 0) exit !EOF
      if(index(atmp,key).ne.0)then
        sgrep = .true.
        exit
      endif
    enddo
    close(ic)
    return
end function sgrep

!============================================================!
! grep the energy from a line of strings
!============================================================!
function grepenergy(line)
    implicit none
    real(wp) :: grepenergy
    character(len=*), intent(in) :: line
    real(wp) :: energy
    character(len=:),allocatable :: atmp
    integer :: ic, i, io
    atmp=trim(line)
    energy=0.0_wp
    do i=1,len_trim(atmp)
      if(len_trim(atmp).lt.1) exit
      read(atmp,*,iostat=io) energy
      if(io > 0)then
        atmp=atmp(2:) 
        atmp=adjustl(atmp)
        cycle
      else
        exit
      endif
    enddo
    grepenergy = energy
    return
end function grepenergy


!=====================================================================================================!
!=====================================================================================================!
! end of the module
!=====================================================================================================!
!=====================================================================================================!
end module strucrd
