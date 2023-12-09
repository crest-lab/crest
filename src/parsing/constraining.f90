!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Philipp Pracht
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

!-----------------------------------------------------------------------------
! sort the .constrains file
!-----------------------------------------------------------------------------
subroutine sort_constrain_file(fname)
  use iomod,only:to_lower
  implicit none
  character(len=256),allocatable :: entire(:),setsave(:)
  character(len=256) :: atmp
  character(len=*) :: fname
  integer :: i
  integer :: io,ich,och,lines
  integer :: setb
  logical :: ex,follow,setset

  inquire (file=fname,exist=ex)
  if (.not. ex) return
  open (newunit=ich,file=fname)
  open (newunit=och,file='.tmpconstrains')
  lines = 0
  do
    read (ich,'(a)',iostat=io) atmp
    if (io .lt. 0) exit
    lines = lines + 1
  end do
  rewind (ich)
  allocate (setsave(lines),entire(lines))
  entire = ''
  setsave = ''
  do i = 1,lines
    read (ich,'(a)',iostat=io) entire(i)
    call To_lower(entire(i))
    entire(i) = adjustl(entire(i))
    if (io .lt. 0) exit
  end do
!---- seperate everything that is written inside of a $set statement
  setb = 1
  follow = .false.
  setset = .false.
  do i = 1,lines
    if (index(entire(i),'$end') .ne. 0) then
      follow = .false.
      entire(i) = ''  !"delete" from memory
    end if
    if ((index(entire(i),'$') .ne. 0) .and. follow) then
      follow = .false.
    end if
    if (follow) then
      setsave(setb) = entire(i)
      setb = setb + 1
      entire(i) = ''  !"delete" from memory
    end if
    if (index(entire(i),'$set') .ne. 0) then
      follow = .true.
      setset = .true.
      entire(i) = ''  !"delete" from memory
    end if
  end do
!---- write the new file in order, ignoring blank lines
!---- first everything not contained within $set-statements
  do i = 1,lines
    if (len_trim(entire(i)) .lt. 1) cycle
    if (index(entire(i),'$') .ne. 0) then
      write (och,'(a)') trim(entire(i))
    else
      write (och,'(2x,a)') trim(entire(i))
    end if
  end do
!---- then all the set stuff
  if (setset) then
    write (och,'(a)') '$set'
    do i = 1,lines
      if (len_trim(setsave(i)) .lt. 1) cycle
      write (och,'(2x,a)') trim(setsave(i))
    end do
  end if
!---- terminate all of it with a $end
  write (och,'(a)') '$end'
!---- deallocate and rename the file
  deallocate (entire,setsave)
  close (ich,status='delete')
  close (och)
  call rename('.tmpconstrains',fname)
end subroutine sort_constrain_file

!-----------------------------------------------------------------------------
! split a set-block line into the keyword and the argument
!-----------------------------------------------------------------------------
subroutine split_set_args(line,argument)
  implicit none
  character(len=*) :: line
  character(len=*) :: argument

  character(len=256) :: dummy
  character(len=1) :: digit
  integer :: i,j

  line = adjustl(line)

  do i = 1,len(trim(line))
    digit = line(i:i)
    if ((digit == ' ') .or. (digit == char(9))) then
      j = i + 1
      argument = adjustl(line(j:))
      j = i - 1
      dummy = trim(line(:j))
      exit
    end if
  end do
  line = trim(dummy)

!        deallocate(strings,floats)

end subroutine split_set_args

!========================================================================================!
! parse the atomlist string and identify specific atoms
!========================================================================================!
subroutine parse_atlist(line,selected,nat,natlist)
  implicit none
  character(len=*) :: line
  integer,intent(out) :: selected

  integer :: i,j,k
  character(len=1) :: digit
  character(len=10),parameter :: numbers = '0123456789'
  character(len=20) :: atom1,atom2
  integer :: at1,at2,nat

  integer,intent(inout) :: natlist(nat)
  logical :: rang,stor

  !>-- initialize
  natlist = 0

  stor = .false.
  rang = .false.

  atom1 = ''
  atom2 = ''
  at1 = 0
  at2 = 0
  selected = 0

  !>-- loop through the input string
  do i = 1,len(trim(line))
    digit = line(i:i)
    !> should exclude tabstops and blanks, 9 is ascii code for  tab
    if (index(numbers,digit) .ne. 0) then 
      if (.not. rang) then
        atom1 = trim(atom1)//digit
        stor = .true.
      else
        atom2 = trim(atom2)//digit
      end if
    end if
    if ((digit == ',') .or. (digit == ' ') .or. (digit == char(9)) &
    &  .or. (i .eq. len(trim(line)))) then  !comma,space,tab or last character
      if (.not. rang .and. stor) then !a single atom
        read (atom1,*) at1
        atom1 = ''
        if (at1 .ge. 1 .and. at1 .le. nat) then
          natlist(at1) = 1
        end if
        stor = .false.
      end if
      if (rang .and. stor) then !> a range of atoms
        read (atom2,*) at2
        atom2 = ''
        stor = .false.
        rang = .false.
        if (at1 .ge. 1 .and. at2 .ge. 1) then
          if (min(at1,at2) .le. nat) then
            k = max(at1,at2)
            do j = min(at1,at2),min(k,nat)   !> order is arbitrary
              natlist(j) = 1
            end do
          end if
        end if
      end if
    end if
    if ((digit == '-') .and. stor) then !> begin to select a range of atoms
      rang = .true.
      read (atom1,*) at1
      atom1 = ''
    end if
  end do
  
  !>-- count selected atoms
  selected = sum(natlist)

  return
end subroutine parse_atlist

!========================================================================================!
! parse atomlist string and identify selected atom types
!========================================================================================!
subroutine parse_atypelist(line,nselect,atlist)
  use strucrd,only:i2e,e2i
  implicit none
  !> Input
  character(len=*) :: line
  !> Output
  integer,intent(out) :: nselect
  logical,intent(inout) :: atlist(118)
  !> Local variables
  integer :: i,j,k
  character(len=1) :: digit
  character(len=10),parameter :: numbers = '0123456789'
  character(len=52),parameter :: letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
  character(len=20) :: atom1,atom2
  integer :: at1,at2
  logical :: rang,stor

  !>-- intitalize
  atlist = .false.

  stor = .false.
  rang = .false.

  atom1 = ''
  atom2 = ''
  at1 = 0
  at2 = 0
  nselect = 0

  !>-- loop through input string
  do i = 1,len(trim(line))
    digit = line(i:i)
    if (index(letters,digit) .ne. 0) then
      if (.not. rang) then
        atom1 = trim(atom1)//digit
        stor = .true.
      else
        atom2 = trim(atom2)//digit
      end if
    end if
    if ((digit == ',') .or. (digit == ' ') .or. (digit == char(9)) &
    &  .or. (i .eq. len(trim(line)))) then  !comma,space,tab or last character
      if (.not. rang .and. stor) then !a single atom type
        atom1 = trim(atom1)
        at1 = e2i(atom1)
        atom1 = ''
        if (at1 > 0) then
          atlist(at1) = .true.
        end if
        stor = .false.
      end if
      if (rang .and. stor) then !a range of atoms
        atom2 = trim(atom2)
        at2 = e2i(atom2)
        atom2 = ''
        stor = .false.
        rang = .false.
        if (at1 > 0 .and. at2 > 0) then
          k = max(at1,at2)
          do j = min(at1,at2),k   !order is arbitrary
            atlist(j) = .true.
          end do
        end if
      end if
    end if
    if ((digit == '-') .and. stor) then !begin to select a range of atoms
      rang = .true.
      atom1 = trim(atom1)
      at1 = e2i(atom1)
      atom1 = ''
    end if
  end do
  do i = 1,118
    if (atlist(i)) nselect = nselect + 1
  end do
  return
end subroutine parse_atypelist

!-----------------------------------------------------------------------------
! modified version of the rdcoord routine to only include selected atoms
!-----------------------------------------------------------------------------
subroutine rdcoord_reduced(fname,n,rn,xyz,iat,rlist)
  use iso_fortran_env,only:wp => real64
  use strucrd
  implicit none
  type(coord) :: struc
  integer :: n,rn
  real*8 :: xyz(3,rn)
  integer :: iat(rn)
  integer :: rlist(n)
  character(len=*) :: fname
  integer :: i,j
  call struc%open(fname) ! read complete structure
  j = 0
  do i = 1,n
    if (rlist(i) == 1) then
      j = j + 1
      xyz(1:3,j) = struc%xyz(1:3,i)
      iat(j) = struc%at(i)
    end if
  end do
  call struc%deallocate()
  return
end subroutine rdcoord_reduced



!-----------------------------------------------------------------------------
! sort the .constrains file
!-----------------------------------------------------------------------------
subroutine sort_constraints(cts)
  use crest_data
  implicit none
  type(constra) :: cts

  character(len=128) :: atmp,btmp
  integer :: i,j
  integer :: bpos
  logical :: follow

!--- some initialization issues
  cts%buff = ''
  bpos = 1
  follow = .false.

!--- collect dublicate blocks
  do i = 1,cts%ndim
    atmp = adjustl(cts%sett(i))
    if (index(atmp(1:1),'$') .ne. 0) then
      if (index(atmp,'$end') .ne. 0) cycle
      cts%buff(bpos) = atmp
      bpos = bpos + 1
      follow = .true.
      do j = 1,cts%ndim
        btmp = adjustl(cts%sett(j))
        if (index(btmp(1:1),'$') .ne. 0) then          !--- if its a new block don't follow ...
          follow = .false.
        end if
        if (index(trim(atmp),trim(btmp)) .ne. 0) then  !--- ... except if its the same keyword
          follow = .true.
          cycle
        end if
        if (follow) then
          if (trim(btmp) .eq. '') cycle                !--- cycle blanks
          write (cts%buff(bpos),'(2x,a)') trim(btmp)
          bpos = bpos + 1
        end if
      end do
    end if
  end do
  !write(cts%buff(bpos),'(''$end'')')

!----- write buffer to settings
  cts%sett = cts%buff

end subroutine sort_constraints

!-----------------------------------------------------------------------------
! read the constraints input file into the buffer
!-----------------------------------------------------------------------------
subroutine read_constrainbuffer(fname,cts)
  use crest_data
  implicit none
  type(constra) :: cts
  character(len=*) :: fname

  character(len=128) :: atmp
  integer :: i
  integer :: io,ich,lines
  integer :: ndim
  logical :: ex

!--- get the dimensions
  inquire (file=fname,exist=ex)
  if (.not. ex) return
  open (newunit=ich,file=fname)
  lines = 0
  do
    read (ich,'(a)',iostat=io) atmp
    if (io < 0) exit
    lines = lines + 1
  end do
  rewind (ich)
!--- allocate
  ndim = max(lines,10)
  call cts%allocate(ndim)

!--- read
  do i = 1,ndim
    read (ich,'(a)',iostat=io) atmp
    if (io < 0) exit
    atmp = adjustl(atmp)
    cts%sett(i) = atmp
  end do

  close (ich)
  return
end subroutine read_constrainbuffer

!-----------------------------------------------------------------------------
! append constraints to opened file
!-----------------------------------------------------------------------------
subroutine write_cts(ich,cts)
  use crest_data
  implicit none
  type(constra) :: cts
  integer :: ich
  integer :: i
!--- not really a "constraint", but convenient for the implementation:
!    change of the gbsa grid
  if (cts%ggrid) then
    if (allocated(cts%gbsagrid)) then
      write (ich,'(a)') '$gbsa'
      write (ich,'(2x,a,a)') 'gbsagrid=',cts%gbsagrid
    end if
  end if
!---- do it only if constaints are given, write them
  if (cts%used) then
    do i = 1,cts%ndim
      if (trim(cts%sett(i)) .ne. '') then
        write (ich,'(a)') trim(cts%sett(i))
      end if
    end do
  end if
!---- apply bondlength constraints globally
  if (cts%cbonds_global) then
    call write_cts_CBONDS(ich,cts)
  end if
  if (cts%dispscal_global) then
    call write_cts_DISP(ich,cts)
  end if
!---- handle the read-in DFT/GFN correction bias

  return
end subroutine write_cts

!-----------------------------------------------------------------------------
! append NCI constraints to opened file
!-----------------------------------------------------------------------------
subroutine write_cts_NCI(ich,cts)
  use crest_data
  implicit none
  type(constra) :: cts
  integer :: ich
  integer :: i
!---- do it only if constaints are given
  if (cts%NCI .and. allocated(cts%pots)) then
    do i = 1,10
      if (trim(cts%pots(i)) .ne. '') then
        write (ich,'(a)') trim(cts%pots(i))
      end if
    end do
    return
  end if
  return
end subroutine write_cts_NCI

!-----------------------------------------------------------------------------
! append NCI constraints to opened file
!-----------------------------------------------------------------------------
subroutine write_cts_NCI_pr(ich,cts)
  use crest_data
  implicit none
  type(constra) :: cts
  integer :: ich
  integer :: i
!---- do it only if constaints are given
  if (cts%NCI .and. allocated(cts%pots)) then
    do i = 1,10
      if (trim(cts%pots(i)) .ne. '') then
        write (ich,'(a,a)') '> ',trim(cts%pots(i))
      end if
    end do
  end if
  return
end subroutine write_cts_NCI_pr

!-----------------------------------------------------------------------------
! append NCI constraints to opened file
!-----------------------------------------------------------------------------
subroutine write_cts_biasext(ich,cts)
  use crest_data
  use iomod
  implicit none
  type(constra) :: cts
  integer :: ich
!---- do it only if constaints are given
  if (cts%usermsdpot) then
    !l = sylnk(cts%rmsdpotfile,'rmsdpot.xyz')
    write (ich,'(a)') '$metadyn'
!        write(ich,'(2x,a,a)') 'bias_input=rmsdpot.xyz'
    write (ich,'(2x,a,a)') 'bias_input=',cts%rmsdpotfile
!       write(ich,'(2x,a,a)') 'bias atoms: ',cts%biasatoms
    write (ich,'(a)') '$end'
    return
  end if
  return
end subroutine write_cts_biasext

!-----------------------------------------------------------------------------
! build a constrainment file for the chosen list of atoms
!-----------------------------------------------------------------------------
subroutine quick_constrain_file(fname,nat,at,atlist)
  use iso_fortran_env,only:output_unit
  use iomod
  implicit none
  !> Input
  integer,intent(in) :: nat
  integer,intent(in) :: at(nat)
  character(len=*) :: atlist
  character(len=*) :: fname
  !> Local variables
  integer,allocatable :: unconstrained(:)
  integer :: ncon

  write (*,'(1x,a,a)') 'Input list of atoms: ',trim(atlist)

  allocate (unconstrained(nat))
  call parse_atlist_new(atlist,ncon,nat,at,unconstrained)
  !> "unconstrained" contains all the *constrained* atoms after parsing
  !> which has to be *inversed* for the next subroutine, therefore:
  unconstrained = unconstrained*(-1)
  unconstrained = unconstrained + 1

  call build_constrain_file(fname,nat,unconstrained)
  write (*,'(1x,i0,a,i0,a)') ncon,' of ',nat,' atoms will be constrained.'
  write (*,'(1x,a,a,a)') 'A reference coord file ',fname,'.ref was created.'
  write (*,'(1x,a,/)') 'The following will be written to <.xcontrol.sample>:'
  call cat_mod(output_unit,' > ','.xcontrol.sample','')
  write (*,*)
  deallocate (unconstrained)
  stop '<.xcontrol.sample> written. exit.'
end subroutine quick_constrain_file

subroutine build_constrain_file(fname,nat,unconstrained)
  use iomod
  implicit none
  integer :: nat
  integer :: unconstrained(nat)  !UNconstrained atoms have value 1, else 0
  character(len=*) :: fname
  character(len=256) :: atstr2
  integer :: ich

  character(len=:),allocatable :: fref

  fref = fname//'.ref'
  call copy(fname,fref)

  open (newunit=ich,file='.xcontrol.sample')
  write (ich,'(a)') '$constrain'
  atstr2 = ''
  call build_atlist(atstr2,nat,unconstrained,.true.)

  write (ich,'(2x,a,a)') 'atoms: ',trim(atstr2)
  write (ich,'(2x,a)') 'force constant=0.5'
  write (ich,'(2x,a,a)') 'reference=',trim(fref)
  write (ich,'(a)') '$metadyn'

  atstr2 = ''
  call build_atlist(atstr2,nat,unconstrained,.false.)

  write (ich,'(2x,a,a)') 'atoms: ',trim(atstr2)
  write (ich,'(a)') '$end'
  close (ich)

  return
end subroutine build_constrain_file

!-----------------------------------------------------------------------------
! Build constraint for first X atoms only
! (e.g. only the heavy atoms in the tautomerization runtype)
!-----------------------------------------------------------------------------
subroutine fix_first_X_atoms(X,forceconst,oname)
  use iso_fortran_env,wp => real64
  implicit none
  integer :: X
  real(wp) :: forceconst
  character(len=*) :: oname
  character(len=:),allocatable :: tmp
  character(len=20) :: dumm
  integer :: ich
  if (oname == '') then
    tmp = '.tmpfix'
  else
    tmp = trim(oname)
  end if
  open (newunit=ich,file=tmp)
  write (ich,'(a)') '$constrain'
  write (ich,'(3x,a,i0)') 'atoms: 1-',X
  write (dumm,'(f16.4)') forceconst
  write (ich,'(3x,a,a)') 'force constant=',adjustl(trim(dumm))
  write (ich,'(a)') '$end'
  close (ich)
  return
end subroutine fix_first_X_atoms

!-----------------------------------------------------------------------------
! build a  atom list in condensed form
!-----------------------------------------------------------------------------
subroutine build_atlist(atstr,nat,atlist,reverse)
  implicit none
  character(len=*) :: atstr
  integer :: nat
  integer :: atlist(nat)
  logical :: reverse
  integer :: i,j,k,l
  character(len=10) :: dum
  integer :: rev

  if (reverse) then
    rev = 0
  else
    rev = 1
  end if
  atstr = ''
  k = 0
  iloop: do i = 1,nat
    if (k .gt. 0) then
      k = k - 1
      cycle iloop
    end if
    if (atlist(i) .eq. rev) then
      write (dum,'(i0)') i
      atstr = trim(atstr)//trim(dum)
      jloop: do j = i + 1,nat
        if (atlist(j) .eq. rev) then
          k = k + 1
        else
          exit jloop
        end if
      end do jloop
      if (k .gt. 1) then
        write (dum,'(i0)') i + k
        atstr = trim(atstr)//'-'//trim(dum)
      else
        k = 0
      end if
      atstr = trim(atstr)//','
    end if
  end do iloop
  l = len_trim(atstr)
  if (atstr(l:l) == ',') atstr(l:l) = ' '
  return
end subroutine build_atlist

!-----------------------------------------------------------------------------
! A control option for static Metadynamics: select only heavy atoms
!-----------------------------------------------------------------------------
subroutine mtdatoms(env)
  use crest_data
  use strucrd
  use zdata
  type(systemdata) :: env
  type(coord) :: mol
  type(zmolecule) :: zmol
!  character(len=*) :: filname
  integer :: i,j
  integer,allocatable :: inc(:)
  character(len=256) :: atstr
  integer :: r,rs
!  call mol%open(trim(filname))
  call env%ref%to(mol)
  call simpletopo(mol%nat,mol%at,mol%xyz,zmol,.false.,.true.,'')
  allocate (inc(env%nat),source=0)
  !-- exclude H atoms
  do i = 1,env%nat
    if (mol%at(i) .ne. 1) then  !heavy atoms
      inc(i) = 1
    end if
    if (env%includeRMSD(i) .ne. 1) then
      inc(i) = 0
    end if
  end do
  !-- exclude rings
  if (zmol%nri .gt. 0) then
    do r = 1,zmol%nri
      rs = zmol%zri(r)%rs
      if (rs <= env%emtd%rmax) then
        do i = 1,rs
          j = zmol%zri(r)%rlist(i)
          inc(j) = 0
        end do
      end if
    end do
  end if
  if (sum(inc) > 0) then
    call build_atlist(atstr,env%nat,inc,.false.)
    env%emtd%atomlist = trim(atstr)
    env%emtd%katoms = sum(inc)
  else
    env%emtd%atomlist = ''
    env%emtd%katoms = mol%nat
  end if
  if(allocated(env%emtd%atomlist2)) deallocate(env%emtd%atomlist2)
  allocate(env%emtd%atomlist2(mol%nat), source =.false.)
  do i=1,mol%nat
    env%emtd%atomlist2(i) = inc(i) > 0
  enddo 
  deallocate (inc)
  call zmol%deallocate()
  call mol%deallocate()
  return
end subroutine mtdatoms

!-----------------------------------------------------------------------------
! read file with REACTOR-settings into memory
!-----------------------------------------------------------------------------
subroutine rdrcontrol(fname,env)
  use crest_data
  use filemod
  implicit none

  type(systemdata) :: env
  character(len=*) :: fname

  type(filetype) :: rc
  integer :: i,n

  if (allocated(env%cts%rctrl)) then
    deallocate (env%cts%rctrl)
  end if

  call rc%open(fname)

  n = rc%nlines
  if (n > 0) then
    allocate (env%cts%rctrl(n))
    env%cts%nrctrl = n

    do i = 1,n
      env%cts%rctrl(i) = trim(rc%line(i))
    end do
    env%cts%ureactor = .true.
  end if
  call rc%close()
  return
end subroutine rdrcontrol
!-----------------------------------------------------------------------------
! append reactor constraints to opened file
!-----------------------------------------------------------------------------
subroutine write_cts_rcontrol(ich,cts)
  use crest_data
  implicit none
  type(constra) :: cts
  integer :: ich,i
  if (cts%ureactor) then
    do i = 1,cts%nrctrl
      write (ich,'(a)') trim(cts%rctrl(i))
    end do
  end if
  return
end subroutine write_cts_rcontrol

!-----------------------------------------------------------------------------
! read file with BONDLENGTH-settings into memory
!-----------------------------------------------------------------------------
subroutine rd_cbonds(fname,env)
  use crest_data
  use filemod
  implicit none

  type(systemdata) :: env
  character(len=*) :: fname

  type(filetype) :: rc
  integer :: i,n

  if (allocated(env%cts%cbonds)) then
    deallocate (env%cts%cbonds)
  end if

  call rc%open(fname)

  n = rc%nlines
  if (n > 0) then
    allocate (env%cts%cbonds(n))
    env%cts%n_cbonds = n
    do i = 1,n
      env%cts%cbonds(i) = trim(rc%line(i))
    end do
  end if
  call rc%close()
  return
end subroutine rd_cbonds

!-----------------------------------------------------------------------------
! append BONDLENGTH constraints to openend file
!-----------------------------------------------------------------------------
subroutine write_cts_CBONDS(ich,cts)
  use crest_data
  implicit none
  type(constra) :: cts
  integer :: ich
  integer :: i
!---- do it only if constaints are given
  if (allocated(cts%cbonds)) then
    do i = 1,cts%n_cbonds
      if (trim(cts%cbonds(i)) .ne. '') then
        write (ich,'(a)') trim(cts%cbonds(i))
      end if
    end do
  end if
  return
end subroutine write_cts_CBONDS

!-----------------------------------------------------------------------------
! append other settings to opened file
!-----------------------------------------------------------------------------
subroutine write_cts_DISP(ich,cts)
  use crest_data
  implicit none
  type(constra) :: cts
  integer :: ich
  character(len=40) :: dum
!---- apply dispersion scaling factor (> xtb 6.4.0)
  write (ich,'(a)') '$gfn'
  write (dum,'(f16.6)') cts%dscal
  write (ich,'(2x,a,a)') 'dispscale=',trim(adjustl(dum))
  return
end subroutine write_cts_DISP

!========================================================================================!
!> parse a list of atoms to exclude from cregen topology check
!========================================================================================!
subroutine parse_topo_excl(env,arg)
  use iso_fortran_env,only:wp => real64
  use crest_data
  use strucrd,only:i2e
  implicit none
  type(systemdata) :: env
  character(len=*) :: arg

  integer :: i,j,k,nat
  integer :: nselect
  integer,allocatable :: natlist(:)
  integer,allocatable :: at(:)
  logical,allocatable :: atypes(:)
  character(len=:),allocatable :: atms

  nat = env%nat
  allocate (natlist(nat),source=0)
  call parse_atlist(arg,nselect,nat,natlist)

  if (nselect > 0) then
    if (allocated(env%excludeTOPO)) deallocate (env%excludeTOPO)
    allocate (env%excludeTOPO(nat),source=.false.)
    do i = 1,nat
      if (natlist(i) > 0) then
        env%excludeTOPO(i) = .true.
        !write(*,*) 'excluding atom',i, ' from topology check'
      else
        env%excludeTOPO(i) = .false.
      end if
    end do
  end if
  deallocate (natlist)

  allocate (at(nat), source=0)
  allocate (atypes(118), source=.false.)
  if (allocated(env%ref%at)) then
    at = env%ref%at
    atypes = .false.
    call parse_atypelist(arg,nselect,atypes)
    if (nselect > 0) then
      if (.not. allocated(env%excludeTOPO)) allocate (env%excludeTOPO(nat),source=.false.)
      do i = 1,118
        if (atypes(i)) then
          do j = 1,nat
            if (at(j) == i) env%excludeTOPO(j) = .true.
          end do
        end if
      end do
    end if
  end if

  !>--- printout construction
  if (allocated(env%excludeTOPO)) then
    nselect = 0
    atypes = .false.
    do i = 1,nat
      if (env%excludeTOPO(i)) then
        nselect = nselect + 1
        if (at(i).ne.0) atypes(at(i)) = .true.
      end if
    end do
    atms = '('
    k = 0
    do i = 1,118
      if (atypes(i)) then
        k = k + 1
        if (k > 1) atms = atms//','
        atms = atms//trim(i2e(i,'nc'))
      end if
    end do
    atms = atms//')'

    write (*,'(2x,a,i0,a,a)') '-notopo <x> : ignoring topology on ',nselect,' atoms ',atms
  end if

  if (allocated(at)) deallocate (at)
  if (allocated(atypes)) deallocate (atypes)
  return
end subroutine parse_topo_excl

!========================================================================================!
!> parse a list of atoms (either by their position in the input file, or by atom type)
!========================================================================================!
subroutine parse_atlist_new(arg,selected,nat,at,natlist)
  use strucrd
  implicit none
  !> Input
  character(len=*) :: arg
  integer,intent(in) :: nat
  integer,intent(in) :: at(nat)
  !> Output
  integer,intent(out) :: selected
  integer,intent(inout) :: natlist(nat)
  !> Local variables
  integer :: i,j
  integer :: nselect
  logical,allocatable :: atypes(:)
  logical,allocatable :: atlist(:)
  !>-- initialize
  natlist = 0

  !>-- first the specific atomlist
  allocate(atlist(nat),source=.false.) 
  call get_atlist(nat,atlist,arg,at) 
  nselect = count(atlist,1)

  !>-- then select atoms
  if (nselect > 0) then
    do j = 1,nat
       if (atlist(j)) natlist(j) = 1
     end do
  end if
  if (allocated(atlist)) deallocate (atlist)

  !>-- count how many have been selected
  selected = sum(natlist)

  return
end subroutine parse_atlist_new

