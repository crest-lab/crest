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

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  append the ".constrains" file to the setblock (ONLY USED IN OLD PARTS OF THE CODE)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine ConstrainsToSET(fname,isplainxyz,constraints)
      use setdgmod
      use iomod
      implicit none
      character(len=*) :: fname
      character(len=*) :: constraints
      character(len=512) :: atmp,dg
      character(len=256),allocatable :: strings(:)
      character(len=128) :: key
      real*8,allocatable  :: floats(:)
      integer :: i,j,k
      integer :: iost,cs,cf,ich
      logical :: ex,isplainxyz

      inquire(file=constraints,exist=ex)
      if(.not.ex)return

      if(isplainxyz)then
        call appendto(constraints,fname)
      else
        open(newunit=ich,file=constraints)
        allocate(floats(3),strings(3))
         outer: do
           read(ich,'(a)',iostat=iost)atmp
           if(iost < 0)exit
      !---- handle $set-blocks in .constrains
           if(index(atmp,'$set').ne.0)then
             setblock: do
               read(ich,'(a)',iostat=iost)atmp
               if(iost < 0)exit outer
               atmp=adjustl(atmp)
               if(index(atmp,'$end').ne.0)exit setblock
               dg=trim(atmp)
               call split_set_args(dg,atmp)
               call setdg(fname,dg,trim(atmp))
             enddo setblock
           endif
      !---- handle all other blocks in .constrains
           if((index(atmp,'$').ne.0).and. &
           &  ((index(atmp,'$set').eq.0).and.(index(atmp,'$end').eq.0)))then
               !write(*,*) trim(atmp)
               key=trim(atmp) !get the argument as the current reference block keyword
               call setdg_block(fname,trim(key),trim(atmp)) !set it to the coord file
             keyword: do
               read(ich,'(a)',iostat=iost)atmp
               if(iost < 0)exit outer
               if(index(atmp,'$').ne.0)then !next keyword found
                 backspace(ich)
                 cycle outer
               endif
               !write(*,*) trim(atmp)
               call setdg_block(fname,trim(key),trim(atmp))
             enddo keyword
           endif
        enddo outer

       deallocate(strings,floats)
       close(ich)
      endif
      return
end subroutine ConstrainsToSET

!-----------------------------------------------------------------------------
! sort the .constrains file
!-----------------------------------------------------------------------------
subroutine sort_constrain_file(fname)
         use iomod, only: to_lower
         implicit none
         character(len=256),allocatable :: entire(:),setsave(:)
         character(len=256) :: atmp
         character(len=*) :: fname
         integer :: i,j,k
         integer :: io,ich,och,lines
         integer :: setb,sete
         logical :: ex,follow,setset
        
         inquire(file=fname,exist=ex)
         if(.not.ex) return
         open(newunit=ich,file=fname)
         open(newunit=och,file='.tmpconstrains')
         lines=0
         do
            read(ich,'(a)',iostat=io)atmp
            if(io.lt.0) exit
            lines=lines+1
         enddo
         rewind(ich)
         allocate(setsave(lines),entire(lines))
         entire=''
         setsave=''
         do i=1,lines
            read(ich,'(a)',iostat=io)entire(i)
            call To_lower(entire(i))
            entire(i)=adjustl(entire(i))
            if(io.lt.0) exit
         enddo
!---- seperate everything that is written inside of a $set statement
         setb=1
         follow=.false.
         setset=.false.
         do i=1,lines
            if(index(entire(i),'$end').ne.0)then
               follow=.false.
               entire(i)=''  !"delete" from memory
            endif
            if((index(entire(i),'$').ne.0).and.follow)then
               follow=.false.
            endif
            if(follow)then
               setsave(setb)=entire(i)
               setb=setb+1
               entire(i)=''  !"delete" from memory
            endif
            if(index(entire(i),'$set').ne.0)then
               follow=.true.
               setset=.true.
               entire(i)=''  !"delete" from memory
            endif
         enddo
!---- write the new file in order, ignoring blank lines
!---- first everything not contained within $set-statements
         do i=1,lines
            if(len_trim(entire(i)).lt.1)cycle
            if(index(entire(i),'$').ne.0)then
               write(och,'(a)')trim(entire(i))
            else
               write(och,'(2x,a)')trim(entire(i))
            endif
         enddo
!---- then all the set stuff
         if(setset)then
          write(och,'(a)')'$set'
          do i=1,lines
            if(len_trim(setsave(i)).lt.1)cycle
            write(och,'(2x,a)')trim(setsave(i))
          enddo
         endif
!---- terminate all of it with a $end
         write(och,'(a)')'$end'
!---- deallocate and rename the file
         deallocate(entire,setsave)
         close(ich,status='delete')
         close(och)
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


        line=adjustl(line)

        do i=1,len(trim(line))
           digit=line(i:i)
           if((digit==' ').or.(digit==char(9)))then
           j=i+1
           argument=adjustl(line(j:))
           j=i-1
           dummy=trim(line(:j))
           exit
           endif
        enddo
        line=trim(dummy)

!        deallocate(strings,floats)

end subroutine split_set_args

!-----------------------------------------------------------------------------
! parse the atom list
!-----------------------------------------------------------------------------
subroutine parse_atlist(line,selected,nat,natlist)
        implicit none
        character(len=*) :: line
        integer,intent(out) :: selected
        
        integer :: i,j,k
        character(len=1) :: digit
        character(len=10),parameter :: numbers='0123456789'
        character(len=20) :: atom1,atom2
        integer :: at1,at2,nat

        integer,intent(inout) :: natlist(nat)
        logical :: rang,stor

        !allocate(natlist(nat))
        natlist=0
 
        stor=.false.
        rang=.false.

        atom1=''
        atom2=''
        at1=0
        at2=0
        selected=0

        do i=1,len(trim(line))
         digit=line(i:i)
         if(index(numbers,digit).ne.0)then  !should exclude tabstops and blanks, 9 is ascii code for tab
            if(.not.rang)then
            atom1=trim(atom1)//digit
            stor=.true.
              else
            atom2=trim(atom2)//digit
            endif
         endif
         if((digit==',').or.(digit==' ').or.(digit==char(9)) &
         &  .or.(i.eq.len(trim(line))))then  !comma,space,tab or last character
            if(.not.rang.and.stor)then !a single atom
               read(atom1,*)at1
               atom1=''
               if(at1.ge.1 .and. at1.le.nat)then
               natlist(at1)=1
               endif
               stor=.false.
            endif
            if(rang.and.stor)then !a range of atoms
              read(atom2,*)at2
              atom2=''
              stor=.false.
              rang=.false.
              if(at1.ge.1 .and. at2.ge.1)then
              if(min(at1,at2).le.nat)then
                k=max(at1,at2)
                do j=min(at1,at2),min(k,nat)   !order is arbitrary
                  natlist(j)=1
                enddo
              endif
              endif
            endif
         endif
         if((digit=='-').and.stor)then !begin to select a range of atoms
            rang=.true.
            read(atom1,*)at1
            atom1=''
         endif
       enddo
       !write(*,*) natlist
       selected=sum(natlist)
       return
end subroutine parse_atlist


!-----------------------------------------------------------------------------
! modified version of the rdcoord routine to only include selected atoms
!-----------------------------------------------------------------------------
subroutine rdcoord_reduced(fname,n,rn,xyz,iat,rlist)
      use iso_fortran_env, only: wp => real64
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
      j=0
      do i=1,n
         if(rlist(i)==1)then
           j=j+1
           xyz(1:3,j)=struc%xyz(1:3,i)
           iat(j)=struc%at(i)  
         endif
      enddo
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
     integer :: i,j,k
     integer :: bpos
     logical :: ex,follow
     
!--- some initialization issues
     cts%buff = ''
     bpos=1
     follow=.false.

!--- collect dublicate blocks
     do i=1,cts%ndim
         atmp=adjustl(cts%sett(i))
         if(index(atmp(1:1),'$').ne.0)then
            if(index(atmp,'$end').ne.0) cycle
            cts%buff(bpos)=atmp
            bpos=bpos+1
            follow=.true.
            do j=1,cts%ndim 
               btmp=adjustl(cts%sett(j))
               if(index(btmp(1:1),'$').ne.0)then          !--- if its a new block don't follow ...
                  follow=.false.
               endif
               if(index(trim(atmp),trim(btmp)).ne.0)then  !--- ... except if its the same keyword
                  follow=.true.
                  cycle
               endif
               if(follow)then
                 if(trim(btmp).eq.'')cycle                !--- cycle blanks
                 write(cts%buff(bpos),'(2x,a)')trim(btmp)
                 bpos=bpos+1
               endif
            enddo
         endif
     enddo
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
     integer :: i,j,k
     integer :: io,ich,och,lines
     integer :: ndim
     logical :: ex,follow,setset

!--- get the dimensions        
     inquire(file=fname,exist=ex)
     if(.not.ex) return
     open(newunit=ich,file=fname)
     lines=0
     do
        read(ich,'(a)',iostat=io)atmp
        if(io < 0) exit
        lines=lines+1
     enddo
     rewind(ich)
!--- allocate
     ndim=max(lines,10)
     call cts%allocate(ndim)

!--- read 
     do i=1,ndim
       read(ich,'(a)',iostat=io)atmp
       if(io < 0) exit
       atmp=adjustl(atmp)
       cts%sett(i)=atmp
     enddo

     close(ich)
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
     integer :: i,j,k
     character(len=40) :: dum
!--- not really a "constraint", but convenient for the implementation:
!    change of the gbsa grid
     if(cts%ggrid)then
     if(allocated(cts%gbsagrid))then
         write(ich,'(a)') '$gbsa'
         write(ich,'(2x,a,a)')'gbsagrid=',cts%gbsagrid
     endif
     endif
!---- do it only if constaints are given, write them
     if(cts%used)then
      do i=1,cts%ndim
         if(trim(cts%sett(i)).ne.'')then
            write(ich,'(a)') trim(cts%sett(i))
         endif
      enddo
     endif
!---- apply bondlength constraints globally
     if(cts%cbonds_global)then
      call write_cts_CBONDS(ich,cts)  
     endif
     if(cts%dispscal_global)then
      call write_cts_DISP(ich,cts)
     endif
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
     integer :: i,j,k
!---- do it only if constaints are given
     if(cts%NCI)then
      do i=1,10
         if(trim(cts%pots(i)).ne.'')then
            write(ich,'(a)') trim(cts%pots(i))
         endif
      enddo
       return
     endif
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
     integer :: i,j,k
!---- do it only if constaints are given
     if(cts%NCI)then
      do i=1,10
       if(trim(cts%pots(i)).ne.'')then
         write(ich,'(a,a)') '> ',trim(cts%pots(i))
       endif
      enddo
     endif
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
     integer :: i,j,k,l
!---- do it only if constaints are given
     if(cts%usermsdpot)then
        !l = sylnk(cts%rmsdpotfile,'rmsdpot.xyz')
        write(ich,'(a)')'$metadyn'
!        write(ich,'(2x,a,a)') 'bias_input=rmsdpot.xyz'
        write(ich,'(2x,a,a)') 'bias_input=',cts%rmsdpotfile
!       write(ich,'(2x,a,a)') 'bias atoms: ',cts%biasatoms        
        write(ich,'(a)')'$end'
       return
     endif
     return
end subroutine write_cts_biasext



!-----------------------------------------------------------------------------
! build a constrainment file for the chosen list of atoms
!-----------------------------------------------------------------------------
subroutine quick_constrain_file(fname,nat,atlist)
     use iso_fortran_env, only: output_unit
     use iomod
     implicit none
     integer :: nat
     character(len=*) :: atlist
     integer,allocatable :: unconstrained(:)
     character(len=*) :: fname
     integer :: ncon,i,j

     write(*,'(1x,a,a)') 'Input list of atoms: ',trim(atlist)

     allocate(unconstrained(nat))
     call parse_atlist(atlist,ncon,nat,unconstrained) !"unconstrained" contains all the constrained atoms
     unconstrained=unconstrained*-1                   !which has to be inversed for the next subroutine 
     unconstrained=unconstrained+1     
     !write(*,*) unconstrained
     call build_constrain_file(fname,nat,unconstrained) 
     write(*,'(1x,i0,a,i0,a)') ncon,' of ',nat,' atoms will be constrained.'
     write(*,'(1x,a,a,a)')'A reference coord file ',fname,'.ref was created.'
     write(*,'(1x,a,/)')'The following will be written to <.xcontrol.sample>:'
     call cat_mod(output_unit,' > ','.xcontrol.sample','')
     write(*,*)
     deallocate(unconstrained)
     stop '<.xcontrol.sample> written. exit.'
end subroutine quick_constrain_file


subroutine build_constrain_file(fname,nat,unconstrained)
     use iomod
     implicit none
     integer :: nat
     integer :: unconstrained(nat)  !UNconstrained atoms have value 1, else 0
     character(len=*) :: fname
     character(len=:),allocatable :: atstr
     character(len=256) :: atstr2
     character(len=10) :: dum
     integer :: ich,i,j,k,l
  
     character(len=:),allocatable :: fref
     
     fref=fname//'.ref'
     call copy(fname,fref)
     
     open(newunit=ich,file='.xcontrol.sample')
     write(ich,'(a)')'$constrain'
     atstr2=''
     call build_atlist(atstr2,nat,unconstrained,.true.)

     write(ich,'(2x,a,a)')'atoms: ',trim(atstr2)
     write(ich,'(2x,a)')'force constant=0.5'
     write(ich,'(2x,a,a)')'reference=',trim(fref)
     write(ich,'(a)')'$metadyn'

     atstr2=''
     call build_atlist(atstr2,nat,unconstrained,.false.)

     write(ich,'(2x,a,a)')'atoms: ',trim(atstr2)
     write(ich,'(a)')'$end'
     close(ich)

     return
end subroutine build_constrain_file

!-----------------------------------------------------------------------------
! Build constraint for first X atoms only 
! (e.g. only the heavy atoms in the tautomerization runtype)
!-----------------------------------------------------------------------------
subroutine fix_first_X_atoms(X,forceconst,oname)
    use iso_fortran_env, wp => real64
    implicit none
    integer :: X
    real(wp) :: forceconst
    character(len=*) :: oname
    character(len=:),allocatable :: tmp
    character(len=20) :: dumm
    integer :: ich
    if(oname=='')then
        tmp='.tmpfix'
    else
        tmp=trim(oname)
    endif
    open(newunit=ich,file=tmp)
     write(ich,'(a)')'$constrain'
     write(ich,'(3x,a,i0)')'atoms: 1-',X
     write(dumm,'(f16.4)') forceconst
     write(ich,'(3x,a,a)')'force constant=',adjustl(trim(dumm))
     write(ich,'(a)')'$end'
    close(ich)
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
       
     if(reverse)then
        rev = 0
     else
        rev = 1
     endif
     atstr=''
     k=0    
     iloop : do i=1,nat
       if(k.gt.0)then
          k=k-1
          cycle iloop
       endif
       if(atlist(i).eq.rev)then
         write(dum,'(i0)')i
         atstr=trim(atstr)//trim(dum)
         jloop: do j=i+1,nat
            if(atlist(j).eq.rev)then
               k=k+1
            else
              exit jloop
            endif
         enddo jloop
         if(k.gt.1)then
           write(dum,'(i0)') i+k
           atstr=trim(atstr)//'-'//trim(dum)
         else
           k=0
         endif
         atstr=trim(atstr)//','
       endif
     enddo iloop
     l=len_trim(atstr)
     if(atstr(l:l)==',') atstr(l:l)=' '
     return
end subroutine build_atlist

!-----------------------------------------------------------------------------
! A control option for static Metadynamics: select only heavy atoms
!-----------------------------------------------------------------------------
subroutine mtdatoms(filname,env)
    use crest_data
    use strucrd
    use zdata
    !type(options) :: opt
    type(systemdata) :: env
    type(coord) :: mol
    type(zmolecule) :: zmol
    character(len=*) :: filname
    integer :: i,j,k,l
    integer,allocatable :: inc(:)
    character(len=256) :: atstr
    integer :: r,rs
    call mol%open(trim(filname))
    call simpletopo(mol%nat,mol%at,mol%xyz,zmol,.false.,'')
    allocate(inc(env%nat), source=0)
    !-- exclude H atoms
    do i=1,env%nat
      if(mol%at(i) .ne. 1)then  !heavy atoms
          inc(i) = 1
      endif
      if(env%includeRMSD(i) .ne. 1)then
          inc(i) = 0
      endif
    enddo
    !-- exclude rings
    if(zmol%nri .gt. 0)then
     do r=1,zmol%nri
        rs = zmol%zri(r)%rs 
        if(rs <= env%emtd%rmax)then
         !write(*,*)'ring',zmol%zri(r)%rlist(1:rs)   
         do i=1,rs
            j =  zmol%zri(r)%rlist(i)
            inc(j) = 0
         enddo 
        endif
     enddo
    endif
    !write(*,*) inc
    if(sum(inc)>0)then
    call build_atlist(atstr,env%nat,inc,.false.)
    env%emtd%atomlist=trim(atstr)
    env%emtd%katoms = sum(inc)
    else
    env%emtd%atomlist=''
    env%emtd%katoms = mol%nat
    endif
    deallocate(inc)
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

    !type(options) :: opt
    type(systemdata) :: env
    character(len=*) :: fname

    type(filetype) :: rc
    integer :: i,n

    if(allocated(env%cts%rctrl))then
        deallocate(env%cts%rctrl)
    endif

    call rc%open(fname)

    n = rc%nlines
    if(n > 0)then
      allocate(env%cts%rctrl(n))
      env%cts%nrctrl = n

      do i=1,n
         env%cts%rctrl(i) = trim(rc%line(i))
      enddo
      env%cts%ureactor = .true.
    endif
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
     if(cts%ureactor)then
      do i=1,cts%nrctrl
          write(ich,'(a)') trim(cts%rctrl(i))
      enddo
     endif
     return
end subroutine write_cts_rcontrol

!-----------------------------------------------------------------------------
! read file with BONDLENGTH-settings into memory
!-----------------------------------------------------------------------------
subroutine rd_cbonds(fname,env)
    use crest_data
    use filemod
    implicit none

    !type(options) :: opt
    type(systemdata) :: env
    character(len=*) :: fname

    type(filetype) :: rc
    integer :: i,n

    if(allocated(env%cts%cbonds))then
        deallocate(env%cts%cbonds)
    endif

    call rc%open(fname)

    n = rc%nlines
    if(n > 0)then
      allocate(env%cts%cbonds(n))
      env%cts%n_cbonds = n
      do i=1,n
         env%cts%cbonds(i) = trim(rc%line(i))
      enddo
    endif
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
     integer :: i,j,k
!---- do it only if constaints are given
     if(allocated(cts%cbonds))then
      do i=1,cts%n_cbonds
         if(trim(cts%cbonds(i)).ne.'')then
            write(ich,'(a)') trim(cts%cbonds(i))
         endif
      enddo
     endif 
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
     integer :: i,j,k
     character(len=40) :: dum
!---- apply dispersion scaling factor (> xtb 6.4.0)
     write(ich,'(a)') '$gfn'
     write(dum,'(f16.6)')cts%dscal
     write(ich,'(2x,a,a)')'dispscale=',trim(adjustl(dum))
     return
end subroutine write_cts_DISP

