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

! Nov 2017, PP 
! Updated Nov 2018
! Updated May 2019
! Updated Jan 2021

!! --------------------------------------------------------------------------------------
!  Sort all H atoms in the coord file to the bottom
!! --------------------------------------------------------------------------------------
subroutine htothebottom(fname,ichrg,n,atmap)
       use crest_parameters
       use strucrd, only: rdnat,rdcoord,i2e
       implicit none

       integer :: i,j
       character(len=*) :: fname
       integer :: n,ichrg
       integer :: atmap(n)
       integer,allocatable :: at(:),atnew(:)
       real(wp),allocatable :: xyz0(:,:)
       integer :: ich,m

       call rdnat(fname,n)
       allocate(xyz0(3,n),at(n),atnew(n))
       call rdcoord(fname,n,at,xyz0)

       m=0
       open (newunit=ich,file='coord')
       write(ich,'(a)')'$coord'
       do i=1,n
          if(at(i).ne.1)then                                     !first all non-H-atoms
             write(ich,'(3F24.14,6x,a2)')xyz0(1:3,i),i2e(at(i))
             m=m+1
             atmap(m) = i
             atnew(m) = at(i)
          endif
       enddo
       do j=1,n
           if(at(j).eq.1)then                                    !then all H-atoms
             write(ich,'(3F24.14,6x,a2)')xyz0(1:3,j),i2e(at(j))
             m=m+1
             atmap(m) = j
             atnew(m) = at(j)
          endif
       enddo
       write(ich,'(a)')'$end'
       write(ich,'(1x,''$chrg'',1x,i0)')ichrg
       write(ich,'(a)')'$end'
       close(ich)

       if(.not.(all(at==atnew)))then
       write(*,*)'Input coordinate lines sorted:'
       write(*,'(1x,a,1x,a5,1x,a5)') 'element','old','new'
       call print_map(6,n,atnew,atmap)
       endif

       deallocate(at,xyz0)
       return
end subroutine htothebottom

subroutine print_map(prch,nat,at,atmap)
    use crest_data
    use strucrd, only: i2e
    implicit none
    integer :: prch
    integer :: nat
    integer :: at(nat)
    integer :: atmap(nat)
    integer :: i
    do i=1,nat
    write(prch,'(1x,a7,1x,i5,1x,i5)') i2e(at(i),'nc'),atmap(i),i
    enddo
    return
end subroutine print_map

!! --------------------------------------------------------------------------------------
!  Sort out all topologically equivalent structures (i.e. conformers)
!! --------------------------------------------------------------------------------------
subroutine cosort(iname,oname,wrscoord,verbose)
      use crest_parameters
      use iomod
      use strucrd, only: wrc0,rdensembleparam,rdensemble,wrxyz
      implicit none

      character(len=*),intent(in) :: iname
      character(len=*),intent(in) :: oname
      logical,intent(in)          :: wrscoord

      integer :: j,k,q,p,r
      real(wp),allocatable :: eread(:),xyz(:,:,:),xyztmp(:,:)
      real(wp),allocatable :: cn(:),bond(:,:)
      integer,allocatable :: at(:),group(:)
      real(wp) :: dE
      character(len=10),allocatable :: itensr(:),itensl(:)   !identifier tensor
      character(len=80) :: str
      integer :: n,nall,nonh,gc,sgc,tgc

      logical :: verbose
      integer :: ochan,ich

      write(*,*)
      write(*,*)'==================================================='
      write(*,'(a)')' Identifying topologically equivalent structures:'

      call rdensembleparam(iname,n,nall)
      allocate(xyz(3,n,nall),eread(nall),group(nall),at(n),xyztmp(3,n),cn(n),bond(n,n))
      call rdensemble(iname,n,nall,at,xyz,eread)

      call countnonh(n,at,nonh)
      !allocate(itens(nonh,nall))
      allocate(itensr(nonh),itensl(nonh))

!---- identifier sorting loops
      gc=1
      group=0
   
      do r=1,nall
        if(group(r).ne.0)cycle
        xyztmp(:,:)=xyz(:,:,r)/bohr
        call get_itens(n,xyztmp,at,nonh,itensr)
        do q=1,nall
          if(r.eq.q)then
            group(r)=gc
            tgc=1
            cycle
          endif
          xyztmp(:,:)=xyz(:,:,q)/bohr
          call get_itens(n,xyztmp,at,nonh,itensl)
          sgc=0
          do p=1,nonh
             !if(itens(p,r)==itens(p,q))then
             if(itensr(p)==itensl(p))then
               sgc=sgc+1
             else
               exit
            endif
          enddo
          if(sgc.eq.nonh)then
             group(q)=gc
             tgc=tgc+1
          endif
        enddo
        if(tgc.gt.1)then
        write(*,'(a,i0,a,i0,a)')' Equivalent to ',r,'. structure: ' &
        &                        ,tgc,' structure(s).'
        endif
        gc=gc+1        
      enddo
      write(*,'(a)')' Done.'
      write(*,'(a,a,a)')' Appending file <',trim(oname),'> with structures.'
      write(*,*)

      open(newunit=ochan,file='tmp')

      write(ochan,*)'==================================================='
      write(ochan,*)'============= ordered structure list =============='
      write(ochan,*)'==================================================='
      if(wrscoord)then
      write(ochan,'(a,a)')' written to <scoord.*> and ',trim(oname)
      write(ochan,*)
      write(ochan,'('' scoord.*     ΔE(kcal/mol)   Etot(Eh)'')')
      else
      write(ochan,'(a,a,a)')' written to file <',trim(oname),'>'
      write(ochan,*)
      write(ochan,'('' structure    ΔE(kcal/mol)   Etot(Eh)'')')
      endif

      open(newunit=ich,file=oname)
      k=0
      do j=1,gc-1
         do r=1,nall
            if(group(r).eq.j)then
              k=k+1
              dE=(eread(r)-eread(1))*627.5095_wp
              write(ochan,'(i5,6x,F10.2,4x,F14.6)') &
              & k,dE,eread(r)
              call wrxyz(ich,n,at,xyz(:,:,r),eread(r))
              where(group.eq.j)group=0
              if(wrscoord)then  ! write new scoord.* if necessary
                write(str,'(''scoord.'',i0)')k
                call wrc0(trim(str),n,at,xyz(:,:,r))
              endif
            else
              cycle
            endif
         enddo
      enddo
      close(ich)
      close(ochan)

      gc=gc-1

      if(nall.ne.gc)then
        write(*,'(a,i0,a,a,a)')' Initial ',nall,' structures from file ',trim(iname),' have'
        write(*,'(a,i0,a)')' been reduced to ',gc,' topologically unique structures.'
      else
        write(*,'(a,i0,a,a,a)')' All initial ',nall,' structures from file ',trim(iname),' are unique.'
      endif

      !deallocate(itens)
      deallocate(itensl,itensr)
      deallocate(bond,cn,xyztmp,at,eread,xyz)
      
      if(verbose) call cat('tmp')
      call remove('tmp')
      !call remove('identify')

end subroutine cosort

!! --------------------------------------------------------------------------------------
!! --------------------------------------------------------------------------------------

subroutine countnonh(nat,iz,n)    !count all non-hydrogen atoms
      implicit none
      integer :: nat,iz(nat),n,i
      n=0
      do i=1,nat
         if(iz(i).eq.1)then
            cycle
         else
            n=n+1
         endif
      enddo
end subroutine countnonh

subroutine counth(nat,iz,n)    !count all hydrogen atoms
      implicit none
      integer :: nat,iz(nat),n,i
      n=0
      do i=1,nat
         if(iz(i).eq.1)then
            cycle
         else
            n=n+1
         endif
      enddo
end subroutine counth


subroutine countnonh2(icn,idarr,n2)
      implicit none
      integer :: icn,idarr(icn),n2,i
      n2=0
      do i=1,icn
         if(idarr(i).eq.1)then
            cycle
         else
            n2=n2+1
         endif
      enddo
end subroutine countnonh2

!! --------------------------------------------------------------------------------------
!  build an identifier tensor for a single molecule
!! --------------------------------------------------------------------------------------
subroutine get_itens(n,xyz,at,ni,itens)
      use crest_parameters
      implicit none
      integer,intent(in)  :: n                    ! number of atoms
      integer,intent(in)  :: ni                   ! number of elements in the identifier
      real(wp),intent(in) :: xyz(3,n)             ! coordinates in Bohrs
      integer,intent(in)  :: at(n)                ! integer atom types
      character(len=10),intent(out) :: itens(ni)  ! Identifier
      integer :: ntopo
      integer,allocatable :: topo(:)
      character(len=10) :: ident
      integer :: p
      itens=''
      ntopo = n*(n+1)/2
      allocate(topo(ntopo))
      call quicktopo(n,at,xyz,ntopo,topo)
      do p=1,ni
        call analyze_neighbours(p,xyz,at,n,ntopo,topo,ident)
        itens(p)=ident
      enddo
      deallocate(topo)
      return
end subroutine get_itens

!! --------------------------------------------------------------------------------------
!  modified version of the neighbours subroutine in select.f
!! --------------------------------------------------------------------------------------
subroutine analyze_neighbours(i,xyz,iat,nat,ntopo,topo,ident)
      use crest_parameters
      implicit none
      integer,intent(in)  :: i
      real(wp),intent(in) :: xyz(3,nat)
      integer,intent(in)  :: iat(nat)
      integer,intent(in)  :: nat
      integer,intent(in) :: ntopo
      integer,intent(in) :: topo(ntopo)
      character(len=10),intent(inout) :: ident
      integer  :: j,k,l,icn,n2
      integer :: lin !this is a function
      integer,allocatable :: idarr(:),neighb(:)
      character(len=10) :: str
      character(len=9) :: str2     
      character(len=1) :: chiral
      logical :: chiralC,chiralX
      chiralC=.false.
      chiralX=.false.
      l=0
      ident=''
      icn=0
      do k=1,nat
         j=lin(i,k)
         if(topo(j)==1)icn=icn+1
      enddo
      allocate(idarr(icn),neighb(icn))                           !identifier array contains atoms
      l=0
      do k=1,nat
         j=lin(i,k)
         if(topo(j)==1)then
             l=l+1
             neighb(l)=k
             idarr(l)=iat(k)
         endif
      enddo

      call idwrite(icn,idarr,str)                   !write the identifier sorted by atom number (hydrogen first, then heavy atoms)
      ident=trim(str)                               !i.e., idarr(1) should be a hydrogen if there are hydrogen bound to atom i
      
      !Check if (pseudo-)chirality has to be determined for atom i
      if((iat(i).eq.6).and.(icn.eq.4))then     !for carbon
         chiralC=.true.
         call countnonh2(icn,idarr,n2)
         !if((idarr(1).eq.idarr(2)).and.(iat(idarr(1)).eq.1))then !methyl and ehtyl group workaround
         !   chiralC=.false.
         if(n2.le.2)then                                          !methyl and ehtyl group workaround
            chiralC=.false.
         endif
      endif
      if((iat(i).ne.6).and.(icn.ge.3))then     !for other heavy atoms
         chiralX=.true.
         call countnonh2(icn,idarr,n2)
         if((icn-n2).gt.1)then                     !if more than one hydrogen are present don't use the chirality identifier
            chiralX=.false.                        !because the results might differ (hydrogen order not fixed in the coordinates)
         endif
         if((icn.eq.3).and.(iat(i).eq.7))chiralX=.false. !special case nitrogen because of the low inversion barrier
         !write(*,*) icn,iat(i)
      endif
 
      if(chiralC.or.chiralX)then                    !add a (pseudo-)chirality identifier (+ or -)
         chiral=''
         call chispat(i,icn,neighb,nat,xyz,chiral)
         str2=trim(ident)
         write(ident,'(a,a)')trim(str2),trim(chiral)
         chiralC=.false.
         chiralX=.false.
      endif
 
      deallocate(neighb,idarr)
end subroutine analyze_neighbours

!! --------------------------------------------------------------------------------------
!  idwrite writing ordered identifier
!! --------------------------------------------------------------------------------------
subroutine idwrite(n,arr,str)
      use strucrd, only: i2e
      implicit none
      integer :: n,arr(n)
      integer :: i
      character(len=*) :: str
      character(len=6) :: dummy
      call quicksort(n,arr)
      str=''
      do i=1,n
         dummy=trim(str)
         write(str,'(a,a)')trim(dummy),trim(i2e(arr(i)))
      enddo      
end subroutine idwrite

!! --------------------------------------------------------------------------------------
!  chispat -  chirality identifier by using a spat product
!! --------------------------------------------------------------------------------------
subroutine chispat(c,icn,neighb,nat,xyz,chiral)
      use crest_parameters
      implicit none
      integer  :: c
      integer  :: nat
      integer  :: icn
      integer  :: neighb(icn)
      real(wp) :: xyz(3,nat)
      character(len=*) :: chiral

      integer :: i,j,k
      real(wp),allocatable :: xyz2(:,:)
      real(wp) :: mat(3,3)
      real(wp) :: sig

      allocate(xyz2(3,icn))

      chiral=''

      !build xyz2 according to atom order of the coordinates. These loops are important since the
      !matrix for the triple product always has to contain the same vectors in the same order.
      xyz2=0
      k=0
      do i=1,nat
         do j=1,icn
            if(i.eq.neighb(j))then
               k=k+1
               xyz2(:,k)=(xyz(:,neighb(j)))-xyz(:,c)         !move coordinate system according to atom c
               xyz2(:,k)=xyz2(:,k)/sqrt(sum(xyz2(1:3,k)**2))
             endif
         enddo
      enddo

      !calculate triple product to determine (pseudo-)chirality
      mat=0
      sig=0
      mat(:,1)=xyz2(:,1)
      mat(:,2)=xyz2(:,2)
      mat(:,3)=xyz2(:,3)
      if(icn.ge.4)then
         do j=4,icn
            mat(:,3)=mat(:,3)-xyz2(:,j)
         enddo
      endif
      call det3x3(mat,sig)
      !write(*,*)sig

      if(sig.lt.0)then
         chiral='-'
      else
         chiral='+'
      endif

      deallocate(xyz2)

end subroutine chispat

!! --------------------------------------------------------------------------------------

subroutine det3x3(m,detval)
      use crest_parameters, only: wp
      implicit none
      real(wp) :: detval
      real(wp) :: m(3,3)
      real(wp) :: vec(3)
      
      vec(1)=(m(2,1)*m(3,2))-(m(3,1)*m(2,2))
      vec(2)=(m(3,1)*m(1,2))-(m(1,1)*m(3,2))
      vec(3)=(m(1,1)*m(2,2))-(m(2,1)*m(1,2))

      detval=(vec(1)*m(1,3))+(vec(2)*m(2,3))+(vec(3)*m(3,3))

end subroutine det3x3

!! ------------------------------------------------------------------
!  identification of methyl groups
!! ------------------------------------------------------------------
subroutine methyl_autocomplete(n,xyz,at,equiv)
      use crest_parameters, only: wp
      implicit none
      integer,intent(in)  :: n                    ! number of atoms
      real(wp),intent(in) :: xyz(3,n)             ! coordinates
      integer,intent(in)  :: at(n)                ! integer atom types
      integer,allocatable :: eqv(:,:)
      integer,intent(inout) :: equiv(n+1,n)

      integer :: i,m
       allocate(eqv(n,3))
       call get_methyl(n,xyz,at,eqv)

       do i=1,n
          if(eqv(i,1).gt.0)then
             m=equiv(1,i)
             if(m.le.1)then
                equiv(1,i)=3
                equiv(2:4,i)=eqv(i,1:3)
             endif
          endif
       enddo
 
       deallocate(eqv)
end subroutine methyl_autocomplete


subroutine get_methyl(n,xyz,at,eqv)
      use crest_parameters, only: wp
      implicit none
      integer,intent(in)  :: n                    ! number of atoms
      real(wp),intent(in) :: xyz(3,n)             ! coordinates
      integer,intent(in)  :: at(n)                ! integer atom types
      integer,intent(inout) :: eqv(n,3)

      real(wp),allocatable :: cn(:)
      real(wp),allocatable :: bond(:,:)

      integer :: p,nonh
      integer :: a,b,c

      logical :: meth
      integer :: nmeth
      integer :: hyd(3)

      allocate(cn(n),bond(n,n))

      call countnonh(n,at,nonh) ! get number of non-Hydrogen atoms

      cn=0.0d0
      bond=0.0d0
      call xcoord(n,at,xyz,cn,bond)
      
      eqv=0

      nmeth=0
      do p=1,n                                 
        if(at(p).ne.6) cycle !only check for carbon
        call methyl(p,at,n,cn,bond,meth,hyd)
        !if(meth)nmeth=nmeth+1
        if(meth)then
          a=hyd(1)
          b=hyd(2)
          c=hyd(3)
          eqv(a,1) = a
          eqv(a,2) = b
          eqv(a,3) = c
          eqv(b,1) = b
          eqv(b,2) = c
          eqv(b,3) = a
          eqv(c,1) = c
          eqv(c,2) = a
          eqv(c,3) = b
        endif
      enddo
      deallocate(bond,cn)

end subroutine get_methyl

!! --------------------------------------------------------------------------------------
!! --------------------------------------------------------------------------------------

subroutine methyl(i,iat,nat,cn,bond,meth,hydrogens)
      use crest_parameters, only: wp
      implicit none
      integer,intent(in)  :: i
      integer,intent(in)  :: iat(nat)
      integer,intent(in)  :: nat
      real(wp),intent(in) :: cn(nat)
      real(wp),intent(inout) :: bond(nat,nat)
      logical,intent(out) :: meth
      integer,intent(out) :: hydrogens(3)
      integer  :: j
      integer  :: icn
      integer  :: k,l
      integer  :: nhyd
      integer,allocatable :: neighb(:)
      meth=.false.
      hydrogens(1:3)=0
      l=0
      icn=floor(cn(i))
      allocate(neighb(icn))                           !identifier array contains atoms
      nhyd=0
      do k=1,icn
         j=maxloc(bond(:,i),1)
         bond(j,i)=0.0d0
         neighb(k)=j
         if(iat(j).eq.1)nhyd=nhyd+1
      enddo
      if((icn.eq.4).and.(nhyd.eq.3))then
        !write(*,*) neighb
        meth=.true.
        do k=1,icn
           j=neighb(k)
           !write(*,*) j,iat(j)
           if(iat(j).eq.1)then
              hydrogens(nhyd)=j
              nhyd=nhyd-1
           endif
        enddo
        !write(*,*) hydrogens
      endif
end subroutine methyl
