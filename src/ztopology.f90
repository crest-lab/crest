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

!================================================================================!
!  This file contains many subroutines related to topology setup of a molecule.
!  It requires the polymorphic datatype "zmolecule", which can be found in the
!  "zdata.f90" module.
!
!  Most routines are based on a representation where the molecule is seen as a
!  "graph" and the topology is analyzed by navigating through the different
!  paths of the graph.
!
!  P. Pracht, 2019/20
!================================================================================!
!================================================================================!
!  SIMPLETOPO: set up topology.
!  All the information is stored in the complex datatype object "zmol"
!  Although the name implies it, there is nothing much "simple" for some
!  of the routines ...
! 
!  simpeltopo_file version is a file wrapper
!  simpeltopo_mol version is a wrapper for a coord object
!================================================================================!
subroutine simpletopo_file(fname,zmol,verbose,getrings,wbofile)
     use iso_fortran_env, wp => real64
     use zdata
     use strucrd, only: rdnat,rdcoord
     implicit none
     character(len=*) :: fname
     type(zmolecule)  :: zmol
     logical          :: verbose
     logical          :: getrings
     character(len=*),optional :: wbofile
     integer :: n
     integer,allocatable :: at(:)
     real(wp),allocatable :: xyz(:,:)
     logical :: ex
     interface
         subroutine simpletopo(n,at,xyz,zmol,verbose,getrings,wbofile)
             import :: zmolecule, wp
             implicit none
             type(zmolecule)  :: zmol       
             logical          :: verbose
             logical          :: getrings
             integer,intent(in)  :: n
             integer,intent(in)  :: at(n)
             real(wp),intent(in) :: xyz(3,n) !in Bohrs
             character(len=*),intent(in),optional :: wbofile
         end subroutine simpletopo
     end interface

     call rdnat(fname,n)
     allocate(at(n),xyz(3,n))
     call rdcoord(fname,n,at,xyz)

     ex = .false.
     if(present(wbofile))then
      inquire(file=wbofile,exist=ex)
     endif
     if(.not.ex)then
        call simpletopo(n,at,xyz,zmol,verbose,getrings,'')    
     else
        call simpletopo(n,at,xyz,zmol,verbose,getrings,wbofile)
     endif
     
     deallocate(xyz,at)
     return
end subroutine simpletopo_file

subroutine simpletopo_mol(mol,zmol,verbose,getrings)
     use iso_fortran_env, wp => real64
     use zdata
     use strucrd
     implicit none
     type(coord)      :: mol    !in 
     type(zmolecule)  :: zmol   !out
     logical          :: verbose
     logical          :: getrings
     interface
         subroutine simpletopo(n,at,xyz,zmol,verbose,getrings,wbofile)
             import :: zmolecule, wp
             implicit none
             type(zmolecule)  :: zmol       
             logical          :: verbose
             logical          :: getrings
             integer,intent(in)  :: n
             integer,intent(in)  :: at(n)
             real(wp),intent(in) :: xyz(3,n) !in Bohrs
             character(len=*),intent(in),optional :: wbofile
         end subroutine simpletopo
     end interface
     call simpletopo(mol%nat,mol%at,mol%xyz,zmol,verbose,getrings,'')
     return
end subroutine simpletopo_mol

!================================================================================!
!  SIMPLETOPO: set up topology.
!  n   - number of atoms
!  at  - atom types (as atom number)
!  xyz - cartesian coordinates (in Bohrs)
!  zmol - polymorphic datatype containing the molecule data
!  verbose - boolean to activate printouts
!  wbofile - (optional) name of the file containing WBOs
!================================================================================!
subroutine simpletopo(n,at,xyz,zmol,verbose,getrings,wbofile)
     use iso_fortran_env, wp => real64
     use zdata
     implicit none
     type(zmolecule)  :: zmol       
     type(zmolecule)  :: zfrag
     type(zring)      :: newring
     logical          :: verbose
     logical          :: getrings
     integer,intent(in)  :: n
     integer,intent(in)  :: at(n)
     real(wp),intent(in) :: xyz(3,n) !in Bohrs
     character(len=*),intent(in),optional :: wbofile
     type(zatom),allocatable :: zat(:)
     real(wp),allocatable :: cn(:)
     real(wp),allocatable :: bond(:,:)
     real(wp),allocatable :: wbo(:,:)
     real(wp),allocatable :: rcov(:)
     real(wp) :: dummy
     integer :: i,j,k,l
     integer :: ntopo
     logical,allocatable :: neighmat(:,:)
     integer,allocatable :: topovec(:)
     integer :: nrings
     character(len=10) :: numchar

     logical :: ex,useWBO

!--- header
     if(verbose)then
       write(*,*)
       call smallhead('TOPOLOGY ANALYSIS') 
     endif

!--- set covalent radii and calculate coordination numbers
     allocate(rcov(94))
     call setrcov(rcov)
     allocate(cn(n),bond(n,n))
     call xcoord2(n,at,xyz,rcov,cn,900.0_wp,bond)
!--- read in WBOs if required
     ex = .false.
     if(present(wbofile))then
      inquire(file=wbofile,exist=ex)
     endif
     if(ex)then
          allocate(wbo(n,n))
          call readwbo("wbo",n, wbo)
          !write(*,*) wbo
          !stop
          useWBO=.true.
      else
        useWBO=.false.
     endif   

     allocate(zat(n))
     if(.not.useWBO)then
!--- get individual neighbour lists and set up molecule object "zat"
       ntopo = n*(n+1)/2
       allocate(neighmat(n,n), source=.false.)
       allocate(topovec(ntopo))
       call bondtotopo(n,at,bond,cn,ntopo,topovec,neighmat)
       call neighbourset(zmol,n,at,xyz,cn,ntopo,topovec)
       deallocate(topovec,neighmat)
     else
!--- neighbour list could also be set up from WBOs, which is not necessarily better
        do i=1,n
          zat(i)%cn = cn(i)
          call wboneighbours(i,n,at,xyz,wbo,0.7d0,zat(i))
          call wborepaireta(i,n,at,xyz,cn,wbo,zat(i))  !handle eta-coordinated atoms
        enddo

     !--- transfer information to zmol class
       zmol%nat=n      !number of atoms
       zmol%zat=zat    !list of atoms, including data
       zmol%at=at      !list of atom types (integer)
     endif
     deallocate(zat)

!--- analyze system for fragments
     call zmol%mrec()

     allocate(zmol%distmat(n,n))
     do i=1,n
        do j=1,n
          zmol%distmat(i,j)=(xyz(1,i)-xyz(1,j))**2 + &
        &                   (xyz(2,i)-xyz(2,j))**2 + &
        &                   (xyz(3,i)-xyz(3,j))**2
          zmol%distmat(i,j)=sqrt(zmol%distmat(i,j))
        enddo
     enddo     

!--- printouts?
     if(verbose)then
      write(*,'(1x,a)')'Coordinates (Bohr):'
      do i=1,n
        call zmol%zat(i)%wrtmline(6)
      enddo
      write(*,*)
      call zmol%wrtable(6)
      write(*,*)

      if(allocated(zmol%molvec))then
        write(*,'(1x,a,i0)') 'Number of fragments in the system:  ',zmol%nfrag
        write(*,*)
      endif
     endif

!--- identify rings
     if(getrings)then
       do i=1,zmol%nfrag
         call zmol%fragment(i,zfrag)
         zfrag%maxring = maxringident !maxringident is a global variable from zdata.f90
         call countrings(zfrag,nrings)
         if(verbose.and.nrings>0)then
             if(zmol%nfrag>1)then
             write(*,'(1x,a,i0)')'Fragment ',i
             write(*,'(1x,a,i0,/)')'Total number of rings in the fragment: ',nrings
             else
             write(*,'(1x,a,i0,/)')'Total number of rings in the system: ',nrings
             endif
         endif
         if(nrings.ge.1)then
          allocate(zfrag%zri(nrings))
          call newgetrings(zfrag,.false.)
          do j=1,zfrag%nri
             newring=zfrag%zri(j)
             do k=1,newring%rs
                newring%rlist(k) = zfrag%revmap(newring%rlist(k))
             enddo
             call zmol%addring(newring)
             call newring%deallocate()
          enddo
         endif
         call zfrag%deallocate()
       enddo
       if(verbose)then
         call zmol%prrings(6)
         if(zmol%nri>0)then
         write(*,'(/,1x,a,i0,/)')'Total number of rings in the system: ',zmol%nri
         endif
       endif
     endif

!--- deallocation of memory
     if(allocated(wbo))deallocate(wbo)
     deallocate(bond,cn)
     deallocate(rcov)
     return
end subroutine simpletopo

!=======================================================================!
!C compute coordination numbers by adding an inverse damping function
!=======================================================================!
subroutine xcoord2(nat,iz,xyz,rcov,cn,cn_thr,bond)
      use iso_fortran_env, wp => real64
      implicit none
      integer,intent(in) :: nat
      integer,intent(in) :: iz(nat)
      real(wp),intent(in) :: xyz(3,nat)
      real(wp),intent(out) :: cn(nat)
      real(wp),intent(in)  :: cn_thr
      real(wp),intent(in)  :: rcov(94)
      real(wp),intent(out) :: bond(nat,nat)
      integer :: i,j,k1
      integer :: iat
      real(wp) :: dx,dy,dz,r,damp,xn,rr,rco,r2,rcovi,rcovj
      k1=16
      bond=0.0d0
      cn=0.0d0
      do i=1,nat
        xn=0.0d0
        rcovi=rcov(iz(i))
        do iat=1,nat
         if(iat.ne.i)then
            dx=xyz(1,iat)-xyz(1,i)
            dy=xyz(2,iat)-xyz(2,i)
            dz=xyz(3,iat)-xyz(3,i)
            r2=dx*dx+dy*dy+dz*dz
            r=sqrt(r2)
            if (r2.gt.cn_thr) cycle
            rcovj=rcov(iz(iat))
! covalent distance in Bohr
            rco=(rcovi+rcovj)*1.0  ! this scaling reduces the size of the clusters
            rr=rco/r
! counting function exponential has a better long-range behavior than MHGs inverse damping
            damp=1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
            bond(iat,i)=damp
            xn=xn+damp
         endif
        enddo
        cn(i)=xn
      enddo
      return
end subroutine xcoord2

!===================================================!
! generate the topo array for a given structure
!===================================================!
subroutine bondtotopo(nat,at,bond,cn,ntopo,topo,neighbourmat)
       use iso_fortran_env, only : wp=>real64
       integer,intent(in)  :: nat
       integer,intent(in) :: at(nat)
       real(wp),intent(inout) :: bond(nat,nat)
       real(wp),intent(in) :: cn(nat)
       integer,intent(in)  :: ntopo
       integer,intent(out) :: topo(ntopo)
       real(wp),allocatable :: cn2(:)
       logical,intent(inout) :: neighbourmat(nat,nat)
       integer :: i,j,k,l
       integer :: icn,rcn
       integer :: lin
       allocate(cn2(nat),source=0.0_wp)
       topo = 0
       neighbourmat=.false.
       !--- some heuristic rules and CN array setup
       do i=1,nat
          cn2(i) = cn(i) 
          rcn = floor(cn(i))
          select case( at(i) ) !additional empirical topology rules
          ! case( 5 ) !B
          !   if( nint(cn(i)) > 4) cn2(i)=4.0_wp
          ! case( 9,17,35,53 ) !F,Cl,Br,I 
          !   cn2(i) = min(cn(i),1.0_wp) 
           case( 6 ) !C
            if((cn(i)-rcn)<0.7_wp)then
            cn2(i) = rcn
            endif
          end select
          !-- extreme CN cases
          if( nint(cn(i))>8)cn2(i) =8.0_wp
          !empirical: rounding down up to .6 is better for topo setup
          if((cn(i)-rcn)<0.6_wp)then 
              cn2(i) = rcn
          endif
       enddo
       !--- build the topology
       do i=1,nat
          icn=nint(cn2(i))
          do k=1,icn
            j=maxloc(bond(:,i),1)
            bond(j,i)=0.0d0
            if (i .eq. j) cycle
            neighbourmat(i,j)=.true. !--important: not automatically (i,j)=(j,i)
          enddo
       enddo
       do i=1,nat
         do j=1,nat
         if(i==j) cycle
         l = lin(i,j)
         !-- only save matching topology --> prevent high CN failures
         if(neighbourmat(i,j).and.neighbourmat(j,i))then
            topo(l) = 1
         else
            ! special case for carbon (because the carbon CN is typically correct)
            ! this helps, e.g. with eta-coordination in ferrocene
            ! (used, except if both are carbon)
            if(.not.(at(i)==6 .and. at(j)==6))then
             if(at(i)==6 .and.neighbourmat(i,j)) topo(l)=1
             if(at(j)==6 .and.neighbourmat(j,i)) topo(l)=1
            endif
         endif
         enddo
       enddo
       deallocate(cn2)
       return
end subroutine bondtotopo

subroutine quicktopo(nat,at,xyz,ntopo,topovec)
    use iso_fortran_env, only : wp => real64
    implicit none
    integer :: nat
    integer :: at(nat)
    real(wp) :: xyz(3,nat) !must be in Bohrs
    integer :: ntopo
    integer :: topovec(ntopo)
    real(wp),allocatable :: rcov(:),cn(:),bond(:,:)
    logical,allocatable :: neighmat(:,:)
    allocate(rcov(94))
    call setrcov(rcov)
    allocate(bond(nat,nat),cn(nat), source=0.0_wp)
    allocate(neighmat(nat,nat), source=.false.)
    cn=0.0d0
    bond=0.0d0
    call xcoord2(nat,at,xyz,rcov,cn,900.0_wp,bond)
    call bondtotopo(nat,at,bond,cn,ntopo,topovec,neighmat)
    deallocate(neighmat,cn,bond,rcov)
    return
end subroutine quicktopo


!-- transfer topology data to the zmol object
subroutine neighbourset(zmol,nat,at,xyz,cn,ntopo,topovec)
      use iso_fortran_env, only : wp => real64
      use zdata
      use strucrd, only: i2e
      implicit none
      type(zmolecule) :: zmol
      integer,intent(in)  :: nat
      real(wp),intent(in) :: xyz(3,nat)
      integer,intent(in)  :: at(nat)
      real(wp),intent(in) :: cn(nat)
      integer,intent(in)     :: ntopo
      integer,intent(in) :: topovec(ntopo)
      type(zatom)         :: zat    !--- "zat" is the complex datatype object for atom i
      integer :: lin
      integer :: i,j,k,l
      integer :: inei

      allocate(zmol%zat(nat))
      zmol%nat=nat    !number of atoms
      zmol%at=at      !list of atom types (integer)
      do i=1,nat
      zmol%zat(i)%atype=at(i)
      zmol%zat(i)%cart(:)=xyz(:,i)
      zmol%zat(i)%el=i2e(at(i),'nc')
      zmol%zat(i)%cn=cn(i)
      zmol%zat(i)%pos=i
        inei = 0
        do j=1,nat
         if(i==j) cycle
         l = lin(i,j)
         if(topovec(l)==1) inei=inei+1
        enddo
        allocate(zmol%zat(i)%ngh(inei))
        allocate(zmol%zat(i)%ngt(inei))
        k=0
        do j=1,nat
         if(i==j)cycle
         l = lin(i,j)
           if(topovec(l)==1)then 
              k=k+1
              zmol%zat(i)%ngh(k)=j      ! the k-th neighbour of atom i is atom j
              zmol%zat(i)%ngt(k)=at(j)  ! atom j has this atom type
           endif
        enddo
        zmol%zat(i)%nei = inei
        call quicksort(inei,zmol%zat(i)%ngh)
        !call quicksort(inei,zmol%zat(i)%ngt)
        do j=1,inei
          zmol%zat(i)%ngt(j)=at(zmol%zat(i)%ngh(j))
        enddo
      enddo

      return
end subroutine neighbourset


!=======================================================================!
!  compute Zmat (internal coordinates) for the zmol struct
!=======================================================================!
subroutine ztopozmat(zmol,pr)
    use iso_fortran_env, wp => real64
    use zdata
    implicit none
    type(zmolecule)  :: zmol       
    logical          :: pr

    integer :: nat
    real(wp),allocatable :: xyz(:,:)
    real(wp),allocatable :: geo(:,:)
    integer,allocatable :: na(:),nb(:),nc(:)
    integer :: i,j,k,l
    real(wp),parameter :: pi =  3.14159265358979D0
    real(wp),parameter :: rad = 180.0d0/pi


    nat = zmol%nat
    allocate(geo(3,nat),xyz(3,nat), source= 0.0d0)
    allocate(na(nat),nb(nat),nc(nat), source = 0)

    do i=1,nat
       xyz(1:3,i) = zmol%zat(i)%cart(1:3)
    enddo   

    !call xyzint(xyz,nat,na,nb,nc,rad,geo)
    call xyzint(xyz,nat,na,nb,nc,1.0d0,geo)
    !call xyzgeo(xyz,nat,na,nb,nc,rad,geo)

    if(pr)then
       call smallhead('INTERNAL COORDINATES (ZMATRIX)') 
       do i=1,nat
       write(*,'(2x,a,4x,3f10.4,1x,3i5)') zmol%zat(i)%el,geo(1:3,i),na(i),nb(i),nc(i)
       enddo
    endif

    call move_alloc(geo,zmol%zmat)
    call move_alloc(na,zmol%zna)
    call move_alloc(nb,zmol%znb)
    call move_alloc(nc,zmol%znc)

    deallocate(xyz)
    return
end subroutine ztopozmat

!=======================================================================!
! With the CN-based neighbourlist it can happen that two atoms
! are shared as a neighbour only in the list of one of the atoms.
! This needs correction and is one of the reasons why WBO-based
! neighbourlists are better.
!=======================================================================!
subroutine crosscheckCNnei(zmol)
      use iso_fortran_env, wp => real64
      use zdata
      implicit none
      type(zmolecule) :: zmol
      integer :: i,j,k,l
      integer :: nei, newnei
      integer,allocatable :: newngh(:)
      integer,allocatable :: newngt(:)
      do i=1,zmol%nat
         nei = zmol%zat(i)%nei
         newnei=nei
         allocate(newngh(nei),source = 0)
         allocate(newngt(nei), source = 0)
         do j=1,nei
            k = zmol%zat(i)%ngh(j)
            if(any(zmol%zat(k)%ngh(:) .eq. i)) then
               newngh(j) = k 
               newngt(j) = zmol%at(k)
               cycle !atom i is also the neighbour of its own neighbour
            else
               !else, the atoms seem to be "artificial" neighbours and
               !has to be removed from the neighbourlist
               newngh(j) = 0 
               newnei=newnei-1
            endif
         enddo
      !--- if the neighbour list was changed we update it
         if(newnei.lt.nei)then
           deallocate(zmol%zat(i)%ngh, zmol%zat(i)%ngt)
           allocate(zmol%zat(i)%ngh(newnei), zmol%zat(i)%ngt(newnei))
           l = 0
           do j=1,nei
              if(newngh(j).ne.0)then
               l = l +1
               zmol%zat(i)%ngh(l) = newngh(j)
               zmol%zat(i)%ngt(l) = newngt(j)
             endif
           enddo
           zmol%zat(i)%nei = newnei
         endif
         deallocate(newngt,newngh)
      enddo

      return
end subroutine crosscheckCNnei

!=======================================================================!
!  get the neighbours of an atom based on the WBOs
!=======================================================================!
subroutine wboneighbours(i,nat,at,xyz,wbo,wbothr,zat)
      use iso_fortran_env, only : wp => real64
      use zdata
      use strucrd, only: i2e
      implicit none
      integer,intent(in)  :: i
      real(wp),intent(in) :: xyz(3,nat)
      real(wp),intent(in) :: wbo(nat,nat)
      real(wp),intent(in) :: wbothr
      integer,intent(in)  :: at(nat)
      integer,intent(in)  :: nat
      type(zatom)         :: zat    !--- "zat" is the complex datatype object for atom i 
      integer :: icn
      integer :: k,l,j
     
      zat%atype=at(i)
      zat%cart(:)=xyz(:,i)
      zat%el=i2e(at(i),'nc')
      zat%pos=i

      icn=0
      do k=1,nat
        if(wbo(i,k).ge.wbothr)then
            icn = icn+1
        endif
      enddo
      zat%nei = icn
      allocate( zat%ngh(icn), zat%ngt(icn) )
      l=0          
      do j=1,nat
          if(wbo(i,j).ge.wbothr)then
          l=l+1
          zat%ngh(l)=j      ! the l-th neighbour of atom i is atom j
          zat%ngt(l)=at(j)  ! atom j has this atom type
          endif
      enddo
      call quicksort(icn,zat%ngh)
      !call quicksort(icn,zat%ngt)
      do j=1,icn
         zat%ngt(j)=at(zat%ngh(j))
      enddo

      return
end subroutine wboneighbours 
subroutine wborepaireta(i,nat,at,xyz,cn,wbo,zat)
      use iso_fortran_env, only : wp => real64
      use zdata
      implicit none
      integer,intent(in)  :: i
      real(wp),intent(in) :: xyz(3,nat)
      real(wp),intent(in) :: cn(nat)
      real(wp),intent(in) :: wbo(nat,nat)
      real(wp) :: wbothr
      integer,intent(in)  :: at(nat)
      integer,intent(in)  :: nat
      type(zatom)         :: zat
 
      if(cn(i).gt.1 .and. zat%nei .lt.1)then
          if(allocated(zat%ngt))deallocate(zat%ngt)
          if(allocated(zat%ngh))deallocate(zat%ngh)
         write(*,*) i,at(i),cn(i)
         call wboneighbours(i,nat,at,xyz,wbo,0.1d0,zat)
      endif
      return
end subroutine wborepaireta

!=================================================================!
! subroutine wbomrec
! molcount: number of total fragments (increased during search)
! nat: overall number of atoms
! wbo: bond orders
! wbothr: threshold for when a bond is present
! molvec: assignment vector of atom to fragment
!=================================================================!
subroutine wbomrec(molcount,nat,wbo,wbothr,molvec)
      use iso_fortran_env, wp => real64
      implicit none
      real(wp) :: wbo(nat,nat)
      real(wp) :: wbothr
      integer :: nat,molvec(nat),i,molcount,at(nat)
      logical :: taken(nat)
      molvec=0
      molcount=1
      taken=.false.
      do i=1,nat
       if(.not.taken(i)) then
         molvec(i)=molcount
         taken(i)=.true.
         call wbofrags(i,nat,wbo,wbothr,taken,molvec,molcount)
         molcount=molcount+1
      endif
      enddo
      molcount=molcount-1
      return
end subroutine wbomrec
!- depth-first search through the WBO list to determine all connected atoms
recursive subroutine wbofrags(i,nat,wbo,wbothr,taken,molvec,molcnt)
      use iso_fortran_env, wp => real64
      implicit none
      integer :: i,nat
      real(wp) :: wbo(nat,nat)
      real(wp) :: wbothr
      integer :: molcnt,molvec(nat),j,iat(nat),k
      logical taken(nat)
      do j=1,nat
         if(i .eq. j) cycle
         if( wbo(i,j) .ge. wbothr)then
         if (.not.taken(j)) then
            molvec(j)=molcnt
            taken(j)=.true.
            call wbofrags(j,nat,wbo,wbothr,taken,molvec,molcnt)
         endif
         endif
      enddo
      return
end subroutine wbofrags

!=======================================================================!
!  get the adjacency and edgelength matrices from the neighbour lists
!=======================================================================!
subroutine zadjacent(zmol,A,E)
      use iso_fortran_env, wp => real64
      use zdata
      implicit none
      type(zmolecule) :: zmol
      integer :: A(zmol%nat,zmol%nat)
      real(wp) :: E(zmol%nat,zmol%nat)
      integer :: i,j,k,l
      integer :: nat
      nat = zmol%nat
      A=0
      E=0.0_wp
      do i=1,nat
        do j=1,nat
          if(any(zmol%zat(i)%ngh(:) .eq. j))then  !only include bonds from the neighbour lists
            E(i,j) = zmol%distmat(i,j)
            A(i,j) = 1
          else
           cycle
          endif
        enddo
      enddo
      return
end subroutine zadjacent

!=======================================================================!
!  compare two arrays and check if their content is identical
!=======================================================================!
subroutine arrcomp(n,narr,m,marr,equi)
      implicit none
      integer :: n
      integer :: narr(*)   
      integer :: m
      integer :: marr(*)
      logical :: equi
      logical,allocatable :: mask(:)
      integer :: i,j,k,l
      if(n.ne.m)then
         equi=.false.
         return
      endif
      allocate( mask(n) )
      mask(:) = narr(1:n) .eq. marr(1:n)
      if(any(.not.mask(:)))then
         equi=.false.
      else
         equi =.true.
      endif
      deallocate(mask)
      return
end subroutine arrcomp

!=======================================================================!
! compare two arrays and check if they have exactly one element in common
!=======================================================================!
function arruniqel(n,narr,m,marr,el)
      implicit none
      logical :: arruniqel
      integer :: n
      integer :: narr(n)
      integer :: m
      integer :: marr(m)
      integer :: el
      integer :: incr
      integer :: i,j
      arruniqel = .false.
      el = 0
      incr = 0
      do i=1,n
        if( any(marr(:).eq.narr(i)) )then
          incr = incr+1
          el = narr(i)
        endif      
      enddo
      if(incr .eq. 1)then
        arruniqel = .true.
        !el = el
      else
        arruniqel = .false.
        el = 0 
      endif
      return
end function arruniqel

!=========================================================================!
! check if all atoms contained in narr are in the neighbour list of atom i
!=========================================================================!
logical function allneigh(zmol,i,n,narr)
      use zdata
      implicit none
      type(zmolecule) :: zmol
      integer :: i
      integer :: n
      integer :: narr(n)
      integer :: j
      logical :: dum
      dum = .true.
      do j=1,n
        dum = dum .and. any(zmol%zat(i)%ngh(:)==narr(j))
      enddo
      allneigh = dum
      return
end function allneigh

!=======================================================================!
! get information about the side chain
!=======================================================================!
subroutine getsideweight(zmol,i,j,mside)
      use iso_fortran_env, only : wp => real64
      use zdata
      use atmasses
      implicit none
      integer,intent(in) :: i  !central atom
      integer,intent(in) :: j  !first atom of the side chain
      type(zmolecule) :: zmol
      real(wp),intent(out) :: mside

      integer :: k,l
      integer :: atm

      logical,allocatable :: taken(:)
      integer,allocatable :: path(:)
  
      mside=0.0_wp

      allocate(path(zmol%nat), source=0)
      l=0
      call recside(zmol,i,j,path,l)

      do k=1,zmol%nat
         if(path(k).eq.0) cycle
         l=path(k)
         atm=zmol%zat(l)%atype
         mside=mside + ams(atm) 
      enddo

      deallocate(path)
      return
end subroutine getsideweight

!=======================================================================!
!  count the number of rings in the system
!  (just a call on the same "newgetrings" routine,
!   but some of the data has to be reset)
!=======================================================================!
subroutine countrings(zmol,nrings)
      use iso_fortran_env, only : wp => real64
      use zdata
      implicit none
      type(zmolecule) :: zmol
      integer,intent(out) :: nrings
      integer :: i
      nrings = 0
      call newgetrings(zmol,.false.)
      nrings=zmol%nri
      do i=1,zmol%nat
         zmol%zat(i)%ring =.false.
      enddo
      return
end subroutine countrings

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!================================================================================================!
! new "getring" routine. 
! very annoying recursive stuff.
! requires a FULLY SET UP TOPOLOGY in the "zmol" datatype 
!
! The search algorithm is some kind of a taylored form of
! an A* search algorithm, fit specially for finding rings.
!
! The following routines belong together:
! newgetrings, startring, recurring2, recside
!
!   DO. NOT. TOUCH.
!
!================================================================================================!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
subroutine newgetrings(zmol,verbose)
      use iso_fortran_env, only : wp => real64
      use zdata
      implicit none
      type(zmolecule) :: zmol
      type(zatom) :: za
      type(zring) :: zri 
      logical :: verbose 
      integer :: k,l
      integer :: ric  
      integer :: atms
      logical :: newring
      logical,parameter :: safety = .true.
        !If safety=.false. the code tries to get ALL rings.
        !But the algo then also counts rings double.
        !I.e. benzene would "have" 6 rings á 6 atoms.
        !With the safety option it stops when all atoms are distributed
        !in rings, which is much more robust
        !E.g. Adamantan will have 3 unique rings instead of correctly 4
        !For now this is the lesser evil.
        !A "best of both worlds"-approach would simply require an
        !on-the-fly evaluation if we already got the identified ring.
        !I marked the position where this would need to be done below.
 
      zmol%nri = 0 !reset number of rings
  
      atms=zmol%nat
      if(zmol%nonterm.lt.1)then
        do k=1,atms
           if(zmol%zat(k)%nei .gt. 1 )then
              zmol%nonterm=zmol%nonterm + 1
           endif
        enddo
      endif

      if(zmol%maxneigh.lt.1)then
        do k=1,atms
         if(zmol%zat(k)%nei .gt. zmol%maxneigh )then
             zmol%maxneigh=zmol%maxneigh + 1
         endif
       enddo
      endif

      ric=0
      do k=1,atms                 !loop over the molecule, i.e., all atoms
         newring=.false.
         za=zmol%zat(k)
         if(za%ring.and.safety) cycle        !cycle atoms already in a ring
         if(za%nei .eq. 1) cycle  !cycle terminal atoms
         if(za%nei .ge. zmol%limitnei) cycle  !cycle highly coordinated atoms (e.g. Cp-->M)
         call startring(zmol,k,newring,zri,verbose)
         if(allocated(zmol%zri).and.newring)then
            ric=ric+1
            zmol%zri(ric) = zri
         endif
         if(newring) zmol%nri = zmol%nri + 1  !add 1 to the ring count
      enddo

      call zri%deallocate !clear space (if allocated)
      return
end subroutine newgetrings

subroutine startring(zmol,k,newring,zri,verbose)
      use iso_fortran_env, only : wp => real64
      use zdata
      implicit none
      type(zmolecule) :: zmol
      type(zatom) :: za,zb
      type(zring) :: zri 
      logical,intent(out) :: newring
      integer,intent(in) :: k          !the starting atom
      logical,intent(in) :: verbose      

      integer :: i,j,l,q,m,n
      integer :: w,v
      integer :: nei

      logical,allocatable :: taken(:)     
      integer :: ntak,tref
      integer,allocatable :: dum(:)
      integer,allocatable :: path(:),path2(:)
      integer,allocatable :: hood(:)
      integer :: npath

      logical,allocatable :: ringpossible(:)
     

      call zri%deallocate !clear space (if allocated)
      newring=.false.     !assume that there is no ring
      za=zmol%zat(k)      !get the atom
      nei=za%nei          !get number of neighbouring atoms
      ntak = 0            !taken counter initialize

      allocate(taken(zmol%nat), source=.false.)
      allocate(path(zmol%nonterm), source=0)
      allocate(path2(zmol%nonterm), source=0)
      taken(k)=.true.     !the start atom obviously has to be taken for the ring,
                          !but must not be included in the counter "ntak"
      npath=0


      !first we have to check if a ring is possible in the first place
      !therefore we have to get the entire neighborhood, i.e., all
      !atoms chained to each of the neighbors of k, and check if two of
      !those neighbors share it.
      allocate(ringpossible(nei), source=.false.)
      allocate(hood(zmol%nat), source=0)
      HOODANALY : do i=1,nei
        if(ringpossible(i))cycle ! we already included this neighboring chain, skip it.
        hood=0
        j=za%ngh(i)
        q=0
        call recside(zmol,k,j,hood,q)
        do w=1,nei
           if(i.eq.w) cycle
           v=za%ngh(w)
           if(any(hood(:)==v))then
             ringpossible(i) = .true. 
             ringpossible(w) = .true.
           endif
        enddo
      enddo HOODANALY
      deallocate(hood)

      !now starts the identification of the rings
      do i=1,nei
         if(.not. ringpossible(i) ) cycle
         if(zmol%maxring.lt.3)then
            tref=zmol%nonterm
         else
            tref=zmol%maxring
         endif
         j=za%ngh(i)
         call recurring2(zmol,k,j,taken,path,path2,ntak,tref)
         if(ntak .gt. 0)then
           exit  !a ring was found! exit the loop
         endif
      enddo
      deallocate(ringpossible)

!---------------
      !how to handle our newfound ring
      if(ntak .gt. 0)then 
        allocate(dum(ntak))
         i=0
         do j=1,tref
            if(path2(j).gt.0)then
            i=i+1
            dum(i)=path2(j)
            l=dum(i)
            zmol%zat(l)%ring=.true.
            endif
        enddo
        call quicksort(ntak,dum)

 !Here we would correctly need to check if we already have
 !identified the ring, and therefore it is a dublicate.
 !but for now we just continue.
        newring =.true.

        !--- if a (truely) new ring was found,show it and save it
        if(newring)then
          !--- save ring data to object (collected outside the routine)
          zri%rs=ntak
          allocate(zri%rlist(ntak), source=0)
          zri%rlist=dum        
          !--- short printout
          if(verbose)then
            call zri%print(6)
          endif
        endif

        deallocate(dum)
      endif
!---------------

      deallocate(path2,path)
      deallocate(taken)
      return
end subroutine startring


recursive subroutine recurring2(zmol,k,j,taken,path,path2,ntak,tref)
      use iso_fortran_env, only : wp => real64
      use zdata
      implicit none                                                       
      type(zmolecule) :: zmol
      integer,intent(in) :: k          !the starting atom
      integer,intent(in) :: j          !the coordinated atom
      integer,intent(inout) :: ntak
      logical,intent(inout) :: taken(zmol%nat)
      integer,intent(inout) :: path(zmol%nonterm)
      integer,intent(inout) :: path2(zmol%nonterm)
      integer,intent(inout) :: tref
    
      integer :: npathold,npathmax
      integer :: nref,iref,nmax

      integer :: i,l,n,m
      integer :: p,q

!      logical,allocatable :: takedum(:)
      integer,allocatable :: ndum(:)
      real(wp),allocatable :: dists(:)

      logical :: ring,dumreset

      ring = .false.

      if(j.eq.k .and. ntak.eq.1) return !for the very first iteration we have to exclude the starting
                                        !atom in this way. (i.e., there are no two-membered rings)

      if(zmol%zat(j)%nei.eq.1) return !for terminal atoms return

      if(zmol%zat(j)%nei.ge.zmol%limitnei) return !for highly coordinated neighbours return (e.g. Cp-->M)
     
      if(k.eq.j)then        !we arrived at the starting atom again! therefore a ring was found.
          ntak=ntak+1
          path(ntak)=k
          if(ntak.lt.tref) tref=ntak
          return
      endif
       
      if(any(path(:)==j)) return !if we already passed the atom, return (i.e., no walking back)

      if(ntak+1 .ge. tref) return !don't go longer paths as the already known ones
      if(ntak+1 .ge. zmol%nonterm) return

      !--- the atom passed the checks, so we consider it for now
      ntak=ntak+1
      path(ntak)=j
      npathold=ntak+1
      npathmax=zmol%nonterm

      nref=ntak

      n=zmol%nat
      m=zmol%zat(j)%nei            !number of next neighbors
      allocate(ndum(m),dists(m))   
      ndum=ntak           
      
      do i=1,m                     !distances to the atom from which the ring was started
         l=zmol%zat(j)%ngh(i)
         dists(i)=zmol%distmat(k,l)
      enddo

      !--- continue with the path along the neighbour list
      do i=1,m
         path(nref+1:npathmax)=0
         !l=zmol%zat(j)%ngh(i)
         p=minloc(dists(:),1)     !|
         dists(p)=1000000.0_wp    !|try to walk the innermost (smallest) cycle by a distance criterium
         l=zmol%zat(j)%ngh(p)     !|
         call recurring2(zmol,k,l,taken,path,path2,ndum(i),tref)
         if(ndum(i).gt.nref) then
            ntak=ndum(i)
            path2=path    !the shortest path is saved
         endif
      enddo
    
      if(any(ndum(:).gt.nref))then
         ring=.true.
         path=path2       !the saved shortest path is returned
      endif

      deallocate(dists,ndum)

      if(.not.ring)then
         ntak=ntak-1      
         taken(j)=.false.
      endif

      return
end subroutine recurring2

!===============================================================!
!  get all atoms of the side chain l of k and write them to path
!===============================================================!
recursive subroutine recside(zmol,k,l,path,j)
      use iso_fortran_env, only : wp => real64
      use zdata
      implicit none
      type(zmolecule) :: zmol
      integer,intent(in) :: k          !the starting atom
      integer,intent(in) :: l          !current atom
      integer,intent(inout) :: j       !postion counter
      integer,intent(inout) :: path(zmol%nat)  !entire path
      integer :: i,p
   
      if( k .eq. l )        return
      if( any(path(:)==l) ) return
      if(j+1 .gt. zmol%nat) return

      j=j+1
      path(j) = l

      if( zmol%zat(l)%nei .gt. 1 )then
      do i=1,zmol%zat(l)%nei
         p=zmol%zat(l)%ngh(i)
         call recside(zmol,k,p,path,j)
      enddo
      endif

      return 
end subroutine recside
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!================================================================================================!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!


!==========================================================!
! "samering" is a function (return value .true./.false.)
! that is used to check is a given list of atoms is part
! of the same ring.
!
! On Input:  zmol - molecule and topology
!            n    - # of atoms to be checked
!            atms - array containing the atoms
! 
! On Output: function value
!            r - the common ring (if function value = .true.)
!==========================================================!
logical function commonring(zmol,n,atms,r)
     use zdata
     implicit none
     type(zmolecule) :: zmol
     integer,intent(in) :: n
     integer,intent(in) :: atms(n)
     integer,intent(out) :: r
     integer :: i,j
     logical :: l1,l2
     commonring = .false. 
     r = 0
  !--- first check if the atoms are in a ring at all.
  !    if not, we can return immediatly.
     l1=.true.
     do i=1,n
       l1 = l1 .and. zmol%zat(atms(i))%ring
     enddo
     if(.not.l1) return
  !--- loop over all rings in the molecule and check if
  !    the atoms are part of it
     OUTER : do i=1,zmol%nri
        l2 = .true.
        INNER : do j=1,n
          l2 = l2 .and. (any(zmol%zri(i)%rlist(:) .eq. atms(j)))
          if(.not.l2) cycle OUTER
        enddo INNER
        if(l2)then
          commonring = .true.
          r = i
        endif
     enddo OUTER
     return
end function commonring

!==========================================================!
! "ringside" is a subroutine that is used to check 
! in which direction of which atom a ring extends,
! starting from atom i.
!
! On Input:  zmol - molecule and topology
!            n    - starting atom
!            ring - which ring to look at
! 
! On Output: atms - array containing the atoms
! 
!==========================================================!
subroutine ringside(zmol,i,ring,atms)
     use zdata
     implicit none
     type(zmolecule) :: zmol
     integer,intent(in) :: i
     integer,intent(in) :: ring
     integer,intent(inout) :: atms(zmol%nat)
     integer :: j,k,n,l
     integer,allocatable :: path(:)
     n=zmol%nat
     allocate(path(n), source = 0)
 !--- loop over neighbour list of i
     do j=1,zmol%zat(i)%nei
        k = zmol%zat(i)%ngh(j)
      !--- if the neighbour and i are in the selected ring, do nothing
        if(any(zmol%zri(ring)%rlist .eq. k)) cycle
      !--- else, get everything attached to  the neighbour
        path = 0
        l = 0
        call recside(zmol,i,k,path,l)
        call mergearr(n,atms,n,path)  !merge path into atms
     enddo
     deallocate(path)
     
     return
end subroutine ringside
