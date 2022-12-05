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

!====================================================================!
! for given input geometry and specified bond identifiy all
! atoms that belong to the ligand
!====================================================================!
subroutine extractligands(infile,centeratom)
    use crest_parameters
    use zdata
    use strucrd
    implicit none
    character(len=*) :: infile
    integer :: centeratom
    type(coord) :: mol
    type(zmolecule) :: zmol
    integer,allocatable :: path(:)
    integer :: i,j,k,l,nb,ich
    character(len=128) :: atmp
    call mol%open(infile)
    call simpletopo_mol(mol,zmol,.false.,.false.)
    allocate(path(zmol%nat), source=0)
    k=centeratom
    do nb=1,zmol%zat(k)%nei
     path = 0
     j=0
     l=zmol%zat(k)%ngh(nb)
     call  recside(zmol,k,l,path,j)
     write(*,'(a,1x,i0)')'Ligand ',nb
     do i=1,zmol%nat
        if(path(i) .ne. 0)then
        write(*,'(1x,a,i0)',advance='no')trim(i2e(zmol%at(path(i)),'nc')),path(i)
        endif
     enddo
     write(*,*)
     write(atmp,'(a,i0,a)') 'ligand-',nb,'.xyz'
     j=0
     do i=1,zmol%nat
       if(path(i).ne.0)j=j+1
     enddo
     open(newunit=ich,file=atmp)
     write(ich,*)j
     write(ich,*)
     l=zmol%zat(k)%ngh(nb)
     write(ich,'(1x,a2,1x,3f20.10)')i2e(mol%at(l),'nc'),mol%xyz(1:3,l)*bohr
     do i=1,zmol%nat
       if(path(i)==l)cycle
       if(path(i) .ne. 0)then
       j=path(i)    
       write(ich,'(1x,a2,1x,3f20.10)')i2e(mol%at(j),'nc'),mol%xyz(1:3,j)*bohr
       endif
     enddo
     close(ich)
    enddo
    call zmol%deallocate()
    call mol%deallocate()
    return
end subroutine extractligands

subroutine exchangeligands(infile,infile2,centeratom,ligandnr)
    use crest_parameters 
    use zdata
    use strucrd
    implicit none
    character(len=*) :: infile
    character(len=*) :: infile2
    integer :: centeratom,ligandnr
    type(coord) :: mol
    type(zmolecule) :: zmol
    integer,allocatable :: path(:)
    integer :: i,j,k,l,m,nb
    integer :: newcenter,newnb
    character(len=128) :: atmp
    interface
        subroutine matchligands(mol1,mol2,metal)
            import :: coord, wp
            implicit none
            type(coord) :: mol1 !reference ligand
            type(coord) :: mol2 !new ligand
            real(wp),optional :: metal(3) !coordinates of the metal center
        end subroutine matchligands
    end interface
    real(wp) :: refdist !distance of the old ligand
    type(coord) :: oldligand
    type(coord) :: newligand
    type(coord) :: newcomplex
    call mol%open(infile)
    call simpletopo_mol(mol,zmol,.false.,.false.)
    call newligand%open(infile2)

    allocate(path(zmol%nat), source=0)

    k=centeratom
    do nb=1,zmol%zat(k)%nei
     if(nb.ne.ligandnr)cycle
     path = 0
     j=0
     l=zmol%zat(k)%ngh(nb)
     call  recside(zmol,k,l,path,j)
     write(*,'(1x,a,i0)')'Exchange ligand ',nb
     do i=1,zmol%nat
        if(path(i) .ne. 0)then
        write(*,'(1x,a,i0)',advance='no')trim(i2e(zmol%at(path(i)),'nc')),path(i)
        endif
     enddo
     write(*,*)
   !--- extract the ligand as ne coord object
     j=0
     do i=1,zmol%nat
       if(path(i).ne.0)j=j+1
     enddo
     oldligand%nat = j
     allocate(oldligand%at(j), source = 0)
     allocate(oldligand%xyz(3,j), source =0.0d0)
     l=zmol%zat(k)%ngh(nb)
     refdist=zmol%dist(k,l)
     m=1
     oldligand%at(m) = mol%at(l)
     oldligand%xyz(1:3,m) = mol%xyz(1:3,l)
     do i=1,oldligand%nat
       if(path(i)==l)cycle
       j=path(i)    
       m=m+1
       oldligand%at(m) = mol%at(j)
       oldligand%xyz(1:3,m) = mol%xyz(1:3,j)
     enddo
     write(atmp,'(a,i0,a)') 'ligand-',nb,'.xyz'
     call oldligand%write(trim(atmp))
   !--- rotate new ligand w.r.t. old ligand
     call matchligands(oldligand,newligand,mol%xyz(1:3,k))
   
   !--- create new structure
     l = mol%nat - oldligand%nat + newligand%nat
     newcomplex%nat = l
     allocate(newcomplex%at(l), source=0)
     allocate(newcomplex%xyz(3,l), source=0.0_wp)  
     m=0
     do i=1,mol%nat
        if(any(path(:)==i))cycle !skip atoms of old ligand
        m=m+1
        newcomplex%at(m) = mol%at(i) 
        newcomplex%xyz(1:3,m) = mol%xyz(1:3,i)
        if(i==k) newcenter = m
     enddo
     do i=1,newligand%nat
        m=m+1
        newnb = m 
        newcomplex%at(m) = newligand%at(i)
        newcomplex%xyz(1:3,m) = newligand%xyz(1:3,i)
     enddo
     call newcomplex%write('newstructure.xyz')


   !--- cleanup
     call newcomplex%deallocate()
     call oldligand%deallocate()
    enddo

    call newligand%deallocate()
    call zmol%deallocate()
    call mol%deallocate()
    return
end subroutine exchangeligands

subroutine matchligands(mol1,mol2,metal)
    use crest_parameters
    use zdata
    use strucrd
    use ls_rmsd
    implicit none
    type(coord) :: mol1 !reference ligand
    type(coord) :: mol2 !new ligand
    real(wp),optional :: metal(3) !coordinates of the metal center
    logical :: centered
    type(zmolecule) :: zmol1 !topology of mol1
    type(zmolecule) :: zmol2 !topology of mol2
    real(wp) :: c0(3),shift(3)
    integer :: i,j,k,nr,nl,pls
    real(wp) :: distref
    real(wp),allocatable :: c1(:,:),c2(:,:)
    real(wp),allocatable :: b1(:,:),b2(:,:)
    integer :: bdim
    real(wp),allocatable :: tmp(:,:)
    real(wp) ::  g(3,3), U(3,3), x_center(3), y_center(3),rmsdval

!--- check if the central metal atom position was given
    centered = .false.
    pls = 1
    if(present(metal))then
        centered = .true.
        pls = 2
    endif

!--- shift to first atoms
    c0=mol1%xyz(1:3,1)
    shift = c0
    do i=1,mol1%nat
    mol1%xyz(1:3,i) = mol1%xyz(1:3,i)-c0(1:3)
    enddo
    c0=mol2%xyz(1:3,1)
    do i=1,mol2%nat
    mol2%xyz(1:3,i) = mol2%xyz(1:3,i)-c0(1:3)
    enddo

!--- identify atoms to base rotation on (from topology)
    call simpletopo_mol(mol1,zmol1,.false.,.false.)
    call simpletopo_mol(mol2,zmol2,.false.,.false.)
    nr = zmol1%zat(1)%nei  !#neighbours of coordination atom in old ligand
    nl = zmol2%zat(1)%nei  !#neighbours of coordination atom in new ligand
    if((nl == 0).or.((nr==0).and.(.not.centered)))then
       do i=1,mol2%nat
       mol2%xyz(1:3,i) = mol2%xyz(1:3,i)+shift(1:3)
       enddo
       return
    endif
    nr = nr + pls
    nl = nl + pls
    allocate(c1(3,nr),c2(3,nl),source =0.0d0)
    !--- get the reference structures
    c1(1:3,i) = mol1%xyz(1:3,1)
    c2(1:3,i) = mol2%xyz(1:3,1)
    i=1
    do j=1,zmol1%zat(1)%nei
        i=i+1
        k = zmol1%zat(1)%ngh(j)
        c1(1:3,i) = mol1%xyz(1:3,k)
    enddo
    i=1
    do j=1,zmol2%zat(1)%nei
        i=i+1
        k = zmol2%zat(1)%ngh(j)
        c2(1:3,i) = mol2%xyz(1:3,k)
    enddo
    if(centered)then
        c1(1:3,nr) = metal(1:3)-shift(1:3)
        distref = getdist(c1(1:3,1),c1(1:3,nr))
        c2(1:3,nl) = 0.0d0
        call estimatecenter(nl,c2,distref)
        c0(1:3) = c2(1:3,nl) 
    endif
    bdim = min(nl,nr)
    if(bdim >= 3)then
      bdim = 3
      allocate(b1(3,bdim),b2(3,bdim),source=0.0d0)
      b1(1:3,1) = c1(1:3,1)
      b1(1:3,2) = c1(1:3,2)
      b1(1:3,3) = c1(1:3,nr)
      b2(1:3,1) = c2(1:3,1)
      b2(1:3,2) = c2(1:3,2)
      b2(1:3,3) = c2(1:3,nl)
    else 
      bdim = 2
      allocate(b1(3,bdim),b2(3,bdim),source=0.0d0)  
      b1(1:3,1) = c1(1:3,1)
      b1(1:3,2) = c1(1:3,nr)
      b2(1:3,1) = c2(1:3,1)
      b2(1:3,2) = c2(1:3,nl)
    endif
    !--- calculate RMSD and rotation matrix
    call rmsd(bdim,b2,b1,1,U,x_center,y_center,rmsdval,.false.,g)
    !--- rotate new ligand 
    allocate(tmp(3,mol2%nat),source=0.0d0)
    tmp = matmul(U(1:3,1:3),mol2%xyz(1:3,1:mol2%nat))
    !--- and shift
    do i=1,mol2%nat
    mol2%xyz(1:3,i) = tmp(1:3,i) + shift(1:3)
    enddo
    deallocate(tmp)
    deallocate(b2,b1)
    deallocate(c2,c1)

!--dummy write    
    !open(newunit=ich,file='ligandalign.xyz')
    !write(ich,*)mol1%nat+mol2%nat+1
    !write(ich,*)
    !do i=1,mol1%nat
    !mol1%xyz(1:3,i) = mol1%xyz(1:3,i)+shift(1:3)
    !write(ich,'(1x,a2,1x,3f20.10)')i2e(mol1%at(i),'nc'),mol1%xyz(1:3,i)*bohr
    !enddo
    !do i=1,mol2%nat
    !write(ich,'(1x,a2,1x,3f20.10)')i2e(mol2%at(i),'nc'),mol2%xyz(1:3,i)*bohr
    !enddo
    !write(ich,'(1x,a2,1x,3f20.10)')i2e(10,'nc'),metal(1:3)*bohr
    !close(ich)

    return
contains
function getdist(r1,r2) result(dist)
    use iso_fortran_env, wp=>real64
    implicit none
    real(wp) :: r1(3),r2(3)
    real(wp) :: dist
    dist=0.0d0
    dist=(r1(1)-r2(1))**2 +(r1(2)-r2(2))**2 + (r1(3)-r2(3))**2
    dist=sqrt(dist)
    return
end function getdist    
subroutine estimatecenter(n,c0,refdist)
    use iso_fortran_env, wp=>real64
    implicit none
    integer :: n
    real(wp) :: c0(3,n)  
    !by convention c0(:,1) should be at (0,0,0)
    !and c0(:,n) is the atom position to be determined
    real(wp) :: refdist
    real(wp) :: center(3),cdist,distrel
    integer :: i,k,n2
    n2=n-2
    center = 0.0_wp
    if(n2>0)then
       k=1
       do i=1,n2
       k=k+1
       center(1:3) = center(1:3) + c0(1:3,k)
       enddo 
       !center = center / float(n2)
       cdist = getdist(c0(1:3,1),center)
       distrel = refdist / cdist
       center = center * (-distrel)
       c0(1:3,n) = center(1:3) 
    else !n2=0 
    !this means that the ligand is an atom and 
    !we can place the center anywhere, e.g. along an axis:
        c0(1,n) = refdist
        c0(2:3,n) = 0.0_wp
    endif    
    return
endsubroutine estimatecenter
end subroutine matchligands    

subroutine ligandtool(infile,newligand,center,oldligand)
    use iso_fortran_env, only: wp=>real64
    use crest_data
    implicit none
    character(len=*) :: infile
    character(len=*) :: newligand
    integer :: center
    integer :: oldligand

    if(oldligand==0)then
        write(*,*) 'No ligand selected in the molecule'
        write(*,'(1x,a,i0,a)') 'Ligands for central atom ',center,':'
        call extractligands(infile,center)
        return
    endif

    call exchangeligands(infile,newligand,center,oldligand)
    write(*,'(1x,a,a,a)')'New geometry written to <','newstructure.xyz','>'

    return
end subroutine ligandtool

!====================================================================!
! A hack to 'flip' hydrogens at OH (technically also SH, NH, etc.)
! currently not used.
!====================================================================!
subroutine ohflip(mol,numstruc)
    use iso_fortran_env, only: wp=>real64
    use crest_data
    use strucrd
    use zdata
    use geo
    implicit none
    type(coord) :: mol,new
    type(zmolecule) :: zmol
    integer :: numstruc
    logical,allocatable :: ohmap(:)
    integer :: i,noh,io
    integer :: theh
    real(wp) :: hpos(3),opos(3),xpos(3)
    real(wp) :: kvec(3),theta,hnew(3)
    real(wp),parameter :: pi = 3.14159265359_wp
    character(len=50) :: atmp

    write(atmp,'(15x,a)') '!hor'

    !-- get topology
    call simpletopo_mol(mol,zmol,.false.,.false.)
    allocate(ohmap(zmol%nat), source = .false.)

    !-- identify OH groups
    do i=1,zmol%nat
    ohmap(i) = isoh(zmol%zat(i))
    enddo
    noh = count(ohmap,1)

    !-- only if there are OH groups
    if(noh > 0)then
     numstruc=noh
    open(newunit=io,file='oh_tmp.xyz')
    !-- loop over all OH
    do i=1,zmol%nat
       if(ohmap(i))then
           new = mol
           new%comment = trim(atmp)
           theh = ohpos(zmol,i,hpos,opos,xpos)
           hpos = hpos - opos
           kvec = xpos - opos
           !-- rotate
           theta = pi
           call rodrot(hpos,kvec,theta)
           !-- shift back
           hnew = hpos + opos 
           new%xyz(:,theh) = hnew
           call new%append(io)
       endif    
    enddo
    close(io)
    else
        numstruc=0
    endif
    deallocate(ohmap)
    call zmol%deallocate()
    call new%deallocate()
    return
contains
function isoh(zatm) result(bool)
    use zdata
    implicit none
    type(zatom) :: zatm
    logical :: bool
    integer :: i,k
    bool =.false.
    k=0
    if(zatm%nei == 2)then
      do i=1,2 
        if(zatm%ngt(i)==1) k=k+1
      enddo  
      if(k==1) bool=.true.
    endif
    return
end function isoh    
function ohpos(zmol,k,hpos,opos,xpos) result(thehatom)
    use iso_fortran_env, only: wp=>real64
    use zdata
    implicit none
    type(zmolecule) :: zmol
    integer :: i,k,j,thehatom
    real(wp) :: hpos(3),opos(3),xpos(3)
    opos = zmol%zat(k)%cart
    do i=1,2
      j=zmol%zat(k)%ngh(i)
      if(zmol%zat(j)%atype==1)then
       hpos = zmol%zat(j)%cart
       thehatom=j
      else
       xpos = zmol%zat(j)%cart
      endif
    enddo
    return
end function ohpos
end subroutine ohflip
!====================================================================!
! file wrapper for the ohflip routine
!====================================================================!
subroutine ohflip_file(infile)
    use iso_fortran_env, only: wp=>real64
    use crest_data
    use strucrd
    use zdata
    implicit none
    character(len=*) :: infile
    type(coord) :: mol
    integer :: k
    call mol%open(infile)
    call ohflip(mol,k)
    call mol%deallocate()
    write(*,*) k,'XH groups in the molecule'
    return
end subroutine ohflip_file

subroutine ohflip_ensemble(infile,maxnew)
    use crest_parameters 
    use strucrd
    use zdata
    use iomod
    implicit none
    character(len=*) :: infile
    type(coord) :: mol
    integer,intent(in) :: maxnew
    integer :: limit,counter,xout
    integer :: nat,nall
    integer :: k
    real(wp),allocatable :: xyz(:,:,:)
    real(wp),allocatable :: er(:)
    integer,allocatable  :: at(:)
    character(len=:),allocatable :: collection

    !--- read in the ensemble parameters
      call rdensembleparam(infile,nat,nall)
    !--- allocate space and read in the ensemble
      allocate(at(nat),er(nall),xyz(3,nat,nall))
      call rdensemble(infile,nat,nall,at,xyz,er)

    !-- loop over ensemble
    limit = maxnew
    collection='oh_ensemble.xyz'
    call remove(collection) 
    k=0
    counter = 0
    do
    k=k+1
    call mol%get(bohr,nat,at,xyz(:,:,k))
    call ohflip(mol,xout)
    call mol%deallocate()
    if(xout==0) then
        exit
    else
     counter = counter + xout   
     call appendto('oh_tmp.xyz',collection)
     call remove('oh_tmp.xyz')
    endif
    if((counter>=limit).or.(k==nall))then
        exit
    endif
    enddo

    deallocate(collection,xyz,er,at)
    return
end subroutine ohflip_ensemble

