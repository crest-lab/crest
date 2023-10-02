!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020 Stefan Grimme, Philipp Pracht
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

!==============================================================!
! the handler function is the sub-program that is called
! by crest, i.e., within this routine the seperate input file is
! read and the calculation is started.
!==============================================================!
subroutine msreact_handler(env,tim)
    use crest_parameters
    use crest_data
    use msmod
    use strucrd
    implicit none

    type(systemdata) :: env
    type(timer) :: tim

    type(msobj) :: mso  !a new msreact object
    type(coord) :: struc

    interface
        subroutine msreact_topowrap(mol,pair,paths,wboname)
            import :: msmol
            implicit none
            type(msmol) :: mol
            integer :: pair(mol%nat*(mol%nat+1)/2)
            integer :: paths(mol%nat*(mol%nat+1)/2,mol%nat)
            character(len=*),optional :: wboname
        end subroutine msreact_topowrap
    end interface
    integer :: nat
    integer,allocatable :: pair(:)
    integer,allocatable :: paths(:,:)
    integer :: k

    call tim%start(1,'MSREACT')

    !-- read the input coord and put it into the
    !   iso-list as Gen 0 structure
    call struc%open('coord')
    struc%xyz=struc%xyz*bohr !to Angström, from this point on by convention!
    call mso%il%append(struc%nat,struc%at,struc%xyz,0.0_wp,env%chrg,0)
    call struc%deallocate()

    !-- additional input file could be read here
    !call msinputreader(mso)


    nat = mso%il%mol(1)%nat
    k=nat*(nat+1)/2
    allocate(pair(k),paths(k,nat))
    call msreact_topowrap(mso%il%mol(1),pair,paths,'wbo')

    !-- setting the threads for correct parallelization
    if(env%autothreads)then
        call ompautoset(env%threads,6,env%omp,env%MAXRUN,0) !set global OMP/MKL variable for xtb jobs
    endif

    !-- do the directory setup and optimizations
    call msreact(mso,mso%il%mol(1),nat,pair,3)


    deallocate(paths,pair)

    call tim%stop(1)
    return
end subroutine msreact_handler

!==============================================================!
! the main implementation of the msreact algo should go here
!==============================================================!
subroutine msreact(mso,mol,nat,pair,nbonds)
      use crest_parameters
      use msmod
      use iomod
      use miscdata, only: rcov
      use utilities, only: lin
      implicit none

      type(msobj) :: mso    !main storage object
      type(msmol) :: mol    ! xyz etc

      integer :: nat
      integer :: pair(nat*(nat+1)/2)
      !integer :: paths(nat*(nat+1)/2,nat)
      integer :: nbonds
      integer :: i,j,k
      integer :: p
      integer :: np
      integer :: io

      character(len=:),allocatable :: subdir
      character(len=40) :: pdir
      character(len=512) :: thisdir
      real(wp)             :: constr_dist

      !-- main subdirectory handling
      call getcwd(thisdir)
      subdir='MSDIR'
      io = makedir(subdir)
      call chdir(subdir)


      !-- get specific pairs
      np=0
      do p=1,nbonds    ! bonds in between
!        write(*,'(1x,a,i0,a)') '1,',p+1,' pairs'
         do i=1,nat
           do j=i,nat
              k=lin(i,j)
              if(p.eq.1.and.pair(k).eq.1.and.(mol%at(i).eq.1.or.mol%at(j).eq.1)) cycle
              if(pair(k)==p)then
                 np = np+1 
                 write(pdir,'(a,i0)')'Pair_',np
                 constr_dist = mso%cdist*(rcov(mol%at(i))+rcov(mol%at(j)))*bohr + float(p)
!                write(*,*) mol%at(i),mol%at(j),constr_dist
                 call isodir(mso,trim(pdir),mol,i,j,constr_dist)
              endif    
           enddo
         enddo
      enddo
             
      write(*,*) '# of distortions',np
      call msreact_jobber(np,'Pair_',.false.)

      call msreact_collect(mol%nat,np,'products.xyz')
      call rename(subdir//'/'//'products.xyz','products.xyz')
      call chdir(thisdir)
      return
end subroutine msreact


!============================================================!
! make a dir for a structure without fragments,
! a controlfile with constraints on atoms A and B (at dist D)
! will be written into the directory
!============================================================!
subroutine isodir(mso,dirname,mol,A,B,D)
    use crest_parameters
    use msmod
    use iomod
    use strucrd, only : wrxyz
    implicit none
    type(msobj) :: mso
    character(len=*) :: dirname
    type(msmol) :: mol
    integer :: A,B
    real(wp) :: D

    character(len=:),allocatable :: fname
    character(len=20) :: dumm
    integer :: io,ich

    io = makedir(dirname) !create the directory
    
    fname = trim(dirname)//'/'//'struc.xyz'
    open(newunit=ich,file=fname)
    call wrxyz(ich,mol%nat,mol%at,mol%xyz)
    close(ich)

    fname = trim(dirname)//'/'//'.CHRG'
    open(newunit=ich,file=fname)
    write(ich,'(i0)') mol%chrg + 1   ! EI +1, DEA -1, CID 0
    close(ich)

    fname = trim(dirname)//'/'//'.xc1'
    open(newunit=ich, file=fname)
    write(ich,'(a)') '$scc'
    write(dumm,'(f16.2)') mso%T
    write(ich,'(1x,a,a)')'temp=',adjustl(trim(dumm))
    write(ich,'(a)')'$constrain'
    write(dumm,'(f16.4)') mso%fc
    write(ich,'(3x,a,a)')'force constant=',adjustl(trim(dumm))
    write(ich,'(3x,a,1x,i0,a,1x,i0,a,1x,f8.5)') 'distance:',A,',',B,',',D
    close(ich)

    fname = trim(dirname)//'/'//'.xc2'
    open(newunit=ich, file=fname)
    write(ich,'(a)') '$scc'
    write(dumm,'(f16.2)') mso%T
    write(ich,'(1x,a,a)')'temp=',adjustl(trim(dumm))
    write(ich,'(a)') '$opt'
    write(ich,'(1x,a)') 'maxcycle=5'
    write(ich,'(a)') '$write'
    write(ich,'(1x,a)') 'wiberg=true'
    close(ich)

    return
end subroutine isodir

!=====================================================================!
! The job construction routine for MSREACT
! (will have to be modified later, for now it is for testing)
!=====================================================================!
subroutine msreact_jobber(ndirs,base,niceprint)
    use crest_parameters
    use msmod
    use iomod
    implicit none
    integer :: ndirs
    character(len=*) :: base
    logical :: niceprint

    character(len=1024) :: jobcall
    character(len=1024) :: jobcall2

    jobcall = ''
    jobcall2 = ''

    write(jobcall,'(a)')  'xtb struc.xyz  --opt loose --input .xc1 > split.out 2>/dev/null'
    write(jobcall2,'(a)') 'xtb xtbopt.xyz --opt crude --input .xc2 > xtb.out 2>/dev/null'
    jobcall = trim(jobcall)//' ; '//trim(jobcall2)

    !-- directories must be numbered consecutively
    call opt_OMP_loop(ndirs,base,jobcall,niceprint)
    write(*,*)
    write(*,*) 'done.'
    return
end subroutine msreact_jobber



!=====================================================================!
! A wrapper to generate the topology for a molecule within the
! MSREACT subprogram
!=====================================================================!
subroutine msreact_topowrap(mol,pair,paths,wboname)
    use crest_parameters
    use msmod
    use zdata
    use adjacency
    use utilities, only: lin
    implicit none
    type(msmol) :: mol
    integer :: pair(mol%nat*(mol%nat+1)/2)
    !integer :: pair(mol%nat,mol%nat)
    integer :: paths(mol%nat*(mol%nat+1)/2,mol%nat)
    character(len=*),optional :: wboname
    type(zmolecule) :: zmol

    integer,allocatable :: A(:,:)
    integer,allocatable :: prev(:,:)
    real(wp),allocatable :: E(:,:)
    real(wp),allocatable :: dist(:,:)

    integer :: lpath,i,j,k
    integer,allocatable :: path(:)
    logical :: ex


    ex=.false.
    if(present(wboname))then
      inquire(file=wboname,exist=ex)
    endif  
    if(ex)then
      call simpletopo(mol%nat,mol%at,mol%xyz,zmol,.false.,.false.,wboname)
    else   
       mol%xyz = mol%xyz / bohr !CN based topo requires Bohrs 
       call simpletopo(mol%nat,mol%at,mol%xyz,zmol,.false.,.false.,'')
       mol%xyz = mol%xyz * bohr
    endif 

    allocate(A(mol%nat,mol%nat),E(mol%nat,mol%nat))
    call zmol%adjacency(A,E)

    allocate(prev(mol%nat,mol%nat),dist(mol%nat,mol%nat))

    call FloydWarshall(mol%nat,A,E,dist,prev)
    allocate(path(mol%nat), source = 0)
    do i=1,mol%nat
      do j=i,mol%nat
       path = 0
       call getPathFW(mol%nat,prev,i,j,path,lpath)
       !write(*,*) path(1:lpath)
       k=lin(i,j)
       pair(k) = lpath - 1 ! number of bonds
       paths(k,:) = path(:)
      enddo
     enddo

    deallocate(dist,prev)
    deallocate(E,A)

    call zmol%deallocate() !clear the zmol memory
    return
end subroutine msreact_topowrap    

!========================================================================!
! collect structures of optimized molecules
! xyz files should still have the same number and order of atoms
!========================================================================!
subroutine msreact_collect(nat,np,outfile)
    use crest_parameters
    use strucrd
    implicit none
    integer :: nat
    integer :: np
    character(len=*) :: outfile
    integer :: ich
    character(len=40) :: pdir
    character(len=:),allocatable :: optfile
    character(len=128) :: newcomment
    integer :: p,p2
    logical :: ex
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    real(wp) :: etot


    allocate(at(nat),xyz(3,nat))
    open(newunit=ich,file=outfile)
    p=0
    do p2=1,np
       write(pdir,'(i0,i0,a,i0)')1,p+1,'Pair_',p2
       write(pdir,'(a,i0)')'Pair_',p2
       optfile=trim(pdir)//'/'//'xtbopt.xyz'
       inquire(file=optfile,exist=ex)
       if(ex)then
          call rdcoord(optfile,nat,at,xyz,etot)
          xyz = xyz*bohr
          write(newcomment,'(1x,f18.8,5x,a)')etot,trim(pdir)
          call wrxyz(ich,nat,at,xyz,newcomment)       
       endif    
    enddo
    close(ich)

    deallocate(xyz,at)
    return
end subroutine msreact_collect
