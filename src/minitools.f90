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

!========================================================!
! A small subroutine to split an ensemble file into
! seperate directories. A max. number of files can be
! be specified.
!========================================================!
subroutine splitfile(fname,up,low)
      use iso_fortran_env, wp => real64
      use iomod
      use strucrd, only: rdensembleparam,rdensemble,wrxyz
      implicit none
      character(len=*) :: fname
      integer :: up,low
      integer :: nref

      character(len=512) :: thispath,tmppath1,tmppath2
      character(len=128) :: atmp,btmp

      real(wp),allocatable :: xyz(:,:,:)
      integer :: nat,nall
      integer :: nc
      integer,allocatable :: at(:)
      integer :: i,r
      logical :: ex

      
      inquire(file=fname,exist=ex)
      if(.not.ex)then
          write(0,'(a,a,a)') "file ",trim(fname)," does not exist. must stop"
          error stop
      endif

      call getcwd(thispath) !current dir= thispath

      call rdensembleparam(fname,nat,nall)
      allocate(xyz(3,nat,nall),at(nat))
      call rdensemble(fname,nat,nall,at,xyz)

      r = makedir("SPLIT")  !create new directory
      call chdir("SPLIT")
      call getcwd(tmppath1)

      if(up.gt.0 .and. low.gt.0)then
         if(low.gt.up)then 
            i   = low
            low = up
            up  = i
         endif    
      endif    
      ! if there are less structures in the file than we need
      if( nall < up ) then
         nc = nall
      else
         nc = up
      endif

      do i=low,nc
         write(tmppath2,'(a,i0)')"STRUC",i
         r = makedir(trim(tmppath2))
         call chdir(tmppath2)
         call wrxyz("struc.xyz",nat,at,xyz(:,:,i))
         call chdir(tmppath1)
      enddo

      call chdir(thispath)
      return
end subroutine splitfile



!========================================================!
! Temporary helper function for 
! trying out the -gfn2@gfnff flag
! Use only x percent of rotamers for each
!========================================================!
subroutine prepentropy(env,fname,percent)
      use iso_fortran_env, wp => real64
      use crest_data
      use iomod
      use strucrd, only: rdensembleparam,rdensemble,wrxyz
      implicit none
      !type(options) :: opt
      type(systemdata) :: env
      character(len=*) :: fname
      real(wp) :: percent
      integer :: nref

      character(len=512) :: thispath,tmppath1,tmppath2
      character(len=128) :: atmp,btmp

      real(wp),allocatable :: xyz(:,:,:)
      integer :: nat,nall
      integer :: nc
      integer,allocatable :: at(:)
      integer,allocatable :: degen(:,:)
      character(len=128),allocatable :: comment(:)
      integer :: ng,ns
      real(wp) :: nsf
      integer :: i,j,k,l,r
      logical,allocatable :: incl(:)
      integer :: ich
      logical :: ex


      inquire(file='cre_members',exist=ex)
      if(ex)then
          open(newunit=ich,file='cre_members')
          read(ich,*)ng
          allocate(degen(3,ng))
          do i=1,ng
             read(ich,*) degen(1:3,i)
          enddo
          close(ich)

          call rename('cre_members','cre_members.backup')
      else
          error stop 'file cre_members not found!'
      endif

      call rdensembleparam(fname,nat,nall)
      allocate(xyz(3,nat,nall),at(nat),incl(nall),comment(nall))
      call rdensemble(fname,nat,nall,at,xyz,comment)

      incl = .false.
      do i=1,ng
         nsf = float(degen(1,i))*percent
         if( nsf .gt. degen(1,i)) nsf = degen(1,i)
         if( nsf .lt. 1.0_wp) nsf = 1.0_wp
         ns = nint(nsf)
         k = degen(2,i)
         l = k + ns - 1
         do j=k,l
            incl(j) = .true.
         enddo
      enddo
      write(*,*) incl

      write(*,*) 'Writing compacted ensemble file'
      open(newunit=ich,file='crest_compact.xyz')
      do i=1,nall
         if(incl(i))then
           call wrxyz(ich,nat,at,xyz(:,:,i),comment(i))
         endif
      enddo


      env%doNMR = .true.
      env%cgf(3)=.true.
      env%confgo = .true.
      env%ensemblename = 'crest_compact.xyz'
      if(env%newcregen)then
       call newcregen(env,0)
      else
       call cregen2(env)
      endif
      env%confgo = .false.
      env%ensemblename = 'crest_compact.xyz.sorted'

      call rename('cre_members.backup','cre_members')
      if(allocated(at))deallocate(at)
      if(allocated(incl))deallocate(incl)
      if(allocated(xyz))deallocate(xyz)
      if(allocated(degen))deallocate(degen)
      return
end subroutine prepentropy      


!- print the anisotropy of the rotational constants
!  for all structures in a given ensemble file
subroutine printaniso(fname,bmin,bmax,bshift)
    use iso_fortran_env, wp => real64
    use strucrd
    use axis_module
    implicit none
    character(len=*) :: fname
    type(ensemble) :: ens

    integer :: nat
    integer :: nall
    real(wp),allocatable :: c1(:,:)
    integer,allocatable :: at(:)

    real(wp),allocatable :: rot(:,:)
    real(wp) :: rotaniso !function
    real(wp),allocatable :: anis(:)

    real(wp) :: bthrerf
    real(wp) :: bmin,bmax,bshift
    real(wp) :: thr
    real(wp) :: dum
    integer :: i,j,k,l

    call ens%open(fname)
    nat=ens%nat
    nall=ens%nall

    allocate(c1(3,nat),at(nat))
    allocate(rot(3,nall))
    allocate(anis(nall))

    at = ens%at

    do i=1,nall
      c1(1:3,:) = ens%xyz(1:3,:,i)
      call axis(nat,at,c1,rot(1:3,i),dum)
      anis(i) = rotaniso(i,nall,rot)
      thr =  bthrerf(bmin,anis(i),bmax,bshift)
      write(*,'(3f10.2,2x,f8.4,2x,f8.4)') rot(1:3,i),anis(i),thr
    enddo

    deallocate(anis,rot,at,c1)

    stop
    return
end subroutine printaniso

!--- read in a file with 1 to 2 columns
!    first column is an energy (in Eh)
!    second column can be a degeneracy
!    a Boltzman weight is calculated for each
!    energy and a average ensemble energy is returned
!    comment lines (#) are ignored
subroutine prbweight(fname,Targ)
    use iso_fortran_env, wp => real64
    implicit none

    character(len=*) :: fname
    character(len=*) :: Targ

    real(wp) :: dum,T
    logical :: ex
    integer :: io,ich

    character(len=128) :: atmp
    integer :: l,j,i
    real(wp),allocatable :: elist(:)
    real(wp),allocatable :: erel(:)
    real(wp),allocatable :: g(:)
    real(wp),allocatable :: p(:)
    real(wp) :: xx(10)
    real(wp) :: eav,elow

    !-- parse temperature from argument
    Targ = trim(adjustl(Targ))
    read(Targ,*,iostat=io) dum
    if(io == 0)then
        T = dum
    else
        T=298.15d0
    endif

    inquire(file=fname,exist=ex)
    if(.not.ex)stop

    open(newunit=ich,file=fname)
    l = 0
    do 
      read(ich,'(a)',iostat=io) atmp
      if(io < 0) exit !EOF
      atmp = trim(adjustl(atmp))
      if(len_trim(atmp).lt.1) cycle
      if(atmp(1:1) == '#') cycle
      l = l + 1
    enddo
    close(ich)

    allocate(elist(l),erel(l), source = 0.0_wp)
    allocate(g(l) , source = 1.0_wp)
    allocate(p(l), source = 0.0_wp)

    open(newunit=ich,file=fname)
    l = 0
    do
      read(ich,'(a)',iostat=io) atmp
      if(io < 0) exit !EOF
      atmp = trim(adjustl(atmp))
      if(len_trim(atmp).lt.1) cycle
      if(atmp(1:1) == '#') cycle
      l = l + 1
      call readl(atmp,xx,j)
      elist(l) = xx(1)
      if(j > 1) g(l) = nint(xx(2))
    enddo
    close(ich)

    dum=minval(elist,1)
    elow=dum
    do i=1,l
      erel(i) = (elist(i) - dum)*627.5905_wp
    enddo
    call entropy_boltz(l,T,erel,g,p)

    atmp=trim(fname)//'.out'
    open(newunit=ich,file=atmp)
    write(ich, '(a,16x,a,2x,a)')'#Etot','degen','pop'
    write(  *, '(a,16x,a,2x,a)')'#Etot','degen','pop'
    do i=1,l
      write(ich,'(f20.8,1x,i0,4x,f10.6)') elist(i),nint(g(i)),p(i)
      write(  *,'(f20.8,1x,i0,4x,f10.6)') elist(i),nint(g(i)),p(i)
    enddo
    close(ich)

    write(*,*)
    write(*,'(1x,a,f6.4)') 'pop sum: ',sum(p)
    do i=1,l
       elist(i) = elist(i)*p(i)
    enddo

    write(*,*) 'E_low:',elow
    write(*,*) 'E_av:',sum(elist)
    stop
end subroutine prbweight

subroutine calceav(fname,T,eav,verbose)
    use iso_fortran_env, wp => real64
    use strucrd
    implicit none
    character(len=*) :: fname
    real(wp) :: eav
    real(wp) :: T
    logical :: verbose
    type(ensemble) :: ens
    real(wp),allocatable :: elist(:)
    real(wp),allocatable :: erel(:)
    real(wp),allocatable :: g(:)
    real(wp),allocatable :: p(:)
    real(wp) :: dum
    integer :: i,j
    call ens%open(fname)
    allocate(elist(ens%nall),erel(ens%nall),g(ens%nall),p(ens%nall))
    elist = ens%er
    erel = 0.0d0
    g = 1.0d0
    p = 0.0d0
    if(verbose)then
        write(*,'(1x,a,1x,a)') 'Energies for file:',trim(fname)
        do i=1,ens%nall
          write(*,'(1x,i6,1x,f20.10)') i,elist(i)
        enddo
    endif
    dum=minval(elist,1)
    do i=1,ens%nall
      erel(i) = (elist(i) - dum)*627.5905_wp
    enddo
    call entropy_boltz(ens%nall,T,erel,g,p)
    do i=1,ens%nall
       elist(i) = elist(i)*p(i)
    enddo
    eav=sum(elist)
    if(verbose)then
        write(*,'(1x,a,1x,f20.10)') 'Weighted energy:',eav
    endif
    deallocate(p,g,erel,elist)
    return
end subroutine calceav

subroutine getelow(fname,elow,verbose)
    use iso_fortran_env, wp => real64
    use strucrd
    implicit none
    character(len=*) :: fname
    real(wp) :: elow
    logical :: verbose
    type(ensemble) :: ens
    real(wp),allocatable :: elist(:)
    real(wp),allocatable :: erel(:)
    real(wp) :: dum
    integer :: i,j

    call ens%open(fname)
    allocate(elist(ens%nall))
    elist = ens%er
    if(verbose)then
        write(*,'(1x,a,1x,a)') 'Energies for file:',trim(fname)
        do i=1,ens%nall
          write(*,'(1x,i6,1x,f20.10)') i,elist(i)
        enddo
    endif
    elow=minval(elist,1)
    if(verbose)then
        write(*,'(1x,a,1x,f20.10)') 'Lowest energy:',elow
    endif
    deallocate(elist)
    call ens%deallocate()
    return
end subroutine getelow

!=================================================================================================!
! set up the topology for a file and analyze it
!=================================================================================================!
subroutine testtopo(fname,env,tmode)
    use iso_fortran_env, wp => real64
    use crest_data
    use iomod
    use atmasses
    use zdata
    implicit none
    !type(options) :: opt
    type(systemdata) :: env
    character(len=*) :: fname
    character(len=:),allocatable :: wbofile
    character(len=*) :: tmode
    character(len=40) :: sumform
    type(zmolecule) :: zmol
    real(wp),allocatable :: xyz(:,:)
    real(wp) :: dum
    integer,allocatable :: inc(:)
    real(wp) :: flex
    real(wp) :: rabc(3),avmom
    real(wp) :: molmass,symnum
    character(len=3) :: symchar
    integer :: nt,i,j,k,l
    logical :: l1
    real(wp),allocatable :: temps(:)
    real(wp),allocatable :: et(:)
    real(wp),allocatable :: ht(:)
    real(wp),allocatable :: gt(:)
    real(wp),allocatable :: stot(:)
 
    call to_lower(tmode)
!--- specify zmol sertup
    select case( tmode )
      case( 'wbo','flexi','all' )
        call xtbsp2(fname,env)  
        call simpletopo_file(fname,zmol,.true.,.true.,'wbo')
      case default
        wbofile = 'none'
        call  simpletopo_file(fname,zmol,.true.,.true.,wbofile)
    end select    
    allocate(xyz(3,zmol%nat))
    call zmol%getxyz(xyz)
    xyz=xyz*bohr   !to angstroem
!--- specify other analysis
     write(*,*)
     select case( tmode )
       case( 'sym','symmetry' )
         call analsym(zmol,dum,.true.)
       case( 'flexi' )
         allocate(inc(zmol%nat), source=1)  
         call flexi(zmol%nat,zmol%nat,inc,flex,dum) 
         write(*,'(1x,a,4x,f6.4)') 'flexibility measure:',flex
         deallocate(inc)
       case( 'zmat' )
         call ztopozmat(zmol,.true.)  
       case( 'formula','sumform' )
         write(*,'(/,1x,a)') trim(sumform(zmol%nat,zmol%at))
         write(*,'(1x,a,i16)') '# atoms: ',zmol%nat
         write(*,'(1x,a,f16.5)') 'Mol.weight: ',molweight(zmol%nat,zmol%at)
       case( 'all' )  
          call ztopozmat(zmol,.true.) 
          write(*,'(/,1x,a)') trim(sumform(zmol%nat,zmol%at))
          write(*,'(1x,a,i16)') '# atoms: ',zmol%nat
          write(*,'(1x,a,f16.5)') 'Mol.weight: ',molweight(zmol%nat,zmol%at)
          call analsym(zmol,dum,.true.)
          allocate(inc(zmol%nat), source=1)
          call flexi(zmol%nat,zmol%nat,inc,flex,dum)
          write(*,'(1x,a,4x,f6.4)') 'flexibility measure:',flex
          deallocate(inc)
       case('thermo')
          !call prepthermo(zmol%nat,zmol%at,xyz,.true., &
          !&    molmass,rabc,avmom,symnum,symchar)    
          if(.not.allocated(env%thermo%temps))then
             call env%thermo%get_temps()
          endif
          nt=env%thermo%ntemps
          allocate(temps(nt),et(nt),ht(nt),gt(nt),stot(nt))
          temps = env%thermo%temps
          call thermo_wrap(env,.true.,zmol%nat,zmol%at,xyz,'', &
          &    nt,temps,et,ht,gt,stot,.false.) 
          deallocate(stot,gt,ht,et,temps)
       case( 'methyl' )
           do i=1,zmol%nat
           l1=zmol%methyl(i)
           !write(*,*) l1
           if(l1) write(*,'(a,i0,a)') 'Atom ',i,' is methyl (or similar)'
           enddo  
       case( 'stereo' )    
        call isstereo(zmol)
       end select
    deallocate(xyz) 
    write(*,*)
    stop
end subroutine testtopo

!===================================================!
! get sumformula as a string from the AT array
!===================================================!
character(len=40) function sumform(nat,at)
   use strucrd, only: i2e
   implicit none
   integer :: nat
   integer :: at(nat)
   integer :: sumat(94)
   integer :: i
   character(len=6) :: str
   sumform=''
   sumat=0
   do i=1,nat
     sumat(at(i))=sumat(at(i))+1
   enddo
   do i=1,94
     if(sumat(i).lt.1)cycle
     write(str,'(a,i0)')trim(adjustl(i2e(i,'nc'))),sumat(i)
     sumform=trim(sumform)//trim(str)
   enddo
   return
end function sumform

!=================================================================!
! read an ensemble and determine the symmetry for all structures
!=================================================================!
subroutine ensemble_analsym(fname,pr)
       use iso_fortran_env, wp => real64
       use crest_data, only: bohr
       use strucrd
       implicit none
       character(len=*) :: fname
       logical :: pr
       integer :: nat,nall
       real(wp),allocatable :: xyz(:,:,:)
       real(wp),allocatable :: c0(:,:)
       real(wp),allocatable :: er(:)
       integer,allocatable  :: at(:)
       integer :: i,j,k,io,ich
       character(len=4) :: sfsym,sfsm
       real(wp),parameter :: desy = 0.1_wp
       integer,parameter  :: maxat = 200
       character(len=80) :: atmp

       call rdensembleparam(fname,nat,nall)
       !--- allocate space and read in the ensemble
       allocate(at(nat),xyz(3,nat,nall),c0(3,nat),er(nall))
       call rdensemble(fname,nat,nall,at,xyz,er)


       if(pr)then
           write(*,*)
           call smallhead('STRUCTURE SYMMETRIES')
           write(*,'(1x,a)') 'Unlisted structures have symmetry C1'
           write(*,'(1x,a)') 'The full list can be found in the file "symmetries"'
           write(*,*)
           write(*,'(12x,a10,2x,a18,2x,a)') 'number','energy/Eh','sym.'
       endif

       xyz = xyz/bohr
       open(file='symmetries',newunit=ich)
       do i=1,nall
          c0(1:3,1:nat) = xyz(:,:,i) 
          call getsymmetry2(.false.,6,nat,at,c0,desy, maxat, sfsym)
          sfsm=sfsym(1:3)
          write(atmp,'(3x,a,i10,2x,f18.8,2x,a)') 'structure',i,er(i),sfsm
          write(ich,'(a)')trim(atmp)
          if(pr)then
            if(trim(sfsm) /= "c1")then
              write(*,'(a)') trim(atmp)
            endif
          endif
       enddo
       close(ich)

       return
end subroutine ensemble_analsym

!=========================================================================!
!
!=========================================================================!
function quick_rmsd(fname,nat,at,xyz,heavy) result(rout)
    use iso_fortran_env, only: wp=>real64,error_unit
    use ls_rmsd
    use strucrd
    use crest_data, only: bohr
    implicit none
    character(len=*) :: fname
    integer :: nat
    integer :: at(nat)
    real(wp) :: xyz(3,nat)
    logical :: heavy
    real(wp) :: rout
    type(coord) :: mol
    real(wp),allocatable :: c0(:,:)
    real(wp),allocatable :: c1(:,:)
    integer :: i,j,k,l
    real(wp),allocatable :: gdum(:,:),Udum(:,:),xdum(:),ydum(:)
    rout = 0.0_wp

    call mol%open(fname)

    if(mol%nat .ne. nat)then
        write(error_unit,*)'dimension mismatch in quick_rmsd()'
        return
    endif
    do i=1,nat
     if(at(i) .ne. mol%at(i))then
        write(error_unit,*)'atom order mismatch in quick_rmsd()'
        return
     endif
    enddo
    if(heavy)then !count heavy atoms
      j=0
      do i=1,nat
       if(at(i).ne.1)then
        j = j +1
       endif
      enddo  
      k=j  
    else
      k=nat
    endif
    allocate(c0(3,k),c1(3,k), source=0.0_wp)
    l=0
    do i=1,nat
      if(heavy .and. at(i)==1)cycle
      l=l+1
      c0(1:3,l) = mol%xyz(1:3,i)*bohr !mol coordinates are in bohr, need conversion
      c1(1:3,l) = xyz(1:3,i) !xyz coordinates are in Ang
    enddo
    allocate(gdum(3,3),Udum(3,3),xdum(3),ydum(3))
    call rmsd(k,c0,c1,0,Udum,xdum,ydum,rout,.false.,gdum)
    deallocate(ydum,xdum,Udum,gdum)
    deallocate(c1,c0)
    return
end function quick_rmsd

subroutine quick_rmsd_tool(fname1,fname2,heavy)
    use iso_fortran_env, only: wp => real64
    use strucrd
    use crest_data, only: bohr
    implicit none
    character(len=*) :: fname1
    character(len=*) :: fname2
    logical :: heavy
    type(coord) :: mol1
    real(wp) :: rmsdval
    real(wp) :: quick_rmsd


    call mol1%open(fname1)

    mol1%xyz =mol1%xyz * bohr !to Angstroem

    rmsdval = quick_rmsd(fname2,mol1%nat,mol1%at,mol1%xyz,heavy)

    if(heavy)then
    write(*,'(1x,a,f16.8)') 'Calculated heavy atom RMSD (Å):',rmsdval
    else
    write(*,'(1x,a,f16.8)') 'Calculated RMSD (Å):',rmsdval
    endif

    return 
end subroutine quick_rmsd_tool

function quick_rmsd2(nat,at,xyz,xyz2,heavy) result(rout)
    use iso_fortran_env, only: wp=>real64,error_unit
    use ls_rmsd
    use crest_data, only: bohr
    implicit none
    integer :: nat
    integer :: at(nat)
    real(wp) :: xyz(3,nat)  !in Angstroem
    real(wp) :: xyz2(3,nat) !in Amgstroem
    logical :: heavy
    real(wp) :: rout
    real(wp),allocatable :: c0(:,:)
    real(wp),allocatable :: c1(:,:)
    integer :: i,j,k,l
    real(wp),allocatable :: gdum(:,:),Udum(:,:),xdum(:),ydum(:)
    rout = 0.0_wp
    if(heavy)then !count heavy atoms
      j=0
      do i=1,nat
       if(at(i).ne.1)then
        j = j +1
       endif
      enddo  
      k=j  
    else
      k=nat
    endif
    allocate(c0(3,k),c1(3,k), source=0.0_wp)
    l=0
    do i=1,nat
      if(heavy .and. at(i)==1)cycle
      l=l+1
      c0(1:3,l) = xyz(1:3,i)
      c1(1:3,l) = xyz2(1:3,i)
    enddo
    allocate(gdum(3,3),Udum(3,3),xdum(3),ydum(3))
    call rmsd(k,c0,c1,0,Udum,xdum,ydum,rout,.false.,gdum)
    deallocate(ydum,xdum,Udum,gdum)
    deallocate(c1,c0)
    return
end function quick_rmsd2
 


!===============================================================================!
! calculate thermostatisitical contributions for a given
! input file and "vibspectrum" file in the turbomole formate
!===============================================================================!
subroutine thermo_mini(env)
    use iso_fortran_env, only: wp=>real64,error_unit
    use strucrd
    use crest_data
    implicit none
    type(systemdata) :: env
    type(coord) :: mol
    real(wp) :: ithr,fscal,sthr
    integer :: nt
    real(wp),allocatable :: temps(:)
    real(wp),allocatable :: et(:)
    real(wp),allocatable :: ht(:)
    real(wp),allocatable :: stot(:)
    real(wp),allocatable :: gt(:)
    integer :: nfreq
    real(wp),allocatable :: freq(:)
    real(wp) :: etot
    logical :: ex

    if(.not.allocated(env%thermo%temps))then
       call env%thermo%get_temps()
    endif
    nt=env%thermo%ntemps
    allocate(temps(nt),et(nt),ht(nt),gt(nt),stot(nt), source=0.0_wp)
    temps = env%thermo%temps

    call mol%open(env%inputcoords)
    etot = mol%energy
    nfreq = 3*mol%nat
    mol%xyz = mol%xyz*bohr !to Angstroem, important!

    allocate(freq(nfreq))
    inquire(file='vibspectrum',exist=ex)
    if(.not.ex)then
        error stop 'vibspectrum file does not exist!'
    endif
    call rdfreq('vibspectrum',nfreq,freq)

    ithr=env%thermo%ithr
    fscal=env%thermo%fscal
    sthr=env%thermo%sthr
    call calcthermo(mol%nat,mol%at,mol%xyz,etot,freq,.true.,ithr,fscal,sthr, &
                    &    nt,temps,et,ht,gt,stot )
    deallocate(freq)
    deallocate(stot,gt,ht,et,temps)
    return
end subroutine thermo_mini



!===============================================================================!
! resort all structures of a given ensemblefile
!===============================================================================!
subroutine resort_ensemble(fname)
    use iso_fortran_env, only: wp=>real64,error_unit
    use strucrd
    use crest_data
    implicit none
    character(len=*) :: fname
    integer :: nat,nall
    real(wp),allocatable :: xyz(:,:,:)
    integer,allocatable :: at(:)
    character(len=128), allocatable :: comm(:)
    logical :: ex
    integer :: i,j,k,ich

    integer,allocatable :: atorder(:)


    inquire(file='.atorder',exist=ex)
    if(.not.ex) error stop 'file .atorder does not exist'

    call rdensembleparam(fname,nat,nall)
    !--- allocate space and read in the ensemble
    allocate(at(nat),xyz(3,nat,nall),comm(nall))
    call rdensemble(fname,nat,nall,at,xyz,comm)

    allocate(atorder(nat))
    open(newunit=ich,file='.atorder')
    do i=1,nat
      read(ich,*) j,k
      atorder(j) = k
    enddo
    close(ich)


    open(newunit=ich, file='ensemble_tmp.xyz')
    do i=1,nall
     write(ich,'(2x,i0)')nat
     write(ich,'(a)') trim(comm(i))
     do j=1,nat
      k = atorder(j) 
      write(ich,'(1x,a2,1x,3f20.10)')i2e(at(k),'nc'),xyz(1:3,k,i)
     enddo
    enddo
    close(ich)

    return
end subroutine resort_ensemble
