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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  property calculations >>> THERMO
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine prepthermo(nat,at,xyz,pr,molmass,rabc,avmom,symnum,symchar)
    use iso_fortran_env, wp => real64
    use atmasses,only : molweight
    use iomod, only: to_lower
    use crest_data, only : bohr
    use axis_module
    implicit none
    integer,intent(in)     :: nat
    integer,intent(in)     :: at(nat)
    real(wp),intent(inout) :: xyz(3,nat) !in Angstroem
    logical,intent(in)     :: pr
    real(wp),intent(out)   :: molmass
    real(wp),intent(inout) :: rabc(3) 
    real(wp),intent(out)   :: avmom
    real(wp),intent(out)   :: symnum 

    real(wp) :: a,b,c
    character(len=4) :: sfsym
    character(len=3) :: sym,symchar
    real(wp),parameter :: desy = 0.1_wp
    integer,parameter  :: maxat = 200

    !--- molecular mass in amu
    molmass = molweight(nat,at)

    if(pr)then
        write(*,'(a,f8.2)')'Mol. weight /amu  : ',molmass
    endif
    !--- rotational constants in cm-1
    rabc = 0.0d0
    call axis(nat,at,xyz,rabc(1:3),avmom)
    if(pr)then
        write(*,'(a,3f8.2)')'Rot. const. /MHz  : ',rabc(1:3)
    endif
    rabc = rabc/2.99792458d+4   ! MHz to cm-1
    a = rabc(1)
    b = rabc(2)
    c = rabc(3)
    if(pr)then
        write(*,'(a,3f8.2)')'Rot. const. /cm-1 : ',rabc(1:3)
    endif
    !--- symmetry number from rotational symmetry
    xyz = xyz/bohr
    call getsymmetry2(.false.,6,nat,at,xyz, desy, maxat, sfsym)
    xyz = xyz*bohr
    sym=sfsym(1:3)
    symchar=sym
    symnum=1.0d0
    if(a.lt.1.d-9.or.b.lt.1.d-9.or.c.lt.1.d-9)then
      if(index(sym,'d').ne.0)symnum=2.0d0
    else
      call to_lower(sym)
      if(index(sym,'c2').ne.0)symnum=2.0d0  
      if(index(sym,'s4').ne.0)symnum=2.0d0
      if(index(sym,'c3').ne.0)symnum=3.0d0
      if(index(sym,'s6').ne.0)symnum=3.0d0
      if(index(sym,'c4').ne.0)symnum=4.0d0
      if(index(sym,'s8').ne.0)symnum=4.0d0
      if(index(sym,'c5').ne.0)symnum=5.0d0
      if(index(sym,'c6').ne.0)symnum=6.0d0
      if(index(sym,'c7').ne.0)symnum=7.0d0
      if(index(sym,'c8').ne.0)symnum=8.0d0
      if(index(sym,'c9').ne.0)symnum=9.0d0
      if(index(sym,'d2').ne.0)symnum=4.0d0
      if(index(sym,'d3').ne.0)symnum=6.0d0
      if(index(sym,'d4').ne.0)symnum=8.0d0
      if(index(sym,'d5').ne.0)symnum=10.0d0 
      if(index(sym,'d6').ne.0)symnum=12.0d0
      if(index(sym,'d7').ne.0)symnum=14.0d0
      if(index(sym,'d8').ne.0)symnum=16.0d0
      if(index(sym,'d9').ne.0)symnum=18.0d0
      if(index(sym,'t' ).ne.0)symnum=12.0d0
      if(index(sym,'td').ne.0)symnum=12.0d0
      if(index(sym,'th').ne.0)symnum=12.0d0
      if(index(sym,'o' ).ne.0)symnum=24.0d0
      if(index(sym,'oh').ne.0)symnum=24.0d0
      if(index(sym,'ih').ne.0)symnum=60.0d0
    endif

    if(pr)then
          write(*,'(a,4x,a)') 'Symmetry:',sym
     endif
   return
end subroutine prepthermo

!====================================================================!
! calcthermo based on xtb's "print_thermo" routine
!====================================================================!
subroutine calcthermo(nat,at,xyz,freq,pr,ithr,fscal,sthr,nt,temps, &
    &      et,ht,gt,stot )
    use iso_fortran_env, wp => real64, iunit=>output_unit
    use crest_thermo
    use atmasses,only : molweight
    use iomod, only: to_lower
    use crest_data, only : bohr
    implicit none
    integer,intent(in)     :: nat
    integer,intent(in)     :: at(nat)
    real(wp),intent(inout) :: xyz(3,nat)  !in Angstroem
    real(wp),intent(inout) :: freq(3*nat) !in cm-1
    logical,intent(in)     :: pr
    real(wp),intent(in) :: ithr     !imag. inv. in cm-1
    real(wp),intent(in) :: fscal    !freq scaling
    real(wp),intent(in) :: sthr     !rotor cut
    integer,intent(in)  :: nt
    real(wp),intent(in) :: temps(nt)
    real(wp) :: et(nt)          !< enthalpy in Eh
    real(wp) :: ht(nt)          !< enthalpy in Eh
    real(wp) :: gt(nt)          !< free energy in Eh
    real(wp) :: stot(nt)        !< entropy in cal/molK
    real(wp) :: ts(nt)          !< entropy*T in Eh
    real(wp) :: rabc(3),a,b,c
    real(wp) :: avmom
    real(wp) :: molmass
    real(wp) :: sym
    real(wp) :: zp
    character(len=3) :: symchar
    logical :: pr2
    logical :: linear =.false.
    logical :: atom=.false.
    integer :: nvib_theo
    integer :: nvib,nimag
    real(wp) :: vibthr
    real(wp),allocatable :: vibs(:)

    integer :: i,j
    integer :: n3,rt
    real(wp) :: adum(nt)
    character(len=64) :: atmp

    character(len=*),parameter :: outfmt = &
      '(9x,"::",1x,a,f24.12,1x,a,1x,"::")'
    character(len=*),parameter :: dblfmt = &
      '(10x,":",2x,a,f24.7,1x,a,1x,":")'
    character(len=*),parameter :: intfmt = &
      '(10x,":",2x,a,i24,       6x,":")'
    character(len=*),parameter :: chrfmt = &
      '(10x,":",2x,a,a24,       6x,":")'

    real(wp),parameter :: autorcm = 219474.63067_wp
    real(wp),parameter :: rcmtoau = 1.0_wp/autorcm
    real(wp),parameter :: autocal = 627.50947428_wp*1000.0_wp


    call prepthermo(nat,at,xyz,pr,molmass,rabc,avmom,sym,symchar)

    n3=3*nat
    allocate( vibs(n3))
    vibthr = 1.0
    a = rabc(1)
    b = rabc(2)
    c = rabc(3)

    nvib_theo=3*nat-6
    if(c.lt.1.d-10) linear=.true.
    if(linear) nvib_theo=3*nat-5
    
    if(a+b+c.lt.1.d-6)then
      atom=.true.
      nvib=0
      nvib_theo=0
    endif

    nvib=0
    vibs = 0.0
    do i=1,n3
      if(abs(freq(i)).gt.vibthr)then
         nvib=nvib+1
         vibs(nvib)=freq(i)
      endif
    enddo
    ! scale
    vibs(1:nvib)=vibs(1:nvib)*fscal
 
    ! invert imaginary modes
    nimag = 0
    do i=1,nvib
      if(vibs(i).lt.0.and.vibs(i).gt.ithr) then
         vibs(i)=-vibs(i)
         if(pr)write(*,*)'inverting freq ',i,vibs(i)  
      endif
      if(vibs(i)<0.0)then
          nimag=nimag+1
      endif
    enddo

    if(pr)then
      write(*,'(a)')
      write(*,'(10x,51("."))')
      write(*,'(10x,":",22x,a,22x,":")') "SETUP"
      write(*,'(10x,":",49("."),":")')
      write(*,intfmt) "# frequencies    ",nvib
      write(*,intfmt) "# imaginary freq.",nimag
      write(atmp,*) linear
      write(*,chrfmt) "linear?          ",trim(atmp)
      write(*,chrfmt) "symmetry         ",adjustr(symchar)
      write(*,intfmt) "rotational number",nint(sym)
      write(*,dblfmt) "scaling factor   ",fscal,"    "
      write(*,dblfmt) "rotor cutoff     ",sthr, "cm⁻¹"
      write(*,dblfmt) "imag. cutoff     ",ithr, "cm⁻¹"
      write(*,'(10x,":",49("."),":")')
    endif

    vibs = vibs * rcmtoau   ! thermodyn needs vibs and zp in Eh

    zp = 0.5_wp * sum(vibs(1:nvib))
    adum=abs(temps-298.15d0)
    rt = minloc(adum,1)  !temperature closest to 298.15 is the ref.
    do j=1,nt
      if( (j==rt) .and. pr)then
          pr2=.true.
      else
          pr2=.false.
      endif
      if(pr2)then
      call print_thermo_sthr_ts(iunit,nvib,vibs,avmom,sthr,temps(j))
      endif
      call thermodyn(iunit,a,b,c,avmom,linear,atom,sym,molmass,vibs,nvib, &
      & temps(j),sthr,et(j),ht(j),gt(j),ts(j),zp,pr2)
      stot(j) = (ts(j)/temps(j))*autocal
    enddo

    if((nt>1) .and. pr)then
       write(iunit,'(a)')
       write(iunit,'(a10)',advance='no') "T/K"
       write(iunit,'(a16)',advance='no') "H(0)-H(T)+PV"
       write(iunit,'(a16)',advance='no') "H(T)/Eh"
       write(iunit,'(a16)',advance='no') "T*S/Eh"
       write(iunit,'(a16)',advance='no') "G(T)/Eh"
       write(iunit,'(a)')
       write(iunit,'(3x,72("-"))')
       do i = 1, nt
          write(iunit,'(3f10.2)',advance='no') temps(i)
          write(iunit,'(3e16.6)',advance='no') ht(i)
          write(iunit,'(3e16.6)',advance='no') et(i)
          write(iunit,'(3e16.6)',advance='no') ts(i)
          write(iunit,'(3e16.6)',advance='no') gt(i)
          if (i == rt) then
            write(iunit,'(1x,"(used)")')
          else
            write(iunit,'(a)')
          endif
       enddo
       write(iunit,'(3x,72("-"))')
    endif

    deallocate(vibs)
    return
end subroutine calcthermo

subroutine thermo_wrap(env,pr,nat,at,xyz,dirname, &
        &  nt,temps,et,ht,gt,stot,bhess)
    use iso_fortran_env, wp => real64
    use crest_data
    use iomod
    use strucrd
    implicit none
    type(systemdata) :: env
    logical,intent(in) :: pr
    integer,intent(in) :: nat
    integer,intent(inout) :: at(nat)
    real(wp),intent(inout) :: xyz(3,nat)  !in Angstroem
    character(len=*) :: dirname
    integer,intent(in)  :: nt
    real(wp),intent(in)  :: temps(nt)
    real(wp),intent(out) :: et(nt)          !< enthalpy in Eh
    real(wp),intent(out) :: ht(nt)          !< enthalpy in Eh
    real(wp),intent(out) :: gt(nt)          !< free energy in Eh
    real(wp),intent(out) :: stot(nt)        !< entropy in cal/molK
    logical,intent(in) :: bhess       !calculate bhess instead?

    logical :: subdir,ex
    integer :: i,io,r,ich
    character(len=1024) :: jobcall
    character(len=*),parameter :: pipe ='2>/dev/null'
    character(len=*),parameter :: xname='freq.xyz'
    character(len=:),allocatable :: optpath
    character(len=:),allocatable :: jobcall2
    character(len=128) :: atmp
    character(len=258) :: thispath
    real(wp) :: etot
    integer :: nfreq
    real(wp),allocatable :: freq(:)
    real(wp) :: ithr,fscal,sthr

    integer :: TID,OMP_GET_THREAD_NUM

!!$OMP PARALLEL PRIVATE(TID)
    TID = OMP_GET_THREAD_NUM()
    !awrite(*,*) '---->',TID
!!$OMP END PARALLEL 
    ich = (TID+1)*1000   ! generate CPU dependent file channel number

    call initsignal()

    optpath=''

    subdir = .false.
    i = len_trim(dirname)
    if( i > 0) subdir = .true.

    write(jobcall,'(a,1x,a,1x,a,'' --ohess '',a,1x,a,1x,a,'' > xtb.out'')') &
    &  trim(env%ProgName),trim(xname),trim(env%gfnver),trim(env%solv), &
    &  '--ceasefiles',trim(pipe)

    if(bhess)then
    write(jobcall,'(a,1x,a,1x,a,'' --bhess loose '',a,1x,a,1x,a,'' > xtb.out'')') &
    &  trim(env%ProgName),trim(xname),trim(env%gfnver),trim(env%solv), &
    &  '--ceasefiles',trim(pipe)
    endif


    if(subdir)then
         call rmrf(trim(dirname))
         r = makedir(trim(dirname))
         optpath=trim(dirname)//'/'
    endif     

    !if(env%chrg.ne.0)then
    !   open(unit=ich,file=trim(optpath)//'.CHRG')
    !   write(ich,*)env%chrg
    !   close(ich)
    !endif
    !if(env%uhf.ne.0)then
    !   open(unit=ich,file=trim(optpath)//'.UHF')
    !   write(ich,*)env%uhf
    !   close(ich)
    !endif
    call env%wrtCHRG(trim(optpath))   
    inquire(file='gfnff_topo',exist=ex)
    if(env%gfnver=='--gff' .and. subdir .and.ex)then
        call getcwd(thispath)
        io = sylnk(trim(thispath)//'/'//'gfnff_topo',trim(optpath)//'gfnff_topo')
    endif
    io = sylnk(trim(thispath)//'/'//env%fixfile,trim(optpath)//env%fixfile)

!$omp critical
    open(unit=ich,file=trim(optpath)//xname)
    call wrxyz(ich,nat,at,xyz)
    if(env%thermo%constrhess)then
        call write_cts(ich,env%cts)
    endif
    close(ich)
!$omp end critical

    if(subdir)then
       jobcall2='cd '//trim(dirname)//' && '//trim(jobcall) 
       !call execute_command_line('cd '//trim(dirname)//' && '//trim(jobcall), exitstat=io)
       call execute_command_line(trim(jobcall2), exitstat=io)
    else
       call execute_command_line(trim(jobcall), exitstat=io)
    endif

    et=0.0_wp
    ht=0.0_wp
    gt=0.0_wp
    stot=0.0_wp

    if(io /= 0)then  !if the calc failed
        return
    endif

!$omp critical
    call rdxmol(trim(optpath)//'xtbopt.xyz',nat,at,xyz,atmp)
    etot = grepenergy(atmp)
    nfreq = 3*nat

    allocate(freq(nfreq))
    call rdfreq(trim(optpath)//'vibspectrum',nfreq,freq)

    ithr=env%thermo%ithr
    fscal=env%thermo%fscal
    sthr=env%thermo%sthr
    call calcthermo(nat,at,xyz,freq,pr,ithr,fscal,sthr, &
    &    nt,temps,et,ht,gt,stot )
    deallocate(freq)
!$omp end critical
    call initsignal()
    return
end subroutine thermo_wrap


!--- read vibspectrum file in TM format
subroutine rdfreq(fname,nmodes,freq)
      use iso_fortran_env, wp => real64
      use crest_data
      use iomod
      implicit none
      character(len=*),intent(in) :: fname
      integer,intent(in)   :: nmodes
      real(wp),intent(out) :: freq(nmodes)    !frequencies
      integer :: k,ich,io,n
      character(len=256) :: atmp
      real(wp) :: floats(10)
      logical :: ex
      integer :: TID,OMP_GET_THREAD_NUM

!!$OMP PARALLEL PRIVATE(TID)
      TID = OMP_GET_THREAD_NUM()
!      write(*,*) '---->',TID
!!$OMP END PARALLEL 
      ich = (TID+1)*1000   ! generate CPU dependent file channel number

      freq=0.0_wp
      inquire(file=fname,exist=ex)
      if(.not.ex) return
      k=1 !modes
      open(file=fname,unit=ich)
      rdfile : do
         read(ich,'(a)',iostat=io) atmp
         if(io<0)exit
         if(index(atmp,'$vibrational spectrum').ne.0)then
           rdblock : do
             read(ich,'(a)',iostat=io) atmp
             if(io<0)exit rdfile
             if(index(atmp,'$end').ne.0)exit rdfile
             if(index(atmp,'#').ne.0)cycle rdblock !skip comment lines
             call readl(atmp,floats,n)
             freq(k)=floats(2)
             k=k+1
           enddo rdblock
         endif
      enddo rdfile
      close(ich)
      return
end subroutine rdfreq



!=====================================================================!
! Calculate S_RRHO averages for a given ensemlbe
!=====================================================================!
subroutine calcSrrhoav(env,ensname)
    use iso_fortran_env, only: wp => real64, ou=>output_unit
    use crest_data
    use strucrd
    use iomod
    implicit none
    type(systemdata) :: env
    character(len=*) :: ensname

    real(wp),allocatable :: cp(:)
    real(wp),allocatable :: hconf(:)
    integer :: nat,nall
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:,:)
    real(wp),allocatable :: er(:)
    real(wp),allocatable :: erel(:)
    real(wp),allocatable :: efree(:,:)
    real(wp),allocatable :: g(:)    !degeneracies, either read from cre_degen2 (if present), or set to unity
    real(wp),allocatable :: p(:,:)  ! populations at different T
    real(wp),allocatable :: gatt(:,:)
    real(wp),allocatable :: satt(:,:)
    real(wp),allocatable :: srrho(:),sav(:)
    real(wp),allocatable :: bsatt(:)
    real(wp),allocatable :: gav(:)
    real(wp),allocatable :: pdum(:)
    real(wp) :: psum,emin,sdum
    real(wp) :: quick_rmsd,rmsdval
    integer :: eloc
    logical :: avbhess

    real(wp),parameter :: Tref = 298.15  !room temperature is reference
    integer :: nt
    real(wp),allocatable :: temps(:)
    real(wp),allocatable :: et(:)
    real(wp),allocatable :: ht(:)
    real(wp),allocatable :: gt(:)
    real(wp),allocatable :: stot(:)
    real(wp),allocatable :: c0(:,:)
    character(len=64) :: atmp
    integer :: i,j,k,ich,io,popf
    logical :: ex
    logical :: niceprint
    real(wp) :: percent
    character(len=52) :: bar
    integer :: ncalc,vz,nlimit,nav
    character(len=512) :: thispath,tmppath
    
    real(wp),parameter :: kcal = 627.50947428_wp

     

!--- read the given ensemble
    call rdensembleparam(trim(ensname),nat,nall)
    allocate(at(nat),xyz(3,nat,nall),er(nall))
    call rdensemble(trim(ensname),nat,nall,at,xyz,er)

    if(any(er(:) > 0.0d0))then
        error stop 'ensemble file must contain energies in Eh! must stop'
    endif

    write(tmppath,'(a)') 'Frequency Calculation and Averages'
    write(*,*)
    call smallhead(trim(tmppath))

!--- temperatures from sys object
    if(.not.allocated(env%thermo%temps))then
        call env%thermo%get_temps()
    endif
    nt=env%thermo%ntemps
    allocate(temps(nt))
    temps = env%thermo%temps

!--- space for populations and degeneracies
    allocate(g(nall), source=1.0_wp)
    allocate(p(nall,nt))
!--- read degeneracies?    
    inquire(file='cre_degen2',exist=ex)
    if(ex)then
      open(newunit=ich,file='cre_degen2')
      read(ich,*) atmp 
      do i=1,nall
         read(ich,*,iostat=io)j,g(i)
         if(io < 0) exit
      enddo
      close(ich)
    else
      g=1
    endif

!--- determine how many hessians must be calculated
    allocate(erel(nall),pdum(nall))
    erel = (er - er(1))*kcal
    call entropy_boltz(nall,Tref,erel,g,pdum)

    ncalc=1  !always the lowest
    psum=pdum(1)
    nlimit = env%thermo%pcap !limit strucs (for VERY large SE)
    do i=2,nall
       psum = psum + pdum(i)
       ncalc = ncalc+1
       if(ncalc == nlimit)then
           exit
       endif
       if(psum > env%thermo%ptot)then
           exit
       endif
    enddo
    deallocate(pdum)
    

!--- print something
    write(*,'(1x,a,i0)') 'Nconf on file      : ',nall
    write(atmp,'(1x,a,f5.2,a)') '(= ',psum*100.0d0,'% total population)'
    write(*,'(1x,a,i0,a)') 'Taken for Hessians : ',ncalc,trim(atmp)
    if(psum < env%thermo%ptot)then
    write(*,'(2x,a,i0,a)') '=> (Limited to ',ncalc,' structures due to amount of calcs.)'
    endif
    write(*,'(1x,a,f8.2,1x,f8.2)') "T range  /K    : ",temps(1),temps(nt)
    write(*,'(1x,a,f17.6,1x,a)')   "scaling factor : ",env%thermo%fscal,"    "
    write(*,'(1x,a,f17.6,1x,a)')   "rotor cutoff   : ",env%thermo%sthr, "cm⁻¹"
    write(*,'(1x,a,f17.6,1x,a)')   "imag. cutoff   : ",env%thermo%ithr, "cm⁻¹"
    write(*,*)

!--- calculate hessians for ncalc lowest structures   
    allocate(gatt(nall,nt),satt(nall,nt), source=0.0_wp)

    io = makedir('HESSIANS')
    call getcwd(thispath)
    inquire(file='gfnff_topo',exist=ex)

    call chdir('HESSIANS')
    if(env%gfnver=='--gff' .and.ex)then
    call getcwd(tmppath)
    io = sylnk(trim(thispath)//'/'//'gfnff_topo',trim(tmppath)//'/'//'gfnff_topo')
    endif
    io = sylnk(trim(thispath)//'/'//env%fixfile,trim(tmppath)//'/'//env%fixfile)

    k=0
    niceprint=env%niceprint

    !--- OMP stuff 
    if(env%autothreads)then
!    call ompautoset(env%threads,7,env%omp,env%MAXRUN,1) !serial with multiple cpus
    call ompautoset(env%threads,7,env%omp,env%MAXRUN,ncalc) !set global OMP/MKL variable for xtb jobs
    endif

    avbhess = env%thermo%avbhess
    write(*,'(1x,a,i0,a)') 'Running ',ncalc,' calculations ...'
!$omp parallel &
!$omp shared( vz,tmppath,ncalc,percent,k,bar,niceprint) &
!$omp shared( nat,at,xyz,c0,et,ht,gt,stot,temps,nt,gatt,satt,avbhess ) 
!$omp single
     allocate(et(nt),ht(nt),gt(nt),stot(nt))
     allocate(c0(3,nat))
     do i=1,ncalc
       call initsignal()
       vz=i
       !$omp task firstprivate( vz ) private( tmppath,et,ht,gt,stot,c0 )
       call initsignal()
       !$omp critical
       write(tmppath,'(''hess'',i0)')vz
       c0(1:3,1:nat) = xyz(1:3,1:nat,vz) 

       !$omp end critical
       call thermo_wrap(env,.false.,nat,at,c0,tmppath, &
       &    nt,temps,et,ht,gt,stot,avbhess)
       !$omp critical
       gatt(vz,1:nt) = gt(1:nt)
       satt(vz,1:nt) = stot(1:nt)  

       !$omp end critical 
       if(.not.env%keepModef) call rmrf(trim(tmppath)) 

       !$omp critical
        k=k+1
        if(niceprint)then
          percent=float(k)/float(ncalc)*100
          call  progbar(percent,bar)
          call printprogbar(percent,bar)
        else
          if(gui)then
             call wrGUIpercent(k,ncalc,100)
          else
            write(6,'(1x,i0)',advance='no')k
            flush(6)
          endif
        endif
      !$omp end critical
      !$omp end task
    enddo
    deallocate(c0)
    deallocate(stot,gt,ht,et)
!$omp taskwait
!$omp end single
!$omp end parallel
    if(niceprint)then
        write(*,'(/)')
    else
        write(*,'(1x,a,/)') 'done.'
    endif
    call chdir(thispath)
    if(.not.env%keepModef) call rmrf('HESSIANS')

!--- process the calculated free energies and entropies into accurate populations
    allocate(srrho(nt),sav(nt),gav(nt),efree(nall,nt))
    srrho=0.0_wp
    sav=0.0_wp
    gav=0.0_wp
    nav=ncalc
    write(ou,'(1x,a)',advance='no') 'calculating averages for G and S ... '
    flush(ou)
    do j=1,nt
      do i=1,ncalc
        if(abs(gatt(i,j)).lt. 1.d-10)then  !failed calcs?
            if(j==1) nav = nav -1
        endif
        gav(j) = gav(j) + gatt(i,j)
        sav(j) = sav(j) + satt(i,j)
      enddo
    enddo
    gav = gav/float(nav)   !get the average G(T)
    sav = sav/float(nav)   !get the avverage S(T)
    write(*,'(a8)')'done.'

!--- get the free energies    
    do j=1,nt
      efree(:,j) = er(:) !all based on etot
      do i=1,nall
        if(i <= ncalc)then
            if(abs(gatt(i,j))<1.d-10)then
            efree(i,j) = efree(i,j) + gav(j) !add |G(T)| (for failed calcs)
            else
            efree(i,j) = efree(i,j) + gatt(i,j) !add G(T)
            endif
        else
!-- for all energies that were not included in the free energy calculation add the average
            efree(i,j) = efree(i,j) + gav(j)
        endif
      enddo
    enddo
    
    !-- make relative energies and calculate Boltzman populations
    write(ou,'(1x,a)',advance='no') 'calculating Boltzmann weights ... '
    flush(ou)
    allocate(pdum(nall))
    do j=1,nt
       emin = minval(efree(:,j),1) !lowest as reference
       erel = (efree(:,j)-emin)*kcal     !to relative energies in kcal/mol
       pdum=0.0d0
       call entropy_boltz(nall,temps(j),erel,g,pdum)
       !call entropy_boltz(ncalc,temps(j),erel,g(1:ncalc),pdum(1:ncalc))
       p(:,j) = pdum(:)
    enddo
    deallocate(pdum)
    write(*,'(a11)')'done.'
    if(env%thermo%printpop)then
    popf = makedir('populations')    
    do j=1,nt
    write(tmppath,'(a,a,a,i0)')'populations','/','.pop_',nint(temps(j))
    open(newunit=popf,file=trim(tmppath))    
      do k=1,nall
      write(popf,'(f16.8)') p(k,j)
      enddo
    close(popf)
    enddo
    endif
!>========================================================================================<!
!>==== after this point p now contains the correct populations based on free energies ====<!
!>========================================================================================<!
!--- S_avRRHO must be calculated relative to the actual DFT reference structure
!--- the corresponding frequencies can be calculated with bhess
     if(env%emtd%bhess)then
         allocate(bsatt(nt))
         allocate(et(nt),ht(nt),gt(nt))
         write(ou,'(1x,a)',advance='no') 'calculating reference S (bhess) ... '
         flush(ou)
         call thermo_wrap(env,.false.,env%emtd%nat,env%emtd%at,  &
       &    env%emtd%xyz,'BHESS', nt,temps,et,ht,gt,bsatt,.true.)
         if(.not.env%keepModef) call rmrf('BHESS')
         deallocate(gt,ht,et)
         write(*,'(a9)')'done.'
     endif
!--- average S_rrho with CORRECT populations
     write(ou,'(1x,a)',advance='no') 'calculating δSrrho ... '
     flush(ou)
     srrho=0.0d0
     do j=1,nt
        !write(*,*) sum(p(1:ncalc,j))
        do i=1,nall
        !do i=1,ncalc
           if(i <= ncalc)then
            if(abs(satt(i,j))<1.d-10)then
            srrho(j) = srrho(j) + p(i,j)*sav(j)! -satt(1,j))    ! (for failed hess calcs)
            else    
            srrho(j) = srrho(j) + p(i,j)*satt(i,j)! -satt(1,j)) ! corrected for different S_rrho
            endif
           else
           srrho(j) = srrho(j) + p(i,j)*sav(j)! -satt(1,j))    
           endif
        enddo
        ! substract the reference value to shift the average
        if(env%emtd%bhess)then
            srrho(j) = srrho(j) - bsatt(j)
        else
            srrho(j) = srrho(j) - satt(1,j)
        endif
     enddo
     write(*,'(a22)')'done.'


     if(env%emtd%bhess)then
      write(*,*)
      call underline('Coordinates for the bhess reference structure (Ångström):')
      call wrxyz(ou,env%emtd%nat,env%emtd%at,env%emtd%xyz)
      write(*,*)'-----------------------------------------------------------'
      write(*,'(1x,a,a,a)') 'as read from <',env%emtd%fromfile,'>'
      inquire(file='crest_best.xyz',exist=ex)
      if(ex)then
      rmsdval = quick_rmsd('crest_best.xyz',env%emtd%nat,env%emtd%at,env%emtd%xyz,.true.)
      write(*,'(1x,a)') 'Heavy-atom RMSD between lowest conformer and this reference :'
      write(*,'(1x,a,f16.6,a)') 'RMSD(heavy) =',rmsdval,' Å'
      write(*,*)
      endif
      write(*,'(1x,a)') 'msRRHO(bhess) reference entropies:'
      do i = 1, nt
          write(*,'(2x,f10.2,2x,f16.6)')temps(i),bsatt(i)
       enddo
     endif

!--- prinout for the average free energy and entropy    
     if((nt>1))then
       write(*,'(a)')
       write(*,'(a10)',advance='no') "T/K"
       write(*,'(a16)',advance='no') "δS/cal/molK"
       write(*,'(a16)',advance='no') "|G(T)|/Eh"
       write(*,'(a16)',advance='no') "G_lowest/Eh"
       write(*,'(a8)',advance='no') "(conf)"
       write(*,'(a)')
       write(*,'(3x,63("-"))')
       do i = 1, nt
          write(*,'(3f10.2)',advance='no') temps(i)
          write(*,'(3e16.6)',advance='no') srrho(i)
          write(*,'(3e16.6)',advance='no') gav(i)
          emin = minval(efree(:,i),1)
          write(*,'(f16.6)',advance='no') emin
          eloc = minloc(efree(:,i),1)
          write(*,'(i8)',advance='no') eloc
          write(*,'(a)')
       enddo
       write(*,'(3x,63("-"))')
       write(*,'(3x,a,a)')'NOTE: if |G(T)| is the averaged ', &
       & 'contributrion to the free energy.'
       write(*,'(3x,a,a,i0,a)')'|G(T)| used only for the higher-energetic ', &
       & 'structures (n > ',ncalc,').'
       write(*,'(3x,a,a)')'All other structures use ', &
       & 'G(T) from the respective Hessian calculations.'
    endif

!--- properties based on free energies
    allocate(cp(nt),hconf(nt))
    do j=1,nt
     emin = minval(efree(:,j),1) !lowest as reference
     erel = (efree(:,j)-emin)*kcal     !to relative energies in kcal/mol
     call entropy_S(nall,temps(j),1.0d0,erel, &
     &    g,sdum,cp(j),hconf(j))  ! g read from file or set to 1
    enddo 
    
    if((nt>1))then
       write(*,*)
       write(*,'(1x,a)') 'Quantities calculated on free energies:'
       write(*,'(a10)',advance='no') "T/K"
       write(*,'(a16)',advance='no') "δSrrho"
       write(*,'(a16)',advance='no') "Cp(T)"
       write(*,'(a16)',advance='no') "[H(T)-H(0)]"
       write(*,'(a)')
       write(*,'(3x,55("-"))')
       do i = 1, nt
          write(*,'(3f10.2)',advance='no') temps(i)
          write(*,'(3e16.6)',advance='no') srrho(i)
          write(*,'(3e16.6)',advance='no') cp(i)
          write(*,'(f16.6)',advance='no') hconf(i)
          write(*,'(a)')
       enddo
       write(*,'(3x,55("-"))')
    endif

    if(allocated(env%emtd%soft))then
       env%emtd%soft(:) = srrho(:)   !this is \overline{S}_{msRRHO}
    endif
    if(allocated(env%emtd%cpoft))then
       env%emtd%cpoft(:) = cp(:)    !this is Cp_conf
    endif
    if(allocated(env%emtd%hoft))then
       env%emtd%hoft(:) = hconf(:)  !this is H_conf
    endif

    if(allocated(bsatt))deallocate(bsatt)
    deallocate(satt,gatt)
    deallocate(efree,gav,sav,srrho)
    deallocate(erel,p,g,temps)
    deallocate(er,xyz,at)
    return
end subroutine calcSrrhoav
