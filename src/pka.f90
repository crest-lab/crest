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

!=========================================================================!
! AUTOMATED pKa CALCULATION SCRIPT
! To use run:
!    crest <input> -pka <acidic H atom>
!
!  <input> should be the Acid and <H atom> specifies 
!  which H atom is removed
!=========================================================================!
subroutine pkaquick(env,tim)
      use iso_fortran_env, wp => real64
      use crest_data
      use strucrd
      implicit none

      type(systemdata) :: env
      type(timer)   :: tim
      type(coord) :: ACID
      type(coord) :: BASE
      type(ensemble) :: ACIDENSEMBLE
      type(ensemble) :: BASEENSEMBLE
      real(wp) :: t
      real(wp) :: GA,GB,dG
      real(wp) :: eatb,gsa,grrhoa,ebtb,gsb,grrhob,dE
      real(wp) :: pKa,dum,pkaref
      interface
          function pKaPolyFER(dG,nc,c,T) result(pka)
              import :: wp
              implicit none
              real(wp) :: dG         !in Eh
              real(wp),optional :: T !in K
              integer :: nc
              real(wp) :: c(nc)
              real(wp) :: pka
          end function pKaPolyFER
      end interface
      integer :: nc
      real(wp),allocatable :: c(:)
      integer :: i,j,h,ich
      integer :: refnat,refchrg,basechrg
      logical :: bhess,ex,rangepka
      real(wp) :: pkamin,pkamax
      integer :: nalla,nallb
      real(wp),allocatable :: gacid(:),pacid(:)
      real(wp),allocatable :: gbase(:),pbase(:)
      real(wp),parameter :: kcal =627.5095d0

      character(len=:),allocatable :: pkaparam
      character(len=128) :: parinfo

      !T=env%tboltz
      T=298.15d0  !FER was fitted with this temperature
      refnat = env%nat
      refchrg = env%chrg
      basechrg = env%chrg-1
      bhess = .true.
      rangepka=.false.
      nc=5
      allocate(c(5), source=0.0_wp)

      !---- FER parameter selection
      parinfo=''
      if(env%ptb%rdCFER)then
         pkaparam = 'external'
      else    
         !pkaparam ='none'
         pkaparam = env%gfnver
      endif  
      select case( pkaparam )
        case( '--gfn1' )
          error stop 'GFN1-xTB not available for pKa calculation via CFER.'   
        case( 'quartic' )
          c(1)= -14501.33156900
          c(2)=    276.23925454
          c(3)=     -1.97288872
          c(4)=      0.00624684 
          c(5)=     -0.00000737
        case('--gfn2' )  
          parinfo='and GFN2 parameters (c0-c3) fitted on PKA74'  
          !c(1)= -1656.74953643 !-1645.21695955
          !c(2)=    23.18527638 ! 23.06345883
          !c(3)=    -0.11102517 !-0.11064255
          !c(4)=     0.00018350 ! 0.00018327
          c(1)=  -1855.025277d0
          c(2)=    26.075982d0
          c(3)=    -0.12496355d0
          c(4)=    0.00020571d0
        case( 'oldparam' )
          c(1)=-1656.7    ! SD = 2.87, NO isorad, 19normal, my COSMOTHERM (out.ccf), adjust for Fabian's out.cosmo, ie -1656.45
          c(2)= 23.185    !    = 2.47 without two extreme outliers (CH3NN+, cycnoform)
          c(3)=-0.11103   !           and checked outlier HClO4
          c(4)= 0.0001835
        case( 'external' )
          call pka_rdparam(env%ptb%cferfile,nc,c,parinfo)
        case default  !do nothing with the values  
          c = 0.0d0
          c(2) = 1.0d0
     end select  


   select case(env%ptb%pka_mode)  
   case( 0 )  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!  
      env%crestver=2
      env%gbsa=.true.
      if((index(env%solv,'--alpb h2o')==0)) then
          error stop 'Error: pKa tool only available for GFN2/ALPB(H2O) calculations'
      endif
      if((env%gfnver.ne.'--gfn2')) then
          error stop 'Error: pKa tool only available for GFN2/ALPB(H2O) calculations'
      endif
      if(env%ptb%h_acidic==0)then
          error stop 'Error: no hydrogen atom selected for acid dissociation'
      endif

!---- very important preperation: resort input so that heavy atoms come first, and then all hydrogens
      allocate(env%ptb%atmap(env%nat))
      call htothebottom('coord',env%chrg,env%nat,env%ptb%atmap)
      !--- check for the selected acidic H atom
      h=0
      do i=1,env%nat
       if(env%ptb%atmap(i)==env%ptb%h_acidic) h=i
      enddo
      if((h.ne.env%ptb%h_acidic).and.(env%ptb%h_acidic>0))then
          write(*,'(1x,a,i0,a,i0)')'The position of the selected acidic hydrogen was updated ', &
          &    env%ptb%h_acidic,' ---> ',h
          env%ptb%h_acidic = h
      endif

!---- first do conformational search in water
      call largehead2('Calculation of the acid')
      !--- because the acid must be the input file, it is saved as "coord"
      call ACID%open('coord')
      h=env%ptb%h_acidic
      if((h>0))then
        if( ACID%at(h).ne.1) error stop 'selected atom in pKa tool is not a hydrogen atom'
      endif
      call ACID%write('acid.xyz')      !--- resave it as 'acid.xyz'

      if(env%preopt)then  !--- optional geometry optimization
          write(6,'(1x,a)',advance='no')'Optimizing acid geometry ... '
          flush(6)
          call miniopt(env,'acid.xyz')
          call xyz2coord('acid.xyz','coord')  !--update coord for possible confsearch
          write(*,*) 'done.'
      endif
      if(env%relax .and. env%preopt)then      !--- optional conformational sampling
      call confscript2i(env,tim) !MTD-GC algo
      call rename('crest_conformers.xyz','crest_pka_acid.xyz')
      call rename('crest_best.xyz','acid.xyz')
      call V2cleanup(.false.)
      endif
      !---finally, update geometry in memory
      call ACID%deallocate()
      call ACID%open('acid.xyz')

!---- then do the (relaxed) deprotonation
      call largehead2('Calculation of the base')
      env%ensemblename='none selected'  !RESET, IMPORTANT!
      env%ptb%threshsort=.true.
    !=================================================================================!      
      !--- write the base input file
      if(env%ptb%h_acidic>0)then            !-- for a manually selected H atom  
        BASE%nat=refnat-1
        allocate(BASE%at(BASE%nat), source=0)
        allocate(BASE%xyz(3,BASE%nat), source=0.0_wp)
        j=0
        do i=1,ACID%nat
           if(i==env%ptb%h_acidic)cycle
           j=j+1
           BASE%at(j) = ACID%at(i)
           BASE%xyz(1:3,j)=ACID%xyz(1:3,i)  
        enddo
      else if(env%ptb%h_acidic==-1)then      !-- "automatic" mode  
        call deprotonate(env,tim)
        call rmrfw('deprotonate_')
        call BASE%open('deprotonated.xyz')
      else if(env%ptb%h_acidic==-2)then      !-- read-in base file mode  
        call BASE%open(env%ptb%pka_baseinp)
      endif
      call BASE%write('base.xyz')
      env%chrg=refchrg-1
      basechrg=refchrg-1
      env%nat = refnat - 1
      env%rednat = refnat - 1
   !=================================================================================!      
      if(env%preopt)then  !--- optional geometry optimization
          write(6,'(1x,a)',advance='no')'Optimizing base geometry ... '
          flush(6)
          call miniopt(env,'base.xyz')
          write(*,*) 'done.'
      endif
      if(env%relax .and. env%preopt)then
      call relaxensemble('base.xyz',env,tim)
      call rename('relax.base.xyz','crest_pka_base.xyz')
      endif
      !-- update geometry in memory
      call BASE%deallocate()
      if(env%relax .and. env%preopt)then
      call BASE%open('crest_pka_base.xyz')
      call BASE%write('base.xyz')
      else
      call BASE%open('base.xyz')
      endif

!--- get the free energies and correction term
      env%chrg=refchrg
      env%nat = refnat
      env%rednat = refnat
      write(*,*)
      call acidbase(env,'acid.xyz','base.xyz',refchrg,.true.,.true.,dE, &
        &  bhess,eatb,gsa,grrhoa,ebtb,gsb,grrhob)
      GA = eatb   !already includes gsa 
      GB = ebtb   !already includes gsb
      if(bhess)then
          GA = GA + grrhoa
          GB = GB + grrhob
      endif
   case( 1 )  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!  
       call largehead2('Read ensemble data for pKa calculation')
       write(*,*)'Note: In order to work properly the ensembles must be in the xyz'
       write(*,*)'format and contain the final free energies in solution (in Eh) as'
       write(*,*)'comment line for each structure.'
       write(*,*)
       call ACIDENSEMBLE%open(env%ptb%pka_acidensemble)
       write(*,'(1x,a,a,a,i0,a)')'Ensemble read for acid: ',env%ptb%pka_acidensemble, &
           & ' with ',ACIDENSEMBLE%nall,' structures'
       call BASEENSEMBLE%open(env%ptb%pka_baseensemble)
       write(*,'(1x,a,a,a,i0,a)')'Ensemble read for base: ',env%ptb%pka_baseensemble, &
           & ' with ',BASEENSEMBLE%nall,' structures'
       write(*,'(1x,a)') 'Calculating population average from read-in (free) energies ...'
       nalla = ACIDENSEMBLE%nall
       allocate(gacid(nalla), source=1.0_wp)
       allocate(pacid(nalla), source=0.0_wp)
       call pka_boltz(nalla,T,ACIDENSEMBLE%er,gacid,pacid)
       GA = 0.0_wp
       do i=1,nalla
         GA = GA + pacid(i)*ACIDENSEMBLE%er(i)
       enddo
       nallb = BASEENSEMBLE%nall
       allocate(gbase(nallb), source=1.0_wp)
       allocate(pbase(nallb), source=0.0_wp)
       call pka_boltz(nallb,T,BASEENSEMBLE%er,gbase,pbase)       
       GB = 0.0_wp
       do i=1,nallb
         GB = GB + pbase(i)*BASEENSEMBLE%er(i)
       enddo
       deallocate(pbase,gbase,pacid,gacid)
       write(*,'(1x,a,f16.8,a)') '|G(AH)| =',GA,' Eh'
       write(*,'(1x,a,f16.8,a)') '|G(A⁻)| =',GB,' Eh'
       dE = 0.0_wp !would be the energy correction Ecorr, ignore for read-in

       !--- for statistics: determine the min and max pKa in ensemble (conformational influence)
       if(nalla>1 .or. nallb>1)then
           rangepka=.true.
           !-- initialize
           dum = GB - GA
           pkamin =  pKaPolyFER(dum,nc,c,T)
           pkamax =  pkamin
           !-- loop
           do i=1,nalla
            do j=1,nallb
              dG = BASEENSEMBLE%er(j) - ACIDENSEMBLE%er(i)
              dum = pKaPolyFER(dG,nc,c,T) 
              if(dum > pkamax ) pkamax = dum
              if(dum < pkamin ) pkamin = dum
            enddo
           enddo
       endif
   end select !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!  
      
!--- calculate the pKa via the free energy relationship
      call largehead2('pKa CALCULATION')
      write(*,*)
      write(*,'(1x,a)') 'For the reaction AH + H₂O --->  A⁻ + H₃O⁺'
      write(*,'(1x,a)') '(Note: H₂O/H₃O⁺ is not included in ΔG)'
      write(*,'(1x,a,f16.2,a)') 'T     =',T,' K'
      write(*,'(1x,a,f16.8,a)') 'G(AH) =',GA,' Eh'

      write(*,'(1x,a,f16.8,a)') 'G(A⁻) =',GB,' Eh' 
      dum = GB - GA
      write(*,'(1x,a,f16.8,a,f8.2,a)') 'ΔG       =',dum,' Eh,',dum*kcal,' kcal/mol'
      dG = dum + dE
      if(env%ptb%pka_mode.ne.1)then
      write(*,'(1x,a,f16.8,a,f8.2,a)') 'ΔG+Ecorr =',dG,' Eh,',dG*kcal,' kcal/mol'
      endif
      write(*,*)

      write(*,'(1x,a)') 'polynomial free energy relationship (FER):'
      write(*,'(2x,a)') 'pKa   = c0 + c1*kdiss + c2*kdiss² + ... + c_n*kdiss^n'
      write(*,'(2x,a,a)') 'with kdiss = ΔG(aq)/ln(10)RT ',trim(parinfo)
      write(*,'(7x,a,f8.2,a)') 'ΔG(aq)=',dG*kcal,' kcal/mol'
      do i=1,nc
        if(c(i)==0.0_wp) cycle
        j=i-1
        write(*,'(7x,a,i0,a,f16.8)') 'c',j,' =',c(i)
      enddo

      pka = pKaPolyFER(dG,nc,c,T)

      write(*,*)
      write(*,'(3x,a)') ' ___________________________________ '
      write(*,'(3x,a)') '|                                   |'
      write(*,'(3x,a,f16.2,a)') '| calculated pKa = ',pka,' |'
      inquire(file='.pkaref',exist=ex)
      if(ex)then
          pkaref=0.0_wp
          open(newunit=ich,file='.pkaref')
          read(ich,*) pkaref
          close(ich)
      write(*,'(3x,a,f16.2,a)') '| exptl. pKa     = ',pkaref,' |'
      dum=pka-pkaref
      write(*,'(3x,a)') '|                       ___________ |'
      write(*,'(3x,a,f16.2,a)') '| ΔpKa           = ',dum,' |'
      endif
      if(rangepka)then
      write(*,'(3x,a)') '|                                   |'
      write(*,'(3x,a,f16.2,a)') '| min.pKa (conf) = ',pkamin,' |'
      write(*,'(3x,a,f16.2,a)') '| max.pKa (conf) = ',pkamax,' |'
      endif
      write(*,'(3x,a)') '|___________________________________|'

!--- cleanup
      call ACID%deallocate()
      call BASE%deallocate()
      deallocate(c)

!--- turn off other property modes via env%npq
      env%npq = 0
      return
contains
    subroutine miniopt(env,fname)
       use iso_fortran_env, wp => real64
       use crest_data
       use iomod
       implicit none
       type(systemdata) :: env
       character(len=*) :: fname
       character(len=1028) :: jobcall
       integer :: io
       logical :: ex
       call env%wrtCHRG('')
       write(jobcall,'(a,1x,a,1x,a,'' --ceasefiles --opt vitght '',a,1x,a,'' >xtb.out'')') &
       &    trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),' 2>/dev/null'
       call execute_command_line(trim(jobcall), exitstat=io)
       inquire(file='xtbopt.xyz',exist=ex)
       if(ex)then
           call rename('xtbopt.xyz',fname)
       endif
       call remove('xtbrestart')
       call remove('xtbtopo.mol')
       call remove('.xtboptok')
       call remove('.CHRG')
       call remove('.UHF')
       return
    end subroutine miniopt
end subroutine pkaquick

subroutine pka_argparse(str,h)
    use iso_fortran_env, only: wp => real64
    implicit none
    character(len=*) :: str
    integer :: h,io
    character(len=:),allocatable :: arg
    real(wp) :: rdum
    logical :: ex
    h=0
    arg=trim(str)
    if(arg(1:1)=='-')then
        h=0
        return
    endif
    if(arg=='auto')then
        h=-1
        return
    endif
    inquire(file=arg,exist=ex)
    if(ex)then
        h=-2
        return
    endif
    read(arg,*,iostat=io) rdum
    if((io==0).and.(rdum>0.0d0)) h = nint(rdum)
    return
end subroutine pka_argparse

subroutine pka_argparse2(env,str1,str2,h)
    use iso_fortran_env, only: wp => real64
    use crest_data
    implicit none
    type(systemdata) :: env
    character(len=*) :: str1,str2
    integer :: h
    logical :: ex,ex2
    h=0
    inquire(file=trim(str1),exist=ex)
    inquire(file=trim(str2),exist=ex2)
    if(ex.and.ex2)then
        h=1
        env%ptb%pka_acidensemble = trim(str1)
        env%ptb%pka_baseensemble = trim(str2)
    endif
    return
end subroutine pka_argparse2

subroutine pka_rdparam(fname,nc,c,parinfo) !c1,c2,c3,c4)
    use iso_fortran_env, only: wp => real64
    use filemod
    implicit none
    character(len=*) :: fname
    integer :: nc
    real(wp) :: c(nc)
    !real(wp),intent(out) :: c1,c2,c3,c4
    character(len=*) :: parinfo
    character(len=:),allocatable :: line,param
    type(filetype) :: f
    real(wp) :: dum
    logical :: ex
    integer :: io,i

    c=0.0_wp

    inquire(file=fname,exist=ex)
    if(ex)then  !>--- read from file
        call f%open(fname)
        line=f%line(1)
        call f%close()

        do i=1,nc
        param=getlarg(line,i)
        read(param,*,iostat=io) dum
        if(io==0)then
            c(i)=dum
        else
            c(i)=0.0_wp
        endif
        enddo
    else  !>--- use method dependent parameters from the paper
        select case( fname )
        case( 'gfn2-tr224' )
          parinfo='and GFN2 parameters (c0-c3) fitted on TR224'  
          c(1)=  6702.3111485d0
          c(2)=  -100.4483504d0
          c(3)=     0.4988740d0
          c(4)=    -0.0008201d0
        case( 'r2scan3c', 'r2scan-3c','r2scan3c-pka74' )
          parinfo='and r2SCAN-3c parameters (c0-c3) fitted on PKA74'  
          c(1)=  -1511.8899792d0
          c(2)=    21.1100681d0
          c(3)=    -0.1011999d0
          c(4)=     0.0001683d0
        case( 'r2scan3c-tr224' )
          parinfo='and r2SCAN-3c parameters (c0-c3) fitted on TR224'            
          c(1)=  5014.2837220d0
          c(2)=   -75.7090708d0
          c(3)=     0.3781462d0
          c(4)=    -0.0006239d0
        case( 'b973c', 'b97-3c', 'b973c-pka74' )
          parinfo='and B97-3c parameters (c0-c3) fitted on PKA74'
          c(1)=  -1835.3033945d0
          c(2)=     25.4169227d0
          c(3)=     -0.1201689d0
          c(4)=      0.0001956d0
        case( 'b973c-tr224' )
          parinfo='and B97-3c parameters (c0-c3) fitted on TR224'  
          c(1)=  3032.1086142d0
          c(2)=   -45.1533848d0
          c(3)=     0.2212100d0
          c(4)=    -0.0003555d0
        case default  
          c = 0.0d0
          c(2) = 1.0d0  
          return
        end select  
     endif

    return
end subroutine pka_rdparam

!========================================================!
! compute Boltzmann populations p for given
! T, level energies e(), and 
! level degeneracies g()  
!========================================================!
subroutine pka_boltz(n,t,e,g,p)
      use iso_fortran_env, wp => real64
      implicit none
      integer n
      real(wp) :: e(n),p(n),g(n)          
      real(wp) :: t,f,esum
      integer :: i
      real(wp),allocatable :: erel(:)
      real(wp) :: emin
      real(wp),parameter :: R = 8.31446261815324_wp/4.184_wp
      real(wp),parameter :: kcal =  627.5095_wp
      !-- first convert into relative kcal/mol energies 
      allocate(erel(n))
      emin=minval(e,1)
      do i=1,n
      erel(i) = e(i) - emin
      erel(i) = erel(i) * kcal
      enddo
      f = 1.0d0 / (T * R * 0.001d0)
      esum=0
      do i=1,n
         esum=esum+g(i)*exp(-erel(i)*f)
      enddo
      do i=1,n
         p(i)=g(i)*exp(-erel(i)*f)/esum
      enddo
      deallocate(erel)
      if(abs(1.0d0-sum(p)).gt.1.d-6) stop 'error in pka_boltz()'
      return
end subroutine pka_boltz


!============================================================================!
! Calculate the pKa value via the linear free Energy relationship (LFER)
!
!
! For the reaction   AH + H₂O --->  A⁻ + H₃O⁺
!
! GA = free energy of species AH 
! GB = free energy of species A⁻
!
!============================================================================!
function pKaLFER(GA,GB,c0,c1,T) result(pka)
    use iso_fortran_env, wp => real64
    implicit none
    real(wp) :: GA,GB
    real(wp) :: dG
    real(wp) :: c0,c1
    real(wp) :: T
    real(wp) :: pka
    real(wp),parameter :: kcal =627.5095d0
    real(wp),parameter :: R = 0.00198720425  !in kcal/K mol
    real(wp) :: factor
    pka = 0.0d0
    factor = log(10.0d0)*R*T
    dG = (GB - GA)*kcal
    pka = c0 * (dG / factor) + c1
    return
end function pKaLFER
!============================================================================!
! Calculate the pKa value via a cubic free Energy relationship (CFER)
!
!
! For the reaction   AH + H₂O --->  A⁻ + H₃O⁺
!
! GA = free energy of species AH 
! GB = free energy of species A⁻
! dG = GB - GA
! 
! CFER:
! pKa = c1 + c2*kdiss + c3(kdiss)^2 + c4(kdiss)^3
!
!============================================================================!
function pKaCFER(dG,c1,c2,c3,c4,T) result(pka)
    use iso_fortran_env, wp => real64
    implicit none
    real(wp) :: dG         !in Eh
    real(wp),optional :: T !in K
    real(wp) :: pka
    real(wp),parameter :: kcal =627.5095d0
    real(wp),parameter :: R = 0.00198720425  !in kcal/K mol
    real(wp) :: logk,logkfix
    real(wp) :: c1,c2,c3,c4
    if(.not.present(T))then
    T=298.15_wp 
    endif
    pka = 0.0d0
    !logk = kcal*dG/0.592452/2.302585
    logk =kcal*dG/(log(10.0d0)*R*T)
    !open(unit=102030, file='.kdiss')
    !write(102030,'(1x,f16.8)') logk
    !close(102030)
    logkfix=kcal*dG/(log(10.0d0)*R*298.15_wp)
    pka = c1 + c2*logk + c3*(logkfix**2) + c4*(logkfix**3)
    return
end function pKaCFER
!============================================================================!
! Calculate the pKa value via a general polynominal free Energy relationship (FER)
!
! For the reaction   AH + H₂O --->  A⁻ + H₃O⁺
!
! GA = free energy of species AH 
! GB = free energy of species A⁻
! dG = GB - GA
! kdiss = dG/ln(10)RT
! 
! Polynomial FER:
! pKa = c1 + c2*kdiss + c3(kdiss)^2 + ... + c_n(kdiss)^(n-1)
!
!============================================================================!
function pKaPolyFER(dG,nc,c,T) result(pka)
    use iso_fortran_env, wp => real64
    implicit none
    real(wp) :: dG         !in Eh
    real(wp),optional :: T !in K
    real(wp) :: pka
    real(wp),parameter :: kcal =627.5095d0
    real(wp),parameter :: R = 0.00198720425  !in kcal/K mol
    real(wp) :: logk,logkfix
    integer :: nc,i,j
    real(wp) :: c(nc)
    if(.not.present(T))then
    T=298.15_wp 
    endif
    pka = 0.0d0
    logk =kcal*dG/(log(10.0d0)*R*T)
    open(unit=102030, file='.kdiss')
    write(102030,'(1x,f16.8)') logk
    close(102030)
    logkfix=kcal*dG/(log(10.0d0)*R*298.15_wp)
    pka = c(1) + c(2)*logk
    do i=3,nc
       j=i-1
       pka=pka + c(i)*(logkfix**j)
    enddo
    return
end function pKaPolyFER
