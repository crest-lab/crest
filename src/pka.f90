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
!    crest <input> -pka
! IMPORTANT NOTE:
!  this is a old version of the pka calculation using the LFER 
!  and currently deactivated.
!  the newer version is the "pkaquick" routine using the CFER.
!=========================================================================!
subroutine pkascript(env,tim)
      use iso_fortran_env, wp => real64
      use crest_data
      use strucrd
      implicit none

      type(systemdata) :: env
      !type(options)    :: opt
      type(timer)   :: tim

      type(timer)   :: timer1
      type(timer)   :: timer2

      real(wp) :: t
      real(wp) :: GA,GB
      real(wp) :: pKa,dum
      real(wp) :: pKaLFER !this is a function
      real(wp) :: c0,c1
      integer :: refnat,refchrg

      character(len=:),allocatable :: anionfile
      character(len=:),allocatable :: pkaparam

      t=298.15d0
      refnat = env%nat
      refchrg = env%chrg

      env%crestver=2
      env%gbsa=.true.
!---- first do conformational search in water
      call largehead('Calculation of the ensemble in water')
      env%solv='--gbsa h2o'
      call confscript2i(env,tim) !MTD-GC algo
      call propcalc(conformerfile,10,env,tim) !hessian calc

      call rename('crest_property.xyz','crest_pka_neutral_G.xyz')
      call rename('crest_rotamers.xyz','crest_pka_neutral_E.xyz')
      call V2cleanup(.false.)

!---- then do the (relaxed) deprotonation
      call largehead('Calculation of deprotonated structure')
      env%ensemblename='none selected'  !RESET, IMPORTANT!
      env%ptb%threshsort=.true.
      call deprotonate(env,tim)
      env%chrg=refchrg-1
      env%nat = refnat - 1
      env%rednat = refnat - 1
      call relaxensemble('deprotonated.xyz',env,tim)
      !stop

      anionfile='crest_pka_deprot_E.xyz'
      call rename('relax.deprotonated.xyz',anionfile)
      call propcalc(anionfile,10,env,tim) !hessian calc
      call rename('crest_property.xyz','crest_pka_deprot_G.xyz')

      call rmrfw('deprotonate') !--- clean up other files

!--- get the lowest free energies
      call getelow('crest_pka_neutral_G.xyz',GA,.false.)
      call getelow('crest_pka_deprot_G.xyz',GB,.false.)

!--- calculate the pKa via the linear free energy relationship
      call largehead('pKa CALCULATION')
      write(*,*)
      write(*,'(1x,a)') 'For the reaction AH + H₂O --->  A⁻ + H₃O⁺'
      write(*,'(1x,a,f16.8,a)') 'G_low(AH) =',GA,' Eh'
      write(*,'(1x,a,f16.8,a)') 'G_low(A⁻) =',GB,' Eh' 
      dum = GB - GA
      write(*,'(1x,a,f16.8,a)') 'ΔG_low    =',dum,' Eh'
      write(*,'(1x,a,f16.2,a)') 'T         =',T,' K'
      write(*,*)
  
      if(env%extLFER)then

      else    
         !pkaparam ='none'
         pkaparam = env%gfnver
      endif  
      select case( pkaparam )
        case( '--gfn1' )
          c0 = 0.397985d0
          c1 = -29.801646d0
        case( '--gfn2' )
          c0 = 0.399969d0
          c1 = -27.294225d0
        case( 'external' )
          c0 = 1.0d0  !placeholder
          c1 = 1.0d0  !placeholder 
        case default  !do nothing with the values
          c0 = 1.0d0
          c1 = 0.0d0
     end select  

      write(*,'(1x,a)') 'LFER :'
      write(*,'(6x,a)') 'pKa = c0*ΔG(aq)/ln(10)RT + c1'
      write(*,'(6x,a,f16.8)') 'c0 =',c0
      write(*,'(6x,a,f16.8)') 'c1 =',c1

      pka = pKaLFER(GA,GB,c0,c1,T)

      write(*,*)
      write(*,'(3x,a)') ' ___________________________ '
      write(*,'(3x,a)') '|                           |'
      write(*,'(3x,a,f8.2,a)') '| calculated pKa = ',pka,' |'
      write(*,'(3x,a)') '|___________________________|'


!--- turn off other property modes via env%npq
      env%npq = 0
      return
end subroutine pkascript

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
      real(wp) :: t
      real(wp) :: GA,GB,dG
      real(wp) :: eatb,gsa,grrhoa,ebtb,gsb,grrhob,dE
      real(wp) :: pKa,dum,pkaref
      real(wp) :: pKaCFER !this is a function
      real(wp) :: c0,c1,c2,c3,c4
      integer :: i,j,k,l,h,ich,io
      integer :: refnat,refchrg
      logical :: bhess,ex
      real(wp),parameter :: kcal =627.5095d0

      character(len=:),allocatable :: anionfile
      character(len=:),allocatable :: pkaparam

      !T=env%tboltz
      T=298.15d0  !CFER was fitted with this temperature
      refnat = env%nat
      refchrg = env%chrg
      bhess = .true.

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
      !--- write the base input file
      if(env%ptb%h_acidic>0)then
        !-- for a manually selected H atom  
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
      else if(env%ptb%h_acidic==-1)then
        !-- "automatic" mode  
        call deprotonate(env,tim)
        call rmrfw('deprotonate_')
        call BASE%open('deprotonated.xyz')
      endif
      call BASE%write('base.xyz')
      env%chrg=refchrg-1
      env%nat = refnat - 1
      env%rednat = refnat - 1
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
      write(*,'(1x,a,f16.8,a,f8.2,a)') 'ΔG+Ecorr =',dG,' Eh,',dG*kcal,' kcal/mol'
      write(*,*)

      if(env%extLFER)then

      else    
         !pkaparam ='none'
         pkaparam = env%gfnver
      endif  
      select case( pkaparam )
        case( '--gfn1' )
          error stop 'GFN1-xTB not available for pKa calculation via CFER.'   
        case( '--gfn2' )
          c1=-1656.7    ! SD = 2.87, NO isorad, 19normal, my COSMOTHERM (out.ccf), adjust for Fabian's out.cosmo, ie -1656.45
          c2= 23.185    !    = 2.47 without two extreme outliers (CH3NN+, cycnoform)
          c3=-0.11103   !           and checked outlier HClO4
          c4= 0.0001835
        case( 'external' )
          c1 = 1.0d0  !placeholder
          c2 = 1.0d0  !placeholder 
          c3 = 1.0d0  !placeholder
          c3 = 1.0d0  !placeholder
        case default  !do nothing with the values  
          c1 = 0.0d0
          c2 = 1.0d0
          c3 = 0.0d0
          c4 = 0.0d0
     end select  

      write(*,'(1x,a)') 'cubic free energy relationship (CFER):'
      write(*,'(2x,a)') 'pKa   = c1 + c2*kdiss + c3*kdiss² + c4*kdiss³'
      write(*,'(2x,a)') 'with kdiss = ΔG(aq)/ln(10)RT'
      write(*,'(7x,a,f8.2,a)') 'ΔG(aq)=',dG*kcal,' kcal/mol'
      write(*,'(7x,a,f16.8)') 'c1 =',c1
      write(*,'(7x,a,f16.8)') 'c2 =',c2
      write(*,'(7x,a,f16.8)') 'c3 =',c3
      write(*,'(7x,a,f16.8)') 'c4 =',c4

      pka = pKaCFER(dG,c1,c2,c3,c4,T)

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
      write(*,'(3x,a)') '|___________________________________|'

!--- cleanup
      call ACID%deallocate()
      call BASE%deallocate()

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
       call env%wrtCHRG('')
       write(jobcall,'(a,1x,a,1x,a,'' --ceasefiles --opt '',a,1x,a,'' >xtb.out'')') &
       &    trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),' 2>/dev/null'
       call execute_command_line(trim(jobcall), exitstat=io)
       inquire(file='xtbopt.xyz',exist=io)
       if(io)then
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
    h=0
    arg=trim(str)
    if(arg=='auto')then
        h=-1
        return
    endif
    read(arg,*,iostat=io) rdum
    if((io==0).and.(rdum>0.0d0)) h = nint(rdum)
    return
end subroutine pka_argparse

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
    real(wp) :: GA,GB
    real(wp) :: dG         !in Eh
    real(wp),optional :: T !in K
    real(wp) :: pka
    real(wp),parameter :: kcal =627.5095d0
    real(wp),parameter :: R = 0.00198720425  !in kcal/K mol
    real(wp) :: logk
    real(wp) :: c1,c2,c3,c4
    if(.not.present(T))then
    T=298.15_wp 
    endif
    pka = 0.0d0
    !logk = kcal*dG/0.592452/2.302585
    logk =kcal*dG/(log(10.0d0)*R*T)
    pka = c1 + c2*logk + c3*(logk**2) + c4*(logk**3)
    return
end function pKaCFER
