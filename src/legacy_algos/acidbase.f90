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

!================================================================================!
! Potential correction for acid base reactions,
! fitted to heterolytic dissociation reactions at DFT level
!================================================================================!
!---------------------------------------------------------!
! For a given acid and base input files, do singlepoints
! and obtain needed data (charges, WBO)
!---------------------------------------------------------!
subroutine acidbase(env,acidfile,basefile,acidchrg,verbose,keepdir,dE, &
        &  bhess,eatb,gsa,grrhoa,ebtb,gsb,grrhob)
      use crest_parameters
      use crest_data
      use iomod
      use strucrd
      use zdata
      implicit none
      type(systemdata) :: env
      character(len=*) :: acidfile
      character(len=*) :: basefile
      integer :: acidchrg
      real(wp) :: dE
      logical :: verbose,keepdir
      logical :: bhess

      type(coord) :: mola
      type(coord) :: molb
      type(zmolecule) :: acid
      type(zmolecule) :: base
      integer :: basechrg
      integer :: i,j,ich,r
      integer :: nata,natb              !number of atoms for acid/bas
      real(wp),allocatable :: wboa(:,:) !WBOs for the acid
      real(wp),allocatable :: wbob(:,:) !WBOs for the base
      real(wp),allocatable :: qa(:)     !atomic charges for the acid
      real(wp),allocatable :: qb(:)     !atomic charges for the base 
      real(wp) :: eatb,gsa,grrhoa
      real(wp) :: ebtb,gsb,grrhob
      real(wp) :: detb
      logical :: ldum
      integer :: X,atX,xH               !the reactive atom and its element
      real(wp) :: qh
      real(wp) :: dgrrho,edft
      real(wp) :: wbosuma,wbosumb
      integer :: ich2
      logical :: ex

      integer :: nqh
      character(len=512) :: thispath
      character(len=8),parameter :: adir = 'ACIDcalc'
      character(len=8),parameter :: bdir = 'BASEcalc'

      call getcwd(thispath)
!--- read structures and check atom numbers      
      call mola%open(trim(acidfile))
      call molb%open(trim(basefile))
      nata = mola%nat
      natb = molb%nat
      if(nata .ne. (natb+1))then
          error stop 'mismatching number of atoms in acid/base correction!'
      endif
      basechrg = acidchrg-1

!--- calculation (singlepoint) to get data for acid      
      eatb=0.0_wp
      gsa=0.0_wp
      grrhoa=0.0_wp
      call simpletopo_mol(mola,acid,.false.,.false.)
      r = makedir(adir)
      call chdir(adir)
      call mola%write('coord')
      call wrshort('.CHRG',acidchrg)
      call ab_singlepoint('coord',env,bhess)
      !-- read energies
      call grepval('xtb.out',"| TOTAL ENERGY",ldum,eatb)  !includes gsa
      call grepval('xtb.out',":: -> Gsolv",ldum,gsa)
      call grepval('xtb.out',":: G(RRHO) contrib.",ldum,grrhoa)
      allocate(wboa(nata,nata),qa(nata), source = 0.0_wp)
      call readwbo("wbo",nata, wboa)
      call ab_rdcharges(nata,qa)
      call chdir(thispath)
      if(.not.keepdir) call rmrf(adir)

!--- calculation (singlepoint) to get data for base
      ebtb=0.0_wp
      gsb=0.0_wp
      grrhob=0.0_wp
      call simpletopo_mol(molb,base,.false.,.false.)
      r = makedir(bdir)
      call chdir(bdir)
      call molb%write('coord')
      call wrshort('.CHRG',basechrg)
      call ab_singlepoint('coord',env,bhess)
      !-- read energies
      call grepval('xtb.out',"| TOTAL ENERGY",ldum,ebtb) !includes gsb
      call grepval('xtb.out',":: -> Gsolv",ldum,gsb)
      call grepval('xtb.out',":: G(RRHO) contrib.",ldum,grrhob)
      allocate(wbob(natb,natb),qb(natb), source = 0.0_wp)
      call readwbo("wbo",natb,wbob)
      call ab_rdcharges(natb,qb)
      call chdir(thispath)
      if(.not.keepdir) call rmrf(bdir)

!--- identify reactive atom X      
     call ab_reactivecenter(acid,base,X)
     if(verbose)then
         if(X > 0)then
         write(*,'(1x,a,i0,3a)')'Atom ',X,' (',i2e(acid%at(X),'nc'), &
         &    ') identified as reactive center.'
         else
         write(*,'(1x,a)')'No unambigous site for the acid/base reaction could be identified.'
         endif
     endif

!--- check is an atom was identified and determine charge of acidic H
     if(X>0)then
        atX = acid%at(X)  
        !--- if X has more than one acidic H, just average their charge
        qh = 0.0_wp
        nqh=0
        do i=1,acid%zat(X)%nei
          j=acid%zat(X)%ngh(i)
          if(acid%at(j)==1)then
             nqh=nqh+1
             qh = qh + qa(j)
             xH=j
          endif
        enddo
        qh = qh/float(nqh)
     else
        atX = 0
        qh = 0.0_wp
     endif

      open(newunit=ich2,file='.ATOM')
      write(ich2,*) xH
      close(ich2)
     
!--- calculate dE
     detb=ebtb-eatb
     if(verbose)then
        write(*,*) 'Gsolv(xtb) A / B    :',gsa,gsb
        write(*,*) 'dE(xtb,uncorr.)     :',detb
        if(bhess)then
        dgrrho=(grrhob-grrhoa)
        write(*,*) 'dG(xtb,uncorr.)     :',detb+(grrhob-grrhoa)
        endif
     endif
     call acidbasepot(dE,X,atX,nata,natb,qh,qa,qb,wboa,wbob)

!--- printout for fit     
    wbosuma = sum(wboa(:,X))
    wbosumb = sum(wbob(:,X))
    inquire(file='.edft',exist=ex)
    if(ex)then
        open(newunit=ich2,file='.edft')
        read(ich2,*) edft
        close(ich2)
    else
        edft=0.0_wp
    endif
    open(newunit=ich,file='.data')
    write(ich,'(2f16.8,4f10.6,f16.8)')edft,detb,(wbosuma-wbosumb),qa(X),qb(X),qh,dgrrho
    close(ich)

     if(verbose)then
        write(*,'(1x,a)')'Energy correction for acid/base reaction:'
        write(*,'(1x,a,f16.10,a)') 'Ecorr =',dE,' Eh'
     endif

     if(allocated(qb)) deallocate(qb)
     if(allocated(wbob))deallocate(wbob)
     if(allocated(qa)) deallocate(qa)
     if(allocated(wboa))deallocate(wboa)
     call acid%deallocate()
     call base%deallocate()
     call mola%deallocate()
     call molb%deallocate()
     return
end subroutine acidbase
subroutine ab_singlepoint(fname,env,bhess)
     use crest_data
     use iomod, only: command
     implicit none
     type(systemdata) :: env
     logical :: bhess
     character(len=*) :: fname
     character(len=1024) :: jobcall
     character(len=80),parameter :: pipe = ' > xtb.out 2>/dev/null'
     integer :: io
     if(.not.bhess)then
     write(jobcall,'(a,1x,a,1x,a,'' --wbo '',a,1x,a,a)') &
     &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),trim(pipe)
     else
     write(jobcall,'(a,1x,a,1x,a,'' --wbo --bhess '',a,1x,a,a)')  &
     &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),trim(pipe)
     endif
     call command(trim(jobcall), io)
     return
end subroutine ab_singlepoint
subroutine ab_reactivecenter(zmola,zmolb,X)
    use zdata
    implicit none
    type(zmolecule) :: zmola  !topology of the acid
    type(zmolecule) :: zmolb  !topology of the base
    integer,intent(out) :: X
    integer :: i,j,n
    integer :: na,nb
    integer :: ha,hb
    X = 0
    n = zmolb%nat  !the base has less atoms
    do i=1,n
      if(zmola%at(i) == 1) cycle
      na = zmola%zat(i)%nei
      nb = zmolb%zat(i)%nei
      if(na == (nb+1))then !Atom i in the acid has one more neighbor
        ha = 0
        hb = 0
        do j=1,na
          if(zmola%zat(i)%ngt(j) == 1) ha=ha+1
        enddo
        do j=1,nb
          if(zmolb%zat(i)%ngt(j) == 1) hb=hb+1
        enddo
        if(ha == (hb+1))then
           X = i 
           exit
        endif
      endif
    enddo
    return
end subroutine ab_reactivecenter
subroutine ab_rdcharges(nat,q)
    use crest_parameters, only: wp
    implicit none
    integer,intent(in) :: nat
    real(wp),intent(out) :: q(nat)
    integer :: i,ich
    logical :: ex
    inquire(file='charges',exist=ex)
    if(.not.ex)then
        q = 0.0_wp
        return
    endif
    open(newunit=ich,file='charges')
    do i=1,nat
      read(ich,*) q(i)
    enddo
    close(ich)
    return
end subroutine ab_rdcharges


!--------------------------------------------------------------!
! Implementation of the potential correction
! dE  - output energy
! X   - (heavy) Atom of the reactive position
! atX - atom type of X
! nata - number of atoms in the acid
! natb - number of atoms in the base
! qh   - charge of the hydrogen atom at the acidic atom
! qa   - atomic charges in the acid
! qb   - atomic charges in the base
! wboa - bond orders in the acid
! wbob - bond orders in the base
!--------------------------------------------------------------!
subroutine acidbasepot(dE,X,atX,nata,natb,qh,qa,qb,wboa,wbob)
    use crest_parameters, only: wp
    implicit none
    real(wp),intent(out) :: dE
    integer :: X
    integer :: atX
    integer :: nata
    integer :: natb
    real(wp) :: qh
    real(wp) :: qa(nata)
    real(wp) :: qb(natb)
    real(wp) :: wboa(nata,nata)
    real(wp) :: wbob(natb,natb)
    real(wp) :: wbosuma,wbosumb
    dE = 0.0_wp
    dE = abparam(1,atX)
    if(X < 1) return
    wbosuma = sum(wboa(X,:))
    wbosumb = sum(wbob(X,:))
    dE = dE + abparam(2,atX)*(wbosuma-wbosumb)
    dE = dE + abparam(3,atX)*qa(X)
    dE = dE + abparam(4,atX)*qb(X)
    dE = dE + abparam(5,atX)*qh
    return
contains
function abparam(i,j) result(val)
    implicit none
    integer :: i,j
    real(wp) :: val
    real(wp) :: dat(6)

    dat = (/0.25_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp/)

    select case( j )
      case ( 6 ) !Carbon C
          dat = (/0.2432399368_wp,0.0068192940_wp,-0.0125805950_wp, &
              &   0.1548006770_wp,0.1430153999_wp,0.0_wp/)
      case( 7 ) !Nitrogen N
          dat = (/0.2570736941_wp,0.0371207152_wp,-0.0557745938_wp,  &
              &   0.1067342982_wp,0.0238365637_wp,0.0_wp/)
      case( 8 ) !Oxygen O
          dat = (/0.1946176893_wp,0.0240456714_wp,0.0282362624_wp, &
              &   0.0148248415_wp,0.1496774031_wp,0.0_wp/)
      case( 9 ) !Fluorine F
          dat = (/0.1222735061_wp,-0.1470616565_wp,0.3283530541_wp, &
              &   -0.3275311076_wp,0.2341983828_wp,0.0_wp/)
      case( 14 ) !Silicon Si
          dat = (/0.2128520556_wp,0.0190478003_wp,-0.0612287051_wp, &
              &   0.0260659207_wp,0.0088543551_wp,0.0_wp/)
      case( 15 ) !Phosporus P
          dat = (/0.2749418789_wp,-0.0234033348_wp,0.0257283686_wp, &
              &   -0.0136988930_wp,-0.0267469314_wp,0.0_wp/)
      case( 16 ) !Sulfur S
          dat = (/0.3010248722_wp,0.0082633594_wp,-0.0806268699_wp, &
              &   0.0939537740_wp,-0.0515873768_wp,0.0_wp/)
      case( 17 ) !Chloride Cl
          dat = (/0.2787361618_wp,-0.0088864104_wp,-0.0083718474_wp, &
              &   0.0140575174_wp,-0.0577563975_wp,0.0_wp/)
      case default
         val = dat(i) 
    end select
    val = dat(i)
    return
end function abparam    
end subroutine acidbasepot

!---------------------------------------------------------!
! A helper routine to generate ensemble files at the
! GFN level containing free energies including the
! acid/base energy correction as comment.
! E_corr will be added to the Base ensemble and 
! reference the lowest Acid structure (w.r.t. WBO and q(X))
!---------------------------------------------------------!
subroutine rewrite_AB_ensemble(env,acensemble,baensemble)
      use crest_parameters, only: wp, bohr
      use crest_data
      use strucrd
      use iomod
      implicit none

      type(systemdata) :: env
      character(len=*) :: acensemble
      character(len=*) :: baensemble
      type(ensemble) :: ACIDENSEMBLE
      type(ensemble) :: BASEENSEMBLE
      real(wp) :: t
      real(wp) :: eatb,gsa,grrhoa,ebtb,gsb,grrhob,dE
      integer :: i,j,l
      logical :: bhess,ldum
      integer :: nalla
      real(wp),allocatable :: ae(:)
      integer :: nallb
      integer :: achrg,bchrg
      real(wp),allocatable :: be(:)
      character(len=512) :: thispath
      character(len=128) :: tmppath
      real(wp),parameter :: kcal =627.5095d0

      T=env%tboltz
      achrg=env%chrg
      bchrg=env%chrg-1
      bhess=.true.

      call new_ompautoset(env,'auto',1,i,j)

!-- first, read the acid ensemble and calculate the free energies
      write(*,*) 'Acid ensemble: ',trim(acensemble)
      call ACIDENSEMBLE%open(acensemble)
      nalla = ACIDENSEMBLE%nall
      allocate(ae(nalla),source=0.0_wp)
      call getcwd(thispath)
      write(*,'(1x,a9,a16,a16,a16,a16)')'structure','E_read','E_calc','E_corr','G_calc'
      do i=1,nalla
      write(tmppath,'(a,i0)')'AC',i
      l = makedir(trim(tmppath))
      call chdir(tmppath)
      call env%wrtCHRG('')
      call wrcoord('coord',ACIDENSEMBLE%nat,ACIDENSEMBLE%at,ACIDENSEMBLE%xyz(:,:,i)/bohr)
      call ab_singlepoint('coord',env,bhess)
      !-- read energies
      call grepval('xtb.out',"| TOTAL ENERGY",ldum,eatb)  !includes gsa
      call grepval('xtb.out',":: -> Gsolv",ldum,gsa)
      call grepval('xtb.out',":: G(RRHO) contrib.",ldum,grrhoa)
      call chdir(thispath)
      call rmrf(trim(tmppath))
      ae(i) = eatb + grrhoa
      write(*,'(1x,i9,2f16.8,a16,f16.8)') i,ACIDENSEMBLE%er(i),eatb,'---',ae(i)
      enddo
      ACIDENSEMBLE%er = ae  
      call ACIDENSEMBLE%write('G_'//trim(acensemble))
      write(*,'(1x,a,a)')'written to: ','G_'//trim(acensemble)


!-- write a reference acid structure (with the lowest energy)
      j=minloc(ae,1)
      write(tmppath,'(f18.8)')ae(j)
      call wrxyz('acid_ref.xyz',ACIDENSEMBLE%nat,ACIDENSEMBLE%at,ACIDENSEMBLE%xyz(:,:,j),tmppath)
      call ACIDENSEMBLE%deallocate()

!--- loop over Base ensemble and calculate G including Ecorr
      write(*,*)
      write(*,*) 'Base ensemble: ',trim(baensemble)
      call BASEENSEMBLE%open(baensemble)
      nallb = BASEENSEMBLE%nall
      allocate(be(nallb),source=0.0_wp)
      call getcwd(thispath)
      write(*,'(1x,a9,a16,a16,a16,a16)')'structure','E_read','E_calc','E_corr','G_calc'
      do i=1,nallb
      call wrxyz('btmp.xyz',BASEENSEMBLE%nat,BASEENSEMBLE%at,BASEENSEMBLE%xyz(:,:,i))
      call acidbase(env,'acid_ref.xyz','btmp.xyz',achrg,.false.,.false.,dE, &
        &  bhess,eatb,gsa,grrhoa,ebtb,gsb,grrhob)
      be(i) = ebtb + grrhob + dE
      write(*,'(1x,i9,4f16.8)') i,BASEENSEMBLE%er(i),ebtb,dE,be(i)
      enddo
      BASEENSEMBLE%er = be
      call BASEENSEMBLE%write('G_'//trim(baensemble))
      call BASEENSEMBLE%deallocate()
      call rmrf('acid_ref.xyz')
      call rmrf('btmp.xyz')
      write(*,'(1x,a,a)')'written to: ','G_'//trim(baensemble)

      return
end subroutine rewrite_AB_ensemble
