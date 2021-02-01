!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020 Stefan Grimme
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
! Subroutine to calculate the molecular entropy of the system
!
! On Input:  nconf     - number of conformers
!            energies  - energy of each conformer
!            nrot      - number of rotamers for each conformer 
!            sym       - the point group of each conformer
!            rotc      - rotational constants A,B,C for each 
!                        conformer in MHz
!            occurence - array to trac which conformer stems
!                        from which timeframe
!            corefac   - degeneracy arising from enantiomers
!            tlen      - total MTD length in ps
!            T         - temperature in K
!            pr        - boolean for printout
!
! On Output: S  - molecular entropy
!            Cp - molecular heat capacity
!
!==============================================================!
subroutine calculateEntropy(nconf_tot,energies,nrot,sym, &
    &                       corefac,T,S,Cp,Hconf,pr,pr2)
     use iso_fortran_env, wp => real64, idp => int64
     implicit none

     integer,intent(in)  :: nconf_tot
     real(wp),intent(in) :: energies(nconf_tot)
     integer(idp),intent(in)  :: nrot(nconf_tot)
     character(len=3)    :: sym(nconf_tot)
     !real(wp),intent(in) :: rotc(3,nconf_tot)
     real(wp),intent(in) :: corefac(nconf_tot)
     logical,intent(in)  :: pr
     logical,intent(in) :: pr2

     real(wp) :: S, T
     real(wp) :: srrho, sconf, hconf, cp
     real(wp) :: dum1,dum2,stest,serr1,serr2,w_extra,x
     real(wp) :: rot(3)
     real(wp),allocatable  :: e (:)
     real(wp),allocatable  :: sr(:)
     real(wp),allocatable  :: g (:)
     real(wp),allocatable  :: p (:)
     integer, allocatable  :: ncslot(:)

     character(len=16) :: nmbr
     character(len=26) :: prnt
     integer :: i,j,k,l,ii
     integer :: ich
     integer :: nconf

     logical :: ex

     real(wp),parameter ::expo=0

      if(pr)then
        write(*,*)
        call smallhead('CALCULATION OF MOLECULAR ENTROPY')
        write(*,*) 'entropies in cal/mol K, (free)energies in kcal/mol'
      endif


      nconf = nconf_tot
      if(nconf.lt.1) stop 'nconf < 1'
      allocate(e(nconf),g(nconf),sr(nconf),p(nconf))

      i = 0
      do ii=1,nconf_tot
        i = i + 1
        g(i) = float(nrot(ii))
        e(i)=627.509541_wp*(energies(ii)-energies(1)) ! in kcal/mol
      enddo
!shift
      e = e - e(1)

! calc conformer populations 
! this includes degeneracies true g()
      call entropy_boltz(nconf,t,e,g,p)

! calc S_conf
      !call entropy_S(nconf,t,expo, e,g,sconf,dum1,dum2) ! g'
      call entropy_S(nconf,t,1.0d0,e,corefac,sconf,dum1,dum2) !"enantiomer" factor
      call entropy_S(nconf,t,1.0d0,e,g,dum1, cp,hconf)  ! g

      S = sconf

      if(pr2)then    
         write(nmbr,'(f16.2)') T
         write(prnt,'(a,a,a)') 'H(',trim(adjustl(nmbr)),'K)-H(0) (conf)'
         write(*,'(1x,a26,2x,f12.6)') adjustr(prnt),hconf
         write(*,'(1x,''                  Cp(conf)  '',f12.6)') cp     ! no contrib from geometry or symmetry change
         write(*,'(1x,''                   S(conf)  '',f12.6)') sconf
       endif

      deallocate(e,g,sr,p)

      if(pr) then
        write(*,'(/,1x,''final cp(total) =            '',f12.6)') cp
        write(*,'(/,1x,''final S (total) =            '',f12.6)') S
        write(*,'(/,1x,''CRE correction to G =        '',f12.6)') -S * T / 1000.0_wp + hconf
      endif

     return
end subroutine calculateEntropy

!========================================================!
! compute Boltzmann populations p for given
! T, level energies e(), and 
! level degeneracies g()  
!========================================================!
subroutine entropy_boltz(n,t,e,g,p)
      use iso_fortran_env, wp => real64
      implicit none
      integer n
      real(wp) :: e(n),p(n),g(n)          
      real(wp) :: t,f,esum
      integer i
      real(wp),parameter :: R = 8.31446261815324/4.184
       
      f = 1.0d0 / (T * R * 0.001d0)

      esum=0
      do i=1,n
         esum=esum+g(i)*exp(-e(i)*f)
      enddo
      do i=1,n
         p(i)=g(i)*exp(-e(i)*f)/esum
      enddo

      if(abs(1.0d0-sum(p)).gt.1.d-6) stop 'error in entropy_boltz()'

end subroutine entropy_boltz

!========================================================!
! compute entropy in cal/mol K for given T
! effective level energies e(), and
! level degeneracies g()  
! returns also heat capacity cp(conf) and thermal 
! correction H-H0(conf)
!========================================================!
subroutine entropy_S(n,T,expo,e,g,s,cp,htmh0)
      use iso_fortran_env, wp => real64
      implicit none
      integer :: n
      real(wp) :: e(n),g(n),s,t,expo,cp,htmh0


      real(wp),parameter :: R  = 8.31446261815324/4.184
      real(wp) :: f,esum,esum2,esum3,gi,pi,sump
      integer :: i
      
      f = 1.0d0 / (T * R * 0.001d0)

      esum =0
      esum2=0
      esum3=0
      do i=1,n
         gi = g(i)**expo      
         pi   =gi*exp(-e(i)*f)
         esum =esum +pi
         esum2=esum2+pi*(e(i)*f)
         esum3=esum3+pi*(e(i)*f)**2
      enddo

!     s    =0
!     sump =0
!     do i=1,n
!        gi = g(i)**expo      
!        pi =gi*exp(-e(i)*f)/esum
!        sump=sump+pi
!        s  =s + pi*log(pi+1.d-14)
!     enddo
!     if(abs(1.0d0-sump).gt.1.d-6) stop 'error in entropy_S()'
!     s = -s * R
!     write(*,*) s  ! identical to eq. below for all gi=1

      s = R * ( log(esum) + esum2/esum ) ! Eq. 30 from https://cccbdb.nist.gov/thermo.asp

      cp= R * ( esum3/esum - (esum2/esum)**2 ) ! Eq. 31 from https://cccbdb.nist.gov/thermo.asp

      htmh0 = T * R * esum2/esum * 0.001d0 ! Eq. 32 from https://cccbdb.nist.gov/thermo.asp

      return
end subroutine entropy_S

!========================================================!
! entropy calculations for documented files
!========================================================!
subroutine newentropyextrapol(env)
    use iso_fortran_env, wp => real64
    use crest_data
    use strucrd
    implicit none
    !type(options) :: opt
    type(systemdata) :: env

    integer :: i,j,k,l
    integer :: iter,nt,itry
    real(wp) :: T
    real(wp),allocatable :: nconflist(:)
    real(wp),allocatable :: cplist(:,:)
    real(wp),allocatable :: slist(:,:)
    real(wp),allocatable :: hlist(:,:)
    real(wp),allocatable :: sfinal(:)
    real(wp),allocatable :: cpfinal(:)
    real(wp),allocatable :: hfinal(:)
    real(wp),allocatable :: srrho(:)

    character(len=128) :: btmp
    character(len=64) :: btmp2
    character(len=16) :: nmbr
    character(len=26) :: prnt
    real(wp) :: s,cp
    real(wp) :: avs,dum1,dum2,val
    real(wp) :: p3,pmin,rmsd,rmin
    logical  :: ok,okmin
    integer  :: nall
    integer  :: rt,ri  
    integer :: ich

    logical :: prt,wrdegen
    real(wp),allocatable :: temps(:)
    real(wp),allocatable :: adum(:)

    !--- determine temperature dependence
    if(.not.allocated(env%thermo%temps))then
      call env%thermo%get_temps()
    endif
    nt=env%thermo%ntemps
    allocate(temps(nt))
    temps = env%thermo%temps
    T = temps(1)

    if(.not.allocated(env%emtd%soft))then
      allocate(env%emtd%soft(nt), source=0.0d0)
    endif
    if(.not.allocated(env%emtd%cpoft))then
      allocate(env%emtd%cpoft(nt), source=0.0d0)
    endif
    if(.not.allocated(env%emtd%hoft))then
      allocate(env%emtd%hoft(nt), source=0.0d0)
    endif

    iter = env%emtd%iterlast + 1   !we always have the zero-th iteration, hence +1
    allocate(nconflist(iter), source=0.0d0)
    allocate(cplist(iter,nt), source=0.0d0)
    allocate(slist(iter,nt), source=0.0d0)
    allocate(hlist(iter,nt), source=0.0d0)

    !--- sorting settings
    env%cgf(2)=.false.
    env%confgo = .true.

    call largehead('ENSEMBLE ENTROPY CALCULATION')
    allocate(adum(nt))
    adum=abs(temps-298.15d0)
    rt = minloc(adum,1)  !temperature closest to 298.15 is the ref.
    deallocate(adum)

    T = 298.15d0
    wrdegen=.false.
    open(unit=111,file='.data')
    open(unit=222,file='.entropydata')
    write(222,*) iter,nt
    ITERATIONLOOP : do i=1,iter
     
     l=i-1  !ITER INCLUDES ZEROTH ITERATATION!
     write(*,'(a)')'======================================================' 
     write(btmp,'(a,i0,a)') 'crest_entropy_rotamer_',l,'.xyz'
     !env%ensemblename = trim(btmp)
 
     !write(*,'(1x,a,1x,a)') 'Analyzing file',trim(btmp)    
     !env%confgo = .true.
     !call newcregen(env,3)
     !call rmrf(trim(btmp)//'.sorted')
     !call rdensembleparam('crest_ensemble.xyz',k,nall)

     !-- new read .dataX files
     write(btmp2,'(a,i0)')'Sdata',l
     open(newunit=ich,file=trim(btmp2))
     read(ich,*) nall
     write(222,*) nall
     do j=1,nt
        read(ich,*)slist(i,j),cplist(i,j),hlist(i,j)
        write(222,'(3F18.10)')slist(i,j),cplist(i,j),hlist(i,j)
     enddo
     close(ich,status='delete')

     write(*,'(1x,a,i0,a)')'Containing ',nall,' conformers :'
     env%confgo = .false.    
     nconflist(i) = float(nall)
    
     !if(i==iter)wrdegen=.true.
     !call entropic(env,.false.,.true.,wrdegen,trim(btmp),T,S,Cp)
     !cplist(i,:)    = env%emtd%Cpoft(:)
     !slist(i,:)     = env%emtd%soft(:)
     !hlist(i,:)     = env%emtd%hoft(:)

     write(nmbr,'(f16.2)') T
     write(prnt,'(a,a,a)') 'H(',trim(adjustl(nmbr)),'K)-H(0) (conf)'
     write(*,'(1x,a26,2x,f12.6)') adjustr(prnt),hlist(i,rt)
     write(*,'(1x,''                  Cp(conf)  '',f12.6)') cplist(i,rt)
     write(*,'(1x,''                   S(conf)  '',f12.6)') slist(i,rt)
       
     write(111,'(i10,2F12.6)') nall,Slist(i,rt),Cplist(i,rt)
    enddo ITERATIONLOOP
    close(222)
    close(111)
    write(*,'(a,/)')'======================================================'

!--- Sconf extrapolation    
    call smallhead('FINAL CALCULATION OF MOLECULAR ENTROPY')
    write(*,*)
    allocate(cpfinal(nt),sfinal(nt),hfinal(nt),srrho(nt), source = 0.0d0)
    TLOOP : do k=1,nt
       T = temps(k)   
       prt =.false.
       if(k.eq.rt) prt=.true.  !printout only for RT
        
       if(iter.gt.3) then ! at least 3 points (sMTDs)
        do itry=iter,5,-1   ! try with decreasing number of points untily fit is valid
          p3   = 1.0d0
          pmin = 1.0d0
          rmin = 1.0d+42
          okmin= .false.
          do ri=1,51  ! constrained search for best p3
            call marqfit(.false.,.false.,iter,itry,slist(1:iter,k),p3,val,rmsd,ok) ! (last) 10 points max
            if(rmsd.lt.rmin .and. ok) then
               pmin = p3
               rmin = rmsd
               okmin= .true.
            endif
            p3 = p3 - 0.01d0 
          enddo
          if(okmin) exit
        enddo
        call marqfit(.false.,prt,iter,itry,slist(1:iter,k),pmin,val,rmsd,ok) ! final with best p3
        if (.not. ok) then
           val = maxval(slist(1:iter,k))
           if(prt) then
              write(*,'(/,''WARNING: large fit error'',/)')
              write(*,'(''FINAL ENTROPY (max conf.)  : '',f12.5)') val
           endif
        else
           if(prt) write(*,'(''FINAL ENTROPY (extrapol.)  : '',f12.5)') val
        endif
       else
        val = maxval(slist(1:iter,k))
        if(prt) then
          write(*,'(/,''WARNING: not enough blocks'',/)')
          write(*,'(''FINAL ENTROPY (max conf.)  : '',f12.5)') val
        endif
       endif
       sfinal(k) = val
    enddo TLOOP


!--- Srrho correction to Sconf
    call calcSrrhoav(env,'crest_conformers.xyz')
    if(env%nreset > 0)then
      write(*,*)  
      write(*,'(1x,a)') 'Attention:'
      write(*,'(1x,a,a)') 'The algorithm has found a better global ', &
      & 'minimum structure and restarted at least once.'
      write(*,'(1x,a)') 'The associated restart can strongly affect the δSrrho values.'
    endif
    srrho(:) = env%emtd%soft(:)
    !--- cpfinal and hfinal come also from calcSrrhoav, but only
    !    to use populations based on the free energies.
    !cpfinal(:) = cplist(iter,:)
    cpfinal(:) = env%emtd%cpoft(:)
    !hfinal(:) = hlist(iter,:)
    hfinal(:) = env%emtd%hoft(:)
    do i=1,nt
       sfinal(i) = sfinal(i) + srrho(i)
    enddo
    
!---- temperature dependence printout
    if(nt > 1)then
       write(*,*)
       call smallhead('TEMPERATURE DEPENDENCE')
       write(*,*)
       write(*,'(1x,a)') 'Final CONFORMATIONAL quantities at given T:'
       write(*,'(a10)',advance='no') "T/K"
       write(*,'(a16)',advance='no') "S(total)"
       write(*,'(a16)',advance='no') "Cp(T)"
       write(*,'(a16)',advance='no') "[H(T)-H(0)]"
       write(*,'(a16)',advance='no') "G(total)"
       write(*,'(a)')
       write(*,'(3x,71("-"))')
       do i = 1, nt
          write(*,'(f10.2)',advance='no') temps(i)
          write(*,'(f16.6)',advance='no') sfinal(i)
          write(*,'(f16.6)',advance='no') cpfinal(i)
          write(*,'(f16.6)',advance='no') hfinal(i)
          dum2 =  -sfinal(i)*temps(i)/1000.0_wp + hfinal(i)
          write(*,'(f16.6)',advance='no') dum2
          write(*,'(a)')
       enddo
       write(*,'(3x,71("-"))')
       write(*,'(3x,a)')'S and Cp in cal/mol*K; H and G in kcal/mol'
       write(*,'(3x,a,a)')'G(total) is the ensemble free energy', &
       &    ' and S(total) = S(conf,extrapol.) + δSrrho'
    endif

!---- last summary print
    T = temps(rt)
    call entropyprintout(T,srrho(rt),sfinal(rt),cpfinal(rt),hfinal(rt))

    deallocate(srrho,hfinal,sfinal,cpfinal)
    deallocate(hlist)
    deallocate(slist)
    deallocate(cplist)
    deallocate(nconflist)

    return
end subroutine newentropyextrapol

subroutine entropyprintout(T,Srrho,S,Cp,H)
    use iso_fortran_env, wp => real64
    implicit none
    real(wp) :: T
    real(wp) :: Srrho
    real(wp) :: S
    real(wp) :: Cp
    real(wp) :: H
    character(len=80) :: atmp,btmp

    write(*,*)
    write(atmp,'(f12.2)')T
    write(btmp,'(a,1x,a,1x,a)')'| FINAL MOLECULAR ENTROPY AT T=',trim(adjustl(atmp)),'K |'
    call smallhead(trim(btmp))

    write(*,'(1x,''  Sconf   =      '',f12.6)') S-Srrho
    write(*,'(1x,''+ δSrrho  =      '',f12.6)') Srrho
    write(*,'(1x,40("-"))')
    write(*,'(1x,''= S(total)  =      '',f12.6)') S
    write(*,'(3x,''G(total)  =      '',f12.6)') -S * T / 1000.0_wp + H
    write(*,'(3x,''H(T)-H(0) =      '',f12.6)') h
    write(*,'(3x,''Cp(total) =      '',f12.6)') cp
    
    write(*,*)

    return
end subroutine entropyprintout


!===============================================================!
! Helper routine to re-do the extrapolation for 
! a read-in file. The two corresponding file formats
! written by crest are .data and .entropydata (see code above)
! Use by:
! crest -redoextrapol .data
!===============================================================!
subroutine redo_extrapol(env,fname,rtin)
    use iso_fortran_env, wp=>real64
    use crest_data
    implicit none
    type(systemdata) :: env
    character(len=*) :: fname
    integer :: ich
    integer :: i,j,k,l
    integer :: iter,nt,rt,rtin
    real(wp) :: T
    real(wp),allocatable :: temps(:)
    real(wp),allocatable :: slist(:,:)
    integer,allocatable  :: nalls(:)
    real(wp),allocatable :: cpfinal(:)
    real(wp),allocatable :: sfinal(:)
    real(wp),allocatable :: hfinal(:)
    real(wp),allocatable :: srrho(:)
    integer :: itry,ri
    real(wp) :: p3,val,pmin,rmin,rmsd
    logical :: okmin,ok,prt
    character(len=80) :: atmp
    integer :: io
    real(wp) :: dum1,dum2

    rt=rtin
    call smallhead('entropy extrapolation from file '//trim(fname))
    write(*,*)

    if(trim(fname)=='.entropydata')then
       open(newunit=ich,file=fname)
       read(ich,*) iter,nt
       allocate(temps(nt),slist(iter,nt), source=0.0_wp)
       allocate(nalls(iter), source=0)
       do i=1,iter
         read(ich,*) nalls(i)
         do j=1,nt
         read(ich,*)slist(i,j),dum1,dum2
         enddo
       enddo
       close(ich)
       if(rtin==0) rt=3 !the third entry for each block is 298K in the crest default
    else if(trim(fname)=='.data')then
        open(newunit=ich,file=fname)
        iter=0
        do
          read(ich,'(a)',iostat=io) atmp
          if(io < 0) exit
          iter=iter+1 
        enddo  
        close(ich)
        nt=1
        rt=1  !.data only has the data for 298K
        allocate(temps(nt),slist(iter,nt), source=0.0_wp)
        allocate(nalls(iter), source=0)
        open(newunit=ich,file=fname)
        do i=1,iter
        read(ich,*)nalls(i),slist(i,nt),dum1
        enddo
        close(ich)
    else
       write(*,*) 'unknown file. must be either ".data" or ".entropydata"'
       error stop
    endif

    allocate(cpfinal(nt),sfinal(nt),hfinal(nt),srrho(nt), source = 0.0d0)
    TLOOP : do k=1,nt
    !      T = temps(k)   !for the read in version we dont need the temperatures
       prt =.false.
       if(k.eq.rt) prt=.true.  !printout only for RT

       if(iter.gt.3) then ! at least 3 points (sMTDs)
        do itry=iter,5,-1   ! try with decreasing number of points untily fit is valid
          p3   = 1.0d0
          pmin = 1.0d0
          rmin = 1.0d+42
          okmin= .false.
          do ri=1,51  ! constrained search for best p3
            call marqfit(.false.,.false.,iter,itry,slist(1:iter,k),p3,val,rmsd,ok) ! (last) 10 points max
            if(rmsd.lt.rmin .and. ok) then
               pmin = p3
               rmin = rmsd
               okmin= .true.
            endif
            p3 = p3 - 0.01d0
          enddo
          if(okmin) exit
        enddo
        call marqfit(prt,prt,iter,itry,slist(1:iter,k),pmin,val,rmsd,ok) ! final with best p3
        if (.not. ok) then
           val = maxval(slist(1:iter,k))
           if(prt) then
              write(*,'(/,''WARNING: large fit error'',/)')
              write(*,'(''FINAL ENTROPY (max conf.)  : '',f12.5)') val
           endif
        else
           if(prt) write(*,'(''FINAL ENTROPY (extrapol.)  : '',f12.5)') val
        endif
       else
        val = maxval(slist(1:iter,k))
        if(prt) then
          write(*,'(/,''WARNING: not enough blocks'',/)')
          write(*,'(''FINAL ENTROPY (max conf.)  : '',f12.5)') val
        endif
       endif
       sfinal(k) = val
    enddo TLOOP

    return
end subroutine redo_extrapol

