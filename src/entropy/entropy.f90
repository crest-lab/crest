!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020-2023 Stefan Grimme, Philipp Pracht
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

subroutine calculateEntropy(nconf_tot,energies,nrot, &
    &                       corefac,T,S,Cp,Hconf,pr,pr2)
!****************************************************************
!* Subroutine to calculate the molecular entropy of the system
!*
!* On Input:  nconf     - number of conformers
!*            energies  - energy of each conformer
!*            nrot      - number of rotamers for each conformer
!*            corefac   - degeneracy arising from enantiomers
!*            tlen      - total MTD length in ps
!*            T         - temperature in K
!*            pr        - boolean for printout
!*
!* On Output: S  - molecular entropy
!*            Cp - molecular heat capacity
!* 
!****************************************************************
  use iso_fortran_env,wp => real64,idp => int64
  implicit none

  integer,intent(in)  :: nconf_tot
  real(wp),intent(in) :: energies(nconf_tot)
  integer(idp),intent(in)  :: nrot(nconf_tot)
  real(wp),intent(in) :: corefac(nconf_tot)
  logical,intent(in)  :: pr
  logical,intent(in) :: pr2

  real(wp) :: S,T
  real(wp) :: sconf,hconf,cp
  real(wp) :: dum1,dum2
  real(wp),allocatable  :: e(:)
  real(wp),allocatable  :: sr(:)
  real(wp),allocatable  :: g(:)
  real(wp),allocatable  :: p(:)

  character(len=16) :: nmbr
  character(len=26) :: prnt
  integer :: i,ii
  integer :: nconf

  real(wp),parameter ::expo = 0

  if (pr) then
    write (*,*)
    call smallhead('CALCULATION OF MOLECULAR ENTROPY')
    write (*,*) 'entropies in cal/mol K, (free)energies in kcal/mol'
  end if

  nconf = nconf_tot
  if (nconf .lt. 1) stop 'nconf < 1'
  allocate (e(nconf),g(nconf),sr(nconf),p(nconf))

  i = 0
  do ii = 1,nconf_tot
    i = i+1
    g(i) = float(nrot(ii))
    e(i) = 627.509541_wp*(energies(ii)-energies(1)) ! in kcal/mol
  end do
!shift
  e = e-e(1)

! calc conformer populations
! this includes degeneracies true g()
  call entropy_boltz(nconf,t,e,g,p)

! calc S_conf
  !call entropy_S(nconf,t,expo, e,g,sconf,dum1,dum2) ! g'
  call entropy_S(nconf,t,1.0d0,e,corefac,sconf,dum1,dum2) !"enantiomer" factor
  call entropy_S(nconf,t,1.0d0,e,g,dum1,cp,hconf)  ! g

  S = sconf

  if (pr2) then
    write (nmbr,'(f16.2)') T
    write (prnt,'(a,a,a)') 'H(',trim(adjustl(nmbr)),'K)-H(0) (conf)'
    write (*,'(1x,a26,2x,f12.6)') adjustr(prnt),hconf
    write (*,'(1x,''                  Cp(conf)  '',f12.6)') cp     ! no contrib from geometry or symmetry change
    write (*,'(1x,''                   S(conf)  '',f12.6)') sconf
  end if

  deallocate (e,g,sr,p)

  if (pr) then
    write (*,'(/,1x,''final cp(total) =            '',f12.6)') cp
    write (*,'(/,1x,''final S (total) =            '',f12.6)') S
    write (*,'(/,1x,''CRE correction to G =        '',f12.6)')-S*T/1000.0_wp+hconf
  end if

  return
end subroutine calculateEntropy

!=========================================================================================!
subroutine entropy_boltz(n,t,e,g,p)
!***************************************************
!* compute Boltzmann populations p for given
!* T, level energies e(), and
!* level degeneracies g()
!***************************************************
  use crest_parameters
  implicit none
  integer n
  real(wp) :: e(n),p(n),g(n)
  real(wp) :: t,f,esum
  integer i
  real(wp),parameter :: R = Rcal

  f = 1.0d0/(T*R*0.001d0)

  esum = 0
  do i = 1,n
    esum = esum+g(i)*exp(-e(i)*f)
  end do
  do i = 1,n
    p(i) = g(i)*exp(-e(i)*f)/esum
  end do

  if (abs(1.0d0-sum(p)) .gt. 1.d-6) stop 'error in entropy_boltz()'

end subroutine entropy_boltz

!=========================================================================================!
subroutine entropy_S(n,T,expo,e,g,s,cp,htmh0)
!***************************************************
!* compute entropy in cal/mol K for given T
!* effective level energies e(), and
!* level degeneracies g()
!* returns also heat capacity cp(conf) and thermal
!* correction H-H0(conf)
!***************************************************
  use crest_parameters
  implicit none
  integer :: n
  real(wp) :: e(n),g(n),s,t,expo,cp,htmh0

  real(wp),parameter :: R = Rcal
  real(wp) :: f,esum,esum2,esum3,gi,p_i
  integer :: i

  f = 1.0d0/(T*R*0.001d0)

  esum = 0
  esum2 = 0
  esum3 = 0
  do i = 1,n
    gi = g(i)**expo
    p_i = gi*exp(-e(i)*f)
    esum = esum+p_i
    esum2 = esum2+p_i*(e(i)*f)
    esum3 = esum3+p_i*(e(i)*f)**2
  end do

  !> Eq. 30 from https://cccbdb.nist.gov/thermo.asp
  s = R*(log(esum)+esum2/esum) 

  !> Eq. 31 from https://cccbdb.nist.gov/thermo.asp
  cp = R*(esum3/esum-(esum2/esum)**2) 

  !> Eq. 32 from https://cccbdb.nist.gov/thermo.asp
  htmh0 = T*R*esum2/esum*0.001d0 

  return
end subroutine entropy_S

!=========================================================================================!
subroutine newentropyextrapol(env)
!*********************************************************
!* Entropy extrapolation.
!*********************************************************
  use crest_parameters
  use crest_data
  use strucrd
  implicit none
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
  real(wp) :: dum2,val
  real(wp) :: p3,pmin,rmsd,rmin
  logical  :: ok,okmin
  integer  :: nall
  integer  :: rt,ri
  integer :: ich

  logical :: prt,wrdegen
  real(wp),allocatable :: temps(:)
  real(wp),allocatable :: adum(:)

!>--- determine temperature dependence
  if (.not.allocated(env%thermo%temps)) then
    call env%thermo%get_temps()
  end if
  nt = env%thermo%ntemps
  allocate (temps(nt))
  temps = env%thermo%temps
  T = temps(1)

  if (.not.allocated(env%emtd%soft)) then
    allocate (env%emtd%soft(nt),source=0.0d0)
  end if
  if (.not.allocated(env%emtd%cpoft)) then
    allocate (env%emtd%cpoft(nt),source=0.0d0)
  end if
  if (.not.allocated(env%emtd%hoft)) then
    allocate (env%emtd%hoft(nt),source=0.0d0)
  end if

  iter = env%emtd%iterlast+1   !> we always have the zero-th iteration, hence +1
  allocate (nconflist(iter),source=0.0d0)
  allocate (cplist(iter,nt),source=0.0d0)
  allocate (slist(iter,nt),source=0.0d0)
  allocate (hlist(iter,nt),source=0.0d0)

  !>--- sorting settings
  env%cgf(2) = .false.
  env%confgo = .true.

  call largehead('ENSEMBLE ENTROPY CALCULATION')
  allocate (adum(nt))
  adum = abs(temps-298.15d0)
  rt = minloc(adum,1)  !temperature closest to 298.15 is the ref.
  deallocate (adum)

  !>--- fallback for diatomics
  if (env%nat .le. 2) then
    T = temps(rt)
    allocate (cpfinal(nt),sfinal(nt),hfinal(nt),srrho(nt),source=0.0d0)
    write (*,*) 'Note: di- or monoatomic molecules have no conformational entropy'
    call entropyprintout(T,srrho(rt),sfinal(rt),cpfinal(rt),hfinal(rt))
    deallocate (srrho,hfinal,sfinal,cpfinal)
    deallocate (hlist,slist,cplist,nconflist)
    deallocate (env%emtd%hoft,env%emtd%cpoft,env%emtd%soft,temps)
    return
  end if

  T = 298.15d0
  wrdegen = .false.
  open (unit=111,file='.data')
  open (unit=222,file='.entropydata')
  write (222,*) iter,nt
  ITERATIONLOOP: do i = 1,iter

    l = i-1  !ITER INCLUDES ZEROTH ITERATATION!
    write (*,'(a)') '======================================================'
    write (btmp,'(a,i0,a)') 'crest_entropy_rotamer_',l,'.xyz'

    !-- new read .dataX files
    write (btmp2,'(a,i0)') 'Sdata',l
    open (newunit=ich,file=trim(btmp2))
    read (ich,*) nall
    write (222,*) nall
    do j = 1,nt
      read (ich,*) slist(i,j),cplist(i,j),hlist(i,j)
      write (222,'(3F18.10)') slist(i,j),cplist(i,j),hlist(i,j)
    end do
    close (ich,status='delete')

    write (*,'(1x,a,i0,a)') 'Containing ',nall,' conformers :'
    env%confgo = .false.
    nconflist(i) = float(nall)

    write (nmbr,'(f16.2)') T
    write (prnt,'(a,a,a)') 'H(',trim(adjustl(nmbr)),'K)-H(0) (conf)'
    write (*,'(1x,a26,2x,f12.6)') adjustr(prnt),hlist(i,rt)
    write (*,'(1x,''                  Cp(conf)  '',f12.6)') cplist(i,rt)
    write (*,'(1x,''                   S(conf)  '',f12.6)') slist(i,rt)

    write (111,'(i10,2F12.6)') nall,Slist(i,rt),Cplist(i,rt)
  end do ITERATIONLOOP
  close (222)
  close (111)
  write (*,'(a,/)') '======================================================'

!>--- Sconf extrapolation
  call smallhead('FINAL CALCULATION OF MOLECULAR ENTROPY')
  write (*,*)
  allocate (cpfinal(nt),sfinal(nt),hfinal(nt),srrho(nt),source=0.0d0)
  TLOOP: do k = 1,nt
    T = temps(k)
    prt = .false.
    if (k .eq. rt) prt = .true.  !printout only for RT

    if (iter .gt. 3) then ! at least 3 points (sMTDs)
      do itry = iter,5,-1   ! try with decreasing number of points untily fit is valid
        p3 = 1.0d0
        pmin = 1.0d0
        rmin = 1.0d+42
        okmin = .false.
        do ri = 1,51  ! constrained search for best p3
          call marqfit(.false.,.false.,iter,itry,slist(1:iter,k),p3,val,rmsd,ok) ! (last) 10 points max
          if (rmsd .lt. rmin.and.ok) then
            pmin = p3
            rmin = rmsd
            okmin = .true.
          end if
          p3 = p3-0.01d0
        end do
        if (okmin) exit
      end do
      call marqfit(.false.,prt,iter,itry,slist(1:iter,k),pmin,val,rmsd,ok) ! final with best p3
      if (.not.ok) then
        val = maxval(slist(1:iter,k))
        if (prt) then
          write (*,'(''WARNING: large fit error expected'',/)')
          write (*,'(''FINAL ENTROPY (max conf.)  : '',f12.5)') val
        end if
      else
        if (prt) write (*,'(''FINAL ENTROPY (extrapol.)  : '',f12.5)') val
      end if
    else
      val = maxval(slist(1:iter,k))
      if (prt) then
        write (*,'(''WARNING: Not enough blocks for extrapolation.'')')
        write (*,'(''         Taking full ensemble entropy instead.'',/)')
        write (*,'(''FINAL ENTROPY (max conf.)  : '',f12.5)') val
      end if
    end if
    sfinal(k) = val
  end do TLOOP

!>--- Srrho correction to Sconf
  call calcSrrhoav(env,'crest_conformers.xyz')
  if (env%nreset > 0) then
    write (*,*)
    write (*,'(1x,a)') 'Attention:'
    write (*,'(1x,a,a)') 'The algorithm has found a better global ', &
    & 'minimum structure and restarted at least once.'
    write (*,'(1x,a)') 'The associated restart can strongly affect the δSrrho values.'
  end if
  srrho(:) = env%emtd%soft(:)
!>--- cpfinal and hfinal come also from calcSrrhoav, but only
!>    to use populations based on the free energies.
  !cpfinal(:) = cplist(iter,:)
  cpfinal(:) = env%emtd%cpoft(:)
  !hfinal(:) = hlist(iter,:)
  hfinal(:) = env%emtd%hoft(:)
  do i = 1,nt
    sfinal(i) = sfinal(i)+srrho(i)
  end do

!>---- temperature dependence printout
  if (nt > 1) then
    write (*,*)
    call smallhead('TEMPERATURE DEPENDENCE')
    write (*,*)
    write (*,'(1x,a)') 'Final CONFORMATIONAL quantities at given T:'
    write (*,'(a10)',advance='no') "T/K"
    write (*,'(a16)',advance='no') "S(total)"
    write (*,'(a16)',advance='no') "Cp(T)"
    write (*,'(a16)',advance='no') "[H(T)-H(0)]"
    write (*,'(a16)',advance='no') "G(total)"
    write (*,'(a)')
    write (*,'(3x,71("-"))')
    do i = 1,nt
      write (*,'(f10.2)',advance='no') temps(i)
      write (*,'(f16.6)',advance='no') sfinal(i)
      write (*,'(f16.6)',advance='no') cpfinal(i)
      write (*,'(f16.6)',advance='no') hfinal(i)
      dum2 = -sfinal(i)*temps(i)/1000.0_wp+hfinal(i)
      write (*,'(f16.6)',advance='no') dum2
      write (*,'(a)')
    end do
    write (*,'(3x,71("-"))')
    write (*,'(3x,a)') 'S and Cp in cal/mol*K; H and G in kcal/mol'
    write (*,'(3x,a,a)') 'G(total) is the ensemble free energy (H-T*S)', &
    &    ' and S(total) = S(conf,extrapol.) + δSrrho'
  end if

!>---- last summary print
  T = temps(rt)
  call entropyprintout(T,srrho(rt),sfinal(rt),cpfinal(rt),hfinal(rt))

  deallocate (srrho,hfinal,sfinal,cpfinal)
  deallocate (hlist)
  deallocate (slist)
  deallocate (cplist)
  deallocate (nconflist)

  return
end subroutine newentropyextrapol

!========================================================================================!
subroutine entropyprintout(T,Srrho,S,Cp,H)
  use crest_parameters, only: wp => real64
  implicit none
  real(wp) :: T
  real(wp) :: Srrho
  real(wp) :: S
  real(wp) :: Cp
  real(wp) :: H
  character(len=80) :: atmp,btmp
  character(len=*),parameter :: sunit = ' (cal mol⁻¹ K⁻¹)'
  character(len=*),parameter :: gunit = ' (kcal mol⁻¹)'

  write (*,*)
  write (atmp,'(f12.2)') T
  write (btmp,'(a,1x,a,1x,a)') '|'//repeat(' ',6)//'FINAL MOLECULAR ENTROPY AT T=', &
  & trim(adjustl(atmp)),'K'//repeat(' ',6)//'|'
  call smallhead(trim(btmp))

  write (*,'(3x,''  Sconf   =      '',f12.6)') S-Srrho
  write (*,'(3x,''+ δSrrho  =      '',f12.6)') Srrho
  write (*,'(1x,53("-"))')
  write (*,'(1x,''= S(total)  =      '',f12.6,1x,a)') S, sunit
  write (*,*)
  write (*,'(3x,''H(T)-H(0) =      '',f12.6)') h
  write (*,'(3x,''G         =      '',f12.6,a)') (-S*T/1000.0_wp),'   (-T*S)'
  write (*,'(1x,53("-"))')
  write (*,'(1x,''= G(total)  =      '',f12.6,a,a)') (-S*T/1000.0_wp+H),'  (H-T*S)',gunit
  write (*,*)
  write (*,'(3x,''Cp(total) =      '',f12.6,1x,a)') cp, sunit
  write (*,*)

  return
end subroutine entropyprintout

!=========================================================================================!
subroutine redo_extrapol(fname,rtin)
!*****************************************************************
!*  Helper routine to re-do the extrapolation for
!*  a read-in file. The two corresponding file formats
!*  written by crest are .data and .entropydata (see code above)
!*  Use by:
!*  crest --redoextrapol .data
!*****************************************************************
  use crest_parameters
  use crest_data
  implicit none
  character(len=*) :: fname
  integer :: ich
  integer :: i,j,k
  integer :: iter,nt,rt,rtin
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

  rt = rtin
  call smallhead('entropy extrapolation from file '//trim(fname))
  write (*,*)

  if (trim(fname) == '.entropydata') then
    open (newunit=ich,file=fname)
    read (ich,*) iter,nt
    allocate (temps(nt),slist(iter,nt),source=0.0_wp)
    allocate (nalls(iter),source=0)
    do i = 1,iter
      read (ich,*) nalls(i)
      do j = 1,nt
        read (ich,*) slist(i,j),dum1,dum2
      end do
    end do
    close (ich)
    if (rtin == 0) rt = 3 !the third entry for each block is 298K in the crest default

  else if (trim(fname) == '.data') then
    open (newunit=ich,file=fname)
    iter = 0
    do
      read (ich,'(a)',iostat=io) atmp
      if (io < 0) exit
      iter = iter+1
    end do
    close (ich)
    nt = 1
    rt = 1  !.data only has the data for 298K
    allocate (temps(nt),slist(iter,nt),source=0.0_wp)
    allocate (nalls(iter),source=0)
    open (newunit=ich,file=fname)
    do i = 1,iter
      read (ich,*) nalls(i),slist(i,nt),dum1
    end do
    close (ich)

  else
    write (*,*) 'unknown file. must be either ".data" or ".entropydata"'
    error stop
  end if

  allocate (cpfinal(nt),sfinal(nt),hfinal(nt),srrho(nt),source=0.0d0)
  TLOOP: do k = 1,nt
    !      T = temps(k)   !for the read in version we dont need the temperatures
    prt = .false.
    if (k .eq. rt) prt = .true.  !printout only for RT

    if (iter .gt. 3) then ! at least 3 points (sMTDs)
      do itry = iter,5,-1   ! try with decreasing number of points untily fit is valid
        p3 = 1.0d0
        pmin = 1.0d0
        rmin = 1.0d+42
        okmin = .false.
        do ri = 1,51  ! constrained search for best p3
          call marqfit(.false.,.false.,iter,itry,slist(1:iter,k),p3,val,rmsd,ok) ! (last) 10 points max
          if (rmsd .lt. rmin.and.ok) then
            pmin = p3
            rmin = rmsd
            okmin = .true.
          end if
          p3 = p3-0.01d0
        end do
        if (okmin) exit
      end do
      call marqfit(prt,prt,iter,itry,slist(1:iter,k),pmin,val,rmsd,ok) ! final with best p3
      if (.not.ok) then
        val = maxval(slist(1:iter,k))
        if (prt) then
          write (*,'(''WARNING: large fit error expected'',/)')
          write (*,'(''FINAL ENTROPY (max conf.)  : '',f12.5)') val
        end if
      else
        if (prt) write (*,'(''FINAL ENTROPY (extrapol.)  : '',f12.5)') val
      end if
    else
      val = maxval(slist(1:iter,k))
      if (prt) then
        write (*,'(''WARNING: Not enough blocks for extrapolation.'')')
        write (*,'(''         Taking full ensemble entropy instead.'',/)')
        write (*,'(''FINAL ENTROPY (max conf.)  : '',f12.5)') val
      end if
    end if
    sfinal(k) = val
  end do TLOOP

  return
end subroutine redo_extrapol

