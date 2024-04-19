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

!=========================================================================================!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!=========================================================================================!

subroutine prepthermo(nat,at,xyz,pr,molmass,rabc,avmom,symnum,symchar)
!***********************************************************************
!* Prepare the calculation of thermodynamic properties of a structure
!* In particular, determine rotational constants and check the symmetry
!***********************************************************************
  use crest_parameters,only:wp,bohr,stdout
  use atmasses,only:molweight
  use iomod,only:to_lower
  use axis_module
  implicit none
  integer,intent(in)     :: nat
  integer,intent(in)     :: at(nat)
  real(wp),intent(inout) :: xyz(3,nat) !> in Angstroem
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

  !>--- molecular mass in amu
  molmass = molweight(nat,at)

  if (pr) then
    write (stdout,'(1x,a,f15.2)') 'Mol. weight /amu  : ',molmass
  end if

  !>--- rotational constants in cm-1
  rabc = 0.0d0
  call axis(nat,at,xyz,rabc(1:3),avmom)
  a = rabc(3)
  b = rabc(2)
  c = rabc(1)
  rabc(1) = a
  rabc(3) = c
  if (pr) then
    write (stdout,'(1x,a,3f15.2)') 'Rot. const. /MHz  : ',rabc(1:3)
  end if
  rabc = rabc/2.99792458d+4   ! MHz to cm-1
  if (pr) then
    write (stdout,'(1x,a,3f15.2)') 'Rot. const. /cm-1 : ',rabc(1:3)
  end if

  !>--- symmetry number from rotational symmetry
  xyz = xyz/bohr
  call getsymmetry2(.false.,6,nat,at,xyz,desy,maxat,sfsym)
  xyz = xyz*bohr
  sym = sfsym(1:3)
  symchar = sym
  symnum = 1.0d0
  if (a .lt. 1.d-9.or.b .lt. 1.d-9.or.c .lt. 1.d-9) then
    if (index(sym,'d') .ne. 0) symnum = 2.0d0
  else
    call to_lower(sym)
    if (index(sym,'c2') .ne. 0) symnum = 2.0d0
    if (index(sym,'s4') .ne. 0) symnum = 2.0d0
    if (index(sym,'c3') .ne. 0) symnum = 3.0d0
    if (index(sym,'s6') .ne. 0) symnum = 3.0d0
    if (index(sym,'c4') .ne. 0) symnum = 4.0d0
    if (index(sym,'s8') .ne. 0) symnum = 4.0d0
    if (index(sym,'c5') .ne. 0) symnum = 5.0d0
    if (index(sym,'c6') .ne. 0) symnum = 6.0d0
    if (index(sym,'c7') .ne. 0) symnum = 7.0d0
    if (index(sym,'c8') .ne. 0) symnum = 8.0d0
    if (index(sym,'c9') .ne. 0) symnum = 9.0d0
    if (index(sym,'d2') .ne. 0) symnum = 4.0d0
    if (index(sym,'d3') .ne. 0) symnum = 6.0d0
    if (index(sym,'d4') .ne. 0) symnum = 8.0d0
    if (index(sym,'d5') .ne. 0) symnum = 10.0d0
    if (index(sym,'d6') .ne. 0) symnum = 12.0d0
    if (index(sym,'d7') .ne. 0) symnum = 14.0d0
    if (index(sym,'d8') .ne. 0) symnum = 16.0d0
    if (index(sym,'d9') .ne. 0) symnum = 18.0d0
    if (index(sym,'t') .ne. 0) symnum = 12.0d0
    if (index(sym,'td') .ne. 0) symnum = 12.0d0
    if (index(sym,'th') .ne. 0) symnum = 12.0d0
    if (index(sym,'o') .ne. 0) symnum = 24.0d0
    if (index(sym,'oh') .ne. 0) symnum = 24.0d0
    if (index(sym,'ih') .ne. 0) symnum = 60.0d0
  end if

  if (pr) then
    write (stdout,'(1x,a,4x,a)') 'Symmetry:',sym
  end if
  return
end subroutine prepthermo

!=========================================================================================!
subroutine calcthermo(nat,at,xyz,freq,pr,ithr,fscal,sthr,nt,temps, &
    &      et,ht,gt,stot)
!**************************************************************
!* Calculate thermodynamic contributions for a given structure
!* from it's frequencies (from second derivatives/the Hessian)
!* Based on xtb's "print_thermo" routine
!**************************************************************
  use crest_parameters,only:wp,bohr,stdout
  use crest_thermo
  use atmasses,only:molweight
  use iomod,only:to_lower
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
  logical :: linear = .false.
  logical :: atom = .false.
  integer :: nvib_theo
  integer :: nvib,nimag
  real(wp) :: vibthr
  real(wp),allocatable :: vibs(:)

  integer :: i,j
  integer :: n3,rt
  real(wp) :: adum(nt)
  character(len=64) :: atmp

  character(len=*),parameter :: outfmt = &
  &  '(9x,"::",1x,a,f24.12,1x,a,1x,"::")'
  character(len=*),parameter :: dblfmt = &
  &  '(10x,":",2x,a,f24.7,1x,a,1x,":")'
  character(len=*),parameter :: intfmt = &
  &  '(10x,":",2x,a,i24,       6x,":")'
  character(len=*),parameter :: chrfmt = &
  &  '(10x,":",2x,a,a24,       6x,":")'

  real(wp),parameter :: autorcm = 219474.63067_wp
  real(wp),parameter :: rcmtoau = 1.0_wp/autorcm
  real(wp),parameter :: autocal = 627.50947428_wp*1000.0_wp

  call prepthermo(nat,at,xyz,pr,molmass,rabc,avmom,sym,symchar)

  n3 = 3*nat
  allocate (vibs(n3))
  vibthr = 1.0
  a = rabc(1)
  b = rabc(2)
  c = rabc(3)

  nvib_theo = 3*nat-6
  if (c .lt. 1.d-10 .or. (symchar=='din')) linear = .true.
  if (linear) nvib_theo = 3*nat-5

  if (a+b+c .lt. 1.d-6) then
    atom = .true.
    nvib = 0
    nvib_theo = 0
  end if

  nvib = 0
  vibs = 0.0
  do i = 1,n3
    if (abs(freq(i)) .gt. vibthr) then
      nvib = nvib+1
      vibs(nvib) = freq(i)
    end if
  end do
  !> scale
  vibs(1:nvib) = vibs(1:nvib)*fscal

  !> invert imaginary modes
  nimag = 0
  do i = 1,nvib
    if (vibs(i) .lt. 0.and.vibs(i) .gt. ithr) then
      vibs(i) = -vibs(i)
      if (pr) write (stdout,'(a,i5," :",f10.2)') 'Inverting frequency',i,vibs(i)
    end if
    if (vibs(i) < 0.0) then
      nimag = nimag+1
    end if
  end do

  if (pr) then
    write (stdout,'(a)')
    write (stdout,'(10x,51("."))')
    write (stdout,'(10x,":",22x,a,22x,":")') "SETUP"
    write (stdout,'(10x,":",49("."),":")')
    write (stdout,intfmt) "# frequencies    ",nvib
    write (stdout,intfmt) "# imaginary freq.",nimag
    write (atmp,*) linear
    write (stdout,chrfmt) "linear?          ",trim(atmp)
    write (stdout,chrfmt) "symmetry         ",adjustr(symchar)
    write (stdout,intfmt) "rotational number",nint(sym)
    write (stdout,dblfmt) "scaling factor   ",fscal,"    "
    write (stdout,dblfmt) "rotor cutoff     ",sthr,"cm⁻¹"
    write (stdout,dblfmt) "imag. cutoff     ",ithr,"cm⁻¹"
    write (stdout,'(10x,":",49("."),":")')
  end if

  vibs = vibs*rcmtoau   ! thermodyn needs vibs and zp in Eh

  zp = 0.5_wp*sum(vibs(1:nvib))
  adum = abs(temps-298.15d0)
  rt = minloc(adum,1)  !temperature closest to 298.15 is the ref.
  do j = 1,nt
    if ((j == rt).and.pr) then
      pr2 = .true.
    else
      pr2 = .false.
    end if
    if (pr2) then
      call print_thermo_sthr_ts(stdout,nvib,vibs,avmom,sthr,temps(j))
    end if
    call thermodyn(stdout,a,b,c,avmom,linear,atom,sym,molmass,vibs,nvib, &
    & temps(j),sthr,et(j),ht(j),gt(j),ts(j),zp,pr2)
    stot(j) = (ts(j)/temps(j))*autocal
  end do

  if ((nt > 1).and.pr) then
    write (stdout,'(a)')
    write (stdout,'(a10)',advance='no') "T/K"
    write (stdout,'(a16)',advance='no') "H(0)-H(T)+PV"
    write (stdout,'(a16)',advance='no') "H(T)/Eh"
    write (stdout,'(a16)',advance='no') "T*S/Eh"
    write (stdout,'(a16)',advance='no') "G(T)/Eh"
    write (stdout,'(a)')
    write (stdout,'(3x,72("-"))')
    do i = 1,nt
      write (stdout,'(3f10.2)',advance='no') temps(i)
      write (stdout,'(3e16.6)',advance='no') ht(i)
      write (stdout,'(3e16.6)',advance='no') et(i)
      write (stdout,'(3e16.6)',advance='no') ts(i)
      write (stdout,'(3e16.6)',advance='no') gt(i)
      if (i == rt) then
        write (stdout,'(1x,"(used)")')
      else
        write (stdout,'(a)')
      end if
    end do
    write (stdout,'(3x,72("-"))')
  end if

  deallocate (vibs)
  return
end subroutine calcthermo

!=========================================================================================!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!=========================================================================================!

subroutine thermo_wrap_legacy(env,pr,nat,at,xyz,dirname, &
        &  nt,temps,et,ht,gt,stot,bhess)
!**********************************************
!* Wrapper for a Hessian calculation to get
!* the thermodynamics of the molecule.
!* Legacy version that uses xtb and reads
!* the frequencies from a vibspectrum file
!*********************************************
  use crest_parameters,only:wp,stdout
  use crest_data
  use iomod
  use strucrd
  implicit none
  !> INPUT
  type(systemdata) :: env
  logical,intent(in) :: pr
  integer,intent(in) :: nat
  integer,intent(inout) :: at(nat)
  real(wp),intent(inout) :: xyz(3,nat)  !> in Angstroem!
  character(len=*) :: dirname
  integer,intent(in)  :: nt
  real(wp),intent(in)  :: temps(nt)
  logical,intent(in) :: bhess       !> calculate bhess instead?
  !> OUTPUT
  real(wp),intent(out) :: et(nt)    !> enthalpy in Eh
  real(wp),intent(out) :: ht(nt)    !> enthalpy in Eh
  real(wp),intent(out) :: gt(nt)    !> free energy in Eh
  real(wp),intent(out) :: stot(nt)  !> entropy in cal/molK
  !> LOCAL
  logical :: subdir,ex
  integer :: i,io,r,ich
  character(len=1024) :: jobcall
  character(len=*),parameter :: pipe = '2>/dev/null'
  character(len=*),parameter :: xname = 'freq.xyz'
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

  optpath = ''

  subdir = .false.
  i = len_trim(dirname)
  if (i > 0) subdir = .true.

  !>-- build the jobcall
  jobcall = ""
  jobcall = trim(jobcall)//trim(env%ProgName)
  if (bhess) then
    jobcall = trim(jobcall)//" "//trim(xname)//' --bhess loose'
  else
    jobcall = trim(jobcall)//" "//trim(xname)//' --ohess'
  end if
  jobcall = trim(jobcall)//" "//trim(env%gfnver)
  jobcall = trim(jobcall)//" "//trim(env%solv)
  if (env%chrg /= 0) then
    jobcall = trim(jobcall)//" --chrg "//to_str(env%chrg)
  end if
  if (env%uhf /= 0) then
    jobcall = trim(jobcall)//" --uhf "//to_str(env%uhf)
  end if
  jobcall = trim(jobcall)//' --ceasefiles > xtb.out '//trim(pipe)

  if (subdir) then
    call rmrf(trim(dirname))
    r = makedir(trim(dirname))
    optpath = trim(dirname)//'/'
  end if

  call env%wrtCHRG(trim(optpath))
  inquire (file='gfnff_topo',exist=ex)
  if (env%gfnver == '--gff'.and.subdir.and.ex) then
    call getcwd(thispath)
    io = sylnk(trim(thispath)//'/'//'gfnff_topo',trim(optpath)//'gfnff_topo')
  end if
  if (index(env%fixfile,'none selected') .eq. 0) then
    io = sylnk(trim(thispath)//'/'//env%fixfile,trim(optpath)//env%fixfile)
  end if

!$omp critical
  open (unit=ich,file=trim(optpath)//xname)
  call wrxyz(ich,nat,at,xyz)
  if (env%thermo%constrhess) then
    call write_cts(ich,env%cts)
  end if
  close (ich)
!$omp end critical

  if (subdir) then
    jobcall2 = 'cd '//trim(dirname)//' && '//trim(jobcall)
    call command( jobcall2, io )
  else
    call command( jobcall, io )
  end if

  et = 0.0_wp
  ht = 0.0_wp
  gt = 0.0_wp
  stot = 0.0_wp

  if (io /= 0) then  !if the calc failed
    return
  end if

!$omp critical
  call rdxmol(trim(optpath)//'xtbopt.xyz',nat,at,xyz,atmp)
  etot = grepenergy(atmp)
  nfreq = 3*nat

  allocate (freq(nfreq))
  call rdfreq(trim(optpath)//'vibspectrum',nfreq,freq)

  ithr = env%thermo%ithr
  fscal = env%thermo%fscal
  sthr = env%thermo%sthr
  call calcthermo(nat,at,xyz,freq,pr,ithr,fscal,sthr, &
  &    nt,temps,et,ht,gt,stot)
  deallocate (freq)
!$omp end critical
  call initsignal()
  return
end subroutine thermo_wrap_legacy

!=========================================================================================!
subroutine rdfreq(fname,nmodes,freq)
!**************************************
!* read vibspectrum file in TM format
!**************************************
  use crest_parameters,only:wp
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

  freq = 0.0_wp
  inquire (file=fname,exist=ex)
  if (.not.ex) return
  k = 1 !modes
  open (file=fname,unit=ich)
  rdfile: do
    read (ich,'(a)',iostat=io) atmp
    if (io < 0) exit
    if (index(atmp,'$vibrational spectrum') .ne. 0) then
      rdblock: do
        read (ich,'(a)',iostat=io) atmp
        if (io < 0) exit rdfile
        if (index(atmp,'$end') .ne. 0) exit rdfile
        if (index(atmp,'#') .ne. 0) cycle rdblock !skip comment lines
        call readl(atmp,floats,n)
        freq(k) = floats(2)
        k = k+1
      end do rdblock
    end if
  end do rdfile
  close (ich)
  return
end subroutine rdfreq

!=========================================================================================!

subroutine thermo_wrap_new(env,pr,nat,at,xyz,dirname, &
        &  nt,temps,et,ht,gt,stot,bhess)
!*********************************************
!* Wrapper for a Hessian calculation to get
!* the thermodynamics of the molecule.
!* Updated version without xtb subprocess
!*********************************************
!*** WARNING: xyz is expected in ANGSTROEM ***
!*********************************************
  use crest_parameters,only:wp,stdout,aatoau
  use crest_data
  use crest_calculator
  use iomod
  use strucrd
  use hessian_tools
  implicit none
  !> INPUT
  type(systemdata) :: env
  logical,intent(in) :: pr
  integer,intent(in) :: nat
  integer,intent(inout) :: at(nat)
  real(wp),intent(inout) :: xyz(3,nat)  !> in Angstroem!
  character(len=*) :: dirname
  integer,intent(in)  :: nt
  real(wp),intent(in)  :: temps(nt)
  logical,intent(in) :: bhess       !> calculate bhess instead?
  !> OUTPUT
  real(wp),intent(out) :: et(nt)    !> enthalpy in Eh
  real(wp),intent(out) :: ht(nt)    !> enthalpy in Eh
  real(wp),intent(out) :: gt(nt)    !> free energy in Eh
  real(wp),intent(out) :: stot(nt)  !> entropy in cal/molK
  !> LOCAL
  type(coord) :: mol
  type(calcdata) :: calctmp
  character(len=10) :: atmp
  logical :: subdir,ex
  integer :: i,io,r,ich
  real(wp) :: etot
  integer :: nfreq
  real(wp),allocatable :: hess(:,:)
  real(wp),allocatable :: freq(:)
  real(wp) :: ithr,fscal,sthr

  integer :: TID,OMP_GET_THREAD_NUM

!!$OMP PARALLEL PRIVATE(TID)
  TID = OMP_GET_THREAD_NUM()
  !awrite(*,*) '---->',TID
!!$OMP END PARALLEL
  ich = (TID+1)*1000   ! generate CPU dependent file channel number

  call initsignal()

  subdir = .false.
  if (len_trim(dirname) > 0) subdir = .true.

!>-- create a calculation object locally, modify calc dir
!$omp critical
  calctmp = env%calc
  calctmp%pr_energies = .false. !> never do that!
  mol%nat = nat
  mol%at = at
  mol%xyz = xyz*aatoau

  do i = 1,calctmp%ncalculations
    write (atmp,'(".",i0)') i
    if (subdir) then
      calctmp%calcs(i)%calcspace = trim(dirname)
    else if (allocated(calctmp%calcs(i)%calcspace)) then
      deallocate (calctmp%calcs(i)%calcspace)
    end if
  end do
!>-- also, allocate frequncy and hessian space
  nfreq = 3*nat
  allocate (freq(nfreq),source=0.0_wp)
  allocate (hess(nfreq,nfreq),source=0.0_wp)
!$omp end critical

!>-- numerical Hessian

  !TODO bhess currently not coded with new calculator
  if (bhess) then
    write (stdout,'("> ",a)') 'bhess not implemented for calculator routines'
  end if
  !else
  call numhess1(mol%nat,mol%at,mol%xyz,calctmp,hess,io)
  !end if

  if (io /= 0) then  !if the calc failed
    return
  end if

!>-- project and get frequencies
!$omp critical
  !>-- Projects and mass-weights the Hessian
  call prj_mw_hess(mol%nat,mol%at,nfreq,mol%xyz,hess)
  !>-- Computes the Frequencies
  call frequencies(mol%nat,mol%at,mol%xyz,nfreq,calctmp,hess,freq,io)
!$omp end critical

!>--- get thermodynamics
!$omp critical
  et = 0.0_wp
  ht = 0.0_wp
  gt = 0.0_wp
  stot = 0.0_wp

  ithr = env%thermo%ithr
  fscal = env%thermo%fscal
  sthr = env%thermo%sthr
  call calcthermo(nat,at,xyz,freq,pr,ithr,fscal,sthr, &
  &    nt,temps,et,ht,gt,stot)
  deallocate (hess,freq)
!$omp end critical
  call initsignal()
  return
end subroutine thermo_wrap_new

!=========================================================================================!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!=========================================================================================!
subroutine calcSrrhoav(env,ensname)
!*******************************************************
!* Calculate S_RRHO averages for a given ensemlbe
!*******************************************************
  use crest_parameters
  use crest_data
  use strucrd
  use iomod
  implicit none
  !> INPUT
  type(systemdata) :: env
  character(len=*) :: ensname
  !> LOCAL
  real(wp),allocatable :: cp(:)
  real(wp),allocatable :: hconf(:)
  integer :: nat,nall
  integer,allocatable :: at(:)
  real(wp),allocatable :: xyz(:,:,:)
  real(wp),allocatable :: er(:)
  real(wp),allocatable :: erel(:)
  real(wp),allocatable :: efree(:,:)
  real(wp),allocatable :: g(:) !> degeneracies, either read from cre_degen2 (if present), or set to unity
  real(wp),allocatable :: p(:,:) !> populations at different T
  integer,allocatable :: pindex(:)
  real(wp),allocatable :: gatt(:,:)
  real(wp),allocatable :: satt(:,:)
  real(wp),allocatable :: srrho(:),sav(:)
  real(wp),allocatable :: bsatt(:)
  real(wp),allocatable :: gav(:)
  real(wp),allocatable :: pdum(:)
  real(wp) :: psum,emin,sdum
  real(wp) :: quick_rmsd,rmsdval
  integer :: eloc,ploc
  logical :: avbhess

  integer :: nt
  real(wp),allocatable :: temps(:)
  real(wp),allocatable :: et(:)
  real(wp),allocatable :: ht(:)
  real(wp),allocatable :: gt(:)
  real(wp),allocatable :: stot(:)
  real(wp),allocatable :: c0(:,:)
  real(wp),allocatable :: sref(:)
  character(len=64) :: atmp
  integer :: i,j,k,ich,io,popf,ii,T,Tn
  logical :: ex
  logical :: niceprint
  real(wp) :: percent
  character(len=52) :: bar
  integer :: ncalc,vz,nlimit,nav
  character(len=512) :: thispath,tmppath

  real(wp),parameter :: Tref = 298.15  !> room temperature is reference
  real(wp),parameter :: kcal = autokcal

!>--- read the given ensemble
  call rdensembleparam(trim(ensname),nat,nall)
  allocate (at(nat),xyz(3,nat,nall),er(nall))
  call rdensemble(trim(ensname),nat,nall,at,xyz,er)

  if (any(er(:) > 0.0d0)) then
    error stop 'ensemble file must contain energies in Eh! must stop'
  end if

  write (tmppath,'(a)') 'Frequency Calculation and Averages'
  write (stdout,*)
  call smallhead(trim(tmppath))

!>--- temperatures from sys object
  if (.not.allocated(env%thermo%temps)) then
    call env%thermo%get_temps()
  end if
  nt = env%thermo%ntemps
  allocate (temps(nt))
  temps = env%thermo%temps

!>--- space for populations and degeneracies
  allocate (g(nall),source=1.0_wp)
  allocate (p(nall,nt))

!>--- read degeneracies?
  inquire (file='cre_degen2',exist=ex)
  if (ex) then
    open (newunit=ich,file='cre_degen2')
    read (ich,*) atmp
    do i = 1,nall
      read (ich,*,iostat=io) j,g(i)
      if (io < 0) exit
    end do
    close (ich)
  else
    g = 1
  end if

!========================================================================================!
!> FREQUENCY CALCULATION AND THERMODYNAMICS
!========================================================================================!
!>--- determine how many hessians must be calculated
!>--- NOTE: this assumes the ensemble is ordered by energy, lowest first.
  allocate (erel(nall),pdum(nall),pindex(nall))
  emin = minval(er(:),1)
  erel = (er-emin)*kcal
  call entropy_boltz(nall,Tref,erel,g,pdum)

!>--- set up index for ensemble that are NOT energy-sorted
  do i = 1,nall
    pindex(i) = i
  end do
  pdum(:) = -pdum(:)             !> Hack because qsort does low-to-high
  call qsort(pdum,1,nall,pindex) !> pindex is what we are after
  pdum(:) = -pdum(:)             !> and switch sign back

!>--- and with the sorted pdum, just count how many calculations we need
  ncalc = 1  !> always take the lowest
  ploc = maxloc(pdum(:),1)
  psum = pdum(ploc)
  nlimit = env%thermo%pcap !> limit strucs (for VERY large SE)
  do i = 2,nall
    psum = psum+pdum(i)
    ncalc = ncalc+1
    if (ncalc == nlimit) then
      exit
    end if
    if (psum > env%thermo%ptot) then
      exit
    end if
  end do
  deallocate (pdum)

!>--- print something
  write (stdout,'(1x,a,i0)') 'Nconf on file      : ',nall
  write (atmp,'(1x,a,f6.2,a)') '(=',psum*100.0d0,'% total population)'
  write (stdout,'(1x,a,i0,a)') 'Taken for Hessians : ',ncalc,trim(atmp)
  if (psum < env%thermo%ptot) then
    write (stdout,'(2x,a,i0,a)') '=> (Limited to ',ncalc,' structures due to amount of calcs.)'
  end if
  write (stdout,'(1x,a,f8.2,1x,f8.2)') "T range  /K    : ",temps(1),temps(nt)
  write (stdout,'(1x,a,f17.6,1x,a)') "scaling factor : ",env%thermo%fscal,"    "
  write (stdout,'(1x,a,f17.6,1x,a)') "rotor cutoff   : ",env%thermo%sthr,"cm⁻¹"
  write (stdout,'(1x,a,f17.6,1x,a)') "imag. cutoff   : ",env%thermo%ithr,"cm⁻¹"
  write (stdout,*)

!>--- calculate Hessians for ncalc lowest structures
  allocate (gatt(nall,nt),satt(nall,nt),source=0.0_wp)

  io = makedir('HESSIANS')
  call getcwd(thispath)
  inquire (file='gfnff_topo',exist=ex)

  call chdir('HESSIANS')
  if (env%legacy) then
    if (env%gfnver == '--gff'.and.ex) then
      call getcwd(tmppath)
      io = sylnk(trim(thispath)//'/'//'gfnff_topo',trim(tmppath)//'/'//'gfnff_topo')
    end if
    if (index(env%fixfile,'none selected') .eq. 0) then
      io = sylnk(trim(thispath)//'/'//env%fixfile,trim(tmppath)//'/'//env%fixfile)
    end if
  end if

  k = 0
  niceprint = env%niceprint

!>--- OMP stuff
  call new_ompautoset(env,'auto',ncalc,T,Tn)

!>--- the parallel loop
  avbhess = env%thermo%avbhess
  write (stdout,'(1x,a,i0,a)') 'Running ',ncalc,' calculations ...'
  call crest_oloop_pr_progress(env,ncalc,0)
!$omp parallel &
!$omp shared( vz,tmppath,ncalc,percent,k,bar,niceprint) &
!$omp shared( env,nat,at,xyz,c0,et,ht,gt,stot,temps,nt,gatt,satt,avbhess,pindex )
!$omp single
  allocate (et(nt),ht(nt),gt(nt),stot(nt))
  allocate (c0(3,nat))
  do i = 1,ncalc
    call initsignal()
    vz = pindex(i) !> restore index
    !$omp task firstprivate( vz ) private( tmppath,et,ht,gt,stot,c0 )
    call initsignal()
    !$omp critical
    write (tmppath,'(''hess'',i0)') vz
    c0(1:3,1:nat) = xyz(1:3,1:nat,vz)
    !$omp end critical

    call thermo_wrap(env,.false.,nat,at,c0,tmppath, &
    &    nt,temps,et,ht,gt,stot,avbhess)

    !$omp critical
    gatt(vz,1:nt) = gt(1:nt)
    satt(vz,1:nt) = stot(1:nt)
    !$omp end critical
    if (.not.env%keepModef) call rmrf(trim(tmppath))

    !$omp critical
    k = k+1
    call crest_oloop_pr_progress(env,ncalc,k)
    !$omp end critical
    !$omp end task
  end do
  deallocate (c0)
  deallocate (stot,gt,ht,et)
!$omp taskwait
!$omp end single
!$omp end parallel
  call crest_oloop_pr_progress(env,ncalc,-1)
  call chdir(thispath)
  if (.not.env%keepModef) call rmrf('HESSIANS')

!========================================================================================!
!>--- process the calculated free energies and entropies into accurate populations
  allocate (srrho(nt),sav(nt),gav(nt),efree(nall,nt))
  srrho = 0.0_wp
  sav = 0.0_wp
  gav = 0.0_wp
  nav = ncalc
  write (stdout,'(1x,a)',advance='no') 'calculating averages for G and S ... '
  flush (stdout)
  do j = 1,nt
    do ii = 1,ncalc
      i = pindex(ii) !> restore index
      if (abs(gatt(i,j)) .lt. 1.d-10) then  !> failed calcs?
        if (j == 1) nav = nav-1
      end if
      gav(j) = gav(j)+gatt(i,j)
      sav(j) = sav(j)+satt(i,j)
    end do
  end do
  gav = gav/float(nav)   !> get the average G(T)
  sav = sav/float(nav)   !> get the avverage S(T)
  write (stdout,'(a8)') 'done.'

!>--- get the free energies
  do j = 1,nt
    efree(:,j) = er(:) !> all based on etot
    do ii = 1,nall
      i = pindex(ii) !> restore index
      if (ii <= ncalc) then
        if (abs(gatt(i,j)) < 1.d-10) then
          efree(i,j) = efree(i,j)+gav(j) !> add |G(T)| (for failed calcs)
        else
          efree(i,j) = efree(i,j)+gatt(i,j) !> add G(T)
        end if
      else
!>-- for all energies that were not included in the free energy calculation add the average
        efree(i,j) = efree(i,j)+gav(j)
      end if
    end do
  end do

!>--- make relative energies and calculate Boltzman populations
  write (stdout,'(1x,a)',advance='no') 'calculating Boltzmann weights ... '
  flush (stdout)
  allocate (pdum(nall))
  do j = 1,nt
    emin = minval(efree(:,j),1)    !> lowest as reference
    erel = (efree(:,j)-emin)*kcal  !> to relative energies in kcal/mol
    pdum = 0.0d0
    call entropy_boltz(nall,temps(j),erel,g,pdum)
    !call entropy_boltz(ncalc,temps(j),erel,g(1:ncalc),pdum(1:ncalc))
    p(:,j) = pdum(:)
  end do
  deallocate (pdum)
  write (stdout,'(a11)') 'done.'
  if (env%thermo%printpop) then
    popf = makedir('populations')
    do j = 1,nt
      write (tmppath,'(a,a,a,i0)') 'populations','/','.pop_',nint(temps(j))
      open (newunit=popf,file=trim(tmppath))
      do k = 1,nall
        write (popf,'(f16.8)') p(k,j)
      end do
      close (popf)
    end do
  end if

!=========================================================================================!
!==== after this point p now contains the correct populations based on free energies =====!
!=========================================================================================!
!>--- S_avRRHO must be calculated relative to the actual DFT reference structure
!>--- the corresponding frequencies can be calculated with bhess
  if (env%emtd%bhess.and.allocated(env%emtd%xyz)) then
    allocate (bsatt(nt))
    allocate (et(nt),ht(nt),gt(nt))
    write (stdout,'(1x,a)',advance='no') 'calculating reference S (bhess) ... '
    flush (stdout)
    call thermo_wrap(env,.false.,env%emtd%nat,env%emtd%at,  &
  &    env%emtd%xyz,'BHESS',nt,temps,et,ht,gt,bsatt,.true.)
    if (.not.env%keepModef) call rmrf('BHESS')
    deallocate (gt,ht,et)
    write (stdout,'(a9)') 'done.'
  end if
!>--- average S_rrho with CORRECT populations
  allocate (sref(nt),source=0.0_wp)
  write (stdout,'(1x,a)',advance='no') 'calculating δSrrho ... '
  flush (stdout)
  srrho = 0.0d0
  do j = 1,nt
    do ii = 1,nall
      i = pindex(ii) !> restore index
      if (ii <= ncalc) then
        if (abs(satt(i,j)) < 1.d-10) then
          srrho(j) = srrho(j)+p(i,j)*sav(j)    !> (for failed hess calcs)
        else
          srrho(j) = srrho(j)+p(i,j)*satt(i,j) !> corrected for different S_rrho
        end if
      else
        srrho(j) = srrho(j)+p(i,j)*sav(j)
      end if
    end do
!>--- substract the reference value to shift the average
    if (env%emtd%bhess.and.allocated(env%emtd%xyz)) then
      sref(j) = bsatt(j) !> if a bhess value is available
    else
      sref(j) = satt(pindex(1),j) !> lowest in ensemble otherwise
    end if
    srrho(j) = srrho(j)-sref(j)
  end do
  write (stdout,'(a22)') 'done.'

  if (env%emtd%bhess.and.allocated(env%emtd%xyz)) then
    write (stdout,*)
    call underline('Coordinates for the bhess reference structure (Ångström):')
    call wrxyz(stdout,env%emtd%nat,env%emtd%at,env%emtd%xyz)
    write (stdout,*) '-----------------------------------------------------------'
    write (stdout,'(1x,a,a,a)') 'as read from <',env%emtd%fromfile,'>'
    inquire (file='crest_best.xyz',exist=ex)
    if (ex) then
      rmsdval = quick_rmsd('crest_best.xyz',env%emtd%nat,env%emtd%at,env%emtd%xyz,.true.)
      write (stdout,'(1x,a)') 'Heavy-atom RMSD between lowest conformer and this reference :'
      write (stdout,'(1x,a,f16.6,a)') 'RMSD(heavy) =',rmsdval,' Å'
      write (stdout,*)
    end if
    write (stdout,'(1x,a)') 'msRRHO(bhess) reference entropies:'
    do i = 1,nt
      write (stdout,'(2x,f10.2,2x,f16.6)') temps(i),bsatt(i)
    end do
  end if

!>--- prinout for the average free energy and entropy
  if ((nt > 1)) then
    write (stdout,'(a)')
    write (stdout,'(a10)',advance='no') "T/K"
    write (stdout,'(a17)',advance='no') "|S(T)|/cal/molK"
    write (stdout,'(a16)',advance='no') "|G(T)|/Eh"
    write (stdout,'(a16)',advance='no') "G_lowest/Eh"
    write (stdout,'(a10)',advance='no') "(conf)"
    write (stdout,'(a)')
    write (stdout,'(3x,65("-"))')
    do i = 1,nt
      write (stdout,'(3f10.2)',advance='no') temps(i)
      write (stdout,'(3e16.6)',advance='no') srrho(i)+sref(i)
      write (stdout,'(3e16.6)',advance='no') gav(i)
      emin = minval(efree(:,i),1)
      write (stdout,'(f16.6)',advance='no') emin
      eloc = minloc(efree(:,i),1)
      write (stdout,'(i10)',advance='no') eloc
      write (stdout,'(a)')
    end do
    write (stdout,'(3x,65("-"))')
    write (stdout,'(3x,a,a)') 'NOTE: if |G(T)| is the averaged ', &
    & 'contributrion to the free energy.'
    write (stdout,'(3x,a,a,i0,a)') '|G(T)| used only for the higher-energetic ', &
    & 'structures (n > ',ncalc,').'
    write (stdout,'(3x,a,a)') 'All other structures use ', &
    & 'G(T) from the respective Hessian calculations.'
  end if

!>--- properties based on free energies
  allocate (cp(nt),hconf(nt))
  do j = 1,nt
    emin = minval(efree(:,j),1)       !> lowest as reference
    erel = (efree(:,j)-emin)*kcal     !> to relative energies in kcal/mol
    call entropy_S(nall,temps(j),1.0d0,erel, &
    &    g,sdum,cp(j),hconf(j))  !> g read from file or set to 1
  end do

  if ((nt > 1)) then
    write (stdout,*)
    write (stdout,'(1x,a)') 'Quantities calculated on free energies:'
    write (stdout,'(a10)',advance='no') "T/K"
    write (stdout,'(a17)',advance='no') "δSrrho"
    write (stdout,'(a16)',advance='no') "Cp(T)"
    write (stdout,'(a16)',advance='no') "[H(T)-H(0)]"
    write (stdout,'(a)')
    write (stdout,'(3x,55("-"))')
    do i = 1,nt
      write (stdout,'(3f10.2)',advance='no') temps(i)
      write (stdout,'(3e16.6)',advance='no') srrho(i)
      write (stdout,'(3e16.6)',advance='no') cp(i)
      write (stdout,'(f16.6)',advance='no') hconf(i)
      write (stdout,'(a)')
    end do
    write (stdout,'(3x,55("-"))')
    write (stdout,'(3x,a,a)') 'NOTE: δSrrho(T) = |S(T)| - Sref(T)'
  end if

  if (allocated(env%emtd%soft)) then
    env%emtd%soft(:) = srrho(:)  !> this is \overline{S}_{msRRHO}
  end if
  if (allocated(env%emtd%cpoft)) then
    env%emtd%cpoft(:) = cp(:)    !> this is Cp_conf
  end if
  if (allocated(env%emtd%hoft)) then
    env%emtd%hoft(:) = hconf(:)  !> this is H_conf
  end if

  if (allocated(sref)) deallocate (sref)
  if (allocated(bsatt)) deallocate (bsatt)
  if (allocated(pindex)) deallocate (pindex)
  deallocate (satt,gatt)
  deallocate (efree,gav,sav,srrho)
  deallocate (erel,p,g,temps)
  deallocate (er,xyz,at)
  return
end subroutine calcSrrhoav
