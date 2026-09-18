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
  use thermochem_module
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
  type(coord) :: mol
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
    call command(jobcall2,io)
  else
    call command(jobcall,io)
  end if

  et = 0.0_wp
  ht = 0.0_wp
  gt = 0.0_wp
  stot = 0.0_wp

  if (io /= 0) then  !if the calc failed
    return
  end if

!$omp critical
  !call rdxmol(trim(optpath)//'xtbopt.xyz',nat,at,xyz,atmp)
  call mol%open(trim(optpath)//'xtbopt.xyz')
  etot = grepenergy(atmp)
  nfreq = 3*mol%nat

  allocate (freq(nfreq))
  call rdfreq(mol,trim(optpath)//'vibspectrum',nfreq,freq)

  ithr = env%thermo%ithr
  fscal = env%thermo%fscal
  sthr = env%thermo%sthr
  call calcthermo(mol%nat,mol%at,mol%xyz,freq,pr,ithr,fscal,sthr, &
  &    nt,temps,et,ht,gt,stot,stdout,emodel=env%thermo%emodel)
  deallocate (freq)
!$omp end critical
  call initsignal()
  return
end subroutine thermo_wrap_legacy

!=========================================================================================!

subroutine rdfreq(mol,fname,nmodes,freq)
!**************************************
!* read vibspectrum file in TM format
!**************************************
  use crest_parameters,only:wp
  use crest_data
  use iomod
  use strucrd
  implicit none
  type(coord),intent(in) :: mol
  character(len=*),intent(in) :: fname
  integer,intent(in)   :: nmodes
  real(wp),intent(out) :: freq(nmodes)    !frequencies

  logical :: ex,ex2,ex3
  type(coord) :: moltmp
  freq(:) = 0.0_wp

  inquire (file=fname,exist=ex)
  if (.not.ex) return

  call minigrep(fname,'$vibrational spectrum',ex)
  if (ex) then
    !> TURBOMOLE "vibspectrum"-style file
    call rdfreq_vibspectrum_file(fname,nmodes,freq)
  end if
  call minigrep(fname,'$orca_hessian_file',ex)
  call minigrep(fname,'$hessian',ex2)
  call minigrep(fname,'$ir_spectrum',ex3)
  !if (ex.and.ex3) then
  !  !> ORCA ".hess" file --> frequencies directly
  !  call rdfreq_orca_ir_spectrum(fname,nmodes,freq)
  !else if (ex.and.ex2) then
  if(ex.and.ex2)then
    !> ORCA ".hess" file --> from Hessian
    moltmp = mol
    call rdfreq_orca_hess(moltmp,fname,nmodes,freq)
  end if

end subroutine rdfreq

subroutine rdfreq_vibspectrum_file(fname,nmodes,freq)
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
end subroutine rdfreq_vibspectrum_file

subroutine rdfreq_orca_ir_spectrum(fname,nmodes,freq)
!**************************************
!* read vibspectrum file in TM format
!**************************************
  use crest_parameters,only:wp,stdout
  use crest_data
  use iomod
  implicit none
  character(len=*),intent(in) :: fname
  integer,intent(in)   :: nmodes
  real(wp),intent(out) :: freq(nmodes)    !frequencies
  integer :: k,ich,io,n,nref
  character(len=256) :: atmp
  real(wp) :: floats(10)
  logical :: ex

  freq = 0.0_wp
  k = 1 !modes
  open (file=fname,unit=ich)
  rdfile: do
    read (ich,'(a)',iostat=io) atmp
    if (io < 0) exit
    if (index(atmp,'$ir_spectrum') .ne. 0) then
      read (ich,'(a)',iostat=io) atmp
      if (io < 0) exit rdfile
      read (atmp,*,iostat=io) nref
      if (io .ne. 0) exit rdfile
      if (nref .ne. nmodes) exit rdfile
      rdblock: do
        read (ich,'(a)',iostat=io) atmp
        if (io < 0) exit rdfile
        if (index(atmp,'$end') .ne. 0) exit rdfile
        if (index(atmp,'#') .ne. 0) cycle rdblock !skip comment lines
        call readl(atmp,floats,n)
        freq(k) = floats(1)
        if (k == nref) exit rdfile
        k = k+1
      end do rdblock
    end if
  end do rdfile
  if (k .ne. nmodes) then
    write (stdout,*) '** WARNING ** error while reading '//trim(fname)
  end if
  close (ich)
  return
end subroutine rdfreq_orca_ir_spectrum

subroutine rdfreq_orca_hess(mol,fname,nmodes,freq)
!**************************************
!* read vibspectrum file in TM format
!**************************************
  use crest_parameters,only:wp,stdout
  use crest_data
  use iomod
  use strucrd
  use thermochem_module
  implicit none
  type(coord),intent(inout) :: mol
  character(len=*),intent(in) :: fname
  integer,intent(in)   :: nmodes
  real(wp),intent(out) :: freq(nmodes)    !frequencies
  integer :: k,ich,io,n,nref
  integer :: ii,jj,kk,iblocks,jblocks,ll
  character(len=256) :: atmp
  real(wp) :: floats(10)
  logical :: ex
  real(wp),allocatable :: hess(:,:)

  freq = 0.0_wp
  allocate (hess(nmodes,nmodes),source=0.0_wp)
  k = 1 !modes
  open (file=fname,unit=ich)
  rdfile: do
    read (ich,'(a)',iostat=io) atmp
    if (io < 0) exit
    if (index(atmp,'$hessian') .ne. 0) then
      read (ich,'(a)',iostat=io) atmp
      if (io < 0) exit rdfile
      read (atmp,*,iostat=io) nref
      if (io .ne. 0) exit rdfile
      if (nref .ne. nmodes) exit rdfile
      iblocks = (floor(real(nref,wp)/5.0_wp))
      jblocks = nref-(iblocks*5)
      rdblock1: do ii = 1,iblocks
        do jj = 0,nref
          read (ich,'(a)',iostat=io) atmp
          if (io < 0) exit rdfile
          if (index(atmp,'$end') .ne. 0) exit rdfile
          if (index(atmp,'#') .ne. 0) cycle rdblock1 !skip comment lines
          call readl(atmp,floats,n)
          if (jj > 0) then
            kk = (ii-1)*5
            do ll = 1,5
              hess(kk+ll,jj) = floats(1+ll)
            end do
          end if
        end do
      end do rdblock1
      if (jblocks > 0) then
        do jj = 0,nref
          read (ich,'(a)',iostat=io) atmp
          if (io < 0) exit rdfile
          if (index(atmp,'$end') .ne. 0) exit rdfile
          if (index(atmp,'#') .ne. 0) cycle
          call readl(atmp,floats,n)
          if (jj > 0) then
            kk = (ii-1)*5
            do ll = 1,jblocks
              hess(kk+ll,jj) = floats(1+ll)
            end do
          end if
        end do
      end if
    end if
  end do rdfile
  if (nref .ne. nmodes) then
    write (stdout,*) '** WARNING ** error while reading '//trim(fname)
  end if
  close (ich)

  write(stdout,'(a)',advance='no') ' Processing (raw) Hessian read from ORCA '//trim(fname)//' ... '
  flush(stdout)
  !!$omp critical
  !>-- Projects and mass-weights the Hessian
  call prj_mw_hess(mol%nat,mol%at,nmodes,mol%xyz,hess)
  !>-- Computes the Frequencies
  call frequencies(mol%nat,mol%at,mol%xyz,nmodes,hess,freq,io)
  !!$omp end critical
  write(stdout,'(a)') 'done.'

  deallocate (hess)
  return
end subroutine rdfreq_orca_hess

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
  use thermochem_module
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
  call prj_mw_hess(mol%nat,mol%at,nfreq,mol%xyz,hess)
  call frequencies(mol%nat,mol%at,mol%xyz,nfreq,hess,freq,io)

!>--- get thermodynamics
  et = 0.0_wp
  ht = 0.0_wp
  gt = 0.0_wp
  stot = 0.0_wp
  ithr = env%thermo%ithr
  fscal = env%thermo%fscal
  sthr = env%thermo%sthr
  call calcthermo(mol%nat,mol%at,mol%xyz,freq,pr,ithr,fscal,sthr, &
  &    nt,temps,et,ht,gt,stot,stdout,emodel=env%thermo%emodel)
  deallocate (hess,freq)
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
  use crest_parameters,only:wp,stdout,autokcal,aatoau
  use crest_data
  use strucrd
  use iomod
  use parallel_interface
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
  integer :: nt
  real(wp),allocatable :: temps(:)
  real(wp),allocatable :: et(:)
  real(wp),allocatable :: ht(:)
  real(wp),allocatable :: gt(:)
  real(wp),allocatable :: sref(:)
  character(len=64) :: atmp
  integer :: i,j,k,ich,io,popf,ii
  logical :: ex
  integer :: ncalc,nlimit,nav
  character(len=512) :: tmppath

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

! ── build coordinate subset and call parallel Hessian loop ───────
  block
    integer :: ii
    real(wp),allocatable :: xyz_calc(:,:,:),er_calc(:)
    real(wp),allocatable :: gt_out(:,:),stot_out(:,:)
    allocate (xyz_calc(3,nat,ncalc),er_calc(ncalc))
    allocate (gt_out(ncalc,nt),stot_out(ncalc,nt))
    do ii = 1,ncalc
      xyz_calc(:,:,ii) = xyz(:,:,pindex(ii))*aatoau  !> Å → Bohr
    end do
    write (stdout,'(1x,a,i0,a)') 'Running ',ncalc,' calculations ...'
    call crest_hessloop(env,nat,ncalc,at,xyz_calc,er_calc,gt_out,stot_out)
    do ii = 1,ncalc
      gatt(pindex(ii),1:nt) = gt_out(ii,1:nt)
      satt(pindex(ii),1:nt) = stot_out(ii,1:nt)
    end do
    deallocate (xyz_calc,er_calc,gt_out,stot_out)
  end block

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
