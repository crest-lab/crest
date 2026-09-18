!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022 - 2023  Philipp Pracht, Gereon Feldmann
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

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> Implementation of a numerical Hessian calculation
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
!> Input/Output:
!>  env  -  crest's systemdata object
!>  tim  -  timer object
!>-----------------------------------------------
subroutine crest_numhess(env,tim)
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
  use thermochem_module
  use gradreader_module
  use xtb_sc
  use oniom_hessian
  use ir_spectrum
  use iomod,only:colorify
  implicit none

  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,io,ich,nat3,n_freqs,i0
  logical :: pr,wr
!========================================================================================!
  type(calcdata) :: calc
  real(wp) :: energy
  real(wp),allocatable :: hess(:,:,:),freq(:,:),grad(:),grad1(:,:),grad2(:,:),heff(:,:)
  real(wp),allocatable :: ohess(:,:),ofreq(:),grad0(:,:),energies0(:)
  real(wp),allocatable :: ir_int(:)
  character(len=60) :: atmp

!========================================================================================!

  write (stdout,*)
  write (stdout,*) " ------------------------------------------------ "
  write (stdout,*) " #    # #    # #    # #    # ######  ####   ###   "
  write (stdout,*) " ##   # #    # ##  ## #    # #      #      #      "
  write (stdout,*) " # #  # #    # # ## # ###### #####   ####   ####  "
  write (stdout,*) " #  # # #    # #    # #    # #           #      # "
  write (stdout,*) " #   ## #    # #    # #    # #      #    # #    # "
  write (stdout,*) " #    #  ####  #    # #    # ######  ####   ####  "
  write (stdout,*) " ------------------------------------------------ "
  write (stdout,*)

!========================================================================================!
  call tim%start(15,'Numerical Hessian')

  call env%ref%to(mol)
  write (stdout,*)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)
!========================================================================================!

  call calc%copy(env%calc)

  !>--- print some info about the calculation
  call calc%info(stdout)

!========================================================================================!

  !>--- start with an initial single point
  write (stdout,'(a)') repeat(":",80)
  write (stdout,'(1x,a)') 'Initial singlepoint calculation ...'
  allocate (grad0(3,mol%nat),source=0.0_wp)
  allocate (energies0(calc%ncalculations),source=0.0_wp)

  call engrad(mol,calc,energy,grad0,io)
  energies0 = calc%etmp

  if (io .ne. 0) then
    write (stdout,*) 'FAILED! Initial singlepoint calculation did not converge.'
    env%iostatus_meta = status_failed
    deallocate (grad0)
    call tim%stop(15)
    return
  end if

  write (atmp,'("Energy = ",f25.15," Eh")') energy
  call smallhead(trim(atmp))
  write (stdout,'(a)') repeat(":",80)
  write (stdout,*)

  deallocate (grad0)

!========================================================================================!

  nat3 = mol%nat*3
  if (calc%ncalculations .eq. 1) then
    write (stdout,'(1x,a)',advance='no') 'Calculating numerical Hessian ...'
  else
    write (stdout,'(1x,a)',advance='no') 'Calculating numerical Hessians ...'
  end if
  flush (stdout)

!========================================================================================!
  if (.not.allocated(calc%ONIOM)) then
!========================================================================================!
!> Regular verision

    if (calc%id .eq. -1) then
      n_freqs = calc%ncalculations+1
    else
      n_freqs = calc%ncalculations
    end if

    allocate (hess(nat3,nat3,calc%ncalculations),source=0.0_wp)
    allocate (freq(nat3,n_freqs),source=0.0_wp)

! ── auto-enable dipole gradient for single tblite level ──────────────────────
    if (calc%ncalculations == 1) then
      if (calc%calcs(1)%id == jobtype%tblite) calc%calcs(1)%rddip = .true.
    end if

!*********************************************************************************
!>--- Computes numerical Hessians and stores them individually for each level
    call numhess2(mol%nat,mol%at,mol%xyz,calc,hess,io)
!*********************************************************************************

    write (stdout,*) 'done.'
    write (stdout,*)

    !> if calc type is set to -1: Computes the effective Hessian of the
    !> first two given calculation levels
    if (calc%id .eq. -1) then

      if (calc%ncalculations .gt. 1) then

        allocate (heff(nat3,nat3),source=0.0_wp)
        allocate (grad1(3,mol%nat),grad2(3,mol%nat),source=0.0_wp)

        grad1 = calc%grdtmp(:,:,1)
        grad2 = calc%grdtmp(:,:,2)

        !>-- Computes effective Hessian
        call effective_hessian(mol%nat,nat3,grad1,grad2,hess(:,:,1),hess(:,:,2),heff)

        !>-- Printout
        call print_hessian(heff,nat3,'','effhess')

        !>-- projection of translation and rotational modes + mass-weighting of Hessian
        call prj_mw_hess(mol%nat,mol%at,nat3,mol%xyz,heff)

        !>-- Comp. of Frequencies
        call frequencies(mol%nat,mol%at,mol%xyz,nat3,heff,freq(:,n_freqs),io)

        !>-- Printout of vibspectrum
        call print_vib_spectrum(mol%nat,mol%at,nat3,mol%xyz,freq(:,n_freqs),'','vibspectrum')

        !>-- Printout of g98 format file
        call print_g98_fake(mol%nat,mol%at,nat3,mol%xyz,freq(:,n_freqs),heff,'','g98.out')

        deallocate (heff)
        deallocate (grad1,grad2)

      else

        write (stdout,*) 'At least two calculation level must be'
        write (stdout,*) 'given for the calculation of the effective Hessian.'
        write (stdout,*)
        env%iostatus_meta = status_config
        return

      end if

    end if

    !> Prints hessian of numerical hessian and does a prj and mass weighting as well for the
    !> computation of normal modes and frequencies
    do i = 1,calc%ncalculations

      if (io .ne. 0) then
        write (stdout,*) 'FAILED!'
        env%iostatus_meta = status_failed

      else

        !> omit numeric suffix when there is only one calculation level
        if (calc%ncalculations == 1) then
          atmp = ''
        else
          write (atmp,*) i
          atmp = trim(adjustl(atmp))
        end if

        !>-- Prints Hessian
        call print_hessian(hess(:,:,i),nat3,'','numhess'//trim(atmp))

        !>--- Print dipole gradients (if they exist)
        call calc%calcs(i)%dumpdipgrad('dipgrad'//trim(atmp))

        !>-- Projects and mass-weights the Hessian
        call prj_mw_hess(mol%nat,mol%at,nat3,mol%xyz,hess(:,:,i))

        !>-- Computes the Frequencies
        call frequencies(mol%nat,mol%at,mol%xyz,nat3,hess(:,:,i),freq(:,i),io)

        if (io .ne. 0) then
          write (stdout,*) 'FAILED!'
          env%iostatus_meta = status_failed
        else

! ── project dipole gradient onto normal modes → IR intensities ───────────
          if (allocated(calc%calcs(i)%dipgrad)) then
            allocate (ir_int(nat3),source=0.0_wp)
            call ir_intensities(mol%nat,mol%at,nat3,hess(:,:,i), &
            &                   calc%calcs(i)%dipgrad,ir_int)
          end if

          !>-- Prints vibspectrum
          call print_vib_spectrum(mol%nat,mol%at,nat3,mol%xyz,freq(:,i), &
          &    '','vibspectrum'//trim(atmp),ir_int=ir_int)

          !>-- Prints g98.out format file (suffix prevents overwrite for empty calcspace)
          call print_g98_fake(mol%nat,mol%at,nat3,mol%xyz,freq(:,i),hess(:,:,i), &
          &    calc%calcs(i)%calcspace,'g98'//trim(atmp)//'.out',ir_int=ir_int)

          if (allocated(ir_int)) deallocate (ir_int)

          call smallhead("Thermo contributions for [[calculation.level]] "//trim(atmp))
          call numhess_thermostat(env,mol,nat3,hess(:,:,i),freq(:,i),energies0(i))

        end if

      end if
    end do
!========================================================================================!
  else
!========================================================================================!
!> ONIOM version
    write (stdout,*)
    flush (stdout)

    allocate (ohess(nat3,nat3),source=0.0_wp)
    allocate (ofreq(nat3),source=0.0_wp)
    call ONIOM_calc_hessians(mol,calc,ohess)

    !>-- Prints Hessian (pure 2nd derivatives in atomic units)
    call print_hessian(ohess(:,:),nat3,'','numhess')

    write (stdout,*) 'Note: The Hessian is printed as nat*3 blocks of nat*3 entries.'
    write (stdout,*)

    !>-- Projects and mass-weights the Hessian (M^1/2*H*M^1/2)
    call prj_mw_hess(mol%nat,mol%at,nat3,mol%xyz,ohess(:,:))

    !>-- Computes the Frequencies (in cm^-1)
    call frequencies(mol%nat,mol%at,mol%xyz,nat3,ohess(:,:),ofreq(:),io)
    if (io .ne. 0) then
      write (stdout,*) 'FAILED!'
      env%iostatus_meta = status_failed
    end if

    !>-- Prints vibspectrum (cm^-1) with artifical intensities
    call print_vib_spectrum(mol%nat,mol%at,nat3,mol%xyz,ofreq(:), &
    &    '','vibspectrum')

    !>-- Prints g98.out format file
    call print_g98_fake(mol%nat,mol%at,nat3,mol%xyz,ofreq(:),ohess(:,:), &
    &    '','g98.out')

    !>--- thermostatistical contributions
    call numhess_thermostat(env,mol,nat3,ohess,ofreq,energy)

!========================================================================================!
  end if
!========================================================================================!

!========================================================================================!
  if (allocated(hess)) deallocate (hess)
  if (allocated(freq)) deallocate (freq)
  if (allocated(ohess)) deallocate (ohess)
  if (allocated(ofreq)) deallocate (ofreq)
!========================================================================================!
  call tim%stop(15)

  return
!========================================================================================!
!========================================================================================!
end subroutine crest_numhess

!========================================================================================!

subroutine numhess_thermostat(env,mol,nat3,hess,freq,etot)
!*****************************************
!* A minimal wrapper of thermo to obtain
!* free energy contributions
!*****************************************
  use crest_parameters
  use crest_data
  use strucrd
  use thermochem_module
  implicit none
  !> INPUT
  type(systemdata) :: env
  type(coord),intent(inout) :: mol
  integer,intent(in) :: nat3
  real(wp),intent(in) :: hess(nat3,nat3)
  real(wp),intent(inout) :: freq(nat3)
  real(wp),intent(in) :: etot
  !> LOCAL
  real(wp) :: ithr,fscal,sthr
  character(len=:),allocatable :: emodel
  integer :: nt,nrt
  real(wp),allocatable :: temps(:),et(:),ht(:),stot(:),gt(:)
  real(wp) :: zpve

  !> inversion threshold
  ithr = env%thermo%ithr
  !> frequency scaling factor
  fscal = env%thermo%fscal
  !> RR-HO interpolation (or cut-off)
  sthr = env%thermo%sthr
  !> Svib model
  emodel = env%thermo%emodel

  if (.not.allocated(env%thermo%temps)) then
    call env%thermo%get_temps()
  end if
  nt = env%thermo%ntemps
  allocate (temps(nt),et(nt),ht(nt),gt(nt),stot(nt),source=0.0_wp)
  temps = abs(env%thermo%temps-298.15_wp)
  !write(*,*) temps
  nrt = minloc(temps(:),1)
  !write(*,*) nrt
  temps = env%thermo%temps

  !> calcthermo wants input in Bohr
  call calcthermo(mol%nat,mol%at,mol%xyz,freq,.true., &
  & ithr,fscal,sthr,nt,temps,et,ht,gt,stot,emodel=emodel)

  zpve = et(nrt)-ht(nrt)
  call print_thermo_summary(stdout,temps(nrt),etot,zpve,gt(nrt))

  deallocate (stot,gt,ht,et,temps)
end subroutine numhess_thermostat

!========================================================================================!

subroutine thermo_standalone(env)
!************************************************************
!* A minimal wrapper of thermo to obtain
!* re-do the thermostatistics either from a read-in hessian
!* or a vibspecturm file
!************************************************************
  use crest_parameters
  use crest_data
  use strucrd
  use thermochem_module
  implicit none
  !> INPUT
  type(systemdata) :: env
  !> LOCAL
  type(coord) :: mol
  integer :: nat3
  real(wp),allocatable :: hess(:,:)
  real(wp),allocatable :: freq(:)
  real(wp) :: etot
  real(wp) :: ithr,fscal,sthr
  character(len=:),allocatable :: emodel
  integer :: nt,nrt
  real(wp),allocatable :: temps(:),et(:),ht(:),stot(:),gt(:)
  real(wp) :: zpve
  integer :: ich,i

  !> header
  write (stdout,'(t10,a)') " _   _                               "
  write (stdout,'(t10,a)') "| |_| |__   ___ _ __ _ __ ___   ___  "
  write (stdout,'(t10,a)') "| __| '_ \ / _ \ '__| '_ ` _ \ / _ \ "
  write (stdout,'(t10,a)') "| |_| | | |  __/ |  | | | | | | (_) |"
  write (stdout,'(t10,a)') " \__|_| |_|\___|_|  |_| |_| |_|\___/ "
  write (stdout,'(t10,a)') "                                     "
  write (stdout,*) "Molecular thermodynamics from the modified and scaled"
  write (stdout,*) "rigid-rotor harmonic-oscillator approximation (msRRHO)"
  write (stdout,*) "See:"
  write (stdout,*) " • S.Grimme, Chem. Eur. J. 2012, 18, 9955–9964."
  write (stdout,*) " • P.Pracht, S.Grimme, Chem. Sci., 2021, 12, 6551-6568."
  write (stdout,*)

  !> input coords
  write (stdout,'(1x,a,t30)',advance='no') 'Reading input coords:'
  if (allocated(env%thermo%coords)) then
    call mol%open(env%thermo%coords)
    write (stdout,'(a)') trim(env%thermo%coords)
  else
    call mol%open(env%inputcoords)
    write (stdout,'(a)') trim(env%inputcoords)
  end if
  nat3 = mol%nat*3
  allocate (freq(nat3),source=0.0_wp)

  !> input frequencies or hessian
  if (allocated(env%thermo%vibfile)) then
    write (stdout,'(1x,a,t30,a)') 'Reading frequencies from:',trim(env%thermo%vibfile)
    call rdfreq(mol,env%thermo%vibfile,nat3,freq)
  else
    write (stdout,'(1x,a)') 'No Hessian or vibspectrum file allocated for thermo routine!'
    call creststop(status_input)
  end if
  write (stdout,*)

  !> energy (maybe read from comment line of xyz)
  etot = mol%energy
  !> inversion threshold
  ithr = env%thermo%ithr
  !> frequency scaling factor
  fscal = env%thermo%fscal
  !> RR-HO interpolation (or cut-off)
  sthr = env%thermo%sthr
  !> Svib model
  emodel = env%thermo%emodel

  if (.not.allocated(env%thermo%temps)) then
    call env%thermo%get_temps()
  end if
  nt = env%thermo%ntemps
  allocate (temps(nt),et(nt),ht(nt),gt(nt),stot(nt),source=0.0_wp)
  temps = abs(env%thermo%temps-298.15_wp)
  !write(*,*) temps
  nrt = minloc(temps(:),1)
  !write(*,*) nrt
  temps = env%thermo%temps

  !> calcthermo wants input in Bohr
  call calcthermo(mol%nat,mol%at,mol%xyz,freq,.true., &
  & ithr,fscal,sthr,nt,temps,et,ht,gt,stot,stdout,emodel=emodel)

  zpve = et(nrt)-ht(nrt)
  call print_thermo_summary(stdout,temps(nrt),etot,zpve,gt(nrt))

  !> for plotting temperature dependencies etc.
  write (stdout,*)
  write (stdout,*) 'Some output will be written to thermo.dump'
  open (newunit=ich,file='thermo.dump')
  do i = 1,nt
    write (ich,'(f12.4,4F20.10)') temps(i),gt(i)+etot,gt(i),ht(i),stot(i)
  end do
  close (ich)

  deallocate (stot,gt,ht,et,temps)
end subroutine thermo_standalone

!========================================================================================!

subroutine print_thermo_summary(iunit,temp,etot,zpve,grrho)
!**********************************************************
!* Print the standard THERMODYNAMICS summary box.
!* Called from numhess_thermostat and thermo_standalone.
!**********************************************************
  use crest_parameters,only:wp,stdout
  implicit none
  integer,intent(in) :: iunit
  real(wp),intent(in) :: temp,etot,zpve,grrho
  character(len=*),parameter :: outfmt = '(10x,"::",1x,a,f24.12,1x,a,1x,"::")'
  write (iunit,*)
  write (iunit,'(10x,a)') repeat(':',50)
  write (iunit,'(10x,"::",7x,a,f12.2,1x,a,8x,"::")') "THERMODYNAMICS at",temp,'K'
  write (iunit,'(10x,a)') repeat(':',50)
  write (iunit,outfmt) 'TOTAL FREE ENERGY',etot+grrho,'Eh'
  write (iunit,'(10x,a)') '::'//repeat('-',46)//'::'
  write (iunit,outfmt) 'total energy     ',etot,'Eh'
  write (iunit,outfmt) 'ZPVE             ',zpve,'Eh'
  write (iunit,outfmt) 'G(RRHO) w/o ZPVE ',grrho-zpve,'Eh'
  write (iunit,outfmt) 'G(RRHO) total    ',grrho,'Eh'
  write (iunit,'(10x,a)') repeat(':',50)
end subroutine print_thermo_summary

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!

subroutine crest_ensemble_hessians(env,tim)
!***********************************************
!* subroutine crest_ensemble_hessians
!* Standalone runtype: numerical Hessian +
!* thermochemistry for every structure in an
!* ensemble file. eread is overwritten with
!* Gibbs free energies (default temperature).
!* Input file : env%ensemblename
!* Output file: crest_ensemble.xyz (ensemblefile)
!***********************************************
  use crest_parameters,only:wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use parallel_interface
  use utilities,only:dumpenergies
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  character(len=:),allocatable :: ensnam
  integer :: nat,nall,T,Tn
  real(wp),allocatable :: xyz(:,:,:),eread(:),etmp(:)
  integer,allocatable  :: at(:)
  logical :: ex
!========================================================================================!
  write (stdout,*)
  inquire (file=env%ensemblename,exist=ex)
  if (ex) then
    ensnam = env%ensemblename
  else
    write (stdout,*) '**ERROR** no ensemble file provided.'
    env%iostatus_meta = status_input
    return
  end if

  call tim%start(14,'Ensemble Hessians')

  call rdensembleparam(ensnam,nat,nall)
  if (nall < 1) then
    write (stdout,*) '**ERROR** empty ensemble file.'
    env%iostatus_meta = status_input
    return
  end if
  allocate (xyz(3,nat,nall),at(nat),eread(nall),etmp(nall))
  call rdensemble(ensnam,nat,nall,at,xyz,eread)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: crest_hessloop requires coordinates in Bohr
  xyz = xyz/bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

  call new_ompautoset(env,'auto',nall,T,Tn)

!========================================================================================!
  write (stdout,*)
  write (stdout,'(10x,"┍",49("━"),"┑")')
  write (stdout,'(10x,"│",16x,a,16x,"│")') "ENSEMBLE HESSIANS"
  write (stdout,'(10x,"┕",49("━"),"┙")')
  write (stdout,*)
  write (stdout,'(1x,a,i0,a,1x,a)') 'Evaluating all ',nall,' structures of file ',trim(ensnam)

  call crest_hessloop(env,nat,nall,at,xyz,etmp)
  eread(:) = eread(:)+etmp(:)

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Back to Angstrom
  xyz = xyz*bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
  call wrensemble(ensemblefile,nat,nall,at,xyz,eread)
  write (stdout,'(/,a,a,a)') 'Ensemble with Gibbs free energies written to <',ensemblefile,'>'

  call dumpenergies('crest.energies',eread)
  write (stdout,'(/,a,a,a)') 'List of free energies written to <','crest.energies','>'

  deallocate (eread,at,xyz)
!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_ensemble_hessians

