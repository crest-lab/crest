!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2026 Philipp Pracht
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
!
! Routines were adapted from the xtb code (github.com/grimme-lab/xtb)
! under the Open-source software LGPL-3.0 Licencse.
!================================================================================!

module dynamics_module
  use crest_parameters
  use crest_calculator
  use iomod
  use strucrd
  use atmasses
  use shake_module
  use metadynamics_module
!$ use omp_lib
  implicit none

  !======================================================================================!
  !--- private module variables and parameters
  private

  !>--- some constants and name mappings
  real(wp),parameter :: amutoau = amutokg*kgtome
  !real(wp),parameter :: fstoau = 41.3413733365614_wp

  !>-- filetypes as integers
  integer,parameter,public :: type_md = 1
  integer,parameter,public :: type_mtd = 2

  !>-- REEXPORTS from metadynamics_module
  public :: mtdpot,mtd_ini,cv_dump,calc_mtd
  public :: cv_std_mtd,cv_rmsd,cv_rmsd_static

  !>-- REEXPORTS from shake_module
  public :: shakedata,init_shake

  public :: mddata
  !======================================================================================!
  !data object that contains settings for a molecular dynamics simulation.
  type :: mddata

    logical :: requested = .false.

    integer :: md_index = 0      !> some index for parallelization
    integer :: simtype = type_md !> type of the molecular dynamics simulation
    logical :: restart = .false.
    character(len=:),allocatable :: restartfile
    logical :: wrtrj = .true.
    character(len=:),allocatable :: trajectoryfile
    !>--- data
    real(wp) :: length_ps = 0.0_wp !20.0_wp  !> total simulation length in ps
    integer :: length_steps = 0 !20000       !> total simulation length in steps
    real(wp) :: tstep = 0.0_wp !1.0_wp       !> timestep in fs
    real(wp) :: dumpstep = 0.0_wp !1000.0_wp !> snapshot dump step in fs
    integer :: sdump = 0 !1000 !> snapshot dump to trajectory every x steps
    integer :: dumped = 0      !> count how many structures have been written
    integer :: printstep = 200 !> control how often (in steps) to print to stdout

    real(wp) :: md_hmass = 0.0_wp !1.00794075_wp !> hydrogen mass

    logical :: shake = .true. !> use SHAKE algorithm
    integer :: nshake = 0  !> number of bonds considered in SHAKE
    type(shakedata) :: shk !> SHAKE bond information

    real(wp) :: tsoll = 0.0_wp !298.15_wp  !> wanted temperature
    real(wp) :: Tavg = 0.0_wp  !> trajectory-average temperature (set after dynamics())
    real(wp) :: Tvar = 0.0_wp  !> trajectory temperature variance (set after dynamics())
    integer  :: Ndf = 0       !> degrees of freedom used in the run (set after dynamics())
    logical :: thermostat = .true. !> apply thermostat?
    character(len=64) :: thermotype = 'berendsen'
    integer :: thermotype_i = 2 !> integer mapping thermostats, 2=berendsen
    real(wp) :: thermo_damp = 500.0_wp !> thermostat damping parameter
    logical :: samerand = .false.

    integer :: blockl = 0    !> block length in MD steps
    integer :: iblock = 0    !> block counter
    integer :: nblock = 0    !> continuous block counter
    integer :: blocknreg = 0 !> block regression points
    integer :: maxblock
    real(wp),allocatable :: blockrege(:)
    real(wp),allocatable :: blocke(:)
    real(wp),allocatable :: blockt(:)

    !>--- the collection of MTD potentials
    integer :: npot = 0
    type(mtdpot),allocatable :: mtd(:)
    integer,allocatable :: cvtype(:)

    !>--- on-the-fly multi-scale level choosing
    integer,allocatable :: active_potentials(:)

  contains
    generic,public :: add => md_add_mtd
    procedure,private :: md_add_mtd
    procedure :: defaults => md_defaults_fallback
  end type mddata

  public :: dynamics
  public :: mdautoset

  character(len=30),parameter,private :: thermostattype(4) = [ &
     & 'None                         ', &
     & 'Berendsen                    ', &
     & 'Langevin                     ', &
     & 'Bussi-Donadio-Parrinello     ']

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine dynamics(mol,dat,calc,pr,term)
!*************************************************************
!* subroutine dynamics
!* perform a molecular dynamics simulation
!* the coordinate propagation is made with an
!* Leap-Frog algorithm (Verlet-type algo)
!*************************************************************
    implicit none

    type(coord) :: mol          !> molecule data (should be in Bohr)
    type(mddata) :: dat         !> MD data
    type(calcdata),target :: calc      !> calculation control
    logical,intent(in) :: pr    !> printout control
    integer,intent(out) :: term !> termination status

    integer :: t,nfreedom
    real(wp) :: tstep_au
    real(wp) :: epot,ekin,edum
    real(wp) :: temp,thermoscal
!>--- averages & errors
    real(wp) :: Tav,Tav2,Epav,Ekav,Eerror
!>--- block statistics (computed before printout)
    real(wp) :: eblk_mean,eblk_std,eblk_drift

    real(wp),allocatable :: grd(:,:)
    real(wp),allocatable :: velo(:,:)
    real(wp),allocatable :: vel(:,:)
    real(wp),allocatable :: veln(:,:)
    real(wp),allocatable :: acc(:,:)
    real(wp),allocatable :: mass(:)
    real(wp),allocatable :: backupweights(:)
    type(coord) :: molo
    real(wp) :: f,rt,rtshift
    real(wp) :: molmass,tmass
    real(wp),allocatable :: cn(:)
    character(len=:),allocatable :: trajectory
    integer :: trj
    character(len=256) :: commentline
    integer :: i,j,k,l,ich,och,io
    integer :: dcount,printcount
    integer,allocatable :: iseed(:)
    logical :: ex,fail,bdump,shakefallback

    call initsignal()

!>--- pre-settings and calculations
    !$omp critical
    call dat%defaults() !> check for unset parameters
!>--- seed RNG before any stochastic use (velocity init or thermostat steps)
    if (dat%samerand) then
      call random_seed(size=i)
      allocate (iseed(i),source=1)
      call random_seed(put=iseed)
      deallocate (iseed)
    else
      call random_seed()
    end if
    term = 0
    tstep_au = dat%tstep*fstoau
    nfreedom = 3*mol%nat
    shakefallback = .false.
    if (calc%nfreeze > 0) then
      nfreedom = nfreedom-3*calc%nfreeze
    end if
    if (dat%shake) then
      if (calc%nfreeze > 0) then
        dat%shk%freezeptr => calc%freezelist
      else
        nullify (dat%shk%freezeptr)
      end if
      if (.not.allocated(dat%shk%wbo)) then
        !> fake bondorder fallback for shake
        shakefallback = .true.
        call mol%cn_to_bond(cn,dat%shk%wbo)
        deallocate (cn)
      end if
      call init_shake(mol%nat,mol%at,mol%xyz,dat%shk,pr)
      dat%nshake = dat%shk%ncons
      nfreedom = nfreedom-dat%nshake
    end if
!>--- averages
    tav = 0.0_wp
    tav2 = 0.0_wp
    eerror = 0.0_wp
    ekav = 0.0_wp
    epav = 0.0_wp
    temp = 0.0_wp

!>--- on-the-fly multiscale definition
    if (allocated(dat%active_potentials)) then
      call calc%active(dat%active_potentials)
    end if

!>--- allocate data fields
    allocate (molo%at(mol%nat),molo%xyz(3,mol%nat))
    allocate (grd(3,mol%nat),vel(3,mol%nat),velo(3,mol%nat),source=0.0_wp)
    allocate (veln(3,mol%nat),acc(3,mol%nat),mass(mol%nat),source=0.0_wp)
    allocate (dat%blocke(dat%blockl),dat%blockt(dat%blockl))
    allocate (dat%blockrege(dat%maxblock))
    !$omp end critical

!>--- settings printout
    if (pr) then
      write (stdout,*)
      write (stdout,'(1x,15("─"),1x,a,1x,14("─"))') 'Molecular Dynamics Settings'
      write (stdout,'("  Simulation type",t25,       ":",1x)',advance='no')
      if (dat%thermostat) then
        write (stdout,'(a9)') 'NVT'
      else
        write (stdout,'(a9)') 'NVE'
      end if
      write (stdout,'("  MD time (length)",t25,       ":",f10.2,a)') dat%length_ps,' ps'
      write (stdout,'("  dt (timetep)",t25,            ":",f10.2,a)') dat%tstep,' fs'
      write (stdout,'("  temperature",t25,    ":",f10.2,a)') dat%tsoll,' K'
      write (stdout,'("  max steps",t25,         ":",i10  )') dat%length_steps
      write (stdout,'("  block length (av.)",t25,":",i10  )') dat%blockl
      write (stdout,'("  dumpstep(trj) /fs",t25, ":",f10.2,1x,"(",i0,")")') dat%dumpstep,dat%sdump
      write (stdout,'("  # deg. of freedom",t25, ":",i10  )') nfreedom
      if (calc%nfreeze > 0) then
        write (stdout,'("  # frozen atoms",t25,     ":",i10  )') calc%nfreeze
      end if
      call thermostatprint(dat,pr)
      write (stdout,'("  SHAKE constraint",t25,   ":",1x,a9)') to_str(dat%shake)
      if (dat%shake) then
        if (shakefallback) then
          write (stdout,'("  SHAKE using CN fallback",t25,":",9x,l)') shakefallback
        end if
        if (dat%shk%shake_mode == 2) then
          write (stdout,'("  # SHAKE bonds",t25,":",i10,a)') dat%nshake,' (all bonds)'
        elseif (dat%shk%shake_mode == 1) then
          write (stdout,'("  # SHAKE bonds",t25,":",i10,a)') dat%nshake,' (H only)'
        end if
      end if
      write (stdout,'("  hydrogen mass /u",t25,":",f10.5 )') dat%md_hmass
      if (allocated(dat%active_potentials)) then
        write (stdout,'("  active potentials",t25,":",i10)') size(dat%active_potentials,1)
      end if

    end if

!>--- set atom masses
    !$omp critical
    molmass = 0.0_wp
    do i = 1,mol%nat
      molmass = molmass+ams(mol%at(i))
      mass(i) = ams(mol%at(i))*amutoau !>-- ams from module atmasses
    end do
    tmass = molmass*amutoau
    do i = 1,mol%nat
      if (mol%at(i) .eq. 1.and.dat%md_hmass .gt. 0.0_wp) then
        mass(i) = dat%md_hmass*amutoau
      end if
    end do
    molmass = molmass*amutokg
    !$omp end critical

!>--- initialize velocities (or read from restart file)
    if (dat%thermostat) then
      f = 1.0_wp
    else
      f = 2.0_wp
    end if
    edum = f*dat%tsoll*0.5_wp*kB*float(nfreedom)
    rtshift = 0.0_wp
    if (.not.dat%restart.or..not.allocated(dat%restartfile)) then
      call mdinitu(mol,dat,velo,mass,edum,pr)
    else
      call rdmdrestart(mol,dat,velo,fail,rtshift,pr)
      if (fail) then
        call mdinitu(mol,dat,velo,mass,edum,pr)
      else
        call ekinet(mol%nat,velo,mass,ekin)
        temp = 2.0_wp*ekin/float(nfreedom)/kB
        tav = temp
        tav2 = temp**2
      end if
    end if
    call ekinet(mol%nat,velo,mass,ekin)
    if (calc%nfreeze > 0) then
      do i = 1,mol%nat
        if (calc%freezelist(i)) then
          acc(:,i) = 0.0_wp
          grd(:,i) = 0.0_wp
          velo(:,i) = 0.0_wp
        end if
      end do
    end if

!>--- initialize MTDs (if required)
    !$omp critical
    if (dat%simtype == type_mtd) then
      call md_init_mtd(mol,dat,pr)
    end if
    !$omp end critical

    !>--- initialize trajectory file
    if (allocated(dat%trajectoryfile)) then
      trajectory = dat%trajectoryfile
    else
      write (commentline,'(a,i0,a)') 'crest_',dat%md_index,'.trj'
      trajectory = trim(commentline)
    end if
    !$omp critical
    if (dat%wrtrj) open (newunit=trj,file=trajectory)
    !>--- under PBC force extxyz trajectory output so the lattice (Lattice=...)
    !>    survives into any downstream re-read/re-optimization of the frames
    mol%wrextxyz = calc%logextxyz .or. allocated(mol%lat)
    !$omp end critical

!>--- begin printout
    if (pr) then
      write (stdout,'(1x,58("─"))')

      write (stdout,'(/,"> ",a)') 'Starting simulation'
      if (.not.dat%thermostat) then
        write (stdout,'(/,11x,"time (ps)",7x,"<Epot>",8x,"Ekin",5x,"<T>",7x,"T",12x, &
           &         "Etot",7x,"error")')
      else
        write (stdout,'(/,11x,"time (ps)",7x,"<Epot>",8x,"Ekin",5x,"<T>",7x,"T",12x, &
           &         "Etot")')
      end if
    end if

    dcount = 0
    printcount = 1
!===============================================================!
!===============================================================!
!>--- begin MD loop
    MD: do t = 1,dat%length_steps
      call initsignal()

      !>>-- STEP 1: calculate energy and forces
      !>--- singlepoint calculation
      epot = 0.0_wp
      grd = 0.0_wp
      call engrad(mol,calc,epot,grd,io)
      mol%energy = epot

      if (io /= 0) then
        if (dat%dumped > 0) then
          term = 2  !> termination during MD
        else
          term = 1  !> termination upon first engrad call
        end if
        exit MD
      end if
      if (t == 1) then
        edum = epot+ekin
      end if

      !>>-- STEP 1.5: calculate metadynamics bias
      if (dat%simtype == type_mtd) then
        !> MTD energy and gradient are added to epot and grd, respectively.
        call md_calc_mtd(mol,dat,epot,grd,pr)
        !$omp critical
        call md_update_mtd(mol,dat,calc,pr)
        !$omp end critical
      end if

      !>--- block data printouts
      call u_block(mol,dat,epot,temp,pr,bdump)
      !>--- MD restart files written for each block rather than at each timestep to reduce I/O
      if (bdump) then
        rt = float(t)*dat%tstep+rtshift
        !$omp critical
        call wrmdrestart(mol,dat,velo,rt)
        !$omp end critical
      end if

      !===========================================!
      !>>-- write to trajectory and printout
      if (dcount == dat%sdump) then
        dcount = 0
        dat%dumped = dat%dumped+1
        !$omp critical
        if (dat%wrtrj) call mol%append(trj)
        !$omp end critical
      end if
      if ((printcount == dat%printstep).or.(t == 1)) then
        if (t > 1) printcount = 0
        if (pr) then
          rt = float(t)*dat%tstep+rtshift
          if (.not.dat%thermostat) then
            write (stdout,'(i7,f10.2,F16.5,F12.4,2F8.1,F16.5,4F10.4)') &
               &   t,0.001_wp*rt, (Epav+Epot)/float(t), &
               &   Ekin,Tav/float(t),temp,Epot+Ekin, &
               &   Edum/float(t)-Epot-Ekin
          else
            write (stdout,'(i7,f10.2,F16.5,F12.4,2F8.1,F16.5)') &
               &   t,0.001_wp*rt, (Epav+epot)/float(t), &
               &   Ekin,Tav/float(t),temp,Epot+Ekin
          end if
        end if
      end if

      !>--- compute the acceleration at t
      do i = 1,mol%nat
        acc(:,i) = -grd(:,i)/mass(i)
      end do

      !>--- special setup for frozen atoms
      if (calc%nfreeze > 0) then
        do i = 1,mol%nat
          if (calc%freezelist(i)) then
            acc(:,i) = 0.0_wp
            grd(:,i) = 0.0_wp
          end if
        end do
      end if

      !>--- store positions (at t); velocities are at t-1/2dt
      !$omp critical
      molo%nat = mol%nat
      molo%at = mol%at
      molo%xyz = mol%xyz
      !$omp end critical

      !>>-- STEP 2: temperature and pressure/density control
      !>--- estimate(!) velocities at t
      veln = velo+0.5_wp*tstep_au*acc

      !>--- compute kinetic energy and termperature
      call ekinet(mol%nat,veln,mass,ekin)
      temp = 2.0_wp*ekin/float(nfreedom)/kB

      !>--- THERMOSTATING and velocity update
      if (dat%thermotype_i == 3) then !'langevin') then
        if (dat%thermostat) then
          call langevin_step(mol%nat,dat,mass,velo,acc,tstep_au,vel)
        else
          vel = velo+acc*tstep_au
        end if
      else
        !>>-- STEP 3: velocity update for scaling thermostats
        !>--- compute trial velocity (unscaled leapfrog step)
        vel = velo+acc*tstep_au
        !>--- Bussi/CSVR: recompute ekin from the actual trial velocity so that
        !>--- the scaling scal=sqrt(K_new/ekin) gives Ekin(vel_scaled)=K_new exactly.
        !>--- Using the half-step estimate veln would give Ekin(vel_scaled)=K_new*(Efull/Ehalf)
        !>--- and systematically overshoot the target temperature by ~Efull/Ehalf.
        if (dat%thermostat.and. &
           & dat%thermotype_i == 4) then
          !& (trim(dat%thermotype) == 'bussi'.or.trim(dat%thermotype) == 'csvr')) then
          call ekinet(mol%nat,vel,mass,ekin)
          temp = 2.0_wp*ekin/float(nfreedom)/kB
        end if
        call thermostating(mol,dat,temp,ekin,nfreedom,thermoscal)
        vel = thermoscal*vel
      end if

      !>--- update positions to t+dt, except for frozen atoms, and not at the final step
      if (t < dat%length_steps) then
        if (calc%nfreeze > 0) then
          do i = 1,mol%nat
            if (.not.calc%freezelist(i)) then
              mol%xyz(:,i) = molo%xyz(:,i)+vel(:,i)*tstep_au
            else
              vel(:,i) = 0.0_wp
            end if
          end do
        else
          mol%xyz = molo%xyz+vel*tstep_au
        end if
      end if

      !>--- estimate new velocities at t
      veln = 0.5_wp*(velo+vel)

      !>--- compute kinetic energy and temperature for average tracking
      call ekinet(mol%nat,veln,mass,ekin)
      temp = 2.0_wp*ekin/float(nfreedom)/kB

      !>--- apply SHAKE at t+dt?
      if (dat%shake.and.dat%nshake > 0) then
        call do_shake(mol%nat,molo%xyz,mol%xyz,vel,acc,mass,tstep_au,dat%shk,pr,io)
        if (io /= 0) then
          if (dat%dumped > 0) then
            term = 2  !> termination during MD
          else
            term = 1  !> termination upon first engrad call
          end if
          exit MD
        end if
      end if

      !>--- update velocities
      velo = vel

      !>--- remove translational and rotational componetnts of the velocity
      !>    (under PBC the fixed cell forbids rotation removal -> translation only)
      call rmrottr(mol%nat,mass,velo,mol%xyz,allocated(mol%lat))

      !>>-- Update averages and counter
      edum = edum+epot+ekin
      eerror = edum/float(t)-epot-ekin
      tav = tav+temp
      tav2 = tav2+temp**2
      epav = epav+epot
      ekav = ekav+ekin
      dcount = dcount+1
      printcount = printcount+1

    end do MD
!>--- finish MD loop
!===============================================================!
!===============================================================!
!>--- close trajectory file
    !$omp critical
    if (dat%wrtrj) close (trj)
    !$omp end critical

!>--- block energy statistics (blockrege still allocated here)
    eblk_mean = 0.0_wp
    eblk_std = 0.0_wp
    eblk_drift = 0.0_wp
    if (dat%nblock >= 1) then
      eblk_mean = sum(dat%blockrege(1:dat%nblock))/real(dat%nblock,wp)
    end if
    if (dat%nblock >= 2) then
      eblk_std = sqrt(sum((dat%blockrege(1:dat%nblock)-eblk_mean)**2) &
               &      /real(dat%nblock-1,wp))
      !> drift in Eh/ps: endpoint slope over all completed blocks
      eblk_drift = (dat%blockrege(dat%nblock)-dat%blockrege(1)) &
                 & /real(dat%nblock-1,wp) &
                 & /(dat%blockl*dat%tstep*1.0e-3_wp)
    end if

!>--- averages printout
    if (pr) then
      write (stdout,*)
      write (stdout,'(1x,a)') 'average properties'
      write (stdout,'(1x,a)') repeat('-',42)
      write (stdout,'(1x,a,t22,f16.8,a)') '<Epot>',Epav/float(t),' Eh'
      write (stdout,'(1x,a,t22,f16.8,a)') '<Ekin>',Ekav/float(t),' Eh'
      write (stdout,'(1x,a,t22,f16.8,a)') '<Etot>', (Ekav+Epav)/float(t),' Eh'
      write (stdout,'(1x,a,t28,f10.2,a)') '<T>',Tav/float(t),' K'
      write (stdout,'(1x,a,t28,f10.2,a)') 'Tvar',tav2/float(t)-(tav/float(t))**2,' K²'
      write (stdout,'(1x,a,t29,f10.2,a)') 'σ(T)',sqrt(max(tav2/float(t)-(tav/float(t))**2,0.0_wp)),' K'
      if (dat%nblock >= 2) then
        write (stdout,'(1x,a)') repeat('-',42)
        write (stdout,'(1x,a,t28,i10,a)') 'blocks',dat%nblock,' '
        write (stdout,'(1x,a,t23,f16.8,a)') 'block σ(Epot)',eblk_std,' Eh'
        write (stdout,'(1x,a,t28,es10.2,a)') 'drift',eblk_drift,' Eh/ps'
      end if
      write (stdout,'(1x,a)') repeat('-',42)
    end if

!>--- write restart file
    rt = float(dat%length_steps)*dat%tstep+rtshift
    call wrmdrestart(mol,dat,velo,rt)

!>--- termination printout
    if (pr) then
      select case (term)
      case (0)
        write (stdout,*) 'normal MD termination'
      case (1)
        write (stderr,*) 'error in MD calculation'
      case (2)
        write (stdout,*) 'MD terminated, but still taking as converged.'
      end select
    end if

!>--- store trajectory-average temperature for callers
    dat%Tavg = tav/float(t)
    dat%Tvar = tav2/float(t)-(tav/float(t))**2
    dat%Ndf = nfreedom

!>--- deallocate data
    deallocate (dat%blockrege,dat%blockt,dat%blocke)
    deallocate (mass,acc,veln)
    deallocate (vel,velo,grd)
    deallocate (molo%xyz,molo%at)

!>--- restore weights if necessary
    if (allocated(dat%active_potentials)) then
      call calc%active_restore()
    end if

    return
  end subroutine dynamics

!========================================================================================!
  subroutine mdautoset(dat,iostatus)
!************************************************
!* subroutine mdautoset
!* convert real-time settings (ps,fs) to steps
!* Intended to restore settings once reset
!************************************************
    implicit none
    type(mddata) :: dat
    integer,intent(out) :: iostatus
    real(wp) :: dum
    integer :: idum
    integer :: i,j,k,l,ich,och,io
    logical :: ex

    iostatus = 0

    if (dat%length_ps .le. 0.0_wp.or. &
    &  dat%tstep .le. 0.0_wp) then
      write (stderr,*) 'need valid simulation length and time step!'
      write (stderr,*) 'abort MD.'
      iostatus = -1
      return
    end if

    !>--- MD length to steps
    dum = (dat%length_ps*1000.0_wp)/dat%tstep
    dat%length_steps = nint(dum)

    !>--- adjust structure dump to trajectory file
    dum = max(1.0_wp, (dat%dumpstep/dat%tstep))
    dat%sdump = nint(dum)

    return
  end subroutine mdautoset

!========================================================================================!
! subroutine ekinet
! calculate kinetic energy from velocities and masses
  subroutine ekinet(n,velo,mass,E)
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: velo(3,n),mass(n)
    real(wp),intent(out) :: e
    integer :: i
    e = 0.0_wp
    do i = 1,n
      e = e+mass(i)*(velo(1,i)**2+velo(2,i)**2+velo(3,i)**2)
    end do
    e = e*0.5_wp
    return
  end subroutine ekinet

!========================================================================================!
! subroutine random_gauss
! sample one standard normal N(0,1) via Box-Muller transform
  subroutine random_gauss(z)
    !*************************************
    !* Sample one standard normal N(0,1) *
    !* via Box-Muller transform.         *
    !*************************************
    implicit none
    real(wp),intent(out) :: z
    real(wp) :: u1,u2
    real(wp),parameter :: twopi = 6.28318530717958647693_wp
    do
      call random_number(u1)
      if (u1 > 0.0_wp) exit
    end do
    call random_number(u2)
    z = sqrt(-2.0_wp*log(u1))*cos(twopi*u2)
    return
  end subroutine random_gauss

!========================================================================================!
! subroutine random_gamma
! sample x ~ Gamma(a,1) using the Marsaglia-Tsang (2000) algorithm
  subroutine random_gamma(a,x)
    !*****************************************************
    !* Sample x ~ Gamma(a,1) using Marsaglia-Tsang 2000. *
    !* For a < 1: apply the relation Gamma(a) = Gamma(a+1)*U^(1/a). *
    !*****************************************************
    implicit none
    real(wp),intent(in)  :: a
    real(wp),intent(out) :: x
    real(wp) :: d,c,z,v,u,aa
    logical  :: small
    small = (a < 1.0_wp)
    aa = merge(a+1.0_wp,a,small)
    d = aa-1.0_wp/3.0_wp
    c = 1.0_wp/sqrt(9.0_wp*d)
    do
      do
        call random_gauss(z)
        v = (1.0_wp+c*z)**3
        if (v > 0.0_wp) exit
      end do
      call random_number(u)
      if (u < 1.0_wp-0.0331_wp*z**2*z**2) then
        x = d*v; exit
      end if
      if (log(u) < 0.5_wp*z**2+d*(1.0_wp-v+log(v))) then
        x = d*v; exit
      end if
    end do
    if (small) then
      call random_number(u)
      x = x*u**(1.0_wp/a)
    end if
    return
  end subroutine random_gamma

!========================================================================================!
! subroutine u_block
! update block data and printout
  subroutine u_block(mol,dat,epot,temp,pr,bdump)
    implicit none
    type(coord) :: mol
    type(mddata) :: dat
    real(wp),intent(in) :: epot
    real(wp),intent(in) :: temp
    logical,intent(in) :: pr
    logical,intent(out) :: bdump
    integer :: i,j,k,l,ich,och,io
    logical :: ex
    integer :: nreg
    real(wp) :: bave,bavt,slope

    bdump = .false.
    if (dat%iblock == dat%blockl) then
      dat%nblock = dat%nblock+1
      dat%iblock = 0
      bdump = .true.
      call blocksd(mol%nat,dat%blockl,dat%blocke,dat%blockt,bave,bavt)
      dat%blocknreg = dat%blocknreg+1
      nreg = dat%blocknreg
      dat%blockrege(nreg) = bave
      if (nreg .ge. 4) then
        call regress(nreg-3,nreg,dat%blockrege,slope)
      else
        slope = 99.0_wp
      end if
      if (pr) then
        write (stdout,'("block <Epot> / <T> :",f14.5,f7.1,4x, &
           &             "drift:",d10.2,3x,"Tbath :",f6.1)')  &
           &             bave,bavt,slope,dat%tsoll
      end if
    else
      dat%iblock = dat%iblock+1
      dat%blocke(dat%iblock) = epot
      dat%blockt(dat%iblock) = temp
    end if

    return

  contains
    subroutine regress(n1,n2,rege,slope)
      implicit none
      real(wp) :: rege(*),slope
      integer :: n1,n2,n
      real(wp) :: sx,sy,sxx,sxy,x
      integer :: i,j,k,l,ich,och,io

      n = n2-n1+1
      sx = 0.0_wp
      sy = 0.0_wp
      sxy = 0.0_wp
      sxx = 0.0_wp
      x = 0.0_wp
      do i = n1,n2
        x = x+1.0_wp
        sx = sx+x
        sy = sy+rege(i)
        sxx = sxx+x**2
        sxy = sxy+x*rege(i)
      end do

      slope = (dble(n)*sxy-sx*sy)/(dble(n)*sxx-sx*sx)
      return
    end subroutine regress

    subroutine blocksd(n,nbl,ebl,tbl,esd,tsd)
      implicit none
      integer :: n,nbl
      real(wp) :: ebl(nbl),tbl(nbl),dum,av,esd,tsd
      integer :: i

      dum = 0.0_wp
      do i = 1,nbl
        dum = dum+ebl(i)
      end do
      av = dum/dble(nbl)
      esd = av

      dum = 0.0_wp
      do i = 1,nbl
        dum = dum+tbl(i)
      end do
      av = dum/dble(nbl)
      tsd = av

      return
    end subroutine blocksd
  end subroutine u_block

!========================================================================================!
! subroutines wrmdrestart & rdmdrestart
! write a file containing coordinates and velocities to restart the simulation
  subroutine wrmdrestart(mol,dat,velo,realtime_fs)
    implicit none
    type(coord) :: mol
    type(mddata) :: dat
    real(wp),intent(in) :: velo(3,mol%nat)
    real(wp),intent(in) :: realtime_fs
    integer :: i,j,k,l,ich,och,io
    logical :: ex
    character(len=256) :: atmp
    if (.not.allocated(dat%restartfile).or.dat%restart) then
      !>--- we must not overwrite the user-provided restart file!
      write (atmp,'(a,i0,a)') 'crest_',dat%md_index,'.mdrestart'
    else
      atmp = dat%restartfile
    end if
    open (newunit=ich,file=trim(atmp)) !dat%restartfile)
    write (ich,*) realtime_fs
    do i = 1,mol%nat
      write (ich,'(6D22.14)') mol%xyz(1:3,i),velo(1:3,i)
    end do
    close (ich)
    return
  end subroutine wrmdrestart

  subroutine rdmdrestart(mol,dat,velo,fail,rtshift,pr)
    implicit none
    type(coord) :: mol
    type(mddata) :: dat
    real(wp),intent(inout) :: velo(3,mol%nat)
    logical,intent(out) :: fail
    real(wp),intent(out) :: rtshift
    logical,intent(in) :: pr
    real(wp) :: dum
    character(len=256) :: atmp
    integer :: i,j,k,l,ich,och,io
    logical :: ex

    fail = .false.
    rtshift = 0.0_wp
    if (allocated(dat%restartfile)) then
      inquire (file=dat%restartfile,exist=ex)
    end if
    if (dat%restart.and.ex) then
      open (newunit=ich,file=dat%restartfile)
      do
        read (ich,*,iostat=io) dum
        if (io < 0) exit
        if (dum > 0.0_wp) rtshift = dum
        do i = 1,mol%nat
          read (ich,'(a)',iostat=io) atmp
          if (io < 0) exit
          read (atmp,'(6D22.14)',iostat=io) mol%xyz(1:3,i),velo(1:3,i)
          if (io /= 0) exit
        end do
        exit
      end do
      if (io /= 0) then
        write (0,*) 'failed to read MD restart file.'
        fail = .true.
      end if
      close (ich)
    else
      fail = .true.
    end if
    if (.not.fail.and.pr) then
      write (stdout,'("  read RESTART file",t25,":",1x,a9)') to_str(.not.fail)
      write (stdout,'("  restart file",t25,":",1x,a)') dat%restartfile
    end if

    return
  end subroutine rdmdrestart

!========================================================================================!
! subroutine mdinitu
! initialize MD velocities, either by reading a restart file
! or by setting them randomly
  subroutine mdinitu(mol,dat,velo,mass,Ekin,pr)
    implicit none
    type(coord) :: mol
    type(mddata) :: dat
    real(wp),intent(inout) :: velo(3,mol%nat)
    real(wp),intent(in) :: mass(mol%nat)
    real(wp),intent(in) :: Ekin
    logical,intent(in) :: pr
    real :: x(3),ranf
    real(wp) :: eperat,v,f,t,edum,f2
    logical :: newvelos
    integer :: i

    newvelos = .true.

    !>--- newly initialized
    if (newvelos) then
      eperat = Ekin/(3.0_wp*float(mol%nat))
      do i = 1,mol%nat
        call random_number(x)
        f2 = 1.0_wp
        if (mol%at(i) .eq. 1) f2 = 2.0_wp
        v = sqrt(2*eperat/mass(i))
        f = 1.0_wp
        if (x(1) .gt. 0.5_wp) f = -1.0_wp
        velo(1,i) = v*f*f2
        f = 1.0_wp
        if (x(2) .gt. 0.5_wp) f = -1.0_wp
        velo(2,i) = v*f*f2
        f = 1.0d0
        if (x(3) .gt. 0.5_wp) f = -1.0_wp
        velo(3,i) = v*f*f2
      end do
    end if
    call ekinet(mol%nat,velo,mass,edum)
    t = edum/(0.5_wp*3.0_wp*float(mol%nat)*0.316681534524639E-05)
    return
  end subroutine mdinitu

!========================================================================================!
! subroutine thermostating
! helper routine to re-scale velocities,
! i.e., thermostating

  subroutine thermostating(mol,dat,t,ekin,nfreedom,scal)
    !**************************************************************
    !* Apply velocity-scaling thermostat and return scale factor. *
    !* Supports: berendsen, bussi/csvr                            *
    !**************************************************************
    implicit none
    type(coord) :: mol
    type(mddata) :: dat
    real(wp),intent(in) :: t
    real(wp),intent(in) :: ekin
    integer,intent(in)  :: nfreedom
    real(wp),intent(out) :: scal
    real(wp) :: c1,K_ref,K_new,gam,xi,chi2
    integer  :: Nf

    scal = 1.0_wp

    if (.not.dat%thermostat) return

    select case (dat%thermotype_i) !(trim(dat%thermotype))
    case (2) !('berendsen')
      scal = sqrt(1.0d0+(dat%tstep/dat%thermo_damp) &
                 &     *(dat%tsoll/t-1.0_wp))
    case (4) !('bussi','csvr')
      !> Bussi-Donadio-Parrinello CSVR (J. Chem. Phys. 126, 014101, 2007)
      !> Stochastic velocity rescaling for correct NVT ensemble sampling.
      !>
      !> Exact discrete update (paper Eq. A7):
      !>   K_new = K*c1 + K_ref*(1-c1)/Nf * sum_i(z_i^2) + 2*sqrt(K*K_ref*c1*(1-c1)/Nf)*z1
      !>   E[K_new] = K*c1 + K_ref*(1-c1)  →  converges to K_ref = Nf/2*kB*T
      !>
      !> Setting A = sqrt(K_ref*(1-c1)/Nf), B = sqrt(K*c1):
      !>   sum term:  A^2*(z1^2 + chi2(Nf-1))  with chi2(Nf-1) = 2*Gamma((Nf-1)/2)
      !>   cross term: 2*A*B*z1
      !>   combined: (A*z1 + B)^2 + A^2*chi2(Nf-1)   <- always non-negative
      Nf = nfreedom
      c1 = exp(-dat%tstep/dat%thermo_damp)
      K_ref = 0.5_wp*real(Nf,wp)*kB*dat%tsoll
      xi = sqrt(K_ref*(1.0_wp-c1)/real(Nf,wp))   ! A = sqrt(K_ref*(1-c1)/Nf)
      call random_gauss(chi2)                         ! z1 ~ N(0,1)
      K_new = (xi*chi2+sqrt(ekin*c1))**2             ! (A*z1 + B)^2, B = sqrt(K*c1)
      if (Nf > 1) then
        call random_gamma(0.5_wp*real(Nf-1,wp),gam)  ! gam ~ Gamma((Nf-1)/2)
        K_new = K_new+xi**2*2.0_wp*gam               ! + A^2 * chi2(Nf-1)
      end if
      scal = sqrt(K_new/max(ekin,1.0e-30_wp))
    case default
      !>-- (no scaling; langevin uses a separate integration path)
      scal = 1.0_wp
    end select

    return
  end subroutine thermostating

!========================================================================================!
! subroutine langevin_step
! BBK Langevin velocity update for leapfrog MD
  subroutine langevin_step(nat,dat,mass,velo,acc,tstep_au,vel)
    !*************************************************************
    !* Langevin BBK velocity update for leapfrog MD.             *
    !* velo = v(t-dt/2), acc = a(t), vel = v(t+dt/2) on output.  *
    !* Uses dat%thermo_damp as relaxation time tau (in fs).      *
    !*************************************************************
    implicit none
    integer,intent(in)      :: nat
    type(mddata),intent(in) :: dat
    real(wp),intent(in)     :: mass(nat)
    real(wp),intent(in)     :: velo(3,nat),acc(3,nat)
    real(wp),intent(in)     :: tstep_au
    real(wp),intent(out)    :: vel(3,nat)
    real(wp) :: tau_au,gamma_dt,c1,c2,sigma
    real(wp) :: xi
    integer  :: i,k

    tau_au = dat%thermo_damp*fstoau
    gamma_dt = tstep_au/tau_au
    c1 = (1.0_wp-0.5_wp*gamma_dt)/(1.0_wp+0.5_wp*gamma_dt)
    c2 = tstep_au/(1.0_wp+0.5_wp*gamma_dt)

    do i = 1,nat
      sigma = sqrt(2.0_wp*kB*dat%tsoll*tstep_au/(mass(i)*tau_au)) &
            & /(1.0_wp+0.5_wp*gamma_dt)
      do k = 1,3
        call random_gauss(xi)
        vel(k,i) = c1*velo(k,i)+c2*acc(k,i)+sigma*xi
      end do
    end do
    return
  end subroutine langevin_step

  subroutine thermostatprint(dat,pr)
    implicit none
    type(mddata) :: dat
    logical,intent(in) :: pr
    integer :: i,j,k,l,ich,och,io

    if (.not.pr) return
    write (stdout,'("  thermostat",t25,":")',advance='no')
    if (dat%thermostat) then
      select case (trim(dat%thermotype))
      case ('berendsen')
        write (stdout,'(1x,a  )') trim(thermostattype(2))
      case ('bussi','csvr')
        write (stdout,'(1x,a )') trim(thermostattype(4))
      case ('langevin','bbk')
        write (stdout,'(1x,a  )') trim(thermostattype(3))
      case default
        write (stdout,'(1x,a  )') trim(thermostattype(1))
      end select
    else
      write (stdout,'(1x,a  )') trim(thermostattype(1))
    end if

    return
  end subroutine thermostatprint

  subroutine thermostat2int(dat)
    implicit none
    type(mddata) :: dat
    integer :: i,j,k,l,ich,och,io

    if (dat%thermostat) then
      select case (trim(dat%thermotype))
      case ('berendsen')
        dat%thermotype_i = 2
      case ('bussi','csvr')
        dat%thermotype_i = 4
      case ('langevin','bbk')
        dat%thermotype_i = 3
      case default
        dat%thermotype_i = 1
        dat%thermostat = .false.
      end select
    else
      dat%thermotype_i = 1
      dat%thermostat = .false.
    end if

    return
  end subroutine thermostat2int

!========================================================================================!
! subroutine zeroz
! remove z-directional acceleration of selected atoms
  subroutine zeroz(nat,acc,apply)
    implicit none
    integer,intent(in) :: nat
    real(wp),intent(inout) :: acc(3,nat)
    logical,intent(in) :: apply(nat)
    integer :: i
    do i = 1,nat
      if (apply(i)) acc(3,i) = 0.0_wp
    end do
    return
  end subroutine zeroz

!========================================================================================!
! subroutine rmrottr
! some MATHs to remove translational and rotational velocities
  subroutine rmrottr(natoms,mass,vel_atom,c,pbc)
    !***********************************************************
    !* Remove the net translational (and, in free space, the
    !* net rotational) component of the atomic velocities.
    !* Under periodic boundary conditions (pbc=.true.) the cell
    !* is fixed, so rigid-body rotation is not a true symmetry;
    !* only the translational drift is projected out then.
    !***********************************************************
    implicit none
    integer :: natoms
    real(wp) :: vel_atom(3,natoms),c(3,natoms),mass(natoms)
    logical,intent(in),optional :: pbc !> periodic cell -> translation only
    real(wp) :: rlm(3),ram(3),omega(3)
    real(wp) :: ixx,iyy,izz,ixy,ixz,iyz,dummy
    real(wp) :: fixx,fiyy,fizz,fixy,fixz,fiyz,COM(3)
    real(wp) :: inertia(3,3),angmom(3)
    real(wp) :: tmass
    integer :: i
    logical :: rmrot

    !> under PBC remove translation only (rotation is invalid vs. a fixed cell)
    rmrot = .true.
    if (present(pbc)) rmrot = .not.pbc

    rlm = 0.0
    ram = 0.0
    omega = 0.0
    call centerofmass(natoms,c,mass,tmass,COM)
    angmom = 0.0
    if (rmrot) then
      do i = 1,natoms
        c(1,i) = c(1,i)-COM(1)
        c(2,i) = c(2,i)-COM(2)
        c(3,i) = c(3,i)-COM(3)
        angmom(1) = angmom(1)+mass(i)*(c(2,i)*vel_atom(3,i)-&
           &                                      c(3,i)*vel_atom(2,i))
        angmom(2) = angmom(2)+mass(i)*(c(3,i)*vel_atom(1,i)-&
           &                                      c(1,i)*vel_atom(3,i))
        angmom(3) = angmom(3)+mass(i)*(c(1,i)*vel_atom(2,i)-&
           &                                      c(2,i)*vel_atom(1,i))
      end do
      ixx = 0.0
      iyy = 0.0
      izz = 0.0
      ixy = 0.0
      ixz = 0.0
      iyz = 0.0
      do i = 1,natoms
        ixx = ixx+mass(i)*(c(2,i)*c(2,i)+c(3,i)*c(3,i))
        iyy = iyy+mass(i)*(c(3,i)*c(3,i)+c(1,i)*c(1,i))
        izz = izz+mass(i)*(c(1,i)*c(1,i)+c(2,i)*c(2,i))
        ixy = ixy-mass(i)*c(1,i)*c(2,i)
        ixz = ixz-mass(i)*c(1,i)*c(3,i)
        iyz = iyz-mass(i)*c(2,i)*c(3,i)
      end do
      inertia(1,1) = ixx
      inertia(2,2) = iyy
      inertia(3,3) = izz
      inertia(1,2) = ixy
      inertia(2,1) = ixy
      inertia(1,3) = ixz
      inertia(3,1) = ixz
      inertia(2,3) = iyz
      inertia(3,2) = iyz
      call dmatinv(inertia,3,3,dummy)
      omega = matmul(inertia,angmom)
    end if
    do i = 1,natoms
      rlm(1) = rlm(1)+mass(i)*vel_atom(1,i)
      rlm(2) = rlm(2)+mass(i)*vel_atom(2,i)
      rlm(3) = rlm(3)+mass(i)*vel_atom(3,i)
    end do
    do i = 1,natoms
      ram(1) = (omega(2)*c(3,i)-omega(3)*c(2,i))
      ram(2) = (omega(3)*c(1,i)-omega(1)*c(3,i))
      ram(3) = (omega(1)*c(2,i)-omega(2)*c(1,i))

      vel_atom(1,i) = vel_atom(1,i)-rlm(1)/tmass-ram(1)
      vel_atom(2,i) = vel_atom(2,i)-rlm(2)/tmass-ram(2)
      vel_atom(3,i) = vel_atom(3,i)-rlm(3)/tmass-ram(3)

      !> restore the COM offset only if we centered for the rotation part
      if (rmrot) then
        c(1,i) = c(1,i)+COM(1)
        c(2,i) = c(2,i)+COM(2)
        c(3,i) = c(3,i)+COM(3)
      end if
    end do

  contains
    subroutine centerofmass(natoms,c,mass,totmass,COM)
      implicit none
      integer natoms
      real(wp) :: c(3,natoms),totmass,mass(natoms),COM(3)
      integer :: i,j

      COM(1) = 0.0
      COM(2) = 0.0
      COM(3) = 0.0

      totmass = 0
      do i = 1,natoms
        totmass = totmass+mass(i)
        COM(1) = COM(1)+mass(i)*c(1,i)
        COM(2) = COM(2)+mass(i)*c(2,i)
        COM(3) = COM(3)+mass(i)*c(3,i)
      end do

      COM(1) = COM(1)/totmass
      COM(2) = COM(2)/totmass
      COM(3) = COM(3)/totmass
      return
    end subroutine centerofmass

    subroutine dmatinv(a,ldm,n,d)
      implicit none
      integer,intent(in) :: ldm,n
      real(wp),intent(out) :: d
      real(wp),intent(inout) :: a(ldm,*)
      integer :: i,j,k,l(n),m(n)
      real(wp) :: biga,temp
      real(wp),parameter :: tol = 1.0d-12
      !
      d = 1.0_wp
      !
      do k = 1,n
        l(k) = k
        m(k) = k
        biga = a(k,k)
        do j = k,n
          do i = k,n
            if (abs(biga) .lt. abs(a(j,i))) then
              biga = a(j,i)
              l(k) = i
              m(k) = j
            end if
          end do
        end do
        j = l(k)
        if (j .gt. k) then
          do i = 1,n
            temp = -a(i,k)
            a(i,k) = a(i,j)
            a(i,j) = temp
          end do
        end if
        i = m(k)
        if (i .gt. k) then
          do j = 1,n
            temp = -a(k,j)
            a(k,j) = a(i,j)
            a(i,j) = temp
          end do
        end if
        if (abs(biga) .lt. tol) then
          d = 0.0_wp
          return
        end if
        do i = 1,n
          if (i .ne. k) a(k,i) = a(k,i)/(-biga)
        end do
        do i = 1,n
          do j = 1,n
            if (i .ne. k) then
              if (j .ne. k) a(j,i) = a(k,i)*a(j,k)+a(j,i)
            end if
          end do
        end do
        do j = 1,n
          if (j .ne. k) a(j,k) = a(j,k)/biga
        end do
        d = max(-1.0d25,min(1.0d25,d))
        d = d*biga
        a(k,k) = 1.0_wp/biga
      end do
      !
      k = n
      do
        !
        k = k-1
        if (k .le. 0) exit
        i = l(k)
        if (i .gt. k) then
          do j = 1,n
            temp = a(k,j)
            a(k,j) = -a(i,j)
            a(i,j) = temp
          end do
        end if
        j = m(k)
        if (j .gt. k) then
          do i = 1,n
            temp = a(i,k)
            a(i,k) = -a(i,j)
            a(i,j) = temp
          end do
        end if
      end do
      !
      return
    end subroutine dmatinv
  end subroutine rmrottr

!========================================================================================!
  subroutine md_init_mtd(mol,dat,pr)
    implicit none
    type(coord) :: mol
    type(mddata) :: dat
    logical :: pr
    integer :: i

    if (dat%simtype .ne. type_mtd) return
    if (dat%npot < 1) return

    do i = 1,dat%npot
      call mtd_ini(mol,dat%mtd(i),dat%tstep,dat%length_ps,pr)
    end do

    return
  end subroutine md_init_mtd

!========================================================================================!
  subroutine md_add_mtd(self,mtd)
    implicit none
    class(mddata) :: self
    type(mtdpot) :: mtd
    type(mtdpot),allocatable :: mtdtmp(:)
    integer,allocatable :: cvtmp(:)
    integer :: i,j,k

    k = self%npot+1
    allocate (mtdtmp(k))
    allocate (cvtmp(k))
    if (k > 1) then
      do i = 1,k-1
        mtdtmp(i) = self%mtd(i)
        cvtmp(i) = self%cvtype(i)
      end do
    end if
    mtdtmp(k) = mtd
    cvtmp(k) = mtd%mtdtype
    call move_alloc(mtdtmp,self%mtd)
    call move_alloc(cvtmp,self%cvtype)
    self%npot = k

    !> if a metadynamics potential was added for the first time
    !> change MD runtype to MTD accordingly
    if (self%simtype == type_md) then
      self%simtype = type_mtd
    end if

    return
  end subroutine md_add_mtd

!========================================================================================!
  subroutine md_update_mtd(mol,dat,calc,pr)
!*********************************************
!* Update the collective variables for each
!* metadynamics potential saved in mddata
!*********************************************
    implicit none
    type(coord) :: mol
    type(mddata) :: dat
    type(calcdata) :: calc
    logical :: pr
    integer :: i

    if (dat%simtype .ne. type_mtd) return
    if (dat%npot < 1) return

    do i = 1,dat%npot
      select case (dat%cvtype(i))
      case (cv_rmsd,cv_rmsd_static)
        call cv_dump(mol,dat%mtd(i),0.0_wp,pr)
      case default
        cycle
      end select
    end do

    return
  end subroutine md_update_mtd

!========================================================================================!
  subroutine md_calc_mtd(mol,dat,epot,grd,pr)
!***********************************************
!* Calculate energy and gradient contributions
!* from metadynamics potentials and add them to
!* the current total energy and gradient
!***********************************************
!$  use omp_lib
    implicit none
    type(coord) :: mol
    type(mddata) :: dat
    real(wp),intent(inout) :: epot
    real(wp),intent(inout) :: grd(3,mol%nat)
    logical :: pr
    integer :: i
    real(wp) :: emtd
    real(wp),allocatable :: grdmtd(:,:)

    if (dat%simtype .ne. type_mtd) return
    if (dat%npot < 1) return

    !$omp critical
    allocate (grdmtd(3,mol%nat),source=0.0_wp)
    !$omp end critical
    do i = 1,dat%npot
      call calc_mtd(mol,dat%mtd(i),emtd,grdmtd)
      epot = epot+emtd
      grd = grd+grdmtd
    end do
    !$omp critical
    deallocate (grdmtd)
    !$omp end critical

    return
  end subroutine md_calc_mtd

!========================================================================================!

  subroutine md_defaults_fallback(self)
!********************************************
!* Check if selected parameters that are
!* necessary to be able to run an MD
!* have been set, and restore default
!* values if not so
!********************************************
    implicit none
    class(mddata) :: self
    real(wp) :: dum

    if (self%length_ps <= 0.0_wp) then
      !> total runtime in ps
      self%length_ps = 20.0_wp
    end if
    if (self%tstep <= 0.0_wp) then
      !> time step in fs
      self%tstep = 1.0_wp
    end if
    if (self%length_steps <= 0) then
      !> simulation steps
      self%length_steps = nint(self%length_ps*1000.0_wp/self%tstep)
    end if
    if (self%tsoll <= 0.0_wp) then
      !> target temperature
      self%tsoll = 298.15_wp
    end if

    if (self%dumpstep <= 0.0_wp) then
      !> dump frequency in fs
      self%dumpstep = 1000.0_wp
    end if
    if (self%sdump <= 0) then
      !> trajectory structure dump every x steps
      dum = max(1.0_wp, (self%dumpstep/self%tstep))
      self%sdump = nint(dum)
    end if

    if (self%md_hmass <= 0.0_wp) then
      !> hydrogen mass
      self%md_hmass = 1.00794075_wp
    end if

    if (self%shake) then
      !> SHAKE, if turned on but no mode selected
      if (self%shk%shake_mode == 0) then
        self%shk%shake_mode = 2 !> all bonds
      end if
    end if

    !> block length (for average analysis and restart dump)
    if (self%blockl <= 0) then
      self%blockl = min(5000,idint(5000.0_wp/self%tstep))
    end if
    self%maxblock = nint(self%length_steps/float(self%blockl))

    call thermostat2int(self)

  end subroutine md_defaults_fallback
!========================================================================================!
!========================================================================================!
end module dynamics_module
