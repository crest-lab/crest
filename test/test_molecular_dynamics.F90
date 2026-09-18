module test_molecular_dynamics
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use crest_parameters
  use crest_calculator
  use strucrd
  use crest_testmol
  use dynamics_module
  use iomod,only:remove
  implicit none
  private

  public :: collect_mol_dynamics

  real(wp),parameter :: thr = 5e+6_wp*epsilon(1.0_wp)

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for using molecular dynamics routines in CREST
!========================================================================================!
!========================================================================================!

!> Collect all exported unit tests
  subroutine collect_mol_dynamics(testsuite)
    !> Collection of tests
    type(unittest_type),allocatable,intent(out) :: testsuite(:)

!&<
    testsuite = [ &
#ifdef WITH_GFNFF
    new_unittest("Compiled gfnff subproject     ",test_compiled_gfnff), &
    new_unittest("molecular dynamics (SHAKE off)",test_md_shake_off), &
    new_unittest("molecular dynamics (SHAKE on) ",test_md_shake_on), &
    new_unittest("molecular dynamics (SHAKE H)  ",test_md_shake_honly), &
    new_unittest("thermostat: berendsen         ",test_md_thermostat_berendsen), &
    new_unittest("thermostat: bussi (CSVR)      ",test_md_thermostat_bussi), &
    new_unittest("thermostat: langevin (BBK)    ",test_md_thermostat_langevin) &
#else
    new_unittest("Compiled gfnff subproject",test_compiled_gfnff,should_fail=.true.) &
#endif
    ]
!&>

  end subroutine collect_mol_dynamics

!========================================================================================!

  subroutine test_compiled_gfnff(error)
    type(error_type),allocatable,intent(out) :: error
#ifndef WITH_GFNFF
    write (*,'("       ...")') 'gfnff not compiled, expecting fail.'
    allocate (error)
#endif
  end subroutine test_compiled_gfnff

!========================================================================================!
!  The three MD tests below intentionally only set up the shared infrastructure (calcdata
!  and a test molecule) and provide placeholders for the MD-specific calls/checks.
!  Fill the marked sections with your MD driver + assertions.
!========================================================================================!

  subroutine test_md_shake_off(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    type(mddata) :: mdyn
    integer :: io
    logical :: pr

    !> setup calculator backend
    call sett%create('gfnff')
    call calc%add(sett)

    !> test molecule
    call get_testmol('methane',mol)

    !> MD setup
    pr = .false.
    io = 0
    mdyn%length_ps = 200.0_wp
    call mdyn%defaults()
    mdyn%shake = .false.
    mdyn%restart = .true. !> turn on restart reading (for determinic results)
    mdyn%wrtrj = .false. !> turn off trajectory dump
    call write_fake_restart(mol,mdyn%restartfile)

    !> run
    call dynamics(mol,mdyn,calc,pr,io)

    !> cleanup
    call remove(mdyn%restartfile)
    call remove('crest_0.mdrestart')

    !> checks
    call check(error,io,0)
    if (allocated(error)) return
    !> Average temperature must be within ±50 K of thermostat target (compiler-portable)
    call check(error,mdyn%Tavg,mdyn%tsoll,thr=50.0_wp)
    if (allocated(error)) return
  end subroutine test_md_shake_off

!========================================================================================!

  subroutine test_md_shake_on(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    type(mddata) :: mdyn
    integer :: io,i,ia,ib
    logical :: wr,pr
    real(wp) :: d

    !> setup calculator backend
    call sett%create('gfnff')
    call calc%add(sett)

    !> get test molecule
    call get_testmol('methane',mol)

    !> MD setup
    pr = .false.
    io = 0
    mdyn%length_ps = 50.0_wp
    mdyn%Tsoll = 450.0_wp
    call mdyn%defaults()
    mdyn%shake = .true.
    mdyn%restart = .true. !> turn on restart reading (for determinic results)
    mdyn%wrtrj = .false. !> turn off trajectory dump
    call write_fake_restart(mol,mdyn%restartfile)

    !> run
    call dynamics(mol,mdyn,calc,pr,io)

    !> cleanup
    call remove(mdyn%restartfile)
    call remove('crest_0.mdrestart')

    !> checks
    call check(error,io,0)
    if (allocated(error)) return
    !> Average temperature must be within ±50 K of thermostat target (compiler-portable)
    call check(error,mdyn%Tavg,mdyn%tsoll,thr=50.0_wp)
    if (allocated(error)) return
    !> SHAKE: all constrained bonds must satisfy their target lengths
    do i = 1,mdyn%shk%ncons
      ia = mdyn%shk%conslist(1,i)
      ib = mdyn%shk%conslist(2,i)
      d = norm2(mol%xyz(:,ia)-mol%xyz(:,ib))
      call check(error,d**2,mdyn%shk%distcons(i),thr=1e-4_wp)
      if (allocated(error)) return
    end do
  end subroutine test_md_shake_on

!========================================================================================!

  subroutine test_md_shake_honly(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol

    type(mddata) :: mdyn
    integer :: io,i,ia,ib
    logical :: wr,pr
    real(wp) :: d

    !> setup calculator backend
    call sett%create('gfnff')
    call calc%add(sett)

    !> get test molecule
    call get_testmol('caffeine',mol)

    !> MD setup
    pr = .false.
    io = 0
    mdyn%length_ps = 10.0_wp  !> shorter runtime because the mol is larger
    call mdyn%defaults()
    mdyn%shake = .true.
    mdyn%shk%shake_mode = 1
    mdyn%restart = .true. !> turn on restart reading (for determinic results)
    mdyn%wrtrj = .false. !> turn off trajectory dump
    call write_fake_restart(mol,mdyn%restartfile)

    !> run
    call dynamics(mol,mdyn,calc,pr,io)

    !> cleanup
    call remove(mdyn%restartfile)
    call remove('crest_0.mdrestart')

    !> checks
    call check(error,io,0)
    if (allocated(error)) return
    !> Average temperature must be within ±50 K of thermostat target (compiler-portable)
    call check(error,mdyn%Tavg,mdyn%tsoll,thr=50.0_wp)
    if (allocated(error)) return
    !> SHAKE: all constrained bonds must satisfy their target lengths
    do i = 1,mdyn%shk%ncons
      ia = mdyn%shk%conslist(1,i)
      ib = mdyn%shk%conslist(2,i)
      d = norm2(mol%xyz(:,ia)-mol%xyz(:,ib))
      call check(error,d**2,mdyn%shk%distcons(i),thr=1e-4_wp)
      if (allocated(error)) return
    end do
  end subroutine test_md_shake_honly

!========================================================================================!

  subroutine test_md_thermostat_berendsen(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    type(mddata) :: mdyn
    integer :: io
    logical :: pr

    call sett%create('gfnff')
    call calc%add(sett)
    call get_testmol('methane',mol)

    pr = .false.
    io = 0
    mdyn%length_ps = 50.0_wp
    mdyn%tstep = 1.0_wp
    call mdyn%defaults()
    mdyn%shake = .false.
    mdyn%samerand = .true.   !> deterministic RNG seed
    mdyn%thermotype = 'berendsen'
    mdyn%thermo_damp = 500.0_wp
    mdyn%restart = .true.
    mdyn%wrtrj = .false.
    call write_fake_restart(mol,mdyn%restartfile)

    call dynamics(mol,mdyn,calc,pr,io)

    call remove(mdyn%restartfile)
    call remove('crest_0.mdrestart')

    !> Simulation must complete without error
    call check(error,io,0)
    if (allocated(error)) return
    !> Average temperature within ±50 K of target.
    !> Berendsen's distinguishing property is deterministic (non-stochastic) convergence;
    !> with τ=500 fs >> dt, scal≈1 each step so variance stays near canonical — no variance check.
    call check(error,mdyn%Tavg,mdyn%tsoll,thr=50.0_wp)
    if (allocated(error)) return
  end subroutine test_md_thermostat_berendsen

!========================================================================================!

  subroutine test_md_thermostat_bussi(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    type(mddata) :: mdyn
    integer :: io
    logical :: pr
    real(wp) :: tvar

    call sett%create('gfnff')
    call calc%add(sett)
    call get_testmol('caffeine',mol)

    pr = .false.
    io = 0
    mdyn%length_ps = 25.0_wp
    call mdyn%defaults()
    mdyn%shake = .false.
    mdyn%samerand = .true.   !> deterministic RNG seed
    mdyn%thermotype = 'bussi'
    mdyn%thermo_damp = 500.0_wp
    mdyn%restart = .true.
    mdyn%wrtrj = .false.
    call write_fake_restart(mol,mdyn%restartfile,tstart=300.0_wp)

    call dynamics(mol,mdyn,calc,pr,io)

    call remove(mdyn%restartfile)
    call remove('crest_0.mdrestart')

    !> Simulation must complete without error
    call check(error,io,0)
    if (allocated(error)) return
    !> Average temperature within ±50 K of target
    call check(error,mdyn%Tavg,mdyn%tsoll,thr=50.0_wp)
    if (allocated(error)) return
    !> Bussi samples the canonical NVT ensemble: Tvar must be close to the
    !> canonical reference 2*<T>^2/Ndf. We check Tvar > 10% of canonical —
    !> 10x margin below the expected ~100%, and far above Berendsen's ~0.1%.
    tvar = mdyn%Tvar
    if (tvar < 0.1_wp*(2.0_wp*mdyn%Tavg**2/real(mdyn%Ndf,wp))) then
      call test_failed(error,'Bussi thermostat does not show canonical T fluctuations: '// &
        & 'expected Tvar ~ 2*Tavg^2/Ndf')
      return
    end if
  end subroutine test_md_thermostat_bussi

!========================================================================================!

  subroutine test_md_thermostat_langevin(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    type(mddata) :: mdyn
    integer :: io
    logical :: pr
    real(wp) :: tvar

    call sett%create('gfnff')
    call calc%add(sett)
    call get_testmol('caffeine',mol)

    pr = .false.
    io = 0
    mdyn%length_ps = 25.0_wp
    call mdyn%defaults()
    mdyn%shake = .false.
    mdyn%samerand = .true.   !> deterministic RNG seed
    mdyn%thermotype = 'langevin'
    mdyn%thermo_damp = 298.15_wp
    mdyn%restart = .true.
    mdyn%wrtrj = .false.
    call write_fake_restart(mol,mdyn%restartfile)

    call dynamics(mol,mdyn,calc,pr,io)

    call remove(mdyn%restartfile)
    call remove('crest_0.mdrestart')

    !> Simulation must complete without error
    call check(error,io,0)
    if (allocated(error)) return
    !> Average temperature within ±50 K of target (equipartition theorem)
    call check(error,mdyn%Tavg,mdyn%tsoll,thr=50.0_wp)
    if (allocated(error)) return
    !> Langevin (BBK) samples the canonical NVT ensemble via per-atom stochastic
    !> friction: Tvar must be close to the canonical reference 2*<T>^2/Ndf.
    !> We check Tvar > 10% of canonical — same criterion as for Bussi.
    tvar = mdyn%Tvar
    if (tvar < 0.1_wp*(2.0_wp*mdyn%Tavg**2/real(mdyn%Ndf,wp))) then
      call test_failed(error,'Langevin thermostat does not show canonical T fluctuations: '// &
        & 'expected Tvar ~ 2*Tavg^2/Ndf')
      return
    end if
  end subroutine test_md_thermostat_langevin

!========================================================================================!

  subroutine write_fake_restart(mol,restartfile,tstart)
    implicit none
    type(coord),intent(in) :: mol
    character(len=:),allocatable,intent(out) :: restartfile
    real(wp),intent(in),optional :: tstart
    integer :: ich,ii
    restartfile = 'crest_test.mdrestart'
    open (newunit=ich,file=restartfile)
    if (present(tstart)) then
      write (ich,*) tstart
    else
      write (ich,*) 500.0_wp
    end if
    do ii = 1,mol%nat
      write (ich,'(6D22.14)') mol%xyz(1:3,ii),mol%xyz(1:3,ii)*0.00005_wp
    end do
    close (ich)
  end subroutine write_fake_restart

!========================================================================================!
!========================================================================================!
end module test_molecular_dynamics
