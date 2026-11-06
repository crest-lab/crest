module test_ddx
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use crest_parameters
  use strucrd
  use crest_testmol
#ifdef WITH_DDX
  use crest_ddx_pc,only:ddx_pc_engrad
  use crest_electrostatic,only:electrostatic_core
  use crest_surface,only:surface_engrad
  use crest_solvation,only:solvation_data,solvation_setup,solvation_core
#endif
  implicit none
  private

  public :: collect_ddx

  real(wp),parameter :: eps_water = 78.4_wp

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for the standalone ddX point-charge solvation engine
!========================================================================================!
!========================================================================================!

  subroutine collect_ddx(testsuite)
    type(unittest_type),allocatable,intent(out) :: testsuite(:)
!&<
    testsuite = [ &
#ifdef WITH_DDX
    new_unittest("ddX point-charge energy < 0  ",test_ddx_energy), &
    new_unittest("ddX explicit grad vs finite-diff",test_ddx_fd), &
    new_unittest("EEQ-BC grad vs finite-diff   ",test_eeq_fd), &
    new_unittest("SASA nonpolar grad vs FD     ",test_surface_fd), &
    new_unittest("ddX full grad (chain) vs FD  ",test_solv_full_fd), &
    new_unittest("solvation composite grad vs FD",test_composite_fd), &
    new_unittest("GFN2/ALPB params load + run  ",test_gfn2_params) &
#else
    new_unittest("ddX not compiled",test_ddx_nocompile,should_fail=.true.) &
#endif
    ]
!&>
  end subroutine collect_ddx

#ifndef WITH_DDX
  subroutine test_ddx_nocompile(error)
    type(error_type),allocatable,intent(out) :: error
    allocate (error)
  end subroutine test_ddx_nocompile
#endif

#ifdef WITH_DDX
!> Deterministic, charge-neutral set of fabricated point charges
  subroutine fake_charges(mol,q)
    type(coord),intent(in) :: mol
    real(wp),allocatable,intent(out) :: q(:)
    integer :: i
    allocate (q(mol%nat))
    do i = 1,mol%nat
      q(i) = 0.1_wp*sin(real(i,wp))
    end do
    q(:) = q(:)-sum(q)/real(mol%nat,wp)   !> enforce neutrality
  end subroutine fake_charges

  subroutine test_ddx_energy(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    real(wp),allocatable :: q(:),grad(:,:)
    real(wp) :: energy
    integer :: io
    call get_testmol('cytosine',mol)
    call fake_charges(mol,q)
    allocate (grad(3,mol%nat))
    call ddx_pc_engrad(mol,q,'cpcm',eps_water,energy,grad,io)
    call check(error,io,0)
    if (allocated(error)) return
    !> polarization of a charge distribution always lowers the energy
    call check(error,energy < 0.0_wp)
  end subroutine test_ddx_energy

  subroutine test_ddx_fd(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    real(wp),allocatable :: q(:),grad(:,:)
    real(wp) :: energy,el,er,h,dev,rms
    integer :: io,iat,k,n
    real(wp),parameter :: thr = 1.0e-5_wp
    h = 5.0e-4_wp
    call get_testmol('cytosine',mol)
    call fake_charges(mol,q)
    allocate (grad(3,mol%nat))
    !> analytic explicit gradient (charges held fixed -> no dqdr term)
    call ddx_pc_engrad(mol,q,'cpcm',eps_water,energy,grad,io)
    call check(error,io,0)
    if (allocated(error)) return
    !> central finite differences at fixed charges
    rms = 0.0_wp; n = 0
    do iat = 1,mol%nat
      do k = 1,3
        mol%xyz(k,iat) = mol%xyz(k,iat)+h
        call ddx_e(mol,q,er,io)
        mol%xyz(k,iat) = mol%xyz(k,iat)-2.0_wp*h
        call ddx_e(mol,q,el,io)
        mol%xyz(k,iat) = mol%xyz(k,iat)+h
        dev = (er-el)/(2.0_wp*h)-grad(k,iat)
        rms = rms+dev*dev; n = n+1
      end do
    end do
    rms = sqrt(rms/real(n,wp))
    write (*,'("       ... ddX grad RMS dev =",es12.4)') rms
    call check(error,rms < thr)
  end subroutine test_ddx_fd

!> energy-only helper for the finite-difference loop
  subroutine ddx_e(mol,q,energy,io)
    type(coord),intent(in) :: mol
    real(wp),intent(in)    :: q(:)
    real(wp),intent(out)   :: energy
    integer,intent(out)    :: io
    real(wp),allocatable   :: g(:,:)
    allocate (g(3,mol%nat))
    call ddx_pc_engrad(mol,q,'cpcm',eps_water,energy,g,io)
  end subroutine ddx_e

!> EEQ-BC electrostatic gradient vs central finite differences
  subroutine test_eeq_fd(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    real(wp),allocatable :: grad(:,:),qat(:),g(:,:)
    real(wp) :: energy,el,er,h,dev,rms
    integer :: io,iat,k,n
    real(wp),parameter :: thr = 1.0e-6_wp
    h = 1.0e-4_wp
    call get_testmol('cytosine',mol)
    allocate (grad(3,mol%nat),qat(mol%nat),g(3,mol%nat))
    call electrostatic_core(mol,0,'eeqbc',energy,grad,qat,io)
    call check(error,io,0)
    if (allocated(error)) return
    rms = 0.0_wp; n = 0
    do iat = 1,mol%nat
      do k = 1,3
        mol%xyz(k,iat) = mol%xyz(k,iat)+h
        call electrostatic_core(mol,0,'eeqbc',er,g,qat,io)
        mol%xyz(k,iat) = mol%xyz(k,iat)-2.0_wp*h
        call electrostatic_core(mol,0,'eeqbc',el,g,qat,io)
        mol%xyz(k,iat) = mol%xyz(k,iat)+h
        dev = (er-el)/(2.0_wp*h)-grad(k,iat)
        rms = rms+dev*dev; n = n+1
      end do
    end do
    rms = sqrt(rms/real(n,wp))
    write (*,'("       ... EEQ-BC grad RMS dev =",es12.4)') rms
    call check(error,rms < thr)
  end subroutine test_eeq_fd

!> Nonpolar SASA gradient vs central finite differences
  subroutine test_surface_fd(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    real(wp),allocatable :: grad(:,:),tension(:),sasa(:),g(:,:)
    real(wp) :: energy,el,er,h,dev,rms
    integer :: io,iat,k,n
    !> SASA gradients are limited by the Lebedev grid + neighbour-list cutoff,
    !> so the finite-difference consistency floors out around 1e-5
    real(wp),parameter :: thr = 1.0e-4_wp
    h = 1.0e-4_wp
    call get_testmol('cytosine',mol)
    allocate (grad(3,mol%nat),tension(mol%nat),sasa(mol%nat),g(3,mol%nat))
    tension(:) = 0.01_wp   !> arbitrary uniform surface tension
    call surface_engrad(mol,tension,energy,grad,io,sasa=sasa)
    call check(error,io,0)
    if (allocated(error)) return
    rms = 0.0_wp; n = 0
    do iat = 1,mol%nat
      do k = 1,3
        mol%xyz(k,iat) = mol%xyz(k,iat)+h
        call surface_engrad(mol,tension,er,g,io)
        mol%xyz(k,iat) = mol%xyz(k,iat)-2.0_wp*h
        call surface_engrad(mol,tension,el,g,io)
        mol%xyz(k,iat) = mol%xyz(k,iat)+h
        dev = (er-el)/(2.0_wp*h)-grad(k,iat)
        rms = rms+dev*dev; n = n+1
      end do
    end do
    rms = sqrt(rms/real(n,wp))
    write (*,'("       ... SASA grad RMS dev =",es12.4)') rms
    call check(error,rms < thr)
  end subroutine test_surface_fd

!> Full ddX gradient (explicit + dq/dR chain term) vs finite differences with
!> geometry-dependent EEQ-BC charges -- validates the charge-response term.
  subroutine test_solv_full_fd(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    real(wp),allocatable :: grad(:,:),geeq(:,:),qat(:),dqdr(:,:,:)
    real(wp) :: energy,eeq,el,er,h,dev,rms
    integer :: io,iat,k,n
    real(wp),parameter :: thr = 1.0e-5_wp
    h = 5.0e-4_wp
    call get_testmol('cytosine',mol)
    allocate (grad(3,mol%nat),geeq(3,mol%nat),qat(mol%nat),dqdr(3,mol%nat,mol%nat))
    !> analytic full gradient: charges + dq/dR from EEQ-BC, fed into ddX
    call electrostatic_core(mol,0,'eeqbc',eeq,geeq,qat,io,dqdr=dqdr)
    call ddx_pc_engrad(mol,qat,'cpcm',eps_water,energy,grad,io,dqdr=dqdr)
    call check(error,io,0)
    if (allocated(error)) return
    !> finite differences with charges recomputed at each perturbed geometry
    rms = 0.0_wp; n = 0
    do iat = 1,mol%nat
      do k = 1,3
        mol%xyz(k,iat) = mol%xyz(k,iat)+h
        call solv_e(mol,er,io)
        mol%xyz(k,iat) = mol%xyz(k,iat)-2.0_wp*h
        call solv_e(mol,el,io)
        mol%xyz(k,iat) = mol%xyz(k,iat)+h
        dev = (er-el)/(2.0_wp*h)-grad(k,iat)
        rms = rms+dev*dev; n = n+1
      end do
    end do
    rms = sqrt(rms/real(n,wp))
    write (*,'("       ... ddX full grad RMS dev =",es12.4)') rms
    call check(error,rms < thr)
  end subroutine test_solv_full_fd

!> ddX solvation energy with self-consistently recomputed EEQ-BC charges
  subroutine solv_e(mol,energy,io)
    type(coord),intent(in) :: mol
    real(wp),intent(out)   :: energy
    integer,intent(out)    :: io
    real(wp),allocatable   :: q(:),g(:,:)
    real(wp) :: eeq
    allocate (q(mol%nat),g(3,mol%nat))
    call electrostatic_core(mol,0,'eeqbc',eeq,g,q,io)
    call ddx_e(mol,q,energy,io)
  end subroutine solv_e

!> Full solvation composite (polar + tension + hbond) gradient vs finite
!> differences -- the end-to-end check with both charges and SASA varying
  subroutine test_composite_fd(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    type(solvation_data) :: solv
    real(wp),allocatable :: grad(:,:),g(:,:)
    real(wp) :: energy,el,er,h,dev,rms
    integer :: io,iat,k,n
    real(wp),parameter :: thr = 1.0e-4_wp   !> floored by the SASA integrator
    h = 5.0e-4_wp
    call get_testmol('cytosine',mol)
    allocate (grad(3,mol%nat),g(3,mol%nat))
    !> inject amplified dummy CDS parameters to stress the nonpolar gradient
    solv%charge_model = 'eeqbc'; solv%smodel = 'cpcm'; solv%do_hbond = .true.
    solv%eps = eps_water; solv%probe = 1.0_wp*aatoau; solv%loaded = .true.
    allocate (solv%tension(mol%nat),source=0.01_wp)
    allocate (solv%hbond(mol%nat),source=0.005_wp)
    !> analytic gradient of the full composite
    call solvation_core(mol,0,solv,energy,grad,io)
    call check(error,io,0)
    if (allocated(error)) return
    !> finite differences (charges + SASA recomputed at each geometry)
    rms = 0.0_wp; n = 0
    do iat = 1,mol%nat
      do k = 1,3
        mol%xyz(k,iat) = mol%xyz(k,iat)+h
        call solvation_core(mol,0,solv,er,g,io)
        mol%xyz(k,iat) = mol%xyz(k,iat)-2.0_wp*h
        call solvation_core(mol,0,solv,el,g,io)
        mol%xyz(k,iat) = mol%xyz(k,iat)+h
        dev = (er-el)/(2.0_wp*h)-grad(k,iat)
        rms = rms+dev*dev; n = n+1
      end do
    end do
    rms = sqrt(rms/real(n,wp))
    write (*,'("       ... composite grad RMS dev =",es12.4)') rms
    call check(error,rms < thr)
  end subroutine test_composite_fd

!> Load real GFN2/ALPB water parameters via the data object and run the
!> composite, checking the gradient against finite differences end to end
  subroutine test_gfn2_params(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    type(solvation_data) :: solv
    real(wp),allocatable :: grad(:,:),g(:,:)
    real(wp) :: energy,el,er,h,dev,rms
    integer :: io,iat,k,n
    real(wp),parameter :: thr = 1.0e-4_wp
    h = 5.0e-4_wp
    call get_testmol('cytosine',mol)
    allocate (grad(3,mol%nat),g(3,mol%nat))
    solv%charge_model = 'eeqbc'; solv%smodel = 'cpcm'
    solv%solvent = 'water'; solv%do_hbond = .true.
    call solvation_setup(mol,solv,io)
    call check(error,io,0)
    if (allocated(error)) return
    call solvation_core(mol,0,solv,energy,grad,io)
    call check(error,io,0)
    if (allocated(error)) return
    call check(error,energy < 0.0_wp)   !> aqueous solvation should stabilize
    if (allocated(error)) return
    rms = 0.0_wp; n = 0
    do iat = 1,mol%nat
      do k = 1,3
        mol%xyz(k,iat) = mol%xyz(k,iat)+h
        call solvation_core(mol,0,solv,er,g,io)
        mol%xyz(k,iat) = mol%xyz(k,iat)-2.0_wp*h
        call solvation_core(mol,0,solv,el,g,io)
        mol%xyz(k,iat) = mol%xyz(k,iat)+h
        dev = (er-el)/(2.0_wp*h)-grad(k,iat)
        rms = rms+dev*dev; n = n+1
      end do
    end do
    rms = sqrt(rms/real(n,wp))
    write (*,'("       ... GFN2 composite E =",f12.6," grad RMS =",es12.4)') energy,rms
    call check(error,rms < thr)
  end subroutine test_gfn2_params
#endif

!========================================================================================!
end module test_ddx
