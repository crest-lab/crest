module test_oniom
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use crest_parameters
  use crest_calculator
  use oniom_hessian
  use strucrd
  use crest_testmol
  use lwoniom_module
  implicit none
  private

  public :: collect_oniom

  !> consistency of identical code paths (should agree almost exactly)
  real(wp),parameter :: thr_e = 1.0e-8_wp
  !> finite-difference comparisons (step 0.005 bohr, central differences)
  real(wp),parameter :: thr_fd = 1.0e-5_wp
  real(wp),parameter :: fdstep = 0.005_wp

  !> the test system: caffeine, with the methyl group on N9 (atoms 14,22-24)
  !> as the high-level region. The single N9-C14 bond cut generates one
  !> link H atom, so the model system has 5 atoms.
  integer,parameter :: methyl(4) = [14,22,23,24]

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for ONIOM calculations (via the lwONIOM library) in crest
!========================================================================================!
!========================================================================================!

!> Collect all exported unit tests
  subroutine collect_oniom(testsuite)
    !> Collection of tests
    type(unittest_type),allocatable,intent(out) :: testsuite(:)

!&<
    testsuite = [ &
#if defined(WITH_LWONIOM) && defined(WITH_GFNFF) && defined(WITH_GFN0)
    new_unittest("Compiled lwONIOM subproject   ",test_compiled_lwoniom), &
    new_unittest("ONIOM singlepoint consistency ",test_oniom_sp), &
    new_unittest("ONIOM gradient vs. finite diff",test_oniom_grad_fd), &
    new_unittest("ONIOM numerical Hessian       ",test_oniom_hessian) &
#else
    new_unittest("Compiled lwONIOM subproject",test_compiled_lwoniom,should_fail=.true.) &
#endif
    ]
!&>
  end subroutine collect_oniom

  subroutine test_compiled_lwoniom(error)
    type(error_type),allocatable,intent(out) :: error
#if !(defined(WITH_LWONIOM) && defined(WITH_GFNFF) && defined(WITH_GFN0))
    write (*,'("       ...")') 'lwoniom (or gfnff/gfn0) not compiled, expecting fail.'
    allocate (error)
#endif
  end subroutine test_compiled_lwoniom

#if defined(WITH_LWONIOM) && defined(WITH_GFNFF) && defined(WITH_GFN0)

!========================================================================================!

  subroutine setup_oniom_calc(mol,calc)
!**********************************************************************
!* Set up the 2-layer ONIOM test system: caffeine with the N9-methyl
!* group at the GFN0-xTB level ("high"), embedded in the GFN-FF
!* description of the full molecule ("low").
!* Uses the direct lwoniom_initialize() API (no input file), so the
!* theory level ids default to the layer numbers: layer 1 -> GFN-FF
!* (calculator 1), layer 2 -> GFN0-xTB (calculator 2).
!**********************************************************************
    implicit none
    type(coord),intent(out) :: mol
    type(calcdata),intent(out) :: calc
    type(calculation_settings) :: sett
    logical,allocatable :: layer(:,:)

    call get_testmol('caffeine',mol)

    call sett%create('gfnff')
    call calc%add(sett)
    call sett%deallocate()
    call sett%create('gfn0')
    call calc%add(sett)
    call sett%deallocate()

! ── ONIOM partitioning: layer 1 = all atoms, layer 2 = methyl group ───────────
    allocate (layer(mol%nat,2),source=.false.)
    layer(:,1) = .true.
    layer(methyl,2) = .true.
    allocate (calc%ONIOM)
    call lwoniom_initialize(mol%nat,mol%at,mol%xyz,calc%ONIOM,layer)
    call calc%ONIOMexpand()
    deallocate (layer)
  end subroutine setup_oniom_calc

!========================================================================================!

  subroutine test_oniom_sp(error)
!**********************************************************************
!* Full ONIOM singlepoint through the crest calculator, verified
!* against a manual reconstruction: independent single-level
!* calculations on the model geometries from lwoniom_get_jobgeo(),
!* recombined as E = E_low(real) + E_high(model) - E_low(model).
!* Also checks translational invariance of the projected gradient.
!**********************************************************************
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol,fmol
    type(calcdata) :: calc,fcalc
    type(calculation_settings) :: fsett
    type(lwoniom_job),allocatable :: jobs(:)
    real(wp) :: energy,eref,ej
    real(wp),allocatable :: grad(:,:),fgrad(:,:)
    integer,allocatable :: atj(:)
    real(wp),allocatable :: xyzj(:,:)
    integer :: io,j,natj

    call setup_oniom_calc(mol,calc)
    allocate (grad(3,mol%nat),source=0.0_wp)

    call engrad(mol,calc,energy,grad,io)
    call check(error,io,0)
    if (allocated(error)) return

! ── manual ONIOM reconstruction from independent single-level runs ────────────
    call lwoniom_get_jobs(calc%ONIOM,jobs)
    call check(error,size(jobs),3)
    if (allocated(error)) return

    eref = 0.0_wp
    do j = 1,size(jobs)
      call lwoniom_get_jobgeo(calc%ONIOM,jobs(j),natj,atj,xyzj)
      call fmol%deallocate()
      fmol%nat = natj
      fmol%at = atj
      fmol%xyz = xyzj
      call fcalc%reset()
      call fsett%deallocate()
      if (jobs(j)%theoryid == 1) then
        call fsett%create('gfnff')
      else
        call fsett%create('gfn0')
      end if
      call fcalc%add(fsett)
      if (allocated(fgrad)) deallocate (fgrad)
      allocate (fgrad(3,natj),source=0.0_wp)
      call engrad(fmol,fcalc,ej,fgrad,io)
      call check(error,io,0)
      if (allocated(error)) return
      if (jobs(j)%level == oniom_low) then
        eref = eref-ej
      else
        eref = eref+ej
      end if
    end do

    call check(error,energy,eref,thr=thr_e)
    if (allocated(error)) return

! ── projected ONIOM gradient must remain translationally invariant ────────────
    if (norm2(sum(grad,dim=2)) > 1.0e-7_wp) then
      call test_failed(error,"ONIOM gradient violates translational invariance")
      print '(3es21.14)',sum(grad,dim=2)
    end if
  end subroutine test_oniom_sp

!========================================================================================!

  subroutine test_oniom_grad_fd(error)
!**********************************************************************
!* Verify the Jacobian-projected ONIOM gradient against central
!* finite differences of the ONIOM energy, for one atom of the
!* high-level region, the cut-bond host atom (whose displacement
!* also moves the link atom), and one unaffected outer atom.
!**********************************************************************
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    type(calcdata) :: calc
    real(wp) :: energy,ep,em,fd
    real(wp),allocatable :: grad(:,:)
    integer :: io,t
    !> tested components: methyl C (x), cut host N9 (y), carbonyl O8 (z)
    integer,parameter :: tatom(3) = [14,9,8]
    integer,parameter :: tdir(3) = [1,2,3]

    call setup_oniom_calc(mol,calc)
    allocate (grad(3,mol%nat),source=0.0_wp)

    call engrad(mol,calc,energy,grad,io)
    call check(error,io,0)
    if (allocated(error)) return

    do t = 1,3
      mol%xyz(tdir(t),tatom(t)) = mol%xyz(tdir(t),tatom(t))+fdstep
      call engrad(mol,calc,ep,grad,io)
      call check(error,io,0)
      if (allocated(error)) return
      mol%xyz(tdir(t),tatom(t)) = mol%xyz(tdir(t),tatom(t))-2.0_wp*fdstep
      call engrad(mol,calc,em,grad,io)
      call check(error,io,0)
      if (allocated(error)) return
      mol%xyz(tdir(t),tatom(t)) = mol%xyz(tdir(t),tatom(t))+fdstep
      fd = (ep-em)*0.5_wp/fdstep

      !> re-evaluate the analytical gradient at the original geometry
      call engrad(mol,calc,energy,grad,io)
      call check(error,io,0)
      if (allocated(error)) return
      call check(error,grad(tdir(t),tatom(t)),fd,thr=thr_fd)
      if (allocated(error)) then
        print '(a,i0,a,i0)','component: atom ',tatom(t),' dir ',tdir(t)
        print '(2es21.14)',grad(tdir(t),tatom(t)),fd
        return
      end if
    end do
  end subroutine test_oniom_grad_fd

!========================================================================================!

  subroutine test_oniom_hessian(error)
!**********************************************************************
!* ONIOM numerical Hessian: check exact symmetry and verify one
!* column against central finite differences of the analytical
!* (Jacobian-projected) ONIOM gradient. Since the fragment
!* coordinates are exactly linear in the real-system coordinates,
!* sum(J^T H_frag J) and d(g_ONIOM)/dx must agree to FD accuracy.
!**********************************************************************
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    type(calcdata) :: calc
    real(wp) :: energy
    real(wp),allocatable :: hess(:,:),grad(:,:),gp(:,:),gm(:,:),fdcol(:)
    integer :: io,nat3,k
    !> tested column: x component of the methyl C (atom 14)
    integer,parameter :: tatom = 14,tdir = 1

    call setup_oniom_calc(mol,calc)
    nat3 = 3*mol%nat
    allocate (hess(nat3,nat3),source=0.0_wp)

    call ONIOM_calc_hessians(mol,calc,hess)

! ── the reconstructed ONIOM Hessian must be symmetric ─────────────────────────
    if (maxval(abs(hess-transpose(hess))) > 1.0e-10_wp) then
      call test_failed(error,"ONIOM Hessian is not symmetric")
      return
    end if

! ── compare one column against FD of the analytical ONIOM gradient ────────────
    allocate (grad(3,mol%nat),gp(3,mol%nat),gm(3,mol%nat),source=0.0_wp)
    k = (tatom-1)*3+tdir

    mol%xyz(tdir,tatom) = mol%xyz(tdir,tatom)+fdstep
    call engrad(mol,calc,energy,gp,io)
    call check(error,io,0)
    if (allocated(error)) return
    mol%xyz(tdir,tatom) = mol%xyz(tdir,tatom)-2.0_wp*fdstep
    call engrad(mol,calc,energy,gm,io)
    call check(error,io,0)
    if (allocated(error)) return
    mol%xyz(tdir,tatom) = mol%xyz(tdir,tatom)+fdstep

    fdcol = reshape((gp-gm)*0.5_wp/fdstep, [nat3])
    if (maxval(abs(hess(:,k)-fdcol)) > 5.0e-5_wp) then
      call test_failed(error,"ONIOM Hessian column does not match FD gradient")
      print '(a,es21.14)','max deviation: ',maxval(abs(hess(:,k)-fdcol))
      return
    end if
  end subroutine test_oniom_hessian

#endif

!========================================================================================!
!========================================================================================!
end module test_oniom
