module test_aiss
!***************************************************************
!* Unit tests for the aISS docking port (src/qcg/aiss): geometry
!* helpers, IFF parameters, wall energies vs. CREST constraints,
!* the IFF-mode localization, IFF energies and a short docking.
!***************************************************************
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use crest_parameters
  use crest_calculator
  use strucrd
  use aiss_param,pi_aiss => pi
  implicit none
  private

  public :: collect_aiss

  real(wp),parameter :: thr = 1.0e-9_wp

!> water monomers (Bohr)
  integer,parameter :: wat_at(3) = [8,1,1]
  real(wp),parameter :: wat_xyz(3,3) = reshape([ &
  &  0.000000_wp, 0.000000_wp, 0.221664_wp, &
  &  0.000000_wp, 1.430900_wp,-0.886655_wp, &
  &  0.000000_wp,-1.430900_wp,-0.886655_wp],[3,3])

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for the aISS docking
!========================================================================================!
!========================================================================================!

  subroutine collect_aiss(testsuite)
    type(unittest_type),allocatable,intent(out) :: testsuite(:)
!&<
    testsuite = [ &
    new_unittest("rotation matrix orthonormal   ",test_rotmat), &
    new_unittest("rigid move keeps distances    ",test_move2), &
    new_unittest("docking CMA                   ",test_cmadock), &
    new_unittest("packed eigensolver            ",test_rsp_packed), &
    new_unittest("IFF parameters                ",test_iff_param), &
    new_unittest("polynomial wall energy        ",test_wall_poly), &
    new_unittest("logfermi wall energy          ",test_wall_logfermi), &
    new_unittest("walls = CREST constraints     ",test_wall_vs_constraint), &
    new_unittest("gene sorting                  ",test_sort6) &
#ifdef WITH_TBLITE
    ,new_unittest("IFF water dimer (GFN2 LMO)    ",test_iff_water_gfn2) &
#ifdef WITH_GFN0
    ,new_unittest("IFF water dimer (GFN0 LMO)    ",test_iff_water_gfn0) &
#endif
#ifdef WITH_GFNFF
    ,new_unittest("docking water dimer (GFN-FF)  ",test_dock_water_gfnff) &
#endif
#endif
    ]
!&>
  end subroutine collect_aiss

!========================================================================================!

  subroutine test_rotmat(error)
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: rot(3,3),id(3,3),det
    integer :: i
    call rotmat([0.0_wp,0.0_wp,0.0_wp,0.3_wp,-1.1_wp,2.4_wp],rot)
    id = matmul(transpose(rot),rot)
    do i = 1,3
      id(i,i) = id(i,i)-1.0_wp
    end do
    call check(error,maxval(abs(id)) < thr,.true.)
    if (allocated(error)) return
    det = rot(1,1)*(rot(2,2)*rot(3,3)-rot(2,3)*rot(3,2)) &
    &    -rot(1,2)*(rot(2,1)*rot(3,3)-rot(2,3)*rot(3,1)) &
    &    +rot(1,3)*(rot(2,1)*rot(3,2)-rot(2,2)*rot(3,1))
    call check(error,det,1.0_wp,thr=thr)
  end subroutine test_rotmat

  subroutine test_move2(error)
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: xnew(3,3),gene(6)
    gene = [1.0_wp,-2.0_wp,0.5_wp,0.7_wp,0.2_wp,-1.3_wp]
    call move2(3,wat_xyz,xnew,gene)
    call check(error,norm2(xnew(:,1)-xnew(:,2)),norm2(wat_xyz(:,1)-wat_xyz(:,2)),thr=thr)
    if (allocated(error)) return
    call check(error,norm2(xnew(:,2)-xnew(:,3)),norm2(wat_xyz(:,2)-wat_xyz(:,3)),thr=thr)
    if (allocated(error)) return
    !> the molecule is rotated about the origin, then shifted by gene(1:3)
    call check(error,norm2(xnew(:,1)-gene(1:3)),norm2(wat_xyz(:,1)),thr=thr)
  end subroutine test_move2

  subroutine test_cmadock(error)
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: r(3)
    call cmadock(3,3,wat_at,wat_xyz,r)
    call check(error,r(2),0.0_wp,thr=thr)
    if (allocated(error)) return
    !> 15.9994*0.221664 + 2*1.0079*(-0.886655) over the total mass
    call check(error,r(3),(15.9994_wp*0.221664_wp-2.0_wp*1.00790_wp*0.886655_wp) &
    &          /(15.9994_wp+2.0_wp*1.00790_wp),thr=1.0e-6_wp)
    if (allocated(error)) return
    !> only the first atom
    call cmadock(1,3,wat_at,wat_xyz,r)
    call check(error,r(3),0.221664_wp,thr=1.0e-6_wp)
  end subroutine test_cmadock

  subroutine test_rsp_packed(error)
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: a(6),w(3),z(3,3)
    !> [[2,1,0],[1,2,0],[0,0,5]] packed upper = [a11,a12,a22,a13,a23,a33]
    a = [2.0_wp,1.0_wp,2.0_wp,0.0_wp,0.0_wp,5.0_wp]
    call rsp_packed(a,3,w,z)
    call check(error,w(1),1.0_wp,thr=thr)
    if (allocated(error)) return
    call check(error,w(2),3.0_wp,thr=thr)
    if (allocated(error)) return
    call check(error,w(3),5.0_wp,thr=thr)
  end subroutine test_rsp_packed

  subroutine test_iff_param(error)
    type(error_type),allocatable,intent(out) :: error
    call set_iff_param(0)
    call check(error,par_rep_scal,0.1232398_wp,thr=thr)
    if (allocated(error)) return
    call check(error,par_pos_sig,par_pos1_lp,thr=thr)
    if (allocated(error)) return
    call check(error,val_e(8),6.0_wp,thr=thr)
    if (allocated(error)) return
    call check(error,val_e(29),11.0_wp,thr=thr)
    if (allocated(error)) return
    call check(error,lpatom(8),1)
    if (allocated(error)) return
    call check(error,sigatom(17),1)
  end subroutine test_iff_param

!========================================================================================!

  subroutine test_wall_poly(error)
    type(error_type),allocatable,intent(out) :: error
    type(aiss_wall) :: w(1)
    real(wp) :: xyz(3,1),e,s
    w(1)%radius = [6.0_wp,4.0_wp,3.0_wp]
    xyz(:,1) = [3.0_wp,2.0_wp,1.0_wp]
    s = (3.0_wp/6.0_wp)**2+(2.0_wp/4.0_wp)**2+(1.0_wp/3.0_wp)**2
    e = 0.0_wp
    call cavity_egrad(1,xyz,1,w,e)
    call check(error,e,s**15,thr=thr)
    if (allocated(error)) return
    !> QCG-mode ellipsoid is the same function for alpha=30
    e = 0.0_wp
    call cavitye(1,xyz,w(1)%radius,e)
    call check(error,e,s**15,thr=thr)
  end subroutine test_wall_poly

  subroutine test_wall_logfermi(error)
    type(error_type),allocatable,intent(out) :: error
    type(aiss_wall) :: w(1)
    real(wp) :: xyz(3,2),e,d
    w(1)%radius = [5.0_wp,5.0_wp,5.0_wp]
    w(1)%logfermi = .true.
    w(1)%temp = 400.0_wp
    w(1)%beta = 5.0_wp
    w(1)%list = [2]
    xyz(:,1) = [9.0_wp,0.0_wp,0.0_wp]  !> not in the list
    xyz(:,2) = [0.0_wp,6.0_wp,0.0_wp]
    d = 6.0_wp
    e = 0.0_wp
    call cavity_egrad(2,xyz,1,w,e)
    call check(error,e,kB*400.0_wp*log(1.0_wp+exp(5.0_wp*(d-5.0_wp))),thr=thr)
  end subroutine test_wall_logfermi

  subroutine test_wall_vs_constraint(error)
    !> screening walls and the constraints used in the optimizations agree
    use aiss_dock,only:aiss_add_walls
    type(error_type),allocatable,intent(out) :: error
    type(aiss_wall) :: w(2)
    type(calcdata) :: calc
    real(wp) :: xyz(3,3),e_iff,e_con,etmp,grd(3,3)
    integer :: i
    xyz = wat_xyz*2.5_wp
    w(1)%radius = [4.0_wp,3.5_wp,3.0_wp]
    w(2)%radius = [3.0_wp,3.0_wp,3.0_wp]
    w(2)%list = [2,3]
    do i = 1,2
      w(:)%logfermi = (i == 2)
      e_iff = 0.0_wp
      call cavity_egrad(3,xyz,2,w,e_iff)
      call calc%reset()
      call aiss_add_walls(calc,3,w)
      e_con = 0.0_wp
      block
        integer :: k
        do k = 1,calc%nconstraints
          call calc_constraint(3,xyz,calc%cons(k),etmp,grd)
          e_con = e_con+etmp
        end do
      end block
      call check(error,e_con,e_iff,thr=1.0e-8_wp*max(1.0_wp,abs(e_iff)))
      if (allocated(error)) return
    end do
  end subroutine test_wall_vs_constraint

  subroutine test_sort6(error)
    type(error_type),allocatable,intent(out) :: error
    real(wp) :: g(7,4)
    g = 0.0_wp
    g(7,:) = [3.0_wp,-1.0_wp,2.0_wp,-5.0_wp]
    g(1,:) = [1.0_wp,2.0_wp,3.0_wp,4.0_wp]
    call sort6(4,g)
    call check(error,g(7,1),-5.0_wp,thr=thr)
    if (allocated(error)) return
    call check(error,g(1,1),4.0_wp,thr=thr)
    if (allocated(error)) return
    call check(error,g(7,4),3.0_wp,thr=thr)
  end subroutine test_sort6

!========================================================================================!

#ifdef WITH_TBLITE
  subroutine water_dimer_iff(lmolvl,d,e_cold,e_disp,io)
    !> IFF setup of two waters and energies at two placements of B
    use aiss_iff_prepare,only:iff_precomp
    use aiss_iff_ini,only:init_iff
    use aiss_iff_energy,only:iff_e
    character(len=*),intent(in) :: lmolvl
    type(iff_data),intent(out) :: d
    real(wp),intent(out) :: e_cold,e_disp
    integer,intent(out) :: io
    type(coord) :: mol
    real(wp) :: icoord(6),icoord0(6),etmp
    mol%nat = 3
    mol%at = wat_at
    mol%xyz = wat_xyz
    mol%chrg = 0
    mol%uhf = 0
    call set_iff_param(0)
    call d%allocate(3,3)
    call iff_precomp(mol,lmolvl,300.0_wp,d,1,etmp,io)
    if (io /= 0) return
    call iff_precomp(mol,lmolvl,300.0_wp,d,2,etmp,io)
    if (io /= 0) return
    icoord = 0.0_wp
    icoord0 = 0.0_wp
    call init_iff(d,0,icoord,icoord0,0)
    call iff_e(d,.false.,0,e_cold,icoord)
    !> B 5.5 Bohr above A, its hydrogens pointing at the oxygen of A
    call iff_e(d,.false.,-1,e_disp,[0.0_wp,0.0_wp,5.5_wp,0.0_wp,0.0_wp,0.0_wp])
  end subroutine water_dimer_iff

  subroutine test_iff_water_gfn2(error)
    type(error_type),allocatable,intent(out) :: error
    type(iff_data) :: d
    real(wp) :: e_cold,e_disp
    integer :: io
    call water_dimer_iff('gfn2',d,e_cold,e_disp,io)
    call check(error,io,0)
    if (allocated(error)) return
    !> 2 LP + 2 sigma LMOs; the complete overlap of the cold fusion is huge
    call check(error,d%nlmo1,4)
    if (allocated(error)) return
    call check(error,abs(sum(d%q1)) < 1.0e-6_wp,.true.)
    if (allocated(error)) return
    call check(error,e_cold > 1.0e10_wp,.true.)
    if (allocated(error)) return
    call check(error,e_disp,-3.809342069172195e-03_wp,thr=1.0e-7_wp)
  end subroutine test_iff_water_gfn2

#ifdef WITH_GFN0
  subroutine test_iff_water_gfn0(error)
    type(error_type),allocatable,intent(out) :: error
    type(iff_data) :: d
    real(wp) :: e_cold,e_disp
    integer :: io
    call water_dimer_iff('gfn0',d,e_cold,e_disp,io)
    call check(error,io,0)
    if (allocated(error)) return
    call check(error,d%nlmo1,4)
    if (allocated(error)) return
    call check(error,e_disp,-6.125897441531194e-03_wp,thr=1.0e-7_wp)
  end subroutine test_iff_water_gfn0
#endif

#ifdef WITH_GFNFF
  subroutine test_dock_water_gfnff(error)
    use aiss_dock,only:aiss_dock_run
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: molA,molB,best
    type(aiss_config) :: cfg
    real(wp) :: eint
    integer :: io
    molA%nat = 3
    molA%at = wat_at
    molA%xyz = wat_xyz
    molB = molA
    call cfg%set_fast()
    cfg%optlvl = 'gfnff'
    cfg%lmolvl = 'gfn2'
    cfg%n_opt = 2
    cfg%samerand = .true.
    cfg%wrfiles = .false.
    call aiss_dock_run(molA,molB,cfg,best,eint,io)
    call check(error,io,0)
    if (allocated(error)) return
    call check(error,best%nat,6)
    if (allocated(error)) return
    !> GFN-FF water dimer binding is about -4 to -6 kcal/mol
    call check(error,eint < -2.0_wp.and.eint > -10.0_wp,.true.)
  end subroutine test_dock_water_gfnff
#endif
#endif

!========================================================================================!
end module test_aiss
