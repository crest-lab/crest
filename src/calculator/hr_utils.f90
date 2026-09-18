module hr_utils
  use iso_fortran_env,only:wp => real64
  use crest_calculator
  use crest_parameters
  use optimize_maths
  use modelhessian_module
  use axis_module
  use strucrd
  implicit none
  private

  public initialize_hessian

contains

  subroutine initialize_hessian(calc,type,xyz,nat,at,hess,hguess,pr) !>Matrix is forced to be positive definite
    type(calcdata),intent(inout) :: calc
    integer,intent(in) :: type
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in) :: at(nat),nat
    real(wp),intent(inout) :: hess(:)
    real(wp),optional,intent(in) :: hguess
    logical,intent(in) :: pr
    type(calcdata),allocatable :: newcalc
    type(mhparam) :: mhset
    integer :: k,i,j,idx,io,nat3

    real(wp),allocatable :: hess_full(:,:)

    real(wp),allocatable :: pmode(:,:),grad(:,:)
    real(wp) :: rot(3),dumi
    logical :: linear
    type(coord) :: mol

    nat3 = 3*nat

    !!$omp critical
    !allocate (pmode(nat3,1)) ! dummy allocated
    !!$omp end critical

    !$omp critical
    allocate (newcalc)
    allocate (hess_full(nat3,nat3),source=0.0_wp)
    !$omp end critical

    select case (type)
    case (0) !>Initialize as a scaled identity
      if (present(hguess)) then
        k = 0
        do i = 1,nat3
          do j = 1,i
            k = k+1
            if (i /= j) then
              hess(k) = 0.0_wp
            else
              hess(k) = hguess
            end if
          end do
        end do
      else
        write (stdout,*) "No hguess provided"
      end if
    case (1)
      !$omp critical
      !write(stdout,*) calc%calcs(1)%chrg
      call newcalc%create('gfnff',chrg=calc%calcs(1)%chrg,uhf=calc%calcs(1)%uhf) !> Different levels?? and what happens to solvent??
      !$omp end critical
      call numhess1(nat,at,xyz,newcalc,hess_full(:,:),io)
      call dsqtoh(nat3,hess_full(:,:),hess(:)) !>Pack Hessian
    case (2)
        !$omp critical
      call newcalc%create('gfn0',chrg=calc%calcs(1)%chrg,uhf=calc%calcs(1)%uhf) !> Different levels?? and what happens to solvent??
      !$omp end critical
      call numhess1(nat,at,xyz,newcalc,hess_full(:,:),io)
      call dsqtoh(nat3,hess_full(:,:),hess(:))
    case (3)
        !$omp critical
      call newcalc%create('gfn1',chrg=calc%calcs(1)%chrg,uhf=calc%calcs(1)%uhf) !> Different levels?? and what happens to solvent??
      !$omp end critical
      call numhess1(nat,at,xyz,newcalc,hess_full(:,:),io)
      call dsqtoh(nat3,hess_full(:,:),hess(:))
    case (4)
        !$omp critical
      call newcalc%create('gfn2',chrg=calc%calcs(1)%chrg,uhf=calc%calcs(1)%uhf) !> Different levels?? and what happens to solvent??
      !$omp end critical
      call numhess1(nat,at,xyz,newcalc,hess_full(:,:),io)
      call dsqtoh(nat3,hess_full(:,:),hess(:))
    case (5)
      !$omp critical
      mhset%model = calc%mh_type
      call modhes(calc,mhset,nat,xyz,at,hess(:),pr)
      !$omp end critical
    end select

    !call axis(nat,at,xyz,rot,dumi)
    !linear = (rot(3) .lt. 1.d-10).or.(nat == 2)

    !if (.not.linear) then
    !    if (calc%nfreeze == 0) then
    !      call trproj(nat,nat3,xyz,hess,.false.,0,pmode,1)  !> normal
    !    else
    !      call trproj(nat,nat3,xyz,hess,.false.,calc%freezelist) !> fozen atoms
    !    end if
    !end if

    call force_positive_definiteness(hess,nat3)

  end subroutine initialize_hessian

  subroutine force_positive_definiteness(hess,nat3)
    real(wp),intent(inout) :: hess(:)
    integer,intent(in) :: nat3
    real(wp),allocatable :: eigvec(:,:),eigval(:)
    real(wp),allocatable :: work(:)
    integer,allocatable :: iwork(:)
    integer :: lwork,liwork,info,i,j,k,l
    real(wp) :: elow,damp

    allocate (eigvec(nat3,nat3),eigval(nat3))
    lwork = 1+6*nat3+2*nat3*nat3
    liwork = 8*nat3
    allocate (work(lwork),iwork(liwork))

    call dspevd('V','U',nat3,hess(:),eigval,eigvec,nat3, &
                work,lwork,iwork,liwork,info)

    if (info /= 0) then
      write (*,*) "dspevd failed, info = ",info
      stop
    end if

    elow = minval(eigval)
    damp = max(1.0e-4_wp-elow,0.0_wp)
    eigval = eigval+damp

    hess(:) = 0.0_wp
    k = 0
    do j = 1,nat3
      do i = 1,j
        k = k+1
        hess(k) = 0.0_wp
        do l = 1,nat3
          hess(k) = hess(k)+eigval(l)*eigvec(i,l)*eigvec(j,l)
        end do
      end do
    end do

    deallocate (eigvec,eigval,work,iwork)

  end subroutine force_positive_definiteness

end module hr_utils
