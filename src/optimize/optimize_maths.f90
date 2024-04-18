!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2022 Philipp Pracht
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

!> Implementation of linear algebra & co. needed for optimization routines
module optimize_maths
  use iso_fortran_env,only:wp => real64,sp => real32
  implicit none
  private

  public :: trproj
  interface trproj
    module procedure :: trproj_normal
    module procedure :: trproj_frozen
  end interface
  public :: detrotra8
  public :: solver_sdavidson
  public :: solver_sspevx
  public :: solver_ssyevx

!========================================================================================!
!========================================================================================!
contains !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine detrotra8(linear,n,xyz,h,eig)
    implicit none
    integer,intent(in)     :: n
    real(wp),intent(in)    :: xyz(3,n)    !> values from projected Lindh diag
    real(wp),intent(in)    :: h(3*n,3*n)  !> values from projected Lindh diag
    real(wp),intent(inout) :: eig(3*n)    !> eigenvectors from projected Lindh diag
    logical,intent(in)     :: linear
    integer               :: i,j,k,kk,ii,nn,n3,nend
    integer,allocatable   :: ind(:)
    real(wp),allocatable  :: tmp(:,:)
    real(wp),allocatable  :: e(:)
    real(wp)              :: a0,b0,c0

    n3 = 3*n
    allocate (tmp(3,n),e(n3),ind(n3))
    nn = 0
    do ii = 1,n3
      if (eig(ii) .gt. 0.05_wp) cycle  !> only lowest checked
      kk = 0
      do j = 1,n
        do k = 1,3
          kk = kk+1
          tmp(k,j) = xyz(k,j)+h(kk,ii) !> distort along mode ii
        end do
      end do
      c0 = 0            !> compared all interatomic distances of original and distortet geom.
      do i = 2,n
        do j = 1,i-1
          a0 = sqrt((xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2+(xyz(3,i)-xyz(3,j))**2)
          b0 = sqrt((tmp(1,i)-tmp(1,j))**2+(tmp(2,i)-tmp(2,j))**2+(tmp(3,i)-tmp(3,j))**2)
          c0 = c0+(a0-b0)**2
        end do
      end do
      nn = nn+1
      e(nn) = sqrt(c0/n)*abs(eig(ii))  !> weight by Lindh eigenvalue
      ind(nn) = nn
    end do

    call detrotra_qsort(e,1,nn,ind)   !> sort

    nend = 6
    if (linear) nend = 5

    do i = 1,nend
      eig(ind(i)) = 0.0_wp  !> identifier for rot/tra
    end do
    deallocate (ind,e,tmp)
    return
  contains
    recursive subroutine detrotra_qsort(a,first,last,ind)
      implicit none
      real(wp) :: a(*),x,t
      integer  :: ind(*)
      integer  :: first,last
      integer  :: i,j,ii
      x = a((first+last)/2)
      i = first
      j = last
      do
        do while (a(i) < x)
          i = i+1
        end do
        do while (x < a(j))
          j = j-1
        end do
        if (i >= j) exit
        t = a(i); a(i) = a(j); a(j) = t
        ii = ind(i); ind(i) = ind(j); ind(j) = ii
        i = i+1
        j = j-1
      end do
      if (first < i-1) call qsort(a,first,i-1,ind)
      if (j+1 < last) call qsort(a,j+1,last,ind)
    end subroutine detrotra_qsort
  end subroutine detrotra8

!========================================================================================!

  subroutine trproj_normal(natoms,nat3,xyz,hess,ldebug,nmode,mode,ndim)
!**********************************************************************
!*  subroutine trproj drives projection of hessian out of the
!*  space of translational and rotational motions:
!*  first get xyz c.m.; second get transl. and rot. projection matrix
!*
!*  get center of mass coordinates with unit mass
!*
!* Input
!*   natoms  = number of atoms
!*   nat3    = 3*natoms
!*   xyz     = cartesian coordinates
!*   ldebug  = debug flag = .true. for debugging
!*
!* Ouput
!*       xyzucm  = temporary c.m. coordinates
!*
!*       hess    = projected hessian out of space of transl. and rot.
!*                 motion
!**********************************************************************
    implicit none
    !> Input
    logical,intent(in) :: ldebug
    integer,intent(in) :: natoms,nat3,nmode,ndim
    real(wp),dimension(3,natoms) :: xyz
    real(wp),dimension(nat3,ndim) :: mode
    !> Ouput
    real(wp),dimension(nat3*(nat3+1)/2) :: hess
    !> Local
    integer :: i
    real(wp) :: xm,ym,zm
    real(wp),dimension(3,natoms) :: xyzucm

    !> temporary coordinates, c.m. shifted
    xyzucm(:,:) = xyz(:,:)
    xm = 0.0_wp
    ym = 0.0_wp
    zm = 0.0_wp
    do i = 1,natoms
      xm = xm+xyzucm(1,i)
      ym = ym+xyzucm(2,i)
      zm = zm+xyzucm(3,i)
    end do
    xm = xm/natoms
    ym = ym/natoms
    zm = zm/natoms
    do i = 1,natoms
      xyzucm(1,i) = xyzucm(1,i)-xm
      xyzucm(2,i) = xyzucm(2,i)-ym
      xyzucm(3,i) = xyzucm(3,i)-zm
    end do

    !> get translational and rotational projection matrix
    call gtrprojm(natoms,nat3,xyzucm,hess,ldebug,nmode,mode,ndim)

    return
  end subroutine trproj_normal

!========================================================================================!

  subroutine gtrprojm(natoms,nat3,xyzucm,hess,ldebug,nmode,mode,ndim)
!*******************************************************************
!* calculating the translational-rotational projection matrix
!*
!* Input
!*   natoms  = number of atoms
!*   nat3    = 3*natoms
!*   xyzucm  = coords c.m. from gxyzucm.f
!*   hess    = hessian
!*   ldebug  = debug flag = .true. for debugging
!*
!* Ouput
!*   fmat    = F-matrix with translational and rotational vectors
!*   pmat    = projection matrix P = (1-FFt)
!*   hess    = projected hessian
!*******************************************************************
    implicit none
    !> Input
    logical,intent(in) :: ldebug
    integer,intent(in) :: natoms,nat3,nmode,ndim
    real(wp),dimension(3,natoms) :: xyzucm
    real(wp),dimension(nat3,ndim) :: mode
    !> Ouput
    real(wp),dimension(nat3*(nat3+1)/2) :: hess
    !> Local
    integer :: i,ii,iii
    real(wp),allocatable :: fmat(:,:)
    integer :: nprj,fnrozen

    nprj = 6
    if (nmode .gt. 0) nprj = nprj+nmode
    !if (nmode .lt. 0) nprj = nprj + fixset%n * 3
    allocate (fmat(nat3,nprj))
    fmat(:,:) = 0.0_wp

    if (nmode .ge. 0) then
      do i = 1,natoms
        !> translation vectors
        do ii = 1,3
          fmat(3*(i-1)+ii,ii) = 1.0_wp
        end do
        !> rotational vectors
        fmat(3*(i-1)+1,4) = 0.0_wp
        fmat(3*(i-1)+2,4) = -xyzucm(3,i)
        fmat(3*(i-1)+3,4) = xyzucm(2,i)
        fmat(3*(i-1)+1,5) = xyzucm(3,i)
        fmat(3*(i-1)+2,5) = 0.0_wp
        fmat(3*(i-1)+3,5) = -xyzucm(1,i)
        fmat(3*(i-1)+1,6) = -xyzucm(2,i)
        fmat(3*(i-1)+2,6) = xyzucm(1,i)
        fmat(3*(i-1)+3,6) = 0.0_wp
      end do
    end if

    !> NMF
    if (nmode .gt. 0) then
      do i = 1,nmode
        fmat(1:nat3,6+i) = mode(1:nat3,i)
      end do
    end if

    if (ldebug) then
      write (*,'(a)')
      write (*,'(a)') ' Basis vectors before orthonormalization'
      write (*,'(3e22.14)') fmat
    end if

    !> do Gramm-Schmidt orthogonalization
    call dblckmgs(nat3,nprj,nat3,fmat)

    !> do projection
    call dsyprj(nat3,nprj,fmat,nat3,hess)

    deallocate (fmat)
    return
  end subroutine gtrprojm

!========================================================================================!

  subroutine trproj_frozen(natoms,nat3,xyz,hess,ldebug,freezelist)
!*************************************************************************
!*  subroutine trproj drives projection of hessian out of the
!*  space of translational and rotational motions:
!*  first get xyz c.m.; second get transl. and rot. projection matrix
!*  This is the frozen atom version
!*  get center of mass coordinates with unit mass
!*
!* Input
!*   natoms  = number of atoms
!*   nat3    = 3*natoms
!*   xyz     = cartesian coordinates
!*   ldebug  = debug flag = .true. for debugging
!*   freezelist = an array of booleans that are .true. for frozen atoms
!*
!* Ouput
!*       xyzucm  = temporary c.m. coordinates
!*
!*       hess    = projected hessian out of space of transl. and rot.
!*                 motion
!************************************************************************
    implicit none
    !> Input
    logical,intent(in) :: ldebug
    integer,intent(in) :: natoms,nat3
    logical,intent(in) :: freezelist(natoms)
    real(wp),dimension(3,natoms) :: xyz
    !> Ouput
    real(wp),dimension(nat3*(nat3+1)/2) :: hess
    !> Local
    integer :: i
    real(wp) :: xm,ym,zm
    real(wp),dimension(3,natoms) :: xyzucm

    !> temporary coordinates, c.m. shifted
    xyzucm(:,:) = xyz(:,:)
    xm = 0.0_wp
    ym = 0.0_wp
    zm = 0.0_wp
    do i = 1,natoms
      xm = xm+xyzucm(1,i)
      ym = ym+xyzucm(2,i)
      zm = zm+xyzucm(3,i)
    end do
    xm = xm/natoms
    ym = ym/natoms
    zm = zm/natoms
    do i = 1,natoms
      xyzucm(1,i) = xyzucm(1,i)-xm
      xyzucm(2,i) = xyzucm(2,i)-ym
      xyzucm(3,i) = xyzucm(3,i)-zm
    end do

    !> get translational and rotational projection matrix
    call gtrprojm_frozen(natoms,nat3,xyzucm,hess,ldebug,freezelist)

    return
  end subroutine trproj_frozen

!========================================================================================!

  subroutine gtrprojm_frozen(natoms,nat3,xyzucm,hess,ldebug,freezelist)
!*******************************************************************
!* calculating the translational-rotational projection matrix
!* taking into account frozen atoms
!*
!* Input
!*   natoms  = number of atoms
!*   nat3    = 3*natoms
!*   xyzucm  = coords c.m. from gxyzucm.f
!*   hess    = hessian
!*   ldebug  = debug flag = .true. for debugging
!*
!* Ouput
!*   fmat    = F-matrix with translational and rotational vectors
!*   pmat    = projection matrix P = (1-FFt)
!*   hess    = projected hessian
!******************************************************************
    implicit none
    !> Input
    logical,intent(in) :: ldebug
    integer,intent(in) :: natoms,nat3
    logical,intent(in) :: freezelist(natoms)
    real(wp),dimension(3,natoms) :: xyzucm
    !> Ouput
    real(wp),dimension(nat3*(nat3+1)/2) :: hess
    !> Local
    integer :: i,ii,iii,f
    real(wp),allocatable :: fmat(:,:)
    integer :: nprj,nfrozen

    nfrozen = count(freezelist,1)
    nprj = 6
    nprj = nprj+nfrozen*3  !> correct dimension with frozen deg. of freedom
    allocate (fmat(nat3,nprj))
    fmat(:,:) = 0.0_wp

    !> exact fixing version
    do i = 1,natoms
      !> rotational vectors
      fmat(3*(i-1)+1,1) = 0.0_wp
      fmat(3*(i-1)+2,1) = -xyzucm(3,i)
      fmat(3*(i-1)+3,1) = xyzucm(2,i)
      fmat(3*(i-1)+1,2) = xyzucm(3,i)
      fmat(3*(i-1)+2,2) = 0.0_wp
      fmat(3*(i-1)+3,2) = -xyzucm(1,i)
      fmat(3*(i-1)+1,3) = -xyzucm(2,i)
      fmat(3*(i-1)+2,3) = xyzucm(1,i)
      fmat(3*(i-1)+3,3) = 0.0_wp
    end do
    i = 0
    do f = 1,natoms
      if (freezelist(f)) then
        i = i+1
        iii = f
        do ii = 1,3
          fmat(3*(iii-1)+ii,3+(i-1)*3+ii) = 1.0_wp
        end do
      end if
    end do

    if (ldebug) then
      write (*,'(a)')
      write (*,'(a)') ' Basis vectors before orthonormalization'
      write (*,'(3e22.14)') fmat
    end if

    !> do Gramm-Schmidt orthogonalization
    call dblckmgs(nat3,nprj,nat3,fmat)

    !> do projection
    call dsyprj(nat3,nprj,fmat,nat3,hess)

    deallocate (fmat)
    return
  end subroutine gtrprojm_frozen

!========================================================================================!

  subroutine dhtosq(n,a,b)
!********************************************
!* expand trigonal matrix b to full matrix a
!* a and b may refer to the same address
!********************************************
    implicit none
    integer,intent(in)  :: n
    real(wp),intent(out) :: a(n,n)
    real(wp),intent(in)  :: b(n*(n+1)/2)
    integer :: i,j,ioff

    do i = n,1,-1
      ioff = i*(i-1)/2
      do j = i,1,-1
        a(j,i) = b(ioff+j)
      end do
    end do

    do i = 1,n
      do j = 1,i-1
        a(i,j) = a(j,i)
      end do
    end do
    return
  end subroutine dhtosq

!========================================================================================!

  subroutine dsyprj(nbdim,m,bmat,n,asym)
!**********************************************************************
!* subroutine dsyprj
!* Performs projection of a real symmetric matrix asym packed in an
!* upper triangular form:
!*
!*     asym = (1 - bmat*bmat')*asym*(1 - bmat*bmat')
!*
!* where (1 - bmat*bmat') is a projector constructed from bmat
!* The matrix asym is overwritten on output; bmat is left unchanged.
!* Operation count is proportional to n*n*m
!*
!* Input/Output:
!* nbdim - First dimension of bmat as declared in the calling routine
!* m     - Number of columns in the matrix bmat
!* bmat  - n*m matrix used to build projector
!* n     - Dimension of the matrix asym
!* asym  - Symmetric matrix in packed upper triangular form
!*********************************************************************
    implicit none
!> Input/Output:
    integer,intent(in)    :: nbdim,m,n
    real(wp),intent(in)    :: bmat(nbdim,m)
    real(wp),intent(inout) :: asym(n*(n+1)/2)
!> Local:
    integer  :: i,j,ij
    real(wp),allocatable :: scrb(:,:)
    real(wp),allocatable :: scra(:,:)
!> external BLAS routines
    external :: dsymm
    external :: dgemm
!-----------------------------------------------------------------
    allocate (scrb(n,m),scra(n,n))

!> Expand trigonal matrix asym to full matrix on scra
    call dhtosq(n,scra,asym)
!> Calculate scrb = asym*bmat (BLAS)
    call dsymm('l','u',n,m,1.0_wp,scra,n,bmat,nbdim,0.0_wp,scrb,n)
!> Calculate scra = scrb*bmat' (BLAS)
    call dgemm('n','t',n,n,m,1.0_wp,scrb,n,bmat,nbdim,0.0_wp,scra,n)

!> Calculate asym = asym - scra - scra'
    do i = 1,n
      do j = 1,i
        ij = i*(i-1)/2+j
        asym(ij) = asym(ij)-scra(i,j)-scra(j,i)
      end do
    end do

!> Calculate scrb' = scra'*bmat (BLAS)
    call dgemm('t','n',n,m,n,1.0_wp,scra,n,bmat,nbdim,0.0_wp,scrb,n)
!> Calculate scra = bmat*scrb' (BLAS)
    call dgemm('n','t',n,n,m,1.0_wp,bmat,nbdim,scrb,n,0.0_wp,scra,n)

!> Calculate asym = asym + scra
    do i = 1,n
      do j = 1,i
        ij = i*(i-1)/2+j
        asym(ij) = asym(ij)+scra(i,j)
      end do
    end do

    deallocate (scra,scrb)
    return
  end subroutine dsyprj

!========================================================================================!

  subroutine dblckmgs(m,n,ndim,darray)
!**********************************************************************
!* subroutine dblckmgs
!* Subroutine performs modified Gramm-Schmidt orthonormalization of a
!* real matrix. Orthonormalization is done in-place, so the darray is
!* overwritten on exit. Linearly dependent vectors are set to zero.
!*
!* Input:
!* m      - Number of rows in the matrix darray
!* n      - Number of columns in the matrix darray
!* ndim   - First array dimension as declared in the calling routine
!* darray - Array to be orthonormalized
!*
!* Output:
!* darray - Orthonormalized array
!*********************************************************************
    implicit none
!> Input/Output:
    integer,intent(in) :: m,n,ndim
    real(wp),dimension(ndim,n),intent(out) :: darray
!> Local:
    integer :: ii,jj,kk,ll,ibsize,nblcks,istrt,jstrt,iend,ncol,ierr
    real(wp) :: tmp
    real(wp),dimension(:,:),allocatable :: smat
    real(wp) :: thr
!> external BLAS routines
    external :: dgemm
    real(wp),external :: ddot
    external :: dscal
    external :: daxpy
!-----------------------------------------------------------------
    thr = epsilon(1.0_wp) !Threshold for zero vectors
!-----------------------------------------------------------------
! Block size optimized for Athlon 1200 MHz with 2.0GB memory for
! matrices up to 5000x5000
    ibsize = 60
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Allocate overlap matrix
!-----------------------------------------------------------------
    allocate (smat(ibsize,ibsize))

!-----------------------------------------------------------------
! Calculate the number of blocks
!-----------------------------------------------------------------
    nblcks = (n+ibsize-1)/ibsize
    ibsize = min(n,ibsize)

!-----------------------------------------------------------------
! Orthogonalize the first block using modified Schmidt
!-----------------------------------------------------------------
    do ii = 1,ibsize
      !> Linear dependence check by
      !> dot product of the column vector (BLAS)
      tmp = ddot(m,darray(:,ii),1,darray(:,ii),1)
      if (tmp < thr) then
        darray(1:m,ii) = 0.0_wp
        cycle
      end if
      !> scale column vectors
      tmp = 1.0_wp/sqrt(tmp)
      call dscal(m,tmp,darray(:,ii),1)
      do jj = ii+1,ibsize
        tmp = ddot(m,darray(:,ii),1,darray(:,jj),1)
        call daxpy(m,-tmp,darray(:,ii),1,darray(:,jj),1)
      end do
    end do
!-----------------------------------------------------------------
! Loop over remaining blocks
!-----------------------------------------------------------------
    do ii = 1,nblcks-1
!> Initial and final column and number of columns in the block ii+1
      istrt = ii*ibsize+1
      iend = (ii+1)*ibsize
      iend = min(n,iend)
      ncol = iend-istrt+1
!> Orthogonalize the block ii+1 against the previous ones
      do jj = 1,ii
!> Initial index of the block jj
        jstrt = (jj-1)*ibsize+1
        call dgemm('t','n',ibsize,ncol,m,1.0_wp,darray(1,jstrt), &
        &          ndim,darray(1,istrt),ndim,0.0_wp,smat,ibsize)
        call dgemm('n','n',m,ncol,ibsize,-1.0_wp,darray(1,jstrt), &
        &          ndim,smat,ibsize,1.0_wp,darray(1,istrt),ndim)
      end do

!> Othogonalize vectors on the block ii+1 among themselves using modified Schmidt
      do kk = istrt,iend
        tmp = ddot(m,darray(:,kk),1,darray(:,kk),1)
!> Linear dependence case
        if (tmp < thr) then
          darray(1:m,kk) = 0.0_wp
          cycle
        end if
        tmp = 1.0_wp/sqrt(tmp)
        call dscal(m,tmp,darray(:,kk),1)
        do ll = kk+1,iend
          tmp = ddot(m,darray(:,kk),1,darray(:,ll),1)
          call daxpy(m,-tmp,darray(:,kk),1,darray(:,ll),1)
        end do
      end do
    end do

!> Clean up and return
    deallocate (smat)
    return
  end subroutine dblckmgs

!========================================================================================!

  subroutine solver_sdavidson(n,crite,Hp,C,e,fail,pr)
!*****************************************************************
!* subroutine solver_sdavidson
!*
!* Davidson method to iteratively diagonalize
!* a subspace of a matrix to provide its first
!* few lowest (or highest) eigenvalues.
!* In this version it is hard-coded to the lowest eigenvalue.
!*
!* Input:
!*      n -  dimension of the matrix to be diagonalized
!*  crite - eigenvalue convergence threshold
!*     Hp - the matrix to be diagonalized in packed form
!*      C - eigenvevtor(s)
!*     pr - print statement
!* Output:
!*      e - eigenvalues
!*   fail - exit status boolean
!*
!* Note the SINGLE PRECISION!
!*****************************************************************
    implicit none
    logical,intent(in) :: pr
    logical,parameter :: ini = .false.
    integer :: n       ! dimension
    integer,parameter :: nr = 1
    real(sp) :: crite   ! eigenvalue convergence threshold
    real(sp) :: Hp(n*(n+1)/2)  ! matrix to be diagonalized
    real(sp) :: C(n,nr) ! eigenvectors
    real(sp) :: e(nr)   ! eigenvalues
    logical,intent(out) :: fail
    !> Local
    integer,parameter :: maxiter = 100        ! maximum # of iterations
    integer :: iter,ineue(1),janf!,lun1,lun2
    integer :: iideks(maxiter),idum,j,jalt,ilauf,jneu
    integer :: l1,l2,k,LWORK,LIWORK,INFO,i,ien,ialt,memlun2
    integer,allocatable :: iwork(:)
    logical :: lconf
    real(sp),allocatable :: lun1(:,:),lun2(:,:)
    integer,parameter :: initial_dyn_array_size = 10
    real(sp) :: valn(1),uim,s,denerg
    real(sp),allocatable :: adiag(:),vecf1(:),vecf2(:),w(:)
    real(sp),allocatable :: Uaug(:,:),d(:),aux(:)
    real(sp),allocatable :: AB(:,:),av(:),tmpav(:,:)
    !> LAPACK & BLAS
    external :: sspmv
    real(sp),external :: sdot
    external :: ssyevd
    external :: saxpy

    fail = .true.

    if (pr) then
      write (*,'(/,10x,''******************************************'')')
      write (*,'(10x,''*            multi-root davidson (R4)    *'')')
      write (*,'(10x,''******************************************'',/)')
      write (*,*) 'dim ',n,' # roots ',1
    end if

    allocate (adiag(n),vecf1(n),vecf2(n),w(n),av(maxiter*(maxiter+1)/2))

    allocate (lun1(n,initial_dyn_array_size),lun2(n,initial_dyn_array_size), &
       &     source=0.0_sp)

    !> H * C for initialization
    call smwrite(n,lun1,C(:,1),1)
    call sspmv('u',n,1.0_sp,Hp,C(:,1),1,0.0_sp,vecf2,1)
    call smwrite(n,lun2,vecf2,1)

    !> make array iideks
    iideks(1) = 1
    do idum = 2,maxiter
      iideks(idum) = iideks(idum-1)+idum
    end do
    valn = 0
    lconf = .false.
    e = 0
    do i = 1,n
      adiag(i) = HP(i*(i+1)/2)
    end do
    av(1) = sdot(n,C(:,1),1,vecf2,1)

!>--------------------------------
!> Davidson algo loop
!>--------------------------------
    j = 1
    DAVIDSON: do iter = 1,maxiter-1
      lwork = 1+6*j+2*j**2
      liwork = 8*j
      allocate (Uaug(j,j),d(j),iwork(liwork),aux(lwork))
      k = 0
      do l1 = 1,j
        do l2 = 1,l1
          k = k+1
          Uaug(l2,l1) = av(k)
          Uaug(l1,l2) = av(k)
        end do
      end do
      call ssyevd('V','U',j,Uaug,j,d,aux,LWORK,IWORK,LIWORK,INFO)
      valn(1:1) = d(1:1)

      !> create and save vectors on vecf1
      vecf1 = 0.0_sp
      do i = 1,j
        call smread(n,lun1,w,i)
        uim = Uaug(i,1)
        call saxpy(n,uim,w,1,vecf1,1)
      end do

      !> calculate E*bi
      vecf2 = -valn(1)*vecf1
      !> calculate h*bi-e*bi (overwrites vecf2)
      do i = 1,j
        call smread(n,lun2,w,i)
        memlun2 = i
        uim = Uaug(i,1)
        call saxpy(n,uim,w,1,vecf2,1)
      end do
      deallocate (aux,iwork,d,Uaug)
      C(1:n,1) = vecf1

      !> calculate (h*bi - e*bi)/(e - haa); (saved as vecf2)
      vecf1 = vecf2/(valn(1)-adiag)

      !> check for convergence of Davidson algo
      denerg = abs(valn(1)-e(1))
      lconf = denerg .lt. crite
      if (pr) write (*,*) iter,lconf,denerg,valn(1:1)
      if (lconf) then
        if (pr) write (*,*) 'all roots converged'
        fail = .false.
        exit DAVIDSON
      end if

      if (j .gt. 0) then
        ialt = j
        !> orthogonalize
        do jalt = 1,ialt
          call smread(n,lun1,w,jalt)
          s = -sdot(n,w,1,vecf1,1)
          call saxpy(n,s,w,1,vecf1,1)
        end do
        !> normalize remaining
        s = sdot(n,vecf1,1,vecf1,1)
        if (s .gt. 0.00000001) then
          s = 1.0_sp/sqrt(s)
          vecf1 = vecf1*s
          ialt = ialt+1
          call smwrite(n,lun1,vecf1,jalt)
        else
          fail = .false.
          exit DAVIDSON
        end if
      end if

      !> H * C
      call sspmv('u',n,1.0_sp,Hp,vecf1,1,0.0_sp,vecf2,1)
      call smwrite(n,lun2,vecf2,memlun2+1)

      !> calculate matrix elements for next iteration
      do jalt = 1,j
        call smread(n,lun1,w,jalt)
        ilauf = iideks(j)+jalt
        av(ilauf) = sdot(n,w,1,vecf2,1)
        ilauf = ilauf+1+j
      end do
      av(iideks(j+1)) = sdot(n,vecf2,1,vecf1,1)
      !> increase expansion space and iterate further
      e = valn
      j = j+1
    end do DAVIDSON
!>--------------------------------
!> end algo loop
!>--------------------------------
    if (pr.and.fail) write (*,*) 'Warning: davidson not properly converged'

    deallocate (adiag,vecf1,vecf2,w,av,lun1,lun2)
    return
  contains
    !> write array v onto iwo
    subroutine smwrite(n,iwo,v,irec)
      implicit none
      real(sp),intent(inout),allocatable :: iwo(:,:)
      real(sp),intent(in)  :: v(n)
      integer,intent(in)  :: n,irec
      real(sp),allocatable :: tmp(:,:)
      integer :: d2,dn
      d2 = size(iwo,2)
      if (irec > d2) then
        dn = d2+d2/2+1
        allocate (tmp(n,dn))
        tmp(:,:d2) = iwo
        deallocate (iwo)
        call move_alloc(tmp,iwo)
      end if
      iwo(:,irec) = v
      return
    end subroutine smwrite
    !> read array v from iwo
    subroutine smread(n,iwo,v,irec)
      implicit none
      real(sp),intent(out) :: v(n)
      real(sp),intent(in)  :: iwo(:,:)
      integer,intent(in)  :: n,irec
      v = iwo(:,irec)
      return
    end subroutine smread
  end subroutine solver_sdavidson

!========================================================================================!

  subroutine solver_ssyevx(n,thr,A,U,e,fail)
!****************************************************************
!* subroutine solver_ssyevx
!* Wrapper for LAPACK's ssyevx routine:
!* SSYEVD computes all eigenvalues and eigenvectors of a
!* real symmetric matrix A using a divide and conquer algorithm.
!****************************************************************
    implicit none
    integer,intent(in) :: n
    real(sp),intent(in) :: thr
    real(sp),intent(inout) :: A(:,:)
    real(sp),intent(inout) :: U(:,:)
    real(sp),intent(inout) :: e(:)
    logical,intent(out) :: fail
    integer :: i,j,k
    integer :: lwork,info
    real(sp),allocatable :: work(:)
    integer,allocatable :: iwork(:)
    integer,allocatable :: ifail(:)
    real(sp) :: dum
    !> LAPACK
    external :: ssyevx
    fail = .false.
    lwork = 1+8*n+n**2
    allocate (iwork(5*n),work(lwork),ifail(n))
    j = 1
    call ssyevx('V','I','U',n,A,n,dum,dum,j,j,thr, &
    &           i,e,U,n,work,lwork,iwork,ifail,info)
    if (info .ne. 0) fail = .true.
    deallocate (iwork,work,ifail)
  end subroutine solver_ssyevx

!========================================================================================!

  subroutine solver_sspevx(n,thr,A,U,e,fail)
!*********************************************************
!* subroutine solver_sspevx
!* wrapper for LAPACK's sspevx routine:
!* SSPEVX computes all eigenvalues and eigenvectors of a
!* real symmetric matrix A in packed storage using a
!* divide and conquer algorithm.
!*********************************************************
    implicit none
    integer,intent(in) :: n
    real(sp),intent(in) :: thr
    real(sp),intent(inout) :: A(:)
    real(sp),intent(inout) :: U(:,:)
    real(sp),intent(inout) :: e(:)
    logical,intent(out) :: fail
    integer :: i,j,k
    integer :: info
    real(sp),allocatable :: work(:)
    integer,allocatable :: iwork(:)
    integer,allocatable :: ifail(:)
    real(sp) :: dum
    !> LAPACK
    external :: sspevx
    fail = .false.
    allocate (iwork(5*n),work(8*n),ifail(n))
    j = 1
    call sspevx('V','I','U',n,A,dum,dum,j,j,thr,i,e,U,n,work,iwork,ifail,info)
    if (info .ne. 0) fail = .true.
    deallocate (iwork,work,ifail)
  end subroutine solver_sspevx

!========================================================================================!
!========================================================================================!
end module optimize_maths
