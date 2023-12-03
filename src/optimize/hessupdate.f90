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
module hessupdate_module
  use iso_fortran_env,only:wp => real64

  public :: bfgs
  public :: powell
  public :: sr1
  public :: bofill
  public :: schlegel

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine bfgs(nat3,gnorm,grad,grado,dx,hess)
!**************************************************************************
!* subroutine bfgs
!* Performs BFGS update of Hessian matrix
!*
!*   Hₖ = Hₖ₋₁ + ΔHᵇᶠᵍˢ
!*
!* with
!*
!*              ΔgΔgᵀ     Hₖ₋₁ ΔxΔxᵀ Hₖ₋₁
!*   ΔHᵇᶠᵍˢ =   ─────  -  ───────────────
!*              ΔgᵀΔx       Δxᵀ Hₖ₋₁ Δx
!*
!* Input:
!* nat3        = dimension parameter as declared in the calling routine
!* grad        = gradient ( gₖ )
!* grado = gradient one cycle before ( gₖ₋₁ )
!* dx    = displ = Δx = xₖ - xₖ₋₁  ; k=cycle
!* hess        = Hessian matrix (Hₖ₋₁) on In- and updated Hessian (Hₖ) on Output
!**************************************************************************
    implicit none

    !> Input:
    integer,intent(in) :: nat3
    real(wp),intent(in) :: grad(nat3)
    real(wp),intent(in) :: grado(nat3)
    real(wp),intent(in) :: dx(nat3)
    real(wp),intent(in) :: gnorm
    !> Output:
    real(wp),intent(inout) :: hess(nat3*(nat3+1)/2)
    !> Local:
    integer  :: i,j,ij,ii
    real(wp),allocatable :: dg(:),Hdx(:)
    real(wp) :: dHBFGS
    real(wp) :: dxHdx,idxHdx,HdxdxHdx
    real(wp) :: dgdx,idgdx,dgdgdx
    real(wp),parameter :: thrs = 1.d-12
    real(wp),parameter :: thr = 1d-2
    !> BLAS:
    external :: dspmv
    real(wp),external :: ddot
    !---------------------------------------------------------------------
    allocate (dg(nat3),Hdx(nat3),source=0.0_wp)

    !> calculate Δg = gₖ - gₖ₋₁
    dg(1:nat3) = grad(1:nat3)-grado(1:nat3)

    !> calculate Hdx = Hₖ₋₁ * Δx
    call dspmv('u',nat3,1.0_wp,hess,dx,1,0.0_wp,Hdx,1)

    !> calculate dxHdx = Δxᵀ * Hₖ₋₁ * Δx
    dxHdx = ddot(nat3,Hdx,1,dx,1)

    !> calculate dgdx = Δgᵀ * Δx
    dgdx = ddot(nat3,dg,1,dx,1)

    !> inverse dxHdx and dgdx
    idxHdx = 1.0_wp/dxHdx
    idgdx = 1.0_wp/dgdx

    if (dgdx > thrs.and.dxHdx > thrs) then
      !$omp parallel default(none) &
      !$omp shared(nat3,idgdx,idxHdx,dg,Hdx) &
      !$omp private(i,j,ii,ij,dgdgdx,HdxdxHdx,dHBFGS) &
      !$omp shared(hess)
      !$omp do
      do i = 1,nat3
        !> Hessian index mapping
        ii = i*(i-1)/2

        !> calculate dgdgdx = Δg / (Δgᵀ * Δx)   (first index)
        dgdgdx = dg(i)*idgdx

        !> calculate HdxdxHdx = (Hₖ₋₁ * Δx) / (Δxᵀ * Hₖ₋₁ * Δx) (first index)
        HdxdxHdx = Hdx(i)*idxHdx

        do j = 1,i
          !> Hessian index mapping
          ij = ii+j

          !> calculate Hessian elements of ΔHᵇᶠᵍˢ (second indices)
          dHBFGS = dg(j)*dgdgdx-Hdx(j)*HdxdxHdx

          !> calculate updated Hessian Hₖ
          hess(ij) = hess(ij)+dHBFGS
        end do
      end do
      !$omp end do
      !$omp end parallel
    end if

    !> limit diagonal to (thr=0.01 slightly better than thr=0.001)
    ij = 0
    do i = 1,nat3
      ij = ij+i
      if (abs(hess(ij)) .lt. thr) hess(ij) = thr
    end do

    !> deallocate
    if (allocated(Hdx)) deallocate (Hdx)
    if (allocated(dg)) deallocate (dg)

    return
  end subroutine bfgs

!========================================================================================!

  subroutine powell(nat3,gnorm,grad,grado,dx,hess)
!**************************************************************************
!* subroutine powell:
!* Performs Powell-symmetric-Broyden update of Hessian matrix
!*
!*   Hₖ = Hₖ₋₁ + ΔHᵖˢᵇ
!*
!* with
!*
!*              ⎛ (Δg - Hₖ₋₁ Δx) Δxᵀ + Δx (Δg - Hₖ₋₁Δx)ᵀ⎞
!*   ΔHᵖˢᵇ  =   ⎜ ──────────────────────────────────────⎟
!*              ⎝                 Δxᵀ Δx                ⎠
!*
!*              ⎛(Δxᵀ (Δg - Hₖ₋₁Δx)) Δx Δxᵀ ⎞
!*            - ⎜────────────────────────── ⎟
!*              ⎝        (Δxᵀ Δx)²          ⎠
!*
!* Input:
!* nat3   = dimension parameter as declared in the calling routine=
!*         3*natoms
!* grad = gradient ( gₖ )
!* grado = gradient one cycle before ( gₖ₋₁ )
!* dx    = displ = Δx = xₖ - xₖ₋₁  ; k=cycle
!* hess = Hessian matrix (Hₖ₋₁) on In- and updated Hessian (Hₖ) on Output
!**************************************************************************
    implicit none
    !> Input:
    integer,intent(in) :: nat3
    real(wp),intent(in) :: grad(nat3),grado(nat3),dx(nat3),gnorm
    !> Output:
    real(wp),intent(inout) :: hess(nat3*(nat3+1)/2)
    !> Local:
    integer :: i,j,ij
    real(wp) :: dHPSB
    real(wp) :: dxdx,dxdgHdx
    real(wp),allocatable :: dgHdx(:)
    real(wp) :: thrs,scal,damp,dampO,dampD
    !> BLAS
    real(wp),external :: ddot
    !---------------------------------------------------------------------
    allocate (dgHdx(nat3),source=0.0_wp)
    thrs = 1.d-14

    !> calculate Hdx = Hₖ₋₁ * Δx
    call dspmv('u',nat3,1.0d0,hess,dx,1,0.0d0,dgHdx,1)

    !> calculate (Δg - Hₖ₋₁ Δx)
    dgHdx(1:nat3) = grad(1:nat3)-grado(1:nat3)-dgHdx(1:nat3)

    !> calculate scalar dxdx = ΔxᵀΔx
    dxdx = ddot(nat3,dx,1,dx,1)

    if (dxdx > thrs) then
      !> calculate (Δxᵀ (Δg - Hₖ₋₁Δx)) / (ΔxᵀΔx)
      dxdgHdx = ddot(nat3,dgHdx,1,dx,1)/dxdx

      do i = 1,nat3
        do j = 1,i
          !> Hessian index mapping
          ij = i*(i-1)/2+j

          !> calculate [(Δg - Hₖ₋₁ Δx) Δxᵀ + Δx (Δg - Hₖ₋₁Δx)ᵀ]
          dHPSB = dgHdx(i)*dx(j)+dx(i)*dgHdx(j)
          !> substract dxdgHdx*(ΔxΔxᵀ)
          dHPSB = dHPSB-dx(i)*dxdgHdx*dx(j)
          !> divide everything by scalar (ΔxᵀΔx)
          dHPSB = dHPSB/dxdx

          !> update Hessian by ΔHᵖˢᵇ
          hess(ij) = hess(ij)+dHPSB
        end do
      end do
    end if

    if (allocated(dgHdx)) deallocate (dgHdx)
    return
  end subroutine powell

!========================================================================================!

  subroutine sr1(nat3,gnorm,grad,grado,dx,hess)
!*************************************************************************
!* subroutine sr1
!* Performs symmetric rank one (SR1) update of Hessian matrix
!*
!*   Hₖ = Hₖ₋₁ + ΔHˢʳ¹
!*
!* with
!*
!*             (Δg - Hₖ₋₁ Δx)(Δg - Hₖ₋₁Δx)ᵀ
!*   ΔHˢʳ¹  =  ────────────────────────────
!*                  (Δg - Hₖ₋₁Δx)ᵀ Δx
!*
!* Input:
!* nat3        = dimension parameter as declared in the calling routine
!* grad        = gradient ( gₖ )
!* grado = gradient one cycle before ( gₖ₋₁ )
!* dx    = displ = Δx = xₖ - xₖ₋₁  ; k=cycle
!* hess        = Hessian matrix (Hₖ₋₁) on In- and updated Hessian (Hₖ) on Output
!*************************************************************************
    implicit none

    !> Input:
    integer,intent(in) :: nat3
    real(wp),intent(in) :: grad(nat3)
    real(wp),intent(in) :: grado(nat3)
    real(wp),intent(in) :: dx(nat3)
    real(wp),intent(in) :: gnorm
    !> Output:
    real(wp),intent(inout) :: hess(nat3*(nat3+1)/2)
    !> Local:
    integer  :: i,j,ij,ii
    real(wp),allocatable :: dg(:),dgHdx(:)
    real(wp) :: dHSR1
    real(wp) :: tddtd
    real(wp) :: dgHdxdx
    real(wp),parameter :: thrs = 1.d-12
    real(wp),parameter :: thr = 1d-2
    !> BLAS:
    external :: dspmv
    real(wp),external :: ddot
    !---------------------------------------------------------------------
    allocate (dg(nat3),dgHdx(nat3),source=0.0_wp)

    !> calculate Δg = gₖ - gₖ₋₁
    dg(1:nat3) = grad(1:nat3)-grado(1:nat3)

    !> calculate Hdx = Hₖ₋₁ * Δx
    call dspmv('u',nat3,1.0_wp,hess,dx,1,0.0_wp,dgHdx,1)
    !> calculate (Δg - Hₖ₋₁ Δx)
    dgHdx = dg-dgHdx

    !> calculate dgHdxdx = (Δg - Hₖ₋₁Δx)ᵀ Δx
    dgHdxdx = ddot(nat3,dgHdx,1,dx,1)

    if (abs(dgHdxdx) .gt. thrs) then
      !$omp parallel default(none) &
      !$omp shared(nat3,dgHdxdx,dg,dgHdx) &
      !$omp private(i,j,ii,ij,tddtd,dHSR1) &
      !$omp shared(hess)
      !$omp do
      do i = 1,nat3
        !> Hessian index mapping
        ii = i*(i-1)/2

        !> calculate tddtd = (Δg - Hₖ₋₁Δx)/ ((Δg - Hₖ₋₁Δx)ᵀ Δx) (first index)
        tddtd = dgHdx(i)/dgHdxdx

        do j = 1,i
          !> Hessian index mapping
          ij = ii+j

          !> calculate Hessian elements of ΔHˢʳ¹ (second index)
          dHSR1 = dgHdx(j)*tddtd

          !> calculate updated Hessian Hₖ
          hess(ij) = hess(ij)+dHSR1
        end do
      end do
      !$omp end do
      !$omp end parallel
    end if

    !> deallocate
    if (allocated(dgHdx)) deallocate (dgHdx)
    if (allocated(dg)) deallocate (dg)

    return
  end subroutine sr1

!========================================================================================!

  subroutine bofill(nat3,gnorm,grad,grado,dx,hess)
!****************************************************************************
!* subroutine bofill:
!* Performs Bofill update of Hessian matrix
!*
!*   Hₖ = Hₖ₋₁ + ΔHᵇᵒᶠⁱˡˡ
!*
!* with
!*
!*   ΔHᵇᵒᶠⁱˡˡ = φ ΔHˢʳ¹ + (1-φ) ΔHᵖˢᵇ
!*
!* and
!*         ⎛((Δg - Hₖ₋₁Δx)ᵀ Δx)²⎞
!*    φ =  ⎜────────────────────⎟
!*         ⎝|Δg - Hₖ₋₁Δx|² |Δx|²⎠
!*
!* Input:
!* nat3   = dimension parameter as declared in the calling routine=
!*         3*natoms
!* grad = gradient ( gₖ )
!* grado = gradient one cycle before ( gₖ₋₁ )
!* dx    = displ = Δx = xₖ - xₖ₋₁  ; k=cycle
!* hess = Hessian matrix (Hₖ₋₁) on In- and updated Hessian (Hₖ) on Output
!****************************************************************************
    implicit none
    !> Input:
    integer,intent(in) :: nat3
    real(wp),intent(in) :: grad(nat3),grado(nat3),dx(nat3),gnorm
    !> Output:
    real(wp),intent(inout) :: hess(nat3*(nat3+1)/2)
    !> Local:
    integer :: i,j,ij
    real(wp) :: dds,ddsd,ddtd,iddtd,dHPSB,dHSR1,dHBOFILL
    real(wp) :: ddh,ddhx,tddtd
    real(wp),allocatable :: dgHdx(:)
    real(wp) :: phi
    real(wp) :: thrs,scal,damp,dampO,dampD
    !> BLAS
    real(wp),external :: ddot
    !---------------------------------------------------------------------
    allocate (dgHdx(nat3),source=0.0_wp)
    thrs = 1.d-14

    !> calculate Hdx = Hₖ₋₁ * Δx
    call dspmv('u',nat3,1.0d0,hess,dx,1,0.0d0,dgHdx,1)

    !> calculate (Δg - Hₖ₋₁ Δx)
    dgHdx(1:nat3) = grad(1:nat3)-grado(1:nat3)-dgHdx(1:nat3)

    !> calculate scalar dds = ΔxᵀΔx = |Δx|²
    dds = ddot(nat3,dx,1,dx,1)

    !> calculate ddtd = (Δg - Hₖ₋₁Δx)ᵀ Δx
    ddtd = ddot(nat3,dgHdx,1,dx,1)
    !> inverse ddtd
    iddtd = 1.0_wp/ddtd

    !> calculate ((Δg - Hₖ₋₁Δx)ᵀ Δx)²
    ddhx = ddtd*ddtd
    !> calculate scalar |Δg - Hₖ₋₁Δx|²
    ddh = ddot(nat3,dgHdx,1,dgHdx,1)

    if (dds > thrs) then
      !> calculate (Δxᵀ (Δg - Hₖ₋₁Δx)) / (ΔxᵀΔx)
      ddsd = ddot(nat3,dgHdx,1,dx,1)/dds

      !> calculate φ
      phi = ddhx/(ddh*dds)

      do i = 1,nat3

        !> calculate tddtd = (Δg - Hₖ₋₁Δx)/ ((Δg - Hₖ₋₁Δx)ᵀ Δx) (first index)
        tddtd = dgHdx(i)*iddtd

        do j = 1,i
          !> Hessian index mapping
          ij = i*(i-1)/2+j

          !> calculate [(Δg - Hₖ₋₁ Δx) Δxᵀ + Δx (Δg - Hₖ₋₁Δx)ᵀ]
          dHPSB = dgHdx(i)*dx(j)+dx(i)*dgHdx(j)
          !> substract ddsd*(ΔxΔxᵀ)
          dHPSB = dHPSB-dx(i)*ddsd*dx(j)
          !> divide everything by scalar (ΔxᵀΔx)
          dHPSB = dHPSB/dds

          !> calculate Hessian elements of ΔHˢʳ¹ (second index)
          dHSR1 = dgHdx(j)*tddtd

          !> calculate ΔHᵇᵒᶠⁱˡˡ
          dHBOFILL = phi*dHSR1+(1.0_wp-phi)*dHPSB

          !> update Hessian by ΔHᵇᵒᶠⁱˡˡ
          hess(ij) = hess(ij)+dHBOFILL
        end do
      end do
    end if

    if (allocated(dgHdx)) deallocate (dgHdx)
    return
  end subroutine bofill

!========================================================================================!

  subroutine schlegel(nat3,gnorm,grad,grado,dx,hess)
!*************************************************************************
!* subroutine schlegel:
!* Performs Farkas-Schlegel update of Hessian matrix
!*
!*   Hₖ = Hₖ₋₁ + ΔHᶠˢ
!*
!* with
!*
!*   ΔHᶠˢ = φ^½ ΔHˢʳ¹ + (1-φ^½) ΔHᵇᶠᵍˢ
!*
!* and
!*           ⎛((Δg - Hₖ₋₁Δx)ᵀ Δx)²⎞½
!*    φ^½ =  ⎜────────────────────⎟
!*           ⎝|Δg - Hₖ₋₁Δx|² |Δx|²⎠
!*
!* Input:
!* nat3   = dimension parameter as declared in the calling routine=
!*         3*natoms
!* grad = gradient ( gₖ )
!* grado = gradient one cycle before ( gₖ₋₁ )
!* dx    = displ = Δx = xₖ - xₖ₋₁  ; k=cycle
!* hess = Hessian matrix (Hₖ₋₁) on In- and updated Hessian (Hₖ) on Output
!*************************************************************************
    implicit none
    !> Input:
    integer,intent(in) :: nat3
    real(wp),intent(in) :: grad(nat3),grado(nat3),dx(nat3),gnorm
    !> Output:
    real(wp),intent(inout) :: hess(nat3*(nat3+1)/2)
    !> Local:
    integer :: i,j,ij
    real(wp) :: dds,ddsd,ddtd,iddtd,dHBFGS,dHSR1,dHFS
    real(wp) :: ddh,ddhx,tddtd
    real(wp) :: dxHdx,idxHdx,dgdx,idgdx,dgdgdx,HdxdxHdx
    real(wp),allocatable :: dgHdx(:)
    real(wp),allocatable :: dg(:),Hdx(:)
    real(wp) :: phi
    real(wp) :: thrs,scal,damp,dampO,dampD
    !> BLAS
    real(wp),external :: ddot
    !---------------------------------------------------------------------
    allocate (dg(nat3),dgHdx(nat3),Hdx(nat3),source=0.0_wp)
    thrs = 1.d-14

    !> calculate Δg = gₖ - gₖ₋₁
    dg(1:nat3) = grad(1:nat3)-grado(1:nat3)

    !> calculate Hdx = Hₖ₋₁ * Δx
    call dspmv('u',nat3,1.0d0,hess,dx,1,0.0d0,Hdx,1)

    !> calculate (Δg - Hₖ₋₁ Δx)
    dgHdx(1:nat3) = dg(1:nat3)-Hdx(1:nat3)

    !> calculate scalar dds = ΔxᵀΔx = |Δx|²
    dds = ddot(nat3,dx,1,dx,1)

    !> calculate ddtd = (Δg - Hₖ₋₁Δx)ᵀ Δx
    ddtd = ddot(nat3,dgHdx,1,dx,1)
    !> inverse ddtd
    iddtd = 1.0_wp/ddtd

    !> calculate ddtd = Δxᵀ * Hₖ₋₁ * Δx
    dxHdx = ddot(nat3,Hdx,1,dx,1)

    !> calculate dds = Δgᵀ * Δx
    dgdx = ddot(nat3,dg,1,dx,1)

    !> inverse dxHdx and dgdx
    idxHdx = 1.0_wp/dxHdx
    idgdx = 1.0_wp/dgdx

    !> calculate ((Δg - Hₖ₋₁Δx)ᵀ Δx)²
    ddhx = ddtd*ddtd
    !> calculate scalar |Δg - Hₖ₋₁Δx|²
    ddh = ddot(nat3,dgHdx,1,dgHdx,1)

    if (dds > thrs) then
      !> calculate φ^½
      phi = ddhx/(ddh*dds)
      phi = sqrt(phi)

      do i = 1,nat3

        !> calculate tddtd = (Δg - Hₖ₋₁Δx)/ ((Δg - Hₖ₋₁Δx)ᵀ Δx) (SR1 first index)
        tddtd = dgHdx(i)*iddtd

        !> calculate dgdgdx = Δg / (Δgᵀ * Δx)   (BFGS first index)
        dgdgdx = dg(i)*idgdx
        !> calculate HdxdxHdx = (Hₖ₋₁ * Δx) / (Δxᵀ * Hₖ₋₁ * Δx) (BFGS first index)
        HdxdxHdx = Hdx(i)*idxHdx

        do j = 1,i
          !> Hessian index mapping
          ij = i*(i-1)/2+j

          !> calculate Hessian elements of ΔHᵇᶠᵍˢ (BFGS second indices)
          dHBFGS = dg(j)*dgdgdx-Hdx(j)*HdxdxHdx

          !> calculate Hessian elements of ΔHˢʳ¹ (SR1 second index)
          dHSR1 = dgHdx(j)*tddtd

          !> calculate ΔHᶠˢ
          dHFS = phi*dHSR1+(1.0_wp-phi)*dHBFGS

          !> update Hessian by ΔHᶠˢ
          hess(ij) = hess(ij)+dHFS
        end do
      end do
    end if

    if (allocated(Hdx)) deallocate (Hdx)
    if (allocated(dgHdx)) deallocate (dgHdx)
    if (allocated(dg)) deallocate (dg)
    return
  end subroutine schlegel

!========================================================================================!
!========================================================================================!
end module hessupdate_module
