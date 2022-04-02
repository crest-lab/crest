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
   use iso_fortran_env, only: wp=>real64

   public :: bfgs
   public :: powell

contains
!========================================================================================!
!> subroutine bfgs
!> Performs BFGS update of Hessian matrix
!>
!>   Hₖ = Hₖ₋₁ + ΔHᵇᶠᵍˢ
!>
!> with 
!> 
!>              ΔgΔgᵀ     Hₖ₋₁ ΔxΔxᵀ Hₖ₋₁ 
!>   ΔHᵇᶠᵍˢ =   ─────  -  ───────────────
!>              ΔgᵀΔx       Δxᵀ Hₖ₋₁ Δx
!>
!> Input:
!> nat3	= dimension parameter as declared in the calling routine
!> grad	= gradient ( gₖ )
!> grado = gradient one cycle before ( gₖ₋₁ )
!> dx    = displ = Δx = xₖ - xₖ₋₁  ; k=cycle
!> hess	= Hessian matrix (Hₖ₋₁) on In- and updated Hessian (Hₖ) on Output
!>------------------------------------------------------------------------
subroutine bfgs(nat3,gnorm,grad,grado,dx,hess)
   implicit none

   !> Input:
   integer, intent(in) :: nat3
   real(wp),intent(in) :: grad(nat3)
   real(wp),intent(in) :: grado(nat3)
   real(wp),intent(in) :: dx(nat3)
   real(wp),intent(in) :: gnorm
   !> Output:
   real(wp),intent(inout) :: hess(nat3*(nat3+1)/2)
   !> Local:
   integer  :: i,j,ij,ii
   real(wp),allocatable :: dg(:), Hdx(:)
   real(wp) :: ddtd, dds, dHBFGS
   real(wp) :: iddtd, idds, sdds, tddtd
   real(wp),parameter :: thrs=1.d-12
   real(wp),parameter :: thr=1d-2
   !> BLAS:
   external :: dspmv
   real(wp),external :: ddot
   !---------------------------------------------------------------------
   allocate( dg(nat3), Hdx(nat3), source = 0.0_wp )

   !> calculate Δg = gₖ - gₖ₋₁
   dg(1:nat3) = grad(1:nat3) - grado(1:nat3)

   !> calculate Hdx = Hₖ₋₁ * Δx
   call dspmv('u',nat3,1.0_wp,hess,dx,1,0.0_wp,Hdx,1)

   !> calculate ddtd = Δxᵀ * Hₖ₋₁ * Δx
   ddtd = ddot(nat3,Hdx,1,dx,1)

   !> calculate dds = Δgᵀ * Δx
   dds  = ddot(nat3,dg,1,dx,1)

   !> inverse ddtd and dds 
   iddtd = 1.0_wp / ddtd
   idds  = 1.0_wp / dds

   if(dds > thrs .and. ddtd > thrs) then
   !$omp parallel default(none) &
   !$omp shared(nat3,idds,iddtd,dg,Hdx) &
   !$omp private(i,j,ii,ij,sdds,tddtd,dHBFGS) &
   !$omp shared(hess)
   !$omp do
      do i=1,nat3
         !> Hessian index mapping
         ii = i*(i-1)/2     
      
         !> calculate ssds = Δg / (Δgᵀ * Δx)   (first index)
         sdds  = dg(i)*idds

         !> calculate tddtd = (Hₖ₋₁ * Δx) / (Δxᵀ * Hₖ₋₁ * Δx) (first index)
         tddtd = Hdx(i)*iddtd

         do j=1,i
            !> Hessian index mapping
            ij = ii + j

            !> calculate Hessian elements of ΔHᵇᶠᵍˢ (second indices)
            dHBFGS = dg(j)*sdds - Hdx(j)*tddtd

            !> calculate updated Hessian Hₖ
            hess(ij) = hess(ij) + dHBFGS
         end do
      end do
   !$omp end do
   !$omp end parallel
   endif

   !> limit diagonal to (thr=0.01 slightly better than thr=0.001)
   ij=0
   do i=1,nat3
      ij=ij+i
      if(abs(hess(ij)).lt.thr)hess(ij)=thr
   enddo

   !> deallocate
   if(allocated(Hdx))deallocate(Hdx)
   if(allocated(dg))deallocate(dg)

   return
end subroutine bfgs

!========================================================================================!
!> subroutine powell:
!> Performs Powell-symmetric-Broyden update of Hessian matrix
!>
!>   Hₖ = Hₖ₋₁ + ΔHᵖˢᵇ
!>
!> with 
!> 
!>              ⎛ (Δg - Hₖ₋₁ Δx) Δxᵀ + Δx (Δg - Hₖ₋₁Δx)ᵀ⎞
!>   ΔHᵖˢᵇ  =   ⎜ ──────────────────────────────────────⎟
!>              ⎝                 Δxᵀ Δx                ⎠
!>             
!>              ⎛(Δxᵀ (Δg - Hₖ₋₁Δx)) Δx Δxᵀ ⎞
!>            - ⎜────────────────────────── ⎟
!>              ⎝        (Δxᵀ Δx)²          ⎠
!>
!> Input:
!> nat3   = dimension parameter as declared in the calling routine=
!>         3*natoms
!> grad = gradient ( gₖ )
!> grado = gradient one cycle before ( gₖ₋₁ )
!> dx    = displ = Δx = xₖ - xₖ₋₁  ; k=cycle
!> hess = Hessian matrix (Hₖ₋₁) on In- and updated Hessian (Hₖ) on Output
!>--------------------------------------------------------------------
subroutine powell(nat3,gnorm,grad,grado,dx,hess)
   implicit none
   !> Input:
   integer, intent(in) :: nat3
   real(wp),intent(in) :: grad(nat3),grado(nat3),dx(nat3),gnorm
   !> Output:
   real(wp),intent(inout) :: hess(nat3*(nat3+1)/2)
   !> Local:
   integer :: i,j,ij
   real(wp) :: dds,ddtd,dHPSB
   real(wp),allocatable :: dgHdx(:)
   real(wp) :: thrs, scal, damp, dampO,dampD
   !> BLAS
   real(wp),external :: ddot
   !---------------------------------------------------------------------
   allocate( dgHdx(nat3), source=0.0_wp)
   thrs=1.d-14

   !> calculate Hdx = Hₖ₋₁ * Δx
   call dspmv('u',nat3,1.0d0,hess,dx,1,0.0d0,dgHdx,1)

   !> calculate (Δg - Hₖ₋₁ Δx)
   dgHdx(1:nat3) = grad(1:nat3) - grado(1:nat3) - dgHdx(1:nat3)

   !> calculate scalar dds = ΔxᵀΔx
   dds  = ddot(nat3,dx,1,dx,1)

   if(dds > thrs) then
      !> calculate (Δxᵀ (Δg - Hₖ₋₁Δx)) / (ΔxᵀΔx)
      ddtd = ddot(nat3,dgHdx,1,dx,1)/dds

      do i=1,nat3
         do j=1,i
            !> Hessian index mapping
            ij = i*(i-1)/2 + j

            !> calculate [(Δg - Hₖ₋₁ Δx) Δxᵀ + Δx (Δg - Hₖ₋₁Δx)ᵀ]
            dHPSB=dgHdx(i)*dx(j) + dx(i)*dgHdx(j)
            !> substract ddtd*(ΔxΔxᵀ)
            dHPSB=dHPSB - dx(i)*ddtd*dx(j)
            !> divide everything by scalar (ΔxᵀΔx)
            dHPSB = dHPSB/dds

            !> update Hessian by ΔHᵖˢᵇ 
            hess(ij) = hess(ij) + dHPSB
         end do
      end do
   endif
   
   if(allocated(dgHdx))deallocate(dgHdx)
   return
end subroutine powell

subroutine hdamp(gnorm,dampO,dampD)
   implicit none
   real(wp) gnorm,dampO,dampD
   dampD = 1.0
   dampO = 1.0
   return
   ! damping of H update
   if(gnorm > 0.5)then
      dampD = 1.0
      dampO = 0.0
   endif
   if(gnorm < 0.5 .and. gnorm > 0.2)then
      dampD = 1.0
      dampO = 0.2
   endif
   if(gnorm < 0.2 .and. gnorm > 0.05) then
      dampD = 1.0
      dampO = 0.5
   endif
end subroutine hdamp


!========================================================================================!
!> subroutine sr1
!> Performs symmetric rank one (SR1) update of Hessian matrix
!>
!>   Hₖ = Hₖ₋₁ + ΔHˢʳ¹
!>
!> with 
!> 
!>             (Δg - Hₖ₋₁ Δx)(Δg - Hₖ₋₁Δx)ᵀ 
!>   ΔHˢʳ¹  =  ────────────────────────────
!>                  (Δg - Hₖ₋₁Δx)ᵀ Δx
!>
!> Input:
!> nat3	= dimension parameter as declared in the calling routine
!> grad	= gradient ( gₖ )
!> grado = gradient one cycle before ( gₖ₋₁ )
!> dx    = displ = Δx = xₖ - xₖ₋₁  ; k=cycle
!> hess	= Hessian matrix (Hₖ₋₁) on In- and updated Hessian (Hₖ) on Output
!>------------------------------------------------------------------------
subroutine sr1(nat3,gnorm,grad,grado,dx,hess)
   implicit none

   !> Input:
   integer, intent(in) :: nat3
   real(wp),intent(in) :: grad(nat3)
   real(wp),intent(in) :: grado(nat3)
   real(wp),intent(in) :: dx(nat3)
   real(wp),intent(in) :: gnorm
   !> Output:
   real(wp),intent(inout) :: hess(nat3*(nat3+1)/2)
   !> Local:
   integer  :: i,j,ij,ii
   real(wp),allocatable :: dg(:), dgHdx(:)
   real(wp) :: ddtd, dHSR1
   real(wp) :: iddtd, tddtd
   real(wp),parameter :: thrs=1.d-8
   real(wp),parameter :: thr=1d-3
   !> BLAS:
   external :: dspmv
   real(wp),external :: ddot
   !---------------------------------------------------------------------
   allocate( dg(nat3), dgHdx(nat3), source = 0.0_wp )

   !> calculate Δg = gₖ - gₖ₋₁
   dg(1:nat3) = grad(1:nat3) - grado(1:nat3)

   !> calculate Hdx = Hₖ₋₁ * Δx
   call dspmv('u',nat3,1.0_wp,hess,dx,1,0.0_wp,dgHdx,1)
   !> calculate (Δg - Hₖ₋₁ Δx)
   dgHdx = dg - dgHdx

   !> calculate ddtd = (Δg - Hₖ₋₁Δx)ᵀ Δx
   ddtd = ddot(nat3,dgHdx,1,dx,1)
   !> inverse ddtd  
   iddtd = 1.0_wp / ddtd

!   if(abs(ddtd).gt.thrs)then
   !$omp parallel default(none) &
   !$omp shared(nat3,iddtd,dg,dgHdx) &
   !$omp private(i,j,ii,ij,tddtd,dHSR1) &
   !$omp shared(hess)
   !$omp do
      do i=1,nat3
         !> Hessian index mapping
         ii = i*(i-1)/2     
      
         !> calculate tddtd = (Δg - Hₖ₋₁Δx)/ ((Δg - Hₖ₋₁Δx)ᵀ Δx) (first index)
         tddtd = dgHdx(i)*iddtd

         do j=1,i
            !> Hessian index mapping
            ij = ii + j

            !> calculate Hessian elements of ΔHˢʳ¹ (second index)
            dHSR1 = dgHdx(j)*tddtd

            !> calculate updated Hessian Hₖ
            hess(ij) = hess(ij) + dHSR1
         end do
      end do
   !$omp end do
   !$omp end parallel
!   endif

   !> limit diagonal
!   ij=0
!   do i=1,nat3
!      ij=ij+i
!      if(abs(hess(ij)).lt.thr)hess(ij)=thr
!   enddo

   !> deallocate
   if(allocated(dgHdx))deallocate(dgHdx)
   if(allocated(dg))deallocate(dg)

   return
end subroutine sr1


!========================================================================================!
!> subroutine bofill:
!> Performs Bofill update of Hessian matrix
!>
!>   Hₖ = Hₖ₋₁ + ΔHᵇᵒᶠⁱˡˡ
!>
!> with 
!> 
!>   ΔHᵇᵒᶠⁱˡˡ = φ ΔHˢʳ¹ + (1-φ) ΔHᵖˢᵇ 
!>
!> and             
!>         ⎛((Δg - Hₖ₋₁Δx)ᵀ Δx)²⎞
!>    φ =  ⎜────────────────────⎟
!>         ⎝|Δg - Hₖ₋₁Δx|² |Δx|²⎠
!>
!> Input:
!> nat3   = dimension parameter as declared in the calling routine=
!>         3*natoms
!> grad = gradient ( gₖ )
!> grado = gradient one cycle before ( gₖ₋₁ )
!> dx    = displ = Δx = xₖ - xₖ₋₁  ; k=cycle
!> hess = Hessian matrix (Hₖ₋₁) on In- and updated Hessian (Hₖ) on Output
!>--------------------------------------------------------------------
subroutine bofill(nat3,gnorm,grad,grado,dx,hess)
   implicit none
   !> Input:
   integer, intent(in) :: nat3
   real(wp),intent(in) :: grad(nat3),grado(nat3),dx(nat3),gnorm
   !> Output:
   real(wp),intent(inout) :: hess(nat3*(nat3+1)/2)
   !> Local:
   integer :: i,j,ij
   real(wp) :: dds,ddsd,ddtd,iddtd,dHPSB,dHSR1,dHBOFILL
   real(wp) :: ddh,ddhx,tddtd
   real(wp),allocatable :: dgHdx(:)
   real(wp) :: phi
   real(wp) :: thrs, scal, damp, dampO,dampD
   !> BLAS
   real(wp),external :: ddot
   !---------------------------------------------------------------------
   allocate( dgHdx(nat3), source=0.0_wp)
   thrs=1.d-14

   !> calculate Hdx = Hₖ₋₁ * Δx
   call dspmv('u',nat3,1.0d0,hess,dx,1,0.0d0,dgHdx,1)

   !> calculate (Δg - Hₖ₋₁ Δx)
   dgHdx(1:nat3) = grad(1:nat3) - grado(1:nat3) - dgHdx(1:nat3)

   !> calculate scalar dds = ΔxᵀΔx = |Δx|²
   dds  = ddot(nat3,dx,1,dx,1)

   !> calculate ddtd = (Δg - Hₖ₋₁Δx)ᵀ Δx
   ddtd = ddot(nat3,dgHdx,1,dx,1)
   !> inverse ddtd  
   iddtd = 1.0_wp / ddtd

   !> calculate ((Δg - Hₖ₋₁Δx)ᵀ Δx)²
   ddhx = ddtd * ddtd
   !> calculate scalar |Δg - Hₖ₋₁Δx|²
   ddh = ddot(nat3,dgHdx,1,dgHdx,1)

   if(dds > thrs) then
      !> calculate (Δxᵀ (Δg - Hₖ₋₁Δx)) / (ΔxᵀΔx)
      ddsd = ddot(nat3,dgHdx,1,dx,1)/dds
   
      !> calculate φ
      phi = ddhx / ( ddh * dds )

      do i=1,nat3

         !> calculate tddtd = (Δg - Hₖ₋₁Δx)/ ((Δg - Hₖ₋₁Δx)ᵀ Δx) (first index)
         tddtd = dgHdx(i)*iddtd

         do j=1,i
            !> Hessian index mapping
            ij = i*(i-1)/2 + j

            !> calculate [(Δg - Hₖ₋₁ Δx) Δxᵀ + Δx (Δg - Hₖ₋₁Δx)ᵀ]
            dHPSB=dgHdx(i)*dx(j) + dx(i)*dgHdx(j)
            !> substract ddsd*(ΔxΔxᵀ)
            dHPSB=dHPSB - dx(i)*ddsd*dx(j)
            !> divide everything by scalar (ΔxᵀΔx)
            dHPSB = dHPSB/dds

            !> calculate Hessian elements of ΔHˢʳ¹ (second index)
            dHSR1 = dgHdx(j)*tddtd

            !> calculate ΔHᵇᵒᶠⁱˡˡ
            dHBOFILL = phi*dHSR1 + (1.0_wp - phi)*dHPSB  

            !> update Hessian by ΔHᵇᵒᶠⁱˡˡ
            hess(ij) = hess(ij) + dHBOFILL
         end do
      end do
   endif
   
   if(allocated(dgHdx))deallocate(dgHdx)
   return
end subroutine bofill

!========================================================================================!
end module hessupdate_module
