!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2017-2020 Stefan Grimme, Sebastian Ehlert (xtb)
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
module wiberg_mayer
  use iso_fortran_env,only:wp => real64
  implicit none

  public :: get_wbo_uhf,get_wbo_rhf,density_matrix
  public :: split_foccab,write_wbo

  public :: get_wbo
  interface get_wbo
    module procedure get_wbo_rhf
    module procedure get_wbo_uhf
  end interface get_wbo

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine get_wbo_rhf(nat,nao,P,S,ao2at,wbo)
!******************************************************************
!* Calculate the Wiberg/Mayer type bond order for
!* a given density matrix and overlap matrix.
!* Also requires the information which AO belongs to what atom.
!* The (Mayer) bond order between atoms A and B is calculated as
!*
!*   Bᴬᴮ = ΣₛΣₜ(PS)ₛₜ(PS)ₜₛ  with s∈ A, t∈ B
!*
!* for the Wiberg BO, S=I, so Bᴬᴮ=ΣΣP²
!* this version is only valid for closed shell cases
!******************************************************************
    implicit none
    !> Input
    integer,intent(in)  :: nat        !> number of atoms
    integer,intent(in)  :: nao        !> number of AOs
    real(wp),intent(in) :: P(nao,nao) !> Density matrix
    real(wp),intent(in) :: S(nao,nao) !> Overlap matrix
    integer,intent(in)  :: ao2at(nao) !> mapping: which AO belongs to what atom
    !> Output
    real(wp),intent(out) :: wbo(nat,nat) !> Bond order matrix
    !> Locals
    integer :: iao,jao,iat,jat
    real(wp) :: pao
    real(wp),allocatable :: PS(:,:)
    !> BLAS
    external :: dgemm

    allocate (PS(nao,nao),source=0.0_wp)

    call dgemm('N','N',nao,nao,nao,1.0d0,P,nao,S,nao,0.0d0,PS,nao)

    !$omp parallel do default(none) collapse(2) &
    !$omp shared(nao, ao2at, PS, wbo) private(iao, jao, iat, jat, pao)
    do iao = 1,nao
      do jao = 1,nao
        iat = ao2at(iao)
        jat = ao2at(jao)
        pao = merge(PS(iao,jao)*PS(jao,iao),0.0_wp,iat /= jat)
        !$omp atomic
        wbo(jat,iat) = wbo(jat,iat)+pao
      end do
    end do

    deallocate (PS)
  end subroutine get_wbo_rhf

!========================================================================================!

  subroutine get_wbo_uhf(nat,nao,Pa,Pb,S,ao2at,wbo)
!*******************************************************************
!* Calculate the Wiberg/Mayer type bond order for
!* a given density matrix and overlap matrix.
!* Also requires the information which AO belongs to what atom.
!* The (Mayer) bond order between atoms A and B is calculated as
!*
!*   Bᴬᴮ = 2 ΣₛΣₜ(PaS)ₛₜ(PaS)ₜₛ+ (PbS)ₛₜ(PbS)ₜₛ  with s∈ A, t∈ B
!*
!* for the Wiberg BO, S=I, so Bᴬᴮ=ΣΣP²
!* this version is valid for closed shell and open cases
!******************************************************************
    implicit none
    !> Input
    integer,intent(in)  :: nat         !> number of atoms
    integer,intent(in)  :: nao         !> number of AOs
    real(wp),intent(in) :: Pa(nao,nao) !> Density matrix (alpha)
    real(wp),intent(in) :: Pb(nao,nao) !> Density matrix (beta)
    real(wp),intent(in) :: S(nao,nao)  !> Overlap matrix
    integer,intent(in)  :: ao2at(nao)  !> mapping: which AO belongs to what atom
    !> Output
    real(wp),intent(out) :: wbo(nat,nat) !> Bond order matrix
    !> Locals
    integer :: iao,jao,iat,jat
    real(wp) :: pao
    real(wp),allocatable :: PS(:,:)
    !> BLAS
    external :: dgemm

    allocate (PS(nao,nao),source=0.0_wp)

    wbo = 0.0_wp
    ! alpha
    call dgemm('N','N',nao,nao,nao,1.0d0,Pa,nao,S,nao,0.0d0,PS,nao)

    !$omp parallel do default(none) collapse(2) &
    !$omp shared(nao, ao2at, PS, wbo) private(iao, jao, iat, jat, pao)
    do iao = 1,nao
      do jao = 1,nao
        iat = ao2at(iao)
        jat = ao2at(jao)
        pao = merge(PS(iao,jao)*PS(jao,iao),0.0_wp,iat /= jat)
        !$omp atomic
        wbo(jat,iat) = wbo(jat,iat)+pao
      end do
    end do

    ! beta
    PS = 0.0_wp
    call dgemm('N','N',nao,nao,nao,1.0d0,Pb,nao,S,nao,0.0d0,PS,nao)

    !$omp parallel do default(none) collapse(2) &
    !$omp shared(nao, ao2at, PS, wbo) private(iao, jao, iat, jat, pao)
    do iao = 1,nao
      do jao = 1,nao
        iat = ao2at(iao)
        jat = ao2at(jao)
        pao = merge(PS(iao,jao)*PS(jao,iao),0.0_wp,iat /= jat)
        !$omp atomic
        wbo(jat,iat) = wbo(jat,iat)+pao
      end do
    end do

    wbo = wbo*2.0_wp

    deallocate (PS)
  end subroutine get_wbo_uhf

!========================================================================================!

  subroutine write_wbo(wbo,cutoff,filename)
!**********************************************************
!* write WBO file
!* optional arguments are a file name and a cut-off value
!**********************************************************
    implicit none
    !> INPUT
    real(wp),intent(in) :: wbo(:,:)
    real(wp),intent(in),optional :: cutoff
    character(len=*),intent(in),optional :: filename
    !> LOCAL
    integer :: iunit,i,j,nat

    if (present(filename)) then
      open (newunit=iunit,file=trim(filename))
    else
      open (newunit=iunit,file='wbo')
    end if

    nat = size(wbo,1)
    if (present(cutoff)) then
      do i = 1,nat
        do j = 1,i-1
          if (wbo(j,i) > cutoff) then
            write (iunit,*) i,j,wbo(j,i)
          end if
        end do
      end do
    else
      do i = 1,nat
        do j = 1,i-1
          write (iunit,*) i,j,wbo(j,i)
        end do
      end do
    end if

    close (iunit)
  end subroutine write_wbo

!========================================================================================!

  subroutine density_matrix(nao,focc,C,P)
!****************************
!* compute density matrix
!*  C: MO coefficient matix
!*  focc: occupation matrix
!*  P  density matrix
!***************************
    integer,intent(in)  :: nao
    real(wp),intent(in)  :: focc(:)
    real(wp),intent(in)  :: C(:,:)
    real(wp),intent(out) :: P(:,:)
    integer :: i,m
    real(wp),allocatable :: Ptmp(:,:)
    !> LAPACK
    external :: dgemm,daxpy

    allocate (Ptmp(nao,nao))
    Ptmp = 0.0_wp
    ! Ptmp = C * n
    do i = 1,nao
      call daxpy(nao,focc(i),C(:,i),1,Ptmp(:,i),1)
    end do
    ! P = C * Ptmp^T = C * n * C^T
    call dgemm('N','T',nao,nao,nao,1.0d0,C,nao,Ptmp,nao,0.0d0,P,nao)

    deallocate (Ptmp)

  end subroutine density_matrix

!========================================================================================!

  subroutine split_foccab(nao,focc,nela,nelb,focca,foccb)
!*************************************************************
!* Split the occupation tensor focc into 
!* the alpha and beta occupations focca and foccb, xtb-style
!* nao: number of AOs
!* nela,nelb: total number of alpha and beta electrons
!* focc: occupation tensor
!* focca,foccb: alpha and beta occupations
!*************************************************************
    implicit none
    integer,intent(in) :: nao
    real(wp),intent(in) :: focc(nao)
    real(wp),intent(in) :: nela,nelb
    real(wp),intent(out),allocatable :: focca(:),foccb(:)
    integer :: i,j
    real(wp) :: nel
    allocate(focca(nao),foccb(nao), source=0.0_wp)
    nel = nela 
    do i=1,nao
      if(nel > 1.0_wp)then
       focca(i) = 1.0_wp
       nel = nel - 1.0_wp
      else
       focca(i) = nel
       exit 
      endif
    enddo
    nel = nelb
    do i=1,nao
      if(nel > 1.0_wp)then
       foccb(i) = 1.0_wp
       nel = nel - 1.0_wp
      else
       foccb(i) = nel
       exit
      endif
    enddo
  end subroutine split_foccab

!========================================================================================!
!========================================================================================!
end module wiberg_mayer
