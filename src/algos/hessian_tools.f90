!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Gereon Feldmann
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
! Routines were adapted from the xtb code (github.com/grimme-lab/xtb)
! under the Open-source software LGPL-3.0 Licencse.
!================================================================================!

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> Routines for the computation of a projected mass-weighted Hessian
!> Routines for the computation of frequencies from the Hessian
!> Rotuines for the computation of the effective Hessian at an MECP based on: DOI:10.1039/A907723E
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!

module hessian_tools
  use crest_parameters, only:wp,stdout
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
  use optimize_maths
  use iomod

  public :: frequencies

!=========================================================================================!
!=========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!=========================================================================================!
!=========================================================================================!


  subroutine frequencies(nat,at,xyz,nat3,calc,prj_mw_hess,freq,io)
!*************************************************
!* Returns the Frequencies from a Hessian in cm-1
!*************************************************
    implicit none

    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(calcdata) :: calc
    real(wp) :: prj_mw_hess(nat3,nat3)

    integer :: io,nat3
    logical :: pr
    real(wp) :: energy
    real(wp) :: freq(nat3)
    real(wp),allocatable :: pmode(:,:)

    integer,allocatable :: iwork(:)
    real(wp),allocatable :: work(:)

    integer :: lwork,liwork,info,i
    !>LAPCK
    external :: dsyevd

    nat3 = nat*3

    !Parameters for diagonalization
    lwork = 1+6*nat3+2*nat3**2
    liwork = 3+5*nat3

    allocate (work(lwork),iwork(liwork))

    !Diagonalization
    call dsyevd('V','U',nat3,prj_mw_hess,nat3,freq,work,lwork,iwork,liwork,info)

    deallocate (work,iwork)

    !Convert eigenvalues to frequencies
    do i = 1,nat3
      if (freq(i) .gt. 0.0_wp) then
        freq(i) = sqrt(freq(i))*219474.63_wp
      else
        freq(i) = -sqrt(abs(freq(i)))*219474.63_wp
      end if
    end do

    return

  end subroutine frequencies


  subroutine mass_weight_hess(nat,at,nat3,hess)
    use atmasses
    implicit none

    !> Mass weighting the Hessian
    integer,intent(in) :: nat                   !Number of atoms
    integer,intent(in) :: at(nat)               !atomic number of all atoms

    real(wp),intent(inout) :: hess(nat3,nat3)   !Hessian matrix
    real(wp) :: mass_in_au             !Masses of all atoms of the periodic table
    integer :: i,j,nat3,i3,i33,j3,j33

    mass_in_au = (1.66054e-27_wp/9.1094e-31_wp)**2

    !amv = ams(1:118)

    do i = 1,nat
      do j = i,nat

        i3 = 3*(i-1)+1
        i33 = 3*(i-1)+3
        j3 = 3*(j-1)+1
        j33 = 3*(j-1)+3

        hess(i3:i33,j3:j33) = 1/sqrt(ams(at(i))*ams(at(j))*mass_in_au)*hess(i3:i33,j3:j33)
        !Hessian is symmetric hence upper triangular can be copied
        hess(j3:j33,i3:i33) = hess(i3:i33,j3:j33)

      end do
    end do

    return
  end subroutine mass_weight_hess

  subroutine dsqtoh(n,a,b)
!****************************************************
!* converts upper triangle of a matrix into a vector
!****************************************************
    implicit none
    integer,intent(in)  :: n
    real(wp),intent(in) :: a(n,n)
    real(wp),intent(out)  :: b(n*(n+1)/2)
    integer :: i,j,k

    k = 0
    do i = 1,n
      do j = 1,i
        k = k+1
        b(k) = a(i,j)
      end do
    end do

  end subroutine dsqtoh

  subroutine dhtosq(n,a,b)
!*********************************************************
!* converts upper triangle vector into a symmetric matrix
!*********************************************************
    implicit none
    integer,intent(in)  :: n
    real(wp),intent(out) :: a(n,n)
    real(wp),intent(in)  :: b(n*(n+1)/2)
    integer :: i,j,k

    k = 0
    do i = 1,n
      do j = 1,i
        k = k+1
        a(j,i) = b(k)
        a(i,j) = b(k)
      end do
    end do

    return
  end subroutine dhtosq

!=========================================================================================!

  subroutine prj_mw_hess(nat,at,nat3,xyz,hess)
!***************************************************************
!* Projection of the translational and rotational DOF out of 
!* the numerical Hessian plus the mass-weighting of the Hessian
!***************************************************************
    implicit none

    integer,intent(in) :: nat,nat3
    integer :: at(nat)
    real(wp),intent(inout) :: hess(nat3,nat3)
    real(wp) ::  xyz(3,nat)
    !real(wp) ::  hess_ut(nat3*(nat3+1)/2),pmode(nat3,1)
    real(wp),allocatable ::  hess_ut(:),pmode(:,:)

    allocate(hess_ut(nat3*(nat3+1)/2), source=0.0_wp)
    allocate(pmode(nat3,1), source=0.0_wp)

    !> Transforms matrix of the upper triangle vector
    call dsqtoh(nat3,hess,hess_ut)

    !> Projection
    call trproj(nat,nat3,xyz,hess_ut,.false.,0,pmode,1)

    !> Transforms vector of the upper triangle into matrix
    call dhtosq(nat3,hess,hess_ut)

    !> Mass weighting
    call mass_weight_hess(nat,at,nat3,hess)

    deallocate(pmode,hess_ut)
  end subroutine prj_mw_hess

!=========================================================================================!

  subroutine print_vib_spectrum(nat,at,nat3,xyz,freq,dir,fname)
!*********************************************************************
!* Prints the frequencies in Turbomoles "vibspectrum" format 
!* The intensity is only artficially set to 1000 for every vibration!!
!**********************************************************************
    integer,intent(in) :: nat,nat3
    integer :: at(nat),i
    real(wp) ::  xyz(3,nat)
    real(wp) ::  freq(nat3),thr
    character(len=*) :: fname
    character(len=*) :: dir

    thr = 0.01_wp
    if (len_trim(dir) .eq. 0) then
      open (newunit=ich,file=fname)
    else
      if(directory_exist(dir))then
      open (newunit=ich,file=dir//'/'//fname)
      else
      open (newunit=ich,file=fname)
      endif
    end if


    write (ich,'("$vibrational spectrum")')
    write (ich,'("#  mode    symmetry    wave number    IR intensity    selection rules")')
    write (ich,'("#                       1/cm              km/mol         IR    RAMAN")')

    do i = 1,nat3
      if (abs(freq(i)) .lt. thr) then
        write (ich,'(i6,9x,    f18.2,f16.5,7x," - ",5x," - ")') &
          i,freq(i),0.0_wp
      else
        write (ich,'(i6,8x,"a",f18.2,f16.5,7x,"YES",5x,"YES")') &
          i,freq(i),1000.0_wp
      end if
    end do

    write (ich,'("$end")')

    close (ich)

  end subroutine print_vib_spectrum

!=========================================================================================!

  subroutine print_g98_fake(nat,at,nat3,xyz,freq,hess,dir,fname)
!****************************************************************
!* Prints the vibration spectrum of the a system as a g98.out.
!* Routine is adapted from the xtb code.
!****************************************************************
    integer,intent(in) :: nat,nat3
    integer :: at(nat)
    integer  :: gu,i,j,ka,kb,kc,la,lb,k

    real(wp) ::  xyz(3,nat)
    real(wp),intent(in) :: hess(nat3,nat3)
    real(wp) ::  freq(nat3),red_mass(nat3),force(nat3),ir_int(nat3),zero(1),f2(nat3),u(nat3,nat3)

    character(len=2) :: irrep
    character(len=*) :: fname
    character(len=*) :: dir

    irrep = 'a'

    red_mass = 99.0
    force = 99.0
    ir_int = 99.0
    zero = 0.0

    k = 0

    do i = 1,nat3
      if (abs(freq(i)) .gt. 1.d-1) then
        k = k+1
        u(1:nat3,k) = hess(1:nat3,i)
        f2(k) = freq(i)
      end if
    end do

    if (len_trim(dir) .eq. 0) then
      open (newunit=gu,file=fname)
    else
      if(directory_exist(dir))then
      open (newunit=gu,file=dir//'/'//fname)
      else
      open (newunit=gu,file=fname) 
      endif
    end if

    write (gu,'('' Entering Gaussian System'')')
    write (gu,'('' *********************************************'')')
    write (gu,'('' Gaussian 98:'')')
    write (gu,'('' frequency output generated by the crest code'')')
    write (gu,'('' *********************************************'')')

    write (gu,*) '                        Standard orientation:'
    write (gu,*) '---------------------------------------------', &
        & '-----------------------'
    write (gu,*) ' Center     Atomic     Atomic', &
        & '              Coordinates (Angstroms)'
    write (gu,*) ' Number     Number      Type ', &
        & '             X           Y           Z'
    write (gu,*) '-----------------------', &
        & '---------------------------------------------'
    j = 0
    do i = 1,nat
      write (gu,111) i,at(i),j,xyz(1:3,i)*0.52917726
    end do
    write (gu,*) '----------------------', &
        & '----------------------------------------------'
    write (gu,*) '    1 basis functions        1 primitive gaussians'
    write (gu,*) '    1 alpha electrons        1 beta electrons'
    write (gu,*)
111 format(i5,i11,i14,4x,3f12.6)

    write (gu,*) 'Harmonic frequencies (cm**-1), IR intensities',' (km*mol⁻¹),'
    write (gu,*) 'Raman scattering activities (A**4/amu),', &
        & ' Raman depolarization ratios,'
    write (gu,*) 'reduced masses (AMU), force constants (mDyne/A)', &
        & ' and normal coordinates:'

    ka = 1
    kc = 3

60  kb = min0(kc,k)
    write (gu,100) (j,j=ka,kb)
    write (gu,105) (irrep,j=ka,kb)
    write (gu,110) ' Frequencies --', (f2(j),j=ka,kb)
    write (gu,110) ' Red. masses --', (red_mass(j),j=ka,kb)
    write (gu,110) ' Frc consts  --', (force(j),j=ka,kb)
    write (gu,110) ' IR Inten    --', (ir_int(j),j=ka,kb)
    write (gu,110) ' Raman Activ --', (zero,j=ka,kb)
    write (gu,110) ' Depolar     --', (zero,j=ka,kb)
    write (gu,*) 'Atom AN      X      Y      Z        X      Y', &
        & '      Z        X      Y      Z'
    la = 1
70  lb = nat
    do i = la,lb
      write (gu,130) i,at(i), (u(i*3-2,j),u(i*3-1,j),u(i*3,j),j=ka,kb)
    end do
    if (lb .eq. nat) go to 90
    go to 70
90  if (kb .eq. k) then
      return
    end if

    ka = kc+1
    kc = kc+3
    go to 60

100 format(3(20x,i3))
105 format(3x,3(18x,a5))
110 format(a15,f11.4,12x,f11.4,12x,f11.4)
130 format(2i4,3(2x,3f7.2))

    write (gu,'(''end of file'')')
    close (gu)

  end subroutine print_g98_fake

!=========================================================================================!


  subroutine print_hessian(hess,nat3,dir,fname)
!*******************************
!* Prints the numerical hessian
!*******************************
    integer :: nat3,i,j,k
    real(wp) :: hess(nat3,nat3)
    character(len=*) :: fname
    character(len=*) :: dir

    if (len_trim(dir) .eq. 0) then
      open (newunit=ich,file=fname)
      write (stdout,'(1x,a)',advance='no') 'Will be written to file "'//fname//'" ...'
    else
      if(directory_exist(dir))then
      open (newunit=ich,file=dir//'/'//fname)
      write (stdout,'(1x,a)',advance='no') 'Will be written to file "'//dir//'/'//fname//'" ...'
      else
      open (newunit=ich,file=fname)
      write (stdout,'(1x,a)',advance='no') 'Will be written to file "'//fname//'" ...'
      endif
    end if
    flush(stdout)

    write (ich,'(1x,a)') '$hessian'
    do i = 1,nat3
      k = 0
      do j = 1,nat3
        k = k+1
        if (k .le. 4) then
          write (ich,'(f16.8)',advance='no') hess(i,j)
        else
          write (ich,'(f16.8)') hess(i,j)
          k = 0
        end if
      end do
      if (k .ne. 0) then
        write (ich,*)
      end if
    end do
    write (ich,'(1x,a)') '$end'
    close (ich)

    write (stdout,*) 'done.'
    write (stdout,*)

  end subroutine print_hessian

!=========================================================================================!

  subroutine effective_hessian(nat,nat3,grad1_i,grad2_i,hess1,hess2,heff)
!******************************************************************
!* Effective Hessian at an MECP is computed via Eq. 27 and Eq. 28 
!* in https://doi.org/10.1002/qua.25124
!******************************************************************
    implicit none
    integer,intent(in) :: nat,nat3
    integer :: i,j,ii
    real(wp),intent(in) :: grad1_i(3,nat3),grad2_i(3,nat3)
    real(wp) :: grad1(nat3),grad2(nat3),dot

    real(wp),intent(in) :: hess1(nat3,nat3),hess2(nat3,nat3)

    real(wp) :: gnorm1,gnorm2,grad_diff_norm
    real(wp) :: grad_diff(nat3),heff_temp(nat3,nat3)

    real(wp),intent(inout) :: heff(nat3,nat3)
    real(wp),allocatable :: proj_vec(:,:)

    real(wp) :: freq(nat3)

    integer,allocatable :: iwork(:)
    real(wp),allocatable :: work(:)

    integer :: lwork,liwork,info

    allocate (proj_vec(nat3,nat3),source=0.0_wp)

    grad1 = reshape(grad1_i, (/nat3/))
    grad2 = reshape(grad2_i, (/nat3/))

    gnorm1 = norm2(grad1)

    gnorm2 = norm2(grad2)

    grad_diff = grad1-grad2

    grad_diff_norm = norm2(grad_diff)

    dot = dot_product(grad1,grad2)

    if (dot .gt. 0.0_wp) then !sloped: dot > 0.0 --> -  | peaked: dot <= 0.0 --> +

      write (stdout,*) 'MECI is considered as a sloped CI'
      write (stdout,*)

      heff = (gnorm1*hess2-gnorm2*hess1)/grad_diff_norm

    else

      write (stdout,*) 'MECI is considered as a peaked CI'
      write (stdout,*)

      heff = (gnorm1*hess2+gnorm2*hess1)/grad_diff_norm

    end if

    !Outer Product of grad_diff

    !Building projection matrix

    !proj_vec = 1 - (dg/|dg| o dg.T/|dg|) = 1 - (dg o dg.T)/|dg|**2

    grad_diff_norm = grad_diff_norm**2

    do i = 1,nat3
      proj_vec(i,:) = -grad_diff(i)*grad_diff/grad_diff_norm
      proj_vec(i,i) = proj_vec(i,i)+1
    end do

    !Projection
    heff = matmul(matmul(proj_vec,heff),proj_vec)

    !Check if hess1 and hess2 are assigned correctly, otherwise change
    lwork = 1+6*nat3+2*nat3**2
    liwork = 3+5*nat3
    allocate (work(lwork),iwork(liwork))

    heff_temp = heff

    call dsyevd('V','U',nat3,heff_temp,nat3,freq,work,lwork,iwork,liwork,info)

    deallocate (work,iwork)

    if (0 .gt. sum(freq)) then
      heff = -heff
    end if

  end subroutine effective_hessian

!=========================================================================================!

 subroutine calculate_frequencies(calc,nat,at,xyz,freq,io,constraints)
!*******************************************************
!* Bundels several routines from this module to
!* calculate the vib. frequencies for a given structure
!* The output frequencies are in cm-1 
!*******************************************************
     implicit none
     !> INPUT
     type(calcdata),intent(inout) :: calc
     integer,intent(in)  :: nat
     integer,intent(in)  :: at(nat)
     real(wp),intent(in) :: xyz(3,nat)
     logical,intent(in),optional :: constraints
     !> OUTPUT
     integer,intent(out)  :: io
     real(wp),intent(out) :: freq(3*nat)
     !> LOCAL
     real(wp),allocatable :: hessian(:,:)
     real(wp),allocatable :: chess(:,:)
     type(calcdata) :: dummycalc
     integer :: nat3,ncalc,i

     io = 0
     freq = 0.0_wp
     nat3 = nat*3
     ncalc = calc%ncalculations

     allocate(hessian(nat3,nat3), source = 0.0_wp)

     !>--- Hessian from combined energy and gradient
     call numhess1(nat,at,xyz, calc,hessian,io) 
     if( io /= 0 ) return

     !>--- do we consider contributions from the constraints?
     !>    (yes, by default, they are in the hessian from numhess1,
     !>    if we DO NOT want them, we need to take them out again)
     if(present(constraints))then
       if(.not.constraints)then
         dummycalc = calc !> new dummy calculation
         dummycalc%id = 0  !> set to zero so that ONLY constraints are considered
         dummycalc%ncalculations = 0
         dummycalc%pr_energies = .false.
         allocate (chess(nat3,nat3),source=0.0_wp)
         call numhess1(nat,at,xyz,dummycalc,chess,io)
         hessian(:,:) = hessian(:,:) - chess(:,:)  
         deallocate( chess )
       endif
     endif
     
     do i = 1,calc%ncalculations

      !>-- Projects and mass-weights the Hessian
      call prj_mw_hess(nat,at,nat3,xyz, hessian(:,:))

      !>-- Computes the Frequencies
      call frequencies(nat,at,xyz,nat3, calc, hessian(:,:), freq(:),io)
     end do

     deallocate( hessian )
     return
 end subroutine calculate_frequencies

!=========================================================================================!
end module hessian_tools
