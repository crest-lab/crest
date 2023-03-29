!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2023 Philipp Pracht, Stefan Grimme
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
!================================================================================!

!=========================================================================================!
! estimate conformational flexibility from geom. and WBO
! the value of flex is normalized i.e. between 0(rigid) and 1(long alkane)
!=========================================================================================!
subroutine flexi(mol,rednat,includeAtom,flex)
  use iso_fortran_env,only:wp => real64
  use strucrd
  use zdata,only:readwbo
  implicit none
  !> INPUT
  type(coord),intent(in) :: mol                  !> input molecule
  integer,intent(in)     :: rednat               !> should equal sum(includeAtom(:))
  integer,intent(in)     :: includeAtom(mol%nat) !> 0 if not considered, 1 if considered
  !> OUTPUT
  real(wp),intent(out) :: flex
  !> LOCAL
  real(wp),allocatable::xyz(:,:),rcov(:),cn(:),cring(:,:),wbo(:,:),wbofull(:,:)
  integer,allocatable::at(:),map(:),b(:,:),sring(:),map2(:)
  real(wp) :: dx,dy,dz,r,r2,rco,val,ringf,doublef,branch,effectivNat,hybf
  real(wp) :: av2
  integer :: i,j,k,l,m,n,rn
  logical :: ex

  !> CN less than this is considerd as terminating atom (H, F, ...)
  real(wp),parameter :: thr = 1.2_wp

  flex = 0.0_wp

  n = mol%nat !> local modification
  rn = rednat !> Number of atoms selected by atomlist+/-
  !> rn=n if atomlist+/- is not used

  allocate (xyz(3,rn),at(rn),rcov(94),map(rn),cn(rn),b(rn,rn))
  allocate (sring(rn),cring(12,rn),wbo(rn,rn),wbofull(n,n),map2(n))
  call setrcov(rcov)

!>--- map given structure to considered substructure
  if (rn .ne. n) then
    call rdcoord_reduced('coord',n,rn,xyz,at,includeAtom)
  else
    !call rdcoord('coord',rn,at,xyz)
    at = mol%at
    xyz = mol%xyz
  end if
  call ncoord(rn,rcov,at,xyz,cn,400.0d0)

!>--- map the new (reduced) coordinate order to the original
  j = 0
  map2 = 0   !> every EXcluded atom should get a zero
  do i = 1,n
    if (includeAtom(i) .gt. 0) then
      j = j+1
      map2(i) = j
    end if
  end do

!>--- read the WBO file and map it to the (reduced) coordinates
!>--- if it does not exist, single and double bonds are not distinguished
  wbo = 0
  wbofull = 0
  inquire (file='wbo',exist=ex)
  if (ex) then
    call readwbo('wbo',n,wbofull)

!>--- map wbofull to wbo (only if atoms are excluded from the rmsd)
    if (rn .ne. n) then
      do i = 1,n-1
        if (includeAtom(i) .gt. 0) then
          k = map2(i)
          do j = i+1,n
            if (includeAtom(j) .gt. 0) then
              l = map2(j)
              wbo(k,l) = wbofull(i,j)
              wbo(l,k) = wbo(k,l)
            else
              cycle !> j cycle
            end if
          end do
        else
          cycle !> i cycle
        end if
      end do
    else !> rn == n
      wbo = wbofull
    end if !>-- if(rn.ne.n) end
  end if !>-- if(ex) end
  deallocate (map2,wbofull)

!+++++++++++++++++++++++++++++++++++++++!
  n = rn !<---- so I don't have to change the following original code
!+++++++++++++++++++++++++++++++++++++++!

!>--- adjacency matrix b setup
  map = 0
  b(:,:) = 0
  do i = 1,n
    if (cn(i) .lt. thr) cycle
    do j = 1,n
      if (i .eq. j) cycle
      dx = xyz(1,j)-xyz(1,i)
      dy = xyz(2,j)-xyz(2,i)
      dz = xyz(3,j)-xyz(3,i)
      r2 = dx*dx+dy*dy+dz*dz
      r = sqrt(r2)
      rco = rcov(at(i))+rcov(at(j))
      if (r .lt. rco) then
        if (cn(j) .gt. thr) then
          b(j,i) = 1
          map(i) = map(i)+1
        end if
      end if
    end do
  end do

!>--- setup ring data for all atoms, i.e., in which ring and how large
  call minringsizes(n,at,xyz,sring)

!>--- the actual flexibility measure calculation
  m = 0
  av2 = 0.0d0
  do i = 1,n
    do j = 1,i-1
      !> count all "bonds"
      if (b(j,i) .gt. 0) m = m+1

      !> no branch on terminating bonds (e.g. Me)
      if (map(i) .eq. 1.or.map(j) .eq. 1) cycle

      !> adjacent pair ij
      if (b(j,i) .gt. 0) then

        hybf = 1.0

        !> sp2 C are less flexible
        if (at(i) .eq. 6.and.cn(i) .lt. 3.3) hybf = hybf*0.5_wp
        if (at(j) .eq. 6.and.cn(j) .lt. 3.3) hybf = hybf*0.5_wp

        !> double bond term
        doublef = 1.0-exp(-4.0*(wbo(j,i)-2.0)**6)

        !> branching measure
        branch = 2.0/sqrt(dble(map(i))*dble(map(j)))

        !> flex. reduction for atoms in rings
        ringf = 1.0
        k = min(sring(i),sring(j))
        !> a ring is (even at infinite size) a factor of 2 less flexible, factors are empirical
        if (k .gt. 0) ringf = 0.5*(1.0-exp(-0.06*float(k)))
        !> (adjusted such that c-C20 is converged)

        !>  put it together
        val = branch*ringf*doublef*hybf
        !> quadratic av
        av2 = av2+val**2
      end if
    end do
  end do
  if (m .gt. 0) then
    av2 = sqrt(av2/dble(m))
  end if

  flex = av2
  effectivNat = av2*dble(n)
  deallocate (xyz,at,rcov,map,cn,b,sring,cring,wbo)

  return
end subroutine flexi
subroutine minringsizes(nat,at,xyz,sring)
  use crest_parameters
  use zdata
  implicit none
  integer :: nat
  integer :: at(nat)
  real(wp) :: xyz(3,nat)
  integer :: sring(nat)
  type(zmolecule) :: zmol
  integer :: i,j

  call simpletopo(nat,at,xyz,zmol,.false.,.true.,'')
  do i = 1,zmol%nat
    sring(i) = 0
    do j = 1,zmol%nri
      if (any(zmol%zri(j)%rlist == i)) then
        if (sring(i) == 0) then
          sring(i) = zmol%zri(j)%rs
        else if (zmol%zri(j)%rs < sring(i)) then
          sring(i) = zmol%zri(j)%rs
        end if
      end if
    end do
  end do
  call zmol%deallocate()
  return
end subroutine minringsizes
!========================================================================================!

!> new version of the non-covalent flexibility measure.
!> the nci-flexibility is estimated RELATIVE to a bioorganic
!> molecule (crambin) and its hydrogen-bond and dispersion
!> energy per atom.
subroutine nciflexi_gfnff(mol,flexval)
  use crest_parameters
  use strucrd
  use calc_type
  use gfnff_api
  implicit none
  !> INPUT
  type(coord),intent(in) :: mol
  !> OUTPUT
  real(wp),intent(out)   :: flexval
  !> LOCAL
  integer :: io
  logical :: ex
  real(wp) :: ehb,edisp
  real(wp) :: energy
  real(wp),allocatable :: grad(:,:)
  type(gfnff_data),allocatable :: ff_dat
!>--- reset
  flexval = 0.0_wp
  ehb = 0.0_wp
  edisp = 0.0_wp

!>--- preprocessor statement (only executed when compiled with GFN-FF support)
#ifdef WITH_GFNFF
  allocate (ff_dat)
  allocate (grad(3,mol%nat),source=0.0_wp)

!>--- set up a GFN-FF calculation
  write (stdout,'(1x,a)',advance='no') 'Calculating NCI flexibility...'
  flush (stdout)
  call gfnff_api_setup(mol,mol%chrg,ff_dat,io,.false.)
  if (io /= 0) then
    write (stdout,'(a)') ' failed.'
    deallocate (grad,ff_dat)
    return
  end if

!>--- perform a singlepoint calculation
  call gfnff_sp(mol,ff_dat,energy,grad,io)
  if (io /= 0) then
    write (stdout,'(a)') ' failed.'
    deallocate (grad,ff_dat)
    return
  end if
  ehb = ff_dat%res%e_hb
  edisp = ff_dat%res%e_disp
  write (stdout,'(a)') ' done.'

!>--- normalize by number of atoms
  ehb = ehb/mol%nat
  edisp = edisp/mol%nat

!>--- NCI flexi is determined RELATIVE to a reference molecule (Crambin)
  flexval = 0.5_wp*(1.0_wp-(ehb/(-0.00043374_wp)))
  flexval = flexval+0.5_wp*(1.0_wp-(edisp/(-0.00163029_wp)))

  deallocate (grad,ff_dat)
#else

  write (stdout,'(a)') 'Caclulation of NCI flexibility skipped.'
  write (stdout,'(a)') 'Set up compilation with -DWITH_GFNFF=true to enable this feature.'
  flexval = 0.0_wp

#endif

  return
end subroutine nciflexi_gfnff
!========================================================================================!
