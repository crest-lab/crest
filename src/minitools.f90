!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020 Philipp Pracht
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

subroutine splitfile(fname,up,low)
!********************************************************
!* A small subroutine to split an ensemble file into
!* seperate directories. A max. number of files can be
!* be specified.
!********************************************************
  use crest_parameters
  use iomod
  use strucrd,only:rdensemble,coord
  implicit none
  character(len=*) :: fname
  integer :: up,low

  character(len=512) :: thispath,tmppath1,tmppath2

  real(wp),allocatable :: xyz(:,:,:)
  integer :: nat,nall
  integer :: nc
  integer,allocatable :: at(:)
  integer :: i,r
  logical :: ex
  type(coord),allocatable :: structures(:)

  inquire (file=fname,exist=ex)
  if (.not.ex) then
    write (0,'(a,a,a)') "file ",trim(fname)," does not exist. must stop"
    error stop
  end if

  call getcwd(thispath) !current dir= thispath

  call rdensemble(fname,nall,structures)

  r = makedir("SPLIT")  !create new directory
  call chdir("SPLIT")
  call getcwd(tmppath1)

  if (up .gt. 0.and.low .gt. 0) then
    if (low .gt. up) then
      i = low
      low = up
      up = i
    end if
  end if
  ! if there are less structures in the file than we need
  if (nall < up) then
    nc = nall
  else
    nc = up
  end if

  do i = low,nc
    write (tmppath2,'(a,i0)') "STRUC",i
    r = makedir(trim(tmppath2))
    call chdir(tmppath2)
    call structures(i)%write("struc.xyz")
    call chdir(tmppath1)
  end do

  call chdir(thispath)
  write (stdout,*)
  write (stdout,'(a,i0,a,i0,a)') '>  Created subdirectories SPLIT/STRUC{',low,'-',nc,'}/'
  write (stdout,'(a)') '>  All directories contain a "struc.xyz" with molecular coordinates'
  write (stdout,'(a)') '>  The order of SPLIT/STRUC*/ is the same as in '//trim(fname)
  write (stdout,*)
  write (stdout,'(a)') 'exit.'
  return
end subroutine splitfile

!=========================================================================================!

subroutine printaniso(fname,bmin,bmax,bshift)
!****************************************************
!* print the anisotropy of the rotational constants
!* for all structures in a given ensemble file
!****************************************************
  use crest_parameters
  use strucrd
  use axis_module
  use rotaniso_mod
  implicit none
  character(len=*) :: fname
  type(coord),allocatable :: structures(:)

  integer :: nat
  integer :: nall
  real(wp),allocatable :: c1(:,:)
  integer,allocatable :: at(:)

  real(wp),allocatable :: rot(:,:)
  real(wp),allocatable :: anis(:)
  real(wp) :: evec(3,3),evecavg(3,3)

  real(wp) :: bmin,bmax,bshift
  real(wp) :: thr
  real(wp) :: dum
  integer :: i

  call rdensemble(fname,nall,structures)
  nat = structures(1)%nat

  allocate (c1(3,nat),at(nat))
  allocate (rot(3,nall))
  allocate (anis(nall))

  at(:) = structures(1)%at(:)
  evecavg(:,:) = 0.0_wp
  do i = 1,nall
    c1(1:3,:) = structures(i)%xyz(1:3,:)*autoaa
    call axis(nat,at,c1,rot(1:3,i),dum,evec)
    evecavg(:,:) = evecavg(:,:)+evec(:,:)
  end do
  evecavg(:,:) = evecavg(:,:)/real(nall)

  do i = 1,nall
    c1(1:3,:) = structures(i)%xyz(1:3,:)*autoaa
    call axis(nat,at,c1,rot(1:3,i),dum,evec)
    anis(i) = rotaniso(rot(1:3,i))
    thr = bthrerf(bmin,anis(i),bmax,bshift)
    write (*,'(3f10.2,2x,f8.4,2x,f8.4)') rot(1:3,i),anis(i),thr
  end do

  deallocate (anis,rot,at,c1)

  stop
  return
end subroutine printaniso

!=========================================================================================!

subroutine rotalign_tool(fname)
!****************************************************
!* print the anisotropy of the rotational constants
!* for all structures in a given ensemble file
!****************************************************
  use crest_parameters
  use strucrd
  use axis_module
  implicit none
  character(len=*) :: fname
  type(coord),allocatable :: structures(:)

  integer :: nat
  integer :: nall
  real(wp),allocatable :: c1(:,:),c2(:,:)
  integer,allocatable :: at(:)

  real(wp),allocatable :: rot(:,:)
  real(wp) :: rotaniso !function
  real(wp),allocatable :: anis(:)
  real(wp) :: evec(3,3),evecavg(3,3)

  real(wp) :: bthrerf
  real(wp) :: bmin,bmax,bshift
  real(wp) :: thr
  real(wp) :: dum
  integer :: i

  real(wp),parameter :: Ry90(3,3) = reshape([ &
                  &    0.0_wp,0.0_wp,-1.0_wp, &
                  &    0.0_wp,1.0_wp,0.0_wp,  &
                  &    1.0_wp,0.0_wp,0.0_wp   &
                  &    ], [3,3])

  call rdensemble(fname,nall,structures)
  nat = structures(1)%nat

  allocate (c1(3,nat),c2(3,nat),at(nat))
  allocate (rot(3,nall))

  at(:) = structures(1)%at(:)
  write (*,'(3a10)') 'A/MHz','B/MHz','C/MHz'
  do i = 1,nall
    c1(1:3,:) = structures(i)%xyz(1:3,:)*autoaa
    call axis(nat,at,c1,rot(1:3,i))
    c2 = c1
    structures(i)%xyz = c2*aatoau
    write (*,'(3f10.2)') rot(1:3,i)
  end do

  deallocate (rot,at,c2,c1)

  call wrensemble('rotalign.xyz',nall,structures)

  stop
  return
end subroutine rotalign_tool

!=========================================================================================!

subroutine prbweight(fname,Targ)
!*****************************************************
!* read in a file with 1 to 2 columns
!* first column is an energy (in Eh)
!* second column can be a degeneracy
!* a Boltzman weight is calculated for each
!* energy and a average ensemble energy is returned
!* comment lines (#) are ignored
!*****************************************************
  use crest_parameters
  use iomod, only: readl
  implicit none

  character(len=*) :: fname
  character(len=*) :: Targ

  real(wp) :: dum,T
  logical :: ex
  integer :: io,ich

  character(len=128) :: atmp
  integer :: l,j,i
  real(wp),allocatable :: elist(:)
  real(wp),allocatable :: erel(:)
  real(wp),allocatable :: g(:)
  real(wp),allocatable :: p(:)
  real(wp) :: xx(10)
  real(wp) :: elow

  !-- parse temperature from argument
  Targ = trim(adjustl(Targ))
  read (Targ,*,iostat=io) dum
  if (io == 0) then
    T = dum
  else
    T = 298.15d0
  end if

  inquire (file=fname,exist=ex)
  if (.not.ex) stop

  open (newunit=ich,file=fname)
  l = 0
  do
    read (ich,'(a)',iostat=io) atmp
    if (io < 0) exit !EOF
    atmp = trim(adjustl(atmp))
    if (len_trim(atmp) .lt. 1) cycle
    if (atmp(1:1) == '#') cycle
    l = l+1
  end do
  close (ich)

  allocate (elist(l),erel(l),source=0.0_wp)
  allocate (g(l),source=1.0_wp)
  allocate (p(l),source=0.0_wp)

  open (newunit=ich,file=fname)
  l = 0
  do
    read (ich,'(a)',iostat=io) atmp
    if (io < 0) exit !EOF
    atmp = trim(adjustl(atmp))
    if (len_trim(atmp) .lt. 1) cycle
    if (atmp(1:1) == '#') cycle
    l = l+1
    call readl(atmp,xx,j)
    elist(l) = xx(1)
    if (j > 1) g(l) = nint(xx(2))
  end do
  close (ich)

  dum = minval(elist,1)
  elow = dum
  do i = 1,l
    erel(i) = (elist(i)-dum)*627.5905_wp
  end do
  call entropy_boltz(l,T,erel,g,p)

  atmp = trim(fname)//'.out'
  open (newunit=ich,file=atmp)
  write (ich,'(a,16x,a,2x,a)') '#Etot','degen','pop'
  write (*,'(a,16x,a,2x,a)') '#Etot','degen','pop'
  do i = 1,l
    write (ich,'(f20.8,1x,i0,4x,f10.6)') elist(i),nint(g(i)),p(i)
    write (*,'(f20.8,1x,i0,4x,f10.6)') elist(i),nint(g(i)),p(i)
  end do
  close (ich)

  write (*,*)
  write (*,'(1x,a,f6.4)') 'pop sum: ',sum(p)
  do i = 1,l
    elist(i) = elist(i)*p(i)
  end do

  write (*,*) 'E_low:',elow
  write (*,*) 'E_av:',sum(elist)
  stop
end subroutine prbweight

!=========================================================================================!

subroutine calceav(fname,T,eav,verbose)
  use crest_parameters
  use strucrd
  implicit none
  character(len=*) :: fname
  real(wp) :: eav
  real(wp) :: T
  logical :: verbose
  type(ensemble) :: ens
  real(wp),allocatable :: elist(:)
  real(wp),allocatable :: erel(:)
  real(wp),allocatable :: g(:)
  real(wp),allocatable :: p(:)
  real(wp) :: dum
  integer :: i
  call ens%open(fname)
  allocate (elist(ens%nall),erel(ens%nall),g(ens%nall),p(ens%nall))
  elist = ens%er
  erel = 0.0d0
  g = 1.0d0
  p = 0.0d0
  if (verbose) then
    write (*,'(1x,a,1x,a)') 'Energies for file:',trim(fname)
    do i = 1,ens%nall
      write (*,'(1x,i6,1x,f20.10)') i,elist(i)
    end do
  end if
  dum = minval(elist,1)
  do i = 1,ens%nall
    erel(i) = (elist(i)-dum)*627.5905_wp
  end do
  call entropy_boltz(ens%nall,T,erel,g,p)
  do i = 1,ens%nall
    elist(i) = elist(i)*p(i)
  end do
  eav = sum(elist)
  if (verbose) then
    write (*,'(1x,a,1x,f20.10)') 'Weighted energy:',eav
  end if
  deallocate (p,g,erel,elist)
  return
end subroutine calceav

!=========================================================================================!

subroutine getelow(fname,elow,verbose)
  use crest_parameters
  use strucrd
  implicit none
  character(len=*) :: fname
  real(wp) :: elow
  logical :: verbose
  type(ensemble) :: ens
  real(wp),allocatable :: elist(:)
  integer :: i

  call ens%open(fname)
  allocate (elist(ens%nall))
  elist = ens%er
  if (verbose) then
    write (*,'(1x,a,1x,a)') 'Energies for file:',trim(fname)
    do i = 1,ens%nall
      write (*,'(1x,i6,1x,f20.10)') i,elist(i)
    end do
  end if
  elow = minval(elist,1)
  if (verbose) then
    write (*,'(1x,a,1x,f20.10)') 'Lowest energy:',elow
  end if
  deallocate (elist)
  call ens%deallocate()
  return
end subroutine getelow

!=========================================================================================!

subroutine testtopo(fname,env,tmode)
!************************************************
!* set up the topology for a file and analyze it
!************************************************
  use crest_parameters
  use crest_data
  use iomod
  use atmasses
  use zdata
  use strucrd
  use molbuilder_classify
  implicit none
  type(systemdata) :: env
  character(len=*) :: fname
  character(len=:),allocatable :: wbofile
  character(len=*) :: tmode
  type(zmolecule) :: zmol
  type(coord) :: mol
  type(coord_classify) :: molc
  real(wp),allocatable :: xyz(:,:)
  real(wp) :: dum
  integer,allocatable :: inc(:)
  real(wp) :: flex
  integer :: nt,i,maxl,j,ich
  logical :: l1
  real(wp),allocatable :: temps(:)
  real(wp),allocatable :: et(:)
  real(wp),allocatable :: ht(:)
  real(wp),allocatable :: gt(:)
  real(wp),allocatable :: stot(:)
  integer,allocatable :: amat(:,:)
  integer,allocatable :: path(:)

  call to_lower(tmode)
!--- specify zmol sertup
  select case (tmode)
  case ('wbo','flexi','all')
    call xtbsp2(fname,env)
    call simpletopo_file(fname,zmol,.true.,.true.,'wbo')
  case default
    wbofile = 'none'
    call simpletopo_file(fname,zmol,.true.,.true.,wbofile)
  end select
  allocate (xyz(3,zmol%nat))
  call zmol%getxyz(xyz)
  mol%nat = zmol%nat
  mol%at = zmol%at
  mol%xyz = xyz
  xyz = xyz*bohr   !to angstroem
!--- specify other analysis
  write (*,*)
  select case (tmode)
  case ('sym','symmetry')
    call analsym(zmol,dum,.true.)

  case ('flexi')
    allocate (inc(zmol%nat),source=1)
    call flexi(mol,zmol%nat,inc,flex)
    write (*,'(1x,a,4x,f6.4)') 'flexibility measure:',flex
    deallocate (inc)

  case ('zmat')
    call ztopozmat(zmol,.true.)

  case ('formula','sumform')
    write (*,'(/,1x,a)') trim(sumform(zmol%nat,zmol%at))
    write (*,'(1x,a,i16)') '# atoms: ',zmol%nat
    write (*,'(1x,a,f16.5)') 'Mol.weight: ',molweight(zmol%nat,zmol%at)

  case ('all')
    call ztopozmat(zmol,.true.)
    write (*,'(/,1x,a)') trim(sumform(zmol%nat,zmol%at))
    write (*,'(1x,a,i16)') '# atoms: ',zmol%nat
    write (*,'(1x,a,f16.5)') 'Mol.weight: ',molweight(zmol%nat,zmol%at)
    call analsym(zmol,dum,.true.)
    allocate (inc(zmol%nat),source=1)
    call flexi(mol,zmol%nat,inc,flex)
    write (*,'(1x,a,4x,f6.4)') 'flexibility measure:',flex
    deallocate (inc)

  case ('thermo')
    if (.not.allocated(env%thermo%temps)) then
      call env%thermo%get_temps()
    end if
    nt = env%thermo%ntemps
    allocate (temps(nt),et(nt),ht(nt),gt(nt),stot(nt))
    temps = env%thermo%temps

    if (.not.env%legacy.and.env%calc%ncalculations == 0) then
      call env2calc_setup(env)
    end if
    call thermo_wrap(env,.true.,zmol%nat,zmol%at,xyz,'', &
    &    nt,temps,et,ht,gt,stot,.false.)
    deallocate (stot,gt,ht,et,temps)

  case ('methyl')
    do i = 1,zmol%nat
      l1 = zmol%methyl(i)
      if (l1) write (*,'(a,i0,a)') 'Atom ',i,' is methyl (or similar)'
    end do

  case ('stereo')
    call isstereo(zmol)

  case ('back','backchain')
    if (.not.allocated(Amat)) then
      allocate (Amat(zmol%nat,zmol%nat),source=0)
    end if
    if (.not.allocated(path)) allocate (path(zmol%nat),source=0)
    call zmol%adjacency(Amat)
    call backchain(zmol%nat,Amat,zmol%at,maxl,path)
    open (newunit=ich,file='backchain.xyz')
    write (ich,*) maxl
    write (ich,*)
    do i = 1,maxl
      j = path(i)
      write (ich,'(a,3F25.15)') i2e(mol%at(j)),mol%xyz(:,j)*autoaa
    end do
    close (ich)

  case ('func')
    call setup_classify(mol,molc)
    call functional_group_classify(molc)
    call molc%print_funcgroups(stdout)


  end select
  deallocate (xyz)
  write (*,*)
  call creststop(status_normal)
end subroutine testtopo

!========================================================================================!

subroutine ensemble_analsym(fname,pr)
!*****************************************************************
!* read an ensemble and determine the symmetry for all structures
!*****************************************************************
  use crest_parameters
  use strucrd
  use symmetry_i,only: getsym
  implicit none
  character(len=*) :: fname
  logical :: pr
  integer :: nat,nall
  real(wp),allocatable :: xyz(:,:,:)
  real(wp),allocatable :: c0(:,:)
  real(wp),allocatable :: er(:)
  integer,allocatable  :: at(:)
  integer :: i,ich
  character(len=4) :: sfsym,sfsm
  character(len=80) :: atmp

  call rdensembleparam(fname,nat,nall)
  !--- allocate space and read in the ensemble
  allocate (at(nat),xyz(3,nat,nall),c0(3,nat),er(nall))
  call rdensemble(fname,nat,nall,at,xyz,er)

  if (pr) then
    write (*,*)
    call smallhead('STRUCTURE SYMMETRIES')
    write (*,'(1x,a)') 'Unlisted structures have symmetry C1'
    write (*,'(1x,a)') 'The full list can be found in the file "symmetries"'
    write (*,*)
    write (*,'(12x,a10,2x,a18,2x,a)') 'number','energy/Eh','sym.'
  end if

  xyz = xyz/bohr
  open (file='symmetries',newunit=ich)
  do i = 1,nall
    c0(1:3,1:nat) = xyz(:,:,i)
    call getsym(.false.,6,nat,at,c0,sfsym)
    sfsm = sfsym(1:3)
    write (atmp,'(3x,a,i10,2x,f18.8,2x,a)') 'structure',i,er(i),sfsm
    write (ich,'(a)') trim(atmp)
    if (pr) then
      if (trim(sfsm) /= "c1") then
        write (*,'(a)') trim(atmp)
      end if
    end if
  end do
  close (ich)

  return
end subroutine ensemble_analsym

!=========================================================================================!

function quick_rmsd(fname,nat,at,xyz,heavy) result(rout)
  use crest_parameters
  use ls_rmsd
  use strucrd
  implicit none
  character(len=*) :: fname
  integer :: nat
  integer :: at(nat)
  real(wp) :: xyz(3,nat)
  logical :: heavy
  real(wp) :: rout
  type(coord) :: mol
  real(wp),allocatable :: c0(:,:)
  real(wp),allocatable :: c1(:,:)
  integer :: i,j,k,l
  real(wp),allocatable :: gdum(:,:),Udum(:,:),xdum(:),ydum(:)
  rout = 0.0_wp

  call mol%open(fname)

  if (mol%nat .ne. nat) then
    write (stderr,*) 'dimension mismatch in quick_rmsd()'
    return
  end if
  do i = 1,nat
    if (at(i) .ne. mol%at(i)) then
      write (stderr,*) 'atom order mismatch in quick_rmsd()'
      return
    end if
  end do
  if (heavy) then !count heavy atoms
    j = 0
    do i = 1,nat
      if (at(i) .ne. 1) then
        j = j+1
      end if
    end do
    k = j
  else
    k = nat
  end if
  allocate (c0(3,k),c1(3,k),source=0.0_wp)
  l = 0
  do i = 1,nat
    if (heavy.and.at(i) == 1) cycle
    l = l+1
    c0(1:3,l) = mol%xyz(1:3,i)*bohr !mol coordinates are in bohr, need conversion
    c1(1:3,l) = xyz(1:3,i) !xyz coordinates are in Ang
  end do
  allocate (gdum(3,3),Udum(3,3),xdum(3),ydum(3))
  call rmsd(k,c0,c1,0,Udum,xdum,ydum,rout,.false.,gdum)
  deallocate (ydum,xdum,Udum,gdum)
  deallocate (c1,c0)
  return
end function quick_rmsd

subroutine quick_rmsd_tool(fname1,fname2,heavy)
  use crest_parameters
  use strucrd
  use irmsd_module
  implicit none
  character(len=*),intent(in) :: fname1
  character(len=*),intent(in) :: fname2
  logical,intent(in) :: heavy
  type(coord) :: mol,ref
  real(wp) :: rmsdval
  integer :: i
  logical,allocatable :: mask(:)

  call ref%open(fname1)
  call mol%open(fname2)

  if (heavy) then
    allocate (mask(ref%nat),source=.false.)
    do i = 1,ref%nat
      if (ref%at(i) > 1) then
        mask(i) = .true.
      end if
    end do
    rmsdval = rmsd(ref,mol,mask=mask)
  else
    rmsdval = rmsd(ref,mol)
  end if

  rmsdval = rmsdval*autoaa
  if (heavy) then
    write (*,'(1x,a,f16.8)') 'Calculated heavy atom RMSD (Å):',rmsdval
  else
    write (*,'(1x,a,f16.8)') 'Calculated RMSD (Å):',rmsdval
  end if

  return
end subroutine quick_rmsd_tool

function quick_rmsd2(nat,at,xyz,xyz2,heavy) result(rout)
  use crest_parameters
  use ls_rmsd
  implicit none
  integer :: nat
  integer :: at(nat)
  real(wp) :: xyz(3,nat)  !in Angstroem
  real(wp) :: xyz2(3,nat) !in Amgstroem
  logical :: heavy
  real(wp) :: rout
  real(wp),allocatable :: c0(:,:)
  real(wp),allocatable :: c1(:,:)
  integer :: i,j,k,l
  real(wp),allocatable :: gdum(:,:),Udum(:,:),xdum(:),ydum(:)
  rout = 0.0_wp
  if (heavy) then !count heavy atoms
    j = 0
    do i = 1,nat
      if (at(i) .ne. 1) then
        j = j+1
      end if
    end do
    k = j
  else
    k = nat
  end if
  allocate (c0(3,k),c1(3,k),source=0.0_wp)
  l = 0
  do i = 1,nat
    if (heavy.and.at(i) == 1) cycle
    l = l+1
    c0(1:3,l) = xyz(1:3,i)
    c1(1:3,l) = xyz2(1:3,i)
  end do
  allocate (gdum(3,3),Udum(3,3),xdum(3),ydum(3))
  call rmsd(k,c0,c1,0,Udum,xdum,ydum,rout,.false.,gdum)
  deallocate (ydum,xdum,Udum,gdum)
  deallocate (c1,c0)
  return
end function quick_rmsd2

!=========================================================================================!

subroutine quick_hungarian_match(fname1,fname2,heavy)
  use crest_parameters
  use strucrd
  use hungarian_module
  use axis_module,only:axis
  implicit none
  character(len=*),intent(in) :: fname1
  character(len=*),intent(in) :: fname2
  logical,intent(in) :: heavy
  type(coord) :: mol,ref
  real(wp) :: rmsdval
  integer :: i,ii,jj,nat,io,ich
  logical,allocatable :: mask(:)
  real(wp),allocatable :: C(:,:)
  real(wp),allocatable :: answers(:)
  integer,allocatable :: mapping(:)
  integer,allocatable :: hmap(:),rhmap(:)
  integer,allocatable :: a(:),b(:)
  real(wp) :: dists(3)

  call ref%open(fname1)
  call mol%open(fname2)

  !> align to rotational axes and shift to center of mass
  call axis(ref%nat,ref%at,ref%xyz)
  call axis(mol%nat,mol%at,mol%xyz)

  if (heavy) then
    allocate (mask(ref%nat),source=.false.)
    allocate (hmap(ref%nat),rhmap(ref%nat),source=0)
    nat = count((ref%at(:) > 1))
    ii = 0
    do i = 1,ref%nat
      if (ref%at(i) > 1) then
        mask(i) = .true.
        ii = ii+1
        hmap(i) = ii
        rhmap(ii) = i
      end if
    end do
  else
    allocate (mask(ref%nat),source=.true.)
    nat = ref%nat
  end if

  allocate (C(nat,nat),answers(nat))
  allocate (mapping(nat+1))
  do ii = 1,nat
    if (.not.mask(ii)) cycle
    do jj = 1,nat
      if (.not.mask(jj)) cycle
      dists(:) = (ref%xyz(:,ii)-mol%xyz(:,jj))**2
      if (heavy) then
        C(hmap(jj),hmap(ii)) = sqrt(sum(dists))
      else
        C(jj,ii) = sqrt(sum(dists))
      end if
    end do
  end do
  allocate (a(nat),b(nat))
  call lsap(C,nat,nat,a,b)

  write (*,'(a,3(1x,a))') 'Assignment:',fname2,'-->',fname1
  do i = 1,nat
    if (heavy) then
      write (*,'(i6," --> ",i6)') rhmap(a(i)),rhmap(b(i))
    else
      write (*,'(i6," --> ",i6)') a(i),b(i)
    end if
  end do
  write (*,*)
  !> write the rotated and shifted coordinates to one file
  open (newunit=ich,file='lsap.xyz')
  call ref%append(ich)
  call mol%append(ich)
  close (ich)

  !> reconstruct RMSD from assignment (since our costs are already distances!)
  rmsdval = 0.0_wp
  do i = 1,nat
    rmsdval = rmsdval+C(a(i),b(i))/real(nat,wp)
  end do
  rmsdval = sqrt(abs(rmsdval))*autoaa
  if (heavy) then
    write (*,'(1x,a,f16.8)') 'Calculated heavy atom RMSD (Å):',rmsdval
  else
    write (*,'(1x,a,f16.8)') 'Calculated RMSD (Å):',rmsdval
  end if

  return
end subroutine quick_hungarian_match

!=========================================================================================!

subroutine resort_ensemble(fname)
!************************************************
!* resort all structures of a given ensemblefile
!************************************************
  use crest_parameters
  use strucrd
  use crest_data
  implicit none
  character(len=*) :: fname
  integer :: nat,nall
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable :: at(:)
  character(len=128),allocatable :: comm(:)
  logical :: ex
  integer :: i,j,k,ich

  integer,allocatable :: atorder(:)

  inquire (file='.atorder',exist=ex)
  if (.not.ex) error stop 'file .atorder does not exist'

  call rdensembleparam(fname,nat,nall)
  !--- allocate space and read in the ensemble
  allocate (at(nat),xyz(3,nat,nall),comm(nall))
  call rdensemble(fname,nat,nall,at,xyz,comm)

  allocate (atorder(nat))
  open (newunit=ich,file='.atorder')
  do i = 1,nat
    read (ich,*) j,k
    atorder(j) = k
  end do
  close (ich)

  open (newunit=ich,file='ensemble_tmp.xyz')
  do i = 1,nall
    write (ich,'(2x,i0)') nat
    write (ich,'(a)') trim(comm(i))
    do j = 1,nat
      k = atorder(j)
      write (ich,'(1x,a2,1x,3f20.10)') i2e(at(k),'nc'),xyz(1:3,k,i)
    end do
  end do
  close (ich)

  return
end subroutine resort_ensemble

!=========================================================================================!

subroutine etoerel(n,er,erel,factor)
!*****************************************************************
!* convert an array of absolute energies (in Eh) to relative ones
!*****************************************************************
  use crest_parameters
  implicit none
  integer,intent(in)  :: n
  real(wp),intent(in) :: er(n)
  real(wp),intent(in) :: factor
  real(wp),intent(out) :: erel(n)
  erel(:) = (er(:)-minval(er,1))*factor
end subroutine etoerel

!=========================================================================================!

subroutine backchain(nat,A,at,maxl,path)
!*****************************************************
!* identify the "backchain" of a molecular structure
!* The backchain here will be defined as the longest
!* path between two heavy (i.e. non-Hydrogen) atoms
!* that is walkable in the molecular graph
!*****************************************************
  use crest_parameters
  use adjacency
  implicit none
  integer,intent(in) :: nat
  integer,intent(in) :: A(nat,nat)
  integer,intent(in) :: at(nat)
  integer,intent(out) :: maxl
  integer,intent(out) :: path(nat)
  !> LOCAL
  real(wp),allocatable :: dist(:,:)
  integer,allocatable :: prev(:,:)

  integer :: i,j,k,l
  integer :: maxi,maxj

  !call wbo2adjacency(nat,wbo,A,0.02_wp)

  allocate (dist(nat,nat),source=0.0_wp)
  allocate (prev(nat,nat),source=0)
  call FloydWarshall_simple(nat,A,dist,prev)

  !allocate(path(nat), source=0)
  maxl = 0
  maxi = 0
  maxj = 0
  do i = 1,nat
    if (at(i) == 1) cycle !> no H
    do j = 1,i-1
      if (at(j) == 1) cycle !> no H
      path = 0
      call getPathFW(nat,prev,i,j,path,l)
      if (l > maxl) then
        maxi = i
        maxj = j
        maxl = l
      end if
    end do
  end do

  write (*,'(a,i0,a,i0,a,i0)') 'Longest chain identified between atoms ',maxi,' and ', &
  & maxj,' with a length of ',maxl
  path = 0
  call getPathFW(nat,prev,maxi,maxj,path,maxl)

  do i = maxl,1,-1
    write (*,*) path(i)
  end do

end subroutine backchain

!=========================================================================================!
!=========================================================================================!
