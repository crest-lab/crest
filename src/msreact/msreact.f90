!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020 Stefan Grimme, Philipp Pracht
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

!==============================================================!
! the handler function is the sub-program that is called
! by crest, i.e., within this routine the seperate input file is
! read and the calculation is started.
!==============================================================!
subroutine msreact_handler(env,tim)
  use crest_parameters
  use crest_data
  use msmod
  use strucrd
  use atmasses
  use iomod
  use zdata
  implicit none

  type(systemdata) :: env
  type(timer) :: tim

  type(msobj) :: mso  !a new msreact object
  type(coord) :: struc

  interface
    subroutine msreact_topowrap(mol,pair,paths,wboname)
      import :: msmol
      implicit none
      type(msmol) :: mol
      integer :: pair(mol%nat*(mol%nat+1)/2)
      integer :: paths(mol%nat*(mol%nat+1)/2,mol%nat)
      character(len=*),optional :: wboname
    end subroutine msreact_topowrap
  end interface
  integer :: nat
  integer,allocatable :: pair(:)
  integer,allocatable :: paths(:,:)
  integer :: k,T,Tn
  integer :: nisomers,nfragpairs
  real(wp) :: molmass
  real(wp) :: estart ! energy of input structure

  call tim%start(1,'MSREACT')

  call msreact_head(env,mso)

  ! get starting energy of lowest fragment
  call get_input_energy(env,mso%T,estart)
  !-- read the input coord and put it into the
  !   product-list as first structure
  call struc%open('coord')
  struc%xyz = struc%xyz*bohr !to Angstrom, from this point on by convention!
  molmass = molweight(struc%nat,struc%at)
  call mso%pl%append(struc%nat,struc%at,struc%xyz,estart,env%chrg,1,molmass)
  if (env%mslargeprint) call wrshort_real("mass",mso%pl%mol(1)%molmass) ! write mass
  call struc%deallocate()

  !-- additional input file can be read here
  call msinputreader(mso,env%msinput)

  nat = mso%pl%mol(1)%nat
  k = nat*(nat+1)/2
  allocate (pair(k),paths(k,nat))
  call msreact_topowrap(mso%pl%mol(1),pair,paths,'wbo')

  !-- setting the threads for correct parallelization
  call new_ompautoset(env,'auto',0,T,Tn)

  !-- do the directory setup and optimizations
  call msreact(env,mso,nat,pair)
  deallocate (paths,pair)

  if (env%mslargeprint) then
    write (*,*) "All structures written to file <crest_allproducts.xyz>"
    call wrplist(mso,'crest_allproducts.xyz')
  end if
  ! detect fragmented structures
  call detect_fragments(mso,env%mslargeprint,nisomers,nfragpairs)
  call increase_fragdist(mso,env%mslargeprint)
  write (*,*) "Remaining structures:",mso%pl%nmol-1

  ! call write_products(env,mso) print out structures if lprint is chosen?

  ! sort out topologically identical structures
  call sortoutduplicates(env,mso)

  write (*,'(a,f10.4,a)') "Calling CREGEN routine and sort out structures above energy window of",env%ewin," kcal/mol"
  ! write plist to file to use cregen routine
  call wrplist(mso,'crest_unique_products.xyz')
  if (env%mslargeprint) then
    write (*,*) "All unique structures written to file <crest_unique_products.xyz>"
  end if
  ! sort according to energies and RMSD, remove highe energy structures
  call newcregen(env,13,'unique_products.xyz')
  call mso%pl%dealloc()
  ! read in again plist after cregen and remove starting structure
  call rdplist(mso,estart,'crest_msreact_products.xyz')
  write (*,*) "Remaining structures after cregen:",mso%pl%nmol
  call detect_fragments(mso,env%mslargeprint,nisomers,nfragpairs)

  ! write final structures to file and directories
  call write_fragments(env,mso,estart,nisomers,nfragpairs,'crest_msreact_products.xyz',env%mslargeprint)

  !cleanup files
  if (.not.env%mslargeprint) then
    call rmrf('MSDIR')
    call rmrf('crest_allproducts.xyz')
    call rmrf('crest_unique_products.xyz')
  end if
  call rmrf('crest_unique_products.sorted')
  call rmrf('scoord.1')
  call rmrf('xtbscreen.xyz xtblmoinfo')
  call tim%stop(1)
  return
end subroutine msreact_handler

!==============================================================!
! the main implementation of the msreact algo
!==============================================================!
subroutine msreact(env,mso,nat,pair)
  use crest_parameters
  use msmod
  use iomod
  use miscdata,only:rcov
  use utilities,only:lin
  implicit none

  type(msobj) :: mso    !main storage object
  type(msmol) :: inmol  ! xyz of molecule
  type(msmol) :: inmoldisp  ! xyz of molecule with displacement
  type(systemdata) :: env ! system data

  integer :: nat
  integer :: pair(nat*(nat+1)/2) ! gives number of bonds between two atoms lin(i,j)
  integer :: nbonds
  integer :: i,j,k,l
  integer :: p
  integer :: np
  integer :: io,T,Tn
  integer ::  nbaseat ! number of lewis basic atoms according to a xTB LMO analysis
  integer,allocatable :: basicatlist(:) ! list of lewis basic atoms according to a xTB LMO analysis

  character(len=:),allocatable :: subdir
  character(len=40) :: pdir
  character(len=30) :: base
  character(len=512) :: thisdir
  real(wp) :: constr_dist

  inmol = mso%pl%mol(1)

  !-- main subdirectory handling
  call getcwd(thisdir)
  subdir = 'MSDIR'
  io = makedir(subdir)
  call chdir(subdir)

  !-- get specific pairs
  base = 'Pair_'
  np = 0
  do i = 1,nat
    do j = i,nat
      k = lin(i,j)
      p = pair(k)
      if (p .eq. 1.and.(inmol%at(i) .eq. 1.or.inmol%at(j) .eq. 1)) cycle ! do not distort all X-H bonds
      if (p .gt. 0.and.p .le. env%msnbonds) then
        !write(*,*) " distort bond ",i,"-",j," with ",p," bonds"
        np = np+1
        write (pdir,'(a,i0)') trim(base),np
        constr_dist = mso%cdist*(rcov(inmol%at(i))+rcov(inmol%at(j)))*bohr+float(p)
        ! write(*,*) inmol%at(i),inmol%at(j),constr_dist
        call isodir(mso,trim(pdir),inmol,i,j,constr_dist,mso%fc)
      end if
    end do
  end do

  !attractive potential for H-shifts
  !for H-shifts: xtb -lmo , read xtb.out -> get list of atoms with pi or lone-pair -> H-allowed to shift there
  if (env%msattrh) then
    call chdir(thisdir)
    call xtblmo(env,.true.)
    call readbasicpos(env,nbaseat,basicatlist)
    call chdir(subdir)
    write (*,*) "Add attractive potential for H-shifts:"
    !-- setting the threads for correct parallelization (were changed in xtblmo call)
    call new_ompautoset(env,'auto',0,T,Tn)
    do i = 1,nat
      do j = i,nat
        k = lin(i,j)
        p = pair(k)
        if (p .ge. 2) then ! only add attractive potential for pairs with at least 2 bonds distance
          if ((inmol%at(i) .eq. 1.and.(findloc(basicatlist,j,1) /= 0)).or. & ! pair of H and basic atom?
          &  (inmol%at(j) .eq. 1.and.(findloc(basicatlist,i,1) /= 0))) then
            if (msmoldist(inmol,i,j) .lt. mso%distthr_att) then
              np = np+1
              !write(*,*) "add attractive potential for pair ", np," between ",i," and ",j
              write (pdir,'(a,i0)') 'Pair_',np
              constr_dist = mso%cdist_att*(rcov(inmol%at(i))+rcov(inmol%at(j)))*bohr+float(p)
              call isodir(mso,trim(pdir),inmol,i,j,constr_dist,mso%fc_attr)
            end if
          end if
        end if
      end do
    end do
  end if

  write (*,*) '# of distortions',np
  call msreact_jobber(env,np,base,.true.,.false.)
  call msreact_collect(nat,np,base,mso)

  !-- do additional shifting of atoms and subsequent optimization
  if (env%msnshifts .gt. 0) then
    write (*,'(a,i0,a)') "Shift atoms ",env%msnshifts," times and reoptimize"
    base = 'Distortion_'
    np = 0
    do l = 1,env%msnshifts
      call shiftatoms(mso,inmol,inmoldisp)
      np = np+1
      write (pdir,'(a,i0)') trim(base),np
      call isodiropt(mso,trim(pdir),inmoldisp)
    end do
    call msreact_jobber(env,np,base,.false.,.false.)
    ! collect everything in isomer list
    call msreact_collect(nat,np,base,mso)
  end if

  ! do additional shifting of atoms and subsequent optimization with constrained bonds
  if (env%msnshifts2 .gt. 2) then
    write (*,'(a,i0,a)') "Shift atoms",env%msnshifts2," times and apply repulsive potential on bonds of distorted structure"
    base = 'Distortedpair_'
    np = 0
    do l = 1,env%msnshifts2
      call shiftatoms(mso,inmol,inmoldisp)
      do i = 1,nat
        do j = i,nat
          k = lin(i,j)
          p = pair(k)
          if (p .eq. 1.and.(inmol%at(i) .eq. 1.or.inmol%at(j) .eq. 1)) cycle
          if (pair(k) .le. env%msnbonds) then
            np = np+1
            write (pdir,'(a,i0)') trim(base),np
            constr_dist = mso%cdist*(rcov(inmoldisp%at(i))+rcov(inmoldisp%at(j)))*bohr+float(p)
            call isodir(mso,trim(pdir),inmoldisp,i,j,constr_dist,mso%fc)
          end if
        end do
      end do
    end do
    call msreact_jobber(env,np,base,.true.,.false.)
    ! collect everything in isomer list
    call msreact_collect(nat,np,base,mso)
  end if

  call chdir(thisdir)
  return
end subroutine msreact

!============================================================!
! make a dir for a constrained optimization of a structure,
! a controlfile with constraints on atoms A and B (at dist D)
! will be written into the directory
!============================================================!
subroutine isodir(mso,dirname,mol,A,B,D,fc)
  use crest_parameters
  use msmod
  use iomod
  use strucrd,only:wrxyz
  implicit none
  type(msobj) :: mso
  character(len=*) :: dirname
  type(msmol),intent(in) :: mol
  integer :: A,B
  real(wp) :: D
  real(wp),intent(in) :: fc ! force constant for bias potential

  character(len=:),allocatable :: fname
  character(len=20) :: dumm
  integer :: io,ich

  io = makedir(dirname) !create the directory

  fname = trim(dirname)//'/'//'struc.xyz'
  open (newunit=ich,file=fname)
  call wrxyz(ich,mol%nat,mol%at,mol%xyz)
  close (ich)

  fname = trim(dirname)//'/'//'.CHRG'
  open (newunit=ich,file=fname)
  write (ich,'(i0)') mol%chrg    ! EI +1, DEA -1, CID 0
  close (ich)

  ! write control file for bond length constraint

  fname = trim(dirname)//'/'//'.xc1'
  open (newunit=ich,file=fname)
  write (ich,'(a)') '$scc'
  write (dumm,'(f16.2)') mso%T
  write (ich,'(1x,a,a)') 'temp=',adjustl(trim(dumm))
  write (ich,'(a)') '$constrain'
  write (dumm,'(f16.4)') fc
  write (ich,'(3x,a,a)') 'force constant=',adjustl(trim(dumm))
  write (ich,'(3x,a,1x,i0,a,1x,i0,a,1x,f8.5)') 'distance:',A,',',B,',',D
  close (ich)

  fname = trim(dirname)//'/'//'.xc2'
  open (newunit=ich,file=fname)
  write (ich,'(a)') '$scc'
  write (dumm,'(f16.2)') mso%T
  write (ich,'(1x,a,a)') 'temp=',adjustl(trim(dumm))
  write (ich,'(a)') '$opt'
  write (ich,'(1x,a)') 'maxcycle=',mso%maxc ! only allow 15 cycles to avoid recombination of fragments
  write (ich,'(a)') '$write'
  write (ich,'(1x,a)') 'wiberg=true'
  close (ich)

  return
end subroutine isodir

!============================================================!
! make a dir for a optimization of a displaced structure,
! a controlfile for an optimization after the atom shifting
! will be written into the directory
!============================================================!
subroutine isodiropt(mso,dirname,mol)
  use iso_fortran_env,only:wp => real64
  use msmod
  use iomod
  use strucrd,only:wrxyz
  implicit none
  type(msobj) :: mso
  character(len=*) :: dirname
  type(msmol),intent(in) :: mol
  integer :: A,B
  real(wp) :: D

  character(len=:),allocatable :: fname
  character(len=20) :: dumm
  integer :: io,ich

  io = makedir(dirname) !create the directory

  fname = trim(dirname)//'/'//'struc.xyz'
  open (newunit=ich,file=fname)
  call wrxyz(ich,mol%nat,mol%at,mol%xyz)
  close (ich)

  fname = trim(dirname)//'/'//'.CHRG'
  open (newunit=ich,file=fname)
  write (ich,'(i0)') mol%chrg   ! EI +1, DEA -1, CID 0, give in crest call
  close (ich)

  fname = trim(dirname)//'/'//'.xc1'
  open (newunit=ich,file=fname)
  write (ich,'(a)') '$scc'
  write (dumm,'(f16.2)') mso%T
  write (ich,'(1x,a,a)') 'temp=',adjustl(trim(dumm))
  close (ich)
  return
end subroutine isodiropt

!=====================================================================!
! The job construction routine for MSREACT
! (will have to be modified later, for now it is for testing)
!=====================================================================!
subroutine msreact_jobber(env,ndirs,base,constr,niceprint)
  use crest_parameters
  use msmod
  use iomod
  implicit none
  integer :: ndirs
  character(len=*) :: base
  logical :: niceprint,constr
  character(len=1024) :: jobcall,jobcall2
  type(systemdata) :: env ! system data

  jobcall = ''
  jobcall2 = ''

  if (constr) then
    write (jobcall,'(a)') 'xtb struc.xyz  --opt loose --input .xc1 '//trim(env%gfnver)//' > split.out 2>/dev/null'
    write (jobcall2,'(a)') 'xtb xtbopt.xyz --opt crude --input .xc2 '//trim(env%gfnver)//' > xtb.out 2>/dev/null'
    jobcall = trim(jobcall)//' ; '//trim(jobcall2)
  else
    write (jobcall,'(a)') 'xtb struc.xyz --opt crude --input .xc1 '//trim(env%gfnver)//' > xtb.out 2>/dev/null'
  end if

  !-- directories must be numbered consecutively
  call opt_OMP_loop(ndirs,base,jobcall,niceprint)
  write (*,*)
  write (*,*) 'done.'
  return
end subroutine msreact_jobber

!=====================================================================!
! A wrapper to generate the topology for a molecule within the
! MSREACT subprogram
!=====================================================================!
subroutine msreact_topowrap(mol,pair,paths,wboname)
  use crest_parameters
  use msmod
  use zdata
  use adjacency
  use utilities,only:lin
  implicit none
  type(msmol) :: mol
  integer :: pair(mol%nat*(mol%nat+1)/2)
  !integer :: pair(mol%nat,mol%nat)
  integer :: paths(mol%nat*(mol%nat+1)/2,mol%nat)
  character(len=*),optional :: wboname
  type(zmolecule) :: zmol

  integer,allocatable :: A(:,:)
  integer,allocatable :: prev(:,:)
  real(wp),allocatable :: E(:,:)
  real(wp),allocatable :: dist(:,:)

  integer :: lpath,i,j,k
  integer,allocatable :: path(:)
  logical :: ex

  ex = .false.
  if (present(wboname)) then
    inquire (file=wboname,exist=ex)
  end if
  if (ex) then
    call simpletopo(mol%nat,mol%at,mol%xyz,zmol,.false.,.false.,wboname)
  else
    mol%xyz = mol%xyz/bohr !CN based topo requires Bohrs
    call simpletopo(mol%nat,mol%at,mol%xyz,zmol,.false.,.false.,'')
    mol%xyz = mol%xyz*bohr
  end if

  allocate (A(mol%nat,mol%nat),E(mol%nat,mol%nat))
  call zmol%adjacency(A,E)

  allocate (prev(mol%nat,mol%nat),dist(mol%nat,mol%nat))

  call FloydWarshall(mol%nat,A,E,dist,prev)
  allocate (path(mol%nat),source=0)
  do i = 1,mol%nat
    do j = i,mol%nat
      path = 0
      call getPathFW(mol%nat,prev,i,j,path,lpath)
      !write(*,*) path(1:lpath)
      k = lin(i,j)
      pair(k) = lpath-1 ! number of bonds
      paths(k,:) = path(:)
    end do
  end do

  deallocate (dist,prev)
  deallocate (E,A)

  call zmol%deallocate() !clear the zmol memory
  return
end subroutine msreact_topowrap

!========================================================================!
! collect structures of optimized molecules
! xyz files should still have the same number and order of atoms
!========================================================================!
subroutine msreact_collect(nat,np,base,mso)
  use crest_parameters
  use strucrd
  use msmod
  use atmasses
  implicit none
  integer :: nat
  integer :: np
  integer :: ich
  integer :: chrg ! charge of the molecule
  type(msobj) :: mso
  character(len=40) :: pdir
  character(len=*) :: base
  character(len=:),allocatable :: optfile
  character(len=128) :: newcomment
  integer :: p
  logical :: ex
  integer,allocatable :: at(:)
  real(wp),allocatable :: xyz(:,:)
  real(wp) :: etot
  real(wp) :: molmass

  allocate (at(nat),xyz(3,nat))

  do p = 1,np
    write (pdir,'(a,i0)') trim(base),p
    optfile = trim(pdir)//'/'//'xtbopt.xyz'
    inquire (file=optfile,exist=ex)
    if (ex) then
      call rdcoord(optfile,nat,at,xyz,etot)
      xyz = xyz*bohr
      call rdshort(trim(pdir)//'/'//'.CHRG',chrg)
      molmass = molweight(nat,at)
      call mso%pl%append(nat,at,xyz,etot,chrg,0,molmass)
    end if
  end do
  close (ich)

  deallocate (xyz,at)
  return
end subroutine msreact_collect

!============================================================!
! shift atoms of input structure randomly
! to generate more structures
! for planar molecules important
!============================================================!
subroutine shiftatoms(mso,molin,molout)
  use iso_fortran_env,only:wp => real64
  use msmod
  use iomod
  use strucrd,only:wrxyz
  implicit none
  type(msobj) :: mso
  type(msmol) :: molin
  type(msmol) :: molout
  real(wp) :: D
  real(wp) :: x,shift
  real(wp) :: dist ! distance to shift atoms
  character(len=80) :: fname
  character(len=20) :: dumm
  integer :: io,ich,i,j,count,incr

  dist = mso%atomshift ! distance of displacement
  molout = molin
  do i = 1,molin%nat
    do j = 1,3
      call Random_Number(x)
      ! positive and negative direction possible
      x = 2.0_wp*x-1.0_wp
      shift = x*dist
      molout%xyz(j,i) = molin%xyz(j,i)+shift
    end do
  end do
  return
end subroutine shiftatoms

