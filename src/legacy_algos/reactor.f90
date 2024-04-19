!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020 Stefan Grimme
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
!===================================================================!
! This file contains routines related to the xtb nano-reactor, 
! first published in S. Grimme, JCTC, 2019, 15, 2847-2862.
!===================================================================!

!======================================================!
! main routine for the nanoreactor
!======================================================!
subroutine reactor(env, tim)
  use crest_parameters, only: wp, bohr
  use crest_data
  use iomod
  use zdata
  use strucrd, only: rdensembleparam,rdensemble
  implicit none

  interface
    subroutine analyse_fragments(molvec, molvecSize, totalFragCount, atomsPerFragCount)
      implicit none
      integer, intent(in):: molvecSize
      integer, dimension(molvecSize), intent(in):: molvec
      integer, intent(out):: totalFragCount
      integer, dimension(:), allocatable, intent(out):: atomsPerFragCount
    end subroutine analyse_fragments
  end interface

  !-- Forward declarations
  logical:: are_zequals_different

  type(systemdata):: env    ! MAIN STORAGE OS SYSTEM DATA
  type(timer):: tim

  !--- declarations by Tim

  ! DECLARATIONS
  ! for purpose of variables look at the lines where they are used
  character(len=100):: trjFile = "xtb.trj"
  integer:: nat              ! <-- nat
  integer:: snapCount              ! <-- nstr

  ! ensemble variables
  real(wp), dimension(:,:,:), allocatable:: file_xyz
  real(wp), dimension(:), allocatable    :: energy
  integer, dimension(:), allocatable    :: iat   ! <-- Ordnungszahlen der Atome
  integer:: allocStat

  ! mrec variables
  integer:: molcount, i
  integer, dimension(:), allocatable  :: molvec  ! <-- assing atom number to fragment

  ! zgrp and zequal variables
  integer:: j, k, l
  integer, dimension(:,:), allocatable:: molvecs  ! dimesion(nat, nall)
  type(zequal), dimension(:), allocatable:: zequals
  integer:: fragCount
  integer, dimension(:), allocatable:: atomsPerFragCount
  type(zgrp), dimension(:), allocatable:: zgrps

  ! variables for comparison of all zequals
  logical, dimension(:), allocatable:: taken
  type(zequal):: zequal_ref


  real(wp),allocatable :: bonds(:,:)
  integer,allocatable :: nghref(:),ngh(:)
  integer :: nghdim

  integer :: ndirs

  !--- start timer
  call tim%start(1, 'nano reactor')
  call prreactorhd()



  !--- 1. set up reactor metadynamics
  !    wall-potential, mass density, simulation length


  !--- 2. run the metadynamics


  !--- 3. read the trajectory and analyze it 
  !    get fragments, fragment charges, etc.
  write(*,*)"Trajectory file: ", trim(trjFile)

  ! get counts of atoms and structures
  call rdensembleparam(trjFile, nat, snapCount)
  write(*,*)
  write(*,*)"Number of atoms      = ", nat
  write(*,*)"Number of snapshots  = ", snapCount

  ! prepare array for next step
  allocate(file_xyz(3, nat, snapCount), stat = allocStat)
  allocate(energy(snapCount), stat = allocStat)
  allocate(iat(nat), stat = allocStat)

  ! fill those arrays with their respective values
  call rdensemble(trjFile, nat, snapCount, iat, file_xyz, energy)
  !print *, ""
  !print *, "File xyz      = ", file_xyz(:, 1, 1:3)
  !print *, "Atomic number = ", iat
  !print *, "Energy        = ", energy(1:3)

  ! prepare further arrays
  ! molvecs will hold all molvec for all snapshots
  allocate(molvecs(nat, snapCount), stat = allocStat)
  allocate(molvec(nat), stat = allocStat)
  allocate(zequals(snapCount), stat = allocStat)

  ! get molvec
  file_xyz = file_xyz/bohr
  do i = 1, snapCount
    call mrec(molcount, file_xyz(1:3, 1:nat, i), nat, iat, molvec)
    molvecs(:, i) = molvec(:)

    ! counts total number of fragments in a snapshot and how many atoms each fragment has (latter outputed as vector)
    call analyse_fragments(molvec, nat, fragCount, atomsPerFragCount) 

    ! zgrps is array; holds all zgrps for one snapshot
    ! write atomnumbers in their respective zgrp
    allocate(zgrps(fragCount))
    do j = 1, fragCount
      call zgrps(j)%allocate(atomsPerFragCount(j))
      l = 1
      do k = 1, nat
        if (molvec(k) == j) then
          zgrps(j)%mem(l) = k
          l = l+1
        end if
      end do
    end do

    ! sort zgrps, so that the zgrps are sorted by their first atoms position number in ascending order
    call sort_zgrps(zgrps, fragCount)

    ! zqeuals is array with so many zqeual objects as there are snapshots
    ! filling all parameters of the zequal object
    call zequals(i)%allocate(nat)
    zequals(i)%ng = fragCount
    zequals(i)%grp(:fragCount) = zgrps
    zequals(i)%ord = molvec
    call zequals(i)%geteng

    ! some printing to check, that all variables are doing fine
    !print *, ""
    !print *, "i = ", i
    !print *, "Molcount            = ", molcount
    !! print *, "Molvec              = ", molvec
    !print *, "Molvecs(i)          = ", molvecs(:, i)
    !print *, "Number of fragments = ", fragCount
    !print *, "Atoms per fragment  = ", atomsPerFragCount

    !print *, ""
    !print *, "<ZQEUAL VARS>"
    !print *, "Number of groups        = ", zequals(i)%ng
    !print *, "Grp more than 1 nuclei  = ", zequals(i)%eng
    !print *, "Group order = ", zequals(i)%ord
    !print *, "Groups"
    !do j = 1, fragCount
    !  print *, "j      = ", j
    !  print *, "Groups = ", zgrps(j)%mem
    !  print *, ""
    !end do

    deallocate(atomsPerFragCount)
    deallocate(zgrps)
  end do


  ! <Comparing all Zeuqal objects>
  allocate(taken(snapCount), source=.FALSE.)
  taken(1) = .TRUE.
  call update_zequal_ref(zequals(1), zequal_ref)
  nghdim=nat*(nat+1)/2
  allocate(bonds(nat,nat))
  allocate(ngh(nghdim),nghref(nghdim))
  call reactorneighbours(nat,iat,file_xyz(1:3,1:nat,1),bonds,nghdim,nghref)

  do j = 2, snapCount
    taken(j) = are_zequals_different(zequal_ref, zequals(j))
    if (taken(j)) then
      call update_zequal_ref(zequals(j), zequal_ref)
    else
      call reactorneighbours(nat,iat,file_xyz(1:3,1:nat,j),bonds,nghdim,ngh)
      call arrcomp(nghdim,nghref,nghdim,ngh,taken(j))
      taken(j) = .not.taken(j)
      if(taken(j)) nghref = ngh
    end if
  End do

  print *, ""
  print *, "<Comparing structures>"
  print *, "Taken = ", taken

  write(*,'(1x,i0,a,i0,a)')count(taken(:)),' of ',snapCount,' taken.'

  ! </Comparing all Zeuqal objects>
  !--- 4. optimize the fragments
  if(env%restartopt)then
    call reactorreopt(env,nat,iat,snapCount,file_xyz,taken,zequals,ndirs)
    call collectproducts('OPTIM','TMPFRG',ndirs,'crest_products.xyz',env%riso)
  endif

  deallocate(nghref,ngh,bonds)
  deallocate(taken)
  deallocate(zequals)
  deallocate(molvec)
  deallocate(iat)
  deallocate(energy)
  deallocate(file_xyz)

  !--- stop timer
  call tim%stop(1)   

  return
end subroutine reactor

!======================================================!
! subroutine reactor_setup
!   read in coordinates from "coord" file and calculate
!   a standard ellipsoide potential.
!   The ellipsoide axes then have to be scaled to 
!   fit the required reactor density
!======================================================!
subroutine reactor_setup(env)
  use crest_parameters, only: wp
  use crest_data
  use iomod
  use filemod
  use atmasses
  use strucrd, only: rdnat,rdcoord,wrc0
  use axis_module 
  implicit none
  type(systemdata):: env    ! MAIN STORAGE OS SYSTEM DATA
  integer:: nat                   !# of atoms
  integer, allocatable  :: at(:)    !atom types
  real(wp), allocatable:: xyz(:,:)  ! coordinates
  real(wp):: rabc(3)  ! axes of the ellipsoide potential
  real(wp) :: mass
  real(wp) :: dens       !density
  real(wp) :: dum
  integer :: j
  character(len = 256):: atmp,btmp
  type(filetype):: f

  call remove("rcontrol")
  call f%open("rcontrol")

  !---------------------------
  if(env%preactormtd)then
    write(*,'(/,1x,a)') 'Metadynamics settings:'
    call f%write("$md")
    write(btmp,'(f12.2)')env%mdtime
    write(*,'(1x,a,a,a)')'Simulation time: ',trim(adjustl(btmp)),' ps'
    write(atmp, '("time=",a)') trim(adjustl(btmp))
    call f%write(trim(atmp))
    write(btmp,'(f12.1)')1.0d0
    write(atmp, '("step=",a)')trim(adjustl(btmp))
    call f%write(trim(atmp))
    write(atmp, '("shake=",i0)')0
    call f%write(trim(atmp))

    call f%write('$set')
    write(atmp,'(2x,a,2x,i0)')'mddump',2000
    call f%write(trim(atmp))

    call f%write('$metadyn')
    write(atmp, '("save=",i0)') env%metadlist(1)
    call f%write(trim(atmp))
    write(btmp,'(f12.6)')env%metadfac(1)*env%nat
    write(*,'(1x,a,a,a)') 'Vbias (k): ',trim(adjustl(btmp)),' Eh'
    write(atmp, '("kpush=",a)') adjustl(trim(btmp))
    call f%write(trim(atmp))
    write(btmp,'(f12.6)')env%metadexp(1)
    write(*,'(1x,a,a,a)') 'Vbias (α): ',trim(adjustl(btmp)),' Bohr⁻²'
    write(atmp, '("alp=",a)')adjustl(trim(btmp))
    call f%write(trim(atmp))

  endif

  !-------------------------
  if(env%preactorpot)then
    write(*,'(/,1x,a)') "Generating spherical logfermi potential:"
    !-- read coord
    call rdnat('coord',nat)
    allocate(at(nat), xyz(3, nat))
    call rdcoord('coord',nat, at, xyz)

    !--- CMA trafo
    call axistrf(nat,nat,at,xyz)

    !--- get the "box"
    do j=1,3
       rabc(j)=maxval(xyz(j,1:nat))-minval(xyz(j,1:nat))
    enddo
    dum=maxval(rabc)
    rabc(1:3)=dum

    !-- scale wall potential for given input density
    mass = molweight(nat, at)
    call reactor_pot(rabc, mass, env%rdens, dens)

    write(*,'(1x,a,f9.3,a)') 'Spherical cavity radius : ',rabc(1),' Bohr' 
    write(*,'(1x,a,f9.1,a)') 'Logfermi temperature    : ',env%tempfermi,' K'

    !-- write CMA coords
    call wrc0('coord',nat, at, xyz)

    !-- write xtb input file
    call f%write("$wall")
    call f%write("  potential=logfermi")
    write(atmp, '("  sphere:",1x,g0,",",1x,"all")') rabc(1)
    call f%write(trim(atmp))
    write(atmp, '(f16.1)')env%tempfermi
    atmp = '  temp='//trim(adjustl(atmp))
    call f%write(trim(atmp))
    deallocate(at, xyz)
  endif

  if(f%nlines .gt.1)then
    call f%write('$end')
    write(*,'(/,1x,a)') 'Base settings written to file <rcontrol>'
  endif    
  call f%clearblanks
  call f%flushclose
  
  return
end subroutine reactor_setup

!======================================================!
! subroutine reactor_pot
!   scale a standard-ellipsoide potential to the
!   the required density
!======================================================!
subroutine reactor_pot(ax, wei, densref, dens)
  use crest_parameters, only: wp, bohr
  implicit none
  real(wp), intent(inout):: ax(3)
  real(wp), intent(in)    :: wei
  real(wp), intent(in)    :: densref
  real(wp), intent(out)   :: dens

  real(wp):: scal
  real(wp) :: mass

  real(wp),parameter :: pi43  = 3.1415926540_wp*(4.0_wp/3.0_wp)
  real(wp),parameter :: pi    = 3.1415926540_wp
  real(wp),parameter :: third = 1.0_wp/3.0_wp

  dens = 0.0_wp

  !-- Volume of ellipsoide: 4/3 * π * a*b*c
  !vol = pi43*ax(1) * ax(2) * ax(3)

  !-- weight
  mass = wei *1.66053886E-27

  !-- unscaled density
  dens = mass/(ax(1)*ax(2)*ax(3))

  !-- transform into g/cm³   
  dens = dens * ( 0.001d0 / (1.0d-30*bohr**3.0d0))
  !dens=0.001*mass/(1.0d-30*ax(1)*ax(2)*ax(3)*bohr**3)

   write(*, '(1x, a, f8.2,a)')'Reactor density (unscaled):',dens,' g/cm³'

!   return

  if(densref .lt. 0.0d0)then
    !-- return unscaled potential
    ax = ax*1.0d0      
    return
  endif

  !-- get scaling factor of the axes = (ρ/ρ_ref)^(1/3)
  scal = (dens/densref) ** third

  !-- apply scaling
  ax = ax*scal
  !vol = pi43*ax(1) * ax(2) * ax(3)
  dens=0.001*mass/(1.0d-30*ax(1)*ax(2)*ax(3)*bohr**3)
    
  write(*, '(1x, a, f8.2,a)')'Reactor density (from input):',dens,' g/cm³'


  return
end subroutine reactor_pot


! ### TIM'S PROCEDURES ###

subroutine analyse_fragments(molvec, molvecSize, totalFragCount, atomsPerFragCount)
  !!! counts total number of fragments in one snapshot and counts how many atoms each fragment has
  implicit none

  integer, intent(in):: molvecSize
  integer, dimension(molvecSize), intent(in):: molvec
  integer, intent(out):: totalFragCount
  integer, dimension(:), allocatable, intent(out):: atomsPerFragCount
  integer:: i, ii

  totalFragCount = maxval(molvec)
  allocate(atomsPerFragCount(totalFragCount), source = 0)

  do i = 1, molvecSize
    ii = molvec(i)
    atomsPerFragCount(ii) = atomsPerFragCount(ii) + 1
  end do

  return
end subroutine analyse_fragments

function are_zequals_different(zequal_1, zequal_2) result(different)
  ! TRUE : zequal objects are different
  ! FALSE: zequal objects are equal

  use zdata

  implicit none

  !-- Forward declarations
  logical:: are_1d_int_arrays_different

  logical:: different
  type(zequal), intent(in):: zequal_1, zequal_2

  different = .TRUE.

  !-- Compare total number of zgrps
  ! if (zequal_1%ng /= zequal_2%ng) then
  !   return
  ! End if

  !-- Compare number of zgrps with only one atom
  ! if (zequal_1%eng /= zequal_2%eng) then
  !   return
  ! End if

  !-- Compare mapping of atoms onto zgrps
  if (are_1d_int_arrays_different(zequal_1%ord, zequal_1%nat, zequal_2%ord, zequal_2%nat)) then
    return
  End if

  !-- Compare grps
  ! if (size(zequal_1%grp) /= size(zequal_2%grp)) then
  !   return
  ! else
  !   do i = 1, zequal_1%ng
  !     if (are_zgrps_different(zequal_1%grp(i), zequal_2%grp(i))) then
  !       return
  !     End if
  !   End do
  ! end if

  different = .FALSE.
  return
End function are_zequals_different

subroutine update_zequal_ref(src, dest)
  use zdata

  implicit none

  type(zequal), intent(in):: src
  type(zequal), intent(out):: dest

  call dest%allocate(src%nat)
  dest%ng  = src%ng
  dest%grp = src%grp
  dest%ord = src%ord
  call dest%geteng
End subroutine update_zequal_ref

function are_1d_int_arrays_different(arr_1, arr_1_dim, arr_2, arr_2_dim) result(different)
  implicit none

  logical:: different
  integer, intent(in):: arr_1_dim, arr_2_dim
  integer, dimension(arr_1_dim), intent(in):: arr_1
  integer, dimension(arr_2_dim), intent(in):: arr_2

  different = .TRUE.

  if (size(arr_1) == size(arr_2)) then
    if (all(arr_1 .eq. arr_2)) then
      different = .FALSE.
    end if
  end if

  return
end function are_1d_int_arrays_different

function are_zgrps_different(zgrp_1, zgrp_2) result(different)
  use zdata

  implicit none

  ! Forward declarations
  logical:: are_1d_int_arrays_different

  logical:: different
  type(zgrp), intent(in):: zgrp_1, zgrp_2

  different = .TRUE.

  if (.NOT. are_1d_int_arrays_different(zgrp_1%mem, zgrp_1%nm, zgrp_2%mem, zgrp_2%nm)) then
    different = .FALSE.
  End if

  return
end function are_zgrps_different

subroutine sort_zgrps(zgrps, zgrps_dim)
  use zdata

  implicit none

  integer, intent(in):: zgrps_dim
  type(zgrp), dimension(zgrps_dim), intent(inout):: zgrps
  integer:: i, j, comp, comp_pos

  do i = 1, zgrps_dim-1
    comp = zgrps(i)%mem(1)
    comp_pos = i
    do j = 2, zgrps_dim
      if (comp > zgrps(j)%mem(1)) then
        comp = zgrps(j)%mem(1)
        comp_pos = j
      End if
    End do
    if (comp_pos /= i) then
      call swap_pair(zgrps, zgrps_dim, i, comp_pos)
    End if
  End do
End subroutine sort_zgrps

subroutine swap_pair(zgrps, zgrps_dim, pos1, pos2)
  use zdata

  implicit none

  integer, intent(in):: zgrps_dim
  type(zgrp), dimension(zgrps_dim), intent(inout):: zgrps
  integer, intent(in):: pos1, pos2
  type(zgrp):: temp_zgrp

  temp_zgrp = zgrps(pos1)
  zgrps(pos1) = zgrps(pos2)
  zgrps(pos2) = temp_zgrp
  return
End subroutine swap_pair


!===============================================================!
! get a neighbour list as an array
!===============================================================!
subroutine reactorneighbours(nat,at,xyz,bond,ndim,ngh)
    use crest_parameters, only: wp, sp
    use miscdata, only: rcov
    use utilities
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in) :: ndim
    real(wp),intent(inout) :: bond(nat,nat)
    integer,intent(out) :: ngh(ndim)  !neighbour list, one dimensional

    real(wp),allocatable :: cn(:)
    integer :: i,j,k,l
    integer :: icn
    real(wp) :: rcn,rcn2

    ngh=0 !reset
    allocate(cn(nat))
    call xcoord2(nat,at,xyz,rcov,cn,800.0_wp,bond) 

    do i=1,nat
      rcn=floor(cn(i))
      rcn2= cn(i) - rcn
      if(rcn2.gt. 0.6d0)then
         icn=nint(cn(i))
      else
         icn=nint(rcn)
      endif
      if(icn.ge.6)cycle !skip highly coordinated molecules
      do k=1,icn
         j=maxloc(bond(:,i),1)
         bond(j,i)=0.0d0
         l=lin(j,i)
         ngh(l) = 1  !atoms i and j are neighbours
      enddo
    enddo

    deallocate(cn)
    return
end subroutine reactorneighbours

!===============================================================!
! Write directories for post optimization
!===============================================================!
subroutine reactorreopt(env,nat,at,nall,xyz,taken,frags,ndirs)
    use crest_parameters, only: wp, bohr
    use crest_data
    use zdata
    use iomod
    use strucrd, only: wrxyz
    implicit none
    type(systemdata) :: env
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    integer,intent(in) :: nall
    real(wp),intent(in) :: xyz(3,nat,nall) !should be in Bohr
    logical,intent(in) :: taken(nall)
    type(zequal),intent(in) :: frags(nall)
    integer,intent(out) :: ndirs

    !real(wp),allocatable :: chrgs(:,:)
    integer :: i,j,k,l,m,n
    integer :: io
    logical,allocatable :: mask(:)

    character(len=512) :: thispath,optpath
    character(len=256) :: ctmp

    character(len=:),allocatable :: xnam
    character(len=:),allocatable :: xnam2
    character(len=1024) :: jobcall

    write(*,*)
   !---- get current path
    call getcwd(thispath)

    xnam='struc.xyz'
    ndirs = 0
   !--- maybe for a future version we want atom charces from SPs
    !call smallhead('singlepoints on taken structures')
    !optpath='CHARGESP'
    !call rmrf(optpath)
    !io = makedir(trim(optpath))
    !call chdir(trim(optpath))
    !call getcwd(optpath)
    !do i=1,nall
    ! if(taken(i))then
    !   ndirs=ndirs+1
    !   write(ctmp,'(a,i0)')'TMPCHRG',ndirs
    !   io = makedir(trim(ctmp))
    !   if(env%chrg .ne. 0)then
    !     open(newunit=ich,file=trim(ctmp)//'/'//'.CHRG')
    !      write(ich,*)env%chrg
    !     close(ich)
    !   endif
    !   if(env%uhf .ne. 0)then
    !     open(newunit=ich,file=trim(ctmp)//'/'//'.UHF')
    !      write(ich,*)env%uhf
    !    close(ich)
    !   endif    
    !   open(newunit=ich,file=trim(ctmp)//'/'//xnam)
    !   call wrxyz(ich,nat,at,xyz(:,:i)*bohr) 
    !   close(ich)
    ! else  
    !   cycle
    ! endif
    !enddo
    !--- SPs would go here
    !--- read 'charges' file afterwards
    !    allocate(chrgs(nat,nall)) ...etc.
    call chdir(trim(thispath)) !go back

   !--- Otherwise, do the fragments
    write(*,*)
    call smallhead('optimization of fragments')
    allocate(mask(nat))
    ndirs = 0 !reset
    optpath='OPTIM'
    call rmrf(optpath)
    io = makedir(trim(optpath))
    call chdir(trim(optpath))
    call getcwd(optpath)
    ILOOP : do i=1,nall
      if(taken(i))then
          k=frags(i)%ng
          JLOOP : do j=1,k
            l=frags(i)%grp(j)%nm
            if(l.le.1)cycle !cycle single atoms
            ndirs=ndirs+1
            write(ctmp,'(a,i0)')'TMPFRG',ndirs
            io = makedir(trim(ctmp))
    !--- currently the charge for all systems is assumed = 0 (neutral)
           mask=.false.
           do m=1,l
            n=frags(i)%grp(j)%mem(m)
            mask(n)=.true.
           enddo
           xnam2=trim(ctmp)//'/'//xnam
           call wrxyz(xnam2,nat,at,xyz(:,:,i)*bohr,mask)
          enddo JLOOP
      else
          cycle
      endif
    enddo ILOOP
    deallocate(mask)

    ! setting the threads for correct parallelization
    call new_ompautoset(env,'auto',ndirs,i,j)
    write(*,'(1x,''Starting optimization of reactor products'')')
    write(*,'(1x,i0,'' jobs to do.'')')ndirs

    write(jobcall,'(a,1x,a,1x,a,1x,a)')trim(env%ProgName),xnam, &
    &  trim(env%gfnver),'--opt >xtb.out 2>>xtb.out'
    call opt_OMP_loop(ndirs,'TMPFRG',jobcall,env%niceprint)
    write(*,'(/,1x,a)') 'done.'
    call chdir(trim(thispath)) !return

    return
end subroutine reactorreopt


subroutine collectproducts(optdir,base,ndirs,oname,iso)
    use crest_parameters, only: wp, bohr
    use iomod
    use strucrd, only: wrxyz,rdnat,rdcoord
    implicit none
    character(len=*),intent(in) :: optdir
    character(len=*),intent(in) :: base
    integer,intent(in) :: ndirs
    character(len=*),intent(in) :: oname
    character(len=:),allocatable :: path
    logical :: iso

    logical :: ex
    integer :: i,k
    integer :: ich
    integer :: natiso
    character(len=256) :: atmp

    real(wp),allocatable :: xyz(:,:)
    integer,allocatable :: at(:)
    character(len=40),allocatable :: sumformulas(:)
    character(len=40) :: sumform
    real(wp),allocatable :: energies(:)
    integer,allocatable :: nats(:)
    logical,allocatable :: taken(:)

    real(wp) :: de
    real(wp),parameter :: ethr = 0.05_wp  !energy threshold for same structure
    real(wp),parameter :: autokcal = 627.509541_wp

    allocate(sumformulas(ndirs))
    allocate(energies(ndirs), source= 0.0_wp)
    allocate(nats(ndirs), source = 0)
    allocate(taken(ndirs), source=.false.)

    inquire(file=oname,exist=ex)
    if(ex) call remove(oname)
    open(newunit=ich,file=oname)
    
    path=''
    do i=1,ndirs
       write(atmp,'(a,i0)')base,i
       path=optdir//'/'//trim(atmp)//'/'//'xtbopt.xyz'
       inquire(file=path,exist=ex)
       if(.not.ex) cycle !maybe some optimizations failed. they are cycled through.
       call rdnat(path,nats(i))
       if(i==1)natiso=nats(i)
       allocate(xyz(3,nats(i)),at(nats(i)))
       call rdcoord(path,nats(i),at,xyz,energies(i))

       sumformulas(i) = sumform(nats(i),at)

       !--- check for duplicates based on composition and total energy difference
       taken(i) = .true.
       do k=1,i-1
          if(.not.taken(k))cycle
          if(nats(i).ne.nats(k)) cycle
          if(trim(sumformulas(i))==trim(sumformulas(k)))then
            de=(energies(i)-energies(k))*autokcal
            if(abs(de).le.ethr) taken(i)=.false.
          else
            cycle
          endif
       enddo
       !--- exclude non-isomers?
       if(iso)then
           if(nats(i).ne.natiso) taken(i)=.false. 
       endif

       !--- if the structure was taken, write it to the collective ensemble
       if(taken(i))then
        call wrxyz(ich,nats(i),at,xyz*bohr,energies(i))
       endif    

       deallocate(at,xyz)
    enddo
    close(ich)


    write(*,*)
    call smallhead('reactor products summary')
    k=0
    write(*,'(1x,a,2x,a,4x,a,8x,a)')'structure','#atoms','Etot','composition'
    do i=1,ndirs
       if(taken(i))then
       k=k+1 
       write(*,'(1x,i6,3x,i6,1x,f16.8,3x,a)') k,nats(i),energies(i), &
      & trim(sumformulas(i))     
       endif
    enddo

    write(*,*)
    write(*,'(1x,a,a,a)')'Structures written to file "',trim(oname),'"'

    deallocate(taken,nats,energies,sumformulas)
    return
end subroutine collectproducts
