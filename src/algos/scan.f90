!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022 Philipp Pracht
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

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> Implementation of a scans
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
!> Input/Output:
!>  env  -  crest's systemdata object
!>  tim  -  timer object
!>-----------------------------------------------
subroutine crest_scan(env,tim)
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
  use geo
  use iomod, only: makedir
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich
  logical :: pr,wr
!========================================================================================!
  type(calcdata) :: calc
  type(calcdata) :: calcclean
  real(wp) :: accuracy,etemp

  real(wp) :: energy
  real(wp),allocatable :: grad(:,:)

!========================================================================================!
  write (*,*)
  !call system('figlet scan')
  write (stdout,*) " ___  ___ __ _ _ __  "
  write (stdout,*) "/ __|/ __/ _` | '_ \ "
  write (stdout,*) "\__ \ (_| (_| | | | |"
  write (stdout,*) "|___/\___\__,_|_| |_|"
  write (stdout,*) "                     "
!========================================================================================!
  call new_ompautoset(env,'max',0,i,j)
  call ompprint_intern()
  call tim%start(14,'Coordinate scan')
!========================================================================================!
  call env%ref%to(mol)
  write (*,*)
  write (*,*) 'Input structure:'
  call mol%append(stdout)
  write (*,*)
!========================================================================================!

  if (env%calc%nscans < 1) then
    stop 'no scans to run!'
  end if

!========================================================================================!
  allocate (grad(3,mol%nat),source=0.0_wp)
  calc = env%calc
  calcclean = env%calc

  !>--- initialize scanning
  io = makedir('scanfiles') 
  call initscans(mol,calc)

  !>--- runscan
  i = 1
  !calcclean%elog = 'scan.elog'
  !open(newunit=calcclean%eout_unit , file=calcclean%elog)
  call runscan(mol,calc,calcclean,i)
  !close(calcclean%eout_unit) 

  deallocate (grad)
!========================================================================================!
  call tim%stop(14)
  return

!========================================================================================!
!========================================================================================!
contains
!========================================================================================!
!========================================================================================!

  subroutine initscans(mol,calc)
    implicit none
    !> INPUT
    type(coord)    :: mol
    type(calcdata) :: calc
    !> LOCAL
    integer :: i,j,k,l,kloc(1)
    integer :: atms(4)
    integer :: nsteps
    real(wp),allocatable :: tmppoints(:)
    real(wp) :: d,delx,mi,ma,dref
    character(len=10) :: scannmbr
    type(constraint) :: constr

    !>--- allocate points for each scan
    do i = 1,calc%nscans
      if (allocated(calc%scans(i)%points)) deallocate (calc%scans(i)%points)
      nsteps = calc%scans(i)%steps
      allocate (calc%scans(i)%points(nsteps),source=0.0_wp)

      allocate (tmppoints(nsteps),source=0.0_wp)
      select case (calc%scans(i)%type)
      case (1) !>-- distance
        mi = min(calc%scans(i)%maxval,calc%scans(i)%minval)
        ma = max(calc%scans(i)%maxval,calc%scans(i)%minval)
        d = ma - mi
        tmppoints(1) = mi
        if (nsteps > 1) then
          delx = abs(d) / float(nsteps - 1)
          do j = 2,nsteps
            tmppoints(j) = tmppoints(j - 1) + delx
          end do
        end if
        write (stdout,*) 'Adding scan of interatomic distance:'
        write (stdout,'(2x,a,i5,i5)') 'atoms:',calc%scans(i)%atms(1:2)
        write (stdout,'(2x,a,f10.4,a,f10.4,a)') 'distance:',mi,' to',ma,' a.u.'
        write (stdout,'(2x,a,1x,i0)') 'steps:',nsteps
        call constr%deallocate()
        atms(1:2) = calc%scans(i)%atms(1:2)
        dref = mol%dist(atms(1),atms(2))
        call constr%bondconstraint(atms(1),atms(2),dref,calc%scansforce)
        !>-- add the constraint initialized to the starting geometry
        !>-- and track which constraint corresponds to the scanning coordinate
        call calc%add(constr)
        calc%scans(i)%constrnmbr = calc%nconstraints
        k = minloc(abs(tmppoints(:) - dref),1)
        !call shiftpoints(nsteps,k,tmppoints,calc%scans(i)%restore,.false.)
        calc%scans(i)%points = tmppoints
        !write(*,*) calc%scans(i)%points 
      case (3) !>-- dihedral
        mi = calc%scans(i)%minval*degtorad
        ma = calc%scans(i)%maxval*degtorad
        if(mi==0.0_wp .and. ma==0.0_wp)then
          mi=-180.000_wp*degtorad  + 1.0d-10  !> avoid *exactly* 180° !!!
          ma= 180.000_wp*degtorad  - 1.0d-10  !> avoid *exactly* 180° !!!
        endif
        d = ma - mi
        tmppoints(1) = mi
        if (nsteps > 1) then
          delx = d / float(nsteps - 1)
          do j = 2,nsteps
            tmppoints(j) = tmppoints(j - 1) + delx
          end do
        end if
        write (stdout,*) 'Adding scan of dihedral angle:'
        write (stdout,'(2x,a,4i5)') 'atoms:',calc%scans(i)%atms(1:4)
        write (stdout,'(2x,a,f10.4,a,f10.4,a)') 'angles:',mi*radtodeg,' to',ma*radtodeg,' deg'
        write (stdout,'(2x,a,1x,i0)') 'steps:',nsteps
        call constr%deallocate()
        atms(1:4) = calc%scans(i)%atms(1:4)
        dref = mol%dihedral(atms(1),atms(2),atms(3),atms(4))
        call constr%dihedralconstraint(atms(1),atms(2),atms(3),atms(4),dref*radtodeg,calc%scansforce)
        !>-- add the constraint initialized to the starting geometry
        !>-- and track which constraint corresponds to the scanning coordinate
        call calc%add(constr)
        calc%scans(i)%constrnmbr = calc%nconstraints
        k = polar_minloc(tmppoints,nsteps,dref)
        call shiftpoints(nsteps,k,tmppoints,calc%scans(i)%restore, .true.)
        calc%scans(i)%points = tmppoints !> reference dihedrals in rad
        write(scannmbr,'(i0)')i
        open(newunit=j,file='scanfiles/grid.'//trim(scannmbr)//'.grid')
        write(j,*) "dihedral scan "//trim(scannmbr)
        do k=1,nsteps
          write(j,*) calc%scans(i)%points(k)*radtodeg + 180.0_wp
        enddo
        close(j) 
        !write(*,*) calc%scans(i)%points * radtodeg
      end select
      deallocate (tmppoints)
    end do

    write(*,*) 'constraints',calc%nconstraints

    !>--- set calculations to 1 for the geometry generation
    calc%ncalculations = 1

  end subroutine initscans

  subroutine shiftpoints(n,k,points,k2,circular)
    implicit none
    integer,intent(in) :: n
    integer,intent(in) :: k
    real(wp),intent(inout) :: points(n)
    integer,intent(out) :: k2
    logical,intent(in) :: circular
    integer :: i,j,l
    real(wp) :: tmp
    real(wp),allocatable :: ptmp(:)
    k2 = 1
    if (k > 1) then
      l = k - 1
      k2 = n - l + 1
      allocate(ptmp(l), source=0.0d0)
      if(circular)then
       ptmp(1:l) = points(1:l)
      else
      j = l
      do i=1,l
        ptmp(i) = points(j)
        j = j - 1
      enddo
       
      endif
      do i = 1,l
        do j = 2,n
          points(j - 1) = points(j)
        end do
      end do
      points(k2:n) = ptmp(1:l)
      if(circular) k2 = 0
    end if
  end subroutine shiftpoints

  function polar_minloc(philist,n,phi) result(pos)
     implicit none
     integer :: pos
     real(wp),intent(in) :: phi !> in rad
     integer,intent(in)  :: n
     real(wp),intent(in) :: philist(n) !> in rad
     real(wp),allocatable :: cd(:)
     integer :: i
     allocate(cd(n),source = 0.0_wp)
     do i=1,n
       cd(i) = &
     &    ( cos(philist(i)) - cos(phi) )**2 &
     & +  ( sin(philist(i)) - sin(phi) )**2 
       cd(i) = sqrt(cd(i))
     enddo
     pos = minloc(cd,1)
     deallocate(cd)
  end function polar_minloc 

!========================================================================================!

  recursive subroutine runscan(mol,calc,calcclean,current)
    implicit none
    !> INPUT
    type(coord)    :: mol
    type(calcdata) :: calc
    type(calcdata) :: calcclean
    integer,intent(in) :: current
    !> LOCAL
    type(coord) :: molnew,molbackup
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    real(wp) :: val
    character(len=500) :: scantrj
    character(len=20)  :: scannmbr
    integer :: dumpid,dumpid2,dumpid3,k,j,j2,next, io


    if( current > calc%nscans ) then !> recursive termination
      return
    else if( current == calc%nscans )then
      scantrj = 'scanfiles'//sep//'scan'
      do j=1,calc%nscans-1
        write(scannmbr,'(a,i0)') '_',calc%scans(j)%currentstep 
        scantrj = trim(scantrj)//trim(scannmbr)
      enddo
      open(newunit=dumpid, file=trim(scantrj)//'.xyz')
      calcclean%elog = trim(scantrj)//'.elog'
      open(newunit=calcclean%eout_unit , file=calcclean%elog) 
      open(newunit=dumpid2, file=trim(scantrj)//'.grid')
    endif    
    
    allocate (grad(3,mol%nat),source=0.0_wp)
    !>-- geometry optimization
    pr = .false. !> stdout printout
    wr = .true. !> write crestopt.log
    molbackup = mol
    do j=1,calc%scans(current)%steps 
      !write(*,*) current, calc%scans(current)%steps, j

      calc%scans(current)%currentstep = j
      !>-- restore starting point?
      if( j == calc%scans(current)%restore)then
       mol = molbackup
      endif
      !>-- get the associated constraint
      k = calc%scans(current)%constrnmbr
      !>-- get the new constraint value
      val = calc%scans(current)%points(j)
      !>-- set it to this value
      calc%cons(k)%ref(1) = val
      !>-- run the optimization
      call optimize_geometry(mol,molnew,calc,energy,grad,pr,wr,io)
      mol = molnew

      !>-- down the rabbit hole
      next = current + 1
      call runscan(mol,calc,calcclean,next) 

      !>--- for the lowest level iterations do the printouts
      if( current == calc%nscans )then
       
       write(stdout,'(10i5)') calc%scans(1:calc%nscans)%currentstep
       do j2=1,calc%nscans
         k = calc%scans(j2)%currentstep
         val = calc%scans(j2)%points(k) * radtodeg
         write(dumpid2,'(f9.4)',advance='no') val + 180.0_wp
       enddo
       write(dumpid2,*)
       energy = 0.0_wp
       call engrad(mol,calcclean,energy,grad,io)
       call molnew%appendlog(dumpid, energy)
      endif
    enddo
    deallocate (grad)
    mol = molbackup

    if( current == calc%nscans )then
    close(dumpid)
    close(calcclean%eout_unit)
    close(dumpid2)
    endif
  end subroutine runscan

!========================================================================================!
!========================================================================================!
end subroutine crest_scan
