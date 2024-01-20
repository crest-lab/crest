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
!> subroutine crest_crossing
!> Read in an ensemble file and optimize all structures
!> perform a crossing of z-matrices to generate new structures
!> The procedure is loosely inspired by genetic algorithms
!>------------------------------------------------------
subroutine crest_crossing(env,maxgen,fname,maxpairs)
  use crest_parameters,only:wp,stdout,bohr,autokcal
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env
  integer,intent(in)             :: maxgen
  character(len=*),intent(in),optional   :: fname
  real(wp),intent(in),optional :: maxpairs
  !> LOCAL
  integer :: i,j,k,l,io,ich,c
  logical :: pr,wr,ex
!========================================================================================!
  real(wp) :: rthr,ewin,cthr

  character(len=:),allocatable :: ensnam
  integer :: nat,nall,nalltmp,maxgen2
  real(wp),allocatable :: eread(:),erel(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  real(wp) :: percent,maxtmp
!========================================================================================!
!>--- check for the ensemble file
  if (present(fname)) then
    ensnam = trim(fname)
  else
    ensnam = env%ensemblename
  end if
  inquire (file=ensnam,exist=ex)
  if (.not.ex) then
    write (stdout,*) 'no ensemble file found by crest_crossing()'
    return
  end if

!>---- read the input ensemble
  call rdensembleparam(ensnam,nat,nall)
  if (nall .lt. 1) return
  allocate (xyz(3,nat,nall),at(nat),eread(nall))
  call rdensemble(ensnam,nat,nall,at,xyz,eread)
  if (nall .lt. 2) then
    write (stdout,*) 'Not enough structures to perform GC!'
    return
  end if

!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: crest_crossing requires coordinates in Bohrs
  xyz = xyz/bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

!>--- set OMP parallelization, use max number of threads
  call ompautoset(env%threads,4,env%omp,env%MAXRUN,0) 

!========================================================================================!
!>--- thresholds
  ewin = env%ewin           !> energy window
  rthr = env%rthr*2.0d0     !> standard RMSD threshold
  cthr = 0.3d0              !> CN clash threshold
  if(present(maxpairs))then !> dynamically adjust ewin in order to fit maxpairs
    allocate(erel(nall))
    call etoerel(nall,eread,erel,autokcal)
    nalltmp = count( erel(:) < ewin)
    maxtmp = (float(nalltmp)*(float(nalltmp)-1.0_wp))/2.0_wp
    if( maxtmp > maxpairs)then
      do while (maxtmp > maxpairs)
        ewin = ewin - 0.01_wp
        nalltmp = count( erel(:) < ewin)
        maxtmp = (float(nalltmp)*(float(nalltmp)-1.0_wp))/2.0_wp
      enddo
    endif
    deallocate(erel)
  endif
!========================================================================================!
  !>--- printout header
  write (stdout,*)
  write (stdout,'(5x,''========================================'')')
  write (stdout,'(5x,''|        Structure Crossing (GC)       |'')')
  write (stdout,'(5x,''========================================'')')
  write (stdout,*)
  call ompprint_intern
  write (stdout,'(a,a)')  'input  file name : ',trim(ensnam)
  write (stdout,'(a,i8)') 'number of atoms                :',nat
  write (stdout,'(a,i8)') 'number of points on xyz files  :',nall
  if(present(maxpairs).and.(ewin.ne.env%ewin))then
  write (stdout,'(a,es8.1)') 'max. # of parent structures    :',maxpairs
  write (stdout,'(a,f8.2)') 'adjusted energy window  /kcal  :',ewin
  else
  write (stdout,'(a,f8.2)') 'conformer energy window  /kcal :',ewin
  endif
  write (stdout,'(a,f8.4)') 'CN per atom difference cut-off :',cthr
  write (stdout,'(a,2f8.4)') 'RMSD threshold (Ang, Bohr)     :',rthr,rthr/bohr
  write (stdout,'(a,1x,i8)') 'max. # of generated structures :',maxgen

  maxgen2 = maxgen
  call crossing(nat,nall,at,xyz,eread,ewin,rthr,cthr,maxgen2)


  deallocate (eread,at,xyz)
!========================================================================================!
  return
end subroutine crest_crossing

!========================================================================================!
!========================================================================================!
subroutine crossing(nat,nall,at,xyz,er,ewin,rthr,cthr,maxgen)
  use crest_parameters
  use ls_rmsd
  use strucrd
  use miscdata, only: rcov
  implicit none
  !> INPUT
  integer,intent(in)  :: nat,nall           !> number of atoms, number of structures
  integer,intent(in)  :: at(nat)            !> atomic numbers
  real(wp),intent(inout) :: xyz(3,nat,nall) !> cartesian coordinates in Bohr(!!!) for ensemble
  real(wp),intent(in) :: er(nall)           !> energies in Eh for each structure
  integer,intent(inout) :: maxgen           !> max. number of generated structures
  real(wp),intent(in) :: ewin               !> energy window
  real(wp),intent(in) :: rthr               !> rmsd threshold in Ang
  real(wp),intent(in) :: cthr               !> CN clash threshold
  !> LOCAL
  logical :: pr
  real(sp),allocatable :: zmat(:,:,:),zref(:,:) !> zmat in single precision
  integer,allocatable  :: na(:),nb(:),nc(:)
  real(wp),allocatable :: cdum(:,:),zdum(:,:),xyzref(:,:)
  real(wp),allocatable :: xyzgen(:,:,:),rms(:)
  real(wp),allocatable :: erel(:)
  real(wp),allocatable :: cnref(:)
  real(wp),allocatable :: gdum(:,:),Udum(:,:),xdum(:),ydum(:)  !> rmsd dummy stuff
  integer,allocatable  :: ind(:)
  real(wp) :: rthrbohr,rthrbohr100,rthrzmat,emin,rmsdavg,rval,rval2
  real(wp) :: nmaxref,ncount,ncheck,ncheckstep
  integer :: ierr,ntaken,ident,rcount,maxgen2
  integer :: i,j,m,k,nwin,minpos
  character(len=80) :: atmp
  integer :: dumpio
  logical :: fail,rmsdcheck,stop_crossing,huge_number
  integer :: sdselect !> this is a function

  pr = .true.
  rmsdcheck = .true.

  !>--- molecule storage
  allocate (cdum(3,nat),zdum(3,nat),source=0.0_wp)
  allocate (xyzref(3,nat),source=0.0_wp)
  allocate (zmat(3,nat,nall),zref(3,nat),source=0.0_sp) !> zmat in single precision
  allocate (na(nat),nb(nat),nc(nat),source=0)
  allocate (gdum(3,3),Udum(3,3),xdum(3),ydum(3))

  !>--- parameter setup
  rthrbohr = rthr/bohr  !> Ang to Bohr
  rthrbohr100 = rthrbohr/100.0_wp
  allocate (erel(nall),cnref(nat),source=0.0_wp)
  ! call setrcov(rcov)
  call etoerel(nall,er,erel,autokcal)
  nwin = count((erel(:) < ewin),1)
  minpos = minloc(erel,1)
  if (pr) write (stdout,*) '# in E window',nwin
  nmaxref = (float(nwin)*(float(nwin)-1.0_wp))/2.0_wp

  !>--- lowest conformer to provide reference values
  xyzref(1:3,1:nat) = xyz(1:3,1:nat,minpos)  !> reference Cartesians
  call ycoord(nat,rcov,at,xyzref,cnref,100.0d0) !> refernce CNs
  call XYZINT(xyzref,nat,na,nb,nc,1.0d0,zdum)   !> z-mat connectivity
  zref(:,:) = real(zdum,sp)                    !> reference z-mat

  !>--- convert to internal coordinates
  !$OMP PARALLEL PRIVATE(i,zdum) &
  !$OMP SHARED(nall,nat,na,nb,nc,xyz,zmat)
  !$OMP DO
  do i = 1,nall
    call XYZGEO(xyz(:,:,i),nat,na,nb,nc,1.0d0,zdum)
    zmat(:,:,i) = real(zdum(:,:),sp) !> zmat in single precision
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  !>--- generation loop setup
  ierr = 0          !> discarded
  rmsdavg = 0.0_wp     !> avg. rmsd of generated structures
  rcount = 0          !> counting for avg
  ident = 0          !> counting of discarded identical
  ntaken = 0          !> counting of selected structures
  maxgen2 = maxgen*2   !> check double the amount of structures needed
  ncount = 0.0_wp      !> continous counter
  ncheck = 0.1_wp     !> checkpoint for printout
  ncheckstep = 0.25_wp !> steps for checkpoint printout
  stop_crossing = .false. !> parameter for early termination
  huge_number = .false.
  if(nmaxref > 10*maxgen2)then
    huge_number = .true.
    maxgen2 = maxgen*5
  endif
  allocate (xyzgen(3,nat,maxgen2),rms(maxgen2),source=0.0_wp)
  

!$OMP PARALLEL PRIVATE(i,j,k,m, zdum,cdum, rval,rval2, fail) &
!$OMP PRIVATE(Udum,xdum,ydum,gdum) &
!$OMP SHARED(erel,ewin,nat,at,zref,zmat,na,nb,nc,cnref ) &
!$OMP SHARED(cthr,rthrbohr,rthrbohr100,rms,xyzgen,maxgen,maxgen2) &
!$OMP SHARED(ncount,nmaxref,ncheck,ncheckstep) &
!$OMP SHARED(ierr,rcount,ident,rmsdavg,ntaken,rmsdcheck,stop_crossing,huge_number )
!$OMP DO
  do i = 1,nall
    if (erel(i) > ewin) cycle
    do j = 1,i-1
      if(stop_crossing) cycle !> since we must not jump out of an OMP loop, use this workaround
      if (erel(j) > ewin) cycle
      ncount = ncount+1.0_wp
      !>-- new structure generation
      do m = 1,nat
        !>-- new internal coords
        do k = 1,3
          zdum(k,m) = real(zref(k,m)+zmat(k,m,j)-zmat(k,m,i),wp)
        end do
        !>-- prevent linear bends
        if (pi-zdum(2,m) .lt. 0.001_wp) zdum(2,m) = pi-0.005_wp
        if (zdum(2,m) .lt. 0.001_wp) zdum(2,m) = 0.005_wp
      end do
      call GMETRY(nat,zdum,cdum,na,nb,nc,fail)
      !>--- checks
      if (fail) then
        !$omp atomic
        ierr = ierr+1
        cycle
      end if
      call ycoord2(nat,rcov,at,cdum,cnref,100.d0,cthr,fail) !> CN clashes
      if (fail) then
        !$omp atomic
        ierr = ierr+1
        cycle
      end if
      call rmsd(nat,cdum,xyzref,0,Udum,xdum,ydum,rval,.false.,gdum)
      !call poor_mans_rmsd(nat,xyzref,cdum,rval)
      if (rval <= rthrbohr) then
        !$omp atomic
        ierr = ierr+1
        cycle
      else
        !$omp atomic
        rcount = rcount+1
        !$omp atomic
        rmsdavg = rmsdavg+rval
      end if
      !$omp critical
      if (ntaken < 1) then
        ntaken = ntaken+1
        xyzgen(:,:,ntaken) = cdum
        rms(ntaken) = rval
      else
        !> check if we have the generated structure already
        fail = .false.
        do m = 1,ntaken
          !> check only those with a similar rmsd to the reference
          if (abs(rms(m)-rval) < rthrbohr100 ) then
            if (rmsdcheck) then !> variant 1: Cartesian RMSDs
              call rmsd(nat,cdum,xyzgen(:,:,m),0,Udum,xdum,ydum,rval2,.false.,gdum)
              if (rval2 < rthrbohr) then
                ident = ident+1
                fail = .true.
                exit
              end if
            else !> variant 2: crude approximation, discard it
              ident = ident+1
              fail = .true.
              exit
            end if
          end if
        end do
        if (fail) then
          continue
        else
          if (ntaken < maxgen2) then
            ntaken = ntaken+1
            xyzgen(:,:,ntaken) = cdum
            rms(ntaken) = rval
          elseif(huge_number) then
            stop_crossing = .true.
          else
            k = sdselect(ntaken,rms,rval)
            if (k > 0) then
              xyzgen(:,:,k) = cdum
              rms(k) = rval
            end if
          end if
        end if
      end if
      !$omp end critical
    end do
    !>-- checkpoint printout
    !$omp critical
    if(huge_number)then
    if (pr .and.((float(ntaken)/float(maxgen2)) >= ncheck).and. .not.stop_crossing) then
      ncheck = ncheck+ncheckstep
      if (ncheck > 1.0_wp) ncheck = 1.0_wp
      write (stdout,'(f6.1," % done")') (float(ntaken)/float(maxgen2))*100.0_wp
    end if
    else
    if (pr .and.((ncount/nmaxref) >= ncheck)) then
      ncheck = ncheck+ncheckstep
      if (ncheck > 1.0_wp) ncheck = 1.0_wp
      write (stdout,'(f6.1," % done")') (ncount/nmaxref)*100.0_wp
    end if
    endif
    !$omp end critical
  end do
!$OMP END DO
!$OMP END PARALLEL
  if(pr.and.(ntaken > 0)) write (stdout,'(" finished.")')
  rmsdavg = rmsdavg/float(rcount)

  if (pr) then
    if (ntaken > 0) then
      write (stdout,'(1x,a,f12.5)') 'average rmsd w.r.t input     :',rmsdavg
      write (stdout,'(1x,a,i12)') 'number of clash discarded    :',ierr
      write (stdout,'(1x,a,i12)') 'removed identical structures :',ident
    else
      write (stdout,'(1x,a)') 'no new structures generated'
    end if
  end if

  if (ntaken > 0) then
    !>--- determin order for dumping
    rmsdavg = sum(rms(1:ntaken))/float(ntaken)
    allocate (ind(ntaken),source=0)
    do i = 1,ntaken
      ind(i) = i
      rms(i) = -(rms(i)-rmsdavg)**2
    end do
    call qsort(rms(1:ntaken),1,ntaken,ind)
    maxgen = min(maxgen,ntaken)
    !>--- dump the new structures to file
    !xyzgen = xyzgen*bohr !> to Angström
    open (newunit=dumpio,file='confcross.xyz')
    do i = 1,maxgen
      k = ind(i)
      xyzgen(:,:,k) = xyzgen(:,:,k)*bohr !> to Angström
      write (atmp,'(1x,f12.6,1x,a)') rms(i),'!GC'
      call wrxyz(dumpio,nat,at,xyzgen(:,:,k),trim(atmp))
    end do
    close (dumpio)
    if (pr) then
      write (stdout,'(/,1x,i0,1x,a,/)') maxgen,'structures written to confcross.xyz'
    end if
  end if

  !>--- cleanup
  if (allocated(ind)) deallocate (ind)
  if (allocated(ydum)) deallocate (ydum)
  if (allocated(xdum)) deallocate (xdum)
  if (allocated(Udum)) deallocate (Udum)
  if (allocated(gdum)) deallocate (gdum)
  if (allocated(zmat)) deallocate (zmat)
  if (allocated(zdum)) deallocate (zdum)
  if (allocated(zref)) deallocate (zref)
  if (allocated(xyzref)) deallocate (xyzref)
  if (allocated(na)) deallocate (na)
  if (allocated(nb)) deallocate (nb)
  if (allocated(nc)) deallocate (nc)
  if (allocated(cnref)) deallocate (cnref)
  if (allocated(erel)) deallocate (erel)
  !if (allocated(rcov)) deallocate (rcov)
  if (allocated(cdum)) deallocate (cdum)
  if (allocated(xyzgen)) deallocate (xyzgen)

end subroutine crossing

!========================================================================================!
  function sdselect(n,rms,rval) result(pos)
    use crest_parameters
!> calculates SD of rms()
!> and checks which entry to replace by rval in order to
!> maximise the SD. returns this entry's number (or zero)
    implicit none
    integer :: pos
    integer,intent(in)  :: n
    real(wp),intent(in) :: rms(n)
    real(wp),intent(in) :: rval
    real(wp) :: avg,sdref,sdtmp
    integer :: i,j

    pos = 0
    avg = sum(rms)
    avg = avg/float(n)
    sdref = sum((rms(:)-avg)**2)
    sdref = sqrt(sdref/float(n))

    do i = 1,n
      avg = sum(rms)-rms(i)+rval
      avg = avg/float(n)
      sdtmp = sum((rms(:)-avg)**2)
      sdtmp = sdtmp-(rms(i)-avg)**2
      sdtmp = sdtmp+(rval-avg)**2
      sdtmp = sqrt(sdtmp/float(n))
      if (sdtmp > sdref) then
        sdref = sdtmp
        pos = i
      end if
    end do
  end function sdselect
!========================================================================================!

subroutine poor_mans_rmsd(nat,xyz_A,xyz_B,rmsd)
  use crest_parameters
  implicit none
  integer,intent(in) :: nat
  real(wp),intent(in) :: xyz_A(3,nat)
  real(wp),intent(in) :: xyz_B(3,nat)
  real(wp),intent(out) :: rmsd
  integer :: i,j,k,l
  real(wp) :: dum
  rmsd = 0.0_wp
  do i=1,nat
     dum = sum((xyz_A(1:3,i)-xyz_B(1:3,i))**2.0_wp)
     dum = dum/float(nat)  
     rmsd = rmsd + dum
  enddo
  rmsd = sqrt(rmsd) 

end subroutine poor_mans_rmsd
