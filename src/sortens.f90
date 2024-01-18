!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Philipp Pracht
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

!! ------------------------------------------------------------------
!  Sort ensemble file according to different criteria
!! ------------------------------------------------------------------
subroutine sort_ens(sort,infile,verbose)
  use crest_parameters
  use crest_data
  use strucrd,only:rdensembleparam,rdensemble,wrxyz
  implicit none
  type(protobj)    :: sort

  character(len=*),intent(in) :: infile

  integer  :: n
  integer  :: nall
  real(wp),allocatable :: xyz(:,:,:)
  real(wp),allocatable :: xyz_new(:,:,:)
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: e_new(:)
  integer,allocatable  :: at(:)
  real(wp) :: e_rel

  real(wp),allocatable :: xdum(:,:)
  integer,allocatable  :: molvec(:)
  integer  :: frag
  integer  :: minl,minll(1)
  integer  :: i,k
  integer  :: ich,ochan
  real(wp) :: dE
  logical :: verbose
  character(len=128) :: atmp

  real(wp),parameter :: kcal = 627.5095_wp

  ochan = 6 !stdout

  !sort=env%ptb

  associate (popthr => sort%popthr, &
  & allowFrag => sort%allowFrag,ewin => sort%ewin,     &
  & protdeprot => sort%protdeprot,deprotprot => sort%deprotprot, &
  & nfrag => sort%nfrag,threshsort => sort%threshsort)

    call rdensembleparam(infile,n,nall)
    allocate (xyz(3,n,nall),eread(nall),at(n),xyz_new(3,n,nall),e_new(nall))
    call rdensemble(infile,n,nall,at,xyz,eread)

!--- sort ensemble for energy
    do i = 1,nall
      minll = minloc(eread)
      minl = minll(1)
      e_new(i) = eread(minl)
      xyz_new(:,:,i) = xyz(:,:,minl)
      eread(minl) = 0.0_wp
    end do

!--- sort out structures that dissociated in the process (if requested)
    if (.not.allowFrag) then
      k = 0
      allocate (xdum(3,n),molvec(n))
      do i = 1,nall
        xdum(:,:) = xyz_new(:,:,i)/bohr
        call mrec(frag,xdum,n,at,molvec)
        if (frag .gt. nfrag) then
          e_new(i) = 0.0_wp      ! structures with E=0.0 are not written to file
          k = k+1
        end if
      end do
      deallocate (molvec,xdum)
      if (k .gt. 0) then
        write (*,'(1x,''Structures sorted out due to dissociation:'',i5)') k
      end if
    end if

!--- apply threshold within tautomerization sorting
    if (threshsort) then
      do i = 1,nall
        if (e_new(i) .lt. 0.0_wp) then
          e_rel = (e_new(i)-minval(e_new))*kcal
          if (e_rel .gt. ewin) e_new(i) = 0.0_wp       !will not be included
        end if
      end do
    end if

!--- write the ensemble file sorted
    open (file=infile,newunit=ich)
    do i = 1,nall
      if (e_new(i) .lt. 0.0_wp) then
        write (atmp,*) e_new(i)
        call wrxyz(ich,n,at,xyz_new(:,:,i),atmp)
      end if
    end do
    close (ich)

!--- printout?
    if (verbose) then
      write (ochan,*)
      write (ochan,'(a)') '==================================================='
      write (ochan,'(a)') '============= ordered structure list =============='
      write (ochan,'(a)') '==================================================='
      write (ochan,'(a,a,a)') ' written to file <',trim(infile),'>'
      write (ochan,*)
      write (ochan,'('' structure    ΔE(kcal/mol)   Etot(Eh)'')')
      do i = 1,nall
        dE = (e_new(i)-e_new(1))*kcal
        write (ochan,'(i5,6x,F10.2,4x,F14.6)') i,dE,e_new(i)
      end do
      write (ochan,*)
    end if
!--- deallocate
    deallocate (e_new,xyz_new,at,eread,xyz)

  end associate
end subroutine sort_ens
