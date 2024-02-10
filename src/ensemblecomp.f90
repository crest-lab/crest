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

!================================================================================!
! ENSEMBLE COMPARISON FUNCTION.
! To use execute:
!   crest --compare <ensemble1> <ensemble2>
!================================================================================!
subroutine compare_ensembles(env)
  use iso_fortran_env,only:wp => real64
  use ls_rmsd
  use iomod
  use crest_data
  use strucrd,only:rdensembleparam,rdensemble
  use cregen_interface
  implicit none

  type(systemdata) :: env

  integer :: i,j,k,l,kk,ll
  integer :: be,ed
  integer :: tr1,tr2
  integer :: nat
  integer :: dum
  integer :: nat1,nat2
  integer :: nall1,nall2
  integer :: nconf1,nconf2
  integer :: nrot1,nrot2
  integer :: rcount
  integer :: ich,ich2

  real(wp) :: rthr
  real(wp) :: rval
  real(wp) :: min_1,min_2,min_tot

  real(wp),allocatable :: gdum(:,:),Udum(:,:),xdum(:),ydum(:)  ! rmsd dummy stuff
  real(wp),allocatable :: rmat(:,:),rmat2(:)
  real(wp),allocatable :: c1(:,:),c2(:,:)
  real(wp),allocatable :: xyz1(:,:,:),xyz2(:,:,:)
  real(wp),allocatable :: conf1(:,:,:),conf2(:,:,:)
  real(wp),allocatable :: eread1(:),eread2(:),en1(:),en2(:)
  integer,allocatable  :: at1(:),at2(:)
  integer,allocatable  :: b1(:),b2(:)
  integer,allocatable  :: e1(:),e2(:)
  integer,allocatable  :: rots1(:),rots2(:)

  character(len=128) :: ensname1,ensname1backup
  character(len=128) :: ensname2
  character(len=60)  :: track1,track2
  character(len=6)   :: connect

  logical :: matching,ex

  real(wp),parameter :: autokcal = 627.509541d0

  associate (ensemblename => env%ensemblename,ensemblename2 => env%ensemblename2, &
  & thresholds => env%thresholds,cgf => env%cgf,maxcompare => env%maxcompare)
    nat = env%nat

!----
    call compens_cleanup()

!---- check if two valid files were given
    inquire (file=ensemblename,exist=ex)
    if (.not.ex) then
      write (*,'(2x,a,a,a)') 'File <',trim(ensemblename),'> does not exist!'
      error stop
    end if
    inquire (file=ensemblename2,exist=ex)
    if (.not.ex) then
      write (*,'(2x,a,a,a)') 'File <',trim(ensemblename2),'> does not exist!'
      error stop
    end if

!---- settings
    allocate (gdum(3,3),Udum(3,3),xdum(3),ydum(3))
    rthr = env%rthr      ! RMSD threshold

    ensname1 = 'ensemble1.inp.xyz'
    ensname2 = 'ensemble2.inp.xyz'
    track1 = 'track.1'
    track2 = 'track.2'

!---- sort the input files with CREGEN
    env%confgo = .true. !needs to be activated here
!----- first file
    call smallhead('Sorting file <'//trim(ensemblename)//'>')
    ensname1backup = trim(ensemblename)

    call newcregen(env,0)

    call rename(trim(ensemblename)//'.sorted',trim(ensname1))
    call rename('.cretrack',track1)
    call rdensembleparam(trim(ensname1),nat1,nall1)
    allocate (xyz1(3,nat1,nall1),eread1(nall1),at1(nat1))
    call rdensemble(trim(ensname1),nat1,nall1,at1,xyz1,eread1)
    !---- read tracking information (which conformer has which rotamers from the files?)
    open (newunit=tr1,file=track1)
    read (tr1,*) nconf1  !ensemble 1 has "nconf1" conformers
    !to compare a limited number "maxcompare" conformers, "nconf1/2" is normalized to maxcompare
    !i.e. if nconf > maxcompare, nconf is set to maxcompare; if nconf < maxcompare, it is not changed
    dum = nconf1
    nconf1 = min(nconf1,maxcompare)
    !---- set up the mapping
    allocate (b1(nconf1),e1(nconf1),rots1(nconf1))
    do i = 1,nconf1
      read (tr1,*) l,b1(i),e1(i)
      rots1(i) = e1(i)-b1(i)+1
    end do
    close (tr1)
    write (*,'(2x,a,a,a,i0,a)') 'File <',trim(ensemblename),'> contains ',dum,' conformers.'
    if (dum .lt. maxcompare) then
      write (*,'(2x,a,i0,a)') 'All of the ',dum,' conformers will be taken for the comparison:'
    else
      write (*,'(2x,a,i0,a)') 'The ',nconf1,' lowest conformers will be taken for the comparison:'
    end if
    write (*,'(2x,a9,2x,a9)') 'conformer','#rotamers'
    do i = 1,nconf1
      write (*,'(i9,2x,i9)') i,rots1(i)
    end do
    write (*,*)

!----- second file
    ensemblename = ensemblename2 !done so that CREGEn sorts the correct file
    call smallhead('Sorting file <'//trim(ensemblename2)//'>')

    call newcregen(env,0)

    ensemblename = trim(ensname1backup)
    call rename(trim(ensemblename2)//'.sorted',trim(ensname2))
    call rename('.cretrack',track2)
    call rdensembleparam(trim(ensname2),nat2,nall2)
    allocate (xyz2(3,nat2,nall2),eread2(nall2),at2(nat2))
    call rdensemble(trim(ensname2),nat2,nall2,at2,xyz2,eread2)
!---- read tracking information (which conformer has which rotamers from the files?)
    open (newunit=tr2,file=track2)
    read (tr2,*) nconf2  !ensemble 2 has "nconf2" conformers
    !to compare a limited number "maxcompare" conformers, "nconf1/2" is normalized to maxcompare
    !i.e. if nconf > maxcompare, nconf is set to maxcompare; if nconf < maxcompare, it is not changed
    dum = nconf2
    nconf2 = min(nconf2,maxcompare)
!-----set up the mapping
    allocate (b2(nconf2),e2(nconf2),rots2(nconf2))
    do i = 1,nconf2
      read (tr2,*) l,b2(i),e2(i)
      rots2(i) = e2(i)-b2(i)+1
    end do
    close (tr2)
    write (*,'(2x,a,a,a,i0,a)') 'File <',trim(ensemblename2),'> contains ',dum,' conformers.'
    if (dum .lt. maxcompare) then
      write (*,'(2x,a,i0,a)') 'All of the ',dum,' conformers will be taken for the comparison:'
    else
      write (*,'(2x,a,i0,a)') 'The ',nconf2,' lowest conformers will be taken for the comparison:'
    end if
    write (*,'(2x,a9,2x,a9)') 'conformer','#rotamers'
    do i = 1,nconf2
      write (*,'(i9,2x,i9)') i,rots2(i)
    end do
    write (*,*)

!---- run some checks that are mandatory in order for the compariso to work properly
    if (nat1 /= nat2) then
      write (*,'(a)') "Nat1 /= Nat2 : Number of atoms of the two ensembles don't match!"
      write (*,'(a)') "You are trying to compare two different molecules!"
      error stop "exit."
    else
      nat = nat1
    end if
    do i = 1,nat
      if (at1(i) /= at2(i)) then
        write (*,'(a)') "The ordering of atoms apparently is different between the two ensembles!"
        write (*,'(a)') "This way it is impossible to calculate RMSDs!"
        error stop "exit."
      end if
    end do

!---- set the threads for the RMSD calculation (OMP parallel)
    if (env%autothreads) then
      call ompautoset(env%threads,4,env%omp,env%MAXRUN,0) !mode=4 --> Program intern Threads max
    end if

!---- printout
    call smallhead('Comparing the Ensembles')

!----- set up some more arrays
    allocate (rmat(nconf1,nconf2),c1(3,nat),c2(3,nat),en1(nconf1),en2(nconf2))
!---- get energies of lowest conformers
    do i = 1,nconf1
      en1(i) = eread1(b1(i))
    end do
    do i = 1,nconf2
      en2(i) = eread2(b2(i))
    end do
    deallocate (eread2,eread1)

    write (*,'(1x,a)',advance='no') 'Calculating RMSDs between conformers...'
!-----compare conformer-block-wise, i.e. all rotamers of each conformer pair of ensemble 1 and 2
!----- THIS IS THE MAIN LOOP OF THE SUBROUTINE
    do i = 1,nconf1
      be = b1(i)
      ed = e1(i)
      nrot1 = ed-be+1
      !write(*,*) b,e,nrot1
      allocate (conf1(3,nat,nrot1))
      conf1(1:3,1:nat,1:nrot1) = xyz1(1:3,1:nat,be:ed)
      do j = 1,nconf2
        be = b2(j)
        ed = e2(j)
        nrot2 = ed-be+1
        allocate (conf2(3,nat,nrot2))
        conf2(1:3,1:nat,1:nrot2) = xyz2(1:3,1:nat,be:ed)
        rcount = nrot1*nrot2 !total number of RMSDs between
        allocate (rmat2(rcount))
        rmat2 = 100 !<--only done so we can easily get the minimum rmsd later

        do k = 1,nrot1
          c1(1:3,1:nat) = conf1(1:3,1:nat,k)
!$omp parallel private (l,rcount,c2,ydum,xdum,Udum,gdum) &
!$omp shared (k,c1,rmat2,nat,conf2,nrot2)
!$omp do
          do l = 1,nrot2
            rcount = ((k-1)*nrot2)+l
            c2(1:3,1:nat) = conf2(1:3,1:nat,l)
            call rmsd(nat,c1,c2,0,Udum,xdum,ydum,rmat2(rcount),.false.,gdum) ! all atoms
            !write(*,*)rmat2(rcount)
          end do
!$omp end do
!$omp end parallel
        end do
        rmat(i,j) = minval(rmat2)
        deallocate (rmat2,conf2)
      end do
      deallocate (conf1)
    end do
    write (*,'(1x,a)') 'done.'
!---------
    !--- print the RMSD matrix
    write (*,'(1x,a,f8.4,a)') 'RMSD threshold:',rthr,' Å'
    call pr_rmat(rmat,nconf1,nconf2)
!--------
    write (*,*)
    call smallhead('Correlation between Conformers :')
    write (*,'(4x,a,5x,a10,13x,a,4x,a10)') '#','Ensemble A','#','Ensemble B'
    dum = nconf1+nconf2
    kk = nconf1
    ll = nconf2
    connect = '<---->'
    do
      matching = .false.
      if (kk .gt. 0.and.ll .gt. 0) then
        min_1 = en1(kk)
        min_2 = en2(ll)
        rval = rmat(kk,ll)
        if (rval .le. rthr) then  !<--- two matching conformers from ensembles
          write (*,'(1x,i4,1x,f14.5,3x,a,1x,i4,f14.5)') kk,min_1,connect,ll,min_2
          kk = kk-1
          ll = ll-1
          matching = .true.
        end if
      end if
      if (.not.matching) then
        if (kk .gt. 0.and.ll .gt. 0) then
          min_1 = en1(kk)
          min_2 = en2(ll)
          if (min_1 .ge. min_2) then
            write (*,'(1x,i4,1x,f14.5)') kk,min_1
            kk = kk-1
          else
            write (*,'(30x,i4,f14.5)') ll,min_2
            ll = ll-1
          end if
        end if
        if (kk .gt. 0.and.ll .eq. 0) then
          min_1 = en1(kk)
          write (*,'(1x,i4,1x,f14.5)') kk,min_1
          kk = kk-1
        end if
        if (kk .eq. 0.and.ll .gt. 0) then
          min_2 = en2(ll)
          write (*,'(30x,i4,f14.5)') ll,min_2
          ll = ll-1
        end if
      end if
      if (kk .eq. 0.and.ll .eq. 0) exit
    end do

!----- plotfiles
    !--- convert to relative energies
    min_1 = minval(en1)
    min_2 = minval(en2)
    min_tot = min(min_1,min_2)
    !en1=(en1-min_tot)*autokcal
    !en2=(en2-min_tot)*autokcal
    !--- energies of ensemble 1
    open (newunit=ich,file='energy_1.dat')
    do i = 1,nconf1
      write (ich,'(2x,f10.5,2x,f14.8)') (en1(i)-min_tot)*autokcal,en1(i)
    end do
    close (ich)
    !--- energies of ensemble 2
    open (newunit=ich,file='energy_2.dat')
    do i = 1,nconf2
      write (ich,'(2x,f10.5,2x,f14.8)') (en2(i)-min_tot)*autokcal,en2(i)
    end do
    close (ich)
    !--- correlation between conformers of the two ensembles
    open (newunit=ich2,file='rmsdmatch.dat')
    do i = 1,nconf1
      do j = 1,nconf2
        rval = rmat(i,j)
        if (rval .le. rthr) then
          write (ich2,'(2x,i6,i6)') i,j
        end if
      end do
    end do
    close (ich2)

!-----
    deallocate (rmat,c1,c2,en1,en2)
    deallocate (b1,e1,b2,e2)
    deallocate (xyz2,xyz1,at2,at1)
    deallocate (ydum,xdum,Udum,gdum)
!-----
    call compens_cleanup()

  end associate

end subroutine compare_ensembles

!---------------------------------------------------------------------------------------
subroutine compens_cleanup()
  use iso_fortran_env,only:wp => real64
  use iomod
  implicit none
  call remove('scoord.1')
  call remove('ensemble1.inp.xyz')
  call remove('ensemble2.inp.xyz')
  call remove('track.1')
  call remove('track.2')
  call remove('crest_best.xyz')
  call remove('cregen.out.tmp')
end subroutine compens_cleanup

!--------------------------------------------------------------------------------------
subroutine pr_rmat(rmat,acon,bcon)
  use iso_fortran_env,only:wp => real64
  implicit none
  integer :: acon,bcon
  real(wp) :: rmat(acon,bcon)
  integer :: i,j
  write (*,*)
  write (*,*) 'RMSD matrix:'
  write (*,'(2x,a9)',advance='no') 'conformer'
  do j = 1,bcon
    write (*,'(1x,i10)',advance='no') j
  end do
  write (*,*)
  do i = 1,acon
    write (*,'(1x,i5,5x)',advance='no') i
    do j = 1,bcon
      write (*,'(1x,f10.5)',advance='no') rmat(i,j)
    end do
    write (*,*)
  end do
end subroutine pr_rmat

