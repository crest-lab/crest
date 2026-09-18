!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2024 Philipp Pracht, David Wales
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

module bh_mc_module
  use crest_parameters
  use iomod
  use strucrd,only:coord
  use crest_calculator
  use optimize_module
  use axis_module
  use irmsd_module
  use canonical_mod
  use quicksort_interface,only:ensemble_qsort
  use bh_class_module
  use bh_step_module
  implicit none
  private

!  logical,parameter :: debug = .true.
  logical,parameter :: debug = .false.

  public :: mc

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine mc(calc,mol,bh,iostat,verbosity)
!********************************************************************
!* A thread-safe single basin-hopping MC run
!* Parameters and quenched structures are saved within the bh object
!********************************************************************
    implicit none
    !> IN/OUTPUT
    type(calcdata),intent(inout) :: calc  !> potential settings
    type(coord),intent(inout)    :: mol   !> molecular system
    type(bh_class),intent(inout) :: bh    !> BH settings
    integer,intent(out) :: iostat
    integer,intent(in),optional  :: verbosity  !> printout parameter
    !> LOCAL
    type(coord) :: tmpmol    !> copy to take steps
    type(coord) :: optmol    !> quenched structure
    integer :: iter,iostatus,accepted,discarded,broke
    real(wp) :: etot,ratio
    real(wp),allocatable :: grd(:,:)
    logical :: accept,dupe,broken
    integer :: printlvl,first,last,dynamicseed
    character(len=20) :: tag

    !$omp critical
    write (tag,'("BH[Runner ",i3,"]>")') bh%id
    !$omp end critical

    if (present(verbosity)) then
      printlvl = verbosity
    else
      printlvl = 0
    end if

    iostat = 0

!>--- Add input energy to Markov chain after an initial quench
    !$omp critical
    allocate (grd(3,mol%nat),source=0.0_wp)
    !$omp end critical

    if (printlvl > 0) then
      !$omp critical
      write (stdout,'(a,1x,a)') trim(tag),'Performing '//colorify('initial quench','gold')//"."
      !$omp end critical
    end if

    tmpmol = mol
    call mcquench(calc,bh,tmpmol,optmol,etot,grd,iostatus)
    if (iostatus .ne. 0) then
      !$omp critical
      write (stdout,'(a,1x,a)') trim(tag),colorify('** WARNING **','red')// &
        & ' initial quench failed. Returning.'
      !$omp end critical
      iostat = iostatus
      return
    end if
    !$omp critical
    mol = optmol
    bh%emin = mol%energy
    !$omp end critical
    call bh%add(mol)

!>--- print information about the run?
    if (printlvl > 0) then
      !$omp critical
      call mcheader(bh)
      !$omp end critical
    end if

!>--- seed the RNG?
    if (allocated(bh%seed)) then
      dynamicseed = bh%seed+(bh%iteration-1)+bh%id*1000
      if (printlvl > 1) then
        !$omp critical
        write (stdout,'(a,1x,2(a,i0),a)') trim(tag), &
        & 'Seeding current RNG instance with: ',bh%seed,' (',dynamicseed,')'
        !$omp end critical
      end if
      call RNG_seed(bh%seed)
    end if

!=======================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
    accepted = 0
    discarded = 0
    broke = 0
    MonteCarlo: do iter = 1,bh%maxsteps
      broken = .false.
      dupe = .false.

!>--- Take the step (mol --> tmpmol)
      call takestep(mol,calc,bh,tmpmol)

!>--- Quench it (tmpmol --> optmol)
      call mcquench(calc,bh,tmpmol,optmol,etot,grd,iostatus)

!>--- Accept/reject
      if (iostatus == 0) then  !> successfull optimization

        if (printlvl > 1) then
          !$omp critical
          write (stdout,'(a,1x,a,es17.8,a,es17.8,a)') trim(tag),'Quench E=',etot, &
          &  ' Eh, Markov E=',bh%emin,' Eh'
          !$omp end critical
        end if

        accept = mcaccept(optmol,bh)
        if (accept) then
          accepted = accepted+1

          call axis(optmol%nat,optmol%at,optmol%xyz)

          !> check duplicates here
          call mcduplicate(mol,bh,dupe,broken)

          !$omp critical
          if (printlvl > 1) then
            write (stdout,'(a)',advance='no') repeat(' ',len_trim(tag)+1)// &
            & "Quench "//colorify('ACCEPTED','green')
          end if

          if (broken) then
            broke = broke+1
            if (printlvl > 1) write (stdout,'(a)',advance='no') &
            & ', but REJECTED due to topology mismatch!'
          else if (dupe) then
            discarded = discarded+1
            if (printlvl > 1) write (stdout,'(a)',advance='no') &
            & ', but '//colorify('NOT SAVED','yellow')//' due to duplicate detection!'
          else if (printlvl == 1) then
            write (stdout,'(a,1x,a,a,es17.8,a)') trim(tag),"Quench "//colorify('ACCEPTED','green'), &
              & ', NEW Markov E=',bh%emin,' Eh'
          end if

          if (printlvl > 1) write (stdout,'(/)')
          !$omp end critical
        else
          if (printlvl > 1) then
            !$omp critical
            write (stdout,'(a,a,/)') repeat(' ',len_trim(tag)+1), &
            &                 'Quench '//colorify('REJECTED','red')//', does not fulfill MC criterion'
            !$omp end critical
          end if
          cycle MonteCarlo
        end if
      else
        if (printlvl > 1)then
          !$omp critical
          write (stdout,'(a,1x,a,/)') trim(tag),"Quench "//colorify("FAILED","red")
          !$omp end critical
        endif
        cycle MonteCarlo
      end if

!>--- Update structures
      if (.not.broken) then
        !> continue Markov chain
        !$omp critical
        mol = optmol
        !$omp end critical

        if (.not.dupe) then
          !> Save new unique structures
          call bh%add(mol)
        end if
      end if

    end do MonteCarlo
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!=======================================================================================!

!>--- Post-processing
    first = 1
    last = bh%saved
    call ensemble_qsort(bh%maxsave,bh%structures,first,last)

!>--- Stats
    if (printlvl > 0) then
      !$omp critical
      call mcstats(bh,accepted,discarded,broke)
      !$omp end critical
    end if

    deallocate (grd)
  end subroutine mc

!=========================================================================================!

  subroutine mcheader(bh)
    implicit none
    type(bh_class),intent(in) :: bh
    character(len=80) :: atmp
    integer :: n

    write (stdout,'(t8,a)') '┌'//repeat('─',63)//'┐'
    write (stdout,'(t8,a,1x)',advance='no') '│'
    write (stdout,'(a,3x)',advance='no') 'Starting Basin-Hopping Global Optimization'
    write (stdout,'(a,i3,a)',advance='no') '[Thread/ID ',bh%id,']'
    write (stdout,'(2x,"│")')
    write (stdout,'(t8,a)') '╞'//repeat('═',63)//'╡'

    write (stdout,'(t8,a,1x)',advance='no') '│'
    write (stdout,'(a,f20.10,a)',advance='no') 'Initial energy:',bh%emin,' Eh'
    write (stdout,'(24x,"│")')

    write (stdout,'(t8,a,1x)',advance='no') '│'
    write (stdout,'(a,es10.3,2x)',advance='no') 'T/K: ',bh%temp
    write (stdout,'(a,i5,3x)',advance='no') 'steps: ',bh%maxsteps
    write (stdout,'(a,i5,3x)',advance='no') 'max save: ',bh%maxsave
    write (stdout,'(12x,"│")')

    if (allocated(bh%seed)) then
      write (stdout,'(t8,a,1x)',advance='no') '│'
      write (atmp,'(a,i0)') 'Random number generation (reference) seed: ',bh%seed
      write (stdout,'(a,1x)',advance='no') trim(atmp)
      n = 61-len_trim(atmp)
      write (stdout,'(a)',advance='no') repeat(' ',n)
      write (stdout,'("│")')
    end if

    write (stdout,'(t8,a,1x)',advance='no') '│'
    write (stdout,'(a,a,2x)',advance='no') 'step type: ',steptypestr(bh%steptype)
    write (stdout,'(a,3f9.5)',advance='no') 'step size:',bh%stepsize(1:3)
    write (stdout,'(3x,"│")')

    write (stdout,'(t8,a,1x)',advance='no') '│'
    write (stdout,'(a,f9.5,a)',advance='no') 'Thresholds   ΔRMSD:',bh%rthr,' Å,  '
    write (stdout,'(a,es11.4,a)',advance='no') 'ΔE: ',bh%ethr,' kcal/mol'
    write (stdout,'(5x,"│")')

    write (stdout,'(t8,a)') '└'//repeat('─',63)//'┘'
  end subroutine mcheader

  subroutine mcstats(bh,accepted,discarded,broke)
    implicit none
    type(bh_class),intent(in) :: bh
    integer,intent(in) :: accepted,discarded,broke
    real(wp) :: ratio

    write (stdout,'(t8,a)') colorify('┌'//repeat('─',63)//'┐','green')
    write (stdout,'(t8,a,1x)',advance='no') colorify("│","green")
    write (stdout,'(a,21x)',advance='no') 'Basin-Hopping Statistics'
    write (stdout,'(a,i3,a)',advance='no') '[Thread/ID ',bh%id,']'
    write (stdout,'(2x,a)') colorify("│","green")

    write (stdout,'(t8,a,1x)',advance='no') colorify("│","green")
    write (stdout,'(a,12x,f17.8,a)',advance='no') 'Latest Markov chain energy: ',bh%emin,' Eh'
    write (stdout,'(2x,a)') colorify("│","green")

    write (stdout,'(t8,a,1x)',advance='no') colorify("│","green")
    ratio = real(accepted,wp)/real(bh%maxsteps,wp)
    write (stdout,'(a,f6.2,a)',advance='no') 'MC acceptance ratio ',ratio*100.0_wp,' %, '
    ratio = real(discarded,wp)/real(accepted,wp)
    write (stdout,'(a,f6.2,a)',advance='no') 'similarity rejection  ',ratio*100.0_wp,' %'
    write (stdout,'(2x,a)') colorify("│","green")

    write (stdout,'(t8,a,1x)',advance='no') colorify("│","green")
    ratio = real(broke,wp)/real(accepted,wp)
    write (stdout,'(a,f6.2,a)',advance='no') 'topology rejection  ',ratio*100.0_wp,' %, '
    ratio = real(accepted-discarded-broke,wp)/real(bh%maxsteps,wp)
    write (stdout,'(a,f6.2,a)',advance='no') 'TOTAL ACCEPT ratio    ',ratio*100.0_wp,' %'
    write (stdout,'(2x,a)') colorify("│","green")
    write (stdout,'(t8,a)') colorify('└'//repeat('─',63)//'┘',"green")
  end subroutine mcstats

!=========================================================================================!

  subroutine RNG_seed(iseed)
!*************************************
!* seed the RNG to get a reproducible
!* sequence of random numbers
!*************************************
    integer,intent(in) :: iseed
    integer :: n
    integer,allocatable :: seedArray(:)
    !> 1) Query how many integers are needed to set the seed (compiler dependent!)
    call random_seed(size=n)
    !> 2) Allocate and assign a known pattern
    allocate (seedArray(n))
    seedArray(:) = iseed
    !> 3) Set the seed explicitly
    call random_seed(put=seedArray)
    deallocate (seedArray)
  end subroutine RNG_seed

!=========================================================================================!

  function mcaccept(mol,bh) result(accept)
!**************************************
!* The regular MC acceptance condition
!**************************************
    implicit none
    logical :: accept
    type(coord),intent(in) :: mol
    type(bh_class),intent(in) :: bh
    real(wp) :: eold,enew,temp
    real(wp) :: random,fact
    accept = .false.
    eold = bh%emin
    enew = mol%energy
    temp = bh%temp*kB  !> Kelvin to a.u.

    if (enew .lt. eold) then
      accept = .true.
    else
      call random_number(random)
      fact = exp(-(enew-eold)/temp)
      if (fact .gt. random) accept = .true.
    end if

  end function mcaccept

!=========================================================================================!

  subroutine mcduplicate(mol,bh,dupe,broken)
!*****************************************************
!* Check if a new structure (mol) is already in the
!* list of saved structures (bh%structures)
!*****************************************************
    implicit none
    type(coord),intent(in)       :: mol
    type(bh_class),intent(inout) :: bh
    real(wp) :: rthr,ethr
    logical,intent(out) :: dupe,broken
    !> LOCAL
    integer :: i,j,k,l,nat
    type(canonical_sorter) :: newsort
    real(wp) :: rmsdval,deltaE
    logical :: topocheck

    dupe = .false.
    broken = .false.
    ethr = bh%ethr
    rthr = bh%rthr
    nat = mol%nat
    topocheck = .true.

    if (debug) write (*,*)

    if (.not.allocated(bh%rcache)) then
      if (debug) write (*,*) "allocating RCACHE"
      !$omp critical
      allocate (bh%rcache)
      call bh%rcache%allocate(nat)
      !$omp end critical
    end if

    !$omp critical
    call newsort%init(mol,invtype='apsp+',heavy=.false.)
    !$omp end critical

    COMPARE: do i = 1,bh%saved

      !> Energy difference
      deltaE = (mol%energy-bh%structures(i)%energy)*autokcal

      !> Geometry difference (permutation-invariant RMSD)
      if (topocheck) then
        bh%rcache%rank(1:nat,1) = newsort%rank(1:nat)
        bh%rcache%rank(1:nat,2) = bh%sorters(i)%rank(1:nat)
      end if
      call min_rmsd(mol,bh%structures(i), &
      &        rcache=bh%rcache,rmsdout=rmsdval)

      if (debug) write (*,'(a,es15.4,a,es15.4,a)') 'RMSD=',rmsdval*autoaa, &
      &          ' Å, delta E=',deltaE,' kcal/mol'

      !> Check
      if (abs(deltaE) .lt. ethr.and.rmsdval*autoaa .lt. rthr) then
        dupe = .true.
        if (deltaE < 0.0_wp) then
          !> if the energy is lower, we replace the molecule (better conformation)
          bh%structures(i) = mol
        end if
        exit COMPARE
      end if
    end do COMPARE

    !$omp critical
    call newsort%deallocate()
    !$omp end critical
  end subroutine mcduplicate

  !========================================================================================!

  subroutine mcquench(calc,bh,tmpmol,optmol,etot,grd,iostat)
    implicit none
    !> Input
    type(calcdata),intent(inout) :: calc      !> potential settings
    type(bh_class),intent(inout) :: bh        !> BH settings
    type(coord),intent(inout)    :: tmpmol    !> molecular system
    real(wp),intent(inout)       :: etot      !> quenechd energy
    real(wp),intent(inout)       :: grd(:,:)  !> gradient (temp storage)
    !> Output
    type(coord),intent(out)      :: optmol    !> molecular system output
    integer,intent(out)          :: iostat    !> return status

    integer :: nrefine,ii
    real(wp) :: etmp
    iostat = 1

    !> initial proper quench (refine_lvl = 0)
    call optimize_geometry(tmpmol,optmol,calc,etot,grd, &
     &                      .false.,.false.,iostat)

    !> Special Post-processing via refinement queue
    if (allocated(bh%refine_queue).and.iostat == 0) then

      nrefine = size(bh%refine_queue,1)

      do ii = 1,nrefine
        if (iostat .ne. 0) exit
        calc%refine_stage = bh%refine_queue(ii)
        select case (calc%refine_stage)
        case (1) !> singlepoint (rerank)
          call engrad(optmol,calc,etot,grd,iostat)

        case (2) !> singlepoint (add)
          call engrad(optmol,calc,etmp,grd,iostat)
          etot = etot+etmp

        case (3) !> geometry opt (requench)
          tmpmol = optmol
          call optimize_geometry(tmpmol,optmol,calc,etot,grd, &
          &                      .false.,.false.,iostat)

        case default
          continue
        end select
      end do

      !> RESET refine level for next quench
      calc%refine_stage = 0

      !> Important: last energy must be stored in the optmol
      optmol%energy = etot
    end if

  end subroutine mcquench

!=========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!=========================================================================================!
end module bh_mc_module
