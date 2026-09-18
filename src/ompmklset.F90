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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!c OMP and MKL parallelization settings
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!

subroutine ompmklset(threads)
  use omp_lib
  implicit none
  integer,intent(in) :: threads

  call OMP_Set_Num_Threads(threads)
#ifdef WITH_MKL
  call MKL_Set_Num_Threads(threads)
  call mkl_set_dynamic(0)
#endif
! call openblasset(threads)
end subroutine ompmklset

subroutine openblasset(threads)
  implicit none
  integer,intent(in) :: threads
#ifdef WITH_OPENBLAS
  call openblas_set_num_threads(threads)
#endif
  return
end subroutine openblasset

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!c OMP and MKL parallelization settings (short routine)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!

subroutine ompenvset(omp)
  use iomod
  implicit none
  integer,intent(in) :: omp
  integer :: io

  io = setenv('OMP_NUM_THREADS',omp)
  io = setenv('MKL_NUM_THREADS',omp)
  io = setenv('OPENBLAS_NUM_THREADS',omp) 

end subroutine ompenvset

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!c OMP and MKL autoset switchcase routine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!

subroutine new_ompautoset(env,modus,maxjobs,parallel_jobs,cores_per_job)
!***********************************************************************
!* Determine the OMP thread split (parallel jobs x cores per job) for a
!* given work distribution mode.
!*
!* For the 'auto'/'auto_nested' modes the per-job core reservation Treq
!* is taken from env%calc%maxthreads() (the largest thread count any
!* active calculation level requests). Concurrent jobs are capped so
!* that parallel_jobs*Treq <= env%threads, guaranteeing every job can
!* host its heaviest level without oversubscribing the machine. With all
!* levels at the default (threads unset) Treq=1 and the historical
!* behavior is reproduced exactly.
!***********************************************************************
  use omp_lib
  use crest_data
  use crest_parameters,only:wp,stdout
  implicit none
  type(systemdata),intent(inout) :: env
  character(len=*),intent(in)    :: modus
  integer,intent(in)  :: maxjobs
  integer,intent(out) :: parallel_jobs
  integer,intent(out) :: cores_per_job
  integer :: T,Tdiff,Treq,Tcap,idle
  real(wp) :: Tfrac,Tfloor

  !> The default, all threads allocated to CREST
  T = env%threads
  parallel_jobs = max(1,T)
  cores_per_job = 1
  !> More settings, nested parallelization reset
  call omp_set_max_active_levels(1)

  !> per-job core reservation (heaviest active level), default 1
  Treq = env%calc%maxthreads()

  select case (modus)
  case ('auto','auto_nested')
    !> distribute jobs automatically, reserving Treq cores per job so that
    !> parallel_jobs*Treq <= T (no oversubscription). Treq is a hard cap:
    !> levels are never grown to soak leftover cores, so if T is not a
    !> multiple of Treq the remainder stays idle (warned about below).
    if (Treq > T) then
      !> a single level requests more cores than the whole budget
      write (stdout,'(1x,a,i0,a,i0,a)') &
        & '**WARNING** a calculation level requests ',Treq, &
        & ' cores but only ',T,' are available; running a single job'
      parallel_jobs = 1
      cores_per_job = T
    else
      !> jobs are core-bound at T/Treq, and additionally capped by the
      !> number of available jobs (maxjobs); maxjobs<=0 means "unbounded"
      Tdiff = T/Treq
      if (maxjobs > 0) Tdiff = min(maxjobs,Tdiff)
      parallel_jobs = max(1,Tdiff)
      cores_per_job = max(Treq,T/parallel_jobs)
    end if
    !> inform about idle cores: only hard-capped levels (ORCA %pal) leave
    !> cores unused; internal calcs and generic subprocesses grow into
    !> cores_per_job via OMP_NUM_THREADS / nested OpenMP and soak the rest
    Tcap = env%calc%maxthreads_capped()
    idle = T-parallel_jobs*Tcap
    if (Tcap > 1 .and. idle > 0) then
      write (stdout,'(1x,a,i0,a,i0,a,i0,a,i0,a)') &
        & '**NOTE** capped subprocess levels reserve ',Tcap, &
        & ' core(s) each; with ',parallel_jobs,' parallel job(s) on ',T, &
        & ' threads, ',idle,' core(s) stay idle during those evaluations'
    end if
    if (index(modus,'_nested') .ne. 0 .and. cores_per_job > 1) then
      if (env%omp_allow_nested) then
        !> We should never need more than two active nested layers
        call omp_set_max_active_levels(2)
      endif
    else
#ifdef WITH_MKL
      !call mkl_free_buffers()
      call mkl_set_dynamic(0)
#endif
    end if
    call openblasset(1)

  case ('max')
    !> Both intern and environment variable threads to max
    parallel_jobs = T
    cores_per_job = T

  case ('min','serial')
    !> Both intern and environment variable threads to one (like a serial program)
    parallel_jobs = 1
    cores_per_job = 1

  case ('subprocess','la-focus')
    !> CREST itself uses one thread, and but the environment variable is set to max
    !> which is useful when driving a single subprocess/systemcall
    parallel_jobs = 1
    cores_per_job = T 

    !> the setting may also be used for linear-algebra focused runs, in which case 
    !> nested parallelism should be active
    if (env%omp_allow_nested) then
      !> We should never need more than two active nested layers
      call omp_set_max_active_levels(2)
    endif

  end select

  !> apply the calculated settings
  call ompmklset(parallel_jobs)
  call ompenvset(cores_per_job)
#ifdef WITH_OPENBLAS
!  if(modus.ne.'auto'.and.modus.ne.'auto_nested')then
!      call openblasset(cores_per_job)
!  endif
   call openblasset(1)
#endif
end subroutine new_ompautoset

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!c report the parallelization split applied by new_ompautoset
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!

subroutine ompautoset_summary(env,label,parallel_jobs,cores_per_job)
!***********************************************************************
!* Print a one-line summary of the parallelization split that
!* new_ompautoset just applied: how many jobs run concurrently and
!* how many cores each job may use, plus whether a subprocess level is
!* hard-capped (ORCA %pal) or nested OpenMP is active. Bare subroutine
!* with no optional args -- callers that want silence just don't call it.
!*
!*  env           - system data (threads budget + calc levels)
!*  label         - short task name shown in the line (e.g. 'optimizations')
!*  parallel_jobs - number of concurrently running jobs (T)
!*  cores_per_job - cores reserved per job (Tn)
!***********************************************************************
  use crest_data
  use crest_parameters,only:stdout
  implicit none
  type(systemdata),intent(in) :: env
  character(len=*),intent(in) :: label
  integer,intent(in) :: parallel_jobs,cores_per_job
  integer :: Tcap
  character(len=:),allocatable :: jobs_word,cores_word,extra

  ! ── singular/plural wording ──────────────────────────────────
  jobs_word = ' parallel job'
  if (parallel_jobs > 1) jobs_word = ' parallel jobs'
  cores_word = ' core/job'
  if (cores_per_job > 1) cores_word = ' cores/job'

  ! ── flag hard-capped subprocesses vs. nested OpenMP ──────────
  Tcap = env%calc%maxthreads_capped()
  extra = ''
  if (Tcap > 1) then
    extra = ', subprocess hard-capped'
  else if (env%omp_allow_nested .and. cores_per_job > 1) then
    extra = ', nested OpenMP'
  end if

  write (stdout,'(1x,"↳ ",a,": ",i0,a," × ",i0,a,"  (",i0," threads",a,")")') &
    & trim(label),parallel_jobs,trim(jobs_word),cores_per_job,trim(cores_word), &
    & env%threads,trim(extra)

end subroutine ompautoset_summary

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!c get omp/mkl automatically from the global variables
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!

subroutine ompgetauto(threads,omp,maxrun)
  use omp_lib
  use iomod
  implicit none
  integer,intent(inout) :: threads,omp,maxrun
  integer :: nproc
  integer :: r
  character(len=256) :: val

  call getenv('OMP_NUM_THREADS',val)
  read (val,*,iostat=r) nproc
  if (r .ne. 0) then
    nproc = 1
  end if
  threads = nproc
  maxrun = 1
  omp = nproc

end subroutine ompgetauto

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!c print omp/mkl threads that are used at the moment
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
subroutine ompprint_intern()
  use omp_lib
  implicit none
  integer :: nproc,TID

!$OMP PARALLEL PRIVATE(TID)
  TID = OMP_GET_THREAD_NUM()
  IF (TID .EQ. 0) THEN
    nproc = OMP_GET_NUM_THREADS()
    write (*,*) '============================='
    write (*,*) ' # threads =',nproc
    write (*,*) '============================='
  END IF
!$OMP END PARALLEL
end subroutine ompprint_intern

