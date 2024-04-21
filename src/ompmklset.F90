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
#endif
end subroutine ompmklset

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

end subroutine ompenvset

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!c OMP and MKL autoset switchcase routine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!

subroutine new_ompautoset(env,modus,maxjobs,parallel_jobs,cores_per_job)
  use omp_lib
  use crest_data
  use crest_parameters,only:wp
  implicit none
  type(systemdata),intent(inout) :: env
  character(len=*),intent(in)    :: modus
  integer,intent(in)  :: maxjobs
  integer,intent(out) :: parallel_jobs
  integer,intent(out) :: cores_per_job
  integer :: T,Tdiff
  real(wp) :: Tfrac,Tfloor

  !> The default, all threads allocated to CREST
  T = env%threads
  parallel_jobs = T
  cores_per_job = 1
  !> More settings, nested parallelization reset
  call omp_set_max_active_levels(1)
  !call omp_set_dynamic(.true.)

  select case (modus)
  case ('auto','auto_nested')
    !> distribute jobs automatically:
    !> if more cores are available than maxjobs, try to distribute remaining
    !> threads EVENLY for each job
    if (maxjobs > 0.and.T > maxjobs) then
      parallel_jobs = maxjobs
      Tfrac = real(T)/real(maxjobs)
      Tfloor = floor(Tfrac)
      cores_per_job = max(nint(Tfloor),1)
    end if
    if (index(modus,'_nested') .ne. 0 .and. cores_per_job > 1) then
      if (env%omp_allow_nested) then
        !> We should never need more than two active nested layers
        call omp_set_max_active_levels(2)
      end if
    end if

  case ('max')
    !> Both intern and environment variable threads to max
    parallel_jobs = T
    cores_per_job = T

  case ('min','serial')
    !> Both intern and environment variable threads to one (like a serial program)
    parallel_jobs = 1
    cores_per_job = 1

  case ('subprocess')
    !> CREST itself uses one thread, and but the environment variable is set to max
    !> which is useful when driving a single subprocess/systemcall
    parallel_jobs = 1
    cores_per_job = T

  end select

  !> apply the calculated settings
  call ompmklset(parallel_jobs)
  call ompenvset(cores_per_job)

end subroutine new_ompautoset

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

