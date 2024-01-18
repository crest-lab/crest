!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Philipp Pracht
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
!--------------------------------------------------------------------------------!
!> The original (unmodified) source code can be found under the GNU LGPL 3.0 license
!> Copyright (C) 2019-2020 Sebastian Ehlert, Stefan Grimme
!> at https://github.com/grimme-lab/xtb
!================================================================================!
module crest_type_timer
  use iso_fortran_env,only:wp => real64,int64
  implicit none

  public :: timer
  private

!=========================================================!
  type :: timer

    !> number of timers
    integer,private :: n = 0

    !> printlevel
    logical,private :: verbose = .false.

    real(wp),private :: totwall = 0.0_wp
    real(wp),private :: totcpu = 0.0_wp
    logical,private,allocatable :: running(:)
    real(wp),private,allocatable :: twall(:)
    real(wp),private,allocatable :: tcpu(:)
    character(len=128),private,allocatable :: tag(:)
    integer :: ltag !> max tag length

  contains

    procedure :: new => allocate_timer
    procedure :: allocate => allocate_timer
    procedure :: deallocate => deallocate_timer
    procedure :: measure => timer_measure
    procedure :: write_timing
    procedure :: write => write_all_timings
    procedure :: get => get_timer
    procedure,private :: start_timing
    procedure,private :: stop_timing

    procedure :: init => allocate_timer
    procedure :: clear => deallocate_timer
    procedure :: start => timer_measure
    procedure :: stop => timer_measure

  end type timer
!========================================================!

  integer,parameter,private :: lmax_default = 26

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

!> To initialize timer
  subroutine allocate_timer(self,n,verbose)

    implicit none

    !> instance of timer
    class(timer),intent(inout) :: self

    !> number of timers
    integer,intent(in)           :: n

    !> if verbose
    logical,intent(in),optional :: verbose

    real(wp) :: time_cpu
    real(wp) :: time_wall

    call self%deallocate

    ! capture negative values !
    if (n < 1) return

    self%n = n
    if (present(verbose)) self%verbose = verbose
    allocate (self%twall(0:n),source=0.0_wp)
    allocate (self%tcpu(0:n),source=0.0_wp)
    allocate (self%running(n),source=.false.)
    allocate (self%tag(n)); self%tag = ' '
    self%ltag = 0

    ! launch timer !
    call self%start_timing(0)

  end subroutine allocate_timer

!========================================================================================!
  subroutine deallocate_timer(self)
!***********************
!* To deallocate memory
!***********************
    implicit none

    !> instance of timer
    class(timer),intent(inout) :: self

    self%n = 0
    self%totwall = 0
    self%totcpu = 0
    self%verbose = .false.
    if (allocated(self%twall)) deallocate (self%twall)
    if (allocated(self%tcpu)) deallocate (self%tcpu)
    if (allocated(self%running)) deallocate (self%running)

  end subroutine deallocate_timer

!========================================================================================!
  function get_timer(self,i) result(time)
!********************************************
!* To obtain current elapsed time in seconds
!********************************************
    !> instance of timer
    class(timer),intent(inout) :: self

    !> if specific timer
    integer,intent(in),optional :: i

    integer  :: it
    real(wp) :: tcpu,twall
    real(wp) :: time
    logical  :: running

    ! if i is not given, calculate overall elapsed time !
    if (present(i)) then
      it = i
    else
      it = 0
    end if

    if (it > 0) then
      running = self%running(it)
    else
      running = .true.
    end if

    if (running) then
      call timing(tcpu,twall)
      time = self%twall(it)+twall
    else
      time = self%twall(it)
    end if

  end function get_timer

!========================================================================================!
  subroutine write_timing(self,iunit,i,inmsg,verbose)
!*****************************************
!* Printout function for a specific timer
!*****************************************
    implicit none

    !> instance of timer
    class(timer),intent(inout) :: self

    !> I/O unit
    integer,intent(in) :: iunit

    !> index
    integer,intent(in) :: i

    !> raw message text
    character(len=*),intent(in),optional :: inmsg

    !> if verbose
    logical,intent(in),optional :: verbose

    character(len=128) :: msg
    real(wp) :: cputime,walltime
    integer(int64) ::  cpudays,cpuhours,cpumins
    integer(int64) :: walldays,wallhours,wallmins
    logical :: lverbose
    integer :: lmax

    ! check if tag should be added !
    if (present(inmsg)) then
      msg = inmsg
      lmax = len_trim(inmsg)
    else
      msg = ''
      msg = self%tag(i)
      lmax = self%ltag
    end if
    lmax = max(lmax,lmax_default)

    ! verbosity settings !
    if (present(verbose)) then
      lverbose = verbose
    else
      lverbose = self%verbose
    end if

    !           DAYS   HOURS   MINUTES   SECONDS
    ! DAYS        1     1/24    1/1440   1/86400
    ! HOURS      24      1       1/60     1/3600
    ! MINUTES   1440    60        1        1/60
    ! SECONDS  86400   3600      60         1

    ! convert elapsed CPU time into days, hours, minutes !
    cputime = self%tcpu(i)
    cpudays = int(cputime/86400._wp)
    cputime = cputime-cpudays*86400._wp
    cpuhours = int(cputime/3600._wp)
    cputime = cputime-cpuhours*3600._wp
    cpumins = int(cputime/60._wp)
    cputime = cputime-cpumins*60._wp

    ! convert elapsed wall time into days, hours, minutes !
    walltime = self%twall(i)
    walldays = int(walltime/86400._wp)
    walltime = walltime-walldays*86400._wp
    wallhours = int(walltime/3600._wp)
    walltime = walltime-wallhours*3600._wp
    wallmins = int(walltime/60._wp)
    walltime = walltime-wallmins*60._wp

    !----------!
    ! printout !
    !----------!
    if (lverbose) then
      write (iunit,'(a)') msg(1:lmax)
      write (iunit,'(" * wall-time: ",i5," d, ",i2," h, ",i2," min, ",f6.3," sec")') &
        walldays,wallhours,wallmins,walltime
      write (iunit,'(" *  cpu-time: ",i5," d, ",i2," h, ",i2," min, ",f6.3," sec")') &
        cpudays,cpuhours,cpumins,cputime
      write (iunit,'(1x,"*",1x,"ratio c/w:",1x,f9.3,1x,"speedup")') self%tcpu(i)/self%twall(i)
    else
      write (iunit,'(a,1x,"...",i9," min, ",f6.3," sec")') &
        msg(1:lmax),wallmins,walltime
    end if

  end subroutine write_timing

!========================================================================================!
  subroutine write_all_timings(self,iunit,inmsg,verbose)
!************************************
!* Printout of all saved timings
!************************************
    implicit none

    !> instance of timer
    class(timer),intent(inout) :: self

    !> I/O unit
    integer,intent(in) :: iunit

    !> raw message
    character(len=*),intent(in),optional :: inmsg

    !> verbose?
    logical,intent(in),optional :: verbose

    character(len=128) :: msg
    character(len=256) :: atmp
    real(wp) :: cputime,walltime,partsum,iowall
    real(wp) :: cpusecs,wallsecs,wallsecs_abs
    integer  :: i
    integer(int64) ::  cpudays,cpuhours,cpumins
    integer(int64) :: walldays,wallhours,wallmins
    integer(int64) :: walldays_abs,wallhours_abs,wallmins_abs
    integer :: lmax,barlen,ll
    logical :: verbose_local

    call self%stop_timing(0)

    ! verbose ?
    if (present(verbose)) then
      verbose_local = verbose
    else
      verbose_local = self%verbose
    end if

    ! check if an external message should be added !
    if (present(inmsg)) then
      msg = inmsg//" (total)"
    else
      msg = "total time"
    end if
    lmax = len_trim(msg)
    lmax = max(lmax,self%ltag)
    lmax = max(lmax,lmax_default)
    do i = 1,self%n
      ll = len_trim(self%tag(i))
      if (ll > lmax) lmax = ll
    end do

    !           DAYS   HOURS   MINUTES   SECONDS
    ! DAYS        1     1/24    1/1440   1/86400
    ! HOURS      24      1       1/60     1/3600
    ! MINUTES   1440    60        1        1/60
    ! SECONDS  86400   3600      60         1

    ! convert overall elapsed CPU time into days, hours, minutes !
    cputime = self%tcpu(0)
    cpudays = int(cputime/86400._wp)
    cputime = cputime-cpudays*86400._wp
    cpuhours = int(cputime/3600._wp)
    cputime = cputime-cpuhours*3600._wp
    cpumins = int(cputime/60._wp)
    cputime = cputime-cpumins*60._wp
    cpusecs = cputime

    ! convert overall elapsed wall time into days, hours, minutes !
    walltime = self%twall(0)
    walldays = int(walltime/86400._wp)
    walldays_abs = walldays
    walltime = walltime-walldays*86400._wp
    wallhours = int(walltime/3600._wp)
    wallhours_abs = wallhours
    walltime = walltime-wallhours*3600._wp
    wallmins = int(walltime/60._wp)
    wallmins_abs = wallmins
    walltime = walltime-wallmins*60._wp
    wallsecs = walltime
    wallsecs_abs = wallsecs

    !----------!
    ! printout !
    !----------!
    write (atmp,'(1x,a,i5," d, ",i2," h, ",i2," min, ",f6.3," sec")') &
      msg(1:lmax+6),walldays,wallhours,wallmins,walltime
    write (iunit,'(a)') trim(atmp)
    if (verbose_local) then
      barlen = len_trim(atmp)
      write (iunit,'(1x,a)') repeat('-',barlen)
    end if

    ! printout every timer and corresponding speedup !
    partsum = 0.0_wp
    do i = 1,self%n
      walltime = self%twall(i)
      if (walltime <= 0.0_wp) cycle
      partsum = partsum+walltime
      wallmins = int(walltime/60._wp)
      walltime = walltime-wallmins*60._wp
      msg = self%tag(i)
      write (iunit,'(1x,a,1x,"...",i9," min, ",f6.3," sec (",f7.3,"%)")') &
        msg(1:lmax),wallmins,walltime,100*self%twall(i)/self%twall(0)
    end do
    !> everything else is I/O or setup (e.g. loading parametrizations)
    walltime = self%twall(0)-partsum
    iowall = walltime
    wallmins = int(walltime/60._wp)
    walltime = walltime-wallmins*60._wp
    msg = 'I/O and setup'
    write (iunit,'(1x,a,1x,"...",i9," min, ",f6.3," sec (",f7.3,"%)")') &
      msg(1:lmax),wallmins,walltime,100*iowall/self%twall(0)
    !> and finally, again the total cpu and wall time
    if (verbose_local) then
      write (iunit,'(1x,a)') repeat('-',barlen)

      write (iunit,'(" * wall-time: ",i5," d, ",i2," h, ",i2," min, ",f6.3," sec")') &
        walldays_abs,wallhours_abs,wallmins_abs,wallsecs_abs
      write (iunit,'(" *  cpu-time: ",i5," d, ",i2," h, ",i2," min, ",f6.3," sec")') &
        cpudays,cpuhours,cpumins,cpusecs
      write (iunit,'(1x,"*",1x,"ratio c/w:",1x,f9.3,1x,"speedup")') self%tcpu(0)/self%twall(0)
      write (iunit,'(1x,a)') repeat('-',barlen)
    end if

  end subroutine write_all_timings

!========================================================================================!
  subroutine timer_measure(self,i,inmsg)
!***********************************
!* Automatic start/stop of timer i
!***********************************
    implicit none

    !> instance of timer
    class(timer),intent(inout) :: self

    !> index
    integer,intent(in) :: i

    !> raw message text
    character(len=*),intent(in),optional :: inmsg
    integer :: l

    ! check if appropriate index is given !
    if (i > self%n.or.i < 1) return

    ! switcher between start/stop status !
    if (self%running(i)) then
      call self%stop_timing(i)
    else
      call self%start_timing(i)
    end if

    ! update status !
    self%running(i) = .not.self%running(i)

    ! assign tag to specific timer !
    if (present(inmsg)) then
      self%tag(i) = trim(inmsg)
      l = len_trim(inmsg)
    end if

  end subroutine timer_measure

!========================================================================================!
  subroutine start_timing(self,i)
!********************
!* Start timer i
!********************
    implicit none

    !> instance of timer
    class(timer),intent(inout) :: self

    !> index
    integer,intent(in) :: i

    real(wp) :: time_cpu
    real(wp) :: time_wall

    call timing(time_cpu,time_wall)
    self%tcpu(i) = self%tcpu(i)-time_cpu
    self%twall(i) = self%twall(i)-time_wall

  end subroutine start_timing

!========================================================================================!
  subroutine stop_timing(self,i)
!*******************
!* Stop timer i
!*******************
    implicit none

    !> instance of timer
    class(timer),intent(inout) :: self

    !> index
    integer,intent(in) :: i

    real(wp) :: time_cpu
    real(wp) :: time_wall

    call timing(time_cpu,time_wall)
    self%tcpu(i) = self%tcpu(i)+time_cpu
    self%twall(i) = self%twall(i)+time_wall

  end subroutine stop_timing

!========================================================================================!
  subroutine timing(time_cpu,time_wall)
!********************************************
!* To retrieve the current CPU and wall time
!********************************************
    implicit none

    real(wp),intent(out) :: time_cpu
    real(wp),intent(out) :: time_wall

    !> current value of system clock (time passed from arbitary point)
    integer(int64) :: time_count

    !> number of clock ticks per second (conversion factor b/n ticks and seconds)
    integer(int64) :: time_rate
    integer(int64) :: time_max

    call system_clock(time_count,time_rate,time_max)
    call cpu_time(time_cpu)

    ! elapsed time in seconds !
    time_wall = real(time_count,wp)/real(time_rate,wp)

  end subroutine timing

!========================================================================================!
!========================================================================================!
end module crest_type_timer
