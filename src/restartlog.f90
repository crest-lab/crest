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
!================================================================================!

!> Lightweight restart checkpoint for conformational search runtypes.
!> Records only which stage completed and which file was last written —
!> no ensemble data is stored.

module crest_restartlog
  use crest_parameters,only:wp,stdout
  implicit none
  private

  character(len=*),parameter,public :: restart_file = 'crest.restart'

  !> All state needed to resume a conformational search.
  type,public :: restart_data
    integer  :: version   = 1
    integer  :: runtype   = 0    !> crestver (crest_imtd=2 or crest_imtd2=22)
    integer  :: main_iter = 0    !> env%nreset at checkpoint
    integer  :: mtd_iter  = 0    !> last completed MTD iteration index
    integer  :: nmetadyn  = 0    !> env%nmetadyn (trimmed after first MTD pass)
    character(len=64)  :: stage        = 'none'
    character(len=512) :: last_file    = ''   !> last CREGEN-sorted file written
    real(wp) :: elowest   = 0.0_wp
    real(wp) :: eprivious = 0.0_wp
  end type restart_data

  public :: write_restart_log
  public :: read_restart_log
  public :: restart_file_exists
  public :: print_restart_info
  public :: delete_restart_log

!========================================================================================!
!========================================================================================!
contains !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  logical function restart_file_exists()
!*************************************
!* Returns .true. if crest.restart
!* exists in the current directory.
!*************************************
    implicit none
    inquire(file=restart_file,exist=restart_file_exists)
  end function restart_file_exists

!========================================================================================!

  subroutine write_restart_log(runtype,stage,main_iter,mtd_iter,nmetadyn, &
    &                          elowest,eprivious,last_file_in)
!*************************************************************
!* Write a text-based checkpoint to crest.restart.
!* Called after each MTD iteration and after collectcre.
!*
!* Arguments:
!*   runtype      - crestver constant (crest_imtd or crest_imtd2)
!*   stage        - stage label: 'mtd_loop', 'post_collect', 'done'
!*   main_iter    - env%nreset (MAINLOOP iteration counter)
!*   mtd_iter     - last completed MTD iteration (0 for post-loop stages)
!*   nmetadyn     - env%nmetadyn (may differ from initial after first pass)
!*   elowest      - current lowest energy
!*   eprivious    - previous lowest energy
!*   last_file_in - last CREGEN-sorted file written to disk
!*************************************************************
    implicit none
    integer,intent(in)          :: runtype,main_iter,mtd_iter,nmetadyn
    character(len=*),intent(in) :: stage,last_file_in
    real(wp),intent(in)         :: elowest,eprivious
    integer :: ich,io
    open(newunit=ich,file=restart_file,status='replace',iostat=io)
    if (io /= 0) then
      write(stdout,'(a)') '**WARNING** could not write crest.restart'
      return
    end if
    write(ich,'(a)') '# CREST restart checkpoint - do not edit manually'
    write(ich,'(a,1x,i0)') 'version',   1
    write(ich,'(a,1x,i0)') 'runtype',   runtype
    write(ich,'(a,1x,i0)') 'main_iter', main_iter
    write(ich,'(a,1x,i0)') 'mtd_iter',  mtd_iter
    write(ich,'(a,1x,i0)') 'nmetadyn',  nmetadyn
    write(ich,'(a,1x,a)')  'stage',     trim(stage)
    write(ich,'(a,1x,a)')  'last_file', trim(last_file_in)
    write(ich,'(a,1x,f25.15)') 'elowest',   elowest
    write(ich,'(a,1x,f25.15)') 'eprivious', eprivious
    close(ich)
  end subroutine write_restart_log

!========================================================================================!

  subroutine delete_restart_log()
!*************************************
!* Remove crest.restart from the
!* current directory. Called on
!* successful completion of a run so
!* no stale checkpoint is left behind.
!*************************************
    implicit none
    integer :: ich,io
    logical :: ex
    inquire(file=restart_file,exist=ex)
    if (.not.ex) return
    open(newunit=ich,file=restart_file,status='old',iostat=io)
    if (io == 0) close(ich,status='delete')
  end subroutine delete_restart_log

!========================================================================================!

  subroutine read_restart_log(rdat)
!*************************************************************
!* Read crest.restart into a restart_data object.
!* Unknown keys are silently ignored for forward compatibility.
!*
!* Arguments:
!*   rdat - restart_data object to populate
!*************************************************************
    implicit none
    type(restart_data),intent(out) :: rdat
    integer :: ich,io
    character(len=512) :: line,key,val
    integer :: pos

    rdat = restart_data()  !> initialise with defaults

    open(newunit=ich,file=restart_file,status='old',iostat=io)
    if (io /= 0) then
      write(stdout,'(a)') '**WARNING** could not read crest.restart'
      return
    end if

    do
      read(ich,'(a)',iostat=io) line
      if (io /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle
      pos = index(line,' ')
      if (pos < 2) cycle
      key = line(1:pos-1)
      val = adjustl(line(pos+1:))
      select case(trim(key))
      case('version')
        read(val,*,iostat=io) rdat%version
      case('runtype')
        read(val,*,iostat=io) rdat%runtype
      case('main_iter')
        read(val,*,iostat=io) rdat%main_iter
      case('mtd_iter')
        read(val,*,iostat=io) rdat%mtd_iter
      case('nmetadyn')
        read(val,*,iostat=io) rdat%nmetadyn
      case('stage')
        rdat%stage = trim(val)
      case('last_file')
        rdat%last_file = trim(val)
      case('elowest')
        read(val,*,iostat=io) rdat%elowest
      case('eprivious')
        read(val,*,iostat=io) rdat%eprivious
      end select
    end do
    close(ich)
  end subroutine read_restart_log

!========================================================================================!

  subroutine print_restart_info(rdat)
!*****************************************************
!* Print a summary of the restart checkpoint to stdout.
!*****************************************************
    implicit none
    type(restart_data),intent(in) :: rdat
    character(len=64) :: rtname
    integer :: w
    w = 57

    select case(rdat%runtype)
    case(2)
      rtname = 'iMTD-GC'
    case(22)
      rtname = 'sMTD-iMTD (entropy)'
    case default
      write(rtname,'(a,i0)') 'runtype ',rdat%runtype
    end select

    write(stdout,*)
    write(stdout,'(1x,a)') repeat(':',w)
    write(stdout,'(1x,a,a,a)') ' RESTART DETECTED (',trim(restart_file),')'
    write(stdout,'(1x,a,a)')   '  runtype  : ',trim(rtname)
    write(stdout,'(1x,a,a)')   '  stage    : ',trim(rdat%stage)
    if (rdat%stage == 'mtd_loop') then
      write(stdout,'(1x,a,i0,a,i0,a)') '  MTD iter : ',rdat%mtd_iter, &
        &  ' (MAINLOOP ',rdat%main_iter,')'
    end if
    if (len_trim(rdat%last_file) > 0) then
      write(stdout,'(1x,a,a)')   '  last file: ',trim(rdat%last_file)
    end if
    write(stdout,'(1x,a,f20.10)') '  elowest  : ',rdat%elowest
    write(stdout,'(1x,a)') repeat(':',w)
    write(stdout,*)
  end subroutine print_restart_info

!========================================================================================!
!========================================================================================!
end module crest_restartlog
