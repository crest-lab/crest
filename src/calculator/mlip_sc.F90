!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2026 Philipp Pracht
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

!> module mlip_sc
!> A module containing routines for calling MLIPs though persistent python instances
!> enabled through the fmlip_relay submodule

!=========================================================================================!
module mlip_sc
  use crest_parameters
  use strucrd
  use iomod
#ifdef WITH_FMLIP_RELAY
  use fmlip_relay_client
#endif
  implicit none
  !>--- private module variables and parameters
  private

  character(len=*),parameter :: basebin = 'fmlip-relay-server'

  public :: mlip_params
  type :: mlip_params
    integer :: BASE_PORT = 54320
    integer :: TIMEOUT_SEC = 120
    character(len=:),allocatable :: backend
    character(len=:),allocatable :: modelpath
    character(len=:),allocatable :: modelsize
    !> shared neural-network options (mace* / uma backends)
    character(len=:),allocatable :: device    !> torch device: cpu | cuda | cuda:0
    !> FairChem UMA backend options (--backend uma)
    character(len=:),allocatable :: umamodel  !> checkpoint, e.g. uma-s-1p2 (default), uma-m-1
    character(len=:),allocatable :: umatask   !> task head: omol | omat | omc | oc20 | odac
    integer :: max_threads = 0 !> CPU thread cap per server (--max-threads), synced from
                               !> the level's threads setting; <1 = unset (inherit env)
    integer :: iid = 0
  end type mlip_params

  public :: mlip_engrad_core,fmlip_relay_init,mlips_shutdown

  integer,parameter  :: nopbc(3) = (/0,0,0/)
  integer,parameter  :: allpbc(3) = (/1,1,1/)
  real(wp),parameter :: bigcell(3,3) = reshape( &
    & (/10000.0_wp,0.0_wp,0.0_wp, &
    &   0.0_wp,10000.0_wp,0.0_wp, &
    &   0.0_wp,0.0_wp,10000.0_wp/), [3,3])

  external creststop
!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine fmlip_relay_init(MPAR,iid)
    type(mlip_params),intent(inout) :: MPAR
    integer,intent(in) :: iid
    integer :: io,tmpport
    character(len=256) :: cmd,cmd_0,cmd_1
#ifdef WITH_FMLIP_RELAY
    if (.not.allocated(MPAR%backend)) then
      write (stdout,*)
      write (stdout,*) '** ERROR ** No model backend selected for MLIP'
      write (stdout,*)
      call creststop(20)
    end if

    ! ── fast path: instance already running, nothing to do ───────────────────
    call mlip_ping(iid,io)
    if (io == MLIP_OK) then
      MPAR%iid = iid
      return
    end if

    call checkprog_silent(basebin,verbose=.false.,iostat=io)
    if (io .ne. 0) then
      write (stdout,*)
      write (stdout,*) '** ERROR ** can not find socket server for MLIPs '//basebin
      write (stdout,*) ' Make sure you install it from the fmlip_relay subproject via pip'
      write (stdout,*)
      call creststop(20)
    end if

    !> check if we have limitations for parallelity
    if (iid > MLIP_MAX_INSTANCES) then
      write (stdout,*)
      write (stdout,*) '** ERROR ** exeeding the max number of parallel socket servers for MLIPs '
      write (stdout,*) ' Please request fewer than '//to_str(MLIP_MAX_INSTANCES)
      write (stdout,*)
      call creststop(20)
    end if

    !> options prepping
    tmpport = MPAR%BASE_PORT+iid
    write (cmd_1,'("--dtype float64")')
    !> cap the server's CPU inference threads (torch/BLAS pools) to the
    !> per-level core reservation so parallel instances don't oversubscribe
    if (MPAR%max_threads > 0) then
      write (cmd_1,'(a,1x,a,1x,i0)') trim(cmd_1),'--max-threads',MPAR%max_threads
    end if

    select case (MPAR%backend)
    case ('mace_off','mace_mp')

      if (allocated(MPAR%modelpath)) then
        if (.not.file_exists(MPAR%modelpath)) then
          write (stdout,*)
          write (stdout,*) '** ERROR ** model path allocated but can not find '//trim(MPAR%modelpath)
          write (stdout,*)
          call creststop(20)
        end if
        !> a user-provided checkpoint is served through the generic 'mace' backend
        cmd_0 = '--backend mace --model '//trim(MPAR%modelpath)
      else
        cmd_0 = '--backend '//trim(MPAR%backend)
        if (allocated(MPAR%modelsize)) cmd_0 = trim(cmd_0)//' --mace-model '//trim(MPAR%modelsize)
      end if
      if (allocated(MPAR%device)) cmd_0 = trim(cmd_0)//' --device '//trim(MPAR%device)
      write (cmd,'(a,1x,a,1x,i0,1x,a,1x,a)') basebin,'--port',tmpport,trim(adjustl(cmd_0)),trim(cmd_1)

    case ('uma')
      !> FairChem UMA foundation model (fairchem-core v2). Charge and spin
      !> multiplicity are forwarded per-call via the relay protocol; only the
      !> checkpoint, task head and torch device are fixed at server startup.
      cmd_0 = ''
      if (allocated(MPAR%umamodel)) cmd_0 = trim(cmd_0)//' --uma-model '//trim(MPAR%umamodel)
      if (allocated(MPAR%umatask))  cmd_0 = trim(cmd_0)//' --uma-task '//trim(MPAR%umatask)
      if (allocated(MPAR%device))   cmd_0 = trim(cmd_0)//' --device '//trim(MPAR%device)
      write (cmd,'(a,1x,a,1x,i0,1x,a,1x,a,1x,a)') basebin,'--port',tmpport, &
        & '--backend uma',trim(adjustl(cmd_0)),trim(cmd_1)

    case default

      if (allocated(MPAR%modelpath)) then
        if (.not.file_exists(MPAR%modelpath)) then
          write (stdout,*)
          write (stdout,*) '** ERROR ** model path allocated but can not find '//trim(MPAR%modelpath)
          write (stdout,*)
          call creststop(20)
        end if
        write (cmd,'(a,1x,a,1x,i0,2(1x,a,1x,a),1x,a)') basebin,'--port',tmpport,'--backend', &
          & trim(MPAR%backend),'--model',trim(MPAR%modelpath),trim(cmd_1)
      else
        !> no model path (e.g. lj, dummy backends that need no model file)
        write (cmd,'(a,1x,a,1x,i0,1x,a,1x,a,1x,a)') basebin,'--port',tmpport, &
          & '--backend',trim(MPAR%backend),trim(cmd_1)
      end if
    end select

    !> spawn the server and verify
    call mlip_init(iid,tmpport,trim(cmd)//' 2>/dev/null',MPAR%TIMEOUT_SEC,io)
    if (io /= MLIP_OK) then
      write (stdout,*)
      write (stdout,*) '** ERROR ** failed to initialize MLIP server'
      write (stdout,*)
      call creststop(1)
    end if
    call mlip_ping(iid,io)
    if (io /= MLIP_OK) then
      write (stdout,*)
      write (stdout,*) '** ERROR ** failed to ping MLIP server'
      call creststop(1)
    end if

    MPAR%iid = iid

#else /* WITH_FMLIP_RELAY */
    write (stdout,*) 'Error: Compiled without fmlip-relay support!'
    write (stdout,*) 'Use -DWITH_FMLIP_RELAY=true in the setup to enable this function'
    write (stdout,*)
    call creststop(20)
#endif
  end subroutine fmlip_relay_init

  subroutine mlip_engrad_core(mol,MPAR,energy,gradient,iostatus, &
      &                       charge,spin,iid)
    type(coord),intent(in) :: mol
    type(mlip_params),intent(in)    :: MPAR
    integer,intent(in),optional :: charge
    integer,intent(in),optional :: spin
    integer,intent(in),optional :: iid
    real(wp),intent(out)   :: energy
    real(wp),intent(out)   :: gradient(3,mol%nat)
    integer,intent(out)    :: iostatus

    integer :: chrg,spn,instance_id
    real(wp) :: stress(3,3)

    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
    iostatus = 1

    chrg = 0
    spn = 1
    if (present(charge)) chrg = charge
    if (present(spin)) spn = spin
    instance_id = MPAR%iid
    if (present(iid)) instance_id = iid

#ifdef WITH_FMLIP_RELAY
    if (allocated(mol%lat)) then
      call mlip_compute(instance_id,mol%nat,mol%at,mol%xyz*autoaa,mol%lat,allpbc,0,chrg,spn, &
      &                 energy,gradient,stress,iostatus)
    else
      call mlip_compute(instance_id,mol%nat,mol%at,mol%xyz*autoaa,bigcell,nopbc,0,chrg,spn, &
      &                 energy,gradient,stress,iostatus)
    end if

    !> CREST always works with atomic units, convert from eV and Angstroem:
    energy = energy/autoev
    gradient(:,:) = -gradient(:,:)*(1.0_wp/(autoev*aatoau))
#endif
  end subroutine mlip_engrad_core

!========================================================================================!

  subroutine mlips_shutdown()
    integer :: io
#ifdef WITH_FMLIP_RELAY
    call mlip_finalize_all(io)
#endif
  end subroutine mlips_shutdown

!========================================================================================!
end module mlip_sc
