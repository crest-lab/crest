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

!> A data storage module for hosting ConfSolv as a http server with requests
!> We need to know the server's PORT to create the request
!> CREST can either try to host the server on its own, or the user
!> must provide the PORT.
!> All of this assumes that the ConfSolv submodule was loaded

module ConfSolv_module
  use crest_parameters
  use crest_data
  use strucrd
  use iomod
  implicit none
  public

  !> ConfSolv helper script PID
  integer,allocatable :: cs_pid
  !> ConfSolv helper script name
  character(len=:),allocatable :: cs_bin
  !> ConfSolv port server port
  integer,allocatable :: cs_port
  !> ConfSolv teardown instruction
  logical :: cs_teardown = .false.
  !> Keeping track of setup. Has it been called already?
  logical :: cs_setup = .false. 

  !> ConfSolv parameter location
  character(len=:),allocatable :: cs_param
  !> ConfSolv solvent & smiles
  character(len=:),allocatable :: cs_solvent
  character(len=:),allocatable :: cs_smiles

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine cs_deallocate()
    if (allocated(cs_pid)) deallocate (cs_pid)
    if (allocated(cs_bin)) deallocate (cs_bin)
    if (allocated(cs_port)) deallocate (cs_port)
    if (allocated(cs_param)) deallocate (cs_param)
    if (allocated(cs_solvent)) deallocate (cs_solvent)
  end subroutine cs_deallocate

!=================================================!

  function cs_running() result(running)
    implicit none
    logical :: running
    integer :: io
    character(len=:),allocatable :: job
    running = .false.

    io = 0
    if (.not.allocated(cs_bin).or. &
    &  .not.allocated(cs_pid).or. &
    &  .not.allocated(cs_port)) then
      running = .false.
      return
    end if

    job = trim(cs_bin)//' --test '//to_str(cs_port)

    call command(job,io)

  end function cs_running

!=================================================!
  subroutine cs_shutdown(io)
    implicit none
    integer,intent(out) :: io

    io = 0
    if (cs_teardown.and.allocated(cs_pid).and.allocated(cs_port)) then
      write (stdout,'(/,a,i0)') 'Shutting down http://localhost/',cs_port
      call kill(cs_pid,9,io)
      deallocate (cs_pid)
    end if

  end subroutine cs_shutdown

!=================================================!

  subroutine cs_deploy()
    implicit none
    character(len=:),allocatable :: job
    integer :: io,ich
    character(len=50) :: atmp
    logical :: ex

    if(.not.allocated(cs_bin)) cs_bin = 'confsolvserver'
    call remove('confsolv.out')
    call remove('config_template.toml')

    job = 'nohup '//trim(cs_bin)//' -l '//'> confsolv.out 2>/dev/null'//' &'

    write(stdout,'(2x,a,a)') 'Hosting command: ',trim(job)
    call command(job,io)

    if(io /= 0) error stop '**ERROR** failed to host ConfSolv server'
    !call sleep(3)
    do 
      call sleep(1)
      inquire(file='config_template.toml',exist=ex)
      if(ex) exit
    enddo

    !> read port and pid
    open(newunit=ich,file='confsolv.out')
    read(ich,*) atmp,cs_pid
    read(ich,*) atmp,cs_port
    close(ich) 
    cs_teardown = .true.
    cs_setup = .false.
    write(stdout,'(2x,2(a,i0))') 'ConfSolv server will be running at http://localhost:',cs_port,' with PID ',cs_pid
  end subroutine cs_deploy

!=================================================!

  subroutine cs_write_config(ensname,threads)
    implicit none
    character(len=*) :: ensname
    integer :: threads
    integer :: i,j,k,l,ich,io
    character(len=1024) :: atmp

    call getcwd(atmp)
    open (newunit=ich,file='config.toml')
    call wrtoml_int(ich,'port',cs_port)
    call wrtoml_int(ich,'pid',cs_pid)
    call wrtoml_int(ich,'num_cores',threads)
    if (allocated(cs_param)) then
      call wrtoml(ich,'model_path',cs_param)
    end if
    
    call wrtoml(ich,'xyz_file',trim(atmp)//'/'//trim(ensname))
    call wrtoml(ich,'solvent_file',trim(atmp)//'/'//'crest_solvents.csv')

    call wrtoml_int(ich,'n_threshold_mols',1)
    close (ich)
  contains
    subroutine wrtoml(ch,key,val)
      integer :: ch
      character(len=*) :: key
      character(len=*) :: val
      write (ch,'(a,a,a,a)') trim(key),' =  "',trim(val),'"'
    end subroutine wrtoml
    subroutine wrtoml_int(ch,key,val)
      integer :: ch
      character(len=*) :: key
      integer :: val
      write (ch,'(a,a,i0)') trim(key),' =  ',val
    end subroutine wrtoml_int
  end subroutine cs_write_config

!========================================================================================!

  subroutine cs_write_solvent_csv(solvent,smiles,ch)
!**************************************************************
!* From CREST's side it makes only sense to define ONE solvent
!* despite ConfSolv being able to handle multiple.
!* ConfSolv will read the solvents from a CSV file with the
!* columns SOLVENT_NAME and SMILES
!**************************************************************
    implicit none
    character(len=*),intent(in) :: solvent
    character(len=*),intent(in),optional :: smiles
    integer,intent(in),optional :: ch
    integer :: ich
    if (.not.present(ch)) then
      open (newunit=ich,file='crest_solvents.csv')
    else
      ich = ch
    end if
    !> column names
    write (ich,'(a,",",a)') 'SOLVENT_NAME','SMILES'
    if (present(smiles)) then
      write (ich,'(a,",",a)') solvent,smiles
    else
      !> switch case for available solvents, if no smiles was given
      select case (lowercase(solvent))
      case ('acetate')
        write (ich,'(a,",",a)') solvent,'CC(=O)[O-]'
      case ('acetic acid')
        write (ich,'(a,",",a)') solvent,'CC(=O)O'
      case ('acetone')
        write (ich,'(a,",",a)') solvent,'CC(=O)C'
      case ('acetonitrile')
        write (ich,'(a,",",a)') solvent,'CC#N'
      case ('ammonia')
        write (ich,'(a,",",a)') solvent,'N'
      case ('ammonium')
        write (ich,'(a,",",a)') solvent,'[NH4+]'
      case ('benzene')
        write (ich,'(a,",",a)') solvent,'c1ccccc1'
      case ('benzoate')
        write (ich,'(a,",",a)') solvent,'[O-]C(=O)c1ccccc1'
      case ('benzylacetate')
        write (ich,'(a,",",a)') solvent,'CC(=O)OCc1ccccc1'
      case ('butanone','2-butanone')
        write (ich,'(a,",",a)') solvent,'CCC(=O)C'
      case ('chloride')
        write (ich,'(a,",",a)') solvent,'[Cl-]'
      case ('trichlormethane')
        write (ich,'(a,",",a)') solvent,'C(Cl)(Cl)Cl'
      case ('cyclohexane')
        write (ich,'(a,",",a)') solvent,'C1CCCCC1'
      case ('dibutylamine')
        write (ich,'(a,",",a)') solvent,'CC[C@H](C)N[C@H](C)CC'
      case ('dichlormethane')
        write (ich,'(a,",",a)') solvent,'C(Cl)Cl'
      case ('diethanolamine')
        write (ich,'(a,",",a)') solvent,'OCCNCCO'
      case ('diethanolammonium')
        write (ich,'(a,",",a)') solvent,'OCC[NH2+]CCO'
      case ('diethylamine')
        write (ich,'(a,",",a)') solvent,'CCNCC'
      case ('diethylammonium')
        write (ich,'(a,",",a)') solvent,'CC[NH2+]CC'
      case ('diethylether')
        write (ich,'(a,",",a)') solvent,'CCOCC'
      case ('heptyloctylether')
        write (ich,'(a,",",a)') solvent,'CCCCCCCCOCCCCCCC'
      case ('acetamide')
        write (ich,'(a,",",a)') solvent,'CC(=O)N(C)C'
      case ('diethylformamide')
        write (ich,'(a,",",a)') solvent,'CN(C)C=O'
      case ('dmso')
        write (ich,'(a,",",a)') solvent,'CS(=O)C'
      case ('dioxolone','2-dioxolone')
        write (ich,'(a,",",a)') solvent,'C1COC(=O)O1'
      case ('ethylmethylester')
        write (ich,'(a,",",a)') solvent,'CCOC(=O)OC'
      case ('ethanol')
        write (ich,'(a,",",a)') solvent,'CCO'
      case ('ethylacetate')
        write (ich,'(a,",",a)') solvent,'CCOC(=O)C'
      case ('ethylamine')
        write (ich,'(a,",",a)') solvent,'CCN'
      case ('ethylaminium')
        write (ich,'(a,",",a)') solvent,'CC[NH3+]'
      case ('glycol')
        write (ich,'(a,",",a)') solvent,'OCCO'
      case ('formate')
        write (ich,'(a,",",a)') solvent,'C(=O)[O-]'
      case ('formic acid')
        write (ich,'(a,",",a)') solvent,'C(=O)O'
      case ('butyrolacetone')
        write (ich,'(a,",",a)') solvent,'O=C1CCCO1'
      case ('glycerin')
        write (ich,'(a,",",a)') solvent,'OCC(O)CO'
      case ('water','h2o')
        write (ich,'(a,",",a)') solvent,'O'
      case ('sulfuric acid')
        write (ich,'(a,",",a)') solvent,'O=S(=O)(O)O'
      case ('hexafluorobenzene')
        write (ich,'(a,",",a)') solvent,'Fc1c(F)c(F)c(F)c(F)c1F'
      case ('isooctane')
        write (ich,'(a,",",a)') solvent,'CC(C)CC(C)(C)C'
      case ('isopropanol')
        write (ich,'(a,",",a)') solvent,'CC(O)C'
      case ('methoxide')
        write (ich,'(a,",",a)') solvent,'C[O-]'
      case ('hexane','n-hexane')
        write (ich,'(a,",",a)') solvent,'CCCCCC'
      case ('1-nonadecanol','nonadecanol')
        write (ich,'(a,",",a)') solvent,'CCCCCCCCCCCCCCCCCCCO'
      case ('1-octanol','octanol')
        write (ich,'(a,",",a)') solvent,'OCCCCCCCC'
      case ('p-dichlorobenzene','dichlorobenzene')
        write (ich,'(a,",",a)') solvent,'Clc1ccccc1Cl'
      case ('perfluorohexane')
        write (ich,'(a,",",a)') solvent,'C(C(C(C(F)(F)F)(F)F)(F)F)(C(C(F)(F)F)(F)F)(F)F'
      case ('propanediol')
        write (ich,'(a,",",a)') solvent,'C[C@@H](O)CO'
      case ('tetraethylammoniom')
        write (ich,'(a,",",a)') solvent,'CC[N+](CC)(CC)CC'
      case ('thf','tetrahydrofuran')
        write (ich,'(a,",",a)') solvent,'O1CCCC1'
      case ('toluene')
        write (ich,'(a,",",a)') solvent,'Cc1ccccc1'
      case ('tributylphosphate')
        write (ich,'(a,",",a)') solvent,'O=P(OCCCC)(OCCCC)OCCCC'
      case ('triethanolamine','trolamine')
        write (ich,'(a,",",a)') solvent,'OCCN(CCO)CCO'
      case ('triethanolammonium')
        write (ich,'(a,",",a)') solvent,'OCC[NH+](CCO)CCO'
      case ('triethylamine','net3')
        write (ich,'(a,",",a)') solvent,'CCN(CC)CC'
      case ('triethylammonium')
        write (ich,'(a,",",a)') solvent,'CC[NH+](CC)CC'
      case ('triglyme')
        write (ich,'(a,",",a)') solvent,'COCCOCCOCCOC'
      case ('urea')
        write (ich,'(a,",",a)') solvent,'NC(=O)N'
      case default
        write(stderr,'(2a)') '**ERROR** failed to find matching solvent SMILES for: ',solvent
        error stop
      end select
    end if
    close (ich)
  end subroutine cs_write_solvent_csv

!========================================================================================!
!========================================================================================!
end module ConfSolv_module

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!

subroutine confsolv_request(ensname,nall,ncpus,gsoln,io)
!***********************************************
!* Interface to ConfSolv via https requests
!*
!* Input/Output:
!*  ensname - ensemble file
!*  nall    - number of structures in ensemble
!*  ncpus   - number on cores to run on
!*  gsoln   - ΔΔGsoln for each structure
!*  io      - exit status 
!***********************************************
  use crest_parameters
  use crest_data
  use strucrd
  use iomod
  use ConfSolv_module
  use parse_csv
  implicit none
  !> INPUT
  character(len=*),intent(in) :: ensname
  integer,intent(in) :: nall
  integer,intent(in) :: ncpus
  !> OUTPUT
  real(wp),intent(out) :: gsoln(nall)
  integer,intent(out)  :: io
  !> LOCAL
  integer :: i,j,k,l,ich
  logical :: pr,wr
  character(len=:),allocatable :: job
  real(wp),allocatable :: column(:)
  real(wp) :: avg

  io = 0
  gsoln(:) = 0.0_wp

!>--- setup
  if (allocated(cs_pid).and.allocated(cs_port)) then
  !> user-provided PID and port (no automatic teardown)
    write (stdout,'(2x,a,i0,a,i0)') 'Looking for ConfSolv server (PID: ',cs_pid,') running at '//&
    & 'http://localhost:',cs_port
    !cs_teardown = .false. 
  else
  !> fallback: automatic host (not recommended)
    allocate(cs_pid,cs_port)
    call cs_deploy() 
  end if
  if (allocated(cs_param)) then
    write (stdout,'(2x,a,/,3x,a)') 'pyTorch checkpoint files located at ',cs_param
  else
    write (stderr,*) '**ERROR** cannot run ConfSolv without defining checkpoint file location!'
    error stop
  end if
  if (allocated(cs_solvent).and.allocated(cs_smiles)) then
    write (stdout,'(2x,a,a,3a)') 'Requested ΔΔGsoln for ',cs_solvent,' (SMILES: ',trim(cs_smiles),')'
    call cs_write_solvent_csv(cs_solvent,smiles=cs_smiles)
  else if(allocated(cs_solvent))then
    write (stdout,'(2x,a,a,a)') 'Requested ΔΔGsoln for ',cs_solvent,' (trying to find SMILES ...)'
    call cs_write_solvent_csv(cs_solvent)
  end if
  write (stdout,'(2x,a,a)') 'Processing ensemble file ',trim(ensname)

!>---- creating the request configuration
  write (stdout,'(2x,a)',advance='no') 'Writing config.toml file              ...'
  flush(stdout) 
  call cs_write_config(ensname,ncpus)
  write(stdout,*) 'done.'


  job = ''
  job = trim(job)//' '//cs_bin//' -c config.toml'
!>----- this should only be called once:
  if(.not.cs_setup)then
    write(stdout,'(2x,a)',advance='no') 'Instructing ConfSolv model setup      ...'
    flush(stdout)
    call command(trim(job)//' -s >> confsolv.out 2>/dev/null',io)
    if(io/=0)then
      write(stdout,*)
      write(stderr,'(a)')"**ERROR** failed request to ConfSolv server"
      call cs_shutdown(io)
      error stop
    endif 
    cs_setup = .true.
    write(stdout,*) 'done.'
  endif 

!>---- and then the actual evaluation
  call remove('confsolv.csv')
  call remove('confsolv_uncertainty.csv')
  write(stdout,'(2x,a)',advance='no') 'Evaluation of ConfSolv D-MPNN         ...'
  flush(stdout) 
  call command(trim(job)//' -r >> confsolv.out 2>/dev/null',io)
  write(stdout,*) 'done.'
  if(io/=0)then
    write(stdout,*)
    write(stderr,'(a)')"**ERROR** failed request to ConfSolv server"
    call cs_shutdown(io)
    error stop
  endif

!>--- read ΔΔGsoln
  write(stdout,'(2x,a)',advance='no') 'Reading confsolv.csv                  ...'
  flush(stdout)
  call parse_csv_column_real('confsolv.csv',cs_solvent,column)
  write(stdout,*) 'done.'
  if(size(column,1) == nall)then
    gsoln(:) = column(:)
  else
    write(stdout,'(a)') '**ERROR** dimension mismatch in confsolv_request'
    call cs_shutdown(io)
    error stop
  endif

!>--- read uncertainty
  write(stdout,'(2x,a)',advance='no') 'Reading confsolv_uncertainty.csv      ...'
  flush(stdout)
  call parse_csv_column_real('confsolv_uncertainty.csv',cs_solvent,column)
  write(stdout,*) 'done.'
  if(size(column,1) == nall)then
    avg = sum(column(:))/real(nall,wp)
    write(stdout,'(2x,a,f25.15)') 'Average uncertainty of ConfSolv prediction:',avg
  else
    write(stdout,'(a)') '**ERROR** dimension mismatch in confsolv_request'
    call cs_shutdown(io)
    error stop
  endif

  if(allocated(column)) deallocate(column)
  return
end subroutine confsolv_request

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!
