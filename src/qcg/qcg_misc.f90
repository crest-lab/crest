!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 Christoph Plett, Sebastian Spicher, Philipp Pracht
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

subroutine xtb_sp_qcg(env,fname,success,eout)
!********************************************************
!*  xtb_sp_qcg
!*  A quick single point xtb calculation without wbo
!********************************************************
  use crest_parameters
  use iomod
  use crest_data
  use crest_calculator
  use strucrd
  implicit none
  type(systemdata) :: env
  character(len=*),intent(in) :: fname
  logical,intent(out) :: success
  real(wp),intent(out) :: eout

  character(len=512) :: jobcall
  character(len=*),parameter :: pipe = ' > xtb.out 2> /dev/null'
  logical,parameter :: debug = .false.
  integer :: io,T,Tn

  success = .false.
  eout = 0.0_wp

  if (env%legacy) then
!>----------------------------------------------
!> The original implementation with systemcall
    call remove('gfnff_topo')
    call remove('energy')
    call remove('charges')
    call remove('xtbrestart')

!---- setting threads
    call new_ompautoset(env,'auto',1,T,Tn)

!---- jobcall
    write (jobcall,'(a,1x,a,1x,a,'' --sp '',a,1x,a)') &
    &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),trim(pipe)

    if (debug) write (stdout,*) trim(jobcall)
    call command(trim(jobcall),io)
    call grepval('xtb.out','| TOTAL ENERGY',success,eout)
!---- cleanup
    call remove('energy')
    call remove('charges')
    call remove('xtbrestart')
    call remove('xtbtopo.mol')
    call remove('gfnff_topo')
  else
!>---------------------------------------------
!> New implementation with calculator and api
    block
      type(calcdata) :: calc
      type(coord) :: mol
      real(wp),allocatable :: gradtmp(:,:)

      call mol%open(fname)
      allocate (gradtmp(3,mol%nat))
      call env2calc(env,calc,mol)
      if (debug) call calc%info(stdout)
      call engrad(mol,calc,eout,gradtmp,io)
      success = (io == 0)
    end block
  end if
end subroutine xtb_sp_qcg

!--------------------------------------------------------------------------------------------
! A quick single xtb optimization gets mol and overwrites it with optimized stuff
!--------------------------------------------------------------------------------------------
subroutine xtb_opt_qcg(env,mol,constrain)
  use crest_parameters
  use iomod
  use crest_data
  use qcg_coord_type
  use strucrd
  use optimize_module

  implicit none
  type(systemdata),intent(in) :: env
  type(coord_qcg),intent(inout) :: mol

  character(:),allocatable :: fname
  character(len=512) :: jobcall
  logical :: constrain,const
  real(wp) :: energy
  integer :: io,T,Tn
  character(stdout),parameter :: pipe = ' > xtb_opt.out 2> /dev/null'
  logical,parameter :: debug = .false.

  if (env%legacy) then
    !> LEGACY version with syscall

    !--- Write coordinated
    fname = 'coord'
    call mol%write(fname)

    !---- setting threads
    call new_ompautoset(env,'auto',1,T,Tn)

    !---- jobcall & Handling constraints
    if (constrain.AND.env%cts%used) then
      call write_constraint(env,fname,'xcontrol')
      call mol%write('coord.ref')
      write (jobcall,'(a,1x,a,1x,a,'' --opt --input xcontrol '',a,1x,a)') &
      &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),trim(pipe)
    else
      write (jobcall,'(a,1x,a,1x,a,'' --opt '',a,1x,a)') &
      &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),trim(pipe)
    end if

    call command(trim(jobcall),io)
    !---- cleanup
    call rdcoord('xtbopt.coord',mol%nat,mol%at,mol%xyz)
    call remove('energy')
    call remove('charges')
    call remove('xtbrestart')
    call remove('xtbtopo.mol')
    call remove('gfnff_topo')

  else

    !> NEW version with calculator
    call new_ompautoset(env,'max',1,T,Tn)
    block
      use crest_calculator
      use optimize_module
      type(calcdata) :: calc
      type(coord) :: molin,molout
      real(wp),allocatable :: gradtmp(:,:)

      allocate (gradtmp(3,mol%nat))
      molin = mol%as_coord()
      call env2calc(env,calc,molin)
      if (debug) call calc%info(stdout)

      call optimize_geometry(molin,molout,calc,energy,gradtmp,debug,.false.,io)

      deallocate (gradtmp)
      if (io == 0) then
        call mol%from_coord(molout)
      else
        write (stdout,*) 'FAILURE in QCG optimization!'
        write (stdout,*) 'Stopping run to avoid unecessary compuations'
        call creststop(status_safety)
      end if
    end block

  end if
end subroutine xtb_opt_qcg

subroutine xtb_md_ensemble_qcg(env,solu,solv,clus,resultspath)
  use crest_parameters
  use crest_data
  use qcg_coord_type
  use strucrd
  use iomod
  use qcg_printouts
  use crest_calculator
  use dynamics_module
  implicit none
  !> IN/OUTPUTS
  type(systemdata),intent(inout) :: env
  type(coord_qcg),intent(inout)  :: solu
  type(coord_qcg),intent(inout)  :: solv
  type(coord_qcg),intent(inout)  :: clus
  character(len=*),intent(in) :: resultspath
  !> LOCAL
  integer  :: T,Tn,i,j,k,l,ich,r,io
  real(wp) :: newtemp,newmdtime,newmdstep,newhmass
  real(wp) :: newmetadfac,newmetadexp,newmetadlist
  character(len=:),allocatable :: fname
  character(len=512) :: tmppath,tmppath2
  character(len=1024) :: jobcall
  logical :: ex,mdfail
  type(ensemble) :: dum
  type(calcdata),target :: calc
  type(mddata) :: mddat
  type(coord)  :: mol
  type(mtdpot) :: mtd
  real(wp),allocatable :: cn(:),fakewbo(:,:)
  character(len=*),parameter :: pipe = ' > xtb.out 2>/dev/null'

  call getcwd(tmppath2)

  !---- Setting threads
  call new_ompautoset(env,'auto',1,T,Tn)

  !--- Setting new defaults for MD/MTD in qcg
  if (env%mdtemp .lt. 0.0d0) then
    newtemp = 400.00d0
  else if (.not.env%user_temp) then
    newtemp = 298.0
  else
    newtemp = env%mdtemp
  end if
  env%mdtemp = newtemp

  if (.not.env%user_mdtime) then
    newmdtime = 100.0 !100.0
  else
    newmdtime = env%mdtime
  end if
  env%mdtime = newmdtime

  if (.not.env%user_dumxyz) then
    env%mddumpxyz = 1000
  end if

  if (.not.env%user_mdstep) then
    if (env%ensemble_opt .ne. '--gff') then
      newmdstep = 4.0d0
    else
      newmdstep = 1.5d0
    end if
  else
    newmdstep = env%mdstep
  end if
  env%mdstep = newmdstep

  if (env%ensemble_opt .ne. '--gff') then
    newhmass = 4.0
  else
    newhmass = 5.0
  end if
  env%hmass = newhmass

  if (.not.allocated(env%metadfac)) then
    allocate (env%metadfac(1),source=0.02_wp)
    allocate (env%metadexp(1),source=0.1_wp)
    allocate (env%metadlist(1),source=10)
  end if
  newmetadfac = 0.02_wp
  newmetadexp = 0.1_wp
  newmetadlist = 10.0_wp

  fname = 'coord'

!> --------------------------------------------------------------------
!> Internal calculator version
!> --------------------------------------------------------------------
  if (.not.env%legacy) then
    call mol%open(fname)

    call env_to_mddat(env)
    mddat = env%mddat
    calc = env%calc
    if (.not.env%solv_md) then
      call calc%set_freeze(mol%nat,1,solu%nat)
    end if

    if (env%shake .ne. 0) then
      call mol%cn_to_bond(cn,fakewbo)
      call move_alloc(fakewbo,mddat%shk%wbo)
      if (calc%nfreeze > 0) then
        mddat%shk%freezeptr => calc%freezelist
      end if
    end if

    !> for MTD runtype add the corresponding potential
    if (env%ensemble_method .EQ. 2) then
      mtd%kpush = newmetadfac
      mtd%alpha = newmetadexp
      mtd%cvdump_fs = real(env%mddump)
      mtd%mtdtype = cv_rmsd
      allocate (mtd%atinclude(mol%nat),source=.true.)
      mtd%atinclude(1:clus%nat) = .false. !> only include solvent
      call mddat%add(mtd)
    end if

    !> set output file and run
    mddat%trajectoryfile = 'xtb.trj.xyz'
    call dynamics(mol,mddat,calc,.true.,io)

    if (io .ne. 0) then
      write (stdout,*) 'WARNING: MD run terminated ABNORMALLY'
      call creststop(status_failed)
    end if

    call rename(mddat%trajectoryfile,'crest_rotamers_0.xyz')

!> --------------------------------------------------------------------
!> xtb-syscall version
!> --------------------------------------------------------------------
  else

    !--- Writing constraining file xcontrol
    !--- Providing xcontrol overwrites constraints in coord file

    open (newunit=ich,file='xcontrol')
    if (env%cts%NCI) then
      do i = 1,10
        if (trim(env%cts%pots(i)) .ne. '') then
          write (ich,'(a)') trim(env%cts%pots(i))
        end if
      end do
    end if

    if (.not.env%solv_md) then
      write (ich,'(a)') '$constrain'
      write (ich,'(2x,a,i0)') 'atoms: 1-',solu%nat
      write (ich,'(2x,a)') 'force constant=0.5'
      write (ich,'(2x,a,a)') 'reference=ref.coord'
    end if

    write (ich,'(a)') '$md'
    write (ich,'(2x,a,f10.2)') 'hmass=',newhmass
    write (ich,'(2x,a,f10.2)') 'time=',newmdtime
    write (ich,'(2x,a,f10.2)') 'temp=',newtemp
    write (ich,'(2x,a,f10.2)') 'step=',newmdstep
    write (ich,'(2x,a,i0)') 'shake=',env%shake
    write (ich,'(2x,a,i0)') 'dump=',env%mddumpxyz
    write (ich,'(2x,a)') 'dumpxyz=500.0'

    if (env%ensemble_method .EQ. 2) then
      write (ich,'(a)') '$metadyn'
      write (ich,'(2x,a,i0,a,i0)') 'atoms: ',solu%nat+1,'-',clus%nat
      write (ich,'(2x,a,f10.2)') 'save=',newmetadlist
      write (ich,'(2x,a,f10.2)') 'kpush=',newmetadfac
      write (ich,'(2x,a,f10.2)') 'alp=',newmetadexp
    end if

    if (env%cts%cbonds_md) call write_cts_CBONDS(ich,env%cts)

    close (ich)

!--- Writing jobcall
    write (jobcall,'(a,1x,a,1x,a,'' --md --input xcontrol '',a,1x,a,a)') &
          &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),pipe
!--- slightly different jobcall for QMDFF usage
    if (env%useqmdff) then
      write (jobcall,'(a,1x,a,1x,a,'' --md --input xcontrol --qmdff'',a,1x,a,a)') &
         &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),pipe
    end if
!--- MD
    if (env%ensemble_method .EQ. 1) then
      call normalMD(fname,env,1,newtemp,newmdtime)
      write (stdout,*) 'Starting MD with the settings:'
      write (stdout,'(''     MD time /ps        :'',f8.1)') newmdtime
      write (stdout,'(''     MD Temperature /K  :'',f8.1)') newtemp
      write (stdout,'(''     dt /fs             :'',f8.1)') newmdstep
      write (tmppath,'(a,i0)') 'NORMMD1'

      r = makedir(tmppath)
      call copysub('xcontrol',tmppath)
      call chdirdbug(tmppath)
      call copy('coord','ref.coord')
      call chdirdbug(tmppath2)

      call command('cd '//trim(tmppath)//' && '//trim(jobcall),io)

      inquire (file=trim(tmppath)//'/'//'xtb.trj',exist=ex)
      if (.not.ex.or.io .ne. 0) then
        write (stdout,'(a,i0,a)') '*Warning: MD seemingly failed (no xtb.trj)*'
      else
        write (stdout,*) '*MD finished*'
      end if

      if (env%trackorigin) then
        call set_trj_origins('NORMMD','md')
      end if
      call chdirdbug('NORMMD1')
    end if

!--- MTD

    if (env%ensemble_method .EQ. 2) then
      call MetaMD(env,1,newmdtime,env%metadfac(1),env%metadexp(1), &
         &               env%metadlist(1))
      write (stdout,'(a,i4,a)') 'Starting Meta-MD with the settings:'
      write (stdout,'(''     MTD time /ps       :'',f8.1)') newmdtime
      write (stdout,'(''     dt /fs             :'',f8.1)') newmdstep
      write (stdout,'(''     MTD Temperature /K  :'',f8.1)') newtemp
      write (stdout,'(''     dumpstep(trj) /fs  :'',i8)') env%mddumpxyz
      write (stdout,'(''     Vbias factor k /Eh :'',f8.4)') newmetadfac
      write (stdout,'(''     Vbias exp α /bohr⁻²:'',f8.2)') newmetadexp

      write (tmppath,'(a,i0)') 'METADYN1'
      r = makedir(tmppath)
      call copysub('xcontrol',tmppath)
      call chdirdbug(tmppath)
      call copy('coord','ref.coord')

      call chdirdbug(tmppath2)

      call command('cd '//trim(tmppath)//' && '//trim(jobcall),io)

      inquire (file=trim(tmppath)//'/'//'xtb.trj',exist=ex)
      if (.not.ex.or.io .ne. 0) then
        write (stdout,'(a,i0,a)') '*Warning: Meta-MTD seemingly failed (no xtb.trj)*'
      else
        write (stdout,*) '*MTD finished*'
      end if

      if (env%trackorigin) then
        call set_trj_origins('METADYN','mtd')
      end if

      call chdirdbug('METADYN1')

    end if

    call rename('xtb.trj','crest_rotamers_0.xyz')
    call copysub('crest_rotamers_0.xyz',tmppath2)
    call dum%open('crest_rotamers_0.xyz')

!--- M(T)D stability check
    call minigrep('xtb.out','M(T)D is unstable, emergency exit',mdfail)
    if (dum%nall .eq. 1) then
      call copysub('xtb.out',resultspath)
      write (stdout,*) 'ERROR : M(T)D results only in one structure'
      if (mdfail) then
        write (stdout,*) '        It was unstable'
      else
        write (stdout,*) '        The M(T)D time step might be too large or the M(T)D time too short.'
      end if
      call copysub('xtb.out',resultspath)
      error stop '         Please check the xtb.out file in the ensemble folder'
    end if
    if (mdfail) then
      write (stdout,*)
      write (stdout,*) '   WARNING: The M(T)D was unstable.'
      write (stdout,*) '            Please check the xtb.out file in the ensemble folder.'
      write (stdout,*)
      call copysub('xtb.out',resultspath)
    end if
    call dum%deallocate
    call chdirdbug(tmppath2)
    call clus%write('coord')
    call inputcoords(env,'coord') !Necessary

!--- Optimization
    call print_qcg_opt
    call multilevel_opt(env,99)
  end if

end subroutine xtb_md_ensemble_qcg

!___________________________________________________________________________________
!
! An xTB docking on all available threads
!___________________________________________________________________________________

subroutine xtb_dock(env,fnameA,fnameB,solu,clus)
  use crest_parameters
  use iomod
  use crest_data
  use qcg_coord_type
  implicit none

  type(systemdata)               :: env
  type(coord_qcg),intent(in)     :: solu,clus
  character(len=*),intent(in)    :: fnameA,fnameB
  character(len=80)              :: pipe
  character(len=512)             :: jobcall
  integer                        :: i,ich,T,Tn

  call remove('xtb_dock.out')
  call remove('xcontrol')

  pipe = ' 2>/dev/null'

!---- writing wall pot in xcontrol
  call write_wall(env,solu%nat,solu%ell_abc,clus%ell_abc,'xcontrol')

!---- Write directed stuff, if requested
  if (allocated(env%directed_file)) then
    do i = 1,size(env%directed_number)
      if &
      & ((i == 1.and.env%directed_number(i) >= clus%nmol).OR. &
      & (env%directed_number(i) >= clus%nmol.and.env%directed_number(i-1) < clus%nmol)) &
      & then
        open (newunit=ich,file='xcontrol',status='old',position='append',action='write')
        write (ich,'("$directed")')
        write (ich,'(a,1x,a)') 'atoms:',trim(env%directed_list(i,1))
        write (ich,'("$end")')
      end if
    end do
  end if

!--- Setting threads
  call new_ompautoset(env,'auto',1,T,Tn)

!--- Jobcall docking
  write (jobcall,'(a,1x,''dock'',1x,a,1x,a,1x,a,1x,f4.2,1x,''--nfrag1'',1x,i0,1x,a,1x, &
          & ''--input xcontrol > xtb_dock.out'',a)') &
          &     trim(env%ProgName),trim(fnameA),trim(fnameB),trim(env%gfnver),&
          &     env%optlev,solu%nat,trim(env%docking_qcg_flag),trim(pipe)
  call command(trim(jobcall))

! cleanup
  call remove('wbo')
  call remove('charges')
  call remove('xtbrestart')

end subroutine xtb_dock

!___________________________________________________________________________________
!
! An xTB optimization on all available threads
!___________________________________________________________________________________

subroutine opt_cluster(env,solu,clus,fname,without_pot)
  use crest_parameters
  use iomod
  use crest_data
  use qcg_coord_type

  implicit none

  type(systemdata)              :: env
  type(coord_qcg),intent(in)    :: solu,clus
  character(len=*),intent(in)   :: fname
  logical,optional,intent(in)   :: without_pot
  character(len=*),parameter    :: pipe = ' 2>/dev/null'
  character(len=:),allocatable  :: jobcall
  integer :: T,Tn

  call remove('xtb.out')

!---- writing wall pot in xcontrol
  if (.not.without_pot) then
    call write_wall(env,solu%nat,solu%ell_abc,clus%ell_abc,'xcontrol')
  end if

!--- Setting threads
  call new_ompautoset(env,'subprocess',1,T,Tn)

!--- Jobcall optimization
  jobcall = trim(env%ProgName)//' '//trim(fname)//' --opt '//optlevflag(env%optlev)
  jobcall = trim(jobcall)//' '//trim(env%gfnver)
  if (without_pot) then
    jobcall = trim(jobcall)//' '//trim(env%solv)
  end if
  jobcall = trim(jobcall)//' > xtb_opt.out 2>/dev/null'
  call command(trim(jobcall))

! cleanup
  call remove('wbo')
  call remove('charges')
  call remove('xtbrestart')

!--- Jobcall SP for gbsa model
  if (.not.without_pot) then
    jobcall = trim(env%ProgName)//' xtbopt.coord --sp '//trim(env%gfnver)
    jobcall = trim(jobcall)//' '//trim(env%solv)
    jobcall = trim(jobcall)//' > xtb_sp.out 2>/dev/null'
    call command(trim(jobcall))
  end if

! cleanup
  call remove('wbo')
  call remove('charges')
  call remove('xtbrestart')

end subroutine opt_cluster

!___________________________________________________________________________________
!
! xTB docking calculation performed in parallel
!___________________________________________________________________________________

subroutine ensemble_dock(env,outer_ell_abc,nfrag1,frag1_file,frag2_file,n_shell&
        &,n_solvent,NTMP,TMPdir,conv)
  !*****************************************************************
  !* Dock one solvent molecule (frag2_file) to the cluster
  !* (frag1_file) in each open directory TMPdir<conv(i)>, i=1..NTMP,
  !* via parallel "xtb dock" system calls.
  !*****************************************************************
  use crest_parameters
  use iomod
  use crest_data
  use qcg_coord_type

  implicit none
  type(systemdata)                :: env

  character(len=*),intent(in)     :: TMPdir     !directory name
  integer,intent(in)              :: NTMP       !number of structures to be optimized
  integer,intent(in)              :: nfrag1     !#atoms of larger fragment
  integer,intent(in)              :: conv(env%nqcgclust+1)
  real(wp),intent(in)             :: outer_ell_abc(env%nqcgclust,3)
  integer,intent(in)              :: n_shell,n_solvent

  integer                         :: i,k
  integer                         :: vz,T,Tn
  character(len=20)               :: pipe
  character(len=1024)             :: jobcall
  character(len=512)              :: thispath,tmppath
  character(len=1536),allocatable :: jobs(:)
  character(len=*),intent(in)    :: frag1_file
  character(len=*),intent(in)    :: frag2_file
  character(len=64)               :: frag1
  character(len=64)               :: frag2
  real(wp)                        :: percent
  character(len=2)                :: flag
  integer                         :: funit

! some options
  pipe = '2>/dev/null'
  frag1 = 'solvent_cluster.coord'
  frag2 = 'solvent'
  call getcwd(thispath)

! setting the threads for correct parallelization
  call new_ompautoset(env,'auto',NTMP,T,Tn)

  write (jobcall,'(a,1x,''dock'',1x,a,1x,a,1x,a,1x,f4.2,1x,''--nfrag1'',1x,i0,1x,&
          & ''--input xcontrol --fast > xtb_dock.out '',a)') &
          & trim(env%ProgName),trim(frag1_file),trim(frag2_file),&
          & trim(env%gfnver),env%optlev,nfrag1,trim(pipe)

  flag = '$'
  do i = 1,NTMP
    vz = i
    write (tmppath,'(a,i0)') trim(TMPdir),conv(i)
    call chdirdbug(trim(tmppath))
    open (newunit=funit,file='xcontrol')
    write (funit,'(a,"fix")') trim(flag)
    write (funit,'(3x,"atoms: 1-",i0)') n_shell !Initial number of atoms (starting solvent shell)
    write (funit,'(a,"wall")') trim(flag)
    write (funit,'(3x,"potential=polynomial")')
    write (funit,'(3x,"ellipsoid:",1x,3(g0,",",1x),i0,"-",i0)') outer_ell_abc(conv(vz),:), &
            & n_shell+1,n_shell+n_solvent !Initial number of atoms (starting solvent shell)
    close (funit)
    call chdirdbug(trim(thispath))
  end do

  !> full job strings are set up here and not inside the tasks below: with ifx
  !> (2025.3) the internal write on the dummy arguments TMPdir/conv within an
  !> OpenMP task yields corrupted directory names or a segfault
  allocate (jobs(NTMP))
  do i = 1,NTMP
    write (tmppath,'(a,i0)') trim(TMPdir),conv(i)
    jobs(i) = 'cd '//trim(tmppath)//' && '//trim(jobcall)
  end do

  k = 0 !counting the finished jobs

!___________________________________________________________________________________

!$omp parallel &
!$omp shared( vz,NTMP,percent,k,jobs )
!$omp single
  do i = 1,NTMP
    vz = i
    !$omp task firstprivate( vz )
    call command(trim(jobs(vz)))
    !$omp critical
    k = k+1
    percent = real(k)/real(NTMP)*100
    !$omp end critical
    !$omp end task
  end do

!$omp taskwait
!$omp end single
!$omp end parallel

!___________________________________________________________________________________
  deallocate (jobs)
  call chdirdbug(trim(thispath))

end subroutine ensemble_dock

!___________________________________________________________________________________
!
! xTB CFF optimization performed in parallel
!___________________________________________________________________________________

subroutine cff_opt(pr,env,fname,n12,NTMP,TMPdir,conv,nothing_added,eread)
  !*****************************************************************
  !* Constrained optimization + SP of the CFF clusters in TMPdir<n>.
  !* Redirects to cff_opt_calculator unless env%legacy is set.
  !*
  !* pr            : final run, all env%nqcgclust directories are
  !*                 processed (NTMP and conv are reset accordingly)
  !* n12           : atoms 1-n12 are frozen (0 = no freezing)
  !* NTMP          : number of open jobs
  !* conv          : directory indices of the open jobs
  !* nothing_added : per cluster, no freezing in the final run
  !* eread         : energies, indexed by cluster (directory) number
  !*****************************************************************
  use crest_parameters
  use iomod
  use crest_data
  use strucrd
  use term_ui,only:progress_init,progress_update,progress_finish
  implicit none

  type(systemdata)                :: env
  character(len=*),intent(in)     :: fname      !file base name
  character(len=*),intent(in)     :: TMPdir     !directory name
  integer,intent(inout)           :: NTMP       !number of structures to be optimized
  integer,intent(inout)           :: conv(env%nqcgclust+1)
  logical,intent(in)              :: pr
  logical,intent(in)              :: nothing_added(env%nqcgclust)
  integer,intent(in)              :: n12
  real(wp),intent(inout)          :: eread(env%nqcgclust)
  integer                         :: i,k
  integer                         :: vz,T,Tn
  integer                         :: funit
  character(len=512)              :: thispath,tmppath
  character(len=1024)             :: jobcall
  character(len=2)                :: flag
  real(wp)                        :: percent
  logical :: ex,e_there
  character(len=*),parameter      :: pipe = '2>/dev/null'

  !> redirect to calculator version in new implementation
  if (.not.env%legacy) then
    call cff_opt_calculator(pr,env,fname,n12,NTMP,TMPdir,conv,nothing_added,eread)
    return
  end if

! pr eq true => post opt run, which has to be performed in every directory !!!
  if (pr) then
    k = 0
    NTMP = env%nqcgclust
    do i = 1,env%nqcgclust
      k = k+1
      conv(k) = i
      conv(env%nqcgclust+1) = k
    end do
  end if

! setting the threads for correct parallelization
  call new_ompautoset(env,'auto',NTMP,T,Tn)

  if (pr) then
    write (stdout,'(2x,''Starting optimizations + SP  of structures'')')
    write (stdout,'(2x,i0,'' jobs to do.'')') NTMP
  end if

  call getcwd(thispath)
  do i = 1,NTMP
    write (tmppath,'(a,i0)') trim(TMPdir),conv(i)
    call chdirdbug(trim(tmppath))
    open (newunit=funit,file='xcontrol')
    if (n12 .ne. 0) then
      flag = '$'
      write (funit,'(a,"fix")') trim(flag)
      write (funit,'(3x,"atoms: 1-",i0)') n12 !Initial number of atoms (starting solvent shell)
    end if
    close (funit)
    if (pr.and.nothing_added(conv(i))) call remove('xcontrol')
    call chdirdbug(trim(thispath))
  end do

!--- Jobcall
  write (jobcall,'(a,1x,a,1x,a,'' --input xcontrol --opt '',i0,1x,a,'' >xtb.out'')') &
  &    trim(env%ProgName),trim(fname),trim(env%gfnver),nint(env%optlev),trim(pipe)

  if (NTMP .lt. 1) then
    write (stdout,'(2x,"No structures to be optimized")')
    return
  end if

  k = 0 !counting the finished jobs
  if (pr) then
    call progress_init(env%ps,NTMP,width=50,prefix=" ↳ ", &
      &                suffix="",show_time=.true.,show_eta=.false.)
    call progress_update(env%ps,0,NTMP)
  end if
!___________________________________________________________________________________

!$omp parallel &
!$omp shared( env, vz,jobcall,NTMP,percent,k,TMPdir,conv )
!$omp single
  do i = 1,NTMP
    vz = i
    !$omp task firstprivate( vz ) private( tmppath )
    write (tmppath,'(a,i0)') trim(TMPdir),conv(vz)
    call command('cd '//trim(tmppath)//' && '//trim(jobcall))
    !$omp critical
    k = k+1
    percent = real(k)/real(NTMP)*100
    if (pr) call progress_update(env%ps,k,NTMP)
    !$omp end critical
    !$omp end task
  end do
!$omp taskwait
!$omp end single
!$omp end parallel

!__________________________________________________________________________________

  do i = 1,NTMP
    write (tmppath,'(a,i0)') trim(TMPdir),conv(i)
    call chdirdbug(trim(tmppath))
    call remove('xtbrestart')
    call chdirdbug(trim(thispath))
  end do

  !create the system call for sp (needed for gbsa model)
  write (jobcall,'(a,1x,a,1x,a,'' --sp '',a,1x,a,'' >xtb_sp.out'')') &
  &    trim(env%ProgName),'xtbopt.coord',trim(env%gfnver),trim(env%solv),trim(pipe)

  if (NTMP .lt. 1) then
    write (stdout,'(2x,"Nothing to do")')
    return
  end if

  k = 0 !counting the finished jobs
  if (pr) then
    call progress_init(env%ps,NTMP,width=50,prefix=" ↳ ", &
      &                suffix="",show_time=.true.,show_eta=.false.)
    call progress_update(env%ps,0,NTMP)
  end if
!___________________________________________________________________________________

!$omp parallel &
!$omp shared( vz,jobcall,NTMP,percent,k,TMPdir,conv )
!$omp single
  do i = 1,NTMP
    vz = i
    !$omp task firstprivate( vz ) private( tmppath )
    write (tmppath,'(a,i0)') trim(TMPdir),conv(vz)
    call command('cd '//trim(tmppath)//' && '//trim(jobcall))
    !$omp critical
    k = k+1
    percent = real(k)/real(NTMP)*100
    if (pr) call progress_update(env%ps,k,NTMP)
    !$omp end critical
    !$omp end task
  end do
!$omp taskwait
!$omp end single
!$omp end parallel

!___________________________________________________________________________________

  do i = 1,NTMP
    write (tmppath,'(a,i0)') trim(TMPdir),conv(i)
    call chdirdbug(trim(tmppath))

    call grepval('xtb_sp.out','| TOTAL ENERGY',e_there,eread(conv(i)))

    call remove('xtbrestart')
    !call remove('xcontrol')
    call chdirdbug(trim(thispath))
  end do

  if (pr) call progress_finish(env%ps)

end subroutine cff_opt

subroutine cff_opt_calculator(pr,env,fname,n12,NTMP,TMPdir, &
    & conv,nothing_added,eread)
  !*****************************************************************
  !* Calculator version of cff_opt: gas-phase optimization with
  !* atoms 1-n12 frozen, followed by an SP with implicit solvation.
  !* Arguments as in cff_opt.
  !*****************************************************************
  use crest_parameters
  use iomod
  use crest_data
  use strucrd
  use crest_calculator
  use optimize_module
  use term_ui,only:progress_init,progress_update,progress_finish
  implicit none

  type(systemdata)                :: env
  character(len=*),intent(in)     :: fname      !file base name
  character(len=*),intent(in)     :: TMPdir     !directory name
  integer,intent(inout)           :: NTMP       !number of structures to be optimized
  integer,intent(inout)           :: conv(env%nqcgclust+1)
  logical,intent(in)              :: pr
  logical,intent(in)              :: nothing_added(env%nqcgclust)
  integer,intent(in)              :: n12
  real(wp),intent(inout)          :: eread(env%nqcgclust)
  integer                         :: i,k
  integer                         :: vz,T,Tn
  integer                         :: funit
  character(len=512)              :: thispath,tmppath
  character(len=1024)             :: jobcall
  character(len=2)                :: flag
  real(wp)                        :: percent
  logical :: gbsa_tmp,opt_nofreeze,all_frozen
  integer :: io,nconstraints_tmp
  character(len=40)  :: solv_tmp
  real(wp) :: etot
  type(calcdata),allocatable :: newcalcs(:)
  type(calculation_settings) :: clevel
  type(coord) :: mol,molopt
  type(coord),allocatable :: structures(:)
  real(wp),allocatable :: grd(:,:)
  character(len=*),parameter      :: pipe = '2>/dev/null'

! pr eq true => post opt run, which has to be performed in every directory !!!
  if (pr) then
    k = 0
    NTMP = env%nqcgclust
    do i = 1,env%nqcgclust
      k = k+1
      conv(k) = i
      conv(env%nqcgclust+1) = k
    end do
  end if

  !> setting the threads to accelerate individual energy calculations
  call new_ompautoset(env,'max',NTMP,T,Tn)

  if (pr) then
    write (stdout,'(2x,"Starting optimizations + SP  of structures")')
    write (stdout,'(2x,i0,'' jobs to do.'')') NTMP
  end if
  if (NTMP .lt. 1) then
    write (stdout,'(2x,"No structures to be optimized")')
    return
  end if

  !> Local storage for structures.
  allocate (structures(NTMP))

  call getcwd(thispath)
  do i = 1,NTMP
    write (tmppath,'(a,i0)') trim(TMPdir),conv(i)
    call chdirdbug(trim(tmppath))
    call structures(i)%open(fname)
    call chdirdbug(trim(thispath))
  end do

!--- Jobcall WITHOUT GBSA, back up data
  solv_tmp = env%solv
  gbsa_tmp = env%gbsa
  env%solv = ''
  env%gbsa = .false.

  !> keep everything nice and separate, we have different molecules after all...
  !> also, this will be gas-phase
  call clevel%create(env%gfnver,chrg=env%chrg,uhf=env%uhf)
  allocate (newcalcs(NTMP))
  do i = 1,NTMP
    call newcalcs(i)%add(clevel)
    !> other important settings from env
    newcalcs(i)%optlev = int(env%optlev)
  end do

  !> -------------------------------------------------
  !> OPTIMIZATIONS
  !> -------------------------------------------------
  !> The structures may have different numbers
  !> and it should not not be many structures anyways.
  !> therefore, just optimize them serially.

  k = 0 !counting the finished jobs
  if (pr) then
    call progress_init(env%ps,NTMP,width=50,prefix=" ↳ ", &
      &                suffix="",show_time=.true.,show_eta=.false.)
    call progress_update(env%ps,0,NTMP)
  end if
!___________________________________________________________________________________

  do i = 1,NTMP
    mol = structures(i)
    opt_nofreeze = (pr.and.nothing_added(conv(i)))
    !> safeguard against freezing all, the SP below provides the energy
    all_frozen = (.not.opt_nofreeze.and.n12 >= mol%nat)

    if (.not.all_frozen) then
      if (.not.opt_nofreeze.and.n12 > 0) then
        call newcalcs(i)%set_freeze(mol%nat,1,n12)
      end if
      allocate (grd(3,mol%nat))
      call optimize_geometry(mol,molopt,newcalcs(i),etot,grd, &
        & .false.,.false.,io)
      structures(i) = molopt
      structures(i)%energy = etot
      deallocate (grd)
    end if

    k = k+1
    percent = real(k)/real(NTMP)*100.0_wp
    if (pr) call progress_update(env%ps,k,NTMP)
  end do
  if (pr) call progress_finish(env%ps)
  !> clear up space
  deallocate (newcalcs)
!__________________________________________________________________________________

!> -------------------------------------------------
!> SINGLEPOINTS
!> -------------------------------------------------
!> same as for optimizations, turn on impl. solv again
  env%solv = solv_tmp
  env%gbsa = gbsa_tmp

  !> again, keep everything nice and separate (now with impl. solv)
  call clevel%create(env%gfnver,chrg=env%chrg,uhf=env%uhf,&
    &  solvmodel=env%solv,solvent=env%solvent &
  )
  allocate (newcalcs(NTMP))
  do i = 1,NTMP
    call newcalcs(i)%add(clevel)
  end do

  k = 0 !counting the finished jobs
  if (pr) then
    call progress_init(env%ps,NTMP,width=50,prefix=" ↳ ", &
      &                suffix="",show_time=.true.,show_eta=.false.)
    call progress_update(env%ps,0,NTMP)
  end if
!___________________________________________________________________________________

  do i = 1,NTMP
    mol = structures(i)
    allocate (grd(3,mol%nat))

    call engrad(mol,newcalcs(i),etot,grd,io)
    structures(i)%energy = etot

    k = k+1
    percent = real(k)/real(NTMP)*100.0_wp
    if (pr) call progress_update(env%ps,k,NTMP)

    deallocate (grd)
  end do
  if (pr) call progress_finish(env%ps)
!___________________________________________________________________________________

  !> for compatibility reasons, let's write the optimized geometries
  !> and pass on energies
  call getcwd(thispath)
  do i = 1,NTMP
    eread(conv(i)) = structures(i)%energy
    write (tmppath,'(a,i0)') trim(TMPdir),conv(i)
    call chdirdbug(trim(tmppath))
    call structures(i)%write('xtbopt.coord')
    call chdirdbug(trim(thispath))
  end do

  if (allocated(newcalcs)) deallocate (newcalcs)
  if (allocated(structures)) deallocate (structures)
end subroutine cff_opt_calculator

!___________________________________________________________________________________
!
! xTB SP performed in parallel
!___________________________________________________________________________________
subroutine ens_sp_with_io(env,ens,clus,resultspath)
  use crest_parameters
  use crest_data
  use qcg_coord_type
  use strucrd
  use iomod
  implicit none
  type(systemdata),intent(inout) :: env
  type(ensemble),intent(inout)   :: ens
  type(coord_qcg),intent(inout)  :: clus
  character(len=512),intent(in)  :: resultspath
  character(len=512) :: tmppath2,to,comment
  integer :: i,io,minpos
  logical :: e_there,not_param,ex

  call getcwd(tmppath2)

  !--- Write folder with xyz-coordinates
  do i = 1,ens%nall
    call rdxmolselec('ensemble.xyz',i,clus%nat,clus%at,clus%xyz)
    write (to,'("TMPSP",i0)') i
    io = makedir(trim(to))
    call copysub('.UHF',to)
    call copysub('.CHRG',to)
    call chdirdbug(to)
    call wrxyz('cluster.xyz',clus%nat,clus%at,clus%xyz*bohr)
    call chdirdbug(tmppath2)
  end do
  !--- SP
  write (stdout,*)
  call ens_sp(env,'cluster.xyz',ens%nall,'TMPSP')
  !--- Getting energy
  do i = 1,ens%nall
    call rdxmolselec('ensemble.xyz',i,clus%nat,clus%at,clus%xyz)
    write (to,'("TMPSP",i0)') i
    call chdirdbug(to)
    call grepval('xtb_sp.out','| TOTAL ENERGY',e_there,ens%er(i))
    call chdirdbug(tmppath2)
  end do

  if (.not.e_there) then
    write (stdout,*)
    write (stdout,*) 'Energy not found. Error in xTB computations occured'
    call chdirdbug(to)
    call minigrep('xtb_sp.out','solv_model_loadInternalParam',not_param)
    call chdirdbug(tmppath2)
    if (not_param) then
      write (stdout,*) '  !!!WARNIG: CHOSEN SOLVENT NOT PARAMETERIZED &
      & FOR IMPLICIT SOLVATION MODEL!!!'
      write (stdout,'(''  CHECK IF '',A,'' IS AVAILABLE IN xTB'')') env%solv
      write (stdout,*) '  PLEASE RESTART THE ENSEMBLE GENERATION WITH AVAILABLE&
             &  PARAMETERIZATION IF YOU NEED ENERGIES'
      call copysub('crest_conformers.xyz',resultspath)
      write (stdout,*) '  The enesemble can be found in the <ensemble> directory&
             & as <crest_conformers.xyz>'
      error stop
    end if
  end if
  call ens%write('full_ensemble.xyz')

!--- crest_best structure
  minpos = minloc(ens%er,dim=1)
  write (to,'("TMPSP",i0)') minpos
  call chdirdbug(to)
  call rdxmol('cluster.xyz',clus%nat,clus%at,clus%xyz)
  call chdirdbug(tmppath2)
  write (comment,'(F20.8)') ens%er(minpos)
  inquire (file='crest_best.xyz',exist=ex)
  if (ex) then
    call rmrf('crest_best.xyz') !remove crest_best from
  end if
  call wrxyz('crest_best.xyz',clus%nat,clus%at,clus%xyz,trim(comment))

end subroutine ens_sp_with_io

subroutine ens_sp(env,fname,NTMP,TMPdir)
  use crest_parameters
  use iomod
  use crest_data
  use strucrd
  use term_ui,only:progress_init,progress_update,progress_finish
  implicit none

  type(systemdata)                :: env
  character(len=*),intent(in)     :: fname      !file base name
  character(len=*),intent(in)     :: TMPdir     !directory name
  integer,intent(inout)           :: NTMP       !number of structures to be optimized

  integer                         :: i,k
  integer                         :: vz,T,Tn
  character(len=20)               :: pipe
  character(len=512)              :: thispath,tmppath
  character(len=1024)             :: jobcall
  real(wp)                        :: percent

! setting the threads for correct parallelization
  call new_ompautoset(env,'auto',NTMP,T,Tn)

  write (stdout,'(2x,''---------------------------------------------'')')
  write (stdout,'(2x,''Single point computation with GBSA/ALPB model'')')
  write (stdout,'(2x,''---------------------------------------------'')')
  write (stdout,'(2x,i0,'' jobs to do.'')') NTMP

  pipe = '2>/dev/null'

  call getcwd(thispath)

  if (NTMP .lt. 1) then
    write (stdout,'(2x,"No structures to be optimized")')
    return
  end if

!--- Jobcall
  write (jobcall,'(a,1x,a,1x,a,'' --sp '',a,1x,a,'' > xtb_sp.out'')') &
  &    trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),trim(pipe)

  k = 0 !counting the finished jobs
  call progress_init(env%ps,NTMP,width=50,prefix=" ↳ ", &
    &                suffix="",show_time=.true.,show_eta=.false.)
  call progress_update(env%ps,0,NTMP)

!___________________________________________________________________________________

!$omp parallel &
!$omp shared( env,vz,NTMP,percent,k,TMPdir,jobcall )
!$omp single
  do i = 1,NTMP
    vz = i
    !$omp task firstprivate( vz ) private( tmppath )
    call initsignal()
    !$omp critical
    write (tmppath,'(a,i0)') trim(TMPdir),vz
    !$omp end critical
    call command('cd '//trim(tmppath)//' && '//trim(jobcall))
    !$omp critical
    k = k+1
    percent = real(k)/real(NTMP)*100
    call progress_update(env%ps,k,NTMP)
    !$omp end critical
    !$omp end task
  end do
!$omp taskwait
!$omp end single
!$omp end parallel

!__________________________________________________________________________________

  do i = 1,NTMP
    write (tmppath,'(a,i0)') trim(TMPdir),i
    call chdirdbug(trim(tmppath))
    call remove('xtbrestart')
    call chdirdbug(trim(thispath))
  end do
  call progress_finish(env%ps)

end subroutine ens_sp

!___________________________________________________________________________________
!
! xTB Freq compuatation performed in parallel
!___________________________________________________________________________________

subroutine ens_freq(env,fname,NTMP,TMPdir,opt)
  use crest_parameters
  use iomod
  use crest_data
  use strucrd
  use term_ui,only:progress_init,progress_update,progress_finish
  implicit none

  type(systemdata)                :: env
  character(len=*),intent(in)     :: fname      !file base name
  character(len=*),intent(in)     :: TMPdir     !directory name
  integer,intent(inout)           :: NTMP       !number of structures to be optimized

  integer                         :: i,k
  integer                         :: vz,T,Tn
  character(len=512)              :: thispath,tmppath
  character(len=1024)             :: jobcall
  real(wp)                        :: percent
  logical                         :: opt
  character(len=*),parameter      :: pipe = '2>/dev/null'

  !> redirect to new calculator version if available
  if (.not.env%legacy) then
    call ens_freq_calculator(env,fname,NTMP,TMPdir,opt)
    return
  end if

! setting the threads for correct parallelization
  call new_ompautoset(env,'auto',NTMP,T,Tn)

  write (stdout,'(2x,''Starting reoptimizations + Frequency computation of ensemble'')')
  write (stdout,'(2x,i0,'' jobs to do.'')') NTMP

  call getcwd(thispath)

  if (NTMP .lt. 1) then
    write (stdout,'(2x,"No structures to be optimized")')
    return
  end if

  k = 0 !counting the finished jobs
  call progress_init(env%ps,NTMP,width=50,prefix=" ↳ ", &
    &                suffix="",show_time=.true.,show_eta=.false.)
  call progress_update(env%ps,0,NTMP)

!--- Jobcall
  if (.not.opt) then
    write (jobcall,'(a,1x,a,1x,a,'' --hess '',a,'' >xtb_freq.out'')') &
     &    trim(env%ProgName),trim(fname),trim(env%gfnver),trim(pipe)
  else
    write (jobcall,'(a,1x,a,1x,a,'' --ohess '',a,'' >xtb_freq.out'')') &
    &    trim(env%ProgName),trim(fname),trim(env%gfnver),trim(pipe)
  end if

!___________________________________________________________________________________

!$omp parallel &
!$omp shared( env,vz,NTMP,percent,k,TMPdir,jobcall )
!$omp single
  do i = 1,NTMP
    vz = i
    !$omp task firstprivate( vz ) private( tmppath )
    write (tmppath,'(a,i0)') trim(TMPdir),i
    call command('cd '//trim(tmppath)//' && '//trim(jobcall))
    !$omp critical
    k = k+1
    percent = real(k)/real(NTMP)*100
    call progress_update(env%ps,k,NTMP)
    !$omp end critical
    !$omp end task
  end do
!$omp taskwait
!$omp end single
!$omp end parallel

!__________________________________________________________________________________

  do i = 1,NTMP
    write (tmppath,'(a,i0)') trim(TMPdir),i
    call chdirdbug(trim(tmppath))
    call remove('xtbrestart')
    call chdirdbug(trim(thispath))
  end do
  call progress_finish(env%ps)

end subroutine ens_freq

subroutine ens_freq_calculator(env,fname,NTMP,TMPdir,opt)
  use crest_parameters
  use iomod
  use crest_data
  use strucrd
  use crest_calculator
  use thermochem_module
  use optimize_module
  use term_ui,only:progress_init,progress_update,progress_finish
  implicit none

  type(systemdata)                :: env
  character(len=*),intent(in)     :: fname      !file base name
  character(len=*),intent(in)     :: TMPdir     !directory name
  integer,intent(inout)           :: NTMP       !number of structures to be optimized

  integer                         :: i,k
  integer                         :: vz,T,Tn
  character(len=*),parameter      :: pipe = '2>/dev/null'
  character(len=512)              :: thispath,tmppath
  character(len=1024)             :: jobcall
  real(wp)                        :: percent
  logical                         :: opt

  !> local calculation setup
  type(coord) :: tmpmol,mol
  type(calculation_settings) :: clevel
  type(calcdata),allocatable :: newcalcs(:)
  real(wp),allocatable :: tmpgrd(:,:),hess(:,:),freq(:)
  real(wp) :: etmp
  integer :: n3,io,ich

  real(wp) :: ithr,fscal,sthr
  character(len=:),allocatable :: emodel
  integer :: nt,nfreq,nrt
  real(wp),allocatable :: temps(:),et(:),ht(:),stot(:),gt(:)
  real(wp) :: zpve

! setting the threads for correct parallelization
  call new_ompautoset(env,'max',NTMP,T,Tn)

  write (stdout,'(2x,''Starting reoptimizations + Frequency computation of ensemble'')')
  write (stdout,'(2x,i0,'' jobs to do.'')') NTMP

  call getcwd(thispath)

  if (NTMP .lt. 1) then
    write (stdout,'(2x,"No structures to be optimized")')
    return
  end if

  k = 0 !counting the finished jobs
  call progress_init(env%ps,NTMP,width=50,prefix=" ↳ ", &
    &                suffix="",show_time=.true.,show_eta=.false.)
  call progress_update(env%ps,0,NTMP)

!--- Jobcall
  if (.not.opt) then
    write (jobcall,'(a,1x,a,1x,a,'' --hess '',a,'' >xtb_freq.out'')') &
     &    trim(env%ProgName),trim(fname),trim(env%gfnver),trim(pipe)
  else
    write (jobcall,'(a,1x,a,1x,a,'' --ohess '',a,'' >xtb_freq.out'')') &
    &    trim(env%ProgName),trim(fname),trim(env%gfnver),trim(pipe)
  end if

!--- prepare calcs
  call clevel%create(env%gfnver,chrg=env%chrg,uhf=env%uhf)
  allocate (newcalcs(NTMP))
  do i = 1,NTMP
    call newcalcs(i)%add(clevel)
    newcalcs(i)%optlev = int(env%optlev)
    !call newcalcs(i)%info(stdout)
  end do

!--- prepare thermo
  !> inversion threshold
  ithr = env%thermo%ithr
  !> frequency scaling factor
  fscal = env%thermo%fscal
  !> RR-HO interpolation
  sthr = env%thermo%sthr
  !> Svib model
  emodel = env%thermo%emodel

  !> we just need one temperature
  nt = 1
  allocate (temps(nt),et(nt),ht(nt),gt(nt),stot(nt),source=0.0_wp)
  temps(:) = 298.15_wp

!___________________________________________________________________________________
!> serial runtype
  do i = 1,NTMP
    write (tmppath,'(a,i0)') trim(TMPdir),i
    call chdirdbug(trim(tmppath))
    call tmpmol%open(fname)

    allocate (tmpgrd(3,tmpmol%nat),source=0.0_wp)
    if (opt) then
      call optimize_geometry(tmpmol,mol,newcalcs(i),etmp,tmpgrd, &
      & .false.,.false.,io)
    else
      mol = tmpmol
    end if

    n3 = mol%nat*3
    allocate (hess(n3,n3),source=0.0_wp)
    allocate (freq(n3),source=0.0_wp)

    !>-- compute Hessian
    call numhess2(mol%nat,mol%at,mol%xyz,newcalcs(i),hess,io)

    !>-- Projects and mass-weights the Hessian
    call prj_mw_hess(mol%nat,mol%at,n3,mol%xyz,hess(:,:))

    !>-- Computes the Frequencies
    call frequencies(mol%nat,mol%at,mol%xyz,n3,hess(:,:),freq(:),io)

    !> write dummy "xtb_freq.out"
    open (newunit=ich,file="xtb_freq.out")
    !> calcthermo wants input in Bohr
    call calcthermo(mol%nat,mol%at,mol%xyz,freq,.true., &
    & ithr,fscal,sthr,nt,temps,et,ht,gt,stot,ich,emodel=emodel)
    close (ich)

    deallocate (freq,hess,tmpgrd)
    k = k+1
    call progress_update(env%ps,k,NTMP)
    call chdirdbug(trim(thispath))
  end do

!__________________________________________________________________________________

  call progress_finish(env%ps)
end subroutine ens_freq_calculator

!============================================================!
! subroutine wr_cluster_cut
! Cuts a cluster file and and writes the parts
!
! On Input: fname          - name of the coord file
!           n1             - number of atoms fragment1
!           n2             - number of atmos fragment2
!           iter           - number of solvent molecules
!           fname_solu_cut - name of outputfile fragment1
!           fname_solv_cut - name of outputfile fragment2
!
!============================================================!

subroutine wr_cluster_cut(fname_cluster,n1,n2,iter,fname_solu_cut,fname_solv_cut)
  use crest_parameters
  use strucrd
  implicit none
  integer,intent(in)         :: n1,n2,iter
  real(wp)                    :: xyz1(3,n1)
  real(wp)                    :: xyz2(3,n2*iter)
  integer                     :: at1(n1),at2(n2*iter)
  character(len=*),intent(in) :: fname_cluster,fname_solu_cut,fname_solv_cut
  character(len=256)         :: atmp
  character(len=2)           :: a2
  integer                     :: ich,i,k,stat,io,io2

  !ich = 142
  open (newunit=ich,file=fname_cluster,iostat=stat)
  read (ich,'(a)') atmp
  k = 1
  do i = 1,n1
    read (ich,'(a)',iostat=io) atmp
    if (io < 0) exit
    atmp = adjustl(atmp)
    call coordline(atmp,a2,xyz1(1:3,k),io2)
    at1(k) = e2i(a2)
    k = k+1
  end do
  k = 1
  do i = 1,n2*iter
    read (ich,'(a)',iostat=io) atmp
    if (io < 0) exit
    atmp = adjustl(atmp)
    call coordline(atmp,a2,xyz2(1:3,k),io2)
    at2(k) = e2i(a2)
    k = k+1
  end do

  call wrc0(fname_solu_cut,n1,at1,xyz1)
  call wrc0(fname_solv_cut,n2*iter,at2,xyz2)
  close (ich)

end subroutine wr_cluster_cut

!----------------------------------------------------------------------------
! write a wall potential in a file used as xtb input

subroutine write_wall(env,n1,rabc1,rabc12,fname)
  use crest_parameters
  use crest_data

  implicit none

  type(systemdata)     :: env
  integer,intent(in)  :: n1
  real(wp),intent(in)  :: rabc1(3),rabc12(3)
  character(len=8)    :: flag
  character(len=*)    :: fname
  integer :: funit

  open (newunit=funit,file=fname)
  flag = '$'
  write (funit,'(a,"wall")') trim(flag)
  write (funit,'(3x,"potential=polynomial")')
  write (funit,'(3x,"ellipsoid:",1x,3(g0,",",1x),"all")') rabc12
  write (funit,'(3x,"ellipsoid:",1x,3(g0,",",1x),"1-",i0)') rabc1,n1
  if (env%constrain_solu) then
    write (funit,'("$fix")')
    write (funit,'(3x,"atoms: 1-",i0)') n1
  end if
  call write_cts(funit,env%cts)
  call write_cts_biasext(funit,env%cts)
  if (env%cts%used) then !Only, if user set constrians is an $end written
    write (funit,'(a)') '$end'
  end if
  close (funit)

end subroutine write_wall

subroutine check_dock(neg_E)
  use crest_parameters
  use crest_data
  use iomod,only:minigrep,grepval

  implicit none
  real(wp)             :: int_E
  logical,intent(out) :: neg_E
  logical :: ex
  character(len=*),parameter :: filename = 'xtbscreen.xyz'

  neg_E = .false.
  int_E = 0.0_wp

  call minigrep('xtb_dock.out','  Lowest Interaction Energy: ********** kcal/mol',ex)
  if (ex) return

  call grepval('xtb_dock.out','Lowest Interaction Energy:',ex,int_E)

  if (ex.and.int_E < 0.0_wp) neg_E = .true.

end subroutine check_dock

subroutine write_constraint(env,coord_name,fname)
  use crest_parameters
  use crest_data
  use iomod

  implicit none

  type(systemdata)     :: env
  character(len=*),intent(in)    :: fname,coord_name
  integer :: funit

  call copysub(coord_name,'coord.ref')
  open (newunit=funit,file=fname)
  call write_cts(funit,env%cts)
  call write_cts_biasext(funit,env%cts)
  if (env%cts%used) then !Only, if user set constrians is an $end written
    write (funit,'(a)') '$end'
  end if
  close (funit)

end subroutine write_constraint

!==============================================================================!

subroutine get_interaction_E(env,solu,solv,clus,iter,E_inter)
  use iso_fortran_env,wp => real64
  use crest_data
  use iomod
  use qcg_coord_type
  use strucrd
  implicit none

  type(systemdata)           :: env
  type(coord_qcg),intent(in) :: solu,solv,clus
  real(wp)                   :: e_cluster,e_solute,e_solvent
  real(wp)                   :: E_inter(env%nsolv)           ! interaction energy
  integer                    :: iter
  logical                    :: e_there

  call remove('cluster.coord')

!--- Prepare input coordinate files
  call clus%write('cluster.coord')
  call wr_cluster_cut('cluster.coord',solu%nat,solv%nat,iter,'solute_cut.coord','solvent_cut.coord')

!--- Perform single point calculations and recieve energies
  call xtb_sp_qcg(env,'solute_cut.coord',e_there,e_solute)
  if (.not.e_there) write (stdout,*) 'Solute energy not found'

  call xtb_sp_qcg(env,'solvent_cut.coord',e_there,e_solvent)
  if (.not.e_there) write (stdout,*) 'Solvent energy not found'

  call xtb_sp_qcg(env,'cluster.coord',e_there,e_cluster)
  if (.not.e_there) write (stdout,*) 'Cluster energy not found'

  E_inter(iter) = e_cluster-e_solute-e_solvent

end subroutine get_interaction_E

!===============================================================================!
subroutine chdirdbug(path)
  implicit none
  character(len=*),intent(in) :: path
  logical,parameter :: debug = .false.
  character(len=500) :: debugpath
  call chdir(path)
  if (debug) then
    call getcwd(debugpath)
    write (*,'(a,a)') '>>>>>>> NOW IN ',trim(debugpath)
  end if
end subroutine chdirdbug

!===============================================================================!
subroutine qcg_envcalc_reinit(env,mol,addconstraints,printinfo)
  use crest_parameters
  use crest_data
  use qcg_coord_type
  use strucrd
  use parse_xtbinput
  implicit none
  type(systemdata),intent(inout) :: env
  type(coord),intent(inout) :: mol
  logical,intent(in) :: addconstraints
  logical,intent(in) :: printinfo

  !> clear old data
  if (allocated(env%calc%calcs)) deallocate (env%calc%calcs)
  env%calc%ncalculations = 0

  !> and re-initialize
  call env2calc(env,env%calc,mol)

  !> add constraints from 'cts' if we want this
  if (addconstraints) then
    if (allocated(env%calc%cons)) deallocate (env%calc%cons)
    env%calc%nconstraints = 0
    if (printinfo) then
      write (stdout,'(a,a,a)') 'Parsing xtb-type constraints from internal backup to set up calculators ...'
    end if
    call parse_constraints_from_cts(env%calc,mol,env%cts)
  end if

  !> do a printout to stdout, if selected
  if (printinfo) then
    call env%calc%info(stdout)
  end if

end subroutine qcg_envcalc_reinit

