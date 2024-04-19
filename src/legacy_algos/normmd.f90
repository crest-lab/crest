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

!--------------------------------------------------------------------------------------------
! A single NORMMD run (bzw. its setup)
!--------------------------------------------------------------------------------------------
subroutine normalMD(fname,env,nr,newtemp,newtime)
  use crest_parameters
  use crest_data
  use iomod
  use strucrd,only:wrc0
  use utilities
  implicit none

  type(systemdata) :: env
  real(wp) :: newtemp
  real(wp) :: newtime
  character(len=*)  :: fname
  character(len=256) :: basename,dirname
  character(len=512) :: tmppath,thispath

  integer :: r,nr

  real(wp) :: k,alpha

!---- some settings
  basename = 'NORMMD'  !base name of the directories

  call getcwd(thispath)              !current dir= thispath
  write (dirname,'(a,i0)') trim(basename),nr
  tmppath = trim(dirname)
  call rmrf(tmppath)            !clear old directory
  r = makedir(trim(tmppath))    !make new directory

  call rename(trim(fname),trim(tmppath)//'/'//'coord')
  call env%wrtCHRG(trim(tmppath))
  call copysub(env%fixfile,trim(tmppath))
  if (env%useqmdff) then
    call copysub('solvent',trim(tmppath))
  end if

  if (env%staticmtd) then
    r = sylnk(trim(thispath)//'/'//env%mtdstaticfile,trim(tmppath)//'/'//env%mtdstaticfile)
  end if

  call chdir(trim(tmppath))  !switch to working directory
!---- do stuff here
  call setMDrun2('coord',env%hmass,newtime,newtemp,env%mdstep,env%shake, &
  &    env%mddumpxyz,env%mdskip,env%mddump,-1,env%cts)

  if (env%staticmtd) then
    k = env%rednat*0.0025d0
    alpha = 1.0d0
    call setstaticmtd('coord',k,alpha,env%emtd%mtdramp, &
    &    env%nstatic,env%mtdstaticfile,'')
  end if

  call chdir(thispath)  !go back to orgiginal directory

end subroutine normalMD

!--------------------------------------------------------------------------------------------
! Run several normal MDs on the lowermost conformers
!--------------------------------------------------------------------------------------------
subroutine normalMD_para_OMP(env,lconf,ntemps)
  use crest_parameters
  use crest_data
  use iomod
  use strucrd,only:wrc0,rdensembleparam,rdensemble
  use utilities
  implicit none

  type(systemdata) :: env

  integer :: lconf,ntemps,tot

  integer :: i,j,k
  integer :: TOTAL,vz,io
  real(wp) :: time
  character(len=512) :: thispath,tmppath
  character(len=512) :: jobcall
  character(len=80)  :: fname,pipe,inpnam,outnam,sname

  integer :: nat
  integer :: iz2
  integer,allocatable :: at(:)
  real(wp),allocatable :: xyz(:,:,:),eread(:),scoords(:,:,:)
  real(wp),allocatable :: temperatures(:),temperatures2(:)

  integer :: nclustbackup

  real(wp) :: newtemp

  nat = env%nat

  if (.not.env%entropic) then
    !--- by default the MDs are started on the (so far) lowest conformers.
    call rdensembleparam(conformerfile,iz2,TOTAL)
    allocate (xyz(3,env%nat,TOTAL),at(nat),eread(TOTAL))
    call rdensemble(conformerfile,iz2,TOTAL,at,xyz,eread)
    if (env%staticmtd) then
      env%mtdstaticfile = conformerfile
      env%nstatic = lconf
    end if
  else
    if (env%staticmtd) then
      call rdensembleparam(conformerfile,iz2,TOTAL)
      lconf = min(nint(float(TOTAL)*0.1d0),lconf)
    end if
    !--- for more diversity use PCA/K-means clusters as seeds
    nclustbackup = env%maxcluster
    env%nclust = lconf
    call create_anmr_dummy(nat)
    call CCEGEN(env,.true.,conformerfile)
    call remove('anmr_nucinfo')
    env%nclust = nclustbackup
    call rdensembleparam(clusterfile,iz2,TOTAL)
    allocate (xyz(3,env%nat,TOTAL),at(nat),eread(TOTAL))
    call rdensemble(clusterfile,iz2,TOTAL,at,xyz,eread)
    if (env%staticmtd) then
      env%mtdstaticfile = clusterfile
      env%nstatic = TOTAL
    end if
  end if
  if (TOTAL .lt. lconf) then
    lconf = TOTAL                 !if there are less conformers than lconf, use only that much
  end if
  allocate (scoords(3,nat,lconf))
  scoords = xyz(:,:,1:lconf)/bohr
  deallocate (eread,xyz)

!---- some settings
  tot = lconf*ntemps  !total number of MDs
  call new_ompautoset(env,'auto',tot,i,j)

  if (env%nmdtemp .lt. 0.0d0) then
    newtemp = 400.00d0
  else
    newtemp = env%nmdtemp
  end if

  allocate (temperatures2(tot),temperatures(ntemps))
  temperatures(1) = newtemp
  do i = 2,ntemps
    j = i-1
    if (env%entropic) then
      temperatures(i) = temperatures(j)+200.0_wp
    else
      temperatures(i) = temperatures(j)+100.0_wp   !each MD 100K higher
    end if
  end do

  if (env%entropic) then
    time = env%mdtime*0.25d0 !shorter for entropy mode
  else
    time = env%mdtime*0.5d0 !half the time of each MTD run
  end if

  call getcwd(thispath)

  fname = 'coord'
  pipe = ' > xtb.out 2>/dev/null'

  write (jobcall,'(a,1x,a,1x,a,'' --md '',a,1x,a,a)') &
  &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),pipe
  !--- slightly different jobcall for QMDFF usage
  if (env%useqmdff) then
    write (jobcall,'(a,1x,a,1x,a,'' --md --qmdff'',a,1x,a,a)') &
    &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),pipe
  end if

!---- Small Header
  write (*,'(''-----------------------------------------------'')')
  if (env%staticmtd) then
    write (*,'(''Additional static MTDs on lowest '',i0,'' conformer(s)'')') lconf
  else
    if (.not.env%entropic) then
      write (*,'(''Additional regular MDs on lowest '',i0,'' conformer(s)'')') lconf
    else
      write (*,'(''Additional regular MDs on  '',i0,'' conformer(s)'')') lconf
    end if
  end if
  write (*,'(''-----------------------------------------------'')')
  if (env%staticmtd) then
    write (*,'(/,"Using ",i0," RMSD bias potentials")') env%nstatic
  end if

!---- set up directories
  call rmrfw('scoord.')
  k = 0
  do i = 1,lconf
    write (sname,'(a,i0)') 'scoord.',i
    do j = 1,ntemps
      call wrc0(trim(sname),nat,at,scoords(:,:,i)) !has to be in the inner loop
      k = k+1
      newtemp = temperatures(j)
      temperatures2(k) = newtemp
      call normalMD(trim(sname),env,k,newtemp,time)
    end do
  end do
  deallocate (at,scoords)

!$omp parallel &
!$omp shared( vz,jobcall,tot,time,env,newtemp,temperatures2 )
!$omp single
  do i = 1,tot
    vz = i
    !$omp task firstprivate( vz ) private( tmppath,io )
    call initsignal()
    !$omp critical
    write (*,'(a,i4,a)') 'Starting MD',vz,' with the settings:'
    write (*,'(''     MD time /ps        :'',f8.1)') time
    write (*,'(''     MD Temperature /K  :'',f8.1)') temperatures2(vz)
    write (*,'(''     dt /fs             :'',f8.1)') env%mdstep
    write (*,'(''     dumpstep(trj) /fs  :'',i8)') env%mddumpxyz
    !$omp end critical
    write (tmppath,'(a,i0)') 'NORMMD',vz
    call command('cd '//trim(tmppath)//' && '//trim(jobcall),io)
    write (*,'(a,i0,a)') '*MD ',vz,' finished*'
    call rmrf(trim(tmppath)//'/scoord.*')
    !$omp end task
  end do
!$omp taskwait
!$omp end single
!$omp end parallel

  if (env%trackorigin) then
    call set_trj_origins('NORMMD','md')
  end if
  call checkname_xyz(crefile,inpnam,outnam)
  write (*,'('' Appending file '',a,'' with new structures'')') trim(inpnam)
  call collect_trj_skipfirst('NORMMD',trim(inpnam))  !collect all 'xtb.trj' from the NORMMD directories, but skip the first point! (input structure)
  write (*,*)

  if (.not.env%keepModef) then
    call cleanMTD
  end if

  return
end subroutine normalMD_para_OMP

!-------------------------------------------------------------------------------!
! Statical MTD addition to the normal MDs
! can be used to partially prevent the MD to
! converge to already known regions of the PES
!-------------------------------------------------------------------------------!
subroutine setstaticmtd(infile,k,alpha,ramp,nset,statfile,atomlist)
  use crest_parameters
  use crest_data
  implicit none
  character(len=*) :: infile
  real(wp) :: k
  real(wp) :: alpha
  real(wp) :: ramp
  integer :: nset
  character(len=*) :: statfile
  character(len=*) :: atomlist
  integer :: atlen

  character(len=256) :: atmp
  integer :: ich,ich2,iost

  open (newunit=ich,file=infile)
  open (newunit=ich2,file='tmpcoordfile')
  do
    read (ich,'(a)',iostat=iost) atmp
    if (iost < 0) exit
    if (index(atmp,'$end') .ne. 0) then
      cycle
    else
      write (ich2,'(a)') trim(atmp)
    end if
  end do

  atlen = len_trim(atomlist)

  write (ich2,'(a)') '$metadyn'
  write (ich2,'(2x,a,i0)') 'save=',nset
  write (atmp,'(f12.6)') k
  write (ich2,'(2x,a,a)') 'kpush=',adjustl(trim(atmp))
  write (atmp,'(f12.6)') alpha
  write (ich2,'(2x,a,a)') 'alp=',adjustl(trim(atmp))
  write (ich2,'(2x,a,a)') 'coord=',trim(statfile)
  write (ich2,'(2x,a)') 'static=true'
  if ((ramp-0.03_wp) .ne. 0.0_wp) then
    write (atmp,'(f12.6)') ramp
    write (ich2,'(2x,a,a)') 'ramp=',adjustl(trim(atmp))
  end if
  if (atlen > 1) then
    write (ich2,'(2x,a,1x,a)') 'atoms:',trim(atomlist)
  end if
  write (ich2,'(a)') '$end'

  close (ich)
  close (ich2)

  call rename('tmpcoordfile',trim(infile))

  return
end subroutine setstaticmtd

!--------------------------------------------------------------------------------------------
! Run several static MTDs on the lowermost conformers
!--------------------------------------------------------------------------------------------
subroutine entropyMD_para_OMP(env)
  use crest_parameters
  use crest_data
  use iomod
  use strucrd,only:wrc0,rdensembleparam,rdensemble
  use utilities
  implicit none
  type(systemdata) :: env
  integer :: lconf,tot
  integer :: i,j,k
  integer :: TOTAL,vz,io
  real(wp) :: time
  character(len=512) :: thispath,tmppath
  character(len=512) :: jobcall
  character(len=80)  :: fname,pipe,inpnam,outnam,sname
  integer :: iz2
  integer,allocatable :: at(:)
  real(wp),allocatable :: xyz(:,:,:),eread(:),scoords(:,:,:)
  real(wp),allocatable :: kprint(:)
  integer :: nclustbackup
  real(wp) :: newtemp
  real(wp) :: kpush,alpha
  logical :: ex
  associate (nat => env%nat)
    !--- use PCA/K-means clusters as seeds
    nclustbackup = env%maxcluster

    !--- first the structures that are used as bias
    env%nclust = env%emtd%nbias
    call create_anmr_dummy(nat)
    call smallhead('determining bias structures')
    call CCEGEN(env,.false.,conformerfile)
    call rdensembleparam(clusterfile,iz2,TOTAL)
    if (TOTAL < 1) then
      call copy('crest_best.xyz',clusterfile)
      TOTAL = 1
    end if
    write (*,'(1x,i0,a)') TOTAL,' structures were selected'
    write (*,'(1x,a,/)') 'done.'
    env%mtdstaticfile = "crest_bias.xyz"
    env%nstatic = TOTAL
    call rename(clusterfile,env%mtdstaticfile)

    !--- then get the input structures
    env%nclust = env%emtd%nMDs
    call smallhead('determining MTD seed structures')
    call CCEGEN(env,.false.,conformerfile)
    call rdensembleparam(clusterfile,iz2,TOTAL)
    write (*,'(1x,i0,a)') TOTAL,' structures were selected'
    write (*,'(1x,a,/)') 'done.'
    !--- cleanup
    call remove('anmr_nucinfo')
    env%nclust = nclustbackup

    allocate (xyz(3,env%nat,TOTAL),at(nat),eread(TOTAL))
    call rdensemble(clusterfile,iz2,TOTAL,at,xyz,eread)
    lconf = TOTAL
    tot = TOTAL
    allocate (scoords(3,nat,lconf))
    scoords = xyz(:,:,1:lconf)/bohr
    deallocate (eread,xyz)

!---- some settings
    call new_ompautoset(env,'auto',tot,i,j)

    !--- Temperature
    newtemp = env%nmdtemp
    !--- Runtime
    time = env%mdtime*env%emtd%lenfac
    !--- kpush & alpha
    !kpush = env%rednat*env%emtd%kpush
    kpush = env%emtd%katoms*env%emtd%kpush
    alpha = env%emtd%alpha

    call getcwd(thispath)

    fname = 'coord'
    pipe = ' > xtb.out 2>/dev/null'

    write (jobcall,'(a,1x,a,1x,a,'' --md '',a,1x,a,a)') &
    &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),pipe

!---- Small Header
    write (*,'(''-----------------------------------------------'')')
    write (*,'(''Additional static MTDs on lowest '',i0,'' conformer(s)'')') lconf
    write (*,'(''-----------------------------------------------'')')
    if (env%staticmtd) then
      write (*,'(/,"Using ",i0," RMSD bias potentials")') env%nstatic
    end if

!---- track kpush for printout
    j = tot*env%emtd%nklist
    allocate (kprint(j))

!---- set up directories
    call rmrfw('scoord.')
    k = 0
    do i = 1,lconf
      if (env%emtd%nklist > 1) then
        do j = 1,env%emtd%nklist
          write (sname,'(a,i0)') 'scoord.',i
          call wrc0(trim(sname),nat,at,scoords(:,:,i)) !has to be in the inner loop
          kpush = env%emtd%klist(j)*env%emtd%katoms
          k = k+1
          call entropyMD(trim(sname),env,k,newtemp,time,kpush,alpha)
          kprint(k) = kpush
        end do
      else
        write (sname,'(a,i0)') 'scoord.',i
        call wrc0(trim(sname),nat,at,scoords(:,:,i)) !has to be in the inner loop
        k = k+1
        call entropyMD(trim(sname),env,k,newtemp,time,kpush,alpha)
        kprint(k) = kpush
      end if
    end do
    deallocate (at,scoords)

    tot = k

!$omp parallel &
!$omp shared( vz,jobcall,tot,time,env,newtemp )
!$omp single
    do i = 1,tot
      vz = i
      !$omp task firstprivate( vz ) private( tmppath,io,ex )
      call initsignal()
      !$omp critical
      write (*,'(a,i4,a)') 'Starting MTD',vz,' with the settings:'
      write (*,'(''     MD time /ps        :'',f8.1)') time
      write (*,'(''     MD Temperature /K  :'',f8.1)') newtemp
      write (*,'(''     dt /fs             :'',f8.1)') env%mdstep
      write (*,'(''     dumpstep(trj) /fs  :'',i8)') env%mddumpxyz
      write (*,'(''     Vbias factor k /Eh :'',f8.4)') kprint(vz)
      write (*,'(''     Vbias exp α /bohr⁻²:'',f8.2)') alpha
      !$omp end critical
      write (tmppath,'(a,i0)') 'STATICMTD',vz
      call command('cd '//trim(tmppath)//' && '//trim(jobcall),io)
      inquire (file=trim(tmppath)//'/'//'xtb.trj',exist=ex)
      if (.not.ex.or.io .ne. 0) then
        write (*,'(a,i0,a)') '*Warning: static MTD ',vz,' seemingly failed (no xtb.trj)*'
        call command('cp -r '//trim(tmppath)//' FAILEDMTD')
      else
        write (*,'(a,i0,a)') '*static MTD ',vz,' finished*'
      end if
      call rmrf(trim(tmppath)//'/'//'scoord.*')
      !$omp end task
    end do
!$omp taskwait
!$omp end single
!$omp end parallel

    deallocate (kprint)

    if (env%trackorigin) then
      call set_trj_origins('STATICMTD','smtd')
    end if
    call checkname_xyz(crefile,inpnam,outnam)
    write (*,'('' Appending file '',a,'' with new structures'')') trim(inpnam)
    call collect_trj_skipfirst('STATICMTD',trim(inpnam))  !collect all 'xtb.trj' from the NORMMD directories, but skip the first point! (input structure)
    write (*,*)

    if (.not.env%keepModef) then
      call cleanMTD
    end if

  end associate

end subroutine entropyMD_para_OMP

!--------------------------------------------------------------------------------------------
! A single NORMMD run (bzw. its setup)
!--------------------------------------------------------------------------------------------
subroutine entropyMD(fname,env,nr,newtemp,newtime,k,alpha)
  use crest_parameters
  use crest_data
  use iomod
  use strucrd,only:wrc0
  use utilities
  implicit none

  type(systemdata) :: env
  real(wp) :: newtemp
  real(wp) :: newtime
  character(len=*)  :: fname
  character(len=256) :: basename,dirname
  character(len=512) :: tmppath,thispath

  integer :: r,nr

  real(wp) :: k,alpha

!---- some settings
  basename = 'STATICMTD'  !base name of the directories

  call getcwd(thispath)              !current dir= thispath
  write (dirname,'(a,i0)') trim(basename),nr
  tmppath = trim(dirname)
  call rmrf(tmppath)            !clear old directory
  r = makedir(trim(tmppath))    !make new directory

  call rename(trim(fname),trim(tmppath)//'/'//'coord')
  call env%wrtCHRG(trim(tmppath))
  call copysub(env%fixfile,trim(tmppath))
  if (env%gfnver == '--gff') then
!            r = sylnk(trim(thispath)//'/'//'gfnff_topo',trim(tmppath)//'/'//'gfnff_topo')
  end if

  r = sylnk(trim(thispath)//'/'//env%mtdstaticfile,trim(tmppath)//'/'//env%mtdstaticfile)

  call chdir(trim(tmppath))  !switch to working directory
!---- do stuff here
  call setMDrun2('coord',env%hmass,newtime,newtemp,env%mdstep,env%shake, &
  &    env%mddumpxyz,env%mdskip,env%mddump,-1,env%cts)

!--- static MTD settings
  call setstaticmtd('coord',k,alpha,env%emtd%mtdramp,env%nstatic, &
  &    env%mtdstaticfile,env%emtd%atomlist)

!--- back to original DIR
  call chdir(thispath)

  return
end subroutine entropyMD
