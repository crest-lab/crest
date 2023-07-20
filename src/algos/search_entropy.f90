!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022 Philipp Pracht
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

subroutine crest_search_entropy(env,tim)
!*******************************************************************
!* This is the re-implementation of CREST's sMTD-iMTD workflow
!* from https://doi.org/10.1039/d1sc00621e
!* with calculation of conformational entropy
!* This is a TODO
!*******************************************************************
  use crest_parameters, only: wp,stdout
  use crest_data
  use strucrd
  use calc_type
  use calc_module
  use dynamics_module
  use shake_module
  use iomod
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,m
  logical :: pr,wr
!===========================================================!
  type(calcdata) :: calc
  type(mddata) :: mddat
  type(shakedata) :: shk

  type(mddata),allocatable :: mddats(:)
  integer :: nsim,nallout

  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)
  character(len=:),allocatable :: ensnam
  integer :: nat,nall
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  logical :: dump,ex
  character(len=80) :: atmp,btmp,str
  logical :: multilevel(6)
  logical :: start,lower
!===========================================================!
!>--- printout header
  write (stdout,*)
  write (stdout,'(10x,"┍",49("━"),"┑")')
  write (stdout,'(10x,"│",14x,a,13x,"│")') "CREST ENTROPY SAMPLING"
  write (stdout,'(10x,"┕",49("━"),"┙")')
  write (stdout,*)
  write (stdout,'(1x,a)') 'please cite:'
  write (stdout,'(1x,a)') '• P.Pracht, S.Grimme, Chem. Sci., 2021, 12, 6551-6568.'
  write (stdout,'(1x,a)') '• J.Gorges, S.Grimme, A.Hansen, P.Pracht, PCCP, 2022,24, 12249-12259.'
  write (stdout,*)

!===========================================================!
!>--- setup
  call env%ref%to(mol)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)

!>--- sets the MD length according to a flexibility measure
  call md_length_setup(env) 
  call env_to_mddat(env)

!===========================================================!
!>--- Start mainloop 
  env%nreset = 0
  start = .true.
  MAINLOOP : do
    call printiter
    if (.not. start) then
!>--- clean Dir for new iterations, but leave iteration backup files
      call clean_V2i 
      env%nreset = env%nreset + 1
    else 
!>--- at the beginning, wipe directory clean
      call V2cleanup(.false.)
    end if
!===========================================================!
!>--- Meta-dynamics loop
  mtdloop: do i = 1,env%Maxrestart

    write(stdout,*)
    write(stdout,'(1x,a)') '------------------------------'
    write(stdout,'(1x,a,i0)') 'Meta-Dynamics Iteration ',i
    write(stdout,'(1x,a)') '------------------------------'

    nsim = -1 !>--- enambles automatic MTD setup in init routines
    call crest_search_multimd_init(env,mol,mddat,nsim)
    allocate (mddats(nsim), source=mddat)
    call crest_search_multimd_init2(env,mddats,nsim)

    call tim%start(2,'MTD simulations')
    call crest_search_multimd(env,mol,mddats,nsim)
    call tim%stop(2)
!>--- a file called crest_dynamics.trj should have been written
    ensnam = 'crest_dynamics.trj'
!>--- deallocate for next iteration
    if(allocated(mddats))deallocate(mddats)

!==========================================================!
!>--- Reoptimization of trajectories
    call tim%start(3,'geom. optimization')
    multilevel = (/.true.,.false.,.false.,.false.,.true.,.false./)
    call crest_multilevel_oloop(env,ensnam,multilevel)
    call tim%stop(3)

!>--- save the CRE under a backup name
    call checkname_xyz(crefile,atmp,str)
    call checkname_xyz('.cre',str,btmp)
    call rename(atmp,btmp)
!>--- save cregen output
    call checkname_tmp('cregen',atmp,btmp)
    call rename('cregen.out.tmp',btmp)

!=========================================================!
!>--- cleanup after first iteration and prepare next
    if (i .eq. 1 .and. start) then
      start = .false.
!>-- obtain a first lowest energy as reference
      env%eprivious = env%elowest
!>-- remove the two extreme-value MTDs
      if (.not. env%readbias .and.  env%runver .ne. 33 .and. &
      &   env%runver .ne. 787878 ) then
        env%nmetadyn = env%nmetadyn - 2
      end if
!>-- the cleanup 
      call clean_V2i   
!>-- and always do two cycles of MTDs
      cycle mtdloop 
    endif
!=========================================================!
!>--- Check for lowest energy
    call elowcheck(lower,env)
    if (.not. lower) then
      exit mtdloop
    end if
  enddo mtdloop
!=========================================================!
!>--- collect all ensembles from mtdloop and merge
  write(stdout,*)
  write (stdout,'(''========================================'')')
  write (stdout,'(''           MTD Simulations done         '')')
  write (stdout,'(''========================================'')')
  write (stdout,'(1x,''Collecting ensmbles.'')')
!>-- collecting all ensembles saved as ".cre_*.xyz"
  call collectcre(env)                      
  call newcregen(env,0)
  call checkname_xyz(crefile,atmp,btmp)
!>--- remaining number of structures
  call remaining_in(atmp,env%ewin,nallout) 

!=========================================================!
!>--- (optional) Perform additional MDs on the lowest conformers
  if (env%rotamermds) then
    call tim%start(4,'MD simulations')
    call crest_rotamermds(env,conformerfile)
    call tim%stop(4)

!>--- Reoptimization of trajectories
    call checkname_xyz(crefile,atmp,btmp)
    write(stdout,'('' Appending file '',a,'' with new structures'')')trim(atmp)
    ensnam = 'crest_dynamics.trj'
    call appendto(ensnam,trim(atmp))
    call tim%start(3,'geom. optimization')
    call crest_multilevel_wrap(env,trim(atmp),5)
    call tim%stop(3)

    call elowcheck(lower,env)
    if (lower) then
      call checkname_xyz(crefile,atmp,str)
      call checkname_xyz('.cre',str,btmp)
      call rename(atmp,btmp)
      cycle MAINLOOP
    end if
  end if

!=========================================================!
!!>---- Entropy mode iterative statically biased MDs
!    if (env%entropymd) then
!      call mtdatoms(env)
!      call tim%start(6,'static MTD')
!      call emtdcopy(env,0,stopiter,fail)
!      bref = env%emtd%nbias
!
!      ENTROPYITER: do eit = 1,env%emtd%iter
!        dum = nint(float(env%emtd%nbias) * env%emtd%nbiasgrow)
!        !env%emtd%nbias = nint(float(env%emtd%nbias) * env%emtd%nbiasgrow)
!        env%emtd%nbias = max(env%emtd%nbias + 1,dum)
!        fail = .false.
!        EFALLBACK: do k = 1,env%emtd%maxfallback
!          call printiter2(eit)
!          call tim%start(6,'static MTD')
!          call entropyMD_para_OMP(env)
!          call tim%stop(6)
!          call emtdcheckempty(env,fail,env%emtd%nbias)
!          if (fail) then
!            if (k == env%emtd%maxfallback) then
!              stopiter = .true.
!            else
!              cycle EFALLBACK
!            end if
!          else
!            call tim%start(3,'multilevel OPT')
!            if (env%optlev >= -1.0d0) then
!              call multilevel_opt(env,2)
!            end if
!            call multilevel_opt(env,99)
!
!            call tim%stop(3)
!            !--- if in the entropy mode a lower structure was found
!            !    --> cycle, required for extrapolation
!            call elowcheck(lower,env)
!            if (lower .and. env%entropic) then
!              env%emtd%nbias = bref  !IMPORTANT, reset for restart
!              cycle MAINLOOP
!            end if
!            !--- file handling
!            eit2 = eit
!            call emtdcopy(env,eit2,stopiter,fail)
!            env%emtd%iterlast = eit2
!          end if
!          if (.not. lower .and. fail .and. .not. stopiter) then
!            cycle EFALLBACK
!          end if
!          exit EFALLBACK  !fallback loop is exited on first opportuinity
!        end do EFALLBACK
!        if (stopiter) then
!          exit ENTROPYITER
!        end if
!      end do ENTROPYITER
!    end if

!==========================================================!
!>--- exit mainloop
   exit MAINLOOP
  enddo MAINLOOP

!==========================================================!
!!>--- final ensemble optimization
!    write (stdout,'(/)')
!    write (stdout,'(3x,''================================================'')')
!    write (stdout,'(3x,''|           Final Geometry Optimization        |'')')
!    write (stdout,'(3x,''================================================'')')
!    call tim%start(3,'geom. optimization')
!    call checkname_xyz(crefile,atmp,str)
!    call crest_multilevel_wrap(env,trim(atmp),0) 
!    call tim%stop(3)                 

!==========================================================!
!>--- print CREGEN results and clean up Directory a bit
    write (stdout,'(/)')
    call smallhead('Final Ensemble Information')
    call V2terminating()

!==========================================================!
  return
end subroutine crest_search_entropy

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!

!========================================================================================!
subroutine crest_smtd_mds(env,ensnam)
!***********************************************************
!* set up and perform several sMTD's on a number of 
!* conformers obtained from clustering
!* This is a TODO
!***********************************************************
  use crest_parameters, only: wp,stdout,bohr
  use crest_data
  use strucrd
  use calc_type
  use calc_module
  use dynamics_module
  use shake_module
  implicit none
  type(systemdata),intent(inout) :: env
  character(len=*),intent(in) :: ensnam

  integer :: nsim
  type(mddata) :: mddat
  type(mddata),allocatable :: mddats(:)
  type(coord) :: mol 
  type(coord),allocatable :: mols(:)
  integer :: nat,nall
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  integer :: nstrucs,i,j,k,io
  real(wp) :: temp,newtemp
  character(len=80) :: atmp
  
!>--- coord setup
  call env%ref%to(mol)
  call rdensembleparam(ensnam,nat,nall)
  if (nall .lt. 1) then
    write(stdout,*) 'empty ensemble file',trim(ensnam)
    return
  endif

!>--- determine how many MDs need to be run and setup
  call adjustnormmd(env)
  nstrucs = min(nall, env%nrotammds)
  nsim = nstrucs * env%temps 
  call crest_search_multimd_init(env,mol,mddat,nsim)
  allocate (mddats(nsim), source=mddat)
  call crest_search_multimd_init2(env,mddats,nsim)
!>--- adjust T's and runtimes
  k = 0
  do i=1,env%temps
    !> each T block 100K higher
    temp = env%nmdtemp + (i-1)*100.0_wp
    do j=1,nstrucs
      k= k + 1
      mddats(k)%tsoll = temp
      !> reduce runtime by 50% compared to MTDs
      mddats(k)%length_ps = mddats(k)%length_ps * 0.5_wp
      call mdautoset(mddats(k),io)
    enddo
  enddo 

!>--- read ensemble and prepare mols
  allocate (xyz(3,nat,nall),at(nat),eread(nall))
  call rdensemble(ensnam,nat,nall,at,xyz,eread)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: mols must be in Bohrs
  xyz = xyz / bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
  allocate(mols(nsim), source=mol) 
  k = 0
  do i=1,env%temps
    do j=1,nstrucs
      k = k + 1
      mols(k)%at = at
      mols(k)%xyz(:,:) = xyz(:,:,j)
    enddo
  enddo
  deallocate(eread,at,xyz)
  
!>--- print what we are doing
  write(stdout,*)
  write(atmp,'(''Additional regular MDs on lowest '',i0,'' conformer(s)'')')nstrucs
  call smallheadline(trim(atmp)) 

!>--- and finally, run the MDs
  call crest_search_multimd2(env,mols,mddats,nsim)

  if(allocated(mols))deallocate(mols)
  if(allocated(mddats))deallocate(mddats) 
  return
end subroutine crest_smtd_mds
