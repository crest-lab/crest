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

!> QCG front-ends for the in-process aISS docking and the cluster optimizations
module qcg_aiss
  use crest_parameters,only:wp,stdout
  use crest_data
  use strucrd,only:coord,get_atlist
  use qcg_coord_type
  use crest_calculator,only:calcdata,calculation_settings,engrad
  use optimize_module,only:optimize_geometry
  use aiss_dock
  implicit none
  private

  public :: qcg_aiss_dock
  public :: qcg_aiss_ensemble_dock
  public :: qcg_opt_cluster
  public :: qcg_relax_in_walls

  character(len=*),parameter :: dockout = 'aiss_dock.out'

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!


!========================================================================================!

  subroutine qcg_aiss_dock(env,fnameA,fnameB,solu,clus,neg_E)
    !******************************************************************
    !* Grow step: dock one solvent (fnameB) to the solute or cluster
    !* (fnameA) inside the two QCG wall potentials, as "xtb dock
    !* --nfrag1 <solu%nat> --qcg --input xcontrol" did. The lowest
    !* structure is written to best.xyz.
    !* neg_E : a negative interaction energy was found
    !******************************************************************
    type(systemdata),intent(inout) :: env
    character(len=*),intent(in)    :: fnameA,fnameB
    type(coord_qcg),intent(in)     :: solu,clus
    logical,intent(out)            :: neg_E
    type(coord) :: molA,molB,best,comb
    type(aiss_config) :: cfg
    type(aiss_wall) :: walls(2)
    type(calcdata) :: ucons
    logical,allocatable :: atl(:),dirsel(:)
    integer,allocatable :: directed(:)
    real(wp) :: eint
    integer :: i,io,ich,nfix,T,Tn

    neg_E = .false.
    call molA%open(fnameA)
    call molB%open(fnameB)
    call new_ompautoset(env,'max',1,T,Tn)

    cfg%optlvl = aiss_level_from_flag(env%gfnver)
    cfg%lmolvl = trim(env%aiss_lmolvl)
    cfg%samerand = env%sameRandomNumber
    if (env%docking_qcg_flag == '--qcg') call cfg%set_qcg()
    cfg%nfrag1 = solu%nat
    cfg%chrg = [env%chrg,0]
    cfg%uhf = [env%uhf,0]

    walls(1)%radius = clus%ell_abc
    walls(2)%radius = solu%ell_abc
    walls(2)%list = [(i,i=1,solu%nat)]
    nfix = 0
    if (env%constrain_solu) nfix = solu%nat

!>--- directed docking: sites whose solvent quota is not yet filled
    allocate (directed(0))
    if (allocated(env%directed_file)) then
      allocate (dirsel(molA%nat),source=.false.)
      do i = 1,size(env%directed_number)
        if ((i == 1.and.env%directed_number(i) >= clus%nmol).or. &
        &   (i > 1.and.env%directed_number(i) >= clus%nmol.and. &
        &    env%directed_number(max(i-1,1)) < clus%nmol)) then
          call get_atlist(molA%nat,atl,trim(env%directed_list(i,1)),molA%at)
          dirsel = dirsel.or.atl
        end if
      end do
      directed = pack([(i,i=1,molA%nat)],dirsel)
    end if

    if (env%cts%used) then
      comb%nat = molA%nat+molB%nat
      comb%at = [molA%at,molB%at]
      allocate (comb%xyz(3,comb%nat))
      comb%xyz(:,1:molA%nat) = molA%xyz
      comb%xyz(:,molA%nat+1:) = molB%xyz
      call cts_to_calc(ucons,comb,env%cts)
    end if

    open (newunit=ich,file=dockout)
    call aiss_dock_run(molA,molB,cfg,best,eint,io,walls=walls,nfix=nfix, &
    &                  directed=directed,ucons=ucons,iunit=ich)
    close (ich)
    if (io /= 0) then
      write (stdout,'(a)') '**ERROR** aISS docking failed, see '//dockout
      call creststop(status_failed)
    end if
    neg_E = (eint < 0.0_wp)
  end subroutine qcg_aiss_dock

!========================================================================================!

  subroutine qcg_aiss_ensemble_dock(env,outer_ell_abc,nfrag1,frag1_file,frag2_file, &
  &                                 n_shell,n_solvent,NTMP,TMPdir,conv)
    !******************************************************************
    !* CFF step: dock one solvent (frag2_file) to the cluster
    !* (frag1_file) in each open directory TMPdir<conv(i)>, as
    !* "xtb dock --fast" with the shell frozen and a wall on the new
    !* solvent. The runs are serial, each uses all threads.
    !******************************************************************
    type(systemdata),intent(inout) :: env
    real(wp),intent(in)            :: outer_ell_abc(env%nqcgclust,3)
    integer,intent(in)             :: nfrag1,n_shell,n_solvent,NTMP
    character(len=*),intent(in)    :: frag1_file,frag2_file,TMPdir
    integer,intent(in)             :: conv(env%nqcgclust+1)
    type(coord) :: molA,molB,best
    type(aiss_config) :: cfg
    type(aiss_wall) :: walls(1)
    character(len=512) :: thispath,tmppath
    real(wp) :: eint
    integer :: i,io,ich,T,Tn

    call new_ompautoset(env,'max',1,T,Tn)
    call getcwd(thispath)
    cfg%optlvl = aiss_level_from_flag(env%gfnver)
    cfg%lmolvl = trim(env%aiss_lmolvl)
    cfg%samerand = env%sameRandomNumber
    call cfg%set_fast()
    cfg%nfrag1 = nfrag1
    walls(1)%list = [(i,i=n_shell+1,n_shell+n_solvent)]

    do i = 1,NTMP
      write (tmppath,'(a,i0)') trim(TMPdir),conv(i)
      call chdir(trim(tmppath))
      call molA%open(frag1_file)
      call molB%open(frag2_file)
      walls(1)%radius = outer_ell_abc(conv(i),:)
      open (newunit=ich,file=dockout)
      call aiss_dock_run(molA,molB,cfg,best,eint,io,walls=walls,nfix=n_shell,iunit=ich)
      close (ich)
      call chdir(trim(thispath))
      if (io /= 0) then
        write (stdout,'(a)') '**ERROR** aISS docking failed in '//trim(tmppath)
        call creststop(status_failed)
      end if
    end do
  end subroutine qcg_aiss_ensemble_dock

!========================================================================================!

  subroutine qcg_opt_cluster(env,fname,without_pot,etot)
    !******************************************************************
    !* Calculator version of opt_cluster. The xtb calls never read
    !* the xcontrol written next to them (no --input), so the walls
    !* are absent here as well:
    !*  without_pot=.false. : gas-phase opt + SP with implicit solvent
    !*  without_pot=.true.  : opt with implicit solvent
    !* The geometry goes to xtbopt.<ext of fname>, etot in Eh.
    !******************************************************************
    type(systemdata),intent(inout) :: env
    character(len=*),intent(in)    :: fname
    logical,intent(in)             :: without_pot
    real(wp),intent(out)           :: etot
    type(coord) :: mol,molopt
    type(calcdata) :: calc,spcalc
    real(wp),allocatable :: grd(:,:)
    integer :: io,T,Tn

    call new_ompautoset(env,'max',1,T,Tn)
    call mol%open(fname)
    allocate (grd(3,mol%nat))

    call calc%add(level(without_pot))
    calc%optlev = nint(env%optlev)
    etot = 0.0_wp
    call optimize_geometry(mol,molopt,calc,etot,grd,.false.,.false.,io)
    if (io /= 0) write (stdout,'(a)') 'Warning: cluster optimization did not converge'

    if (index(fname,'.xyz') > 0) then
      call molopt%write('xtbopt.xyz')
    else
      call molopt%write('xtbopt.coord')
    end if

    if (.not.without_pot) then
      call spcalc%add(level(.true.))
      call engrad(molopt,spcalc,etot,grd,io)
    end if

  contains

    function level(solvated) result(clevel)
      logical,intent(in) :: solvated
      type(calculation_settings) :: clevel
      logical :: have_solv
      have_solv = .false.
      if (allocated(env%solv)) have_solv = (len_trim(env%solv) > 0)
      if (solvated.and.have_solv) then
        call clevel%create(env%gfnver,chrg=env%chrg,uhf=env%uhf, &
        &                  solvmodel=env%solv,solvent=env%solvent)
      else
        call clevel%create(env%gfnver,chrg=env%chrg,uhf=env%uhf)
      end if
    end function level

  end subroutine qcg_opt_cluster

!========================================================================================!

  subroutine qcg_relax_in_walls(env,nfreeze,clus)
    !******************************************************************
    !* Relax the ensemble start cluster with the ensemble calculator
    !* (walls, bond constraints), atoms 1..nfreeze frozen as in the MD.
    !* The grown cluster was optimized without walls and can protrude
    !* from them, which heats up the first MD/MTD steps explosively.
    !******************************************************************
    type(systemdata),intent(inout) :: env
    integer,intent(in)             :: nfreeze
    type(coord_qcg),intent(inout)  :: clus
    type(calcdata) :: calc
    type(coord) :: mol,molopt
    real(wp),allocatable :: grd(:,:)
    real(wp) :: e0,e1
    integer :: io

    mol = clus%as_coord()
    calc = env%calc
    if (nfreeze > 0.and.nfreeze < mol%nat) call calc%set_freeze(mol%nat,1,nfreeze)
    allocate (grd(3,mol%nat))
    call engrad(mol,calc,e0,grd,io)
    if (io /= 0) return
    e1 = e0
    call optimize_geometry(mol,molopt,calc,e1,grd,.false.,.false.,io)
    if (io < 0.or.e1 >= e0) return
    call clus%from_coord(molopt)
    write (stdout,'(2x,a,f12.6,a,f12.6,a)') 'Start cluster relaxed in the wall potentials:', &
    &  e0,' ->',e1,' Eh'
  end subroutine qcg_relax_in_walls

!========================================================================================!

  subroutine cts_to_calc(calc,mol,cts)
    !> user constraints from the xtb-style backup (the $-blocks of write_cts)
    use parse_xtbinput,only:parse_constraints_from_cts
    type(calcdata),intent(inout) :: calc
    type(coord),intent(inout)    :: mol
    type(legacy_constraints),intent(in) :: cts
    call parse_constraints_from_cts(calc,mol,cts)
  end subroutine cts_to_calc

!========================================================================================!
end module qcg_aiss
