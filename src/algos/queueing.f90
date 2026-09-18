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

subroutine crest_queue_setup(env,iterate)
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use molbuilder_construct_list
  use molbuilder_construct_mod
  implicit none
  type(systemdata),intent(inout) :: env
  logical,intent(out) :: iterate

  integer :: splitlayers
  integer :: ii,jj,nn,kk,ich
  type(coord),pointer :: reference_mol
  type(coord),target :: mol
  integer,allocatable :: splitatms(:)
  integer :: parentlayer,parentnode
  character(len=1024) :: thispath
  real(wp),allocatable :: qat(:)
  integer,allocatable :: lq(:)

  iterate = .true.

  if (allocated(env%splitqueue)) then

    !> check for incompatible runtypes (or rather, whitelist a few)
    if (.not.any(env%crestver == [crest_imtd,crest_imtd2,crest_sp, &
      & crest_optimize,crest_moldyn,crest_rigcon,crest_trialopt,crest_bh,crest_test])) then
      write (stdout,'(a)') '** ERROR ** Selected CREST runtype incompatible with substructure builder'
      call creststop(status_config)
    end if
    if (allocated(env%ONIOM_input).or.allocated(env%ONIOM_toml)) then
      write (stdout,'(a)') '** ERROR ** ONIOM incompatible with substructure builder'
      call creststop(status_config)
    end if

    !> if the program sees no problem, set the global boolean
    env%substructure_queue = .true.
    splitlayers = size(env%splitqueue,1)

    !> we may need to calculate charges to distribute them:
    if (env%chrg .ne. 0) then
      call calc_charges(env,qat)
    else
      allocate (qat(env%ref%nat),source=0.0_wp)
    end if

    !> start constructing the splitheap
    env%splitheap%nlayer = splitlayers
    allocate (env%splitheap%layer(splitlayers))
    associate (heap => env%splitheap,layer => env%splitheap%layer)

      do ii = 1,heap%nlayer
        layer%id = ii

        nn = env%splitqueue(ii)%natms
        allocate (splitatms(nn))
        splitatms(:) = env%splitqueue(ii)%atms(:)
        if (ii == 1) then
          call env%ref%to(mol)
          reference_mol => mol
        else

          call pick_parent(heap,ii,splitatms,parentlayer,parentnode)
          if (parentlayer == 0) then
            call env%ref%to(mol)
            reference_mol => mol
          else
            mol = heap%layer(parentlayer)%node(parentnode)
            reference_mol => mol
            heap%layer(ii)%parent = parentlayer
            heap%layer(ii)%parentnode = parentnode
          end if
        end if
        layer(ii)%refmol = reference_mol
        call reference_mol%get_cn(layer(ii)%refcn)
        allocate (layer(ii)%reficn(reference_mol%nat))
        layer(ii)%reficn(:) = nint(layer(ii)%refcn(:))
        call binarysplit(reference_mol,splitatms,layer(ii)%node,layer(ii)%alignmap, &
          & ncap=layer(ii)%ncapped,position_mapping=layer(ii)%position_mapping)
        deallocate (splitatms)
        layer(ii)%nnodes = size(layer(ii)%node,1)
        call heap%map_origins_for_layer(ii)
        !> determening charges for fragments
        call sum_charges_layer(env,heap,ii,qat,lq)
        do jj = 1,layer(ii)%nnodes
          layer(ii)%node(jj)%chrg = lq(jj)
        end do
      end do

      call heap%setup_queue()
      call getcwd(thispath)
      !> some backups
      call env%ref%to(heap%originmol)
      heap%originmol%chrg = env%chrg
      heap%origindir = trim(thispath)
      heap%origincalc => env%calc

    end associate
    iterate = .true.
  end if

  return
contains
  subroutine pick_parent(heap,current_layer,splitatms,parentlayer,parentnode)
    use molbuilder_construct_list
    implicit none
    type(construct_heap),intent(inout) :: heap
    integer,intent(inout) :: splitatms(:)
    integer,intent(in) :: current_layer
    integer,intent(out) :: parentlayer,parentnode
    integer :: ii,jj,kk,prev_layer
    logical :: matching

    parentlayer = 0
    parentnode = 0
    if (current_layer .eq. 1) return

    !> iterate through the previous layer and check which node
    !> contains all the split atoms
    prev_layer = current_layer-1
    LAYITER: do while (prev_layer >= 1)
      do ii = 1,heap%layer(prev_layer)%nnodes
        matching = .true.
        do jj = 1,size(splitatms,1)
          matching = matching.and. &
            & any(heap%layer(prev_layer)%origin(ii)%map(:) .eq. splitatms(jj))
        end do
        if (matching) then
          parentlayer = prev_layer
          parentnode = ii
          !> on the first match, exit
          exit LAYITER
        end if
      end do
      !> if no matching parent node was found, try again in one layer further up
      if (parentnode == 0) prev_layer = prev_layer-1
    end do LAYITER

    !> IMPORTANT; we need to update the splitatms with the correctly mapped indices
    !> reflecting their position in the selected parent layer
    if (parentnode .ne. 0) then
      do ii = 1,size(splitatms,1)
        jj = splitatms(ii)
        call heap%find_current_position(jj,parentlayer,parentnode,kk)
        splitatms(ii) = kk
      end do
      !> we also map the current node as a child node of the selected parent
      if (.not.allocated(heap%layer(parentlayer)%childlayer)) then
        ii = heap%layer(parentlayer)%nnodes
        allocate (heap%layer(parentlayer)%childlayer(ii),source=0)
      end if
      heap%layer(parentlayer)%childlayer(parentnode) = current_layer
    end if

  end subroutine pick_parent
  subroutine calc_charges(env,qat)
    use tblite_api,only:tblite_quick_ceh_q
    implicit none
    type(systemdata),intent(inout) :: env
    real(wp),intent(out),allocatable :: qat(:)
    real(wp),allocatable :: qat0(:)
    character(len=256) :: atmp
    type(coord)  :: mol
    integer :: ii
    write (atmp,'(a)') 'Calculating atomic charges under consideration of molecular charge'
    call underline(trim(atmp))
    write (stdout,'(a,i0)') 'Molecular charge : ',env%chrg
    call env%ref%to(mol)
    call tblite_quick_ceh_q(mol,qat, &
      & chrg=env%chrg,uhf=env%uhf,pr=.true.,prch=stdout)
    write (stdout,'(a)') 'Obtained CEH charges for full structure:'
    do ii = 1,mol%nat
      write (stdout,'(3x,a3,2x,f10.6)') i2e(mol%at(ii)),qat(ii)!,qat0(ii),qat(ii)-qat0(ii)
    end do

    write (stdout,'(/,a)') 'NOTE: Total charge for each fragment will be selected automatically by'
    write (stdout,'(a)') '      matching the best atomic charge MAE to these charges.'
  end subroutine calc_charges
  subroutine sum_charges_layer(env,heap,lay,qat,lq)
    use tblite_api,only:tblite_quick_ceh_q
    implicit none
    type(systemdata) :: env
    type(construct_heap) :: heap
    integer,intent(in) :: lay
    real(wp),intent(in) :: qat(:)
    integer,intent(out),allocatable :: lq(:)
    integer :: ii,jj,nat,nnodes,kk,nnat,sign,cc,cc2,chrgs
    integer,allocatable :: ichrgs(:)
    real(wp) :: qtmp0,qtmpc
    real(wp),allocatable :: qtmp(:)
    real(wp),allocatable :: qattmp0(:),qattmpc(:),qattmpref(:),qdum(:)
    real(wp),allocatable :: qattmp(:,:)

    nat = size(qat,1)
    nnodes = heap%layer(lay)%nnodes
    allocate (lq(nnodes),source=0) !> default chrg of 0

    if (env%chrg == 0) return !> return for neutral systems (may need some implementation for zwitter ions)

    write (stdout,'(a,i0,a)') 'Calculating charges for fragments in layer ',lay,' ...'
    sign = 1
    if (env%chrg < 0) sign = -1
    chrgs = abs(env%chrg)+1
    allocate (qtmp(chrgs),source=0.0_wp)
    allocate (ichrgs(chrgs),source=0)
    cc2 = 0
    do cc = 0,env%chrg,sign
      cc2 = cc2+1
      ichrgs(cc2) = cc
    end do

    do ii = 1,nnodes
      qtmp(:) = 0.0_wp
      qtmp0 = 0.0_wp
      qtmpc = 0.0_wp
      nnat = heap%layer(lay)%node(ii)%nat
      allocate (qattmpref(nnat),source=0.0_wp)
      allocate (qattmp(nnat,chrgs))
      !> check different charge settings
      cc2 = 0
      do cc = 0,env%chrg,sign
        cc2 = cc2+1
        call tblite_quick_ceh_q(heap%layer(lay)%node(ii),qdum, &
        & chrg=cc,uhf=env%uhf,pr=.false.,prch=stdout)
        qattmp(:,cc2) = qdum(:)
        do jj = 1,nnat
          kk = heap%layer(lay)%origin(ii)%map(jj)
          if (kk > 0) then
            qattmpref(jj) = qat(kk)
          else
            qattmp(jj,cc2) = 0.0_wp
          end if
          qtmp(cc2) = qtmp(cc2)+abs(qattmp(jj,cc2)-qattmpref(jj))
        end do
      end do
      !> select best charge
      cc = minloc(qtmp,1)
      lq(ii) = ichrgs(cc)
      deallocate (qattmp,qattmpref)
      !write (*,*) 'charge MAEs on frag:',qtmp
      !write (*,*) 'selected charge:',lq(ii)
    end do
    write (stdout,'(2x,a)',advance='no') 'determined charges:'
    write (stdout,*) lq
  end subroutine sum_charges_layer
end subroutine crest_queue_setup

!=============================================================================!
!#############################################################################!
!=============================================================================!

subroutine crest_queue_iter(env,iterate)
  use crest_parameters
  use crest_data
  use strucrd
  use iomod
  use crest_calculator
  implicit none
  type(systemdata),intent(inout),target :: env
  logical,intent(out) :: iterate
  integer :: ii,jj,kk,io,nn,ll,lll,ati,atj
  type(coord) :: mol
  character(len=10) :: atmp
  character(len=*),parameter :: dirname = 'crest_queue_'

  iterate = .false.

  if (allocated(env%splitqueue).and.env%splitheap%nqueue > 0) then
!>--- important restoring to initial calc/dir
    env%calc => env%splitheap%origincalc
    call chdir(env%splitheap%origindir)

    !> next iter
    ii = env%queue_iter+1
    env%queue_iter = ii

    write (stdout,'(/,70("§"))')
    write (stdout,'(a,i0)') "§§§   QUEUE ITERATION ",ii
    write (stdout,'(70("§"))')

    jj = env%splitheap%queue(ii)%layer
    kk = env%splitheap%queue(ii)%node
    associate (heap => env%splitheap,queue => env%splitheap%queue(ii))

      !> create a dedicated work directory
      write (atmp,'(i0)') ii
      queue%workdir = dirname//trim(atmp)
      io = makedir(queue%workdir)
      call chdir(queue%workdir)
      write (stdout,'(a,t28,a,t30,a)') 'Queue work (sub-)directory',':', &
        & trim(queue%workdir)

      !> selecting output file depending on runtype
      select case (env%crestver)
      case (crest_imtd,crest_imtd2)
        queue%file = 'crest_ensemble.xyz'
      case (crest_optimize)
        queue%file = 'crestopt.xyz'
      case (crest_moldyn)
        queue%file = 'crest_dynamics.trj.xyz'
      case (crest_bh)
        queue%file = 'crest_quenched.xyz'
      case default
        queue%file = 'struc.xyz'
      end select
      write (stdout,'(a,t28,a,t30,a)') 'Selected output file',':',queue%file

!>--- new calculator setup section and env update
      call queue%calc%copy(env%calc,ignore_constraints=.true.)
      !> for constraints we must be careful and map them to the new order
      call update_constraints_queue(heap,jj,kk,env%calc,queue%calc)

      mol = env%splitheap%layer(jj)%node(kk)
      call env%ref%load(mol)
      call mol%write('coord')
      call queue%calc%set_charge(mol%chrg) !> the nodes may have different charges saved
      call queue%calc%info(stdout)

      if (allocated(env%ref%wbo)) deallocate (env%ref%wbo)
      env%nat = mol%nat
      env%rednat = mol%nat
      env%chrg = mol%chrg
      if (.not.env%user_mdtime) then
        env%mdtime = -1.0_wp
        env%mddat%length_ps = -1.0_wp
      end if

      env%calc => queue%calc

    end associate
    if (ii < env%splitheap%nqueue) then
      iterate = .true.
    end if

    write (stdout,*)
  end if

contains
  subroutine update_constraints_queue(heap,layer,node,refcalc,newcalc)
    use molbuilder_construct_list
    implicit none
    type(construct_heap) :: heap
    integer :: layer,node
    type(calcdata),intent(in) :: refcalc
    type(calcdata),intent(inout) :: newcalc
    integer :: nn,ll,lll,ati,atj,nn2
    type(constraint),allocatable :: cons(:)
    if (refcalc%nconstraints > 0) then
      nn = refcalc%nconstraints
      allocate (cons(nn))
      do ll = 1,nn
        call cons(ll)%copy(refcalc%cons(ll))
        do lll = 1,cons(ll)%n
          ati = cons(ll)%atms(lll)
          call heap%find_current_position(ati,layer,node,atj)
          cons(ll)%atms(lll) = atj !> overwrite with the current position
        end do
        if (any(cons(ll)%atms(:) .eq. 0)) then
          cons(ll)%active = .false.
        end if
      end do
      !> clean (active) constraints
      nn2 = 0
      do ll = 1,nn
        if (cons(ll)%active) nn2 = nn2+1
      end do
      if (nn2 > 0) then
        newcalc%nconstraints = nn2
        allocate (newcalc%cons(nn2))
        lll = 0
        do ll = 1,nn
          if (cons(ll)%active) then
            lll = lll+1
            call newcalc%cons(lll)%copy(cons(ll))
          end if
        end do
      end if
    end if
  end subroutine update_constraints_queue
end subroutine crest_queue_iter

subroutine crest_queue_iter_resort(env,iterate)
  use crest_parameters
  use crest_data
  use iomod
  use cregen_interface
  implicit none
  type(systemdata),intent(inout) :: env
  logical,intent(in) :: iterate

  character(len=:),allocatable :: file
  logical :: heavytmp,confgotmp,ex

  if (.not. (allocated(env%splitqueue).and.env%splitheap%nqueue > 0)) return

  select case (env%crestver)
  case (crest_imtd,crest_imtd2)

    write (stdout,'(/,75("*"))')
    write (stdout,'(a,i0)') "***  CREGEN heavy-atom resorting for QUEUE iteration ",env%queue_iter
    write (stdout,'(75("*"))')
    ex = .false.
    if (file_exists(crefile//'.xyz')) then
      ex = .true.
      file = crefile//'.xyz'
    else if (file_exists(conformerfile)) then
      ex = .true.
      file = conformerfile
    end if
    heavytmp = env%heavyrmsd
    confgotmp = env%confgo
    env%heavyrmsd = .true.
    env%confgo = .true.
    call newcregen(env,infile=file)
    env%heavyrmsd = heavytmp
    env%confgo = confgotmp
    if (file_exists(file//'.sorted')) then
      call rename(file//'.sorted',ensemblefile)
    end if
  case default
  end select

end subroutine crest_queue_iter_resort

!=============================================================================!
!#############################################################################!
!=============================================================================!

subroutine crest_queue_reconstruct(env,tim)
  use crest_parameters
  use crest_data
  use molbuilder_construct_list
  use molbuilder_construct_mod
  use strucrd
  use iomod
  use crest_calculator
  use utilities,only:checkname_xyz
  use term_ui,only:progress_init,progress_update,progress_finish
  implicit none
  type(systemdata),intent(inout) ::  env
  type(timer),intent(inout) :: tim
  type(coord) :: mol
  integer :: ii,jj,kk,nall
  logical :: ex,multilevel(6)
  type(timer) :: timtmp
  type(coord),allocatable :: structures(:)
  type(calcdata),target :: newcalc
  character(len=256) :: inpnam,outnam
  character(len=*),parameter :: recfile = 'crest_reconstruct.xyz'

  if (.not. (allocated(env%splitqueue).and.env%splitheap%nqueue > 0)) then
    return
  end if

  call tim%start(9,'Queue reconstruction')

  write (stdout,'(/,80("#"))')
  write (stdout,'(3("#"),t25,a,t78,3("#"))') 'QUEUE STRUCTURE RECONSTRUCTION'
  write (stdout,'(80("#"),/)')

  !> reset
  mol = env%splitheap%originmol
  call env%ref%load(mol)
  env%nat = mol%nat
  env%rednat = mol%nat
  env%chrg = mol%chrg
  env%calc => env%splitheap%origincalc
  call chdir(env%splitheap%origindir)

  call env%splitheap%fill_inverse_depth()
  call recusrive_construct(env,env%splitheap,1)
  nall = env%splitheap%layer(1)%nmols
  allocate (structures(nall))
  do ii = 1,nall
    structures(ii) = env%splitheap%layer(1)%mols(ii)
  end do
  !deallocate (env%splitheap%layer(1)%mols)
  deallocate (env%splitheap%layer)
  deallocate (env%splitheap%queue)

  write (stdout,'(/,1x,a)') 'Writing reconstructed structures to: "'//recfile//'"'
  call wrensemble(recfile,nall,structures)
  write (stdout,*)

  call newcalc%copy(env%calc)
  env%calc => newcalc
  call env%calc%info(stdout)

  select case (env%crestver)
  case (crest_optimize)
    call env%ref%load(structures(1))
    call crest_optimization(env,timtmp)
  case default
    call optlev_to_multilev(env%optlev,multilevel)
    call crest_multilevel_oloop(env,recfile,multilevel,0)
    if (env%iostatus_meta .ne. 0) return

    call smallheadline('FINAL GEOMETRY OPTIMIZATION IN QUEUE RECONSTRUCTION')
    call checkname_xyz(crefile,inpnam,outnam)
    call rename(inpnam,recfile)
    call crest_multilevel_wrap(env,recfile,0)

    call V2terminating()
  end select

  if (.not.env%keepmodef) then
    call rmrf('crest_queue_*')
  end if

  call tim%stop(9)

contains
  recursive subroutine recusrive_construct(env,heap,targetlayer)
    use irmsd_module,only:irmsd,rmsd,rmsd_cache,rmsd_core_cache,min_rmsd
    use canonical_mod
    use omp_lib
    implicit none
    type(systemdata),intent(inout) :: env
    type(construct_heap),intent(inout) :: heap
    integer,intent(in) :: targetlayer

    integer :: ii,jj,kk
    integer :: queuepos
    character(len=:),allocatable :: basefile,sidefile
    type(coord),allocatable :: structures_b(:)
    type(coord),allocatable :: structures_s(:)
    type(coord) :: mol,moltmp
    type(coord),allocatable :: moltmp_arr(:)
    integer :: nall_b,nall_s,id_b,id_s
    integer :: rr,io,rg,nregions,max_structs
    integer :: reg_blo(3),reg_bhi(3),reg_slo(3),reg_shi(3)
    integer :: target_bhi,target_shi
    integer :: outer_lo,outer_hi,inner_lo,inner_hi,outer_idx,inner_idx
    logical :: base_is_outer
    integer :: duplicates
    logical :: ex,clash,duplicate
    real(wp) :: RTHR,rmsval,ETHR,deltaE,depthlimit
    real(wp) :: layerfactor_b,layerfactor_s,weight_s,weight_b
    type(rmsd_cache),allocatable :: rcache(:)
    type(rmsd_core_cache),allocatable :: ccache(:)
    type(canonical_sorter) :: canref
    real(wp),allocatable :: xyzscratch(:,:,:,:)
    logical,allocatable :: mask(:)
    integer :: T,Tn,tt
    type(timer) :: profiler

    character(len=*),parameter :: subdir_tmp = 'crest_queue_'
    character(len=:),allocatable :: subdirfile
    character(len=10) :: atmp
    character(len=60) :: btmp

    associate (layer => heap%layer(targetlayer))
      if (layer%nnodes > 2) then
        write (stdout,'(a)') 'currently unhandled edge-case in layer reconstruction:'
        write (stdout,'(a,i0,a)') 'layer ',targetlayer,' was split in more than 2 structures'
        stop
      end if

      layer%inverse_depth = layer%inverse_depth+1.0_wp
      do ii = 1,layer%nnodes
        if (allocated(layer%childlayer)) then
          jj = layer%childlayer(ii)
        else
          jj = 0
        end if
        if (jj == 0.and.ii == 1) then

          do kk = 1,heap%nqueue
            if (heap%queue(kk)%layer == targetlayer.and.heap%queue(kk)%node == ii) then
              basefile = heap%queue(kk)%file
              id_b = kk
            end if
          end do

          write (atmp,'(i0)') id_b
          subdirfile = subdir_tmp//trim(atmp)//'/'//basefile
          inquire (exist=ex,file=subdirfile)
          if (ex) then
            write (stdout,'(a,t26,a,t30,a)',advance='no') &
              & 'Reading fragment(s) from',':',subdirfile
            call rdensemble(subdirfile,nall_b,structures_b)
            write (stdout,'(1x,a,i0,a)') '--> ',nall_b,' structure(s)'

          end if
          layerfactor_b = 1.0_wp

        else if (jj == 0.and.ii == 2) then

          do kk = 1,heap%nqueue
            if (heap%queue(kk)%layer == targetlayer.and.heap%queue(kk)%node == ii) then
              sidefile = heap%queue(kk)%file
              id_s = kk
            end if
          end do

          write (atmp,'(i0)') id_s
          subdirfile = subdir_tmp//trim(atmp)//'/'//sidefile
          inquire (exist=ex,file=subdirfile)
          if (ex) then
            write (stdout,'(a,t26,a,t30,a)',advance='no') &
              & 'Reading fragment(s) from',':',subdirfile
            call rdensemble(subdirfile,nall_s,structures_s)
            write (stdout,'(1x,a,i0,a)') '--> ',nall_s,' structure(s)'
          end if
          layerfactor_s = 1.0_wp

        else
          call recusrive_construct(env,heap,jj)
          if (ii == 1) then
            nall_b = heap%layer(jj)%nmols
            allocate (structures_b(nall_b))
            do kk = 1,nall_b
              structures_b(kk) = heap%layer(jj)%mols(kk)
            end do
            layerfactor_b = heap%layer(jj)%inverse_depth

          else if (ii == 2) then

            nall_s = heap%layer(jj)%nmols
            allocate (structures_s(nall_s))
            do kk = 1,nall_s
              structures_s(kk) = heap%layer(jj)%mols(kk)
            end do
            layerfactor_s = heap%layer(jj)%inverse_depth
            !deallocate (heap%layer(jj)%mols)
          end if
        end if
      end do
      weight_s = layerfactor_s/(layerfactor_s+layerfactor_b)
      weight_b = layerfactor_b/(layerfactor_s+layerfactor_b)

      write (stdout,*)
      write (stdout,'(a,i0)') 'Reconstructing layer : ',targetlayer
      write (stdout,'(2x,a,i0)') 'Base structures         : ',nall_b
      write (stdout,'(2x,a,i0)') 'Side chain structures   : ',nall_s
      write (stdout,'(2x,a,es9.2)') 'Max. combinations       : ',real(nall_b,wp)*real(nall_s,wp)
      write (stdout,'(2x,a,f7.5,a)') 'Similarity threshold    : ',env%rthr,' Å'
      write (stdout,'(2x,a,f7.5,a)') 'ΔE threshold (ETHR)     : ',env%ethr,' kcal/mol'

      layer%nmols = 0
      depthlimit = real(env%queue_maxreconstruct,wp)*(env%queue_depthfac**real(targetlayer-1,wp))
      max_structs = nint(min(real(nall_b,wp)*real(nall_s,wp),depthlimit))
      allocate (layer%mols(max_structs))
      write (stdout,'(2x,a,i0)') 'Capping limit           : ',env%queue_maxreconstruct
      write (stdout,'(2x,a,f4.2,a)') 'Depth factor            : ',env%queue_depthfac,'^(layer-1)'
      write (stdout,'(2x,a,i0)') 'Max. new structs stored : ',max_structs

      RTHR = env%rthr*aatoau   !> RMSD threshold in Bohr
      ETHR = env%ethr/autokcal !> deltaE threshold in hartree
      duplicates = 0
      T = 1
      call new_ompautoset(env,'max',max_structs,T,Tn)
      write (stdout,'(2x,a,i0)') 'OpenMP threads          : ',T
      allocate (ccache(T))
      allocate (rcache(T))
      allocate (moltmp_arr(T))
      allocate (mask(layer%refmol%nat),source=.true.)
      call canref%init(layer%refmol,invtype='apsp+',heavy=.false.)

      do tt = 1,T
        call ccache(tt)%allocate(layer%refmol%nat,scratch=.true.)
        call rcache(tt)%allocate(layer%refmol%nat)
        rcache(tt)%stereocheck = .not. (canref%hasstereo(layer%refmol))
        rcache(tt)%rank(:,1) = canref%rank(:)
        rcache(tt)%rank(:,2) = canref%rank(:)
      end do
      do ii = 1,layer%refmol%nat
        if (layer%refmol%at(ii) == 1) mask(ii) = .false.
      end do
!      write (stdout,'(2x,a)') 'Recombining under heavy-atom RMSD consideration (this may take a while) ... '
      write (stdout,'(2x,a)') 'Recombining under iRMSD consideration (this may take a while) ... '
      call progress_init(env%ps,max_structs,width=50,prefix=" ↳ ", &
        &                suffix="",show_time=.true.,show_eta=.false.)
      call progress_update(env%ps,0,max_structs)

      call profiler%init(1)
      call profiler%start(1)

      ! ── Precompute sampling regions ──────────────────────────────
      !> Region 1 targets max_structs combinations in the correct
      !> weight ratio.  Regions 2–3 expand into remaining structures.
      base_is_outer = (nall_b <= nall_s)
      nregions = 0

      target_bhi = nint(sqrt(real(max_structs,wp)*weight_b/weight_s))
      target_shi = nint(sqrt(real(max_structs,wp)*weight_s/weight_b))

      reg_blo(1) = 1
      reg_slo(1) = 1
      reg_bhi(1) = min(nall_b,target_bhi)
      reg_shi(1) = min(nall_s,target_shi)
      !> reciprocal fill if one dimension was capped
      if (reg_bhi(1) < target_bhi.and.reg_bhi(1) > 0) then
        reg_shi(1) = min(nall_s,nint(real(max_structs,wp)/real(reg_bhi(1),wp)))
      else if (reg_shi(1) < target_shi.and.reg_shi(1) > 0) then
        reg_bhi(1) = min(nall_b,nint(real(max_structs,wp)/real(reg_shi(1),wp)))
      end if
      nregions = 1

      !> Region 2: expand whichever dimension wasn't exhausted
      if (reg_shi(1) == nall_s.and.reg_bhi(1) < nall_b) then
        nregions = 2
        reg_blo(2) = reg_bhi(1)+1
        reg_bhi(2) = nall_b
        reg_slo(2) = 1
        reg_shi(2) = nall_s
      else if (reg_bhi(1) == nall_b.and.reg_shi(1) < nall_s) then
        nregions = 2
        reg_blo(2) = 1
        reg_bhi(2) = nall_b
        reg_slo(2) = reg_shi(1)+1
        reg_shi(2) = nall_s
      else if (reg_bhi(1) < nall_b.and.reg_shi(1) < nall_s) then
        !> Neither exhausted: expand larger dim first, then the other
        nregions = 3
        if (base_is_outer) then
          reg_blo(2) = 1
          reg_bhi(2) = reg_bhi(1)
          reg_slo(2) = reg_shi(1)+1
          reg_shi(2) = nall_s
          reg_blo(3) = reg_bhi(1)+1
          reg_bhi(3) = nall_b
          reg_slo(3) = 1
          reg_shi(3) = nall_s
        else
          reg_blo(2) = reg_bhi(1)+1
          reg_bhi(2) = nall_b
          reg_slo(2) = 1
          reg_shi(2) = reg_shi(1)
          reg_blo(3) = 1
          reg_bhi(3) = nall_b
          reg_slo(3) = reg_shi(1)+1
          reg_shi(3) = nall_s
        end if
      end if

      ! ── Reconstruct by iterating over regions ───────────────────
      regionloop: do rg = 1,nregions
        if (base_is_outer) then
          outer_lo = reg_blo(rg); outer_hi = reg_bhi(rg)
          inner_lo = reg_slo(rg); inner_hi = reg_shi(rg)
        else
          outer_lo = reg_slo(rg); outer_hi = reg_shi(rg)
          inner_lo = reg_blo(rg); inner_hi = reg_bhi(rg)
        end if
        do outer_idx = outer_lo,outer_hi
          do inner_idx = inner_lo,inner_hi
            if (base_is_outer) then
              ii = outer_idx; jj = inner_idx
            else
              ii = inner_idx; jj = outer_idx
            end if

            call attach(structures_b(ii),structures_s(jj),layer%alignmap,mol, &
            & remove_lastx=layer%ncapped,original_map=layer%position_mapping, &
            & clash=clash,reficn=layer%reficn)
            mol%energy = structures_b(ii)%energy+structures_s(jj)%energy
            if (.not.clash) then
              duplicate = .false.

              !$omp parallel &
              !$omp shared(duplicate,duplicates,mol,ccache,rcache,mask,ETHR,moltmp_arr) &
              !$omp private(rr,tt,deltaE,rmsval)
              !$omp do schedule(dynamic)
              do rr = 1,layer%nmols
                if (duplicate) cycle
                tt = omp_get_thread_num()+1
                deltaE = abs(mol%energy-layer%mols(rr)%energy)
                if (deltaE < ETHR) then
                  call moltmp_arr(tt)%copy(layer%mols(rr))
                  call min_rmsd(mol,moltmp_arr(tt),rcache=rcache(tt),rmsdout=rmsval,align=.false.)
                  !$omp critical
                  if (rmsval < RTHR.and..not.duplicate) then
                    duplicate = .true.
                    duplicates = duplicates+1
                  end if
                  !$omp end critical
                end if
              end do
              !$omp end do
              !$omp end parallel

              if (.not.duplicate) then
                layer%nmols = layer%nmols+1
                layer%mols(layer%nmols) = mol
                call progress_update(env%ps,layer%nmols,max_structs)
                if (layer%nmols == max_structs) exit regionloop
              end if
            end if
          end do
        end do
      end do regionloop
      if (layer%nmols < max_structs) then
        call progress_update(env%ps,1,1)
      end if
      call progress_finish(env%ps)
      write (stdout,'(2x,a)') 'done!'
      if (duplicates > 0) then
        write (stdout,'(2x,a,i0)') 'Avoided duplicates       : ',duplicates
      end if
      write (stdout,'(2x,a,i0)') 'Successful combinations  : ',layer%nmols
      call profiler%stop(1)
      write (btmp,'(2x,a)') 'Total runtime for recombination step:'
      call profiler%write_timing(stdout,1,trim(btmp),.true.)
      write (stdout,*)

    end associate
  end subroutine recusrive_construct

end subroutine crest_queue_reconstruct

