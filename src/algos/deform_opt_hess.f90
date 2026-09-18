subroutine deform_opt_hess(calc,mol)  
    use crest_calculator
    use strucrd
    use irmsd_module
    use bh_step_module
    use crest_parameters
    use optimize_module
    use thermochem_module
    use hr_utils
    use optimize_maths
    implicit none
    type(calcdata) :: calc
    type(coord) :: mol
    type(coord) :: molnew,mol_reopt
    real(wp) :: energy,stepsize,rmsdval
    real(wp) :: grad(3,mol%nat) 
    logical :: pr,wr 
    integer :: io, nat3,idx

    real(wp) :: etot 
    real(wp), allocatable :: h_init(:,:)

    if (allocated(calc%chess)) deallocate(calc%chess)

    !allocate (calc%chess)
    !call calc%chess%alloc(mol%nat,calc%hu_steps,calc%hguess,calc%initialize_hr_type, calc%hr_hu_type) !Maybe in future just reset the cash here, or if this works we only call this here

    nat3 = 3*mol%nat
    molnew%nat = mol%nat
    molnew%at = mol%at
    molnew%xyz = mol%xyz

    allocate(h_init(nat3,nat3))

    stepsize = calc%doh_stepsize


    call take_fixed_stepsize_cart(molnew,stepsize,calc)

    pr = .true.
    wr = .true.

    call optimize_geometry(molnew,mol_reopt,calc,energy,grad,pr,wr,io)

    pr = .true.
    wr = .true.

    rmsdval = rmsd(mol,mol_reopt)

    write(stdout,*) 'VALUE OF RMSD FOR REOPTIMISED STRUCTURE',rmsdval

    if (rmsdval .le. 0.1_wp) then
    
        idx = minloc(calc%chess%order,1)
        if (minval(calc%chess%order) .eq. 0) idx = 1

        write(stdout,*)calc%chess%coords(idx,:,:) - molnew%xyz(:,:)
        
        call initialize_hessian(calc,calc%chess%initialize_type,molnew%xyz,molnew%nat,molnew%at,calc%chess%hess(:),calc%chess%hguess,pr)  !> This hguess is set through the hguess variable of the optimizer and needs to be hardcoded/set explicitly before initialization for benchmarking!!
        call dhtosq(nat3,H_init,calc%chess%hess) !> maybe this should all be inside the construct bfgs function later? -> cannot due to circular import!!!
        write(stdout,*)                                                                                                   !> Hessian type (gfnff,mod,identity) is set through input file and is already encoded into the calc object
        write(stdout,*)"THERMO FROM INITIALIZED HESSIAN:"
        write(stdout,*) 
        call calc_thermo_from_hess(molnew,H_init,pr, &
        & calc%nt,calc%temperatures,calc%ithr,calc%fscal,calc%sthr,calc%et, &
        & calc%ht,calc%gt,calc%stot,etot)

        call calc%chess%construct_hessian()

        write (stdout,*)
        write (stdout,*) "THERMO FROM RECONSTRUCTED HESSIAN:" 
        write (stdout,*)

        call calc_thermo_from_hess(molnew,calc%chess%H(:,:),pr, &
        & calc%nt,calc%temperatures,calc%ithr,calc%fscal,calc%sthr,calc%et, &
        & calc%ht,calc%gt,calc%stot,etot)

    else 
        write(stdout,*) "Reoptimised Geometry not equal to initial structure"
    
    endif


end subroutine deform_opt_hess