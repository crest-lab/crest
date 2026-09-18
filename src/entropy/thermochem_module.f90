module thermochem_module
  use crest_parameters
  use symmetry_i
  use optimize_maths
  use atmasses,only:molweight
  use iomod,only:to_lower,directory_exist
  use axis_module
  use strucrd
  use crest_thermo
  implicit none
  private

  public :: frequencies
  public :: effective_hessian
  public :: prj_mw_hess,mass_weight_hess
  public :: calcthermo,calc_thermo_from_hess
  public :: print_vib_spectrum,print_hessian,print_g98_fake

!=============================================================================!
contains  !> MODULE PROCEDURES STARTE HERE
!=============================================================================!

  subroutine frequencies(nat,at,xyz,nat3,prj_mw_hess,freq,io)
!*************************************************
!* Returns the Frequencies from a Hessian in cm-1
!*************************************************
    implicit none

    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    real(wp) :: prj_mw_hess(nat3,nat3)

    integer :: io,nat3
    logical :: pr
    real(wp) :: energy
    real(wp) :: freq(nat3)
    real(wp),allocatable :: pmode(:,:)

    integer,allocatable :: iwork(:)
    real(wp),allocatable :: work(:)

    integer :: lwork,liwork,info,i
    integer :: unit
    !>LAPCK
    external :: dsyevd

    io = 0
    nat3 = nat*3

    !Parameters for diagonalization
    lwork = 1+6*nat3+2*nat3**2
    liwork = 3+5*nat3

    !$omp critical
    allocate (work(lwork),iwork(liwork))
    !$omp end critical

    !Diagonalization
    call dsyevd('V','U',nat3,prj_mw_hess,nat3,freq,work,lwork,iwork,liwork,info)
    io = info
    deallocate (work,iwork)

    !Convert eigenvalues to frequencies
    do i = 1,nat3
      if (freq(i) .gt. 0.0_wp) then
        freq(i) = sqrt(freq(i))*autorcm
      else
        freq(i) = -sqrt(abs(freq(i)))*autorcm
      end if
    end do

    !open (newunit=unit,file="frequencies")
    !write (unit,*) "Frequencies:"
    !do i = 1,size(freq)
    !  write (unit,*) freq(i)
    !end do
    !close (unit)

    return

  end subroutine frequencies

  subroutine mass_weight_hess(nat,at,nat3,hess)
    implicit none

    !> Mass weighting the Hessian
    integer,intent(in) :: nat                   !Number of atoms
    integer,intent(in) :: at(nat)               !atomic number of all atoms

    real(wp),intent(inout) :: hess(nat3,nat3)   !Hessian matrix
    real(wp) :: mass_in_au             !Masses of all atoms of the periodic table
    integer :: i,j,nat3,i3,i33,j3,j33

    !mass_in_au = (1.66054e-27_wp/9.1094e-31_wp)**2
    mass_in_au = (amutokg/metokg)**2

    do i = 1,nat
      do j = i,nat

        i3 = 3*(i-1)+1
        i33 = 3*(i-1)+3
        j3 = 3*(j-1)+1
        j33 = 3*(j-1)+3

        hess(i3:i33,j3:j33) = 1/sqrt(ams(at(i))*ams(at(j))*mass_in_au)*hess(i3:i33,j3:j33)
        !Hessian is symmetric hence upper triangular can be copied
        hess(j3:j33,i3:i33) = hess(i3:i33,j3:j33)

      end do
    end do

    return
  end subroutine mass_weight_hess

!=========================================================================================!

  subroutine prj_mw_hess(nat,at,nat3,xyz,hess,phess_ut)
!***************************************************************
!* Projection of the translational and rotational DOF out of
!* the numerical Hessian plus the mass-weighting of the Hessian
!***************************************************************
    implicit none

    integer,intent(in) :: nat,nat3
    integer,intent(in) :: at(nat)
    real(wp),intent(inout) :: hess(nat3,nat3)
    real(wp),intent(in) ::  xyz(3,nat)
    real(wp),intent(in),optional,target :: phess_ut(:)
    !real(wp) ::  hess_ut(nat3*(nat3+1)/2),pmode(nat3,1)
    real(wp),allocatable,target :: hess_ut(:)
    real(wp),allocatable :: pmode(:,:)
    real(wp),pointer :: phess(:)
    integer :: i

    if (present(phess_ut)) then
      phess => phess_ut
    else
      !$omp critical
      allocate (hess_ut(nat3*(nat3+1)/2),source=0.0_wp)
      phess => hess_ut
      !$omp end critical
    end if
    !$omp critical
    allocate (pmode(nat3,1),source=0.0_wp)
    !$omp end critical

    !> Transforms matrix of the upper triangle vector
    call dsqtoh(nat3,hess,phess)

    !> Projection
    call trproj(nat,nat3,xyz,phess,.false.,0,pmode,1)

    !> Transforms vector of the upper triangle into matrix
    call dhtosq(nat3,hess,phess)

    !> Mass weighting
    call mass_weight_hess(nat,at,nat3,hess)

  end subroutine prj_mw_hess

  !============================================================================!
  !############################################################################!
  !============================================================================!

  subroutine prepthermo(nat,at,xyz,pr,molmass,rabc,avmom,symnum,symchar,iunit)
!***********************************************************************
!* Prepare the calculation of thermodynamic properties of a structure
!* In particular, determine rotational constants and check the symmetry
!***********************************************************************
    implicit none
    integer,intent(in)     :: nat
    integer,intent(in)     :: at(nat)
    real(wp),intent(inout) :: xyz(3,nat) !> in Angstroem
    logical,intent(in)     :: pr
    real(wp),intent(out)   :: molmass
    real(wp),intent(inout) :: rabc(3)
    real(wp),intent(out)   :: avmom
    real(wp),intent(out)   :: symnum
    integer,intent(in)     :: iunit

    real(wp) :: a,b,c
    character(len=4) :: sfsym
    character(len=3) :: sym,symchar
    !>--- molecular mass in amu
    molmass = molweight(nat,at)

    if (pr) then
      write (iunit,'(1x,a,f15.2)') 'Mol. weight /amu  : ',molmass
    end if

    !>--- rotational constants in cm-1
    rabc = 0.0d0
    call axis(nat,at,xyz,rabc(1:3),avmom)
    a = rabc(3)
    b = rabc(2)
    c = rabc(1)
    rabc(1) = a
    rabc(3) = c
    if (pr) then
      write (iunit,'(1x,a,3f15.2)') 'Rot. const. /MHz  : ',rabc(1:3)
    end if
    !rabc = rabc/2.99792458d+4   ! MHz to cm-1
    rabc = rabc*mhztorcm
    if (pr) then
      write (iunit,'(1x,a,3f15.2)') 'Rot. const. /cm-1 : ',rabc(1:3)
    end if

    !>--- symmetry number from rotational symmetry
    xyz = xyz/bohr
    !$omp critical
    call getsym(.false.,iunit,nat,at,xyz,sfsym)
    !$omp end critical
    xyz = xyz*bohr
    sym = sfsym(1:3)
    symchar = sym
    symnum = 1.0d0
    if (a .lt. 1.d-9.or.b .lt. 1.d-9.or.c .lt. 1.d-9) then
      if (index(sym,'d') .ne. 0) symnum = 2.0d0
    else
      call to_lower(sym)
      if (index(sym,'c2') .ne. 0) symnum = 2.0d0
      if (index(sym,'s4') .ne. 0) symnum = 2.0d0
      if (index(sym,'c3') .ne. 0) symnum = 3.0d0
      if (index(sym,'s6') .ne. 0) symnum = 3.0d0
      if (index(sym,'c4') .ne. 0) symnum = 4.0d0
      if (index(sym,'s8') .ne. 0) symnum = 4.0d0
      if (index(sym,'c5') .ne. 0) symnum = 5.0d0
      if (index(sym,'c6') .ne. 0) symnum = 6.0d0
      if (index(sym,'c7') .ne. 0) symnum = 7.0d0
      if (index(sym,'c8') .ne. 0) symnum = 8.0d0
      if (index(sym,'c9') .ne. 0) symnum = 9.0d0
      if (index(sym,'d2') .ne. 0) symnum = 4.0d0
      if (index(sym,'d3') .ne. 0) symnum = 6.0d0
      if (index(sym,'d4') .ne. 0) symnum = 8.0d0
      if (index(sym,'d5') .ne. 0) symnum = 10.0d0
      if (index(sym,'d6') .ne. 0) symnum = 12.0d0
      if (index(sym,'d7') .ne. 0) symnum = 14.0d0
      if (index(sym,'d8') .ne. 0) symnum = 16.0d0
      if (index(sym,'d9') .ne. 0) symnum = 18.0d0
      if (index(sym,'t') .ne. 0) symnum = 12.0d0
      if (index(sym,'td') .ne. 0) symnum = 12.0d0
      if (index(sym,'th') .ne. 0) symnum = 12.0d0
      if (index(sym,'o') .ne. 0) symnum = 24.0d0
      if (index(sym,'oh') .ne. 0) symnum = 24.0d0
      if (index(sym,'ih') .ne. 0) symnum = 60.0d0
    end if

    if (pr) then
      write (iunit,'(1x,a,4x,a)') 'Symmetry:',sym
    end if
    return
  end subroutine prepthermo

  subroutine calcthermo(nat,at,xyz,freq,pr,ithr,fscal,sthr,nt,temps, &
      &      et,ht,gt,stot,iunit_in,emodel)
!**************************************************************
!* Calculate thermodynamic contributions for a given structure
!* from it's frequencies (from second derivatives/the Hessian)
!* Based on xtb's "print_thermo" routine
!**************************************************************

    implicit none
    integer,intent(in)     :: nat
    integer,intent(in)     :: at(nat)
    real(wp),intent(inout) :: xyz(3,nat)  !in Bohr
    real(wp),intent(inout) :: freq(3*nat) !in cm-1
    logical,intent(in)     :: pr
    real(wp),intent(in) :: ithr     !imag. inv. in cm-1
    real(wp),intent(in) :: fscal    !freq scaling
    real(wp),intent(in) :: sthr     !rotor cut
    integer,intent(in)  :: nt
    real(wp),intent(in) :: temps(nt)
    integer,intent(in),optional     :: iunit_in
    character(len=*),intent(in),optional :: emodel
    real(wp) :: et(nt)          !< enthalpy in Eh
    real(wp) :: ht(nt)          !< enthalpy in Eh
    real(wp) :: gt(nt)          !< free energy in Eh
    real(wp) :: stot(nt)        !< entropy in cal/molK
    real(wp) :: ts(nt)          !< entropy*T in Eh
    real(wp) :: rabc(3),a,b,c
    real(wp) :: avmom
    real(wp) :: molmass
    real(wp) :: sym
    real(wp) :: zp
    character(len=3) :: symchar
    logical :: pr2
    logical :: linear = .false.
    logical :: atom = .false.
    integer :: nvib_theo
    integer :: nvib,nimag
    real(wp) :: vibthr
    real(wp),allocatable :: vibs(:)

    integer :: i,j,iunit,emodelunit
    integer :: n3,rt
    real(wp) :: adum(nt)
    character(len=64) :: atmp

    character(len=*),parameter :: outfmt = &
    &  '(9x,"::",1x,a,f24.12,1x,a,1x,"::")'
    character(len=*),parameter :: dblfmt = &
    &  '(10x,":",2x,a,f24.7,1x,a,t63,":")'
    character(len=*),parameter :: intfmt = &
    &  '(10x,":",2x,a,i24,       t63,":")'
    character(len=*),parameter :: chrfmt = &
    &  '(10x,":",2x,a,a24,       t63,":")'

    !real(wp),parameter :: autorcm = 219474.63067_wp
    !real(wp),parameter :: rcmtoau = 1.0_wp/autorcm
    real(wp),parameter :: autocal = autokcal*1000.0_wp

    xyz = xyz*autoaa  !> NOTE: FROM HERE ON WE WORK IN ANGSTRÖM

    if (present(iunit_in)) then
      iunit = iunit_in
    else
      iunit = stdout
    end if

    if (present(emodel)) then
      select case (emodel)
      case ('grimme')
        emodelunit = 1
      case ('truhlar')
        emodelunit = 2
      case default
        emodelunit = 1
      end select
    else
      emodelunit = 1
    end if

    call prepthermo(nat,at,xyz,pr,molmass,rabc,avmom,sym,symchar,iunit)

    n3 = 3*nat
    !$omp critical
    allocate (vibs(n3))
    !$omp end critical
    vibthr = 1.0
    a = rabc(1)
    b = rabc(2)
    c = rabc(3)

    nvib_theo = 3*nat-6
    if (c .lt. 1.d-10.or.(symchar == 'din')) linear = .true.
    if (linear) nvib_theo = 3*nat-5

    if (a+b+c .lt. 1.d-6) then
      atom = .true.
      nvib = 0
      nvib_theo = 0
    end if

    nvib = 0
    vibs = 0.0
    do i = 1,n3
      if (abs(freq(i)) .gt. vibthr) then
        nvib = nvib+1
        vibs(nvib) = freq(i)
      end if
    end do
    !> scale
    vibs(1:nvib) = vibs(1:nvib)*fscal

    !> invert imaginary modes
    nimag = 0
    do i = 1,nvib
      if (vibs(i) .lt. 0.and.vibs(i) .gt. ithr) then
        vibs(i) = -vibs(i)
        if (pr) write (stdout,'(a,i5," :",f10.2)') 'Inverting frequency',i,vibs(i)
      end if
      if (vibs(i) < 0.0) then
        nimag = nimag+1
      end if
    end do

    if (pr) then
      write (iunit,'(a)')
      write (iunit,'(10x,53("."))')
      write (iunit,'(10x,":",23x,a,23x,":")') "SETUP"
      write (iunit,'(10x,":",51("."),":")')
      write (iunit,intfmt) "# frequencies    ",nvib
      write (iunit,intfmt) "# imaginary freq.",nimag
      write (atmp,*) linear
      write (iunit,chrfmt) "linear?          ",trim(atmp)
      write (iunit,chrfmt) "symmetry         ",adjustr(symchar)
      write (iunit,intfmt) "rotational number",nint(sym)
      write (iunit,dblfmt) "scaling factor   ",fscal,"    "
      select case (emodelunit)
      case (1)
        write (iunit,chrfmt) "vib.entropy model      ","Grimme (2012)"
        write (iunit,dblfmt) "rotor cutoff     ",sthr,"cm^-1"
      case (2)
        write (iunit,chrfmt) "vib.entropy model      ","Truhlar (2011)"
        write (iunit,dblfmt) "frequency cutoff ",sthr,"cm^-1"
      end select

      write (iunit,dblfmt) "imag. cutoff     ",ithr,"cm^-1"
      write (iunit,'(10x,":",51("."),":")')
    end if

    vibs = vibs*rcmtoau   ! thermodyn needs vibs and zp in Eh

    zp = 0.5_wp*sum(vibs(1:nvib))
    adum = abs(temps-298.15d0)
    rt = minloc(adum,1)  !temperature closest to 298.15 is the ref.
    do j = 1,nt
      if ((j == rt).and.pr) then
        pr2 = .true.
      else
        pr2 = .false.
      end if
      if (pr2) then
        select case (emodelunit)
        case (1)
          call print_thermo_sthr_ts(iunit,nvib,vibs,avmom,sthr,temps(j))
        case (2)
          call print_thermo_sthr_cut(iunit,nvib,vibs,sthr,temps(j))
        end select
      end if
      call thermodyn(iunit,a,b,c,avmom,linear,atom,sym,molmass,vibs,nvib, &
      & temps(j),sthr,et(j),ht(j),gt(j),ts(j),zp,pr2,emodel=emodelunit)
      stot(j) = (ts(j)/temps(j))*autocal
    end do

    if (pr) then
      write (iunit,'(a)')
      write (iunit,'(a10)',advance='no') "T/K"
      write (iunit,'(a16)',advance='no') "H(0)-H(T)+PV"
      write (iunit,'(a16)',advance='no') "H(T)/Eh"
      write (iunit,'(a16)',advance='no') "T*S/Eh"
      write (iunit,'(a16)',advance='no') "G(T)/Eh"
      write (iunit,'(a)')
      write (iunit,'(3x,72("-"))')
      do i = 1,nt
        write (iunit,'(3f10.2)',advance='no') temps(i)
        write (iunit,'(3e16.6)',advance='no') ht(i)
        write (iunit,'(3e16.6)',advance='no') et(i)
        write (iunit,'(3e16.6)',advance='no') ts(i)
        write (iunit,'(3e16.6)',advance='no') gt(i)
        if (i == rt.and.nt > 1) then
          write (iunit,'(1x,"(used)")')
        else
          write (iunit,'(a)')
        end if
      end do
      write (iunit,'(3x,72("-"))')
    end if

    xyz = xyz*aatoau !> NOTE: BACK TO BOHRS

    !$omp critical
    deallocate (vibs)
    !$omp end critical
    return
  end subroutine calcthermo

  subroutine calc_thermo_from_hess(mol,hess,pr,nt,temps,ithr,&
  & fscal,sthr,et,ht,gt,stot,etot)
    type(coord),intent(inout) :: mol
    integer :: nat3
    integer :: io,iunit
    logical :: pr
    real(wp) :: ithr,fscal,sthr
    real(wp),intent(in) :: temps(nt)
    integer,intent(in) :: nt
    real(wp),allocatable,intent(out) :: et(:),ht(:),gt(:),stot(:)
    real(wp),intent(inout) :: hess(:,:)
    real(wp),allocatable :: freq(:)
    real(wp),intent(in) :: etot
    real(wp) :: zpve
    integer :: nrt
    real(wp),allocatable :: int_temps(:)
    character(len=*),parameter :: outfmt = &
    &  '(10x,"::",1x,a,f24.12,1x,a,1x,"::")'

    nat3 = 3*mol%nat
    !$omp critical
    allocate (freq(nat3))
    allocate (et(nt))
    allocate (ht(nt))
    allocate (gt(nt))
    allocate (stot(nt))
    allocate (int_temps(nt))
    !$omp end critical

    int_temps = abs(temps-298.15_wp)
    nrt = minloc(int_temps(:),1)

    call prj_mw_hess(mol%nat,mol%at,nat3,mol%xyz,hess)

    call frequencies(mol%nat,mol%at,mol%xyz,nat3,hess,freq,io)

    call calcthermo(mol%nat,mol%at,mol%xyz,freq,pr,ithr,fscal,sthr,nt,temps, &
        &      et,ht,gt,stot)

    zpve = et(nrt)-ht(nrt)
    if (pr) then
      write (stdout,*)
      write (stdout,'(10x,a)') repeat(':',50)
      write (stdout,'(10x,"::",7x,a,f12.2,1x,a,8x,"::")') "THERMODYNAMICS at",temps(nrt),'K'
      write (stdout,'(10x,a)') repeat(':',50)
      write (stdout,outfmt) 'TOTAL FREE ENERGY',etot+gt(nrt),'Eh'
      write (stdout,'(10x,a)') '::'//repeat('-',46)//'::'
      write (stdout,outfmt) 'total energy     ',etot,'Eh'
      write (stdout,outfmt) 'ZPVE             ',zpve,'Eh'
      write (stdout,outfmt) 'G(RRHO) w/o ZPVE ',gt(nrt)-zpve,'Eh'
      write (stdout,outfmt) 'G(RRHO) total    ',gt(nrt),'Eh'
      write (stdout,'(10x,a)') repeat(':',50)
    end if

  end subroutine calc_thermo_from_hess

  subroutine effective_hessian(nat,nat3,grad1_i,grad2_i,hess1,hess2,heff)
!******************************************************************
!* Effective Hessian at an MECP is computed via Eq. 27 and Eq. 28
!* in https://doi.org/10.1002/qua.25124
!******************************************************************
    implicit none
    integer,intent(in) :: nat,nat3
    integer :: i,j,ii
    real(wp),intent(in) :: grad1_i(3,nat3),grad2_i(3,nat3)
    real(wp) :: grad1(nat3),grad2(nat3),dot

    real(wp),intent(in) :: hess1(nat3,nat3),hess2(nat3,nat3)

    real(wp) :: gnorm1,gnorm2,grad_diff_norm
    real(wp) :: grad_diff(nat3),heff_temp(nat3,nat3)

    real(wp),intent(inout) :: heff(nat3,nat3)
    real(wp),allocatable :: proj_vec(:,:)

    real(wp) :: freq(nat3)

    integer,allocatable :: iwork(:)
    real(wp),allocatable :: work(:)

    integer :: lwork,liwork,info

    !$omp critical
    allocate (proj_vec(nat3,nat3),source=0.0_wp)
    !$omp end critical

    grad1 = reshape(grad1_i, (/nat3/))
    grad2 = reshape(grad2_i, (/nat3/))

    gnorm1 = norm2(grad1)

    gnorm2 = norm2(grad2)

    grad_diff = grad1-grad2

    grad_diff_norm = norm2(grad_diff)

    dot = dot_product(grad1,grad2)

    if (dot .gt. 0.0_wp) then !sloped: dot > 0.0 --> -  | peaked: dot <= 0.0 --> +

      write (stdout,*) 'MECI is considered as a sloped CI'
      write (stdout,*)

      heff = (gnorm1*hess2-gnorm2*hess1)/grad_diff_norm

    else

      write (stdout,*) 'MECI is considered as a peaked CI'
      write (stdout,*)

      heff = (gnorm1*hess2+gnorm2*hess1)/grad_diff_norm

    end if

    !Outer Product of grad_diff

    !Building projection matrix

    !proj_vec = 1 - (dg/|dg| o dg.T/|dg|) = 1 - (dg o dg.T)/|dg|**2

    grad_diff_norm = grad_diff_norm**2

    do i = 1,nat3
      proj_vec(i,:) = -grad_diff(i)*grad_diff/grad_diff_norm
      proj_vec(i,i) = proj_vec(i,i)+1
    end do

    !Projection
    heff = matmul(matmul(proj_vec,heff),proj_vec)

    !Check if hess1 and hess2 are assigned correctly, otherwise change
    lwork = 1+6*nat3+2*nat3**2
    liwork = 3+5*nat3
    !$omp critical
    allocate (work(lwork),iwork(liwork))
    !$omp end critical

    heff_temp = heff

    call dsyevd('V','U',nat3,heff_temp,nat3,freq,work,lwork,iwork,liwork,info)

    deallocate (work,iwork)

    if (0 .gt. sum(freq)) then
      heff = -heff
    end if

  end subroutine effective_hessian

!============================================================================!
!############################################################################!
!============================================================================!
!> PRINTOUT ROUTINES
  subroutine print_vib_spectrum(nat,at,nat3,xyz,freq,dir,fname,ir_int)
!******************************************************************
!* Prints the frequencies in Turbomole "vibspectrum" format.
!* When ir_int is provided, genuine intensities (km/mol) and
!* a threshold-based IR-activity flag are written; otherwise a
!* dummy intensity of 1000 km/mol is used for all active modes.
!* Input:
!*   ir_int (optional) - IR intensities in km/mol, shape (nat3)
!******************************************************************
    integer,intent(in) :: nat,nat3
    integer :: at(nat),i,ich
    real(wp) ::  xyz(3,nat)
    real(wp) ::  freq(nat3),thr
    real(wp),intent(in),optional :: ir_int(nat3)
    character(len=*) :: fname
    character(len=*) :: dir
    !> IR-activity threshold in km/mol
    real(wp),parameter :: ir_thr = 1.0_wp
    real(wp) :: irint_i
    character(len=3) :: ir_flag

    thr = 0.01_wp
    if (len_trim(dir) .eq. 0) then
      open (newunit=ich,file=fname)
    else
      if (directory_exist(dir)) then
        open (newunit=ich,file=dir//'/'//fname)
      else
        open (newunit=ich,file=fname)
      end if
    end if

    write (ich,'("$vibrational spectrum")')
    write (ich,'("#  mode    symmetry    wave number    IR intensity    selection rules")')
    write (ich,'("#                       1/cm              km/mol         IR    RAMAN")')

    do i = 1,nat3
      if (abs(freq(i)) .lt. thr) then
        write (ich,'(i6,9x,    f18.2,f16.5,7x," - ",5x," - ")') &
          i,freq(i),0.0_wp
      else
        if (present(ir_int)) then
          irint_i = ir_int(i)
          if (irint_i >= ir_thr) then
            ir_flag = 'YES'
          else
            ir_flag = 'NO '
          end if
        else
          irint_i = 1000.0_wp
          ir_flag = 'YES'
        end if
        write (ich,'(i6,8x,"a",f18.2,f16.5,7x,a3,5x,"YES")') &
          i,freq(i),irint_i,ir_flag
      end if
    end do

    write (ich,'("$end")')

    close (ich)

  end subroutine print_vib_spectrum

!=========================================================================================!

  subroutine print_g98_fake(nat,at,nat3,xyz,freq,hess,dir,fname,ir_int)
!******************************************************************
!* Prints the vibration spectrum of a system as a g98.out file.
!* Routine is adapted from the xtb code.
!* When ir_int is provided, genuine IR intensities (km/mol) are
!* written; otherwise a dummy value of 99 km/mol is used.
!* Input:
!*   ir_int (optional) - IR intensities in km/mol, shape (nat3)
!******************************************************************
    integer,intent(in) :: nat,nat3
    integer :: at(nat)
    integer  :: gu,i,j,ka,kb,kc,la,lb,k

    real(wp) ::  xyz(3,nat)
    real(wp),intent(in) :: hess(nat3,nat3)
    real(wp),intent(in),optional :: ir_int(nat3)
    real(wp) ::  freq(nat3),red_mass(nat3),force(nat3),ir_f(nat3),zero(1),f2(nat3),u(nat3,nat3)

    character(len=2) :: irrep
    character(len=*) :: fname
    character(len=*) :: dir

    irrep = 'a'

    red_mass = 99.0
    force = 99.0
    ir_f = 99.0
    zero = 0.0

    k = 0

    do i = 1,nat3
      if (abs(freq(i)) .gt. 1.d-1) then
        k = k+1
        u(1:nat3,k) = hess(1:nat3,i)
        f2(k) = freq(i)
        if (present(ir_int)) ir_f(k) = ir_int(i)
      end if
    end do

    if (len_trim(dir) .eq. 0) then
      open (newunit=gu,file=fname)
    else
      if (directory_exist(dir)) then
        open (newunit=gu,file=dir//'/'//fname)
      else
        open (newunit=gu,file=fname)
      end if
    end if

    write (gu,'('' Entering Gaussian System'')')
    write (gu,'('' *********************************************'')')
    write (gu,'('' Gaussian 98:'')')
    write (gu,'('' frequency output generated by the crest code'')')
    write (gu,'('' *********************************************'')')

    write (gu,*) '                        Standard orientation:'
    write (gu,*) '---------------------------------------------', &
        & '-----------------------'
    write (gu,*) ' Center     Atomic     Atomic', &
        & '              Coordinates (Angstroms)'
    write (gu,*) ' Number     Number      Type ', &
        & '             X           Y           Z'
    write (gu,*) '-----------------------', &
        & '---------------------------------------------'
    j = 0
    do i = 1,nat
      write (gu,111) i,at(i),j,xyz(1:3,i)*0.52917726
    end do
    write (gu,*) '----------------------', &
        & '----------------------------------------------'
    write (gu,*) '    1 basis functions        1 primitive gaussians'
    write (gu,*) '    1 alpha electrons        1 beta electrons'
    write (gu,*)
111 format(i5,i11,i14,4x,3f12.6)

    write (gu,*) 'Harmonic frequencies (cm**-1), IR intensities',' (km*mol⁻¹),'
    write (gu,*) 'Raman scattering activities (A**4/amu),', &
        & ' Raman depolarization ratios,'
    write (gu,*) 'reduced masses (AMU), force constants (mDyne/A)', &
        & ' and normal coordinates:'

    ka = 1
    kc = 3

60  kb = min0(kc,k)
    write (gu,100) (j,j=ka,kb)
    write (gu,105) (irrep,j=ka,kb)
    write (gu,110) ' Frequencies --', (f2(j),j=ka,kb)
    write (gu,110) ' Red. masses --', (red_mass(j),j=ka,kb)
    write (gu,110) ' Frc consts  --', (force(j),j=ka,kb)
    write (gu,110) ' IR Inten    --', (ir_f(j),j=ka,kb)
    write (gu,110) ' Raman Activ --', (zero,j=ka,kb)
    write (gu,110) ' Depolar     --', (zero,j=ka,kb)
    write (gu,*) 'Atom AN      X      Y      Z        X      Y', &
        & '      Z        X      Y      Z'
    la = 1
70  lb = nat
    do i = la,lb
      write (gu,130) i,at(i), (u(i*3-2,j),u(i*3-1,j),u(i*3,j),j=ka,kb)
    end do
    if (lb .eq. nat) go to 90
    go to 70
90  if (kb .eq. k) then
      goto 200
    end if

    ka = kc+1
    kc = kc+3
    go to 60

100 format(3(20x,i3))
105 format(3x,3(18x,a5))
110 format(a15,f11.4,12x,f11.4,12x,f11.4)
130 format(2i4,3(2x,3f7.2))
200 continue
    write (gu,'(''end of file'')')
    close (gu)

  end subroutine print_g98_fake

!=========================================================================================!

  subroutine print_hessian(hess,nat3,dir,fname)
!*******************************
!* Prints the numerical hessian
!*******************************
    integer :: nat3,i,j,k,ich
    real(wp) :: hess(nat3,nat3)
    character(len=*) :: fname
    character(len=*) :: dir

    if (len_trim(dir) .eq. 0) then
      open (newunit=ich,file=fname)
      write (stdout,'(1x,a)',advance='no') 'Will be written to file "'//fname//'" ...'
    else
      if (directory_exist(dir)) then
        open (newunit=ich,file=dir//'/'//fname)
        write (stdout,'(1x,a)',advance='no') 'Will be written to file "'//dir//'/'//fname//'" ...'
      else
        open (newunit=ich,file=fname)
        write (stdout,'(1x,a)',advance='no') 'Will be written to file "'//fname//'" ...'
      end if
    end if
    flush (stdout)

    write (ich,'(1x,a)') '$hessian'
    do i = 1,nat3
      k = 0
      do j = 1,nat3
        k = k+1
        if (k .le. 4) then
          write (ich,'(f16.8)',advance='no') hess(i,j)
        else
          write (ich,'(f16.8)') hess(i,j)
          k = 0
        end if
      end do
      if (k .ne. 0) then
        write (ich,*)
      end if
    end do
    write (ich,'(1x,a)') '$end'
    close (ich)

    write (stdout,*) 'done.'
    write (stdout,*)

  end subroutine print_hessian

!============================================================================!
!############################################################################!
!============================================================================!
end module thermochem_module
