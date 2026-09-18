!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021-2025 Christoph Plett, Philipp Pracht
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

module qcg_utils
  use crest_parameters,only:stdout,wp
  use crest_data
  use iomod
  implicit none
  public

contains

!========================================================================================!
!> Level of theory of all energies entering Esolv/Gsolv
!========================================================================================!
  function qcg_energy_level(env) result(lvl)
    !***********************************************
    !* Return the method flag of the QCG energy level,
    !* i.e. the level of the final ensemble and CFF
    !* single points. The monomer reference energies
    !* must be evaluated at the same level.
    !***********************************************
    type(systemdata),intent(in) :: env
    character(len=20) :: lvl
    if (env%final_gfn2_opt) then
      lvl = '--gfn2'
    else
      lvl = env%ensemble_opt
    end if
  end function qcg_energy_level

!========================================================================================!
!> Convert given QCG coordinate files into (TM format)
!> Write "solute" and "solvent" coordinate files
!========================================================================================!
  subroutine inputcoords_qcg(env,solute,solvent)
    use crest_parameters
    use crest_data
    use strucrd
    use qcg_coord_type
    use iomod
    implicit none

    type(systemdata),intent(inout) :: env
    type(coord_qcg),intent(out) :: solute,solvent

    logical :: ex11,ex21,solu,solv
    type(coord_qcg) :: mol
    integer :: i

!--------------------Checking for input-------------!

    !Solute
    inquire (file=env%solu_file,exist=ex11)
    inquire (file='solute',exist=solu)
    if (solu) call copy('solute','solute.old') !Backup solute file
    if ((.not.ex11).and.(.not.solu)) then
      error stop 'No (valid) solute file! exit.'
    else if ((.not.ex11).and.(solu)) then
      env%solu_file = 'solute'
    end if

    !Solvent
    inquire (file=env%solv_file,exist=ex21)
    inquire (file='solvent',exist=solv)
    if (solu) call copy('solvent','solvent.old') !Backup solvent file
    if ((.not.ex21).and.(.not.solv)) then
      error stop 'No (valid) solvent file! exit.'
    else if ((.not.ex11).and.(solu)) then
      env%solu_file = 'solvent'
    end if

!---------------Handling solute---------------------!
    call mol%open(env%solu_file)
    call mol%write('solute')
    solute%nat = mol%nat
    solute%at = mol%at
    solute%xyz = mol%xyz

    !--- if the input was a SDF file, special handling
    env%sdfformat = .false.
    call checkcoordtype(env%solu_file,i)
    if (i == 31.or.i == 32) then
      !Add sdf stuff here, if somebody needs it
    end if
    !--- Add as ref structure in env
    call env%ref%load(mol)
    call mol%deallocate()
!---------------Handling solvent---------------------!

    call mol%open(env%solv_file)
    call mol%write('solvent')
    solvent%nat = mol%nat
    solvent%at = mol%at
    solvent%xyz = mol%xyz
    call mol%deallocate()

    !--- if the input was a SDF file, special handling
    env%sdfformat = .false.
    call checkcoordtype(env%solv_file,i)
    if (i == 31.or.i == 32) then
      !Add sdf stuff here, if somebody needs it
    end if

    return
  end subroutine inputcoords_qcg

!==============================================================================!

  subroutine write_reference(env,solu,clus)
    use crest_data
    use qcg_coord_type
    use iomod
    use strucrd
    implicit none
    type(systemdata):: env    ! MAIN STORAGE OS SYSTEM DATA
    type(coord_qcg)            :: solu,clus
    type(coord_qcg)            :: ref_mol,ref_clus
    ref_mol = solu
    call rdcoord(env%solu_file,ref_mol%nat,ref_mol%at,ref_mol%xyz) !original solute coordinates
    call remove(env%fixfile)
    ref_clus = clus
    ref_clus%xyz(1:3,1:solu%nat) = solu%xyz
    call wrc0(env%fixfile,ref_clus%nat,ref_clus%at,ref_clus%xyz)
  end subroutine write_reference

!=============================================================================!

  subroutine aver(pr,env,runs,e_tot,S,H,G,sasa,a_present,a_tot)
    use crest_parameters
    use crest_data

    implicit none
!---- Dummy
    type(systemdata),intent(in)   :: env
    integer,intent(in)            :: runs
    real(wp),intent(inout)        :: e_tot
    real(wp),intent(in),optional  :: a_tot
    real(wp),intent(out)          :: S
    real(wp),intent(out)          :: H
    real(wp),intent(out)          :: G
    real(wp),intent(out)          :: sasa
!---- Stack
    logical,intent(in)            :: pr,a_present
    integer                       :: j,jmin
    real(wp)                      :: A
    real(wp)                      :: e0
    real(wp),allocatable          :: de(:)
    real(wp),allocatable          :: p(:)
    real(wp)                      :: pmax
    real(wp)                      :: eav
    real(wp)                      :: area
    real(wp)                      :: beta
    real(wp)                      :: temp
    integer                       :: ich48
    dimension e_tot(runs)
    dimension a_tot(runs)

    temp = env%tboltz
    allocate (de(runs),source=0.0d0)
    allocate (p(runs),source=0.0d0)

    beta = 1./(temp*8.314510/4.184/1000.+1.d-14)
    e0 = e_tot(1)
    de(1:runs) = (e_tot(1:runs)-e0)
    call qcg_boltz(env,runs,de,p)

    A = 0
    eav = 0
    pmax = 0
    area = 0
    do j = 1,runs
      A = A+p(j)*log(p(j)+1.d-12)
      eav = eav+p(j)*e_tot(j)
      if (p(j) .gt. pmax) then
        pmax = p(j)
        jmin = j
      end if
      if (a_present) area = area+p(j)*a_tot(j)
    end do
    sasa = area
    S = (1./beta)*A
    H = eav
    G = eav+S
    if (pr) then
      open (newunit=ich48,file='population.dat')
      write (ich48,'(2x, ''cluster'',2x,''E_norm [Eh]'',2x, ''De [kcal]'', 4x, ''p'')')
      do j = 1,runs
        if (j .lt. 10) then
          write (ich48,'(5x,i0,3x,f11.6,5x,f6.4,3x,f6.4)') j,e_tot(j)/autokcal,de(j),p(j)
        else
          write (ich48,'(5x,i0,2x,f11.6,5x,f6.4,3x,f6.4)') j,e_tot(j)/autokcal,de(j),p(j)
        end if
      end do
      write (ich48,*)
      write (ich48,'(''Ensemble free energy [Eh]:'', f20.10)') G/autokcal
      close (ich48)
    end if

    deallocate (de,p)

  end subroutine aver

  !==============================================================================!

  subroutine get_sphere(pr,mol,r_logical)
    use crest_parameters
    use qcg_coord_type
    use miscdata
    implicit none
    type(coord_qcg),intent(inout) :: mol
    type(coord_qcg) :: dum
    logical        :: pr
    logical        :: r_logical !Determines wether r is overwritten or not
    real(wp),parameter :: pi43 = pi*4.0d0/3.0d0
    real(wp),parameter :: third = 1.0d0/3.0d0
    integer :: i
    real(wp) :: rad(mol%nat),xyz_tmp(3,mol%nat)
    external get_volume

    do i = 1,mol%nat
      rad(i) = bohr*rcov_qcg(mol%at(i))*1.40 ! scale factor adjusted to rough
      xyz_tmp(1:3,i) = bohr*mol%xyz(1:3,i)
    end do

    dum = mol
    dum%xyz = xyz_tmp

    call get_volume(dum,rad)

    mol%atot = dum%atot/bohr**2
    mol%vtot = dum%vtot/bohr**3

    if (r_logical) then
      mol%rtot = mol%vtot*3.0/4.d0/pi
      mol%rtot = mol%rtot**(1.d0/3.d0)
    end if

    if (pr) then
      if (r_logical) then
        write (stdout,'(2x,''molecular radius (Bohr**1):'',F8.2)') mol%rtot
      end if
      write (stdout,'(2x,''molecular area   (Bohr**2):'',F8.2)') mol%atot
      write (stdout,'(2x,''molecular volume (Bohr**3):'',F8.2)') mol%vtot
    end if
  end subroutine get_sphere

  !=============================================================================!
  !
  subroutine cma_shifting(solu,solv)
    use crest_parameters
    use crest_data
    use iomod
    use qcg_coord_type
    use strucrd
    use axis_module,only:cma
    implicit none

    type(coord_qcg)    :: solu,solv

    integer            :: i

    call cma(solu%nat,solu%at,solu%xyz,solu%cma)
    call cma(solv%nat,solv%at,solv%xyz,solv%cma)

    do i = 1,solu%nat
      solu%xyz(1:3,i) = solu%xyz(1:3,i)-solu%cma(1:3)
    end do
    do i = 1,solv%nat
      solv%xyz(1:3,i) = solv%xyz(1:3,i)-solv%cma(1:3)
    end do

  end subroutine cma_shifting

!==============================================================================!
!
  subroutine get_ellipsoid(env,solu,solv,clus,pr1)
    use crest_parameters
    use crest_data
    use iomod
    use qcg_coord_type
    use strucrd
    use axis_module
    implicit none

    type(systemdata)   :: env
    type(coord_qcg)    :: solu,solv,clus
    type(coord_qcg)    :: dummy_solu,dummy_solv
    real(wp)           :: rabc_solu(3),rabc_solv(3)
    real(wp)           :: aniso,sola
    real(wp)           :: rmax_solu,rmax_solv
    real(wp)           :: boxr,roff,r
    character(len=10) :: fname
    logical            :: ex,pr,pr1

    real(wp),parameter :: pi43 = pi*4.0d0/3.0d0
    real(wp),parameter :: third = 1.0d0/3.0d0

    pr = .false. !Outprint deactivated

    fname = 'eaxis.qcg'
    inquire (file=fname,exist=ex)

    if (pr1) then !First time called
!--- Moving all coords to the origin (transformation)
      call axistrf(solu%nat,solu%nat,solu%at,solu%xyz)
!  call axistrf(solv%nat,solv%nat,solv%at,solv%xyz)  !Not done in original QCG code
      call axistrf(clus%nat,solu%nat,clus%at,clus%xyz)

!--- Overwrite solute and solvent coord in original file with transformed and optimized ones
      call solu%write('solute')
      call solv%write('solvent')

!--- Getting axis
      write (stdout,*) 'Solute:'
      call axis(pr1,solu%nat,solu%at,solu%xyz,solu%eax)
      write (stdout,*) 'Solvent:'
      call axis(pr1,solv%nat,solv%at,solv%xyz,solv%eax)
      write (stdout,*)
    end if

!--- Computing anisotropy factor of solute and solvent
    sola = sqrt(1.+(solu%eax(1)-solu%eax(3))/((solu%eax(1)+solu%eax(2)+solu%eax(3))/3.))
    aniso = sqrt(1.+(solv%eax(1)-solv%eax(3))/((solv%eax(1)+solv%eax(2)+solv%eax(3))/3.)) ! =1 for a spherical system

!--- Get maximum intramoleclar distance of solute and solvent
    call getmaxrad(solu%nat,solu%at,solu%xyz,rmax_solu)
    call getmaxrad(solv%nat,solv%at,solv%xyz,rmax_solv)

!--- Getting V and A of dummies
    dummy_solu = solu
    dummy_solv = solv !Why is dummy_solv%vtot different to solv%vtot
    call get_sphere(.false.,dummy_solu,.false.)
    call get_sphere(.false.,dummy_solv,.false.)

!--- Computation of outer Wall
    roff = sola*dummy_solu%vtot/1000
    boxr = ((0.5*aniso*clus%nmol*dummy_solv%vtot+dummy_solu%vtot)/pi43)**third+roff+rmax_solv*0.5 !0.5 both
    r = (boxr**3/(solu%eax(1)*solu%eax(2)*solu%eax(3)))**third       ! volume of ellipsoid = volume of sphere
    rabc_solv = solu%eax*r                              ! outer solvent wall

!--- Computation of inner wall
    roff = sola*dummy_solu%vtot/1000
    boxr = ((sola*dummy_solu%vtot)/pi43)**third+roff+rmax_solu*0.1 !0.1 before
    r = (boxr**3/(solu%eax(1)*solu%eax(2)*solu%eax(3)))**third       ! volume of ellipsoid = volume of sphere
    rabc_solu = solu%eax*r
    dummy_solu%ell_abc(1) = solu%eax(1)**2/sum((solu%eax(1:3))**2)
    dummy_solu%ell_abc(2) = solu%eax(2)**2/sum((solu%eax(1:3))**2)
    dummy_solu%ell_abc(3) = solu%eax(3)**2/sum((solu%eax(1:3))**2)
    rabc_solu = dummy_solu%ell_abc*r

    solu%aniso = sola
    solv%aniso = aniso
    solu%ell_abc = rabc_solu
    clus%ell_abc = rabc_solv*env%potscal

    if (pr1) then
      write (stdout,'(2x,''solvent anisotropy            :'',4f10.3)') aniso
      write (stdout,'(2x,''solute anisotropy             :'',4f10.3)') sola
      write (stdout,'(2x,''roff inner wall               :'',4f10.3)') roff
      write (stdout,'(2x,''solute max dist               :'',4f10.3)') rmax_solu
      write (stdout,'(2x,''solvent max dist              :'',4f10.3)') rmax_solv
      write (stdout,'(2x,''inner unit axis               :'',3f10.3)') dummy_solu%ell_abc(1:3)
      write (stdout,'(2x,''inner ellipsoid/Bohr          :'',3f10.3)') rabc_solu(1:3)
      write (stdout,'(2x,''scaling factor outer ellipsoid:'',3f10.3)') env%potscal
      write (stdout,'(2x,''outer ellipsoid/Bohr          :'',3f10.3)') clus%ell_abc(1:3)
      if (env%potscal .gt. 1.0_wp) write &
           &(stdout,'(2x,''!!!WARNING: A SCALING FACTOR LARGER 1.0 IS ONLY RECOMMENDED FOR MICROSOLVATION'')')
      write (stdout,*)
    end if

  end subroutine get_ellipsoid

!==============================================================================!

  subroutine getmaxrad(n,at,xyz,r)
    use crest_parameters,only:wp
    use miscdata,only:rcov_qcg
    implicit none
    real(wp) :: xyz(3,n),r
    integer :: n,at(n)
    real(wp) :: rx,ry,rz,rr
    integer :: i,j

    r = 0
    do i = 1,n-1
      do j = i+1,n
        rx = xyz(1,i)-xyz(1,j)
        ry = xyz(2,i)-xyz(2,j)
        rz = xyz(3,i)-xyz(3,j)
        rr = sqrt(rx**2+ry**2+rz**2)+rcov_qcg(at(i))+rcov_qcg(at(j))
        if (rr .gt. r) r = rr
      end do
    end do
  end subroutine getmaxrad

!==============================================================================!

  subroutine ellipsout(fname,n,at,xyz,r1)
    use crest_parameters
    use strucrd,only:i2e
    implicit none

    integer            :: i
    integer            :: n,at(n)
    real(wp)           :: xyz(3,n),r1(3)
    real(wp)           :: x,y,z,f,rr
    character(len=*)   :: fname
    integer            :: ich11

    open (newunit=ich11,file=fname)
    write (ich11,'(a)') '$coord'
    do i = 1,n
      write (ich11,'(3F24.14,6x,a)') xyz(1,i),xyz(2,i),xyz(3,i),i2e(at(i))
    end do
    do i = 1,500
      call random_number(x)
      call random_number(f)
      if (f .gt. 0.5) x = -x
      call random_number(y)
      call random_number(f)
      if (f .gt. 0.5) y = -y
      call random_number(z)
      call random_number(f)
      if (f .gt. 0.5) z = -z
      rr = sqrt(x*x+y*y+z*z)
      x = x*r1(1)/rr
      y = y*r1(2)/rr
      z = z*r1(3)/rr
      write (ich11,'(3F24.14,6x,a2)') x,y,z,'he'
    end do
    write (ich11,'(a)') '$end'
    close (ich11)

  end subroutine ellipsout

!==============================================================================!

  subroutine both_ellipsout(fname,n,at,xyz,r1,r2)
    use crest_parameters
    use strucrd,only:i2e
    implicit none

    integer            :: i
    integer            :: n,at(n)
    real(wp)           :: xyz(3,n),r1(3)
    real(wp),optional :: r2(3)
    real(wp)           :: x,y,z,f,rr
    character(len=*)   :: fname
    integer            :: ich11

    open (newunit=ich11,file=fname)
    write (ich11,'(a)') '$coord'
    do i = 1,n
      write (ich11,'(3F24.14,6x,a)') xyz(1,i),xyz(2,i),xyz(3,i),i2e(at(i))
    end do
    do i = 1,500
      call random_number(x)
      call random_number(f)
      if (f .gt. 0.5) x = -x
      call random_number(y)
      call random_number(f)
      if (f .gt. 0.5) y = -y
      call random_number(z)
      call random_number(f)
      if (f .gt. 0.5) z = -z
      rr = sqrt(x*x+y*y+z*z)
      x = x*r1(1)/rr
      y = y*r1(2)/rr
      z = z*r1(3)/rr
      write (ich11,'(3F24.14,6x,a2)') x,y,z,'he'
    end do
    if (present(r2)) then
      do i = 1,100
        call random_number(x)
        call random_number(f)
        if (f .gt. 0.5) x = -x
        call random_number(y)
        call random_number(f)
        if (f .gt. 0.5) y = -y
        call random_number(z)
        call random_number(f)
        if (f .gt. 0.5) z = -z
        rr = sqrt(x*x+y*y+z*z)
        x = x*r2(1)/rr
        y = y*r2(2)/rr
        z = z*r2(3)/rr
        write (ich11,'(3F24.14,6x,a2)') x,y,z,'b'
      end do
    end if
    write (ich11,'(a)') '$end'
    close (ich11)

  end subroutine both_ellipsout

!==============================================================================!

  subroutine analyze_cluster(nsolv,n,nS,nM,xyz,at,av,last)
    use crest_parameters
    use axis_module,only:cma
    implicit none
    real(wp) xyz(3,n)
    real(wp) av,last
    integer n,nS,nM,nsolv,at(n)
    real(wp) xyzM(3,nM)
    integer atm(nM)
    real(wp) xyzS(3,nS)
    integer atS(nS)
    real(wp) x1(3),x2(3),r
    integer i,is,ie

    if (nsolv .eq. 1) return
    xyzS(1:3,1:nS) = xyz(1:3,1:nS)
    atS(1:nS) = at(1:nS)
    call cma(nS,atS,xyzS,x1)

    av = 0
    do i = 1,nsolv
      is = nS+(i-1)*nM+1
      ie = is+nM-1
      xyzM(1:3,1:nM) = xyz(1:3,is:ie)
      atM(1:nM) = at(is:ie)
      call cma(nM,atM,xyzM,x2)
      r = sqrt((x1(1)-x2(1))**2+(x1(2)-x2(2))**2+(x1(3)-x2(3))**2)
      if (i .lt. nsolv) then
        av = av+r
      else
        last = r
      end if
    end do
    av = av/real(nsolv-1,wp)
  end subroutine analyze_cluster

!==============================================================================!
!
  subroutine qcg_boltz(env,n,e,p)
    use crest_parameters
    use crest_data
    implicit none
    type(systemdata),intent(in)    :: env
    integer,intent(in)              :: n
    real(wp),intent(in)             :: e(:)
    real(wp),intent(out)            :: p(:)
    integer                         :: i
    real(wp)                        :: temp
    real(wp)                        :: f,hsum,esum

    temp = env%tboltz
    f = 8.314*temp/4.184d+3
    esum = 0
    do i = 1,n
      esum = esum+exp(-e(i)/f)
    end do
    hsum = 0
    do i = 1,n
      p(i) = exp(-e(i)/f)/esum
    end do
  end subroutine qcg_boltz

!==============================================================================!

!==============================================================================!

  subroutine fill_take(env,n2,n12,rabc,ipos)
    use crest_parameters
    use crest_data
    use strucrd
    use axis_module,only:cma
    implicit none

    type(systemdata)      :: env
    integer,intent(in)   :: n2,n12
    real(wp),intent(in)   :: rabc(3)
    integer,intent(out)  :: ipos
    integer               :: i,m,n21
    integer               :: at2(n2),at12(n12)
    integer               :: counter
    real(wp)              :: xyz2(3,n2),xyz12(3,n12)
    real(wp)              :: etmp(100)
    real(wp)              :: eabc
    real(wp)              :: cma2(3)
    real(wp),allocatable  :: dist(:)

    eabc = 0
    counter = 0
    n21 = n12-n2+1
    if (env%use_xtbiff) then
      call rdxtbiffE('xtbscreen.xyz',m,n12,etmp)
    else
      call rdxtbiffE('best.xyz',m,n12,etmp)
    end if

    allocate (dist(m),source=0.0d0)
    dist = 0.0d0

    do i = 1,m
      if (env%use_xtbiff) then
        call rdxmolselec('xtbscreen.xyz',i,n12,at12,xyz12)
      else
        call rdxmolselec('final_structures.xyz',i,n12,at12,xyz12)
      end if

      at2(1:n2) = at12(n21:n12)
      xyz2(1:3,1:n2) = xyz12(1:3,n21:n12)
      call cma(n2,at2,xyz2,cma2)
      call calc_dist(cma2,rabc,dist(i),eabc)
      if (eabc .gt. 1.0d0) then
        dist(i) = 1.0d42
        counter = counter+1
      end if
    end do

    ipos = minloc(dist(1:m),dim=1)

    if (counter .eq. m) ipos = 0

    deallocate (dist)
  end subroutine fill_take

!==============================================================================!

  subroutine calc_dist(xyz,rabc,dist,eabc)
    use crest_parameters
    implicit none

    real(wp),intent(in)    :: xyz(3)
    real(wp),intent(in)    :: rabc(3)
    real(wp),intent(out)   :: dist
    real(wp),intent(out)   :: eabc
    real(wp)                :: center(3),rc(3)

    center = 0.d0
    rc = (xyz(1:3)-center)
    dist = norm2(rc)
    eabc = sum((xyz(1:3)**2)/(rabc(1:3)**2))
  end subroutine calc_dist

!==============================================================================!

  subroutine sort_min(i,j,col,A)
    use crest_parameters
    implicit none
    integer,intent(in)   :: i,j,col
    real*8,intent(inout) :: A(i,j)
    real*8                :: buf(j)
    integer               :: nsize,irow,krow
    nsize = i

    do irow = 1,nsize
      krow = minloc(A(irow:nsize,col),dim=1)+irow-1
      buf(:) = A(irow,:)
      A(irow,:) = A(krow,:)
      A(krow,:) = buf(:)
    end do
  end subroutine sort_min

!==============================================================================!

  subroutine qcg_dump_sorted_ensemble(ens,e_ens,fname)
    use crest_parameters
    use crest_data
    use strucrd
    implicit none
    type(ensemble),intent(in) :: ens
    real(wp),intent(in)       :: e_ens(ens%nall)
    real(wp),allocatable :: dum(:)
    character(len=*) :: fname
    integer          :: ich,i,e_min
    allocate(dum(ens%nall))
    dum(:) = e_ens(:)
    open (newunit=ich,file=fname)
    do i = 1,ens%nall
      e_min = minloc(dum,dim=1)
      call wrxyz(ich,ens%nat,ens%at,ens%xyz(:,:,e_min),e_ens(e_min))
      dum(e_min) = huge(1.0_wp)
    end do
    close (ich)
    deallocate(dum)
  end subroutine qcg_dump_sorted_ensemble

!==============================================================================!

  subroutine rdtherm(fname,ht,svib,srot,stra,gt)
    use crest_parameters
    use crest_data
    use iomod

    implicit none
! Dummy
    real(wp),intent(out)  :: ht
    real(wp),intent(out)  :: gt
    real(wp),intent(out)  :: svib
    real(wp),intent(out)  :: srot
    real(wp),intent(out)  :: stra
! Stack
    integer :: nn,io,counter,hg_line,ich
    real(wp) :: xx(20)
    logical :: ende
    character(len=*) :: fname
    character(len=128) :: a

    ende = .false.
    counter = 0
    hg_line = 0

    open (newunit=ich,file=fname)
    do while (.not.ende)
      read (ich,'(a)',iostat=io) a
      if (io .lt. 0) then
        ende = .true.
        cycle
      end if
      if (index(a,'G(T)/Eh ') .ne. 0) then
        hg_line = counter
      end if
      if (index(a,'  VIB  ') .ne. 0) then
        call readl(a,xx,nn)
        svib = xx(5)
        if (svib .eq. 0.0d0) then
          call readl(a,xx,nn)
          svib = xx(4)
        end if
      end if
      if (index(a,'  ROT  ') .ne. 0) then
        call readl(a,xx,nn)
        srot = xx(4)
      end if
      if (index(a,'  TR   ') .ne. 0) then
        call readl(a,xx,nn)
        stra = xx(4)
      end if
      if (counter .eq. hg_line+2) then
        call readl(a,xx,nn)
        ht = xx(3)*autokcal
        gt = xx(5)*autokcal
      end if
      counter = counter+1
    end do
    close (ich)
  end subroutine rdtherm

!============================================================!
! Read the Energies from a xtbiff output
!============================================================!

  subroutine rdxtbiffE(fname,m,n,e)
    use crest_parameters
    implicit none
    integer :: m,n
    character(len=*),intent(in) :: fname
    real(wp) :: e(:)

    character(len=128) :: line
    real(wp) :: xx(10)
    integer :: ich,i,j,nn

    open (newunit=ich,file=fname)

    j = 1
10  continue
    read (ich,'(a)',end=999) line
    read (ich,'(a)') line
    call readl(line,xx,nn)
    e(j) = xx(1)
    do i = 1,n
      read (ich,'(a)') line
    end do
    j = j+1
    goto 10

999 close (ich)
    m = j-1
  end subroutine rdxtbiffE

!==============================================================================!
!==============================================================================!
end module qcg_utils
