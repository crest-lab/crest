!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2022 Philipp Pracht
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

!====================================================!
! a small module for constraining potentials
!====================================================!

module constraints
  use crest_parameters
  use strucrd
  use wall_setup
  use iso_fortran_env,only:wp => real64
  implicit none

  !=========================================================================================!
  !--- private module variables and parameters
  private
  integer :: i,j,k,l,ich,och,io
  logical :: ex

  !--- some constants and name mappings
  real(wp),parameter :: deg = 180.0_wp/pi ! 1 rad in degrees
  real(wp),parameter :: fcdefault = 0.01_wp
  real(wp),parameter :: Tdefault = 298.15_wp
  real(wp),parameter :: betadefault = 50.0_wp

  !>--- constrain types
  integer,parameter :: bond = 1
  integer,parameter :: allbonds = 11
  integer,parameter :: angle = 2
  integer,parameter :: dihedral = 3
  integer,parameter :: wall = 4
  integer,parameter :: wall_fermi = 5
  integer,parameter :: box = 6
  integer,parameter :: box_fermi = 7
  integer,parameter :: bondrange = 8

  integer,public,parameter :: na_gapdiff = -1
  integer,public,parameter :: na_gapdiff2 = -2

  integer,parameter :: pharmonic = 1
  integer,parameter :: plogfermi = 2

  !=====================================================!

  public :: constraint
  type :: constraint
    !> quick select
    logical :: active = .true.

    !> required
    integer :: type = 0  !> type of the constraint
    integer :: n = 0  !> number of atoms affected
    integer,allocatable :: atms(:)
    real(wp),allocatable :: ref(:)
    real(wp),allocatable :: fc(:)

    !> special directives
    logical :: frozenatms = .false.
    logical,pointer :: freezeptr(:)

    !> other
    real(wp) :: wscal = 1.0_wp
    integer :: subtype = pharmonic  !> currently unused switch

  contains
    procedure :: print => print_constraint
    procedure :: deallocate => constraint_deallocate
    procedure :: bondconstraint => create_bond_constraint
    generic,public :: sphereconstraint => create_sphere_constraint,create_sphere_constraint_all
    procedure,private :: create_sphere_constraint,create_sphere_constraint_all
    procedure :: sphereupdate => sphere_update_nat
    procedure :: ellipsoid => create_ellips_constraint
    procedure :: angleconstraint => create_angle_constraint
    procedure :: dihedralconstraint => create_dihedral_constraint
    procedure :: gapdiffconstraint => create_gapdiff_constraint
    procedure :: gapdiffconstraint2 => create_gapdiff_constraint2
    procedure :: dummyconstraint => create_dummy_constraint
    procedure :: analyzedummy => analyze_dummy_bond_constraint
    procedure :: rdbondconstraint => analyze_dummy_bond_constraint2
    procedure :: bondrangeconstraint => create_bondrange_constraint
    procedure :: addfreeze => constraint_freezeassoc
    procedure :: complete => complete_defaults
  end type constraint

  !=====================================================!

  public :: scantype
  type :: scantype

    integer :: type = 0  !> 0=nothing, 1=distance, 3=dihedral
    integer :: n = 0
    integer :: steps = 5 !> number of sampling steps for the scan
    real(wp) :: minval = 0.0_wp
    real(wp) :: maxval = 0.0_wp
    integer,allocatable :: atms(:)
    real(wp),allocatable :: points(:)
    integer :: constrnmbr = 0
    integer :: restore = 1
    integer :: currentstep = 0
  contains
    procedure :: deallocate => scantype_deallocate
  end type scantype
  !=====================================================!

  public :: calc_constraint

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine complete_defaults(self,mol)
!****************************************************
!* In case there are some defaults missing, complete
!* with sensible settings
!****************************************************
    class(constraint) :: self
    type(coord) :: mol
    integer :: i,j
    real(wp) :: dummy
    select case (self%type)
!>--- bond constraint
    case (bond)
      if (self%n .ne. 2) error stop '*** ERROR *** wrong number of atoms for bond constraint'
      if (.not.allocated(self%fc)) then
        allocate (self%fc(1))
        self%fc(1) = fcdefault
      end if
      if (.not.allocated(self%ref)) then
        allocate (self%ref(1))
        self%ref(1) = mol%dist(self%atms(1),self%atms(2))
      end if

!>--- bondrange constraint
    case (bondrange)
      if (self%n .ne. 2) error stop '*** ERROR *** wrong number of atoms for bondrange constraint'
      if (.not.allocated(self%fc)) then
        allocate (self%fc(2))
        self%fc(1) = fcdefault/kB !> bondrange doesn't use 300K default! 
        self%fc(2) = betadefault
      else
        if (size(self%fc) < 2) error stop '*** ERROR *** wrong number of parameters for bondrange constraint'
      end if
      if (.not.allocated(self%ref)) then
        allocate (self%ref(2))
        self%ref(1) = mol%dist(self%atms(1),self%atms(2))+0.5_wp
        self%ref(2) = self%ref(1)-1.0_wp
      else
        dummy = minval(self%ref(:))
        self%ref(1) =  maxval(self%ref(:))
        self%ref(2) = dummy
      end if

!>--- angle constraint
    case (angle)
      if (self%n .ne. 3) error stop '*** ERROR *** wrong number of atoms for angle constraint'
      if (.not.allocated(self%fc)) then
        allocate (self%fc(1))
        self%fc(1) = fcdefault
      end if
      if (.not.allocated(self%ref)) then
        allocate (self%ref(1))
        !> reference in rad
        self%ref(1) = mol%angle(self%atms(1),self%atms(2),self%atms(3))
      end if

!>--- dehedral constraint
    case (dihedral)
      if (self%n .ne. 4) error stop '*** ERROR *** wrong number of atoms for dihedral constraint'
      if (.not.allocated(self%fc)) then
        allocate (self%fc(1))
        self%fc(1) = fcdefault
      end if
      if (.not.allocated(self%ref)) then
        allocate (self%ref(1))
        !> reference in rad
        self%ref(1) = mol%dihedral(self%atms(1),self%atms(2),self%atms(3),self%atms(4))
      end if

!>--- wall potential (exponential)
    case (wall)
      if (self%n .eq. 0.or..not.allocated(self%atms)) then
        self%n = mol%nat
        if (allocated(self%atms)) deallocate (self%atms)
        allocate (self%atms(mol%nat))
        do i = 1,mol%nat
          self%atms(i) = i
        end do
      end if
      if (.not.allocated(self%fc)) then
        allocate (self%fc(2))
        self%fc(1) = fcdefault
        self%fc(2) = 6.0_wp
      else
        if (size(self%fc) < 2) error stop '*** ERROR *** wrong number of parameters for wall potential'
      end if
      if (.not.allocated(self%ref)) then
        allocate (self%ref(3))
        call boxpot_core(mol,self%ref,self%wscal)
      else
        if (size(self%ref) < 3) error stop '*** ERROR *** wrong number of reference vaules for wall    potential'
      end if

!>--- wall potential (logfermi)
    case (wall_fermi)
      if (self%n .eq. 0.or..not.allocated(self%atms)) then
        self%n = mol%nat
        if (allocated(self%atms)) deallocate (self%atms)
        allocate (self%atms(mol%nat))
        do i = 1,mol%nat
          self%atms(i) = i
        end do
      end if

      if (.not.allocated(self%fc)) then
        allocate (self%fc(2))
        self%fc(1) = Tdefault
        self%fc(2) = betadefault
      else
        if (size(self%fc) < 2) error stop '*** ERROR *** wrong number of parameters for wall potential'
      end if
      if (.not.allocated(self%ref)) then
        allocate (self%ref(3))
        call boxpot_core(mol,self%ref,self%wscal)
      else
        if (size(self%ref) < 3) error stop '*** ERROR *** wrong number of reference vaules for wall potential'
      end if

!>--- gap difference potential (only needs fc)
    case (na_gapdiff)
      if (.not.allocated(self%fc)) then
        allocate (self%fc(2))
        self%fc(1) = 10.0_wp  !> sigma
        self%fc(2) = 0.005_wp !> alpha
      else
        if (size(self%fc) < 2) error stop '*** ERROR *** wrong number of parameters for gap potential'
      end if

!>--- new gap difference potential (only needs fc)
    case (na_gapdiff2)
      if (.not.allocated(self%fc)) then
        allocate (self%fc(3))
        self%fc(1) = 10.0_wp  !> sigma
        self%fc(2) = 0.005_wp !> alpha
        self%fc(3) = 0.25_wp  !> k
      else
        if (size(self%fc) < 3) error stop '*** ERROR *** wrong number of parameters for gap potential'
      end if

    end select
  end subroutine complete_defaults

!========================================================================================!
!========================================================================================!

  subroutine calc_constraint(n,xyz,constr,energy,grd)
!*****************************************************
!* Calculate energy and gradient contribution for
!* the given constraint.
!* In and outputs are in atomic units (Bohr, Hartree)
!*****************************************************
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    type(constraint) :: constr

    real(wp),intent(out) :: energy
    real(wp),intent(out) :: grd(3,n)

    energy = 0.0_wp
    grd = 0.0_wp

    if(.not.constr%active) return

    select case (constr%type)
    case (bond)
      call bond_constraint(n,xyz,constr,energy,grd)
    case (angle)
      call angle_constraint(n,xyz,constr,energy,grd)
    case (dihedral)
      call dihedral_constraint(n,xyz,constr,energy,grd)
    case (wall,wall_fermi)
      call wall_constraint(n,xyz,constr,energy,grd,constr%type)
    case (box)

    case (box_fermi)

    case (bondrange)
      call bondrange_constraint(n,xyz,constr,energy,grd)
    case default
      return
    end select

    return
  end subroutine calc_constraint

!========================================================================================!
!========================================================================================!

  subroutine print_constraint(self,chnl)
    implicit none
    class(constraint) :: self
    character(len=64) :: art
    character(len=64) :: atoms
    character(len=258) :: values
    character(len=10) :: atm
    integer :: chnl
    logical :: pr
    character(len=*),parameter :: headfmt ='("> constraint: ",a,a)'
    if (self%type == 0) return
    pr = .true.
    select case (self%type)
    case (bond)
      art = 'distance'
      write (atoms,'(1x,"atoms:",1x,i0,",",i0)') self%atms(1:2)
      write (chnl,headfmt) trim(art),trim(atoms)
      write (chnl,'(2x,"d(AA)=",f12.5,1x,"k=",f8.5)') self%ref(1)*autoaa,self%fc(1)
    case (angle)
      art = 'angle'
      write (atoms,'(1x,"atoms:",1x,i0,",",i0,",",i0)') self%atms(1:3)
      write (chnl,headfmt) trim(art),trim(atoms)
      write (chnl,'(2x,"deg=",f6.2,1x,"k=",f8.5)') self%ref(1)*deg,self%fc(1)
    case (dihedral)
      art = 'dihedral'
      write (atoms,'(1x,"atoms:",1x,i0,",",i0,",",i0,",",i0)') self%atms(1:4)
      write (chnl,headfmt) trim(art),trim(atoms)
      write (chnl,'(2x,"deg=",f6.2,1x,"k=",f8.5)') self%ref(1)*deg,self%fc(1)
    case (wall)
      art = 'wall'
      write (atoms,'(1x,"atoms:",1x,i0,a)') self%n,'/all'
      write (chnl,headfmt) trim(art),trim(atoms)
      write (chnl,'(2x,"radii(AA)=",3f12.5)') self%ref(1:3)*autoaa
      write (chnl,'(2x,"kb*T=",f12.5,1x,"exp=",f12.5)') self%fc(1),self%fc(2)
    case (wall_fermi)
      art = 'wall_fermi'
      write (atoms,'(1x,"atoms:",1x,i0,a)') self%n,'/all'
      write (chnl,headfmt) trim(art),trim(atoms)
      write (chnl,'(2x,"radii(AA)=",3f12.5)') self%ref(1:3)*autoaa
      write (chnl,'(2x,"     kb*T=",f12.5,1x," exp=",f12.5)') self%fc(1)*kB,self%fc(2)
    case (na_gapdiff)
      art = 'nonadiabatic gap'
      write (atoms,'(1x,"[",a,"]")') 'σ*ΔE²/(ΔE+α)'
      write (chnl,headfmt) trim(art),trim(atoms)
      write (chnl,'(2x,"σ=",f8.5," α=",f8.5)') self%fc(1:2)
    case (na_gapdiff2)
      art = 'nonadiabatic gap'
      write (atoms,'(1x,"[",a,"]")') 'σ*(exp(-β|ΔE|)+C) * ΔE²/(|ΔE|+α)'
      write (chnl,headfmt) trim(art),trim(atoms)
      write (chnl,'(2x,"σ=",f8.5," α=",f8.5)') self%fc(1:2)
      write (chnl,'(2x,"C=",f8.5," β=",f8.5)') self%fc(3),27.2114_wp
    case (bondrange)
      art = 'bondrange'
      write (atoms,'(1x,"atoms:",1x,i0,",",i0)') self%atms(1:2)
      write (chnl,headfmt) trim(art),trim(atoms)
      write (chnl,'(2x,"upper=",f8.2,1x,"lower=",f8.2)') self%ref(1)*autoaa,self%ref(2)*autoaa
    case default
      art = 'none'
      atoms = 'none'
      values = ' '
      pr = .false.
    end select
    !if (pr) then
    !  write (chnl,'("> ",a,a,a)') 'constraint: ',trim(art),trim(atoms)
    !  write (chnl,'(1x,a)') trim(values)
    !end if

    return
  end subroutine print_constraint

!========================================================================================!
!> subroutine analyze_dummy_bond_constraint
  subroutine analyze_dummy_bond_constraint(self,t,i,rawa)
    implicit none
    class(constraint) :: self
    integer,intent(in) :: i,t
    character(len=*) :: rawa(i)
    real(wp) :: k,dist
    integer :: a1,a2

    call self%deallocate()

    select case (i)
    case (1)
      if (any(rawa(1) == (/character(7)::'all','allauto'/))) then
        call self%dummyconstraint(t)
      end if
    case (2)
      if (any(rawa(1) == (/character(9)::'all','allauto'/))) then
        read (rawa(2),*) k
        call self%dummyconstraint(t,k)
      end if
    case (3)
      read (rawa(1),*) a1
      read (rawa(2),*) a2
      read (rawa(3),*) dist
      dist = abs(dist)
      call self%bondconstraint(a1,a2,dist)
    case (4)
      read (rawa(1),*) a1
      read (rawa(2),*) a2
      read (rawa(3),*) dist
      read (rawa(4),*) k
      dist = abs(dist)
      call self%bondconstraint(a1,a2,dist,k)
    end select

    return
  end subroutine analyze_dummy_bond_constraint

  subroutine analyze_dummy_bond_constraint2(self,i,fa)
    implicit none
    class(constraint) :: self
    integer,intent(in) :: i
    real(wp),intent(in) :: fa(i)
    real(wp) :: k,dist
    integer :: a1,a2

    call self%deallocate()

    select case (i)
    case (3)
      a1 = nint(fa(1))
      a2 = nint(fa(2))
      dist = abs(fa(3))
      call self%bondconstraint(a1,a2,dist)
    case (4)
      a1 = nint(fa(1))
      a2 = nint(fa(2))
      dist = abs(fa(3))
      k = abs(fa(4))
      call self%bondconstraint(a1,a2,dist,k)
    end select

    return
  end subroutine analyze_dummy_bond_constraint2

!========================================================================================!
!> subroutine constraint_deallocate
!> reset and deallocate all data of a given constraint object
  subroutine constraint_deallocate(self)
    implicit none
    class(constraint) :: self
    if (allocated(self%atms)) deallocate (self%atms)
    if (allocated(self%fc)) deallocate (self%fc)
    if (allocated(self%ref)) deallocate (self%ref)
    self%type = 0
    self%subtype = pharmonic
    self%n = 0
    self%frozenatms = .false.
    self%wscal = 1.0_wp
    return
  end subroutine constraint_deallocate

!========================================================================================!
!> subroutine constraint_freezeassoc
!> associate the freezeptr
  subroutine constraint_freezeassoc(self,freezelist)
    implicit none
    class(constraint) :: self
    logical,intent(in),target :: freezelist(:)
    self%freezeptr => freezelist
    self%frozenatms = .true.
    return
  end subroutine constraint_freezeassoc

!========================================================================================!
!> subroutine scantype_deallocate
!> reset and deallocate all data of a given scantype object
  subroutine scantype_deallocate(self)
    implicit none
    class(scantype) :: self
    if (allocated(self%atms)) deallocate (self%atms)
    self%type = 0
    self%steps = 5
    self%n = 0
    self%minval = 0.0_wp
    self%maxval = 0.0_wp
    if (allocated(self%points)) deallocate (self%points)
    self%restore = 1
    return
  end subroutine scantype_deallocate

!========================================================================================!
!> subroutien create_bond_constraint
  subroutine create_dummy_constraint(self,t,k)
    implicit none
    class(constraint) :: self
    integer,intent(in) :: t
    real(wp),optional :: k

    call self%deallocate()
    self%type = t
    self%n = 0
    allocate (self%fc(1),source=fcdefault)
    if (present(k)) then
      self%fc(1) = k
    end if
    return
  end subroutine create_dummy_constraint

!========================================================================================!
!> subroutien create_bond_constraint
  subroutine create_bond_constraint(self,i,j,d,k)
    implicit none
    class(constraint) :: self
    integer,intent(in) :: i,j
    real(wp),intent(in) :: d
    real(wp),optional :: k

    call self%deallocate()
    self%type = bond
    self%n = 2
    allocate (self%atms(2))
    allocate (self%fc(1),source=fcdefault)
    allocate (self%ref(1))
    self%atms(1) = i
    self%atms(2) = j
    self%ref(1) = d
    if (present(k)) then
      self%fc(1) = k
    end if
    return
  end subroutine create_bond_constraint

!========================================================================================!

  subroutine bond_constraint(n,xyz,constr,energy,grd)
!**************************************************
!* constrain the distance between two atoms A...B
!* by an harmonic potential V(r) = 1/2kr²
!* r and k are in atomic units (Bohr, Hartree)
!**************************************************
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    type(constraint) :: constr

    real(wp),intent(out) :: energy
    real(wp),intent(out) :: grd(3,n)
    integer :: iat,jat
    real(wp) :: a,b,c,x
    real(wp) :: dist,ref,k,dum
    real(wp) :: T,beta

    energy = 0.0_wp
    grd = 0.0_wp

    if (constr%n /= 2) return

    iat = constr%atms(1)
    jat = constr%atms(2)
    a = xyz(1,iat)-xyz(1,jat)
    b = xyz(2,iat)-xyz(2,jat)
    c = xyz(3,iat)-xyz(3,jat)
    dist = sqrt(a**2+b**2+c**2)
    ref = constr%ref(1)
    k = constr%fc(1)

    x = dist-ref
    select case (constr%subtype)
    case (pharmonic)
      energy = 0.5_wp*k*(x)**2
      dum = k*x
    case (plogfermi)
      energy = kb*T*log(1.0_wp+exp(beta*x))
      dum = (kb*T*beta*exp(beta*x))/(exp(beta*x)+1.0_wp)
    end select

    grd(1,iat) = dum*(a/dist)
    grd(2,iat) = dum*(b/dist)
    grd(3,iat) = dum*(c/dist)
    grd(1,jat) = -grd(1,iat)
    grd(2,jat) = -grd(2,iat)
    grd(3,jat) = -grd(3,iat)

    return
  end subroutine bond_constraint

!========================================================================================!

  subroutine create_bondrange_constraint(self,i,j,dup,dlow,beta,T)
    implicit none
    class(constraint) :: self
    integer,intent(in) :: i,j
    real(wp),intent(in) :: dup,dlow
    real(wp),optional :: beta
    real(wp),optional :: T

    call self%deallocate()
    self%type = bondrange
    self%n = 2
    allocate (self%atms(2))
    allocate (self%fc(2),source=fcdefault)
    allocate (self%ref(2))
    self%atms(1) = i
    self%atms(2) = j
    self%ref(1) = max(dup,dlow)
    self%ref(2) = min(dup,dlow)
    if (present(T)) then
      self%fc(1) = T/kb
    else
      self%fc(1) = 0.1_wp/kb
    end if
    self%fc(1) = abs(self%fc(1))
    if (present(beta)) then
      self%fc(2) = beta
    else
      self%fc(2) = betadefault
    end if
    return
  end subroutine create_bondrange_constraint

  subroutine bondrange_constraint(n,xyz,constr,energy,grd)
!************************************************************
!* constrain the distance between two atoms A...B
!* via two logfermi potentials V(r) = kb*T*log(1+e^(beta*r))
!* This potential allows to define an upper and lower bound
!* for the AB distance
!************************************************************
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    type(constraint) :: constr

    real(wp),intent(out) :: energy
    real(wp),intent(out) :: grd(3,n)
    integer :: iat,jat
    real(wp) :: a,b,c,x
    real(wp) :: dist,ref,k,dum
    real(wp) :: ref_upper,ref_lower
    real(wp) :: T,beta
    real(wp) :: dr(3)

    energy = 0.0_wp
    grd = 0.0_wp

    if (constr%n /= 2) return

    iat = constr%atms(1)
    jat = constr%atms(2)
    a = xyz(1,iat)-xyz(1,jat)
    b = xyz(2,iat)-xyz(2,jat)
    c = xyz(3,iat)-xyz(3,jat)
    dist = sqrt(a**2+b**2+c**2)
    ref_upper = constr%ref(1)
    ref_lower = constr%ref(2)
    T = constr%fc(1)
    beta = constr%fc(2)

    !> upper bound contribution
    x = dist-ref_upper
    !energy = kb*T*log(1.0_wp+exp(beta*x))
    energy = logfermi(T,beta,x)
    !dum = (kb*T*beta*exp(beta*x))/(exp(beta*x)+1.0_wp)
    dum = dlogfermi(T,beta,x)
    grd(1,iat) = dum*(a/dist)
    grd(2,iat) = dum*(b/dist)
    grd(3,iat) = dum*(c/dist)
    grd(1,jat) = -grd(1,iat)
    grd(2,jat) = -grd(2,iat)
    grd(3,jat) = -grd(3,iat)

    !> lower bound contribution
    x = ref_lower-dist
    !energy = energy +  kb*T*log(1.0_wp+exp(beta*x))
    energy = energy+logfermi(T,beta,x)
    !dum = (kb*T*beta*exp(beta*x))/(exp(beta*x)+1.0_wp)
    dum = dlogfermi(T,beta,x)
    dr(1) = -dum*(a/dist)
    dr(2) = -dum*(b/dist)
    dr(3) = -dum*(c/dist)
    grd(1,iat) = grd(1,iat)+dr(1)
    grd(2,iat) = grd(2,iat)+dr(2)
    grd(3,iat) = grd(3,iat)+dr(3)
    grd(1,jat) = grd(1,jat)-dr(1)
    grd(2,jat) = grd(2,jat)-dr(2)
    grd(3,jat) = grd(3,jat)-dr(3)

    !> energy shift (no gard contribution)
    !> if ref_upper == ref_lower, this will transform the
    !> bondrange constraint into something close to a harmonic potential
    !> and the shift shifts the potential to zero at the minimum
    ref = (ref_upper+ref_lower)/2.0_wp
    x = ref-ref_upper
    dum = logfermi(T,beta,x)
    x = ref_lower-ref
    dum = dum+logfermi(T,beta,x)
    energy = energy-dum

    return
  end subroutine bondrange_constraint

  function logfermi(T,beta,x) result(energy)
    implicit none
    real(wp) :: energy
    real(wp),intent(in) :: T
    real(wp),intent(in) :: beta
    real(wp),intent(in) :: x
    energy = kb*T*log(1.0_wp+exp(beta*x))
  end function logfermi

  function dlogfermi(T,beta,x) result(dldx)
    implicit none
    real(wp) :: dldx
    real(wp),intent(in) :: T
    real(wp),intent(in) :: beta
    real(wp),intent(in) :: x
    dldx = (kb*T*beta*exp(beta*x))/(exp(beta*x)+1.0_wp)
  end function dlogfermi

!========================================================================================!
!> subroutien create_angle_constraint
  subroutine create_angle_constraint(self,a,b,c,d,k)
    implicit none
    class(constraint) :: self
    integer,intent(in) :: a,b,c
    real(wp),intent(in) :: d ! constrain angle in degrees
    real(wp),optional :: k
    real(wp) :: dum,d2
    call self%deallocate()
    self%type = angle
    self%n = 3
    allocate (self%atms(3))
    allocate (self%fc(1),source=fcdefault)
    allocate (self%ref(1))
    self%atms(1) = a
    self%atms(2) = b
    self%atms(3) = c
    d2 = abs(d)
    if (d2 > 360.0_wp) then
      dum = d2
      do
        dum = dum-360.0_wp
        if (dum < 360.0_wp) then
          d2 = dum
          exit
        end if
      end do
    end if
    if (d2 > 180.0_wp) then
      d2 = 360.0_wp-d2
    end if
    self%ref(1) = d2/deg !reference in rad
    if (present(k)) then
      self%fc(1) = k
    end if
    return
  end subroutine create_angle_constraint

!========================================================================================!
!> subroutine angle_constraint
!> constrain angle between atoms A and C, connected via a central atom B:  A-B-C
!> using a harmonic potential
  subroutine angle_constraint(n,xyz,constr,energy,grd)
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    type(constraint) :: constr
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: grd(3,n)
    integer :: i,iat,jat,kat
    real(wp) :: A(3),B(3),C(3)
    real(wp) :: r1(3),r2(3)
    real(wp) :: angle,k,ref,p,d,l1,l2
    real(wp) :: dinv,dum,x,T,beta
    real(wp) :: dadA(3),dadB(3),dadC(3)

    energy = 0.0_wp
    grd = 0.0_wp

    if (constr%n /= 3) return
    if (.not.allocated(constr%atms)) return
    if (.not.allocated(constr%ref)) return
    if (.not.allocated(constr%fc)) return

    ref = constr%ref(1)
    k = constr%fc(1)
    iat = constr%atms(1)
    jat = constr%atms(2)
    kat = constr%atms(3)
    A = xyz(:,iat)
    B = xyz(:,jat)
    C = xyz(:,kat)
    call angle_and_derivatives(A,B,C,angle,dadA,dadB,dadC)

    x = angle-ref
    select case (constr%subtype)
    case (pharmonic) !> harmonic potential
      energy = 0.5_wp*k*(x)**2
      dum = k*(x)
    case (plogfermi) !> logfermi potential
      energy = kb*T*log(1.0_wp+exp(beta*x))
      dum = (kb*T*beta*exp(beta*x))/(exp(beta*x)+1.0_wp)
    end select

    grd(1,iat) = dum*dadA(1)
    grd(2,iat) = dum*dadA(2)
    grd(3,iat) = dum*dadA(3)
    grd(1,jat) = dum*dadB(1)
    grd(2,jat) = dum*dadB(2)
    grd(3,jat) = dum*dadB(3)
    grd(1,kat) = dum*dadC(1)
    grd(2,kat) = dum*dadC(2)
    grd(3,kat) = dum*dadC(3)

    return
  end subroutine angle_constraint

  subroutine angle_and_derivatives(A,B,C,angle,dadA,dadB,dadC)
    implicit none
    real(wp),intent(in) :: A(3),B(3),C(3) !> points spanning the angle
    real(wp),intent(out) :: angle !> the angle in rad
    real(wp),intent(out) :: dadA(3),dadB(3),dadC(3) !> Cartesian derivatives
    real(wp) :: r1(3),r2(3)
    real(wp) :: p,d,l1,l2
    real(wp) :: dinv,dum

    angle = 0.0_wp
    dadA = 0.0_wp
    dadB = 0.0_wp
    dadC = 0.0_wp

    r1 = A-B
    r2 = C-B
    l1 = rlen(r1)
    l2 = rlen(r2)
    p = dot(r1,r2)
    d = p/(l1*l2)
    angle = acos(d)
    if (angle < 1d-6.or.(pi-angle) < 1d-6) then
      dadA(:) = (1.0_wp/l2)*sin(acos(r2(:)/l2))
      dadC(:) = (1.0_wp/l1)*sin(acos(r1(:)/l1))
      if ((pi-angle) < 1d-6) then
        dadA = -dadA
        dadC = -dadC
      end if
    else
      dinv = 1.0_wp/sqrt(1.0_wp-d**2)
      dadA(:) = -dinv*(r2(:)*l1*l2-p*(l2/l1)*r1(:))/(l1**2*l2**2)
      dadC(:) = -dinv*(r1(:)*l1*l2-p*(l1/l2)*r2(:))/(l1**2*l2**2)
    end if
    dadB = -dadA-dadC

    return
  end subroutine angle_and_derivatives

  real(wp) function rlen(r)
    implicit none
    real(wp) :: r(3)
    rlen = 0.0_wp
    rlen = r(1)**2+r(2)**2+r(3)**2
    rlen = sqrt(rlen)
    return
  end function rlen
  real(wp) function dot(r1,r2)
    implicit none
    real(wp) :: r1(3),r2(3)
    dot = 0.0_wp
    dot = r1(1)*r2(1)+r1(2)*r2(2)+r1(3)*r2(3)
    return
  end function dot
  subroutine cross(r1,r2,r3)
    implicit none
    real(wp) :: r1(3),r2(3)
    real(wp) :: r3(3)
    r3 = 0.0_wp
    r3(1) = r1(2)*r2(3)-r1(3)*r2(2)
    r3(2) = r1(3)*r2(1)-r1(1)*r2(3)
    r3(3) = r1(1)*r2(2)-r1(2)*r2(1)
    return
  end subroutine cross

!========================================================================================!
!> subroutien create_dihedral_constraint
!> constrain angle in degrees (°), input should be between -180° and 180°
  subroutine create_dihedral_constraint(self,a,b,c,d,ref,k)
    implicit none
    class(constraint) :: self
    integer,intent(in) :: a,b,c,d
    real(wp),intent(in) :: ref !> constrain angle in degrees
    real(wp),optional :: k
    real(wp) :: dum,d2,sig
    call self%deallocate()
    self%type = dihedral
    self%n = 4
    allocate (self%atms(4))
    allocate (self%fc(1),source=fcdefault)
    allocate (self%ref(1))
    self%atms(1) = a
    self%atms(2) = b
    self%atms(3) = c
    self%atms(4) = d

    d2 = ref
    sig = sign(1.0_wp,ref)
    if (abs(d2) > 360.0_wp) then
      dum = abs(d2)
      do
        dum = dum-360.0_wp
        if (dum < 360.0_wp) then
          d2 = dum
          exit
        end if
      end do
      d2 = d2*sig
    end if
    if (d2 > 180.0_wp) then
      d2 = d2-360.0_wp
    end if
    if (d2 < -180.0_wp) then
      d2 = d2+360.0_wp
    end if
    self%ref(1) = d2/deg !reference in rad
    if (present(k)) then
      self%fc(1) = k
    end if
    return
  end subroutine create_dihedral_constraint
!========================================================================================!
!> subroutine dihedral_constraint
!> constrain dihedral angle spanned by atoms A-B-C-D
!> using a harmonic potential
  subroutine dihedral_constraint(n,xyz,constr,energy,grd)
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    type(constraint) :: constr
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: grd(3,n)
    integer :: i,iat,jat,kat,lat
    real(wp) :: A(3),B(3),C(3),D(3)
    real(wp) :: N1(3),N2(3),Nzero(3)
    real(wp) :: rab(3),rcb(3),rdc(3),na,nb,nc
    real(wp) :: dangle,k,ref,p,l1,l2
    real(wp) :: dinv,dum,x,T,beta
    real(wp) :: dadN1(3),dadN2(3),dad0(3)
    real(wp) :: sig,dDdr(3)
    real(wp) :: dDdA(3),dDdB(3),dDdC(3),dDdD(3)

    energy = 0.0_wp
    grd = 0.0_wp

    if (constr%n /= 4) return
    if (.not.allocated(constr%atms)) return
    if (.not.allocated(constr%ref)) return
    if (.not.allocated(constr%fc)) return

    ref = constr%ref(1)
    k = constr%fc(1)
    iat = constr%atms(1)
    jat = constr%atms(2)
    kat = constr%atms(3)
    lat = constr%atms(4)
    A = xyz(:,iat)
    B = xyz(:,jat)
    C = xyz(:,kat)
    D = xyz(:,lat)
    Nzero = 0.0_wp
    !> vectors spanning the planes (A,B,C) and (D,C,B)
    rab = A-B
    rcb = C-B
    rdc = D-C
    !> get the two normal vectors N1 and N2 for the two planes
    call cross(rab,rcb,N1)
    call cross(rdc,rcb,N2)
    p = dot(N1,rdc)
    sig = -sign(1.0_wp,p)
    call angle_and_derivatives(N1,Nzero,N2,dangle,dadN1,dad0,dadN2)
    dangle = sig*dangle

    x = dangle-ref
    if (x < -180.0_wp/deg) x = x+360.0_wp/deg
    if (x > 180.0_wp/deg) x = x-360.0_wp/deg
    select case (constr%subtype)
    case (pharmonic) !> harmonic potential
      energy = 0.5_wp*k*(x)**2
      dum = k*(x)
    case (plogfermi) !> logfermi potential
      energy = kb*T*log(1.0_wp+exp(beta*x))
      dum = (kb*T*beta*exp(beta*x))/(exp(beta*x)+1.0_wp)
    end select

    call dtorsdr(A,B,C,D,dadN1,dadN2,sig,dDdA,dDdB,dDdC,dDdD)
    grd(1:3,iat) = dum*dDdA(1:3)
    grd(1:3,jat) = dum*dDdB(1:3)
    grd(1:3,kat) = dum*dDdC(1:3)
    grd(1:3,lat) = dum*dDdD(1:3)

    return
  contains
    subroutine dtorsdr(A,B,C,D,dadN1,dadN2,sig,dphidA,dphidB,dphidC,dphidD)
      implicit none
      real(wp),intent(in) :: A(3),B(3),C(3),D(3) !>points spanning the torsion A-B-C-D
      real(wp),intent(in) :: dadN1(3),dadN2(3) !> derivatives of the angle between N1 and N2
      real(wp),intent(in) :: sig !> sign
      real(wp),intent(out) :: dphidA(3),dphidB(3),dphidC(3),dphidD(3) !> Cartesian derivatives

      dphidA(1) = sig*(dadN1(2)*(B(3)-C(3))+dadN1(3)*(C(2)-B(2)))
      dphidA(2) = sig*(dadN1(1)*(C(3)-B(3))+dadN1(3)*(B(1)-C(1)))
      dphidA(3) = sig*(dadN1(1)*(B(2)-C(2))+dadN1(2)*(C(1)-B(1)))

      dphidB(1) = sig*(dadN1(2)*(C(3)-A(3))+dadN1(3)*(A(2)-C(2)) &
      &       +dadN2(2)*(C(3)-D(3))+dadN2(3)*(D(2)-C(2)))
      dphidB(2) = sig*(dadN1(1)*(A(3)-C(3))+dadN1(3)*(C(1)-A(1)) &
      &       +dadN2(1)*(D(3)-C(3))+dadN2(3)*(C(1)-D(1)))
      dphidB(3) = sig*(dadN1(1)*(C(2)-A(2))+dadN1(2)*(A(1)-C(1)) &
      &       +dadN2(1)*(C(2)-D(2))+dadN2(2)*(D(1)-C(1)))

      dphidC(1) = sig*(dadN1(2)*(A(3)-B(3))+dadN1(3)*(B(2)-A(2)) &
      &       +dadN2(2)*(D(3)-B(3))+dadN2(3)*(B(2)-D(2)))
      dphidC(2) = sig*(dadN1(1)*(B(3)-A(3))+dadN1(3)*(A(1)-B(1)) &
      &       +dadN2(1)*(B(3)-D(3))+dadN2(3)*(D(1)-B(1)))
      dphidC(3) = sig*(dadN1(1)*(A(2)-B(2))+dadN1(2)*(B(1)-A(1)) &
      &       +dadN2(1)*(D(2)-B(2))+dadN2(2)*(B(1)-D(1)))

      dphidD(1) = sig*(dadN2(2)*(B(3)-C(3))+dadN2(3)*(C(2)-B(2)))
      dphidD(2) = sig*(dadN2(1)*(C(3)-B(3))+dadN2(3)*(B(1)-C(1)))
      dphidD(3) = sig*(dadN2(1)*(B(2)-C(2))+dadN2(2)*(C(1)-B(1)))

      return
    end subroutine dtorsdr

  end subroutine dihedral_constraint

!========================================================================================!
!> subroutine create_sphere_constraint

  subroutine create_sphere_constraint_all(self,n,r,k,alpha,logfermi)
    implicit none
    class(constraint) :: self
    integer,intent(in) :: n
    logical,allocatable :: atms(:)
    real(wp),intent(in) :: r
    real(wp) :: k,alpha
    logical,intent(in) :: logfermi
    integer :: i,c

    allocate (atms(n),source=.true.)
    call create_sphere_constraint(self,n,atms,r,k,alpha,logfermi)
    deallocate (atms)
    return
  end subroutine create_sphere_constraint_all

  subroutine create_sphere_constraint(self,n,atms,r,k,alpha,logfermi)
    implicit none
    class(constraint) :: self
    integer,intent(in) :: n
    logical,intent(in) :: atms(n)
    real(wp),intent(in) :: r
    real(wp) :: k,alpha
    logical,intent(in) :: logfermi
    integer :: i,c,ii

    call self%deallocate()
    if (logfermi) then
      self%type = wall_fermi
    else
      self%type = wall
    end if
    c = count(atms,1)
    self%n = c
    allocate (self%atms(c))
    allocate (self%fc(2),source=fcdefault)
    allocate (self%ref(3),source=r)
    ii = 0
    do i = 1,n
      if (atms(i))then 
        ii = ii +1
        self%atms(ii) = i
      endif
    end do
    self%ref(:) = r
    self%fc(1) = k
    self%fc(2) = alpha
    return
  end subroutine create_sphere_constraint

  subroutine create_ellips_constraint(self,n,atms,r,k,alpha,logfermi)
    implicit none
    class(constraint) :: self
    integer,intent(in) :: n
    logical,intent(in) :: atms(n)
    real(wp),intent(in) :: r(3)
    real(wp) :: k,alpha
    logical,intent(in) :: logfermi
    integer :: i,c,ii

    call self%deallocate()
    if (logfermi) then
      self%type = wall_fermi
    else
      self%type = wall
    end if
    c = count(atms,1)
    self%n = c
    allocate (self%atms(c))
    allocate (self%fc(2),source=fcdefault)
    allocate (self%ref(3),source=r)
    ii = 0
    do i = 1,n
      if (atms(i))then
        ii = ii +1 
        self%atms(ii) = i
      endif
    end do
    self%ref(:) = r(:)
    self%fc(1) = k
    self%fc(2) = alpha
    return
  end subroutine create_ellips_constraint

  subroutine sphere_update_nat(self,n,atms)
    implicit none
    class(constraint) :: self
    integer,intent(in) :: n
    logical,intent(in) :: atms(n)
    integer :: c,i
    c = count(atms,1)
    self%n = c
    if (allocated(self%atms)) deallocate (self%atms)
    allocate (self%atms(c))
    do i = 1,n
      if (atms(i)) self%atms(i) = i
    end do
    return
  end subroutine sphere_update_nat
!========================================================================================!

  subroutine wall_constraint(n,xyz,constr,energy,grd,subtype)
!*****************************************************************
!* constrain atoms within defined wall potentials
!* the potentials themselves can be polinomial or logfermi type
!*****************************************************************
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    type(constraint) :: constr
    integer,intent(in) :: subtype
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: grd(3,n)
    integer :: i,iat,jat
    real(wp) :: a,b,c,x,y,z,dx,dy,dz
    real(wp) :: dist,ddist,ref
    real(wp) :: k,alpha,dalpha,T,beta
    real(wp) :: fermi,expo,r(3),w(3)

    energy = 0.0_wp
    grd = 0.0_wp

    if (.not.allocated(constr%atms)) return
    if (.not.allocated(constr%ref)) return
    if (.not.allocated(constr%fc)) return

    !>--- xtb defaults are:
    !> sphere_alpha = 30
    !> sphere_beta  = 6.0_wp
    !> sphere_temp  = 300.0_wp

    do i = 1,constr%n
      iat = constr%atms(i)

!> DO NOT APPLY THE WALL POTENTIAL TO FROZEN ATOMS. THE ENERGY WILL JUST EXPLODE
      if (constr%frozenatms) then
        if (constr%freezeptr(iat)) cycle
      end if

      select case (subtype)
      case default
        return
      case (wall)
        !>
        !> V = k*Σ(|R-O|/Rref)^α
        !>
        k = constr%fc(1)
        alpha = constr%fc(2)
        dalpha = alpha-1.0_wp
        x = xyz(1,iat)
        y = xyz(2,iat)
        z = xyz(3,iat)
        a = constr%ref(1)
        b = constr%ref(2)
        c = constr%ref(3)
        dist = (x/a)**2+(y/b)**2+(z/c)**2
        energy = energy+k*(dist**alpha)
        dx = 2.0_wp*(x/(a**2))
        dy = 2.0_wp*(y/(b**2))
        dz = 2.0_wp*(z/(c**2))
        ddist = k*alpha*(dist**dalpha)
        grd(1,iat) = ddist*dx
        grd(2,iat) = ddist*dy
        grd(3,iat) = ddist*dz

      case (wall_fermi)
        !>
        !> V = Σ kT*log{1+exp[β(|R-O|-Rref)]}
        !>
        T = constr%fc(1)
        beta = constr%fc(2)
        ref = maxval(constr%ref(1:3))
        w(1:3) = ref/constr%ref(1:3)
        r = w*(xyz(1:3,iat))
        dist = sqrt(sum(r**2))
        expo = exp(beta*(dist-ref))
        fermi = 1.0_wp/(1.0_wp+expo)
        energy = energy+kB*T*log(1.0_wp+expo)
        grd(:,iat) = grd(:,iat)+kB*T*beta*expo*fermi*(r*w)/(dist+1.0e-14_wp)

      case (box)

      case (box_fermi)

      end select
    end do

    return
  end subroutine wall_constraint

!========================================================================================!
!> subroutien create_gapdiff_constraint
!> (calculation of the constraint in nonadiabatic.f90)
  subroutine create_gapdiff_constraint(self,sigm,alph)
    implicit none
    class(constraint) :: self
    real(wp),intent(in) :: sigm,alph

    call self%deallocate()
    self%type = na_gapdiff
    allocate (self%fc(2),source=fcdefault)
    self%fc(1) = sigm
    self%fc(2) = alph
    return
  end subroutine create_gapdiff_constraint

!========================================================================================!
!> subroutien create_gapdiff_constraint2
!> (calculation of the constraint in nonadiabatic.f90)
  subroutine create_gapdiff_constraint2(self,sigm,alph,c)
    implicit none
    class(constraint) :: self
    real(wp),intent(in) :: sigm,alph,c

    call self%deallocate()
    self%type = na_gapdiff2
    allocate (self%fc(3),source=fcdefault)
    self%fc(1) = sigm
    self%fc(2) = alph
    self%fc(3) = c
    return
  end subroutine create_gapdiff_constraint2

!========================================================================================!
!========================================================================================!
end module constraints
