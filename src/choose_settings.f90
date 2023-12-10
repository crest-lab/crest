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

!> The routines contained here take care of setting MTD runtimes
!> and how to set corresponding bias parameters.
!> Naturally, this is quite imporant for the overall performance of CREST

!=========================================================================================!
!=========================================================================================!
subroutine md_length_setup(env)
!***********************************************************
!* set the total run time according to flexibility measures
!***********************************************************
  use crest_parameters
  use crest_data
  use strucrd
  use zdata, only:readwbo
  implicit none
  !> IN/OUTPUT
  type(systemdata) :: env    !> MAIN STORAGE OS SYSTEM DATA
  !> LOCAL
  real(wp) :: total,minimum,lenthr
  real(wp) :: flex,av1,rfac,nciflex
  type(coord) :: mol 
  logical :: ex
!> get reference geometry  
  call env%ref%to( mol ) 

!> at least 5ps per MTD
  minimum = 5.0d0
!> Maximum of 200 ps, longer runs can only be conducted by user input
  lenthr = 200.0d0

  call smallhead('Generating MTD length from a flexibility measure')

  if ((env%crestver .ne. crest_solv).and..not.env%NCI) then
    write(stdout,'(1x,a)',advance='no') 'Calculating GFN0-xTB WBOs   ...'
!>-- xtb singlepoint to get WBOs (always GFN0)
    call xtbsp(env,0)
    write (stdout,'(1x,a)') 'done.'
!>-- save those WBOs to the reference
    inquire(file='wbo',exist = ex)
    if(ex)then
    if(.not.allocated(env%ref%wbo)) allocate(env%ref%wbo( mol%nat, mol%nat), source=0.0_wp)   
    call readwbo('wbo',mol%nat, env%ref%wbo)
    endif

!>-- covalent flexibility measure based on WBO and structure only
    call flexi( mol, env%rednat, env%includeRMSD, flex)
!>-- NCI flexi based on E(HB)/Nat and E(disp)/Nat
    call nciflexi(env,nciflex)
    write (stdout,'(1x,''    covalent flexibility measure :'',f8.3)') flex
    write (stdout,'(1x,''non-covalent flexibility measure :'',f8.3)') nciflex
!>-- the NCI flex is only relevant if the covalent framework is flexible
    flex = 0.5*flex+0.5*nciflex*sqrt(flex)

    if (env%entropic) then
!>--- special case for entropy mode that depends more strongly on the flexibility than on size
!>--- -8 accounts for small systems having no conf.
      av1 = (flex**1.333333)*max(1,env%rednat-8)
!>--- Maximum of 3000 ps, longer runs can only be conducted by user input
      lenthr = 3000.0d0
!>--- minimum is 5 ps set above
      env%tmtd = 4.50d0*exp(0.165d0*av1)
    else
!>--- normal case
      av1 = (flex**1.000000)*max(1,env%rednat-8)
!>--- Maximum of 500 ps
      lenthr = 500.0d0
      env%tmtd = 3.0d0*exp(0.10d0*av1)
    end if
  else
    write (stdout,'(1x,"System flexiblity is set to 1.0 for NCI mode")')
    flex = 1.0d0
    env%tmtd = 0.10*(env%rednat+0.1*env%rednat*env%rednat)
  end if
  env%flexi = flex
  write (stdout,'(1x,''flexibility measure :'',f8.3)') env%flexi

!>-- rfac is used to scale the total MD length according to special runtypes
  total = max(minimum,env%tmtd)
  select case (env%runver)
  case (2,5,6,33)   !> "-quick","-squick","-mquick"
    rfac = 0.5d0
  case (3)          !> "-qcg"
    rfac = 0.25d0
  case (77)
    rfac = 1.50d0
  case (8)
    rfac = 2.0d0
  case (787878)
    rfac = 14.0d0/float(env%nmetadyn)
  case default  !> everything else 1=default, 4=NCI
    rfac = 1.0d0
  end select
!>-- additional user set scaling
  if (env%scallen) then
    write (stdout,'(1x,''t(MTD) based on flexibility :'',  f8.1)') env%tmtd*rfac
    rfac = env%mdlenfac*rfac
    write (stdout,'(1x,''MTD length is scaled by     :'',  f6.3)') env%mdlenfac
  end if

!>-- ONLY use generated MD length if not already set by the user
  if (env%mdtime .le. 0.0d0) then
    if(env%mddat%length_ps > 0.0_wp)then
      total = env%mddat%length_ps
      write(stdout,'(1x,"t(MTD) / ps set via calculator :",  f8.1)') total
    else if (total .gt. lenthr) then
      total = lenthr
      call mtdwarning(lenthr*rfac)
    end if
    env%mdtime = anint(total)*rfac
  else
    write (stdout,'(1x,''t(MTD) / ps set by command line  :'',  f8.1)') env%mdtime
  end if

  write (stdout,'(1x,''t(MTD) / ps    :'',  f8.1)') env%mdtime
  write (stdout,'(1x,''Σ(t(MTD)) / ps :'',  f8.1,'' ('',i0,'' MTDs)'')') &
  & env%mdtime*float(env%nmetadyn),env%nmetadyn

!> A MTD Vbias snapshot is taken every 1 ps
  env%metadlist(:) = ceiling(env%mdtime)

  return
end subroutine md_length_setup

!========================================================================================!
!========================================================================================!
subroutine defaultGF(env)
!************************************************************
!* Setmetadynamics default Guiding Force Parameter
!* There are different combinations depending on the runtype
!************************************************************
  use crest_parameters 
  use crest_data
  use filemod
  implicit none
  !> IN/OUTPUT
  type(systemdata) :: env
  !> LOCAL
  integer  :: ia,ik,na,nk,m,nmtdyn,nmtdynmax,nrem
  real(wp) :: alp,k
  real(wp) :: kstart,kinc
  real(wp) :: alpinc
  type(filetype) :: biasfile
  logical :: ex
  integer :: i,io
  character(len=:),allocatable :: atmp

  nrem = 0

  if (.not.env%metadynset) then
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    if (env%readbias) then
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      inquire (file='mtdbias',exist=ex)
      if (ex) then
        atmp = ''
        call biasfile%open('mtdbias')
        call biasfile%clearblanks()
        nmtdyn = biasfile%nlines
        call env%allocate(nmtdyn)
        do i = 1,nmtdyn
          atmp = getlarg(biasfile%line(i),1)
          read (atmp,*,iostat=io) k
          if (io == 0) then
            env%metadfac(i) = k*env%rednat
          end if
          atmp = getlarg(biasfile%line(i),2)
          read (atmp,*,iostat=io) alp
          if (io == 0) then
            env%metadexp(i) = alp
          end if
        end do
        !write(*,*) env%metadfac
        !write(*,*) env%metadexp
      else
        error stop "no file 'mtdbias'"
      end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    else
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      select case (env%runver)
!---------- "-quick","-squick"
      case (2,5) 
        na = 3
        nk = 2
        nmtdyn = na*nk
        alp = 1.2d0 ! start value alpha
        kstart = 0.002d0 ! start value k
        alpinc = 2.0 ! increment
        kinc = 2.0     ! increment
!---------- "-mquick"
      case (6) 
        na = 3
        nk = 2
        nmtdyn = na*nk
        alp = 1.0d0 ! start value alpha
        kstart = 0.002d0 ! start value k
        alpinc = 2.0 ! increment
        kinc = 2.0     ! increment
!---------- "-qcg"
      case (3) 
        na = 4
        nk = 3
        nmtdyn = na*nk
        alp = 1.0d0 ! start value alpha
        kstart = 0.00125d0 ! start value k
        alpinc = (3./2.) ! increment
        kinc = (3./2.)   ! increment
!---------- "-nci"
      case (4) 
        na = 3
        nk = 2
        nmtdyn = na*nk
        alp = 1.0d0 ! start value alpha
        kstart = 0.001d0 ! start value k
        alpinc = 2.0 ! increment
        kinc = 2.0     ! increment
!---------- "-singlerun"
      case (45) 
        na = 1
        nk = 1
        nmtdyn = na*nk
        alp = 1.0d0 ! start value alpha
        kstart = 0.001d0 ! start value k
        alpinc = 2.0 ! increment
        kinc = 2.0     ! increment
!---------- "-relax"
      case (33)
        na = 1
        nk = 3
        nmtdyn = na*nk
        alp = 0.8d0 ! start value alpha
        kstart = 0.0030 ! start value k
        alpinc = 2.0 ! increment
        kinc = 2.0     ! increment
!---------- "-compress"
      case (77) 
        na = 3
        nk = 3
        nmtdyn = na*nk
        alp = 1.61803 ! start value alpha
        kstart = 0.005d0 ! start value k
        alpinc = 1.61803 ! increment
        kinc = 2.0     ! increment
!--------- "search_1"
      case (crest_s1,crest_mecp) 
        na = 3
        nk = 3
        nmtdyn = (na*nk)
        alp = 1.61803 ! start value alpha
        kstart = 0.0075d0 ! start value k
        alpinc = 1.61803 ! increment
        kinc = 2.0     ! increment
!---------- "-nmtd <x>"
      case (787878)
        call gfdistribute(env%nmetadyn,nk,na,nrem)
        nmtdyn = na*nk !> nrem will be substracted at the end
        alp = 1.3 ! start value alpha
        kstart = 0.0050 ! start value k
        alpinc = (5./3.) ! increment
        kinc = 1.5d0     ! increment
!---------- "-entropy"
      case (111) 
        na = 6
        nk = 4
        nmtdyn = (na*nk)
        alp = 1.61803 ! start value alpha
        kstart = 0.0075d0 ! start value k
        alpinc = 1.61803 ! increment
        kinc = 2.0     ! increment
!---------- default
      case default
        if (env%iterativeV2) then  !for the default iterative mode
          !=======================================================!
          na = 4
          nk = 3
          nrem = 0
          nmtdyn = (na*nk)+2
          call env%allocate(nmtdyn)   !allocate k(Vbias) and α(Vbias)
          alp = 1.3 ! start value alpha
          kstart = 0.0030 ! start value k
          alpinc = (5./3.) ! increment
          kinc = 2.0     ! increment

          !-- two additional MTDs with extreme values
          env%metadfac(nmtdyn-1) = 0.001*env%rednat
          env%metadexp(nmtdyn-1) = 0.1

          env%metadfac(nmtdyn) = 0.005*env%rednat
          env%metadexp(nmtdyn) = 0.8
          !=======================================================!
        else  !for the non-iterative mode
          !=======================================================!
          na = 6
          nk = 4
          nmtdyn = na*nk
          alp = 1.3 ! start value alpha
          kstart = 0.003d00 ! start value k
          alpinc = (4./3.) ! increment
          kinc = (3./2.)     ! increment
          !======================================================!
        end if
      end select
!>---- settings are generated here
      m = 0
      nmtdynmax = nmtdyn-nrem
      call env%allocate(nmtdynmax)   !allocate k(Vbias) and α(Vbias)
      do ia = 1,na
        k = kstart  ! start value k
        do ik = 1,nk
          m = m+1
          if (m > nmtdynmax) cycle !> skip the last nrem setups
          env%metadfac(m) = k*env%rednat*env%mtd_kscal
          env%metadexp(m) = alp
          k = k/kinc ! increment
        end do
        alp = alp/alpinc ! increment
      end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++!
    end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++!
  end if
  return
!====================================================================!
contains
!====================================================================!

  subroutine gfdistribute(nsoll,k,a,rem)
!**********************************************************
!* subroutine gfdistribute
!* for a given number of requested MTDs (nsoll),
!* specify the number of different kpush (k) and
!* alpha (a) values, and from the product (k*a),
!* determine how many MTDs have to be neglected (nrem).
!* I.e.,
!*        (k*a)-nrem = nsoll
!*
!**********************************************************
    implicit none
    integer,intent(in) :: nsoll
    integer,intent(out) :: k,a,rem
    real(wp) :: nn,nsq
    k = 1
    a = 1
    rem = 0
    if (nsoll < 1) then
      error stop 'MTD setup failure.'
    end if
    nn = float(nsoll)
    nsq = sqrt(nn)
    a = ceiling(nsq)
    k = nint(nsq)
    rem = abs(nsoll-(a*k))
    return
  end subroutine gfdistribute
end subroutine defaultGF

!=========================================================================================!
!=========================================================================================!
subroutine adjustnormmd(env)
!*************************************************************************
!* Dynamically determine the number of normMDs and settings of staticMTDs
!* Set their number and the different temperatures.
!* Defaults for the static MTDs are more lengthy...
!*************************************************************************
  use crest_parameters
  use crest_data
  implicit none
  !> IN/OUTPUT
  type(systemdata) :: env
  !> LOCAL
  integer :: ndum

  if (env%rotamermds) then
!>--- first the number of normMDs on low conformers
    if (env%nrotammds .le. 0) then !> if no user input was set
      !> multiple short MDs, which has a better parallel efficiency
      !> default is 4 
      env%nrotammds = max(1,nint(float(env%nmetadyn)/4.0d0)) 
    end if

!>--- then the temperature range
    if (env%temps .le. 0) then
      !> at how many different temperatures? 
      !> starting at 400k and increasing 100K for each (200 K for -entropy mode)
      env%temps = 2
      if (env%entropic) then
        env%temps = 1
      end if
    end if
!>--- total number of NORMMDs is temps*nrotammds
  end if

!==============================================!
!>--- settings for static MTDS in entropy mode
!==============================================!
  if (env%entropymd) then 
    env%emtd%iter = 20    !> max number of iterations
    env%emtd%nbias = min(150,nint(env%tmtd/4)) !> max number of bias structures
    env%emtd%nbiasgrow = min(1.4d0,1.2d0+env%tmtd*1.d-3) !> increase of nBias in each cycle
    env%emtd%nMDs = 36          !> number of static MTDs
    env%emtd%lenfac = 0.5d0     !> length (relativ to regular MTDs)
    env%emtd%temperature = env%nmdtemp !> sMTD temperature (default 600 K)
    env%emtd%kpush = 1.d-4+env%tmtd*1.0d-6   !> kpush constant PER ATOM, a bit more for flexible systems 1.d-4+env%tmtd*1.d-6 1.5 zu viel, 0.5 zu wenig
    env%emtd%alpha = 1.0d0        !> some alpha
    env%emtd%mtdramp = 0.015d0    !> parameter to control how "fast" bias is applied in MTD
    if (env%crestver == crest_imtd) then
      if (env%emtd%confthr < 0.0d0) then
        env%emtd%confthr = 0.02d0    !> if we gain less than x% NEW conformers, exit
      end if
      if (env%emtd%sconvthr < 0.0d0) then
        env%emtd%sconvthr = 0.005d0   !> if we gain less than x% NEW entropy, exit
      end if
    end if

    if (env%nmdtemp < 0.d0) then   !> if temperature is not set by the user
      env%nmdtemp = 600.0d0
    end if

!>--- for the new alternative iMTD-sMTD runtype, re-adjust settings
    if (env%crestver == crest_imtd2) then
      ndum = 2
      env%emtd%nklist = ndum
      allocate (env%emtd%klist(ndum))
      env%emtd%klist(1) = env%emtd%kpush
      env%emtd%klist(2) = env%emtd%kpush*2.5d0
      env%emtd%nMDs = 12           !> number of static MTDs
      env%emtd%lenfac = 0.5d0      !> half the length because we have 2 kpush
      if (env%emtd%confthr < 0.0d0) then
        env%emtd%confthr = 0.05d0  !> if we gain less than x% NEW conformers, exit
      end if
      if (env%emtd%sconvthr < 0.0d0) then
        env%emtd%sconvthr = 0.01d0  !> if we gain less than x% NEW entropy, exit
      end if
    end if

!>--- Exclude atoms from static MTD bias
    env%emtd%rmax = 0    !> ignore small rings up to this size in bias
    call mtdatoms(env)
  end if
!==============================================!

  return
end subroutine adjustnormmd

!========================================================================================!
!========================================================================================!
subroutine env_to_mddat(env)
!**********************************************
!* Convert CREST's global MD settings to the
!* MD calculator object that will be used as
!* A basis for the newer calculator routines
!**********************************************
  use crest_parameters
  use crest_data
  implicit none
  type(systemdata) :: env
  real(wp) :: dum
!!>--- dont override user-defined settings
!  if(env%mddat%requested) return
!> we will check if any default settings were already set individually, instead
!> the if-statements in the following take care of that

!>--- necessary transfer global settings into mddat object
   if(env%mddat%length_ps <= 0.0_wp)then
   !> total runtime in ps
     env%mddat%length_ps    = env%mdtime
   else
     env%mdtime = env%mddat%length_ps
   endif
   if(env%mddat%tstep <= 0.0_wp)then
   !> time step in fs 
     env%mddat%tstep        = env%mdstep
   endif
   !> simulation steps (would be recovered automatically later, but just to make sure)
   env%mddat%length_steps = nint(env%mddat%length_ps*1000.0_wp / env%mddat%tstep)
   if(env%mddat%tsoll <= 0.0_wp)then
   !> target temperature
     env%mddat%tsoll = env%mdtemp
   endif

   if( env%mddat%dumpstep <= 0.0_wp ) then 
   !> dump frequency in fs
     env%mddat%dumpstep = float(env%mddumpxyz)
   endif
   if(env%mddat%sdump <= 0)then
   !> trajectory structure dump every x steps 
     dum = max(1.0_wp, (env%mddat%dumpstep / env%mddat%tstep))
     env%mddat%sdump = nint(dum)
   endif

   !> The SHAKE setup (special condition referring to the default)
   env%mddat%shake = env%mddat%shake .and.(env%shake > 0) !> SHAKE algorithm?
   if( env%mddat%shake .and. env%mddat%shk%shake_mode == 0)then
   env%mddat%shk%shake_mode = env%shake     !> H-only shake =1, all atom =2
   endif 

   if(env%mddat%md_hmass <= 0.0_wp)then
   !> hydrogen mass (to enable longer timesteps)
     env%mddat%md_hmass = env%hmass 
   endif

   ! TODO: WBO reader if shake is applied and wbo file is present

!>--- set flag to signal present settings
  env%mddat%requested = .true.

end subroutine env_to_mddat


