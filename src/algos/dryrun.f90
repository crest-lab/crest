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

subroutine crest_dry_run(env,tim)
!********************************************************************
!* Dry-run runtype. Prints a formatted summary of all CREST settings
!* and exits cleanly without performing any calculations.
!*
!* Input/Output:
!*  env  -  crest's systemdata object
!*  tim  -  timer object
!********************************************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  use iomod
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  character(len=512) :: dumstr
  character(len=:),allocatable :: ctmp
  logical :: ex
!========================================================================================!

  write (stdout,*)
  call drawbox(stdout,'D R Y   R U N',charset=6,width=80)
  write (stdout,'(1x,a)') 'Dry run was requested.'
  write (stdout,'(1x,a)') 'Running CREST with the chosen arguments would result in the following settings:'
  write (stdout,*)

!========================================================================================!
!> INPUT FILE
!========================================================================================!
  call drawbox(stdout,'Input',charset=4,padl=2,padr=2)
  write (stdout,*)
  ex = file_exists(env%inputcoords)
  write (stdout,'(2x,a,a)',advance='no') 'Input file : ',trim(env%inputcoords)
  if (ex) then
    write (stdout,*)
  else
    write (stdout,'(1x,"( ",a," )")') colorify('NOT FOUND','red')
  end if
  write (stdout,*)

!========================================================================================!
!> RUNTYPE
!========================================================================================!
  call drawbox(stdout,'Job type',charset=4,padl=2,padr=2)
  write (stdout,*)
  select case (env%crestver)
  case (crest_mfmdgc)
    write (stdout,'(2x,a)') 'Conformational search via the MF-MD-GC algorithm ('//colorify("DEPRECATED",'red')//')'
  case (crest_imtd)
    write (stdout,'(2x,a)') 'Conformational search via the iMTD-GC algorithm'
  case (crest_imtd2)
    write (stdout,'(2x,a)') 'Conformational search via the iMTD-sMTD algorithm (-v4)'
  case (crest_mdopt)
    write (stdout,'(2x,a)') 'Ensemble reoptimization (-mdopt)'
  case (crest_mdopt2)
    write (stdout,'(2x,a)') 'Ensemble reoptimization, variant 2 (-mdopt2)'
  case (crest_screen)
    write (stdout,'(2x,a)') 'Ensemble screening and reoptimization (-screen)'
  case (crest_nano)
    write (stdout,'(2x,a)') 'GFNn-xTB nano reactor (-reactor)'
  case (crest_sp)
    write (stdout,'(2x,a)') 'Standalone singlepoint calculation'
  case (crest_optimize)
    if (.not.env%crest_ohess) then
      write (stdout,'(2x,a)') 'Standalone geometry optimization'
    else
      write (stdout,'(2x,a)') 'Standalone geometry optimization followed by numerical Hessian' 
    end if
  case (crest_moldyn)
    write (stdout,'(2x,a)') 'Standalone molecular dynamics simulation'
  case (crest_s1)
    write (stdout,'(2x,a)') 'Conformational search (crest_s1)'
  case (crest_mecp)
    write (stdout,'(2x,a)') 'Minimum energy crossing point (MECP) search'
  case (crest_numhessian)
    write (stdout,'(2x,a)') 'Numerical Hessian calculation'
  case (crest_scanning)
    write (stdout,'(2x,a)') 'Coordinate scan'
  case (crest_rigcon)
    write (stdout,'(2x,a)') 'Rule-based conformer generation'
  case (crest_sorting)
    write (stdout,'(2x,a)') 'Standalone ensemble sorting (CREGEN)'
  case (crest_bh)
    write (stdout,'(2x,a)') 'Basin-hopping conformer search'
  case (crest_bhpt)
    write (stdout,'(2x,a)') 'Basin-hopping with parallel tempering'
  case (crest_none)
    write (stdout,'(2x,a)') '<no runtype specified>'
  case default
    write (stdout,'(2x,a,i0,a)') '<unrecognized runtype: crestver=',env%crestver,'>'
  end select
  write (stdout,*)

!========================================================================================!
!> CALCULATION SETTINGS
!========================================================================================!
  call drawbox(stdout,'Calculation settings',charset=4,padl=2,padr=2)
  write (stdout,*)
  if (associated(env%calc)) then
    if (env%calc%ncalculations > 0) then
      call env%calc%info(stdout,printhdr=.false.)
    else
      write (stdout,'(2x,a)') 'Calculation object associated but no levels defined yet.'
    end if
  else
    write (stdout,'(2x,a)') 'Calculation object not associated (legacy mode or not yet set up).'
  end if
  write (stdout,*)

!========================================================================================!
!> OPTIMIZATION SETTINGS
!========================================================================================!
  call drawbox(stdout,'Optimization settings',charset=4,padl=2,padr=2)
  write (stdout,*)

  write (stdout,'(2x,a,t35,": ",a,1x,"(",i0,")")') 'Optimization level',optlevflag(env%optlev),nint(env%optlev)
  if (associated(env%calc)) then
    block
      use optimize_utils,only:get_optthr
      real(wp) :: ethr,gthr
      integer :: nat,iolev
      nat = env%ref%nat
      iolev = nint(env%optlev)
      call get_optthr(nat,iolev,env%calc,ethr,gthr)
      write (stdout,'(2x,a,t35,": ",i0)') 'Max cycles (calc obj)',env%calc%maxcycle
      write (stdout,'(2x,a,t35,": ",es12.4)') 'Energy convergence  [Eh]',ethr
      write (stdout,'(2x,a,t35,": ",es12.4)') 'Gradient convergence [Eh/a0]',gthr
    end block
  end if
  write (stdout,*)

!========================================================================================!
!> MD / MTD SETTINGS
!========================================================================================!
  call drawbox(stdout,'MD / MTD settings',charset=4,padl=2,padr=2)
  write (stdout,*)
  if (env%mdtime > 0.0_wp) then
    write (stdout,'(2x,a,t35,": ",f10.1,a)') 'Simulation length',env%mdtime,' ps'
  else
    write (stdout,'(2x,a,t35,": ",a)') 'Simulation length','<system dependent>'
  end if
  write (stdout,'(2x,a,t35,": ",f10.1,a)') 'Time step',env%mdstep,' fs'
  write (stdout,'(2x,a,t35,": ",i10)') 'SHAKE mode',env%shake
  write (stdout,'(2x,a,t35,": ",f10.2,a)') 'MD temperature',env%mdtemp,' K'
  write (stdout,'(2x,a,t35,": ",i10,a)') 'Trajectory dump step',env%mddumpxyz,' fs'
  write (stdout,'(2x,a,t35,": ",f10.1,a)') 'MTD Vbias dump',real(env%mddump,wp)/1000.0_wp,' ps'
  if (env%mddat%length_ps > 0.0_wp) then
    write (stdout,*)
    write (stdout,'(2x,a)') 'mddata object (modern MD runtype):'
    write (stdout,'(4x,a,t35,": ",f10.1)') 'length_ps',env%mddat%length_ps
    write (stdout,'(4x,a,t35,": ",f10.4)') 'tstep [fs]',env%mddat%tstep
    write (stdout,'(4x,a,t35,": ",f10.2)') 'T_soll',env%mddat%tsoll
    write (stdout,'(4x,a,t35,": ",l6)') 'SHAKE',env%mddat%shake
    write (stdout,'(4x,a,t35,": ",a)') 'thermostat',trim(env%mddat%thermotype)
  end if
  write (stdout,*)

!========================================================================================!
!> THERMODYNAMICS SETTINGS
!========================================================================================!
  call drawbox(stdout,'Thermodynamics settings',charset=4,padl=2,padr=2)
  write (stdout,*)
  select case (env%thermo%emodel)
  case ('grimme')
    ctmp = 'Grimme (2012)'
  case ('truhlar')
    ctmp = 'Truhlar (2011)'
  end select
  write (stdout,'(2x,a,t35,": ",a15)') 'Vibrational entropy model',ctmp
  write (stdout,'(2x,a,t35,": ",f10.2,a)') 'Imaginary freq. threshold',env%thermo%ithr,' cm^-1'
  write (stdout,'(2x,a,t35,": ",f10.4)') 'Frequency scaling factor',env%thermo%fscal
  write (stdout,'(2x,a,t35,": ",f10.2,a)') 'Rot/vib interpolation threshold',env%thermo%sthr,' cm^-1'
  write (stdout,'(2x,a,t35,": ",f6.2,a,f6.2,a,f6.2)') &
    & 'T range [K] (start/end/step)', &
    & env%thermo%trange(1),'/',env%thermo%trange(2),'/',env%thermo%trange(3)
  write (stdout,'(2x,a,t35,": ",i10)') 'Number of temperature points             : ',env%thermo%ntemps
  write (stdout,*)

!========================================================================================!
!> SORTING / CREGEN SETTINGS
!========================================================================================!
  call drawbox(stdout,'Sorting / CREGEN settings',charset=4,padl=2,padr=2)
  write (stdout,*)
  write (stdout,'(2x,a,t35,": ",f10.4,a)') 'Energy window           ',env%ewin,' kcal/mol'
  write (stdout,'(2x,a,t35,": ",f10.4,a)') 'RTHR (RMSD threshold)  ',env%rthr,' Å'
  write (stdout,'(2x,a,t35,": ",f10.4,a)') 'ETHR (energy threshold)',env%ethr,' kcal/mol'
  write (stdout,'(2x,a,t35,": ",f10.2,a)') 'BTHR (rot. threshold)  ',env%bthr2*100.0d0,' %'
  write (stdout,'(2x,a,t35,": ",f10.2)') 'Boltzmann temperature   ',env%tboltz
  write (stdout,'(2x,a,t35,": ",l6)') 'Heavy-atom RMSD only    ',env%heavyrmsd
  write (stdout,'(2x,a,t35,": ",l6)') 'Topology check in CREGEN',env%checktopo
  write (stdout,*)

!========================================================================================!
!> TECHNICAL SETTINGS
!========================================================================================!
  call drawbox(stdout,'Technical settings',charset=4,padl=2,padr=2)
  write (stdout,*)
  call getcwd(dumstr)
  write (stdout,'(2x,a,t25,": ",a)') 'Working directory',trim(dumstr)
  write (stdout,'(2x,a,t25,": ",i0)') 'CPUs / threads',env%threads
  write (stdout,*)

!========================================================================================!
!> CREST BINARY METADATA  (always last)
!========================================================================================!
  call drawbox(stdout,'CREST binary info',charset=4,padl=2,padr=2)
  write (stdout,*)
  call print_crest_metadata()
  write (stdout,*)

!========================================================================================!
  call creststop(status_normal)
end subroutine crest_dry_run
