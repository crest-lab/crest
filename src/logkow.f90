!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020 Philipp Pracht
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
subroutine logkow(env,tim)
      use iso_fortran_env, wp => real64
      use crest_data
      use strucrd
      implicit none

      type(systemdata) :: env
      !type(options)    :: opt
      type(timer)   :: tim

      type(timer)   :: timer1
      type(timer)   :: timer2

      real(wp) :: t
      real(wp) :: eavw,eavo
      real(wp) :: Kow

      t=298.15d0

      env%crestver=2
      env%gbsa=.true.
!---- first do conformational search in water
      call largehead('Calculation of the ensemble in water')
      env%solv='--gbsa h2o'
      call confscript2i(env,tim) !MTD-GC algo
      call propcalc(conformerfile,10,env,tim) !hessian calc

      call rename('crest_property.xyz','crest_kow_h2o.xyz')
      call rename('crest_rotamers.xyz','crest_h2o.xyz')
      call V2cleanup(.false.)

!---- then do conformational search in octanol
      call largehead('Calculation of the ensemble in 1-octanol')
      env%ensemblename='none selected'  !RESET, IMPORTANT!
      env%solv='--gbsa n-hexane'
      call confscript2i(env,tim) !MTD-GC algo
      call propcalc(conformerfile,10,env,tim) !hessian calc

      call rename('crest_property.xyz','crest_kow_oct.xyz')
      call rename('crest_rotamers.xyz','crest_oct.xyz')
      call V2cleanup(.false.)

!--- calculate logKow
      write(*,'(/,/)')
      call largehead('Calculation of K_ow')
      write(*,*)
      call getelow('crest_kow_h2o.xyz',t,eavw,.true.)
      write(*,*)
      call getelow('crest_kow_oct.xyz',t,eavo,.true.)
      write(*,*)
      write(*,'(1x,a)') 'G(H₂O) = E_gas(H₂O_Geom) + Gsolv(H₂O) + RRHO(H₂O)' 
      write(*,'(1x,a,1x,f20.10)') '       =',eavw
      write(*,'(1x,a)') 'G(1-Oct) = E_gas(1-Oct_Geom) + Gsolv(1-Oct) + RRHO(1-Oct)'
      write(*,'(1x,a,1x,f20.10)') '       =',eavo
      write(*,*)
      write(*,'(1x,a)') 'K_ow = ln( G(H₂O) - G(1-OCt) )'
      kow = eavw - eavo
      kow = log(kow)
      write(*,'(1x,a,f16.8)') 'K_ow =',kow 



!--- turn off other property modes via env%npq
      env%npq = 0

      return
end subroutine logkow
