!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2023 Philipp Pracht, Stefan Grimme
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
subroutine nciflexi_legacy(env,flexval)
  use crest_parameters 
  use crest_data
  use iomod
  use strucrd
  implicit none
  type(systemdata) :: env
  character(len=80) :: fname
  character(len=512) :: jobcall
  integer :: io,T,Tn
  logical :: ex
  real(wp) :: ehb,edisp
  real(wp) :: flexval
  type(coord) :: mol
  character(len=25) :: chrgflag
  character(len=*),parameter ::  pipe = ' > xtb.out 2>/dev/null'

!>--- get coordinates
  call env%ref%to( mol ) 

!>--- some options
  call remove('energy')
  call remove('charges')
  call remove('xtbrestart')

!>--- setting threads
  call new_ompautoset(env,'auto',1,T,Tn)

!>--- new plain coord file
  fname = 'tmp.coord'
  call mol%write(fname)
  write(chrgflag,'(a,i0)') '--chrg ',mol%chrg

!>--- jobcall
  write (stdout,'(1x,a)',advance='no') 'Calculating NCI flexibility...'
  flush(stdout)
  jobcall = ""
  jobcall = trim(jobcall)//trim(env%ProgName)
  jobcall = trim(jobcall)//" "//trim(fname)//" --sp --gfnff"
  jobcall = trim(jobcall)//" "//trim(chrgflag)
  jobcall = trim(jobcall)//" "//trim(env%solv)
  jobcall = trim(jobcall)//pipe
  call command(trim(jobcall), io)

!>--- read E(disp) and E(HB) from output
  call grepval('xtb.out','HB energy',ex,ehb)
  call grepval('xtb.out','dispersion energy',ex,edisp)

  write (stdout,'(a)') ' done.'

!>--- normalize by number of atoms
  ehb = ehb/env%nat
  edisp = edisp/env%nat

!>--- NCI flexi is determined RELATIVE to a reference molecule (Crambin)
  flexval = 0.5_wp*(1.0_wp-(ehb/(-0.00043374_wp)))
  flexval = flexval+0.5_wp*(1.0_wp-(edisp/(-0.00163029_wp)))

!>--- cleanup
  call remove('xtb.out')
  call remove('energy')
  call remove('charges')
  if (env%chargesfile) then
    call env%wrtCHRG('')
  end if
  call remove('xtbrestart')
  call remove('xtbtopo.mol')

  return
end subroutine nciflexi_legacy

