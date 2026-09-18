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

subroutine crest_setup_alkylize(env)
  use crest_parameters
  use crest_data
  use strucrd
  use molbuilder_classify
  implicit none
  type(systemdata),intent(inout) :: env
  type(coord_classify) :: molc
  type(coord) :: mol

  integer :: ii,jj,kk,cc
  integer :: splt(3)

  call env%ref%to(mol)

  call underline("Analyzing Input Structure")

  call setup_classify(mol,molc)
  call functional_group_classify(molc)
  if (molc%nfuncs == 0) then
    write (stdout,'(a)') 'no relevant substructures found'
    return
  else
    write (stdout,'(a)') 'Found the following substructure parts'
    call molc%print_funcgroups(stdout)
  end if

  do ii = 1,molc%nfuncs
    associate (func => molc%funcgroups(ii))
      if (trim(func%name) == 'alkyl') then

        !> only for propane or longer
        if (func%natms > 6) then
          write (stdout,'(a)') 'selected alkyl group for fragment dispatching'
          splt(:) = 0
          splt(1) = func%attached_to
          kk = 1
          do while (kk < 3)
            do jj = 1,molc%nat
              if (molc%Ah(jj,splt(kk)) == 1.and. &
                & .not.any(splt(:) .eq. jj).and. &
                & any(func%ids(:) .eq. jj).and.kk < 3) then
                kk = kk+1
                splt(kk) = jj
              end if
            end do
          end do
          call env%addsplitqueue(splt)
        end if
        write (stdout,'(2x,a,5(1x,i0))') '> shared atoms:',splt(:)
      end if
    end associate
  end do

  write (stdout,*)
end subroutine crest_setup_alkylize

subroutine crest_proxy_nalkane(env,doreturn)
  use crest_parameters
  use crest_data
  use strucrd
  use molbuilder_classify
  use INTERNALS_mod
  implicit none
  type(systemdata),intent(inout) :: env
  logical,intent(out) :: doreturn
  type(coord) :: mol,newmol
  type(coord_classify) :: molc
  integer :: ii,jj
  integer,allocatable  :: na(:),nb(:),nc(:)
  real(wp),allocatable :: zmat(:,:)
  integer :: itmp(4)

  doreturn = .false.

  if (env%alkylize) then
    call env%ref%to(mol)
    call setup_classify(mol,molc)
    call functional_group_classify(molc)

    do ii = 1,molc%nfuncs
      if (molc%funcgroups(ii)%name == 'alkane'.or.  &
      &  (molc%funcgroups(ii)%name == 'alkyl'.and.  &
      &   molc%funcgroups(ii)%natms >= (molc%nat-3)) &
      &  ) then
        write (stdout,'(a)') '> This substructure contains an n-alkane.'
        if (env%alkylizeskip) then
          write (stdout,'(a)') '> SKIPPING sampling and writing LINEAR structure.'
          doreturn = .true.
        else
          write (stdout,'(a)') '> Writing LINEAR structure and sampling independently.' 
        end if

        call molc%get_zmat(.true.)
        call molc%print_zmat(stdout)

        !> ZMAT construction to make the molecule linear
        do jj = 1,molc%nat
          if (molc%ztod(jj) .ne. 0) then
            itmp(1) = molc%at(jj)
            itmp(2) = molc%at(molc%zmap(jj,1))
            itmp(3) = molc%at(molc%zmap(jj,2))
            itmp(4) = molc%at(molc%zmap(jj,3))
            if (all(itmp(:) > 1)) then
              !write(*,*) 'C-C bond:',molc%zmap(jj,1:2)
              molc%zmat(3,jj) = -pi
            end if
          end if
        end do

        call molc%print_zmat(stdout)
        call molc%from_zmat(newmol)
        call newmol%write(conformerfile)

        call env%ref%load(newmol)
        exit
      end if
    end do

  end if

end subroutine crest_proxy_nalkane

