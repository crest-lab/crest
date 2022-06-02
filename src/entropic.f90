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

!======================================================================================================!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!======================================================================================================!
!===============================================================================!
!  main routine for ensemble entropy calculations
!  with CRE "autocomplete"
!
!  On Input:   pr - printout boolean
!              T  - temperature in Kelvin
!
!  On Output:  S  - ensemble entropy
!              Cp - heat capacity
!  
!  following tasks are performed:
!    1. topology for ref. structure (lowest conformer) is set up
!    2. anmr_nucinfo is read for nuclear equivalencies
!    3. groups of equivalent nuclei are sorted into subgroups
!    4. determine the molecule symmetry for each conformer
!    5. Compare all rotamer groups and check for different core structures
!    6. get rotamer prefactors from the subgroups
!    7. determine rotamers for each conformer
!    8. calculate entropy
! 
!==============================================================================!
subroutine entropic(env,pr,pr2,wrdegen,fname,T,S,Cp)
      use iso_fortran_env, wp => real64, idp => int64
      use crest_data
      use zdata
      implicit none

      type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
      !type(options) :: opt       ! MAIN STORAGE OF BOOLEAN SETTINGS

      logical,intent(in)   :: pr
      logical,intent(in)   :: pr2
      logical,intent(in) :: wrdegen
      character(len=*),intent(in) :: fname
      real(wp),intent(in)  :: T
      real(wp),intent(out) :: S
      real(wp),intent(out) :: Cp

      type(zmolecule) :: zmol
      type(zequal) :: groups
      type(zequal) :: subgroups

      real(wp),allocatable :: rotscal(:) !rotamer equivalencies for each conformer
      integer(idp),allocatable :: introtscal(:) !the same, but as integer
      integer,allocatable :: rotfac(:)   !factor for rotatable groups
      real(wp),allocatable :: symfac(:)  !factor for symmetry
      real(wp),allocatable :: corefac(:) !factor for backbones
      real(wp),allocatable :: enantiofac(:) !factor for backbones
      real(wp),allocatable :: rotconst(:,:)  !rotational constants for each conformer
      real(wp),allocatable :: energies(:)    !energy for each conformer
      character(len=3),allocatable :: symsym(:) !symmetry group symbol for each conformer
      integer,allocatable :: inc(:)
      integer :: i,k,nt
      real(wp) :: dum2
      real(wp) :: tdum,sdum,cpdum,hdum

      integer :: maxrotfac !this is a function
      
      type(zensemble) :: zens

   !--- 0. Set defaults, read ensemble
      S = 0.0_wp
      call creread(fname,zens)
      allocate(rotscal(zens%nconf), source = 1.0_wp)
      allocate(energies(zens%nconf), source = 0.0_wp)
      !--- remap energy (from dimension nrotamers -> nconfomers
      do i=1,zens%nconf
        k=zens%nmat(i,1)  !get the conformer
        energies(i) = zens%eread(k)
      enddo



   !--- 1. topology for reference strucuture
      if(env%wbotopo)then
         env%wbofile='wbo'
      else
         env%wbofile='none given'
      endif
      zens%xyz = zens%xyz/bohr
      call simpletopo(zens%nat,zens%at,zens%xyz,zmol,pr,.true.,env%wbofile)
      zens%xyz = zens%xyz*bohr
      allocate(inc(zmol%nat), source = 0)


   !--- 2. read nuclear equivalencies
      if(pr)then
        write(*,*)
        call smallhead('READING NUCLEAR EQUIVALENCIES')
      endif
      call readequals('anmr_nucinfo',zmol,groups)  
      if(pr)then
        call groups%prsum(6) !--- print summary to screen
        write(*,'(1x,a)') 'Unlisted nuclei (groups) are unique.'
      endif

   !--- 3. distribute groups into subgroups basedon topology
      if(pr)then
        write(*,*)
        call smallhead('ANALYZING EQUIVALENCIES')
      endif
      call distsubgr(zmol,groups,subgroups,inc,pr)

   !--- 4. analyze symmetry group for each conformer to get additional prefactors
      if(pr)then
        write(*,*)
        call smallhead('ANALYZING SYMMETRY')
      endif
      !call analsym(zmol,dum,pr)
      allocate(symfac(zens%nconf),source=1.0_wp)
      allocate(symsym(zens%nconf), source='c1 ')
      call allsym(zens,subgroups,symfac,symsym,pr)
      !!write(*,*) symsym

   !--- 5. Run RMSDs for each group of rotamers. Equivalent atoms must be excluded
      if(pr)then
        write(*,*)
        call smallhead('ANALYZING CORE STRUCTURES BY RMSD')
      endif
      allocate(corefac(zens%nconf),source=1.0_wp)
      allocate(enantiofac(zens%nconf),source=1.0_wp)
      call excludeFromRMSD(zmol,inc)
      call intraconfRMSD(zens,inc,symsym,corefac,enantiofac,.false.)   
      !open(newunit=l,file='corefac.list')
      !do i=1,zens%nconf
      !write(l,*) nint(enantiofac(i)),trim(symsym(i))
      !enddo
      !close(l)
      if(pr) write(*,'(1x,a)') 'done.'

   !--- 5.5 calculate rotational constants for all conformers
     ! if(.not.env%confgo .and. .not.env%fullcre)then
     ! allocate(rotconst(3,zens%nconf))
     ! call nconfRot(zens,rotconst,.false.)
     ! endif


   !--- 6. determine multiplicative factors from rotamers
      i=maxrotfac(subgroups)
      allocate(rotfac(i), source = 0 )      
      call countFactors(subgroups,rotfac,i)
     
   !--- 7. calculating total prefactor per conformer
      if(pr)then
        write(*,*)
        call smallhead('TOTAL ROTAMER NUMBERS')
      endif
      allocate(introtscal(zens%nconf), source=int(1,kind=8))
      call degeneracies(zens%nconf,rotscal,introtscal,i,rotfac,corefac, &
      &    symfac,pr,env%xh3)
      !--- for enso and cregen sorting make a stop here   
      if(env%confgo .and. env%fullcre)then
         call writedegen(zens%nconf,introtscal,'cre_degen2')
         return
      else if(wrdegen)then
         call writedegen(zens%nconf,introtscal,'cre_degen2')
      endif


   !--- 8. calculating entropy  
      call calculateEntropy(zens%nconf,energies,introtscal, &
     &   enantiofac,T,S,Cp,Hdum,pr,pr2)

   !--- 9. temperature dependence
      if( env%properties == -45)then
        if(.not.allocated(env%thermo%temps))then
          call env%thermo%get_temps()
        endif
          nt=env%thermo%ntemps
          if(.not.allocated(env%emtd%soft))then
            allocate(env%emtd%soft(nt), source=0.0d0)
          endif
          if(.not.allocated(env%emtd%cpoft))then
            allocate(env%emtd%cpoft(nt), source=0.0d0)
          endif
          if(.not.allocated(env%emtd%hoft))then
            allocate(env%emtd%hoft(nt), source=0.0d0)
          endif
      endif
      if(env%thermo%ntemps > 0)then
          !tdum = env%thermo%trange(1)
          !tstep = env%thermo%trange(3)
          do i =1,env%thermo%ntemps
             tdum = env%thermo%temps(i)
             sdum=0.0d0
             cpdum=0.0d0
             call calculateEntropy(zens%nconf,energies,introtscal, &
             &   enantiofac,tdum,Sdum,Cpdum,Hdum,.false.,.false.)
             !tdum = tdum + tstep
             if(allocated(env%emtd%soft))then
                env%emtd%soft(i) = sdum
             endif
             if(allocated(env%emtd%cpoft))then
                env%emtd%cpoft(i) = cpdum
             endif
             if(allocated(env%emtd%hoft))then
                env%emtd%hoft(i) = hdum
             endif
          enddo
       !---- temperature dependence printout (for standalone version)
       if(env%properties == -45)then
         write(*,*)
         call smallhead('TEMPERATURE DEPENDENCE FOR ENSEMBLE PROPERTIES')
         write(*,*)
         write(*,'(1x,a)') 'CONFORMATIONAL quantities at given T:'
         write(*,'(a10)',advance='no') "T/K"
         write(*,'(a16)',advance='no') "S"
         write(*,'(a16)',advance='no') "Cp(T)"
         write(*,'(a16)',advance='no') "[H(T)-H(0)]"
         write(*,'(a16)',advance='no') "G"
         write(*,'(a)')
         write(*,'(3x,71("-"))')
         !open(unit=98898,file='.cptmp')
         do i = 1, env%thermo%ntemps
         !   write(98898,'(1x,f10.2,1x,f16.6)') env%thermo%temps(i),env%emtd%cpoft(i)
            write(*,'(f10.2)',advance='no') env%thermo%temps(i)
            write(*,'(f16.6)',advance='no') env%emtd%soft(i)
            write(*,'(f16.6)',advance='no') env%emtd%cpoft(i)
            write(*,'(f16.6)',advance='no') env%emtd%hoft(i)
            dum2 =  -env%emtd%soft(i)*env%thermo%temps(i)/1000.0_wp + env%emtd%hoft(i)
            write(*,'(f16.6)',advance='no') dum2
            write(*,'(a)')
         enddo
         !close(98898)
         write(*,'(3x,71("-"))')
         write(*,'(3x,a)')'S and Cp in cal/mol*K; H and G in kcal/mol'
         write(*,'(3x,a,a)')'G(total) is the ensemble free energy', &
         &    ' and S(total) = S(conf,extrapol.) + δSrrho'
       endif
      endif


      if(allocated(enantiofac))deallocate(enantiofac)
      if(allocated(corefac))   deallocate(corefac)
      if(allocated(introtscal))deallocate(introtscal)
      if(allocated(rotfac))    deallocate(rotfac)
      if(allocated(rotconst))  deallocate(rotconst)
      if(allocated(symsym))    deallocate(symsym)
      if(allocated(symfac))    deallocate(symfac)
      if(allocated(inc))       deallocate(inc)
      if(allocated(energies))  deallocate(energies)
      if(allocated(rotscal))   deallocate(rotscal)
      call zens%deallocate()
      return
end subroutine entropic


subroutine degeneracies(nconf,rotscal,introtscal,n,rotfac,corefac,symfac,pr,xh3) 
     use iso_fortran_env, wp => real64, idp => int64
     implicit none
     integer :: nconf
     real(wp) :: rotscal(nconf)  !total factor to be applied to each conformer
     integer(idp),intent(inout) :: introtscal(nconf)  !total factor to be applied to each conformer as integer
     integer :: n
     integer :: rotfac(n)        !factor arising from rotation (same for all conformers)
     real(wp) :: corefac(nconf)  !factor arising from number of core fragements per conformer
     real(wp) :: symfac(nconf)   !factor from the symmetry of the conformer
     logical :: pr
     integer :: i,k,l
     integer(idp) :: factotal !this is a function
     integer(idp) :: bigint
     real(wp) :: rf,sf
     real(wp) :: xh3
     character(len=8) :: dum
     rotscal(:) = 1.0_wp !reset
     introtscal(:) = 1
     if(pr)then
       write(*,'(1x,a)') "Degeneracies arising from single-bond rotations:"
       write(*,'(1x)',advance='no')
       do i=1,n
         if(rotfac(i).gt.0)then
           write(*,'(a)',advance='no') '+---------'
         endif
       enddo
       write(*,'(a)')'+'
       write(*,'(1x)',advance='no')
       do i=1,n
         if(rotfac(i).gt.0)then
           write(dum,'(i0,a,i0)') i,'^',rotfac(i)
           write(*,'(a,1x,a)',advance='no') '|',dum
         endif
       enddo
       write(*,'(a)')'|'
       write(*,'(1x)',advance='no')
       do i=1,n
         if(rotfac(i).gt.0)then
           write(*,'(a)',advance='no') '+---------'
         endif
       enddo
       write(*,'(a)')'+'
     endif    
     !--- save how many XN3 groups we have:
     if(n>2)then
     XH3 = float(rotfac(3))
     endif
     !--- total rot. factor
     l=int(factotal(rotfac,n), 4)
     rf=float(l)
     if(pr) write(*,'(1x,a,i0)') 'Total factor:  *',l


     if(pr)then
       write(*,'(/,1x,a,i0)')'Min. factor from #core structures: ',nint(minval(corefac,1))
       write(*,'(1x,a,i0)')'Max. factor from #core structures: ',nint(maxval(corefac,1))
       write(*,'(/,1x,a,f8.6,/)')'Min. factor from symmetry: ',minval(symfac,1)
     endif

     do i=1,nconf
       rotscal(i) = rotscal(i) * rf * corefac(i) * symfac(i)
       introtscal(i) = introtscal(i) * l
       introtscal(i) = introtscal(i) * nint(corefac(i))
       if(rotscal(i).lt.1.0d0)then
           rotscal(i)=1.0d0
           introtscal(i) = 1
       else
           sf=1.0/symfac(i)
           introtscal(i) = introtscal(i) / nint(sf)
      endif
       !write(*,'(1x,a,i0,a,i0)') 'rotamer degeneracy for conformer ',i,': ',introtscal(i)
     enddo

     if(pr)then
       write(*,'(/,1x,a)') '----------------------------------' 
       !bigint = nint(minval(rotscal,1))
       k=minloc(rotscal,1)
       bigint = introtscal(k)
       write(*,'(1x,a,i0)')'Min. number of rotamers: ',bigint
       !write(*,*) bigint,minval(rotscal,1)
       !bigint = nint(maxval(rotscal,1))
       k=maxloc(rotscal,1)
       bigint=introtscal(k)
       write(*,'(1x,a,i0,/)')'Max. number of rotamers: ',bigint
       !write(*,*) bigint,maxval(rotscal,1)
     endif

     return
end subroutine degeneracies

subroutine writedegen(nall,introtscal,fname)
    use iso_fortran_env
    implicit none
    integer :: nall
    integer(int64) :: introtscal(nall)
    character(len=*) :: fname
    integer :: ch,i
    open(newunit=ch,file=fname)
    write(ch,'(3x,i0)') nall
    do i=1,nall
        write(ch,'(3x,i0,2x,i0)') i,introtscal(i)
    enddo
    close(ch)
    return
end subroutine writedegen

!======================================================================================================!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!======================================================================================================!

!======================================================!
!  read nuclei equivalences from "anmr_nucinfo" file.
!
!  On Input:  fname   - name of the anmr_nucinfo file
!             zmol    - molecule & topology
!
!  On Output: groups  - object that contains all the
!                       groups of nuclear equivalencies  
!======================================================!
subroutine readequals(fname,zmol,groups)
      use iso_fortran_env, wp => real64
      use filemod
      use zdata
      implicit none

      character(len=*) :: fname
      type(filetype)  :: nucinfo
      type(zmolecule) :: zmol
      type(zequal)    :: groups   !"array" of all groups with some additional data
      type(zgrp)      :: newgrp   !a single group. contains the member atoms 
 
      integer :: nat
      character(len=:),allocatable :: line

      integer :: linc
      integer :: i,j,k
      integer :: ng,nm

      logical :: ex

      inquire(file=fname,exist=ex)
      if(.not.ex)then
        error stop 'no "anmr_nucinfo" file!'
      endif


      call nucinfo%open(fname)
      line=getlarg(nucinfo%line(1),1)
      read(line,*)nat
      call nucinfo%close

      !write(*,'(/,a,i0)') 'nucinfo Nat: ',nat
      if(nat.ne.zmol%nat)then
       error stop 'Mismatch in dimension of anmr_nucinfo and number of atoms!'
      endif
   !--- allocate space for quivalencies data
      call groups%allocate(nat)

   !--- open file 'anmr_nucinfo'
      call nucinfo%open(fname)
      linc=2 !line counter
      do i=1,nat
       !-- cycle atom if we got it already sorted
         if(groups%member(i))then
           linc = linc + 2
           cycle
         endif
       !-- otherwise open a new group
         groups%ng = groups%ng + 1  !increment # of groups
         ng = groups%ng
         line=getlarg(nucinfo%line(linc),2)
         read(line,*)nm             !# members in the group
         call newgrp%allocate(nm)
       !-- fill members of group from next line
         linc = linc + 1
         do j=1,nm
          line=getlarg(nucinfo%line(linc),j)
          read(line,*) k
          newgrp%mem(j) = k   !save atom k into the group
          groups%ord(k) = ng  !save the groupnumber for the atom
         enddo
       !-- sort the atoms in the group and then save group
         call quicksort(nm,newgrp%mem)
         groups%grp(ng) = newgrp
       !-- increment line for next atom
         linc = linc + 1
      enddo
   !--- close 'anmr_nucinfo'
      call nucinfo%close

   !--- count groups with more than 1 atom
      call groups%geteng()   !saved in groups%eng


      return
end subroutine readequals


!=====================================================================================================!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!=====================================================================================================!

!======================================================!
!  distribute atoms from groups into subgroups based
!  on the molecular topology.
!  THIS IS THE ROTAMERGROUP IDENTIFICATION ROUTINE.
!
!  On Input:  zmol    - molecule & topology
!             groups  - groups of nuclear equivalencies 
!
!  On Output: rotgr   - more detailed groups of
!                       nuclear equals. (subgroups)
!
!======================================================!
subroutine distsubgr(zmol,groups,rotgr,exclude,pr)
      use iso_fortran_env, wp => real64
      use zdata
      implicit none

      type(zmolecule) :: zmol
      type(zequal)    :: groups   !-- (super-)groups
      type(zequal)    :: rotgr    !-- rotamer sub-groups
      type(zequal)    :: subgr    !-- sub-groups (dummy)
      logical :: pr

      type(zgrp)      :: grp    !a single group
      integer :: i,j,k,l
      integer :: atk,atl,cent,ring
      integer :: ng,nsg,nms
      integer :: R

      integer,allocatable :: gref(:)  !-- tracking of groups
      logical,allocatable :: centeredSubgroup(:)
      integer,allocatable :: centers(:)
      logical :: specialgrpring !this is a function
      logical :: arruniqel      !this is a function
      logical :: allneigh       !this is a function
      logical :: rotagroup      !this is a function
      logical :: commonring     !this is a function
      integer :: auxgref        !this is a function
      integer :: auxgetgroup   !this is a function
      
      integer,allocatable :: alist(:)
      logical,allocatable :: rings(:)
      logical,allocatable :: rings2(:)
      integer,allocatable :: grouptrac(:)

      logical :: dumbool
      integer :: dumarr(1)

      integer :: exclude(zmol%nat)
      integer,allocatable :: dumpath(:)


      call subgr%allocate(zmol%nat)

      allocate(dumpath(zmol%nat), source = 0)
      allocate(centeredSubgroup(zmol%nat), source =.false.)
      allocate(gref(zmol%nat), source = 0)
      allocate(centers(zmol%nat), source = 0)
      if(zmol%nri .gt. 0)then
        allocate( rings(zmol%nri), source = .false.) !track on which ring we've taken a look so far
        allocate( alist(3), source = 0)
      endif
      allocate( rings2(zmol%nat), source = .false.) !tracking if atoms belong to a group of ring-only-atoms
      


      !write(*,*) size(gref)
      !stop
 !first we need to identify all 'simple' rotameric groups.
 !these are the ones that have exactly one topological
 !neighbour in common. E.g. in a methyl group all the hydrogens
 !build a group and have the same carbon atom as a neighbour.
 !In this loop equivalent atoms that cannot rotate are
 !still INCLUDED, e.g., ...-(CH2)-...
 !Unique atoms (i.e., groups with #members=1) are discarded
   !--- loop over existing (super-)groups from anmr_nucinfo
      ng=groups%ng
      SUPERGROUP : do i=1,ng
         if(groups%grp(i)%nm .lt. 2) cycle SUPERGROUP
         grp = groups%grp(i)
   !--- before we analyze anything else, lets check for freely rotatable rings, e.g. Cp
         if(zmol%nri .gt. 0)then !but we only need to to it if we have rings at all.
           if(specialgrpring(zmol,grp,gref,size(gref)))then !gref should be updated in the function
        !--- if the gorup corresponds to a ring we assign it to a subgroup
             do k=1,grp%nm
                atk = grp%mem(k)
                R = gref(atk)
                call subgr%grp(R)%append(atk)
                exclude(atk) = 1                 !exclude these atoms from the RMSD later
                rings2(atk) = .true.  !indicate that this atom is special, w.r.t being ring-rotateable
                dumarr(1) = atk
                dumbool = commonring(zmol,1,dumarr(1),ring)
                dumpath = 0
                call ringside(zmol,atk,ring,dumpath)
                call atmsTobool(zmol%nat,dumpath,rings2)
                call atmsToboolint(zmol%nat,dumpath,exclude)
             enddo
             cycle SUPERGROUP
           endif
         endif
   !--- loop over all member pairs of the group
         do k=1,grp%nm
            atk = grp%mem(k) 
            do l=k+1,grp%nm
              atl = grp%mem(l)
              !write(*,*) atk, atl              
   !--- are the two members of the group bound to a common neighbour?
              if(arruniqel(zmol%zat(atk)%nei,zmol%zat(atk)%ngh, &
      &          zmol%zat(atl)%nei,zmol%zat(atl)%ngh,cent)) then
                 !write(*,*) 'bonded via:',cent
                 R = auxgref(atk,atl,gref,size(gref))
   !--- append these atoms to a new subroup of the initial group i
                 centeredSubgroup(R) = .true.
                 centers(R) = cent
                 call subgr%grp(R)%append(atk)
                 call subgr%grp(R)%append(atl)
                 !write(*,*) subgr%grp(cent)%mem
                 exclude(atk) = 1           !TODO exclude these atoms from the RMSD later
                 exclude(atl) = 1           !TODO exclude these atoms from the RMSD later
              endif
            enddo   
         enddo
      enddo SUPERGROUP
  !It can happen that there are some new subgroups with 
  !a degeneracy of =2, e.g. due to some symmetry.
  !Some of them have to be sorted out because they are not rotateable
  !and will not concern us any further.
  !special case handling of groups with 2, especially important for rings
  !e.g. phenyl groups
   !--- loop over existing (sub-)groups from previous loop
      nsg = maxval(gref,1)
      do i=1,nsg
         if(subgr%grp(i)%nm .ne. 2) cycle ! only look at the groups with 2 atoms
         grp = subgr%grp(i)
   !--- only two  members in the group
         atk = grp%mem(1)
         atl = grp%mem(2)
   !--- are the two members of the group bound to a common neighbour?
         UNIQUE : if(arruniqel(zmol%zat(atk)%nei,zmol%zat(atk)%ngh, &
      &          zmol%zat(atl)%nei,zmol%zat(atl)%ngh,cent)) then
   !--- if both atoms are bound to the same neighbour, select how to proceed. 
   !--- they can already be taken out from the RMSD for later, e.g., the atoms will be CH2 hydrogens
                 exclude(atk) = 1           !TODO exclude these atoms from the RMSD later
                 exclude(atl) = 1           !TODO exclude these atoms from the RMSD later
   !--- First check if the group can rotate
          ROTOR : if(rotagroup(zmol,cent,grp))then
   !--- check if the two atoms and their neighbour are inside the same ring, only these are interesting here
             if(zmol%nri.gt.0)then
               alist(1) = cent
               alist(2) = atk
               alist(3) = atl
               COMRING : if(commonring(zmol,3,alist,ring))then  !if commenring == true, ring=the corresponding ring
                  !write(*,*) 'ring',ring,atk,atl,'center',cent
                  rings2(atk) = .true.
                  rings2(atl) = .true.
            !-- atoms attached at this ring have to be excluded from RMSD, so we got to identify them
                  dumpath = 0
                  call ringside(zmol,atk,ring,dumpath)
                  call atmsTobool(zmol%nat,dumpath,rings2)
                  call atmsToboolint(zmol%nat,dumpath,exclude)
                  dumpath = 0
                  call ringside(zmol,atl,ring,dumpath)
                  call atmsTobool(zmol%nat,dumpath,rings2)
                  call atmsToboolint(zmol%nat,dumpath,exclude)
            !--- continue with the algo
                  R = auxgref(atk,atl,gref,size(gref))
            !--- ensure that each ring is taken only once ...
                  if(.not.rings(ring))then
                    rings(ring) = .true. 
   !--- append these atoms to a new subroup of the initial group i
                    centeredSubgroup(R) = .true.
                    centers(R) = cent
                    call subgr%grp(R)%deallocate
                    call subgr%grp(R)%append(atk)
                    call subgr%grp(R)%append(atl)
            !--- ... otherwise we can deallocate the subgroup
                  else
                    centeredSubgroup(R) = .false.
                    call subgr%grp(R)%deallocate
                  endif
                endif COMRING
             endif
   !--- if the group is not a rotor, empty its space
          else
             centeredSubgroup(i) = .false.
             call subgr%grp(i)%deallocate
          endif ROTOR

         endif UNIQUE
      enddo
  !At this point we should have all free rotor groups and bound rings
  !that can "generate" rotamers, e.g., methyl groups, phenyl groups etc.
  !What is missing are rotamers from free ring rotations with multiple equivalencies
  !e.g. cyclohexane
      !PRELIMINARY VERSION, EXCLUDE THESE CASES
       do i=1,zmol%nat
        if(centeredSubgroup(i))then
           j=centers(i)
           if(.not.allneigh(zmol,j,subgr%grp(i)%nm,subgr%grp(i)%mem))then
             call subgr%grp(i)%deallocate
           endif
        endif
      enddo
  !We might have taken out too many atoms from the RMSD-backbone comparison.
  !Therefore we have to re-evaluate this data at this point, using the
  !information of the sorted groups and the environment (neighbours)
  !of the atom. Iteration is over the RMSD-exclusion-array "exclude"
      allocate(grouptrac(ng), source = 0 )
      CROSSCHECK : do i=1,zmol%nat
      !--- check only for all excluded
         if(exclude(i).eq.1)then
      !--- all "terminal" atoms, e.g., H,F,Br,Cl etc. should be ok and can be cycled through
            if(zmol%zat(i)%nei .eq. 1) cycle CROSSCHECK
      !--- for now, also ignore freely-rotatable ring atoms, or chains attached to rings
            if(rings2(i)) cycle CROSSCHECK
      !--- for all other atoms there is a condition: 
      !    They can only be excluded if they have an even number of attached groups
            !write(*,*) 'atom:',i
            grouptrac(:) = 0
            do j=1,zmol%zat(i)%nei
               k = zmol%zat(i)%ngh(j)
               l = auxgetgroup(k,groups)
               !write(*,*) 'neighbour:',k,'group:',l
               grouptrac(l) = 1
            enddo
            if(sum(grouptrac(:)).gt.2) exclude(i) = 0
         endif
      enddo CROSSCHECK
      deallocate(grouptrac)



      !nsg = count(centeredSubgroup,1)
      nsg = 0
      !write(*,'(1x,''GROUPS'')')
      do i=1,zmol%nat
        if(allocated(subgr%grp(i)%mem))then
        !  do j=1,subgr%grp(i)%nm
        !   write(*,'(1x,i0)',advance='no') subgr%grp(i)%mem(j)
        !  enddo
        !  write(*,*)
         !write(*,*) subgr%grp(i)%prgrp()
         nsg = nsg + 1
        endif
      enddo
      if(pr)then
      write(*,'(1x,a,i0)') 'Number of initial groups:            ',ng
      write(*,'(1x,a,i0)') 'Number of rotamer-equivalent groups: ',nsg
      endif

      subgr%ng = nsg

     deallocate(rings2)
     deallocate(centers)
     deallocate(gref)
     deallocate(centeredSubgroup)
     deallocate(dumpath)
   !--- finally transfer data from the dummy subgroup to the rotgr object
     call rotgr%allocate(nsg)
     rotgr%ng = nsg
     nms=0
     do i=1,zmol%nat
        if(allocated(subgr%grp(i)%mem))then
         nms = nms+1
         rotgr%grp(nms) = subgr%grp(i)
        endif
      enddo

     return
end subroutine distsubgr


integer function auxgetgroup(atm,grps)
      use zdata
      implicit none
      type(zequal) :: grps
      integer :: atm
      integer :: i
      auxgetgroup = 0
      do i=1,grps%ng
        if(any(grps%grp(i)%mem(:) .eq. atm))then
         auxgetgroup = i
         exit
        endif
      enddo
      return
end function auxgetgroup



!==========================================!
! check if elements i or j of array gref
! have a value different from 0.
! If one (or both of them) do, function
! retrun value is this vaule, otherwise
! it is maxval(gref)+1.
! Furthermore elements i and j are given
! the corresponding value
!==========================================!
integer function auxgref(i,j,gref,d)
      implicit none
      integer :: i,j
      integer :: d
      integer :: gref(d)
      logical :: dum1,dum2
      dum1 = gref(i).ne.0
      dum2 = gref(j).ne.0
      if(dum1)then
        auxgref = gref(i)
      endif
      if(dum2)then
        auxgref = gref(j)
      endif
      if(.not.dum1 .and. .not.dum2)then
        auxgref = maxval(gref,1) + 1
      endif
      gref(i) = auxgref
      gref(j) = auxgref
      return
end function auxgref

integer function grefnewindex(gref,d)
     implicit none
     integer :: d
     integer :: gref(d)
     grefnewindex = maxval(gref,1) + 1
     return
end function grefnewindex

integer function maxrotfac(grps)
     use zdata
     implicit none
     type(zequal) :: grps
     integer :: i
     maxrotfac=1
     do i=1,grps%ng
        if(grps%grp(i)%nm.gt.maxrotfac)then
          maxrotfac=grps%grp(i)%nm
        endif
     enddo
     return
end function maxrotfac

subroutine countFactors(grps,rotfac,n)
     use iso_fortran_env, wp=>real64
     use zdata
     implicit none
     type(zequal) :: grps
     integer :: n
     integer :: rotfac(n)
     integer :: i,k
     do i=1,grps%ng
        k=grps%grp(i)%nm
        rotfac(k) = rotfac(k) + 1
     enddo
     return
end subroutine countFactors

function facTotal(rotfac,n)
     use iso_fortran_env, idp => int64
     implicit none
     integer(idp) :: factotal
     integer :: n
     integer :: rotfac(n)
     integer :: i
     facTotal = 1
     do i=1,n
        !write(*,*) i,rotfac(i)
        facTotal = facTotal * (i**rotfac(i))
     enddo
     return
end function facTotal

!======================================================!
!  Analyze a subgroup bound to a common neighbour (i)
!  and estimate if it can freely rotate, i.e., if it can
!  spawn a  rotamer.
!  To do this, there must be one rotatable bond,
!  responsible for the rotameric equivalence of the
!  grouped atoms. For this number it is checked
!  with the number of neighbours of (i).
!  A special case is where #groupmembers=#neighbours(i)
!  This arises from symmerty, e.g. CH4 or neopentane
!
!  On Input: zmol  - molecule & topology
!            i     - atom indicator
!            grp   - group of atoms bound to i
!
!  On Output: function value
!======================================================!
function rotagroup(zmol,i,grp)
     use iso_fortran_env, wp => real64
     use zdata
     implicit none
     type(zmolecule) :: zmol
     integer :: i
     type(zgrp) :: grp
     logical :: rotagroup
     integer :: m,k,km
     rotagroup = .false.
     m = grp%nm
     k = zmol%zat(i)%nei
     km = k - m 
     !if((km.ge.0).and.(km.le.1)) rotagroup = .true.
     if(km.eq.1) rotagroup = .true.
     return
end function rotagroup

!======================================================!
!  Analyze a group w.r.t. rings in the molecule.
!  We want to check if a symmetry equivalent group
!  coincides with a single ring,
!  or if two (or more) rings are entirly included 
!  in the group. The most important example is
!  Ferrocene, where all 10 carbon atoms build a group
!  and we want to seperate them as 2 ring subgroups
!
!  On Input: zmol  - molecule & topology
!            grp   - group of atoms bound to i
!
!  On Output: function value
!======================================================!
logical function  specialgrpring(zmol,grp,gref,d)
     use iso_fortran_env, wp => real64
     use zdata
     implicit none
     type(zmolecule) :: zmol
     type(zgrp) :: grp
     integer :: d
     integer :: gref(d)
     integer :: grefnewindex !this is a function
     type(zring) :: zr
     integer :: i,j
     integer :: m,km,k,at
     integer :: rtot
     logical,allocatable :: mask(:)
     specialgrpring = .false.
   !--- if there are no rings in the system we can return
     if(zmol%nri.lt.1) return

     m = grp%nm
   !first, are all group members in rings?
     allocate( mask(m) )
     do i=1,m
      mask(i) = zmol%zat(grp%mem(i))%ring
     enddo
     km = count(mask,1)
   !-- not all atoms of the group are in rings, so return
   !-- it should be either all or none of the atoms
     if(km.ne.m)then
        deallocate(mask)
        return
     endif
     deallocate(mask)
   !-- check how many different rings are present in the group (and which ones)
     allocate(mask(zmol%nri), source = .false.)
     do j=1,zmol%nri
        zr = zmol%zri(j)
        do i =1,m
           if(any(zr%rlist(:) == grp%mem(i) )) mask(j) = .true.
        enddo
     enddo
     km = count(mask,1)
   !-- count the total number of atoms in all the included rings
     rtot=0
     do j=1,zmol%nri
       if(mask(j))then
         rtot = rtot + zmol%zri(j)%rs
       endif
     enddo

   !--- if #(atoms in the rings) = #(atoms in the group),
   !    we can safely assume a freely rotatable ring (or several equivalent ones, e.g. Cp in ferrocene)
     if(grp%nm == rtot)then
        do j=1,zmol%nri
          if(mask(j))then !make a new reference group for each of the contained rings
            k = grefnewindex(gref,d) 
            do i=1,zmol%zri(j)%rs !set all the ring atoms to this group
              at=zmol%zri(j)%rlist(i)
              if(gref(at).eq.0)then
                  gref(at)=k
              endif
            enddo
          endif
        enddo
        specialgrpring = .true.
     endif 

     deallocate(mask)

     return
end function specialgrpring

!=====================================================================================================!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!=====================================================================================================!
!======================================================!
!  Run RMSDs between all rotamers belonging to the
!  same conformer to identify the underlying number
!  of different "backbone" strucutres.
!  This number is a factor that has to be applied to
!  the rotamer degeneration.
!
!  On Input: zens - object containing geometries
!            inc  - list of atoms included in RMSD (value 1 or 0)
!            pr   - printout .true./.false.
!
!  On Output: fac - #of backbone structures for each
!                   conformer  
!             enantiofac - 1 or 2, if the conf. has an enantiomer 
!======================================================!
subroutine intraconfRMSD(zens,inc,symsym,fac,enantiofac,pr)
      use iso_fortran_env, wp => real64
      use zdata
      use ls_rmsd
      implicit none
      type(zensemble) :: zens
      integer :: inc(zens%nat)
      real(wp) :: fac(zens%nconf)
      real(wp) :: enantiofac(zens%nconf)
      character(len=3) :: symsym(zens%nconf) !symmetry group symbol for each conformer
      logical :: pr
      logical :: enantio

      integer :: i,j,k
      real(wp),allocatable :: gdum(:,:),Udum(:,:),xdum(:), ydum(:)  ! rmsd dummy stuff
      integer :: n,rednat,nall,nrot,nconf

      real(wp),allocatable :: c1(:,:),c1r(:,:),c2(:,:),c2r(:,:)
      real(wp),allocatable :: rmat(:,:),xyz(:,:,:)
      integer :: co,low,up
      !integer :: uniqueCore !this is a function
      integer,allocatable :: groups(:)
      real(wp) :: rmsdenantio
      real(wp),parameter :: rthr = 0.125   !RMSD difference threshold

      character(len=3),allocatable :: symref(:)
      real(wp),allocatable :: symenantio(:)

      nall = zens%nall
      nconf = zens%nconf
      n = zens%nat
      rednat = sum(inc(:))

   !--- we can return if we would only include a single atom in RMSD
      if(rednat.lt.2)then
         fac(:) = 1.0_wp
         return
      endif

      allocate(gdum(3,3),Udum(3,3),xdum(3),ydum(3))
      allocate(c1(3,n),c1r(3,rednat),c2(3,n),c2r(3,rednat))
      allocate(groups(nconf), source = 1)

      allocate(symref(100))
      allocate(symenantio(100), source = 1.0_wp)
      do i=1,100
      symref(i) = ''
      enddo

!---->      
      CONF : do co=1,nconf !loop over conformers
      nrot = zens%nrot(co) 
      low = zens%nmat(co,1)
      up = zens%nmat(co,nrot)
       
      if(nrot.gt.1)then !skip unique conformers
      allocate(rmat(nrot,nrot), source=0.0_wp)
      allocate(xyz(3,n,nrot))
      xyz(1:3,1:n,1:nrot)= zens%xyz(1:3,1:n,low:up)  !extract the right structures from "zens"  

      enantio = .false.
      ROT : do i=1,nrot
         c1(1:3,1:n)=xyz(1:3,1:n,i)
         call cpincluded(n,rednat,c1,c1r,inc)
!$OMP PARALLEL PRIVATE ( j,k,c2,c2r,xdum,ydum,Udum,gdum) &
!$OMP SHARED ( i,c1r,rmat,n,xyz,rednat,inc)
!$OMP DO 
         do j=1,i-1
!           if(heavy)then
!              call heavyrmsd(n,nall,i,j,xyz,at,rmat(k)) ! all heavy + OH
!           else
               c2(1:3,1:n)=xyz(1:3,1:n,j)
               call cpincluded(n,rednat,c2,c2r,inc)
               call rmsd(rednat,c1r,c2r,0,Udum,xdum,ydum,rmat(i,j),.false.,gdum) 
               rmat(j,i) = rmat(i,j)
!              write(*,*) rmat(k),rotdiff(i,j,nall,rot)
!           endif
!---- enantiomer evaluation
            if(rmat(i,j) > rthr)then
              c2(1,1:n) = -xyz(1,1:n,j)    !artifical enantiomer
              c2(2:3,1:n) = xyz(2:3,1:n,j) 
              call cpincluded(n,rednat,c2,c2r,inc)
              call rmsd(rednat,c1r,c2r,0,Udum,xdum,ydum,rmsdenantio,.false.,gdum)
            !$OMP CRITICAL
              if(rmsdenantio < rthr*0.75) enantio = .true.
            !$OMP END CRITICAL  
            endif
         enddo
!$OMP END DO
!$OMP END PARALLEL
      enddo ROT
      deallocate(xyz)
      groups(co) = uniqueCore(rmat,nrot,rthr)
      deallocate(rmat)
      !--have we found any enantiomer for this conformer?
      if(enantio)then
        enantiofac(co) = 2.0d0
      endif
      endif
      k = symmap(symref,symsym(co))
      if(enantiofac(co) > symenantio(k))then
          symenantio(k) = enantiofac(co)
      endif

      if(pr)then
         write(*,'(1x,a,i0,a,i0,a)')'conformer ',co,': ',groups(co),' basic structures'
      endif
      fac(co) = float(groups(co))
!---->      
      enddo CONF

      do co=1,nconf
      k = symmap(symref,symsym(co))
      enantiofac(co) = symenantio(k) 
      enddo
      
      deallocate(symenantio)
      deallocate(symref)
      deallocate(c2r,c2,c1r,c1)
      deallocate(ydum,xdum,Udum,gdum)
      return
contains
integer function uniqueCore(mat,n,thr)
      use iso_fortran_env, wp => real64
      implicit none
      integer :: n
      real(wp) :: mat(n,n)
      real(wp) :: thr
      integer,allocatable :: mp(:)
      integer :: i,j,gr
      allocate(mp(n), source=0)
      gr=0 !group counter
      do i=1,n
         if(mp(i).eq.0)then
            gr = gr+1
            mp(i) = gr
         endif
         do j=1,n
           if(mat(i,j).lt.thr) mp(j)=gr  
        enddo
      enddo      
      !uniqueCore = maxval(mp,1)
      uniqueCore = gr
      deallocate(mp)
      return
end function uniqueCore
integer function symmap(symlist,symflag)
      implicit none
      character(len=3),intent(inout) :: symlist(100)
      character(len=3) :: symflag
      integer :: i
      symmap = 1
      do i=1,100
       if(trim(symlist(i)) == trim(symflag))then
           symmap = i
           exit
       endif
       if(len_trim(symlist(i))<1)then
        symlist(i) = symflag
        symmap = i
        exit
       endif
      enddo
      return
end function symmap
end subroutine intraconfRMSD

!=====================================================================!
! Routines related to tracking array that indicates if atoms
! are to be taken out from RMSD
!=====================================================================!
subroutine excludeFromRMSD(zmol,inc)
      use iso_fortran_env, wp => real64
      use zdata
      implicit none
      type(zmolecule) :: zmol
      integer :: inc(zmol%nat)
      integer :: i      

!>--- invert include array (up to this point it has only the atoms which shall be excluded)
      do i=1,zmol%nat
        if(inc(i).eq.1)then
          inc(i) = 0
        else
          inc(i) = 1
        endif
      enddo

!>--- as a test remove the hydrogen atoms
!      do i=1,zmol%nat
!        if(zmol%at(i).eq.1) inc(i) = 0 
!      enddo

      return
end subroutine excludeFromRMSD

!========================================!
! auxiliary function
!========================================!
subroutine atmsToboolint(n,atms,inc)
      implicit none
      integer :: n
      integer :: atms(n)
      integer :: inc(n)
      integer :: i,j
      do i=1,n
         if(atms(i).eq.0) cycle
         j=atms(i)
         inc(j) = 1
      enddo
      return
end subroutine atmsToboolint

!========================================!
! auxiliary function
!========================================!
subroutine atmsTobool(n,atms,bool)
     implicit none
     integer :: n
     integer :: atms(n)
     logical :: bool(n)
     integer :: i
     do i=1,n
       if(atms(i).ne.0)then
         bool(atms(i)) =.true.
       endif
     enddo
     return
end subroutine atmsTobool


!=====================================================================================================!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!=====================================================================================================!

!======================================================!
!  Analize the symmetry for each conformer, as they are
!  not the same.
!
!  On Input: zens - object containing geometries
!            grps - groups of rotamer-equivalent atoms
!
!  On Output: fac - list of symmetry factors for each
!                   conformer  
!======================================================!
subroutine allsym(zens,grps,fac,symbol,pr)
      use iso_fortran_env, wp => real64
      use crest_data, only: bohr
      use zdata
      implicit none
      type(zensemble) :: zens
      type(zequal) :: grps
      real(wp) :: fac(zens%nconf)
      character(len=3) :: symbol(zens%nconf)
      logical :: pr

      integer :: nat,ich
      integer,allocatable :: at(:)
      real(wp),allocatable :: xyz(:,:)
      real(wp) :: sf
      character(len=4) :: sfsm
      integer :: i,j,l

      if(zens%nat.gt.200)then
         if(pr)then
          write(*,'(1x,a)') 'Skipped since Nat > 200. Assume no symmetry.'
         endif
         fac = 1.0_wp
         return
      endif

      nat = zens%nat
      allocate(xyz(3,nat),at(nat))
      at(1:nat) = zens%at(1:nat)
              

      if(pr)then
        write(*,'(1x,a)') "For large structures/many conformers this might take a moment ..."
        write(*,'(/,1x,a,i0,a)') 'For ',zens%nconf,' conformers the following'
        write(*,'(1x,a)') 'symmetry groups have been found:'
      endif
      open(newunit=ich,file='symmetries')
      do i=1,zens%nconf
         j = zens%nmat(i,1)
         !write(*,*) zens%nmat(i,1)
         do l=1,nat
           xyz(1:3,l) = zens%xyz(1:3,l,j)/bohr ! use this geometry!!!
         enddo
         call analsym_geo(grps,nat,xyz,at,sf,.false.,sfsm)
         fac(i) = sf
         write(ich,'(3x,i0,1x,a10)')i,sfsm
         if(pr)then
            if(trim(sfsm) /= "c1") write(*,'(3x,a,i0,a,a)') 'conformer ',i,' symmetry: ',sfsm
         endif
         symbol(i) = sfsm(1:3)
      enddo
      close(ich)
      if(pr)then
        write(*,'(1x,a)') 'Unlisted conformers have the point group c1.'
      endif


      deallocate(at,xyz)
      return
end subroutine allsym


!======================================================!
!  Analyze the symmetry of a structure and get a factor
!  for the number of rotamers of it
!======================================================!
subroutine analsym(zmol,fac,pr)
      use iso_fortran_env, wp => real64
      use zdata
      implicit none
      type(zmolecule) :: zmol
      real(wp),intent(out) :: fac
      logical :: pr
      integer :: nat
      real(wp),allocatable :: xyz(:,:)
      integer,allocatable  :: at(:)
      integer :: i
      character(len=4) :: sfsym
      real(wp),parameter :: desy = 0.1_wp
      integer,parameter  :: maxat = 200
      fac = 1.0_wp
      nat = zmol%nat
      if(nat .gt. maxat)then
         return
      endif
      allocate(xyz(3,nat),at(nat))
      do i=1,nat
       xyz(:,i) = zmol%zat(i)%cart(:)
      enddo
      at = zmol%at
      call getsymmetry2(.false.,6,nat,at,xyz, desy, maxat, sfsym)

      if(pr)then
          write(*,'(1x,a,4x,a)') 'symmetry:',sfsym(1:3)
      endif

      return
end subroutine analsym

subroutine analsym_geo(grp,nat,xyz,at,fac,pr,sfsm)
      use iso_fortran_env, wp => real64
      use zdata
      implicit none
      type(zequal) :: grp
      real(wp),intent(out) :: fac
      logical :: pr
      integer :: nat
      real(wp) :: xyz(3,nat)
      integer  :: at(nat)
      character(len=4) :: sfsym
      character(len=4),intent(out) :: sfsm
      real(wp),external :: symfactor
      real(wp),parameter :: desy = 0.1_wp
      integer,parameter  :: maxat = 200
      fac = 1.0_wp
      call getsymmetry2(pr,6,nat,at,xyz, desy, maxat, sfsym)

      fac = 1.0_wp / symfactor(grp,sfsym)
      sfsm=sfsym(1:3)
      return
end subroutine analsym_geo

!======================================================!
!  Function that retruns a (HEURISTIC) symmetry factor
!  for a given point group.
!  Rotamer degeneracies are reduced by this factor.
!  Feel free to add other special cases if you can think
!  of any.
!
!  Dieder groups (Dn, Dnd, Dnh) have a factor of n
!  from the principle axis
!  examples: Ethane (D2d), Ferrocene (D5h/D5d)
!  
!  Tetraeder group (Td) has a factor of 4,
!  example: neo-Pentane  
!
!  Special case of C3v: iso-Butane, factor of 3
!  (somewhat ill-defined criteria, but must have
!   a number of rotamer subgroups dividable by 3)
!
!======================================================!
real(wp) function symfactor(grps,sfsym)
      use iso_fortran_env, wp => real64
      use zdata
      type(zequal) :: grps
      character(len=*) :: sfsym !schoenflies group lable
      character(len=1) :: str
      real(wp) :: ax
      integer :: io
      real(wp) :: ngrps,div
      symfactor = 1.0_wp !default
      str = sfsym(1:1)
      select case( str )
       case( "d" )
        read(sfsym(2:2),*,iostat=io) ax
        if(io .eq. 0)then
         if(ax.gt.2)then
         symfactor = ax !factor comes from principal axis
         endif
        endif
       case( "t" )
        symfactor = 4.0_wp
       case( "o" )
        symfactor = 8.0_wp
       case( "c" )
        if(sfsym == "c3v" ) then
         symfactor = 3.0_wp
         ngrps=float(grps%ng - 1)
         div = ngrps - floor(ngrps)
         if(div .gt. 0.0_wp)then
          symfactor = 1.0_wp
         endif
        endif
       case default
        symfactor = 1.0_wp
      end select
      return
end function symfactor
!=====================================================================================================!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!=====================================================================================================!
!======================================================!
!  Calculate rotational constants for all conformers
!  in the zens object
!
!  On Input: zens - object containing geometries
!            pr   - printout .true./.false.
!
!  On Output: rotconst - rotational constants for
!                        all conformers
!
!======================================================!
subroutine nconfRot(zens,rotconst,pr)
      use iso_fortran_env, wp => real64
      use zdata
      use axis_module
      implicit none
      type(zensemble) :: zens
      real(wp) :: rotconst(3,zens%nconf)
      logical :: pr

      integer :: i,k
      integer :: n,nall,nconf

      real(wp),allocatable :: xyz(:,:)
      integer,allocatable :: at(:)

      real(wp) :: abc(3),bdum

      nall = zens%nall
      nconf = zens%nconf
      n = zens%nat

      if(pr)then
        write(*,*)
        call smallhead('ROTATIONAL CONSTANTS (MHz)')
      endif

      allocate(xyz(3,n),at(n))
      at(1:n) = zens%at(1:n)
      do i=1,nconf
        k=zens%nmat(i,1)  !get the conformer
        xyz(1:3,1:n) = zens%xyz(1:3,1:n,k) 
        abc(:) = 0 !reset
        call axis(n,at,xyz,abc(1:3),bdum)
        rotconst(1:3,i) = abc(1:3)
        if(pr)then
          write(*,'(1x,a,i0,a,3f12.2)') 'conf ',i,' rot: ',abc(1:3)
        endif
      enddo

      deallocate(at,xyz)
      return
end subroutine nconfRot


!=========================================================================================!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!=========================================================================================!
!======================================================!
! Read ensemble into memory, and also get 
! 
! On Input: fname  - name of the ensemble file
!
! On Output: zens  - object containing the coordinates
!                    and rotamer degeneracies
!======================================================!
subroutine creread(fname,zens)
      use iso_fortran_env, wp => real64
      use zdata
      use strucrd, only: rdensembleparam,rdensemble
      implicit none
      character(len=*),intent(in) :: fname
      type(zensemble) :: zens
      integer :: nat,nall
      integer :: i,j,k
      integer :: ich,io
      call rdensembleparam(fname,nat,nall)
      call zens%allocate(nat,nall)
      call rdensemble(fname,nat,nall,zens%at,zens%xyz,zens%eread)
      open(newunit=ich,file='cre_degen')
      read(ich,*,iostat=io) zens%nconf
      if(io .ne. 0)then  ! how to handle if "cre_degen" is invalid/not present
        close(ich)
        zens%nconf = nall  
        if(.not.allocated(zens%nrot))then
         allocate(zens%nrot(zens%nconf))
        endif
        zens%nrot(:) = 1 !all structures considered unique, no rotamers
        zens%nrotmax = 1
      else
       if(.not.allocated(zens%nrot))then
         allocate(zens%nrot(zens%nconf))
       endif
       do i=1,zens%nconf
         read(ich,*) j,zens%nrot(i)
       enddo
       close(ich)
       zens%nrotmax = maxval(zens%nrot,1)
      endif
   !--- the following makes only sense if rotamers were already
   !    grouped together in the file
      if(.not.allocated(zens%nmat))then
       allocate(zens%nmat(zens%nconf,zens%nrotmax), source = 0)
      endif
      k=0
      do i=1,zens%nconf     
         do j=1,zens%nrot(i)
           k = k + 1
           zens%nmat(i,j) = k
         enddo
      enddo

      return
end subroutine creread


!======================================================!
! Read the reference structure for bhess in S_avRRHO
! entropy calculation
! 
! On Input: sys - metadataobject
!           fname - file name to be read
!
! On Output: sys - env%emtd%nat/at/xyz will be filled
!                  coordinates are in Angstroem
!======================================================!
subroutine read_bhess_ref(env,fname)
      use iso_fortran_env, wp => real64
      use crest_data
      use strucrd
      implicit none
      type(systemdata) :: env
      character(len=*) :: fname
      type(coord) :: mol
      if(allocated(env%emtd%at)) deallocate(env%emtd%at)
      if(allocated(env%emtd%xyz))deallocate(env%emtd%xyz)
      call mol%open(fname)
      env%emtd%nat = mol%nat
      env%emtd%at = mol%at
      env%emtd%xyz = mol%xyz * bohr
      call mol%deallocate()

      env%emtd%fromfile=trim(fname)
      return
end subroutine read_bhess_ref


