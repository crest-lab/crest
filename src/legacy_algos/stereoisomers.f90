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

!============================================================!
! STEREOISOMERIZATION TOOL                                   !
! first implementation into CREST, Feb.2020                  !
! was previously a standalone binary                         !
!============================================================!
subroutine stereoisomerize(env,tim)
     use iso_fortran_env, wp => real64
     use crest_data
     use zdata
   
     implicit none

     type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
     type(timer)   :: tim

     type(zmolecule) :: zmol    !molecule object for topology storage

     call header_stereo

     call smallheadline('TOPOLOGY SETUP')
     call tim%start(1,'topology setup')
     call identify_stereocenters(env%inputcoords,zmol)
     call tim%stop(1)

     call smallheadline('GENERATION OF STEREOISOMERS')
     call tim%start(2,'generation step')
     call stereoinvert(zmol)
     call tim%stop(2)


     return
end subroutine stereoisomerize
!===================================================================!
! Topology is analyzed and written to "zmol" polymorphic datatype
! The routine is similar to the "simpletopo" routine, but with
! chemoinformatic identifiaction of stereocenters
! based application of (simplified) Cahn-Ingold-Prelog rules
!===================================================================!
subroutine identify_stereocenters(fname,zmol)
     use iso_fortran_env, wp => real64
     use zdata
     use strucrd, only: rdnat,rdcoord
     implicit none
     character(len=*) :: fname
     type(zmolecule) :: zmol
     logical,parameter :: verbose = .true.

     call simpletopo_file(fname,zmol,.true.,.true.,'')

     !---- identify stereocenters
     call isstereo(zmol)

     return
end subroutine identify_stereocenters
!=======================================================================!
!  get the information which atoms are stereo centers
!=======================================================================!
subroutine isstereo(zmol)
      use iso_fortran_env, only : wp => real64
      use zdata
      implicit none
      type(zmolecule) :: zmol

      integer :: i,k
      integer :: nat
      logical :: teststereo
      logical :: ster
      integer,allocatable :: sphere(:)
      logical,allocatable :: stereoabridged(:)

      nat=zmol%nat
      allocate(sphere(nat))

      write(*,*)
      do i=1,nat
        call stereocond1(i,zmol,teststereo,ster)
        if(teststereo)then
          call cip_spherecomp(zmol,i)
          if(zmol%zat(i)%stereo)then
           call RorS(zmol,i)
          endif
        endif
      enddo
      deallocate(sphere)

      !--- print out total number of stereocenters
      allocate(stereoabridged(nat))
      do i=1,nat
        stereoabridged(i)=zmol%zat(i)%stereo
        if(stereoabridged(i))then
         write(*,'(1x,a,a,a,i0,a,a,a)')'Atom ',trim(zmol%zat(i)%el),'(',i,') is a ',trim(zmol%zat(i)%rs),' stereocenter.'
        endif
      enddo
      !--- count total number of stereocenters
      k= count(stereoabridged,1)
      zmol%nstereo = k
      deallocate(stereoabridged)

      return
end subroutine isstereo
!=======================================================================!
!  quick and dirty conditionals for stereoinformation
!=======================================================================!
subroutine stereocond1(i,zmol,teststereo,ster)
      use iso_fortran_env, only : wp => real64
      use zdata
      implicit none
      type(zmolecule) :: zmol
      type(zatom) :: za
      integer,intent(in)  :: i
      logical,intent(out) :: teststereo
      logical,intent(out) :: ster

      integer :: k,l
      real(wp) :: dum

      logical :: con1,con2,con3,con4,con5

      teststereo = .false.
      ster=.false.

      za=zmol%zat(i)

      con1 = za%nei .eq. 4    !for now only with exactly 4 neighbors  

      !1. discard all centers that don't have enough neighbours
      if(con1)then
         teststereo=.true.
      else
         return  !no potential stereocenter       
      endif

      !--- get weight of neighbour chains
         allocate(za%nweig(za%nei), source=0.0_wp)
         do k=1,za%nei
           l=za%ngh(k)
           call getsideweight(zmol,i,l,za%nweig(k))
         enddo
         zmol%zat(i) = za

      !--- check neighbor type
      !2. if all four neighboring atoms are different,
      !   there defenitly is a sterocenter
      con2 = za%ngt(1) .ne. za%ngt(2)
      con3 = za%ngt(2) .ne. za%ngt(3)
      con4 = za%ngt(3) .ne. za%ngt(4)

      if( con2 .and. con3 .and. con4 ) then
         ster = .true.
         ! the priorities are simply determined from the atom numbers
         za%prio = za%ngh
         do k=1,za%nei
           za%prio(k) = -zmol%at(za%prio(k))
         enddo
         do k=1,za%nei
            l=minloc(za%prio(:),1)
            za%prio(l)=k
         enddo
         za%prioset =.true.
         zmol%zat(i)=za

         return
      endif

      !3. We need to know if atom i is in a ring
      con1 = za%ring
      if(con1)then
         !it is, therefore we need to know which of the neighbours are also in the ring
!         allocate(za%nweig(za%nei), source=0.0_wp)
!         do k=1,za%nei
!           l=za%ngh(k)
!           call getsideweight(zmol,i,l,za%nweig(k))
!         enddo
!         zmol%zat(i) = za
      else
         !it is not, therefore try to check via the side-chain weight if we have different
         !side chains. otherwise we will have to do it together with the Cahn-Ingold-Prelog rules
!         allocate(za%nweig(za%nei), source=0.0_wp)
!         do k=1,za%nei
!           l=za%ngh(k)
!           call getsideweight(zmol,i,l,za%nweig(k))
!         enddo
         !write(*,'(1x,a,a,a,i0,a)')'Atom ',trim(zmol%zat(i)%el),'(',i,') weights:'
         !write(*,*) za%nweig
!         zmol%zat(i) = za
         con5 = .true.
         do k=1,(za%nei-1)
            dum = abs( za%nweig(k) - za%nweig(k+1))
            con5 = con5 .and. (dum .gt. 1d-6)
         enddo
         if(con5) ster=.true.
      endif

      return
end subroutine stereocond1
!=======================================================================!
!  get the absolute stereo configuration (only if we have priorities for
!  the respective center from CIP rules)
!=======================================================================!
subroutine RorS(zmol,i)
      use iso_fortran_env, only : wp => real64
      use zdata
      use geo, only: tangle,rodrot
      implicit none
      integer,intent(in) :: i  !central atom
      type(zmolecule) :: zmol
      real(wp),allocatable :: coords(:,:)
      character(len=2),allocatable :: ele(:)
      integer :: k,l,m,n
      real(wp),parameter :: pi = 3.14159265359_wp
      real(wp) :: theta
      real(wp) :: vec(3),uec(3)
      logical :: pr
      pr = .false.
      if(.not.zmol%zat(i)%prioset) return
      k=zmol%zat(i)%nei
      allocate(coords(3,k),ele(k))
      do l=1,k
         do m=1,k
            if(zmol%zat(i)%prio(m) .eq. l)then
               n=zmol%zat(i)%ngh(m)
               coords(:,l)= zmol%zat(n)%cart(:) - zmol%zat(i)%cart(:)
               ele(l)= zmol%zat(n)%el
            endif
         enddo
      enddo
      !rotate the highest prio atom onto z axis (0,0,1)
      !first into into xz plane, therefore get the -x unit vector first
      !and the projection of the neighbour in the xy-plane
      vec=(/ -1.0d0,0.0d0,0.0d0 /)
      uec(1) = coords(1,1)
      uec(2) = coords(2,1)
      uec(3) = 0.0_wp
      theta = tangle(vec,uec)
      !then rotate around z-axis
      vec = (/ 0.0d0,0.0d0,1.0d0 /)
      if(uec(2).lt.0.0_wp) theta= -theta
      do l=1,k
      call rodrot(coords(:,l),vec,theta)
      enddo

      !afterwards angle to the z-axis
      vec = (/ 0.0d0,0.0d0,1.0d0 /)
      uec = coords(:,1)
      theta = tangle(vec,uec)
      !then take this angle and rotate around y
      vec = (/ 0.0d0,1.0d0,0.0d0 /)
      do l=1,k
      call rodrot(coords(:,l),vec,theta)
      enddo

      !as a last step tak the lowest-prio  neighbour and rotate it into the xy-plane
      vec = (/ -1.0d0,0.0d0,0.0d0 /)
      uec(1) = coords(1,k)
      uec(2) = coords(2,k)
      uec(3) = 0.0_wp
      theta = tangle(vec,uec)
      vec = (/ 0.0,0.0,1.0 /)
      if(uec(2).lt.0.0_wp) theta= -theta
      do l=1,k
      call rodrot(coords(:,l),vec,theta)
      enddo

      if(pr)then
       write(*,*)
       do l=1,k
          write(*,'(1x,a,2x,3f15.10,2x,a)') ele(l),coords(:,l)
       enddo
       write(*,*)
      endif

      !determine orientation
      if(coords(2,2) .gt. coords(2,3)) then
          zmol%zat(i)%rs='(R)'
      else
          zmol%zat(i)%rs='(S)'
      endif

      deallocate(coords)
      return
end subroutine RorS

!=======================================================================!
!  Cahn-Ingold-Prelog routines
!=======================================================================!
subroutine cip_spherecomp(zmol,i)
      ! used to compare the spheres for each neighbour of atom i
      use iso_fortran_env, only : wp => real64
      use zdata
      implicit none
      integer,intent(in) :: i  !neighbour i of k
      type(zmolecule) :: zmol
      type(zatom)     :: za
      integer :: k,l,m,n
      integer :: s,nat
      logical :: equi
      integer :: pri
      integer :: p,q

      integer,allocatable :: dumprio(:),dum2(:)
      logical,allocatable :: incr(:)
      integer,allocatable :: sphereA(:)
      integer,allocatable :: sphereB(:)

      za= zmol%zat(i)
      nat=zmol%nat

      allocate(sphereA(nat),sphereB(nat), source=0)

      allocate(dumprio(za%nei),dum2(za%nei), source=0)
      allocate(incr(za%nei), source=.false.)

      PLOOP : do p=1,za%nei
      incr=.false.
      dum2=999
      pri=0
      s=1
      if(p.eq.2)cycle  !the priorities of the first two atoms is
      if(p.eq.1)then   !determined in the first iteration here
         do
           call cip_sphere(zmol,i,1,s,sphereA)
           call cip_sphere(zmol,i,2,s,sphereB)
           m = sum(sphereA)
           n = sum(sphereB)
           call arrcomp(nat,sphereA,nat,sphereB,equi)
           if(.not.equi)then
             call  arrprio(nat,sphereA,sphereB,pri)
             exit  !we got a prio, exit the loop
           endif
           if(m.eq.0 .or. n.eq.0) exit  !if n or m are 0 then there is no s-th sphere
           s=s+1
         enddo
         if(pri.eq.0)then    !we have two neighbours with indistinguishable priorities, exit
            zmol%zat(i)%stereo = .false.
            zmol%zat(i)%prioset = .false.
            exit PLOOP
         else
           if(pri.eq.1)then
              dumprio(1)=1
              dumprio(2)=2
           endif
           if(pri.eq.2)then
              dumprio(1)=2
              dumprio(2)=1
           endif
         endif
      else  !for all other atoms
        QLOOP : do q=1,p-1
           s=1
           pri=0
           do
             call cip_sphere(zmol,i,p,s,sphereA)
             call cip_sphere(zmol,i,q,s,sphereB)
             m = sum(sphereA)
             n = sum(sphereB)
             call arrcomp(nat,sphereA,nat,sphereB,equi)
             if(.not.equi)then
               call  arrprio(nat,sphereA,sphereB,pri)
               exit  !we got a prio, exit the loop
             endif
             if(m.eq.0 .or. n.eq.0) exit  !if n or m are 0 then there is no s-th sphere
             s=s+1
           enddo
           if(pri.eq.0)then !we have two neighbours with indistinguishable priorities, exit
             zmol%zat(i)%stereo = .false.
             zmol%zat(i)%prioset = .false.
             exit PLOOP
           else
             if(pri.eq.1) incr(q)=.true. !p has higher prio than q, and hence q must be incremented
           endif
        enddo QLOOP
        if(any(incr(:)))then
           do l=1,za%nei
              if(incr(l)) dum2(l)=dumprio(l)
           enddo
           k=minval(dum2(:),1)
           dumprio(p)=k
           do l=1,za%nei
             if(incr(l)) dumprio(l)= dumprio(l) + 1
           enddo
        else
          k=maxval(dumprio(:),1)
          dumprio(p)=k+1
        endif
      endif
      enddo PLOOP

      k=1
      do l=1,za%nei
         k=k*dumprio(l)
      enddo
      if(k.ne.0)then
      zmol%zat(i)%stereo = .true.
      zmol%zat(i)%prioset = .true.
      zmol%zat(i)%prio = dumprio
      endif

      deallocate(incr)
      deallocate(dumprio)
      deallocate(sphereB,sphereA)
      return
end subroutine cip_spherecomp
subroutine cip_sphere(zmol,k,i,s,sphere)
      !get the s-th sphere of neighbour i of atom k
      use iso_fortran_env, only : wp => real64
      use zdata
      implicit none
      integer,intent(in) :: k  !central atom
      integer,intent(in) :: i  !neighbour i of k
      integer,intent(in) :: s  !sphere number
      type(zmolecule) :: zmol
      integer,intent(inout) :: sphere(zmol%nat)
      integer,allocatable :: sphereold(:)
      integer :: j,dummy

      allocate(sphereold(zmol%nat), source=0)

      !get the s-th sphere in direction of i
      sphere=0
      dummy=1
      call cip_getsphere(zmol,k,i,dummy,s,sphere,sphereold)

      !sphere contains the atoms, which we have now to convert to
      !integer atom types:
      do j=1,zmol%nat
         if(sphere(j).eq.0) cycle
         sphere(j) = zmol%at(sphere(j))
      enddo
      !use quicksort to sort the sphere in terms of the atom numbers
      sphere = -sphere   !just a dirty trick since quicksort does return a list from lowest-to-highest
                         !but we need the opposite
      call quicksort(zmol%nat,sphere)
      sphere = -sphere

      deallocate(sphereold)
      return
end subroutine cip_sphere

recursive subroutine cip_getsphere(zmol,k,i,s,sref,sphere,sphereold)
      use iso_fortran_env, only : wp => real64
      use zdata
      implicit none
      integer,intent(in) :: k  !central atom
      integer,intent(in) :: i  !neighbour i of k
      integer,intent(inout) :: s  !sphere number
      integer,intent(in) :: sref
      type(zmolecule) :: zmol
      integer,intent(inout) :: sphere(zmol%nat)
      integer,intent(inout) :: sphereold(zmol%nat)
      logical,allocatable :: mask(:)
      integer,allocatable :: sphereref(:)

      integer :: l,y,m,n
      integer :: p,q,r
      integer :: sn

      y=1

      !the innermost sphere is just the neighbor i itself
      if(sref.le.1)then
         sphere=0
         m=zmol%zat(k)%ngh(i)
         sphere(y) = m
         return
      elseif(s.le.1)then
         sphere=0
         m=zmol%zat(k)%ngh(i)
         sphere(y) = m
         sphereold(1) = k
         sn=s+1  !next sphere
         call cip_getsphere(zmol,k,i,sn,sref,sphere,sphereold)
         return
      endif

      allocate(mask(zmol%nat))
      mask = sphere .ne. 0
      m = count(mask,1)
      mask = .false.
      mask = sphereold.ne.0
      r= count(mask,1) +1
      deallocate(mask)

      !update the history list 'sphereold'
      do l=1,m
         q=sphere(l)
         if(any(sphereold(:).eq.q))cycle
         sphereold(r) = q
         r=r+1
      enddo

      allocate(sphereref(zmol%nat))
      sphereref=sphere
      sphere=0
      do l=1,m
         if(zmol%zat(sphereref(l))%nei .eq. 1) cycle
         n=zmol%zat(sphereref(l))%nei
         do p=1,n
            q=zmol%zat(sphereref(l))%ngh(p)
            if(any(sphere(:).eq.q)) cycle    !no doublicates
            if(any(sphereold(:).eq.q)) cycle !no walking back
            sphere(y)=q
            y=y+1
         enddo
      enddo
      deallocate(sphereref)

      if(s.eq.sref)then
         return
      else
         sn=s+1  !next sphere
         call cip_getsphere(zmol,k,i,sn,sref,sphere,sphereold)
      endif

      return
end subroutine cip_getsphere

!=======================================================================!
!  compare two arrays and check which of them has the higher "priority".
!  arrays must be sorted!
!=======================================================================!
subroutine arrprio(n,narr,marr,pri)
      implicit none
      integer :: n
      integer :: narr(*)
      integer :: marr(*)
      integer :: pri
      integer :: i
      pri = 0
      do i=1,n
         if(narr(i) .eq. marr(i)) cycle
         if(narr(i) .gt. marr(i))then
           pri=1
         else
           pri=2
         endif
         exit
      enddo
      return
end subroutine arrprio

!===================================================================!
! after the topology was analyzed and saved to "zmol"
! the generation of 3D structures is done by the following routine
!===================================================================!
subroutine stereoinvert(zmol)
     use iso_fortran_env, wp => real64
     use zdata
     implicit none
    
     type(zmolecule) :: zmol

     integer :: i,j,k,l,m

     logical :: ster
     integer :: nat
     integer :: ich

     real(wp),allocatable :: vector(:)
     integer :: selected(2)
     integer,allocatable :: path1(:),path2(:)

     logical,allocatable :: book(:)
     logical,allocatable :: checkcounter(:)
     integer :: x,xnei,y,ynei,ycount
  
     integer :: nk
    
     nk=zmol%nstereo
     if(nk.ge.1)then
       write(*,'(/,1x,a,i0,a)') 'There are a total of ',nk,' stereocenters in the molecule.'
       write(*,'(1x,a,i0,a)') 'This results in 2ⁿ = ',2**nk,' different stereoisomers.'
     else
       write(*,'(1x,a)') 'There were no stereocenters found in the molecule.'
       return
     endif
     if(zmol%nstereo .ge.14)then
        write(*,'(/,1x,a)') 'WARNING: There are more than 15000 different stereoisomers for this many stereocenters!'
        write(*,'(1x,a)') 'Hence automatic generation was disabled. Please specifiy which stereocenters to consider by hand.'
        return
     endif

     nat=zmol%nat
     !--- allocate additional zmol stuff
     if(.not.allocated(zmol%inverter))then
        allocate( zmol%inverter(nat,nat), source = 0 )
     endif
     if(.not.allocated(zmol%invector))then
        allocate( zmol%invector(3,nat), source = 0 )
     endif
     if(.not.allocated(zmol%stereotrac))then
        allocate( zmol%stereotrac(zmol%nstereo), source = 0 )
     endif

     allocate(vector(3))

      


     !--- loop over all atoms
     write(*,'(1x,a,/)') 'Neighbouring atoms (with everything attached) chosen to be inverted for each stereocenter:'
     m=1
     do i=1,nat
       ster=zmol%zat(i)%stereo  !--- if we arrive at a stereocenter, invert it
       if(ster)then
         zmol%stereotrac(m)=i
         m=m+1

         !--- special case: do we have a stereocenter in a polycyclic ring system (e.g., steroids)
         if(zmol%zat(i)%ring)then
            xnei=zmol%zat(i)%nei
            allocate(checkcounter(xnei), source = .false.)
            allocate(path1(nat), source = 0) 
            xloop : do x=1,xnei
              checkcounter=.false.
              path1=0
              j=zmol%zat(i)%ngh(x)
              l=0
              call recside(zmol,i,j,path1,l)
              yloop: do y=1,xnei
                 ynei=zmol%zat(i)%ngh(y)
                 if(any(path1(:).eq.ynei))then
                    checkcounter(y)=.true.
                 endif
              enddo yloop
              ycount= count(checkcounter,1)
              if(ycount.ge.3)then
                zmol%polycycle = .true.
                zmol%zat(i)%multiring = .true.
                zmol%invector(3,i)=j                                  ! <----
                !--- get the atoms to be rotated in the special polycycle case
                yloop2: do y=1,xnei
                    if(checkcounter(y).eqv..false.)then
                     allocate(path2(nat), source = 0)
                     k=zmol%zat(i)%ngh(y)
                     l=0
                     call recside(zmol,i,k,path2,l)
                     !zmol%inverter(1,i)=i                              ! <----
                     !zmol%inverter(2:nat,i)=path2(1:nat-1)             ! <----
                     call appendarr(zmol%nat,path2,i)
                     zmol%inverter(:,i)=path2(:)
                     deallocate(path2)
                     write(*,'(1x,a,a,a,i0,a,2x,a,1x,a)')'Stereocenter ',trim(zmol%zat(i)%el),'(',i,'):', &
                    & trim(zmol%prsym(k)),'(polycyclic)' 
                     exit yloop2
                    endif
                 enddo yloop2
                 !--- get the atoms for the rotational axis setup in the polycyle case
                 yloop3: do y=1,xnei
                   if(zmol%invector(1,i).eq.0)then
                     if(checkcounter(y).eqv..true.)then
                        k=zmol%zat(i)%ngh(y)
                        if(k.eq.zmol%invector(3,i)) cycle yloop3
                        zmol%invector(1,i)=k                           ! <----
                     endif
                   else if(zmol%invector(2,i).eq.0)then
                     if(checkcounter(y).eqv..true.)then
                        k=zmol%zat(i)%ngh(y)
                        if(k.eq.zmol%invector(3,i)) cycle yloop3
                        zmol%invector(2,i)=k                           ! <----
                        exit yloop3
                     endif
                   endif
                 enddo yloop3
                exit xloop
              endif
            enddo xloop
            deallocate(path1,checkcounter)
         endif

         !--- regular case
         if(.not. zmol%zat(i)%multiring)then
            call whichinvert(zmol,i,selected,vector)  
            allocate(path1(nat),path2(nat), source = 0)         
            j=zmol%zat(i)%ngh(selected(1))
            k=zmol%zat(i)%ngh(selected(2))
            
            write(*,'(1x,a,a,a,i0,a,2x,a,1x,a)')'Stereocenter ',trim(zmol%zat(i)%el),'(',i,'):', &
            & trim(zmol%prsym(j)),trim(zmol%prsym(k))

            l=0
            call recside(zmol,i,j,path1,l)
            l=0
            call recside(zmol,i,k,path2,l)
            call mergearr(nat,path1,nat,path2)
 
            zmol%invector(1,i)=j                                       ! <----
            zmol%invector(2,i)=k                                       ! <----
            zmol%inverter(:,i)=path1(:)                                ! <----

            deallocate(path2,path1)
         endif
       endif
     enddo
     deallocate(vector)



     write(*,'(/,a)') 'Generation process ...'
     call getbmat(zmol,1,0.2_wp)
     !stop
     allocate(book(zmol%nstereo), source = .false. )
     open(newunit=ich,file='stereomers.xyz')
     x = 0
     call recgetisomer(zmol,x,ich,book)
     close(ich)
     deallocate(book)

     return
end subroutine stereoinvert



subroutine whichinvert(zmol,i,these,vec)
     use iso_fortran_env, wp => real64
     use zdata
     implicit none

     type(zmolecule) :: zmol
     type(zatom) :: zat
     integer :: i

     integer,intent(out) :: these(2)  !--- the two neighbours which must be inverted
     real(wp),intent(out) :: vec(3)   !--- vector of the inversion axis

     integer :: sref
     integer :: j,k,l
 
     real(wp),allocatable :: weights(:)
     real(wp),allocatable :: coords(:,:)

     zat = zmol%zat(i)
     
     sref = 0   !--- if sref=2 at some point in the following, we've found the neighbouring chains that we shall invert
     these(:) = 0 !--- non taken at the start
     vec(:) = 0

     do
       !--- in the easiest case we have neighbours that  are terminal atoms
       do j=1,zat%nei
          k=zat%ngh(j)
          if((zmol%zat(k)%nei.eq.1) .and. (sref.lt.2))then
            sref=sref+1
            these(sref) = j
          endif
          if(sref.ge.2)exit
       enddo
       !--- if we found two such neighbours we can exit the identification loop
       if(sref.ge.2)exit
        
        
       !--- otherwise we take the lightest neighbouring chains
       allocate(weights(zat%nei))
       weights=zat%nweig
       do k=1,zat%nei
          if(sref.lt.2)then
            l=minloc(weights,1)
            weights(l)=99999.0_wp
            if(any(these(:).eq.l))cycle  !cycle for already chosen neighbours
            sref=sref+1       
            these(sref)=l
          else
            exit
          endif
       enddo
       deallocate(weights)

       exit !--- exit the identification block
     enddo

     if(sref.lt.2) error stop 'number of atoms to be inverted < 2, something went wrong!'

     !--- get the inversion axis
     allocate(coords(3,3))
     coords(:,1)=zat%cart(:)
     k=zat%ngh(these(1))
     l=zat%ngh(these(2))
     coords(:,2)=zmol%zat(k)%cart(:) - coords(:,1)
     coords(:,3)=zmol%zat(l)%cart(:) - coords(:,1)
 
     vec(:) = coords(:,2) + coords(:,3) 
     deallocate(coords)

     return
end subroutine whichinvert


!==========================================================================================!
!  merge the elements of arr2 into arr1 (or arr1 into arr2, depending on which is larger)
!==========================================================================================!
subroutine mergearr(i,arr1,j,arr2)
     implicit none
     integer :: i,j
     integer :: arr1(i)
     integer :: arr2(j)
     integer :: k,l,m

     if(i.ge.j)then
        do l=1,j
           m = arr2(l)
           if(any(arr1.eq.m))cycle
           do k=1,i
             if(arr1(k).ne.0)cycle
             arr1(k)=m
             exit
           enddo
        enddo
     else
        do l=1,i
           m = arr1(l)
           if(any(arr2.eq.m))cycle
           do k=1,j
             if(arr2(k).ne.0)cycle
             arr2(k)=m
             exit
           enddo
        enddo
     endif     

     return
end subroutine mergearr


subroutine appendarr(i,arr1,atm)
     implicit none
     integer :: i
     integer :: arr1(i)
     integer :: atm
     integer :: k
     if(any(arr1.eq.atm))return
     do k=1,i
        if(arr1(k).ne.0)cycle
        arr1(k)=atm
        exit 
     enddo
     return
end subroutine appendarr

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C  Recursive generation of all stereo isomers.
!C  Should terminate if the "Depth" reaches the number of stereocenters
!C  "Channel" is the output pipe
!C  "Book" is a bookkeeping array for the stereocenters
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
recursive subroutine recgetisomer(zmol,Depth,Channel,Book)
      use iso_fortran_env, wp => real64
      use zdata
      implicit none
      type(zmolecule) :: zmol
      integer :: Depth
      integer :: Channel
      logical :: Book(zmol%nstereo)

      if(Depth.eq.zmol%nstereo)then
       ! if(zmol%polycycle)then
          call buildstruc2(zmol,Book,Channel)
       ! else
       !   call buildstruc(zmol,Book,Channel)
       ! endif
        return
      endif
      Depth = Depth + 1
      !--- proceed with unchanged structure
      call recgetisomer(zmol,Depth,Channel,Book)
      !--- change structure
      Book(Depth) = .true.
      call recgetisomer(zmol,Depth,Channel,Book) 
      Book(Depth) = .false.
      Depth = Depth - 1
   
      return
end subroutine recgetisomer

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C  build the structure
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine buildstruc(zmol,book,ch)
      use iso_fortran_env, wp => real64
      use zdata
      use crest_data, only: bohr
      use geo, only: rodrot,unitv
      implicit none
      type(zmolecule) :: zmol
      logical :: book(zmol%nstereo)
      integer :: ch
      real(wp),allocatable :: xyz(:,:),dum(:)
      integer :: i,j,k     
      real(wp) :: axis(3) 
      integer :: m,p,q,r
 
      real(wp),parameter :: pi = 3.14159265359_wp

      allocate(xyz(3,zmol%nat),dum(3))
      !--- get the original set of coordinates
      do i=1,zmol%nat 
        xyz(:,i)=zmol%zat(i)%cart(:)
      enddo
      !--- loop over all the stereocenters that have been bookkept
      do j=1,zmol%nstereo
         if(book(j))then
           k=zmol%stereotrac(j)
           if(zmol%zat(k)%multiring)then
              r=zmol%invector(3,k)
              dum(:) = xyz(:,r)
           else
              dum(:) = xyz(:,k)
           endif
           !--- shift the molecule to the position of the stereocenter
           do i=1,zmol%nat
              xyz(:,i) = xyz(:,i) - dum(:)
           enddo
           !--- get the rotation axis and normalize its length
           p = zmol%invector(1,k)
           q = zmol%invector(2,k)
           axis(:) = xyz(:,p) + xyz(:,q)
           call unitv(axis)
           !--- rotate everything attached by 180 degree
           do i=1,zmol%nat
              if(zmol%inverter(i,k).eq.0)cycle
              m = zmol%inverter(i,k)
              call rodrot(xyz(1:3,m),axis,pi)
           enddo
           !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
           !C  Preferably there should follow a geometry optimization at
           !C  this position in the code to repair the molecule, but I'll
           !C  skip that for now.
           !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          else
           cycle
         endif
      enddo
      !--- finally write the coordinates to "ch"
      write(ch,'(4x,i0)') zmol%nat
      write(ch,*)
      do i=1,zmol%nat
         write(ch,'(1x,a,1x,f16.8,1x,f16.8,1x,f16.8)') zmol%zat(i)%el,xyz(1:3,i)*bohr
      enddo
      deallocate(dum,xyz)
      return
end subroutine buildstruc

!===============================================================================!
! Quick and Dirty version of the structure generation WITH optimization
!===============================================================================! 
subroutine buildstruc2(zmol,book,ch)
      use iso_fortran_env, wp => real64
      use zdata
      use geo, only: rodrot,unitv
      implicit none
      type(zmolecule) :: zmol
      logical :: book(zmol%nstereo)
      integer :: ch
      real(wp),allocatable :: xyz(:,:),dum(:)
      integer :: i,j,k
      real(wp) :: axis(3)
      integer :: m,p,q,r
      integer,allocatable :: at(:)

      real(wp),parameter :: pi = 3.14159265359_wp
      real(wp),parameter :: bohr=0.529177_wp

      allocate(xyz(3,zmol%nat),dum(3),at(zmol%nat))
      !--- get the original set of coordinates
      do i=1,zmol%nat
        xyz(:,i)=zmol%zat(i)%cart(:)
        at(i)=zmol%zat(i)%atype
      enddo
      !--- loop over all the stereocenters that have been bookkept
      do j=1,zmol%nstereo
         if(book(j))then
           k=zmol%stereotrac(j)
           if(zmol%zat(k)%multiring)then
              r=zmol%invector(3,k)
              dum(:) = xyz(:,r)
           else
              dum(:) = xyz(:,k)
           endif
           !--- shift the molecule to the position of the stereocenter
           do i=1,zmol%nat
              xyz(:,i) = xyz(:,i) - dum(:)
           enddo
           !--- get the rotation axis and normalize its length
           p = zmol%invector(1,k)
           q = zmol%invector(2,k)
           axis(:) = xyz(:,p) + xyz(:,q)
           call unitv(axis)
           !--- rotate everything attached by 180 degree
           do i=1,zmol%nat
              if(zmol%inverter(i,k).eq.0)cycle
              m = zmol%inverter(i,k)
              call rodrot(xyz(1:3,m),axis,pi)
           enddo
           !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
           !C geometry optimization, xyz in-out
           !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
           call optstruc(zmol%nat,at,xyz)

          else
           cycle
         endif
      enddo
      !--- finally write the coordinates to "ch"
      write(ch,'(4x,i0)') zmol%nat
      write(ch,*)
      do i=1,zmol%nat
         write(ch,'(1x,a,1x,f16.8,1x,f16.8,1x,f16.8)') zmol%zat(i)%el,xyz(1:3,i)*bohr
      enddo
      deallocate(at,dum,xyz)
      return
end subroutine buildstruc2
                               

subroutine optstruc(nat,at,xyz)
      use iso_fortran_env, wp => real64
      use iomod
      use strucrd, only: wrc0,rdcoord
      implicit none
      integer :: nat
      integer :: at(nat)
      real(wp) :: xyz(3,nat)
      integer :: r,l,io
      character(len=512) :: thispath
      character(len=:),allocatable :: jobcall

      !call system('rm -r DUMMY 2>/dev/null')
      call execute_command_line('rm -r DUMMY 2>/dev/null', exitstat=io)

      call getcwd(thispath)
      r = makedir('DUMMY')
      l = sylnk(trim(thispath)//'/'//'gfnff_topo','DUMMY'//'/'//'gfnff_topo')
      l = sylnk(trim(thispath)//'/'//'bondlengths','DUMMY'//'/'//'constraints')
      call chdir('DUMMY')
      call wrc0('coord',nat,at,xyz)

      jobcall = 'xtb coord --gff --opt loose --input constraints > xtb.out 2>> xtb.out'
      !call system('xtb coord --gff --opt loose --input constraints > xtb.out 2>> xtb.out')
      call execute_command_line(trim(jobcall), exitstat=io)

      call rename('xtbopt.coord','coord')
      jobcall='xtb coord --gff --opt loose >> xtb.out 2>> xtb.out'
      !call system('xtb coord --gff --opt loose >> xtb.out 2>> xtb.out')
      call execute_command_line(trim(jobcall), exitstat=io)



      call rdcoord('xtbopt.coord',nat,at,xyz)
      call chdir(trim(thispath))

      return
end subroutine optstruc



