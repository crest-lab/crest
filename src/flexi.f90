!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Philipp Pracht
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

!=========================================================================================!
! estimate conformational flexibility from geom. and WBO
! the value of flex is normalized i.e. between 0(rigid) and 1(long alkane)
! SG, Jan. 2019 
! PP, Jan. 2019 modified for use within CREST, include only selected Atoms(atomlist+/-)
!=========================================================================================!
subroutine flexi(nat,rednat,includeRMSD,flex,effectivNat)                      
      use iso_fortran_env, only : wp => real64
      use strucrd, only: rdcoord,i2e
      use zdata, only: readwbo
      implicit none

      integer,intent(in) :: nat
      integer,intent(in) :: rednat
      integer,intent(in) :: includeRMSD(nat)

      real(wp) :: flex

      real(wp),allocatable::xyz(:,:),bond(:,:),rcov(:),cn(:),cring(:,:),wbo(:,:),wbofull(:,:)
      integer,allocatable::at(:),map(:),nb(:,:),b(:,:),sring(:),map2(:)

      real(wp) :: thr
      parameter (thr=1.2)  ! CN less than this is considerd as terminating atom (H, F, ...)

      real(wp) :: dx,dy,dz,r,r2,rco,tmp,val,ringf,doublef,branch,av1,effectivNat,hybf
      real(wp) :: av2
      integer :: i,j,k,l,m,n,rn
      logical samering,ex

      flex = 0

      n=nat     !from module "optionscom"
      rn=rednat !Number of atoms selected by atomlist+/-
                !rn=n if atomlist+/- is not used

      !allocate(xyz(3,n),at(n),rcov(94),map(n),cn(n),nb(20,n),b(n,n),sring(n),cring(10,n),wbo(n,n))
      allocate(xyz(3,rn),at(rn),rcov(94),map(rn),cn(rn),nb(20,rn),b(rn,rn))
      allocate(sring(rn),cring(12,rn),wbo(rn,rn),wbofull(n,n),map2(n))

      call setrcov(rcov) 
     
      !call rdcoord('coord',n,xyz,at)
      if(rn.ne.n)then
        call rdcoord_reduced('coord',n,rn,xyz,at,includeRMSD)
      else
        call rdcoord('coord',rn,at,xyz)
      endif
      call ncoord(rn,rcov,at,xyz,cn,400.0d0)

!--- map the new (reduced) coordinate order to the original
      j=0
      map2=0   !every not-included atom should get a zero
      do i=1,n
         if(includeRMSD(i).gt.0)then
            j=j+1
            map2(i)=j
         endif
      enddo

      !---- debug printouts
      ! write(*,*) map2
      ! do i=1,rn
      !   write(*,'(3f12.6,4x,a2)') xyz(1:3,i),i2e(at(i))
      ! enddo
      !----

!--- read the WBO file and map it to the (reduced) coordinates
      wbo = 0 ! if it does not exist, single and double bonds are not distinguished
      wbofull=0
      inquire(file='wbo',exist=ex)
      if(ex)then
      call readwbo('wbo',n,wbofull)
      !write(*,*) wbofull

!--- map wbo to wbofull (only if atoms are excluded from the rmsd)
      if(rn.ne.n)then
        do i=1,n-1
           if(includeRMSD(i).gt.0)then
              k=map2(i)
              do j=i+1,n
                 if(includeRMSD(j).gt.0)then
                    l=map2(j)
                    wbo(k,l)=wbofull(i,j)
                    wbo(l,k)=wbo(k,l)
                 else
                    cycle  !j cycle
                 endif
              enddo
           else
             cycle         !i cycle
           endif
        enddo
      else !rn=n
        wbo=wbofull
      endif !-- if(rn.ne.n) end
      endif !-- if(ex) end
      deallocate(map2,wbofull)

      !write(*,*) wbo

!---- neighbor list and bond array
      n=rn !<---- so i don't have to change the following original code

      map  = 0
      nb   = 0
      b    = 0
      do i=1,n
         if(cn(i).lt.thr) cycle
         do j=1,n  
            if(i.eq.j) cycle
            dx=xyz(1,j)-xyz(1,i)
            dy=xyz(2,j)-xyz(2,i)
            dz=xyz(3,j)-xyz(3,i)
            r2=dx*dx+dy*dy+dz*dz 
            r=sqrt(r2)    
            rco=rcov(at(i))+rcov(at(j))
            if(r.lt.rco) then
               if(cn(j).gt.thr) then
                  b(j,i)=1
                  map(i)=map(i)+1
                  nb(20,i)=nb(20,i)+1
                  nb(nb(20,i),i)=j
               endif
            endif
         enddo
      enddo

!     call prmati(6,b,n,n,'B')

!---- setup ring data for all atoms, i.e., in which ring and how large

      call minringsizes(n,at,xyz,sring)  !<--- new call, no max ring size

      m = 0
      av2 = 0.0d0
      do i=1,n
         do j=1,i-1
            if(b(j,i).gt.0)  m = m + 1 ! count all bonds! new
            if(map(i).eq.1.or.map(j).eq.1) cycle ! no branch on terminating bonds (e.g. Me)
            if(b(j,i).gt.0) then

               hybf  = 1.0
               if(at(i).eq.6.and.cn(i).lt.3.3) hybf=hybf*1./2. ! sp2 C are less flexible
               if(at(j).eq.6.and.cn(j).lt.3.3) hybf=hybf*1./2.

               doublef = 1.0-exp(-4.0*(wbo(j,i)-2.0)**6 )     ! double bond term 

               branch  = 2.0/sqrt(dble(map(i))*dble(map(j)))  ! branching measure

               ringf = 1.0
               k = min(sring(i),sring(j))
               if(k.gt.0) ringf = 0.5*(1.0-exp(-0.06*float(k))) ! a ring is even at infinite size a factor of 2 less flexible, 0.07 is empirical
                                                                ! (adjusted such that c-C20 is converged)

               val  = branch * ringf * doublef * hybf         ! put it together

               av2  = av2 + val**2                            ! quadratic av

!              write(*,'(2i3,2i2,6f12.4)') i,j,sring(i),sring(j),wbo(j,i),branch,doublef,ringf,hybf,val
            endif
         enddo
      enddo

      if(m.gt.0)then
         av2 = sqrt(av2 / dble(m))
      endif
      
      flex = av2
      effectivNat  = av2 * dble(n)
      deallocate(xyz,at,rcov,map,cn,nb,b,sring,cring,wbo)

      call rmrf('*.zmat')

      return
end subroutine flexi         


subroutine nciflexi(env,flexval)
         use iso_fortran_env, only : wp => real64
         use iomod
         use crest_data
         implicit none
         type(systemdata) :: env
         !type(options)    :: opt
         character(len=80) :: fname,pipe,solv
         character(len=512) :: jobcall
         integer :: io
         logical :: ex
         real(wp) :: ehb,edisp
         real(wp) :: flexval

!---- some options
         pipe=' > xtb.out 2>/dev/null'
         call remove('energy')
         call remove('charges')
         call remove('xtbrestart')
!---- setting threads
         if(env%autothreads)then
            call ompautoset(env%threads,7,env%omp,env%MAXRUN,1) !set the global OMP/MKL variables for the xtb jobs
         endif
!---- new plain coord file
         fname='tmpcoord'
         call copy('coord',fname)
         call clear_setblock(fname)
!---- jobcall
         write(*,'(1x,a)',advance='no')'Calculating NCI flexibility...'
         write(jobcall,'(a,1x,a,1x,a,'' --sp '',a,1x,a,a)') &
         &     trim(env%ProgName),trim(fname),'--gfnff',trim(env%solv),trim(pipe)
         call execute_command_line(trim(jobcall), exitstat=io)

!---- read E(disp) and E(HB) from output
         call grepval('xtb.out','HB energy',ex,ehb)
         call grepval('xtb.out','dispersion energy',ex,edisp)
        
         write(*,'(a)') ' done.'

!---- normalize by number of atoms
         ehb = ehb / env%nat
         edisp = edisp / env%nat

         !--- NCI flexi is determined RELATIVE to a reference molecule (Crambin)
         flexval = 0.5_wp * ( 1.0_wp - ehb / -0.00043374 + 1.0_wp - edisp / -0.00163029 )

!---- cleanup
         call remove(fname)
         call remove('xtb.out')
         call remove('energy')
         call remove('charges')
         if(env%chargesfile)then
             call env%wrtCHRG('')
         endif
         call remove('xtbrestart')
         call remove('xtbtopo.mol')

         return
end subroutine nciflexi

subroutine minringsizes(nat,at,xyz,sring)
    use iso_fortran_env, only: wp=>real64
    use zdata
    implicit none
    integer :: nat
    integer :: at(nat)
    real(wp) :: xyz(3,nat)
    integer :: sring(nat)
    character(len=:),allocatable :: dum
    type(zmolecule) :: zmol
    integer :: i,j,k

    call simpletopo(nat,at,xyz,zmol,.false.,.true.,'')
    do i=1,zmol%nat
       sring(i) = 0
       do j=1,zmol%nri
         if(any(zmol%zri(j)%rlist==i))then
           if( sring(i) == 0 )then
              sring(i) = zmol%zri(j)%rs 
           else if(zmol%zri(j)%rs < sring(i))then
              sring(i) = zmol%zri(j)%rs
           endif
         endif
       enddo
    enddo
    call zmol%deallocate()
    return
end subroutine minringsizes
