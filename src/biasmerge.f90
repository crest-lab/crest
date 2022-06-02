!===============================================================================!
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
!===============================================================================!

!=========================================================================================!
! BIASMERGE is used to determine sufficiently different structures
! that can be used as an static bias. k_i and alpha_i will be determined
! based on RMSDs and the read-in energy difference.
!
! The new ensemble must be in the XYZ ensemble format and contain 2 energies
!    in the comment line: first the DFT and then the TB energy
! The old ensemble is also in the XYZ ensemble format, but in its comment line
!    there must be k_i and alpha_i for the structure. This is NOT modified.
!=========================================================================================!
subroutine biasmerge(env)
      use iso_fortran_env, only: wp => real64, sp => real32, dp => int64, output_unit
      use crest_data
      use strucrd
      use iomod
      implicit none
      type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
      type(timer) :: tim2
      integer :: i,j,k,l,p,ich
      integer :: nat,nall
      integer,allocatable :: at(:)
      real(wp),allocatable :: xyz(:,:,:)
      real(wp),allocatable :: kbiasref(:)
      real(wp),allocatable :: alpref(:)
      real(wp),allocatable :: ecdftref(:)
      real(wp),allocatable :: ectbref(:)
      real(wp),allocatable :: rminmaxref(:,:)
      real(wp) :: rdum
      character(len=128),allocatable :: comments(:)
      character(len=:),allocatable :: reffile
      logical :: refgiven

      integer :: nat2,nall2
      integer,allocatable :: at2(:)
      real(wp),allocatable :: xyz2(:,:,:)
      character(len=128),allocatable :: comments2(:)
      real(wp),allocatable :: ecdft(:)
      real(wp),allocatable :: ectb(:)
      real(wp) :: ecdftmin,ectbmin
      integer,allocatable :: take(:)
      integer :: ntaken
      real(wp) :: xx(10)
      
      real(wp) :: quick_rmsd2  !this is a function
      real(wp),allocatable :: rlist(:)
      real(wp),allocatable :: rminmax(:,:)

      real(wp),allocatable :: kbias(:)
      real(wp),allocatable :: alp(:)
      real(wp) :: mad

      character(len=512) :: newcomment,atmp

      logical :: heavy !use only heavy atoms in rmsd?
      real(wp),parameter :: kcal = 627.509541_wp
      real(wp) :: optlevbackup


      integer :: nalltot
      real(wp),allocatable :: rmat(:,:)
      integer :: numrmsdav
      integer :: atry,ktry,ai,ki
      real(wp),allocatable :: mads(:,:)
      real(wp),allocatable :: aval(:),kvals(:)
      real(wp) :: emintry,emintryTB
      integer :: nat3,nall3,mintry
      integer,allocatable :: at3(:)
      real(wp),allocatable :: xyz3(:,:,:)
      real(wp),allocatable :: eread3(:)
      integer :: bestmad(2)
      real(wp) :: alpscal

      optlevbackup = env%optlev
      env%optlev = env%gescoptlev
  !--- determine some settings
      call tim2%init(20)
      heavy = env%cts%gesc_heavy

      numrmsdav = 3
      atry = 4
      allocate(aval(atry), source=1.0_wp)
      aval(2) = 0.75_wp
      aval(3) = 1.25_wp
      aval(4) = 1.50_wp
      alpscal = 1.0_wp
      ktry = env%kshiftnum
      allocate(kvals(ktry), source=0.0_wp)
      if(ktry > 1)then
        kvals(1) = 1.0_wp  
        do i=2,ktry
        kvals(i) = kvals(i-1) + 1.0_wp
        enddo
      else
        kvals = env%kshift  
      endif
           

  !--- first check preexisting files that shall be built upon
      if(allocated(env%cts%rmsdpotfile))then
          reffile = env%cts%rmsdpotfile
      else
          reffile = 'rmsdpot.xyz'
      endif
      nall = 0
      inquire(file=reffile,exist=refgiven)
      if(refgiven)then
          call rdensembleparam(reffile,nat,nall)
          allocate(at(nat), source=0)
          allocate(xyz(3,nat,nall), source=0.0_wp)
          allocate(comments(nall))
          call rdensemble(reffile,nat,nall,at,xyz,comments)
          allocate(ecdftref(nall),ectbref(nall), source=0.0_wp)
          do i=1,nall
            call readl(comments(i),xx,j)
            ecdftref(i) = xx(3)                    
            ectbref(i)  = xx(4)
          enddo
      else
        nall = 0    
      endif    

  !--- then read the NEW ensemble that shall be converted into a bias
      if(allocated(env%biasfile))then
      !================================================================================!    
          call rdensembleparam(env%biasfile,nat2,nall2)
          allocate(at2(nat2), source=0)
          allocate(xyz2(3,nat2,nall2), source=0.0_wp)
          allocate(comments2(nall2))
          call rdensemble(env%biasfile,nat2,nall2,at2,xyz2,comments2)


          nalltot = nall + nall2
          allocate(rmat(nalltot,nalltot), source=0.0_wp)


          !--- check the new ensemble for too similar structures
          allocate(take(nall2),source=1)
          allocate(rminmax(3,nall2), source=0.0_wp)
          allocate(rlist(nall2))
          do i=1,nall2
            rlist=0.0_wp
            if(take(i)==0)cycle
            do j=1,nall2
            if(i==j)cycle
            rlist(j) = quick_rmsd2(nat2,at2,xyz2(:,:,i),xyz2(:,:,j),heavy)
            if(rlist(j) < env%rthr2)then
                take(j) = 0  !RMSD similarity check is done here
                write(*,'(1x,a,i0,a,i0,a,f10.6,a)') 'structure ',j,' discarded (similar to structure ',i,', RMSD = ',rlist(j),')'
            endif   
            rmat(i,j) = rlist(j)
            enddo
            rminmax(2,i) = maxval(rlist,1) !get max
            do j=1,nall2
            if(take(j)==0)rlist(j) = 999.0_wp
            enddo
            rlist(i) = 999.0_wp
            rminmax(1,i) = minval(rlist,1) !get min
          enddo

          !--- with the remaining structures also check the old bias ensemble (if present)
          if(refgiven)then
            allocate(rminmaxref(3,nall), source=-1.0_wp)  
            do i=1,nall
             k=i+nall2
             rlist=999.0_wp
             do j=1,nall
                if(i==j) cycle
                rdum = quick_rmsd2(nat,at,xyz(:,:,i),xyz(:,:,j),heavy)
                l=j+nall2
                rmat(k,l) = rdum
                if(rminmaxref(1,i)<0.0_wp) rminmaxref(1,i) =rdum
                if(rminmaxref(2,i)<0.0_wp) rminmaxref(2,i) =rdum
                if(rdum < rminmaxref(1,i)) rminmaxref(1,i) =rdum
                if(rdum > rminmaxref(2,i)) rminmaxref(2,i) =rdum
             enddo
             rlist=0.0_wp
             do j=1,nall2
               if(take(j)==0)cycle
               rlist(j) = quick_rmsd2(nat,at,xyz(:,:,i),xyz2(:,:,j),heavy)
               if(rlist(j) < env%rthr2)then !RMSD similarity check is done here
                   take(j) = 0  
               else !maybe update minmax values?   
                 rmat(k,j) = rlist(j)  
                 rmat(j,k) = rlist(j)
                 if(rlist(j) > rminmax(2,j)) rminmax(2,j) = rlist(j)
                 if(rlist(j) < rminmax(1,j)) rminmax(1,j) = rlist(j)
                 if(rlist(j) > rminmaxref(2,i)) rminmaxref(2,i) = rlist(j)
                 if(rlist(j) < rminmaxref(1,i)) rminmaxref(1,i) = rlist(j)
               endif    
             enddo
            enddo
          endif

          !--- get an RMSD average (for the x lowest RMSDs) to calculate alpha later
          if(refgiven)then
              do i=1,nall
                k=i+nall2
                l=nalltot
                do j=1,nalltot
                if(rmat(k,j) < env%rthr2)then
                    rmat(k,j) = 999.0_wp
                    l=l-1
                endif
                enddo
                if(l>1)then
                  rminmaxref(3,i) = 0.0_wp    
                  do j=1,min(numrmsdav,l)
                    p = minloc(rmat(k,:),1)
                    rminmaxref(3,i) = rminmaxref(3,i) + rmat(k,p)
                    rmat(k,p) = 999.0_wp
                  enddo
                  rdum = float(min(numrmsdav,l))
                  rminmaxref(3,i) = rminmaxref(3,i) / rdum
                else !fallback for a single bias structure
                 rminmaxref(3,i)=env%rthr2
                endif   
              enddo
          endif
          do i=1,nall2
                k=i
                l=nalltot
                do j=1,nalltot
                if(rmat(k,j) < env%rthr2)then
                    rmat(k,j) = 999.0_wp
                    l=l-1
                endif
                enddo
                if(l>1)then
                  rminmax(3,i) = 0.0_wp    
                  do j=1,min(numrmsdav,l)
                    p = minloc(rmat(k,:),1)
                    !write(*,*)rmat(k,p) 
                    rminmax(3,i) = rminmax(3,i) + rmat(k,p)
                    rmat(k,p) = 999.0_wp
                  enddo
                  rdum = float(min(numrmsdav,l))
                  rminmax(3,i) = rminmax(3,i) / rdum
                else !fallback for a single bias structure
                 rminmax(3,i)=env%rthr2
               endif   
          enddo

          !--- if we have remaining new structures,calculate k_i and alpha
          ntaken = sum(take)
          if(ntaken > 0)then
            allocate(ecdft(ntaken),ectb(ntaken), source=0.0_wp)
            allocate(kbias(ntaken),alp(ntaken), source=0.0_wp)
            k=0
            do i=1,nall2
               if(take(i)==0)cycle
               k=k+1
               call readl(comments2(i),xx,j)
               ecdft(k) = xx(1)
               ectb(k)  = xx(2)
            enddo   
            k = minloc(ecdft,1)
            ecdftmin = ecdft(k)
            ectbmin  = ectb(k)
            if(refgiven)then
                k = minloc(ecdftref,1)
                if(ecdftref(k) < ecdftmin)then
                    ecdftmin = ecdftref(k)
                    ectbmin = ectbref(k)
                endif
                ecdftref = ecdftref - ecdftmin
                ectbref = ectbref - ectbmin
            endif
            ecdft = ecdft - ecdftmin
            ectb  = ectb - ectbmin

       !#########################################################################################!
            if(atry>1 .or. ktry>1)then
       !#########################################################################################!
            call smallhead('OPTIMIZING kshift AND alpha')
            if(refgiven)allocate(kbiasref(nall),alpref(nall), source=0.0_wp)
            allocate(mads(atry,ktry),source=0.0_wp)
            do ai=1,atry
            do ki=1,ktry
            write(*,*)
            write(*,'(1x,a,f6.1,2x,a,f12.4)') 'kshift (kcal/mol):',kvals(ki),'alpha prefactor:',aval(ai)
            j = 0
            if(refgiven)then
              do i=1,nall
                kbiasref(i) = (ecdftref(i) - ectbref(i)) - (kvals(ki)/kcal)
                alpref(i) = aval(ai) / rminmaxref(3,i)**2.0_wp
                j=j+1
              enddo
            endif
            k=0
            do i=1,nall2
            if(take(i)==0)cycle
            k=k+1
            kbias(k) = (ecdft(k) - ectb(k)) - (kvals(ki)/kcal)
            alp(k) = aval(ai) / rminmax(3,i)**2.0_wp
            j=j+1
            enddo

            env%cts%usermsdpot = .true.
            call getcwd(newcomment)
            write(atmp,'(a,i0,a,i0,a)')'rmsdpot-',ai,'-',ki,'.xyz'
            env%cts%rmsdpotfile=trim(newcomment)//'/'//trim(atmp)
            open(newunit=ich,file=trim(atmp))
            emintry=0.0_wp
            mintry = 1
            l=0
            if(refgiven)then
              do i=1,nall
              write(newcomment,'(2F12.6,2x,2F18.8)') kbiasref(i), alpref(i), &
              &    ecdftref(i)+ecdftmin,ectbref(i)+ectbmin
              l=l+1
              if((ecdftref(i)+ecdftmin) < emintry)then
              emintry=ecdftref(i)+ecdftmin
              mintry = l
              endif
              call wrxyz(ich,nat,at,xyz(:,:,i),newcomment)
              enddo
            endif
            if(ntaken > 0)then
              k=0
              do i=1,nall2
              if(take(i)==0)cycle
              k=k+1
              write(newcomment,'(2F12.6,2x,2F18.8)') kbias(k), alp(k), &
              &    ecdft(k)+ecdftmin, ectb(k)+ectbmin
              l=l+1
              if((ecdft(k)+ecdftmin) < emintry)then
                emintry=ecdft(k)+ecdftmin
                mintry = l
              endif
              call wrxyz(ich,nat2,at2,xyz2(:,:,i),trim(newcomment))
              enddo
            endif
            close(ich)
            env%ensemblename=trim(atmp)
            call remove('crest_ensemble.xyz')
            call mdopt(env,tim2)   !optimization with potential
            call rdensembleparam('crest_ensemble.xyz',nat3,nall3)
            allocate(at3(nat3), source=0)
            allocate(xyz3(3,nat3,nall3), source=0.0_wp)
            allocate(eread3(nall3))
            call rdensemble('crest_ensemble.xyz',nat3,nall3,at3,xyz3,eread3)
            mad = 0.0_wp
            emintryTB=eread3(mintry)
            l=0
            if(refgiven)then
              do i=1,nall
              l=l+1
              mad = mad + abs(ecdftref(i) - (eread3(l)-emintryTB))
              enddo
            endif
            if(ntaken > 0)then
              k=0
              do i=1,nall2
              if(take(i)==0)cycle
              k=k+1
              l=l+1
              mad = mad + abs(ecdft(k) - (eread3(l)-emintryTB))
              enddo
            endif
            write(*,*) 'Erel(TB) / kcal/mol :'
            write(*,*) (eread3-emintryTB)*kcal
            mads(ai,ki) = kcal*mad/float(l)
            write(*,*) 'mad (DFT-TB) ', kcal*mad/float(l)
            deallocate(eread3,xyz3,at3)
            enddo
            enddo

            write(*,'(a)') "================"
            bestmad = minloc(mads)
            write(*,'(1x,a,f8.4)')'best MAD =',mads(bestmad(1),bestmad(2))
            write(*,'(1x,a)') 'obtained with:'
            write(*,'(1x,a,f6.1,2x,a,f12.4)') 'kshift (kcal/mol):',kvals(bestmad(2)), &
            &  'alpha prefactor:',aval(bestmad(1))
            env%kshift=kvals(bestmad(2))
            alpscal=aval(bestmad(1))

            call rmrfw('rmsdpot-')
       !#########################################################################################!
            endif
       !#########################################################################################!

            write(*,*)
            write(*,*) 'gESC overwiev (energies in kcal/mol)'
            write(*,'(4x,a10)',advance='no') 'dE(DFT)'
            write(*,'(a10)',advance='no') 'dE(TB)'
            write(*,'(a10)',advance='no') 'kbias'
            write(*,'(a10)',advance='no') 'alpha'
            write(*,'(a10)',advance='no') 'RMSD(min)'
            write(*,'(a10)',advance='no') 'RMSD(max)'
            write(*,'(a10)',advance='no') 'RMSD(avg)'
            write(*,*)
            mad = 0.0_wp
            j = 0
            if(refgiven)then
              if(.not.allocated(kbiasref))then  
              allocate(kbiasref(nall),alpref(nall), source=0.0_wp)    
              endif
              do i=1,nall   
                kbiasref(i) = (ecdftref(i) - ectbref(i)) - (env%kshift/kcal)
                !alpref(i) = 1.0_wp / rminmaxref(1,i)**2.0_wp
                alpref(i) = alpscal / rminmaxref(3,i)**2.0_wp
                write(*,'(i3,a,7F10.4)') i,'*', kcal*ecdftref(i),kcal*ectbref(i),kcal*kbiasref(i),alpref(i), &
                & rminmaxref(1,i),rminmaxref(2,i),rminmaxref(3,i)
                mad = mad + abs(ecdftref(i) - ectbref(i))
                j=j+1
              enddo
            endif
            k=0
            do i=1,nall2
            if(take(i)==0)cycle
            k=k+1
            kbias(k) = (ecdft(k) - ectb(k)) - (env%kshift/kcal)
            !alp(k) = 1.0_wp / rminmax(1,i)**2.0_wp
            alp(k) = alpscal / rminmax(3,i)**2.0_wp
            write(*,'(i3,1x,7F10.4)') i, kcal*ecdft(k),kcal*ectb(k),kcal*kbias(k),alp(k), &
                & rminmax(1,i),rminmax(2,i),rminmax(3,i)
                
            mad = mad + abs(ecdft(k) - ectb(k))
            j=j+1
            enddo
            write(*,*) 'mad (DFT-TB,unoptimized) ', kcal*mad/float(j)
            write(*,'(1x,a,a,a)') 'Structures marked by "*" were read from <',reffile,'>'
          else  
            write(*,*) 'No new structures were added to the bias!'  
            write(*,'(1x,a,a,a)') 'File <',reffile,'> will remain unchanged.'
          endif
      !================================================================================!    
      else
        error stop 'No bias file specified!'  
      endif

  !--- write new bias file    
       if(ntaken > 0)then
       open(newunit=ich,file="rmsdpot.xyz")
       if(refgiven)then
         ecdftref = ecdftref + ecdftmin
         ectbref = ectbref + ectbmin  
         do i=1,nall
         write(newcomment,'(2F12.6,2x,2F18.8)') kbiasref(i), alpref(i),ecdftref(i),ectbref(i)
         call wrxyz(ich,nat,at,xyz(:,:,i),newcomment)
         enddo    
       endif
       if(ntaken > 0)then
         ectb = ectb + ectbmin
         ecdft = ecdft + ecdftmin  
         k=0  
         do i=1,nall2 
         if(take(i)==0)cycle
         k=k+1
         write(newcomment,'(2F12.6,2x,2F18.8)') kbias(k), alp(k), ecdft(k), ectb(k)
         call wrxyz(ich,nat2,at2,xyz2(:,:,i),trim(newcomment))
         enddo
       endif
       close(ich)
       write(*,'(1x,a)')'New bias ensemble file <rmsdpot.xyz> written.'   
       endif

  !--- some preperation if we want to automatically start a conf. sampling:     
       env%cts%usermsdpot = .true.
       call getcwd(newcomment)
       env%cts%rmsdpotfile=trim(newcomment)//'/'//'rmsdpot.xyz' 
       !write(*,*) trim(env%cts%rmsdpotfile)
       env%optlev = optlevbackup

  !--- clear space
      if(allocated(rmat))deallocate(rmat)
      if(allocated(kbias))deallocate(kbias)
      if(allocated(alp))deallocate(alp)
      if(allocated(alpref))deallocate(alpref)
      if(allocated(kbiasref))deallocate(kbiasref)
      if(allocated(ectb))deallocate(ectb)
      if(allocated(ecdft))deallocate(ecdft)
      if(allocated(rlist))deallocate(rlist)
      if(allocated(rminmax))deallocate(rminmax)
      if(allocated(take))deallocate(take)
      if(allocated(at2))deallocate(at2)
      if(allocated(xyz2))deallocate(xyz2)
      if(allocated(comments2))deallocate(comments2)
      if(allocated(ectbref))deallocate(ectbref)
      if(allocated(ecdftref))deallocate(ecdftref)
      if(allocated(at))deallocate(at)
      if(allocated(xyz))deallocate(xyz)
      if(allocated(comments))deallocate(comments)
      return
end subroutine biasmerge      
