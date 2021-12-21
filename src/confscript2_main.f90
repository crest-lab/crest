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

!============================================================================================!
! iMTD(RMSD)-GC Algorithm (also Entropy mode and iMTD-sMTD Algo --v4)
! This is the algo for the conformational search.
!============================================================================================!
subroutine confscript2i(env,tim)
     use iso_fortran_env, only : wp => real64
     use crest_data
     use iomod
     use strucrd, only: coord2xyz,xyz2coord
     implicit none
     type(systemdata) :: env
     !type(options) :: opt
     type(timer)   :: tim
     integer :: i,j,k,l,m
     integer :: eit,eit2
     real(wp) :: time
     real(wp) :: ethr,ediff,ewin
     integer :: dum,bref
     character(len=256) :: str,atmp,btmp
     integer :: nallout
     logical :: lower
     logical :: start
     logical :: ex
     logical :: stopiter,fail
     integer :: fcount

     real(wp) :: autokcal
     parameter (autokcal=627.509541d0)

     settingLogs: associate( performMTD => env%performMTD, restartopt => env%restartopt, &
     &         rotamermds => env%rotamermds, performCross => env%performCross,           &
     &         doNMR => env%doNMR)

     settingData: associate( nat => env%nat, cgf => env%cgf, thresholds => env%thresholds)


!--- some defaults
     ewin = env%ewin    !EWIN

!--- prepare the directory
     if(performMTD)then
     call V2cleanup(restartopt)
     endif
!--- if we only have 2 atoms, do not do anything, except copying the coords
     if(nat.le.2)then
         call catchdiatomic(env)
         return
     endif


!---- Get do a single trial MTD to test the settings
      call V2mdlength(env)  !set the MD length according to a flexibility measure
!---- set number of regular MDs
      call adjustnormmd(env)
 
      if(performMTD)then
         call tim%start(1,'test MD')
         call trialmd(env)    !calculate a short 1ps test MTD to check settings
         call tim%stop(1)
      endif

!---- copy the original coord
      call copy('coord','coord.original')
      call coord2xyz('coord','.history.xyz')

!=====================================================================================================!
!=====================================================================================================!

     if(env%performMTD)then
     write(*,*)
     write(*,*)'list of Vbias parameters applied:'
     do m=1,env%nmetadyn
        write(*,'(''$metadyn '',f10.5,f8.3,i5)') env%metadfac(m)/env%rednat,env%metadexp(m)
     enddo
     endif

     env%nreset=0
     start=.true.
   MAINLOOP: do
      if(env%iterativeV2)then
      call printiter
      endif
      if(.not.start)then
         call clean_V2i  !--- clean Dir for new iterations
         env%nreset = env%nreset+1
      else !--at the beginning clean existing backup-ensembles
         call rmrfw('.cre_') 
      endif
!-------- iterative loop over MTDs
      mtdloop: do i=1,env%Maxrestart
         !---- Small Header
         write(*,*)
         write(*,'(''========================================'')')
         if(env%Maxrestart .gt. 1)then
         write(*,'(''            MTD Iteration '',i2,13x)')i
         else
         write(*,'(''            MTD Simulations '',15x)')
         endif
         write(*,'(''========================================'')')
       !---- do the MTDs
         call tim%start(2,'MTD')
         call MetaMD_para_OMP(env)
         call tim%stop(2)

         MULTILEVELSKIP : if(env%multilevelopt)then
       !---- Optimize using confopt (parallel)
         call tim%start(3,'multilevel OPT')
         write(*,*)
         write(*,'(''-----------------------'')')
         write(*,'(''Multilevel Optimization'')')
         write(*,'(''-----------------------'')')
         write(*,*)
         if(env%optlev >= -2.0d0)then
          call multilevel_opt(env,4)
         endif
         call append_INPUT_to('coord',nat,'input')   !include the input structure into the last optimization
         if(.not.env%superquick .and. env%optlev >= 1.0d0)then
            if(.not.env%entropic)then
              call multilevel_opt(env,5)
            else
              call multilevel_opt(env,99)   !--- the last CREGEN is done within this subroutine
            endif
         endif
         call tim%stop(3)
       !---- save the CRE under a backup name
         call checkname_xyz(crefile,atmp,str)        
         call checkname_xyz('.cre',str,btmp)
         call rename(atmp,btmp)
       !---- save cregen output
         call checkname_tmp('cregen',atmp,btmp)
         call rename('cregen.out.tmp',btmp)
        else 
         exit mtdloop   
        endif MULTILEVELSKIP

       !---- in the first cycle just save the energy and cycle
         if(i.eq.1 .and. start)then
         !---- clean the dir
          start=.false.       !--- only in the first cycle of MAINLOOP the MTDs are done at least 2 times
          if(.not.env%readbias .and. env%runver.ne.33)then
          dum = env%nmetadyn -2
          env%nmetadyn = dum !--- only in the first cycle two MTDs are done additionally with extreme Vbias settings
          endif
          call clean_V2i
          env%eprivious=env%elowest !in the first cycle only save the energy
         !---- save the new best conformer
           inquire(file='crest_best.xyz',exist=ex)
           if(ex)then
           call XYZappendto('crest_best.xyz','.history.xyz')
           call xyz2coord('crest_best.xyz','coord')             !new reference coord to start the MTDs with
           endif

          cycle  !--- always do at least 2 cycles
         endif

       !---- check elowest
        call elowcheck(lower,env)
        if(.not.lower)then
          exit mtdloop
        endif
      enddo mtdloop

!--------
      write(*,'(''========================================'')')
      if(env%Maxrestart .gt. 1)then
      write(*,'(''            MTD Iterations done         '')')
      else
      write(*,'(''           MTD Simulations done         '')')
      endif
      write(*,'(''========================================'')')
       write(*,'(1x,''Collecting ensmbles.'')')   
       call collectcre(env)                      !--- collecting all ensembles saved as ".cre_*.xyz"
       if(.not.env%newcregen)then
             call cregen2(env)  !Legacy subroutine
       else
           if(.not.env%entropic .and. env%crestver.ne.22)then 
             call newcregen(env,0)
           else
             call newcregen(env,2)
           endif
       endif
       call checkname_xyz(crefile,atmp,btmp)
       call remaining_in(atmp,ewin,nallout) !--- remaining number of structures
       write(*,*)
!====================================================================!
!---- (Optional) sampling of additional XH positions
      if(env%doOHflip)then
          call XHorient(env,conformerfile)
          inquire(file='oh_ensemble.xyz',exist=ex)
          if(ex)then
            call checkname_xyz(crefile,atmp,btmp)
            call appendto('oh_ensemble.xyz',atmp)       
            call remove('oh_ensemble.xyz')
            if(.not.env%entropic .and. env%crestver.ne.22)then
              call newcregen(env,0)
            else
              call newcregen(env,2)
            endif
            call remaining_in(btmp,ewin,nallout)
            write(*,*)
          endif
      endif
!====================================================================!
!---- Perform additional MDs on the lowest conformers
      if(env%rotamermds)then
          call tim%start(4,'MD ')
          call normalMD_para_OMP(env,env%nrotammds,env%temps)
          if(env%multilevelopt)then
           if(env%optlev >= 1.0d0)then   
           call multilevel_opt(env,6)
           else
           call multilevel_opt(env,99)    
           endif
           call elowcheck(lower,env)
          else
           lower = .false.
          endif
          call tim%stop(4)
          if(lower)then
            call checkname_xyz(crefile,atmp,str)  
            call checkname_xyz('.cre',str,btmp)
            call rename(atmp,btmp)
            if(env%iterativeV2)cycle MAINLOOP
          endif
      endif
!====================================================================!
!---- Genetic crossing
     if(env%performCross)then
        call tim%start(5,'GC')
        call cross2(env)
        call tim%stop(5)
        call confg_chk3(env)       
        call elowcheck(lower,env)
        if(lower)then
          call checkname_xyz(crefile,atmp,str)
          call checkname_xyz('.cre',str,btmp)
          call rename(atmp,btmp)
          if(env%iterativeV2) cycle MAINLOOP
        endif
     endif
!====================================================================!
!---- Entropy mode iterative statically biased MDs
     if(env%entropymd)then
       call mtdatoms('coord',env)
       call tim%start(6,'static MTD')
       call emtdcopy(env,0,stopiter,fail)
       bref=env%emtd%nbias

       ENTROPYITER : do eit=1,env%emtd%iter
         dum = nint(float(env%emtd%nbias) * env%emtd%nbiasgrow)
         !env%emtd%nbias = nint(float(env%emtd%nbias) * env%emtd%nbiasgrow)
         env%emtd%nbias = max(env%emtd%nbias+1,dum)
         fail = .false.
         EFALLBACK : do k=1,env%emtd%maxfallback
         call printiter2(eit)
         call tim%start(6,'static MTD')
         call entropyMD_para_OMP(env)
         call tim%stop(6)
         call emtdcheckempty(env,fail,env%emtd%nbias)
         if(fail)then
             if(k==env%emtd%maxfallback)then
               stopiter=.true.
             else
               cycle EFALLBACK
             endif
         else
            call tim%start(3,'multilevel OPT')
            if(env%optlev >= -1.0d0)then
             call multilevel_opt(env,2)
            endif
            call multilevel_opt(env,99)

            call tim%stop(3)
            !--- if in the entropy mode a lower structure was found 
            !    --> cycle, required for extrapolation
            call elowcheck(lower,env)
            if(lower.and.env%entropic)then
                env%emtd%nbias = bref  !IMPORTANT, reset for restart
                cycle MAINLOOP
            endif
            !--- file handling
            eit2=eit
            call emtdcopy(env,eit2,stopiter,fail)
            env%emtd%iterlast = eit2
         endif
         if(.not.lower .and. fail .and. .not.stopiter)then
           cycle EFALLBACK
         endif
         exit EFALLBACK  !fallback loop is exited on first opportuinity
         enddo EFALLBACK
         if(stopiter)then
             exit ENTROPYITER
         endif
       enddo ENTROPYITER
     endif

!-------
     exit MAINLOOP !--- if this point is reached, i.e., there weren't any further restarts, exit the loop
  enddo MAINLOOP

!=====================================================================================================!
!=====================================================================================================!

      !if(.not.env%entropic .and. .not.(env%crestver == 22))then
      if(.not.env%entropymd)then
      !------ last optimization (with user set optlevel)
          write(*,*)
          write(*,*)
          write(*,'(3x,''================================================'')')
          write(*,'(3x,''|           Final Geometry Optimization        |'')')
          write(*,'(3x,''================================================'')')

          if(doNMR) cgf(3)=.true.                        !--- if NMR equivalencies are requested, turn them on here
          call tim%start(7,'')
          call multilevel_opt(env,99)   !--- the last CREGEN is done within this subroutine
          call tim%stop(7)                                 !--- optlevel is userset
      else
      !------ or just sort the last ensemble for entropy mode    
          !if(doNMR) cgf(3)=.true.
          !call confg_chk3(env) 
      endif

!---- print CREGEN results and clean up Directory a bit
    if(env%crestver .ne. crest_solv) then
      call V2terminating()
    end if

      end associate settingData
      end associate settingLogs

end subroutine confscript2i

!=========================================================================================!
!--- set the total run time according to the mRMSD criterium
!=========================================================================================!
subroutine V2mdlength(env)
      use iso_fortran_env, only : wp => real64
      use crest_data
      implicit none
      type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
      !type(options) :: opt       ! MAIN STORAGE OF BOOLEAN SETTINGS
      integer :: mdlenfactor
      real(wp) :: total,k,minimum,fRMSD,tmtd,lenthr
      real(wp) :: flex,av1,rfac,nciflex,flextot

      settingBool: associate( quick => env%quick, QCG => env%QCG, V2i => env%iterativeV2, NCI => env%NCI)

       minimum=5.0d0                    !at least 5ps per MTD
       lenthr=200.0d0     ! Maximum of 200 ps, longer runs can only be conducted by user input      

       write(*,*)
       write(*,'(''------------------------------------------------'')')
       write(*,'(''Generating MTD length from a flexibility measure'')')
       write(*,'(''------------------------------------------------'')')

       if((env%crestver .ne. crest_solv) .and. .not.NCI)then
         write(*,'(1x,a)',advance='no')'Calculating WBOs...'
         call xtbsp(env,0)       !xtb singlepoint to get WBOs (always GFN0)
         write(*,'(1x,a)') 'done.'
         call flexi(env%nat,env%rednat,env%includeRMSD,flex,av1)
         call nciflexi(env,nciflex) !NCI flexi based on E(HB)/Nat and E(disp)/Nat
         write(*,'(1x,''    covalent flexibility measure :'',f8.3)') flex
         write(*,'(1x,''non-covalent flexibility measure :'',f8.3)') nciflex
         flex = 0.5*flex + 0.5*nciflex*sqrt(flex)   ! the NCI flex is only relevant if the covalent framework is flexible
         if(env%entropic)then ! special case for entropy mode
           av1 = (flex**1.333333)*max(1,env%rednat-8) ! depends more strongly on the flexibility than on size, -8 accounts for small systems having no conf.
           lenthr=3000.0d0                            ! Maximum of 3000 ps, longer runs can only be conducted by user input      
           env%tmtd = 4.50d0*exp(0.165d0*av1)         ! minimum is 5 ps set above
         else                 ! normal case
           av1 = (flex**1.000000)*max(1,env%rednat-8) ! 
           lenthr=500.0d0                             ! Maximum of 500 ps
           env%tmtd = 3.0d0*exp(0.10d0*av1)           ! 
         endif    
       else
         write(*,'(1x,"System flexiblity is set to 1.0 for NCI mode")')  
         flex=1.0d0
         env%tmtd = 0.10*(env%rednat+0.1*env%rednat*env%rednat)
       endif
       env%flexi = flex
       write(*,'(1x,''flexibility measure :'',f8.3)') env%flexi

!---- rfac is used to scale the total MD length according to special runtypes
       total=max(minimum,env%tmtd)
       select case( env%runver ) 
         case( 2,5,6,33 )   ! "-quick","-squick","-mquick"
            rfac = 0.5d0
         case( 3 )     ! "-qcg"
            rfac = 0.25d0
         case( 77 ) 
            rfac = 1.50d0  
         case( 8 ) 
            rfac = 2.0d0
         case default  ! everything else 1=default, 4=NCI
            rfac = 1.0d0
       end select

        if(env%scallen)then
            write(*,'(1x,''t(MTD) based on flexibility :'',  f8.1)') env%tmtd*rfac
            rfac = env%mdlenfac * rfac
            write(*,'(1x,''MTD length is scaled by     :'',  f6.3)') env%mdlenfac
        endif

      if(env%mdtime.le.0.0d0)then     !<-- ONLY use generated MD length if not already set by the user
        if(total.gt.lenthr)then
           total=lenthr
           call mtdwarning(lenthr*rfac)
        endif
        env%mdtime = anint(total) * rfac
      else  
         write(*,'(1x,''t(MTD) / ps set by command line  :'',  f8.1)')env%mdtime
      endif

       write(*,'(1x,''t(MTD) / ps    :'',  f8.1)')env%mdtime       
       write(*,'(1x,''Σ(t(MTD)) / ps :'',  f8.1,'' ('',i0,'' MTDs)'')') &
       & env%mdtime*float(env%nmetadyn),env%nmetadyn
      
      env%metadlist(:)=ceiling(env%mdtime)  !each ps of the MTD a Vbias snapshot is taken
      end associate settingBool
      return
end subroutine V2mdlength

!===========================================================================================!
! Set METADYN default Guiding Force Parameter
! There are different combinations depending on the runtype
!===========================================================================================!
subroutine iV2defaultGF(env)
     use iso_fortran_env, only : wp => real64
     use crest_data
     use filemod
     implicit none
     
     type(systemdata) :: env
     !type(options)    :: opt

     integer  :: ia,ik,na,nk,m,nmtdyn
     real(wp) :: alp,k
     real(wp) :: kstart,kinc
     real(wp) :: alpinc

     type(filetype) :: biasfile
     logical :: ex
     integer :: i,j,io
     character(len=:),allocatable :: atmp

     settingBool: associate( quick => env%quick, QCG => env%QCG, V2i => env%iterativeV2)
!----
     if(.not.env%metadynset)then
        if(env%readbias)then
           inquire(file='mtdbias',exist=ex) 
           if(ex)then
           atmp=''    
           call biasfile%open('mtdbias')
           call biasfile%clearblanks()
           nmtdyn=biasfile%nlines
           call env%allocate(nmtdyn)
           do i=1,nmtdyn
             atmp=getlarg(biasfile%line(i),1)
             read(atmp,*,iostat=io) k
             if(io == 0)then
               env%metadfac(i) = k*env%rednat
             endif
             atmp=getlarg(biasfile%line(i),2)
             read(atmp,*,iostat=io) alp
             if(io == 0)then
               env%metadexp(i) = alp
             endif
           enddo           
           !write(*,*) env%metadfac
           !write(*,*) env%metadexp
          else  
           error stop "no file 'mtdbias'"
          endif          
        else    
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++!    

     select case( env%runver )
!----------
       case( 2,5 ) ! "-quick","-squick"
         na=3
         nk=2
         nmtdyn=na*nk
         alp=1.2d0 ! start value alpha
         kstart=0.002d0 ! start value k
         alpinc = 2.0 ! increment
         kinc =  2.0     ! increment
!----------
       case( 6 ) ! "-mquick""
         na=3
         nk=2
         nmtdyn=na*nk
         alp=1.0d0 ! start value alpha
         kstart=0.002d0 ! start value k
         alpinc = 2.0 ! increment
         kinc =  2.0     ! increment
!----------
       case( 3 ) ! "-qcg"
         na=4
         nk=3
         nmtdyn=na*nk
         alp=1.0d0 ! start value alpha
         kstart=0.00125d0 ! start value k
         alpinc = (3./2.) ! increment
         kinc =  (3./2.)   ! increment
!----------
       case( 4 )  ! "-nci"
         na=3
         nk=2
         nmtdyn=na*nk
         alp=1.0d0 ! start value alpha
         kstart=0.001d0 ! start value k
         alpinc = 2.0 ! increment
         kinc =  2.0     ! increment
!----------
       case( 45  )  ! "-singlerun"
         na=1
         nk=1
         nmtdyn=na*nk
         alp=1.0d0 ! start value alpha
         kstart=0.001d0 ! start value k
         alpinc = 2.0 ! increment
         kinc =  2.0     ! increment
!----------
       case( 33 ) ! "-relax"
         na=1
         nk=3
         nmtdyn=na*nk
         alp=0.8d0 ! start value alpha
         kstart=0.0030 ! start value k
         alpinc = 2.0 ! increment
         kinc =  2.0     ! increment
!----------
       case( 77 ) ! "-compress"
         na=3
         nk=3
         nmtdyn=na*nk
         alp=1.61803 ! start value alpha
         kstart=0.005d0 ! start value k
         alpinc = 1.61803 ! increment
         kinc =  2.0     ! increment
!----------
       case( 111 ) ! "-entropy"
         na=6
         nk=4 
         nmtdyn=(na*nk)
         alp=1.61803 ! start value alpha
         kstart=0.0075d0 ! start value k 
         alpinc = 1.61803 ! increment
         kinc =  2.0     ! increment
!----------
       case default
        if(V2i)then  !for the default iterative mode
        !=======================================================!    
         na=4
         nk=3
         nmtdyn=(na*nk)+2
         call env%allocate(nmtdyn)   !allocate k(Vbias) and α(Vbias)
         alp=1.3 ! start value alpha
         kstart=0.0030 ! start value k 
         alpinc = (5./3.) ! increment
         kinc =  2.0     ! increment
 
         !-- two additional MTDs with extreme values
         env%metadfac(nmtdyn-1)=0.001*env%rednat
         env%metadexp(nmtdyn-1)=0.1

         env%metadfac(nmtdyn)=0.005*env%rednat
         env%metadexp(nmtdyn)=0.8
        !=======================================================! 
        else  !for the non-iterative mode
        !=======================================================!    
         na=6
         nk=4
         nmtdyn=na*nk
         alp=1.3 ! start value alpha
         kstart=0.003d00 ! start value k
         alpinc = (4./3.) ! increment
         kinc =  (3./2.)     ! increment
        !======================================================! 
        endif
       end select
!---- settings are generated here
        call env%allocate(nmtdyn)   !allocate k(Vbias) and α(Vbias)
        m=0
        do ia=1,na
            k=kstart  ! start value k
            do ik=1,nk
               m = m +1
               env%metadfac(m) =k*env%rednat
               env%metadexp(m) =alp
               k = k / kinc ! increment
            enddo
            alp = alp / alpinc ! increment
        enddo

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++!
       endif
!----
       end if
     end associate settingBool
     return
end subroutine iV2defaultGF

!=========================================================================!
! Dynamically determine the number of normMDs and settings of staticMTDs
! Set their number and the different temperatures.
! Defaults for the static MTDs are more lengthy...
!=========================================================================!
subroutine adjustnormmd(env)
    use crest_data
    implicit none
    !type(options) :: opt
    type(systemdata) :: env
    integer :: ndum

    if(env%rotamermds)then
      !--- first the number of normMDs on low conformers
        if(env%nrotammds.le.0)then !if no user input was set
            env%nrotammds = max(1,nint(float(env%nmetadyn)/4.0d0)) !more but shorter, which is petter concerning parallel efficeiency, 4 for default
        endif

     !--- then the temperature range
        if(env%temps .le. 0)then
           env%temps = 2  !at how many different temperatures? starting at 400k and increasing 100K for each (200 K for -entropy mode)
           if(env%entropic)then
               env%temps = 1
           endif
        endif
    !--- total number of NORMMDs is temps*nrotammds
    endif

!--- settings for static MTDS in entropy mode
     if(env%entropymd)then !special case for entropy mode
         env%emtd%iter = 20            !max number of iterations
         env%emtd%nbias  = min(150,nint(env%tmtd/4)) !max number of bias structures 
         env%emtd%nbiasgrow = min(1.4d0,1.2d0+env%tmtd*1.d-3)!increase of nBias in each cycle
         env%emtd%nMDs   = 36          !number of static MTDs 
         env%emtd%lenfac = 0.5d0     !length (relativ to regular MTDs) 
         env%emtd%temperature = env%nmdtemp !sMTD temperature (default 600 K)
         env%emtd%kpush = 1.d-4+env%tmtd*1.0d-6   !kpush constant PER ATOM, a bit more for flexible systems 1.d-4+env%tmtd*1.d-6 1.5 zu viel, 0.5 zu wenig
         env%emtd%alpha = 1.0d0        !some alpha
         env%emtd%mtdramp = 0.015d0    !parameter to control how "fast" bias is applied in MTD
         if(env%crestver==crest_imtd)then
            if(env%emtd%confthr<0.0d0)then
            env%emtd%confthr = 0.02d0    !if we gain less than x% NEW conformers, exit
            endif
            if(env%emtd%sconvthr<0.0d0)then
            env%emtd%sconvthr = 0.005d0   !if we gain less than x% NEW entropy, exit
            endif
         endif

         if( env%nmdtemp < 0.d0)then   !if temperature is not set by the user
          env%nmdtemp = 600.0d0
         endif

         !-- for the new alternative iMTD-sMTD runtype, re-adjust settings
         if(env%crestver==crest_imtd2)then  
            ndum=2 
            env%emtd%nklist=ndum
            allocate(env%emtd%klist(ndum))
            env%emtd%klist(1) = env%emtd%kpush
            env%emtd%klist(2) = env%emtd%kpush*2.5d0
            env%emtd%nMDs   = 12         !number of static MTDs
            env%emtd%lenfac = 0.5d0      !half the length because we have 2 kpush
            if(env%emtd%confthr<0.0d0)then
            env%emtd%confthr = 0.05d0    !if we gain less than x% NEW conformers, exit
            endif
            if(env%emtd%sconvthr<0.0d0)then
            env%emtd%sconvthr = 0.01d0   !if we gain less than x% NEW entropy, exit
            endif
         endif

         !--- Exclude atoms from static MTD bias
         env%emtd%rmax = 0    ! ignore small rings up to this size in bias
         call mtdatoms('coord',env)
     endif

    return
end subroutine adjustnormmd

