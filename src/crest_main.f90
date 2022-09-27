!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2022 Philipp Pracht
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
!  This is the code of the Conformer-Rotamer Ensemble Sampling Tool (CREST).
!=========================================================================================!
program CREST
      use iso_fortran_env, wp => real64
      use crest_data
      use iomod

      type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
      type(timer)   :: tim
      
      integer :: i,j,l,args
      character(len=:),allocatable :: arg(:)
      character(len=:),allocatable :: infile
      character(len=512) :: thisdir
      character(len=1024) :: cmd
      real(wp) :: dumfloat,dumfloat2,d3,d4,d5,d6,d7,d8
      logical :: ex,ex1,ex2

      call initsignal() !SIGTERM catcher

!=========================================================================================!
!>  Initialize system clock time
      call tim%init(20)

!=========================================================================================!
!> set defaults and pars flags
      args = iargc()
      l=len_trim(cmd)
      allocate(arg(args), source=repeat(' ',l))
      do i=1,args
         call getarg(i,arg(i))
      enddo
      call parseflags(env,arg,args)
      deallocate(arg)

!=========================================================================================!
!> scratch dir handling

      if(env%scratch)then
      call getcwd(thisdir)
      call scrdir(env)
      endif

!=========================================================================================!
!>   OMP_NUM_THREAD handling
!=========================================================================================!
      if(.not.env%autothreads)then
          call ompquickset(env%omp)
      else
         if(.not.env%threadssetmanual)then
          call ompgetauto(env%threads,env%omp,env%MAXRUN)
          call ompprint()
         endif
         call ompautoset(env%threads,0,env%omp,env%MAXRUN,0)   !<--- default
      endif

!=========================================================================================!
!>   DRY run stop
!=========================================================================================!
      if(env%dryrun)then
       call crest_dry(env)
      endif
     
!=========================================================================================!
!> SOME I/O STUFF
!=========================================================================================!
!>--- check for the coord file in the working directory
      if (env%crestver.eq.crest_solv) then
         inquire(file='solute',exist=ex1)
         inquire(file='solvent',exist=ex2)
         if (.not.ex1) then
            error stop 'No solute file found. Exit.'
         else if (.not.ex2) then
            error stop 'No solvent file found. Exit.'
        end if
      else
         inquire(file='coord',exist=ex)
         if(.not.ex)then
            error stop 'No coord file found. Exit.'
         end if
      end if


!>--- call zsort subroutine?
      if(env%autozsort)then
         write(*,'(''-------------------------'')')
         write(*,'(''Starting z-matrix sorting'')')
         write(*,'(''-------------------------'')')
         call zsort
         if(env%onlyZsort)then
         write(*,*)
         write(*,*) 'The z-matrix of the input coord file has been sorted.'
         write(*,*) 'The sorted file in the TM format is called "zcoord"'
         write(*,*)
         write(*,*) 'exit.'
         stop
         endif       
         call rename('zcoord','coord')
         call rmrf('*.zmat')
      end if

!=========================================================================================!
!>        PRE-CONFSEARCH PROPERTY CALCS
!=========================================================================================!
     select case( env%properties )
  !>--- only CREGEN routine
       case( p_cregen )
           call tim%start(1,'CREGEN')
           write(*,*)'Using only the cregen sorting routine.'
           env%cgf(6)=.true.   !write confg output to file <tmp>
           if(env%doNMR)then
              env%cgf(3)=.true.
              if(.not.env%fullcre)then
              env%cgf(2)=.false.
              endif
           endif
           if(env%newcregen)then
           call newcregen(env,0)
           else
           call cregen2(env)
           endif
           if(env%doNMR .and. env%fullcre)then
             call entropic(env,.true.,.false.,.false.,env%ensemblename, &
             &    env%tboltz,dumfloat,dumfloat2)
           endif
           if(env%cluster)then
             call ccegen(env,.true.,'crest_ensemble.xyz')
           endif
           call tim%stop(1)
           call propquit(tim)
  !>--- only ensemble comparison
       case( p_compare )
        call compare_ensembles(env)         !compare ensembles
        call propquit(tim)
  !>--- protonation tool
       case( p_protonate )
        call protonate(env,tim)
        call propquit(tim)
  !>--- deprotonation
       case( p_deprotonate )
        call deprotonate(env,tim)
        call propquit(tim)
  !>--- tautomerization
       case( p_tautomerize )
        call tautomerize(env,tim)
        call propquit(tim)
  !>--- extended tautomerization
       case( p_tautomerize2 )
        call tautomerize_ext(env%ensemblename,env,tim)
        call propquit(tim)
  !>--- stereoisomerization
       case( p_isomerize )
        call stereoisomerize(env,tim)
        call propquit(tim)    

  !>--- reactor setup
       case( p_reactorset )
        call reactor_setup(env)

        stop
  !>--- enhanched ensemble entropy
       case( p_CREentropy )
        call entropic(env,.true.,.true.,.false.,env%ensemblename, &
        &    env%tboltz,dumfloat,dumfloat2)
        call propquit(tim)
  !>--- calculate hessians and average thermo. contrib
       case( p_rrhoaverage )
         call tim%start(4,'freq+thermo')
         call calcSrrhoav(env,env%ensemblename)  
         call tim%stop(4)
         call propquit(tim)
  !>--- to PCA and k-Means clustering for given file
       case( p_cluster )
        call ccegen(env,.true.,env%ensemblename) 
        call propquit(tim)
  !>--- properties for enesemble file
       case( p_propcalc )
        call propcalc(env%ensemblename,env%properties2,env,tim)        
        call propquit(tim)
  !>--- calculate potential correction for acid/base reaction
       case( p_acidbase )
         call tim%start(4,'acid/base')
         if(env%ptb%pka_mode==0)then
         call acidbase(env,env%ptb%pka_acidensemble,env%ptb%pka_baseensemble,env%chrg,.true., &
             & .false.,dumfloat,.false.,d3,d4,d5,d6,d7,d8)
         else
         call rewrite_AB_ensemble(env,env%ptb%pka_acidensemble,env%ptb%pka_baseensemble)
         endif
         call tim%stop(4)
         call propquit(tim)  
  !>--- calculate potential correction for acid/base reaction
       case( p_ligand )
         call tim%start(4,'')
         call ligandtool(env%ptb%infile,env%ptb%newligand, &
         &    env%ptb%centeratom,env%ptb%ligand)
         call tim%stop(4)
         call propquit(tim)  
  !>--- wrapper for the thermo routine       
       case( p_thermo )
         call tim%start(4,'')
         call thermo_mini(env)
         call tim%stop(4)
         call propquit(tim)  
  !>--- ensemble merging tool
       case( p_gesc1,p_gesc2 )
         call tim%start(9,'')
         call biasmerge(env)
         call tim%stop(9)
        if( env%properties == -9224) call propquit(tim)
  !>--- do nothing here
       case default
         continue
     end select
!=========================================================================================!
!>         PRE-OPTIMIZATION OF THE GEOMETRY
!=========================================================================================!
      if(env%preopt)then
         call xtbopt(env)
      else if(env%presp)then
         call xtbsp(env)
      endif
!=========================================================================================!
!>         CONFORMATION SEARCH CALLS START HERE
!=========================================================================================!
      select case( env%crestver )
        case( crest_mfmdgc )
           call confscript1(env,tim)  !> MF-MD-GC algo
        case( crest_imtd,crest_imtd2 )
           call confscript2i(env,tim) !> MTD-GC algo
        case( crest_mdopt )
           call mdopt(env,tim)        !> MDOPT
        case( crest_screen )
           call screen(env,tim)       !> SCREEN
        case( crest_nano )                      
           call reactor(env,tim)      !> NANO-REACTOR
        case( crest_compr )                                        
           call compress(env,tim)     !> MTD COMPRESS mode   
        case( crest_msreac ) 
           call msreact_handler(env,tim) !> MSREACT sub-program   
        case( crest_pka )
           call pkaquick(env,tim)   
        case( crest_solv )             !> microsolvation tools   
           call crest_solvtool(env, tim) 
        case default
           continue
      end select

      if((any((/crest_mfmdgc,crest_imtd,crest_imtd2/)==env%crestver)) &
      &  .and.( env%outputsdf .or. env%sdfformat) )then    
         !call wrsdfens(env%sdf,conformerfile,conformerfilebase//'.sdf')
         call new_wrsdfens(env,conformerfile,conformerfilebase//'.sdf',.false.)
      endif

!=========================================================================================!
!>        POST-CONFSEARCH PROPERTY CALCS
!=========================================================================================!
     if(env%npq .gt. 0)then
     infile="crest_rotamers.xyz"    
     do i=1,env%npq
       j=env%pqueue(i)
       select case( j )
         case( 1:8,10,20,100 )
           call propcalc(conformerfile,j,env,tim)  
         case( 45 ) 
            call tim%start(15,'entropy eval.') 
            call newentropyextrapol(env)
            call tim%stop(15)
         case( 50:59 )  !hybrid reoptimization (e.g. gfn2@gff)
           call propcalc(infile,j,env,tim)
           infile='crest_reopt.xyz'  
         case( 70  ) !PCA and clustering
           call ccegen(env,.true.,conformerfile)
         case( 555 )
           call tautomerize_ext(infile,env,tim)
         case default
           continue
       end select
     enddo
     endif

!=========================================================================================!
!> go back from scratch directory
     if(env%scratch)then
       call chdir(thisdir)
       call scrend(env)
     endif


!=========================================================================================!
!> Evaluate and print timings
      call eval_timer(tim)
      write(*,*)
      write(*,*)'CREST terminated normally.'
!> end of main program
end program CREST

!=========================================================================================!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!=========================================================================================!
