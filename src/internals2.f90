!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Philipp Pracht
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

module INTERNALS_mod
   use crest_parameters
   use adjacency
   use geo
   use strucrd, only: i2e
   implicit none
   public


!========================================================================================!
!========================================================================================!
contains !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!


  SUBROUTINE BETTER_XYZINT(nat,xyz,amat,NA,NB,NC,geo)
!***********************************************************************
!*
!* BETTER_XYZINT works out the internal coordinates of a molecule
!* based on a given connectivity.
!*
!*  INPUT:
!*       nat    = number of atoms
!*       xyz    = cartesian coordinates
!*       amat   = adjacency matrix
!*
!*  OUTPUT:
!*       NA,NB,NC = mapping of internals
!*       geo      = the zmatrix (angles will always be in radians!)
!*
!***********************************************************************
    use crest_parameters
    use adjacency
    implicit none
    !> INPUT
    integer,intent(in)  :: nat
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in)  :: amat(nat,nat)
    !> OUTPUT
    integer,intent(out)  :: NA(nat)
    integer,intent(out)  :: NB(nat)
    integer,intent(out)  :: NC(nat)
    real(wp),intent(out) :: geo(3,nat)
    !> LOCAL
    integer :: V,i,j,k,l
    real(wp),allocatable :: dist(:,:)
    integer,allocatable  :: prev(:,:)
    logical,allocatable  :: taken(:)
    integer,allocatable  :: tmppath(:),path(:)
    !>-- init
    NA = 0
    NB = 0
    NC = 0
    geo = 0.0_wp

    V = nat
    allocate (dist(V,V),source=0.0_wp)
    allocate (prev(V,V),source=0)
    call FloydWarshall(V,Amat,dist,prev)

    !>-- set up
    allocate (taken(nat),source=.false.)
    call init_first_four(nat,xyz,dist,prev,taken,NA,NB,NC)
    allocate (tmppath(V),path(4),source=0)
    do while (any(.not.taken(:)))
      call walkbonds(nat,xyz,dist,prev,taken,4,tmppath,path)
      k = path(1)
      taken(k) = .true.
      NA(k) = path(2)
      NB(k) = path(3)
      NC(k) = path(4)
    end do
    deallocate (path,tmppath)
    deallocate (prev,dist)

    !>-- finally, generate zmatrix and save to geo
    call XYZGEO2(nat,xyz,NA,NB,NC,1.0_wp,GEO)
    RETURN
  CONTAINS
    subroutine init_first_four(V,xyz,dist,prev,taken,NA,NB,NC)
      implicit none
      integer :: V
      real(wp),intent(in) :: xyz(3,V)
      real(wp) :: dist(V,V)
      integer  :: prev(V,V)
      logical  :: taken(V)
      integer  :: NA(V),NB(V),NC(V)
      integer,allocatable :: tmppath(:),path(:)
      integer :: i,j,k,l
      !> first at origin
      taken(1) = .true.
      !> if we have only 3 or less atoms, don't bother
      if (V <= 3) then
        taken = .true.
        if (V > 1) NA(2) = 1
        if (V > 2) NA(3) = 2
        if (V > 2) NB(3) = 1
        return
      end if
      !> set up zmat entries 2 to 4
      allocate (tmppath(V),path(4),source=0)
      do i = 2,4
        call walkbonds(V,xyz,dist,prev,taken,4,tmppath,path)
        k = path(1)
        taken(k) = .true.
        NA(k) = path(2)
        if (i > 2) NB(k) = path(3)
        if (i > 3) NC(k) = path(4)
      end do
    end subroutine init_first_four
    recursive subroutine walkbonds(V,xyz,dist,prev,taken,nmax,tmppath,path)
      implicit none
      integer,intent(in) :: V
      real(wp),intent(in) :: xyz(3,V)
      real(wp),intent(in) :: dist(V,V)
      integer,intent(in)  :: prev(V,V)
      logical,intent(in)  :: taken(V)
      integer,intent(in) :: nmax
      integer,intent(inout) :: path(nmax)
      integer,intent(inout) :: tmppath(V)
      integer :: newnmax,i,j
      integer :: k,lpath
      logical :: found
      if (nmax <= 1) then !> return condition 1
        do i = 1,V
          if (.not.taken(i)) then
            k = i
            exit
          end if
        end do
        path(1) = k
        return
      end if
      found = .false.
      ILOOP: do i = 1,V
        !> skip to next untaken
        if (taken(i)) cycle
        JLOOP: do j = 1,V
          !> look for path in taken
          if (i == j) cycle
          tmppath = 0
          call getPathFW(V,prev,i,j,tmppath,lpath)
          if (lpath == nmax) then
            if (all(taken(tmppath(2:lpath)))) then
              !> save the first viable path
              found = .true.
              path(1:nmax) = tmppath(1:lpath)
              !> and exit further checks
              exit ILOOP
            end if
          end if
        end do JLOOP
      end do ILOOP
      if (found) then
        return !> retrun condition 2
      else
        !> if no path was found, search a shorter one and
        !> fill up the last spot by a distance criterium.
        !> this is done recursively until nmax=1 (see return condition 1),
        !> so path is ensured to exist!
        newnmax = nmax-1
        call walkbonds(V,xyz,dist,prev,taken,newnmax,tmppath,path(1:newnmax))
        i = path(1)
        call closesttaken(V,xyz,i,taken,path(1:newnmax),j)
        path(nmax) = j
      end if
    end subroutine walkbonds
    subroutine closesttaken(V,xyz,a,taken,path,b)
      use geo,only:distance
      implicit none
      integer,intent(in)  :: V,a
      real(wp),intent(in) :: xyz(3,V)
      logical,intent(in)  :: taken(V)
      integer,intent(in)  :: path(:)
      integer,intent(out) :: b
      integer :: i
      real(wp) :: d,dtmp
      b = 0
      d = huge(d)
      do i = 1,V
        if (taken(i).and..not.any(path .eq. i)) then
          dtmp = distance(xyz(:,a),xyz(:,i))
          if (dtmp < d) then
            b = i
          end if
        end if
      end do
      if (b /= 0) return
      !> if no closest could be selected from taken ∩ path
      !> check all other
      do i = 1,V
        if (i == a) cycle
        if (.not.taken(i)) then
          dtmp = distance(xyz(:,a),xyz(:,i))
          if (dtmp < d) then
            b = i
          end if
        end if
      end do
    end subroutine closesttaken
  END SUBROUTINE BETTER_XYZINT

  SUBROUTINE BANGLE2(XYZ,I,J,K,ANGLE)
    use crest_parameters
    use geo
    implicit none
    real(wp),intent(in) :: xyz(3,*)
    integer,intent(in)  :: i,j,k
    real(wp),intent(out) :: ANGLE
    real(wp) :: u(3),v(3)
    u(1:3) = xyz(1:3,I)-xyz(1:3,J)
    v(1:3) = xyz(1:3,K)-xyz(1:3,J)
    ANGLE = tangle(u,v)
  END SUBROUTINE BANGLE2

  SUBROUTINE DIHED2(XYZ,I,J,K,L,ANGLE)
    use crest_parameters
    use geo
    implicit none
    real(wp),intent(in) :: xyz(3,*)
    integer,intent(in)  :: i,j,k,l
    real(wp),intent(out) :: ANGLE
    real(wp) :: u(3),v(3),w(3)
    u(1:3) = xyz(1:3,J)-xyz(1:3,I)
    v(1:3) = xyz(1:3,K)-xyz(1:3,J)
    w(1:3) = xyz(1:3,L)-xyz(1:3,K)
    ANGLE = dihedral(u,v,w)
  END SUBROUTINE DIHED2

!========================================================================================!

  SUBROUTINE XYZGEO2(nat,xyz,NA,NB,NC,DEGREE,GEO)
    use crest_parameters
!***********************************************************************
!*
!*   XYZGEO CONVERTS COORDINATES FROM CARTESIAN TO INTERNAL.
!*
!*    INPUT:
!*      nat    = number of atoms
!*      xyz    = array of cartesian coordinates
!*      NA     = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
!*               BY DISTANCE
!*      NB     = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
!*               BY ANGLE
!*      NC     = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
!*               BY DIHEDRAL
!*      DEGREE = should be equal 1.0 if output is in Rad,
!*               or 180/pi if output is in degree
!*
!*    OUTPUT:
!*      GEO    = internal coordinates consisting of a
!*               distance and two angles for each atom
!*
!***********************************************************************
    implicit none
    integer,intent(in)   :: nat
    real(wp),intent(in)  :: xyz(3,nat)
    integer,intent(in)   :: NA(nat),NB(nat),NC(nat)
    real(wp),intent(in)  :: DEGREE
    real(wp),intent(out) :: GEO(3,nat)

    integer :: i,j,k,l

    GEO(:,:) = 0.0_wp

    !> first atom at zero
    GEO(1:3,1) = 0.0_wp
    !> select second atom
    do i = 2,nat
      j = NA(i)
      k = NB(i)
      l = NC(i)
      GEO(1,i) = sqrt((xyz(1,i)-xyz(1,j))**2+ &
    &                  (xyz(2,i)-xyz(2,j))**2+ &
    &                  (xyz(3,i)-xyz(3,j))**2)
      if (k /= 0) then
        call BANGLE2(xyz,i,j,k,geo(2,i))
        geo(2,i) = geo(2,i)*DEGREE
        if (l /= 0) then
          call DIHED2(xyz,i,j,k,l,geo(3,i))
          geo(3,i) = geo(3,i)*DEGREE
        end if
      end if
    end do
    return
  END SUBROUTINE XYZGEO2

!========================================================================================!
  SUBROUTINE GMETRY2(nat,geo,coord,na,nb,nc)
!***********************************************************************
!*
!*   GMETRY2 converts coordinates from internal to cartesian
!*
!*    INPUT:
!*      nat    = number of atoms
!*      geo    = array of internal coordinates (distance in Bohr, 2x angle in rad)
!*      NA     = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
!*               BY DISTANCE
!*      NB     = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
!*               BY ANGLE
!*      NC     = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
!*               BY DIHEDRAL
!*
!*    OUTPUT:
!*      COORD  = cartesian coordinates
!*
!***********************************************************************
    use crest_parameters
    implicit none
    integer,intent(in)   :: nat
    real(wp),intent(in)  :: geo(3,nat)
    integer,intent(in)   :: na(nat),nb(nat),nc(nat)
    real(wp),intent(out) :: coord(3,nat)

    real(wp) :: COSC,COSA,XB,YB,ZB,RBC,tmp
    real(wp) :: XA,YA,ZA,XPA,XPB,XYB,YPA,XQA,ZQA,YZA
    real(wp) :: XD,YD,ZD,SINA,SIND,COSD,XRD
    real(wp) :: YPD,ZPD,XPD,ZQD,XQD,YQD
    real(wp) :: COSTH,SINTH,COSPH,SINPH,COSKH,SINKH
    integer :: i,j,K,l,MA,MB,MC
    logical,allocatable :: taken(:)
    real(wp),parameter :: verylarge = 1.0d100

    !COORD(:,:) = 0.0_wp
    COORD(:,:) = huge(tmp)
    !allocate (taken(nat))
    !taken(:) = .false.

!>--- first atom at origin
    COORD(1:3,1) = 0.0_wp
    !taken(1) = .true.
!>--- second atom
    do i = 2,nat
      if (na(i) == 1.and.nb(i) == 0.and.nc(i) == 0) then
        COORD(1,i) = GEO(1,i)
        COORD(2,i) = 0.0_wp
        COORD(3,i) = 0.0_wp
        !taken(i) = .true.
      else
        cycle
      end if
    end do
!>--- third atom
    do i = 2,nat
      if (na(i) /= 0.and.nb(i) /= 0.and.nc(i) == 0) then
        COSC = COS(GEO(2,i))
        j = na(i)
        if (j == 1) then
          COORD(1,i) = COORD(1,j)+GEO(1,i)*COSC
        else
          COORD(1,i) = COORD(1,j)-GEO(1,i)*COSC
        end if
        COORD(2,i) = GEO(1,i)*SIN(GEO(2,i))
        COORD(3,i) = 0.0_wp
        !taken(i) = .true.
      else
        cycle
      end if
    end do

   ! TAKELOOP: do while (any(.not.taken(:)))
   TAKELOOP : do while ( any(COORD(:,:) > verylarge))
      DO I = 2,nat
!>--- CYCLE the atom if already generated
        !if (taken(i)) cycle
        if (COORD(1,i) < verylarge ) cycle
!>--- CYCLE if any of the depending atoms have not been generated
       ! if ((.not.taken(NA(i))).or.(.not.taken(NB(i))) &
       !&   .or.(.not.taken(NC(i)))) cycle
        if ( coord(1,NA(i)) > verylarge .or. &
        &    coord(1,NB(i)) > verylarge .or. &
        &    coord(1,NC(i)) > verylarge) cycle

        COSA = COS(GEO(2,I))
        MB = NB(I)
        MC = NA(I)
        XB = COORD(1,MB)-COORD(1,MC)
        YB = COORD(2,MB)-COORD(2,MC)
        ZB = COORD(3,MB)-COORD(3,MC)
        RBC = 1.0D00/DSQRT(XB*XB+YB*YB+ZB*ZB)

        IF (ABS(COSA) .GE. 0.9999999999D0) THEN
!>--- ATOMS MC, MB, AND (I) ARE COLLINEAR
          RBC = GEO(1,I)*RBC*COSA
          COORD(1,I) = COORD(1,MC)+XB*RBC
          COORD(2,I) = COORD(2,MC)+YB*RBC
          COORD(3,I) = COORD(3,MC)+ZB*RBC
          CYCLE
        END IF

!>--- THE ATOMS ARE NOT COLLINEAR
40      MA = NC(I)
        XA = COORD(1,MA)-COORD(1,MC)
        YA = COORD(2,MA)-COORD(2,MC)
        ZA = COORD(3,MA)-COORD(3,MC)

!>-- ROTATE ABOUT THE Z-AXIS TO MAKE YB=0, AND XB POSITIVE.  IF XYB IS
!>-- TOO SMALL, FIRST ROTATE THE Y-AXIS BY 90 DEGREES.
        XYB = DSQRT(XB*XB+YB*YB)
        K = -1
        !IF (XYB.GT.0.1D00) GO TO 50
        IF (XYB .LE. 0.1D00) THEN
          XPA = ZA
          ZA = -XA
          XA = XPA
          XPB = ZB
          ZB = -XB
          XB = XPB
          XYB = DSQRT(XB*XB+YB*YB)
          K = 1
        END IF

!>--- ROTATE ABOUT THE Y-AXIS TO MAKE ZB VANISH
50      COSTH = XB/XYB
        SINTH = YB/XYB
        XPA = XA*COSTH+YA*SINTH
        YPA = YA*COSTH-XA*SINTH
        SINPH = ZB*RBC
        COSPH = DSQRT(ABS(1.D00-SINPH*SINPH))
        XQA = XPA*COSPH+ZA*SINPH
        ZQA = ZA*COSPH-XPA*SINPH

!>--- ROTATE ABOUT THE X-AXIS TO MAKE ZA=0, AND YA POSITIVE.
        YZA = DSQRT(YPA*YPA+ZQA*ZQA)
        IF (YZA .LT. 1.D-4) THEN
          COSKH = 1.D0
          SINKH = 0.D0
        ELSE
          COSKH = YPA/YZA
          SINKH = ZQA/YZA
        END IF

!>--   COORDINATES :-   A=(XQA,YZA,0),   B=(RBC,0,0),  C=(0,0,0)
!>--   NONE ARE NEGATIVE.
!>--   THE COORDINATES OF i ARE EVALUATED IN THE NEW FRAME.
        SINA = SIN(GEO(2,I))
        SIND = -SIN(GEO(3,I))
        COSD = COS(GEO(3,I))
        XD = GEO(1,I)*COSA
        YD = GEO(1,I)*SINA*COSD
        ZD = GEO(1,I)*SINA*SIND

!>--- TRANSFORM THE COORDINATES BACK TO THE ORIGINAL SYSTEM.
        YPD = YD*COSKH-ZD*SINKH
        ZPD = ZD*COSKH+YD*SINKH
        XPD = XD*COSPH-ZPD*SINPH
        ZQD = ZPD*COSPH+XD*SINPH
        XQD = XPD*COSTH-YPD*SINTH
        YQD = YPD*COSTH+XPD*SINTH
        !IF (K.LT.1) GO TO 90
        IF (K .GE. 1) THEN
          XRD = -ZQD
          ZQD = XQD
          XQD = XRD
        END IF
90      COORD(1,I) = XQD+COORD(1,MC)
        COORD(2,I) = YQD+COORD(2,MC)
        COORD(3,I) = ZQD+COORD(3,MC)

        !taken(I) = .true.
      END DO
    end do TAKELOOP
100 CONTINUE
110 CONTINUE

    !deallocate (taken)
    RETURN
  END SUBROUTINE GMETRY2

!========================================================================================!

  subroutine print_zmat(ch,nat,at,geo,na,nb,nc,nice)
    use crest_parameters
    use strucrd
    implicit none
    integer,intent(in) :: ch,nat
    integer,intent(in) :: na(nat),nb(nat),nc(nat)
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: geo(3,nat)
    logical,intent(in) :: nice
    character(len=120) :: atmp
    integer :: i

    if(nice)then
    write (ch,'(1x,a5,1x,a12,1x,a12,1x,a12,a5,a5,a5)') 'A','d(AB)','θ(ABC)','ϕ(ABCD)','B','C','D'
    do i = 1,nat
      if (na(i) .ne. 0) then
        if (nb(i) .ne. 0) then
          if (nc(i) .ne. 0) then
            write (atmp,'(1x,i5,1x,3f12.4,3i5)') i,geo(1:3,i),na(i),nb(i),nc(i)
          else
            write (atmp,'(1x,i5,1x,2f12.4,a12,2i5)') i,geo(1:2,i),'-',na(i),nb(i)
            atmp = trim(atmp)//'    -'
          end if
        else
          write (atmp,'(1x,i5,1x,f12.4,a12,a12,i5)') i,geo(1,i),'-','-',na(i)
          atmp = trim(atmp)//'    -    -'
        end if
      else
        write (atmp,'(1x,i5,1x,a12,a12,a12)') i,'-','-','-'
        atmp = trim(atmp)//'    -    -    -'
      end if
      write (ch,'(a)') trim(atmp)
    end do
    else
      write(ch,*) nat
     do i = 1,nat
      write (atmp,'(1x,a,1x,3f12.4,3i5)') i2e(at(i),'nc'),geo(1:3,i),na(i),nb(i),nc(i)
      write (ch,'(a)') trim(atmp)
     end do

    endif

  end subroutine print_zmat

!========================================================================================!

  subroutine rd_zmat(fname,nat,at,zmat,na,nb,nc)
    use crest_parameters
    use strucrd
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nat
    integer,intent(out) :: at(nat)
    real(wp),intent(out) :: zmat(3,nat)
    integer,intent(out) :: na(nat),nb(nat),nc(nat)
    integer :: ich,i,j,k,io
    character(len=4) :: sym
    character(len=300) line

    at(:) = 0
    na(:) = 0 
    nb(:) = 0
    nc(:) = 0
    zmat(:,:) = 0.0_wp
    open(newunit=ich,file=trim(fname))
    read(ich,*) j
    if(j /= nat) error stop 'Nat mismatch in rd_zmat()' 
    do i=1,nat
      read(ich,'(a)') line
      call zmatline(line,sym,zmat(:,i),na(i),nb(i),nc(i),io)
      if(io /= 0) error stop 'error while reading zmat'
      at(i) = e2i(sym)  
    enddo 
    close(ich)
   contains
    subroutine zmatline(line,sym,xyz,a,b,c,io)
    implicit none
    character(len=*) :: line
    character(len=*) :: sym
    real(wp) :: xyz(3)
    integer :: a,b,c
    integer,intent(out) :: io

    io = 0
    xyz(1:3) = 0
    a = 0
    b = 0
    c = 0
    read (line,*,iostat=io) sym,xyz(1:3),a,b,c
    if (io .ne. 0) then
      read (line,*,iostat=io) sym,xyz(1:2),a,b
    if(io.ne.0)then
       read (line,*,iostat=io) sym,xyz(1),a
    if(io.ne.0)then
       read (line,*,iostat=io) sym
    endif
    endif
    end if

    return
    end subroutine zmatline

  end subroutine rd_zmat

!========================================================================================!
end module INTERNALS_mod
