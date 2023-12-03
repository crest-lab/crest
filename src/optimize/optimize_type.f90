!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2017-2020 Stefan Grimme, Sebastian Ehlert (xtb)
! Copyright (C) 2021 - 2022 Philipp Pracht
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
!
! Routines were adapted from the xtb code (github.com/grimme-lab/xtb)
! under the Open-source software LGPL-3.0 Licencse.
!================================================================================!
module optimize_type
   use iso_fortran_env, only: wp=>real64
   use optimize_maths, only: detrotra8
   implicit none

   public :: optimizer
   public :: convergence_log
   private

   type optimizer
      integer  :: n    !< number of atoms
      integer  :: n3   !< dimension of hessian
      integer  :: nvar !< actual dimension
      real(wp) :: hlow
      real(wp) :: hmax
      real(wp),allocatable :: hess(:)
      real(wp),allocatable :: B(:,:)
      real(wp),allocatable :: eigv(:)
      real(wp),allocatable :: coord(:)
      real(wp),allocatable :: xyz(:,:)
   contains
   procedure :: allocate => allocate_anc
   procedure :: allocate2 => allocate_opt
   procedure :: deallocate => deallocate_anc
   procedure :: write => write_anc
   procedure :: new => generate_anc_blowup
   procedure :: get_cartesian
   end type optimizer


  type :: convergence_log
    integer :: nlog
    real(wp),allocatable :: elog(:)
    real(wp),allocatable :: glog(:)
  contains
    procedure :: set_eg_log
    procedure :: get_averaged_energy
    procedure :: get_averaged_gradient
  end type convergence_log
  interface convergence_log
    module procedure new_convergence_log
  end interface convergence_log
  

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

subroutine allocate_anc(self,n,nvar,hlow,hmax)
   implicit none
   class(optimizer),intent(inout)  :: self
   integer,      intent(in)     :: n
   integer,      intent(in)     :: nvar
   integer                      :: n3
   real(wp),intent(in),optional :: hlow
   real(wp),intent(in),optional :: hmax
   n3 = 3*n
   call self%deallocate
   self%n    = n
   self%n3   = 3*n
   self%nvar = nvar
   if (present(hlow)) self%hlow = hlow
   if (present(hmax)) self%hmax = hmax
   allocate( self%hess(nvar*(nvar+1)/2), source = 0.0_wp )
   allocate( self%B(n3,n3),              source = 0.0_wp )
   allocate( self%eigv(n3),              source = 0.0_wp )
   allocate( self%coord(nvar),           source = 0.0_wp )
   allocate( self%xyz(3,n),              source = 0.0_wp )
end subroutine allocate_anc


subroutine allocate_opt(self,n)!,nvar)
   implicit none
   class(optimizer),intent(inout)  :: self
   integer,      intent(in)     :: n
   !integer,      intent(in)     :: nvar
   integer :: nvar
   integer                      :: n3
   !real(wp),intent(in),optional :: hlow
   !real(wp),intent(in),optional :: hmax
   n3 = 3*n
   call self%deallocate
   self%n    = n
   self%n3   = 3*n
   nvar = 3*n
   self%nvar = nvar
   allocate( self%hess(nvar*(nvar+1)/2), source = 0.0_wp )
   !allocate( self%B(n3,n3),              source = 0.0_wp )
   allocate( self%eigv(n3),              source = 0.0_wp )
   !allocate( self%coord(nvar),           source = 0.0_wp )
   !allocate( self%xyz(3,n),              source = 0.0_wp )
end subroutine allocate_opt



!========================================================================================!

subroutine deallocate_anc(self)
   implicit none
   class(optimizer),intent(inout) :: self
   self%n    = 0
   self%n3   = 0
   self%nvar = 0
   if (allocated(self%hess )) deallocate( self%hess  )
   if (allocated(self%B    )) deallocate( self%B     )
   if (allocated(self%eigv )) deallocate( self%eigv  )
   if (allocated(self%coord)) deallocate( self%coord )
   if (allocated(self%xyz  )) deallocate( self%xyz   )
end subroutine deallocate_anc

!========================================================================================!
!> @brief print information about current approximate normal coordinates to unit
subroutine write_anc(self,iunit,comment)
   implicit none
   class(optimizer),   intent(in) :: self    !< approximate normal coordinates
   integer,         intent(in) :: iunit   !< file handle
   character(len=*),intent(in) :: comment !< name of the variable
   character(len=*),parameter  :: dfmt = '(1x,a,1x,"=",1x,g0)'

   write(iunit,'(72(">"))')
   write(iunit,'(1x,"*",1x,a)') "Writing 'optimizer' class"
   write(iunit,'(  "->",1x,a)') comment
   write(iunit,'(72("-"))')
   write(iunit,'(1x,"*",1x,a)') "status of the fields"
   write(iunit,dfmt) "integer :: n           ",self%n
   write(iunit,dfmt) "integer :: n3          ",self%n3
   write(iunit,dfmt) "integer :: nvar        ",self%nvar
   write(iunit,dfmt) "real    :: hlow        ",self%hlow
   write(iunit,dfmt) "real    :: hmax        ",self%hmax
   write(iunit,'(72("-"))')
   write(iunit,'(1x,"*",1x,a)') "allocation status"
   write(iunit,dfmt) "allocated? hess(:)     ",allocated(self%hess)
   write(iunit,dfmt) "allocated? B(:)        ",allocated(self%B)
   write(iunit,dfmt) "allocated? eigv(:)     ",allocated(self%eigv)
   write(iunit,dfmt) "allocated? coord(:)    ",allocated(self%coord)
   write(iunit,dfmt) "allocated? xyz(:,:)    ",allocated(self%xyz)
   write(iunit,'(72("-"))')
   write(iunit,'(1x,"*",1x,a)') "size of memory allocation"
   if (allocated(self%hess)) then
   write(iunit,dfmt) "size(1) :: hess(*)     ",size(self%hess,1)
   endif
   if (allocated(self%B)) then
   write(iunit,dfmt) "size(1) :: B(*,:)      ",size(self%B,1)
   write(iunit,dfmt) "size(2) :: B(:,*)      ",size(self%B,2)
   endif
   if (allocated(self%eigv)) then
   write(iunit,dfmt) "size(1) :: eigv(*)     ",size(self%eigv,1)
   endif
   if (allocated(self%coord)) then
   write(iunit,dfmt) "size(1) :: coord(*)    ",size(self%coord,1)
   endif
   if (allocated(self%xyz)) then
   write(iunit,dfmt) "size(1) :: xyz(*,:)    ",size(self%xyz,1)
   write(iunit,dfmt) "size(2) :: xyz(:,*)    ",size(self%xyz,2)
   endif
   write(iunit,'(72("<"))')
end subroutine write_anc

!========================================================================================!

subroutine generate_anc_blowup(self,xyz,hess,pr,linear,fail)
   implicit none
   class(optimizer),intent(inout) :: self
   real(wp),     intent(in)    :: xyz(3,self%n)
   real(wp),     intent(inout) :: hess(self%n3,self%n3)
   logical,      intent(in)    :: pr
   logical,      intent(in)    :: linear
   logical,      intent(out)   :: fail

   real(wp),parameter   :: thr1 = 1.0e-10_wp
   real(wp),parameter   :: thr2 = 1.0e-11_wp
   integer, parameter   :: maxtry = 4
   integer  :: i,itry
   integer  :: nvar
   integer  :: info
   integer  :: lwork
   integer  :: liwork
   integer, allocatable :: iwork(:)
   real(wp) :: elow,damp,thr
   real(wp),allocatable :: aux(:)

   !> LAPACK
   external :: dsyevd

   fail = .false.
   self%xyz = xyz

   thr = thr2
   lwork  = 1 + 6*self%n3 + 2*self%n3**2
   liwork = 8 * self%n3

   allocate(iwork(liwork), source = 0 )
   allocate(aux(lwork), source = 0.0_wp )

   call dsyevd('V','U',self%n3,hess,self%n3,self%eigv, &
      &        aux,lwork,iwork,liwork,info)

   deallocate(aux,iwork)

   call detrotra8(linear,self%n,self%xyz,hess,self%eigv) 

   elow = minval(self%eigv,mask=(abs(self%eigv) > thr1))

   damp = max(self%hlow - elow,0.0_wp)
   where(abs(self%eigv) > thr2) self%eigv = self%eigv + damp
!   do i = 1, self%n3
!      if (abs(self%eigv(i)) > thr2 ) self%eigv(i) = self%eigv(i) + damp
!   enddo

   if(pr)then
      write(*,*) 'Shifting diagonal of input Hessian by ', damp
      write(*,*) 'Lowest  eigenvalues of input Hessian'
      write(*,'(6F12.6)')(self%eigv(i),i=1,min(18,self%n3))
      write(*,*) 'Highest eigenvalues'
      write(*,'(6F12.6)')(self%eigv(i),i=self%n3-5,self%n3)
      write(*,*)
   endif

   fail = .true.
   get_anc: do itry = 1, maxtry
      self%B = 0.0_wp
      self%hess = 0.0_wp
      nvar = 0
      ! take largest (positive) first
      do i = self%n3, 1, -1
         if (abs(self%eigv(i)) > thr .and. nvar < self%nvar) then
            nvar = nvar+1
            self%B(:,nvar) = hess(:,i)
            self%hess(nvar+nvar*(nvar-1)/2) = &
               min(max(self%eigv(i),self%hlow),self%hmax)
         endif
      enddo

      if (nvar.ne.self%nvar) then
         thr = thr * 0.1_wp
         cycle get_anc
      endif

      fail = .false.
      exit get_anc
   enddo get_anc

   if (fail) then
      if(pr) write(*,*) 'nvar, self%nvar',nvar,self%nvar
      return
   end if

   call sort(self%n3,self%nvar,self%hess,self%B)

   self%coord = 0.0_wp
   return
end subroutine generate_anc_blowup

!========================================================================================!

subroutine generate_anc_packed(self,xyz,hess,pr,fail)
   implicit none
   class(optimizer),intent(inout) :: self
   real(wp),     intent(in)    :: xyz(3,self%n)
   real(wp),     intent(inout) :: hess(self%n3*(self%n3+1)/2)
   logical,      intent(in)    :: pr
   logical,      intent(out)   :: fail

   real(wp),parameter   :: thr1 = 1.0e-10_wp
   real(wp),parameter   :: thr2 = 1.0e-11_wp
   integer, parameter   :: maxtry = 4
   integer  :: i,itry
   integer  :: nvar
   integer  :: info
   integer  :: lwork
   integer  :: liwork
   integer, allocatable :: iwork(:)
   real(wp) :: elow,damp,thr
   real(wp),allocatable :: aux(:)
   real(wp),allocatable :: u(:,:)

   !> LAPACK
   external :: dspevd

   self%xyz = xyz

   thr = thr2
   lwork  = 1 + 6*self%n3 + 2*self%n3**2
   liwork = 8 * self%n3

   allocate(iwork(liwork), source = 0 )
   allocate(aux(lwork), source = 0.0_wp )
   allocate(u(self%n3,self%n3), source = 0.0_wp )

   call dspevd('V','U',self%n3,hess,self%eigv,u,self%n3, &
      &        aux,lwork,iwork,liwork,info)

   !elow = 1.0e+99_wp
   elow = minval(self%eigv,mask=(abs(self%eigv) > thr1))
   !do i = 1, self%n3
   !   if (abs(self%eigv(i)) > thr1 ) elow = min(elow,self%eigv(i))
   !enddo

   damp = max(self%hlow - elow,0.0_wp)
   where(abs(self%eigv) > thr2) self%eigv = self%eigv + damp
!   do i = 1, self%n3
!      if (abs(self%eigv(i)) > thr2 ) self%eigv(i) = self%eigv(i) + damp
!   enddo

   if(pr)then
      write(*,*) 'Shifting diagonal of input Hessian by ', damp
      write(*,*) 'Lowest  eigenvalues of input Hessian'
      write(*,'(6F12.6)')(self%eigv(i),i=1,min(18,self%n3))
      write(*,*) 'Highest eigenvalues'
      write(*,'(6F12.6)')(self%eigv(i),i=self%n3-5,self%n3)
      write(*,*)
   endif

   fail = .true.
   get_anc: do itry = 1, maxtry
      self%B = 0.0_wp
      self%hess = 0.0_wp
      nvar = 0
      ! take largest (positive) first
      do i = self%n3, 1, -1
         if (abs(self%eigv(i)) > thr .and. nvar < self%nvar) then
            nvar = nvar+1
            self%B(:,nvar) = u(:,i)
            self%hess(nvar+nvar*(nvar-1)/2) = &
               min(max(self%eigv(i),self%hlow),self%hmax)
         endif
      enddo

      if (nvar.ne.self%nvar) then
         thr = thr * 0.1_wp
         cycle get_anc
      endif

      fail = .false.
      exit get_anc

   enddo get_anc

  if (fail) then
      if(pr) write(*,*) 'nvar, selv%nvar',nvar,self%nvar
      return
   end if

   call sort(self%n3,self%nvar,self%hess,self%B)

   self%coord = 0.0_wp
   return
end subroutine generate_anc_packed

!========================================================================================!

pure subroutine sort(nat3,nvar,hess,b)
   implicit none
   integer :: ii,k,j,m,i
   integer, intent(in)    :: nat3,nvar
   real(wp),intent(inout) :: hess(nvar*(nvar+1)/2)
   real(wp),intent(inout) :: b(nat3,nat3)
   real(wp) :: pp,sc1
   real(wp),allocatable   :: edum(:)
   allocate( edum(nvar), source = 0.0_wp )

   do k=1,nvar
      edum(k)=hess(k+k*(k-1)/2)
   enddo
!  sort
   do ii = 2, nvar
      i = ii - 1
      k = i
      pp= edum(i)
      do j = ii, nvar
         if (edum(j) .gt. pp) cycle
         k = j
         pp= edum(j)
      enddo
      if (k .eq. i) cycle
      edum(k) = edum(i)
      edum(i) = pp
      do m=1,nat3
         sc1=b(m,i)
         b(m,i)=b(m,k)
         b(m,k)=sc1
      enddo
   enddo

   do k=1,nvar
      hess(k+k*(k-1)/2)=edum(k)
   enddo
   return
end subroutine sort

!========================================================================================!

subroutine get_cartesian(self,xyz)
   implicit none
   class(optimizer),intent(in) :: self
   integer :: m,i,j,k
   real(wp),intent(out) :: xyz (3,self%n)
   real(wp) :: dum
    external :: dgemv
!>  generate cartesian displacement vector
   xyz = self%xyz
   call dgemv('n',self%n3,self%nvar,1.0_wp,self%B,self%n3,self%coord,1,1.0_wp,xyz,1)
   return
end subroutine get_cartesian

!========================================================================================!

  pure function new_convergence_log(nmax) result(self)
    integer,intent(in) :: nmax
    type(convergence_log) :: self
    self%nlog = 0
    allocate (self%elog(nmax))
    allocate (self%glog(nmax))
  end function new_convergence_log

  pure function get_averaged_energy(self) result(val)
    class(convergence_log),intent(in) :: self
    real(wp) :: eav,val
    integer :: i,j,low
    integer,parameter :: nav = 5

    ! only apply it if sufficient number of points i.e. a "tail" can exist
    ! with the censo blockl = 8 default, this can first be effective in the second
    if (self%nlog .lt. 3 * nav) then
      val = self%elog(self%nlog)
    else
      low = max(1,self%nlog - nav + 1)
      j = 0
      eav = 0
      do i = self%nlog,low,-1
        j = j + 1
        eav = eav + self%elog(i)
      end do
      val = eav / float(j)
    end if

  end function get_averaged_energy

  pure function get_averaged_gradient(self) result(deriv)
    class(convergence_log),intent(in) :: self
    real(wp) :: gav,deriv
    integer :: i,j,low
    integer,parameter :: nav = 5

    ! only apply it if sufficient number of points i.e. a "tail" can exist
    ! with the censo blockl = 8 default, this can first be effective in the second
    if (self%nlog .lt. 3 * nav) then
      deriv = self%glog(self%nlog)
    else
      low = max(1,self%nlog - nav + 1)
      j = 0
      gav = 0
      do i = self%nlog,low,-1
        j = j + 1
        gav = gav + self%glog(i)
      end do
      ! adjust the gradient norm to xtb "conventions" because e.g. a noisy
      ! DCOSMO-RS gradient for large cases can never (even on average)
      ! become lower than the "-opt normal" thresholds
      deriv = gav / float(j) / 2.d0
    end if

  end function get_averaged_gradient

  pure subroutine set_eg_log(self,e,g)
    class(convergence_log),intent(inout) :: self
    real(wp),intent(in) :: e,g
    real(wp),allocatable :: dum(:)
    integer :: k,k2
    k = size(self%elog)
    if (self%nlog >= k) then
      k2 = k + 1
      allocate (dum(k2))
      dum(1:k) = self%elog(1:k)
      call move_alloc(dum,self%elog)
      allocate (dum(k2))
      dum(1:k) = self%glog(1:k)
      call move_alloc(dum,self%glog)
    end if
    self%nlog = self%nlog + 1
    self%elog(self%nlog) = e
    self%glog(self%nlog) = g
  end subroutine set_eg_log


!========================================================================================!
!========================================================================================!
end module optimize_type
