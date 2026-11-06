module hungarian_module
!************************************************************
!* Implementations of 
!* A) The Hungarian (Kuhn-Munkres) Algorithm
!*    in O(n³) time (Edmons & Karp / Tomizawa).
!*
!* B) A Rectengular linear assignment problem algorithm
!*    (LSAP) accodring to 
!*    D.F. Crouse, IEEE Trans. Aerosp. Electron. Syst.,
!*    2016, 52, 1679-1696, doi: 10.1109/TAES.2016.140952
!*
!* Implemented in single precision with a cache to
!* circumvent repeated memory allocation.
!*
!* Also includes some wrappers for standalone use
!************************************************************
  use iso_fortran_env,sp => real32,wp => real64
  implicit none
  private

  public :: hungarian
  interface hungarian
    module procedure hungarian_cached
    module procedure hungarian_wrap_int
    module procedure hungarian_wrap_sp
    module procedure hungarian_wrap_wp
  end interface hungarian

  public :: lsap
  interface lsap
    module procedure lsap_cached
    module procedure lsap_wrap_int
    module procedure lsap_wrap_sp
    module procedure lsap_wrap_wp
  end interface lsap

  real(sp),parameter,private :: inf = huge(1.0_sp) !> Use huge intrinsic for large numbers
  integer,parameter,private  :: infi = huge(1) !> Use huge intrinsic for large numbers

  public :: assignment_cache
  type :: assignment_cache
    integer :: J,W
    real(sp),allocatable :: Cost(:)  !> Cost(J*W), 1D for more efficient memory access
    !> Hungarian algo related
    real(sp),allocatable :: answers(:) !> answers(J)
    integer,allocatable  :: job(:)     !> job(W+1)
    real(sp),allocatable :: ys(:)      !> ys(J)
    real(sp),allocatable :: yt(:)      !> yt(W+1)
    real(sp),allocatable :: Ct(:,:)    !> Ct(W,J)
    real(sp),allocatable :: min_to(:)  !> min_to(W+1)
    integer,allocatable  :: prv(:)     !> prv(W+1)
    logical,allocatable  :: in_Z(:)    !> in_Z(W+1)
    !> LSAP related
    integer,allocatable  :: a(:),b(:)   !> a(J), b(J)
    real(sp),allocatable :: u(:),v(:)   !> u(J), v(W)
    real(sp),allocatable :: shortestPathCosts(:) !> ...(W)
    integer,allocatable  :: path(:),remaining(:) !> path(W), remaining(W)
    integer,allocatable  :: col4row(:),row4col(:) !> col4row(J), row4col(W)
    logical,allocatable  :: SR(:),SC(:)     !> SR(J), SC(W)
  contains
    procedure :: allocate => allocate_assignment_cache
    procedure :: deallocate => deallocate_assignment_cache
  end type assignment_cache

  interface ckmin
  !  module procedure ckmin_int
    module procedure ckmin_sp
  end interface ckmin

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine allocate_assignment_cache(self,J,W,lsapcache)
    implicit none
    integer,intent(in) :: J
    integer,intent(in) :: W
    class(assignment_cache),intent(inout) :: self
    logical,intent(in),optional :: lsapcache
    logical :: yesno
    yesno = .false.
    if (present(lsapcache)) yesno = lsapcache

    !> Store dimensions
    self%J = J
    self%W = W
    if (J > W) then
      error stop 'linear assignment problems require rectengular matrices!'
    end if
    allocate (self%Cost(J*W))

    !> Allocate arrays based on input dimensions & algo type
    if (.not.yesno) then
      !> Hungarian algo cache:
      allocate (self%answers(J))
      allocate (self%job(W+1))
      !> Allocate workspace arrays
      allocate (self%ys(J))
      allocate (self%yt(W+1))
      allocate (self%Ct(W,J))
      allocate (self%min_to(W+1))
      allocate (self%prv(W+1))
      allocate (self%in_Z(W+1))
    else
      !> LSAP cache
      allocate (self%a(J),self%b(J))
      allocate (self%u(J),self%v(W),self%shortestPathCosts(W))
      allocate (self%path(W),self%col4row(J),self%row4col(W))
      allocate (self%SR(J),self%SC(W),self%remaining(W))
    end if
  end subroutine allocate_assignment_cache

  subroutine deallocate_assignment_cache(self)
    implicit none
    class(assignment_cache),intent(inout) :: self
    ! Deallocate arrays if they are allocated
    if (allocated(self%Cost)) deallocate (self%Cost)
    if (allocated(self%answers)) deallocate (self%answers)
    if (allocated(self%job)) deallocate (self%job)
    if (allocated(self%ys)) deallocate (self%ys)
    if (allocated(self%yt)) deallocate (self%yt)
    if (allocated(self%Ct)) deallocate (self%Ct)
    if (allocated(self%min_to)) deallocate (self%min_to)
    if (allocated(self%prv)) deallocate (self%prv)
    if (allocated(self%in_Z)) deallocate (self%in_Z)
    if (allocated(self%a)) deallocate (self%a)
    if (allocated(self%b)) deallocate (self%b)
    if (allocated(self%u)) deallocate (self%u)
    if (allocated(self%v)) deallocate (self%v)
    if (allocated(self%shortestPathCosts)) deallocate (self%shortestPathCosts)
    if (allocated(self%path)) deallocate (self%path)
    if (allocated(self%col4row)) deallocate (self%col4row)
    if (allocated(self%row4col)) deallocate (self%row4col)
    if (allocated(self%SR)) deallocate (self%SR)
    if (allocated(self%SC)) deallocate (self%SC)
    if (allocated(self%remaining)) deallocate (self%remaining)
  end subroutine deallocate_assignment_cache

!========================================================================================!

  !logical function ckmin_int(a,b) result(yesno)
  !  !> Helper function to compute the minimum and update
  !  integer,intent(inout) :: a
  !  integer,intent(in) :: b
  !  yesno = .false.
  !  if (b < a) then
  !    a = b
  !    yesno = .true.
  !  end if
  !end function ckmin_int

  logical function ckmin_sp(a,b) result(yesno)
    !> Helper function to compute the minimum and update
    real(sp),intent(inout) :: a
    real(sp),intent(in) :: b
    yesno = .false.
    if (b < a) then
      a = b
      yesno = .true.
    end if
  end function ckmin_sp

  subroutine hungarian_cached(cache,J,W)
    !****************************************************************
    !* Hungarian algorithm implementation to solve an assignment
    !* problem in O(n³) time.
    !* This implementation refers to a cache, which is created
    !* to avoid repeated memory allocation.
    !* Passing J and W explicitly enables reuse of memory
    !* for smaller sub-problems (i.e. cache%J >= J, W accoridingly)
    !* Unfortunately, this algorithm has problems with
    !* assignments of equal cost.
    !*
    !* Inputs (all within cache, except J and W):
    !*   C(J, W) - Cost matrix of dimensions J-by-W,
    !*             where C(jj, ww) is the cost to assign
    !*             jj-th job to ww-th worker
    !*   J       - Number of jobs
    !*   W       - Number of workers
    !* Outputs (all within cache):
    !*   answers(J) - Vector of length J, where answers(jj) is
    !*                the minimum cost to assign the first jj
    !*                jobs to distinct workers
    !*   job(W+1)   - Vector where job(ww) is the job assigned to
    !*                the ww-th worker (or -1 if no job is assigned)
    !****************************************************************
    integer,intent(in) :: J
    integer,intent(in) :: W
    type(assignment_cache),intent(inout) :: cache
    integer  :: jj_cur,ww_cur,jj,ww_next,ww
    real(sp) :: delta

    !> IMPORTANT: associate to have shorter variable names
    associate (C => cache%Cost, &
              & answers => cache%answers, &
              & job => cache%job, &
              & ys => cache%ys, &
              & yt => cache%yt, &
              & Ct => cache%Ct, &
              & min_to => cache%min_to, &
              & prv => cache%prv, &
              & in_Z => cache%in_Z)

      job = -1
      ys = 0
      yt = 0
      !Ct = transpose(reshape(C,[J,W]))

      do jj_cur = 1,J   !> O(n¹)
        ww_cur = W+1
        job(ww_cur) = jj_cur
        min_to = inf
        prv = -1
        in_Z = .false.

        do while (job(ww_cur) /= -1)  !> O(n¹) -> O(n²)
          in_Z(ww_cur) = .true.
          jj = job(ww_cur)
          delta = inf
          do ww = 1,W !> O(n²) -> O(n³)
            if (.not.in_Z(ww)) then
              !if (ckmin(min_to(ww),Ct(ww,jj)-ys(jj)-yt(ww))) then
              if (ckmin(min_to(ww),C(jj+(ww-1)*J)-ys(jj)-yt(ww))) then
                prv(ww) = ww_cur
              end if
              if (ckmin(delta,min_to(ww))) then
                ww_next = ww
              end if
            end if
          end do

          do ww = 1,W+1
            if (in_Z(ww)) then
              ys(job(ww)) = ys(job(ww))+delta
              yt(ww) = yt(ww)-delta
            else
              min_to(ww) = min_to(ww)-delta
            end if
          end do
          ww_cur = ww_next
        end do

        !> Update assignments along alternating path
        do while (ww_cur /= W+1)
          job(ww_cur) = job(prv(ww_cur))
          ww_cur = prv(ww_cur)
        end do

        answers(jj_cur) = -yt(W+1)
      end do

    end associate
  end subroutine hungarian_cached

!========================================================================================!

  subroutine hungarian_wrap_int(C,J,W,answers,job)
    !*********************************************
    !* Wrapper for integer precision
    !*********************************************
    implicit none
    integer,intent(in) :: J
    integer,intent(in) :: W
    integer,intent(in) :: C(J,W)
    integer,intent(out) :: answers(J)
    integer,intent(out) :: job(W+1)
    type(assignment_cache) :: cache

    call cache%allocate(J,W)
    cache%Cost(1:J*W) = reshape(real(C(1:J,1:W),sp),[J*W])
    call hungarian_cached(cache,J,W)

    answers(1:J) = nint(cache%answers(1:J))
    job(1:W+1) = cache%job(1:W+1)
    call cache%deallocate()
  end subroutine hungarian_wrap_int

  subroutine hungarian_wrap_sp(C,J,W,answers,job)
    !*********************************************
    !* Wrapper for single precision
    !*********************************************
    implicit none
    integer,intent(in) :: J
    integer,intent(in) :: W
    real(sp),intent(in) :: C(J,W)
    real(sp),intent(out) :: answers(J)
    integer,intent(out) :: job(W+1)
    type(assignment_cache) :: cache

    call cache%allocate(J,W)
    cache%Cost(1:J*W) = reshape(C(1:J,1:W),[J*W])
    call hungarian_cached(cache,J,W)

    answers(1:J) = cache%answers(1:J)
    job(1:W+1) = cache%job(1:W+1)
    call cache%deallocate()
  end subroutine hungarian_wrap_sp

  subroutine hungarian_wrap_wp(C,J,W,answers,job)
    !*********************************************
    !* Wrapper for double precision
    !*********************************************
    implicit none
    integer,intent(in) :: J
    integer,intent(in) :: W
    real(wp),intent(in) :: C(J,W)
    real(wp),intent(out) :: answers(J)
    integer,intent(out) :: job(W+1)
    type(assignment_cache) :: cache

    call cache%allocate(J,W)
    cache%Cost(1:J*W) = reshape(real(C(1:J,1:W),sp),[J*W])
    call hungarian_cached(cache,J,W)

    answers(1:J) = real(cache%answers(1:J),wp)
    job(1:W+1) = cache%job(1:W+1)
    call cache%deallocate()
  end subroutine hungarian_wrap_wp

!========================================================================================!
!========================================================================================!

!****************************************************************
!* The following implements an alternative algorithm capable
!* to better handle assignments with equivalent costs
!* The algorithm follows
!*   D.F. Crouse, IEEE Trans. Aerosp. Electron. Syst.,
!*   2016, 52, 1679-1696, doi: 10.1109/TAES.2016.140952
!*
!* The source code is a free Fortran-adaptation
!* of the C++ lsap algorithm in SciPy
!****************************************************************

  function augmenting_path(nr,nc,cost,u,v,path,row4col, &
                          &  shortestPathCosts,i,SR,SC, &
                          &  remaining,minValue) result(sink)
    implicit none
    integer,intent(in) :: nr              !> Number of columns (jobs)
    integer,intent(in) :: nc              !> Number of columns (workers)
    real(sp),intent(in) :: cost(:)        !> Cost matrix (1D, nr*nc length)
    real(sp),intent(inout) :: u(:)        !> Dual variables for rows (jobs)
    real(sp),intent(inout) :: v(:)        !> Dual variables for columns (workers)
    integer,intent(inout) :: path(:)      !> Path array
    integer,intent(inout) :: row4col(:)   !> Array storing which row is assigned to which column
    real(sp),intent(inout) :: shortestPathCosts(:) !> Array for storing shortest path costs
    integer,intent(inout) :: i            !> Current row being processed
    logical,intent(inout) :: SR(:)        !> Boolean array for rows
    logical,intent(inout) :: SC(:)        !> Boolean array for columns
    integer,intent(inout) :: remaining(:) !> Array of remaining columns to be processed
    real(sp),intent(inout) :: minValue    !> Minimum value of the path cost
    integer :: sink                       !> The resulting sink (column) from the augmenting path
    integer :: num_remaining,indx,j,it
    real(sp) :: lowest,r

    minValue = 0.0_sp
    num_remaining = nc

    !> Initialize the remaining array in reverse order
    do it = 1,nc
      remaining(it) = nc-(it-1)
    end do

    !> Initialize SR, SC, and shortestPathCosts
    SR = .false.
    SC = .false.
    shortestPathCosts = inf  !> Set to a very large value

    !> Start finding the shortest augmenting path
    sink = -1
    do while (sink == -1)
      indx = -1
      lowest = inf
      SR(i) = .true.

      do it = 1,num_remaining
        j = remaining(it)
        r = minValue+cost(i+((j-1)*nr))-u(i)-v(j)
        if (r < shortestPathCosts(j)) then
          path(j) = i
          shortestPathCosts(j) = r
        end if

        !> Choose the smallest cost or a new sink node
        if (shortestPathCosts(j) < lowest.or. &
            (shortestPathCosts(j) == lowest.and.row4col(j) == -1)) then
          lowest = shortestPathCosts(j)
          indx = it
        end if
      end do

      minValue = lowest
      if (minValue == inf) then  !> Infeasible cost matrix
        sink = -1
        return
      end if

      j = remaining(indx)
      if (row4col(j) == -1) then
        sink = j
      else
        i = row4col(j)
      end if

      SC(j) = .true.

      remaining(indx) = remaining(num_remaining)
      num_remaining = num_remaining-1
    end do
  end function augmenting_path

  subroutine swap(x,y)
    implicit none
    integer,intent(inout) :: x,y
    integer :: temp
    temp = x
    x = y
    y = temp
  end subroutine swap

  subroutine lsap_cached(lcache,nr,nc,maximize,iostatus)
    implicit none
    type(assignment_cache),intent(inout) :: lcache
    integer,intent(in) :: nr,nc
    logical,intent(in) :: maximize
    integer :: iostatus
    integer :: curRow,currowtmp,i,j,jj,sink
    real(sp) :: minValue
    !> error codes
    integer,parameter :: RECTANGULAR_LSAP_TRANSPOSED = 1
    integer,parameter :: RECTANGULAR_LSAP_INFEASIBLE = 2

    !> use associates to offload allocation outside the routine
    associate (cost => lcache%Cost, &
             & a => lcache%a,b => lcache%b, &
             & u => lcache%u,v => lcache%v, &
             & shortestPathCosts => lcache%shortestPathCosts, &
             & path => lcache%path,col4row => lcache%col4row,row4col => lcache%row4col, &
             & remaining => lcache%remaining,SR => lcache%SR,SC => lcache%SC)

      !> Handle trivial inputs
      if (nr == 1.or.nc == 1) then
        a(1) = 1
        b(1) = 1
        iostatus = 0
        return
      end if

      !> Determine if we need to transpose the matrix
      !> Let the user handle that outside the call
      if (nc < nr) then
        iostatus = RECTANGULAR_LSAP_TRANSPOSED
        return
      end if

      !> Negate the cost matrix for maximization
      if (maximize) then
        cost = -cost
      end if

      !> Initialize
      u(:) = 0.0_sp
      v(:) = 0.0_sp
      col4row(:) = -1
      row4col(:) = -1
      path(:) = -1

      !> Iteratively build the solution
      do curRow = 1,nr
        curRowtmp = curRow
        !> Call augmenting_path routine
        sink = augmenting_path(nr,nc,cost,u,v,path,row4col,  &
                              & shortestPathCosts,curRowtmp, &
                              & SR,SC,remaining,minValue)
        if (sink < 0) then
          iostatus = RECTANGULAR_LSAP_INFEASIBLE
          return
        end if

        !> Update dual variables
        u(curRow) = u(curRow)+minValue
        do i = 1,nr
          if (SR(i).and.i /= curRow) then
            u(i) = u(i)+minValue-shortestPathCosts(col4row(i))
          end if
        end do

        do j = 1,nc
          if (SC(j)) then
            v(j) = v(j)-minValue+shortestPathCosts(j)
          end if
        end do

        !> Augment previous solution
        j = sink
        do jj=1,nc+1  !> avoid infinite loop
          i = path(j)
          row4col(j) = i
          call swap(col4row(i),j)
          if (i == curRow) exit
        end do
      end do

      !> Finalize the assignment
      do i = 1,nr
        a(i) = i
        b(i) = col4row(i)
      end do

      iostatus = 0
    end associate
  end subroutine lsap_cached

!========================================================================================!

  subroutine lsap_wrap_int(C,J,W,a,b)
    !*********************************************
    !* Wrapper for integer precision
    !*********************************************
    implicit none
    integer,intent(in) :: J
    integer,intent(in) :: W
    integer,intent(in) :: C(J,W)
    integer,intent(out),allocatable :: a(:)
    integer,intent(out),allocatable :: b(:)
    type(assignment_cache) :: cache
    integer :: io

    call cache%allocate(J,W,.true.)
    cache%Cost(1:J*W) = reshape(real(C(1:J,1:W),sp), [J*W])
    call lsap_cached(cache,J,W,.false.,io)

    allocate(a(J), b(J))
    a(1:J) = cache%a(1:J)
    b(1:J) = cache%b(1:J)
    call cache%deallocate()
  end subroutine lsap_wrap_int


  subroutine lsap_wrap_sp(C,J,W,a,b)
    !*********************************************
    !* Wrapper for single precision
    !*********************************************
    implicit none
    integer,intent(in) :: J
    integer,intent(in) :: W
    real(sp),intent(in) :: C(J,W)
    integer,intent(out),allocatable :: a(:)
    integer,intent(out),allocatable :: b(:)
    type(assignment_cache) :: cache
    integer :: io

    call cache%allocate(J,W,.true.)
    cache%Cost(1:J*W) = reshape(C(1:J,1:W), [J*W])
    call lsap_cached(cache,J,W,.false.,io)

    allocate(a(J), b(J))
    a(1:J) = cache%a(1:J)
    b(1:J) = cache%b(1:J)
    call cache%deallocate()
  end subroutine lsap_wrap_sp


  subroutine lsap_wrap_wp(C,J,W,a,b)
    !*********************************************
    !* Wrapper for double precision
    !*********************************************
    implicit none
    integer,intent(in) :: J
    integer,intent(in) :: W
    real(wp),intent(in) :: C(J,W)
    integer,intent(out),allocatable :: a(:)
    integer,intent(out),allocatable :: b(:)
    type(assignment_cache) :: cache
    integer :: io

    call cache%allocate(J,W,.true.)
    cache%Cost(1:J*W) = reshape(real(C(1:J,1:W),sp), [J*W])
    call lsap_cached(cache,J,W,.false.,io)

    allocate(a(J), b(J))
    a(1:J) = cache%a(1:J)
    b(1:J) = cache%b(1:J)
    call cache%deallocate()
  end subroutine lsap_wrap_wp

!========================================================================================!
!========================================================================================!
end module hungarian_module
