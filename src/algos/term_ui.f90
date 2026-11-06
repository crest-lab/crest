module term_ui
  use crest_parameters
  use iso_c_binding,only:c_int
  implicit none
  private

  interface
    integer(c_int) function c_isatty(fd) bind(c,name="isatty")
      use iso_c_binding,only:c_int
      integer(c_int),value :: fd
    end function c_isatty
  end interface

  !> Public API
  public :: ansi_enabled,set_ansi_enabled
  public :: fg,bg,style,reset,strip_ansi
  public :: printc,eprintc
  public :: progress_state
  public :: progress_init,progress_update,progress_finish
  public :: clear_line

  !> ============================================================
  !> Configuration
  !> ============================================================
  logical :: ansi_enabled = .true.

  integer,parameter :: C_RESET = -1

  !> Basic 8-color codes (foreground base 30..37, background base 40..47)
  integer,parameter,public :: &
   &  BLACK = 0,RED = 1,GREEN = 2,YELLOW = 3, &
   &  BLUE = 4,MAGENTA = 5,CYAN = 6,WHITE = 7

  !> Style SGR codes
  integer,parameter,public :: &
   &  S_BOLD = 1, &
   &  S_DIM = 2, &
   &  S_UNDERLINE = 4, &
   &  S_BLINK = 5, &
   &  S_REVERSE = 7, &
   &  S_HIDDEN = 8

  !> ============================================================
  !> Progress state
  !> ============================================================
  type :: progress_state
    real(wp) :: t0 = 0.0_wp
    integer      :: width = 40
    character(:),allocatable :: prefix
    character(:),allocatable :: suffix
    character(:),allocatable :: fill_char
    character(:),allocatable :: empty_char
    character(:),allocatable :: left_cap
    character(:),allocatable :: right_cap
    integer(int64) :: last_draw_ms = -huge(0_int64)
    integer(int64) :: min_interval_ms = 50_int64  !> throttle redraw
    logical :: show_time = .true.
    logical :: show_eta = .true.
    logical :: started = .false.
    logical :: tty = .false.        !> auto-set by progress_init
    real(wp) :: increment = 10.0_wp !> % step for non-TTY mode
    real(wp) :: barrier = 0.0_wp   !> current % threshold for non-TTY mode
  end type progress_state

!========================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================!

  !> ============================================================
  !> Time helper (seconds as wp) - based on system_clock
  !> ============================================================
  function now_seconds() result(t)
    real(wp) :: t
    integer(int64) :: count,rate
    call system_clock(count=count,count_rate=rate)
    if (rate > 0_int64) then
      t = real(count,wp)/real(rate,wp)
    else
      t = 0.0_wp
    end if
  end function now_seconds

  function now_millis_i64() result(ms)
    integer(int64) :: ms
    integer(int64) :: count,rate
    call system_clock(count=count,count_rate=rate)
    if (rate > 0_int64) then
      ms = int((real(count,wp)*1000.0_wp)/real(rate,wp),int64)
    else
      ms = 0_int64
    end if
  end function now_millis_i64

  subroutine set_ansi_enabled(flag)
    logical,intent(in) :: flag
    ansi_enabled = flag
  end subroutine set_ansi_enabled

  !> ============================================================
  !> ANSI builders
  !> ============================================================
  pure function sgr(code) result(s)
    integer,intent(in) :: code
    character(:),allocatable :: s
    if (.not.ansi_enabled) then
      s = ""
    else
      s = achar(27)//"["//itoa(code)//"m"
    end if
  end function sgr

  pure function reset() result(s)
    character(:),allocatable :: s
    s = sgr(0)
  end function reset

  pure function fg(color,bright) result(s)
    integer,intent(in) :: color
    logical,intent(in),optional :: bright
    character(:),allocatable :: s
    integer :: base
    logical :: b
    if (.not.ansi_enabled) then
      s = ""
      return
    end if
    b = .false.; if (present(bright)) b = bright
    if (color == C_RESET) then
      s = achar(27)//"[39m"
    else
      base = 30+max(0,min(7,color))
      if (b) then
        !> "Bright" variant using 90..97 (widely supported)
        base = 90+max(0,min(7,color))
      end if
      s = achar(27)//"["//itoa(base)//"m"
    end if
  end function fg

  pure function bg(color,bright) result(s)
    integer,intent(in) :: color
    logical,intent(in),optional :: bright
    character(:),allocatable :: s
    integer :: base
    logical :: b
    if (.not.ansi_enabled) then
      s = ""
      return
    end if
    b = .false.; if (present(bright)) b = bright
    if (color == C_RESET) then
      s = achar(27)//"[49m"
    else
      base = 40+max(0,min(7,color))
      if (b) then
        base = 100+max(0,min(7,color))
      end if
      s = achar(27)//"["//itoa(base)//"m"
    end if
  end function bg

  pure function style(code) result(s)
    integer,intent(in) :: code
    character(:),allocatable :: s
    s = sgr(code)
  end function style

  !> ============================================================
  !> Printing helpers
  !> ============================================================
  subroutine printc(msg,unit,advance)
    character(len=*),intent(in) :: msg
    integer,intent(in),optional :: unit
    logical,intent(in),optional :: advance
    integer :: u
    logical :: adv
    u = stdout; if (present(unit)) u = unit
    adv = .true.; if (present(advance)) adv = advance
    if (adv) then
      write (u,'(a)') msg
    else
      write (u,'(a)',advance='no') msg
    end if
  end subroutine printc

  subroutine eprintc(msg,advance)
    character(len=*),intent(in) :: msg
    logical,intent(in),optional :: advance
    call printc(msg,unit=stderr,advance=advance)
  end subroutine eprintc

  !> Clears current line and returns carriage to start (no newline).
  pure function clear_line() result(s)
    character(:),allocatable :: s
    if (.not.ansi_enabled) then
      s = achar(13)  !> CR only
    else
      !> CR + "erase line" (EL2) is common; some use EL0. We'll do EL2.
      s = achar(13)//achar(27)//"[2K"
    end if
  end function clear_line

  !> ============================================================
  !> Progress bar API
  !> ============================================================
  subroutine progress_init(ps,total,width,prefix,suffix,min_interval_ms,show_time,show_eta, &
                           fill_char,empty_char,left_cap,right_cap)
    !***********************************************************
    !* Initialise a progress_state object.
    !* Automatically detects whether stdout is a TTY and selects
    !* either the fancy ANSI bar or the sequential |>x% mode.
    !* total: total number of steps (used to calibrate the non-TTY
    !*        print interval).
    !***********************************************************
    implicit none
    type(progress_state),intent(inout) :: ps
    integer,intent(in),optional :: total
    integer,intent(in),optional :: width
    character(len=*),intent(in),optional :: prefix,suffix
    integer(int64),intent(in),optional :: min_interval_ms
    logical,intent(in),optional :: show_time,show_eta
    character(len=*),intent(in),optional :: fill_char,empty_char,left_cap,right_cap

    ps%t0 = now_seconds()
    ps%started = .true.

    if (present(width)) ps%width = max(5,width)
    if (present(min_interval_ms)) ps%min_interval_ms = max(0_int64,min_interval_ms)
    if (present(show_time)) ps%show_time = show_time
    if (present(show_eta)) ps%show_eta = show_eta

    if (present(prefix)) then
      ps%prefix = prefix
    else
      ps%prefix = ""
    end if

    if (present(suffix)) then
      ps%suffix = suffix
    else
      ps%suffix = ""
    end if

    if (present(fill_char)) then
      ps%fill_char = fill_char
    else
      ps%fill_char = "█"
    end if

    if (present(empty_char)) then
      ps%empty_char = empty_char
    else
      ps%empty_char = "░"
    end if

    if (present(left_cap)) then
      ps%left_cap = left_cap
    else
      ps%left_cap = "["
      ps%left_cap = "┃"
    end if

    if (present(right_cap)) then
      ps%right_cap = right_cap
    else
      ps%right_cap = "]"
      ps%right_cap = "┃"
    end if

    ps%last_draw_ms = -huge(0_int64)

    ! ── TTY detection & non-TTY state ────────────────────────────
    ps%tty = (c_isatty(1_c_int) /= 0)
    ps%barrier = 0.0_wp
    ps%increment = 10.0_wp
    if (.not.ps%tty .and. present(total)) then
      if (total > 20000) ps%increment = 1.0_wp
      if (total > 10000 .and. total <= 20000) ps%increment = 2.5_wp
      if (total > 5000  .and. total <= 10000) ps%increment = 5.0_wp
      if (total > 1000  .and. total <= 5000)  ps%increment = 7.5_wp
    end if
  end subroutine progress_init

  subroutine progress_update(ps,curr,tot,unit,force)
    !***********************************************************
    !* Update the progress display.
    !* TTY: renders the in-place ANSI bar.
    !* non-TTY: prints sequential |>x% markers.
    !***********************************************************
    implicit none
    type(progress_state),intent(inout) :: ps
    integer,intent(in) :: curr,tot
    integer,intent(in),optional :: unit
    logical,intent(in),optional :: force

    integer :: u
    integer(int64) :: current,total
    logical :: do_force
    integer(int64) :: ms,elapsed_ms
    real(wp) :: frac,elapsed_s,rate,eta_s
    integer :: filled,nfilled
    character(:),allocatable :: line,pct,time_str,eta_str
    character(:),allocatable :: bar
    real(wp) :: pct_r
    character(len=5) :: atmp

    if (.not.ps%started) call progress_init(ps)

    current = int(curr,int64)
    total = int(tot,int64)
    u = stdout; if (present(unit)) u = unit
    do_force = .false.; if (present(force)) do_force = force

    if (ps%tty) then
      ! ── Fancy in-place ANSI bar ───────────────────────────────────
      if (current == total) do_force = .true.

      ms = now_millis_i64()
      if (.not.do_force) then
        if (ps%last_draw_ms >= 0_int64) then
          if (ms-ps%last_draw_ms < ps%min_interval_ms) return
        end if
      end if
      ps%last_draw_ms = ms

      if (total <= 0_int64) then
        frac = 0.0_wp
      else
        frac = real(max(0_int64,min(current,total)),wp)/real(total,wp)
      end if

      filled = int(frac*real(ps%width,wp))
      filled = max(0,min(ps%width,filled))
      nfilled = max(0,(ps%width-filled))

      pct = fmt_percent(frac)

      elapsed_s = now_seconds()-ps%t0
      elapsed_ms = int(elapsed_s*1000.0_wp,int64)

      if (ps%show_time) then
        time_str = " "//dim_text("("//fmt_hms(elapsed_ms)//")")
      else
        time_str = ""
      end if

      if (ps%show_eta.and.current > 0_int64.and.total > 0_int64.and.current < total) then
        rate = real(current,wp)/max(1.0e-12_wp,elapsed_s)
        eta_s = real(total-current,wp)/max(1.0e-12_wp,rate)
        eta_str = " "//dim_text("ETA "//fmt_hms(int(eta_s*1000.0_wp,int64)))
      else if (ps%show_eta.and.total > 0_int64.and.current >= total) then
        eta_str = " "//dim_text("ETA 00:00")
      else
        eta_str = ""
      end if

      if (ansi_enabled) then
        if (frac >= 0.999_wp) then
          bar = ps%left_cap// &
            & fg(GREEN,bright=.true.)//repeat(ps%fill_char,filled)//reset()// &
            & repeat(ps%empty_char,nfilled)//ps%right_cap
        else
          bar = ps%left_cap// &
            & fg(GREEN,bright=.false.)//repeat(ps%fill_char,filled)//reset()// &
            & fg(YELLOW,bright=.false.)//repeat(ps%empty_char,nfilled)//reset()// &
            & ps%right_cap
        end if
      else
        bar = ps%left_cap//repeat(ps%fill_char,filled)// &
              & repeat(ps%empty_char,nfilled)//ps%right_cap
      end if

      line = clear_line()//ps%prefix//bar//" "//pct//ps%suffix//eta_str//time_str

      call printc(line,unit=u,advance=.false.)
      call flush_unit(u)

    else
      ! ── non-TTY: sequential |>x% printout ────────────────────────
      if (current <= 0_int64) return
      if (total <= 0_int64) return
      pct_r = real(current,wp)/real(total,wp)*100.0_wp
      if (pct_r >= ps%barrier) then
        write (atmp,'(f5.1)') pct_r
        write (u,'(1x,a)',advance='no') '|>'//trim(adjustl(atmp))//'%'
        ps%barrier = min(ps%barrier+ps%increment,100.0_wp)
        call flush_unit(u)
      end if

    end if
  end subroutine progress_update

  subroutine progress_finish(ps,unit,newline)
    !***********************************************************
    !* Finalise the progress display.
    !* TTY: prints a terminating newline.
    !* non-TTY: prints a newline followed by " done." to close the
    !*          inline |>x% sequence.
    !***********************************************************
    implicit none
    type(progress_state),intent(inout) :: ps
    integer,intent(in),optional :: unit
    logical,intent(in),optional :: newline
    integer :: u
    logical :: nl

    u = stdout; if (present(unit)) u = unit
    nl = .true.; if (present(newline)) nl = newline

    if (ps%tty) then
      if (nl) call printc("",unit=u,advance=.true.)
    else
      if (nl) write (u,'(/,1x,a)') 'done.'
    end if
    call flush_unit(u)
    ps%started = .false.
  end subroutine progress_finish

  !> ============================================================
  !> Utilities: formatting, flushing, ANSI stripping
  !> ============================================================
  subroutine flush_unit(u)
    integer,intent(in) :: u
    flush (u)
  end subroutine flush_unit

  pure function itoa(i) result(s)
    integer,intent(in) :: i
    character(:),allocatable :: s
    character(len=64) :: buf
    write (buf,'(i0)') i
    s = trim(buf)
  end function itoa

  pure function itoa_i64(i) result(s)
    integer(int64),intent(in) :: i
    character(:),allocatable :: s
    character(len=64) :: buf
    write (buf,'(i0)') i
    s = trim(buf)
  end function itoa_i64

  pure function fmt_percent(frac) result(s)
    real(wp),intent(in) :: frac
    character(:),allocatable :: s
    integer :: p
    real(wp) :: rp
    character(8) :: buf
    rp = min(100.0_wp*max(0.0_wp,min(1.0_wp,frac))+0.5_wp,100.0_wp)
    write (buf,'(f5.1,a)') rp,"%"
    s = adjustl(buf)
  end function fmt_percent

  pure function fmt_hms(ms) result(s)
    integer(int64),intent(in) :: ms
    character(:),allocatable :: s
    integer(int64) :: t,hh,mm,ss
    character(32) :: buf

    t = max(0_int64,ms/1000_int64)
    hh = t/3600_int64
    mm = (t-hh*3600_int64)/60_int64
    ss = t-hh*3600_int64-mm*60_int64

    if (hh > 0_int64) then
      write (buf,'(i0,":",i2.2,":",i2.2)') hh,int(mm),int(ss)
    else
      write (buf,'(i2.2,":",i2.2)') int(mm),int(ss)
    end if
    s = trim(buf)
  end function fmt_hms

  pure function dim_text(t) result(s)
    character(len=*),intent(in) :: t
    character(:),allocatable :: s
    if (.not.ansi_enabled) then
      s = t
    else
      s = style(S_DIM)//t//reset()
    end if
  end function dim_text

  pure function strip_ansi(s_in) result(s_out)
    character(len=*),intent(in) :: s_in
    character(:),allocatable :: s_out
    integer :: i,n
    character :: c
    logical :: in_esc
    in_esc = .false.
    s_out = ""
    n = len_trim(s_in)

    i = 1
    do while (i <= n)
      c = s_in(i:i)
      if (.not.in_esc) then
        if (c == achar(27)) then
          in_esc = .true.
        else
          s_out = s_out//c
        end if
      else
        !> We are inside ESC[ ... m ; consume until 'm' or end
        if (c == "m") then
          in_esc = .false.
        end if
      end if
      i = i+1
    end do
  end function strip_ansi

!=============================================================================!
!#############################################################################!
!=============================================================================!
end module term_ui

