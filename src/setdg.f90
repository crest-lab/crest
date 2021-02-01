 module setdgmod
       interface setdg
       module procedure setdg_subr_real
       module procedure setdg_subr_int
       module procedure setdg_subr_string
       module procedure setdg_subr_metadyn
       end interface setdg
 
       interface setdg_block
       module procedure setdg_block_plain
       module procedure setdg_block_key
       end interface setdg_block

       contains
!-----------------------------------------------------------------------------
! split a set-block line and get the keyword
!-----------------------------------------------------------------------------
subroutine get_set_arg(line,argument)
        implicit none
        character(len=*) :: line
        character(len=*) :: argument

        character(len=256) :: dummy
        character(len=1) :: digit
        integer :: i,j


        line=adjustl(line)

        do i=1,len(trim(line))
           digit=line(i:i)
           if((digit==' ').or.(digit==char(9)) &
           & .or. (digit=='='))then
           j=i-1
           argument=adjustl(line(1:j))
           exit
           endif
        enddo

end subroutine get_set_arg

!-----------------------------------------------------------------------------
! sets and modifies the $set block in a coord file
! if the keyword exists it is overwritten
!-----------------------------------------------------------------------------
subroutine setdg_subr_real(fnam,dg,tmp2)
      implicit none
      character(len=*),intent(in) :: fnam
      character(len=*),intent(in) :: dg
      character(len=512) :: dg2
      character(len=512) :: tmp
      character(len=64)  :: argum
      integer :: ich,och,io
      integer :: i,j,k,nn
      real*8 :: tmp2

      write(dg2,'(a,f14.4)')trim(dg),tmp2

      open(newunit=ich,file=trim(fnam))
      open(newunit=och,file='.setdgtmp')
      do
        read(ich,'(a)',iostat=io)tmp
        if(io > 0)then
           write(och,'(''$end'')')
           exit
        endif
        if(index(tmp,'$set').ne.0)then
           backspace(ich)
           exit
        else
          write(och,'(a)')trim(tmp)
        endif
        if(index(tmp,'$end').ne.0)exit
      enddo

      outer: do
        read(ich,'(a)',iostat=io)tmp
        if(io < 0)then                   !done if no $set-block is present
           write(och,'(''$set'')')
           write(och,'(a)') trim(dg2)
           write(och,'(''$end'')')
           exit
        endif
        if(trim(tmp) == '') cycle outer  !skip blank lines
        write(och,'(a)')trim(tmp)
        if(index(tmp,'$set') /= 0) then  !done if a $setblock is present
           inner: do
             read(ich,'(a)',iostat=io)tmp
             if((io < 0).or.(index(tmp,'$end') /= 0))then !terminate set block
                write(och,'(a)') trim(dg2)
                write(och,'(''$end'')')
                exit outer              !important exit of outer loop
             endif
             if(trim(tmp) == '') cycle inner !skip blank lines
             call get_set_arg(trim(tmp),argum)     
             if(index(argum,trim(dg)) /= 0)then  !skip a write if the argument "dg" is already present
                cycle inner
             else
                write(och,'(a)')trim(tmp)
             endif
           enddo inner
        endif
      enddo outer

      close(ich,status='delete')
      close(och)

      call rename('.setdgtmp',trim(fnam))

end subroutine setdg_subr_real

!-----------------------------------------------------------------------------
! sets and modifies the $set block in a coord file
! if the keyword exists it is overwritten
!-----------------------------------------------------------------------------
subroutine setdg_subr_int(fnam,dg,tmp2)
      implicit none
      character(len=*),intent(in) :: fnam
      character(len=*),intent(in) :: dg
      character(len=512) :: dg2
      character(len=512) :: tmp
      character(len=64)  :: argum
      integer :: ich,och,io
      integer :: i,j,k,nn
      integer :: tmp2

      write(dg2,'(a,i10)')trim(dg),tmp2

      open(newunit=ich,file=trim(fnam))
      open(newunit=och,file='.setdgtmp')
      do
        read(ich,'(a)',iostat=io)tmp
        if(io < 0)then
           write(och,'(''$end'')')
           exit
        endif
        if(index(tmp,'$set').ne.0)then
           backspace(ich)
           exit
        else
          write(och,'(a)')trim(tmp)
        endif
        if(index(tmp,'$end').ne.0)exit
      enddo

      outer: do
        read(ich,'(a)',iostat=io)tmp
        if(io < 0)then                   !done if no $set-block is present
           write(och,'(''$set'')')
           write(och,'(a)') trim(dg2)
           write(och,'(''$end'')')
           exit
        endif
        if(trim(tmp) == '') cycle outer  !skip blank lines
        write(och,'(a)')trim(tmp)
        if(index(tmp,'$set') /= 0) then  !done if a $setblock is present
           inner: do
             read(ich,'(a)',iostat=io)tmp
             if((io < 0).or.(index(tmp,'$end') /= 0))then !terminate set block
                write(och,'(a)') trim(dg2)
                write(och,'(''$end'')')
                exit outer              !important exit of outer loop
             endif
             if(trim(tmp) == '') cycle inner !skip blank lines
             call get_set_arg(trim(tmp),argum)
             if(index(argum,trim(dg)) /= 0)then  !skip a write if the argument "dg" is already present
                cycle inner
             else
                write(och,'(a)')trim(tmp)
             endif
           enddo inner
        endif
      enddo outer

      close(ich,status='delete')
      close(och)

      call rename('.setdgtmp',trim(fnam))

end subroutine setdg_subr_int

!-----------------------------------------------------------------------------
! sets and modifies the $set block in a coord file
! if the keyword exists it is overwritten
!-----------------------------------------------------------------------------
subroutine setdg_subr_string(fnam,dg,tmp2)
      implicit none
      character(len=*),intent(in) :: fnam
      character(len=*),intent(in) :: dg
      character(len=*),intent(in) :: tmp2
      character(len=512) :: dg2
      character(len=512) :: tmp
      character(len=64)  :: argum
      integer :: ich,och,io
      integer :: i,j,k,nn

      write(dg2,'(a,4x,a)')trim(dg),trim(tmp2)

      open(newunit=ich,file=trim(fnam))
      open(newunit=och,file='.setdgtmp')
      do
        read(ich,'(a)',iostat=io)tmp
        if(io < 0)then
           write(och,'(''$end'')')
           exit
        endif
        if(index(tmp,'$set').ne.0)then
           backspace(ich)
           exit
        else
          write(och,'(a)')trim(tmp)
        endif
        if(index(tmp,'$end').ne.0)exit
      enddo

      outer: do
        read(ich,'(a)',iostat=io)tmp
        if(io < 0)then                   !done if no $set-block is present
           write(och,'(''$set'')')
           write(och,'(a)') trim(dg2)
           write(och,'(''$end'')')
           exit
        endif
        if(trim(tmp) == '') cycle outer  !skip blank lines
        write(och,'(a)')trim(tmp)
        if(index(tmp,'$set') /= 0) then  !done if a $setblock is present
           inner: do
             read(ich,'(a)',iostat=io)tmp
             if((io < 0).or.(index(tmp,'$end') /= 0))then !terminate set block
                write(och,'(a)') trim(dg2)
                write(och,'(''$end'')')
                exit outer              !important exit of outer loop
             endif
             if(trim(tmp) == '') cycle inner !skip blank lines
             call get_set_arg(trim(tmp),argum)
             if(index(argum,trim(dg)) /= 0)then  !skip a write if the argument "dg" is already present
                cycle inner
             else
                write(och,'(a)')trim(tmp)
             endif
           enddo inner
        endif
      enddo outer

      close(ich,status='delete')
      close(och)

      call rename('.setdgtmp',trim(fnam))

end subroutine setdg_subr_string
 
!-----------------------------------------------------------------------------
! sets and modifies the $set block in a coord file
! if the keyword exists it is overwritten
! The same as setdg_subr_string, but looking for the keyword $set instead of $end
!-----------------------------------------------------------------------------
subroutine setdg_constr(fnam,dg,tmp2)
      implicit none
      character(len=*),intent(in) :: fnam
      character(len=*),intent(in) :: dg
      character(len=*),intent(in) :: tmp2
      character(len=512) :: dg2
      character(len=512) :: tmp
      character(len=64)  :: argum
      integer :: ich,och,io
      integer :: i,j,k,nn

      write(dg2,'(a,4x,a)')trim(dg),trim(tmp2)

      open(newunit=ich,file=trim(fnam))
      open(newunit=och,file='.setdgtmp')
      do
        read(ich,'(a)',iostat=io)tmp
        if(io < 0)then
          write(och,'(a)')'$set'
          exit
        endif
        write(och,'(a)')trim(tmp)
        if(index(tmp,'$set').ne.0)exit
      enddo

      outer: do
        read(ich,'(a)',iostat=io)tmp
        if((io < 0).or.(index(tmp,'$end') /= 0))then    !done if EOF is reached or the argument was not present
           write(och,'(a)') trim(dg2)
           write(och,'(''$end'')')
           exit outer
        endif
        if(trim(tmp) == '') cycle outer  !skip blank lines
        call get_set_arg(trim(tmp),argum)
        if(index(argum,trim(dg)) /= 0)then  !skip a write if the argument "dg" is already present
           cycle outer
        else
           write(och,'(a)')trim(tmp)
        endif
      enddo outer

      close(ich,status='delete')
      close(och)

      call rename('.setdgtmp',trim(fnam))

end subroutine setdg_constr

!------------------------------------------------------------------------------------------
! set the metadyn argument
!------------------------------------------------------------------------------------------
subroutine setdg_subr_metadyn(fnam,dg,tmp2,tmp3,tmp4)
      implicit none
      character(len=*),intent(in) :: fnam
      character(len=*),intent(in) :: dg
      real*8,intent(in) :: tmp2
      real*8,intent(in) :: tmp3
      integer,intent(in) :: tmp4
      character(len=512) :: dg2
      character(len=512) :: tmp
      character(len=64)  :: argum
      integer :: ich,och,io
      integer :: i,j,k,nn

      write(dg2,'(a,1x,f10.5,1x,f10.5,3x,i0)')trim(dg),tmp2,tmp3,tmp4

      open(newunit=ich,file=trim(fnam))
      open(newunit=och,file='.setdgtmp')
      do
        read(ich,'(a)',iostat=io)tmp
        if(io < 0)then
           write(och,'(''$end'')')
           exit
        endif
        if(index(tmp,'$set').ne.0)then
           backspace(ich)
           exit
        else
          write(och,'(a)')trim(tmp)
        endif
        if(index(tmp,'$end').ne.0)exit
      enddo

      outer: do
        read(ich,'(a)',iostat=io)tmp
        if(io < 0)then                   !done if no $set-block is present
           write(och,'(''$set'')')
           write(och,'(a)') trim(dg2)
           write(och,'(''$end'')')
           exit
        endif
        if(trim(tmp) == '') cycle outer  !skip blank lines
        write(och,'(a)')trim(tmp)
        if(index(tmp,'$set') /= 0) then  !done if a $setblock is present
           inner: do
             read(ich,'(a)',iostat=io)tmp
             if((io < 0).or.(index(tmp,'$end') /= 0))then !terminate set block
                write(och,'(a)') trim(dg2)
                write(och,'(''$end'')')
                exit outer              !important exit of outer loop
             endif
             if(trim(tmp) == '') cycle inner !skip blank lines
             call get_set_arg(trim(tmp),argum)
             if(index(argum,trim(dg)) /= 0)then  !skip a write if the argument "dg" is already present
                cycle inner
             else
                write(och,'(a)')trim(tmp)
             endif
           enddo inner
        endif
      enddo outer

      close(ich,status='delete')
      close(och)

      call rename('.setdgtmp',trim(fnam))

end subroutine setdg_subr_metadyn

!-----------------------------------------------------------------------------------------
! set the optlev argument
!-----------------------------------------------------------------------------------------
subroutine setdg_optlev(fnam,optlev)
      implicit none
      character(len=*),intent(in) :: fnam
      real*8,intent(in) :: optlev
      character(len=512) :: dg2
      character(len=512) :: tmp
      character(len=64)  :: argum
      integer :: ich,och,io
      integer :: i,j,k,nn

      write(dg2,'(''optlev'',i5)')nint(optlev)

      open(newunit=ich,file=trim(fnam))
      open(newunit=och,file='.setdgtmp')
      do
        read(ich,'(a)',iostat=io)tmp
        if(io < 0)then
           write(och,'(''$end'')')
           exit
        endif
        if(index(tmp,'$set').ne.0)then
           backspace(ich)
           exit
        else
          write(och,'(a)')trim(tmp)
        endif
        if(index(tmp,'$end').ne.0)exit
      enddo

      outer: do
        read(ich,'(a)',iostat=io)tmp
        if(io < 0)then                   !done if no $set-block is present
           write(och,'(''$set'')')
           write(och,'(a)') trim(dg2)
           write(och,'(''$end'')')
           exit
        endif
        if(trim(tmp) == '') cycle outer  !skip blank lines
        write(och,'(a)')trim(tmp)
        if(index(tmp,'$set') /= 0) then  !done if a $setblock is present
           inner: do
             read(ich,'(a)',iostat=io)tmp
             if((io < 0).or.(index(tmp,'$end') /= 0))then !terminate set block
                write(och,'(a)') trim(dg2)
                write(och,'(''$end'')')
                exit outer              !important exit of outer loop
             endif
             if(trim(tmp) == '') cycle inner !skip blank lines
             call get_set_arg(trim(tmp),argum)
             if(index(argum,'optlev') /= 0)then  !skip a write if the argument "dg" is already present
                cycle inner
             else
                write(och,'(a)')trim(tmp)
             endif
           enddo inner
        endif
      enddo outer

      close(ich,status='delete')
      close(och)

      call rename('.setdgtmp',trim(fnam))

end subroutine setdg_optlev

!-----------------------------------------------------------------------------------------
! add a line "dg"  after the coord block, but before the the $set block (if present)
!-----------------------------------------------------------------------------------------
subroutine setdg_block_plain(fnam,dg)
      implicit none
      character(len=*),intent(in) :: fnam
      character(len=*),intent(in) :: dg
      character(len=512) :: tmp
      integer :: ich,och,io
      integer :: i,j,k,nn

      open(newunit=ich,file=trim(fnam))
      open(newunit=och,file='.setdgtmp')
      do
        read(ich,'(a)',iostat=io)tmp
        if(io < 0)then
           write(och,'(''$end'')')
           exit
        endif
        write(och,'(a)')trim(tmp)
        if(index(tmp,'$end').ne.0)exit
      enddo

      outer: do
        read(ich,'(a)',iostat=io)tmp
        if(io < 0)then                   !done if no $set-block is present
           write(och,'(a)')trim(dg)      !just add the line dg at the end of the file
           exit
        endif
        if(trim(tmp) == '') cycle outer  !skip blank lines
        if(index(tmp,'$set') /= 0) then  !done if the $set-block is present
           write(och,'(a)')trim(dg)      !add the line
           write(och,'(a)')'$set'        !then write until the end of the file 
           inner: do
             read(ich,'(a)',iostat=io)tmp
             if(io < 0) exit outer
             write(och,'(a)')trim(tmp)
           enddo inner
       endif
       write(och,'(a)')trim(tmp)
      enddo outer

      close(ich,status='delete')
      close(och)

      call rename('.setdgtmp',trim(fnam))

end subroutine setdg_block_plain

!-----------------------------------------------------------------------------------------
! add a line "dg"  after the coord block, but before the the $set block (if present)
! also looks checks if there is already a corresponding $-block before the $setblock
!-----------------------------------------------------------------------------------------
subroutine setdg_block_key(fnam,key,dg)
      implicit none
      character(len=*),intent(in) :: fnam
      character(len=*),intent(in) :: key
      character(len=*),intent(in) :: dg
      character(len=512) :: tmp
      integer :: ich,och,io
      integer :: i,j,k,nn

      open(newunit=ich,file=trim(fnam))
      open(newunit=och,file='.setdgtmp')
      do
        read(ich,'(a)',iostat=io)tmp
        if(io < 0)then
           write(och,'(''$end'')')
           exit
        endif
        write(och,'(a)')trim(tmp)
        if(index(tmp,'$end').ne.0)exit
      enddo

      outer: do
        read(ich,'(a)',iostat=io)tmp
        if(io < 0)then                   !done if no $set-block is present
           write(och,'(a)')trim(dg)      !just add the line dg at the end of the file
           exit
        endif
        if(trim(tmp) == '') cycle outer  !skip blank lines
        if(index(tmp,trim(key)) /= 0) then  !done if there is already a corresponding block
           write(och,'(a)')trim(key)
           if(index(trim(dg),trim(key)).eq.0)then
              write(och,'(a)')trim(dg)      !add the line, if it is not the key itself
           endif
           inner1: do  !copy the rest of the file
             read(ich,'(a)',iostat=io)tmp
             if(io < 0) exit outer
             write(och,'(a)')trim(tmp)
           enddo inner1
        endif
        if(index(tmp,'$set') /= 0) then  !done if the $set-block is present
           write(och,'(a)')trim(dg)      !add the line
           write(och,'(a)')'$set'        
           inner2: do  !copy the rest of the file
             read(ich,'(a)',iostat=io)tmp
             if(io < 0) exit outer
             write(och,'(a)')trim(tmp)
           enddo inner2
       endif
       write(och,'(a)')trim(tmp)
      enddo outer

      close(ich,status='delete')
      close(och)

      call rename('.setdgtmp',trim(fnam))

end subroutine setdg_block_key


!-----------------------------------------------------------------------------------------
! screen file and correct the $end-statement positions
!-----------------------------------------------------------------------------------------
subroutine clear_end(fnam)
      implicit none
      character(len=*),intent(in) :: fnam
      character(len=512) :: tmp
      integer :: ich,och,io
      integer :: i,j,k,nn

      open(newunit=ich,file=trim(fnam))
      open(newunit=och,file='.setdgtmp')
      do !check for $coord keyword
        read(ich,'(a)',iostat=io)tmp
        if(io < 0)then
           return
           exit
        endif
        write(och,'(a)')trim(tmp)
        if(index(tmp,'$coord').ne.0)exit
      enddo
      do !check for first $end keyword
        read(ich,'(a)',iostat=io)tmp
        if(io < 0)then
           exit
        endif
        write(och,'(a)')trim(tmp)
        if(index(tmp,'$end').ne.0)exit
      enddo
      do !check for the rest of the file, keep no $end 
        read(ich,'(a)',iostat=io)tmp
        if(io < 0)then                   !exit loop
           exit
        endif
        if(trim(tmp) == '') cycle        !skip blank lines
        if(index(tmp,'$end') /= 0) cycle !cycle through $end
        write(och,'(a)')trim(tmp)   !copy everything else
      enddo

      write(och,'(a)')'$end' !write a $end at the end

      close(ich,status='delete')
      close(och)

      call rename('.setdgtmp',trim(fnam))

end subroutine clear_end

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
end module setdgmod
