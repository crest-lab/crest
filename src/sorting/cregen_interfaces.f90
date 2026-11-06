
!=========================================================================================!
!=========================================================================================!
!> Interfaces for use CREGEN (and related)
!=========================================================================================!
!=========================================================================================!

module cregen_interface
!*******************************************************
!* module to load an interface to the newcregen routine
!* mandatory to handle the optional input arguments
!*******************************************************
  use unionize_module
  implicit none
  interface
    subroutine newcregen(env,quickset,infile,structurelist)
      use crest_parameters
      use crest_data
      use crest_restartlog
      use strucrd
      implicit none
      type(systemdata),intent(inout) :: env
      integer,intent(in),optional :: quickset
      character(len=*),intent(in),optional :: infile
      type(coord),allocatable,intent(inout),optional :: structurelist(:)
    end subroutine newcregen

    subroutine cregen_irmsd_all(nall,structures,printlvl,iinversion)
      use strucrd
      implicit none
      !> INPUT
      integer,intent(in) :: nall
      type(coord),intent(inout),target :: structures(nall)
      integer,intent(in),optional :: printlvl
      integer,intent(in),optional :: iinversion
    end subroutine cregen_irmsd_all

    subroutine cregen_irmsd_sort(env,nall,structures,groups,allcanon,printlvl)
      use crest_data
      use strucrd
      implicit none
      !> INPUT
      type(systemdata),intent(inout) :: env
      integer,intent(in) :: nall
      type(coord),intent(inout),target :: structures(nall)
      integer,intent(inout) :: groups(nall)
      logical,intent(in),optional :: allcanon
      integer,intent(in),optional :: printlvl
    end subroutine cregen_irmsd_sort

  end interface
!>--- Additional Related RE-EXPORTS
  public :: unionizeEnsembles
end module cregen_interface

!=========================================================================================!
!=========================================================================================!
!> Interfaces for routines used WITHIN CREGEN
!=========================================================================================!
!=========================================================================================!

module cregen_subroutines
!*************************************
!* interfaces for cregen subroutines
!*************************************
  implicit none
  interface

    subroutine cregen_files(env,fname,oname,cname,simpleset,userinput,ensembleinput,iounit)
      use crest_data
      implicit none
      type(systemdata),intent(inout) :: env
      character(len=:),allocatable,intent(inout) :: fname
      character(len=:),allocatable,intent(inout) :: oname
      character(len=:),allocatable,intent(inout) :: cname
      integer,intent(in) :: simpleset
      logical,intent(in) :: userinput
      logical,intent(in) :: ensembleinput
      integer,intent(out) :: iounit
    end subroutine cregen_files

    subroutine cregen_discardbroken(ch,env,topocheck,structures,newnall)
      use crest_data
      use strucrd
      use cregen_utils
      type(systemdata),intent(in) :: env
      integer,intent(in) :: ch
      logical,intent(in) :: topocheck
      type(coord),intent(inout),allocatable,target :: structures(:)
      integer,intent(out) :: newnall
    end subroutine cregen_discardbroken

    subroutine cregen_topocheck(ch,env,checkez,structures,newnall)
      use crest_data
      use strucrd
      use cregen_utils
      implicit none
      type(systemdata) :: env
      integer,intent(in) :: ch
      logical,intent(in) :: checkez
      type(coord),intent(inout),allocatable,target :: structures(:)
      integer,intent(out) :: newnall
    end subroutine cregen_topocheck

    subroutine cregen_esort(ch,structures,nallout,ewin)
      use crest_parameters
      use strucrd
      use quicksort_interface
      implicit none
      integer,intent(in) :: ch
      type(coord),intent(inout),allocatable :: structures(:)
      integer,intent(out) :: nallout
      real(wp),intent(in),optional :: ewin
    end subroutine cregen_esort

    subroutine cregen_CRE_new(env,nall,structures,groups,rthresh,ethr,bthr, &
    &                         printlvl,ch)
      use crest_parameters
      use crest_data
      use strucrd
      implicit none
      type(systemdata),intent(inout) :: env
      integer,intent(inout) :: nall
      type(coord),intent(inout),allocatable,target :: structures(:)
      integer,intent(out),allocatable :: groups(:)
      real(wp),intent(in) :: RTHRESH
      real(wp),intent(in) :: ETHR
      real(wp),intent(in) :: BTHR
      integer,intent(in),optional :: printlvl
      integer,intent(in),optional :: ch
    end subroutine cregen_CRE_new

    function cregen_majority_periodic(structures) result(periodic)
      use crest_parameters
      use strucrd
      implicit none
      type(coord),intent(in) :: structures(:)
      logical :: periodic
    end function cregen_majority_periodic

    subroutine cregen_discard_periodic(ch,structures)
      use crest_parameters
      use strucrd
      implicit none
      integer,intent(in) :: ch
      type(coord),intent(inout),allocatable :: structures(:)
    end subroutine cregen_discard_periodic

    subroutine cregen_CRE_periodic(env,nall,structures,groups,rthresh,ethr,printlvl,ch)
      use crest_parameters
      use crest_data
      use strucrd
      implicit none
      type(systemdata),intent(inout) :: env
      integer,intent(inout) :: nall
      type(coord),intent(inout),allocatable,target :: structures(:)
      integer,intent(out),allocatable :: groups(:)
      real(wp),intent(in) :: RTHRESH
      real(wp),intent(in) :: ETHR
      integer,intent(in),optional :: printlvl
      integer,intent(in),optional :: ch
    end subroutine cregen_CRE_periodic

    subroutine cregen_rmsdalign(nall,structures)
      use crest_parameters
      use irmsd_module
      use strucrd
      implicit none
      integer,intent(in) :: nall
      type(coord),intent(inout) :: structures(nall)
    end subroutine cregen_rmsdalign

    subroutine cregen_file_wr(env,fname,structures)
      use crest_parameters
      use crest_data
      use strucrd
      use utilities,only:boltz
      implicit none
      type(systemdata),intent(inout) :: env
      character(len=*),intent(in) :: fname
      type(coord),intent(inout) :: structures(:)
    end subroutine cregen_file_wr

    subroutine cregen_conffile(env,cname,structures,ng,degen)
      use crest_parameters
      use crest_data
      use strucrd
      implicit none
      type(systemdata),intent(inout) :: env
      character(len=*),intent(in) :: cname
      type(coord),intent(inout) :: structures(:)
      integer,intent(in) :: ng
      integer,intent(in) :: degen(3,ng)
    end subroutine cregen_conffile
  end interface
end module cregen_subroutines
