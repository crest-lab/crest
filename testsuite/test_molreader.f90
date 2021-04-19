module test_molreader
   use mctc_env_accuracy, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, &
      & test_failed
   use strucrd
   use iomod
   use testsuitedata
   implicit none
   private

   public :: collect_molreader

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: thr3 = 1.0e-8_wp
   real(wp),parameter :: bohr  = 0.52917726_wp
   integer :: io,i

contains


!> Collect all exported unit tests
subroutine collect_molreader(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("uracil_coord", test_coord_uracil), &
      & new_unittest("uracil_xyz", test_xyz_uracil) &  !, &
      & ]

end subroutine collect_molreader


subroutine test_rw_gen(error, filename, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> file name to read
   character(len=*) :: filename

   !> Reference stucture data
   type(coord), intent(in) :: ref

   !> Molecular structure data
   type(coord) :: mol

   call mol%open(filename)

   !> test nat
   if (mol%nat .ne. ref%nat) then
      call test_failed(error, "Number of atoms numbers do not match")
      print'(3es21.14)', mol%nat
   end if

   !> test atom types
   do i=1,mol%nat
     if(mol%at(i) .ne. ref%at(i))then
         call test_failed(error, "Atom types do not match")
     print'(3es21.14)', mol%at(i)
     endif
   enddo

   !> test coordinates
   if (any(abs(mol%xyz - ref%xyz) > thr3)) then
      call test_failed(error, "Atomic coordinates do not match")
      print'(3es21.14)', mol%xyz
   end if

   return
end subroutine test_rw_gen

subroutine test_coord_uracil(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(coord) :: ref
   call get_uracil(ref)
   call ref%write('uracil.coord')
   call test_rw_gen(error, 'uracil.coord', ref)
   !call remove('uracil.coord')

end subroutine test_coord_uracil

subroutine test_xyz_uracil(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(coord) :: ref
   call get_uracil(ref)
   !call ref%write('uracil.xyz')
   open(newunit=io,file='uracil.xyz')
   write(io,'(3x,i0)') ref%nat
   write(io,'(3x,f16.8)') -10.12345_wp
   do i=1,ref%nat
     write(io,'(a4,3f16.8)') i2e(ref%at(i),'nc'),ref%xyz(1:3,i)*bohr
   enddo
   call test_rw_gen(error, 'uracil.xyz', ref)
   !call remove('uracil.xyz')

end subroutine test_xyz_uracil









end module test_molreader
