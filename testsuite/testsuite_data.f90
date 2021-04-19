module testsuitedata
   use mctc_env_accuracy 
   use crest_data 
   use strucrd
   implicit none
   private

   public :: build_env

   !> data routines
   public :: get_uracil

contains


!> build a systemdata object with default information
subroutine build_env(env)

   type(systemdata) :: env
   integer :: nra
   character(len=256) arg(1)
   nra = 1
   call parseflags(env,arg,nra) !defaults are st in the parser routine

   return
end subroutine build_env


!> some presaved data
subroutine get_uracil(mol)
    implicit none
    type(coord) :: mol
    mol%nat = 12
    allocate(mol%xyz(3,12), source=0.0_wp)
    allocate(mol%at(12))

    mol%at = (/8,6,7,6,6,6,8,7,1,1,1,1/)
    mol%xyz = reshape( (/&
       &  1.05396_wp,  0.02675_wp, -0.16406_wp, &  
       &  2.40385_wp, -0.01045_wp, -0.01820_wp, &
       &  2.98978_wp,  1.18949_wp,  0.02887_wp, &
       &  4.32632_wp,  1.15814_wp,  0.17720_wp, &
       &  5.05350_wp, -0.00667_wp,  0.27393_wp, &
       &  4.31437_wp, -1.16540_wp,  0.20726_wp, &
       &  4.91576_wp, -2.38801_wp,  0.28763_wp, &
       &  2.98691_wp, -1.20331_wp,  0.06179_wp, &
       &  0.93344_wp,  0.98938_wp, -0.18618_wp, &
       &  4.81293_wp,  2.12851_wp,  0.21720_wp, &
       &  6.12608_wp,  0.00635_wp,  0.39479_wp, &
       &  5.87206_wp, -2.26141_wp,  0.36107_wp &
    &/) , (/3, 12/) )

    return
end subroutine get_uracil



end module testsuitedata
