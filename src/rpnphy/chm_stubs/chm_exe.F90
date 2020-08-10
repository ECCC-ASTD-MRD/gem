!============================================================================!
!         Environnement Canada         |        Environment Canada           !
!                                      |                                     !
! - Service meteorologique du Canada   | - Meteorological Service of Canada  !
! - Direction generale des sciences    | - Science and Technology Branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http://www.ec.gc.ca                             !
!============================================================================!
!
! Projet/Project : GEM-MACH
! Fichier/File   : chm_exe.ftn90
! Creation       : A. Kallaur, H. Landry, S. Menard - July 2005 and more
! Description    : Computes the chemical tranformation for a
!                  list of species and reactions in the atmosphere.
!
! Extra info     : The object is the Development of chosen chemical scheme that
!                  begins here in this subroutine: chemical transformations,
!                  with this chosen scheme.
!
! Arguments:
!           IN
!             busdyn    -->    busdyn_size
!             busper    -->    busper_size (also refered as fsiz elsewhere)
!             busvol    -->    busvol_size
!
!             busdyn_size        -->    dimension of busdyn   (physics)
!             fsiz      -->    dimension of busper   (physics)
!             busvol    -->    dimension of F_bvol_siz (physics)
!             slab_index-->    slice number
!             step      -->    timestep number
!
!  Modification1: -----> STUB VERSION STUB VERSION STUB VERSION STUB VERSION <-----
!
!=============================================================================

subroutine chm_exe2(busdyn, busper, busvol,  &
     busdyn_size, busper_size, busvol_size, &
     slab_index, step)

   implicit none

   ! Subroutine arguments
   integer, intent(in) :: busdyn_size, busper_size, busvol_size 
   real, intent(inout) :: busdyn(busdyn_size), busper(busper_size), busvol(busvol_size)
   integer, intent(in) :: slab_index, step

   if (slab_index < 0) print *,'Called chm_exe stub with:', busdyn(1), busper(1), busvol(1), &
        busdyn_size, busper_size, busvol_size, &
        slab_index, step


   return
end subroutine chm_exe2
