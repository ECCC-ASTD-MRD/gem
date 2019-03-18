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
!             busent    -->    dimension of busent_size
!             busdyn    -->    busdyn_size
!             busper    -->    busper_size (also refered as fsiz elsewhere)
!             busvol    -->    busvol_size
!
!             busent_size        -->    dimension of busent   (physics)
!             busdyn_size        -->    dimension of busdyn   (physics)
!             fsiz        -->    dimension of busper   (physics)
!             busvol        -->    dimension of F_bvol_siz (physics)
!             dt              -->    timestep (sec.)
!             slab_index               -->    slice number
!             step            -->    timestep number
!             F_cpu             -->    cpu number executing slice "slab_index"
!             ni            -->    horizontal running length (physics)
!             nk              -->    vertical dimension
!
!            IN/OUT
!             F_chmbusent       -->   Entry input field(s)  (chemical)
!             F_chmbusdyn       -->   Dynamics input field  (chemical)
!             F_chmbusper       -->   Historic variables for the chemistry
!             F_chmbusvol       -->   Chemical tendencies and other output
!                                          fields from the physics
!
!  Modification1: -----> STUB VERSION STUB VERSION STUB VERSION STUB VERSION <-----
!
!=============================================================================

subroutine chm_exe(busent, busdyn, busper, busvol, &
     busent_size, busdyn_size, busper_size, busvol_size, &
     timestep, slab_index, step, ni, nk)

   implicit none

   ! Subroutine arguments
   integer, intent(in) :: busent_size, busdyn_size, busper_size, busvol_size 
   real, intent(inout) :: busent(busent_size), busdyn(busdyn_size), busper(busper_size), busvol(busvol_size)
   integer, intent(in) :: timestep, slab_index, step, ni, nk



   return
end subroutine chm_exe
