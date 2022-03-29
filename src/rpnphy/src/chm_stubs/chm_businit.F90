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
! Fichier/File   : chm_businit.ftn90
! Creation       : H. Landry (Janvier 2008)
! Description    : Initialize chemistry buses by calling gesdict for
!                  any declared field in the master array of species
!
!
! Extar info     : This s/r should remain independant of any chemical scheme,
!                  package, emissions or else
!
! Arguments: IN
!              F_ni --> horizontal dimension
!              F_nk       --> vertical dimension
!
! Modification1 : -----> STUB VERSION STUB VERSION STUB VERSION STUB VERSION <-----
!
!==============================================================================

subroutine chm_businit(F_ni, F_nk)

   implicit none

! Subroutine arguments
   integer, intent(in) :: F_ni, F_nk

   if (F_ni < 0) print *,'Called chm_businit stub with:', F_ni, F_nk

   return

end subroutine chm_businit
