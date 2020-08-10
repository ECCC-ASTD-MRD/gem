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
! Fichier/File   : chm_load_emissions.ftn90
! Creation       : A. Kallaur (MSC/ARQI), Octobre 2006
! Description    : Read from file (if need be) with the surface
!                  biogenic, anthropogenic and major point sources.
!
! Modification1   : December 2010 (A. Kallaur)
!                   Adaptation to GEM4.4.0+ new interface.
!                   Chemical buses no longer exist, so all emissions data
!                   is written to the permanent Physics bus.
!
! Modification2   : June 2020 (D. Akingunola)
!                   Added GEM input API arguments for reading area emissions
!                   (GEM >= 5.1.x)
!         -----> STUB VERSION  STUB VERSION  STUB VERSION  STUB VERSION  <-----
! Extra info     :
!
! Arguments: IN
!              F_basedir_S    --> base path for input data file
!              gem_tstep_num  --> current Model time step
!              inputobj       --> INPUTIO object
!              nbvar_input    --> Number of fields in inputobj
!
!==============================================================================
subroutine chm_load_emissions2(F_basedir_S, gem_tstep_num, inputobj, nbvar_input)
   use inputio_mod, only: INPUTIO_T
   implicit none

! Subroutine arguments
   character(len=*)               :: F_basedir_S
   integer, intent(in)            :: gem_tstep_num
   type(INPUTIO_T), intent(inout) :: inputobj
   integer(kind=4), intent(inout) :: nbvar_input

   !#Note: the next line is only to prevent compiler complain for unused var
   if (gem_tstep_num < 0) print *,'Called chm_load_emissions stub with:', F_basedir_S, gem_tstep_num, nbvar_input !, inputobj

   return
end subroutine chm_load_emissions2

