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
! Modification2   : -----> STUB VERSION  STUB VERSION  STUB VERSION  STUB VERSION  <-----
! Extra info     :
!
! Arguments: IN
!              bper_siz       --> bus size
!              bper_len       --> bus length (bus size*p_nj)
!              gem_datimeo    --> GEM date of origin
!              gem_dt_8       --> Model time step length
!              gem_tstep_num  --> current Model time step
!                                 (Higher up; called F_stepno and Lctl_step)
!
!            IN-OUT
!              busper3D       --> Physics permanent bus 
!                                 (all chemical vars are now there)
!
!==============================================================================
subroutine chm_load_emissions (F_basedir_S, gem_datimeo, i_gem_dt, F_phy_lcl_ni, F_phydim_ni, gem_tstep_num)

   implicit none

! Subroutine arguments
   character(len=*)    :: F_basedir_S
   integer, intent(in) :: gem_datimeo
   integer, intent(in) :: gem_tstep_num
   integer, intent(in) :: i_gem_dt, F_phy_lcl_ni, F_phydim_ni

   return
end subroutine chm_load_emissions

