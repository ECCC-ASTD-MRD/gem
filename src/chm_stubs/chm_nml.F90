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
! Creation       : Alexander Kallaur, Sept. 2011
! Description    : Read chemistry namelists and initialize chemistry configuration
!             ------------> STUB VERSION  STUB VERSION  STUB VERSION <------------
!
! Extra info     :
!
! Arguments:
!           IN     :  F_namelist-> name of namelist file
!                     lun_out   -> Unit number of diagnostic message outputs
!
!
!           IN/OUT : Logival value of chm_nml function
!
!=============================================================================

integer function chm_nml (F_namelist,lun_out)
   implicit none
!
!  Input arguments
!
   character(len=*), intent(in) :: F_namelist
   integer,        intent(in) :: lun_out
 
   chm_nml = 1
 
   if (F_namelist == '' .and. lun_out >= 0) print *,'Called chm_nml stub with:', F_namelist, lun_out

end function chm_nml
