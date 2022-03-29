!     ##################
      MODULE YOE_CONDENS
!     ##################

!!****  *YOE_CONDENS * -  constants module
!!
!!    PURPOSE
!!    -------
!       The module contains constants used in routine CONDENS.
!!
!!**  METHOD
!!    ------
!!      The "condensation" constants are set to their numerical values
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (module YOE_CONDENS)
!!
!!
!!    AUTHOR
!!    ------
!!      P. Marquet   *Meteo-France CNRM/GMGEC/EAC*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/06/2002
!!   Last modified  07/06/2002
!-------------------------------------------------------------------------------

!*       0.    DECLARATIONS
!              ------------


IMPLICIT NONE

!-------------------------------------------------------------------------------

!*       1.    Set the thermodynamical and numerical constants for
!              the deep convection parameterization
!              ---------------------------------------------------


real :: XCTFRZ1    ! begin of freezing interval
real :: XCTFRZ2    ! end of freezing interval

real :: XCLSGMAX   ! Max-tropospheric length for Sigma_s 
real :: XCLSGMIN   ! Lower-height Min-length for Sigma_s

real :: XCSIGMA    ! constant in sigma_s parameterization

real :: XCSIG_CONV ! scaling factor for ZSIG_CONV as
                   ! function of mass flux

END MODULE YOE_CONDENS
