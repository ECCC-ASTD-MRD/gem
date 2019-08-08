!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html

!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

subroutine CONVECT_CONDENS1(KLON, &
     &  KICE, PPRES, PTHL, PRW, PRCO, PRIO, PZ, &
     &  PT, PEW, PRC, PRI, PLV, PLS, PCPH)
   !!**** Compute temperature cloud and ice water content from enthalpy and r_w
   !
   !
   !!    PURPOSE
   !!    -------
   !!     The purpose of this routine is to determine cloud condensate
   !!     and to return values for L_v, L_s and C_ph
   !
   !
   !!**  METHOD
   !!    ------
   !!     Condensate is extracted iteratively
   !
   !
   !!    EXTERNAL
   !!    --------
   !!     None
   !
   !
   !!    IMPLICIT ARGUMENTS
   !!    ------------------
   !
   !!      Module YOMCST
   !!          RG                   ! gravity constant
   !!          RALPW, RBETW, RGAMW ! constants for water saturation pressure
   !!          RALPS, RBETS, RGAMS ! constants for ice saturation pressure
   !!          RATM                 ! reference pressure
   !!          RD, RV             ! gaz  constants for dry air and water vapor
   !!          RCPD, RCPV           ! specific heat for dry air and water vapor
   !!          RCW, RCS             ! specific heat for liquid water and ice
   !!          RTT                  ! triple point temperature
   !!          RLVTT, RLSTT         ! vaporization, sublimation heat constant
   !
   !!    IMPLICIT ARGUMENTS
   !!    ------------------
   !!      Module YOE_CONVPAR
   !!          XTFRZ1               ! begin of freezing interval
   !!          XTFRZ2               ! end of freezing interval
   !
   !!    REFERENCE
   !!    ---------
   !
   !!      Book1,2 of documentation ( routine CONVECT_CONDENS)
   !
   !!    AUTHOR
   !!    ------
   !!      P. BECHTOLD       * Laboratoire d'Aerologie *
   !
   !!    MODIFICATIONS
   !!    -------------
   !!      Original    07/11/95
   !!   Last modified  04/10/97
   !-------------------------------------------------------------------------------


   !*       0.    DECLARATIONS
   !              ------------

   use YOMCST
   use YOE_CONVPAR

   implicit none
!!!#include <arch_specific.hf>
#define _ZERO_   0.0
#define _ONE_    1.0

   !*       0.1   Declarations of dummy arguments :

   integer, intent(IN)                :: KLON    ! horizontal loop index
   integer, intent(IN)                :: KICE    ! flag for ice ( 1 = yes,
   !                0 = no ice )
   real, dimension(KLON),   intent(IN) :: PPRES  ! pressure
   real, dimension(KLON),   intent(IN) :: PTHL   ! enthalpy (J/kg)
   real, dimension(KLON),   intent(IN) :: PRW    ! total water mixing ratio
   real, dimension(KLON),   intent(IN) :: PRCO   ! cloud water estimate (kg/kg)
   real, dimension(KLON),   intent(IN) :: PRIO   ! cloud ice   estimate (kg/kg)
   real, dimension(KLON),   intent(IN) :: PZ     ! level height (m)

   real, dimension(KLON),   intent(OUT):: PT     ! temperature
   real, dimension(KLON),   intent(OUT):: PRC    ! cloud water mixing ratio(kg/kg)
   real, dimension(KLON),   intent(OUT):: PRI    ! cloud ice mixing ratio  (kg/kg)
   real, dimension(KLON),   intent(OUT):: PLV    ! latent heat L_v
   real, dimension(KLON),   intent(OUT):: PLS    ! latent heat L_s
   real, dimension(KLON),   intent(OUT):: PCPH   ! specific heat C_ph
   real, dimension(KLON),   intent(OUT):: PEW    ! water saturation mixing ratio

   !*       0.2   Declarations of local variables KLON

   integer :: JITER          ! iteration index
   real    :: ZEPS, ZEPSA    ! R_d / R_v, 1 / ZEPS
   real    :: ZCVOCD         ! RCPV / RCPD
   real    :: ZRDOCP         ! R_d / C_pd

   real, dimension(KLON)    :: ZEI           ! ice saturation mixing ratio
   real, dimension(KLON)    :: ZWORK1, ZWORK2, ZWORK3, ZT ! work arrays


   !-------------------------------------------------------------------------------

   !*       1.     Initialize temperature and Exner function
   !               -----------------------------------------

   ZRDOCP      = RD   / RCPD
   ZEPS        = RD   / RV
   ZEPSA       = _ONE_ / ZEPS
   ZCVOCD      = RCPV  / RCPD


   ! Make a first temperature estimate, based e.g. on values of
   !  r_c and r_i at lower level

   !! Note that the definition of ZCPH is not the same as used in
   !! routine CONVECT_SATMIXRATIO
   PCPH(:)   = RCPD + RCPV * PRW(:)
   ZWORK1(:) = ( _ONE_ + PRW(:) ) * RG * PZ(:)
   PT(:)     = ( PTHL(:) + PRCO(:) * RLVTT + PRIO(:) * RLSTT - ZWORK1(:) )   &
        & / PCPH(:)
   PT(:)     = max(180., min( 330., PT(:) ) ) ! set overflow bounds in
   ! case that PTHL=0


   !*       2.     Enter the iteration loop
   !               ------------------------

   do JITER = 1,6
      PEW(:) = exp( RALPW - RBETW / PT(:) - RGAMW * log( PT(:) ) )
      ZEI(:) = exp( RALPS - RBETS / PT(:) - RGAMS * log( PT(:) ) )
      PEW(:) = ZEPS * PEW(:) / ( PPRES(:) - PEW(:) )
      ZEI(:) = ZEPS * ZEI(:) / ( PPRES(:) - ZEI(:) )

      PLV(:)    = RLVTT + ( RCPV - RCW ) * ( PT(:) - RTT ) ! compute L_v
      PLS(:)    = RLSTT + ( RCPV - RCS ) * ( PT(:) - RTT ) ! compute L_i

      ZWORK2(:) = ( PT(:) - XTFRZ2 ) / ( XTFRZ1 - XTFRZ2 ) ! freezing interval
      ZWORK2(:) = max( _ZERO_, min(_ONE_, ZWORK2(:) ) )
      if ( KICE==0 ) ZWORK2(:) = _ONE_
      ZWORK2(:) = ZWORK2(:) * ZWORK2(:)
      ZWORK3(:) = ( _ONE_ - ZWORK2(:) ) * ZEI(:) + ZWORK2(:) * PEW(:)
      PRC(:)    = max( _ZERO_, ZWORK2(:) * ( PRW(:) - ZWORK3(:) ) )
      PRI(:)    = max( _ZERO_, ( _ONE_ - ZWORK2(:) ) * ( PRW(:) - ZWORK3(:) ) )
      ZT(:)     = ( PTHL(:) + PRC(:) * PLV(:) + PRI(:) * PLS(:) - ZWORK1(:) )   &
           & / PCPH(:)
      PT(:) = PT(:) + ( ZT(:) - PT(:) ) * 0.4  ! force convergence
      PT(:) = max( 175., min( 330., PT(:) ) )
   enddo


end subroutine CONVECT_CONDENS1

