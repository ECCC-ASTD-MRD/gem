!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
! version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

subroutine CONVECT_CHEM_TRANSPORT1(KLON, KLEV, KCH, PCH1, PCH1C,     &
     & KDPL, KPBL, KLCL, KCTL, KLFS, KDBL, &
     & PUMF, PUER, PUDR, PDMF, PDER, PDDR, &
     & PTIMEC, PDXDY, PMIXF, PLMASS, PWSUB)
   !!**** Compute  modified chemical tracer values due to convective event
   !!
   !!
   !!    PURPOSE
   !!    -------
   !!      The purpose of this routine is to determine the final adjusted
   !!      environmental values of the chemical tracers
   !!      The final convective tendencies can then be evaluated in the main
   !!      routine DEEP_CONVECT by (PCH1C-PCH1)/PTIMEC
   !!
   !!
   !!**  METHOD
   !!    ------
   !!      Identical to the computation of the conservative variables in the
   !!      main deep convection code
   !!
   !!    EXTERNAL
   !!    --------
   !!
   !!    IMPLICIT ARGUMENTS
   !!    ------------------
   !!      Module YOMCST
   !!          RG                 ! gravity constant
   !!
   !!     Module YOE_CONVPAREXT
   !!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
   !!
   !!    AUTHOR
   !!    ------
   !!      P. BECHTOLD       * Laboratoire d'Aerologie *
   !!
   !!    MODIFICATIONS
   !!    -------------
   !!
   !!      Original    11/12/97
   !!
   !----------------------------------------------------------------------------

   !*       0.    DECLARATIONS
   !              ------------

   use YOMCST
   use YOE_CONVPAREXT

   implicit none

!!!#include <arch_specific.hf>
#define _ZERO_   0.0
#define _HALF_   0.5
#define _ONE_    1.0

   !*       0.1   Declarations of dummy arguments :

   integer,                intent(IN) :: KLON     ! horizontal dimension
   integer,                intent(IN) :: KLEV     ! vertical dimension
   integer,                intent(IN) :: KCH      ! number of passive tracers

   real,dimension(KLON,KLEV,KCH),intent(IN) :: PCH1 ! grid scale tracer concentr.
   real,dimension(KLON,KLEV,KCH),intent(OUT):: PCH1C! conv adjusted tracer concntr.

   integer, dimension(KLON), intent(IN) :: KDPL   ! index for departure level
   integer, dimension(KLON), intent(IN) :: KPBL   ! index for top of source layer
   integer, dimension(KLON), intent(IN) :: KLCL   ! index lifting condens. level
   integer, dimension(KLON), intent(IN) :: KCTL   ! index for cloud top level
   integer, dimension(KLON), intent(IN) :: KLFS   ! index for level of free sink
   integer, dimension(KLON), intent(IN) :: KDBL   ! index for downdraft base level

   real, dimension(KLON,KLEV), intent(IN) :: PUMF ! updraft mass flux (kg/s)
   real, dimension(KLON,KLEV), intent(IN) :: PUER ! updraft entrainment (kg/s)
   real, dimension(KLON,KLEV), intent(IN) :: PUDR ! updraft detrainment (kg/s)
   real, dimension(KLON,KLEV), intent(IN) :: PDMF ! downdraft mass flux (kg/s)
   real, dimension(KLON,KLEV), intent(IN) :: PDER ! downdraft entrainment (kg/s)
   real, dimension(KLON,KLEV), intent(IN) :: PDDR ! downdraft detrainment (kg/s)

   real, dimension(KLON),     intent(IN) :: PTIMEC! convection time step
   real, dimension(KLON),     intent(IN) :: PDXDY ! grid area (m^2)
   real, dimension(KLON),     intent(IN) :: PMIXF ! mixed fraction at LFS
   real, dimension(KLON,KLEV),intent(IN) :: PLMASS! mass of model layer (kg)
   real, dimension(KLON,KLEV),intent(IN) :: PWSUB ! envir. compensating subsidence(Pa/s)

   !*       0.2   Declarations of local variables :

   integer :: INCH1          ! number of chemical tracers
   integer :: IIE, IKB, IKE  ! horizontal + vertical loop bounds
   integer :: IKS            ! vertical dimension
   integer :: JI             ! horizontal loop index
   integer :: JK, JKP        ! vertical loop index
   integer :: JN             ! chemical tracer loop index
   integer :: JSTEP          ! fractional time loop index
   integer :: JKLD, JKLP     ! loop index for levels
   integer, dimension(KLON) :: ITSTEP !fractional convective time step

   real, dimension(KLON,KLEV)     :: ZOMG ! compensat. subsidence (Pa/s)
   real, dimension(KLON,KLEV,KCH) :: ZUCH1, ZDCH1 ! updraft/downdraft values
   real, dimension(KLON)          :: ZTIMEC! adjust fractional convective time step
   real, dimension(KLON,KLEV)     :: ZTIMC! 2D work array for ZTIMEC
   real, dimension(KLON,KLEV,KCH) :: ZCH1MFIN, ZCH1MFOUT
   ! work arrays for environm. compensat. mass
   real, dimension(KLON,KCH)      :: ZWORK1, ZWORK2, ZWORK3

   !----------------------------------------------------------------------------

   !*       0.3   Compute loop bounds
   !              -------------------

   INCH1  = KCH
   IIE    = KLON
   IKB    = 1 + JCVEXB
   IKS    = KLEV
   IKE    = KLEV - JCVEXT
   !JKMAX  = MAXVAL( KCTL(:) ) -> KCTL(JI)


   !*      2.      Updraft computations
   !               --------------------

   ZUCH1(:,:,:) = _ZERO_

   !*      2.1     Initialization  at LCL
   !               ----------------------------------

   do JI = 1, IIE
      JKLD = KDPL(JI)
      JKLP = KPBL(JI)
      ZWORK1(JI,:) = _HALF_ * ( PCH1(JI,JKLD,:) + PCH1(JI,JKLP,:) )
   enddo

   !*      2.2     Final updraft loop
   !               ------------------

   do JI = 1, IIE
      !DO JK = MINVAL( KDPL(:) ), KCTL(JI)
      do JK = KDPL(JI), KCTL(JI)
         JKP = JK + 1

         do JN = 1, INCH1
            if ( KDPL(JI) <= JK .and. KLCL(JI) > JK ) ZUCH1(JI,JK,JN) = ZWORK1(JI,JN)

            if ( KLCL(JI) - 1 <= JK .and. KCTL(JI) > JK ) then
               !if you have reactive i.e. non-passive tracers
               ! add the corresponding sink term in the following equation
               ZUCH1(JI,JKP,JN) = ( PUMF(JI,JK) * ZUCH1(JI,JK,JN) +  &
                    &   PUER(JI,JKP) * PCH1(JI,JK,JN) )  /           &
                    & ( PUMF(JI,JKP) + PUDR(JI,JKP) + 1.E-7 )
            endif
         enddo
      enddo
   enddo

   !*      3.      Downdraft computations
   !               ----------------------

   ZDCH1(:,:,:) = _ZERO_

   !*      3.1     Initialization at the LFS
   !               -------------------------

   ZWORK1(:,:) = spread( PMIXF(:), DIM=2, NCOPIES=INCH1 )
   do JI = 1, IIE
      JK = KLFS(JI)
      ZDCH1(JI,JK,:) = ZWORK1(JI,:) * PCH1(JI,JK,:) +                          &
           &                  ( _ONE_ - ZWORK1(JI,:) ) * ZUCH1(JI,JK,:)
   enddo

   !*      3.2     Final downdraft loop
   !               --------------------

   !DO JK = MAXVAL( KLFS(:) ), IKB + 1, -1
   do JI = 1, IIE
      do JK = KLFS(JI), IKB + 1, -1
         JKP = JK - 1
         do JN = 1, INCH1
            if ( JK <= KLFS(JI) .and. JKP >= KDBL(JI) ) then
               ZDCH1(JI,JKP,JN) = ( ZDCH1(JI,JK,JN) * PDMF(JI,JK) -  &
                    &   PCH1(JI,JK,JN) *  PDER(JI,JKP) ) /           &
                    & ( PDMF(JI,JKP) - PDDR(JI,JKP) - 1.E-7 )
            endif
         enddo
      enddo
   enddo


   !*      4.      Final closure (environmental) computations
   !               ------------------------------------------

   PCH1C(:,IKB:IKE,:) = PCH1(:,IKB:IKE,:) ! initialize adjusted envir. values

   do JK = IKB, IKE
      ZOMG(:,JK) = PWSUB(:,JK) * PDXDY(:) / RG ! environmental subsidence
   enddo

   !PLEASE CHECK HERE!!!
   !Recalculate the fractional time step like in convect_closure_shal -vlee
   !ZTIMEC(:) = PTIMEC(:) / REAL( KFTSTEPS ) ! adjust  fractional time step
   ! to be an integer multiple of PTIMEC
   ZTIMEC(:) = PTIMEC(:)
   ITSTEP(:) = int( PTIMEC(:) / ZTIMEC(:) ) + 1
   ZTIMEC(:) = PTIMEC(:) / real( ITSTEP(:) )
   where ( PTIMEC(:) < _ONE_ ) ZTIMEC(:) = _ZERO_
   ZTIMC(:,:)= spread( ZTIMEC(:), DIM=2, NCOPIES=IKS )

   !These variables are inside the loop in the convect_uv_transport_shal
   !perhaps due to the fractional time step loop that I have implemented
   !in here with itstep(:) instead of kftsteps -vlee
   !PLEASE CHECK HERE
   !ZCH1MFIN(:,:,:)   = _ZERO_
   !ZCH1MFOUT(:,:,:)  = _ZERO_


   do JI = 1, IIE
      FRACTIONAL_STEPS:do JSTEP = 1, ITSTEP(JI)! Enter the fractional time step loop
         !DO JSTEP = 1, KFTSTEPS ! Enter the fractional time step loop
         ZCH1MFIN(:,:,:)   = _ZERO_
         ZCH1MFOUT(:,:,:)  = _ZERO_

         do JK = IKB + 1, KCTL(JI)
            JKP = max( IKB + 1, JK - 1 )
            ZWORK3(JI,1) = ZOMG(JI,JK)
            ZWORK1(JI,1) = sign( _ONE_, ZWORK3(JI,1) )
            ZWORK2(JI,1) = _HALF_ * ( 1. + ZWORK1(JI,1) )
            ZWORK1(JI,1) = _HALF_ * ( 1. - ZWORK1(JI,1) )
            ZCH1MFIN(JI,JK,:)  = - ZWORK3(JI,1) * PCH1C(JI,JKP,:) * ZWORK1(JI,1)
            ZCH1MFOUT(JI,JK,:) =   ZWORK3(JI,1) * PCH1C(JI,JK,:)  * ZWORK2(JI,1)
            ZCH1MFIN(JI,JKP,:) = ZCH1MFIN(JI,JKP,:) + ZCH1MFOUT(JI,JK,:) * ZWORK2(JI,1)
            ZCH1MFOUT(JI,JKP,:)= ZCH1MFOUT(JI,JKP,:) + ZCH1MFIN(JI,JK,:) * ZWORK1(JI,1)
         end do
         !
         do JK = IKB + 1, KCTL(JI)
            do JN = 1, INCH1
               PCH1C(JI,JK,JN) = PCH1C(JI,JK,JN) + ZTIMC(JI,JK) / PLMASS(JI,JK) *  (    &
                    ZCH1MFIN(JI,JK,JN) + PUDR(JI,JK) * ZUCH1(JI,JK,JN) +      &
                    PDDR(JI,JK) * ZDCH1(JI,JK,JN) - ZCH1MFOUT(JI,JK,JN) -     &
                    ( PUER(JI,JK) + PDER(JI,JK) ) * PCH1(JI,JK,JN)    )
               !  PCH1C(JI,JK,JN) = MAX( 0., PCH1C(JI,JK,JN) )
            end do
         end do

      enddo FRACTIONAL_STEPS! final values
   enddo

end subroutine CONVECT_CHEM_TRANSPORT1

