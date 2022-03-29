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

subroutine CONVECT_UV_TRANSPORT_DEEP1(KLON, KLEV, PU, PV, PUC, PVC, &
     & KDPL, KPBL, KLCL, KCTL, KLFS, KDBL,KDDT, &
     & PUMF, PUER, PUDR, PDMF, PDER, PDDR,      &
     & PDXDY, PMIXF, PLMASS, PWSUB,     &
     & KFTSTEPS, PTIMC, ITSTEP, GWORK1)

   !!**** Compute  modified horizontal wind components due to convective event


   !!    PURPOSE
   !!    -------
   !!      The purpose of this routine is to determine convective adjusted
   !!      horizontal wind components u and v
   !!      The final convective tendencies can then be evaluated in the main
   !!      routine DEEP_CONVECT by (PUC-PU)/PTIMEC


   !!**  METHOD
   !!    ------
   !!      Identical to the computation of the conservative variables in the
   !!      tracer routine but includes pressure term

   !!    EXTERNAL
   !!    --------

   !!    IMPLICIT ARGUMENTS
   !!    ------------------
   !!      Module YOMCST
   !!          RG                 ! gravity constant

   !!     Module YOE_CONVPAREXT
   !!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries

   !!    AUTHOR
   !!    ------
   !!      P. BECHTOLD       * Laboratoire d'Aerologie *

   !!    MODIFICATIONS
   !!    -------------

   !!      Original    11/02/02

   !-------------------------------------------------------------------------------


   !*       0.    DECLARATIONS
   !              ------------

   use YOMCST
   use YOE_CONVPAR
   use YOE_CONVPAREXT

   implicit none
!!!#include <arch_specific.hf>
#define _ZERO_   0.0
#define _HALF_   0.5
#define _ONE_    1.0
#define _TWO_    2.0

   !*       0.1   Declarations of dummy arguments :

   integer,                intent(IN) :: KLON     ! horizontal dimension
   integer,                intent(IN) :: KLEV     ! vertical dimension

   real,dimension(KLON,KLEV),intent(IN) :: PU     ! horizontal wind in x (m/s)
   real,dimension(KLON,KLEV),intent(IN) :: PV     ! horizontal wind in x (m/s)
   real,dimension(KLON,KLEV),intent(OUT):: PUC    ! convective adjusted value of u (m/s)
   real,dimension(KLON,KLEV),intent(OUT):: PVC    ! convective adjusted value of v (m/s)

   integer, dimension(KLON), intent(IN) :: KDPL   ! index for departure level
   integer, dimension(KLON), intent(IN) :: KPBL   ! index for top of source layer
   integer, dimension(KLON), intent(IN) :: KLCL   ! index lifting condens. level
   integer, dimension(KLON), intent(IN) :: KCTL   ! index for cloud top level
   integer, dimension(KLON), intent(IN) :: KLFS   ! index for level of free sink
   integer, dimension(KLON), intent(IN) :: KDBL   ! index for downdraft base level
   integer, dimension(KLON), intent(IN) :: KDDT   ! index for downdraft detrainment top level (60mb)

   real, dimension(KLON,KLEV), intent(IN) :: PUMF ! updraft mass flux (kg/s)
   real, dimension(KLON,KLEV), intent(IN) :: PUER ! updraft entrainment (kg/s)
   real, dimension(KLON,KLEV), intent(IN) :: PUDR ! updraft detrainment (kg/s)
   real, dimension(KLON,KLEV), intent(IN) :: PDMF ! downdraft mass flux (kg/s)
   real, dimension(KLON,KLEV), intent(IN) :: PDER ! downdraft entrainment (kg/s)
   real, dimension(KLON,KLEV), intent(IN) :: PDDR ! downdraft detrainment (kg/s)

   real, dimension(KLON),     intent(IN) :: PDXDY ! grid area (m^2)
   real, dimension(KLON),     intent(IN) :: PMIXF ! mixed fraction at LFS
   real, dimension(KLON,KLEV),intent(IN) :: PLMASS! mass of model layer (kg)
   real, dimension(KLON,KLEV),intent(IN) :: PWSUB ! envir. compensating subsidence(Pa/s)
   integer,                   intent(IN) :: KFTSTEPS  ! maximum fractional time steps
   real, dimension(KLON),     intent(IN) :: PTIMC   ! fractional convective time step
   integer,dimension(KLON),   intent(IN) :: ITSTEP  ! # of fractional convective timesteps
   logical, dimension(KLON),  intent(IN) :: GWORK1  ! flag for newly activated columns including MASS CONSERVATION criteria from closure

   !*       0.2   Declarations of local variables :

   integer :: IIE, IKB, IKE  ! horizontal + vertical loop bounds
   integer :: IKS            ! vertical dimension
   integer :: JI             ! horizontal loop index
   integer :: JK, JKP        ! vertical loop index
   integer :: JSTEP          ! fractional time loop index
   integer :: JKLD, JKLP, JKMAX ! loop index for levels

   integer, parameter             :: IUV = 2    ! for u and v
   real, dimension(KLON,KLEV)     :: ZOMG       ! compensat. subsidence (Pa/s)
   real, dimension(KLON,KLEV,IUV) :: ZUUV, ZDUV ! updraft/downdraft values
   real, dimension(KLON)          :: ZTIMEC     ! fractional convective time step
   real, dimension(KLON,KLEV)     :: ZTIMC      ! 2D work array for ZTIMEC
   real, dimension(KLON,KLEV,IUV) :: ZUVMFIN, ZUVMFOUT
   ! work arrays for environm. compensat. mass
   real, dimension(KLON,IUV)      :: ZWORK1, ZWORK2, ZWORK3

   integer, dimension(KLON)  :: ICOUNT       ! timestep counter
   logical, dimension(KLON)  :: GWORK3

   !-------------------------------------------------------------------------------

   !*       0.3   Compute loop bounds
   !              -------------------

   IIE    = KLON
   IKB    = 1 + JCVEXB
   IKS    = KLEV
   IKE    = KLEV - JCVEXT
   JKMAX  = maxval( KCTL(:) )


   !*      2.      Updraft computations
   !               --------------------

   ZUUV(:,:,:) = _ZERO_

   !*      2.1     Initialization  at LCL
   !               ----------------------------------

   do JI = 1, IIE
      JKLD = KDPL(JI)
      JKLP = KPBL(JI)
      ZWORK1(JI,1) = _HALF_ * ( PU(JI,JKLD) + PU(JI,JKLP) )
      ZWORK1(JI,2) = _HALF_ * ( PV(JI,JKLD) + PV(JI,JKLP) )
   enddo

   !*      2.2     Final updraft loop
   !               ------------------

   do JK = minval( KDPL(:) ), JKMAX
      JKP = JK + 1

      do JI = 1, IIE
         if ( KDPL(JI) <= JK .and. KLCL(JI) > JK ) then
            ZUUV(JI,JK,1) = ZWORK1(JI,1)
            ZUUV(JI,JK,2) = ZWORK1(JI,2)
         end if

         if ( KLCL(JI) - 1 <= JK .and. KCTL(JI) > JK ) then
            ! instead of passive tracers equations
            ! wind equations also include pressure term
            ZUUV(JI,JKP,1) = ( PUMF(JI,JK) * ZUUV(JI,JK,1) +                   &
                 &   PUER(JI,JKP) * PU(JI,JK) )  /                 &
                 & ( PUMF(JI,JKP) + PUDR(JI,JKP) + 1.E-7 ) +  &
                 &   XUVDP * ( PU(JI,JKP) - PU(JI,JK) )
            ZUUV(JI,JKP,2) = ( PUMF(JI,JK) * ZUUV(JI,JK,2) +                   &
                 &   PUER(JI,JKP) * PV(JI,JK) )  /                 &
                 & ( PUMF(JI,JKP) + PUDR(JI,JKP) + 1.E-7 ) +  &
                 &   XUVDP * ( PV(JI,JKP) - PV(JI,JK) )
         endif
      enddo

   enddo

   !*      3.      Downdraft computations
   !               ----------------------

   ZDUV(:,:,:) = _ZERO_

   !*      3.1     Initialization at the LFS
   !               -------------------------

   ZWORK1(:,:) = spread( PMIXF(:), DIM=2, NCOPIES=IUV )
   do JI = 1, IIE
      JK = KLFS(JI)
      JKP=JK+1
      ZDUV(JI,JK,1) = ZWORK1(JI,1) * PU(JI,JK) +                          &
           &            ( _ONE_ - ZWORK1(JI,1) ) * ZUUV(JI,JK,1)
      ZDUV(JI,JK,2) = ZWORK1(JI,2) * PV(JI,JK) +                          &
           &            ( _ONE_ - ZWORK1(JI,2) ) * ZUUV(JI,JK,2)
   enddo

   !*      3.2     Final downdraft loop
   !               --------------------

   do JK = maxval( KLFS(:) ), IKB + 1, -1
      JKP = JK - 1
      do JI = 1, IIE
         if ( JK <= KLFS(JI) .and. JKP >= KDDT(JI) ) then
            ZDUV(JI,JKP,1) = ( ZDUV(JI,JK,1) * PDMF(JI,JK) -                  &
                 &     PU(JI,JK) * PDER(JI,JKP) ) /               &
                 & ( PDMF(JI,JKP) - PDDR(JI,JKP) - 1.E-7 ) + &
                 &   XUVDP * ( PU(JI,JKP) - PU(JI,JK) )
            ZDUV(JI,JKP,2) = ( ZDUV(JI,JK,2) * PDMF(JI,JK) -                  &
                 &     PV(JI,JK) *  PDER(JI,JKP) ) /              &
                 & ( PDMF(JI,JKP) - PDDR(JI,JKP) - 1.E-7 ) + &
                 &   XUVDP * ( PV(JI,JKP) - PV(JI,JK) )

         endif
         !Jing
         if ( JKP < KDDT(JI) .and. JKP >= KDBL(JI) ) then
            ZDUV(JI,JKP,1) =ZDUV(JI,JK,1)
            ZDUV(JI,JKP,2) =ZDUV(JI,JK,2)
         endif

         !Jing



      enddo
   enddo


   !*      4.      Final closure (environmental) computations
   !               ------------------------------------------

   PUC(:,IKB:IKE) = PU(:,IKB:IKE) ! initialize adjusted envir. values
   PVC(:,IKB:IKE) = PV(:,IKB:IKE) ! initialize adjusted envir. values

   do JK = IKB, IKE
      ZOMG(:,JK) = PWSUB(:,JK) * PDXDY(:) / RG ! environmental subsidence
   enddo

   ZTIMEC(:) = PTIMC(:)
   ! to be an integer multiple of PTIMEC
   where ( PTIMC(:) < _ONE_ ) ZTIMEC(:) = _ZERO_
   ZTIMC(:,:)= spread( ZTIMEC(:), DIM=2, NCOPIES=IKS )

   ICOUNT(:) = 0

   do JSTEP = 1, KFTSTEPS ! Enter the fractional time step loop

      ICOUNT(:) = ICOUNT(:) + 1
      GWORK3(:) =  ITSTEP(:) >= ICOUNT(:) .and. GWORK1(:)

      ZUVMFIN(:,:,:)   = _ZERO_
      ZUVMFOUT(:,:,:)  = _ZERO_

      do JK = IKB + 1, JKMAX
         JKP = max( IKB + 1, JK - 1 )
         do JI = 1, IIE
            if ( GWORK3(JI) .and. JK <= KCTL(JI) )  then
               ZWORK3(JI,1) = ZOMG(JI,JK)
               ZWORK1(JI,1) = sign( _ONE_, ZWORK3(JI,1) )
               ZWORK2(JI,1) = _HALF_ * ( _ONE_ + ZWORK1(JI,1) )
               ZWORK1(JI,1) = _HALF_ * ( _ONE_ - ZWORK1(JI,1) )
               ZUVMFIN(JI,JK,1)  = - ZWORK3(JI,1) * PUC(JI,JKP) * ZWORK1(JI,1)
               ZUVMFOUT(JI,JK,1) =   ZWORK3(JI,1) * PUC(JI,JK)  * ZWORK2(JI,1)
               ZUVMFIN(JI,JK,2)  = - ZWORK3(JI,1) * PVC(JI,JKP) * ZWORK1(JI,1)
               ZUVMFOUT(JI,JK,2) =   ZWORK3(JI,1) * PVC(JI,JK)  * ZWORK2(JI,1)
               ZUVMFIN(JI,JKP,1) = ZUVMFIN(JI,JKP,1) + ZUVMFOUT(JI,JK,1) * ZWORK2(JI,1)
               ZUVMFIN(JI,JKP,2) = ZUVMFIN(JI,JKP,2) + ZUVMFOUT(JI,JK,2) * ZWORK2(JI,1)
               ZUVMFOUT(JI,JKP,1)= ZUVMFOUT(JI,JKP,1)+ ZUVMFIN(JI,JK,1)  * ZWORK1(JI,1)
               ZUVMFOUT(JI,JKP,2)= ZUVMFOUT(JI,JKP,2)+ ZUVMFIN(JI,JK,2)  * ZWORK1(JI,1)
            end if
         end do
      end do

      do JK = IKB + 1, JKMAX
         do JI = 1, IIE
            if (  GWORK3(JI) .and. JK <= KCTL(JI) ) then
               PUC(JI,JK) = PUC(JI,JK) + ZTIMC(JI,JK) / PLMASS(JI,JK) *  (       &
                    &   ZUVMFIN(JI,JK,1) + PUDR(JI,JK) * ZUUV(JI,JK,1) +    &
                    &   PDDR(JI,JK) * ZDUV(JI,JK,1) - ZUVMFOUT(JI,JK,1) -   &
                    &   ( PUER(JI,JK) + PDER(JI,JK) ) * PU(JI,JK)    )
               PVC(JI,JK) = PVC(JI,JK) + ZTIMC(JI,JK) / PLMASS(JI,JK) *  (       &
                    &   ZUVMFIN(JI,JK,2) + PUDR(JI,JK) * ZUUV(JI,JK,2) +    &
                    &   PDDR(JI,JK) * ZDUV(JI,JK,2) - ZUVMFOUT(JI,JK,2) -   &
                    &   ( PUER(JI,JK) + PDER(JI,JK) ) * PV(JI,JK)    )
            end if
         end do
      end do

   enddo ! Exit the fractional time step loop


end subroutine CONVECT_UV_TRANSPORT_DEEP1

