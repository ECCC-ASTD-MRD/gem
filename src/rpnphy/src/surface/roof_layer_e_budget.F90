!-------------------------------------- LICENCE BEGIN ------------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer, 
!                     version 3; Last Modified: May 7, 2008.
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
!-------------------------------------- LICENCE END --------------------------------------
!   ##########################################################################
    SUBROUTINE ROOF_LAYER_E_BUDGET(PT_ROOF, PQSAT_ROOF,                        &
                             PTA, PQA, PPS,                                    &
                             PLW, PTSTEP,                                      &
                             PEMIS_ROOF, PHC_ROOF, PTC_ROOF, PD_ROOF,          &
                             PTI_BLD, PAC_BLD, PDELT_ROOF,                     &
                             PDELTSNOW_ROOF, PGSNOW_ROOF,                      &
                             PRHOA, PAC_ROOF, PAC_ROOF_WAT,                    &
                             PABS_SW_ROOF, PABS_LW_ROOF                        )
!   ##########################################################################
!
!!****  *ROOF_LAYER_E_BUDGET*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the evoultion of surface temperature of roofs
!         
!     
!!**  METHOD
!     ------
!
!
!
!
!    5 : equation for evolution of Ts_roof
!        *********************************
!
!     dTt_1(t) / dt = 1/(dt_1*Ct_1) * (  Rn - H - LE
!                                      - 2*Kt_1*(Tt_1-Tt_2)/(dt_1 +dt_2)       )
!
!     dTt_k(t) / dt = 1/(dt_k*Ct_k) * (- 2*Kt_k-1*(Tt_k-Tt_k-1)/(dt_k-1 +dt_k) 
!                                      - 2*Kt_k  *(Tt_k-Tt_k+1)/(dt_k+1 +dt_k) )
!
!       with
!
!       K*_k  = (d*_k+ d*_k+1)/(d*_k/k*_k+ d*_k+1/k*_k+1)
!
!       Rn = (dir_Rg + sca_Rg) (1-a) + emis * ( Rat - sigma Ts**4 (t+dt) )
!
!       H  = rho Cp CH V ( Ts (t+dt) - Tas )
!
!       LE = rho Lv CH V ( qs (t+dt) - qas )
!
!      where the as subscript denotes atmospheric values at ground level
!      (and not at first half level)
!
!      The tridiagonal systel is linearized with
!
!       using      Ts**4 (t+dt) = Ts**4 (t) + 4*Ts**3 (t) * ( Ts(t+dt) - Ts(t) )
!
!       and  qs (t+dt) = Hu(t) * qsat(t) + Hu(t) dqsat/dT * ( Ts(t+dt) - Ts(t) )
!
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODD_CST
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!	V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    23/01/98 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XCPD, XLVTT, XSTEFAN
!
USE MODE_THERMOS
!
USE MODI_TRIDIAG_GROUND
!
implicit none
!!!#include <arch_specific.hf>
!
!*      0.1    declarations of arguments
!
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_ROOF      ! roof layers temperatures
REAL, DIMENSION(:), INTENT(INOUT) :: PQSAT_ROOF     ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PTA            ! air temperature at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PQA            ! air specific humidity
                                                    ! at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PPS            ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PLW            ! atmospheric infrared radiation
REAL,                 INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_ROOF     ! roof emissivity
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_ROOF       ! heat capacity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_ROOF       ! thermal conductivity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_ROOF        ! depth of roof layers
REAL, DIMENSION(:), INTENT(IN)    :: PTI_BLD        ! inside building temp.
REAL, DIMENSION(:), INTENT(IN)    :: PAC_BLD        ! aerodynamical resistance
                                                    ! inside building itself
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_ROOF     ! fraction of water
REAL, DIMENSION(:), INTENT(IN)    :: PDELTSNOW_ROOF ! roof snow fraction
REAL, DIMENSION(:), INTENT(IN)    :: PGSNOW_ROOF    ! roof snow conduction
!                                                   ! heat fluxes at mantel
!                                                   ! base
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA          ! air density
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROOF       ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROOF_WAT   ! aerodynamical conductance (for water)
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROOF   ! absorbed solar radiation
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_ROOF   ! absorbed infra-red rad.
!
!*      0.2    declarations of local variables
!
REAL :: ZIMPL = 0.5        ! implicit coefficient
REAL :: ZEXPL = 0.5        ! explicit coefficient
!
REAL, DIMENSION(SIZE(PTA)) :: ZDELTFREE_ROOF ! snow-free fraction
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROOF,2)) :: ZA ! lower diag.
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROOF,2)) :: ZB ! main  diag.
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROOF,2)) :: ZC ! upper diag.
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROOF,2)) :: ZY ! r.h.s.
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROOF,2)) :: ZX ! solution
!
REAL, DIMENSION(SIZE(PTA)) :: ZRHO_AC_ROOF ! conductance * rho
REAL, DIMENSION(SIZE(PTA)) :: ZRHO_AC_ROOF_WAT ! ! conductance * rho (for water)
REAL, DIMENSION(SIZE(PTA)) :: ZDQSAT_ROOF  ! dq_sat/dTs
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROOF,2)) :: ZMTC_O_D
! mean thermal conductivity over distance between 2 layers
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROOF,2)) :: ZHC_D_ROOF
! thermal capacity times layer depth
REAL, DIMENSION(SIZE(PTA)) :: ZTS_ROOF
! roof surface temperature
!
INTEGER :: IROOF_LAYER           ! number of roof layers
INTEGER :: JLAYER                ! loop counter
!-------------------------------------------------------------------------------
!
PABS_LW_ROOF(:) = 0.
ZDELTFREE_ROOF=0.
!
ZDELTFREE_ROOF=1.-PDELTSNOW_ROOF
!
!*      1.     Layer thermal properties
!              ------------------------
!
IROOF_LAYER = SIZE(PT_ROOF,2)
ZMTC_O_D(:,:) = 0.
ZHC_D_ROOF(:,:) = 0.
!
DO JLAYER=1,IROOF_LAYER-1
  ZMTC_O_D(:,JLAYER) = 2./(   PD_ROOF(:,JLAYER  )/PTC_ROOF(:,JLAYER  ) &
                            + PD_ROOF(:,JLAYER+1)/PTC_ROOF(:,JLAYER+1) )
  ZHC_D_ROOF(:,JLAYER) = PHC_ROOF(:,JLAYER) * PD_ROOF(:,JLAYER)
END DO
!
ZMTC_O_D(:,IROOF_LAYER) = 2. * PTC_ROOF(:,IROOF_LAYER) &
                             / PD_ROOF (:,IROOF_LAYER)
ZMTC_O_D(:,IROOF_LAYER) = 1./(  1./ZMTC_O_D(:,IROOF_LAYER)      &
                              + 1./(XCPD*PRHOA(:)*PAC_BLD(:))   )
!
ZHC_D_ROOF(:,IROOF_LAYER) = PHC_ROOF(:,IROOF_LAYER) &
                          * PD_ROOF (:,IROOF_LAYER)
!
!-------------------------------------------------------------------------------
!
!*      2.     Roof Ts coefficients
!              --------------------
!
!
ZRHO_AC_ROOF    (:) = PRHOA(:) * PAC_ROOF    (:)
ZRHO_AC_ROOF_WAT(:) = PRHOA(:) * PAC_ROOF_WAT(:)
!
ZTS_ROOF(:) = PT_ROOF(:,1)
!
!*      2.1    dqsat/dTs, and humidity for roofs
!              ---------------------------------
!
ZDQSAT_ROOF(:) = DQSAT(ZTS_ROOF(:),PPS(:),PQSAT_ROOF(:))
!
!
!*      2.2    coefficients
!              ------------
!
ZA(:,1) =   0.

ZB(:,1) =   ZHC_D_ROOF(:,1)/PTSTEP                                           &
  + ZIMPL * ( ZDELTFREE_ROOF(:) * (                                          &
               4. * PEMIS_ROOF(:) * XSTEFAN * ZTS_ROOF(:)**3                 &
             + ZRHO_AC_ROOF(:) * XCPD                                        &
             + ZRHO_AC_ROOF_WAT(:) * XLVTT * PDELT_ROOF(:) * ZDQSAT_ROOF(:)  &
                                  )                                          &
             + ZMTC_O_D(:,1)                                                 &
            )

ZC(:,1) =                                                                    &
    ZIMPL * ( - ZMTC_O_D(:,1)                                                &
            )
!
ZY(:,1) =   ZHC_D_ROOF(:,1)/PTSTEP * PT_ROOF(:,1)                            &
          + ZDELTFREE_ROOF(:) * (                                            &
          + PABS_SW_ROOF(:) + PEMIS_ROOF(:)*PLW(:)                           &
          + ZRHO_AC_ROOF(:) * XCPD * PTA(:)                                  &
          - ZRHO_AC_ROOF_WAT(:) * XLVTT * PDELT_ROOF(:) * PQSAT_ROOF(:)      &
          + ZRHO_AC_ROOF_WAT(:) * XLVTT * PDELT_ROOF(:) * PQA(:)             &
                                )                                            &
          + PDELTSNOW_ROOF(:) * PGSNOW_ROOF(:)                               &
  + ZIMPL * (   ZDELTFREE_ROOF(:) * (                                        &
                  3. * PEMIS_ROOF(:) * XSTEFAN * ZTS_ROOF(:)** 4             &
                + ZRHO_AC_ROOF_WAT(:) * XLVTT * PDELT_ROOF(:)                &
                                  * ZDQSAT_ROOF(:) * ZTS_ROOF(:)             &
                                    )                                        &
            )                                                                &
  + ZEXPL * (   ZDELTFREE_ROOF(:) * (                                        &
              - PEMIS_ROOF(:) * XSTEFAN * ZTS_ROOF(:)** 4                    &
              - ZRHO_AC_ROOF(:) * XCPD * ZTS_ROOF(:)                         &
                                    )                                        &
              - ZMTC_O_D(:,1) * PT_ROOF(:,1)                                 &
              + ZMTC_O_D(:,1) * PT_ROOF(:,2)                                 &
            )
!
!-------------------------------------------------------------------------------
!
!*      3.     Other layers coefficients
!              -------------------------
!
DO JLAYER=2,IROOF_LAYER-1
  ZA(:,JLAYER) =                                                       &
         ZIMPL * ( - ZMTC_O_D(:,JLAYER-1)                              &
                 )

  ZB(:,JLAYER) =   ZHC_D_ROOF(:,JLAYER)/PTSTEP                         &
       + ZIMPL * (   ZMTC_O_D(:,JLAYER-1)                              &
                   + ZMTC_O_D(:,JLAYER  )                              &
                 )

  ZC(:,JLAYER) =                                                       &
         ZIMPL * ( - ZMTC_O_D(:,JLAYER  )                              &
                   )
!
  ZY(:,JLAYER) =   ZHC_D_ROOF(:,JLAYER)/PTSTEP * PT_ROOF(:,JLAYER)     &
       + ZEXPL * (   ZMTC_O_D(:,JLAYER-1) * PT_ROOF(:,JLAYER-1)        &
                   - ZMTC_O_D(:,JLAYER-1) * PT_ROOF(:,JLAYER  )        &
                   - ZMTC_O_D(:,JLAYER  ) * PT_ROOF(:,JLAYER  )        &
                   + ZMTC_O_D(:,JLAYER  ) * PT_ROOF(:,JLAYER+1)        &
                 )
END DO
!
!-------------------------------------------------------------------------------
!
!*      4.     Inside layer coefficients
!              -------------------------
!
ZA(:,IROOF_LAYER) =                                                        &
            ZIMPL * ( - ZMTC_O_D(:,IROOF_LAYER-1)                          &
                    )

ZB(:,IROOF_LAYER) =   ZHC_D_ROOF(:,IROOF_LAYER) / PTSTEP                   &
          + ZIMPL * (   ZMTC_O_D(:,IROOF_LAYER-1)                          &
                      + ZMTC_O_D(:,IROOF_LAYER  )                          &
                    )

ZC(:,IROOF_LAYER) =   0.
!
ZY(:,IROOF_LAYER) =     ZHC_D_ROOF(:,IROOF_LAYER)/PTSTEP                   &
                                                * PT_ROOF(:,IROOF_LAYER)   &
                      + ZMTC_O_D(:,IROOF_LAYER) * PTI_BLD(:)               &
          + ZEXPL * (   ZMTC_O_D(:,IROOF_LAYER-1)                          &
                       * PT_ROOF(:,IROOF_LAYER-1)                          &
                      - ZMTC_O_D(:,IROOF_LAYER-1)                          &
                       * PT_ROOF(:,IROOF_LAYER  )                          &
                      - ZMTC_O_D(:,IROOF_LAYER  )                          &
                       * PT_ROOF(:,IROOF_LAYER  )                          &
                    )
!
!-------------------------------------------------------------------------------
!
!*      5.     Tri-diagonal system resolution
!              ------------------------------
!
CALL TRIDIAG_GROUND(ZA,ZB,ZC,ZY,ZX)
!
PT_ROOF(:,:) = ZX(:,:)
!
!-------------------------------------------------------------------------------
!
!*      6.     Infra-red radiation absorbed by roofs
!              -------------------------------------
!
PABS_LW_ROOF(:) = PEMIS_ROOF(:) * (PLW(:) - XSTEFAN * PT_ROOF(:,1)** 4)
!
!-------------------------------------------------------------------------------
!
!*      7.     New saturated specified humidity near the roof surface
!              ------------------------------------------------------
!
!
PQSAT_ROOF(:) =  QSAT(PT_ROOF(:,1),PPS(:))
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ROOF_LAYER_E_BUDGET
