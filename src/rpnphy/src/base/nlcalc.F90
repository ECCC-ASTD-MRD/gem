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
!
      SUBROUTINE NLCALC( FNN, FRAC, FNGAUSS, FNNONLOC, &
                         WB_NG, WTHL_NG, WQW_NG, WCLD_PROF, &
                         U, V, UW_NG, VW_NG, F_CS, &
                         T, QV, QC, FRV, Z0M, FB_SURF, HPAR, TAU, &
                         Z,ZE,PS,S,N,M,NK,NCLD)
      use tdpack
      implicit none
!
      INTEGER M,N,NK,NCLD
      REAL WCLD_PROF(N,NK,NCLD)
      REAL WB_NG(N,NK), WTHL_NG(N,NK), WQW_NG(N,NK)
      REAL U(M,NK), V(M,NK), UW_NG(N,NK), VW_NG(N,NK)
      REAL T(M,NK),QV(M,NK),QC(N,NK),FNN(N,NK),FRAC(N,NK),FNGAUSS(N,NK),FNNONLOC(N,NK)
      REAL Z(N,NK),ZE(N,NK), S(N,NK), PS(N)
      REAL FB_SURF(N), Z0M(N), FRV(N), HPAR(N)
      REAL F_CS(N,NK)
      REAL TAU
!
!Author
!           A.Lock (May 2004)
!
!Object
!           Calculate the terms pertaining to non-local mixing from cumulus clouds
!           Ref: Lock and Mailhot, 2006: Combining non-local scalings
!            with a TKE closure for mixing in boundary-layer clouds. BLM 121, 313-338.
!           (based on Grant and Brown 1999, QJRMS 125, 1913-1935;
!            Grant 2001, QJRMS 127, 407-421; and Grant and Lock 2004,
!            QJRMS 130, 401-422).
!
!Arguments
!
!           - Output -
! WCLD_PROF profiles of cloud-layer velocity scale (added to EN^0.5 in 
!           the LH04 stable length-scale formulation)
!           Element 1 is for LH, 2 for LM, 3 for LEPS 
! WB_NG    non-gradient part of the buoyancy flux
! WTHL_NG  non-gradient part of the theta_l flux
! WQW_NG   non-gradient part of the qw flux
! UW_NG    non-gradient part of the uw flux
! VW_NG    non-gradient part of the vw flux
! F_CS     factor in stable mixing length coefficient
! FRAC     cloud fraction
! FNNONLOC non-local cumulus cloud fraction
!
!          - Input/Output -
! HPAR     Parcel top height (IN: from previous timestep)
!
!          - Input -
! T        temperature
! QV       specific humidity
! QC       condensed water 
! FNN      buoyancy flux enhancement factor
! FB_SURF  surface buoyancy flux
! Z0M      Surface roughness length for momentum
! FRV      friction velocity
! Z        vertical levels
! ZE       staggered vertical levels (for energy)
! S        sigma levels
! PS       surface pressure
! TAU      model timestep
! N        horizontal dimension
! NK       vertical dimension
! FNGAUSS  Gaussian part of FRAC
!
!
!IMPLICITS
!
!***********************************************************************
      REAL A_PARCEL,B_PARCEL,MAX_T_GRAD
      PARAMETER ( &
       A_PARCEL=0.2 &
      ,B_PARCEL=3.26 &
      ,MAX_T_GRAD=1.0E-3)
!***********************************************************************
! local variables
!
       LOGICAL &
        TOPINV(N) &! indicates top of inversion being reached
       ,ABOVE_LCL  ! indicates being above the LCL
!
      INTEGER J,K &
      ,K_PAR(N)   &! level for start of parcel ascent
      ,KTPAR(N)   &! highest theta level below inversion (at HPAR)
      ,KTINV(N)   &! top level of inversion
      ,K_LCL(N)   &! level of lifting condensation
      ,K_TRT      &! top of cloud-base transition zone
      ,INTERP_INV &! flag to interpolated inversion heights (1=yes)
      ,TOPBL(N)    ! 1 => top of bl reached 
                   ! 2 => max allowable height reached
!
      REAL &
        PRESS(N,NK)     &! pressure
       ,EXNER(N,NK)     &! sigma^cappa
       ,TH(N,NK)        &! potential temperature
       ,THL(N,NK)       &! liquid water potential tempeature
       ,QW(N,NK)        &! total water spec humidity
       ,THVL(N,NK)      &! virtual thl
       ,THV(N,NK)       &! virtual th
       ,THV_PAR(N,NK)    ! virtual th for parcel
      REAL &
        THL_PAR(N)      &! parcel thl
       ,QW_PAR(N)       &! parcel qw
       ,TH_REF(N)       &! reference theta for parcel ascent
       ,TH_PAR_KP1(N)   &! parcel theta at level below
       ,P_LCL(N)        &! pressure of LCL
       ,THV_PERT(N)     &! threshold for parcel thv
       ,Z_LCL(N)        &! height of LCL
       ,ZH(N)           &! boundary layer depth (from RI)
       ,HPAR_OLD(N)     &! height of cloud-top on previous timestep
       ,HPAR_MAX(N)     &! Maximum allowed height for cloud-top
!                        ! (to limit growth rate of boundary layer)
       ,W_STAR(N)       &! sub-cloud layer velocity scale (m/s)
       ,DTHV_MAX(N)     &! maximum parcel thv excess in cloud layer
       ,CAPE(N)         &! convective available potential energy (m2/s2)
       ,DBDZ_INV(N)     &! buoyancy gradient across inversion
       ,DZ_INV_CU(N)     ! inversion thickness above Cu
!
      REAL  &
        VIRT_FACTOR     &
       ,Z_SURF          &! height of surface layer
       ,W_S             &! velocity scale
       ,THV_SD          &! standard deviation of thv in surface layer
       ,QSAT,DQSATDT,QSATFAC &! saturation coefficients in buoyancy parameters
       ,T_REF           &! reference temperature
       ,QC_PAR,QC_ENV   &! parcel and environment liquid water
       ,VAP_PRESS       &! Vapour pressure.
       ,T_LCL           &! temperature of LCL
!!$       ,T_PAR           &! T of parcel
       ,TH_PAR          &! theta of parcel
       ,DPAR_BYDZ       &! parcel thv gradient
       ,DENV_BYDZ       &! environement thv gradient
       ,GAMMA_FA        &! free atmospheric lapse rate
!!$       ,GAMMA_INV       &! inversion lapse rate
       ,GAMMA_CLD       &! cloud layer lapse rate
       ,WB_SCALE        &! buoyancy flux scaling (m2/s3)
       ,W_CLD           &! cloud layer velocity scale (m/s)
       ,V_SUM           &! cubic average of velocity scales (m/s)
       ,Z_CLD           &! cloud layer depth (m)
       ,DZ_INV_CU_REC   &! reconstructed inversion thickness above Cu
       ,VSCALSQ_INCLD   &! incloud squared velocity scale
       ,M_BASE          &! cloud base mass flux (m/s)
       ,Z_PR,ZE_PR      &! scaled height
       ,ZPR_TOP         &! inversion top in scaled coordinate
       ,GRAD            &! gradient of F_ng at z'=0.5
       ,A,B,C           &! coefficients of cubic shape function
       ,F_NG,F_NGE      &! non-gradient shape functions
       ,WTHL_LCL,WQW_LCL & ! fluxes at LCL
       ,UW_LCL,VW_LCL   &! cloud base stresses
       ,Z_TRANS_TOP     &! top of cloud base transition zone
       ,THL_TRT,QW_TRT  &! theta_l and QW at Z_TRANS_TOP
       ,U_TRT,V_TRT     &! U and V at Z_TRANS_TOP
       ,Z0,Z1,Z2,Z3     &! heights for polynomial interpolation
       ,D0,D1,D2,D3     &! values for polynomial interpolation 
       ,A2,A3,XI        &! work variables for polynomial interpolation 
       ,A_POLY,B_POLY,C_POLY &! coefficients in polynomial interpolation 
       ,RI              &! Richardson number
       ,C_STAB_H,C_STAB_M,C_STAB_D &! constants for stable mixing length
       ,GRID_INT        &! THV integral over inversion
       ,ZHDISC          &! height of subgrid interpolated inversion 
       ,ZTINV           &! height of top of inversion
       ,THV_TINV,THV_TPAR! theta_v at top and base of inversion
!*********************************************************

 ! Initialize non-gradient terms
      wb_ng = 0.
      wthl_ng = 0.
      wqw_ng = 0.
      uw_ng = 0.
      vw_ng = 0.

      DO J = 1, N
        HPAR_OLD(J) = HPAR(J)
!       ! Set maximum height of these "shallow" cu to 4km
!        HPAR_MAX(J) = 4000.0 ! Z(J,2)  ! if not restricting growth rate
!       ! Limit boundary layer growth rate to 0.14 m/s (approx 500m/hour)
!       ! as a crude way to allow for turbulence spin-up
        HPAR_MAX(J) = MAX( 200., MIN( 4000.0, HPAR_OLD(J)+TAU*0.14 ) )
        ZH(J) = 1000.0   ! default value only used if RI<RIC everywhere!
        TOPBL(J) = 0
      END DO
      DO K = 1, NK
      DO J = 1, N
        FNNONLOC(J,K) = 0.0   ! initialise cumulus cloud fraction to zero
        PRESS(J,K) = PS(J)*S(J,K)
        EXNER(J,K) = S(J,K)**CAPPA
        TH(J,K) = T(J,K)/EXNER(J,K)
        THL(J,K)= TH(J,K)- (CHLC/CPD)*QC(J,K)/EXNER(J,K)
        QW(J,K) = QV(J,K) + QC(J,K)
        THVL(J,K)= THL(J,K) * ( 1. + DELTA*QW(J,K) )
        VIRT_FACTOR = 1. + DELTA*QV(J,K) - QC(J,K)
        THV(J,K) = TH(J,K) * VIRT_FACTOR 
        THV_PAR(J,K) = THV(J,K)   ! default for stable bls
      END DO
      END DO
      DO K = NK-1,1,-1
      DO J = 1, N
          RI = (U(J,K)-U(J,K+1))**2+(V(J,K)-V(J,K+1))**2
          RI = (GRAV*(Z(J,K)-Z(J,K+1))*(THV(J,K)-THV(J,K+1))/THV(J,K)) / &
                MAX( 1.E-14, RI )
        IF ( RI .GT. RIC .AND. TOPBL(J) .EQ. 0 ) THEN
          ZH(J)    = Z(J,K+1)
          TOPBL(J) = 1
        ENDIF
      END DO
      END DO
!-----------------------------------------------------------------------
! 1. Set up parcel 
!-----------------------------------------------------------------------
! Start parcel ascent from grid-level above top of surface layer, taken
! to be at a height, z_surf, given by 0.1*ZH
!-----------------------------------------------------------------------
      DO J = 1, N
        K_PAR(J) = NK
        K_LCL(J) = 1
        HPAR(J) = ZH(J)  ! initialise to bl depth (from RI)
        IF (FB_SURF(J) .GT. 0.0) THEN
         Z_SURF = 0.1 * ZH(J)
         DO WHILE ( Z(J,K_PAR(J)) .LT. Z_SURF .AND. &
!                   ! not reached Z_SURF
                    THVL(J,K_PAR(J)-1) .LE. THVL(J,K_PAR(J)) )
!                   ! not reached inversion
            K_PAR(J) = K_PAR(J) - 1
         ENDDO
         THL_PAR(J) = THL(J,K_PAR(J))
         QW_PAR(J)  = QW(J,K_PAR(J))
        ENDIF   ! test on unstable
      ENDDO
!-----------------------------------------------------------------------
! Calculate temperature and pressure of lifting condensation level
! using approximations from Bolton (1980)
!-----------------------------------------------------------------------
      DO J = 1, N
        IF (FB_SURF(J) .GT. 0.0) THEN
         VAP_PRESS = QV(J,K_PAR(J)) * &
                     PRESS(J,K_PAR(J)) / ( 100.0*EPS1 )
         IF (VAP_PRESS .GT. 0.0) THEN
           T_LCL = 55.0 + 2840.0 / ( 3.5*ALOG(T(J,K_PAR(J))) &
                                     - ALOG(VAP_PRESS) - 4.805 )
           P_LCL(J) =  PRESS(J,K_PAR(J)) * &
                      ( T_LCL / T(J,K_PAR(J)) )**(1.0/CAPPA)
         ELSE
           P_LCL(J) = PS(J)
         ENDIF
!        ! K_LCL is model level BELOW the lifting condensation level
         K_LCL(J) = NK
         DO WHILE( K_LCL(J).GT.2 .AND. PRESS(J,K_LCL(J)-1).GT.P_LCL(J)) 
           K_LCL(J) = K_LCL(J) -1
         ENDDO
!        ! interpolate to find height of LCL
         Z_LCL(J) = Z(J,K_LCL(J)-1) + ( Z(J,K_LCL(J))-Z(J,K_LCL(J)-1) ) &
                            *( PRESS(J,K_LCL(J)-1) - P_LCL(J)) &
                            /( PRESS(J,K_LCL(J)-1) - PRESS(J,K_LCL(J)) )
         Z_LCL(J) = MAX( 100.0, Z_LCL(J) )
!-----------------------------------------------------------------------
! Threshold on parcel buoyancy for ascent, THV_PERT, is related to 
! standard deviation of thv in surface layer
!-----------------------------------------------------------------------
         W_S = 1.E-10 + ( FB_SURF(J)*ZH(J) + FRV(J)**3 )**(1.0/3.0)
         THV_SD = 1.93 * FB_SURF(J) * THV(J,K_PAR(J)) / ( GRAV * W_S )

         THV_PERT(J)= MAX( A_PARCEL, &
               MIN( MAX_T_GRAD*ZH(J), B_PARCEL*THV_SD ) )

         TH_REF(J) = THL_PAR(J)
         TH_PAR_KP1(J) = THL_PAR(J)
        ENDIF   ! test on unstable
      ENDDO
!-----------------------------------------------------------------------
!! 2  Parcel ascent:
!-----------------------------------------------------------------------
! Lift parcel conserving its THL and QW.
! Calculate parcel QC by linearising q_sat about a reference temperature 
! equal to an estimate of the parcel's temperature found by extrapolating 
! the parcel's potential temperature gradient up from the grid-level below.
      DO K = NK, 1, -1
      DO J = 1, N
        IF (FB_SURF(J) .GT. 0.0) THEN
!!$          T_REF   = TH_REF(J)*EXNER(J,K)
          T_REF   = max(TH_REF(J)*EXNER(J,K),100.)
          QSAT    = FOQSA( T_REF, PRESS(J,K) )
          DQSATDT = FODQA( QSAT, T_REF )
          QSATFAC = 1./(1.+(CHLC/CPD)*DQSATDT)
          QC_PAR  = MAX( 0.0, QSATFAC*( QW_PAR(J) - QSAT  &
                    - (THL_PAR(J)-TH_REF(J))*EXNER(J,K)*DQSATDT ) )

!         ! Add on the difference in the environment's qc as calculated by the
!         ! partial condensation scheme and with the above linear approximation.
!         ! This then imitates partial condensation in the parcel.
          QC_ENV  = MAX( 0.0, QSATFAC*( QW(J,K) - QSAT  &
                    - (THL(J,K)-TH_REF(J))*EXNER(J,K)*DQSATDT ) )
          QC_PAR  = QC_PAR + QC(J,K) - QC_ENV
          TH_PAR  = THL_PAR(J) + (CHLC/CPD)*QC_PAR/EXNER(J,K)
          THV_PAR(J,K) = TH_PAR *  &
                         (1.+DELTA*QW_PAR(J)-(1.+DELTA)*QC_PAR)
          IF (K .GT. 1 .AND. K .LT. NK-1) THEN
!         ! extrapolate reference TH gradient up to next grid-level
            Z_PR      = (Z(J,K-1)-Z(J,K))/(Z(J,K)-Z(J,K+1))
            TH_REF(J) = TH_PAR*(1.0+Z_PR) - TH_PAR_KP1(J)*Z_PR
            TH_PAR_KP1(J) = TH_PAR
          ENDIF

        ENDIF   ! test on unstable
      ENDDO
      ENDDO

!-----------------------------------------------------------------------
!! 3 Identify layer boundaries
!-----------------------------------------------------------------------
      DO J = 1, N
        TOPBL(J) = 0
        TOPINV(J)= .FALSE.
        KTPAR(J) = NK
        KTINV(J) = NK
        DTHV_MAX(J) = 0.0
        DBDZ_INV(J) = 0.003  ! impose a weak lower limit inversion lapse rate 
!                            ! (~1.e-4 s^-2, converted from K/m to s^-2 later)
      ENDDO

      DO K = NK-1, 1, -1
      DO J = 1, N

        IF (FB_SURF(J) .GT. 0.0) THEN
!         !------------------------------------------------------------
!         ! Set flag to true when level BELOW is above the lcl 
!         ! and above LCL transition zone
!         !------------------------------------------------------------
           ABOVE_LCL = K+1 .LT. K_LCL(J) .AND. Z(J,K) .GT. 1.1*Z_LCL(J)
!         !-------------------------------------------------------------
!         ! Calculate vertical gradients in parcel and environment THV
!         !-------------------------------------------------------------
          DPAR_BYDZ = (THV_PAR(J,K) - THV_PAR(J,K+1)) / &
                      (Z(J,K) - Z(J,K+1))
          DENV_BYDZ = (THV(J,K) - THV(J,K+1)) / &
                      (Z(J,K) - Z(J,K+1))
!         !-------------------------------------------------------------
!         ! Find top of inversion - where parcel has minimum buoyancy
!         !-------------------------------------------------------------
          IF ( TOPBL(J) .GT. 0 .AND. .NOT. TOPINV(J) ) THEN 
            DBDZ_INV(J) = MAX( DBDZ_INV(J), DENV_BYDZ )
            IF ( K+1 .LE. KTPAR(J)-2 .AND. ( &
!                ! Inversion at least two grid-levels thick
                 DENV_BYDZ .LE. DPAR_BYDZ .OR.  &
!                ! => at a parcel buoyancy minimum
                 ZE(J,K) .GT. HPAR(J)+MIN(1000., 0.5*HPAR(J)) )) THEN
!                ! restrict inversion thickness < 1/2 bl depth and 1km
              TOPINV(J) = .TRUE.
              KTINV(J) = K+1
            ENDIF
          ENDIF
!         !-------------------------------------------------------------
!         ! Find base of inversion - where parcel has maximum buoyancy
!         !                          or is negatively buoyant
!         !-------------------------------------------------------------
          IF ( TOPBL(J).EQ.0 .AND. K .LT. K_PAR(J) .AND. &
               (  ( THV_PAR(J,K)-THV(J,K) .LE. - THV_PERT(J) ) .OR. &
!
!                      plume non buoyant
!
                  ( ABOVE_LCL .AND. (DENV_BYDZ .GT. 1.25*DPAR_BYDZ) ) &
!
!                      or environmental virtual temperature gradient
!                      significantly larger than parcel gradient
!                      above lifting condensation level
!
               ) &
             ) THEN
!
            TOPBL(J) = 1
            KTPAR(J) = K+1   ! marks most buoyant theta-level 
!                            ! (just below inversion)
            HPAR(J)  = ZE(J,KTPAR(J)-1)
            DBDZ_INV(J) = MAX( DBDZ_INV(J), DENV_BYDZ )
          ENDIF

          IF ( TOPBL(J).EQ.0 .AND. Z(J,K+1) .GE. HPAR_MAX(J) ) THEN
!                      gone above maximum allowed height
            TOPBL(J) = 2
            KTPAR(J) = K+2
            HPAR(J)  = ZE(J,KTPAR(J)-1)
            DBDZ_INV(J) = MAX( DBDZ_INV(J), DENV_BYDZ )
          ENDIF

          IF ( TOPBL(J).EQ.0 .AND. K .LT. K_LCL(J) ) THEN
!           within cloud layer
            DTHV_MAX(J) = MAX( DTHV_MAX(J), THV_PAR(J,K)-THV(J,K) )
          ENDIF

        ENDIF   ! test on unstable
      ENDDO
      ENDDO
!--------------------------------------------------------------------------
!! 3.1 Interpolate inversion base (max buoyancy excess) between grid-levels
!--------------------------------------------------------------------------
      DO J = 1, N
        IF ( KTPAR(J) .LT. NK ) THEN
!         !-----------------------------------------------------
!         ! parcel rose successfully
!         !-----------------------------------------------------
          HPAR(J)    = ZE(J,KTPAR(J)-1)
          INTERP_INV=1
          IF (TOPBL(J) .EQ. 2) THEN
!           ! Stopped at max allowable height
            INTERP_INV= 0
            HPAR(J)  = HPAR_MAX(J)
          ENDIF
          IF ( INTERP_INV .EQ. 1 .AND. KTPAR(J) .GT. 2 ) THEN
!            !-----------------------------------------------------------
!            ! interpolate height by fitting a cubic to 3 parcel excesses,
!            ! at Z1 (top of cloud layer) and the two grid-levels above,
!            ! and matching cloud layer gradient (D1-D0) at Z1.
!            ! Polynomial is f(z) = a*z^3 + b*z2 + cz + d
!            !-----------------------------------------------------------
              K = KTPAR(J)-2
              Z3=Z(J,K)  -Z(J,K+2)
              D3=THV_PAR(J,K)-THV(J,K)
              Z2=Z(J,K+1)-Z(J,K+2)
              D2=THV_PAR(J,K+1)-THV(J,K+1)
              Z1=Z(J,K+2)-Z(J,K+2)
              D1=THV_PAR(J,K+2)-THV(J,K+2)
              Z0=Z(J,K+3)-Z(J,K+2)
              D0=THV_PAR(J,K+3)-THV(J,K+3)
              C_POLY = (D1-D0)/(Z1-Z0)
              A2= D2 - D1 - C_POLY*Z2
              A3= D3 - D1 - C_POLY*Z3
              B_POLY = (A3-A2*(Z3*Z3*Z3)/(Z2*Z2*Z2))/(Z3*Z3*(1.0-Z3/Z2))
              A_POLY = (A2-B_POLY*Z2*Z2)/Z2**3

              XI=B_POLY*B_POLY-3.*A_POLY*C_POLY
              IF (A_POLY .NE. 0 .AND. XI .GT. 0.0) THEN
!               ! HPAR is then the height where the above polynomial
!               ! has zero gradient
                HPAR(J) = Z(J,K+2)-(B_POLY+SQRT(XI))/(3.*A_POLY)
                HPAR(J) = MAX( MIN( HPAR(J), Z(J,K) ), Z(J,K+2) )
                IF ( HPAR(J) .GT. Z(J,KTPAR(J)-1) ) THEN
                  KTPAR(J)=KTPAR(J)-1
                ENDIF
              ENDIF
          ENDIF

        ENDIF   ! parcel rose
      ENDDO
!-----------------------------------------------------------------------
!! 4. Integrate parcel excess buoyancy
!-----------------------------------------------------------------------
      DO J = 1, N
        CAPE(J) = 0.0
        IF (K_LCL(J) .GT. 1 .AND. KTPAR(J) .LT. K_LCL(J)  &
!                                 ! parcel rose above LCL
                            .AND. DTHV_MAX(J) .GT. 0.0 ) THEN
!                                 ! parcel became positively buoyant
          K = K_LCL(J)-1
!         ! assume parcel is neutrally buoyant at LCL
          CAPE(J) = 0.5*(Z(J,K)-Z_LCL(J))* &
                         (THV_PAR(J,K)-THV(J,K))
          IF ( THV_PAR(J,K)-THV(J,K) .LT. 0.0 ) THEN
!              ! Parcel wasn't buoyant just above the LCL
!              ! so adjust CAPE for this negative buoyancy 
!              ! (ie. assume parcel is neutrally buoyant at LCL)
            CAPE(J) = CAPE(J) -  &
                     (HPAR(J)-Z_LCL(J))*(THV_PAR(J,K)-THV(J,K))
          ENDIF
!         ! now calculate CAPE at interpolated top of profile
          K = KTPAR(J)
          CAPE(J) = CAPE(J) + MAX( 0.0,  &
                    (HPAR(J)-Z(J,K))*(THV_PAR(J,K)-THV(J,K)) )
        ENDIF
      ENDDO
      DO K = 1, NK-1
      DO J = 1, N
        IF (K .LE. K_LCL(J)-2 .AND. K .GE. KTPAR(J) ) THEN
!         ! within main cloud-layer
          CAPE(J) = CAPE(J) + 0.5*(Z(J,K)-Z(J,K+1))* &
                       (THV_PAR(J,K)-THV(J,K)+THV_PAR(J,K+1)-THV(J,K+1))
        ENDIF
      ENDDO
      ENDDO
!-----------------------------------------------------------------------
!! 5. Calculate non-gradient fluxes and velocity scales
!-----------------------------------------------------------------------
      WCLD_PROF = 0.
      DO J = 1, N
        DZ_INV_CU(J) = 0.0
        W_STAR(J) = 0.0
        IF (FB_SURF(J) .GT. 0.0) THEN
         W_STAR(J) = ( FB_SURF(J)*HPAR(J) )**(1.0/3.0)  ! dry bl scale
         DBDZ_INV(J) = GRAV*DBDZ_INV(J)/THV(J,KTPAR(J)) ! convert to buoyancy units
         DZ_INV_CU(J) = 0.2*HPAR(J)                     ! default for no CAPE
         IF (CAPE(J) .GT. 0.0) THEN
          K = K_LCL(J)
!         ! calculate cumulus scaling parameters
          W_STAR(J) = ( FB_SURF(J)*Z_LCL(J) )**(1.0/3.0)
          M_BASE    = 0.04*W_STAR(J)
          CAPE(J)   = ( GRAV/THV(J,K) ) * CAPE(J)
          W_CLD     = ( M_BASE * CAPE(J) )**(1./3.)
          V_SUM     = (W_CLD**3+W_STAR(J)**3)**(1./3.)
          Z_CLD     = MAX( 50., HPAR(J) - Z_LCL(J) )
          WB_SCALE  = ( W_CLD**3/Z_CLD ) * SQRT( M_BASE/W_CLD ) !Eq. 17 of LM06
!         ! calculate fluxes at LCL
          Z_TRANS_TOP = 1.1*Z_LCL(J)
          K_TRT = MIN( NK-1, K_LCL(J) )
          DO WHILE (Z(J,K_TRT) .LT. Z_TRANS_TOP .AND. K_TRT .GE. 2) 
!           ! K_TKT is first theta-grid-level above Z_TRANS_TOP
            K_TRT = K_TRT - 1
          ENDDO
          Z_PR = (Z(J,K_TRT)-Z_TRANS_TOP)/(Z(J,K_TRT)-Z(J,K_TRT+1))
!           ! Assuming wchi = - mb*dchi_transition_zone
            THL_TRT  = THL(J,K_TRT) - (THL(J,K_TRT)-THL(J,K_TRT+1))*Z_PR
            QW_TRT   =  QW(J,K_TRT) - ( QW(J,K_TRT)- QW(J,K_TRT+1))*Z_PR
          U_TRT  = U(J,K_TRT) - (U(J,K_TRT)-U(J,K_TRT+1))*Z_PR
          V_TRT  = V(J,K_TRT) - (V(J,K_TRT)-V(J,K_TRT+1))*Z_PR
          WTHL_LCL = M_BASE *( THL(J,K_PAR(J)) - THL_TRT )
          WQW_LCL  = M_BASE *(  QW(J,K_PAR(J)) - QW_TRT )
          UW_LCL   = M_BASE *(   U(J,K_LCL(J)) - U_TRT )  ! use K_LCL for momentum as 
          VW_LCL   = M_BASE *(   V(J,K_LCL(J)) - V_TRT )  ! mixed layer not so well-mixed
!         !----------------------------------------------------------
!         ! Estimate inversion thickness using energy conservation 
!         ! for cumulus parcel with KE scale VSCAL_INCLD encountering 
!         ! a region with stability DBDZ_INV (=maximum stability across
!         ! the model's inversion)
!         !----------------------------------------------------------
          VSCALSQ_INCLD = 2.0*CAPE(J)
          DZ_INV_CU(J)  = SQRT( VSCALSQ_INCLD/DBDZ_INV(J) )
!         ! 
!         ! If inversion is unresolved (less than 3 grid-levels thick)
!         ! then use profile reconstruction 
!         ! 
          IF ( KTPAR(J) .GE. 5 ) THEN
          IF ( DZ_INV_CU(J) .LT. Z(J,KTPAR(J)-3) - Z(J,KTPAR(J)) ) THEN
!           !
!           ! First interpolate to find height of discontinuous inversion
!           !
            K = KTPAR(J)
            GAMMA_CLD = (THV(J,K)-THV(J,K+1))/(Z(J,K)-Z(J,K+1))
            IF (K+2 .LT. K_LCL(J)) THEN
              GAMMA_CLD = MIN( GAMMA_CLD, ( THV(J,K+1)-THV(J,K+2) ) &
                                         /(   Z(J,K+1)-  Z(J,K+2) ) )
            ENDIF
            GAMMA_CLD = MAX(0.0, GAMMA_CLD)
            GAMMA_FA = (THV(J,K-4)-THV(J,K-3))/(Z(J,K-4)-Z(J,K-3))
            GAMMA_FA = MAX(0.0, GAMMA_FA)
!           ! Integrate thv over the inversion grid-levels
            GRID_INT =  (THV(J,K-1)-THV(J,K))*(ZE(J,K-2)-ZE(J,K-1)) &
                      + (THV(J,K-2)-THV(J,K))*(ZE(J,K-3)-ZE(J,K-2)) &
                      + (THV(J,K-3)-THV(J,K))*( Z(J,K-3)-ZE(J,K-3)) 
!           ! Solve quadratic equation for Dz=ZHDISC-Z(K)
            C_POLY = (THV(J,K-3)-THV(J,K))*(Z(J,K-3)-Z(J,K)) &
                     - 0.5*GAMMA_FA*(Z(J,K-3)-Z(J,K))**2 - GRID_INT 
            B_POLY = -(THV(J,K-3)-THV(J,K)-GAMMA_FA*(Z(J,K-3)-Z(J,K)))
            A_POLY = 0.5*(GAMMA_CLD+GAMMA_FA)
            XI     = B_POLY*B_POLY-4.*A_POLY*C_POLY

            IF (XI.GE.0.0.AND.( A_POLY.NE.0.0.OR.B_POLY.NE.0.0 )) THEN
              IF (A_POLY .EQ. 0.0) then 
                DZ_INV_CU_REC = -C_POLY/B_POLY
              ELSE
                DZ_INV_CU_REC = (-B_POLY-SQRT(XI))/(2.*A_POLY)
              ENDIF
              ZHDISC  = Z(J,K)+DZ_INV_CU_REC
!             !
!             ! Now calculate inversion stability given Dz=V^2/DB
!             ! Solve quadratic equation for Dzi
!             !
              C_POLY = -VSCALSQ_INCLD*THV(J,K-1)/GRAV
              B_POLY = THV(J,K-3)-GAMMA_FA *(Z(J,K-3)-ZHDISC) &
                      -THV(J,K)  -GAMMA_CLD*(ZHDISC  -Z(J,K))
              A_POLY = 0.5*(GAMMA_CLD+GAMMA_FA)
              XI=B_POLY*B_POLY-4.*A_POLY*C_POLY

              IF (XI.GE.0.0.AND.( A_POLY.NE.0.0.OR.B_POLY.NE.0.0 )) THEN
                IF (A_POLY .EQ. 0.0) then 
                  DZ_INV_CU_REC = -C_POLY/B_POLY
                ELSE
                  DZ_INV_CU_REC = (-B_POLY+SQRT(XI))/(2.*A_POLY)
                ENDIF
                DZ_INV_CU_REC = MAX(MIN(DZ_INV_CU_REC,2.0*(ZHDISC-Z(J,KTPAR(J)))),1.)
                IF (DZ_INV_CU_REC .LE. DZ_INV_CU(J)) THEN
!                 ! Use reconstructed inversion thickness if it is 
!                 ! less than the estimate from grid-level lapse rate
!                 ! and recalculate inversion stability for F_CS(RI)
                  DZ_INV_CU(J) = DZ_INV_CU_REC
                  ZHDISC  = ZHDISC - 0.5*DZ_INV_CU_REC
                  ZTINV   = ZHDISC + DZ_INV_CU_REC
                  THV_TINV = THV(J,K-3)-GAMMA_FA *(Z(J,K-3)-ZTINV)
                  THV_TPAR = THV(J,K)  +GAMMA_CLD*(ZHDISC -Z(J,K))
                  DBDZ_INV(J) = (GRAV/THV(J,K-1)) *  &
                                (THV_TINV-THV_TPAR)/DZ_INV_CU_REC
                  DBDZ_INV(J) = MAX(0.0, DBDZ_INV(J) )
                ENDIF
              ENDIF  ! interpolation for DZ_INV_CU successful
            ENDIF    ! interpolation for ZHDISC successful

          ENDIF  ! inversion not resolved
          ENDIF  ! KTPAR ge 5
          DZ_INV_CU(J)= MAX( 1., DZ_INV_CU(J) )

          ZPR_TOP     = 1. + MIN(1., DZ_INV_CU(J)/Z_CLD )
          DO K=1,NK-1

!           ! Z_PR=0 at cloud-base, 1 at cloud-top
            Z_PR = ( ZE(J,K) - Z_LCL(J) )/ Z_CLD
            IF ( Z_PR.LE.0.0 ) THEN
! ----------------------------------------------------------------
!   Sub-cloud layer profiles
! ----------------------------------------------------------------
!             ! ZE_PR=0 at surface, 1 at cloud-base
              ZE_PR = ZE(J,K) / Z_LCL(J)
!
!             ! non-gradient conserved variable fluxes
!
              F_NG = ZE_PR
              WTHL_NG(J,K) = WTHL_LCL* F_NG
              WQW_NG(J,K)  = WQW_LCL * F_NG
              UW_NG(J,K)   = UW_LCL * F_NG
              VW_NG(J,K)   = VW_LCL * F_NG
            ELSE 
! ----------------------------------------------------------------
!   Cloud layer profiles
! ----------------------------------------------------------------
              IF (Z_PR .LE. ZPR_TOP) THEN
!
!               ! non-gradient conserved variable fluxes
!
                IF (Z_PR .LE. 1.0) THEN
!                 ! cubic function with df/dz=grad at z=0.5
!                 !                     f=1,0.5,0 at z=0,0.5,1
                  GRAD=-0.5
                  A = -4.*(GRAD+1.)
                  B =  6.*(GRAD+1.)
                  C = -2.*(GRAD+1.) - 1.
                  F_NG = A*Z_PR**3+B*Z_PR**2+C*Z_PR+1.
                  WTHL_NG(J,K) = WTHL_LCL* F_NG 
                  WQW_NG(J,K)  = WQW_LCL * F_NG
                  UW_NG(J,K)   = UW_LCL * F_NG
                  VW_NG(J,K)   = VW_LCL * F_NG
                ENDIF
!
!               ! Velocity scales for stable mixing length calculation
!               ! First WCLD_PROF(1) for heat, humidity
!
                IF (Z_PR .LE. 0.9) THEN
                  F_NG = 1.8*(Z_PR**1.4)*(1.543-Z_PR)  ! function fitted to LES l_h
                ELSE
                  ZE_PR = (Z_PR-0.9)/(ZPR_TOP-0.9)  ! from 0 to 1
                  F_NG = 0.5*(1.+COS(PI*ZE_PR))
                ENDIF
                C_STAB_H = 0.072
                WCLD_PROF(J,K,1) = C_STAB_H * 1.3 * V_SUM * F_NG
                C_STAB_M = 0.15   ! 0.05*Min(3,2*Ri) and assume Ri large
                WCLD_PROF(J,K,2) = C_STAB_M * 0.8 * V_SUM * F_NG
!               ! Second WCLD_PROF(3) for dissipation
                IF (Z_PR .LE. 0.1) THEN
                  F_NGE = Z_PR/0.1
                ELSE IF (Z_PR .LE. 1.0) THEN
                  F_NGE = ( (1.-Z_PR) + 0.8*(Z_PR-0.1) )/0.9
                ELSE
                  F_NGE = 0.8*(ZPR_TOP-Z_PR)/(ZPR_TOP-1.0)
                ENDIF
                C_STAB_D = 0.3   ! 0.1*Min(3,5*Ri) and assume Ri large
                WCLD_PROF(J,K,3) = C_STAB_D * 0.7 * V_SUM * F_NGE 
!
!               ! Non-gradient function for WB (Eq. A3 of LM06)
                IF ( Z_PR .LE. 0.9 ) THEN
!                 ! function with gradient=0 at z=0.9
!                 !                    f=0,1 at z=0,0.9
                  ZE_PR = Z_PR/0.9
                  F_NG  = 0.5 * SQRT(ZE_PR) * (3.0-ZE_PR)
                ELSE
                  ZE_PR = (Z_PR-0.9)/(ZPR_TOP-0.9)  ! from 0 to 1
                  F_NG  = 0.5 * (1.+COS(PI*ZE_PR))
                ENDIF
                
                ! Non-gradient buoyancy flux component (Eq. 16 of LM06)
                WB_NG(J,K) = (1.-FNN(J,K))*3.7*F_NG*WB_SCALE

              ENDIF ! Z_PR le ZPR_TOP

            ENDIF   ! Z_PR gt 0
!           !--------------------------------------
!           ! Cloud fraction enhancement 
!           ! (on Z rather than ZE levels)
!           !--------------------------------------
            Z_PR = ( Z(J,K) - Z_LCL(J) )/ Z_CLD
!           ! Z_PR=0 at cloud-base, 1 at cloud-top
!
            IF (Z_PR .GT. 0.0) THEN
              F_NG = 0.0
              IF ( Z_PR .LE. 0.9 ) THEN
                 F_NG = 1.0+3.0*EXP(-5.0*Z_PR)     ! =4 at cloud-base
              ELSE IF ( Z_PR .LT. ZPR_TOP ) THEN
                 ZE_PR = (Z_PR-0.9)/(ZPR_TOP-0.9)  ! from 0 to 1
                 F_NG  = 0.5*(1.+COS(PI*ZE_PR))
              ENDIF
              FNNONLOC(J,K) = 0.5*F_NG*MIN(0.5,M_BASE/W_CLD)
!             ! FRAC_GAUSS is cloud fraction calculated by Bechtold scheme.
!             ! Combine cloud fractions (using maximum) for Ri calculation 
!             ! in blcloud 
              FRAC(J,K) = MAX( FNGAUSS(J,K), FNNONLOC(J,K) )
            ENDIF   ! Z_PR le 0

          ENDDO   ! loop over K
         ENDIF ! Test on CAPE
        ENDIF ! Test on FB_SURF
      ENDDO
!----------------------------------------------------------------------------
!! 6. Calculate factor for c_s in stable length scale formulation:
!       0 = Mailhot and Lock "small" values
!       1 = Lenderink and Holtslag (larger) originals
!----------------------------------------------------------------------------
      DO J = 1, N
!       ! For stable/neutral, try F_CS a function of surface inhomogenity.
!       ! Using Z0M as a proxy.  Transition on log scale 
!       ! from f_cs~0 for z0<0.002 to f_cs~2 for z0>0.5
!       ! with f_cs = 1 at z0=0.1m, 0.5 at z0=0.05m

        F_CS(J,NK) = 1.+TANH(ALOG(Z0M(J))+2.3)  ! ALOG(0.1)=-2.3

!       ! Add stability dependence to remove effect for bls with zi/L <-2
        F_NG = KARMAN*W_STAR(J)**3/(1.E-10+FRV(J)**3)
        F_CS(J,NK) = F_CS(J,NK)*MAX(0.0, 1.- 0.5*F_NG )

!        F_CS(J,NK) = 0.0          ! Or use small values for c_s

        DO K=1,NK-1
          Z_PR = EXP(-ZE(J,K)/500.)
          F_CS(J,K) = Z_PR*F_CS(J,NK)
        ENDDO
        IF (FB_SURF(J) .GT. 0.0) THEN
!         ! Make XI ~1 for weak inversions, ~0 for more stable (Sc)
          RI = HPAR(J)*SQRT(DBDZ_INV(J))/ &
               (MAX(0.0,2.0*CAPE(J))**1.5 +W_STAR(J)**3)**(1./3.)
          XI = 0.5*(1.0-TANH((RI-12.)/4.))
          DO K=1,NK
!           ! Linear increase/decrease with height across inversion 
             Z_PR = 0.0
             IF ( ZE(J,K) .LE. HPAR(J) ) THEN
               Z_PR = ZE(J,K)/HPAR(J)
             ELSE IF (ZE(J,K) .LE. HPAR(J)+DZ_INV_CU(J) ) THEN
               Z_PR = 1.
             ELSE IF (ZE(J,K) .LT. HPAR(J)+2.*DZ_INV_CU(J) ) THEN
               Z_PR = (HPAR(J)+2.*DZ_INV_CU(J)-ZE(J,K))/DZ_INV_CU(J)
             ENDIF
             F_CS(J,K) = MAX( F_CS(J,K), &   ! factor of 2 to give extra boost to bl top mixing
                              2.0*Z_PR*XI )  ! particularly for cloud-free surface-heated bls
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END
