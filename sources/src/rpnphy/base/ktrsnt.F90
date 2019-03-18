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
!**S/P KTRSNT3
!
      SUBROUTINE KTRSNT3 ( CTT,CQT,ilab,CCF,QCKTL,QCKTI,DBDT,CRRL,CRRI, &
                          QCN,TP,TM,QP,QM,GZM,TQDF,PSP,PSM, &
                          SIGMA, TAU, KSHAL, NI, NK)
      use tdpack
      implicit none
#include <arch_specific.hf>
!
      INTEGER NI,NK
      REAL CTT(NI,NK),CQT(NI,NK)
      INTEGER ilab(NI,NK)
      REAL CCF(NI,NK),QCKTL(NI,NK),QCKTI(NI,NK),DBDT(NI),TQDF(NI,NK)
      REAL TP(NI,NK),TM(NI,NK),QP(NI,NK),QM(NI,NK),GZM(NI,NK)
      REAL PSP(NI),PSM(NI),SIGMA(NI,NK),CRRL(NI),CRRI(NI),QCN(NI,NK)
      REAL TAU
      REAL KSHAL(NI)
!
!Authors
!          Claude Girard and Gerard Pellerin 1995
!
!Revision
! 001      G.Pellerin (Nov 98) Added kuo65 option
!                      and change accession closure
! 002      C.Girard   (Dec 98) Added detrainement
!
! 003      A-M.Leduc  (March 2002) Automatic arrays
!
! 004      S.Belair, A-M. Leduc (nov 2002) added convective counter kshal
!                         as argument. ktrsnt--->ktrsnt2
! 005      G.Pellerin (Avril 2003) -  CVMG... Replacements
! 006      G. Pellerin (Mai 03) - Conversion IBM
!                     - calls to vexp routine (from massvp4 library)
!                     - calls to optimized routine MFOQST
!
! 005      S.Belair   (April 2003) Outputs for in-cloud water (QCKT)
!
! 006      S.Belair, A-M. Leduc (Mai 2003) calculation of ice fraction (ZFICE)
!                    and of saturation with respect to water and ice in
!                    temperature interval MAXFRZ and MINFRZ. Output
!                    QCKTL and QCKTI.
! 007      A. Plante, C. Girard (Feb 2004) KTRSNT2 -> KTRSNT3 :
!                    Add precipitation and evaporation of precipitation under
!                    clouds. Output CRRL, CRRI and QCN.
! 008      B. Bilodeau (May 2005) - QCKTL and QCKTI average values
!                                   instead of in-cloud values
!
!
!Object
!          To calculate the convective tendencies of T and Q
!          using a scheme with a "Kuo65-type closure".
!          Geleyn's method is used to obtain the cloud profiles.
!
!Arguments
!
!            - Outputs -
! CTT      convective temperature tendency
! CQT      convective specific humidity tendency
! ilab     flag array: an indication of convective activity
! CCF      estimated cumulus cloud fraction
! QCKTL    estimated cumulus cloud liquid water for radiation only
! QCKTI    estimated cumulus cloud solid water for radiation only
! DBDT     estimated averaged cloud fraction growth rate
! CRRL     Precipitation rate (liquid)
! CRRI     Precipitation rate (solid)
! QCN      estimated cumulus total cloud liquid for precipitation computation only
!          note that QCN is not equal to QCKTL + QCKTI
!            - Inputs -
! TP       temperature at (t+dt)
! TM       temperature at (t-dt)
! QP       specific humidity at (t+dt)
! QM       specific humidity at (t-dt)
! GZM      geopotential
! TQDF     tendance diffusive de couche limite (t+dt)
! PSP      surface pressure at (t+dt)
! PSM      surface pressure at (t-dt)
! SIGMA    sigma levels
! TAU      effective timestep (2*dt)
! NI       horizontal dimension
! NK       vertical dimension
!
!Notes
!          The routine is divided into 5 parts:
!           1)allocation and position for work space
!           2)preliminary computations
!           3)cloud ascent and flagging
!           4)total moisture accession calculations
!           5)cloud heating and moistening (drying) calculations
!
!           STPEVP      EVAPORATION COEFFICIENT FOR STRATIFORM PRECIPITATION
!           HDPMX       MAXIMUM PRECIPITATION CHANGE DUE TO EVAPORATION
!           PFRAC       FRACTION OF TILE THAT IS PRECIPITATING
!*
      LOGICAL LO
      INTEGER IS,IKA,IKB,jk,jkm1,jl,MODP,stat
      INTEGER, DIMENSION(NI) :: KMIN
      INTEGER, external :: neark
      REAL ZTVC,rgrav3,rcpd,rcpv
      REAL ENTRM,TAUCU,DELTA2,DZETR
      REAL ZCOR,ZQSC,DETRN,ZQ0,ZT0,ZN0
      REAL ZQCD,ZK,ZDH,temp1,temp2, temp3
      REAL MAXFRZ, MINFRZ, DGZFIX
      REAL PRCP, XN, STPEVP, HKPEVP, HCIMPR, XX, YY, HDPMX, XP, XDP
      REAL PFRAC, ZPC
!
!*********************************************************
!     AUTOMATIC ARRAYS
!*********************************************************
!
      REAL, dimension(NI,NK) :: ZPP
      REAL, dimension(NI   ) :: ZPT
      REAL, dimension(NI,NK) :: ZDSG
      REAL, dimension(NI,NK) :: ZDP
      REAL, dimension(NI,NK) :: ZSDP
      REAL, dimension(NI,NK) :: ZQAC
      REAL, dimension(NI,NK) :: ZSQAC
      REAL, dimension(NI,NK) :: ZLDCP
      REAL, dimension(NI,NK) :: ZQSE
      REAL, dimension(NI,NK) :: ZTC
      REAL, dimension(NI,NK) :: ZQC
      REAL, dimension(NI   ) :: ZQCNT
      REAL, dimension(NI,NK) :: ZTE
      REAL, dimension(NI,NK) :: ZQE
      REAL, dimension(NI,NK) :: ZTVE
      REAL, dimension(NI,NK) :: ZDQ
      REAL, dimension(NI,NK) :: ZSDQ
      REAL, dimension(NI,NK) :: ZDT
      REAL, dimension(NI,NK) :: ZSDT
      REAL, dimension(NI,NK) :: ZSDH
      REAL, dimension(NI,NK) :: ZFICE
      REAL, dimension(NI,NK) :: DFMX
      REAL, dimension(NI   ) :: ZCP
      REAL, dimension(NI   ) :: ZLDCP0
      REAL, dimension(NI   ) :: CPR
      LOGICAL, dimension(NI   ) :: LO1

!
!****************************************************
!
!*    PHYSICAL CONSTANTS.
!     -------- ----------

      rcpd = 1./CPD
      rgrav3 = 1./(GRAV*1.E3)
!     typiquement entraine selon labda=1/GH (H=2km)
      DZETR = 2.E+03
      ENTRM = 1./(DZETR*GRAV)
      TAUCU = 900.
      DELTA2 = CPV/CPD - 1.
      DGZFIX=200.*GRAV
!      HKPEVP = 4.15E-4
      HKPEVP = 2.E-4
      STPEVP = 2. * GRAV * HKPEVP
!
!     ------------------------------------------------------------------
!
!*         1.     ALLOCATION AND POSITION FOR WORK SPACE.
!                 ---------- --- -------- --- ---- ------
!

!***
!
!     METHOD.
!     -------
!
!          IN (3) A NEARLY ADIABATIC ASCENT IS ATTEMPTED FOR A CLOUD
!     PARCEL STARTING FROM THE LOWEST MODEL LAYER. THIS CLOUD ASCENT
!     IS COMPUTED IN TERMS OF TEMPERATURE AND SPECIFIC HUMIDITY.
!     ENTRAINMENT IS SIMULATED VIA AN ENTRAINMENT PARAMETER.
!     THE LAYERS ARE FLAGGED ACCORDING TO THE FOLLOWING CODE:
!     0 = STABLE OR INACTIVE LAYER,
!     1 = PART OF THE WELL MIXED BOUNDARY LAYER OR DRY UNSTABLE LAYER,
!     2 = MOIST UNSTABLE OR ACTIVE OR CLOUD LAYER.
!     THE 1-FLAGS ARE RESET TO 0-FLAGS FOR THE NEXT SECTION.
!          IN (4) THE INTEGRATED MOIST AND DRY ENTHALPY ACCESSIONS
!     FOR EACH CLOUD LAYER ARE STORED INTO ALL THE CORRESPONDING
!     LAYERS IF THE FIRST IS POSITIVE WHILE THE SECOND IS NEGATIVE,
!     OTHERWISE, THE 2-FLAGS ARE ALSO RESET TO 0-FLAGS.
!          IN (5) THE ACTUAL MODIFICATIONS OF TEMPERATURE AND SPECIFIC
!     HUMIDITY ARE COMPUTED. A CLOUD-COVER VALUE IS ESTIMATED BY
!     COMPARING THE TIME AT WHICH THE ENVIRONMENT WOULD REACH
!     EQUILIBRIUM WITH THE CLOUD TO A PRESCRIBED CLOUD LIFE-TIME.
!
!     ------------------------------------------------------------------
!
!*         2.     PRELIMINARY COMPUTATIONS.
!                 ----------- -------------
!
!
!          2.0     FRACTION OF TOTAL CONDENSATE THAT IS SOLID (ZFICE)
!
!                              Calculation of ice fraction with linear
!                              variation between maxfrz and minfrz.
!                              DFMX is the value of the derivative w/r to T.
!                              Saturation is w/r to liquid for T > MAXFRZ and w/r
!                              to solid for T < MINFRZ and mixed phased in between.
!
      MAXFRZ = 268.16
      MINFRZ = 258.16
!
!
      ZFICE(:,:) =  ( MAXFRZ - TP(:,:)  ) / ( MAXFRZ - MINFRZ )
      ZFICE(:,:) = MIN( MAX( ZFICE(:,:) , 0. ) , 1. )
!
      WHERE( TP(:,:) < MINFRZ .OR. TP(:,:) > MAXFRZ )
        DFMX(:,:) = 0.
      ELSEWHERE
        DFMX(:,:) = - 1. / ( MAXFRZ - MINFRZ )
      END WHERE

!     Find the level nearest to a prescribed distance (in Pa) above the surface
      stat = neark(sigma,psp,1000.,ni,nk,kmin)
!
!
!
!*         2.1     ENVIRONMENTAL PROFILES AND PARAMETERS,
!*                 DRY AND MOIST ENTHALPY ACCESSIONS (divided by cp)
!*                 AND INITIALIZATIONS.
!
      DO jl=1,NI
         ZDSG(jl,1)=0.5*(SIGMA(jl,2)-SIGMA(jl,1))
         ZDSG(jl,NK)=0.5*(1.-SIGMA(jl,NK-1))+0.5*(1.-SIGMA(jl,NK))
         DBDT(jl) = 0.
         ZLDCP0(jl) =  ( CHLC + ZFICE(jl,NK)*CHLF ) * rCPD
         ZSDT(jl,NK)=0.
         CRRI(jl)=0.
         CRRL(jl)=0.
      END DO
!
      DO jk=2,NK-1
         DO jl=1,NI
            ZDSG(jl,jk)=0.5*(SIGMA(jl,jk+1)-SIGMA(jl,jk-1))
         END DO
      END DO
!
      DO jk=1,NK
         DO jl=1,NI
            ZPP(jl,jk)=SIGMA(jl,jk)*PSP(jl)
            ZDP(jl,jk)=ZDSG(jl,jk)*PSP(jl)
            ZTE(jl,jk)=TP(jl,jk)
            ZQSE(jl,jk)=FQSMX( ZTE(jl,jk), ZPP(jl,jk), ZFICE(jl,jk) )
            ZQE(jl,jk)=amin1(ZQSE(jl,jk),QM(jl,jk))
            ZTVE(jl,jk) = FOTVT( ZTE(jl,jk), ZQE(jl,jk) )
            ZLDCP(jl,jk) = ( CHLC + ZFICE(jl,jk)*CHLF ) &
                     / ( CPD*(1.+DELTA2*ZQE(jl,jk)) )
!
            ZQAC(jl,jk)=TQDF(jl,jk)*ZDP(jl,jk)
!
            if ( ilab(jl,jk) .ge. 1 ) ZQAC(jl,jk)=-1.
            ilab(jl,jk) = 0
            CTT(jl,jk) = 0.0
            CQT(jl,jk) = 0.0
            CCF(jl,jk) = 0.0
            QCN(jl,jk) = 0.0
         END DO
      END DO
!
!*         2.2     SPECIFY TC AND QC AT THE LOWEST LAYER TO START THE
!*                 CLOUD ASCENT. CHECK FOR POSITIVE ACCESSION
!*                 BETWEEN SURFACE AND CLOUD BASE.
!*                 ZQC=0 INDICATES STABLE CONDITIONS.
!
      DO jl=1,NI
         CPR(jl) = 0.
         ZTC(jl,NK)=ZTE(jl,KMIN(jl))
         ZQC(jl,NK)=0.
          LO=ZQAC(jl,NK).GT.0.
         IF (LO) THEN
            ZQC(jl,NK)=ZQE(jl,KMIN(jl))
            ilab(jl,NK) = 1
         ENDIF
      END DO
!
!     ------------------------------------------------------------------
!
!*         3.     CLOUD ASCENT AND FLAGGING.
!                 ----- ------ --- ---------
!
!*         3.1     CALCULATE TC AND QC AT UPPER LEVELS BY DRY ADIABATIC
!*                 LIFTING FOLLOWED BY LATENT HEAT RELEASE WHEN REQUIRED.
!*                 CONDENSATION CALCULATIONS ARE DONE WITH TWO ITERATIONS.
!***
      DO jk=NK-1,1,-1
!***
         DO jl=1,NI
            ZCP(jl)=CPD*(1.+DELTA2*ZQC(jl,jk+1))
            ZTC(jl,jk)=ZTC(jl,jk+1)+(GZM(jl,jk+1)-GZM(jl,jk))* &
               (1./ZCP(jl)+ENTRM*MAX(0.,ZTC(jl,jk+1)-ZTE(jl,jk+1)))
            ZQC(jl,jk)=ZQC(jl,jk+1)+(GZM(jl,jk+1)-GZM(jl,jk))* &
               (            ENTRM*MAX(0.,ZQC(jl,jk+1)-ZQE(jl,jk+1)))
            ZTVC = FOTVT( ZTC(jl,jk), ZQC(jl,jk) )
               LO= ZTVC.GT.ZTVE(jl,jk) .AND. ZQC(jl,jk).NE.0.
            IF (LO) ilab(jl,jk) = 1
         END DO
!
         DO jl=1,NI
             temp1 = ZTC(jl,jk)
             temp2 = ZPP(jl,jk)
            ZQSC=FQSMX( temp1, temp2, ZFICE(jl,jk) )
            temp3 = FDLESMX( temp1, ZFICE(jl,jk), DFMX(jl,jk) )
            ZCOR=ZLDCP(jl,jk)*FDQSMX( ZQSC, temp3 )
            ZQCD=AMAX1(0.,(ZQC(jl,jk)-ZQSC)/(1.+ZCOR))
            QCKTL(jl,jk) = ( 1.-ZFICE(jl,jk) ) * ZQCD
            QCKTI(jl,jk) = ZFICE(jl,jk) * ZQCD
            ZQC(jl,jk)=ZQC(jl,jk)-ZQCD
            ZTC(jl,jk)=ZTC(jl,jk)+ZQCD*ZLDCP(jl,jk)
               LO1(jl)=ZQCD.NE.0.
         END DO
!
         LO=.FALSE.
         DO jl=1,NI
            LO=LO.OR.LO1(jl)
         END DO
!
         IF (LO) THEN
            DO jl=1,NI
                temp1 = ZTC(jl,jk)
                temp2 = ZPP(jl,jk)
               ZQSC=FQSMX( temp1, temp2, ZFICE(jl,jk) )
               temp3 = FDLESMX( temp1, ZFICE(jl,jk), DFMX(jl,jk) )
               ZCOR=ZLDCP(jl,jk)*FDQSMX( ZQSC, temp3 )
               ZQCD=(ZQC(jl,jk)-ZQSC)/(1.+ZCOR)
               if (.not. LO1(jl)) ZQCD = 0.
               QCKTL(jl,jk) = ( 1.-ZFICE(jl,jk) ) * ZQCD   + QCKTL(jl,jk)
               QCKTI(jl,jk) = ZFICE(jl,jk) * ZQCD   +  QCKTI(jl,jk)
               ZQC(jl,jk)=ZQC(jl,jk)-ZQCD
               ZTC(jl,jk)=ZTC(jl,jk)+ZQCD*ZLDCP(jl,jk)
            END DO
         ENDIF
!
         DO jl=1,NI
             temp1 = ZTC(jl,jk)
             temp2 = ZQC(jl,jk)
            ZTVC = FOTVT( temp1, temp2 )
               LO= ZTVC.GT.ZTVE(jl,jk)  .AND. LO1(jl)
            IF (LO) ilab(jl,jk) = 2
               LO1(jl)=ilab(jl,jk).EQ.0
            if (LO1(jl)) ZTC(jl,jk) = ZTE(jl,jk)
            if (LO1(jl)) ZQC(jl,jk) = 0.
         END DO
!
!*         3.2     IF NOT AT THE TOP CHECK FOR NEW LIFTING LEVEL, I.E.
!*                 MOISTURE ACCESSION IN A STABLE LAYER.
!***
         IF (jk.NE.1) THEN
            DO jl=1,NI
                  LO=LO1(jl).AND.(ZQAC(jl,jk).GT.0.)
               if (LO)  ZTC(jl,jk) = ZTE(jl,min(jk,KMIN(jl)))
               if (LO) ZQC(jl,jk) = ZQE(jl,min(jk,KMIN(jl)))
            END DO
         ENDIF
!***
      END DO
!***
!*         3.3     ilab=0 UNLESS ilab=2
!*                 IKA INDICATES THE HIGHEST TOP OF A CLOUD
!*                 (TO AVOID UNNECESSARY COMPUTATIONS LATER).
!
      IKA=NK+1
!
      DO jk=1,NK
!
         DO jl=1,NI
               LO=(ilab(jl,jk).EQ.1)
            IF (LO) ilab(jl,jk) = 0
         END DO
!
         IF (IKA.EQ.NK+1) THEN
            IS=0
            DO jl=1,NI
               IS=IS+ilab(jl,jk)
            END DO
            IF (IS.NE.0) IKA=jk
         ENDIF
!
      END DO

!***
      IF (IKA.EQ.NK+1) GO TO 600
!***
!     ------------------------------------------------------------------
!
!*         4.     TOTAL MOISTURE ACCESSION
!                 ----- ------ ---------
!*                 TOTAL MOISTURE ACCESSION BE > 0
!*                 IKB IS AN UPDATE OF IKA.
!
      DO jl=1,NI
         ZSQAC(jl,NK) = 0.0
         ZSDP(jl,NK) = 0.0
      END DO
!
      DO jk=NK-1,IKA,-1
         DO jl=1,NI
            LO=ilab(jl,jk).eq.2
            if (LO) then
               ZSQAC(jl,jk) = ZSQAC(jl,jk+1)+ZQAC(jl,jk)
               ZSDP(jl,jk) = ZSDP(jl,jk+1)+ZDP(jl,jk)
            else
               ZSQAC(jl,jk) = 0.
               ZSDP(jl,jk) = 0.
            endif
         END DO
      END DO
!
      IKB=NK+1
!
      DO jk=IKA,NK-1
      jkm1=max0(jk-1,1)
!
         DO jl=1,NI
            LO=(ilab(jl,jk).EQ.2).AND.(ilab(jl,jkm1).EQ.2)
            if (LO) then
               ZSQAC(jl,jk) = ZSQAC(jl,jkm1)
               ZSDP(jl,jk) = ZSDP(jl,jkm1)
            endif
               LO = ZSQAC(jl,jk).gt.0. .and. ZSDP(jl,jk).gt.0.
            IF (.not.LO) ilab(jl,jk) = 0
         END DO
!
         IF (IKB.EQ.NK+1) THEN
            IS=0
            DO jl=1,NI
               IS=IS+ilab(jl,jk)
            END DO
            IF (IS.NE.0) IKB=jk
         ENDIF
!
      END DO
!***
      IF (IKB.EQ.NK+1) GO TO 600
!
!     Estimation of total cloud water content for precipitation.
!     Method : from the lowest level in the cloud we rise a particule
!              for the constant height DGZFIX and compute the condensation
!              like in 3.1. This becomes the total cloud water which is
!              then supposed constant through the cloud. Therefore there is
!              not dependency on the modele level distribution. The value is
!              then adjusted if needed in order to avoid negative
!              precipitation. An optinal linear profile is also available.

      DO jk=NK-1,IKB,-1
         jkm1=max0(jk-1,1)
         DO jl=1,NI
            LO=ilab(jl,jk+1).eq.0.and.ilab(jl,jk).eq.2
            if(LO)then
!              First cloud level.
!              Compute tentative total cloud water.
               ZN0=0.
               ZQ0=FQSMX(ZTC(jl,jk),ZPP(jl,jk), ZFICE(jl,jk) )
               ZT0=ZTC(jl,jk)-DGZFIX/CPD
!               print*,ZT0+ZLDCP0(jl)*ZQ0,ZQ0+ZN0

               temp2=exp(alog(ZPP(jl,jk))-DGZFIX/ &
                    (RGASD*.5*(ZTC(jl,jk)+ZT0)))
               ZQSC=FQSMX( ZT0, temp2, ZFICE(jl,jk) )
               temp3 = FDLESMX( ZT0, ZFICE(jl,jk), DFMX(jl,jk) )
               ZCOR=ZLDCP(jl,jk)*FDQSMX( ZQSC, temp3 )
               ZQCD=(ZQ0-ZQSC)/(1.+ZCOR)
               ZT0=ZT0+ZQCD*ZLDCP0(jl)
               ZQ0=ZQ0-ZQCD
               ZN0=ZN0+ZQCD
!
!              Second iteration (out of loop for vectorization)
               temp2=exp(alog(ZPP(jl,jk))-DGZFIX/ &
                    (RGASD*.5*(ZTC(jl,jk)+ZT0)))
               ZQSC=FQSMX( ZT0, temp2, ZFICE(jl,jk) )
               temp3 = FDLESMX( ZT0, ZFICE(jl,jk), DFMX(jl,jk) )
               ZCOR=ZLDCP(jl,jk)*FDQSMX( ZQSC, temp3 )
               ZQCD=(ZQ0-ZQSC)/(1.+ZCOR)
               ZT0=ZT0+ZQCD*ZLDCP0(jl)
               ZQ0=ZQ0-ZQCD
               ZN0=ZN0+ZQCD
!               print*,ZT0+ZLDCP0(jl)*ZQ0,ZQ0+ZN0
!               print*,'T,P,QS,QC=',ZT0-tcdk,temp2*.01,ZQ0,ZN0
               QCN(jl,jk)=ZN0
            else
               QCN(jl,jk)=0.
            endif
            LO=ilab(jl,jk+1).eq.2.and.ilab(jl,jk).eq.2
            if(LO)then
!              In the cloud but not the first level.
               QCN(jl,jk)=QCN(jl,jk+1)
            endif
         ENDDO
      ENDDO

!     Check if tentative total cloud water is too large.
      DO jk=NK-1,IKB,-1
         jkm1=max0(jk-1,1)
         DO jl=1,NI
            LO=ilab(jl,jk).EQ.2
            IF(LO)THEN
               temp1 = ZTC(jl,jk)
               temp2 = ZQC(jl,jk)
               ZTVC = FOTVT( temp1, temp2 )
               ZSDT(jl,jk)=ZSDT(jl,jk+1)+(ZTVC-ZTVE(jl,jk))*ZDP(jl,jk)
            ELSE
               ZSDT(jl,jk)=0.
            ENDIF
            LO=ilab(jl,jk).eq.2.and.ilab(jl,jkm1).eq.0
            IF(LO)THEN
!              Cloud top, time to check
               QCN(jl,jk)=MIN(QCN(jl,jk),ZSDT(jl,jk)/ &
                    (ZSDP(jl,jk)*ZLDCP0(jl)))
            ENDIF
         ENDDO
      ENDDO
!     Update cloud water
!      if(.false.)then
!        Constant profile
!         DO jk=IKB,NK-1
!            DO jl=1,NI
!               LO=ilab(jl,jk).eq.2.and.ilab(jl,jk-1).eq.2
!               IF(LO)THEN
!                 In cloud but not cloud top.
!                  QCN(jl,jk)=QCN(jl,jk-1)
!               ENDIF
!           ENDDO
!         ENDDO
!      else
!        Linear profile
         DO jk=IKB,NK-1
            jkm1=max0(jk-1,1)
            DO jl=1,NI
               LO=ilab(jl,jk).eq.2.and.ilab(jl,jkm1).eq.0
               IF(LO)THEN
!                 Top of cloud
                  ZPT(jl)=.5*(ZPP(jl,jkm1)+ZPP(jl,jk))
                  ZQCNT(jl)=2.*QCN(jl,jk)
               ENDIF
               LO=ilab(jl,jk).eq.2
               IF(LO)THEN
!                 In cloud
                  ZPC=.25*(ZPP(jl,jkm1)+2.*ZPP(jl,jk)+ZPP(jl,jk+1))
                  QCN(jl,jk)=ZQCNT(jl)/ZSDP(jl,jk)* &
                       (ZSDP(jl,jk)+ZPT(jl)-ZPC)
               ENDIF
            ENDDO
         ENDDO
!      endif
!***

!***
!     ------------------------------------------------------------------
!
!*         5.     HEATING AND MOISTENING
!                 ----------------------
!
!*         5.1     COMPUTE THE TOTAL CLOUD-ENVIRONMENT ENTHALPY
!*                 DIFFERENCE IN CLOUD LAYERS.
!
      DO jl=1,NI
         ZSDH(jl,NK)=0.
         ZSDQ(jl,NK)=0.
      END DO
!
      DO jk=NK-1,IKB,-1
         DO jl=1,NI
             temp1 = ZTC(jl,jk)
             temp2 = ZQC(jl,jk)
            ZTVC = FOTVT( temp1, temp2 )
            ZDQ(jl,jk)=(ZQSE(jl,jk)+QCN(jl,jk)-ZQE(jl,jk))*ZDP(jl,jk)
            ZDT(jl,jk)=(ZTVC-ZTVE(jl,jk)-ZLDCP0(jl)*QCN(jl,jk)) &
                 *ZDP(jl,jk)
!            ZDQ(jl,jk)=(ZQSE(jl,jk)-ZQE(jl,jk))*ZDP(jl,jk)
!            ZDT(jl,jk)=(ZTVC-ZTVE(jl,jk))*ZDP(jl,jk)
            ZDH = ZDT(jl,jk)+ZLDCP0(jl)*ZDQ(jl,jk)
               LO=ilab(jl,jk).EQ.2
            if (LO) then
               ZSDH(jl,jk) = ZSDH(jl,jk+1)+ZDH
               ZSDQ(jl,jk) = ZSDQ(jl,jk+1)+ZDQ(jl,jk)
            else
               ZSDH(jl,jk) = 0.
               ZSDQ(jl,jk) = 0.
            endif
         END DO
      END DO
!
      DO jk=IKB+1,NK-1
         DO jl=1,NI
               LO=(ilab(jl,jk).EQ.2).AND.(ilab(jl,jk-1).EQ.2)
            if (LO) then
               ZSDH(jl,jk) = ZSDH(jl,jk-1)
               ZSDQ(jl,jk) = ZSDQ(jl,jk-1)
            endif
         END DO
      END DO
!
!*         5.2     COMPUTE CONVECTIVE HEATING AND MOISTENING.
!*                 ESTIMATE CONVECTIVE CLOUD FRACTION.
!
      DO jk=IKB,NK-1
         DO jl=1,NI
!
            LO=ilab(jl,jk).eq.0
            if (LO) then
               ZQAC(jl,jk) = 0.
               ZSQAC(jl,jk) = 0.
            endif

               LO=ZSDH(jl,jk).GT.0.
            if (.not. LO) ZSDH(jl,jk) = -1.
               LO=ZSDQ(jl,jk).GT.0.
            if (.not. LO) ZSDQ(jl,jk) = -1.
!
            ZK = ZLDCP0(jl)*ZSQAC(jl,jk)/ZSDH(jl,jk)
!
            CQT(jl,jk) = (ZK*ZDQ(jl,jk)-ZQAC(jl,jk))/ZDP(jl,jk)
            CTT(jl,jk) = (ZK*ZDT(jl,jk)            )/ZDP(jl,jk)
!*  ajouter du detrainement
            DETRN=0.0
            CQT(jl,jk) = CQT(jl,jk)+DETRN*CTT(jl,jk)/ZLDCP0(jl)
            CTT(jl,jk) = CTT(jl,jk)-DETRN*CTT(jl,jk)
!
            CPR(jl) = CPR(jl) + CTT(jl,jk)/ZLDCP0(jl)*ZDP(jl,jk)
!
            DBDT(jl) = AMAX1(DBDT(jl),ZK)
!
!           Precipitation computation
            CRRL(jl)=CRRL(jl)+CTT(jl,jk)/(GRAV*ZLDCP0(jl))*ZDP(jl,jk)
!
!           Evaporation of precipitation (See Appendix 5 in Scientific
!           Description of RPN Physics Library -Version 3.6- )
!            print*,CRRL(jl)
            LO=ilab(jl,jk).eq.0.and.CRRL(jl).gt.0..and.DBDT(jl).ne.0.
            if(LO)then
!              Fraction of tile that is precipitating.
               PFRAC=DBDT(jl)*TAUCU
!               print*,'Fraction',PFRAC
!              Precipitation flux in that fraction.
               PRCP=CRRL(jl)/PFRAC
!              Compute (1+beta).
               HCIMPR = FDLESMX(ZTE(jl,jk),ZFICE(jl,jk),DFMX(jl,jk))
               HCIMPR = 1.+ZLDCP(jl,jk)*HCIMPR
!              Compute N (A5.2).
               XN = STPEVP*TAU*HCIMPR
!              Avoid division by zero.
               XX = MAX(PRCP,1.E-16)
               XX = SQRT(XX)
!              Compute Delta Pmax (A5.4)
               HDPMX = 2.*HKPEVP* &
                    MAX(0.,ZQSE(jl,jk)-ZQE(jl,jk))*ZDP(jl,jk)/XN
!              Explicit estimate.
               YY = 0.5*XN*HDPMX/XX
!              Implicit correction.
               YY = YY / ( 1. + 0.5*XN*XX)
!              YY=1 => complete evaporation.
               YY = MIN(YY ,1.)
               XP = PRCP*( 1.-YY)**2
!              Avoid over saturartion.
               XDP = -MIN(PRCP-XP,HDPMX)*PFRAC
!               temp1=CRRL(jl)
!               temp2=CTT(jl,jk)
!               temp3=CQT(jl,jk)
!              Update precip. flux temp. and specific hum. trend du to evaporation.
               CRRL(jl)=CRRL(jl)+XDP
               CTT(jl,jk) = CTT(jl,jk)+XDP*GRAV*ZLDCP0(jl)/ZDP(jl,jk)
               CQT(jl,jk) = CQT(jl,jk)-XDP*GRAV/ZDP(jl,jk)
!               print*,'DP,DT,DQ=',
!     $              CRRL(jl)-temp1,
!     $              CTT(jl,jk)-temp2,CQT(jl,jk)-temp3
            endif
         END DO
      END DO
!
      DO jl=1,NI
         CRRL(jl)=max(0.,CRRL(jl))
         CRRI(jl)=ZFICE(jl,NK)*CRRL(jl)
         CRRL(jl)=max(0.,CRRL(jl)-CRRI(jl))
         CPR(jl) = max( 1.E-12, CPR(jl)*rGRAV3 )
      END DO
         call vslog (cpr,cpr,ni)
      DO jl=1,NI
         CPR(jl) = 2.5 + .125 * CPR(jl)
         CPR(jl) = min( max( DBDT(jl) * TAUCU / ( 1. + DBDT(jl)*TAUCU ) , &
                                               CPR(jl) ) , 0.8 )
      END DO
!
      DO jk=IKB,NK-1
         DO jl=1,NI
             LO=ilab(jl,jk).ne.2
             if (LO) then
                CCF(jl,jk) = 0.
             else
                CCF(jl,jk) = CPR(jl)
             endif
!!        CCF(jl,jk) = CCF(jl,jk)* min((SIGMA(jl,jk)/0.8)**2, 1.0 )
             temp1=(SIGMA(jl,jk)*1.25)*(SIGMA(jl,jk)*1.25)
             CCF(jl,jk) = CCF(jl,jk)* min(temp1, 1.0 )
         END DO
      END DO
!
!
!       consistency check between the cloud fraction and cloud water
!       convert QCKTL and QCKTI from in-cloud values to average values
!
      DO jk=1,NK
      DO jl=1,NI
        lo=CCF(jl,jk).GE.0.01
        if ( .not. lo) QCKTL(jl,jk) = 0.
        if ( .not. lo) QCKTI(jl,jk) = 0.
      END DO
      END DO
!
!
!       tendency check

          DO jk=1,NK
             DO jl=1,NI
               if (cqt(jl,jk).gt.1.E-10) KSHAL(jl)=1.
             END DO
          END DO

!***
!     ------------------------------------------------------------------
!
!*         6.     RETURN WORKSPACE.
!                 ------ ----------
  600 CONTINUE
!
!
!       conversion from in-cloud values to average grid value
      DO jk=1,NK
      DO jl=1,NI
        QCKTL(jl,jk) = QCKTL(jl,jk) * CCF(jl,jk)
        QCKTI(jl,jk) = QCKTI(jl,jk) * CCF(jl,jk)
      END DO
      END DO
!
      RETURN
      END
