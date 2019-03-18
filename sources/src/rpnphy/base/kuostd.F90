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
!**S/P KUOSTD
!
      SUBROUTINE KUOSTD ( CTT,CQT,ilab,CCF,DBDT, &
                          TP,TM,QP,QM,GZM,PSP,PSM, &
                          SIGMA, TAU, NI, NK )
      use tdpack
      implicit none
#include <arch_specific.hf>
!
!
      INTEGER NI,NK
      REAL CTT(NI,NK),CQT(NI,NK)
      INTEGER ilab(NI,NK)
      REAL CCF(NI,NK),DBDT(NI)
      REAL TP(NI,NK),TM(NI,NK),QP(NI,NK),QM(NI,NK),GZM(NI,NK)
      REAL PSP(NI),PSM(NI),SIGMA(NI,NK)
      REAL TAU
!
!Authors
!          Claude Girard and Gerard Pellerin 1995
!Revisions
!001       G. Pellerin (Mai 03) - CVMG... Replacements
!002       G. Pellerin (Mai 03) - Conversion IBM
!                  - calls to vexp routine (from massvp4 library)
!                  - calls to optimized routine mfoqst3
!
!Object
!          To calculate the convective tendencies of T and Q
!          according to the assumptions of Kuo (1974).
!          Geleyn's method is used to obtain the cloud profiles.
!
!Arguments
!
!            - Outputs -
! CTT      convective temperature tendency
! CQT      convective specific humidity tendency
! ilab     flag array: an indication of convective activity
! CCF      estimated cumulus cloud fraction
! DBDT     estimated averaged cloud fraction growth rate
!
!            - Inputs -
! TP       temperature at (t+dt)
! TM       temperature at (t-dt)
! QP       specific humidity at (t+dt)
! QM       specific humidity at (t-dt)
! GZM      geopotential
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
!           4)total moisture convergence and mean beta-parameter
!             calculations
!           5)cloud heating and moistening (drying) calculations
!
!*
      LOGICAL LO
      INTEGER IS,IKA,IKB,jk,jl
      REAL ZQCD,ZK,temp1,temp2
      REAL ZTVC
      REAL ENTRM,TAUCU,CHLS,DELTA2
      REAL ZCOR,ZEPSDP,ZKQ
      real rCPD,rCPDv,rGRAV3
!
!***********************************************************************
!     AUTOMATIC ARRAYS
!***********************************************************************
!
      logical, dimension(ni) :: lo1
      real,    dimension(ni) :: zcp   ,   zldcp0,   zqsc  ,   cpr
!
      real, dimension(ni,nk) :: zpp   ,   zdsg  ,   zdp   ,   zsdp  ,&
                                zqac  ,   zldcp ,   ztac  ,   zstac ,&
                                zqse  ,   ztc   ,   zqc   ,   zte   ,&
                                zqe   ,   ztve  ,   zdq   ,   zdt   ,& 
                                zsqac ,   zbeta ,   zsdq  ,   zsdt
!
!***********************************************************************
!
!*    PHYSICAL CONSTANTS.
!     -------- ----------

      ENTRM = 5.E-6
      TAUCU = 1800.
      DELTA2 = CPV/CPD - 1.
      CHLS   = CHLC + CHLF
      rCPD   = 1./CPD
      rgrav3 = 1./(GRAV*1.E3)
!
!*    SECURITY PARAMETER.
!     --------------------
!
!     *ZEPSDP* AVOIDS DIVIDING BY ZERO IN THE ABSENCE OF CLOUD
!
      ZEPSDP=1.E-12
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
!     3 = LIFTING LEVEL IF DIFFERENT FROM THE GROUND.
!          IN (4) THE TOTAL MOISTURE CONVERGENCE FOR EACH SLAB OF
!     NON-0-FLAGS IS STORED INTO ALL THE CORRESPONDING LAYERS IF IT IS
!     POSITIVE AND THE SAME IS DONE FOR THE MEAN OF 1-RELATIVE HUMIDITY
!     FOR EACH CORRESPONDING SLAB OF 2-FLAGS.
!          IN (5) THE ACTUAL MODIFICATIONS OF TEMPERATURE AND SPECIFIC
!     HUMIDITY ARE COMPUTED. FIRST THE ENVIRONMENTAL MOISTENING IS
!     TAKEN PROPORTIONAL TO THE MOISTURE CONVERGENCE WEIGHTED BY "BETA"
!     AND TO THE SATURATION DEFICIT OF SPECIFIC HUMIDITY.
!     SECOND THE FORMATION OF PRECIPITATION IS TAKEN PROPORTIONAL TO THE
!     MOISTURE CONVERGENCE WEIGHTED BY "1 - BETA"  AND TO THE
!     TEMPERATURE DIFFERENCE BETWEEN THE CLOUD AND THE ENVIRONMENT.
!     "BETA" , THE MOISTENING PARAMETER, IS A FUNCTION OF MEAN REL. HUM.
!     A CLOUD-COVER VALUE IS OBTAINED BY COMPARING THE TIME AT WHICH THE
!     ENVIRONMENT WOULD REACH EQUILIBRIUM WITH THE CLOUD TO A PRESCRIBED
!     LIFE-TIME VALUE FOR THE CLOUD ITSELF.
!
!     ------------------------------------------------------------------
!
!*         2.     PRELIMINARY COMPUTATIONS.
!                 ----------- -------------
!
!*         2.1     ENVIRONMENTAL PROFILES AND PARAMETERS,
!*                 TEMPERATURE (times L/cp) AND MOISTURE ACCESSIONS
!*                 AND INITIALIZATIONS.
!
      DO jl=1,NI
         DBDT(jl) = 0.
            LO = TP(jl,NK).LT.TRPL
         if (LO) then
            ZLDCP0(jl) = CHLS * rCPD
         else
            ZLDCP0(jl) = CHLC * rCPD
         endif

         ZDSG(jl,1)=0.5*(SIGMA(jl,2)-SIGMA(jl,1))
         ZDSG(jl,NK)=0.5*(1.-SIGMA(jl,NK-1))+0.5*(1.-SIGMA(jl,NK))
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
         END DO
      END DO
!
       CALL mfoqst3(ZQSE,ZTE,ZPP,NI,NK,NI)

      DO jk=1,NK
         DO jl=1,NI
!            ZQSE(jl,jk)=FOQST( ZTE(jl,jk), ZPP(jl,jk) )
            ZQE(jl,jk)=QM(jl,jk)
            ZTVE(jl,jk) = FOTVT( ZTE(jl,jk), ZQE(jl,jk) )
               LO = ZTE(jl,jk).LT.TRPL
!            rCPDv=1./( CPD*(1.+DELTA2*ZQE(jl,jk)) )
            if (LO) then
               ZLDCP(jl,jk) = CHLS / ( CPD*(1.+DELTA2*ZQE(jl,jk)) )
            else
               ZLDCP(jl,jk) = CHLC / ( CPD*(1.+DELTA2*ZQE(jl,jk)) )
            endif
!
            ZTAC(jl,jk)=(TP(jl,jk)-TM(jl,jk))*ZDP(jl,jk)/TAU
            ZQAC(jl,jk)=(QP(jl,jk)-ZQE(jl,jk))*ZDP(jl,jk)/TAU
!
            ZQAC(jl,jk)= ZQAC(jl,jk)*ilab(jl,jk)
!
            ilab(jl,jk) = 0
            CTT(jl,jk) = 0.0
            CQT(jl,jk) = 0.0
            CCF(jl,jk) = 0.0
         END DO
      END DO
!
!*         2.2     SPECIFY TC AND QC AT THE LOWEST LAYER TO START THE
!*                 CLOUD ASCENT. CHECK FOR POSITIVE MOISTURE ACCESSION
!*                 BETWEEN SURFACE AND CLOUD BASE.
!*                 ZQC=0 INDICATES STABLE CONDITIONS.
!
      DO jl=1,NI
         CPR(jl) = 0.
         ZTC(jl,NK)=ZTE(jl,NK)
         ZQC(jl,NK)=0.
         IF (ZQAC(jl,NK).GT.0.) THEN
            ZQC(jl,NK)=ZQE(jl,NK)
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
        CALL mfoqst3(ZQSC,ZTC(1,jk),ZPP(1,jk),NI,1,NI)

         DO jl=1,NI
!            ZQSC=FOQST( ZTC(jl,jk),  ZPP(jl,jk) )
            ZCOR=ZLDCP(jl,jk)*FODQS( ZQSC(jl), ZTC(jl,jk) )
            ZQCD=AMAX1(0.,(ZQC(jl,jk)-ZQSC(jl))/(1.+ZCOR))
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
          CALL mfoqst3(ZQSC,ZTC(1,jk),ZPP(1,jk),NI,1,NI)

            DO jl=1,NI
!               ZQSC=FOQST( ZTC(jl,jk), ZPP(jl,jk) )
               ZCOR=ZLDCP(jl,jk)*FODQS( ZQSC(jl), ZTC(jl,jk) )
               ZQCD=(ZQC(jl,jk)-ZQSC(jl))/(1.+ZCOR)
               IF(.not.LO1(jl)) ZQCD = 0.
               ZQC(jl,jk)=ZQC(jl,jk)-ZQCD
               ZTC(jl,jk)=ZTC(jl,jk)+ZQCD*ZLDCP(jl,jk)
            END DO
         ENDIF
!
         DO jl=1,NI
            ZTVC = FOTVT( ZTC(jl,jk), ZQC(jl,jk) )
               LO= ZTVC.GT.ZTVE(jl,jk)  .AND. LO1(jl)
            IF (LO) ilab(jl,jk) = 2
               LO1(jl)=ilab(jl,jk).EQ.0
            if (LO1(jl)) THEN
               ZTC(jl,jk) = ZTE(jl,jk)
               ZQC(jl,jk) = 0.
            endif
         END DO
!
!*         3.2     IF NOT AT THE TOP CHECK FOR NEW LIFTING LEVEL, I.E.
!*                 MOISTURE CONVERGENCE IN A STABLE LAYER.
!***
         IF (jk.NE.1) THEN
            DO jl=1,NI
                  LO=LO1(jl).AND.(ZQAC(jl,jk).GT.0.)
               if (LO) then
                  ZTC(jl,jk) = ZTE(jl,jk)
                  ZQC(jl,jk) = ZQE(jl,jk)
               endif
            END DO
         ENDIF
!***
      END DO
!***
!*         3.3     ilab=0 FOR DRY UNSTABLE LAYERS IF NO CLOUD IS ABOVE
!*                 ilab=3 FOR LIFTING LEVEL IF LAYER ABOVE IS UNSTABLE.
!*                 IKA INDICATES THE HIGHEST TOP OF A CLOUD
!*                 (TO AVOID UNNECESSARY COMPUTATIONS LATER).
!
      IKA=NK+1
      DO jl=1,NI
            LO=ilab(jl,1).EQ.1
         IF (LO) ilab(jl,1) = 0
      END DO
!
      IS=0
      DO jl=1,NI
         IS=IS+ilab(jl,1)
      END DO
      IF (IS.NE.0) IKA=1
!
      DO jk=2,NK
!
         DO jl=1,NI
               LO=(ilab(jl,jk).EQ.1).AND.(ilab(jl,jk-1).EQ.0)
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
      DO jk=NK,IKA+1,-1
         DO jl=1,NI
               LO=(ilab(jl,jk).EQ.0).AND.(ilab(jl,jk-1).NE.0)
            IF (LO) ilab(jl,jk) = 3
         END DO
      END DO
!
!     ------------------------------------------------------------------
!
!*         4.     TOTAL MOISTURE CONVERGENCE AND MEAN BETA-PARAMETER.
!                 ----- -------- ----------- --- ---- ---------------
!
!*         4.1     CALCULATE TOTAL MOISTURE ACCESSION FOR UNSTABLE
!*                 LAYERS AND THE PARTITION PARAMETER BETA
!*                 AVERAGED OVER CLOUD LAYERS.
!*                 IKB IS AN UPDATE OF IKA.
!
      DO jl=1,NI
            LO=ilab(jl,NK).NE.0
          if (LO) then
            ZSQAC(jl,NK) = ZQAC(jl,NK)
         else
            ZSQAC(jl,NK) = 0.
         endif
         ZSTAC(jl,NK) = 0.0
         ZQSE(jl,NK) = AMAX1(ZQSE(jl,NK),ZQE(jl,NK))
         ZBETA(jl,NK) = 0.0
         ZSDP(jl,NK) = 0.0
      END DO
!
      DO jk=NK-1,IKA,-1
         DO jl=1,NI
             LO=ilab(jl,jk).NE.0
            if (LO) then
               ZSQAC(jl,jk) = ZQAC(jl,jk)
            else
               ZSQAC(jl,jk) = 0.
            endif
               LO=LO.AND.(ilab(jl,jk).NE.3)
            if (LO) then
               ZSQAC(jl,jk) = ZSQAC(jl,jk+1)+ZSQAC(jl,jk)
            endif
               LO=ilab(jl,jk).eq.2
            if (LO) then
               ZSTAC(jl,jk) = ZSTAC(jl,jk+1)+ZTAC(jl,jk)
                ZSDP(jl,jk) = ZSDP(jl,jk+1)+ZDP(jl,jk)
            else
               ZSTAC(jl,jk) = 0.
                ZSDP(jl,jk) = 0.
            endif
            ZQSE(jl,jk) = AMAX1(ZQSE(jl,jk),ZQE(jl,jk))
            ZBETA(jl,jk) = ZQE(jl,jk) / ZQSE(jl,jk)
            if (LO) then
             ZBETA(jl,jk) = (ZSDP(jl,jk+1)*ZBETA(jl,jk+1)+ZDP(jl,jk) &
                          *ZBETA(jl,jk))/AMAX1(ZSDP(jl,jk),ZEPSDP)
            else
             ZBETA(jl,jk) = 0.
            endif
         END DO
      END DO
!
      DO jl=1,NI
!gp       ZBETA(jl,IKA) = MIN( 1.0, 4.0*(1-ZBETA(jl,IKA))**3 )
         temp1=(1-ZBETA(JL,IKA))*(1-ZBETA(JL,IKA))
          temp1=temp1*(1-ZBETA(JL,IKA))
         ZBETA(JL,IKA) = MIN( 1.0, 4.0*temp1 )
            LO = ( ZSQAC(jl,IKA).LE.0. .AND. ilab(jl,IKA).EQ.2 ) &
            .or. ( ZSTAC(jl,IKA).GT.0. .AND. ilab(jl,IKA).EQ.2 ) &
            .or. ( ilab(jl,IKA).NE.2 )
         IF (LO) ilab(jl,IKA) = 0
      END DO
!
      IKB=IKA
      IS=0
      DO jl=1,NI
            if (LO) then
               ZSTAC(jl,jk) = ZSTAC(jl,jk+1)+ZTAC(jl,jk)
                ZSDP(jl,jk) = ZSDP(jl,jk+1)+ZDP(jl,jk)
            else
               ZSTAC(jl,jk) = 0.
                ZSDP(jl,jk) = 0.
            endif
         IS=IS+ilab(jl,IKA)
      END DO
      IF (IS.EQ.0) IKB=NK+1
!
      DO jk=IKA+1,NK
!
         DO jl=1,NI
            ZBETA(jl,jk) = MIN( 1.0, 4.0*(1-ZBETA(jl,jk))**3 )
               LO=(ilab(jl,jk).EQ.2).AND.(ilab(jl,jk-1).EQ.2)
            if (LO) then
               ZSQAC(jl,jk) = ZSQAC(jl,jk-1)
               ZSTAC(jl,jk) = ZSTAC(jl,jk-1)
               ZBETA(jl,jk) = ZBETA(jl,jk-1)
            endif
               LO = ( ZSQAC(jl,jk).LE.0. .AND. ilab(jl,jk).EQ.2 ) &
               .or. ( ZSTAC(jl,jk).GT.0. .AND. ilab(jl,jk).EQ.2 ) &
               .or. ( ilab(jl,jk).NE.2 .AND. ilab(jl,jk-1).EQ.0 )
            IF (LO) ilab(jl,jk) = 0
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
!***
!     ------------------------------------------------------------------
!
!*         5.     HEATING AND MOISTENING
!                 ----------------------
!
!*         5.1     COMPUTE THE TOTAL CLOUD-ENVIRONMENT
!*                 MOISTURE AND TEMPERATURE DIFFERENCES.
!
      DO jl=1,NI
         ZTC(jl,NK)=ZTE(jl,NK)
         ZQC(jl,NK)=ZQE(jl,NK)
         ZDQ(jl,NK)=0.
         ZDT(jl,NK)=0.
         ZSDQ(jl,NK)=0.
         ZSDT(jl,NK)=0.
      END DO
!
      DO jk=NK-1,IKB,-1
         DO jl=1,NI
               LO=ilab(jl,jk).EQ.2
           if (.not.LO) then
              ZTC(jl,jk) = ZTE(jl,jk)
              ZQC(jl,jk) = ZQE(jl,jk)
            endif
            ZTVC = FOTVT( ZTC(jl,jk), ZQC(jl,jk) )
            ZDQ(jl,jk) = (ZQSE(jl,jk)-ZQE(jl,jk))*ZDP(jl,jk)
            ZDT(jl,jk) = (ZTVC-ZTVE(jl,jk))*ZDP(jl,jk)
            if (LO) then
               ZSDQ(jl,jk) = ZSDQ(jl,jk+1)+ZDQ(jl,jk)
               ZSDT(jl,jk) = ZSDT(jl,jk+1)+ZDT(jl,jk)
            else
               ZSDQ(jl,jk) = 0.
               ZSDT(jl,jk) = 0.
            endif
         END DO
      END DO
!
      DO jk=IKB+1,NK
         DO jl=1,NI
               LO=(ilab(jl,jk).EQ.2).AND.(ilab(jl,jk-1).EQ.2)
            if (LO) then
               ZSDQ(jl,jk) = ZSDQ(jl,jk-1)
               ZSDT(jl,jk) = ZSDT(jl,jk-1)
            endif
         END DO
      END DO
!
!*         5.2     COMPUTE CONVECTIVE HEATING AND MOISTENING.
!*                 ESTIMATE CONVECTIVE CLOUD FRACTION.
!
      DO jk=IKB,NK
         DO jl=1,NI
!
             LO=ilab(jl,jk).eq.0
            if (LO) ZQAC(jl,jk) = 0.
               LO=ilab(jl,jk).ne.2
            if (LO) ZSQAC(jl,jk) = 0.
               LO=ZSDQ(jl,jk).GT.0.
            if (.not. LO) ZSDQ(jl,jk) = -1.
               LO=ZSDT(jl,jk).GT.0.
            if (.not. LO) ZSDT(jl,jk) = -1.
!
            ZKQ =                ZBETA(jl,jk)*ZSQAC(jl,jk)/ZSDQ(jl,jk)
            ZK = ZLDCP0(jl)*(1.-ZBETA(jl,jk))*ZSQAC(jl,jk)/ZSDT(jl,jk)
!
            CQT(jl,jk) = (ZKQ*ZDQ(jl,jk)-ZQAC(jl,jk))/ZDP(jl,jk)
            CTT(jl,jk) = ZK*ZDT(jl,jk)/ZDP(jl,jk)
!
            CPR(jl) = CPR(jl) + CTT(jl,jk)/ZLDCP0(jl)*ZDP(jl,jk)
!
            DBDT(jl) = AMAX1(DBDT(jl),ZKQ)
!
         END DO
      END DO
!
!      DO jl=1,NI
!         CPR(jl) = CPR(jl) / ( GRAV * 1.E3 )
!         CPR(jl) = 2.5 + .125 * alog( max( 1.E-12, CPR(jl) ) )
!         CPR(jl) = min( max( DBDT(jl) * TAU , CPR(jl) ) , 0.8 )
!      END DO
!
      DO jl=1,NI
            cpr(jl) = max( 1.E-12, CPR(jl)*rGRAV3 )
      END DO
         call vslog (cpr,cpr,ni)
      DO jl=1,NI
         CPR(jl) = 2.5 + .125 * CPR(jl)
         CPR(jl) = min( max( DBDT(jl) * TAU , CPR(jl) ) , 0.8 )
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
             temp1=(SIGMA(jl,jk)*1.25)*(SIGMA(jl,jk)*1.25)
             CCF(jl,jk) = CCF(jl,jk)* min(temp1, 1.0 )
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
      RETURN
      END
