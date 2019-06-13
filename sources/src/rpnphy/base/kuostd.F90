!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END ---------------------------

subroutine KUOSTD2(CTT,CQT,ilab,CCF,DBDT, &
     TP,TM,QP,QM,GZM,PSP, &
     SIGMA, TAU, NI, NK)
   use tdpack
   implicit none
!!!#include <arch_specific.hf>

   integer NI,NK
   real CTT(NI,NK),CQT(NI,NK)
   integer ilab(NI,NK)
   real CCF(NI,NK),DBDT(NI)
   real TP(NI,NK),TM(NI,NK),QP(NI,NK),QM(NI,NK),GZM(NI,NK)
   real PSP(NI),SIGMA(NI,NK)
   real TAU

   !@Authors Claude Girard and Gerard Pellerin 1995
   !@Revisions
   !001       G. Pellerin (Mai 03) - CVMG... Replacements
   !002       G. Pellerin (Mai 03) - Conversion IBM
   !                  - calls to vexp routine (from massvp4 library)
   !                  - calls to optimized routine mfoqst3
   !@Object
   !          To calculate the convective tendencies of T and Q
   !          according to the assumptions of Kuo (1974).
   !          Geleyn's method is used to obtain the cloud profiles.
   !@Arguments
   !            - Outputs -
   ! CTT      convective temperature tendency
   ! CQT      convective specific humidity tendency
   ! ilab     flag array: an indication of convective activity
   ! CCF      estimated cumulus cloud fraction
   ! DBDT     estimated averaged cloud fraction growth rate
   !            - Inputs -
   ! TP       temperature at (t+dt)
   ! TM       temperature at (t-dt)
   ! QP       specific humidity at (t+dt)
   ! QM       specific humidity at (t-dt)
   ! GZM      geopotential
   ! PSP      surface pressure at (t+dt)
   ! SIGMA    sigma levels
   ! TAU      effective timestep (2*dt)
   ! NI       horizontal dimension
   ! NK       vertical dimension
   !@Notes
   !          The routine is divided into 5 parts:
   !           1)allocation and position for work space
   !           2)preliminary computations
   !           3)cloud ascent and flagging
   !           4)total moisture convergence and mean beta-parameter
   !             calculations
   !           5)cloud heating and moistening (drying) calculations

   logical LO
   integer IS,IKA,IKB,jk,jl
   real ZQCD,ZK,temp1
   real ZTVC
   real ENTRM,TAUCU,CHLS,DELTA2
   real ZCOR,ZEPSDP,ZKQ
   real rCPD,rGRAV3


   logical, dimension(ni) :: lo1
   real,    dimension(ni) :: zcp   ,   zldcp0,   zqsc  ,   cpr

   real, dimension(ni,nk) :: zpp   ,   zdsg  ,   zdp   ,   zsdp  ,&
        zqac  ,   zldcp ,   ztac  ,   zstac ,&
        zqse  ,   ztc   ,   zqc   ,   zte   ,&
        zqe   ,   ztve  ,   zdq   ,   zdt   ,&
        zsqac ,   zbeta ,   zsdq  ,   zsdt

   !***********************************************************************

   !*    PHYSICAL CONSTANTS.
   !     -------- ----------

   ENTRM = 5.E-6
   TAUCU = 1800.
   DELTA2 = CPV/CPD - 1.
   CHLS   = CHLC + CHLF
   rCPD   = 1./CPD
   rgrav3 = 1./(GRAV*1.E3)

   !*    SECURITY PARAMETER.
   !     --------------------

   !     *ZEPSDP* AVOIDS DIVIDING BY ZERO IN THE ABSENCE OF CLOUD

   ZEPSDP=1.E-12

   !     ------------------------------------------------------------------

   !*         1.     ALLOCATION AND POSITION FOR WORK SPACE.
   !                 ---------- --- -------- --- ---- ------


   !     METHOD.
   !     -------

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

   !     ------------------------------------------------------------------

   !*         2.     PRELIMINARY COMPUTATIONS.
   !                 ----------- -------------

   !*         2.1     ENVIRONMENTAL PROFILES AND PARAMETERS,
   !*                 TEMPERATURE (times L/cp) AND MOISTURE ACCESSIONS
   !*                 AND INITIALIZATIONS.

   do jl=1,NI
      DBDT(jl) = 0.
      LO = TP(jl,NK).lt.TRPL
      if (LO) then
         ZLDCP0(jl) = CHLS * rCPD
      else
         ZLDCP0(jl) = CHLC * rCPD
      endif

      ZDSG(jl,1)=0.5*(SIGMA(jl,2)-SIGMA(jl,1))
      ZDSG(jl,NK)=0.5*(1.-SIGMA(jl,NK-1))+0.5*(1.-SIGMA(jl,NK))
   end do

   do jk=2,NK-1
      do jl=1,NI
         ZDSG(jl,jk)=0.5*(SIGMA(jl,jk+1)-SIGMA(jl,jk-1))
      end do
   end do

   do jk=1,NK
      do jl=1,NI
         ZPP(jl,jk)=SIGMA(jl,jk)*PSP(jl)
         ZDP(jl,jk)=ZDSG(jl,jk)*PSP(jl)
         ZTE(jl,jk)=TP(jl,jk)
      end do
   end do

   call mfoqst3(ZQSE,ZTE,ZPP,NI,NK,NI)

   do jk=1,NK
      do jl=1,NI
         !            ZQSE(jl,jk)=FOQST( ZTE(jl,jk), ZPP(jl,jk) )
         ZQE(jl,jk)=QM(jl,jk)
         ZTVE(jl,jk) = FOTVT( ZTE(jl,jk), ZQE(jl,jk) )
         LO = ZTE(jl,jk).lt.TRPL
         !            rCPDv=1./( CPD*(1.+DELTA2*ZQE(jl,jk)) )
         if (LO) then
            ZLDCP(jl,jk) = CHLS / ( CPD*(1.+DELTA2*ZQE(jl,jk)) )
         else
            ZLDCP(jl,jk) = CHLC / ( CPD*(1.+DELTA2*ZQE(jl,jk)) )
         endif

         ZTAC(jl,jk)=(TP(jl,jk)-TM(jl,jk))*ZDP(jl,jk)/TAU
         ZQAC(jl,jk)=(QP(jl,jk)-ZQE(jl,jk))*ZDP(jl,jk)/TAU

         ZQAC(jl,jk)= ZQAC(jl,jk)*ilab(jl,jk)

         ilab(jl,jk) = 0
         CTT(jl,jk) = 0.0
         CQT(jl,jk) = 0.0
         CCF(jl,jk) = 0.0
      end do
   end do

   !*         2.2     SPECIFY TC AND QC AT THE LOWEST LAYER TO START THE
   !*                 CLOUD ASCENT. CHECK FOR POSITIVE MOISTURE ACCESSION
   !*                 BETWEEN SURFACE AND CLOUD BASE.
   !*                 ZQC=0 INDICATES STABLE CONDITIONS.

   do jl=1,NI
      CPR(jl) = 0.
      ZTC(jl,NK)=ZTE(jl,NK)
      ZQC(jl,NK)=0.
      if (ZQAC(jl,NK).gt.0.) then
         ZQC(jl,NK)=ZQE(jl,NK)
         ilab(jl,NK) = 1
      endif
   end do

   !     ------------------------------------------------------------------

   !*         3.     CLOUD ASCENT AND FLAGGING.
   !                 ----- ------ --- ---------

   !*         3.1     CALCULATE TC AND QC AT UPPER LEVELS BY DRY ADIABATIC
   !*                 LIFTING FOLLOWED BY LATENT HEAT RELEASE WHEN REQUIRED.
   !*                 CONDENSATION CALCULATIONS ARE DONE WITH TWO ITERATIONS.
   !***
   do jk=NK-1,1,-1
      !***
      do jl=1,NI
         ZCP(jl)=CPD*(1.+DELTA2*ZQC(jl,jk+1))
         ZTC(jl,jk)=ZTC(jl,jk+1)+(GZM(jl,jk+1)-GZM(jl,jk))* &
              (1./ZCP(jl)+ENTRM*max(0.,ZTC(jl,jk+1)-ZTE(jl,jk+1)))
         ZQC(jl,jk)=ZQC(jl,jk+1)+(GZM(jl,jk+1)-GZM(jl,jk))* &
              (            ENTRM*max(0.,ZQC(jl,jk+1)-ZQE(jl,jk+1)))
         ZTVC = FOTVT( ZTC(jl,jk), ZQC(jl,jk) )
         LO= ZTVC.gt.ZTVE(jl,jk) .and. ZQC(jl,jk).ne.0.
         if (LO) ilab(jl,jk) = 1
      end do

      call mfoqst3(ZQSC,ZTC(1,jk),ZPP(1,jk),NI,1,NI)

      do jl=1,NI
         !            ZQSC=FOQST( ZTC(jl,jk),  ZPP(jl,jk) )
         ZCOR=ZLDCP(jl,jk)*FODQS( ZQSC(jl), ZTC(jl,jk) )
         ZQCD=AMAX1(0.,(ZQC(jl,jk)-ZQSC(jl))/(1.+ZCOR))
         ZQC(jl,jk)=ZQC(jl,jk)-ZQCD
         ZTC(jl,jk)=ZTC(jl,jk)+ZQCD*ZLDCP(jl,jk)
         LO1(jl)=ZQCD.ne.0.
      end do

      LO=.false.
      do jl=1,NI
         LO=LO.or.LO1(jl)
      end do

      if (LO) then
         call mfoqst3(ZQSC,ZTC(1,jk),ZPP(1,jk),NI,1,NI)

         do jl=1,NI
            !               ZQSC=FOQST( ZTC(jl,jk), ZPP(jl,jk) )
            ZCOR=ZLDCP(jl,jk)*FODQS( ZQSC(jl), ZTC(jl,jk) )
            ZQCD=(ZQC(jl,jk)-ZQSC(jl))/(1.+ZCOR)
            if(.not.LO1(jl)) ZQCD = 0.
            ZQC(jl,jk)=ZQC(jl,jk)-ZQCD
            ZTC(jl,jk)=ZTC(jl,jk)+ZQCD*ZLDCP(jl,jk)
         end do
      endif

      do jl=1,NI
         ZTVC = FOTVT( ZTC(jl,jk), ZQC(jl,jk) )
         LO= ZTVC.gt.ZTVE(jl,jk)  .and. LO1(jl)
         if (LO) ilab(jl,jk) = 2
         LO1(jl)=ilab(jl,jk).eq.0
         if (LO1(jl)) then
            ZTC(jl,jk) = ZTE(jl,jk)
            ZQC(jl,jk) = 0.
         endif
      end do

      !*         3.2     IF NOT AT THE TOP CHECK FOR NEW LIFTING LEVEL, I.E.
      !*                 MOISTURE CONVERGENCE IN A STABLE LAYER.
      !***
      if (jk.ne.1) then
         do jl=1,NI
            LO=LO1(jl).and.(ZQAC(jl,jk).gt.0.)
            if (LO) then
               ZTC(jl,jk) = ZTE(jl,jk)
               ZQC(jl,jk) = ZQE(jl,jk)
            endif
         end do
      endif
      !***
   end do
   !***
   !*         3.3     ilab=0 FOR DRY UNSTABLE LAYERS IF NO CLOUD IS ABOVE
   !*                 ilab=3 FOR LIFTING LEVEL IF LAYER ABOVE IS UNSTABLE.
   !*                 IKA INDICATES THE HIGHEST TOP OF A CLOUD
   !*                 (TO AVOID UNNECESSARY COMPUTATIONS LATER).

   IKA=NK+1
   do jl=1,NI
      LO=ilab(jl,1).eq.1
      if (LO) ilab(jl,1) = 0
   end do

   IS=0
   do jl=1,NI
      IS=IS+ilab(jl,1)
   end do
   if (IS.ne.0) IKA=1

   do jk=2,NK

      do jl=1,NI
         LO=(ilab(jl,jk).eq.1).and.(ilab(jl,jk-1).eq.0)
         if (LO) ilab(jl,jk) = 0
      end do

      if (IKA.eq.NK+1) then
         IS=0
         do jl=1,NI
            IS=IS+ilab(jl,jk)
         end do
         if (IS.ne.0) IKA=jk
      endif

   end do
   !***
   if (IKA.eq.NK+1) GO TO 600
   !***
   do jk=NK,IKA+1,-1
      do jl=1,NI
         LO=(ilab(jl,jk).eq.0).and.(ilab(jl,jk-1).ne.0)
         if (LO) ilab(jl,jk) = 3
      end do
   end do

   !     ------------------------------------------------------------------

   !*         4.     TOTAL MOISTURE CONVERGENCE AND MEAN BETA-PARAMETER.
   !                 ----- -------- ----------- --- ---- ---------------

   !*         4.1     CALCULATE TOTAL MOISTURE ACCESSION FOR UNSTABLE
   !*                 LAYERS AND THE PARTITION PARAMETER BETA
   !*                 AVERAGED OVER CLOUD LAYERS.
   !*                 IKB IS AN UPDATE OF IKA.

   do jl=1,NI
      LO=ilab(jl,NK).ne.0
      if (LO) then
         ZSQAC(jl,NK) = ZQAC(jl,NK)
      else
         ZSQAC(jl,NK) = 0.
      endif
      ZSTAC(jl,NK) = 0.0
      ZQSE(jl,NK) = AMAX1(ZQSE(jl,NK),ZQE(jl,NK))
      ZBETA(jl,NK) = 0.0
      ZSDP(jl,NK) = 0.0
   end do

   do jk=NK-1,IKA,-1
      do jl=1,NI
         LO=ilab(jl,jk).ne.0
         if (LO) then
            ZSQAC(jl,jk) = ZQAC(jl,jk)
         else
            ZSQAC(jl,jk) = 0.
         endif
         LO=LO.and.(ilab(jl,jk).ne.3)
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
      end do
   end do

   do jl=1,NI
      !gp       ZBETA(jl,IKA) = MIN( 1.0, 4.0*(1-ZBETA(jl,IKA))**3 )
      temp1=(1-ZBETA(JL,IKA))*(1-ZBETA(JL,IKA))
      temp1=temp1*(1-ZBETA(JL,IKA))
      ZBETA(JL,IKA) = min( 1.0, 4.0*temp1 )
      LO = ( ZSQAC(jl,IKA).le.0. .and. ilab(jl,IKA).eq.2 ) &
           .or. ( ZSTAC(jl,IKA).gt.0. .and. ilab(jl,IKA).eq.2 ) &
           .or. ( ilab(jl,IKA).ne.2 )
      if (LO) ilab(jl,IKA) = 0
   end do

   IKB=IKA
   IS=0
   do jl=1,NI
      if (LO) then
         ZSTAC(jl,jk) = ZSTAC(jl,jk+1)+ZTAC(jl,jk)
         ZSDP(jl,jk) = ZSDP(jl,jk+1)+ZDP(jl,jk)
      else
         ZSTAC(jl,jk) = 0.
         ZSDP(jl,jk) = 0.
      endif
      IS=IS+ilab(jl,IKA)
   end do
   if (IS.eq.0) IKB=NK+1

   do jk=IKA+1,NK

      do jl=1,NI
         ZBETA(jl,jk) = min( 1.0, 4.0*(1-ZBETA(jl,jk))**3 )
         LO=(ilab(jl,jk).eq.2).and.(ilab(jl,jk-1).eq.2)
         if (LO) then
            ZSQAC(jl,jk) = ZSQAC(jl,jk-1)
            ZSTAC(jl,jk) = ZSTAC(jl,jk-1)
            ZBETA(jl,jk) = ZBETA(jl,jk-1)
         endif
         LO = ( ZSQAC(jl,jk).le.0. .and. ilab(jl,jk).eq.2 ) &
              .or. ( ZSTAC(jl,jk).gt.0. .and. ilab(jl,jk).eq.2 ) &
              .or. ( ilab(jl,jk).ne.2 .and. ilab(jl,jk-1).eq.0 )
         if (LO) ilab(jl,jk) = 0
      end do

      if (IKB.eq.NK+1) then
         IS=0
         do jl=1,NI
            IS=IS+ilab(jl,jk)
         end do
         if (IS.ne.0) IKB=jk
      endif

   end do
   !***
   if (IKB.eq.NK+1) GO TO 600
   !***
   !     ------------------------------------------------------------------

   !*         5.     HEATING AND MOISTENING
   !                 ----------------------

   !*         5.1     COMPUTE THE TOTAL CLOUD-ENVIRONMENT
   !*                 MOISTURE AND TEMPERATURE DIFFERENCES.

   do jl=1,NI
      ZTC(jl,NK)=ZTE(jl,NK)
      ZQC(jl,NK)=ZQE(jl,NK)
      ZDQ(jl,NK)=0.
      ZDT(jl,NK)=0.
      ZSDQ(jl,NK)=0.
      ZSDT(jl,NK)=0.
   end do

   do jk=NK-1,IKB,-1
      do jl=1,NI
         LO=ilab(jl,jk).eq.2
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
      end do
   end do

   do jk=IKB+1,NK
      do jl=1,NI
         LO=(ilab(jl,jk).eq.2).and.(ilab(jl,jk-1).eq.2)
         if (LO) then
            ZSDQ(jl,jk) = ZSDQ(jl,jk-1)
            ZSDT(jl,jk) = ZSDT(jl,jk-1)
         endif
      end do
   end do

   !*         5.2     COMPUTE CONVECTIVE HEATING AND MOISTENING.
   !*                 ESTIMATE CONVECTIVE CLOUD FRACTION.

   do jk=IKB,NK
      do jl=1,NI

         LO=ilab(jl,jk).eq.0
         if (LO) ZQAC(jl,jk) = 0.
         LO=ilab(jl,jk).ne.2
         if (LO) ZSQAC(jl,jk) = 0.
         LO=ZSDQ(jl,jk).gt.0.
         if (.not. LO) ZSDQ(jl,jk) = -1.
         LO=ZSDT(jl,jk).gt.0.
         if (.not. LO) ZSDT(jl,jk) = -1.

         ZKQ =                ZBETA(jl,jk)*ZSQAC(jl,jk)/ZSDQ(jl,jk)
         ZK = ZLDCP0(jl)*(1.-ZBETA(jl,jk))*ZSQAC(jl,jk)/ZSDT(jl,jk)

         CQT(jl,jk) = (ZKQ*ZDQ(jl,jk)-ZQAC(jl,jk))/ZDP(jl,jk)
         CTT(jl,jk) = ZK*ZDT(jl,jk)/ZDP(jl,jk)

         CPR(jl) = CPR(jl) + CTT(jl,jk)/ZLDCP0(jl)*ZDP(jl,jk)

         DBDT(jl) = AMAX1(DBDT(jl),ZKQ)

      end do
   end do

   !      DO jl=1,NI
   !         CPR(jl) = CPR(jl) / ( GRAV * 1.E3 )
   !         CPR(jl) = 2.5 + .125 * alog( max( 1.E-12, CPR(jl) ) )
   !         CPR(jl) = min( max( DBDT(jl) * TAU , CPR(jl) ) , 0.8 )
   !      END DO

   do jl=1,NI
      cpr(jl) = max( 1.E-12, CPR(jl)*rGRAV3 )
   end do
   call vslog (cpr,cpr,ni)
   do jl=1,NI
      CPR(jl) = 2.5 + .125 * CPR(jl)
      CPR(jl) = min( max( DBDT(jl) * TAU , CPR(jl) ) , 0.8 )
   end do

   do jk=IKB,NK-1
      do jl=1,NI
         LO=ilab(jl,jk).ne.2
         if (LO) then
            CCF(jl,jk) = 0.
         else
            CCF(jl,jk) = CPR(jl)
         endif
         temp1=(SIGMA(jl,jk)*1.25)*(SIGMA(jl,jk)*1.25)
         CCF(jl,jk) = CCF(jl,jk)* min(temp1, 1.0 )
      end do
   end do
   !***
   !     ------------------------------------------------------------------

   !*         6.     RETURN WORKSPACE.
   !                 ------ ----------
600 continue

   return
end subroutine KUOSTD2
