!-------------------------------------- LICENCE BEGIN -------------------------
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

subroutine KTRSNT_MG2(CTT,CQT,ilab,CCF,QCKTL,QCKTI,DBDT, &
     TP,QM,GZM,TQDF,PSP, &
     SIGMA, KSHAL, NI, NK)
   use tdpack
   use phy_options
   implicit none
!!!#include <arch_specific.hf>

   integer NI,NK
   real CTT(NI,NK),CQT(NI,NK)
   integer ilab(NI,NK)
   real CCF(NI,NK),QCKTL(NI,NK),QCKTI(NI,NK),DBDT(NI),TQDF(NI,NK)
   real TP(NI,NK),QM(NI,NK),GZM(NI,NK)
   real PSP(NI),SIGMA(NI,NK)
   real KSHAL(NI)

   !@Authors Claude Girard and Gerard Pellerin 1995
   !@Revision
   ! 001      G.Pellerin (Nov 98) Added kuo65 option
   !                      and change accession closure
   ! 002      C.Girard   (Dec 98) Added detrainement
   ! 003      A-M.Leduc  (March 2002) Automatic arrays
   ! 004      S.Belair, A-M. Leduc (nov 2002) added convective counter kshal
   !                         as argument. ktrsnt--->ktrsnt2
   ! 005      G.Pellerin (Avril 2003) -  CVMG... Replacements
   ! 006      G. Pellerin (Mai 03) - Conversion IBM
   !                     - calls to vexp routine (from massvp4 library)
   !                     - calls to optimized routine MFOQST
   ! 005      S.Belair   (April 2003) Outputs for in-cloud water (QCKT)
   ! 006      S.Belair, A-M. Leduc (Mai 2003) calculation of ice fraction (ZFICE)
   !                    and of saturation with respect to water and ice in
   !                    temperature interval MAXFRZ and MINFRZ. Output
   !                    QCKTL and QCKTI.
   ! 007      A-M.Leduc  (Feb 2004) Add key KTICEFRAC. Key true by default.
   !                     if set to false, zfice=0.0 and qcktl is total cloud water.
   ! 008      B. Bilodeau (Jul 2004) - Rename ktrsnt_mg
   ! 009      B. Bilodeau (May 2005) - QCKTL and QCKTI average values
   !                                   instead of in-cloud values
   ! 010      B. Bilodeau (May 2005) - New comdeck fintern
   !@Object
   !          To calculate the convective tendencies of T and Q
   !          using a scheme with a "Kuo65-type closure".
   !          Geleyn's method is used to obtain the cloud profiles.
   !@Arguments
   !            - Outputs -
   ! CTT      convective temperature tendency
   ! CQT      convective specific humidity tendency
   ! ilab     flag array: an indication of convective activity
   ! CCF      estimated cumulus cloud fraction
   ! QCKTL    estimated cumulus cloud liquid water (average over grid)
   ! QCKTI    estimated cumulus cloud solid water (average over grid)
   ! DBDT     estimated averaged cloud fraction growth rate
   !            - Inputs -
   ! TP       temperature at (t+dt)
   ! QM       specific humidity at (t-dt)
   ! GZM      geopotential
   ! TQDF     tendance diffusive de couche limite (t+dt)
   ! PSP      surface pressure at (t+dt)
   ! SIGMA    sigma levels
   ! NI       horizontal dimension
   ! NK       vertical dimension
   !@Notes
   !          The routine is divided into 5 parts:
   !           1)allocation and position for work space
   !           2)preliminary computations
   !           3)cloud ascent and flagging
   !           4)total moisture accession calculations
   !           5)cloud heating and moistening (drying) calculations

   logical LO
   integer IS,IKA,IKB,jk,jkm1,jl,stat
   integer, dimension(NI) :: KMIN
   integer, external :: neark
   real ZTVC,rgrav3,rcpd
   real ENTRM,TAUCU,DELTA2,DZETR,CHLS
   real ZCOR,ZQSC,DETRN
   real ZQCD,ZK,ZDH,temp1,temp2, temp3
   real MAXFRZ, MINFRZ

   real, dimension(NI,NK) :: ZPP
   real, dimension(NI,NK) :: ZDSG
   real, dimension(NI,NK) :: ZDP
   real, dimension(NI,NK) :: ZSDP
   real, dimension(NI,NK) :: ZQAC
   real, dimension(NI,NK) :: ZSQAC
   real, dimension(NI,NK) :: ZLDCP
   real, dimension(NI,NK) :: ZQSE
   real, dimension(NI,NK) :: ZTC
   real, dimension(NI,NK) :: ZQC
   real, dimension(NI,NK) :: ZTE
   real, dimension(NI,NK) :: ZQE
   real, dimension(NI,NK) :: ZTVE
   real, dimension(NI,NK) :: ZDQ
   real, dimension(NI,NK) :: ZSDQ
   real, dimension(NI,NK) :: ZDT
   real, dimension(NI,NK) :: ZSDH
   real, dimension(NI,NK) :: ZFICE
   real, dimension(NI,NK) :: DFMX
   real, dimension(NI   ) :: ZCP
   real, dimension(NI   ) :: ZLDCP0
   real, dimension(NI   ) :: CPR
   logical, dimension(NI   ) :: LO1


   !****************************************************

   !*    PHYSICAL CONSTANTS.
   !     -------- ----------

   rcpd = 1./CPD
   rgrav3 = 1./(GRAV*1.E3)
   !     typiquement entraine selon labda=1/GH (H=2km)
   DZETR = 2.E+03
   ENTRM = 1./(DZETR*GRAV)
   TAUCU = 900.
   DELTA2 = CPV/CPD - 1.
   CHLS   = CHLC + CHLF

   !     ------------------------------------------------------------------

   !*         1.     ALLOCATION AND POSITION FOR WORK SPACE.
   !                 ---------- --- -------- --- ---- ------


   !***

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
   !     THE 1-FLAGS ARE RESET TO 0-FLAGS FOR THE NEXT SECTION.
   !          IN (4) THE INTEGRATED MOIST AND DRY ENTHALPY ACCESSIONS
   !     FOR EACH CLOUD LAYER ARE STORED INTO ALL THE CORRESPONDING
   !     LAYERS IF THE FIRST IS POSITIVE WHILE THE SECOND IS NEGATIVE,
   !     OTHERWISE, THE 2-FLAGS ARE ALSO RESET TO 0-FLAGS.
   !          IN (5) THE ACTUAL MODIFICATIONS OF TEMPERATURE AND SPECIFIC
   !     HUMIDITY ARE COMPUTED. A CLOUD-COVER VALUE IS ESTIMATED BY
   !     COMPARING THE TIME AT WHICH THE ENVIRONMENT WOULD REACH
   !     EQUILIBRIUM WITH THE CLOUD TO A PRESCRIBED CLOUD LIFE-TIME.

   !     ------------------------------------------------------------------

   !*         2.     PRELIMINARY COMPUTATIONS.
   !                 ----------- -------------


   !          2.0     FRACTION OF TOTAL CONDENSATE THAT IS SOLID (ZFICE)

   !                              Calculation of ice fraction with linear
   !                              variation between maxfrz and minfrz.
   !                              DFMX is the value of the derivative w/r to T.
   !                              Saturation is w/r to liquid for T > MAXFRZ and w/r
   !                              to solid for T < MINFRZ and mixed phased in between.

   MAXFRZ = 268.16
   MINFRZ = 258.16


   if(kticefrac)then
      ZFICE(:,:) =  ( MAXFRZ - TP(:,:)  ) / ( MAXFRZ - MINFRZ )
      ZFICE(:,:) = min( max( ZFICE(:,:) , 0. ) , 1. )

      where( TP(:,:) < MINFRZ .or. TP(:,:) > MAXFRZ )
         DFMX(:,:) = 0.
      elsewhere
         DFMX(:,:) = - 1. / ( MAXFRZ - MINFRZ )
      end where
   else
      ZFICE(:,:) = 0.0
   endif

   !     Find the level nearest to a prescribed distance (in Pa) above the surface
   stat = neark(sigma,psp,1000.,ni,nk,kmin)



   !*         2.1     ENVIRONMENTAL PROFILES AND PARAMETERS,
   !*                 DRY AND MOIST ENTHALPY ACCESSIONS (divided by cp)
   !*                 AND INITIALIZATIONS.

   do jl=1,NI
      ZDSG(jl,1)=0.5*(SIGMA(jl,2)-SIGMA(jl,1))
      ZDSG(jl,NK)=0.5*(1.-SIGMA(jl,NK-1))+0.5*(1.-SIGMA(jl,NK))
      DBDT(jl) = 0.
      if (kticefrac) then
         ZLDCP0(jl) =  ( CHLC + ZFICE(jl,NK)*CHLF ) * rCPD
      else
         LO = TP(jl,NK).lt.TRPL
         if (LO)then
            ZLDCP0(jl) =  CHLS * rCPD
         else
            ZLDCP0(jl) = CHLC * rCPD
         endif
      endif
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
         if(kticefrac)then
            ZQSE(jl,jk)=FQSMX( ZTE(jl,jk), ZPP(jl,jk), ZFICE(jl,jk) )
         else
            ZQSE(jl,jk)=FOQST( ZTE(jl,jk), ZPP(jl,jk) )
         endif
         ZQE(jl,jk)=amin1(ZQSE(jl,jk),QM(jl,jk))
         ZTVE(jl,jk) = FOTVT( ZTE(jl,jk), ZQE(jl,jk) )
         if(kticefrac)then
            ZLDCP(jl,jk) = ( CHLC + ZFICE(jl,jk)*CHLF ) &
                 / ( CPD*(1.+DELTA2*ZQE(jl,jk)) )
         else
            LO=ZTE(jl,jk).lt.TRPL
            if (LO) then
               ZLDCP(jl,jk) = CHLS / ( CPD*(1.+DELTA2*ZQE(jl,jk)))
            else
               ZLDCP(jl,jk) = CHLC / ( CPD*(1.+DELTA2*ZQE(jl,jk)) )
            endif
         endif


         ZQAC(jl,jk)=TQDF(jl,jk)*ZDP(jl,jk)

         if ( ilab(jl,jk) .ge. 1 ) ZQAC(jl,jk)=-1.
         ilab(jl,jk) = 0
         CTT(jl,jk) = 0.0
         CQT(jl,jk) = 0.0
         CCF(jl,jk) = 0.0
      end do
   end do

   !*         2.2     SPECIFY TC AND QC AT THE LOWEST LAYER TO START THE
   !*                 CLOUD ASCENT. CHECK FOR POSITIVE ACCESSION
   !*                 BETWEEN SURFACE AND CLOUD BASE.
   !*                 ZQC=0 INDICATES STABLE CONDITIONS.

   do jl=1,NI
      CPR(jl) = 0.
      ZTC(jl,NK)=ZTE(jl,KMIN(jl))
      ZQC(jl,NK)=0.
      LO=ZQAC(jl,NK).gt.0.
      if (LO) then
         ZQC(jl,NK)=ZQE(jl,KMIN(jl))
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

      do jl=1,NI
         temp1 = ZTC(jl,jk)
         temp2 = ZPP(jl,jk)
         if(kticefrac)then
            ZQSC=FQSMX( temp1, temp2, ZFICE(jl,jk) )
            temp3 = FDLESMX( temp1, ZFICE(jl,jk), DFMX(jl,jk) )
            ZCOR=ZLDCP(jl,jk)*FDQSMX( ZQSC, temp3 )
         else
            ZQSC=FOQST( temp1, temp2 )
            ZCOR=ZLDCP(jl,jk)*FODQS( ZQSC, temp1 )
         endif
         ZQCD=AMAX1(0.,(ZQC(jl,jk)-ZQSC)/(1.+ZCOR))
         QCKTL(jl,jk) = ( 1.-ZFICE(jl,jk) ) * ZQCD
         QCKTI(jl,jk) = ZFICE(jl,jk) * ZQCD
         ZQC(jl,jk)=ZQC(jl,jk)-ZQCD
         ZTC(jl,jk)=ZTC(jl,jk)+ZQCD*ZLDCP(jl,jk)
         LO1(jl)=ZQCD.ne.0.
      end do

      LO=.false.
      do jl=1,NI
         LO=LO.or.LO1(jl)
      end do

      if (LO) then
         do jl=1,NI
            temp1 = ZTC(jl,jk)
            temp2 = ZPP(jl,jk)
            if(kticefrac)then
               ZQSC=FQSMX( temp1, temp2, ZFICE(jl,jk) )
               temp3 = FDLESMX( temp1, ZFICE(jl,jk), DFMX(jl,jk) )
               ZCOR=ZLDCP(jl,jk)*FDQSMX( ZQSC, temp3 )
            else
               ZQSC=FOQST( temp1, temp2 )
               ZCOR=ZLDCP(jl,jk)*FODQS( ZQSC, temp1 )
            endif
            ZQCD=(ZQC(jl,jk)-ZQSC)/(1.+ZCOR)
            if (.not. LO1(jl)) ZQCD = 0.
            QCKTL(jl,jk) = ( 1.-ZFICE(jl,jk) ) * ZQCD   + QCKTL(jl,jk)
            QCKTI(jl,jk) = ZFICE(jl,jk) * ZQCD   +  QCKTI(jl,jk)
            ZQC(jl,jk)=ZQC(jl,jk)-ZQCD
            ZTC(jl,jk)=ZTC(jl,jk)+ZQCD*ZLDCP(jl,jk)
         end do
      endif

      do jl=1,NI
         temp1 = ZTC(jl,jk)
         temp2 = ZQC(jl,jk)
         ZTVC = FOTVT( temp1, temp2 )
         LO= ZTVC.gt.ZTVE(jl,jk)  .and. LO1(jl)
         if (LO) ilab(jl,jk) = 2
         LO1(jl)=ilab(jl,jk).eq.0
         if (LO1(jl)) ZTC(jl,jk) = ZTE(jl,jk)
         if (LO1(jl)) ZQC(jl,jk) = 0.
      end do

      !*         3.2     IF NOT AT THE TOP CHECK FOR NEW LIFTING LEVEL, I.E.
      !*                 MOISTURE ACCESSION IN A STABLE LAYER.
      !***
      if (jk.ne.1) then
         do jl=1,NI
            LO=LO1(jl).and.(ZQAC(jl,jk).gt.0.)
            if (LO)  ZTC(jl,jk) = ZTE(jl,min(jk,KMIN(jl)))
            if (LO) ZQC(jl,jk) = ZQE(jl,min(jk,KMIN(jl)))
         end do
      endif
      !***
   end do
   !***
   !*         3.3     ilab=0 UNLESS ilab=2
   !*                 IKA INDICATES THE HIGHEST TOP OF A CLOUD
   !*                 (TO AVOID UNNECESSARY COMPUTATIONS LATER).

   IKA=NK+1

   do jk=1,NK

      do jl=1,NI
         LO=(ilab(jl,jk).eq.1)
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
   !     ------------------------------------------------------------------

   !*         4.     TOTAL MOISTURE ACCESSION
   !                 ----- ------ ---------
   !*                 TOTAL MOISTURE ACCESSION BE > 0
   !*                 IKB IS AN UPDATE OF IKA.

   do jl=1,NI
      ZSQAC(jl,NK) = 0.0
      ZSDP(jl,NK) = 0.0
   end do

   do jk=NK-1,IKA,-1
      do jl=1,NI
         LO=ilab(jl,jk).eq.2
         if (LO) then
            ZSQAC(jl,jk) = ZSQAC(jl,jk+1)+ZQAC(jl,jk)
            ZSDP(jl,jk) = ZSDP(jl,jk+1)+ZDP(jl,jk)
         else
            ZSQAC(jl,jk) = 0.
            ZSDP(jl,jk) = 0.
         endif
      end do
   end do

   IKB=NK+1

   do jk=IKA,NK-1
      jkm1=max0(jk-1,1)

      do jl=1,NI
         LO=(ilab(jl,jk).eq.2).and.(ilab(jl,jkm1).eq.2)
         if (LO) then
            ZSQAC(jl,jk) = ZSQAC(jl,jkm1)
            ZSDP(jl,jk) = ZSDP(jl,jkm1)
         endif
         LO = ZSQAC(jl,jk).gt.0. .and. ZSDP(jl,jk).gt.0.
         if (.not.LO) ilab(jl,jk) = 0
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

   !*         5.1     COMPUTE THE TOTAL CLOUD-ENVIRONMENT ENTHALPY
   !*                 DIFFERENCE IN CLOUD LAYERS.

   do jl=1,NI
      ZSDH(jl,NK)=0.
      ZSDQ(jl,NK)=0.
   end do

   do jk=NK-1,IKB,-1
      do jl=1,NI
         temp1 = ZTC(jl,jk)
         temp2 = ZQC(jl,jk)
         ZTVC = FOTVT( temp1, temp2 )
         ZDQ(jl,jk) = (ZQSE(jl,jk)-ZQE(jl,jk))*ZDP(jl,jk)
         ZDT(jl,jk) = (ZTVC-ZTVE(jl,jk))*ZDP(jl,jk)
         ZDH = ZDT(jl,jk)+ZLDCP0(jl)*ZDQ(jl,jk)
         LO=ilab(jl,jk).eq.2
         if (LO) then
            ZSDH(jl,jk) = ZSDH(jl,jk+1)+ZDH
            ZSDQ(jl,jk) = ZSDQ(jl,jk+1)+ZDQ(jl,jk)
         else
            ZSDH(jl,jk) = 0.
            ZSDQ(jl,jk) = 0.
         endif
      end do
   end do

   do jk=IKB+1,NK-1
      do jl=1,NI
         LO=(ilab(jl,jk).eq.2).and.(ilab(jl,jk-1).eq.2)
         if (LO) then
            ZSDH(jl,jk) = ZSDH(jl,jk-1)
            ZSDQ(jl,jk) = ZSDQ(jl,jk-1)
         endif
      end do
   end do

   !*         5.2     COMPUTE CONVECTIVE HEATING AND MOISTENING.
   !*                 ESTIMATE CONVECTIVE CLOUD FRACTION.

   do jk=IKB,NK-1
      do jl=1,NI

         LO=ilab(jl,jk).eq.0
         if (LO) then
            ZQAC(jl,jk) = 0.
            ZSQAC(jl,jk) = 0.
         endif

         LO=ZSDH(jl,jk).gt.0.
         if (.not. LO) ZSDH(jl,jk) = -1.
         LO=ZSDQ(jl,jk).gt.0.
         if (.not. LO) ZSDQ(jl,jk) = -1.

         ZK = ZLDCP0(jl)*ZSQAC(jl,jk)/ZSDH(jl,jk)

         CQT(jl,jk) = (ZK*ZDQ(jl,jk)-ZQAC(jl,jk))/ZDP(jl,jk)
         CTT(jl,jk) = (ZK*ZDT(jl,jk)            )/ZDP(jl,jk)
         !*  ajouter du detrainement
         DETRN=0.0
         CQT(jl,jk) = CQT(jl,jk)+DETRN*CTT(jl,jk)/ZLDCP0(jl)
         CTT(jl,jk) = CTT(jl,jk)-DETRN*CTT(jl,jk)

         CPR(jl) = CPR(jl) + CTT(jl,jk)/ZLDCP0(jl)*ZDP(jl,jk)

         DBDT(jl) = AMAX1(DBDT(jl),ZK)

      end do
   end do

   do jl=1,NI
      CPR(jl) = max( 1.E-12, CPR(jl)*rGRAV3 )
   end do
   call vslog (cpr,cpr,ni)
   do jl=1,NI
      CPR(jl) = 2.5 + .125 * CPR(jl)
      CPR(jl) = min( max( DBDT(jl) * TAUCU / ( 1. + DBDT(jl)*TAUCU ) , &
           CPR(jl) ) , 0.8 )
   end do

   do jk=IKB,NK-1
      do jl=1,NI
         LO=ilab(jl,jk).ne.2
         if (LO) then
            CCF(jl,jk) = 0.
         else
            CCF(jl,jk) = CPR(jl)
         endif
         !!        CCF(jl,jk) = CCF(jl,jk)* min((SIGMA(jl,jk)/0.8)**2, 1.0 )
         temp1=(SIGMA(jl,jk)*1.25)*(SIGMA(jl,jk)*1.25)
         CCF(jl,jk) = CCF(jl,jk)* min(temp1, 1.0 )
      end do
   end do


   !       consistency check between the cloud fraction and cloud water

   do jk=1,NK
      do jl=1,NI
         lo=CCF(jl,jk).ge.0.01
         if ( .not. lo) QCKTL(jl,jk) = 0.
         if ( .not. lo) QCKTI(jl,jk) = 0.
      end do
   end do


   !       tendency check

   do jk=1,NK
      do jl=1,NI
         if (cqt(jl,jk).gt.1.E-10) KSHAL(jl)=1.
      end do
   end do

   !***
   !     ------------------------------------------------------------------

   !*         6.     RETURN WORKSPACE.
   !                 ------ ----------
600 continue

   !       conversion from in-cloud values to average grid value
   do jk=1,NK
      do jl=1,NI
         QCKTL(jl,jk) = QCKTL(jl,jk) * CCF(jl,jk)
         QCKTI(jl,jk) = QCKTI(jl,jk) * CCF(jl,jk)
      end do
   end do

   return
end subroutine KTRSNT_MG2
