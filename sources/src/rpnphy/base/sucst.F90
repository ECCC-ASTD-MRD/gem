SUBROUTINE SUCST(KULOUT,KDAT,KSSS,KPRINTLEV)

!**** *SUCST * - Routine to initialize the constants of the model.

!     Purpose.
!     --------
!           Initialize and print the common YOMCST + initialize
!         date and time of YOMRIP.

!**   Interface.
!     ----------
!        *CALL* *SUCST (..)

!        Explicit arguments :
!        --------------------

!        KULOUT  - logical unit for the output
!        KDAT    - date in the form AAAAMMDD
!        KSSS    - number of seconds in the day
!        KPRINTLEV - printing level

!        Implicit arguments :
!        --------------------
!        COMMON YOMCST
!        COMMON YOMRIP

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-10-15
!        Additions : 90-07-30 (J.-F. Geleyn)
!                    91-11-15 (M. Deque)
!                    96-08-12 M.Hamrud - Reduce printing
!     ------------------------------------------------------------------


USE YOMCST   , ONLY : RPI      ,RCLUM    ,RHPLA    ,RKBOL    ,&
            &RNAVO    ,RDAY     ,REA      ,REPSM    ,RSIYEA   ,&
            &RSIDAY   ,ROMEGA   ,RA       ,RG       ,R1SA     ,&
            &RSIGMA   ,RI0      ,R        ,RMD      ,RMV      ,&
            &RMO3     ,RD       ,RV       ,RCPD     ,RCPV     ,&
            &RCVD     ,RCVV     ,RKAPPA   ,RETV     ,RCW      ,&
            &RCS      ,RLVTT    ,RLSTT    ,RLVZER   ,RLSZER   ,&
            &RLMLT    ,RTT      ,RATM     ,RDT      ,RESTT    ,&
            &RALPW    ,RBETW    ,RGAMW    ,RALPS    ,RBETS    ,&
            &RGAMS    ,RALPD    ,RBETD    ,RGAMD
USE YOMRIP   , ONLY : RTIMST   ,RTIMTR

IMPLICIT NONE
!!!#include <arch_specific.hf>

#define _ZERO_   0.0
#define _HALF_   0.5
#define _ONE_    1.0
#define _TWO_    2.0


!     DUMMY INTEGER SCALARS
integer :: KDAT
integer :: KPRINTLEV
integer :: KSSS
integer :: KULOUT


!     LOCAL INTEGER SCALARS
integer :: IA, ID, IDAT, IM, ISSS, J

!     LOCAL REAL SCALARS
real :: ZDE, ZET, ZJU, ZRS, ZRSREL, ZTETA, ZTI


#include "fctast.cdk"
#include "fcttrm.cdk"
#include "fcttim.cdk"
!      -----------------------------------------------------------------

!*       1.    DEFINE FUNDAMENTAL CONSTANTS.
!              -----------------------------


RPI=_TWO_*ASIN(_ONE_)
RCLUM=299792458.
RHPLA=6.6260755E-34
RKBOL=1.380658E-23
RNAVO=6.0221367E+23

!     ------------------------------------------------------------------

!*       2.    DEFINE ASTRONOMICAL CONSTANTS.
!              ------------------------------

RDAY=86400.
REA=149597870000.
REPSM=0.409093

RSIYEA=365.25*RDAY*_TWO_*RPI/6.283076
RSIDAY=RDAY/(_ONE_+RDAY/RSIYEA)
ROMEGA=_TWO_*RPI/RSIDAY

IDAT=KDAT
ISSS=KSSS
ID=NDD(IDAT)
IM=NMM(IDAT)
IA=NCCAA(IDAT)
ZJU=RJUDAT(IA,IM,ID)
ZTI=RTIME(IA,IM,ID,ISSS)
RTIMST=ZTI
RTIMTR=ZTI
ZTETA=RTETA(ZTI)
ZRS=RRS(ZTETA)
ZDE=RDS(ZTETA)
ZET=RET(ZTETA)
ZRSREL=ZRS/REA

!     ------------------------------------------------------------------

!*       3.    DEFINE GEOIDE.
!              --------------

RG=9.80665
RA=6371229.
R1SA=REAL(_ONE_/REAL(RA,KIND(_ONE_)),KIND(R1SA))

!     ------------------------------------------------------------------

!*       4.    DEFINE RADIATION CONSTANTS.
!              ---------------------------

!RSIGMA=_TWO_ * RPI**5 * RKBOL**4 /(15.* RCLUM**2 * RHPLA**3)
!PV-feb2015: The unnecessary calculation of the stephan-boltzmann constant above was causing a crash
!   since the passage to the intel compiler ; specify the value explicitly
RSIGMA=0.5669800000000e-07
RI0=1370.

!     ------------------------------------------------------------------

!*       5.    DEFINE THERMODYNAMIC CONSTANTS, GAS PHASE.
!              ------------------------------------------

R=RNAVO*RKBOL
RMD=28.9644
RMV=18.0153
RMO3=47.9942
RD=1000.*R/RMD
RV=1000.*R/RMV
RCPD=3.5*RD
RCVD=RCPD-RD
RCPV=4. *RV
RCVV=RCPV-RV
RKAPPA=RD/RCPD
RETV=RV/RD-_ONE_

!     ------------------------------------------------------------------

!*       6.    DEFINE THERMODYNAMIC CONSTANTS, LIQUID PHASE.
!              ---------------------------------------------

RCW=4218.

!     ------------------------------------------------------------------

!*       7.    DEFINE THERMODYNAMIC CONSTANTS, SOLID PHASE.
!              --------------------------------------------

RCS=2106.

!     ------------------------------------------------------------------

!*       8.    DEFINE THERMODYNAMIC CONSTANTS, TRANSITION OF PHASE.
!              ----------------------------------------------------

RTT=273.16
RDT=11.82
RLVTT=2.5008E+6
RLSTT=2.8345E+6
RLVZER=RLVTT+RTT*(RCW-RCPV)
RLSZER=RLSTT+RTT*(RCS-RCPV)
RLMLT=RLSTT-RLVTT
RATM=100000.

!     ------------------------------------------------------------------

!*       9.    SATURATED VAPOUR PRESSURE.
!              --------------------------

RESTT=611.14
RGAMW=(RCW-RCPV)/RV
RBETW=RLVTT/RV+RGAMW*RTT
RALPW=LOG(RESTT)+RBETW/RTT+RGAMW*LOG(RTT)
RGAMS=(RCS-RCPV)/RV
RBETS=RLSTT/RV+RGAMS*RTT
RALPS=LOG(RESTT)+RBETS/RTT+RGAMS*LOG(RTT)
RGAMD=RGAMS-RGAMW
RBETD=RBETS-RBETW
RALPD=RALPS-RALPW

!     ------------------------------------------------------------------

!*      10.    PRINTS

IF (KPRINTLEV >= 1) THEN
  WRITE(KULOUT,'(''0*** Constants of the ICM   ***'')')
  WRITE(KULOUT,'('' *** Fundamental constants ***'')')
  WRITE(KULOUT,'(''           PI = '',E14.7,'' -'')')RPI
  WRITE(KULOUT,'(''            c = '',E14.7,''m s-1'')')RCLUM
  WRITE(KULOUT,'(''            h = '',E14.7,''J s'')')RHPLA
  WRITE(KULOUT,'(''            K = '',E14.7,''J K-1'')')RKBOL
  WRITE(KULOUT,'(''            N = '',E14.7,''mol-1'')')RNAVO
  WRITE(KULOUT,'('' *** Astronomical constants ***'')')
  WRITE(KULOUT,'(''          day = '',E14.7,'' s'')')RDAY
  WRITE(KULOUT,'('' half g. axis = '',E14.7,'' m'')')REA
  WRITE(KULOUT,'('' mean anomaly = '',E14.7,'' -'')')REPSM
  WRITE(KULOUT,'('' sideral year = '',E14.7,'' s'')')RSIYEA
  WRITE(KULOUT,'(''  sideral day = '',E14.7,'' s'')')RSIDAY
  WRITE(KULOUT,'(''        omega = '',E14.7,'' s-1'')')ROMEGA

  WRITE(KULOUT,'('' The initial date of the run is :'')')
  WRITE(KULOUT,'(1X,I8,1X,I5,5X,I4,1X,I2,1X,I2)')IDAT,ISSS,IA,IM,ID
  WRITE(KULOUT,'('' The Julian date is : '',F11.2)') ZJU
  WRITE(KULOUT,'('' Time of the model  : '',F15.2,'' s'')')ZTI
  WRITE(KULOUT,'('' Distance Earth-Sun : '',E14.7,'' m'')')ZRS
  WRITE(KULOUT,'('' Relative Dist. E-S : '',E14.7,'' m'')')ZRSREL
  WRITE(KULOUT,'('' Declination        : '',F12.5)') ZDE
  WRITE(KULOUT,'('' Eq. of time        : '',F12.5,'' s'')')ZET
  WRITE(KULOUT,'('' ***         Geoide         ***'')')
  WRITE(KULOUT,'(''      Gravity = '',E14.7,'' m s-2'')')RG
  WRITE(KULOUT,'('' Earth radius = '',E14.7,'' m'')')RA
  WRITE(KULOUT,'('' Inverse E.R. = '',E14.7,'' m'')')R1SA
  WRITE(KULOUT,'('' ***        Radiation       ***'')')
  WRITE(KULOUT,'('' Stefan-Bol.  = '',E14.7,'' W m-2 K-4'')')  RSIGMA
  WRITE(KULOUT,'('' Solar const. = '',E14.7,'' W m-2'')')RI0
  WRITE(KULOUT,'('' *** Thermodynamic, gas     ***'')')
  WRITE(KULOUT,'('' Perfect gas  = '',E14.7)') R
  WRITE(KULOUT,'('' Dry air mass = '',E14.7)') RMD
  WRITE(KULOUT,'('' Vapour  mass = '',E14.7)') RMV
  WRITE(KULOUT,'('' Ozone   mass = '',E14.7)') RMO3
  WRITE(KULOUT,'('' Dry air cst. = '',E14.7)') RD
  WRITE(KULOUT,'('' Vapour  cst. = '',E14.7)') RV
  WRITE(KULOUT,'(''         Cpd  = '',E14.7)') RCPD
  WRITE(KULOUT,'(''         Cvd  = '',E14.7)') RCVD
  WRITE(KULOUT,'(''         Cpv  = '',E14.7)') RCPV
  WRITE(KULOUT,'(''         Cvv  = '',E14.7)') RCVV
  WRITE(KULOUT,'(''      Rd/Cpd  = '',E14.7)') RKAPPA
  WRITE(KULOUT,'(''     Rv/Rd-1  = '',E14.7)') RETV
  WRITE(KULOUT,'('' *** Thermodynamic, liquid  ***'')')
  WRITE(KULOUT,'(''         Cw   = '',E14.7)') RCW
  WRITE(KULOUT,'('' *** thermodynamic, solid   ***'')')
  WRITE(KULOUT,'(''         Cs   = '',E14.7)') RCS
  WRITE(KULOUT,'('' *** Thermodynamic, trans.  ***'')')
  WRITE(KULOUT,'('' Fusion point  = '',E14.7)') RTT
  WRITE(KULOUT,'('' RTT-Tx(ew-ei) = '',E14.7)') RDT
  WRITE(KULOUT,'(''        RLvTt  = '',E14.7)') RLVTT
  WRITE(KULOUT,'(''        RLsTt  = '',E14.7)') RLSTT
  WRITE(KULOUT,'(''        RLv0   = '',E14.7)') RLVZER
  WRITE(KULOUT,'(''        RLs0   = '',E14.7)') RLSZER
  WRITE(KULOUT,'(''        RLMlt  = '',E14.7)') RLMLT
  WRITE(KULOUT,'('' Normal press. = '',E14.7)') RATM
  WRITE(KULOUT,'('' Latent heat :  '')')
  WRITE(KULOUT,'(10(1X,E11.4))') (10.*J,J=-4,4)
  WRITE(KULOUT,'(10(1X,E11.4))') (RLV(RTT+10.*J),J=-4,4)
  WRITE(KULOUT,'(10(1X,E11.4))') (RLS(RTT+10.*J),J=-4,4)
  WRITE(KULOUT,'('' *** Thermodynamic, satur.  ***'')')
  WRITE(KULOUT,'('' Fusion point = '',E14.7)') RTT
  WRITE(KULOUT,'(''      es(Tt)  = '',E14.7)') RESTT
  WRITE(KULOUT,'('' es(T) :  '')')
  WRITE(KULOUT,'(10(1X,E11.4))') (10.*J,J=-4,4)
  WRITE(KULOUT,'(10(1X,E11.4))') (ESW(RTT+10.*J),J=-4,4)
  WRITE(KULOUT,'(10(1X,E11.4))') (ESS(RTT+10.*J),J=-4,4)
  WRITE(KULOUT,'(10(1X,E11.4))') (ES (RTT+10.*J),J=-4,4)
ENDIF

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE SUCST






