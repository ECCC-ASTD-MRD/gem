
subroutine serfin2(VS,VV,SURFACE,PROFILS,NT,NK,NSURF,NPROF, &
     DGRW, LAT, LON, &
     DEGRAD, MODELE, IG)
   use tdpack_const, only: PI
   implicit none
!!!#include <arch_specific.hf>

   integer NSURF,NK,NT,NPROF
   real VS(NT,NSURF),VV(NK,NT,NPROF)
   character(len=*) SURFACE(*),PROFILS(*)
   real LAT,LON
   real DGRW, DEGRAD
   character(len=3) MODELE
   integer IG(4)

   !@Author R. Benoit
   !@Revision
   ! 001      V.Alex.(Feb 87) Documentation
   ! 002      N. Brunet  (May90)
   !                Standardization of thermodynamic functions
   ! 003      N. Brunet  (May91)
   !                New version of thermodynamic functions
   !                and file of constants
   ! 004      B. Bilodeau  (July 1991)- Adaptation to UNIX
   ! 005      B. Bilodeau  (Jan94) - Add rotation of GU and GV
   ! 006      B. Bilodeau  (Feb94) - Remove calculations involving
   !                                 TSTAR, QSTAR and PSTAR
   ! 007      B. Bilodeau (Feb 1997) - Rotation of winds from GEF grid
   ! 008      B. Bilodeau (Jan 2001) - Automatic arrays
   ! 009      B. Dugas    (Dec 2005) - Correct erroneous wind components
   !                rotation by calling PLLWFGFW
   !@Object normalize units of the extracted time-series data
   !        normalisation des unites des diverses series
   !@Arguments
   !          - Input/Output -
   ! VS       time-serie values of surface variables requested
   ! VV       time-serie values of profile variables requested
   !          - Input -
   ! SURFACE  names of time series surface variables requested
   ! PROFILS  names of time series profile variables requested
   ! NT       timestep number
   ! NK       vertical dimension
   ! NSURF    number of surface variables requested
   ! NPROF    number of profile variables requested
   ! DGRW     angle between the Greenwich meridian and the X
   !          axis of the model
   ! LAT      latitude of the station
   ! LON      longitude of the station
   ! DEGRAD   conversion constant from degrees to radians
   ! MODELE   model name (EFR or GEF)
   ! IG       IG1,IG2,IG3 and IG4 of the GEF grid descriptors

   integer, external :: INDSERI

   integer IELAU,IELAV,IELATU
   integer IELATV,IELAGU,IELAGV
   integer IELAUD,IELAVD,IELAWS,IELAWD
   save    IELAU,IELAV,IELATU
   save    IELATV,IELAGU,IELAGV

   !***********************************************************************

   !  UU , VV , TU , TV , VE
   IELAU=INDSERI('UUWE',PROFILS,NPROF)
   IELAV=INDSERI('VVSN',PROFILS,NPROF)
   IELATU=INDSERI('TU',PROFILS,NPROF)
   IELATV=INDSERI('TV',PROFILS,NPROF)
   IELAGU=INDSERI('GU',PROFILS,NPROF)
   IELAGV=INDSERI('GV',PROFILS,NPROF)

   if (MODELE == 'EFR' .or. MODELE == 'GEF') then
      print *, 'ERROR: (serfin) model no longer supported: '//trim(MODELE)
      print *, '---- ABORT ----'
      call flush(6)
      stop
   endif

   !    Time-serie values of surface variables requested : either the pair
   !    UDWE, VDSN or WSPD, WD

   IELAUD=INDSERI('UDWE',SURFACE,NSURF)
   IELAVD=INDSERI('VDSN',SURFACE,NSURF)
   IELAWS=INDSERI('WSPD',SURFACE,NSURF)
   IELAWD=INDSERI('WD',  SURFACE,NSURF)

   return
end subroutine serfin2
