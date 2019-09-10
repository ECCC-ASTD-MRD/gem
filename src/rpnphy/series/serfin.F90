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
!-------------------------------------- LICENCE END --------------------------

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

   integer I,J
   integer IELAU,IELAV,IELATU
   integer IELATV,IELAGU,IELAGV
   integer IELAUD,IELAVD,IELAWS,IELAWD
   save    IELAU,IELAV,IELATU
   save    IELATV,IELAGU,IELAGV
   real EPSIL,X,Y,PSIMLON,HNOR
   real DEG2RAD, THETA

   real, dimension(1,nt)  :: lat1, lon1, x1, y1
   real, dimension(nk,nt) :: lat2, lon2, x2, y2

   !***********************************************************************

   DEG2RAD = PI/180.

   !  UU , VV , TU , TV , VE
   IELAU=INDSERI('UUWE',PROFILS,NPROF)
   IELAV=INDSERI('VVSN',PROFILS,NPROF)
   IELATU=INDSERI('TU',PROFILS,NPROF)
   IELATV=INDSERI('TV',PROFILS,NPROF)
   IELAGU=INDSERI('GU',PROFILS,NPROF)
   IELAGV=INDSERI('GV',PROFILS,NPROF)

   if (MODELE == 'EFR') then

      !  FORMULES PRISES DANS 'LONGITUDE LATITUDE GRIDS'(D.ROBERTSON,'78)

      EPSIL=1.E-30

      if ( IELAU*IELAV.gt.0 ) then
         do J=1,NT
            do I=1,NK
               X=VV(I,J,IELAU)
               Y=VV(I,J,IELAV)
               PSIMLON=-(DGRW+LON)*DEGRAD+atan2(Y+EPSIL,X+EPSIL)
               HNOR=sqrt(X**2+Y**2)
               VV(I,J,IELAU)=HNOR*sin(PSIMLON)
               VV(I,J,IELAV)=-HNOR*cos(PSIMLON)
            enddo
         enddo
      endif
      if ( IELATU*IELATV.gt.0 ) then
         do J=1,NT
            do I=1,NK
               X=VV(I,J,IELATU)
               Y=VV(I,J,IELATV)
               PSIMLON=-(DGRW+LON)*DEGRAD+atan2(Y+EPSIL,X+EPSIL)
               HNOR=sqrt(X**2+Y**2)
               VV(I,J,IELATU)=HNOR*sin(PSIMLON)
               VV(I,J,IELATV)=-HNOR*cos(PSIMLON)
            enddo
         enddo
      endif
      if ( IELAGU*IELAGV.gt.0 ) then
         do J=1,NT
            do I=1,NK
               X=VV(I,J,IELAGU)
               Y=VV(I,J,IELAGV)
               PSIMLON=-(DGRW+LON)*DEGRAD+atan2(Y+EPSIL,X+EPSIL)
               HNOR=sqrt(X**2+Y**2)
               VV(I,J,IELAGU)=HNOR*sin(PSIMLON)
               VV(I,J,IELAGV)=-HNOR*cos(PSIMLON)
            enddo
         enddo
      endif

   else if (MODELE == 'GEF') then

      if (IELAU*IELAV.gt.0) then

         do J=1,NT
            do I=1,NK
               LAT2(I,J) = LAT
               LON2(I,J) = LON
               X2  (I,J) = VV(I,J,IELAU)
               Y2  (I,J) = VV(I,J,IELAV)
            end do
         end do

         call PLLWFGFW(X2,Y2,LAT2,LON2,NT*NK,1, &
              'E',IG(1),IG(2),IG(3),IG(4))

         do J=1,NT
            do I=1,NK
               THETA = PI/2 - Y2(I,J)*DEG2RAD
               VV(I,J,IELAU) = -X2(I,J)*cos(THETA)
               VV(I,J,IELAV) = -X2(I,J)*sin(THETA)
            end do
         end do

      endif

      if (IELATU*IELATV.gt.0) then


         do J=1,NT
            do I=1,NK
               LAT2(I,J) = LAT
               LON2(I,J) = LON
               X2  (I,J) = VV(I,J,IELATU)
               Y2  (I,J) = VV(I,J,IELATV)
            end do
         end do

         call PLLWFGFW(X2,Y2,LAT2,LON2,NT*NK,1, &
              'E',IG(1),IG(2),IG(3),IG(4))

         do J=1,NT
            do I=1,NK
               THETA = PI/2 - Y2(I,J)*DEG2RAD
               VV(I,J,IELATU) = -X2(I,J)*cos(THETA)
               VV(I,J,IELATV) = -X2(I,J)*sin(THETA)
            end do
         end do

      endif

      if (IELAGU*IELAGV.gt.0) then

         do J=1,NT
            do I=1,NK
               LAT2(I,J) = LAT
               LON2(I,J) = LON
               X2  (I,J) = VV(I,J,IELAGU)
               Y2  (I,J) = VV(I,J,IELAGV)
            end do
         end do

         call PLLWFGFW(X2,Y2,LAT2,LON2,NT*NK,1, &
              'E',IG(1),IG(2),IG(3),IG(4))

         do J=1,NT
            do I=1,NK
               THETA = PI/2 - Y2(I,J)*DEG2RAD
               VV(I,J,IELAGU) = -X2(I,J)*cos(THETA)
               VV(I,J,IELAGV) = -X2(I,J)*sin(THETA)
            end do
         end do

      endif

   endif

   !    Time-serie values of surface variables requested : either the pair
   !    UDWE, VDSN or WSPD, WD

   IELAUD=INDSERI('UDWE',SURFACE,NSURF)
   IELAVD=INDSERI('VDSN',SURFACE,NSURF)
   IELAWS=INDSERI('WSPD',SURFACE,NSURF)
   IELAWD=INDSERI('WD',  SURFACE,NSURF)

   if (MODELE.eq.'GEF') then

      if (IELAUD*IELAVD.gt.0) then
         do J=1,NT
            LAT1(1,J) = LAT
            LON1(1,J) = LON
            X1  (1,J) = VS(J,IELAUD)
            Y1  (1,J) = VS(J,IELAVD)
         end do

         call PLLWFGFW(X1,Y1,LAT1,LON1,1,NT, &
              'E',IG(1),IG(2),IG(3),IG(4))

         do J=1,NT
            THETA = PI/2 - Y1(1,J)*DEG2RAD
            VS(J,IELAUD) = -X1(1,J)*cos(THETA)
            VS(J,IELAVD) = -X1(1,J)*sin(THETA)
         end do

      endif
      if ( IELAWS*IELAWD.gt.0 ) then
         do J=1,NT
            LAT1(1,J) = LAT
            LON1(1,J) = LON
            X1  (1,J) = VS(J,IELAUD)
            Y1  (1,J) = VS(J,IELAVD)
         end do

         call PLLWFGFW(X1,Y1,LAT1,LON1,1,NT, &
              'E',IG(1),IG(2),IG(3),IG(4))

         do J=1,NT
            VS(J,IELAWS) = X1(1,J)
            VS(J,IELAWD) = Y1(1,J)
         end do

      endif

   endif

   return
end subroutine serfin2
