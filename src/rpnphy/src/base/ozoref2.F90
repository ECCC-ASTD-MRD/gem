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

!*@/
subroutine OZOREF3(O3F,LREF,DLAT,NP,NMAX,LBL,NLAT,ALAT,F)
   use tdpack, only: PI
   implicit none
!!!#include <arch_specific.hf>
   integer LREF,NP,NMAX,IB,I,J,K,L,NLAT
   real SLOPE,XLATI
   real DL,B,F1,F2,A1,A2,OZON,FAC,CONV
   real O3F(NMAX,LREF)
   integer LBL(NP)
   real DLAT(NP)
   real F(NLAT,LREF),ALAT(NLAT)

   !@Author
   !          L.Garand (1997)
   !          rewritten from original CCRN code
   !@Revisions
   ! 001      B. Bilodeau (Jan 2000) - Exit if latitudes out of bounds
   !@Object
   !          to calculate the ozone mixing ratio (kg/kg) at ozone
   !          climatological levels for desired array of latitudes
   !@Arguments
   !          - Output -
   ! O3F      ozone (kg O3/kg air) at  each reference
   !          climatological level for DLAT latitudes
   !          - Input -
   ! LREF     number of climatological ozone level
   ! DLAT     latitude in radians of model points
   ! NP       number of points to process
   ! NMAX     number of maximum points permitted
   ! LBL      work field
   ! NLAT     number of latitude climatological bands
   ! ALAT     climatological ozone latitudes in degrees
   ! F        climatological ozone field in PPMV
   !*@/

   FAC = 180./PI

   do J=1,NP
      IB=0
      XLATI= FLOAT( nint(DLAT(J)*FAC) )

      do I=1,NLAT
         if( XLATI.lt.ALAT(I) .and. IB.eq.0 ) IB=I
      enddo

      if ( XLATI.eq. 90.0 ) IB=NLAT

      if(IB.le.1) then
         write(6,6030) XLATI
         write(6,*) (ALAT(I),i=1,NLAT)
         call physeterror('ozoref', 'O3 Inpol out bounds in latitude')
         return
6030     format(1X,' O3 INPOL OUT BOUNDS IN LATITUDE:',E12.4)
      endif
      LBL(J)=IB-1
   enddo

   !  interpolate to desired latitudes
   !  and transform into kg/kg using CONV:
   !  1.E-6 converts PPMV to PPV , 48.0=M(O3), 28.964= M(dry air)

   CONV = 1.E-6*48./28.964

   do L=1,LREF
      do J=1,NP
         K=LBL(J)
         A1=ALAT(K)
         A2=ALAT(K+1)
         F1=F(K,L)
         F2=F(K+1,L)
         DL=DLAT(J)  *FAC
         SLOPE = (F2-F1)/(A2-A1)
         B = F1 - SLOPE*A1
         OZON = SLOPE*DL + B
         O3F(J,L)=OZON*CONV
      enddo
   enddo

   return
end subroutine OZOREF3
