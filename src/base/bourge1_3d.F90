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
!-------------------------------------- LICENCE END --------------------------

subroutine bourge1_3d(fneige,fip,t,s,ps,ni,nk)
   use tdpack_const, only: TCDK, RGASD
   implicit none
!!!#include <arch_specific.hf>
   integer, intent(in) :: ni,nk
   real fneige(ni,nk+1),fip(ni,nk+1),t(ni,nk),s(ni,nk),ps(ni)

!AUTHOR
!          Andre Plante, April 2006, based on :
!          Andre Methot, Andre Plante (2002)
!
!Revision
! 001      Andre Plante (April 2006) - Initial version
! 002      Anne-Marie Leduc . Andre Plante (November 2007) - add parameter ps and change name to bourge1_3d.
!                             Do calculation of fice and fip when pressure is above 200 mb.
!
!Arguments
!
!          - Output -
! fneige   Fraction of precipitation that is in form of snow.
! fip      Fraction of liquid precipitation that is in the form of ice pellets.
!
!          - Input -
! t        temperature
! s        sigma (p/ps) or pressure levels
! ps       surface pressure
! ni       horizontal dimension
! nk       vertical dimension
!
!Object
!          bourge_3d differs from bourge1 only for the output of 3D snow and ice pellets fractions.
!
!          This routine produces a diagnostic for discriminating precipitation types.
!          Based on temperature and pressure (sigma) profiles, two fractions are computed for
!          all levels. The snow fraction, "fneige", gives the fraction of the precipitation
!          that should be falling as snow. Of the remaining portion (1-fneige) there is a
!          fraction "fip" that is refrozen (ice pellels), the rest is rain or freezing rain
!          depending of the temperature. Note that the actual precipitation types are not
!          computed here but only the two fractions necessary to compute these 4 types.
!
!Notes:
!          This diagnostic is a 3-D extension of Bourgouin's method
!          (Wea. and Forcasting. 2000,vol15).
!          There is some discrepencies between this 3-D method and the original.
!
!          Some mixtures of precipitation types not allowed in original method are allowed here
!          ( Freezing rain and snow,  Ice pellets and snow, etc).   This is then an approximation
!          that is going outside of the statistical signifiance of the original work.
!
!          Also, statistical thresholds for snow melting in a warm layer are always the same in this
!          3-D method, regardless of the profil below the warm layer.   On the other hand, Bourgouin
!          2000 suggest different thresholds for snow melting in the warm layer according to the full
!          sounding categorical situation ( rain-snow vs freezing rain for instance ).
!          The statistical threshold used here correspond to the rain-snow findings from Bourgouin
!          in all cases.  Strictly speaking, this is going outside of statistical singifiance from
!          the original work from Bourgouin.
!
!
!     Physical constants & parameters

!     Bourgouin's parametres (Wea. Forecasting 2000, 15, pp 583-592)
      real m1,m2,f1,f2,fslope
      parameter(m1=5.6,m2=13.2,f1=46.,f2=66.,fslope=.66)

      integer i,k,warm(ni)
      real wa(ni),ca(ni),wrk1(ni,nk),area,seuil1,seuil2,mfactor
      real fice(ni),ficetop(ni),riptop(ni),rip(ni)
      real wa3d(ni,nk)
      real press

!*

!     Factor for linear interpolation in the context of partial melting.
      mfactor=1./(m2-m1)

      do k=2,nk-1
         do i=1,ni
            wrk1(i,k)=log( (s(i,k+1)+s(i,k)) / (s(i,k)+s(i,k-1)) )
         enddo
      enddo

      do i=1,ni
         wrk1(i,nk)=log( (1.+s(i,nk)) / (s(i,nk)+s(i,nk-1)) )
         fneige(i,1)= 1.
         fneige(i,2)= 1.
         fip(i,1)= 0.
         fip(i,2)= 0.
         fice(i)= 1.
         wa(i)  = 0.
         ca(i)  = 0.
         rip(i) = 0.
         warm(i)= 0
      enddo

      do k=2,nk
         do i=1,ni


!        do calculation of fneige and fip only for points where press is greater or equal to 200 mb or 20000 pa.
!        otherwise problems with the meso strato model where above zero degree layers in the
!        stratosphere give liquid precipitation at the surface.

          press = s(i,k)*ps(i)
!___________________________________________________________________________

          if (press.lt.20000) then

            fneige(i,k+1)= 1.
            fip   (i,k+1)= 0.

          else

            if ( t(i,k) .gt. tcdk ) then
!------------------------------------------WARM LAYER-------------------
               if (warm(i).eq.0)then
!                 Debut de la couche chaude
                  warm(i)=1
                  ficetop(i)=fice(i)
                  riptop(i)=rip(i)
               endif

!              Reduction of total warm area due to cold area above.
               if ( ca(i) .gt. 0. ) &
                    wa(i) = max(0.,wa(i) - rip(i)*ca(i))

               area= rgasd * (t(i,k)-tcdk) * wrk1(i,k)
               wa(i)  = wa(i) + area

               wa3d(i,k)=wa(i)

!              if wa(i) < m1 then warm area is too small for metling, therefore fice do not change.
!              if wa(i) > m2 area large enough for complet melting
!              else partial melting (lineair interpolation between m1 and m2)
!
!                            Partial melting
!                 ^                 |
!              1 -|------------.    v   .
!                 |            .\       .
!              f  |            . \      .
!              a  |            .  \     .
!              c  | fice do    .   \    .  fice = 0
!              t  | no change  .    \   .  (complete melting)
!              o  |            .     \  .
!              r  |            .      \ .
!                 |            .       \.
!              0 -------------------------------------------->
!                              |        |                wa
!                              m1       m2
!
!              The melting factor is computed as fallow:

               area=min( 1. , max( 0. , 1. - (wa(i)-m1)*mfactor ) )

!              We melt snow (fice) and ice pellets (rip)

               fice(i)= ficetop(i)* area
               rip(i) = riptop(i) * area

!              Forget about above cold layer.
               ca(i) = 0.

            else if ( wa(i) .gt. .00001 ) then
!--------------------------------------------COLD LAYER BELOW WARM LAYER
               if (warm(i).eq.1)then
!                 Debut de la couche froide
                  warm(i)=0
                  riptop(i)=rip(i)
               endif

               area= rgasd *(tcdk-t(i,k))* wrk1(i,k)
               ca(i)  = ca(i) + area

               wa3d(i,k)=-ca(i)

               seuil1= f1 + fslope * wa(i)
               seuil2= f2 + fslope * wa(i)

!                            Partial freezing
!                 ^                |
!              1 -|            .   v    .-------------------
!                 |            .       /.
!              f  |            .      / .
!              a  |            .     /  .
!              c  |            .    /   .
!              t  | no freezing.   /    .  complete freezing
!              o  |            .  /     .
!              r  |            . /      .
!                 |            ./       .
!              0 -------------------------------------------->
!                              |        |               ca
!                           seuil1   seuil2
!
!              The factor is computed as fallow:

               area=min( 1. , &
                    max( 0. , (ca(i)-seuil1)/(seuil2-seuil1) ) )

!              We make ice pellets from available liquid water fraction.

               rip(i) = riptop(i) + (1.-fice(i)-riptop(i))*area

            endif

            fneige(i,k+1)=fice(i)
            if ( fice(i) .gt. 0.9999 ) then
               fip(i,k+1) = 0.
            else
               fip(i,k+1) = rip(i) / ( 1. - fice(i) )
            endif

       endif
!________________________________________________________________________________________________

         enddo
      enddo

      end
