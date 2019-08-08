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

subroutine compute_re1(REI, REW, &
     lwc, iwc, aird, &
     mg, ml, &
     lmx, nk)
   use tdpack_const
   implicit none
!!!#include <arch_specific.hf>
#include "nbsnbl.cdk"

   integer lmx, nk, ioptrew
   real mg(lmx), ml(lmx)
   real lwc(lmx,nk), iwc(lmx,nk)
   real aird(lmx,nk)
   real rew(lmx,nk), rei(lmx,nk)
   parameter(ioptrew=1)

   !@Author Paul Vaillancourt / J. Li (Oct 2003)
   !@Revisions
   ! 001      Jason Cole (April 2006) - Modified to just compute the
   !                    liquid and ice effective radii
   !@Object
   !          A new scheme to calculate cloud optical properties,
   !          cloud microphysical input based on cldoptx4
   !@Arguments
   !          - Output -
   ! rei      effective radius for ice clouds
   ! rew      effective radius for liquid water clouds
   !
   !          - Input -
   ! LWC      in-cloud liquid water content in g/m^3
   ! IWC      in-cloud ice water content in  g/m^3
   ! MG       ground cover (ocean=0.0,land <= 1.)  (LMX)
   ! ML       fraction of lakes (0.-1.) (LMX)
   ! LMX      number of profiles to compute
   ! NK       number of layers
   !
   !***********************************************************************

   real, save :: third

   integer i, k
   real rec_grav
   real zrieff(lmx,nk), zrieff_vs(lmx,nk)
   real epsilon,epsilon2,betan,betad
   real rec_cdd(lmx,nk), vs1(lmx,nk)

#include "cldop.cdk"

   data third /0.3333333/

   rec_grav=1./grav


   ! COMPUTE THE EFFECTIVE RADIUS FOR LIQUID AND ICE CLOUD PARTICLES


   do K = 1, NK
      do I =1, LMX

         !...        EFFECTIVE RADIUS FOR WATER CLOUDS, SET NUMBER OF DROPS PER
         !           CM^3, 100 FOR WATER AND 500 FOR LAND

         if (MG(I) .le. 0.5 .and. ML(I) .le. 0.5)                then
            !              CDD=50.
            REC_CDD(I,K) =  0.01
         else
            !              CDD=250.
            REC_CDD(I,K) =  0.002
         endif

         !           calcul prealable au vspown de ZRIEFF
         !           Units of IWC must be g/m3
         !           for parameterization of REI

         ZRIEFF_VS(I,K) = IWC(I,K) !1000. * IWC(I,K) * AIRD(I,K)
      end do
   end do


   if (IOPTREW .eq. 1)                                           then
      do K = 1, NK
         do I = 1, LMX
            VS1(I,K) = (LWC(I,K)/1000.0) * REC_CDD(I,K) !* AIRD(I,K) * REC_CDD(I,K)
         enddo
      enddo
      call VSPOWN1(REW, VS1, THIRD, NK * LMX)

      !...    THIS PARAMETERIZATION FROM H. BARKER, BASED ON AIRCRAFT DATA
      !       RANGE 4-17 MICRON IS THAT SPECIFIED BY SLINGO FOR
      !       PARAMETERIZATIONS

      do K = 1, NK
         do I = 1, LMX
            REW(I,K) = min ( max (4., 754.6 * REW(I,K)), 17.0)
         end do
      end do
   elseif (IOPTREW .eq. 2)                                       then
      do K = 1, NK
         do I = 1, LMX
            VS1(I,K) = (1.0 + LWC(I,K)/AIRD(I,K) * 1.E4) &
                 * (LWC(I,K)/1000.0) * REC_CDD(I,K) !* AIRD(I,K) * REC_CDD(I,K)
         enddo
      enddo
      call VSPOWN1(REW, VS1,THIRD, NK * LMX)

      do K = 1, NK
         do I = 1, LMX
            REW(I,K) =  3000. * REW(I,K)
            REW(I,K) =  min (max (2.5, REW(I,K)), 50.0)
         end do
      end do
   elseif (IOPTREW .eq. 3)                                       then
      do K = 1, NK
         do I = 1, LMX
            EPSILON =  1.0 - 0.7 * exp(- 0.001 / REC_CDD(I,K))
            EPSILON2 =  EPSILON * EPSILON
            BETAD =  1.0 + EPSILON2
            BETAN =  BETAD + EPSILON2
            !            REW(I,K) = 620.3504944*((BETAN*BETAN*LWC(i,k)*aird(i,k))
            !     1                  / (BETAD / REC_CDD(I,K)) )**THIRD
            REW(I,K) = 620.3504944*((BETAN*BETAN*LWC(i,k)) &
                 / (BETAD / REC_CDD(I,K)) )**THIRD
            REW(I,K) =  min (max (2.5, REW(I,K)), 17.0)
            !            REW(I,K) =  MIN (MAX (4.0, REW(I,K)), 17.0)
         end do
      end do
   endif

   !...        EFFECTIVE RADIUS FOR ICE CLOUDS


   call VSPOWN1(ZRIEFF, ZRIEFF_VS, 0.216, NK * LMX)
   do K = 1, NK
      do I = 1, LMX
         if (IWC(I,K) .ge. 1.E-9)                              then
            ZRIEFF(I,K) = 83.8 * ZRIEFF(I,K)
         else
            ZRIEFF(I,K) = 20.
         endif
         REI(I,K) =  max (min (ZRIEFF(I,K), 50.0), 20.0)
      end do
   end do

   return
end subroutine compute_re1
