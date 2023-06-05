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

subroutine CLDWIN1(BM,FM,T,HU,PS,SIGMA,NI,NK,SATUCO)
   use, intrinsic :: iso_fortran_env, only: INT64
   use tdpack, only: FOQST, FOQSA
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer NI,NK
   real FM(NI,NK),BM(NI,NK)
   real T(NI,NK),HU(NI,NK)
   real PS(NI)
   real TRAV(NI,NK)
   real SIGMA(NI,NK)

   logical SATUCO

   !@Author Janusz Pudykiewicz (Nov. 1985)
   !@Revision
   ! 001      J. Mailhot (March 1993) 0 <= RH <= 1
   ! 002      R. Sarrazin (May 1994) Bugs correction
   !          and function funb0 func0
   ! 003      B. Bilodeau (Aug 1994) New physics interface
   ! 004      R. Sarrazin (June 95) Corrections; add cloud fraction
   ! 005      B. Bilodeau (Jan 2001) Automatic arrays
   ! 006      M. Lepine  (March 2003) -  CVMG... Replacements
   !@Object initialize the cloud water field
   !@Arguments
   !          - Output -
   ! BM       stratiform cloud fraction
   ! FM       cloud water field
   !          - Input -
   ! T        temperature
   ! HU       specific humidity
   ! PS       surface pressure
   ! TRAV     work space
   ! SIGMA    vertical discretization
   ! NI       X-(horizontal grid dimension)
   ! NK       number of vertical levels
   ! SATUCO   .TRUE. if water/ice phase for saturation
   !          .FALSE. if water phase only for saturation

   !  STATEMENT FUNCTION:
   !----------------------
   !    >      FQSAT(XT,XP)=
   !    >*
   !    >     X0.622/((XP/100.)*EXP(5418./XT - 21.656) - 0.378)

   !     USE COMDECKS FROM PHYSICS LIBRARY   ----------------
   !                                                        *
   real TEMP1,WMR

   real, dimension(NI) :: WPU0

   !***********************************************************************

   real FUNB0,FUNC0
   real YPU0,XWFLO,XSIGMA,YPSGBT,YPUSG1,YPUS,YPTCI,XTTM,XW1

   integer NIK,K,I

   real, parameter :: PU0 = 0.85
   real, parameter :: PUS = 1.0
   real, parameter :: PMR = 0.2E-3
   real, parameter :: PUSBG1 = 0.99
   real, parameter :: PSGBT = 0.8
   real, parameter :: PTCI = 238.

#include "comphy.cdk"
   include "phyinput.inc"

   !____FUNCTION USED FOR CORRECTING OF THE U00 FOR SIGMA>0.9

   FUNB0(YPU0,XWFLO,XSIGMA,YPSGBT,YPUSG1)= &
        YPU0+XWFLO*((XSIGMA-YPSGBT)/(1.-YPSGBT))* &
        (YPUSG1-YPU0)

   !____FUNCTION USED FOR CORRECTING OF THE U00 FOR CIRRUS LEVEL

   FUNC0(YPUS,YPTCI,XW1,XTTM)= &
        YPUS-((YPUS-XW1)/(1.+0.15*(YPTCI-XTTM)))



   !______COMMON VALUES

   NIK = NI * NK

   !-----------------------------------------------------------

   !      INITIALISATION OF THE CLOUD WATER CONTAIN
   !      USING EMPIRICAL RELATION BETWEEN RELATIVE
   !      HUMIDITY AND PARTIAL CLOUD COVER.

   !___1) COMPUTE RELATIVE HUMIDITY


   !     -----------------------------------------
   if (SATUCO) then
      do K=1,NK
         do I=1,NI

            TRAV(I,K) = HU(I,K) / &
                 FOQST( T(I,K) , PS(I)*SIGMA(I,K))

            TRAV(I,K) = max( 0.0 , min(TRAV(I,K),1.0) )

         enddo
      enddo
   else
      do K=1,NK
         do I=1,NI

            TRAV(I,K) = HU(I,K) / &
                 FOQSA( T(I,K) , PS(I)*SIGMA(I,K))

            TRAV(I,K) = max( 0.0 , min(TRAV(I,K),1.0) )

         enddo
      enddo
   endif

   DOK: do K=1,NK

      ! correction to pu0 near ground

      TEMP1 = 0.75
      do I=1,NI
         WPU0(I)=PU0
         ! 0.5 on continent, 1.0 over ocean
         if( SIGMA(I,K) .ge. PSGBT ) &
              WPU0(I)=FUNB0(PU0,TEMP1,SIGMA(I,K),PSGBT,PUSBG1)
      enddo

      ! correction to pu0 for low temperatures

      do I=1,NI
         if (T(I,K).le.PTCI) then
            TEMP1 = WPU0(I)
            WPU0(I) = FUNC0(PUS,PTCI,TEMP1,T(I,K))
         endif
      enddo


      !___2) COMPUTE PARTIAL CLOUD COVER



      DIAG_FTOT: if (.not.any('ftot'==phyinread_list_s(1:phyinread_n))) then
         do I=1,NI

            if ((TRAV(I,K).ge.WPU0(I)).and.(TRAV(I,K).le.PUS)) then
               BM(I,K)=((TRAV(I,K)-WPU0(I))/(PUS-WPU0(I)))
            else
               BM(I,K) = 0.
            endif

         enddo


         do I=1,NI

            BM(I,K) = 1. - sqrt( 1. - BM(I,K) )

         enddo
      endif DIAG_FTOT


      !___3) ESTIMATE CLOUD WATER CONTAIN FROM PARTIAL CLOUD
      !      COVER.

      DIAG_LWC: if (.not.any('lwc'==phyinread_list_s(1:phyinread_n))) then
         do I=1,NI

            ! reduce pmr for cold temperature

            if( T(I,K) .lt. 268. ) then
               WMR=PMR/(1.0+0.5*(268.-T(I,K))**0.5)
            else
               WMR=PMR
            endif

            if(TRAV(I,K).le.PUS) then
               FM(I,K)=.9*BM(I,K)*WMR
            else
               FM(I,K)=   BM(I,K)*WMR
            endif

         enddo
      endif DIAG_LWC

   enddo DOK
   !--------------------------------------------------------
   return
end subroutine CLDWIN1
