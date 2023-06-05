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
!-------------------------------------- LICENCE END ----------------------------

subroutine SERALC2(NI,NJ,NK)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
   integer NI,NJ,NK

   !@Author M. Lepine  -  RFE model code revision project (Feb 87)
   !@Revision
   ! 001      B. Bilodeau  (July 1991)- Adaptation to UNIX
   ! 002      B. Bilodeau  (July 2002)- Initialize SURFACE(M,2) and
   !                                    PROFILS(M,2)
   !@Object INITIALISATION DE CERTAINES VARIABLES DES SERIES TEMPORELLES
   !        to initialize certain variables for time-series
   !@Arguments
   !          - Input -
   ! NI       1st horizontal dimension
   ! NJ       2nd horizontal dimension
   ! NK       vertical dimension

#include "series.cdk"

   integer i,m,l,k

      NINJNK(1) = NI
      NINJNK(2) = NJ
      NINJNK(3) = NK
      NSTAT   = 0
      nstat_g = 0
      NSURF   = 0
      NPROF   = 0
      do 10 I = 1,MXSTT
        IJSTAT(I,1) = 0
        JSTAT(I)    = 0
        istat_g(i)  = 0
        jstat_g(i)  = 0
        statnum(i)  = 0
        name(i)     = ''
  10  continue

      do M=1,MXSRF
         SURFACE(M,1) = '        '
         SURFACE(M,2) = '        '
         do L=1,MXSTT
            SERS(L,M) = 0.0
         end do
      end do

      do M=1,MXPRF
         PROFILS(M,1) = '        '
         PROFILS(M,2) = '        '
         do L=1,MXSTT
         do K=1,MXNVO
            SERP(K,L,M) = 0.0
         end do
         end do
      end do

      return
      end
