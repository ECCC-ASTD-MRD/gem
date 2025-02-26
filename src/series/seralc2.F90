
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
