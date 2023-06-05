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

subroutine rigrad1b(RI,GAMA,GAMAQ,TBL,DUDZ2,T,TVE,Q,QE, &
     SIGMA, WW,  N, M, NK)
   use tdpack
   implicit none
!!!#include <arch_specific.hf>
   integer N,M,NK
   real RI(N,NK),GAMA(N,NK),GAMAQ(N,NK),TBL(N)
   real DUDZ2(N,NK),T(M,NK),TVE(N,NK),Q(M,NK),QE(N,NK)
   real SIGMA(N,NK)
   real WW(N)
   real TE,VIRCOR,DZ,TVK,TVKP,FAC,DLNP
   integer j,k

   !@Author J. Mailhot RPN(Feb 1990)

   !@Revision
   ! 001      J.Mailhot RPN(Feb 1990)shallow convection (GELEYN)
   ! 002      G.Pellerin(August90)Adaptation to thermo functions
   ! 003      Y. Delage (Nov 1990)Options of shallow convection
   ! 004      Y. Delage  (Nov 1990)- Removal of BETA and PRI
   ! 005      N. Brunet  (May91)
   !                New version of thermodynamic functions
   !                and file of constants
   ! 006      C. Girard (Nov92) New parameterization of the
   !          shallow convection
   ! 007      G. Pellerin(May95) Extract level of BLH
   ! 008      C. Girard (Nov95) Added calculations of GAMAQ
   ! 009      A. Plante (June 2003) - IBM conversion
   !             - calls to vslog routine (from massvp4 library)
   ! 010      B. Bilodeau (May 2005) - New comdeck fintern
   ! 011      L. Spacek (Dec 2007) - Add calculation of gama at nk

   !@Object calculate the gradient Richardson number

   !@Arguments
   !          - Outputs -
   ! RI       gradient Richardson number
   ! GAMA     equilibrium gradient term for temperature
   ! GAMAQ    equilibrium gradient term for moisture
   ! TBL      Level corresponding to top of Unstable boundary layer
   !          - Inputs -
   ! DUDZ2    vertical shear of the wind squared
   ! T        temperature
   ! TVE      virtual temperature at 'E' levels
   ! Q        specific humidity
   ! QE       specific humidity at 'E' levels
   ! SIGMA    sigma level values
   ! WW       work field
   ! N        horizontal dimension
   ! M        1st dimension of T and Q
   ! NK       vertical dimension


   do j = 1, N
      WW(j) = FOTVT( T(j,1), Q(j,1) )
   end do

   do k = 1, NK - 1
      do j = 1, N

         RI(j,k) = log(SIGMA(j,k+1)/SIGMA(j,k))
         
         !           TEMPERATURES VIRTUELLES

         TVK = WW(j)
         WW(j) = FOTVT( T(j,k+1), Q(j,k+1) )
         TVKP = WW(j)
         TE = FOTTV( TVE(j,k), QE(j,k) )
         VIRCOR = TVE(j,k) / TE

         DLNP =  RI(j,k)
         DZ = RGASD * TVE(j,k) * DLNP / GRAV

         !           RI

         FAC = GRAV / ( TVE(j,k) * (max(DUDZ2(j,k),1.e-10)) )
         RI(j,k) = FAC * ( ( TVK - TVKP ) / DZ + GRAV/CPD )

         !           GAMA

         GAMA(j,k) = - GRAV / ( CPD * VIRCOR )

         GAMAQ(j,k) = 0.

         !           TOP OF THE unstable BOUNDARY LAYER

         if( RI(j,k).gt.0. ) TBL(j) = k

      end do
   end do

   do j = 1, N

      RI(j,NK) = RI(j,NK-1)
      TE = FOTTV( TVE(j,NK), QE(j,NK) )
      VIRCOR = TVE(j,NK) / TE
      GAMA(j,NK) = - GRAV / ( CPD * VIRCOR )

   end do

   return
end subroutine rigrad1b
