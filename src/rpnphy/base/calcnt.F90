!-------------------------------------- LICENCE BEGIN -------------------------                     version 3; Last Modified: May 7, 2008.
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

subroutine calcNT(liqwpin,icewpin,cloud,ecc,ni,nk,nkp)
   implicit none
!!!#include <arch_specific.hf>

   integer, intent(in) :: ni,nk,nkp
   real liqwpin(ni,nk), icewpin(ni,nk), cloud(ni,nk)

   !@AUTHOR
   !     P. Vaillancourt  (Dec 2008)
   !
   !@OBJECT
   !     Reproduce variable NT - effective cloud cover - as it is when using the newrad radiative transfer scheme
   !     See code in cldoptx4 from L. Garand for more details
   !
   !@ARGUMENTS
   !          - Output -
   ! ECC     effective cloud amount
   !
   !          -Input -
   !     liqwpin  in-cloud liquid water path (g/m^2)
   !     icewpin  in-cloud ice    water path (g/m^2)
   !     cloud    layer cloud amount (0. to 1.) (LMX,NK)
   !     ni  horizontal dimension
   !     nk  number of layers
   !     nkp  number of layers +1

   !***********************************************************************
   real, dimension(ni,nk  ) :: ew
   real, dimension(ni,nk  ) :: ei
   real, dimension(ni,nk  ) :: eneb
   real, dimension(ni) :: trmin
   real, dimension(ni) :: tmem
   real, dimension(ni) :: ecc
   real, dimension(ni,nkp) :: ff

   integer :: i,k,kind
   real :: emiss,xnu

   ! diffusivity factor of Elsasser
   real, parameter :: ELSA = 1.66
   real, parameter :: REI = 15.
   real, parameter :: REC_REI = 1. / REI
   real, parameter :: KI = .0003 + 1.290 * REC_REI
 
   do k=1,nk
      do I=1,ni
         EW(i,k) = 1. - exp(-0.087 * ELSA * liqwpin(i,k))
         EI(i,k) = 1. - exp(-ELSA * KI * icewpin(i,k))
         EMISS   = 1. - (1.-EI(i,k)) * (1.-EW(i,k))
         ENEB(i,k) = cloud(i,k)*emiss
      end do
   end do
   
   !...  maximum random overlap

   do i=1,ni
      ff(i,1) = 1.
      tmem(i) = 1.
      trmin(i) = 1.
   enddo
   do k=2,nkp
      kind = k-2
      kind = max0(kind,1)
      do i=1,ni
         xnu = 1.-eneb(i,k-1)
         if (cloud(i,kind) < 0.01) then
            tmem(i) = ff(i,k-1)
            trmin(i) = xnu
         else
            trmin(i) = min(trmin(i), xnu)
         endif
         ff(i,k) = tmem(i) * trmin(i)
      enddo
   enddo

   do  i=1,ni
      ecc(i) = 1.-ff(i,nkp)
   enddo

   return
end subroutine calcNT
