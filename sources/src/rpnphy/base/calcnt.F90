!-------------------------------------- LICENCE BEGIN ------------------------------------
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
!-------------------------------------- LICENCE END --------------------------------------
!** S/P CALCNT
!
      subroutine calcNT(liqwpin,icewpin,cloud,ecc,ni,nk,nkp)
!
      implicit none
#include <arch_specific.hf>
!
      integer ni,nk,nkp
      real liqwpin(ni,nk), icewpin(ni,nk), cloud(ni,nk)
!
!AUTHOR
!     P. Vaillancourt  (Dec 2008)
!
!REVISION
!
! 001
!
!OBJECT
!     Reproduce variable NT - effective cloud cover - as it is when using the newrad radiative transfer scheme
!     See code in cldoptx4 from L. Garand for more details
!
!ARGUMENTS
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
!
!MODULES
!
!
!***********************************************************************
!     AUTOMATIC ARRAYS
!***********************************************************************
!
      real, dimension(ni,nk  ) :: ew
      real, dimension(ni,nk  ) :: ei
      real, dimension(ni,nk  ) :: eneb
      real, dimension(ni) :: trmin
      real, dimension(ni) :: tmem
      real, dimension(ni) :: ecc
      real, dimension(ni,nkp) :: ff

!
!***********************************************************************
!
        integer i,k,kind
        real rei, rec_rei, ki
        real elsa, emiss,xnu
!
#include "tdpack_const.hf"
!
! diffusivity factor of Elsasser
        data elsa/1.66/
!
!
      REI = 15.
      REC_REI = 1. / 15.
      KI = .0003 + 1.290 * REC_REI
!
      CALL VSEXP (EW,-0.087*elsa*liqwpin,nk*ni)
      CALL VSEXP (EI,-elsa*ki*icewpin,nk*ni)
!
      do k=1,nk
      do I=1,ni
            EW(i,k) = 1. - EW(i,k)
            EI(i,k) = 1. - EI(i,k)
            EMISS = 1. - (1.-EI(i,k))* (1.-EW(i,k))
            ENEB(i,k)= cloud(i,k)*emiss
      end do
      end do
!
!
!...  maximum random overlap

      do I=1,ni
            ff(i,1)=1.
            tmem(i)=1.
            trmin(i)=1.
      enddo
      do k=2,nkp
            kind=k-2
            kind=max0(kind,1)
            do I=1,ni
               xnu=1.-eneb(i,k-1)
               if(cloud(i,kind).lt.0.01) then
                 tmem(i)= ff(i,k-1)
                 trmin(i)= xnu
               else
                 trmin(i)=min(trmin(i),xnu)
               endif
               ff(i,k)= tmem(i) * trmin(i)
            enddo
      enddo
!
      do  i=1,ni
         ecc(i)=1.-ff(i,nkp)
      enddo


      return
      end
