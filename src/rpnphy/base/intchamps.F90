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
!**s/r intchamps
      subroutine intchamps( ch, chrd, s, srd, n, nk, nkrd)
      implicit none
!!!#include <arch_specific.hf>

      integer n, nk, nkrd
      real chrd(n,nkrd), ch(n,nk), s(n,nk), srd(n,nkrd)
!
!Author:
!          Marc Gagnon
!
!Arguments:
!          - output -
! ch       Interpolated field.
!
!          - input -
! chrd     Reduced field.
! s        Full field sigma levels.
! srd      Reduced field sigma levels.
! n        Number of points.
! nk       Number of levels in ch and s.
! nkrd     Number of reduced levels in chrd and srd.
!
!Object:
!          Interpolate a work field from a reduced set of sigma
!          levels to a full set of sigma levels.
!*

      integer nimax
      parameter (nimax=1000)

      integer i, k, kk(nimax)
      real delt(nimax)

!     La position courante dans les niveaux reduits en dans kk.
!     Chaque point est independant, comme dans s et srd.

      do i=1,n
        kk(i)=1
        delt(i)=(chrd(i,2)-chrd(i,1))/(srd(i,2)-srd(i,1))
      enddo

!     Interpoler les points a tous les niveaux a partir des niveaux
!     reduits.  Quand le niveau courant depasse le niveau reduit en
!     direction du sol, on avance au prochain niveau reduit.  Rendu
!     au sol, on avance pas mais la valeur est calculee correctement
!     car le dernier niveau est le meme pour ch et chrd.

      do k=1,nk
        do i=1,n

          if( s(i,k) .ge. srd(i,kk(i)+1) .and. kk(i) .lt. nkrd-1 ) then
            kk(i)=kk(i)+1
            delt(i)=(chrd(i,kk(i)+1)-chrd(i,kk(i)))/ &
                    (srd(i,kk(i)+1)-srd(i,kk(i)))
          endif

!         La valeur est proportionnelle a la difference de niveau.

          ch(i,k)=chrd(i,kk(i))+delt(i)*(s(i,k)-srd(i,kk(i)))
        enddo
      enddo

      end
