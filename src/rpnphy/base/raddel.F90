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
!**s/r raddel
      subroutine raddel( del, s, sh, n, nk, nkp)
      implicit none
!!!#include <arch_specific.hf>

      integer n, nk, nkp
      real del(n,nk), s(n,nkp), sh(n,nk)
!
!Author:
!          Marc Gagnon
!
!Arguments:
!          - output -
! del      Thickness between levels.
! s        Flux sigma levels.
!
!          - input -
! sh       Sigma levels.
! n        Number of points.
! nk       Number of levels.
! nkp      Number of levels including ground.
!
!Object:
!          To compute the thickness between levels and flux levels.
!
!*
      integer k,i

      do i=1,n
         s(i,1)=2.*sh(i,1)-((sh(i,1)+sh(i,2))*0.5)

!        s(i,1)=amax1(s(i,1),0.0003)

         s(i,1)=amax1(s(i,1),sh(i,1)/2.)
         s(i,nkp)=1.0
      enddo

      do k=2,nk
         do i=1,n
            s(i,k)=(sh(i,k)+sh(i,k-1))*0.5
            del(i,k-1)=s(i,k)-s(i,k-1)
         enddo
      enddo

      do i=1,n
         del(i,nk)=1.-s(i,nk)
      enddo

      end
