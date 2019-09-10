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
!**s/p rdmoy
      subroutine rdmoy( chrd, ch, vtot, inrd, ni, nk, nkrd)
      implicit none
!!!#include <arch_specific.hf>

      integer ni, nk, nkrd
      integer inrd(nkrd)
      real chrd(ni,nkrd), ch(ni,nk), vtot(ni)

!Author:
!          Marc Gagnon
!
!Arguments:
!
!          - Output -
! chrd     array after reduction
!
!          - Input -
! ch       original array before reduction
! vtot     work field
! inrd     list of reduced levels
! ni       horizontal dimension
! nk       number of levels
! nkrd     reduced number of levels
!
!Object:
!          select a subset of the levels of a 2D (ni by nk) array,
!          keeping the average values of the rejected levels
!
!Notes:
!          ni .gt. 0
!          nkrd .le. nk
!          inrd(*) is included in [1..nk]
!*

      integer krd, kc, ks, k, i, n

      do krd=1,nkrd
        kc = inrd(krd)

        if( krd .ne. nkrd) then
          ks = inrd(krd+1)
        else
          ks = kc
        endif

        do i=1,ni
          vtot(i) = ch(i,kc)
        enddo

        do k=kc+1,ks-1
          do i=1,ni
            vtot(i) = vtot(i) + ch(i,k)
          enddo
        enddo

        n = (ks-1) - (kc+1) + 2

        if( n .lt. 1 ) then
          n = 1
        endif

        do i=1,ni
          chrd(i,krd) = vtot(i) / n
        enddo

      enddo

      end
