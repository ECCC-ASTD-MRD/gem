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
!**S/R MFSLESMX  -  MIXED PHASE SATURATION VAPOR PRESSURE CALCULATION
!
      SUBROUTINE MFDLESMX(RES,TT,FF,DF,NI,NK)
      use tdpack_const
      implicit none
#include <arch_specific.hf>
!
      INTEGER NI,NK
      REAL RES(NI,NK),TT(NI,NK),FF(NI,NK),DF(NI,NK)
!
!Author
!          A. Plante (May 2003), based on FDLESMX from J. Mailhot
!
!Revision
!
!Object
!          To calculate mixed phase saturation vapor pressure
!
!Arguments
!
!          - Output -
! RES      mixed phase saturation vapor pressure
!
!          - Input -
! TT       temperature
! FF       ice fraction
! DF       value of derivative w/r to T for saturation w/r to ice or liquid
! NI       horizontal dimension
! NK       vertical dimension
!*
      integer i,k
      real*8, dimension(ni,nk) :: work8
      real*8, dimension(ni,nk) :: foewad
      real*8, dimension(ni,nk) :: fesid
      real*8, dimension(ni,nk) :: fesmxd

#define __FORTRAN__
#include "tdpack_func.h"

!***********************************************************************
      do k=1,nk
         do i=1,ni
            foewad(i,k)=FOEWAF(tt(i,k))
            fesid(i,k) =FESIF(tt(i,k))
         enddo
      enddo
      call vexp(foewad,foewad,ni*nk)
      call vexp(fesid,fesid,ni*nk)
      do k=1,nk
         do i=1,ni
            foewad(i,k)=aerk1w*foewad(i,k)
            fesid(i,k) =aerk1i*fesid(i,k)
            fesmxd(i,k)=FESMXX(ff(i,k),fesid(i,k),foewad(i,k))
         enddo
      enddo
      do k=1,nk
      do i=1,ni
      res(i,k)=&
      FDLESMXX(tt(i,k),ff(i,k),df(i,k),foewad(i,k),fesid(i,k),fesmxd(i,k))
      enddo
      enddo
!
      RETURN
      END
