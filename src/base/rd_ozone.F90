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
!**s/r rd_ozone  -- Perform the actual reading and organization of
!                   the data in the ozone file
!
      subroutine rd_ozone (iun, rbuf, dim, status)
      implicit none
!!!#include <arch_specific.hf>
!
      integer iun,dim(*),status
      real rbuf(*)
!
!Author
!          M. Desgagne (Spring 2008)
!
!Arguments
!          - Input -
! iun      fortran unit
!          - Input/Output -
! rbuf     read buffer
!          - Output -
! status   exit status for the routine (=0 if OK,  =-1 otherwise)
!
#include "radiation.cdk"
#include "ozopnt.cdk"
!
      integer i,ilir,mi,mj,mk,m,NLP,code
      integer fstinf,fstluk
!
!-----------------------------------------------------------------
!
      code   = status
      status = -1

      if (code.eq.200) then
         if (fstinf (iun,NLACL,NPCL,mk,-1,' ',-1,-1,1,' ','O3').lt.0) &
         return
         dim(1) = 3
         dim(2) = NLACL*NPCL*13 + NLACL + NPCL
         dim(3) = NLACL
         dim(4) = NPCL
         status = 0
         return
      endif

      NLP = NLACL*NPCL

      if ( (code.gt.200) .and. (code.le.300) )then
         ilir = fstinf (iun,mi, mj,mk,-1,' ',-1,-1,-1,' ','ZLAT')
         if (ilir.lt.0) return
         ilir = fstluk (rbuf, ilir, mi, mj, mk)
         ilir = fstinf (iun,mi, mj,mk,-1,' ',-1,-1,-1,' ','PREF')
         if (ilir.lt.0) return
         ilir = fstluk (rbuf(NLACL+1), ilir, mi, mj, mk)
         do m=1,12
            ilir = fstinf (iun,mi,mj,mk,-1,' ',-1,-1,m,' ','O3')
            if (ilir.lt.0) return
            ilir = fstluk (rbuf(NLACL+NPCL+(m-1)*NLP+1),ilir,mi,mj,mk)
         end do
         status = 0
      endif

      if (code.ge.300) then
         NLACL = dim(3)
         NPCL  = dim(4)
         NLP   = NLACL*NPCL

         allocate (gozon12(NLP,12), goz(NLP+NLACL+NPCL))

         DO i=1,NLACL
            goz(NLP      +i)=  rbuf(i)
         ENDDO

         DO i=1,NPCL
            goz(NLP+NLACL+i)= rbuf(NLACL+i)
         ENDDO
         do m=1,12
         do i=1,NLP
            gozon12(i,m) = rbuf(NLACL+NPCL+(m-1)*NLP+i)
         end do
         end do
         status = 0
      endif
!
!-----------------------------------------------------------------
!
      return
      end
