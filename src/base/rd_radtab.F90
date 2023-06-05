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
!**s/r rd_radtab  -- Perform the actual reading and organization of
!                    the data in the radiation file
!
      subroutine rd_radtab (iun, rbuf, dim, status)
      implicit none
!!!#include <arch_specific.hf>
!
      integer iun,dim(*),status
      real rbuf(*)
!
!Author
!          M. Desgagne (Spring 2008)
!
!Object
!          Reads and organizes the data in cible according to
!          flag ozotable:   .true.  ===> for the ozone data
!                           .false. ===> for the radiation table
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
      integer ilir,inbr,mi,mj,mk,code
      integer fstinf,fstlir
!
!-----------------------------------------------------------------
!
      code   = status
      status = -1

      if (code.eq.200) then
         ilir   = fstinf (iun,mi,mj,mk,-1,' ',-1,-1,-1,'C','G1')
         if (.not.(( mi.eq.mxx.and.mj.eq.ntt) .and. (ilir.ge.0))) return
         dim(1) = 1
         dim(2) = ntotal
         status = 0
         return
      endif

      if ( (code.gt.200) .and. (code.le.300) )then
         inbr=fstlir(rbuf(g1),  iun,mi,mj,mk,-1,' ',-1,-1,-1,' ','G1')
         inbr=fstlir(rbuf(g2),  iun,mi,mj,mk,-1,' ',-1,-1,-1,' ','G2')
         inbr=fstlir(rbuf(g3),  iun,mi,mj,mk,-1,' ',-1,-1,-1,' ','G3')
         inbr=fstlir(rbuf(th2o),iun,mi,mj,mk,-1,' ',-1,-1,-1,' ','2O')
         inbr=fstlir(rbuf(tro3),iun,mi,mj,mk,-1,' ',-1,-1,-1,' ','T3')
         inbr=fstlir(rbuf(yg3), iun,mi,mj,mk,-1,' ',-1,-1,-1,' ','Y3')
         inbr=fstlir(rbuf(bcn), iun,mi,mj,mk,-1,' ',-1,-1,-1,' ','BN')
         inbr=fstlir(rbuf(dbcn),iun,mi,mj,mk,-1,' ',-1,-1,-1,' ','DN')
         inbr=fstlir(rbuf(bo3), iun,mi,mj,mk,-1,' ',-1,-1,-1,' ','B3')
         inbr=fstlir(rbuf(dbo3),iun,mi,mj,mk,-1,' ',-1,-1,-1,' ','D3')
         inbr=fstlir(rbuf(to3), iun,mi,mj,mk,-1,' ',-1,-1,-1,' ','3O')
         inbr=fstlir(rbuf(uu),  iun,mi,mj,mk,-1,' ',-1,-1,-1,' ','2U')
         inbr=fstlir(rbuf(tt),  iun,mi,mj,mk,-1,' ',-1,-1,-1,' ','2T')
         status = 0
      endif

      if (code.ge.300) then
         g      = rbuf(1:ntotal)
         status = 0
      endif
!
!-----------------------------------------------------------------
!
      return
      end
