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

!/@*
subroutine fillagg(bus,bussiz,ptsurf,ptsurfsiz,indx_sfc,surflen)
   use sfcbus_mod
   implicit none
!!!#include <arch_specific.hf>
   !@Object
   !@Arguments
   !             - Input/Output -
   ! BUS         bus for the given surface scheme
   !             - Input -
   ! BUSSIZ      dimension of bus
   ! PTSURF      surface pointers
   ! PTSURFSIZ   dimension of ptsurf
   ! INDX_SFC    integer value (1-4) corresponding to each surface type
   ! SURFLEN     horizontal dimension (row length) of the gathered points
   !             for the given surface type on the full model row
   !             To copy the fields FC, FV, TSURF, etc.
   !             into the arrays that subroutine AGREGE
   !             will use to calculate the aggregated
   !             values. This subroutine needs to be called
   !             4 times

   integer bussiz,indx_sfc, ptsurfsiz, surflen
   integer ptsurf(ptsurfsiz)
   real, target :: bus(bussiz)

   !@Author B. Bilodeau  Spring 2001
   !@Revisions
   ! 001         B. Bilodeau (Feb 2004) - Revised logic to facilitate the
   !                                      addition of new types of surface
   !*@/

   real, pointer, dimension(:) :: champ_agg, champ_sfc

   integer i,ik,j,k

   !# fonction-formule
   integer :: fn_busidx, fn_ptr, fn_j, fn_k
   fn_busidx(fn_ptr,fn_j,fn_k) = ptsurf(fn_ptr)+(fn_k-1)*surflen+fn_j-1
   !--------------------------------------------------------------
   do j=1,nvarsurf
      if (vl(j)%agg.eq.indx_agrege.and.vl(j)%mul.gt.1) then
         champ_agg => bus(fn_busidx(j,1,indx_agrege):)
         champ_sfc => bus(fn_busidx(j,1,indx_sfc   ):)
         do k=1,vl(j)%niveaux
            do i=1,surflen
               ik = (k-1)*surflen + i
               champ_agg(ik) = champ_sfc(ik)
            end do
         end do
      endif
   end do
   !--------------------------------------------------------------   
   return
end subroutine fillagg
