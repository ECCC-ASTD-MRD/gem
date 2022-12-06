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

!/@*
subroutine intozon2(jour, mois)
   implicit none
!!!#include <arch_specific.hf>
   integer, intent(in) :: jour,mois
   !@Author B. Dugas (Winter 2001) - From litozon2
   !@Object
   !          produces a new ozone field valid at jour/mois from
   !          the previously read climatology table gozon12.
   !          Allows bit-reproducibility for time integrations
   !          running a multiple of 1 day per "clone".
   !@Arguments
   !          - Input -
   ! mois     month of ozone record
   ! jour     day of the month
   !@Notes
   !          Monthly climatological values are supposed to be valid
   !          on the 15th. No interpolation is needed (done) when
   !          the input jour is 15
   !*@/

#include <rmn/msg.h>
#include "radiation.cdk"
#include "ozopnt.cdk"

   integer, parameter :: ANNEE(12) = &
        (/31,28,31,30,31,30,31,31,30,31,30,31/)

   real :: total, ecoule
   integer :: basem, destm, j, nlp
   character(len=256) ::tmp_S
   !-----------------------------------------------------------------
   nlp = nlacl*npcl

   !# doit-on interpoler ?
   if (jour < 15) then
      if (mois == 1) then
         destm = 1
         basem = 12
      else
         destm = mois
         basem = destm-1
      endif
      ecoule = jour+ANNEE(basem)-15
   else if (jour > 15) then
      if (mois == 12) then
         basem = 12
         destm = 1
      else
         basem = mois
         destm = basem+1
      endif
      ecoule = jour-15
   else
      basem = mois
      destm = basem
   endif

   if (destm /= basem) then
      total = 1./real(ANNEE(basem))
      !# interpoler pour le jour courant.
      do j=1,nlp
         goz(j) = gozon12(j,basem) + &
              (gozon12(j,destm)-gozon12(j,basem))*total*ecoule
      enddo
      write(tmp_S, '(a,1x,i0,1x,a,1x,i0,1x,a,1x,f0.3,1x,a,1x,i0,1x,a,1x,i0)') 'day=', jour, ', month=', mois, ' [',total*ecoule*100.,'%]',basem, ', ', destm
      call msg(MSG_INFO,'intozon: ozone interpolated to '//tmp_S)
   else
      do j=1,nlp
         goz(j) = gozon12(j,destm)
      enddo
      write(tmp_S, '(a,1x,i0,1x,a,1x,i0)') 'day=', jour, ', month=', mois
      call msg(MSG_INFO,'intozon: ozone interpolated to '//tmp_S)
   endif
   !-----------------------------------------------------------------
   return
end subroutine intozon2
