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
!**s/r fonction schal - chaleur latente selon TT et SWPH
!
      Function schal(tt, ti, swph)
      use tdpack, only: chlc, trpl, folv, fols
      implicit none
!!!#include <arch_specific.hf>
!
      Real tt, ti, schal
!
      Logical swph
!
!author
!       N. Brunet (septembre 2000)
!
!revision
!
!object
!       to calculate the latent heat (J Kg-1) according
!       to temperature and swph
!
!arguments
!       tt - temperature in deg K
!       ti - temperature (K) at which we start calculating
!            latent heat of sublimation
!            if swph=false, ti is n/a
!            ti must be .LE. trpl
!
!note
!      if(not.swph) then schal is equal to folv
!      if(swph) at tt .le. ti, schal is equal to fols
!               at tt .ge. trpl, schal is folv
!               at tt between ti and trpl, schal is a linear
!               interpolation .
!*
!---------------------------------------------------------------
!
      Real lsti, lvt0
!--------------------------------------------------------------------
!
      If(swph)Then
!
         If(tt.Ge.trpl) Then
            schal = folv(tt)
         End If
         If(tt.Le.ti) Then
            schal = fols(tt)
         End If
         If(tt.Lt.trpl .And. tt.Gt.ti .And. ti.Ne.trpl)Then
            lvt0 = chlc
            lsti = fols(ti)
            schal = (lsti-lvt0)*(tt-ti)/(ti-trpl)+lsti
            schal = Max(Min(schal,lsti), lvt0)
         End If
!
      Else
!
         schal = folv(tt)
!
      End If
!
      Return
      End
