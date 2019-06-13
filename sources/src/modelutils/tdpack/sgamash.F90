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
!**s/r fonction sgamash - pente pseudo-adiabat sature (dt/dz)
!
      Function sgamash(tt, pr, swph, ti)
      use tdpack, only: grav, cpd, eps1, rgasd, schal, foqst, foqsa
      implicit none
!!!#include <arch_specific.hf>
!
      Real tt, pr, ti, sgamash
!
      Logical swph
!
!author
!       N. Brunet (septembre 2000)
!
!revision
!
!object
!       to calculate saturated pseudo-adiabatic lapse rate
!       dt/dz - deg K / m
!       --------- dt/dz is < 0 ---------
!
!arguments
!       tt - temperature (K) at which we calculate the lapse rate
!       pr - pressure (pa) at which we calculate the lapse rate
!       swph - .true., to consider water and ice phase
!              .false, to consider water phase only
!       ti - temperature (K) at which we start calculating
!            latent heat of sublimation
!            if swph=false, ti is n/a
!            ti must be .LE. trpl
!*
!---------------------------------------------------------------
      Real latheat, x, z
      Real gammad
!--------------------------------------------------------------------
!
      gammad = grav / cpd
!
!     calcule la chaleur latente
      latheat = schal(tt, ti, swph)
!
      x = latheat / rgasd
      z = eps1 * latheat**2 / (rgasd * cpd)
!
      If(swph)Then
         sgamash = (-Dble(gammad))*(1.d0+Dble(x)*foqst(tt,pr)/tt) / &
                   (1.d0+Dble(z)*foqst(tt,pr)/tt**2)
      Else
         sgamash = (-Dble(gammad))*(1.d0+Dble(x)*foqsa(tt,pr)/tt) / &
                   (1.d0+Dble(z)*foqsa(tt,pr)/tt**2)
      End If
!
      End Function sgamash
