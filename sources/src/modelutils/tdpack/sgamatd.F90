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
!**s/r fonction sgamatd - pente de TD lors de VV non-sature
!
      Function sgamatd(td, tt, ti, pr, typv, swph)
      use tdpack, only: schal, rgasd, eps1, grav
      implicit none
!!!#include <arch_specific.hf>
!
      Real sgamatd
      Real td, tt, ti, pr
!
      Integer typv
!
      Logical swph
!
!author
!       N. Brunet (septembre 2000)
!
!object
!       to compute the td lapse rate during an adiabatic
!       and unsaturated ascent
!       2 derivatives are offered: dtd/dp and dtd/dz; the choice is
!       controlled via "typv"
!
!       dtd/dp is > 0 and in K/pa
!       dtd/dz is < 0 and in K/m
!
!arguments
!       td - dew point temperature (K)
!       tt - temperature (K)
!       ti - temperature (K) at which we start calculating
!            latent heat of sublimation
!            if swph=false, ti is n/a
!            ti must be .le. trpl
!       pr  - pressure (pa)
!       typv - if = 1; compute dtd/dp
!                 = 2: compute dtd/dz
!       swph - if .true.: phase ice and water are considered
!                 .false.: phase water for all temperatures with
!                          computation of saturation
!*
!----------------------------------------------------------------
      Real latheat
!
!----------------------------------------------------------------
!
!      ---calcule chaleur latente
      latheat = schal(td, ti, swph)
!
      If(typv .Eq. 1) sgamatd = (rgasd*(td**2))/(eps1*latheat*pr)
!
      If(typv .Eq. 2) sgamatd = -grav*(td**2)/(eps1*latheat*tt)
!
      End Function sgamatd
