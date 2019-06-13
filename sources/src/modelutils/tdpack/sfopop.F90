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
!**fonction sfopop - resout l'eqn de Poisson pour la pression
!
      Function sfopop(t0, tf2, p0)
      use tdpack, only: fopoip
      implicit none
!!!#include <arch_specific.hf>
!
      Real sfopop, t0, tf2, p0
!
!author
!      N. Brunet (sept 2000)
!
!object
!      solve Poisson equation for pressure
!
!arguments
!
!      t0  - temperature (K) at a given pressure p0
!      p0  - pressure (Pa)
!      tf2 - temperature (K) at the target level
!
!      result is in Pa
!*
!--------------------------------------------------------------------
!
      sfopop = fopoip(t0, tf2, p0)
!
      End Function sfopop
