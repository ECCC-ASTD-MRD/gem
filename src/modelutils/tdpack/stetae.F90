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
!**s/r fonction stetae - calcule thetae
!
      Function stetae(tt, td, pr)
      use tdpack, only: sttlcl, fopoip, fopoit, foqsa, schal, sesahu3, foefq, cpd, tcdk
      implicit none
!!!#include <arch_specific.hf>
!
      Real stetae
      Real tt, td, pr
!
!author
!     N. Brunet (septembre 2000)
!
!revision
!
!object
!     to compute thetae, the equivalent potential temperature
!
!arguments
!      tt - temperature in deg K
!      td - dew point temperature in deg K
!      pr - pressure (Pa) at the level of tt
!
!note
!      the saturation computations are done with respect to
!      water only
!*
!---------------------------------------------------------------
      Real tl, pl
      Real qsat, chal, teta
      Real es, hu, e, prd
      Real trm3
      Real tetaea
      Real y, ttc, prm
!--------------------------------------------------------------------
!
      If ((tt-td).Gt.0.)Then
!
         tl = sttlcl(td, tt)
         pl = fopoip(tt, tl, pr)
!
      Else
!
         tl = tt
         pl = pr
!
      End If
!
      qsat = foqsa(tl, pl)
!
      chal = schal(tl, -1., .False.)
!
!     il faut calculer "e" et passer "p-e" a fopoit
!
      es=0.
      hu = sesahu3(es,tl,pl,.False.)
      e = foefq(hu,pl)
      prd = pl - e
!
      teta = fopoit(tl,prd,100000.)
!
!     --- ici tetaea est un tetae temporaire
      tetaea = teta*Exp((chal*qsat)/(cpd*tl))
!
!     --- maintenant on calcule le 3e terme ---
!
      es = tt - td
      ttc = tt - tcdk
      prm = pr/100.
!
      y = (-0.1017 + 0.0005*ttc)*es + 0.08*ttc - 7.06 + &
          (1000. - prm)*0.00178
      trm3 = Exp(Exp(y))
!
      stetae = tetaea * trm3
!
      End Function stetae
