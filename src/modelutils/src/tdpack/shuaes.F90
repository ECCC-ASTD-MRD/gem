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
!**fonction shuaes3  -  passage de hu a es
!
      Function shuaes3(hu,tt,ps,swph)
      use tdpack, only: aerk1w, aerk1i, aerk2w, aerk2i, aerk3w, aerk3i, foefq, trpl
      implicit none
!!!#include <arch_specific.hf>
      Real shuaes3, hu, tt, ps
      Logical swph
!
!Author
!          N. Brunet  (Jan91)
!
!Object
!          to return dew point depression (Celsius) calculated from
!          specific humidity, temperature and pressure
!
!Arguments
!
!          - Input -
! hu       specific humidity in kg/kg
! tt       temperature in K
! ps       pressure in Pa
! swph     .true. to consider water and ice phase
!          .false. to consider water phase only
!
!Notes
!     If hu <=0, we don't change the value of hu, but we take
!     the maximum value of: max(hu,0.0000000001).
!     to avoid the occurence of taking the log of a negative
!     number
!*
!--------------------------------------------------------------------
      Real e, cte, td, petit, alpha
!--------------------------------------------------------------------
!
      petit = 0.0000000001
      alpha = log(aerk1w/aerk1i)
!note: cte_ice=cte+alpha
!
      e = foefq(Max(petit,hu),ps)
      cte = alog(e/Real(aerk1w))
      td = (Real(aerk3w)*cte - Real(aerk2w)*trpl)/(cte - Real(aerk2w))
      If(td.Lt.trpl.And.swph)Then
         td = (Real(aerk3i)*(cte + alpha) - Real(aerk2i)*trpl)/(cte + alpha - Real(aerk2i))
      End If
!
      shuaes3 = tt-td
!
      End Function shuaes3
