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
!**s/r mhuaes3  -  passage de hu a es
!
      Subroutine mhuaes3(es,hu,tt,ps,swph,ni,nk,n)
      use tdpack
      implicit none
!!!#include <arch_specific.hf>
!
      Integer ni, nk, n
      Real es(ni,nk), hu(ni,nk), tt(ni,nk)
      Real ps(ni,*)
      Logical swph
!
!Author
!          N. Brunet  (Jan91)
!
!Revision
! 001      B. Bilodeau  (August 1991)- Adaptation to UNIX
! 002      B. Bilodeau (January 2001) - Automatic arrays
! 003      B. Bilodeau (September 2003) - IBM conversion
!                       - call to vslog from massvp4 library
!
!Object
!          to calculate the dew point depression from specific
!          humidity, temperature and pressure
!
!Arguments
!
!          - Output -
! es       dew point depressions in degrees Celsius
!
!          - Input -
! hu       specific humidity in kg/kg
! tt       temperature in K
! ps       pressure in Pa
! swph     .true. to consider water and ice phase
!          .false. to consider water phase only
! ni       horizontal dimension
! nk       vertical dimension
! n        number of treated points
!
!Notes
!          if hu <= 0, the value of hu is not changed but the
!          function max(hu,0.0000000001) will prevent the log
!          of a negative number.
!*
!--------------------------------------------------------------------
      Real, Dimension(n,nk) :: cte
      Real td, petit , alpha
      Integer  k, i
!--------------------------------------------------------------------
!
!
      petit = 0.0000000001
      alpha = log(aerk1w/aerk1i)
!note: cte_ice=cte+alpha
!
      Do k=1,nk
      Do i=1,n
          cte(i,k) = (foefq(Max(petit,hu(i,k)),ps(i,k)))/Real(aerk1w)
      Enddo
      Enddo
!
      Call vslog(cte,cte,n*nk)
!
      Do k=1,nk
      Do i=1,n
         td = (Real(aerk3w)*cte(i,k) - Real(aerk2w)*trpl)/ &
              (cte(i,k) - Real(aerk2w))
!
         If(td.Lt.trpl.And.swph) &
            td = (Real(aerk3i)*(cte(i,k)+alpha) - Real(aerk2i)*trpl) &
                  /(cte(i,k) + alpha - Real(aerk2i))
!
         es(i,k) = tt(i,k) - td
      End Do
      End Do
!
      End Subroutine mhuaes3
