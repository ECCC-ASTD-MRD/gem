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
!**
!*
!
      REAL FUNCTION INTERPVEG( JULIEN, TABLE )
!
      implicit none
#include <arch_specific.hf>
!
!
      REAL JULIEN
      REAL TABLE(13)
!
!
!Author
!           S. Belair (February 1999)
!
!
!Revision
! 001
!
!
!Object
!           Interpolate the characteristics of the vegetation
!           (vegetation fraction and leaf area index mainly)
!           for the day of the year, using a table with monthly
!           values
!
!
!Arguments
!
!           - Input -
! JULIEN    Julian day
! TABLE     Table of monthly values (the first month is repeated
!           at the end, so the dimension is 13 instead of 12)
!
!
!*
!**
!*
!
      INTEGER MONTH
      REAL MONTHL, HALFM, YEARL
      REAL DAY, DAYOFMONTH
!
!
      yearl  = 366.
      monthl = yearl / 12.
      halfm  = yearl / 12. / 2.
!
!
      day        = julien - halfm
      if (day.lt.0.) day = day + yearl
!
      month      = INT( day / monthl ) + 1
      dayofmonth = day - float(month-1)*monthl
!
!
      interpveg  = table(month) &
                 + dayofmonth / monthl * (table(month+1)-table(month))
!
!
      RETURN
      END

