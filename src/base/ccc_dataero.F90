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
      block data dataero
!
      implicit none
#include "nbsnbl.cdk"
#include "ccc_aeros.cdk"
!
      integer i,j
!
!     optical parameters for the standard aerosols of the radiation
!     commission
!
!     model 1= continental, 2=maritime, 3=urban, 4=volcanic, 5=stratospheric
!
!     for background aerosol based on rpn old rad, broad band results
!     for solar, the effect for longwave is small and neglected
!
!
      data ((extab(i,j), j = 1, 5), i = 1, nbs)                        / &
            .730719, .912819, .725059, .745405, .682188, &
            .730719, .912819, .725059, .745405, .682188, &
            .730719, .912819, .725059, .745405, .682188, &
            .730719, .912819, .725059, .745405, .682188                /
      data ((omab(i,j), j = 1, 5), i = 1, nbs)                         / &
            .872212, .982545, .623143, .944887, .997975, &
            .872212, .982545, .623143, .944887, .997975, &
            .872212, .982545, .623143, .944887, .997975, &
            .872212, .982545, .623143, .944887, .997975                /
      data ((gab(i,j), j = 1, 5), i = 1, nbs)                          / &
            .647596, .739002, .580845, .662657, .624246, &
            .647596, .739002, .580845, .662657, .624246, &
            .647596, .739002, .580845, .662657, .624246, &
            .647596, .739002, .580845, .662657, .624246                /
      data ((absab(i,j), j = 1, 5), i = 1, nbl)                        / &
            .0, .0, .0, .0, .0, &
            .0, .0, .0, .0, .0, &
            .0, .0, .0, .0, .0, &
            .0, .0, .0, .0, .0, &
            .0, .0, .0, .0, .0, &
            .0, .0, .0, .0, .0, &
            .0, .0, .0, .0, .0, &
            .0, .0, .0, .0, .0, &
            .0, .0, .0, .0, .0                                         /
!
      end block data dataero
