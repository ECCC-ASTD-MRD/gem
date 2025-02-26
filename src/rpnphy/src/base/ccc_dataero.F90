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
