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
!**S/P PREINTP - DETERMINES THE PRESSURE INTERPOLATION POINTS
!
      subroutine ccc1_preintp (inpt, inptm, dip, dip0, pp, il1, il2, ilg,lay)
      implicit none
!!!#include <arch_specific.hf>
!
      integer ilg, lay, il1, il2, jends, k, i, j, inpdif, m, n
      real p0(ilg)  !#, pm
      real dip(ilg,lay), dip0(ilg)
      real pp(ilg,lay), standp(28)
      integer inpt(ilg,lay), inptm(ilg,lay)
!
!Authors
!        J. Li, M. Lazare, CCCMA, rt code for gcm4
!        (Ref: J. Li, H. W. Barker, 2005:
!        JAS Vol. 62, no. 2, pp. 286\226309)
!        P. Vaillancourt, D. Talbot, RPN/CMC;
!        adapted for CMC/RPN physics (May 2006)
!
!Revisions
! 001
!Object
!        This subroutine determines the pressure interpolation points
!
!Arguments
!
! inpt   number of the level for the standard input data pressures
!        (for 28 interpolation levels)
! inptm  number of the level for the standard input data pressures
!        (for 18 interpolation levels below 1 mb)
! p      pressure at middle of each layer
! dip    interpolation factor for pressure between two
!        neighboring standard input data pressure levels
! dip0   interpolation factor for pressure above model top level
!----------------------------------------------------------------------
!*
      data standp / 5.0000e-04, 1.4604e-03, 2.9621e-03, 6.0080e-03, &
                    1.2186e-02, 2.4717e-02, 5.0134e-02, 1.0169e-01, &
                    2.0625e-01, 4.1834e-01, &
                    1.2180, 1.8075, 2.6824, 3.9806, 5.9072, 8.7662, &
                    13.0091, 19.3054, 28.6491, 42.5151, 63.0922, &
                    93.6284, 138.9440, 206.1920, 305.9876, 454.0837, &
                    673.8573, 1000.0000 /
!
      jends = 27
      do 500 k = 1, lay
        inpdif =  0
        do 200 i = il1, il2
          inpt(i,k)   =  0
          do 100 j = 1, jends
            if (pp(i,k) .gt. standp(j))                               then
              inpt(i,k) =  inpt(i,k) + 1
            endif
  100   continue
!
!----------------------------------------------------------------------
!     calculate arrays dip and dit required later for gasopt routines.
!     also, set values of inpt for a given level to be negative if all
!     longitude values are the same. this is also used in the gasopt
!     routines to improve performance by eliminating the unnecessary
!     indirect-addressing if inpt is negative for a given level.
!     note that for inpt=0, it is assumed that levels are more or
!     less horizontal in pressure, so scaling by -1 still preserves
!     the value of zero and no indirect-addressing is done in the
!     gasopt routines.
!----------------------------------------------------------------------
!
          if(inpt(i,k) .ne. inpt(1,k) )  inpdif = 1
          m  =  inpt(i,k)
          n  =  m + 1
          if (m .gt. 0)                                             then
            dip(i,k)  = (pp(i,k) - standp(m)) / (standp(n) - standp(m))
          else
            dip(i,k)  =  pp(i,k) / standp(1)
          endif
  200   continue
!
!       when all values along i are the same
!       we add 1000
!
        if(inpdif .eq. 0)                                           then
          do 250 i = il1, il2
            inpt(i,k) =  inpt(i,k) + 1000
  250     continue
        endif
!
        do 300 i = il1, il2
          inptm(i,k)  =  inpt(i,k) - 10
  300   continue
!
  500 continue
!
!----------------------------------------------------------------------
!     interpolation factor for lattenu and sattenu (attenuation above
!     model top
!     note : remove commented lines below if top is less than .0005
!----------------------------------------------------------------------
!       pm =  pp(1,1)
!       do 700 i = il1, il2
!         pm          =  min (pm, pp(i,1))
! 700   continue
!
!       if (pm .le. 0.0005)                                         then
!         do 800 i = il1, il2
!           dip0(i)   =  0.0
! 800     continue
!       else
          do i = il1, il2
            p0(i)   = pp(i,1)*pp(i,1) / pp(i,2)
            dip0(i) = sqrt(p0(i) * pp(i,1))
            dip0(i) = (dip0(i) - pp(i,1)) / (p0(i) - pp(i,1))
          enddo
!       endif
!
      return
      end
