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
!**S/P PLANCK - PLANCK FUNCTION
!
      subroutine ccc2_planck2 (bf, bs, urbf, bf0, urbf0, dbf, tfull, gt, ib, &
                         il1, il2, ilg, lay, lev)
      implicit none
!!!#include <arch_specific.hf>
!
      integer ilg, lay, lev, ib, il1, il2
      real  bf(ilg,lev), bs(ilg), bf0(ilg), urbf(ilg,lay), urbf0(ilg), &
            dbf(ilg,lay), tfull(ilg,lev), gt(ilg)
!
!Authors
!
!        J. Li, M. Lazare, CCCMA, rt code for gcm4
!        (Ref: J. Li, H. W. Barker, 2005:
!        JAS Vol. 62, no. 2, pp. 286\226309)
!        P. Vaillancourt, D. Talbot, RPN/CMC;
!        adapted for CMC/RPN physics (May 2006)
!
!Revisions
!
! 001   M.Lazare  (Dec 05,2007) -  new version for gcm15g:  - xx now internal work array.
!
! 002   P. vaillancourt, A-M Leduc  (Apr 2010) - Integrating the Li/Lazare 2009 version. (planck2)
!
! 003    P. Vaillancourt (Feb 12) - assume temperature is isothermal above model top
!
!Object
!
!       Calculation of planck function in valid range 120 - 360 k
!
!Arguments
!
! bf     blackbody intensity integrated over each band at each
!        level in units w / m^2 / sr.
! bs     the blackbody intensity at the surface.
! bf0    the blackbody intensity at the toa (assume 210 k).
! tfull  temperature at each level
! gt     temperature at ground
! uu     1 / diffusivity factor
! urbf   uu times the difference of log(bf) for two neighbor levels
!        used for exponential source function (li, 2002 jas p3302)
! urbf0  uu times the difference of log(bf) for toa and level 1
! dbf    difference of bf for two neighbor levels used for linear
!        source function (li, 2002 jas p3302)
!
!*
!---------------------------------------------------------------------
!     0.0040816327 = 1 / 245 (245 the standard temperature for poly.
!     fit)
!----------------------------------------------------------------------
!
      integer i, j, k, km1, kp1
      real rtstand, uu
      real dt, xxt, xx0
      real  xp(6,9),xx(ilg,lay)


      data  uu / 0.60653066 /, rtstand / 0.0040816327 /
      data ((xp(i,j), i = 1, 6), j = 1, 9) / &
        -2.9876423e+00,    1.3660089e+01,   -1.2944461e+01, &
         1.1775748e+01,   -1.9236798e+01,    2.3584435e+01, &
        -1.6414103e+00,    1.1898535e+01,   -1.1262182e+01, &
         1.0236863e+01,   -1.6677772e+01,    2.0423136e+01, &
         6.5215205e-01,    9.2657366e+00,   -8.5872301e+00, &
         7.6765044e+00,   -1.2287254e+01,    1.4990547e+01, &
         1.5442143e+00,    7.2253228e+00,   -6.7811515e+00, &
         6.1572299e+00,   -9.8725011e+00,    1.1997278e+01, &
         1.2777580e+00,    6.1257638e+00,   -5.7906013e+00, &
         5.3296782e+00,   -8.7529282e+00,    1.0741367e+01, &
         2.1005257e+00,    5.2376301e+00,   -4.8915631e+00, &
         4.5030997e+00,   -7.3199981e+00,    8.9204038e+00, &
         2.9091223e+00,    3.9860795e+00,   -3.5829565e+00, &
         3.2692193e+00,   -5.1799711e+00,    6.2157752e+00, &
         2.7856424e+00,    2.8179582e+00,   -2.3780464e+00, &
         2.1432949e+00,   -3.4540206e+00,    4.1814100e+00, &
         2.4623332e+00,    1.8731841e+00,   -1.3659538e+00, &
         1.1484948e+00,   -1.5975564e+00,    1.7791135e+00 /
!
      do 100 i = il1, il2
        dt          =  gt(i) * rtstand - 1.0
        bs(i)       =  exp( xp(1,ib) + &
                                  dt * (xp(2,ib) + dt * (xp(3,ib) + &
                                  dt * (xp(4,ib) + dt * (xp(5,ib) + &
                                  dt *  xp(6,ib) )))) )
!
! The following line extrapolates the temperature above model top for moon layer temperature
!        dt          = (2. * tfull(i,1) - tfull(i,2)) * rtstand - 1.0
! The following line assumes isothermal temperature above model top
        dt          =  tfull(i,1) * rtstand - 1.0

        xxt         =  xp(1,ib) + dt * (xp(2,ib) + dt * (xp(3,ib) + &
                                  dt * (xp(4,ib) + dt * (xp(5,ib) + &
                                  dt *  xp(6,ib) ))))
!
        dt          =  tfull(i,1) * rtstand - 1.0
        xx0         =  xp(1,ib) + dt * (xp(2,ib) + dt * (xp(3,ib) + &
                                  dt * (xp(4,ib) + dt * (xp(5,ib) + &
                                  dt *  xp(6,ib) ))))
        dt          =  tfull(i,2) * rtstand - 1.
        xx(i,1)     =  xp(1,ib) + dt * (xp(2,ib) + dt * (xp(3,ib) + &
                                  dt * (xp(4,ib) + dt * (xp(5,ib) + &
                                  dt *  xp(6,ib) ))))

        bf0(i)      =  exp(xxt)
        urbf0(i)    =  uu * (xx0 - xxt)
        bf(i,1)     =  exp(xx0)
        bf(i,2)     =  exp(xx(i,1))
        urbf(i,1)   =  uu * (xx(i,1) - xx0)
        dbf(i,1)    =  bf(i,2) - bf(i,1)
  100 continue

      do k = 2, lay
        km1 = k - 1
        kp1 = k + 1
        do i = il1, il2
          dt        =  tfull(i,kp1) * rtstand - 1.0
          xx(i,k)   =  xp(1,ib) + dt * (xp(2,ib) + dt * (xp(3,ib) + &
                                  dt * (xp(4,ib) + dt * (xp(5,ib) + &
                                  dt *  xp(6,ib) ))))

          bf(i,kp1) =  exp(xx(i,k))
          urbf(i,k) =  uu * (xx(i,k) - xx(i,km1))
          dbf(i,k)  =  bf(i,kp1) - bf(i,k)
        enddo
      enddo

      return
      end
