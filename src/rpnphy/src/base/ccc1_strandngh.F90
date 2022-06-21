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
!**S/P STRANDNGH - CALCULATION OF THE DOWNWARD SOLAR FLUX
!
      subroutine ccc1_strandngh (tran, gwgh, atten, taua, tauoma, &
                            taucs, tauomc, cldfrac, rmu, dp, &
                            o3, qq, ib, ig, inpt, &
                            dip, dt, lev1, gh, cut, &
                            il1, il2, ilg, lay, lev, &
                            taug, s)
!
      implicit none
!!!#include <arch_specific.hf>
!
      integer ilg, lay, lev, ib, ig, lev1, il1, il2, i, k, kp1, ng
      integer im, ng2
      real gwgh, cut, absc
      real tran(ilg,2,lev), atten(ilg), taua(ilg,lay), tauoma(ilg,lay), &
           taucs(ilg,lay), tauomc(ilg,lay), cldfrac(ilg,lay), rmu(ilg), &
           dp(ilg,lay), o3(ilg,lay), qq(ilg,lay), dip(ilg,lay), &
           dt(ilg,lay), taug(ilg,lay), s(ilg,lay)
      integer inpt(ilg,lay)
      logical gh
      real tau,dtr1
      integer init

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
! 001    M.Lazarre,K.Winger,P.Vaillancourt   (Apr 08) - use integer variables instead of actual integers
!
!Object
!
!        Calculation of the downward solar flux under the condition that
!        the extinction coefficient of gas is very large, the scattering
!        effects can be neglected. the cloud optical depth is much smaller
!        than the gaseous optical depth, the cloud effect is very small
!        and be treated simply
!
! hrs      solar heating rate (k / sec)
! tran     downward flux
! atten    attenuation factor for downward flux from toa to the
!          model top level
! taucs    cloud optical depth
! cldfrac  cloud fraction
! rmu      cos of solar zenith angle
! dp       air mass path for a model layer (explained in raddriv)
! o3       o3 mass mixing ratio
! qq       water vapor mass mixing ratio
! inpt     number of the level for the standard input data pressures
! dip      interpolation factor for pressure between two
!          neighboring standard input data pressure levels
! dt       layer temperature - 250 k
!
!Implicites
!
#include "ccc_tracegases.cdk"
#include "ccc1_bandsh.cdk"
!
!*
      do 10 i = il1, il2
        tran(i,1,1)           =  atten(i)
        tran(i,2,1)           =  atten(i)
   10 continue
!
      if (ib .eq. 1)                                                then
!
!----------------------------------------------------------------------
!     band1 for uvc (35700 - 50000 cm^-1), nongray gaseous absorption
!     of o2  and o3. solar energy  6.9015 w m^-2
!----------------------------------------------------------------------
!
        if (ig .eq. 3)                                              then
          do 105 k = 1, lay
            kp1 = k + 1
            do i = il1, il2
              tau             = (cs1o3gh(ig) * o3(i,k) + cs1o2gh3 * &
                                 rmo2) * dp(i,k) + taua(i,k)
              dtr1            = exp(- (tau    - tauoma(i,k)) / rmu(i))
              tran(i,1,kp1)   =  tran(i,1,k) * dtr1   
!
              if (cldfrac(i,k) .lt. cut)                            then
                tran(i,2,kp1) =  tran(i,2,k) * dtr1   
              else
                absc          = (1.0-cldfrac(i,k))*dtr1   +cldfrac(i,k)* &
                                 exp( -(tau   +taucs(i,k)-tauomc(i,k)) &
                                 / rmu(i))
                tran(i,2,kp1) =  tran(i,2,k) * absc
              endif
           enddo
  105     continue
        else
          do 115 k = 1, lay
            kp1 = k + 1
            do i = il1, il2
              tau             =  cs1o3gh(ig) * o3(i,k) * dp(i,k) + &
                                 taua(i,k)
              dtr1            =  exp(- (tau    - tauoma(i,k)) / rmu(i))
              tran(i,1,kp1)   =  tran(i,1,k) * dtr1   
!
              if (cldfrac(i,k) .lt. cut)                            then
                tran(i,2,kp1) =  tran(i,2,k) * dtr1   
              else
                absc          =(1.0-cldfrac(i,k))*dtr1   +cldfrac(i,k) * &
                                 exp( - (tau   +taucs(i,k)-tauomc(i,k)) &
                                 / rmu(i))
                tran(i,2,kp1) =  tran(i,2,k) * absc
              endif
           enddo
  115     continue
        endif
      gwgh =  gws1gh(ig)
!
      else if (ib .eq. 2)                                           then
!
!----------------------------------------------------------------------
!     band (8400 - 14500 cm^-1), nongray gaseous absorption of o2
!     and o3. solar energy 8.72450 w m^-2
!----------------------------------------------------------------------
!
        if (ig .eq. 1)                                              then
          ng =  1
          init=2
          call ccc1_tline1 (taug, cs2h2ogh(1,1), qq, ng, dp, dip, &
                       dt, inpt, lev1, gh, ntl, init, &
                       il1, il2, ilg, lay, s)
        else
          im =  ig - 1
          ng =  6
          init=2
          call ccc1_tline1 (taug, cs2o2gh(1,1,im), qq, ng, dp, dip, &
                       dt, inpt, lev1, gh, ntl, init, &
                       il1, il2, ilg, lay, s)
        endif
!
        do 205 k = 1, lay
          kp1 = k + 1
          do i = il1, il2
            tau               =  taug(i,k) + taua(i,k)
            dtr1              =  exp(- (tau    - tauoma(i,k)) / rmu(i))
            tran(i,1,kp1)     =  tran(i,1,k) * dtr1   
!
            if (cldfrac(i,k) .lt. cut)                              then
              tran(i,2,kp1)   =  tran(i,2,k) * dtr1   
            else
              absc            = (1.0-cldfrac(i,k))*dtr1   +cldfrac(i,k)* &
                                 exp( -(tau   +taucs(i,k)-tauomc(i,k)) &
                                 / rmu(i))
              tran(i,2,kp1)   =  tran(i,2,k) * absc
            endif
         enddo
  205   continue
!
        gwgh =  gws2gh(ig)
!
      else if (ib .eq. 3)                                           then
!
!----------------------------------------------------------------------
!     band (4200 - 8400 cm^-1), nongray gaseous absorption of h2o and
!     co2. solar energy 4.0330 w m^-2
!----------------------------------------------------------------------
!
        if (ig .le. 2)                                              then
          ng2 =  3
          call ccc1_tline2 (taug, cs3h2ogh(1,1,ig), cs3co2gh(1,1,ig), qq, o3, &
                       ng2, dp, dip, dt, inpt, &
                       lev1, gh, ntl, il1, il2, ilg, lay, s)

        else
          ng =  3
          init=2
          call ccc1_tline1 (taug, cs3co2gh(1,1,ig), qq, ng, dp, dip, &
                       dt, inpt, lev1, gh, ntl, init, &
                       il1, il2, ilg, lay, s)
        endif
!
        do 305 k = 1, lay
          kp1 = k + 1
          do i = il1, il2
            tau               =  taug(i,k) + taua(i,k)
            dtr1              =  exp(- (tau    - tauoma(i,k)) / rmu(i))
            tran(i,1,kp1)     =  tran(i,1,k) * dtr1
!
            if (cldfrac(i,k) .lt. cut)                              then
              tran(i,2,kp1)   =  tran(i,2,k) * dtr1
            else
              absc            = (1.0-cldfrac(i,k))*dtr1   +cldfrac(i,k)* &
                                 exp( - (tau   +taucs(i,k)-tauomc(i,k)) &
                                 / rmu(i))
              tran(i,2,kp1)   =  tran(i,2,k) * absc
            endif
         enddo
  305   continue
!
        gwgh =  gws3gh(ig)
!
      else if (ib .eq. 4)                                           then
!
!----------------------------------------------------------------------
!     band (2500 - 4200 cm^-1), nongray gaseous absorption of h2o
!     and co2. solar energy 9.6020 w m^-2
!----------------------------------------------------------------------
!
        if (ig .le. 3)                                              then
          ng2 =  3
          call ccc1_tline2 (taug, cs4h2ogh(1,1,ig), cs4co2gh(1,1,ig), qq, o3, &
                       ng2, dp, dip, dt, inpt, &
                       lev1, gh, ntl, il1, il2, ilg, lay, s)
        else if (ig .eq. 4 .or. ig .eq. 6 .or. ig .eq. 8)           then
          ng =  1
          if (ig .eq. 4)  im = 4
          if (ig .eq. 6)  im = 5
          if (ig .eq. 8)  im = 6
          init=2
          call ccc1_tline1 (taug, cs4h2ogh(1,1,im), qq, ng, dp, dip, &
                       dt, inpt, lev1, gh, ntl, init, &
                       il1, il2, ilg, lay, s)
        else
          ng =  3
          if (ig .eq. 5)  im = 4
          if (ig .eq. 7)  im = 5
          if (ig .eq. 9)  im = 6
          init=2
          call ccc1_tline1 (taug, cs4co2gh(1,1,im), qq, ng, dp, dip, &
                       dt, inpt, lev1, gh, ntl, init, &
                       il1, il2, ilg, lay, s)
        endif
!
        do 405 k = 1, lay
          kp1 = k + 1
          do i = il1, il2
            tau               =  taug(i,k) + taua(i,k)
            dtr1              =  exp(- (tau    - tauoma(i,k)) / rmu(i))
            tran(i,1,kp1)     =  tran(i,1,k) * dtr1
!
            if (cldfrac(i,k) .lt. cut)                              then
              tran(i,2,kp1)   =  tran(i,2,k) * dtr1
            else
              absc            = (1.0-cldfrac(i,k))*dtr1   +cldfrac(i,k)* &
                                 exp(-(tau   +taucs(i,k)-tauomc(i,k)) &
                                 / rmu(i))
              tran(i,2,kp1)   =  tran(i,2,k) * absc
            endif
         enddo
  405   continue
!
        gwgh =  gws4gh(ig)
!
      endif
!
      return
      end
