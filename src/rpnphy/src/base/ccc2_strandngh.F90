!-------------------------------------- LICENCE BEGIN
!------------------------------------
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
      subroutine ccc2_strandngh4 (tran, gwgh, atten, taua, tauoma, &
                            taucs, tauomc, cldfrac, rmu, dp, &
                            o3, qq, co2, ch4, o2, ib, ig, inpt, &
                            dip, dt, lev1, gh, cut, &
                            il1, il2, ilg, lay, lev)
!
      implicit none
!!!#include <arch_specific.hf>
!
      integer ilg, lay, lev, ib, ig, lev1, il1, il2
      real gwgh, cut
      real tran(ilg,2,lev), atten(ilg), taua(ilg,lay), tauoma(ilg,lay), &
           taucs(ilg,lay), tauomc(ilg,lay), cldfrac(ilg,lay), rmu(ilg), &
           dp(ilg,lay), o3(ilg,lay),  o2(ilg,lay), co2(ilg,lay), ch4(ilg,lay), &
           qq(ilg,lay), dip(ilg,lay), dt(ilg,lay)
       integer inpt(ilg,lay)
       logical gh

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
! 001  M. Lazare (May 05,2006) -  previous version strandngh2 for gcm15e/f:
!                              - pass integer variables "init" and
!                                 "nit" instead of actual integer
!                                 values, to "tline_" routines.
! 002  M.Lazare, L.Solheim, J. Li (Apr 18,2008) - previous version strandngh3 for gcm15g:
!                             - cosmetic change to add threadprivate
!                            for common block "trace", in support
!                            of "radforce" model option.
!                            - calls tline{1,2}y instead of tline{1,2}x.
!                            - updating o3, adding ch4 and using
!                              kurucz solar function.
! 003  M.Lazarre,K.Winger,P.Vaillancourt   (Apr 08) - use integer variables instead of actual integers
! 004  J. Li (Feb 09,2009)      - new version for gcm15h:
!                               - 3d ghg implemented, thus no need
!                                 for "trace" common block or
!                                 temporary work arrays to hold
!                                 mixing ratios of ghg depending on
!                                 a passed, specified option.
! 005  P. vaillancourt, A-M. Leduc (Apr 2010) - Integrating the Li/Lazare 2009 version.(strandngh4)
!                                add arguments: co2,ch4,o2
!                                remove ng, ng2 , gw1 ,cs1o3, cs1o21, and remove arguments , s, taug



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
!
!             - Output -
! tran     downward flux
! gwgh                               A DEFINIR   AML
!
!            - Input -
! atten    attenuation factor for downward flux from toa to the
!          model top level
! taua     arerosol optical depth              definition  TO CHECK  AML
! tauoma   arerosol optical depth              definition  TO CHECK  AML
! taucs    cloud optical depth
! tauomc                                       definition  TO CHECK  AML
! cldfrac  cloud fraction
! rmu      cos of solar zenith angle
! dp       air mass path for a model layer (explained in raddriv)
! o3       o3 mass mixing ratio
! qq       water vapor mass mixing ratio
! co2,     ch4, o2 are mass mixing ratio
! ib                                 A DEFINIR   AML
! ig                                 A DEFINIR   AML
! inpt     number of the level for the standard input data pressures
! dip      interpolation factor for pressure between two
!          neighboring standard input data pressure levels
! dt       layer temperature - 250 k
! lev1     a level close to 1 mb determined in cldifm  TO CHECK   AML
! gh       logical key      A DEFINIR   AML
! cut                              A DEFINIR   AML
! il1      1
! il2      horizontal dimension
! ilg      horizontal dimension
! lay      number of model levels

!
!Implicites
!
#include "ccc_tracegases.cdk"
#include "ccc2_bandsh.cdk"
!
!*
      integer i, k, kp1, im, init
      real absc,dto3
      real tau,dtr1,taug(ilg,lay)

      do 10 i = il1, il2
        tran(i,1,1)           =  atten(i)
        tran(i,2,1)           =  atten(i)
   10 continue
!
      if (ib .eq. 1)                                                then
!
!----------------------------------------------------------------------
!     band1 for uvc (35700 - 50000 cm^-1), nongray gaseous absorption
!     of o2  and o3. solar energy  7.5917 w m^-2
!----------------------------------------------------------------------
!
        if (ig .eq. 3)                                              then
          do 105 k = 1, lay
            kp1 = k + 1
            do i = il1, il2
              dto3            =  dt(i,k) + 23.13
              tau             = ((cs1o3gh(1,ig) + &
                                 dto3 * (cs1o3gh(2,ig) + &
                                 dto3 * cs1o3gh(3,ig))) * o3(i,k) + &
                                 cs1o2gh3 * o2(i,k)) * dp(i,k) + &
                                 taua(i,k)
               dtr1            =  exp( - (tau    - tauoma(i,k)) / rmu(i) )
!              dto3            =  dt(i,k) + 23.13
!              tau             = ((cs1o3gh(1,ig) +
!     1                           dto3 * (cs1o3gh(2,ig) +
!     2                           dto3 * cs1o3gh(3,ig))) * o3(i,k) +
!     3                           cs1o2gh3 * o2(i,k)) * dp(i,k) +
!     4                           taua(i,k)
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
              dto3            =  dt(i,k) + 23.13
              tau             = (cs1o3gh(1,ig) + dto3 * (cs1o3gh(2,ig) + &
                                 dto3 * cs1o3gh(3,ig))) * o3(i,k) * &
                                 dp(i,k) + taua(i,k)
              dtr1            =   exp( - (tau    - tauoma(i,k)) / rmu(i) )
!             tau             = (cs1o3gh(1,ig) + dto3 * (cs1o3gh(2,ig) +
!     1                           dto3 * cs1o3gh(3,ig))) * o3(i,k) *
!     2                           dp(i,k) + taua(i,k)
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
!     and o3. solar energy 8.9036 w m^-2
!----------------------------------------------------------------------
!
        if (ig .eq. 1)                                              then
          init=2
          call ccc2_tline1z (taug, cs2h2ogh(1,1), qq, dp, dip, &
                       dt, inpt, lev1, gh, ntl, init, &
                       il1, il2, ilg, lay)
        else
          im =  ig - 1
          init=2
          call ccc2_tline1z (taug, cs2o2gh(1,1,im), o2, dp, dip, &
                       dt, inpt, lev1, gh, ntl, init, &
                       il1, il2, ilg, lay)
        endif
!
        do 205 k = 1, lay
          kp1 = k + 1
          do i = il1, il2
            tau               =  taug(i,k) + taua(i,k)
            dtr1              =  exp( - (tau    - tauoma(i,k)) / rmu(i) )
!           tau               =  taug(i,k) + taua(i,k)
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
!     co2. solar energy 7.4453 w m^-2
!----------------------------------------------------------------------
!
        if (ig .le. 2)                                              then
          call ccc2_tline2z (taug, cs3h2ogh(1,1,ig), cs3co2gh(1,1,ig), qq, &
                       co2, dp, dip, dt, inpt, &
                       lev1, gh, ntl, il1, il2, ilg, lay)

        else
          init=2
          call ccc2_tline1z (taug, cs3co2gh(1,1,ig), co2, dp, dip, &
                       dt, inpt, lev1, gh, ntl, init, &
                       il1, il2, ilg, lay)
        endif
!
        do 305 k = 1, lay
          kp1 = k + 1
          do i = il1, il2
            tau               =  taug(i,k) + taua(i,k)
            dtr1              =  exp( - (tau    - tauoma(i,k)) / rmu(i))
!           tau               =  taug(i,k) + taua(i,k)
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
!     and co2. solar energy 7.0384 w m^-2
!----------------------------------------------------------------------
!
        if (ig .le. 3)                                              then
          call ccc2_tline2z (taug, cs4h2ogh(1,1,ig), cs4co2gh(1,1,ig), qq, &
                       co2, dp, dip, dt, inpt, &
                       lev1, gh, ntl, il1, il2, ilg, lay)
        else if (ig .eq. 6 .or. ig .eq. 8)           then
          if (ig .eq. 6)  im = 5
          if (ig .eq. 8)  im = 6
          init=2
          call ccc2_tline1z (taug, cs4h2ogh(1,1,im), qq, dp, dip, &
                       dt, inpt, lev1, gh, ntl, init, &
                       il1, il2, ilg, lay)
        else if (ig .eq. 5 .or. ig .eq. 7 .or. ig .eq. 9)           then
          if (ig .eq. 5)  im = 4
          if (ig .eq. 7)  im = 5
          if (ig .eq. 9)  im = 6
          init=2
          call ccc2_tline1z (taug, cs4co2gh(1,1,im), co2, dp, dip, &
                       dt, inpt, lev1, gh, ntl, init, &
                       il1, il2, ilg, lay)
        else
         im = 4
          call ccc2_tline2z (taug, cs4h2ogh(1,1,im), cs4ch4gh, qq, ch4, &
                        dp, dip, dt, inpt, lev1, gh, ntl, &
                        il1, il2, ilg, lay)


        endif
!
        do 405 k = 1, lay
          kp1 = k + 1
          do i = il1, il2
            tau               =  taug(i,k) + taua(i,k)
            dtr1              =  exp( - (tau    - tauoma(i,k)) / rmu(i) )
!           tau               =  taug(i,k) + taua(i,k)
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
