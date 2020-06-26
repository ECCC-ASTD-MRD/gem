!-------------------------------------- LICENCE BEGIN -------------------------
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
!-------------------------------------- LICENCE END ---------------------------

subroutine ccc2_lwtran4 (fu, fd, slwf, tauci, omci, &
     gci, fl, taual, taug, bf, &
     bs, urbf, dbf, em0, cldfrac, &
     cldm, anu, nct, ncd, ncu, &
     ncum, ncdm, lev1, cut, maxc, &
     il1, il2, ilg, lay, lev)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
   integer ilg, lay, lev, lev1, il1, il2, l1, l2, maxc, i
   integer k, km1, km2, kp1, kp2, kx, kxm, kk, km3, nanu
   real cut, ubeta, epsd, epsu, zeta, taul2, ssalb
   real sf, ww, cow, tau2, sanu, xx, yy, dtr2, embk
   real ru, wt, sx, sy, p1, p2, zz, x2, y2, x3, y3, wgrcow
   real taudtr, bkins, anutau, dtrgw, fwins, fmbk
   real fu(ilg,2,lev), fd(ilg,2,lev)
   real slwf(ilg), tauci(ilg,lay), omci(ilg,lay), gci(ilg,lay), &
        fl(ilg,lay), taual(ilg,lay), taug(ilg,lay), bf(ilg,lev), &
        bs(ilg),urbf(ilg,lay),dbf(ilg,lay),em0(ilg),cldfrac(ilg,lay), &
        cldm(ilg,lay), anu(ilg,lay), taul1(ilg,lay), rtaul1(ilg,lay)
   real taucsg(ilg,lay), term1(ilg), emisw(ilg,lay), scatbk(ilg,lay), &
        scatfw(ilg,lay), scatsm(ilg,4,lay), taum(ilg,4,lay), &
        xd(ilg,4,lay), xu(ilg,4,lay), dtr(ilg,4,lay), fx(ilg,4,lev), &
        fy(ilg,4,lev), fw(ilg,4,lev), s1(ilg), c2(ilg)
   real embs,abse0
   real(REAL64) :: dtr_vs(ilg,lay)
   integer nct(ilg), ncd(ilg,lay), ncu(ilg,lay), ncum(lay), ncdm(lay)

   data  ru / 1.6487213 /

   !@Authors
   !        J. Li, M. Lazare, CCCMA, rt code for gcm4
   !        (Ref: J. Li, H. W. Barker, 2005:
   !        JAS Vol. 62, no. 2, pp. 286\226309)
   !        P. Vaillancourt, D. Talbot, RPN/CMC;
   !        adapted for CMC/RPN physics (May 2006)
   !@Revisions
   ! 001     P.Vaillancourt, K.Winger, M.Lazarre (Sep 2006) :
   !         replace int par nint, remplace 1-cld-eps par max(1-cld,eps)
   ! 002     M.LAZARE (Dec 2007) . previous version lwtran3 for gcm15g:
   !                              - s,emisw,scatbk,scatfw,term1,s1,t
   !                                now internal work arrays.
   !                              - support added for emissivity<1.
   ! 003     J. Cole (Feb 2009) new version for gcm15h:
   !                            - correct bug in specification of abse0.
   ! 004     P. Vaillancourt (May 2010) adapted gcm15h code for CMC/RPN physics
   !@Object LONGWAVE RADIATIVE TRANSFER
   !         Calculation of longwave radiative transfer using absorption
   !         approximation. The finite cloud effect is properly considered
   !         with random and full overlap assumption. Cloud subgrid
   !         variability is included (based on Li, 2002 JAS p3302; Li and
   !         Barker JAS p3321).
   !@Arguments
   !          - Output -
   ! fu       upward infrared flux
   ! fd       downward infrared flux
   ! taucsg   total scattering
   ! term1    surface albedo
   ! emisw
   ! scatbk   backward scattering
   ! scatfw   forward scattering
   ! scatsm   internal scattering
   ! taum     taum(1) a factor related to tau in zeta factor for
   !          linear source term; taum(2) the cumulated taum(1) for
   !          subgrid variability calculation
   ! xd       the emission part in the downward flux transmission
   ! xu       the emission part in the upward flux transmission
   !          (Li, 2002 JAS p3302)
   ! dtr      direct transmission
   ! fx       the same as fy but for the downward flux
   ! fy       upward flux for pure clear portion (1) and pure cloud
   !          portion (2)
   ! fw       a term for transfer within cldm
   ! s1       rmug
   ! c2
   !          - Intput -
   ! slwf     input solar flux at model top level for each band
   ! tauci    cloud optical depth for the infrared
   ! omci     cloud single scattering albedo times optical depth
   ! gci      cloud asymmetry factor times omci
   ! fl       square of cloud asymmetry factor
   ! taual    aerosol optical depth for the infrared
   ! taug     gaseous optical depth for the infrared
   ! bf       blackbody intensity integrated over each band at each
   !          level in units w / m^2 / sr. therefor a pi factor needed
   !          for flux
   ! bs       the blackbody intensity at the surface.
   ! urbf     u times the difference of log(bf) for two neighbor
   !          levels used for exponential source function
   ! dbf      difference of bf for two neighbor levels used for
   !          linear source function
   ! em0      surface emission
   ! cldfrac  cloud fraction
   ! cldm     maximum portion in each cloud block, in which the exact
   !          solution for subgrid variability is applied li and
   !          Barker JAS p3321).
   ! anu      nu factor for cloud subgrid variability
   ! nct      the highest cloud top level for the longitude and
   !          latitude loop (ilg)
   ! ncd      layer inside a cloud block accounted from cloud top
   ! ncu      layer inside a cloud block accounted from cloud bottom
   ! ncum     maximum loop number cloud vertical correlation accounted
   !          from lower level to higher level
   ! ncdm     maximum loop number cloud vertical correlation accounted
   !          from higher level to lower level
   !----------------------------------------------------------------------

   !----------------------------------------------------------------------
   !     initialization for first layer. calculate the downward flux in
   !     the second layer
   !     combine the optical properties for the infrared,
   !     1, aerosol + gas; 2, cloud + aerosol + gas.
   !     fd (fu) is down (upward) flux, fx (fy) is the incident flux
   !     above (below) the considered layer.
   !     gaussian integration and diffusivity factor, ru (li jas 2000)
   !     above maxc, exponential source function is used
   !     below maxc, linear source function is used
   !----------------------------------------------------------------------

   l1 =  lev1
   l2 =  lev1 + 1
   do i = il1, il2
      fd(i,1,lev1)              =  slwf(i)
      fd(i,2,lev1)              =  slwf(i)
      fx(i,1,lev1)              =  fd(i,1,lev1)
      fx(i,2,lev1)              =  fd(i,2,lev1)
   enddo

   do k = l2, lev
      km1 = k - 1
      do i = il1, il2
         taul1(i,km1)            =  taual(i,km1) + taug(i,km1)
         rtaul1(i,km1)           =  taul1(i,km1) * ru
         dtr_vs(i,k-l2+1)        =  exp(dble(-rtaul1(i,km1)))
      enddo
   enddo

   do k = l2, maxc
      km1 = k - 1
      do i = il1, il2
         dtr(i,1,km1)            =  dtr_vs(i,k-l2+1)
         ubeta                   =  urbf(i,km1) / (taul1(i,km1)+1.e-20)
         epsd                    =  ubeta + 1.0
         epsu                    =  ubeta - 1.0

         if (abs(epsd) .gt. 0.001) then
            xd(i,1,km1)           = (bf(i,k) - bf(i,km1) * &
                 dtr(i,1,km1)) / epsd
         else
            xd(i,1,km1)           = rtaul1(i,km1)*bf(i,km1)*dtr(i,1,km1)
         endif

         if (abs(epsu) .gt. 0.001) then
            xu(i,1,km1)           = (bf(i,k) * dtr(i,1,km1) - &
                 bf(i,km1)) / epsu
         else
            xu(i,1,km1)           = rtaul1(i,km1)*bf(i,k)*dtr(i,1,km1)
         endif

         fd(i,1,k)               =  fd(i,1,km1) * dtr(i,1,km1) + &
              xd(i,1,km1)
         fd(i,2,k)               =  fd(i,2,km1) * dtr(i,1,km1) + &
              xd(i,1,km1)
         fx(i,1,k)               =  fd(i,2,k)
         fx(i,2,k)               =  fd(i,2,k)
      enddo
   enddo

   !----------------------------------------------------------------------
   !     add the layers downward from the second layer to the surface.
   !     determine the xu for the upward path.
   !     using exponential source function for clr flux calculation and
   !     also for all sky flux in cloud free layers.
   !----------------------------------------------------------------------

   if (maxc .lt. lev) then
      DO250: do k = maxc + 1, lev
         km1 = k - 1
         km2 = km1 - 1
         DO225: do i = il1, il2
            dtr(i,1,km1)          =  dtr_vs(i,k-l2+1)
            ubeta                 =  urbf(i,km1)/(taul1(i,km1) + 1.e-20)
            epsd                  =  ubeta + 1.0
            epsu                  =  ubeta - 1.0

            if (abs(epsd) .gt. 0.001) then
               xd(i,1,km1)         = (bf(i,k) - bf(i,km1) * &
                    dtr(i,1,km1)) / epsd
            else
               xd(i,1,km1)         = rtaul1(i,km1)*bf(i,km1)*dtr(i,1,km1)
            endif

            if (abs(epsu) .gt. 0.001) then
               xu(i,1,km1)         = (bf(i,k) * dtr(i,1,km1) - &
                    bf(i,km1)) / epsu
            else
               xu(i,1,km1)         = rtaul1(i,km1)*bf(i,k)*dtr(i,1,km1)
            endif

            fd(i,1,k)             =  fd(i,1,km1) * dtr(i,1,km1) + &
                 xd(i,1,km1)
            if (cldfrac(i,km1) .lt. cut) then
               fd(i,2,k)           =  fd(i,2,km1) * dtr(i,1,km1) + &
                    xd(i,1,km1)
               fx(i,1,k)           =  fd(i,2,k)
               fx(i,2,k)           =  fd(i,2,k)
            else
               taul2               =  tauci(i,km1) + taug(i,km1)
               ssalb               =  omci(i,km1) / (taul2 + 1.e-20)
               sf                  =  ssalb * fl(i,km1)
               ww                  = (ssalb - sf) / (1.0 - sf)
               cow                 =  1.0 - ww
               taum(i,1,km1)       = (cow * taul2 * (1.0 - sf) + &
                    taual(i,km1)) * ru
               zeta                =  dbf(i,km1) / taum(i,1,km1)
               tau2                =  taum(i,1,km1) + taum(i,1,km1)

               sanu                =  anu(i,km1)
               xx                  =  sanu / (sanu + taum(i,1,km1))
               yy                  =  sanu / (sanu + tau2)

               if (sanu .le. 0.50) then
                  dtr(i,2,km1)      =  sqrt(xx)
                  dtr2              =  sqrt(yy)
                  emisw(i,km1)      =  zeta * (sqrt(1.0 + tau2) - 1.0)
                  embk              =  1.0 - sqrt(1.0 + tau2 + tau2)
               else if (sanu .gt. 0.50 .and. sanu .le. 1.0) then
                  wt                =  2.0 * sanu - 1.0
                  sx                =  sqrt(xx)
                  sy                =  sqrt(yy)
                  dtr(i,2,km1)      =  sx + (xx - sx) * wt
                  dtr2              =  sy + (yy - sy) * wt
                  p1                =  sqrt(1.0 + tau2) - 1.0
                  emisw(i,km1)      =  zeta * (p1 + (log(1.0 + &
                       taum(i,1,km1)) - p1) * wt)
                  p2                =  1.0 - sqrt(1.0 + tau2 + tau2)
                  embk              =  p2 - (log(1.0 + tau2) + p2) * wt
               else if (sanu .gt. 1.0 .and. sanu .le. 2.0) then
                  wt                =  sanu - 1.0
                  dtr(i,2,km1)      =  xx + (xx * xx - xx) * wt
                  dtr2              =  yy + (yy * yy - yy) * wt
                  zz                =  sanu / (sanu - 1.0)
                  p1                =  log(1.0 + taum(i,1,km1))
                  emisw(i,km1)      =  zeta * (p1 + (zz* (1.0 - xx) - p1)* &
                       wt)
                  p2                =  - log(1.0 + tau2)
                  embk              =  p2 + (zz * (yy - 1.0) - p2) * wt
               else if (sanu .gt. 2.0 .and. sanu .le. 3.0) then
                  x2                =  xx * xx
                  y2                =  yy * yy
                  wt                =  sanu - 2.0
                  dtr(i,2,km1)      =  x2 + (xx * x2 - x2) * wt
                  dtr2              =  y2 + (yy * y2 - y2) * wt
                  zz                =  sanu / (sanu - 1.0)
                  emisw(i,km1)      =  zz * zeta * &
                       (1.0 - xx + (xx - x2) * wt)
                  embk              =  zz * (yy - 1.0 + (y2 - yy) * wt)
               else if (sanu .gt. 3.0 .and. sanu .le. 4.0) then
                  x2                =  xx * xx
                  y2                =  yy * yy
                  x3                =  x2 * xx
                  y3                =  y2 * yy
                  wt                =  sanu - 3.0
                  dtr(i,2,km1)      =  x3 + (x2 * x2 - x3) * wt
                  dtr2              =  y3 + (y2 * y2 - y3) * wt
                  zz                =  sanu / (sanu - 1.0)
                  emisw(i,km1)      =  zz * zeta * &
                       (1.0 - x2 + (x2 - x3) * wt)
                  embk              =  zz * (y2 - 1.0 + (y3 - y2) * wt)

                  !----------------------------------------------------------------------
                  !     for anu > 4, the inhomoeneity effect is very weak, for saving
                  !     the integer anu is assumed. for anu > 20, homogenous is assumed
                  !----------------------------------------------------------------------

               else if (sanu .gt. 4.0 .and. sanu .le. 20.0) then
                  nanu              =  nint(sanu)
                  dtr(i,2,km1)      =  xx ** nanu
                  dtr2              =  yy ** nanu
                  zz                =  sanu / (sanu - 1.0)
                  emisw(i,km1)      =  zz* zeta * (1.0 - dtr(i,2,km1) /xx)
                  embk              =  zz* (dtr2 / yy - 1.0)
               else
                  dtr(i,2,km1)      =  exp(- taum(i,1,km1))
                  dtr2              =  exp(- tau2)
                  zz                =  sanu / (sanu - 1.0)
                  emisw(i,km1)      =  zz * zeta * (1.0 - dtr(i,2,km1) /xx)
                  embk              =  zz * (dtr2 / yy - 1.0)

               endif

               xd(i,2,km1)         =  bf(i,k) - bf(i,km1) * &
                    dtr(i,2,km1) - emisw(i,km1)
               xu(i,2,km1)         =  bf(i,km1) - bf(i,k) * &
                    dtr(i,2,km1) + emisw(i,km1)

               wgrcow              =  ww * gci(i,km1) / cow
               taudtr              =  taum(i,1,km1) * dtr(i,2,km1)

               scatfw(i,km1)       =  wgrcow * xx * taudtr
               scatbk(i,km1)       =  0.5 * wgrcow * (dtr2 - 1.0)

               xx                  =  wgrcow * (2.0 * emisw(i,km1) + &
                    (0.5 * embk - taudtr) * zeta)
               scatsm(i,1,km1)     =  - scatbk(i,km1) * bf(i,k) - &
                    scatfw(i,km1) * bf(i,km1) - xx
               scatsm(i,2,km1)     =  - scatbk(i,km1) * bf(i,km1) - &
                    scatfw(i,km1) * bf(i,k) + xx

               if (k .eq. l2) then
                  fx(i,1,k)         =  fx(i,1,km1) * dtr(i,1,km1) + &
                       xd(i,1,km1)
                  fx(i,2,k)         =  fx(i,2,km1) * dtr(i,2,km1) + &
                       xd(i,2,km1)
               else if (cldfrac(i,km1) .le. cldfrac(i,km2)) then
                  fx(i,1,k)         = (fx(i,2,km1)+(1.0-cldfrac(i,km2)) / &
                       max(1.0 - cldfrac(i,km1),1.e-10) * &
                       (fx(i,1,km1) - fx(i,2,km1)) ) * &
                       dtr(i,1,km1) + xd(i,1,km1)
                  fx(i,2,k)         =  fx(i,2,km1) * dtr(i,2,km1) + &
                       xd(i,2,km1)
               else if (cldfrac(i,km1) .gt. cldfrac(i,km2)) then
                  fx(i,1,k)         =  fx(i,1,km1) * dtr(i,1,km1) + &
                       xd(i,1,km1)
                  fx(i,2,k)         = (fx(i,1,km1) + &
                       cldfrac(i,km2) / cldfrac(i,km1) * &
                       (fx(i,2,km1) - fx(i,1,km1))) * &
                       dtr(i,2,km1) + xd(i,2,km1)
               endif

               fd(i,2,k)           =  fx(i,1,k) + cldfrac(i,km1) * &
                    (fx(i,2,k) - fx(i,1,k))
            endif
         enddo DO225
      enddo DO250
   endif

   !----------------------------------------------------------------------
   !     initialization for surface
   !----------------------------------------------------------------------

   k = lev - 1
   do i = il1, il2
      embs                      =  em0(i) * bs(i)
      abse0                     =  1.0 - em0(i)
      fu(i,1,lev)               =  embs + abse0 * fd(i,1,lev)
      fy(i,1,lev)               =  embs + abse0 * fx(i,1,lev)
      fy(i,2,lev)               =  embs + abse0 * fx(i,2,lev)
      fu(i,2,lev)               =  fy(i,1,lev) + cldfrac(i,k) * &
           (fy(i,2,lev) - fy(i,1,lev))
      fw(i,2,lev)               =  fy(i,2,lev)

      !----------------------------------------------------------------------
      !     determining the upward flux for the first lay above surface
      !----------------------------------------------------------------------

      fu(i,1,k)                 =  fu(i,1,lev) * dtr(i,1,k) + &
           xu(i,1,k)
      fy(i,1,k)                 =  fy(i,1,lev) * dtr(i,1,k) + &
           xu(i,1,k)


      if (cldfrac(i,k) .lt. cut) then
         taucsg(i,k)             =  0.0
         fy(i,2,k)               =  fy(i,2,lev) * dtr(i,1,k) + &
              xu(i,1,k)
         fu(i,2,k)               =  fy(i,1,k)
         fw(i,2,k)               =  fy(i,2,k)
         taum(i,2,k)             =  0.0
      else
         taucsg(i,k)             =  scatbk(i,k) * fx(i,2,k) + &
              scatfw(i,k) * fy(i,2,lev) + &
              scatsm(i,2,k)

         fy(i,2,k)               =  fy(i,2,lev) * dtr(i,2,k) + &
              xu(i,2,k) + taucsg(i,k)
         fu(i,2,k)               =  fy(i,1,k) + cldfrac(i,k) * &
              (fy(i,2,k) - fy(i,1,k))
         fw(i,2,k)               =  fy(i,2,k)
         taum(i,2,k)             =  taum(i,1,k)
      endif
   enddo

   !----------------------------------------------------------------------
   !     add the layers upward from the second layer to maxc
   !     scattering effect for upward path is included
   !----------------------------------------------------------------------

   DO450: do k = lev - 2, maxc, - 1
      kp1 = k + 1
      kp2 = k + 2
      DO400: do i = il1, il2
         if (k .ge. nct(i)) then
            fu(i,1,k)             =  fu(i,1,kp1) * dtr(i,1,k) + &
                 xu(i,1,k)

            if (cldfrac(i,k) .lt. cut) then
               fu(i,2,k)           =  fu(i,2,kp1) * dtr(i,1,k) + &
                    xu(i,1,k)
               fy(i,1,k)           =  fu(i,2,k)
               fy(i,2,k)           =  fu(i,2,k)
               fw(i,2,k)           =  fu(i,2,k)
               taum(i,2,k)         =  0.0
            else

               !----------------------------------------------------------------------
               !     fy(i,2,k) contains unperturbed + backward scattering effect +
               !     forward scattering effect + internal scattering effect
               !    (li and fu, jas 2000)
               !----------------------------------------------------------------------

               if (cldfrac(i,k) .le. cldfrac(i,kp1) .or. &
                    cldfrac(i,k) - cldm(i,k) .lt. cut) then

                  fy(i,1,k)         = ( fy(i,2,kp1)+(1.0-cldfrac(i,kp1)) / &
                       max(1.0 - cldfrac(i,k),1.e-10) * &
                       (fy(i,1,kp1) - fy(i,2,kp1)) ) * &
                       dtr(i,1,k) + xu(i,1,k)
                  c2(i)             =  fy(i,2,kp1)
               else
                  fy(i,1,k)         =  fy(i,1,kp1) * dtr(i,1,k) + &
                       xu(i,1,k)
                  c2(i)             =  fy(i,1,kp1) + &
                       (cldfrac(i,kp1) - cldm(i,kp1)) / &
                       (cldfrac(i,k) - cldm(i,k)) * &
                       (fy(i,2,kp1) - fy(i,1,kp1))
               endif

               bkins               =  scatbk(i,k) * fx(i,2,k) + &
                    scatsm(i,2,k)
               fy(i,2,k)           =  c2(i) * (dtr(i,2,k) + scatfw(i,k))+ &
                    xu(i,2,k) + bkins
               taum(i,2,k)         =  taum(i,2,kp1) + taum(i,1,k)
               s1(i)               =  0.0
               taucsg(i,k)         =  bkins + scatfw(i,k) * fy(i,2,kp1)
               term1(i)            =  0.0

               if (ncu(i,k) .gt. 1) then
                  kx = k + ncu(i,k)
                  kxm = kx - 1

                  sanu              =  anu(i,kxm)
                  anutau            =  sanu / (sanu + taum(i,2,k))
                  if (sanu .le. 0.50) then
                     dtrgw           =  sqrt(anutau)
                  else if (sanu .gt. 0.50 .and. sanu .le. 1.0) then
                     xx              =  sqrt(anutau)
                     dtrgw           =  xx + 2.0 * (sanu - 0.50) * &
                          (anutau - xx)
                  else if (sanu .gt. 1.0 .and. sanu .le. 2.0) then
                     dtrgw           =  anutau + (sanu - 1.0) * anutau * &
                          (anutau - 1.0)
                  else if (sanu .gt. 2.0 .and. sanu .le. 3.0) then
                     xx              =  anutau * anutau
                     dtrgw           =  xx + (sanu - 2.0) * xx * &
                          (anutau - 1.0)
                  else if (sanu .gt. 3.0 .and. sanu .le. 4.0) then
                     xx              =  anutau * anutau * anutau
                     dtrgw           =  xx + (sanu - 3.0) * xx * &
                          (anutau - 1.0)
                  else if (sanu .gt. 4.0 .and. sanu .le. 20.0) then
                     dtrgw           =  anutau ** (nint(sanu))
                  else
                     dtrgw           =  exp(- taum(i,2,k))
                  endif

                  term1(i)          = (fw(i,2,kx) - bf(i,kx)) * dtrgw
                  s1(i)             = (emisw(i,kp1) + taucsg(i,kp1)) * &
                       dtr(i,2,k)
               endif
            endif
         endif
      enddo DO400

      !----------------------------------------------------------------------
      !     determining the terms going into the correlation calculations
      !     for subgrid variability for cldm portion.
      !----------------------------------------------------------------------

      if (ncum(k) .gt. 2) then
         do kk = kp2, k + ncum(k) - 1
            do i = il1, il2
               if (k .ge. nct(i) .and. cldfrac(i,k) .ge. cut .and. &
                    ncu(i,k) .gt. 2 .and. kk .le. k + ncu(i,k) - 1) then

                  sanu                =  anu(i,kk)
                  anutau              =  sanu / (sanu + &
                       taum(i,2,k) - taum(i,2,kk))
                  if (sanu .le. 0.50) then
                     dtrgw             =  sqrt(anutau)
                  else if (sanu .gt. 0.50 .and. sanu .le. 1.0) then
                     xx                =  sqrt(anutau)
                     dtrgw             =  xx + 2.0 * (sanu - 0.50) * &
                          (anutau - xx)
                  else if (sanu .gt. 1.0 .and. sanu .le. 2.0) then
                     dtrgw             =  anutau + (sanu - 1.0) * anutau * &
                          (anutau - 1.0)
                  else if (sanu .gt. 2.0 .and. sanu .le. 3.0) then
                     xx                =  anutau * anutau
                     dtrgw             =  xx + (sanu - 2.0) * xx * &
                          (anutau - 1.0)
                  else if (sanu .gt. 3.0 .and. sanu .le. 4.0) then
                     xx                =  anutau * anutau * anutau
                     dtrgw             =  xx + (sanu - 3.0) * xx * &
                          (anutau - 1.0)
                  else if (sanu .gt. 4.0 .and. sanu .le. 20.0) then
                     dtrgw             =  anutau ** (nint(sanu))
                  else
                     dtrgw             =  exp(-(taum(i,2,k) - taum(i,2,kk)))
                  endif

                  s1(i)               =  s1(i) + &
                       (emisw(i,kk) + taucsg(i,kk)) * dtrgw
               endif
            enddo
         enddo
      endif

      !----------------------------------------------------------------------
      !     in cldm region consider the correlation between different layers
      !----------------------------------------------------------------------

      do i = il1, il2
         if (k .ge. nct(i)) then
            if (cldfrac(i,k) .ge. cut) then
               if (ncu(i,k) .eq. 1) then
                  fw(i,2,k)         =  fy(i,2,k)
                  fu(i,2,k)         =  fy(i,1,k)+cldfrac(i,k)*(fy(i,2,k) - &
                       fy(i,1,k))
               else
                  fw(i,2,k)         =  term1(i) + s1(i) + bf(i,k) + &
                       emisw(i,k) + taucsg(i,k)
                  fu(i,2,k)         =  cldm(i,k) * (fw(i,2,k) - &
                       fy(i,2,k)) + fy(i,1,k) + &
                       cldfrac(i,k)*(fy(i,2,k)-fy(i,1,k))
               endif
            endif
         endif
      enddo
   enddo DO450

   !----------------------------------------------------------------------
   !     add the layers upward above the highest cloud  to the toa, no
   !     scattering
   !----------------------------------------------------------------------

   do k = lev - 1, l1, - 1
      kp1 = k + 1

      do i = il1, il2
         if (kp1 .le. nct(i)) then
            fu(i,1,k)             =  fu(i,1,kp1) * dtr(i,1,k) + &
                 xu(i,1,k)
            fu(i,2,k)             =  fu(i,2,kp1) * dtr(i,1,k) + &
                 xu(i,1,k)
         endif

         !----------------------------------------------------------------------
         !     scattering effect for downward path at the top layer of the
         !     highest cloud
         !----------------------------------------------------------------------

         if (k .eq. nct(i)) then
            fw(i,1,k)             =  fx(i,1,k)
            fwins                 =  scatsm(i,1,k) + &
                 scatfw(i,k) * fx(i,2,k)
            fmbk                  =  fx(i,2,k) * dtr(i,2,k) + &
                 xd(i,2,k) + fwins
            fx(i,2,kp1)           =  fmbk + scatbk(i,k) * fy(i,2,kp1)
            taum(i,2,k)           =  taum(i,1,k)
            taucsg(i,k)           =  scatbk(i,k) * fw(i,2,kp1) + fwins

            fw(i,1,kp1)           =  fmbk + scatbk(i,k) * fw(i,2,kp1)
            fd(i,2,kp1)           =  fx(i,1,kp1) + cldfrac(i,k) * &
                 (fx(i,2,kp1) - fx(i,1,kp1))
         endif
      enddo
   enddo

   !----------------------------------------------------------------------
   !     scattering effect for downward path in from maxc to the surface
   !----------------------------------------------------------------------

   DO750: do k = maxc + 2, lev
      km1 = k - 1
      km2 = k - 2
      km3 = k - 3
      DO700: do i = il1, il2
         if (km2 .ge. nct(i)) then
            if (cldfrac(i,km1) .lt. cut) then
               fd(i,2,k)           =  fd(i,2,km1) * dtr(i,1,km1) + &
                    xd(i,1,km1)
               fx(i,1,k)           =  fd(i,2,k)
               fx(i,2,k)           =  fd(i,2,k)
               fw(i,1,k)           =  fd(i,2,k)
               taum(i,2,km1)       =  0.0
            else
               if (cldfrac(i,km1) .le. cldfrac(i,km2) .or. &
                    cldfrac(i,km1) - cldm(i,km1) .lt. cut) then

                  fx(i,1,k)         = (fx(i,2,km1)+(1.0-cldfrac(i,km2)) / &
                       max(1.0 - cldfrac(i,km1),1.e-10) * &
                       (fx(i,1,km1) - fx(i,2,km1))) * &
                       dtr(i,1,km1) + xd(i,1,km1)
                  c2(i)             =  fx(i,2,km1)
               else
                  fx(i,1,k)         =  fx(i,1,km1) * dtr(i,1,km1) + &
                       xd(i,1,km1)
                  c2(i)             =  fx(i,1,km1) + &
                       (cldfrac(i,km2) - cldm(i,km2)) / &
                       (cldfrac(i,km1) - cldm(i,km1)) * &
                       (fx(i,2,km1) -  fx(i,1,km1))
               endif

               fx(i,2,k)           =  c2(i) * dtr(i,2,km1) + xd(i,2,km1)+ &
                    scatbk(i,km1) * fy(i,2,k) + &
                    scatfw(i,km1) * c2(i) + &
                    scatsm(i,1,km1)

               taum(i,2,km1)       =  taum(i,2,km2) + taum(i,1,km1)
               s1(i)               =  0.0
               taucsg(i,km1)       =  scatbk(i,km1) * fw(i,2,k) + &
                    scatfw(i,km1) * fw(i,1,km1) + &
                    scatsm(i,1,km1)
               term1(i)            =  0.0

               if (ncd(i,km1) .gt. 1) then
                  kx = k - ncd(i,km1)
                  sanu              =  anu(i,kx)
                  anutau            =  sanu / (sanu + taum(i,2,km1))
                  if (sanu .le. 0.50) then
                     dtrgw           =  sqrt(anutau)
                  else if (sanu .gt. 0.50 .and. sanu .le. 1.0) then
                     xx              =  sqrt(anutau)
                     dtrgw           =  xx + 2.0 * (sanu - 0.50) * &
                          (anutau - xx)
                  else if (sanu .gt. 1.0 .and. sanu .le. 2.0) then
                     dtrgw           =  anutau + (sanu - 1.0) * anutau * &
                          (anutau - 1.0)
                  else if (sanu .gt. 2.0 .and. sanu .le. 3.0) then
                     xx              =  anutau * anutau
                     dtrgw           =  xx + (sanu - 2.0) * xx * &
                          (anutau - 1.0)
                  else if (sanu .gt. 3.0 .and. sanu .le. 4.0) then
                     xx              =  anutau * anutau * anutau
                     dtrgw           =  xx + (sanu - 3.0) * xx * &
                          (anutau - 1.0)
                  else if (sanu .gt. 4.0 .and. sanu .le. 20.0) then
                     dtrgw           =  anutau ** (nint(sanu))
                  else
                     dtrgw           =  exp(- taum(i,2,km1))
                  endif

                  term1(i)          = (fw(i,1,kx) - bf(i,kx)) * dtrgw
                  s1(i)             = (taucsg(i,km2) - emisw(i,km2)) * &
                       dtr(i,2,km1)
               endif
            endif
         endif
      enddo DO700
      !----------------------------------------------------------------------
      !     determining the terms going into the correlation calculations
      !     for cldm portion.
      !----------------------------------------------------------------------

      if (ncdm(km1) .gt. 2) then

         !----------------------------------------------------------------------
         !     note that in the following loop, "km1" is actually the
         !     representative variable, so that k-ncd(i,km1) is actually
         !     km1-ncd(i,km1)+1. the simpler form is used only for
         !     computational efficiency.
         !----------------------------------------------------------------------

         do kk = km3, k - ncdm(km1), - 1
            do i = il1, il2
               if (km2 .ge. nct(i) .and. cldfrac(i,km1) .ge. cut .and. &
                    ncd(i,km1) .gt. 2 .and. kk .ge. k - ncd(i,km1)) then

                  sanu                =  anu(i,kk)
                  anutau              =  sanu / (sanu + &
                       taum(i,2,km1) - taum(i,2,kk))
                  if (sanu .le. 0.50) then
                     dtrgw             =  sqrt(anutau)
                  else if (sanu .gt. 0.50 .and. sanu .le. 1.0) then
                     xx                =  sqrt(anutau)
                     dtrgw             =  xx + 2.0 * (sanu - 0.50) * &
                          (anutau - xx)
                  else if (sanu .gt. 1.0 .and. sanu .le. 2.0) then
                     dtrgw             =  anutau + (sanu - 1.0) * anutau * &
                          (anutau - 1.0)
                  else if (sanu .gt. 2.0 .and. sanu .le. 3.0) then
                     xx                =  anutau * anutau
                     dtrgw             =  xx + (sanu - 2.0) * xx * &
                          (anutau - 1.0)
                  else if (sanu .gt. 3.0 .and. sanu .le. 4.0) then
                     xx                =  anutau * anutau * anutau
                     dtrgw             =  xx + (sanu - 3.0) * xx * &
                          (anutau - 1.0)
                  else if (sanu .gt. 4.0 .and. sanu .le. 20.0) then
                     dtrgw             =  anutau ** (nint(sanu))
                  else
                     dtrgw             =  exp(-(taum(i,2,km1)-taum(i,2,kk)))
                  endif

                  s1(i)               =  s1(i) - &
                       (emisw(i,kk) - taucsg(i,kk)) * dtrgw
               endif
            enddo
         enddo
      endif

      do i = il1, il2
         if (km2 .ge. nct(i)) then
            if (cldfrac(i,km1) .ge. cut) then
               if (ncd(i,km1) .eq. 1) then
                  fw(i,1,k)         =  fx(i,2,k)
                  fd(i,2,k)         =  fx(i,1,k) + cldfrac(i,km1) * &
                       (fx(i,2,k) - fx(i,1,k))
               else
                  fw(i,1,k)         =  term1(i) + s1(i) + bf(i,k) - &
                       emisw(i,km1) + taucsg(i,km1)
                  fd(i,2,k)         =  cldm(i,km1) * &
                       (fw(i,1,k) - fx(i,2,k)) + &
                       fx(i,1,k) + cldfrac(i,km1) * &
                       (fx(i,2,k) - fx(i,1,k))
               endif
            endif
         endif
      enddo
   enddo DO750

   return
end subroutine ccc2_lwtran4
