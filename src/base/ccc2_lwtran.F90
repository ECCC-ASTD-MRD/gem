
module ccc2_lwtran
   implicit none
   private
   public :: ccc2_lwtran1
contains
   
subroutine ccc2_lwtran1(fu, fd, &
     slwf, tauci, omci, gci, fl, taual, taug, bf, &
     bs, urbf, dbf, em0, cldfrac, cldm, anu, &
     nct, ncd, ncu, ncum, ncdm, &
     lev1, cut, maxc, ni, lay, lev)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer, intent(in) :: lev1, maxc, ni, lay, lev
   integer, intent(in) :: nct(ni), ncd(ni,lay), ncu(ni,lay), ncum(lay), ncdm(lay)
   real, intent(in) :: cut
   real, intent(in) :: slwf(ni), tauci(ni,lay), omci(ni,lay), gci(ni,lay), &
        fl(ni,lay), taual(ni,lay), taug(ni,lay), bf(ni,lev), &
        bs(ni), urbf(ni,lay), dbf(ni,lay), em0(ni), cldfrac(ni,lay), &
        cldm(ni,lay), anu(ni,lay)
   real, intent(out) :: fu(ni,2,lev), fd(ni,2,lev)
   
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

   !          - ?internal? -   
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
   !          latitude loop (ni)
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
   
   integer :: l1, l2, i
   integer :: k, km1, km2, kp1, kp2, kx, kxm, kk, km3, nanu
   real :: ubeta, epsd, epsu, zeta, taul2, ssalb
   real :: sf, ww, cow, tau2, sanu, xx, yy, dtr2, embk
   real :: wt, sx, sy, p1, p2, zz, x2, y2, x3, y3, wgrcow
   real :: taudtr, bkins, anutau, dtrgw, fwins, fmbk
   real :: taul1(ni,lay), rtaul1(ni,lay)
   real :: taucsg(ni,lay), term1(ni), emisw(ni,lay), scatbk(ni,lay), &
        scatfw(ni,lay), scatsm(ni,4,lay), taum(ni,4,lay), &
        xd(ni,4,lay), xu(ni,4,lay), dtr(ni,4,lay), fx(ni,4,lev), &
        fy(ni,4,lev), fw(ni,4,lev), s1(ni), c2(ni)
   real :: embs,abse0
   real(REAL64) :: dtr_vs(ni,lay)

   real, parameter :: ru = 1.6487213
   !----------------------------------------------------------------------
   
   l1 = lev1
   l2 = lev1 + 1
   do i = 1, ni
      fd(i,1,lev1) = slwf(i)
      fd(i,2,lev1) = slwf(i)
      fx(i,1,lev1) = fd(i,1,lev1)
      fx(i,2,lev1) = fd(i,2,lev1)
   enddo

   do k = l2, lev
      km1 = k - 1
      do i = 1, ni
         taul1(i,km1)     = taual(i,km1) + taug(i,km1)
         rtaul1(i,km1)    = taul1(i,km1) * ru
         dtr_vs(i,k-l2+1) = exp(dble(-rtaul1(i,km1)))  !#TODO: R8 needed? always copied to a R4
      enddo
   enddo

   do k = l2, maxc
      km1 = k - 1
      do i = 1, ni
         dtr(i,1,km1) = dtr_vs(i,k-l2+1)
         ubeta        = urbf(i,km1) / (taul1(i,km1)+1.e-20)
         epsd         = ubeta + 1.0
         epsu         = ubeta - 1.0

         if (abs(epsd) > 0.001) then
            xd(i,1,km1) = (bf(i,k) - bf(i,km1) * dtr(i,1,km1)) / epsd
         else
            xd(i,1,km1) = rtaul1(i,km1)*bf(i,km1)*dtr(i,1,km1)
         endif

         if (abs(epsu) > 0.001) then
            xu(i,1,km1) = (bf(i,k) * dtr(i,1,km1) - bf(i,km1)) / epsu
         else
            xu(i,1,km1) = rtaul1(i,km1)*bf(i,k)*dtr(i,1,km1)
         endif

         fd(i,1,k) = fd(i,1,km1) * dtr(i,1,km1) + xd(i,1,km1)
         fd(i,2,k) = fd(i,2,km1) * dtr(i,1,km1) + xd(i,1,km1)
         fx(i,1,k) = fd(i,2,k)
         fx(i,2,k) = fd(i,2,k)
      enddo
   enddo

   !----------------------------------------------------------------------
   !     add the layers downward from the second layer to the surface.
   !     determine the xu for the upward path.
   !     using exponential source function for clr flux calculation and
   !     also for all sky flux in cloud free layers.
   !----------------------------------------------------------------------

   IF_LEV: if (maxc < lev) then
      DO_K_250: do k = maxc + 1, lev
         km1 = k - 1
         km2 = km1 - 1
         !#TODO: 3rd most expensive loop
         DO_I_225: do i = 1, ni
            dtr(i,1,km1) = dtr_vs(i,k-l2+1)
            ubeta        = urbf(i,km1)/(taul1(i,km1) + 1.e-20)
            epsd         = ubeta + 1.0
            epsu         = ubeta - 1.0

            if (abs(epsd) > 0.001) then
               xd(i,1,km1) = (bf(i,k) - bf(i,km1) * dtr(i,1,km1)) / epsd
            else
               xd(i,1,km1) = rtaul1(i,km1)*bf(i,km1)*dtr(i,1,km1)
            endif

            if (abs(epsu) > 0.001) then
               xu(i,1,km1) = (bf(i,k) * dtr(i,1,km1) - bf(i,km1)) / epsu
            else
               xu(i,1,km1) = rtaul1(i,km1)*bf(i,k)*dtr(i,1,km1)
            endif

            fd(i,1,k)      = fd(i,1,km1) * dtr(i,1,km1) + xd(i,1,km1)
            
            IF_CUT: if (cldfrac(i,km1) < cut) then
               fd(i,2,k)   = fd(i,2,km1) * dtr(i,1,km1) + xd(i,1,km1)
               fx(i,1,k)   = fd(i,2,k)
               fx(i,2,k)   = fd(i,2,k)
            else
               taul2       = tauci(i,km1) + taug(i,km1)
               ssalb       = omci(i,km1) / (taul2 + 1.e-20)
               sf          = ssalb * fl(i,km1)
               ww          = (ssalb - sf) / (1.0 - sf)
               cow         = 1.0 - ww
               taum(i,1,km1) = (cow * taul2 * (1.0 - sf) + taual(i,km1)) * ru
               zeta        = dbf(i,km1) / taum(i,1,km1)
               tau2        = taum(i,1,km1) + taum(i,1,km1)

               sanu        = anu(i,km1)
               xx          = sanu / (sanu + taum(i,1,km1))
               yy          = sanu / (sanu + tau2)

               IF_SANU: if (sanu <= 0.50) then
                  dtr(i,2,km1)  = sqrt(xx)
                  dtr2          = sqrt(yy)
                  emisw(i,km1)  = zeta * (sqrt(1.0 + tau2) - 1.0)
                  embk          = 1.0 - sqrt(1.0 + tau2 + tau2)
               else if (sanu > 0.50 .and. sanu <= 1.0) then
                  wt            = 2.0 * sanu - 1.0
                  sx            = sqrt(xx)
                  sy            = sqrt(yy)
                  dtr(i,2,km1)  = sx + (xx - sx) * wt
                  dtr2          = sy + (yy - sy) * wt
                  p1            = sqrt(1.0 + tau2) - 1.0
                  emisw(i,km1)  = zeta * (p1 + (log(1.0 + &
                       taum(i,1,km1)) - p1) * wt)
                  p2            = 1.0 - sqrt(1.0 + tau2 + tau2)
                  embk          = p2 - (log(1.0 + tau2) + p2) * wt
               else if (sanu > 1.0 .and. sanu <= 2.0) then
                  wt            = sanu - 1.0
                  dtr(i,2,km1)  = xx + (xx * xx - xx) * wt
                  dtr2          = yy + (yy * yy - yy) * wt
                  zz            = sanu / (sanu - 1.0)
                  p1            = log(1.0 + taum(i,1,km1))
                  emisw(i,km1)  = zeta * (p1 + (zz* (1.0 - xx) - p1)* wt)
                  p2            = - log(1.0 + tau2)
                  embk          = p2 + (zz * (yy - 1.0) - p2) * wt
               else if (sanu > 2.0 .and. sanu <= 3.0) then
                  x2            = xx * xx
                  y2            = yy * yy
                  wt            = sanu - 2.0
                  dtr(i,2,km1)  = x2 + (xx * x2 - x2) * wt
                  dtr2          = y2 + (yy * y2 - y2) * wt
                  zz            = sanu / (sanu - 1.0)
                  emisw(i,km1)  = zz * zeta * (1.0 - xx + (xx - x2) * wt)
                  embk          = zz * (yy - 1.0 + (y2 - yy) * wt)
               else if (sanu > 3.0 .and. sanu <= 4.0) then
                  x2            = xx * xx
                  y2            = yy * yy
                  x3            = x2 * xx
                  y3            = y2 * yy
                  wt            = sanu - 3.0
                  dtr(i,2,km1)  = x3 + (x2 * x2 - x3) * wt
                  dtr2          = y3 + (y2 * y2 - y3) * wt
                  zz            = sanu / (sanu - 1.0)
                  emisw(i,km1)  = zz * zeta * (1.0 - x2 + (x2 - x3) * wt)
                  embk          = zz * (y2 - 1.0 + (y3 - y2) * wt)

                  !----------------------------------------------------------------------
                  !     for anu > 4, the inhomoeneity effect is very weak, for saving
                  !     the integer anu is assumed. for anu > 20, homogenous is assumed
                  !----------------------------------------------------------------------

               else if (sanu > 4.0 .and. sanu <= 20.0) then
                  nanu          = nint(sanu)
                  dtr(i,2,km1)  = xx ** nanu
                  dtr2          = yy ** nanu
                  zz            = sanu / (sanu - 1.0)
                  emisw(i,km1)  = zz* zeta * (1.0 - dtr(i,2,km1) /xx)
                  embk          = zz* (dtr2 / yy - 1.0)
               else
                  dtr(i,2,km1)  = exp(- taum(i,1,km1))
                  dtr2          = exp(- tau2)
                  zz            = sanu / (sanu - 1.0)
                  emisw(i,km1)  = zz * zeta * (1.0 - dtr(i,2,km1) /xx)
                  embk          = zz * (dtr2 / yy - 1.0)
               endif IF_SANU

               xd(i,2,km1)     = bf(i,k) - bf(i,km1) * &
                    dtr(i,2,km1) - emisw(i,km1)
               xu(i,2,km1)     = bf(i,km1) - bf(i,k) * &
                    dtr(i,2,km1) + emisw(i,km1)

               wgrcow          = ww * gci(i,km1) / cow
               taudtr          = taum(i,1,km1) * dtr(i,2,km1)

               scatfw(i,km1)   = wgrcow * xx * taudtr
               scatbk(i,km1)   = 0.5 * wgrcow * (dtr2 - 1.0)

               xx              = wgrcow * (2.0 * emisw(i,km1) + &
                    (0.5 * embk - taudtr) * zeta)
               scatsm(i,1,km1) = - scatbk(i,km1) * bf(i,k) - &
                    scatfw(i,km1) * bf(i,km1) - xx
               scatsm(i,2,km1) = - scatbk(i,km1) * bf(i,km1) - &
                    scatfw(i,km1) * bf(i,k) + xx

               if (k == l2) then
                  fx(i,1,k)     = fx(i,1,km1) * dtr(i,1,km1) + xd(i,1,km1)
                  fx(i,2,k)     = fx(i,2,km1) * dtr(i,2,km1) + xd(i,2,km1)
               else if (cldfrac(i,km1) <= cldfrac(i,km2)) then
                  fx(i,1,k)     = (fx(i,2,km1)+(1.0-cldfrac(i,km2)) / &
                       max(1.0 - cldfrac(i,km1),1.e-10) * &
                       (fx(i,1,km1) - fx(i,2,km1)) ) * &
                       dtr(i,1,km1) + xd(i,1,km1)
                  fx(i,2,k)     = fx(i,2,km1) * dtr(i,2,km1) + xd(i,2,km1)
               else if (cldfrac(i,km1) > cldfrac(i,km2)) then
                  fx(i,1,k)     = fx(i,1,km1) * dtr(i,1,km1) + xd(i,1,km1)
                  fx(i,2,k)     = (fx(i,1,km1) + &
                       cldfrac(i,km2) / cldfrac(i,km1) * &
                       (fx(i,2,km1) - fx(i,1,km1))) * &
                       dtr(i,2,km1) + xd(i,2,km1)
               endif

               fd(i,2,k)       = fx(i,1,k) + cldfrac(i,km1) * &
                    (fx(i,2,k) - fx(i,1,k))
            endif IF_CUT
         enddo DO_I_225
      enddo DO_K_250
   endif IF_LEV

   !----------------------------------------------------------------------
   !     initialization for surface
   !----------------------------------------------------------------------

   k = lev - 1
   do i = 1, ni
      embs        = em0(i) * bs(i)
      abse0       = 1.0 - em0(i)
      fu(i,1,lev) = embs + abse0 * fd(i,1,lev)
      fy(i,1,lev) = embs + abse0 * fx(i,1,lev)
      fy(i,2,lev) = embs + abse0 * fx(i,2,lev)
      fu(i,2,lev) = fy(i,1,lev) + cldfrac(i,k) * (fy(i,2,lev) - fy(i,1,lev))
      fw(i,2,lev) = fy(i,2,lev)

      !----------------------------------------------------------------------
      !     determining the upward flux for the first lay above surface
      !----------------------------------------------------------------------

      fu(i,1,k) = fu(i,1,lev) * dtr(i,1,k) + xu(i,1,k)
      fy(i,1,k) = fy(i,1,lev) * dtr(i,1,k) + xu(i,1,k)


      if (cldfrac(i,k) < cut) then
         taucsg(i,k) = 0.0
         fy(i,2,k)   = fy(i,2,lev) * dtr(i,1,k) + xu(i,1,k)
         fu(i,2,k)   = fy(i,1,k)
         fw(i,2,k)   = fy(i,2,k)
         taum(i,2,k) = 0.0
      else
         taucsg(i,k) = scatbk(i,k) * fx(i,2,k) + &
              scatfw(i,k) * fy(i,2,lev) + scatsm(i,2,k)

         fy(i,2,k)   = fy(i,2,lev) * dtr(i,2,k) + xu(i,2,k) + taucsg(i,k)
         fu(i,2,k)   = fy(i,1,k) + cldfrac(i,k) * (fy(i,2,k) - fy(i,1,k))
         fw(i,2,k)   = fy(i,2,k)
         taum(i,2,k) = taum(i,1,k)
      endif
   enddo

   !----------------------------------------------------------------------
   !     add the layers upward from the second layer to maxc
   !     scattering effect for upward path is included
   !----------------------------------------------------------------------

   DO_K_450: do k = lev - 2, maxc, - 1
      kp1 = k + 1
      kp2 = k + 2
      DO_I_400: do i = 1, ni
         if (k >= nct(i)) then
            fu(i,1,k) = fu(i,1,kp1) * dtr(i,1,k) + xu(i,1,k)

            if (cldfrac(i,k) < cut) then
               fu(i,2,k) = fu(i,2,kp1) * dtr(i,1,k) + xu(i,1,k)
               fy(i,1,k) = fu(i,2,k)
               fy(i,2,k) = fu(i,2,k)
               fw(i,2,k) = fu(i,2,k)
               taum(i,2,k) = 0.0
            else

               !----------------------------------------------------------------------
               !     fy(i,2,k) contains unperturbed + backward scattering effect +
               !     forward scattering effect + internal scattering effect
               !    (li and fu, jas 2000)
               !----------------------------------------------------------------------

               if (cldfrac(i,k) <= cldfrac(i,kp1) .or. &
                    cldfrac(i,k) - cldm(i,k) < cut) then

                  fy(i,1,k) = ( fy(i,2,kp1)+(1.0-cldfrac(i,kp1)) / &
                       max(1.0 - cldfrac(i,k),1.e-10) * &
                       (fy(i,1,kp1) - fy(i,2,kp1)) ) * dtr(i,1,k) + xu(i,1,k)
                  c2(i)     = fy(i,2,kp1)
               else
                  fy(i,1,k) = fy(i,1,kp1) * dtr(i,1,k) + xu(i,1,k)
                  c2(i)     = fy(i,1,kp1) + &
                       (cldfrac(i,kp1) - cldm(i,kp1)) / &
                       (cldfrac(i,k) - cldm(i,k)) * &
                       (fy(i,2,kp1) - fy(i,1,kp1))
               endif

               bkins       = scatbk(i,k) * fx(i,2,k) + scatsm(i,2,k)
               fy(i,2,k)   = c2(i) * (dtr(i,2,k) + scatfw(i,k))+ &
                    xu(i,2,k) + bkins
               taum(i,2,k) = taum(i,2,kp1) + taum(i,1,k)
               s1(i)       = 0.0
               taucsg(i,k) = bkins + scatfw(i,k) * fy(i,2,kp1)
               term1(i)    = 0.0

               if (ncu(i,k) > 1) then
                  kx = k + ncu(i,k)
                  kxm = kx - 1

                  sanu   = anu(i,kxm)
                  anutau = sanu / (sanu + taum(i,2,k))
                  if (sanu <= 0.50) then
                     dtrgw = sqrt(anutau)
                  else if (sanu > 0.50 .and. sanu <= 1.0) then
                     xx    = sqrt(anutau)
                     dtrgw = xx + 2.0 * (sanu - 0.50) * (anutau - xx)
                  else if (sanu > 1.0 .and. sanu <= 2.0) then
                     dtrgw = anutau + (sanu - 1.0) * anutau * (anutau - 1.0)
                  else if (sanu > 2.0 .and. sanu <= 3.0) then
                     xx    = anutau * anutau
                     dtrgw = xx + (sanu - 2.0) * xx * (anutau - 1.0)
                  else if (sanu > 3.0 .and. sanu <= 4.0) then
                     xx    = anutau * anutau * anutau
                     dtrgw = xx + (sanu - 3.0) * xx * (anutau - 1.0)
                  else if (sanu > 4.0 .and. sanu <= 20.0) then
                     dtrgw = anutau ** (nint(sanu))
                  else
                     dtrgw = exp(- taum(i,2,k))
                  endif

                  term1(i) = (fw(i,2,kx) - bf(i,kx)) * dtrgw
                  s1(i)    = (emisw(i,kp1) + taucsg(i,kp1)) * dtr(i,2,k)
               endif
            endif
         endif
      enddo DO_I_400

      !----------------------------------------------------------------------
      !     determining the terms going into the correlation calculations
      !     for subgrid variability for cldm portion.
      !----------------------------------------------------------------------

      if (ncum(k) > 2) then
         DO_KK: do kk = kp2, k + ncum(k) - 1
            !#TODO: 2nd most expensive loop
            do i = 1, ni
               if (k >= nct(i) .and. cldfrac(i,k) >= cut .and. &
                    ncu(i,k) > 2 .and. kk <= k + ncu(i,k) - 1) then

                  sanu     = anu(i,kk)
                  anutau   = sanu / (sanu + taum(i,2,k) - taum(i,2,kk))
                  if (sanu <= 0.50) then
                     dtrgw = sqrt(anutau)
                  else if (sanu > 0.50 .and. sanu <= 1.0) then
                     xx    = sqrt(anutau)
                     dtrgw = xx + 2.0 * (sanu - 0.50) * (anutau - xx)
                  else if (sanu > 1.0 .and. sanu <= 2.0) then
                     dtrgw = anutau + (sanu - 1.0) * anutau * (anutau - 1.0)
                  else if (sanu > 2.0 .and. sanu <= 3.0) then
                     xx    = anutau * anutau
                     dtrgw = xx + (sanu - 2.0) * xx * (anutau - 1.0)
                  else if (sanu > 3.0 .and. sanu <= 4.0) then
                     xx    = anutau * anutau * anutau
                     dtrgw = xx + (sanu - 3.0) * xx * (anutau - 1.0)
                  else if (sanu > 4.0 .and. sanu <= 20.0) then
                     dtrgw = anutau ** (nint(sanu))
                  else
                     dtrgw = exp(-(taum(i,2,k) - taum(i,2,kk)))
                  endif

                  s1(i)    = s1(i) + (emisw(i,kk) + taucsg(i,kk)) * dtrgw
               endif
            enddo
         enddo DO_KK
      endif

      !----------------------------------------------------------------------
      !     in cldm region consider the correlation between different layers
      !----------------------------------------------------------------------

      do i = 1, ni
         if (k >= nct(i)) then
            if (cldfrac(i,k) >= cut) then
               if (ncu(i,k) == 1) then
                  fw(i,2,k) = fy(i,2,k)
                  fu(i,2,k) = fy(i,1,k)+cldfrac(i,k)*(fy(i,2,k) - fy(i,1,k))
               else
                  fw(i,2,k) = term1(i) + s1(i) + bf(i,k) + emisw(i,k) + taucsg(i,k)
                  fu(i,2,k) = cldm(i,k) * (fw(i,2,k) - &
                       fy(i,2,k)) + fy(i,1,k) + &
                       cldfrac(i,k)*(fy(i,2,k)-fy(i,1,k))
               endif
            endif
         endif
      enddo
   enddo DO_K_450

   !----------------------------------------------------------------------
   !     add the layers upward above the highest cloud  to the toa, no
   !     scattering
   !----------------------------------------------------------------------

   do k = lev - 1, l1, - 1
      kp1 = k + 1

      do i = 1, ni
         if (kp1 <= nct(i)) then
            fu(i,1,k) = fu(i,1,kp1) * dtr(i,1,k) + xu(i,1,k)
            fu(i,2,k) = fu(i,2,kp1) * dtr(i,1,k) + xu(i,1,k)
         endif

         !----------------------------------------------------------------------
         !     scattering effect for downward path at the top layer of the
         !     highest cloud
         !----------------------------------------------------------------------

         if (k == nct(i)) then
            fw(i,1,k)   = fx(i,1,k)
            fwins       = scatsm(i,1,k) + scatfw(i,k) * fx(i,2,k)
            fmbk        = fx(i,2,k) * dtr(i,2,k) + xd(i,2,k) + fwins
            fx(i,2,kp1) = fmbk + scatbk(i,k) * fy(i,2,kp1)
            taum(i,2,k) = taum(i,1,k)
            taucsg(i,k) = scatbk(i,k) * fw(i,2,kp1) + fwins

            fw(i,1,kp1) = fmbk + scatbk(i,k) * fw(i,2,kp1)
            fd(i,2,kp1) = fx(i,1,kp1) + cldfrac(i,k) * (fx(i,2,kp1) - fx(i,1,kp1))
         endif
      enddo
   enddo

   !----------------------------------------------------------------------
   !     scattering effect for downward path in from maxc to the surface
   !----------------------------------------------------------------------

   DO_K_750: do k = maxc + 2, lev
      km1 = k - 1
      km2 = k - 2
      km3 = k - 3
      DO_I_700: do i = 1, ni
         if (km2 >= nct(i)) then
            if (cldfrac(i,km1) < cut) then
               fd(i,2,k) = fd(i,2,km1) * dtr(i,1,km1) + xd(i,1,km1)
               fx(i,1,k) = fd(i,2,k)
               fx(i,2,k) = fd(i,2,k)
               fw(i,1,k) = fd(i,2,k)
               taum(i,2,km1) = 0.0
            else
               if (cldfrac(i,km1) <= cldfrac(i,km2) .or. &
                    cldfrac(i,km1) - cldm(i,km1) < cut) then

                  fx(i,1,k) = (fx(i,2,km1)+(1.0-cldfrac(i,km2)) / &
                       max(1.0 - cldfrac(i,km1),1.e-10) * &
                       (fx(i,1,km1) - fx(i,2,km1))) * &
                       dtr(i,1,km1) + xd(i,1,km1)
                  c2(i)     = fx(i,2,km1)
               else
                  fx(i,1,k) = fx(i,1,km1) * dtr(i,1,km1) + xd(i,1,km1)
                  c2(i)     = fx(i,1,km1) + &
                       (cldfrac(i,km2) - cldm(i,km2)) / &
                       (cldfrac(i,km1) - cldm(i,km1)) * &
                       (fx(i,2,km1) -  fx(i,1,km1))
               endif

               fx(i,2,k)    = c2(i) * dtr(i,2,km1) + xd(i,2,km1)+ &
                    scatbk(i,km1) * fy(i,2,k) + &
                    scatfw(i,km1) * c2(i) + &
                    scatsm(i,1,km1)

               taum(i,2,km1) = taum(i,2,km2) + taum(i,1,km1)
               s1(i)         = 0.0
               taucsg(i,km1) = scatbk(i,km1) * fw(i,2,k) + &
                    scatfw(i,km1) * fw(i,1,km1) + scatsm(i,1,km1)
               term1(i)      = 0.0

               if (ncd(i,km1) > 1) then
                  kx = k - ncd(i,km1)
                  sanu   = anu(i,kx)
                  anutau = sanu / (sanu + taum(i,2,km1))
                  if (sanu <= 0.50) then
                     dtrgw = sqrt(anutau)
                  else if (sanu > 0.50 .and. sanu <= 1.0) then
                     xx    = sqrt(anutau)
                     dtrgw = xx + 2.0 * (sanu - 0.50) * (anutau - xx)
                  else if (sanu > 1.0 .and. sanu <= 2.0) then
                     dtrgw = anutau + (sanu - 1.0) * anutau * (anutau - 1.0)
                  else if (sanu > 2.0 .and. sanu <= 3.0) then
                     xx    = anutau * anutau
                     dtrgw = xx + (sanu - 2.0) * xx * (anutau - 1.0)
                  else if (sanu > 3.0 .and. sanu <= 4.0) then
                     xx    = anutau * anutau * anutau
                     dtrgw = xx + (sanu - 3.0) * xx * (anutau - 1.0)
                  else if (sanu > 4.0 .and. sanu <= 20.0) then
                     dtrgw = anutau ** (nint(sanu))
                  else
                     dtrgw = exp(- taum(i,2,km1))
                  endif

                  term1(i) = (fw(i,1,kx) - bf(i,kx)) * dtrgw
                  s1(i)    = (taucsg(i,km2) - emisw(i,km2)) * dtr(i,2,km1)
               endif
            endif
         endif
      enddo DO_I_700
      !----------------------------------------------------------------------
      !     determining the terms going into the correlation calculations
      !     for cldm portion.
      !----------------------------------------------------------------------

      if (ncdm(km1) > 2) then

         !----------------------------------------------------------------------
         !     note that in the following loop, "km1" is actually the
         !     representative variable, so that k-ncd(i,km1) is actually
         !     km1-ncd(i,km1)+1. the simpler form is used only for
         !     computational efficiency.
         !----------------------------------------------------------------------

         DO_KK2: do kk = km3, k - ncdm(km1), - 1
            !#TODO: 1st expensive loop
            do i = 1, ni
               if (km2 >= nct(i) .and. cldfrac(i,km1) >= cut .and. &
                    ncd(i,km1) > 2 .and. kk >= k - ncd(i,km1)) then

                  sanu     = anu(i,kk)
                  anutau   = sanu / (sanu + taum(i,2,km1) - taum(i,2,kk))
                  if (sanu <= 0.50) then
                     dtrgw  = sqrt(anutau)
                  else if (sanu > 0.50 .and. sanu <= 1.0) then
                     xx     = sqrt(anutau)
                     dtrgw  = xx + 2.0 * (sanu - 0.50) * (anutau - xx)
                  else if (sanu > 1.0 .and. sanu <= 2.0) then
                     dtrgw  = anutau + (sanu - 1.0) * anutau * (anutau - 1.0)
                  else if (sanu > 2.0 .and. sanu <= 3.0) then
                     xx     = anutau * anutau
                     dtrgw  = xx + (sanu - 2.0) * xx * (anutau - 1.0)
                  else if (sanu > 3.0 .and. sanu <= 4.0) then
                     xx     = anutau * anutau * anutau
                     dtrgw  = xx + (sanu - 3.0) * xx * (anutau - 1.0)
                  else if (sanu > 4.0 .and. sanu <= 20.0) then
                     dtrgw  = anutau ** (nint(sanu))
                  else
                     dtrgw  = exp(-(taum(i,2,km1)-taum(i,2,kk)))
                  endif

                  s1(i)     = s1(i) - (emisw(i,kk) - taucsg(i,kk)) * dtrgw
               endif
            enddo
         enddo DO_KK2
      endif

      do i = 1, ni
         if (km2 >= nct(i) .and. cldfrac(i,km1) >= cut) then
!!$            if (cldfrac(i,km1) >= cut) then
               if (ncd(i,km1) == 1) then
                  fw(i,1,k) = fx(i,2,k)
                  fd(i,2,k) = fx(i,1,k) + cldfrac(i,km1) * &
                       (fx(i,2,k) - fx(i,1,k))
               else
                  fw(i,1,k) = term1(i) + s1(i) + bf(i,k) - &
                       emisw(i,km1) + taucsg(i,km1)
                  fd(i,2,k) = cldm(i,km1) * &
                       (fw(i,1,k) - fx(i,2,k)) + &
                       fx(i,1,k) + cldfrac(i,km1) * &
                       (fx(i,2,k) - fx(i,1,k))
               endif
!!$            endif
         endif
      enddo
   enddo DO_K_750

   return
end subroutine ccc2_lwtran1

end module ccc2_lwtran
