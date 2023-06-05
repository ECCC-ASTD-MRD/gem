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

subroutine ccc1_lwtragh(fu, fd, slwf, tauci, omci, &
     taual, taug, bf, urbf, cldfrac, &
     em0, bs, cut, il1, il2, &
     ilg, lay, lev)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer ilg, lay, lev, il1, il2, i, k, km1, km2, kp1
   real ru, ubeta, epsd, epsu, taul2, cow, cut
   real ctaul2, crtaul2
   real fu(ilg,2,lev), fd(ilg,2,lev)
   real slwf(ilg), tauci(ilg,lay), omci(ilg,lay), &
        taual(ilg,lay), taug(ilg,lay), bf(ilg,lev), urbf(ilg,lay), &
        cldfrac(ilg,lay), em0(ilg), bs(ilg), &
        taul1(ilg,lay), rtaul1(ilg,lay)
   real xu(ilg,2,lay), xd(ilg,2,lay), dtr(ilg,2,lay), fy(ilg,2,lev), &
        fx(ilg,2,lev)
   real(REAL64) :: dtr_vs(ilg,lay)

   !@Authors
   !        J. Li, M. Lazare, CCCMA, rt code for gcm4
   !        (Ref: J. Li, H. W. Barker, 2005:
   !        JAS Vol. 62, no. 2, pp. 286\226309)
   !        P. Vaillancourt, D. Talbot, RPN/CMC;
   !        adapted for CMC/RPN physics (May 2006)
   !@Revisions
   ! 001    P. Vaillancourt, K. Winger (Sep 2006) :
   !        remplace 1-cld-eps par max(1-cld,eps)
   ! 002    M. Desgagne, P. Vaillancourt, M. Lazarre (Dec 2008) :
   !        Give appropriate dimensions to xu,xd,dtr,fy and fx;
   !        Initialize fx(I,1,1) and fx(I,2,1) correctly
   !@Object  OPTICAL DEPTHS CALCULATION
   !        In the g space with interval close 1 (very large optical depth)
   !        or in the case with cloud absorption is very small or the weight
   !        of flux and cooling rate are very small. The cloud radiative
   !        process can be highly simplified. The absorption approximation
   !        method is used and cloud random and maximum overlap is
   !        considered, but cloud scattering and inhomogeneity are ignored.
   !        The exponential source planck function is used which is more
   !        accurate in the region above 200 mb in comparison with linear
   !        source function.
   !@Arguments
   ! fu       upward infrared flux
   ! fd       downward infrared flux
   ! slwf     input solar flux at model top level for each band
   ! tauci    cloud optical depth for the infrared
   ! omci     cloud single scattering albedo times optical depth
   ! taual    aerosol optical depth for the infrared
   ! taug     gaseous optical depth for the infrared
   ! bf       blackbody intensity integrated over each band at each
   !          level in units w / m^2 / sr.
   ! bs       the blackbody intensity at the surface.
   ! urbf     u times the difference of log(bf) for two neighbor levels
   !          used for exponential source function (li, 2002 jas p3302)
   ! cldfrac  cloud fraction
   ! em0      surface emission
   ! xu       the emission part in the upward flux transmission
   !          (Li, 2002 JAS p3302)
   ! xd       the emission part in the downward flux transmission
   ! dtr      direct transmission
   ! fy       upward flux for pure clear portion (1) and pure cloud  portion (2)
   ! fx       the same as fy but for the downward flux
   !----------------------------------------------------------------------

   data  ru / 1.6487213 /

   !----------------------------------------------------------------------
   !     initialization for first layer. calculate the downward flux in
   !     the second layer
   !     combine the optical properties for the infrared,
   !     1, aerosol + gas; 2, cloud + aerosol + gas.
   !     fd (fu) is down (upward) flux
   !     the overlap between solar and infrared in 4 - 10 um is
   !     considered, slwf is the incoming solar flux
   !     singularity for xd and xu has been considered as li jas 2002
   !----------------------------------------------------------------------

   do k = 2, lev
      km1 = k - 1
      do i = il1, il2
         taul1(i,km1)    =  taual(i,km1) + taug(i,km1)
         rtaul1(i,km1)   =  taul1(i,km1) * ru
         dtr_vs(i,km1)   =  exp(dble(- rtaul1(i,km1)))
      enddo
   enddo

   DO100: do i = il1, il2
      fd(i,1,1)         =  slwf(i)
      fd(i,2,1)         =  slwf(i)
      fx(i,1,1)         =  slwf(i)
      fx(i,2,1)         =  slwf(i)

      dtr(i,1,1)        =  dtr_vs(i,1)
      ubeta             =  urbf(i,1) / (taul1(i,1) + 1.e-20)
      epsd              =  ubeta + 1.0
      epsu              =  ubeta - 1.0

      if (abs(epsd) .gt. 0.001)                                   then
         xd(i,1,1)       = (bf(i,2) - bf(i,1) * dtr(i,1,1)) / epsd
      else
         xd(i,1,1)       =  rtaul1(i,1) * bf(i,1) * dtr(i,1,1)
      endif
      if (abs(epsu) .gt. 0.001)                                   then
         xu(i,1,1)       = (bf(i,2) * dtr(i,1,1) - bf(i,1)) / epsu
      else
         xu(i,1,1)       =  rtaul1(i,1) * bf(i,2) * dtr(i,1,1)
      endif

      fd(i,1,2)         =  fd(i,1,1) * dtr(i,1,1) + xd(i,1,1)

      if (cldfrac(i,1) .lt. cut)                                  then
         fx(i,1,2)       =  fd(i,1,2)
         fx(i,2,2)       =  fd(i,1,2)
         fd(i,2,2)       =  fd(i,1,2)
      else
         taul2           =  tauci(i,1) + taul1(i,1)
         cow             =  1.0 - omci(i,1) / taul2
         ctaul2          =  cow * taul2
         crtaul2         =  ctaul2 * ru
         dtr(i,2,1)      =  exp (- crtaul2)
         ubeta           =  urbf(i,1) / (ctaul2)
         epsd            =  ubeta + 1.0
         epsu            =  ubeta - 1.0

         if (abs(epsd) .gt. 0.001)                                 then
            xd(i,2,1)     = (bf(i,2) - bf(i,1) * dtr(i,2,1)) / epsd
         else
            xd(i,2,1)     =  crtaul2 * bf(i,1) * dtr(i,2,1)
         endif
         if (abs(epsu) .gt. 0.001)                                 then
            xu(i,2,1)     = (bf(i,2) * dtr(i,2,1) - bf(i,1)) / epsu
         else
            xu(i,2,1)     =  crtaul2 * bf(i,2) * dtr(i,2,1)
         endif

         fx(i,1,2)       =  fx(i,1,1) * dtr(i,1,1) + xd(i,1,1)
         fx(i,2,2)       =  fx(i,2,1) * dtr(i,2,1) + xd(i,2,1)
         fd(i,2,2)       =  fx(i,1,2) + &
              cldfrac(i,1) * (fx(i,2,2) - fx(i,1,2))
      endif
   enddo DO100

   DO250: do k = 3, lev
      km1 = k - 1
      km2 = km1 - 1
      do i = il1, il2
         dtr(i,1,km1)    =  dtr_vs(i,km1)
         ubeta           =  urbf(i,km1) / (taul1(i,km1) + 1.e-20)
         epsd            =  ubeta + 1.0
         epsu            =  ubeta - 1.0

         if (abs(epsd) .gt. 0.001)                                 then
            xd(i,1,km1)   = (bf(i,k) - bf(i,km1) * dtr(i,1,km1)) / epsd
         else
            xd(i,1,km1)   =  rtaul1(i,km1) * bf(i,km1) * dtr(i,1,km1)
         endif
         if (abs(epsu) .gt. 0.001)                                 then
            xu(i,1,km1)   = (bf(i,k) * dtr(i,1,km1) - bf(i,km1)) / epsu
         else
            xu(i,1,km1)   =  rtaul1(i,km1) * bf(i,k) * dtr(i,1,km1)
         endif

         fd(i,1,k)       =  fd(i,1,km1) * dtr(i,1,km1) + xd(i,1,km1)

         if (cldfrac(i,km1) .lt. cut)                              then
            fd(i,2,k)     =  fd(i,2,km1) * dtr(i,1,km1) + xd(i,1,km1)
            fx(i,1,k)     =  fd(i,2,k)
            fx(i,2,k)     =  fd(i,2,k)
         else
            taul2         =  tauci(i,km1) + taul1(i,km1)
            cow           =  1.0 - omci(i,km1) / taul2
            ctaul2        =  cow * taul2
            crtaul2       =  ctaul2 * ru
            dtr(i,2,km1)  =  exp (- crtaul2)
            ubeta         =  urbf(i,km1) / (ctaul2)
            epsd          =  ubeta + 1.0
            epsu          =  ubeta - 1.0

            if (abs(epsd) .gt. 0.001)                               then
               xd(i,2,km1) = (bf(i,k) - bf(i,km1) * dtr(i,2,km1)) / epsd
            else
               xd(i,2,km1) =  crtaul2 * bf(i,km1) * dtr(i,2,km1)
            endif
            if (abs(epsu) .gt. 0.001)                               then
               xu(i,2,km1) = (bf(i,k) * dtr(i,2,km1) - bf(i,km1)) / epsu
            else
               xu(i,2,km1) =  crtaul2 * bf(i,k) * dtr(i,2,km1)
            endif

            if (cldfrac(i,km1) .le. cldfrac(i,km2))                 then
               fx(i,1,k)   = ( fx(i,2,km1) + (1.0 - cldfrac(i,km2)) / &
                    max(1.0 - cldfrac(i,km1),1.e-10) * &
                    (fx(i,1,km1) - fx(i,2,km1)) ) * &
                    dtr(i,1,km1) + xd(i,1,km1)
               fx(i,2,k)   =  fx(i,2,km1) * dtr(i,2,km1) + xd(i,2,km1)
            else if (cldfrac(i,km1) .gt. cldfrac(i,km2))            then
               fx(i,1,k)   =  fx(i,1,km1) * dtr(i,1,km1) + xd(i,1,km1)
               fx(i,2,k)   = (fx(i,1,km1)+cldfrac(i,km2)/cldfrac(i,km1) * &
                    (fx(i,2,km1) - fx(i,1,km1))) * &
                    dtr(i,2,km1) + xd(i,2,km1)
            endif

            fd(i,2,k)     =  fx(i,1,k) + cldfrac(i,km1) * (fx(i,2,k) - &
                 fx(i,1,k))
         endif
      enddo
   enddo DO250

   do i = il1, il2
      fu(i,1,lev)      =  fd(i,1,lev) + em0(i) * (bs(i) - fd(i,1,lev))
      fy(i,1,lev)      =  fx(i,1,lev) + em0(i) * (bs(i) - fx(i,1,lev))
      fy(i,2,lev)      =  fx(i,2,lev) + em0(i) * (bs(i) - fx(i,2,lev))

      if (cldfrac(i,lay) .gt. cut)                                then
         fu(i,2,lev)    =  fy(i,1,lev) + &
              cldfrac(i,lay) * (fy(i,2,lev) - fy(i,1,lev))
      else
         fu(i,2,lev)    =  fy(i,2,lev)
      endif

      fu(i,1,lay)      =  fu(i,1,lev) * dtr(i,1,lay) + xu(i,1,lay)

      if (cldfrac(i,lay) .lt. cut)                                then
         fu(i,2,lay)    =  fu(i,2,lev) * dtr(i,1,lay) + xu(i,1,lay)
         fy(i,1,lay)    =  fu(i,2,lay)
         fy(i,2,lay)    =  fu(i,2,lay)
      else
         fy(i,1,lay)    =  fy(i,1,lev) * dtr(i,1,lay) + xu(i,1,lay)
         fy(i,2,lay)    =  fy(i,2,lev) * dtr(i,2,lay) + xu(i,2,lay)
         fu(i,2,lay)    =  fy(i,1,lay) + &
              cldfrac(i,lay) * (fy(i,2,lev) - fy(i,1,lev))
      endif
   enddo

   do k = lev - 2, 1, - 1
      kp1 = k + 1
      do i = il1, il2
         fu(i,1,k)      =  fu(i,1,kp1) * dtr(i,1,k) + xu(i,1,k)

         if (cldfrac(i,k) .lt. cut)                                then
            fu(i,2,k)    =  fu(i,2,kp1) * dtr(i,1,k) + xu(i,1,k)
            fy(i,1,k)    =  fu(i,2,k)
            fy(i,2,k)    =  fu(i,2,k)
         else
            if (cldfrac(i,k) .lt. cldfrac(i,kp1))                   then
               fy(i,1,k)  = ( fy(i,2,kp1) + (1.0 - cldfrac(i,kp1)) / &
                    (1.0 - cldfrac(i,k)) * (fy(i,1,kp1) - &
                    fy(i,2,kp1)) ) * dtr(i,1,k) + xu(i,1,k)
               fy(i,2,k)  =  fy(i,2,kp1) * dtr(i,2,k) + xu(i,2,k)
            else
               fy(i,1,k)  =  fy(i,1,kp1) * dtr(i,1,k) + xu(i,1,k)
               fy(i,2,k)  = ( fy(i,1,kp1) + cldfrac(i,kp1)/cldfrac(i,k) * &
                    (fy(i,2,kp1) - fy(i,1,kp1)) ) * dtr(i,2,k) + &
                    xu(i,2,k)
            endif

            fu(i,2,k)    =  fy(i,1,k) + &
                 cldfrac(i,k) * (fy(i,2,k) - fy(i,1,k))
         endif
      enddo
   enddo

   return
end subroutine ccc1_lwtragh
