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

subroutine hines_heat(heat, diffco,                                 &
     &                alt, bvfreq, density, sigma_t, sigma_alpha,   &
     &                flux, visc_mol, kstar, f1, f2, f3, f5, f6,    &
     &                naz, il1, il2, lev1, lev2, nlons, nlevs,      &
     &                nazmth, losigma_t )
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   integer :: naz, il1, il2, lev1, lev2, nlons, nlevs, nazmth
   real(REAL64) :: kstar, f1, f2, f3, f5, f6
   real(REAL64) :: heat(nlons,nlevs), diffco(nlons,nlevs)
   real(REAL64) :: alt(nlons,nlevs), bvfreq(nlons,nlevs), density(nlons,nlevs)
   real(REAL64) :: sigma_t(nlons,nlevs),  sigma_alpha(nlons,nlevs,nazmth)
   real(REAL64) :: flux(nlons,nlevs,nazmth), visc_mol(nlons,nlevs)
   logical :: losigma_t(nlons,nlevs)

   !@Author
   !  aug. 6/95 - c. mclandress
   !  2001      - m. charron
   !@Object
   !  This routine calculates the gravity wave induced heating and
   !  diffusion coefficient on a longitude by altitude grid for
   !  the hines' doppler spread gravity wave drag parameterization scheme.
   !  This routine can be used for nonzero minimum cutoff wavenumber (m_min)
   !  only in the case of spectral slope=1, in which case m_min is not needed
   !  since its vertical derivative is zero.
   !@arguments
   !                 - output -
   ! heat            gravity wave heating (k/sec).
   ! diffco          diffusion coefficient (m^2/sec)
   !                 - input -
   ! bvfreq          background brunt vassala frequency (rad/sec).
   ! density         background density (kg/m^3).
   ! sigma_t         total rms horizontal wind (m/s).
   ! visc_mol        molecular viscosity (m^2/s).
   ! kstar           typical gravity wave horizontal wavenumber (1/m).
   ! slope           slope of incident vertical wavenumber spectrum.
   ! f1,f2,f3,f5,f6  hines's fudge factors.
   ! il1             first longitudinal index to use (il1 >= 1).
   ! il2             last longitudinal index to use (il1 <= il2 <= nlons).
   ! lev1            first altitude level to use (lev1 >=1).
   ! lev2            last altitude level to use (lev1 < lev2 <= nlevs).
   ! nlons           number of longitudes.
   ! nlevs           number of vertical levels.
   ! nazmth          azimuthal array dimension (nazmth >= naz).
   ! losigma_t       .true. for total sigma not zero


   ! internal variables.

   integer  :: ii,i, l, n, lev1p, lev2m
   real(REAL64) :: m_sub_m_turb, m_sub_m_mol, m_sub_m, dendz2
   real(REAL64) :: heatng1(il1:il2)
   real(REAL64) :: visc, visc_min

   real(REAL64) :: dfdz(nlons,nlevs,nazmth)

   real(REAL64)  :: cpd
   real(REAL64) :: zero
   !-----------------------------------------------------------------------

   zero=0.
   cpd=1004.
   visc_min = 1.e-10

   lev1p = lev1 + 1
   lev2m = lev2 - 1

   do l = lev1p,lev2m
      do i = il1,il2
         if (losigma_t(i,l)) then
            dendz2 = density(i,l) * ( alt(i,l-1) - alt(i,l) )
            visc    = max ( visc_mol(i,l), visc_min )
            m_sub_m_turb = bvfreq(i,l) / ( f2 * sigma_t(i,l) )
            m_sub_m_mol  = (bvfreq(i,l)*kstar/visc)**0.33333333/f3
            m_sub_m      = min ( m_sub_m_turb, m_sub_m_mol )
            do ii=1,8
               !ts          dfdz(i,l,:) = ( flux(i,l-1,:) - flux(i,l,:) ) / dendz2 &
               !ts               & * ( f1*sigma_alpha(i,l,:) + bvfreq(i,l)/m_sub_m )
               dfdz(i,l,ii) = ( flux(i,l-1,ii) - flux(i,l,ii) ) / dendz2 &
                    & * ( f1*sigma_alpha(i,l,ii) + bvfreq(i,l)/m_sub_m )
            enddo
         endif
      end do
   end do

   do i = il1,il2
      if (losigma_t(i,lev1)) then
         dendz2 = density(i,lev1) * ( alt(i,lev1) - alt(i,lev1p) )
         visc    = max ( visc_mol(i,lev1), visc_min )
         m_sub_m_turb = bvfreq(i,lev1) / ( f2 * sigma_t(i,lev1) )
         m_sub_m_mol  = (bvfreq(i,lev1)*kstar/visc)**0.33333333/f3
         m_sub_m      = min ( m_sub_m_turb, m_sub_m_mol )
         do ii=1,8
            !ts       dfdz(i,lev1,:) = -flux(i,lev1,:) / dendz2 &
            !ts            & * ( f1*sigma_alpha(i,lev1,:) + bvfreq(i,lev1)/m_sub_m )
            dfdz(i,lev1,ii) = -flux(i,lev1,ii) / dendz2 &
                 & * ( f1*sigma_alpha(i,lev1,ii) + bvfreq(i,lev1)/m_sub_m )
         enddo
      endif
   end do

   do i = il1,il2
      if (losigma_t(i,lev2)) then
         dendz2 = density(i,lev2) * ( alt(i,lev2m) - alt(i,lev2) )
         visc    = max ( visc_mol(i,lev2), visc_min )
         m_sub_m_turb = bvfreq(i,lev2) / ( f2 * sigma_t(i,lev2) )
         m_sub_m_mol  = (bvfreq(i,lev2)*kstar/visc)**0.33333333/f3
         m_sub_m      = min ( m_sub_m_turb, m_sub_m_mol )
         do ii=1,8
            !ts       dfdz(i,lev2,:) = ( flux(i,lev2m,:) - flux(i,lev2,:) ) / dendz2 &
            !ts            & * ( f1*sigma_alpha(i,lev2,:) + bvfreq(i,lev2)/m_sub_m )
            dfdz(i,lev2,ii) = ( flux(i,lev2m,ii) - flux(i,lev2,ii) ) / dendz2 &
                 & * ( f1*sigma_alpha(i,lev2,ii) + bvfreq(i,lev2)/m_sub_m )
         enddo
      endif
   end do

   !  heating and diffusion.


   !  maximum permissible value of cutoff wavenumber is the smaller
   !  of the instability-induced wavenumber (m_sub_m_turb) and
   !  that imposed by molecular viscosity (m_sub_m_mol).


   do l = lev1,lev2
      heatng1=0
      do n=1,naz
         do i = il1,il2
            if (losigma_t(i,l)) then
               heatng1(i) = heatng1(i) - f5 * dfdz(i,l,n)
            endif
         enddo
      enddo
      do i = il1,il2
         if (losigma_t(i,l)) then
            visc    = max ( visc_mol(i,l), visc_min )
            m_sub_m_turb = bvfreq(i,l) / ( f2 * sigma_t(i,l) )
            m_sub_m_mol  = (bvfreq(i,l)*kstar/visc)**0.33333333/f3
            m_sub_m      = min ( m_sub_m_turb, m_sub_m_mol )
            !ts          heatng = 0.
            !ts          DO n=1,naz
            !ts             heatng = heatng - f5 * dfdz(i,l,n)
            !ts          ENDDO
            diffco(i,l) = f6 * heatng1(i)**0.33333333 / m_sub_m**1.33333333
            ! The turubulent diffusion is limited by the molecular viscosity following
            ! Akmaev (JGR, 2001) and Akmaev et al. (Ann. Geoph., 1997)
            diffco(i,l) = max ( diffco(i,l) - visc_mol(i,l), zero )

            heat(i,l)   = heatng1(i) / cpd
         endif
      end do
   end do

   return
   !-----------------------------------------------------------------------
end subroutine hines_heat
