!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END --------------------------
!/@*
SUBROUTINE hines_extro4( nlons, nlevs, nazmth,                          &
                           drag_u, drag_v, heat, diffco, flux_u, flux_v,  &
                           vel_u, vel_v, bvfreq, density, visc_mol, alt,  &
                           rmswind, anis, k_alpha, sigsqmcw,              &
                           m_alpha,  mmin_alpha, sigma_t, sigmatm,        &
                           kount, trnch,                                  &
                           lev2, hflt, lorms, iheatcal)

    USE mo_gwspectrum,     ONLY: kstar, m_min,            &
         &                       naz, slope, f1, f2, f3, f5, f6,  &
         &                       icutoff, alt_cutoff, smco, nsmax
    use phy_status, only: phy_error_L
    implicit none
#include <arch_specific.hf>

    INTEGER :: nlons, nlevs, nazmth, lev2, hflt, kount, trnch

    REAL*8    :: drag_u(nlons,nlevs),   drag_v(nlons,nlevs) 
    REAL*8    :: heat(nlons,nlevs),     diffco(nlons,nlevs)
    REAL*8    :: flux_u(nlons,nlevs),   flux_v(nlons,nlevs)
    REAL*8    :: flux(nlons,nlevs,nazmth)
    REAL*8    :: vel_u(nlons,nlevs),    vel_v(nlons,nlevs)
    REAL*8    :: bvfreq(nlons,nlevs),   density(nlons,nlevs)
    REAL*8    :: visc_mol(nlons,nlevs), alt(nlons,nlevs)
    REAL*8    :: rmswind(nlons),      bvfb(nlons),   densb(nlons)
    REAL*8    :: anis(nlons,nazmth)
    REAL*8    :: sigma_t(nlons,nlevs), sigsqmcw(nlons,nlevs,nazmth)
    REAL*8    :: sigma_alpha(nlons,nlevs,nazmth), sigmatm(nlons,nlevs)

    REAL*8    :: m_alpha(nlons,nlevs,nazmth), v_alpha(nlons,nlevs,nazmth)
    REAL*8    :: ak_alpha(nlons,nazmth),      k_alpha(nlons,nazmth)
    REAL*8    :: mmin_alpha(nlons,nazmth)    
    REAL*8    :: smoothr1(nlons,nlevs), smoothr2(nlons,nlevs)

    LOGICAL :: lorms(nlons), losigma_t(nlons,nlevs)

!@Authors
!  aug. 13/95 - c. mclandress
!  sept. /95  - n. mcfarlane
!  1995- 2002 - e. manzini
!@Revision
!  001   A. Zadra & R. McTaggart-Cowan (June 2011) - provide upper boundary index through interface
!  001   A. Plante - Remove toplev (notop). Set lev1 to 1 since our model to is well below 120km.
!@Object
!  main routine for hines' "extrowave" gravity wave parameterization based
!  on hines' doppler spread theory. this routine calculates zonal
!  and meridional components of gravity wave drag, heating rates
!  and diffusion coefficient on a longitude by altitude grid.
!  no "mythical" lower boundary region calculation is made. 
!@Argumentsxs
!              - Output -
! drag_u       zonal component of gravity wave drag (m/s^2).
! drag_v       meridional component of gravity wave drag (m/s^2).
! heat         gravity wave heating (k/sec).
! diffco       diffusion coefficient (m^2/sec)
! flux_u       zonal component of vertical momentum flux (pascals)
! flux_v       meridional component of vertical momentum flux (pascals)
!              - Input -
! vel_u        background zonal wind component (m/s).
! vel_v        background meridional wind component (m/s).
! bvfreq       background brunt vassala frequency (radians/sec).
! densit       background density (kg/m^3) 
! visc_mol     molecular viscosity (m^2/s)
! alt          altitude of momentum, density, buoyancy levels (m)
!              (note: levels ordered so that alt(i,1) > alt(i,2), etc.)
! rmswind      root mean square gravity wave wind at lowest level (m/s).
! anis         anisotropy factor (sum over azimuths is one)
! lorms        .true. for drag computation (column selector)
! k_alpha      horizontal wavenumber of each azimuth (1/m).
! lev2         index of last level (eg bottom) for drag calculation 
!              (i.e., lev1 < lev2 <= nlevs).
! nlons        number of longitudes.
! nlevs        number of vertical levels.
! nazmth       azimuthal array dimension (nazmth >= naz).
!              - ouput diagnostics -
! m_alpha      cutoff vertical wavenumber (1/m).
! mmin_alpha   minimum value of cutoff wavenumber.
! sigma_t      total rms horizontal wind (m/s).
!              - work arrays -
! v_alpha      wind component at each azimuth (m/s) and if iheatcal=1
!              holds vertical derivative of cutoff wavenumber.
! sigma_alpha  total rms wind in each azimuth (m/s).
! ak_alpha     spectral amplitude factor at each azimuth 
!              (i.e.,{ajkj}) in m^4/s^2.
! densb        background density at bottom level.
! bvfb         buoyancy frequency at bottom level and
!              work array for icutoff = 1.
! losigma_t    .true. for total sigma not zero
!*@/

    !
    !  internal variables.
    !
    INTEGER :: i, n, l, lev1, il1, il2, iprint, iheatcal

    !----------------------------------------------------------------------- 
    !

    ! range of longitude index:
    il1 = 1      
    il2 = nlons     

    lev1=1              ! top level index

    iprint = 0       !     * iprint     = 1 to print out various arrays.

    !
    !  buoyancy and density at bottom level.
    !
    DO i = il1,il2
       bvfb(i)  = bvfreq(i,lev2)
       densb(i) = density(i,lev2)
    END DO
    !
    !  initialize some variables
    !
    DO n = 1,naz
       DO l=1,lev2
          DO i=il1,il2
             m_alpha(i,l,n) =  m_min
          END DO
       END DO
    END DO
    !
    !  compute azimuthal wind components from zonal and meridional winds.
    !

    CALL hines_wind ( v_alpha,   & 
         &                  vel_u, vel_v, naz,   &
         &                  il1, il2, lev1, lev2, nlons, nlevs, nazmth )

    !  calculate cutoff vertical wavenumber and velocity variances.
    !
    CALL hines_wavnum ( m_alpha, sigma_t, sigma_alpha, ak_alpha,   &
         &              mmin_alpha, losigma_t,                     &
         &              v_alpha, visc_mol, density, densb,         &
         &              bvfreq, bvfb, rmswind, anis, lorms,        &
         &              sigsqmcw, sigmatm,                         &
         &              il1, il2, lev1, lev2, nlons, nlevs, nazmth)
    if (phy_error_L) return

    !  smooth cutoff wavenumbers and total rms velocity in the vertical 
    !  direction nsmax times, using flux_u as temporary work array.
    !   
    IF (nsmax.GT.0)  THEN
       DO n = 1,naz
          DO l=1,lev2
             DO i=il1,il2
                smoothr1(i,l) = m_alpha(i,l,n)
             END DO
          END DO
          CALL vert_smooth (smoothr1,smoothr2, smco, nsmax,  &
               &                       il1, il2, lev1, lev2, nlons, nlevs )
          DO l=lev1,lev2
             DO i=il1,il2
                m_alpha(i,l,n) = smoothr1(i,l)
             END DO
          END DO
       END DO
       CALL vert_smooth ( sigma_t, smoothr2, smco, nsmax, &
            &                     il1, il2, lev1, lev2, nlons, nlevs )
    END IF
    !
    !  calculate zonal and meridional components of the
    !  momentum flux and drag.
    !
    CALL hines_flux3( flux_u, flux_v, flux, drag_u, drag_v,       &
         &            alt, density, densb,                        &
         &            m_alpha,  ak_alpha, k_alpha,                &
         &            m_min, slope, naz,                          &
         &            il1, il2, lev1, lev2, nlons, nlevs, nazmth, &
         &            kount, trnch, lorms, hflt )

    !  cutoff drag above alt_cutoff, using bvfb as temporary work array.
    !
    IF (icutoff.EQ.1)  THEN
       CALL hines_exp ( drag_u, bvfb, alt, alt_cutoff,  &
            &                   il1, il2, lev1, lev2, nlons, nlevs )
       CALL hines_exp ( drag_v, bvfb, alt, alt_cutoff,  &
            &                   il1, il2, lev1, lev2, nlons, nlevs )
    END IF
    if (phy_error_L) return

    !  print out various arrays for diagnostic purposes.
    !
    IF (iprint.EQ.1)  THEN
       CALL hines_print ( flux_u, flux_v, drag_u, drag_v, alt,     &
            &             sigma_t, sigma_alpha, v_alpha, m_alpha,  &
            &             1, 1, il1, il2, lev1, lev2,              &
            &             naz, nlons, nlevs, nazmth)
    END IF
    !
    !  if not calculating heating rate and diffusion coefficient then finished.
    !
    IF (iheatcal.NE.1)  RETURN

    call physeterror('hines_extro4', 'You are not supposed to call hines_heat')
    return

    !TODO: WTF: this would never be called!
    !  heating rate and diffusion coefficient.
    CALL hines_heat ( heat, diffco,                                 &
         &            alt, bvfreq, density, sigma_t, sigma_alpha,   &
         &            flux, visc_mol, kstar, f1, f2, f3, f5, f6,    &
         &            naz, il1, il2, lev1, lev2, nlons, nlevs,      &
         &            nazmth, losigma_t )

    RETURN
    !-----------------------------------------------------------------------
 END SUBROUTINE hines_extro4
