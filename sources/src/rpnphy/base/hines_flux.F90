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
!*** S/P HINES_FLUX
  SUBROUTINE hines_flux3( flux_u, flux_v, flux, drag_u, drag_v,        & 
       &                  alt, density, densb,                         &
       &                  m_alpha, ak_alpha, k_alpha,                  &
       &                  m_min, slope, naz,                           &
       &                  il1, il2, lev1, lev2, nlons, nlevs, nazmth,  &
       &                  kount, trnch, lorms, hflt )
    implicit none
#include <arch_specific.hf>

    INTEGER  naz, il1, il2, lev1, lev2, hflt
    INTEGER  nlons, nlevs, nazmth, kount, trnch
    REAL*8  slope, m_min
    REAL*8  flux_u(nlons,nlevs), flux_v(nlons,nlevs)
    REAL*8  flux(nlons,nlevs,nazmth)
    REAL*8  drag_u(nlons,nlevs), drag_v(nlons,nlevs)
    REAL*8  alt(nlons,nlevs),    density(nlons,nlevs), densb(nlons)
    REAL*8  m_alpha(nlons,nlevs,nazmth)
    REAL*8  ak_alpha(nlons,nazmth), k_alpha(nlons,nazmth)

    LOGICAL lorms(nlons)
!
!Authors
!  aug. 6/95 - c. mclandress
!       2001 - m. charron
!
!Revision
! 001   L. Spacek (Sep 2008) - calculate drag using centered differences
! 002   M. Charron & A. Zadra (June 2011) - use closest levels for centered differences
!
!Modules
!
!Object
!  Calculate zonal and meridional components of the vertical flux 
!  of horizontal momentum and corresponding wave drag (force per unit mass)
!  on a longitude by altitude grid for the hines' doppler spread 
!  gwd parameterization scheme.
!  note: only 4 or 8 azimuths can be used.
!
!           - Output -
! flux_u    zonal component of vertical momentum flux (pascals)
! flux_v    meridional component of vertical momentum flux (pascals)
! drag_u    zonal component of drag (m/s^2).
! drag_v    meridional component of drag (m/s^2).
!
!            - Input -
! alt        altitudes (m).
! density    background density (kg/m^3).
! densb      background density at bottom level (kg/m^3).
! m_alpha    cutoff vertical wavenumber (1/m).
! ak_alpha   spectral amplitude factor (i.e., {ajkj} in m^4/s^2).
! k_alpha    horizontal wavenumber (1/m).
! slope      slope of incident vertical wavenumber spectrum.
! m_min      minimum allowable cutoff wavenumber (1/m)
!            for spectral slope of one.
! naz        actual number of horizontal azimuths used (must be 4 or 8).
! il1        first longitudinal index to use (il1 >= 1).
! il2        last longitudinal index to use (il1 <= il2 <= nlons).
! lev1       first altitude level to use (lev1 >=1). 
! lev2       last altitude level to use (lev1 < lev2 <= nlevs).
! nlons      number of longitudes.
! nlevs      number of vertical levels.
! nazmth     azimuthal array dimension (nazmth >= naz).
! lorms      .true. for drag computation (column selector)
!
!  constant in data statement.
!
! cos45      cosine of 45 degrees. 		
!

    !
    !  internal variables.
    !
    INTEGER  i, l, lev1p, lev2p, lev2m, k, it
    REAL*8  cos45, dendz, dendz2
    REAL :: work_u(nlons,nlevs), work_v(nlons,nlevs)
    !-----------------------------------------------------------------------
    cos45 = 0.7071068   
    !
    lev1p = lev1 + 1
    lev2m = lev2 - 1
    lev2p = lev2 + 1
    !
    !  sum over azimuths for case where slope = 1.
    !

    IF ( ABS(slope-1.) .LT. EPSILON(1.) )  THEN
       !
       !  case with 4 azimuths.
       !
       IF (naz.EQ.4)  THEN
          DO l = lev1,lev2
             DO i = il1,il2
                flux(i,l,:) = ak_alpha(i,:)*k_alpha(i,:)*(m_alpha(i,l,:)-m_min)
                flux_u(i,l) = flux(i,l,1) - flux(i,l,3)
                flux_v(i,l) = flux(i,l,2) - flux(i,l,4)
             END DO
          END DO
       END IF
       !
       !  case with 8 azimuths.
       !
       IF (naz.EQ.8)  THEN
          DO l = lev1,lev2
             DO k = 1, nazmth
                DO i = il1,il2
                  flux(i,l,k) = ak_alpha(i,k)*k_alpha(i,k)*(m_alpha(i,l,k)-m_min)
                END DO
             END DO
             DO i = il1,il2
                flux_u(i,l) = flux(i,l,1) - flux(i,l,5) + cos45 *     &
                     ( flux(i,l,2) - flux(i,l,4) - flux(i,l,6) + flux(i,l,8) )
                flux_v(i,l) = flux(i,l,3) - flux(i,l,7) + cos45 *     &
                     ( flux(i,l,2) + flux(i,l,4) - flux(i,l,6) - flux(i,l,8) )
             END DO
          END DO
       END IF

    END IF
    !
    !  sum over azimuths for case where slope not equal to 1.
    !
    IF ( ABS(slope-1.) .GT. EPSILON(1.) )  THEN
       !
       !  case with 4 azimuths.
       !
       IF (naz.EQ.4)  THEN
          DO l = lev1,lev2
             DO i = il1,il2
                flux(i,l,:) = ak_alpha(i,:)*k_alpha(i,:)*m_alpha(i,l,:)**slope
                flux_u(i,l) = flux(i,l,1) - flux(i,l,3)
                flux_v(i,l) = flux(i,l,2) - flux(i,l,4)
             END DO
          END DO
       END IF
       !
       !  case with 8 azimuths.
       !
       IF (naz.EQ.8)  THEN
          DO l = lev1,lev2
             DO k = 1, nazmth
                DO i = il1,il2
                  flux(i,l,k) = ak_alpha(i,k)*k_alpha(i,k)*m_alpha(i,l,k)**slope
                END DO
             END DO
             DO i = il1,il2
                flux_u(i,l) = flux(i,l,1) - flux(i,l,5) + cos45 *     &
                     ( flux(i,l,2) - flux(i,l,4) - flux(i,l,6) + flux(i,l,8) )
                flux_v(i,l) = flux(i,l,3) - flux(i,l,7) + cos45 *     &
                     ( flux(i,l,2) + flux(i,l,4) - flux(i,l,6) - flux(i,l,8) )
             END DO
          END DO
       END IF

    END IF
    !
    !  calculate flux from sum.
    !
    DO l = lev1,lev2
       DO i = il1,il2
          flux_u(i,l) = flux_u(i,l) * densb(i) / slope
          flux_v(i,l) = flux_v(i,l) * densb(i) / slope
       END DO
       DO k = 1, nazmth
          DO i = il1,il2
            flux(i,l,k) = flux(i,l,k) * densb(i) / slope
          END DO
       END DO
    END DO
    !
    !  filter fluxes (hflt iterations)
    !
    work_u = 0.
    work_v = 0.
    
    DO it=1,hflt
       DO i = il1,il2
          work_u(i,lev1) = 0.25*(2.*flux_u(i,lev1  ) + &
                                    flux_u(i,lev1+1) )
          work_v(i,lev1) = 0.25*(2.*flux_v(i,lev1  ) + &
                                    flux_v(i,lev1+1) )   
          work_u(i,lev2) = 0.25*(   flux_u(i,lev2-1) + &
                                 2.*flux_u(i,lev2  ) )
          work_v(i,lev2) = 0.25*(   flux_v(i,lev2-1) + &
                                 2.*flux_v(i,lev2  ) )                                     
       END DO
       DO l = lev1+1,lev2-1
          DO i = il1,il2     
             work_u(i,l) = 0.25*(   flux_u(i,l-1) + &
                                 2.*flux_u(i,l  ) + &
                                    flux_u(i,l+1) )
             work_v(i,l) = 0.25*(   flux_v(i,l-1) + &
                                 2.*flux_v(i,l  ) + &
                                    flux_v(i,l+1) )
          END DO
       END DO
       DO l = lev1,lev2
          DO i = il1,il2     
             flux_u(i,l) = work_u(i,l) 
             flux_v(i,l) = work_v(i,l) 
          END DO
       END DO      
    END DO
    !
    !  calculate drag at full levels using centered differences
    !      
    DO l = lev1p,lev2m
       DO i = il1,il2
          IF (lorms(i)) THEN
             dendz2 = density(i,l) * ( alt(i,l-1) - alt(i,l) ) 
             drag_u(i,l) = - ( flux_u(i,l-1) - flux_u(i,l) ) / dendz2
             drag_v(i,l) = - ( flux_v(i,l-1) - flux_v(i,l) ) / dendz2
         ENDIF
       END DO
    END DO

    !
    !  drag at first and last levels using one-side differences.
    ! 
    DO i = il1,il2
       IF (lorms(i)) THEN
          dendz = density(i,lev1) * ( alt(i,lev1) - alt(i,lev1p) ) 
          drag_u(i,lev1) =  flux_u(i,lev1)  / dendz
          drag_v(i,lev1) =  flux_v(i,lev1)  / dendz
       ENDIF
    END DO

    DO i = il1,il2
       IF (lorms(i)) THEN
          dendz = density(i,lev2) * ( alt(i,lev2m) - alt(i,lev2) )
          drag_u(i,lev2) = - ( flux_u(i,lev2m) - flux_u(i,lev2) ) / dendz
          drag_v(i,lev2) = - ( flux_v(i,lev2m) - flux_v(i,lev2) ) / dendz
       ENDIF
    END DO
    IF (nlevs .GT. lev2) THEN
       DO i = il1,il2
          IF (lorms(i)) THEN
             dendz = density(i,lev2p) * ( alt(i,lev2) - alt(i,lev2p) )
             drag_u(i,lev2p) = -  flux_u(i,lev2)  / dendz
             drag_v(i,lev2p) = - flux_v(i,lev2)  / dendz
          ENDIF
       END DO
    ENDIF

    RETURN
    !-----------------------------------------------------------------------
 END SUBROUTINE hines_flux3
