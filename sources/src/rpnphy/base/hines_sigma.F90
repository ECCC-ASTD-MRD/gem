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
!*** S/P HINES_SIGMA
  SUBROUTINE hines_sigma (sigma_t,sigma_alpha,sigsqh_alpha,  &
       &                        naz,lev,il1,il2,nlons,nlevs,nazmth)
    implicit none
#include <arch_specific.hf>

    INTEGER  lev, naz, il1, il2
    INTEGER  nlons, nlevs, nazmth
    REAL*8  sigma_t(nlons,nlevs)
    REAL*8  sigma_alpha(nlons,nlevs,nazmth)
    REAL*8  sigsqh_alpha(nlons,nlevs,nazmth)
!
!Author
!
!  aug. 7/95 - c. mclandress
!
!Object
!
!  This routine calculates the total rms and azimuthal rms horizontal 
!  velocities at a given level on a longitude by altitude grid for 
!  the hines' doppler spread gwd parameterization scheme.
!  note: only four or eight azimuths can be used.
!  output arguements:
!
!Revision
!
!Modules
!
!Arguments
!
!               - Output -
! sigma_t       total rms horizontal wind (m/s).
! sigma_alpha   total rms wind in each azimuth (m/s).
!
!               - Input -
! sigsqh_alpha  portion of wind variance from waves having wave
!               normals in the alpha azimuth (m/s).
! naz           actual number of horizontal azimuths used (must be 4 or 8).
! lev           altitude level to process.
! il1           first longitudinal index to use (il1 >= 1).
! il2           last longitudinal index to use (il1 <= il2 <= nlons).
! nlons         number of longitudes.
! nlevs         number of vertical levels.
! nazmth        azimuthal array dimension (nazmth >= naz).
!
!  subroutine arguements.
!
    !
    !  internal variables.
    !
    INTEGER  i, n
    REAL*8  sum_even, sum_odd, above_zero
    !-----------------------------------------------------------------------     
    !
    !  calculate azimuthal rms velocity for the 4 azimuth case.
    !
    IF (naz.EQ.4)  THEN
       DO i = il1,il2
          above_zero= max(0.d0,sigsqh_alpha(i,lev,1)+sigsqh_alpha(i,lev,3))
          sigma_alpha(i,lev,1) = SQRT(above_zero)
          above_zero= max(0.d0,sigsqh_alpha(i,lev,2)+sigsqh_alpha(i,lev,4))
          sigma_alpha(i,lev,2) = SQRT(above_zero)
          sigma_alpha(i,lev,3) = sigma_alpha(i,lev,1)
          sigma_alpha(i,lev,4) = sigma_alpha(i,lev,2)
       END DO
    END IF
    !
    !  calculate azimuthal rms velocity for the 8 azimuth case.
    !
    IF (naz.EQ.8)  THEN
       DO i = il1,il2
          sum_odd  = ( sigsqh_alpha(i,lev,1) + sigsqh_alpha(i,lev,3)   &
               &             + sigsqh_alpha(i,lev,5) + sigsqh_alpha(i,lev,7) ) / 2.
          sum_even = ( sigsqh_alpha(i,lev,2) + sigsqh_alpha(i,lev,4)   &
               &             + sigsqh_alpha(i,lev,6) + sigsqh_alpha(i,lev,8) ) / 2.
          above_zero= max(0.d0, sigsqh_alpha(i,lev,1)   &
               &                           + sigsqh_alpha(i,lev,5) + sum_even )
          sigma_alpha(i,lev,1) = SQRT(above_zero)
! this next sqrt was indicated to produce FP invalid operation by xlf13
! on at least one point on processor oe-0000-0385 out of 480 in a GY25 grid
          above_zero= max(0.d0, sigsqh_alpha(i,lev,2)   &
               &                           + sigsqh_alpha(i,lev,6) + sum_odd )
          sigma_alpha(i,lev,2) = SQRT(above_zero)
          above_zero= max(0.d0, sigsqh_alpha(i,lev,3)   &
               &                           + sigsqh_alpha(i,lev,7) + sum_even )
          sigma_alpha(i,lev,3) = SQRT(above_zero)
          above_zero= max(0.d0, sigsqh_alpha(i,lev,4)   &
               &                           + sigsqh_alpha(i,lev,8) + sum_odd )
          sigma_alpha(i,lev,4) = SQRT(above_zero)
          sigma_alpha(i,lev,5) = sigma_alpha(i,lev,1)
          sigma_alpha(i,lev,6) = sigma_alpha(i,lev,2)
          sigma_alpha(i,lev,7) = sigma_alpha(i,lev,3)
          sigma_alpha(i,lev,8) = sigma_alpha(i,lev,4)
       END DO
    END IF
    !
    !  calculate total rms velocity.
    !
    DO i = il1,il2
       sigma_t(i,lev) = 0.
    END DO
    DO n = 1,naz
       DO i = il1,il2
          sigma_t(i,lev) = sigma_t(i,lev) + sigsqh_alpha(i,lev,n)
       END DO
    END DO
    DO i = il1,il2
       above_zero= max(0.d0, sigma_t(i,lev))
       sigma_t(i,lev) = SQRT(above_zero)
    END DO

    RETURN
    !-----------------------------------------------------------------------     
  END SUBROUTINE hines_sigma
