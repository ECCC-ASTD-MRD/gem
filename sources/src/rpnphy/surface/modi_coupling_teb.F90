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
!   ########################
    MODULE MODI_COUPLING_TEB2
!   ########################
!
INTERFACE
!
SUBROUTINE COUPLING_TEB2(PTSTEP, KYEAR, KMONTH, KDAY, PTIME, PTSUN, PZENITH, PAZIM,         &
               PZREF, PUREF, PZS, PU, PV, PQA, PTA, PRHOA, PRAIN, PSNOW, PLW, PDIR_SW,     &
               PSCA_SW, PSW_BANDS, PPS, PPA,                                               &
               PSFTQ, PSFTH, PSFU, PSFV, PTRAD, PDIR_ALB, PSCA_ALB, PEMIS, PLAT            )
!
!*      0.1    declarations of arguments

INTEGER,            INTENT(IN)  :: KYEAR     ! current year (UTC)
INTEGER,            INTENT(IN)  :: KMONTH    ! current month (UTC)
INTEGER,            INTENT(IN)  :: KDAY      ! current day (UTC)
REAL,               INTENT(IN)  :: PTIME     ! current time since midnight (UTC, s)
REAL, DIMENSION(:), INTENT(IN)  :: PTSUN     ! solar time                    (s from midnight)
REAL,               INTENT(IN)  :: PTSTEP    ! atmospheric time-step                 (s)
REAL, DIMENSION(:), INTENT(IN)  :: PZREF     ! height of T,q forcing                 (m)
REAL, DIMENSION(:), INTENT(IN)  :: PUREF     ! height of wind forcing                (m)
!
REAL, DIMENSION(:), INTENT(IN)  :: PTA       ! air temperature forcing               (K)
REAL, DIMENSION(:), INTENT(IN)  :: PQA       ! air humidity forcing                  (kg/m3)
REAL, DIMENSION(:), INTENT(IN)  :: PRHOA     ! air density                           (kg/m3)
REAL, DIMENSION(:), INTENT(IN)  :: PU        ! zonal wind                            (m/s)
REAL, DIMENSION(:), INTENT(IN)  :: PV        ! meridian wind                         (m/s)
REAL, DIMENSION(:,:),INTENT(IN) :: PDIR_SW   ! direct  solar radiation (on horizontal surf.)
!                                            !                                       (W/m2)
REAL, DIMENSION(:,:),INTENT(IN) :: PSCA_SW   ! diffuse solar radiation (on horizontal surf.)
!                                            !                                       (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PSW_BANDS ! mean wavelength of each shortwave band (m)
REAL, DIMENSION(:), INTENT(IN)  :: PZENITH   ! zenithal angle       (radian from the vertical)
REAL, DIMENSION(:), INTENT(IN)  :: PAZIM     ! azimuthal angle      (radian from North, clockwise)
REAL, DIMENSION(:), INTENT(IN)  :: PLW       ! longwave radiation (on horizontal surf.)
!                                            !                                       (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PPS       ! pressure at atmospheric model surface (Pa)
REAL, DIMENSION(:), INTENT(IN)  :: PPA       ! pressure at forcing level             (Pa)
REAL, DIMENSION(:), INTENT(IN)  :: PZS       ! atmospheric model orography           (m)
REAL, DIMENSION(:), INTENT(IN)  :: PSNOW     ! snow precipitation                    (kg/m2/s)
REAL, DIMENSION(:), INTENT(IN)  :: PRAIN     ! liquid precipitation                  (kg/m2/s)
!
REAL, DIMENSION(:), INTENT(OUT) :: PSFTH     ! flux of heat                          (W/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PSFTQ     ! flux of water vapor                   (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT) :: PSFU      ! zonal momentum flux                   (Pa)
REAL, DIMENSION(:), INTENT(OUT) :: PSFV      ! meridian momentum flux                (Pa)
!
REAL, DIMENSION(:), INTENT(OUT) :: PTRAD     ! radiative temperature                 (K)
REAL, DIMENSION(:,:),INTENT(OUT):: PDIR_ALB  ! direct albedo for each spectral band  (-)
REAL, DIMENSION(:,:),INTENT(OUT):: PSCA_ALB  ! diffuse albedo for each spectral band (-)
REAL, DIMENSION(:), INTENT(OUT) :: PEMIS     ! emissivity                            (-)
REAL, DIMENSION(:), INTENT(IN)  :: PLAT      ! latitude                              (deg)
!
END SUBROUTINE COUPLING_TEB2
!
!
END INTERFACE
!
!
END MODULE MODI_COUPLING_TEB2
