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
    MODULE MODI_SNOW_COVER_1LAYER
!
!
!
INTERFACE
!
!
    SUBROUTINE SNOW_COVER_1LAYER(PTSTEP, PANSMIN, PANSMAX, PTODRY,         &
                                 PRHOSMIN, PRHOSMAX, PRHOFOLD, OALL_MELT,  &
                                 PDRAIN_TIME, PWCRN, PZ0SN, PZ0HSN,        &
                                 PTSNOW, PASNOW, PRSNOW, PWSNOW, PTS_SNOW, &
                                 PESNOW,                                   &
                                 PTG, PABS_SW, PLW1, PLW2,                 &
                                 PTA, PQA, PVMOD, PPS, PRHOA, PSR,         &
                                 PZREF, PUREF,                             &
                                 PRNSNOW, PHSNOW, PLESNOW, PGSNOW, PMELT   )
!
!
!*      0.1    declarations of arguments 
!
!
REAL,                 INTENT(IN)    :: PTSTEP   ! time step
REAL,                 INTENT(IN)    :: PANSMIN  ! minimum snow albedo
REAL,                 INTENT(IN)    :: PANSMAX  ! maximum snow albedo
REAL,                 INTENT(IN)    :: PTODRY   ! snow albedo decreasing constant
REAL,                 INTENT(IN)    :: PRHOSMIN ! minimum snow density
REAL,                 INTENT(IN)    :: PRHOSMAX ! maximum snow density
REAL,                 INTENT(IN)    :: PRHOFOLD ! snow density increasing constant
LOGICAL,              INTENT(IN)    :: OALL_MELT! T --> all snow runs off if
                                                ! lower surf. temperature is
                                                ! positive
REAL,                 INTENT(IN)    :: PDRAIN_TIME ! drainage folding time (days)
REAL,                 INTENT(IN)    :: PWCRN    ! critical snow amount necessary
                                                ! to cover the considered surface
REAL,                 INTENT(IN)    :: PZ0SN    ! snow roughness length for momentum
REAL,                 INTENT(IN)    :: PZ0HSN   ! snow roughness length for heat
REAL, DIMENSION(:), INTENT(INOUT) :: PWSNOW   ! snow reservoir (kg/m2)
REAL, DIMENSION(:), INTENT(INOUT) :: PTSNOW   ! snow temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PASNOW   ! snow albedo
REAL, DIMENSION(:), INTENT(INOUT) :: PRSNOW   ! snow density
REAL, DIMENSION(:), INTENT(INOUT) :: PTS_SNOW ! snow surface temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PESNOW   ! snow emissivity
REAL, DIMENSION(:), INTENT(IN)    :: PTG      ! underlying ground temperature
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW  ! absorbed SW energy (Wm-2)
REAL, DIMENSION(:), INTENT(IN)    :: PLW1     ! LW coef independant of TSNOW
                                              ! (Wm-2)     usually equal to:
                                              !      emis_snow * LW_down
                                              !
REAL, DIMENSION(:), INTENT(IN)    :: PLW2     ! LW coef dependant   of TSNOW
                                              ! (Wm-2 K-4) usually equal to:
                                              ! -1 * emis_snow * stefan_constant
                                              !
REAL, DIMENSION(:), INTENT(IN)    :: PTA      ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA      ! specific humidity
                                                ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD    ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PPS      ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA    ! air density
                                                ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PSR      ! snow rate
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
                                              ! atmospheric level (temperature)
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the first
                                              ! atmospheric level (wind)
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSNOW ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSNOW  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESNOW ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSNOW  ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT   ! snow melting rate (kg/m2/s)
!
!
END SUBROUTINE SNOW_COVER_1LAYER
!
!
! 
END INTERFACE
!
!
!
END MODULE MODI_SNOW_COVER_1LAYER
