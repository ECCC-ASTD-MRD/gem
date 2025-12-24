!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.

!auto_modi:spll_snowcro.D
MODULE MODI_SNOW3L
INTERFACE


      SUBROUTINE SNOW3L(HSNOWRES, OMEB, HIMPLICIT_WIND,                   &
                PPEW_A_COEF, PPEW_B_COEF,                                 &
                PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,       &
                PSNOWSWE,PSNOWRHO,PSNOWHEAT,PSNOWALB,                     &
                  PSNOWAGE, PTSTEP, PPS, PSR, PUNLOAD, PRR, PPSN3L,PRESA_SV,&
                  PTAR,PTAC,PTG,PSW_RAD,PQA,PVMOD,PWIND_DRIFT,              &
                  PLW_RAD, PRHOA,                                           &
                PUREF,PEXNS,PEXNA,PDIRCOSZW,                              &
                PZREF,PZ0,PZ0EFF,PZ0H,PALB,                               &
                PSOILCOND,PD_G,PLVTT,PLSTT,                               &
                PSNOWLIQ,PSNOWTEMP,PSNOWDZ,                               &
                  PTHRUFAL,PSNOWMELT,PSNREFREEZ,                            &
                  PGRNDFLUX,PEVAPCOR,PSOILCOR,                              &
                PGFLXCOR,PSNOWSFCH, PDELHEATN, PDELHEATN_SFC,             &
                  PDELPHASEN, PDELPHASEN_SFC,                               &
                  PSWNETSNOW,PSWNETSNOWS,PLWNETSNOW,PRESTOREN,              &
                PRNSNOW,PHSNOW,PGFLUXSNOW,                                &
                PHPSNOW,PLES3L,PLEL3L,PEVAP,PSNDRIFT,PRI,                 &
                PEMISNOW,PCDSNOW,PUSTAR,PCHSNOW,PSNOWHMASS,PQS,           &
                PPERMSNOWFRAC,PFORESTFRAC,PZENITH,                        &
                HSNOWDRIFT,OSNOWDRIFT_SUBLIM                     )  
!     ##########################################################################        

IMPLICIT NONE
REAL, INTENT(IN)                       :: PTSTEP
 CHARACTER(LEN=*), INTENT(IN)          :: HSNOWRES
LOGICAL, INTENT(IN)                    :: OMEB       ! True = coupled to MEB. This means surface fluxes ae IMPOSED
 CHARACTER(LEN=*), INTENT(IN)          :: HIMPLICIT_WIND   ! wind implicitation option
REAL, DIMENSION(:), INTENT(IN)         :: PPS, PTAR,PTAC, PSW_RAD, PQA, PVMOD, PWIND_DRIFT,PLW_RAD, PSR, PRR, PUNLOAD 
REAL, DIMENSION(:), INTENT(IN)         :: PTG, PSOILCOND, PD_G, PPSN3L
REAL, DIMENSION(:), INTENT(IN)         :: PZREF, PUREF, PEXNS, PEXNA, PDIRCOSZW, PRHOA, PZ0, PZ0EFF, &
                                       PALB, PZ0H, PPERMSNOWFRAC, PFORESTFRAC 
REAL, DIMENSION(:), INTENT(IN)         :: PPEW_A_COEF, PPEW_B_COEF,                   &
                                        PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF,      &
                                        PPEQ_B_COEF  
REAL, DIMENSION(:), INTENT(IN)    :: PLVTT, PLSTT ! = latent heats for hydrology
REAL, DIMENSION(:), INTENT(OUT)      :: PSNOWALB
REAL, DIMENSION(:,:), INTENT(INOUT)    :: PSNOWHEAT, PSNOWRHO, PSNOWSWE
REAL, DIMENSION(:,:), INTENT(INOUT)    :: PSNOWAGE  ! Snow grain age
REAL, DIMENSION(:,:), INTENT(OUT)    :: PSNOWTEMP
REAL, DIMENSION(:,:), INTENT(OUT)      :: PSNOWLIQ, PSNOWDZ
REAL, DIMENSION(:), INTENT(OUT)     :: PTHRUFAL, PEVAPCOR, PSOILCOR, PGFLXCOR, &
                                       PRESTOREN, PSNOWSFCH, PDELHEATN, PDELHEATN_SFC, &
                                       PDELPHASEN, PDELPHASEN_SFC, PSNOWMELT, PSNREFREEZ
REAL, DIMENSION(:), INTENT(INOUT)      :: PGRNDFLUX
REAL, DIMENSION(:), INTENT(INOUT)      :: PRNSNOW, PHSNOW, PGFLUXSNOW, PLES3L, PLEL3L, &
                                          PHPSNOW, PCDSNOW, PUSTAR, PEVAP
REAL, DIMENSION(:), INTENT(OUT)        :: PSNDRIFT
REAL, DIMENSION(:), INTENT(INOUT)      :: PSWNETSNOW, PLWNETSNOW, PSWNETSNOWS
REAL, DIMENSION(:), INTENT(INOUT)      :: PCHSNOW,PRI
REAL, DIMENSION(:), INTENT(OUT)        :: PEMISNOW, PSNOWHMASS
REAL, DIMENSION(:), INTENT(OUT)        ::  PQS
REAL, DIMENSION(:), INTENT(IN)         :: PZENITH ! solar zenith angle
 CHARACTER(4), INTENT(IN)              :: HSNOWDRIFT        ! Snowdrift scheme :
LOGICAL, INTENT(IN)                    :: OSNOWDRIFT_SUBLIM ! activate sublimation during drift
REAL, DIMENSION(:)                :: PRESA_SV     ! Additional aerodynamic surface resistance if Ta and HR above canopy
END SUBROUTINE SNOW3L
END INTERFACE
END MODULE MODI_SNOW3L
