!copyright (C) 2001  MSC-RPN COMM  %%%RPNPHY%%%
module soil_frzdepth_svs2_mod
  implicit none
  public
contains
      SUBROUTINE SOIL_FRZDEPTH_SVS2( &
        WSOIL, ISOIL, DL, DZ, FRZDEPTH, FRZTHICK, THWDEPTH, N, NLEVELS)

      USE SFC_OPTIONS
      USE SVS_CONFIGS

      IMPLICIT NONE

      ! Input
      INTEGER N ! Number of grid points
      INTEGER NLEVELS ! Number of levels in soil variables and properties

      REAL, DIMENSION(N, NLEVELS) :: WSOIL, ISOIL
      REAL, DIMENSION(NLEVELS) :: DL, DZ

      ! Output
      REAL, DIMENSION(N) :: FRZDEPTH, FRZTHICK, THWDEPTH

      !
      !Author
      !          B. Bouchard (April 2026)
      !Revisions
      !
      !Object
      ! Compute the variables related to frozen ground (depth and thickness of frozen ground; depth of thawed ground)

      !
      !Arguments
      !
      !          - INPUT -
      !
      !          ---  Soil layer properties   ---
      !
      ! DL       depth of the soil layer [m]
      ! DZ       thickness of the soil layer [m]
      !
      !          --- Prognostic variables of SVS 
      !
      ! WSOIL    soil volumetric water content [m3/m3]
      ! ISOIL    soil volumetric ice content [m3/m3]
      ! 
      !          - OUTPUT -
      !
      !          ---  depth and thickness of frozeng round, depth of thawed ground  ---
      !
      ! FRZDEPTH depth of frozen ground [m]
      ! FRZTHICK thickness of frozen ground [m]
      ! THWDEPTH depth of thawed ground [m]
      !
      !          -  DIMENSIONS  -
      !
      ! N        number of grid cells
      ! NLLEVELS number of levels in soil variables and properties.

      ! Local Variable and arrays
      INTEGER I, K
      INTEGER FRZLAY !Deepest frozen soil layer
      REAL, DIMENSION(N, NLEVELS)   :: FRZFRAC, FRZMASK, ZTOPVEC

      DO I=1,N
    
         !Initialize thickness of frozen ground to 0
         FRZTHICK(I) = 0.0
         FRZLAY = 1
         DO K=1,NLEVELS
            
            !Vectorize top boundary of each layer over space (to compute freezing depth, thickness and thawing depth)
            ZTOPVEC(I,K) = DL(K) - DZ(K)
    
            !Compute the fraction of frozen liquid water (WSOIL including the residual unfrozen water so FRZFRAC does not reach 1.0) - to be corrected!
            FRZFRAC(I,K) = MAX(0.0,MIN(ISOIL(I,K)/(WSOIL(I,K) + ISOIL(I,K)),1.0))
            IF (FRZFRAC(I,K) .GT. 0.0) THEN
               !Mask to estimate the thickness of frozen ground
               FRZMASK(I,K) = 1.0
               !Deepest layer with ice
               FRZLAY = K
            ELSE
               FRZMASK(I,K) = 0.0
            END IF
            
            !Cumulate the thickness of each frozen layer
            FRZTHICK(I) = FRZTHICK(I) + FRZMASK(I,K)*DZ(K)
         END DO
         
         !Frozen depth: top of the deepest frozenlayer + ice fraction of the deepest frozen layer X thickness of that layer
         FRZDEPTH(I) = ZTOPVEC(I,FRZLAY) + FRZFRAC(I,FRZLAY)*DZ(FRZLAY)
         
         !Freezing thickness: consider only the frozen fraction of the deepest frozen layer
         IF(FRZTHICK(I) .GT. 0.0) THEN
             FRZTHICK(I) = FRZTHICK(I) - DZ(FRZLAY)*(1 - FRZFRAC(I,FRZLAY))
         END IF
    
         !Thawing depth: difference between the freezing depth and the frozen thickness
         THWDEPTH(I) = MAX(FRZDEPTH(I)-FRZTHICK(I),0.0)
      END DO


    END SUBROUTINE SOIL_FRZDEPTH_SVS2
  end module soil_frzdepth_svs2_mod
