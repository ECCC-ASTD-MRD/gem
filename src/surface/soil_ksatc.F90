!copyright (C) 2001  MSC-RPN COMM  %%%RPNPHY%%%

      SUBROUTINE SOIL_KSATC( &
        WSOIL, ISOIL, WSAT, KSAT, KSATC, N, NLEVELS)

      USE SFC_OPTIONS
      USE SVS_CONFIGS

      IMPLICIT NONE

      ! Input
      INTEGER N ! Number of grid points
      INTEGER NLEVELS ! Number of levels in soil variables and properties

      REAL, DIMENSION(N, NLEVELS) :: WSOIL, ISOIL, WSAT, KSAT

      ! Output
      REAL, DIMENSION(N, NLEVELS) :: KSATC

      !
      !Author
      !          V. Fortin (May 2024)
      !Revisions
      !
      !Object
      ! Adjust hydraulic conductivity at saturation based on ice content

      !
      !Arguments
      !
      !          - INPUT -
      !
      !          ---  Soil properties   ---
      !
      ! WSAT     soil water content at saturation for ice-free soil [m/s]
      ! KSAT     soil hydraulic conductivity for ice-free soil [m/s]
      !
      !          --- Prognostic variables of SVS not modified by SOIL_WSAT_KSATC ---
      !
      ! WSOL     soil volumetric water content [m3/m3]
      ! ISOL     frozen soil volumetric water [m3/m3]
      ! 
      !          - OUTPUT -
      !
      !          ---  Soil properties modified to take into account frozen soil content  ---
      !
      ! KSATC    soil hydraulic conductivity modified based on ice content of soil [m/s]
      !
      !          -  DIMENSIONS  -
      !
      ! N        number of grid cells
      ! NLLEVELS number of levels in soil variables and properties. Not necessarily equal to
      !          NL_SVS because function also used to compute WSATC and KSATC at layer boundaries

      ! Local Variable and arrays
      INTEGER I, K

      IF (lsoil_freezing_svs1 .AND. soil_ksat_ice .EQ. 'ZHANGGRAY97') THEN
         ! Correction factor taken from Zhang and Gray (1997). Same as CLASS 3.6
         DO I=1,N
            DO K=1,NLEVELS
               KSATC(I,K) = KSAT(I,K) * (1.0-MAX(0.0,MIN((WSAT(I,K)-CRITWATER)/WSAT(I,K),ISOIL(I,K)/WSAT(I,K))))**2.
            END DO
         END DO
      ELSE IF (lsoil_freezing_svs1 .AND. soil_ksat_ice .EQ. 'BOONE2000') THEN
         ! Impedance factor taken from SURFEX (Boone et al., 2000)
         DO I=1,N
            DO K=1,NLEVELS
               KSATC(I,K) = KSAT(I,K) * EXP(LOG(10.0)*(-6*ISOIL(I,K)/(ISOIL(I,K)+WSOIL(I,K))))
            END DO
         END DO
      ELSE
         ! No change to KSAT as a function of ice content (not very realistic!)
         KSATC(:,:) = KSAT(:,:)
      ENDIF 

      END SUBROUTINE SOIL_KSATC
