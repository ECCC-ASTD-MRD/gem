!/* EDITFST - Collection of useful routines in C and FORTRAN
! * Copyright (C) 1975-2014  Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
!** S/P ZAP - CHANGE LES ARGUMENTS DU LABEL A LA SORTIE
!
      SUBROUTINE ZAP(TV, NV, LBL, DATE, IP1, IP2, IP3)
      use configuration
      IMPLICIT NONE 
  
      INTEGER, intent(IN) ::  TV, NV, LBL(*), DATE, IP1, IP2, IP3
!
!AUTEURS
!VERSION ORIGINALE Y. BOURASSA JUL 91
!REVISION      001 "      "    JUL 92 COUPE DW DE DATE
!              002 "      "    OCT 92 PADING DU LABEL
!              003 B. Dugas    fev 12 valider date avec newdate
!              004 M. Valin    fev 14 initialisation de LIS, intent 
!LANGUAGA FTN77
!  
!ARGUMENTS
!ENTRE   TV   -  TYPEVAR 
!  "     NV   -  NOMVAR  
!  "     LBL  -  ETIKET  
!  "     DATE -  DATE
!  "     IP1  -  IP1 
!  "     IP2  -  IP2
!  "     IP3  -  IP3
!
!MODULE
!  
      EXTERNAL     FSTCVT, ARGDOPE, HOLACAR
      INTEGER,     EXTERNAL :: NEWDATE
!*
      INTEGER      IER,DAT1,DAT2,NDATE
      INTEGER      FSTCVT, ARGDOPE, I, LIS(10)
      CHARACTER *1 G
      CHARACTER *2 T
      CHARACTER *4 N
      CHARACTER *12 E

      ZA = .FALSE.  ! initialize to do nothing values
      ZD = -1
      Z1 = -1
      Z2 = -1
      Z3 = -1
      ZT = '??'
      ZN = '????'
      ZE = '????????????' 

      GO TO(70, 60, 50, 40, 30, 20, 10) NP
   10 IF(IP3 .NE. -1) THEN
         ZA = .TRUE.
         Z3 = IP3
         IF( DEBUG ) WRITE(6,*)' ZIP3 = ',Z3
      ENDIF
   20 IF(IP2 .NE. -1) THEN
         ZA = .TRUE.
         Z2 = IP2
         IF( DEBUG ) WRITE(6,*)' ZIP2 = ',Z2
      ENDIF
   30 IF(IP1 .NE. -1) THEN
         ZA = .TRUE.
         Z1 = IP1
         IF( DEBUG ) WRITE(6,*)' ZIP1 = ',Z1
      ENDIF
   40 IF(DATE .NE. -1) THEN
         ! ZD = DATE - DATE/1000000000*1000000000 !
         ! VERIFIER QUE CE NOUVEAU STAMP EST VALIDE
         IER = NEWDATE( DATE, DAT1,DAT2, -3 )
         IF (IER == 0) THEN
            IER = NEWDATE( NDATE, DAT1,DAT2, +3 )
            IF (IER == 0 .AND. NDATE == DATE) THEN
               ZA = .TRUE.
               ZD = DATE
               IF( DEBUG ) WRITE(6,*)' ZDAT = ',ZD
            ENDIF
         ENDIF
      ENDIF
   50 IF(LBL(1) .NE. -1) THEN
         lis = 0
         I = ARGDOPE(3, LIS, 10)  ! nombre de strings + table de localisation de readlx
         CALL HOLACAR(ZE, LIS, I, LBL, 12)
         IF(ZE .NE. '????????????') THEN
            ZA = .TRUE.
            IF( DEBUG ) WRITE(6,*)' ZETIKET= -',ZE,'-'
         ENDIF
      ENDIF
   60 IF(NV .NE. -1) THEN
         I = FSTCVT(NV, -1, -1, -1, ZN, T, E, G, .TRUE.)
         IF(ZN .NE. '????') THEN
            ZA = .TRUE.
            IF( DEBUG ) WRITE(6,*)' ZNOM= -',ZN,'-'
         ENDIF
      ENDIF
   70 IF(TV .NE. -1) THEN
         I = FSTCVT(-1, TV, -1, -1, N, ZT, E, G, .TRUE.)
            IF(ZT .NE. '??') THEN
            ZA = .TRUE.
            IF( DEBUG ) WRITE(6,*)' ZTYP= -',ZT,' -'
         ENDIF
      ENDIF
     
      RETURN 
      END 




