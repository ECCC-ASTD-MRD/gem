*/* RMNLIB - Library of useful routines for C and FORTRAN programming
* * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
* *                          Environnement Canada
* *
* * This library is free software; you can redistribute it and/or
* * modify it under the terms of the GNU Lesser General Public
* * License as published by the Free Software Foundation,
* * version 2.1 of the License.
* *
* * This library is distributed in the hope that it will be useful,
* * but WITHOUT ANY WARRANTY; without even the implied warranty of
* * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* * Lesser General Public License for more details.
* *
* * You should have received a copy of the GNU Lesser General Public
* * License along with this library; if not, write to the
* * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
* * Boston, MA 02111-1307, USA.
* */
*
***S/P NPSPIN - COMPUTE DIRECTION AND SPEED AT POLE NORTH OR SOUTH
*              USING ALLPTS ON THE 2 CLOSESTS LATITUDE CIRCLE


      SUBROUTINE NPSPIN (SPEED,PSI,U,V,LI,LJ,KASE,PHI,ARBPSI,DELLAM)
      REAL U(LI,LJ), V(LI,LJ)
*
*AUTEUR   - P.SARRAZIN DEC 78 DRPN DORVAL QUEBEC CANADA
*
*REVISION 001  C. THIBEAULT  -  MAR 83  CONVERSION AU CODE CRAY
*
*LANGAGE  - FORTRAN
*
*OBJET(NPSPIN)
*         - USING ALL POINTS ON THE 2 LATITUDES CIRCLE COMPUTE MEAN
*
*APPEL    - CALL NPSPIN(SPEED, PSI, U, V, LI, LJ, KASE, PHI, CLAMDA)
*
*ARGUMENTS
*   OUT   - SPEED  - WIND SPEED IN KNOTS
*         - PSI    - ANGLE BETWEEN THE WIND AND GREEWICH
*   IN    - U,V    - VECTOR FIELD ON LAT LON
*         - LI,LJ  - SIZE OF FIELD U,V
*         - KASE   - 1=GLOBAL 2=NORTH 3=SOUTH
*         - PHI    - LATITUDE OF PT U OR V
*         - ARBPSI - DIRECTION ARBITRAIRE DU VENT QUAND CELUI-CI EST
*                    ZERO. (SI SPEED=0, PSI=ARBPSI)
*         - DELLAM - DIFFERENCE DE LONGITUDE ENTRE DEUX POINTS DE LA
*                    GRILLE.
*
*LIBRAIRIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=NPSPIN
*         - OBJET   RMNLIB,ID=RMNP
*
* ------------------------------------------------------------------------------
*




#include "kases.cdk"
#include "pi.cdk"


* POLIN FUNCTION(9/8 OF UV1 LESS 1/8 OF UV2) / NO OF PTS


      POLIN (UV1 , UV2 ,I) = (UV1 * 1.125 - UV2*0.125) /FLOAT(I)
*
* POSITION OF COORDINATE OF SOUTH OR NORTH POLE


      J1 =1
      J2 =2
       IF(PHI.GT. 90.0 .AND.KASE.EQ.GLOBAL) J1 = LJ
       IF(PHI.GT. 90.0 .AND.KASE.EQ.GLOBAL) J2 = LJ-1


       UV2J1 = 0.0
       UV2J2 = 0.0
      VV2J1 = 0.0
      VV2J2 = 0.0


* ACCUMULATION OF SCALAR QUANTITE
* ------------------------------------------------------------------------------


       DO 20 I=1,LI
*
* ACCUMULATE SUMS OF P-S COMPONENTS.
*
        CLAMBI = (I-1)*DELLAM*DGTORD
        S = -SIN(CLAMBI)
        C =  COS(CLAMBI)
        UI = S*U(I,J1) - C*V(I,J1)
        VI = C*U(I,J1) + S*V(I,J1)
        UV2J1 = UV2J1 + UI
        VV2J1 = VV2J1 + VI
        UI = S*U(I,J2) - C*V(I,J2)
        VI = C*U(I,J2) + S*V(I,J2)
        UV2J2 = UV2J2 + UI
        VV2J2 = VV2J2 + VI
   20   CONTINUE
*
      UPOLE = POLIN (UV2J1,UV2J2,LI)
      VPOLE = POLIN (VV2J1,VV2J2,LI)
      SPEED = SQRT(UPOLE**2 + VPOLE**2)
                        PSI = ARBPSI
      IF (SPEED.NE.0.0) PSI = RDTODG*ATAN2(VPOLE,UPOLE)
      IF (PSI.LT.0.0) PSI = PSI + 360.0


       RETURN
       END
