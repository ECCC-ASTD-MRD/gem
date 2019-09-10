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
***S/R DATEC - CONVERT JULIAN DAY TO YEAR, MONTH AND DAY.
*
      SUBROUTINE DATEC(JD,I,J,K)
*
*AUTHOR   - FLIEGEL AND FLANDERN
*
*REVISION 001   C. THIBEAULT  -  NOV 79  DOCUMENTATION
*REVISION 002   C. THIBEAULT  -  MAR 83  DATEC (ANCIENNEMENT DATE)
*                                CHANGEMENT DE NOM A CAUSE DE LA ROUTINE
*                                DU SYSTEME QUI PORTE LE MEME NOM
*
*LANGUAGE - fortran
*
*OBJECT(DATEC)
*         - DATEC CONVERTS A JULIAN DAY TO THREE INTEGERS, YEAR(I),
*           MONTH(J) AND DAY(K).
*
*LIBRARIES
*         - SOURCE   RMNSOURCELIB,ID=RMNP     DECK=DATEC
*         - OBJECT   RMNLIB,ID=RMNP
*
*USAGE    - CALL DATEC(JD,I,J,K)
*
*ARGUMENTS
*   IN    - JD - JULIAN DAY, A UNIQUE INTEGER WHICH MAPS ONE-TO-ONE ONTO
*                TRIPLES OF INTEGERS REPRESENTING YEAR, MONTH,
*                DAY OF MONTH.
*   OUT   - I  - YEAR
*         - J  - MONTH
*         - K  - DAY OF THE MONTH
*
*NOTES    - COPIED FROM "COMMUNICATIONS OF THE ACM" (1968), PAGE 657.
*         - IT COVERS A PERIOD OF 7980 YEARS WITH DAY 1 STARTING
*           AT YEAR=-4713, MONTH=11, DAY=25
*
*-------------------------------------------------------------------------------
*
      L= JD+68569
      N= 4*L/146097
      L= L-(146097*N+3)/4
      I= 4000*(L+1)/1461001
      L= L-1461*I/4+31
      J= 80*L/2447
      K= L-2447*J/80
      L= J/11
      J= J+2-12*L
      I= 100*(N-49)+I+L
*
*-------------------------------------------------------------------
*
      RETURN
      END
