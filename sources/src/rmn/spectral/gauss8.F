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
      SUBROUTINE GAUSS8(NRACP,RACP,PG,SIA,RAD,PGSSIN2,SINM1,SINM2,SIN2)

*     ***************************************************************** 
*     CALCULE LES NRACP RACINES POSITIVES DU POLYNOME DE LEGENDRE DE
*     DEGRE 2*NRACP (ICI-APRES NOTE PN) DEFINI SUR L INTERVALLE DES 
*     COLATITUDES ALLANT DE 0 (POLE NORD) A PI (POLE SUD). ON SAIT QUE
*     LES 2*NRACP RACINES SONT ANTI-SYMETRIQUES P/R A L EQUATEUR PI/2,
*     ETANT POSITIVES ENTRE COLAT=0 ET COLAT =PI/2. 
*     ON CALCULE ENSUITE LES POIDS DE GAUSS ASSOCIES AUX COLATITUDES
*     GAUSSIENNES (ICI APRES NOTEES CG), AINSI QU UN CERTAIN NOMBRE DE
*     FONCTIONS DE CG DEFINIES PLUS LOIN. ON RAPPELLE ENFIN QUE LA LATI-
*     TUDE LAT=COLAT-PI/2, ET DONC QUE SIN(LAT)=COS(COLAT). 
*     NRACP        : NOMBRE DE RACINES POSITIVES DU POLYNOME DE LEGENDRE
*                  : DE DEGRE 2*NRACP.
*     RACP(I)      : RACINES DE PN, =SIN(LG)=COS(CG). 
*     PG(I)        : POIDS DE GAUSS CORRESPONDANTS. 
*     SIA(I)       : SIN(CG)=COS(LG). 
*     RAD(I)       : COLATITUDE CG EN RADIANS.
*     PGSSIN2(I)   : POIDS DE GAUSS / (SIN(CG))**2. 
*     SINM1(I)     : (SIN(CG))**-1. 
*     SINM2(I)     : (SIN(CG))**-2. 
*     VOIR NST 8, CHAP. A, PP.1-7, ET APPENDICE D12, PP. 26-27. 
*     VERSION REVISEE PAR MICHEL BELAND, 9 DECEMBRE 1980. 
*     Version codee en REAL*8 par Bernard Dugas, 4 janvier 1994.
*     ***************************************************************** 

*     -----------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
*     ----------------------------------------------------------------- 
      INTEGER  NRACP
      REAL*8   RACP(NRACP),  PG(NRACP),      SIA(NRACP)
      REAL*8   SINM1(NRACP), SINM2(NRACP),   SIN2(NRACP)
      REAL*8   RAD(NRACP),   PGSSIN2(NRACP)
      EXTERNAL ORDLEG8
*     --------------------------------------------------------------

*     ON DEMANDE UNE PRECISION DE 1.E-13 POUR LES RACINES DE PN.

      XLIM = 1.E-13 
      PIS2 = 2.0*ATAN(1.D0)
      IR   = 2*NRACP
      FI   = IR
      FI1  = FI+1. 
      FN   = NRACP 

*     ON UTILISE UNE FORMULE ASYMPTOTIQUE POUR OBTENIR LES VALEURS
*     APPROXIMATIVES DES COLATITUDES GAUSSIENNES
*     CG(I) = (PI/2) * (2*I-1)/(2*NRACP). 
*     VOIR ABRAMOWITZ AND STEGUN, P. 787, EQU. 22.16.6 .

      DO 20 I=1,NRACP 
         DOT     = I-1
         RACP(I) =-PIS2*(DOT+.5)/FN + PIS2
         RACP(I) = SIN(RACP(I)) 
   20 CONTINUE

*     ON CALCULE ENSUITE LES CONSTANTES FACTEURS DE P(N+1) ET P(N-1)
*     DANS L EXPRESSION DE LA PSEUDO-DERIVEE DE PN. 

      DN  = FI/SQRT(4.*FI*FI-1.) 
      DN1 = FI1/SQRT(4.*FI1*FI1-1.) 
      A   = DN1*FI
      B   = DN*FI1
      IRP = IR + 1
      IRM = IR -1 

*     ON EMPLOIE ENSUITE UNE METHODE DE NEWTON POUR AUGMENTER LA PREC.
*     SI RACTEMP EST UNE SOL. APPROXIMATIVE  DE PN(RACP)=0., ALORS LA 
*     SEQUENCE RACTEMP(I+1)=RACTEMP(I)-PN(RACTEMP(I))/DER.PN(RACTEMP(I))
*     CONVERGE VERS RACP DE FACON QUADRATIQUE.
*     VOIR ABRAMOWITZ AND STEGUN, P.18, EQU. 3.9.5. 
*     ORDLEG CALCULE LA VALEUR DE PN (RACP) , NORMALISE.

      DO 50 I=1,NRACP 
   42    CALL ORDLEG8(G,RACP(I),IR) 
            CALL ORDLEG8( GM,RACP(I),IRM ) 
            CALL ORDLEG8( GP,RACP(I),IRP ) 
            GT      = (A*GP-B*GM)/(RACP(I)*RACP(I)-1.) 
            RACTEMP = RACP(I) - G/GT
            GTEMP   = RACP(I) - RACTEMP 
            RACP(I) = RACTEMP 
         IF( ABS(GTEMP).GT.XLIM) GOTO 42
   50 CONTINUE

*     ON CALCULE ENSUITE LES POIDS DE GAUSS SELON L ALGORITHME
*     PG(I) = 2./[(1.-RACP(I)**2)*(DER.PN(RACP(I)))**2].
*     VOIR ABRAMOWITZ AND STEGUN, P.887, EQU. 25.4.29.
*     NOTE: ON DOIT MULTIPLIER LA PRECEDENTE FORMULE PAR UN FACTEUR 
*     DE DENORMALISATION, LES PN DONNES PAR ORDLEG ETANT NORMALISES.
*     ON SE SERT D UNE FORMULE DE RECURRENCE POUR LA DERIVEE DE PN. 

      DO 60 I=1,NRACP 
         A          = 2.*(1.-RACP(I)**2)
         CALL ORDLEG8( B,RACP(I),IRM )
         B          = B*B*FI*FI 
         PG(I)      = A*(FI-.5)/B 
         RAD(I)     = ACOS(RACP(I))
         SIA(I)     = SIN(RAD(I)) 
         C          = (SIA(I))**2 
         SINM1(I)   = 1./SIA(I)
         SINM2(I)   = 1./C 
         PGSSIN2(I) = PG(I)/C 
         SIN2(I)    = C 
   60 CONTINUE

      RETURN
      END 
