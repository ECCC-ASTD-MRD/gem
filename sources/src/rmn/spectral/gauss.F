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
      SUBROUTINE GAUSS(NRACP,RACP,PG,SIA,RAD,PGSSIN2,SINM1,SINM2,SIN2)
C 
C     ***************************************************************** 
C     CALCULE LES NRACP RACINES POSITIVES DU POLYNOME DE LEGENDRE DE
C     DEGRE 2*NRACP (ICI-APRES NOTE PN) DEFINI SUR L INTERVALLE DES 
C     COLATITUDES ALLANT DE 0 (POLE NORD) A PI (POLE SUD). ON SAIT QUE
C     LES 2*NRACP RACINES SONT ANTI-SYMETRIQUES P/R A L EQUATEUR PI/2,
C     ETANT POSITIVES ENTRE COLAT=0 ET COLAT =PI/2. 
C     ON CALCULE ENSUITE LES POIDS DE GAUSS ASSOCIES AUX COLATITUDES
C     GAUSSIENNES (ICI APRES NOTEES CG), AINSI QU UN CERTAIN NOMBRE DE
C     FONCTIONS DE CG DEFINIES PLUS LOIN. ON RAPPELLE ENFIN QUE LA LATI-
C     TUDE LAT=COLAT-PI/2, ET DONC QUE SIN(LAT)=COS(COLAT). 
C     NRACP        : NOMBRE DE RACINES POSITIVES DU POLYNOME DE LEGENDRE
C                  : DE DEGRE 2*NRACP.
C     RACP(I)      : RACINES DE PN, =SIN(LG)=COS(CG). 
C     PG(I)        : POIDS DE GAUSS CORRESPONDANTS. 
C     SIA(I)       : SIN(CG)=COS(LG). 
C     RAD(I)       : COLATITUDE CG EN RADIANS.
C     PGSSIN2(I)   : POIDS DE GAUSS / (SIN(CG))**2. 
C     SINM1(I)     : (SIN(CG))**-1. 
C     SINM2(I)     : (SIN(CG))**-2. 
C     VOIR NST 8, CHAP. A, PP.1-7, ET APPENDICE D12, PP. 26-27. 
C     VERSION REVISEE PAR MICHEL BELAND, 9 DECEMBRE 1980. 
C     ***************************************************************** 
C 
C 
C     ----------------------------------------------------------------- 
      DIMENSION RACP(1),PG(1),SIA(1),RAD(1),PGSSIN2(1),SINM1(1),SINM2(1)
     1         ,SIN2(1) 
C     --------------------------------------------------------------
C 
C     ON DEMANDE UNE PRECISION DE 1.E-13 POUR LES RACINES DE PN.
C 
#if defined (CRAY) || defined (NEC64)
      XLIM=1.E-13 
#else
      XLIM=1.E-6
#endif
      PI = 3.1415926535898
      IR = 2*NRACP
      FI=FLOAT(IR)
      FI1=FI+1. 
      FN=FLOAT(NRACP) 
C 
C     ON UTILISE UNE FORMULE ASYMPTOTIQUE POUR OBTENIR LES VALEURS
C     APPROXIMATIVES DES COLATITUDES GAUSSIENNES
C     CG(I) = (PI/2) * (2*I-1)/(2*NRACP). 
C     VOIR ABRAMOWITZ AND STEGUN, P. 787, EQU. 22.16.6 .
C 
      DO 20 I=1,NRACP 
      DOT=FLOAT(I-1)
      RACP(I)=-PI*.5*(DOT+.5)/FN + PI*.5
      RACP(I) =  SIN(RACP(I)) 
   20 CONTINUE
C 
C     ON CALCULE ENSUITE LES CONSTANTES FACTEURS DE P(N+1) ET P(N-1)
C     DANS L EXPRESSION DE LA PSEUDO-DERIVEE DE PN. 
C 
      DN = FI/SQRT(4.*FI*FI-1.) 
      DN1=FI1/SQRT(4.*FI1*FI1-1.) 
      A = DN1*FI
      B = DN*FI1
      IRP = IR + 1
      IRM = IR -1 
C 
C     ON EMPLOIE ENSUITE UNE METHODE DE NEWTON POUR AUGMENTER LA PREC.
C     SI RACTEMP EST UNE SOL. APPROXIMATIVE  DE PN(RACP)=0., ALORS LA 
C     SEQUENCE RACTEMP(I+1)=RACTEMP(I)-PN(RACTEMP(I))/DER.PN(RACTEMP(I))
C     CONVERGE VERS RACP DE FACON QUADRATIQUE.
C     VOIR ABRAMOWITZ AND STEGUN, P.18, EQU. 3.9.5. 
C     ORDLEG CALCULE LA VALEUR DE PN (RACP) , NORMALISE.
C 
      DO 50 I=1,NRACP 
   42 CALL ORDLEG(G,RACP(I),IR) 
      CALL ORDLEG(GM,RACP(I),IRM) 
      CALL ORDLEG(GP,RACP(I),IRP) 
      GT = (A*GP-B*GM)/(RACP(I)*RACP(I)-1.) 
      RACTEMP = RACP(I) - G/GT
      GTEMP = RACP(I) - RACTEMP 
      RACP(I) = RACTEMP 
      IF( ABS(GTEMP).GT.XLIM) GO TO 42
   50 CONTINUE
C 
C     ON CALCULE ENSUITE LES POIDS DE GAUSS SELON L ALGORITHME
C     PG(I) = 2./[(1.-RACP(I)**2)*(DER.PN(RACP(I)))**2].
C     VOIR ABRAMOWITZ AND STEGUN, P.887, EQU. 25.4.29.
C     NOTE: ON DOIT MULTIPLIER LA PRECEDENTE FORMULE PAR UN FACTEUR 
C     DE DENORMALISATION, LES PN DONNES PAR ORDLEG ETANT NORMALISES.
C     ON SE SERT D UNE FORMULE DE RECURRENCE POUR LA DERIVEE DE PN. 
C 
      DO 60 I=1,NRACP 
      A=2.*(1.-RACP(I)**2)
      CALL ORDLEG(B,RACP(I),IRM)
      B = B*B*FI*FI 
      PG(I)=A*(FI-.5)/B 
      RAD(I) =   ACOS(RACP(I))
      SIA(I) =  SIN(RAD(I)) 
      C=(SIA(I))**2 
      SINM1(I) = 1./SIA(I)
      SINM2(I) = 1./C 
      PGSSIN2(I) =PG(I)/C 
      SIN2(I)=C 
   60 CONTINUE
C 
      RETURN
      END 
