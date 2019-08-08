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
      SUBROUTINE FASP8(SC,LR,LM,CFC, ALP,LALP,WEIGHT) 
***S/P FASP8        FOURIER A SPECTRAL

*  ENTREE FASP8 

*AUTEUR   M.VALIN   8-5-78

*REVISION 
*         001  J.D.HENDERSON    VERSION 6  (FORTRAN)
*         002  M.VALIN   1977   VERSION 7 (COMPASS) 
*         003  C.GIRARD  1977   VERSION 8 (FORTRAN) 
*         004  M.VALIN   1977   VERSION 8 (COMPASS) 
*         005  B.Dugas   1994   VERSION CODEE REAL/REAL*8

*LANGAGE  COMPASS  7600 

*OBJET(FASP8) 
*         CONTRIBUTION DES ONDES DE FOURIER COMPLEXES  CFC  POUR UNE
*         LATITUDE, A AJOUTER AUX COEFFICIENTS SPECTRAUX SC 
*         CE SOUS-PROGRAMME EST APPELE A CHAQUE LATITUDE, ET,AVANT LE 
*         PREMIER APPEL, SC EST MIS A  (0.0,0.0)

*PRECISION
*         COMPLEXE SIMPLE 

*DIRECTIVES(FASP8)

*APPEL    CALL FASP8(SC,LR,LM,CFC,ALP,LALP,WEIGHT)

*ARGUMENTS
*       S-SC     COEFFICIENTS SPECTRAUX 
*                CAS GLOBAL          ONDES  0,1,2,...,(LR-1)
*                CAS SYMETRIQUE      ONDES  0,2,4,...,2*(LR-1)
*                CAS ANTISYMETRIQUE  ONDES  1,3,5,...,2*(LR-1)+1
*         LR     INFORMATION DE CONTROLE DE TRONCATURE
*         LM     DIMENSION DE LR ET DE CFC
*         CFC    COEFFICIENTS DE FOURIER
*         ALP    COEFFICIENTS DES POLYNOMES DE LEGENDRE POUR UNE
*                LATITUDE 
*         LALP   DIMENSION DE ALP ET DE SC
*         WEIGHT COEFFICIENT DE PONDERATION DES CONTRIBUTIONS POUR
*                CETTE LATITUDE 

*IMPLICITES 
*         AUCUN 

**INTERNE 
*MEMOIRE
*         LES TABLEAUX SC,CFC, ET ALP SONT  DANS LE LCM 
*         LR,LM,ET LALP SONT DANS LE SCM

*NOTES
*         LES RANGEES DE SC ET ALP DOIVENT AVOIR LA MEME STRUCTURE

      INTEGER LALP,LM,LR(LM)
      REAL*8  ALP(LALP),WEIGHT
      COMPLEX SC(LALP),CFC(LM)

      INTEGER MR,MR1,LMA,MN,M
      REAL*8  CFCR,CFCI

      MR  = 0
      LMA = LM

      IF (LR(LM).EQ.0) 
     +LMA=LM-1

      DO 30  M=1,LMA 

         CFCR = REAL( CFC(M) )*WEIGHT
         CFCI = AIMAG( CFC(M) )*WEIGHT 
         MR1  = MR+1
         MR   = MR+LR(M) 

         DO 20 MN=MR1,MR

            SC(MN) = SC(MN)
     +     +CMPLX( SNGL( CFCR*ALP(MN) ),SNGL( CFCI*ALP(MN) ) )

   20    CONTINUE

   30 CONTINUE

      RETURN
      END 
