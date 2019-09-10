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
      SUBROUTINE FASP2(SC,LR,LM,CFC, ALP,LALP,WEIGHT) 
***S/P FASP2        FOURIER A SPECTRAL
* 
*  ENTREE FASP2 
* 
*AUTEUR   M.VALIN   8-5-78
* 
*REVISION 
*         001  J.D.HENDERSON    VERSION 6  (FORTRAN)
*         002  M.VALIN   1977   VERSION 7 (COMPASS) 
*         003  C.GIRARD  1977   VERSION 8 (FORTRAN) 
*         004  M.VALIN   1977   VERSION 8 (COMPASS) 
* 
*LANGAGE  COMPASS  7600 
* 
*OBJET(FASP2) 
*         CONTRIBUTION DES ONDES DE FOURIER COMPLEXES  CFC  POUR UNE
*         LATITUDE, A AJOUTER AUX COEFFICIENTS SPECTRAUX SC 
*         CE SOUS-PROGRAMME EST APPELE A CHAQUE LATITUDE, ET,AVANT LE 
*         PREMIER APPEL, SC EST MIS A  (0.0,0.0)
* 
*PRECISION
*         COMPLEXE SIMPLE 
* 
*DIRECTIVES(FASP2)
* 
*APPEL    CALL FASP2(SC,LR,LM,CFC,ALP,LALP,WEIGHT)
* 
*ARGUMENTS
*       S-SC     COEFFICIENTS SPECTRAUX 
*                CAS GLOBAL          ONDES  0,1,2,...,(LR-1)
*                CAS SYMETRIQUE      ONDES  0,2,4,...,2*(LR-1)
*                CAS ANTISYMETRIQUE  ONDES  1,3,5,...,2*(LR-1)+1
*         LR     INFORMATION DE CONTROLE DE TRONCATURE
*         LM     DIMENSION DE LR
*         CFC    COEFFICIENTS DE FOURIER
*         ALP    COEFFICIENTS DES POLYNOMES DE LEGENDRE POUR UNE
*                LATITUDE 
*         LALP   INUTILE (RELIQUAT DE LA VERSION 7) 
*         WEIGHT COEFFICIENT DE PONDERATION DES CONTRIBUTIONS POUR
*                CETTE LATITUDE 
* 
*IMPLICITES 
*         AUCUN 
* 
**INTERNE 
*MEMOIRE
*         LES TABLEAUX SC,CFC, ET ALP SONT  DANS LE LCM 
*         LR,LM,ET LALP SONT DANS LE SCM
* 
*NOTES
*         LES RANGEES DE SC ET ALP DOIVENT AVOIR LA MEME STRUCTURE
* 
      COMPLEX SC(1),CFC(1)
      DIMENSION ALP(1),LR(1)
  
      MR=0
      LMA=LM
      IF(LR(LM).EQ.0) LMA=LM-1
  
      DO  30  M=1,LMA 
      CFCR=REAL(CFC(M))*WEIGHT
      CFCI=AIMAG(CFC(M))*WEIGHT 
      MR1=MR+1
      MR=MR+LR(M) 
  
      DO 20  MN=MR1,MR
  
      SC(MN)=SC(MN)+CMPLX(CFCR*ALP(MN),CFCI*ALP(MN))
   20 CONTINUE
  
   30 CONTINUE
  
      RETURN
      END 
