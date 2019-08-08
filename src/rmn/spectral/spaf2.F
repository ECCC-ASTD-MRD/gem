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
      SUBROUTINE SPAF2(CFC,SC,LR,LM,ALP,LALP) 
***S/P  SPAF2       SPECTRAL A FOURIER
*  ENTREE SPAF2 
* 
* 
*AUTEUR   M.VALIN RPN  8-5-78 
*REVISION 001  J.D.HENDERSON   OCT 1975 (FORTRAN) 
*         002  M.VALIN   1977   POUR SPECTRAL V-7 
*         004  M.VALIN   8-5-78 POUR SPECTRAL V-8 
*         003  C.GIRARD  1977  (FORTRAN POUR SPECTRAL V-8 TRONCATURE VAR
*LANGAGE  COMPASS 7600
* 
*OBJET(SPAF2) 
*         CALCULER LES COEFFIGIENTS DE FOURIER A PARTIR DES 
*         COEFFICIENTS SPECTRAUX
* 
*PRECISION
*         COMPLEXE SIMPLE 
* 
*DIRECTIVES(SPAF2)
* 
*APPEL
*         CALL SPAF2(CFC,SC,LR,LM,ALP,LALP) 
* 
*ARGUMENTS
*       S-CFC    COEFFICIENTS DU SPECTRE DE FOURIER 
*                ONDES 0 A (LM-1) 
*         SC     CAS GLOBAL  ONDES  0,1,2,...,(LR-1)
*                CAS SYMETRIQUE  ONDES  0,2,4,...,2*(LR-1)
*                CAS ANTISYMETRIQUE  ONDES  1,3,5,...,2*(LR-1)+1
*         LR     TABLEAU CONTENANT LES INDICATIONS DE TRONCATURE
*         LM     NOMBRE D'ONDES A TRAITER 
*         ALP    COEFFICIENTS DES POLYNOMES DE LEGENDRE POUR UNE
*                LATITUDE 
*         LALP   NE SERT A RIEN 
* 
*              CFC,SC,LR  SONT DES TABLEAUX 
* 
*IMPLICITES 
*         AUCUN 
* 
**INTERNE 
*MEMOIRE
*         SC,CFC,ET ALP SONT DANS LE LCM
*         LR,LM,LALP SONT DANS LE SCM 
* 
*NOTES
*         LES RANGEES DE ALP ET DE SC DOIVENT AVOIR LA MEME STRUCTURE 
* 
**   CODE FORTRAN EQUIVALENT (REVISION 003) 
**
  
      COMPLEX CFC(1),SC(1)
      DIMENSION ALP(1),LR(1)
  
      MR=0
      CFC(LM)=(0.,0.) 
      LMA=LM
      IF(LR(LM).EQ.0) LMA=LM-1
  
      DO  30  M=1,LMA 
      FCR=0.
      FCI=0.
      MR1=MR+1
      MR=MR+LR(M) 
  
      DO  20  MN=MR1,MR 
      FCR=FCR  +  ALP(MN)*REAL(SC(MN))
      FCI=FCI  +  ALP(MN)*AIMAG(SC(MN)) 
   20 CONTINUE
  
      CFC(M)=CMPLX(FCR,FCI) 
   30 CONTINUE
  
      RETURN
      END 
