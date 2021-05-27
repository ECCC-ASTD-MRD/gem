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
      SUBROUTINE EWFDR2(AF,BF,CF,ILGH,ILEV,LM,Q,IGN)
C 
C     ***** NOV 1977 - CLAUDE GIRARD *****
C     * CALCULATE EAST-WEST DERIVATIVE OF FOURIER COEFF  CF 
C     * MULTIPLIED BY  Q(A REAL CONSTANT) 
C     * IN OPTION ADD(IGN=+1) OR SUBSTRACT(IGN=-1)  BF TO RESULT
C     * STORE FINAL RESULT IN  AF 
C 
      COMPLEX AF(ILGH,ILEV),BF(ILGH,ILEV),CF(ILGH,ILEV) 
C---------------------------------------------------------------------- 
C 
      IF(IGN.EQ. 0)  GO TO 100
      IF(IGN.EQ.+1)  GO TO 200
      IF(IGN.EQ.-1)  GO TO 300
C 
      WRITE(6,1000)  IGN
 1000 FORMAT(5H IGN=,I5,13H IS NOT VALID) 
      GO TO 400 
C 
  100 CONTINUE
C 
      DO 150 L=1,ILEV 
      DO 150 M=1,LM 
      FMSQ=FLOAT(M-1)*Q 
      AF(M,L)=         CMPLX(-FMSQ*AIMAG(CF(M,L)),+FMSQ*REAL(CF(M,L)))
  150 CONTINUE
      GO TO 400 
C 
  200 CONTINUE
C 
      DO 250 L=1,ILEV 
      DO 250 M=1,LM 
      FMSQ=FLOAT(M-1)*Q 
      AF(M,L)= BF(M,L)+CMPLX(-FMSQ*AIMAG(CF(M,L)),+FMSQ*REAL(CF(M,L)))
  250 CONTINUE
      GO TO 400 
C 
  300 CONTINUE
C 
      DO 350 L=1,ILEV 
      DO 350 M=1,LM 
      FMSQ=FLOAT(M-1)*Q 
      AF(M,L)=-BF(M,L)+CMPLX(-FMSQ*AIMAG(CF(M,L)),+FMSQ*REAL(CF(M,L)))
  350 CONTINUE
C 
  400 RETURN
      END 
