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
      SUBROUTINE DIMCAL(LRS,LRA,LM,LA,IR,KHEM,NTRUNC) 
C 
C     *****   NOV 1975  -  JOHN D. HENDERSON  ****
C     * COMPUTES ROW LENGTHS AND DIMENSIONS FOR SPECTRAL ARRAYS 
C     * GIVEN THE RESOLUTION (IR,IEW) AND REGION OF INTEGRATION (KHEM). 
C 
C     * NTRUNC = 1 RHOMBOIDAL,  NTRUNC=2 TRIANGULAR.
C 
C     * IF THE MODEL IS HEMISPHERIC ALL FIELDS ARE EITHER SYMMETRIC 
C     *  OR ANTISYMMETRIC, AND THE ZEROES ARE SQUEEZED OUT. 
C 
C     * LRS = ROW LENGTH OF SYMMETRIC FIELDS
C     * LRA = ROW LENGTH OF ANTISYMMETRIC FIELDS
C     * LM  = NUMBER OF ROWS IN EACH FIELD. 
C     * LA  = SEPARATION BETWEEN SUCCESSIVE LEVELS OF ALL OTHER FIELDS. 
C 
      DIMENSION LRS(1),LRA(1) 
C-----------------------------------------------------------------------
      IF(NTRUNC.EQ.1) GO TO 50
      IF(NTRUNC.EQ.2) GO TO 150 
C 
C     * RHOMBOIDAL CASE 
C 
   50 LA = 0
      DO 100 M=1,LM 
C 
      LRS(M) = IR+1 
       IF(KHEM.NE.0) LRS(M) = IR/2 + 1
C 
      LRA(M) = IR+1 
      IF(KHEM.NE.0) LRA(M) = (IR+1)/2 
C 
      LA = LA + LRS(M)
  100 CONTINUE
      GO TO 300 
C 
C     * TRIANGULAR CASE 
C 
  150 IF(LM-1.LE.IR) GO TO 151
      WRITE(6,6000) 
      STOP
C 
  151 LA = 0
      DO 200 M=1,LM 
C 
      LRS(M) = IR-M+2 
      IF(KHEM.NE.0) LRS(M) = (IR-M+1)/2 + 1 
C 
      LRA(M) = IR-M+2 
      IF(KHEM.NE.0) LRA(M) = (IR-M+2)/2 
C 
      LA = LRS(M) + LA
  200 CONTINUE
C 
  300 RETURN
 6000 FORMAT(37H0 INCORRECT TRIANGULAR SPECIFICATION )
      END 
