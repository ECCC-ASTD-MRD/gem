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
***S/R AFIX8 - SETS A REAL*8 ARRAY TO A CONSTANT
*
      SUBROUTINE AFIX8 (R,CON,NA)
      INTEGER NA
      REAL*8 CON, R(NA)
*VERSION REAL*8  DE LA ROUTINE AFIX (J. CAVEEN - SEPT 1993)
*
*AUTHOR   - HERSH MITCHELL  -  DEC 74
*
*REVISION 001:  C. THIBEAULT  -  JUL 79  DOCUMENTATION
*                                        CALL TO VECLIB
*REVISION 002:  C. THIBEAULT  -  MAR 83  CONVERSION AU CODE CRAY
*
*LANGUAGE  - fortran
*
*OBJECT(AFIX8)
*         - SETS ARRAY R(NA) TO CONSTANT CON
*
*LIBRARIES
*         - SOURCE   RMNSOURCELIB,ID=RMNP     DECK=AFIX8
*         - OBJECT   RMNLIB,ID=RMNP
*
*USAGE    - CALL AFIX8(R,CON,NA)
*
*ARGUMENTS
*  OUT    - R   - RESULTS
*  IN     - CON - CONSTANT
*         - NA  - LENGTH OF ARRAY R
*
*------------------------------------------------------------------
*
*
*------------------------------------------------------------------
*
      DO 10 I=1,NA
      R(I) = CON
   10 CONTINUE
*
*-------------------------------------------------------------------
*
      RETURN
      END
