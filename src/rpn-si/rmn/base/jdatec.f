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
***S/R JDATEC - COMPUTES THE JULIAN CALENDAR DAY, GIVEN A YEAR,
*              A MONTH AND A DAY.
      SUBROUTINE JDATEC(JD,I,J,K)
*
*
*AUTHOR   - C. THIBEAULT  -  JAN 1980
*
*REVISION 001    C. THIBEAULT - JAN 80  DOCUMENTATION
*
*LANGUAGE - fortran
*
*OBJECT(JDATEC)
*         - COMPUTES THE JULIAN CALENDAR DAY, GIVEN A YEAR,
*           A MONTH AND A DAY.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=JDATEC
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL JDATEC(JD,I,J,K)
*
*ARGUMENTS
*         - JD - JULIAN DAY, A UNIQUE INTEGER WHICH MAPS ONE-TO-ONE
*                ONTO TRIPLES OF INTEGERS REPRESENTING YEAR, MONTH,
*                DAY OF THE MONTH.
*         - I  - YEAR
*         - J  - MONTH
*         - K  - DAY OF THE MONTH.
*
*NOTES    - COPIED FROM "COMMUNICATIONS OF THE ACM" (1968), PAGE 657.
*         - IT COVERS A PERIOD OF 7980 YEARS WITH DAY 1 STARTING
*           AT YEAR=-4713, MONTH=11, DAY=25.
*
*-----------------------------------------------------------------------------
*
*  CALCULATE JULIAN CALENDAR DAY
*
      JD        = K-32075+1461*(I+4800+(J-14)/12)/4
     1             +367*(J-2-(J-14)/12*12)/12-3
     2             *((I+4900+(J-14)/12)/100)/4
*
      RETURN
      END
