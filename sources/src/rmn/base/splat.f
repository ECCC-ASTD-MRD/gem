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
***S/R SPLAT - DRAWS UP TO FIFTEEN CURVES ON THE SAME GRAPH
*
      SUBROUTINE SPLAT(A,LI,NC,II,INCR,GMIN,GMAX)
      REAL A(LI,NC)

*
*AUTHOR   - J.D. HENDERSON  -  MAR 75
*
*REVISION 001   C. THIBEAULT  -  NOV 79  DOCUMENTATION
*REVISION 002   C. THIBEAULT  -  MAR 83  CONVERSION AU CODE CRAY
*
*LANGUAGE - fortran
*
*OBJECT(SPLAT)
*         - DRAWS UP TO FIFTEEN CURVES ON THE SAME GRAPH, THE Y-AXIS EXTENDS
*           HORIZONTALLY ACROSS THE PAGE AND THE X-AXIS EXTENDS DOWN THE
*           PAGE FOR AS MANY PAGES AS REQUIRED
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=SPLAT
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL SPLAT(A,LI,NC,II,INCR,GMIN,GMAX)
*
*ARGUMENTS
*         - A    - ARRAY TO BE PLOTTED
*         - LI   - X-DIMENSION
*         - NC   - NUMBER OF CURVES TO BE DRAWN (MAXIMUM 15)
*         - II   - NUMBER OF ELEMENTS IN THE X-DIRECTION TO BE DRAWN
*                  FROM EACH ROW OF A, MUST BE LESS OR EQUAL TO LI
*         - INCR - X-AXIS INCREMENT, ABSOLUTE VALUE IS USED.
*                - INCR=1 DRAWS EVERY POINT
*                - INCR=2 DRAWS EVERY SECOND POINT ETC...
*                - IF INCR>0, THE VALUES OF CURVE NC ARE WRITTEN ALONG
*                  THE X-AXIS
*                - IF INCR<0, THESE VALUES ARE OMITTED.
*         - GMIN - MINIMUM VALUE ON Y-AXIS
*         - GMAX - MAXIMUM VALUE ON Y-AXIS
*                - VALUES OF THE CURVE FALLING OUTSIDE THESE LIMITS ARE
*                  SIMPLY OMITTED, IF GMIN=GMAX, THE SUBROUTINE SELECTS
*                  THE SMALLEST AND LARGEST OF ALL THE VALUES TO BE
*                  PLOTTED AS THE GRAPH LIMITS.
*
*NOTES    - IF II OR NC IS ZERO, THE SUBROUTINE SIMPLY RETURN WITH
*           NO GRAPH DRAWN.
*               DEPENDENT ***
*
*-------------------------------------------------------------------------------
*
      INTEGER NG(101),S(4),IGR(15)
      DATA S/1H ,1HI,1H-,1H*/
      DATA IGR/1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9,
     1         1HA,1HB,1HC,1HD,1HE,1HF/
*
*-----------------------------------------------------------------------
*
      IF(INCR.EQ.0.OR.NC.EQ.0) GO TO 99
      INC=IABS(INCR)
      AMIN=GMIN
      AMAX=GMAX
      IF(AMIN.LT.AMAX) GO TO 14
      AMIN=A(1,1)
      AMAX=AMIN
*
      DO 12 J=1,NC
      DO 12 I=1,II,INC
      T=A(I,J)
      IF(T.LT.AMIN) AMIN=T
      IF(T.GT.AMAX) AMAX=T
   12 CONTINUE
*
      IF(AMIN.EQ.AMAX) AMAX=AMIN+1.
   14 RANGE=AMAX-AMIN
      DA100=100./RANGE
      AX=(AMAX+AMIN)*.5
      WRITE(6,605) IGR(NC),RANGE
      WRITE(6,610) AMIN,AX,AMAX
      WRITE(6,615)
*
      DO 24 I=1,II,INC
      DO 15 K=1,101
   15 NG(K)=S(1)
      NG(51)=S(2)
      IF(I/10*10.EQ.I) NG(51)=S(3)
*
      DO 20 J=1,NC
      Y=(A(I,J)-AMIN)*DA100+1.5
      K=INT(Y)
      IF(K.GE.1.AND.K.LE.101) NG(K)=IGR(J)
   20 CONTINUE
*
      IF(INCR.GT.0) WRITE(6,620) I,A(I,NC),NG
      IF(INCR.LT.0) WRITE(6,625) NG
   24 CONTINUE
*
      WRITE(6,615)
      WRITE(6,610) AMIN,AX,AMAX
   99 RETURN
*
*-----------------------------------------------------------------------
*
  605 FORMAT(8H1 CURVE ,A1,T87,7HRANGE =,1PE12.5)
  610 FORMAT(1H ,16X,1PE12.5,T68,E12.5,T118,E12.5)
  615 FORMAT(1H ,21X,10(10HI---------),1HI)
  620 FORMAT(1H ,I4,1PE15.6,2X,101A1)
  625 FORMAT(1H ,21X,101A1)
      END
