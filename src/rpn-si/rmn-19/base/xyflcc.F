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
      SUBROUTINE XYFLCC(X,Y,DLAT,DLON,PHI1,PHI2,LAM0,DELTA,NHEM)
      REAL X,Y,DLAT,DLON,PHI1,PHI2,LAM0,DELTA,NHEM


#include "piprm.cdk"


      REAL XPHI1,XPHI2,XLAM0,XDELTA,R0,ALPHA,C1,S1
      SAVE XPHI1,XPHI2,XLAM0,XDELTA,R0,ALPHA,C1,S1
      REAL C2,S2,R,A,B


      LOGICAL XYALL


      DATA XPHI1,XPHI2,XDELTA/-999999.,-999999.,0./
      XYALL = .FALSE.
      GO TO 10


      ENTRY LCCFXY(X,Y,DLAT,DLON,PHI1,PHI2,LAM0,DELTA,NHEM)


      XYALL = .TRUE.
10    CONTINUE


      IF(XPHI1.NE.PHI1 .OR. XPHI2.NE.PHI2) THEN
         XPHI1 = PHI1
         XPHI2 = PHI2
         C1 = COS(DGTORD*XPHI1)
         S1 = SIN(DGTORD*XPHI1)
         S2 = SIN(DGTORD*XPHI2)
         C2 = COS(DGTORD*XPHI2)
         A = ALOG(C1/C2)
         B = ALOG((1+S2)/(1+S1))
         ALPHA = A/(A+B)
         XDELTA = 0.0
         PRINT *,' ALPHA =',ALPHA
      ENDIF
      IF(XDELTA.NE.DELTA) THEN
         XDELTA = DELTA
         R0 = 6371100.*C1*((1+S1)/C1)**ALPHA/ALPHA/XDELTA
         PRINT *,' R0 =',R0
      ENDIF
      IF(XYALL)THEN
         R = SQRT(X*X+Y*Y)
         DLON = LAM0
         IF(X.NE.0.) DLON = LAM0 + RDTODG*ATAN2(X,-Y)/ALPHA
         DLAT = 90. - RDTODG*2.*ATAN((R/R0)**(1./ALPHA))
      ELSE
         R = R0*(COS(DGTORD*DLAT)/(1.+SIN(DGTORD*DLAT)))**ALPHA
         X = R*SIN(DGTORD*(DLON-LAM0))
         Y = -R*COS(DGTORD*(DLON-LAM0))
      ENDIF
      RETURN
      END
