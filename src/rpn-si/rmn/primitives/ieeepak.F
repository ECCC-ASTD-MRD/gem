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
#include <ftnmacros.hf>
#if !defined (_FLOAT1)
***S/P IEEEPAK - COMPACTION / DECOMPACTION DE REELS ("format ieee")
*
      SUBROUTINE IEEEPAK(IFLD,IPK,NI,NJ,NPAK,serpas,KIND)
      IMPNONE
      INTEGER NI,NJ,IFLD(0:*),IPK(*),NPAK,KIND,serpas
*
*AUTHOR    - M. LEPINE  - Nov 91
*
*LANGUAGE - FORTRAN
*
*OBJET(IEEEPAK)
*         - COMPACTER OU DECOMPACTER DES NOMBRES REELS 
*           en format ieee  (1 signe 8 exposant 32-nbits mantisse)
*
*PRECISION
*         - L'ALGORITHME DE COMPACTION GARDE LES NBITS BITS LES
*           MOINS SIGNIFICATIFS.
*
*LIBRAIRIES
*
*ARGUMENTS
* IN/OUT  - IFLD - INTEGER ARRAY, UNPACKED (NI,NJ)
* IN/OUT  - IPK  - THE SAME ARRAY PACKED
*   IN    - NI   - X-DIMENSION
*   IN    - NJ   - Y-DIMENSION
*   IN    - NPAK - NPAK > 0, DENSITE DE COMPACTION
*                  NPAK = 0 OU 1, PAS DE COMPACTION
*                  NPAK < 0, NOMBRE DE BITS A GARDER PAR NOMBRE
*   IN    - kind - 1 = PACK
*                  2 = UNPACK
*   in    - serpas sert pas !!
**
*
#include "xxpakc.cdk"
*
*
      INTEGER I,INDX,ITEMP,LEFT,ITEM,NINJ,INEXT
      INTEGER RBITS, RLEFT, ISAVE
*
      XXPERR = 0
      IF (NPAK .GT. 1) THEN
        NBITS = MAX(1,BITMOT/NPAK)
      ELSE IF (NPAK .LT. 0) THEN
        NBITS = MIN(64,-NPAK)
      ELSE
        NBITS = BITMOT
      ENDIF
*
      NINJ = NI*NJ
*
      RBITS = RMASK(NBITS)
      NWR = (NINJ*NBITS+BITMOT-1) / BITMOT
*

#if defined (ALL64)
      IF (NBITS .NE. 32) THEN
         print *,'IEEEPAK ERROR: NBITS MUST BE 32, NBITS=',nbits
         stop 'IEEEPAK' 
      ENDIF
      IF (KIND .EQ. 1) THEN
*                                               COMPACTER
         CALL MOVR8R4(IFLD,IPK,NINJ)
      ELSE
*                                               DECOMPACTER
         CALL MOVR4R8(IPK,IFLD,NINJ)
      ENDIF
#else

      IF (KIND .EQ. 1) THEN
*                                               COMPACTER
        IF (NBITS .GE. BITMOT) THEN
           IF (NBITS .EQ. BITMOT) THEN
              DO 10 I = 0,NWR-1
                 IPK(I+1) = IFLD(I)
 10           CONTINUE
           ELSEIF (NBITS .EQ. 64) THEN
              CALL MOVE6432(IFLD,IPK,NINJ)
           ELSE
              print *,'IEEEPAK ERROR: NBITS MUST BE <= 32 OR 64, NBITS='
     %         ,nbits
              stop 'IEEEPAK' 
           ENDIF
        ELSE
          INDX = 0
          LEFT = BITMOT
          ITEMP = 0
          DO 100 I=0,NINJ-1
            ITEM = IAND(RSHIFT(IFLD(I),32-NBITS) , RBITS)
            IF (NBITS .LE. LEFT) THEN
              ITEMP = IOR(LYSHIFT(ITEMP,NBITS) , ITEM)
              LEFT = LEFT - NBITS
            ELSE
              INDX = INDX+1
              IPK(INDX) = IOR(LYSHIFT(ITEMP,LEFT) ,
     %                        RSHIFT(ITEM,NBITS-LEFT))
              ITEMP = IAND(ITEM , RMASK(NBITS-LEFT))
              LEFT = LEFT - NBITS + BITMOT
            ENDIF
 100      CONTINUE

          IF (LEFT .LT. BITMOT) THEN
            INDX = INDX + 1
            IPK(INDX) = IOR(LYSHIFT(ITEMP,LEFT),
     %                  IAND(IPK(INDX),RMASK(LEFT)))
          ENDIF
        ENDIF

      ELSE
*                                               DECOMPACTER
        IF (NBITS .GE. BITMOT) THEN
           IF (NBITS .EQ. BITMOT) THEN
              DO 110 I = 0,NWR-1
                 IFLD(I) = IPK(I+1)
 110          CONTINUE
           ELSEIF (NBITS .EQ. 64) THEN
              CALL MOVE6432(IPK,IFLD,NINJ)
           ELSE
              print *,'IEEEPAK ERROR: NBITS MUST BE <= 32 OR 64, NBITS='
     %               ,nbits
              stop 'IEEEPAK' 
           ENDIF
        ELSE
          ITEMP = 0
          ISAVE = IPK(1+ITEMP)
          DO 120 I=NWR,2,-1
            IFLD(I-NWR+NINJ-1) = IPK(I+ITEMP)
 120      CONTINUE
          INDX = (NINJ-NWR) + 1
          LEFT = BITMOT
          ITEMP = LXSHIFT(ISAVE,BITMOT-LEFT)
          INEXT = IFLD(INDX)

          DO 200 I=0,NINJ-1
            IF (NBITS .LE. LEFT) THEN
              IFLD(I) = IAND(RYSHIFT(ITEMP,BITMOT-NBITS),RBITS)
              ITEMP = LYSHIFT(ITEMP,NBITS)
              LEFT = LEFT - NBITS
            ELSE
              RLEFT = RMASK(NBITS-LEFT)
              IFLD(I)=IOR(IAND(IAND(RYSHIFT(ITEMP,BITMOT-NBITS),RBITS),
     %            COMPL(RLEFT)),IAND(RYSHIFT(INEXT,BITMOT-(NBITS-LEFT)),
     %            RLEFT))
              INDX = INDX+1
              ITEMP = LYSHIFT(INEXT,NBITS-LEFT)
              INEXT = IFLD(INDX)
              LEFT = LEFT - NBITS + BITMOT
            ENDIF
            IFLD(I) = LSHIFT(IFLD(I),32-NBITS)
 200      CONTINUE
        ENDIF
      ENDIF
*
#endif
      RETURN
      END
#endif
#if defined (_FLOAT1)
***S/P IEEEPAK - COMPACTION / DECOMPACTION DE REELS ("format ieee")
*
      SUBROUTINE IEEEPAK(IFF,IPK,NI,NJ,NPAK,serpas,KIND)
      IMPNONE
      INTEGER NI,NJ,IFF(0:*),IPK(*),NPAK,KIND,serpas
*
*AUTHOR    - M. LEPINE  - Nov 92
*
*LANGUAGE - RATFOR
*
*OBJET(IEEEPAK)
*         - COMPACTER OU DECOMPACTER DES NOMBRES REELS (version NEC)
*           en format ieee  (1 signe 8 exposant 32-nbits mantisse)
*
*PRECISION
*         - L'ALGORITHME DE COMPACTION GARDE LES NBITS BITS LES
*           MOINS SIGNIFICATIFS.
*
*LIBRAIRIES
*
*ARGUMENTS
* IN/OUT  - IFLD - INTEGER ARRAY, UNPACKED (NI,NJ)
* IN/OUT  - IPK  - THE SAME ARRAY PACKED
*   IN    - NI   - X-DIMENSION
*   IN    - NJ   - Y-DIMENSION
*   IN    - NPAK - NPAK > 0, DENSITE DE COMPACTION
*                  NPAK = 0 OU 1, PAS DE COMPACTION
*                  NPAK < 0, NOMBRE DE BITS A GARDER PAR NOMBRE
*   IN    - kind - 1 = PACK
*                  2 = UNPACK
*   in    - serpas sert pas !!
**
*
#include "xxpakc.cdk"
*
*
      INTEGER I,INDX,ITEMP,LEFT,ITEM,NINJ,INEXT
      INTEGER IFLD,IERR,IMANT,IEXP,INFINITY,BIT01
      POINTER (IPFLD,IFLD(0:NI*NJ))
      INTEGER RBITS, RLEFT, ISAVE
      INTEGER  Z800000,Z400000,Z200000
      INTEGER EXPADJ(0:3)
      DATA Z800000 /HEX(800000)/
      DATA Z400000 /HEX(400000)/
      DATA Z200000 /HEX(200000)/
      DATA INFINITY/HEX(FF800000)/
      DATA EXPADJ /0,3,2,1/
*
      CALL HPALLOC(IPFLD,NI*NJ,IERR,1)
      XXPERR = 0
      IF (NPAK .GT. 1) THEN
        NBITS = MAX(1,BITMOT/NPAK)
      ELSE IF (NPAK .LT. 0) THEN
        NBITS = MIN(BITMOT,-NPAK)
      ELSE
        NBITS = BITMOT
      ENDIF
*
      NBITS = MIN(32,NBITS)
      NINJ = NI*NJ
*
      RBITS = RMASK(NBITS)
      NWR = (NINJ*NBITS+BITMOT-1) / BITMOT
*
      IF (KIND .EQ. 1) THEN
*                                                COMPACTER
        DO 100 I = 0,NINJ-1
#if defined (ALL64)
             IFLD(I) = RSHIFT(IAND(IFF(I),MASK(1)),32)
             IMANT = IAND(RSHIFT(IFF(I),32),RMASK(24))
             IEXP = LSHIFT(IAND(RXSHIFT(RSHIFT(IFF(I),32),24),127)-64,2)
     %              + 126
#else
             IFLD(I) = IAND(IFF(I),MASK(1))
             IMANT = IAND(IFF(I),RMASK(24))
             IEXP = LSHIFT(IAND(RXSHIFT(IFF(I),24),127)-64,2) + 126
#endif
             IF (IAND(IMANT,Z800000) .NE. 0) THEN
               IFLD(I) = IOR(IFLD(I),IAND(IMANT,RMASK(23)))
             ELSE IF (IAND(IMANT,Z400000) .NE. 0) THEN
               IFLD(I) = IOR(IFLD(I),LXSHIFT(IAND(IMANT,RMASK(22)),1))
               IEXP = IEXP - 1
             ELSE IF (IAND(IMANT,Z200000) .NE. 0) THEN
               IFLD(I) = IOR(IFLD(I),LXSHIFT(IAND(IMANT,RMASK(21)),2))
               IEXP = IEXP - 2
             ELSE
               IFLD(I) = IOR(IFLD(I),LXSHIFT(IAND(IMANT,RMASK(20)),3))
               IEXP = IEXP - 3
             ENDIF
             IFLD(I) = IOR(IFLD(I),LSHIFT(IEXP,23))
             IF (IEXP .LT. 0) IFLD(I) = 0
             IF (IEXP .GT. 254) IFLD(I) = INFINITY
 100    CONTINUE
        IF (NBITS .EQ. BITMOT) THEN
          DO 110 I = 0,NWR-1
            IPK(I+1) = IFLD(I)
 110      CONTINUE
*         CALL MOVLEV(IFLD,IPK,NWR)
        ELSE
          INDX = 0
          LEFT = BITMOT
          ITEMP = 0
          DO 200 I=0,NINJ-1 
            ITEM = IAND(RSHIFT(IFLD(I),32-NBITS) , RBITS)
            IF (NBITS .LE. LEFT) THEN
              ITEMP = IOR(LYSHIFT(ITEMP,NBITS) , ITEM)
              LEFT = LEFT - NBITS
            ELSE
              INDX = INDX+1
              IPK(INDX) = IOR(LYSHIFT(ITEMP,LEFT) ,
     %                        RSHIFT(ITEM,NBITS-LEFT))
              ITEMP = IAND(ITEM , RMASK(NBITS-LEFT))
              LEFT = LEFT - NBITS + BITMOT
            ENDIF
 200      CONTINUE

          IF (LEFT .LT. BITMOT) THEN
            INDX = INDX + 1
            IPK(INDX) = IOR(LYSHIFT(ITEMP,LEFT),
     %                      IAND(IPK(INDX),RMASK(LEFT)))
          ENDIF
        ENDIF

      ELSE
*                                                DECOMPACTER
        IF (NBITS .EQ. BITMOT) THEN
          DO 220 I = 0,NWR-1
            IFLD(I) = IPK(I+1)
 220      CONTINUE
        ELSE
          ITEMP = 0
          ISAVE = IPK(1+ITEMP)
#if defined (ALL64)
*CDIR NOVECTOR
#endif
          DO 230 I=NWR,2,-1
            IFLD(I-NWR+NINJ-1) = IPK(I+ITEMP)
 230      CONTINUE
          INDX = (NINJ-NWR) + 1
          LEFT = BITMOT
          ITEMP = LXSHIFT(ISAVE,BITMOT-LEFT)
          INEXT = IFLD(INDX)

          DO 300 I=0,NINJ-1
            IF (NBITS .LE. LEFT) THEN
              IFLD(I) = IAND(RYSHIFT(ITEMP,BITMOT-NBITS),RBITS)
              ITEMP = LYSHIFT(ITEMP,NBITS)
              LEFT = LEFT - NBITS
            ELSE 
              RLEFT = RMASK(NBITS-LEFT)
              IFLD(I) =  IOR(
     %         IAND(IAND(RYSHIFT(ITEMP,BITMOT-NBITS),RBITS),COMPL(RLEFT)),
     %         IAND(RYSHIFT(INEXT,BITMOT-(NBITS-LEFT)),RLEFT))
              INDX = INDX+1
              ITEMP = LYSHIFT(INEXT,NBITS-LEFT)
              INEXT = IFLD(INDX)
              LEFT = LEFT - NBITS + BITMOT
            ENDIF
            IFLD(I) = LSHIFT(IFLD(I),32-NBITS)
 300      CONTINUE
        ENDIF
        DO 400 I = 0,NINJ-1
            IFF(I) = IAND(IFLD(I),LSHIFT(1,31))
            IMANT = IOR(IAND(IFLD(I),RMASK(23)),Z800000)
            IEXP = IAND(RSHIFT(IFLD(I),23),RMASK(8)) - 126
            BIT01 = IAND(IEXP,RMASK(2))
            IEXP = (IEXP+EXPADJ(BIT01) + 256) / 4
            IMANT = RSHIFT(IMANT,EXPADJ(BIT01))
#if defined (ALL64)
            IFF(I) = LXSHIFT(IOR(IOR(IFF(I),LSHIFT(IEXP,24)),
     %                      IAND(IMANT,RMASK(24))),32)
#else
            IFF(I) = IOR(IOR(IFF(I),LSHIFT(IEXP,24)),IAND(IMANT,RMASK(24)))
#endif
 400    CONTINUE
      ENDIF
*
      RETURN
      END
#endif
