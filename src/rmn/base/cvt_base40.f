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
***FONCTION   B40TOC  BASE40 TO CHARACTER CONVERSION
      CHARACTER*8  FUNCTION  B40TOC  ( IWRD )
*
*AUTEUR       GASTON BOISVERT - DDO - DORVAL - 683-8221
*                             - NOV. 1985 -
*
*LANGAGE      fortran 5
*
*OBJET        ROUTINE TO CONVERT A BASE 40 INTEGER TO A CHARACTER STRING
*
*ARGUMENTS    IWRD   - INPUT  - INTEGER   - THE BASE 40 INTEGER VALUE TO
*                                           CONVERT TO A CHARACTER STRING
*
*             VAL    - OUTPUT - CHARACTER - THE CHARACTER STRING
*
**
      PARAMETER    ( IBASE=40 )
      CHARACTER*1  BSE40( 0:IBASE-1 ), ICAR
*
      INTEGER  WORD
*
*   INITIALIZE BSE40 UPPER CASE CHARACTER SET
*

      DATA  BSE40  / '0', '1', '2', '3', '4', '5', '6', '7',
     *               '8', '9', 'A', 'B', 'C', 'D', 'E', 'F',
     *               'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
     *               'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
     *               'W', 'X', 'Y', 'Z', '/', '?', ':', ' ' /
*
      WORD = IWRD
*
*   INITIALIZE STRING TO BLANKS
*
      B40TOC(1:) = ' '
*
*   LOOP TO PROCESS ALL CHARACTERS IN WORD
*
      DO  300  IORD=8,1,-1
          IDVDN = WORD / IBASE
          IRMDR = WORD - ( IDVDN*IBASE )
*
*   GET EQUIVALENT CHARACTER FROM BSE40 ARRAY
*
          ICAR  = BSE40 (IRMDR)
          B40TOC(IORD:IORD) = ICAR
          WORD  = IDVDN
  300 CONTINUE
*
      RETURN
      END
***FONCTION   CTOB40   CHARACTER TO BASE 40 CONVERSION
      INTEGER  FUNCTION  CTOB40  ( STRING )
*
*AUTEUR       GASTON BOISVERT - DDO - DORVAL - 683-8221
*                             - NOV. 1985 -
*
*LANGAGE      fortran 5
*
*OBJET        ROUTINE TO CONVERT A STRING TO AN INTEGER VALUE IN BASE 40
*
*ARGUMENTS
*           STRING - INPUT  - CHARACTER - THE CHARACTER VALUE TO CONVERT
*                                         TO BASE 40
*
*             VAL  - OUTPUT -   INTEGER - INTEGER IN BASE 40
*
**
      PARAMETER      ( IBASE=40 )
      CHARACTER*1    BSE40( 0:IBASE-1 ), ICAR
      CHARACTER*(*)  STRING
*
      INTEGER  POWER, VALUE
*
*   INITIALIZE BSE40 UPPER CASE CHARACTER SET
*

      DATA  BSE40  / '0', '1', '2', '3', '4', '5', '6', '7',
     *               '8', '9', 'A', 'B', 'C', 'D', 'E', 'F',
     *               'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
     *               'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
     *               'W', 'X', 'Y', 'Z', '/', '?', ':', ' ' /
*
      CTOB40 = 0
      VALUE  = 0
*
*   GET THE NUMBER OF CHARACTERS IN STRING
*
      NCARS = LEN ( STRING )
*
*   LOOP ON THE NUMBER OF CHARACTERS IN STRING
*
      DO  300  II=1,NCARS
          ICAR = STRING (II:II)
          POWER = NCARS-II
*
*   LOOP IN BSE40 ARRAY TO FIND ORDINAL OF CHARACTER
*
          DO   200  IORD=0,IBASE-1
               IF ( ICAR .EQ. BSE40(IORD) ) THEN
                  VALUE = VALUE + ( IORD * ( IBASE**POWER ) )
                  GO TO 300
               END IF
  200     CONTINUE
*
*   PRINT MESSAGE FOR ILLEGAL CHARACTER
*
          VALUE = VALUE + ( 36 * ( IBASE**POWER ) )
          WRITE ( 6,'(" ILLEGAL CHARACTER IN WORD ", A)' ) ICAR
          WRITE ( 6,'(" CHARACTER REPLACED BY A SLASH")' )
  300 CONTINUE
*
      CTOB40 = VALUE
*
      RETURN
      END
