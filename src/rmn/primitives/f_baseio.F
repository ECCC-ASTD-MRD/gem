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
        INTEGER FUNCTION qqqf7op(iun,name,lrec,rndflag,unfflag,lmult)
        integer iun,lrec
        character *(*) name
        integer rndflag,unfflag
        
        integer lng
        qqqf7op=0
        lng = len(name)
*        print *,'opening file ',name(1:lng),' as unit ',iun
*        print *,'lrec=',lrec,' flags = ',rndflag,unfflag
*        print *,'len=',lng
        if (rndflag.eq.1) then
          if (unfflag.eq.1) then
*            print *,'ACCESS=DIRECT,RECL=',lrec
            OPEN(iun,FILE=name(1:lng),ACCESS='DIRECT',RECL=lrec*lmult,
     %           ERR=77)

          else
*            print *,'ACCESS=DIRECT,RECL=',lrec*4
            OPEN(iun,FILE=name(1:lng),ACCESS='DIRECT',
     %           FORM='FORMATTED',RECL=lrec*4,ERR=77)
          endif
        else
         if ((name(1:lng).EQ.'input') .OR. (name(1:lng).EQ.'$input')
     %  .OR. (name(1:lng).EQ.'output') .OR. (name(1:lng).EQ.'$output')
     %  .OR. (name(1:lng).EQ.'$in') .OR. (name(1:lng).EQ.'$out'))
     %    then
*            print *,' STDIN or STDOUT'
            return
          else
            if (unfflag.eq.1) then
*              print *,'UNFORMATTED open'
              OPEN(iun,FILE=name(1:lng),FORM='UNFORMATTED',ERR=77)
            else
*              print *,'FORMATTED open'
#if !defined (HP)
              OPEN(iun,FILE=name(1:lng),FORM='FORMATTED',DELIM='QUOTE',
     %             ERR=77)        
#else
              OPEN(iun,FILE=name(1:lng),FORM='FORMATTED',ERR=77)        
#endif
            endif
          endif
        endif
        return
 77     continue
        qqqf7op = -1
        return
        end

        integer FUNCTION ftnclos(iun)
        integer iun
        
        ftnclos = 0
        CLOSE(iun)
        return
        end
      INTEGER FUNCTION LONGUEUR(NOM)
      implicit none
      CHARACTER * (*) NOM
*
      INTEGER LNG,I
*
      LNG = LEN(NOM)
      DO 10 I = LNG,1,-1
         IF (NOM(I:I) .EQ. ' ') THEN
            LNG = LNG - 1
         ELSE
            GOTO 20
         ENDIF
 10   CONTINUE
 20   CONTINUE
      LONGUEUR = LNG
      RETURN
      END

