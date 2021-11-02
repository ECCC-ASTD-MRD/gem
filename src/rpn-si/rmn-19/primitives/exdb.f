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
****FONCTION EXDB    IMPRESSION DE BOITE DE DEBUT D'EXECUTION
*
*
      INTEGER FUNCTION EXDB(in_TITRE,REVIS,FLAG)
      IMPLICIT NONE
      CHARACTER*(*) in_TITRE,REVIS,FLAG
*
      INTEGER EXDBPLUS
      EXTERNAL EXDBPLUS

      CHARACTER *90 UNUSEDSTRING
      EXDB=EXDBPLUS(in_TITRE,REVIS,FLAG,UNUSEDSTRING,0)
      RETURN
      END

**FONCTION EXDBPLUS    IMPRESSION DE BOITE DE DEBUT D'EXECUTION
*
*
      INTEGER FUNCTION EXDBPLUS(in_TITRE,REVIS,FLAG,SUPP,NSUP)
      IMPLICIT NONE
      INTEGER NSUP
      CHARACTER*(*) in_TITRE,REVIS,FLAG,SUPP(NSUP)
*
*AUTEUR  M.VALIN RPN MARS 1983
*Revision 002 M. Lepine Octobre 1998 - Ajout de la signature RMNLIB
*Revision 003 M. Lepine Octobre 2008 - Anciennement exdb, ajout de lignes d'impression
*                                      supplementaires en option
*
*LANGAGE RATFOR
*        IL Y A DES DEPENDANCES CDC NOS/SCOPE 2
*
*OBJET(EXDB)
*        IMPRIMER UNE BOITE INDIQUANT LE DEBUT DE L'EXECUTION
*        D'UN PROGRAMME FORTRAN. EXDB RENVOIE EGALEMENT LE
*        DATE TIME STAMP (NOMBRE ENTIER DE 10 CHIFFRES) DE
*        LA PASSE OPERATIONNELLE COURANTE.
*
*        EXFIN IMPRIME UNE BOITE INDIQUANT LA FIN DE L'EXECUTION
*        D'UN PROGRAMME FORTRAN. LA VALEUR RENVOYEE EST ZERO.
*
*CALL STDLIB
*
*ARGUMENTS
*  E     TITRE     CHAINE DE CARACTERES. SEULS LES 7 PREMIERS
*                  SERONT IMPRIMES.
*  E     REVIS     CHAINE DE CARACTERES INDIQUANT LA VERSION DU
*                  PRODUIT (EXDB) OU UN MESSAGE DE FIN (EXFIN). SEULS LES
*                  DIX PREMIERS CARACTERES SERONT IMPRIMES.
*  E     FLAG      CHAINE DE CARACTERES (7), FLAG PEUT ETRE
*                  - 'NON' AUQUEL CAS IOPDATM RENVOIE 2010101011.
*                  - UNE CLE DE 7 CARACTERES OU MOINS CORRESPONDANT
*                    A D'AUTRES "DATE TIME STAMP".
*                  - UN NOMBRE DE 7 CHIFFRES YYJJJZZ OU YY EST L'ANNEE,
*                    JJJ LE JOUR DANS L'ANNEE (1 A 366) ET ZZ L'HEURE
*                    (00 - 24).
*                  - IOPDATM RENVOIE ALORS LE DATE TIME STAMP QUI
*                    CORRESPOND A CETTE DATE ET CETTE HEURE.
*                  - UN NOM DE FICHIER. IOPDATM LIRA LE PREMIER MOT DE CE
*                    FICHIER EN LE CONSIDERANT COMME UN DATE TIME STAMP.
*                    SI IOPDATM RENCONTRE UNE ERREUR, LE "STAMP" RENVOYE
*                    SERA 1010101011.
*  E     SUPP      LIGNES D'INFORMATION SUPPLEMENTAIRE A IMPRIMER
*  E     NSUP      NOMBRE DE LIGNE SUPPLEMENTAIRES A IMPRIMER
*
*MESSAGES
*        EXDB ET EXFIN ECRIVENT SUR L'UNITE FORTRAN 6. L'USAGER
*        EST RESPONSABLE D'OUVRIR  CETTE UNITE CORRECTEMENT.
*
*MODULES
*        IOPDATM   POUR OBTENIR LE DATE TIME STAMP  (RMNLIB5)
*        DATE      POUR OBTENIR LA DATE DU SYSTEME  (RMNLIB5)
*        DATMGP2   POUR RECONSTITUER LE DATE TIME STAMP  (RMNLIB5)
*        TIME      POUR OBTENIR L'HEURE  (LIBRAIRIE FORTRAN)
*        SECOND    POUR OBTENIR LE TEMPS CPU  (LIBRAIRIE FORTRAN)
*
*NOTES
*        EXDB VA CHERCHER DANS RA+56 A RA+63 L'IMAGE DE LA CARTE
*        DE CONTROLE QUI A SERVI A APPELER LE PROGRAMME PRINCIPAL.
*        CECI EST UNE DEPENDANCE CDC  (NOS OU SCOPE 2)
*        EXDB SE SERT DU FICHIER FORTRAN 88, L'OUVRE ET LE
*        RETOURNE.
*
**
*


      CHARACTER *24 CDATIM
      CHARACTER *105 VERSION,titre,tempstring
      INTEGER EXFIN,IOPDATM,I,IDATIM(14)
      REAL T1,SECOND
      external flush_stdout
      SAVE T1

      titre = ' '
      titre(1:min(len(in_titre),90)) = in_titre(1:min(len(in_titre),90))
      IDATIM(14) = IOPDATM(FLAG)
      CALL DATMGP2(IDATIM)
      CALL FDATE(CDATIM)
* Obtenir la version de rmnlib utilisee
      CALL RMNLIB_version(VERSION,.false.)
      WRITE(6,500) TITRE,REVIS,VERSION,CDATIM
      DO I=1,NSUP
        tempstring=SUPP(I)
        WRITE(6,460) tempstring
      ENDDO
      IF (FLAG.NE.'NON') THEN
         write(6,450) flag
         WRITE(6,700) (IDATIM(I),I=7,14)
      ENDIF
      WRITE(6,800) 'BEGIN  EXECUTION     '
 450  format(3X,'*',107X,'*',/3x,'*',8x,a8,t91,'*')
 460  format(3x,'*',107X,'*',/3x,'*',2x,a105,'*')
 500  FORMAT(1H1,
     %     /,3X,'*',107('*'),'*',
     %     /,3X,'*',107X,'*',
     %     /,3X,'*',3X,A57,10X,A10,27X,'*',
     %     /,3X,'*',107X,'*',
     %     /,3X,'*',1X,A105,1X,'*',
     %     /,3X,'*',107X,'*',
     %     /,3X,'*',107X,'*',
     %     /,3X,'*',3X,A24,46X,34X,'*')
 600  FORMAT(1H1,
     %     /,3X,'*',107('*'),'*',
     %     /,3X,'*',107X,'*',
     %     /,3X,'*',4X,A7,53X,A10,17X,'*',
     %     /,3X,'*',107X,'*',
     %     /,3X,'*',107X,'*',
     %     /,3X,'*',3X,A10,56X,A10,16X,'*')
 700  FORMAT(3X,'*',107X,'*',
     %     /,3X,'*',3X,7A4,I12,10X,'*')
 800  FORMAT(3X,'*',107X,'*',
     %     /,3X,'*',3X,A20,84X,'*',
     %     /,3X,'*',107X,'*',
     %     /,3X,'*',107('*'),'*')

      T1 = SECOND( )
      EXDBPLUS = IDATIM(14)


      RETURN
* ENTREE EXFIN   BOITE DE FIN D'EXECUTION
      ENTRY EXFIN(in_TITRE,REVIS,FLAG)

      call flush_stdout()
      titre = ' '
      titre(1:min(len(in_titre),90)) = in_titre(1:min(len(in_titre),90))
      CALL FDATE(CDATIM)
      WRITE(6,501) TITRE,REVIS,CDATIM,'END EXECUTION       ',
     %     SECOND( )-T1
 501  FORMAT(
     %     /,3X,'*',107('*'),'*',
     %     /,3X,'*',107X,'*',
     %     /,3X,'*',3X,A57,3X,A10,34X,'*',
     %     /,3X,'*',107X,'*',
     %     /,3X,'*',3X,A24,46X,34X,'*',
     %     /,3X,'*',107X,'*',
     %     /,3X,'*',3X,A20,84X,'*',
     %     /,3X,'*',107X,'*',
     %     /,3X,'*',3X,'CP SECS = ',F10.3,84X,'*',
     %     /,3X,'*',107X,'*',
     %     /,3X,'*',107('*'),'*')
 601  FORMAT(
     %     /,3X,'*',107('*'),'*',
     %     /,3X,'*',107X,'*',
     %     /,3X,'*',3X,A7,53X,A10,27X,'*',
     %     /,3X,'*',107X,'*',
     %     /,3X,'*',3X,A10,50X,A10,27X,'*',
     %     /,3X,'*',107X,'*',
     %     /,3X,'*',3X,A20,84X,'*',
     %     /,3X,'*',107X,'*',
     %     /,3X,'*',3X,'CP SECS = ',F10.3,70X,'*',
     %     /,3X,'*',107X,'*',
     %     /,3X,'*',107('*'),'*')

*      call memoirc(0)

      EXFIN = 0
      

      RETURN
      END
