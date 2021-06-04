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
      Real Function jjulien( DEET,NPAS,IDATE )
*
      Implicit    none
*
      Integer     NPAS,IDATE
      Real        DEET
*
*     Auteur: B.Dugas, RPN - 17 mars 1993.
*
C     Version ajustee des nouvelles fonctions ne
C     pouvant produire le nouveau type de STAMP.
C     Les modifications sont en fait des ajouts
C     de commentaires (desactivation de code).
*
*     Description:
*      Cette fonction calcule le nombre de jour (ordinal) 
*      depuis le debut de l'annee en cours, utilisant la 
*      date d'analyse, le nombre de pas de temps effectues 
*      depuis cette date et la taille du pas de temps.
*
*     Parametres d'entree:
*      DEET  = Taille du pas-de-temps en seconde.
*      NPAS  = Nombre de pas de temps depuis IDATE.
*      IDATE = DATE TIME STAMP CMC valide.
*
*      Si DEET et/ou NPAS sont nuls, le jour correspondant
*      a IDATE est tout de meme retourne.
*
      Integer     jour,mois,annee
      Integer     jdebut,jfin,datim(14),is1,is2
      Real*8      heures
*
      External    datmgp,incdatr,jdatec
*
*----------------------------------------------------------------
***    Calculer le nombre d'heures depuis 
***    le debut de l'integration.
*
      If (     DEET .le. 0
     +    .or. NPAS .le. 0 )                                   Then
          heures = 0.0
      Else
          heures = dble(npas)/( 3600./dble( deet ) )
      End If
*
***    Determiner le date-time-stamp correspondant.
*
      Call incdatr( datim(14), IDATE,  heures )
*
***    Extraire l'annee, le mois, le jour et l'heure correspondante.
*
*     Call datmgp( datim )
      call newdate(datim(14),is1,is2,-3)
*
*     heures = datim(5)
      heures = is2/1000000
*     annee  = datim(4)
      annee  = is1/10000
*     jour   = datim(3)
      jour   = mod(is1,100)
*     mois   = datim(2)
      mois   = mod(is1/100,100)
*
***    Verifier s'il y a aussi des minutes.
*
C     If (datim(6).gt.7) heures = heures + ( datim(6)-8 )/120.0
*
***    Trouver le jour julien du debut de cette 
***    annee de meme que le jour julien final.
*
      Call jdatec( jdebut, annee,  01,   01  )
      Call jdatec( jfin,   annee, mois, jour )
*
***    jjulien est la difference de ces deux jours.
*
      heures  = heures / 24.0
      jjulien = jfin-jdebut+1 + heures
*
      Return
      End
