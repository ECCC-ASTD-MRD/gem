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

       function qvmcks(adresse,nbelem,mode)
      
        integer qvmcks
        integer nbelem, mode
        integer adresse(nbelem)
***s/p qvmcks
*
*objet(qvmcks)
*          Calculer et retourner le check-sum du bloc
*          memoire associe a la clef inkey
*
*auteur J. Caveen  -  novembre 1993
*
*arguments
*         in  adresse-  adresse memoire du premier element du tableau
*         in  nbelem -  nombre d'elements pour lesquels on fait le check-sum
*         in  mode   -  mode de calcul
*                       pour cette version, mode=1 seulement => xor du champ
*
**
          integer cks, i

          qvmcks = -1

          cks = 0
*
*         check - sum en mode 1 : summ de toutes les entrees de adresse
*
          if (mode .eq. 1) then
              do i = 1, nbelem
                  cks = cks + adresse(i)
              enddo
          endif

          qvmcks = cks
          return
          end

