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
      program litht
      integer unf,key1
      integer date0,deet,npas,ni,nj,nk,nbits,datyp,ip1,ip2,ip3
      integer ig1,ig2,ig3,ig4,swa,lng,dltf,ubc
      integer extra1,extra2,extra3
      integer fstinf,fstprm
      external fstinf,fstprm

      character *8 etiket
      character *2 nomvar
      character *1 grtyp, typvar

      unf = 10
      ier=fnom(unf,'/users/dor/armn/mid/data4/tcm/g10/19900914.000000',
     %           'std+rnd',0)
      print *,'Debug fnom=',ier
      ier = fstouv(unf,'RND')
      key1 = fstinf (unf,ni1,nj1,nk1,-1,' ',0,0,0,' ','HT')
      print*, key1
      key1 = fstinf (unf,ni1,nj1,nk1,-1,' ',-1,-1,-1,' ','HT')
      print*, key1

      ier = fstprm(key1,DATEO,DEET, NPAS, NI, NJ, NK, NBITS, DATYP, IP1,
     %         IP2, IP3, TYPVAR, NOMVAR, ETIKET, GRTYP, IG1, IG2, IG3,
     %         IG4, SWA, LNG, DLTF, UBC, EXTRA1, EXTRA2, EXTRA3)
      print *,'Debug apres fstprm ip1,ip2,ip3 =',ip1,ip2,ip3

      ier = fstfrm(10)

      stop
      end
