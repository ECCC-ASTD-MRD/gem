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
#if !defined(NEC)
      subroutine qqqfasp(sc,cfcs,cfca,lp,lr,na,lmc,fcd2,lpd1,n,lb)
      integer lr,na,lmc,fcd2,lpd1,n,lb
      real sc(na,n),cfcs(lmc,fcd2,lb),cfca(lmc,fcd2,lb)
      real*8 lp(lpd1,lb),tempo

      integer mn,k,ilb
      
*     print *,' lr,na,lmc,fcd2,lpd1,n,lb = ',lr,na,lmc,fcd2,lpd1,n,lb
      do 21 mn=1,lr,2
      do 21 k=1,n
      tempo = sc(mn,k)
      do 20 ilb=1,lb
20    tempo=tempo+lp(mn,ilb)*cfcs(1,k,ilb)
21    sc(mn,k) = tempo
      do 41 mn=2,lr,2
      do 41 k=1,n
      tempo = sc(mn,k)
      do 40 ilb=1,lb
40    tempo=tempo+lp(mn,ilb)*cfca(1,k,ilb)
41    sc(mn,k) = tempo

      return

      entry qqqspaf(sc,cfcs,cfca,lp,lr,na,lmc,fcd2,lpd1,n,lb)
      
*     print *,' lr,na,lmc,fcd2,lpd1,n,lb = ',lr,na,lmc,fcd2,lpd1,n,lb
      do 50 ilb=1,lb
      do 50 k=1,n
         cfcs(1,k,ilb) = 0.0
      	 cfca(1,k,ilb) = 0.0
50    continue

      do 61 ilb=1,lb
      do 61 k=1,n
      tempo = 0
      do 60 mn=1,lr,2
60    tempo=tempo+lp(mn,ilb)*sc(mn,k)
      cfcs(1,k,ilb)=tempo
61    continue

      do 71 ilb=1,lb
      do 71 k=1,n
      tempo = 0
      do 70 mn=2,lr,2
70    tempo=tempo+lp(mn,ilb)*sc(mn,k)
      cfca(1,k,ilb)=tempo
71    continue

      return
      end
#endif
