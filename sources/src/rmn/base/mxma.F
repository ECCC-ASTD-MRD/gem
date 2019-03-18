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
#if !defined (NEC)
      subroutine mxma(a,na,iad,b,nb,ibd,c,nc,icd,nar,nac,nbc)
C
C     Computes a matrix, A, times a matrix, B, and allows for
C          arbitrary spacing of matrix elements.
C     Martin J. McBride.  8/21/85.
C     General Electric CRD, Information System Operation.
C
      integer na,iad,nb,ibd,nc,icd,nar,nac,nbc
      integer ac,ar,bc,br,cc,cr,i,j,k
      real a(1,1),b(1,1),c(1,1)

      bc = 1
      cc = 1
      do 30 j = 1,nbc
         cr = 1
         ar = 1
         do 20 i = 1,nar
            c(cr,cc) = 0.0
            ac = 1
            br = 1
            do 10 k = 1,nac
               c(cr,cc) = c(cr,cc) + a(ar,ac)*b(br,bc)
               ac = ac + iad
               br = br + nb
   10       continue
            cr = cr + nc
            ar = ar + na
   20    continue
         cc = cc + icd
         bc = bc + ibd
   30 continue
      return
      end

      subroutine mxma8(a,na,iad,b,nb,ibd,c,nc,icd,nar,nac,nbc)
C
C     Computes a matrix, A, times a matrix, B, and allows for
C          arbitrary spacing of matrix elements.
C     Martin J. McBride.  8/21/85.
C     General Electric CRD, Information System Operation.
C
      integer na,iad,nb,ibd,nc,icd,nar,nac,nbc
      integer ac,ar,bc,br,cc,cr,i,j,k
      real*8 a(1,1),b(1,1),c(1,1)

      bc = 1
      cc = 1
      do 30 j = 1,nbc
         cr = 1
         ar = 1
         do 20 i = 1,nar
            c(cr,cc) = 0.0
            ac = 1
            br = 1
            do 10 k = 1,nac
               c(cr,cc) = c(cr,cc) + a(ar,ac)*b(br,bc)
               ac = ac + iad
               br = br + nb
   10       continue
            cr = cr + nc
            ar = ar + na
   20    continue
         cc = cc + icd
         bc = bc + ibd
   30 continue
      return
      end
#endif
