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
      subroutine mxm(a,nar,b,nac,c,nbc)
C
C    Routine on VAX to emulate CRAY SCILIB routine, MXM p. 4-22
C         LIBRARY manual, implemented in VAX single precision.
C
      integer nar,nac,nbc,i,j,k
      real a(nar,1),b(nac,1),c(nar,1)

      do 30 j=1,nbc
         do 20 i = 1,nar
            c(i,j) = 0.0
            do 10 k = 1,nac
               c(i,j) = c(i,j) + a(i,k)*b(k,j)
   10       continue
   20    continue
   30 continue
      return
      end
