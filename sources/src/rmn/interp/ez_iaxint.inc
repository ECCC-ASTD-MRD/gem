!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
      subroutine ez_iaxint(zo,px,npts,ax,z,ni, i1,i2,ordint)
      implicit none
!*******
!Auteur: Y.Chartier, drpn
!        Fevrier 1991
!
!Objet:  Interpolation bi-cubique de points a partir 
!        d'une grille source irreguliere.
!*******

      integer npts,ni,i1,i2,ordint
      real zo(npts),px(npts)
      real ax(i1:i2)
      real z(ni)
!
!  npts   : nombre de points a interpoler
!  i1:i2  : dimension de la grille source selon x
!  zo     : vecteur de sortie contenant les valeurs interpolees
!  px     : vecteur contenant la position x des points que l'on
!           veut interpoler
!  ax     : vecteur contenant la pos. des points sur l'axe des X.
!  cx     : vecteur contenant une table de differences selon X.
!  z      : valeurs de la grille source.
!
!***************************************************************************
!
!  *   *   *   *
!  
!  *   *   *   *
!        #        .eq.>   pt (x,y)
!  *  (=)  *   *  .eq.> = pt (i, j)
!
!  *   *   *   *
!
!
!     

!  structure identique pour cy(j,1..6)

      real a11,a12,a13,a14
      real b1
      real x1,x2,x3
      integer i, n
      
      real a1,a2,a3,a4,x,c1,c2,c3,c4,c5,c6
      real fa, fa2, fa3, fa4
      real dx,dy
      
!     definition des fonctions in-line

#include "ez_qqqxtrp.cdk"
#include "zlin.cdk"

      fa(a1,a2,a3,a4,x,x1,x2,x3)=a1+(x-x1)*(a2+(x-x2)*(a3+a4*(x-x3)))
      fa2(c1,a1,a2)=c1*(a2-a1)
      fa3(c1,c2,c3,a1,a2,a3)=c2*(c3*(a3-a2)-c1*(a2-a1))
      fa4(c1,c2,c3,c4,c5,c6,a1,a2,a3,a4)=c4*(c5*(c6*(a4-a3)-      c3*(a3-a2)) - c2*(c3*(a3-a2)-c1*(a2-a1)))

      if (ordint.eq.cubique) then
         do 10 n=1,npts
            i = min(ni-2,max(2,ifix(px(n))))
            
            x = ax(i) + (ax(i+1)-ax(i))*(px(n)-i)
            
            a11 = z(i-1)
            a12 = fa2((1.0/(ax(i)-ax(i-1))),z(i-1),z(i))
            a13 = fa3((1.0/(ax(i)-ax(i-1))),1.0/(ax(i+1)-ax(i-1)), 1.0/(ax(i+1)-ax(i)),z(i-1),z(i),z(i+1))
            a14 = fa4((1.0/(ax(i)-ax(i-1))),1.0/(ax(i+1)-ax(i-1)), 1.0/(ax(i+1)-ax(i)),1.0/(ax(i+2)-ax(i-1)),1.0/(ax(i+2)-ax(i)),1.0/(ax(i+2)-ax(i+1)),z(i-1),z(i),z(i+1),z(i+2))
            zo(n)  = fa(a11,a12,a13,a14,x,ax(i-1),ax(i),ax(i+1))
            
            
 10      continue
         endif

         if (ordint.eq.lineair) then
            do 20 n=1,npts
               i = min(i2-1,max(i1,ifix(px(n))))
               
               x = ax(i) + (ax(i+1)-ax(i))*(px(n)-i)
               
               dx = (x - ax(i))/(ax(i+1)-ax(i))
               
               zo(n) = zlin((z(i)),z(i+1),dx)
 20         continue
         endif
         
         if (ordint.eq.voisin) then
            do 30 n=1,npts
               i = min(i2,max(i1,nint(px(n))))
               zo(n) = z(i)
 30         continue
         endif
         
         return
         end

