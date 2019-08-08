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
      subroutine ez_gggdint_nw(zo,px,py,npts,ay,z,i1,i2,j1,j2)
      implicit none
!*******
!Auteur: Y.Chartier, drpn
!        Fevrier 1991
!
!Objet:  Interpolation bi-cubique de points a partir 
!        d'une grille gaussienne.
!*******

      integer npts,i1,i2,j1,j2,degree
      real zo(npts),px(npts),py(npts)
      real ay(j1:j2),cy(6)
      real z(i1:i2,j1:j2)
!
!  npts   : nombre de points a interpoler
!  i1:i2  : dimension de la grille source selon x
!  j1:j2  : dimension de la grille source selon y
!  zo     : vecteur de sortie contenant les valeurs interpolees
!  px     : vecteur contenant la position x des points que l'on
!           veut interpoler
!  py     : vecteur contenant la position y des points que l'on
!           veut interpoler
!  ay     : vecteur contenant la pos. des points sur l'axe des Y.
!  cy     : vecteur contenant une table de differences selon Y.
!  z      : valeurs de la grille source.
!
!**********************************************************************
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
      integer i, j, m, n,stride
      real*8 x, x1, x2, x3, x4
      real*8 b1, b2,  b3,  b4
      real*8 b11, b12, b13, b14
      real*8 y,y1,y2,y3,y4
      real*8 y11, y12, y13, y14
      real*8 ay1, ay2, ay3, ay4
      real*8 fa, fa2, fa3, fa4
      real*8 a1,a2,a3,a4,c1,c2,c3,c4,c5,c6
      
#include "cubic8.cdk"

      fa(a1,a2,a3,a4,x,x1,x2,x3)=a1+(x-x1)*(a2+(x-x2)*(a3+a4*(x-x3)))
      fa2(c1,a1,a2)=c1*(a2-a1)
      fa3(c1,c2,c3,a1,a2,a3)=c2*(c3*(a3-a2)-c1*(a2-a1))
      fa4(c1,c2,c3,c4,c5,c6,a1,a2,a3,a4)=c4*(c5*(c6*(a4-a3)-c3*      (a3-a2)) - c2*(c3*(a3-a2)-c1*(a2-a1)))
      
      do n=1,npts
         i = min(i2-2,max(i1+1,ifix(px(n))))
         j = min(j2-2,max(j1+1,ifix(py(n))))
         
         dx = px(n) - i
         
         y1=cubic(dble(z(i-1,j-1)),dble(z(i  ,j-1)),dble(z(i+1,j-1)),dble(z(i+2,j-1)),dx)
         y2=cubic(dble(z(i-1,j  )),dble(z(i  ,j  )),dble(z(i+1,j  )),dble(z(i+2,j  )),dx)
         y3=cubic(dble(z(i-1,j+1)),dble(z(i  ,j+1)),dble(z(i+1,j+1)),dble(z(i+2,j+1)),dx)
         y4=cubic(dble(z(i-1,j+2)),dble(z(i  ,j+2)),dble(z(i+1,j+2)),dble(z(i+2,j+2)),dx)
         
         y = ay(j) + (ay(j+1)-ay(j))*(py(n)-j)
          
!     interpolation finale selon y
         
         ay1=ay(j-1)
         ay2=ay(j)
         ay3=ay(j+1)
         ay4=ay(j+2)
         
         cy(1) = 1.0 / (ay2 - ay1)
         cy(2) = 1.0 / (ay3 - ay1)
         cy(3) = 1.0 / (ay3 - ay2)
         cy(4) = 1.0 / (ay4 - ay1)
         cy(5) = 1.0 / (ay4 - ay2)
         cy(6) = 1.0 / (ay4 - ay3)

         y11 = y1
         y12 = fa2(dble(cy(1)),y1,y2)
         y13 = fa3(dble(cy(1)),dble(cy(2)),dble(cy(3)),y1,y2,y3)
         y14 = fa4(dble(cy(1)),dble(cy(2)),dble(cy(3)),dble(cy(4)),         dble(cy(5)),dble(cy(6)),y1,y2,y3,y4)
         zo(n) = fa(y11,y12,y13,y14,y,ay1,ay2,ay3)
      enddo
      
      return
      end
      
