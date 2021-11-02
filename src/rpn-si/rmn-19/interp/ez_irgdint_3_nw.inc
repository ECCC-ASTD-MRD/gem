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
      subroutine ez_irgdint_3_nw(zo,px,py,npts,ax,ay,cx,cy,z,i1,i2,j1,j2)
      implicit none
!*******
!     Auteur: Y.Chartier, drpn
!     Fevrier 1991
!     
!     Objet:  Interpolation bi-cubique de points a partir 
!     d'une grille source irreguliere.
!*******
      
      integer npts,i1,i2,j1,j2,degree
      real zo(npts),px(npts),py(npts)
      real fa, fa2, fa3, fa4
      real ax(i1:i2),ay(j1:j2),cx(i1:i2,6),cy(j1:j2,6)
      real z(i1:i2,j1:j2)
!     
!     npts   : nombre de points a interpoler
!     1:ni  : dimension de la grille source selon x
!     1:nj  : dimension de la grille source selon y
!     zo     : vecteur de sortie contenant les valeurs interpolees
!     px     : vecteur contenant la position x des points que l'on
!              veut interpoler
!     py     : vecteur contenant la position y des points que l'on
!              veut interpoler
!     ax     : vecteur contenant la pos. des points sur l'axe des X.
!     ay     : vecteur contenant la pos. des points sur l'axe des Y.
!     cx     : vecteur contenant une table de differences selon X.
!     cy     : vecteur contenant une table de differences selon Y.
!     z      : valeurs de la grille source.
!     
!***************************************************************************
!     
!     *   *   *   *
!     
!     *   *   *   *
!           #        .eq.>   pt (x,y)
!     *  (=)  *   *  .eq.> = pt (i, j)
!     
!     *   *   *   *
!     
!     
!     
!     cx(i,1) = 1.0 / (x2-x1)
!     cx(i,2) = 1.0 / (x3-x1)
!     cx(i,3) = 1.0 / (x3-x2)
!     cx(i,4) = 1.0 / (x4-x1)
!     cx(i,5) = 1.0 / (x4-x2)
!     cx(i,6) = 1.0 / (x4-x3)
!     
!  structure identique pour cy(j,1..6)

      integer i, j, m, n,stride
      real*8 a11,a12,a13,a14,a21,a22,a23,a24
      real*8 a31,a32,a33,a34,a41,a42,a43,a44
      real*8 b1,b2,b3,b4,b11,b12,b13,b14
      real*8 x1,x2,x3,x4,y1,y2,y3,y4
      real*8 a1,a2,a3,a4,x,y,c1,c2,c3,c4,c5,c6
      real*8 z1,z2,z3,z4
      
#include "ez_qqqxtrp.cdk"

!     definition des fonctions in-line


      fa(a1,a2,a3,a4,x,x1,x2,x3)=a1+(x-x1)*(a2+(x-x2)*(a3+a4*(x-x3)))
      fa2(c1,a1,a2)=c1*(a2-a1)
      fa3(c1,c2,c3,a1,a2,a3)=c2*(c3*(a3-a2)-c1*(a2-a1))
      fa4(c1,c2,c3,c4,c5,c6,a1,a2,a3,a4)=c4*(c5*(c6*(a4-a3)-      c3*(a3-a2)) - c2*(c3*(a3-a2)-c1*(a2-a1)))
      
      do n=1,npts
         i = min(i2-2,max(i1+1,ifix(px(n))))
         j = min(j2-2,max(j1+1,ifix(py(n))))
         
         x = ax(i) + (ax(i+1)-ax(i))*(px(n)-i)
         y = ay(j) + (ay(j+1)-ay(j))*(py(n)-j)
         
         x1=ax(i-1)
         x2=ax(i)
         x3=ax(i+1)
         x4=ax(i+2)
         
         y1=ay(j-1)
         y2=ay(j)
         y3=ay(j+1)
         y4=ay(j+2)
         
!     interpolation 1ere rangee selon x
         
         z1=z(i-1,j-1)
         z2=z(i  ,j-1)
         z3=z(i+1,j-1)
         z4=z(i+2,j-1)
         
         a11 = z1
         a12 = fa2(dble(cx(i,1)),z1,z2)
         a13 = fa3(dble(cx(i,1)),dble(cx(i,2)),dble(cx(i,3)),z1,z2,z3)
         a14 = fa4(dble(cx(i,1)),dble(cx(i,2)),dble(cx(i,3)),dble(cx(i,4)),dble(cx(i,5)),dble(cx(i,6)),         z1,z2,z3,z4)
         b1  = fa(a11,a12,a13,a14,x,x1,x2,x3)
         
!     interpolation 2eme rangee selon x
         
         z1=z(i-1,j)
         z2=z(i  ,j)
         z3=z(i+1,j)
         z4=z(i+2,j)
         
         a21 = z1
         a22 = fa2(dble(cx(i,1)),z1,z2)
         a23 = fa3(dble(cx(i,1)),dble(cx(i,2)),dble(cx(i,3)),z1,z2,z3)
         a24 = fa4(dble(cx(i,1)),dble(cx(i,2)),dble(cx(i,3)),dble(cx(i,4)),         dble(cx(i,5)),dble(cx(i,6)),z1,z2,z3,z4)
         b2  = fa(a21,a22,a23,a24,x,x1,x2,x3)
         
!     interpolation 3eme rangee selon x
         
         z1=z(i-1,j+1)
         z2=z(i  ,j+1)
         z3=z(i+1,j+1)
         z4=z(i+2,j+1)
         
         a31 = z1
         a32 = fa2(dble(cx(i,1)),z1,z2)
         a33 = fa3(dble(cx(i,1)),dble(cx(i,2)),dble(cx(i,3)),z1,z2,z3)
         a34 = fa4(dble(cx(i,1)),dble(cx(i,2)),dble(cx(i,3)),dble(cx(i,4)),         dble(cx(i,5)),dble(cx(i,6)),z1,z2,z3,z4)
         b3  = fa(a31,a32,a33,a34,x,x1,x2,x3)
         
!     interpolation 4eme rangee selon x
         
         z1=z(i-1,j+2)
         z2=z(i  ,j+2)
         z3=z(i+1,j+2)
         z4=z(i+2,j+2)
         
         a41 = z1
         a42 = fa2(dble(cx(i,1)),z1,z2)
         a43 = fa3(dble(cx(i,1)),dble(cx(i,2)),dble(cx(i,3)),z1,z2,z3)
         a44 = fa4(dble(cx(i,1)),dble(cx(i,2)),dble(cx(i,3)),dble(cx(i,4)),         dble(cx(i,5)),dble(cx(i,6)),z1,z2,z3,z4)
         b4  = fa(a41,a42,a43,a44,x,x1,x2,x3)
         
!     interpolation finale selon y
         
         b11 = b1
         b12 = fa2(dble(cy(j,1)),b1,b2)
         b13 = fa3(dble(cy(j,1)),dble(cy(j,2)),dble(cy(j,3)),b1,b2,b3)
         b14 = fa4(dble(cy(j,1)),dble(cy(j,2)),dble(cy(j,3)),dble(cy(j,4)),         dble(cy(j,5)),dble(cy(j,6)),b1,b2,b3,b4)
         zo(n) = fa(b11,b12,b13,b14,y,y1,y2,y3)
      enddo
      
      return
      end
