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
      subroutine ez_rgdint_3_wnnc(zo,px,py,npts,z,ni,j1,j2,wrap)
!*******
!Auteur: Y.Chartier, drpn
!        Fevrier 1991
!
!Objet:  Interpolation bi-cubique de points a partir d'une grille
!        source reguliere.
!        
!*******
      implicit none


      integer npts,ni,j1,j2,degree,wrap
      real zo(npts),px(npts),py(npts)
      real z(ni,j1:j2)
!
!  npts   : nombre de points a interpoler
!  i1:i2  : dimension de la grille source selon x
!  j1:nj  : dimension de la grille source selon y
!  zo     : vecteur de sortie contenant les valeurs interpolees
!  px     : vecteur contenant la position x des points que l'on
!         : veut interpoler
!  py     : vecteur contenant la position y des points que l'on
!         : veut interpoler
!  z      : valeurs de la grille source.
!
!  wrap est est le facteur de "wrap around" dans le cadre d'une grille globale
!  pour une grille de type 'A' ou 'G', wrap = 2
!  pour une grille de type 'B', wrap = 1
!  dans tous les autres case wrap = 0
!
!===========================================
!
!     *   *   *   *
!     
!     *   *   *   *
!           #        ==>   pt (x,y)
!     *  (=)  *   *  ==> = pt (iind, jind)
!
!     *   *   *   *
!
!===========================================

      real*8 y1,y2,y3,y4
      integer m,n,i,j,stride
      integer imoins1, iplus1, iplus2, limite

#include "ez_qqqxtrp.cdk"
#include "cubic8.cdk"

      limite = ni+2-wrap
      
      do n=1,npts
         i = min(ni-2+wrap,max(1,max(2-wrap,ifix(px(n)))))
         j = min(j2-2,max(j1+1,ifix(py(n))))
         
         if (wrap.gt.0.and.(i.le.1).or.i.ge.(ni-1)) then
            imoins1 = mod(limite+i-1,limite)
            iplus1  = mod(limite+i+1,limite)
            iplus2  = mod(limite+i+2,limite)
            
            if (imoins1.eq.0) imoins1 = ni
            if (i.eq.0) i = ni
            if (iplus1.eq.0) iplus1 = ni
            if (iplus2.eq.0) iplus2 = ni
            
            if (wrap.eq.1) then
               if (iplus2.eq.ni) iplus2 = 2
               if (imoins1.eq.ni) imoins1=ni-1
            endif
         else
            imoins1 = i-1
            iplus1  = i+1
            iplus2  = i+2
         endif
         dx = px(n) - i
         dy = py(n) - j
         
         y1=cubic(dble(z(imoins1,j-1)),dble(z(i ,j-1)),dble(z(iplus1,j-1)),dble(z(iplus2,j-1)),dx)
         y2=cubic(dble(z(imoins1,j  )),dble(z(i ,j  )),dble(z(iplus1,j  )),dble(z(iplus2,j  )),dx)
         y3=cubic(dble(z(imoins1,j+1)),dble(z(i ,j+1)),dble(z(iplus1,j+1)),dble(z(iplus2,j+1)),dx)
         y4=cubic(dble(z(imoins1,j+2)),dble(z(i ,j+2)),dble(z(iplus1,j+2)),dble(z(iplus2,j+2)),dx)
         
         zo(n)=cubic(y1,y2,y3,y4,dy)
      enddo
      
      return
      end

