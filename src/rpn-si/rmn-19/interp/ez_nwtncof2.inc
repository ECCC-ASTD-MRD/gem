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
      subroutine ez_nwtncof2(cx,cy,ax,ay,i1,i2,j1,j2)
!*******
!     Auteur: Y. Chartier, drpn
!     Fevrier 1991
!     
!Objet:  Calcul de coefficients uti1ises dans la forme newtonienne
!        de l'interpolation de Lagrange.
!
!*************************************************************
!   -----*-------------*------#------*----------*------->
!        x1            x2     x      x3         x4
!
!*************************************************************
!     cx(i,1) = 1.0 / (x2-x1)
!     cx(i,2) = 1.0 / (x3-x1)
!     cx(i,3) = 1.0 / (x3-x2)
!     cx(i,4) = 1.0 / (x4-x1)
!     cx(i,5) = 1.0 / (x4-x2)
!     cx(i,6) = 1.0 / (x4-x3)
!
!     structure identique pour cy(j,1..6)

!*******
      implicit none
          
      integer ni,nj,i1,j1,i2,j2
      real cx(i1:i2,6),cy(j1:j2,6),ax(i1:i2),ay(j1:j2)
      
      integer i,j
      
      print *, ax
      print *
      print *, ay
      do 10 i=i1+1,i2-2
         cx(i,1) = 1. / (ax(i  ) - ax(i-1))
         cx(i,2) = 1. / (ax(i+1) - ax(i-1))
         cx(i,3) = 1. / (ax(i+1) - ax(i  ))
         cx(i,4) = 1. / (ax(i+2) - ax(i-1))
         cx(i,5) = 1. / (ax(i+2) - ax(i  ))
         cx(i,6) = 1. / (ax(i+2) - ax(i+1))
 10   continue
      
      do 20 j=j1+1,j2-2
         cy(j,1) = 1. / (ay(j  ) - ay(j-1))
         cy(j,2) = 1. / (ay(j+1) - ay(j-1))
         cy(j,3) = 1. / (ay(j+1) - ay(j  ))
         cy(j,4) = 1. / (ay(j+2) - ay(j-1))
         cy(j,5) = 1. / (ay(j+2) - ay(j  ))
         cy(j,6) = 1. / (ay(j+2) - ay(j+1))
         print *, (cy(j,i),i=1,6)
 20   continue
      
      return
      end

!********************************************************************
!**
!**

