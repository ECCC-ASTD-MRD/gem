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
!**s/r ez_vrotf - rotation of the components of the wind once interpolated
!              to the rotated grid of the model
!
      subroutine ez_vrotf(u, v, lonp, latp, lon, lat, ro, xyz, uvcart, ni, nj)
      implicit none
      integer ni, nj 
      real    u(ni,nj), lonp(ni,nj), lon(ni,nj), uvcart(3,ni*nj),          ro(3,3), v(ni,nj),  latp(ni,nj), lat(ni,nj), xyz(3,ni*nj)
!
!author michel roch - april 1990
!
!revision
!	001 - yvon bourrassa - mai/juin 1993 - remove dynamic allocation 
!	                       change calling sequence
!	002 - michel roch - documentation
!
!arguments
!   in/out  u         composante u du vent sur grille non tournee en entree
!                     composante u du vent sur grille tournee en sortie
!           v         composante v du vent sur grille non tournee en entree
!                     composante v du vent sur grille tournee en sortie
!    in     lonp      longitudes d'origine dans le systeme non tourne 
!           latp      latitudes d'origine dans le systeme non tourne
!           lon       longitudes de grille variable dans le systeme tourne
!           lat       latitudes de grille variable dans le systeme tourne
!           ro        matrice de transformation du systeme non tourne 
!                     au systeme de coordonnees tourne
!           ni        dimension e-o de la grille a sortir
!           nj        dimension n-s de la grille a sortir
!    out    uvcart    champ de travail
!           xyz       champ de travail
!
!
      external ez_uvacart, mxm, ez_cartauv

!     calcul des vent en espace cartesiennes
      call  ez_uvacart(xyz, u, v, lonp, latp, ni, nj)

!     calcul des vents dans l'espace cartesien avec rotation
      call mxm(ro, 3, xyz, 3, uvcart, ni*nj)

!     calcul des composantes de vent dans le systeme de
!                coordonnees tourne par rapport a la geographie reelle
      call ez_cartauv(u, v, uvcart, lon, lat, ni, nj)

      return
      end 
