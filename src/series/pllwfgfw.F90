!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

!/@ RMNLIB - Library of useful routines for C and FORTRAN programming
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
! @/

subroutine pllwfgfw(spd,dir,xlat,xlon,li,lj,grtyp,ig1,ig2,ig3,ig4)
   implicit none
!!!#include <arch_specific.hf>

   integer   li,lj
   real      spd(li,lj), dir(li,lj), xlat(li,lj),xlon(li,lj)
   character grtyp
   integer   ig1,ig2,ig3,ig4

   !@auteur y. chartier - april 94

   !@revision
   ! B.Dugas (dec 2005): ajouter un appel a EZ_CROT et corriger la doc

   !@objet
   !     - passe de vent de grille (composantes u et v)
   !     - a vitesse et direction.
   ! converts from u-v (grid) wind components to
   !                standard meteorological speed and direction

   !@arguments
   !     in/out - spd   - a l'entree contient la composante u
   !                      a la sortie la vitesse.
   !     in/out - dir   - a l'entree contient la composante v
   !                      a la sortie la direction
   !     in     - xlat  - latitudes (vraies) ou sont ces vents
   !     in     - xlon  - longitudes (vraies) ou sont ces vents
   !     in     - li    - premiere dimension des champs spd et dir
   !     in     - lj    - deuxieme dimension des champs spd et dir
   !     in     - grtyp - type de grille (utilise au decodage des igx)
   !     in     - igx   - pour x=1,2,3,4: descripteurs de grille codes
   !-------------------------------------------------------------


   real r(3,3),ri(3,3)

   common /qqqmrot/ r,ri

   integer ier
   real    xlat1,xlon1,xlat2,xlon2
   real    xlatgf,xlongf,uvcart,xyz

   pointer( xlatgp  , xlatgf(li,*) )
   pointer( xlongp  , xlongf(li,*) )
   pointer( uvcartp , uvcart( 3,*) )
   pointer( xyzp    ,    xyz( 3,*) )

   call hpalloc( xlatgp,   li*lj,ier,0 )
   call hpalloc( xlongp,   li*lj,ier,0 )
   call hpalloc( uvcartp,3*li*lj,ier,0 )
   call hpalloc( xyzp,   3*li*lj,ier,0 )

   call cigaxg(grtyp, xlat1, xlon1, xlat2, xlon2, &
        ig1  , ig2  , ig3  , ig4)
   call ez_crot(r,ri,  xlon1, xlat1, xlon2, xlat2)

   call ez_gfxyfll(xlon,xlat,xlongf,xlatgf, li*lj, &
        xlat1,xlon1, xlat2, xlon2)

   call ez_vrotf2(spd,dir, xlon,xlat,xlongf,xlatgf, &
        ri,xyz,uvcart,li,lj)

   call ez_llwfgdw(spd,dir, xlongf,li,lj, 'L',0,0,0,0)

   call hpdeallc(xlatgp ,ier,0)
   call hpdeallc(xlongp ,ier,0)
   call hpdeallc(uvcartp,ier,0)
   call hpdeallc(xyzp   ,ier,0)

   return
end subroutine pllwfgfw
