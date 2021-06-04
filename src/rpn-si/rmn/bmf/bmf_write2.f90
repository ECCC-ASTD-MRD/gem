!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2005  Environnement Canada
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
function bmf_write2(no_pe,nom,ni,istart,iend,nj,jstart,jend,nk,kstart, &
              kend,time1,time2,hgrid,vgrid,dtyp,scat,ndata,vecteur) result(ierr)
!
!AUTEUR       Luc Corbeil (bmf_write2)
!
!REVISION
!
!ARGUMENTS
!
!______________________________________________________________________ 
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! nom                | nom de la variable a ecrire                     |
! NI                 | 1ere dimension des donnes                       |
! Istart             | indice de debut                                 |
! Iend               | indice de fin                                   |
! NJ                 | 2eme dimension des donnes                       |
! Jstart             | indice de debut                                 |
! Jend               | indice de fin                                   |
! NK                 | 3eme dimension des donnes                       |
! Kstart             | indice de debut                                 |
! Kend               | indice de fin                                   |
! time1              | AAAAMMJJ (annee/mois/jour)                      |
! time2              | HHMMSSCC (heure/minute/secondes/centiemes)      |
! hgrid              | Descripteur de la grille horizontale            |
! vgrid              | Descripteur de la grille verticale              |
! dtyp               | Type de donnes                                  |
! scat               | scatter list                                    |
! vecteur            | tableau des donnees                             |
! ierr               | error status                                    |
!----------------------------------------------------------------------+
!
! DESCRIPTION
!
! Routine which write in the specicified file an array and his specifications.
!
! The array is of size ndata*nb_words where nb_words depends ofthe
! type of variable written.
! Here are the recognized variable types (dtyp):
!        integer  =>  40
!        real*4   =>  41
!        integer*8=>  80  (not fully supported yet)
!        real*8   =>  81
!        complex  =>  82  (not fully supported yet)
! ni,nj,nk, (ijk)start and (ijk)end are attributes which will permit
! to place each slice of data correctly in the field. Timestamps with
! the variable name are used for unicity. At reading time, an array of
! size at least (ni,nj,nk) will be used to place each slice, limited
! by (istart:iend, jstart:jend, kstart:kend).
!
! We must have:
! 1 <= imin < imax <= ni
! 1 <= jmin < jmax <= nj
! 1 <= kmin < kmax <= nk
!
! For writing a field with size from 1-a to nx+a, an appropriate change
! to imin, imax and ni must be performed (here, imin=1, imax=nx+2*a and
! ni=nx+2*a). See documentation of bmf_get in order to retreive such a field.
! Notice: at this point , attributes hgrid, vgrid and scat are not used.
!

use bmf_modsplit

! now standard in f90...

  implicit none

! Entry

  integer, intent(IN)::  no_pe
  character*4, intent(IN) :: nom
  integer, intent(IN) :: ni,istart,iend,nj,jstart,jend,nk,kstart,kend
  integer, intent(IN) :: time1,time2
  integer, intent(IN) :: hgrid,vgrid,dtyp,scat
  integer, intent(IN) :: ndata
  integer, intent(IN) :: vecteur(ndata*(dtyp/10/4))
  integer :: ierr

! Declarations

  integer :: size,head_size,data_size,iun
  real rvecteur

  integer  bmf_connect
  external bmf_connect

  ierr=0

! Check if iun is open

  ierr= bmf_connect(no_pe)
  iun=split_unit(no_pe)
  call bmf_write(iun,nom,ni,istart,iend,nj,jstart,jend,nk,kstart, &
              kend,time1,time2,hgrid,vgrid,dtyp,scat,ndata,vecteur)
end function bmf_write2
