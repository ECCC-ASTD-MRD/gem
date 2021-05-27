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
subroutine bmf_write(iun,nom,ni,istart,iend,nj,jstart,jend,nk,kstart, &
              kend,time1,time2,hgrid,vgrid,dtyp,scat,ndata,vecteur)
!
!AUTEUR       Luc Corbeil (bmf_write)
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
!----------------------------------------------------------------------+
!
! DESCRIPTION
!
! Routine which write in the specicified file an array and his specifications.
! We strongly recommend that the file is opened as following: 
! (routine FNOM, part of rmnlib)  
!
!  ierr = FNOM(iun,filename,'SEQ/UNF',0)
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
! now standard in f90...
 use bmf_mod
  implicit none
! Entry
  integer::  iun
  character*4 :: nom
  integer :: ni,istart,iend,nj,jstart,jend,nk,kstart,kend
  integer :: time1,time2
  integer :: hgrid,vgrid,dtyp,scat
  integer :: ndata
  integer :: vecteur(max(ndata*(dtyp/10/4),1))
! Declarations
  integer :: size,head_size,data_size,sizeofint,toto(2)
  integer, allocatable, dimension(:) :: cdata
  real rvecteur
  integer bmf_char2i, err
!
  size=dtyp/10/4
  head_size=4*(13)+4
  data_size=ndata*size
!
  sizeofint=loc(toto(2))-loc(toto(1))
  if(dtyp.eq.bmf_character) then
     allocate(cdata(2+(ndata-1)/sizeofint))
     err = bmf_char2i(vecteur,ndata,cdata,2+(ndata-1)/sizeofint)
  endif

  write(iun) head_size,nom
  write(iun) ni,istart,iend
  write(iun) nj,jstart,jend
  write(iun) nk,kstart,kend
  write(iun) time1,time2
  write(iun) hgrid,vgrid
  write(iun) dtyp,scat
  write(iun) ndata
  write(iun) head_size
  write(iun) data_size
  if(dtyp.eq.bmf_character) then
     write(iun) cdata
  else
     write(iun) vecteur  
  endif
  write(iun) data_size
end subroutine bmf_write
