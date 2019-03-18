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
subroutine bmf_catalog(nom,ni,istart,iend,nj,jstart,jend,nk, &   
                      kstart,kend,time1,time2,hgrid,vgrid,dtyp,scat,ndata)
!
!AUTEUR       Luc Corbeil (bmf_catalog)
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
! ndata              | nombre de blocs                                 |
!----------------------------------------------------------------------+
!
!DESCRIPTION
!
! Subroutine which lists all fields read with bmf_gobe and their attributes.
! Those fields can be retrieved via bmf_get. The main goal of this subroutine
! is to find the size needed to stock a field in an array. 
!
  use bmf_mod
! now standard in f90...
  implicit none
        integer :: length
        character*4 :: nom(bmf_length)
        integer :: ni(bmf_length),nj(bmf_length),nk(bmf_length)
        integer :: istart(bmf_length),iend(bmf_length),jstart(bmf_length)
        integer :: jend(bmf_length),kstart(bmf_length),kend(bmf_length)     
        integer :: ndata(bmf_length),time1(bmf_length),time2(bmf_length)
        integer :: hgrid(bmf_length),vgrid(bmf_length),dtyp(bmf_length)
        integer :: scat(bmf_length) 
        type(bmf_liste), pointer :: champ_courant
        integer :: indice,i,ttemp,ttemp2,icourant
        character*4 nom_tempo
        logical trouve
!      write(*,*) 'DEBUT'
   length=bmf_length
   champ_courant=> liste
   indice=1
   icourant=1
   do i=1,length
      ndata(i)=0
   enddo
   do
      if(.not.associated(champ_courant)) EXIT
      nom_tempo=champ_courant%bmf_champ%nom
      ttemp=champ_courant%bmf_champ%time1
      ttemp2=champ_courant%bmf_champ%time2
      trouve=.false.
      do i=1,indice-1
         if((nom(i).eq.nom_tempo).and.(ttemp.eq.time1(i)) &
           .and.(ttemp2.eq.time2(i))) then
            icourant=i
            trouve=.true.
         else
            icourant=indice
         endif
      enddo
      nom(icourant)=nom_tempo
!   write(*,*) bmf_length, 'GLAM1', g_lam, nom(icourant)
      ni(icourant)=champ_courant%bmf_champ%ni
      nj(icourant)=champ_courant%bmf_champ%nj
      nk(icourant)=champ_courant%bmf_champ%nk
      istart(icourant)=champ_courant%bmf_champ%istart
      iend(icourant)=champ_courant%bmf_champ%iend
!    write(*,*) bmf_length, 'GLAM1', g_lam, iend(icourant)
      jstart(icourant)=champ_courant%bmf_champ%jstart
      jend(icourant)=champ_courant%bmf_champ%jend
      kstart(icourant)=champ_courant%bmf_champ%kstart
      kend(icourant)=champ_courant%bmf_champ%kend
!   write(*,*) bmf_length, 'GLAM1', g_lam, kend(icourant),icourant,loc(G_lam)
      dtyp(icourant)=champ_courant%bmf_champ%dtyp
!   write(*,*) bmf_length, 'GLAM1', g_lam, dtyp(icourant),loc(dtyp(icourant))
      ndata(icourant)=champ_courant%bmf_champ%ndata
!   write(*,*) bmf_length, 'GLAM1ndata', g_lam, ndata(icourant)
      time1(icourant)=ttemp
      time2(icourant)=ttemp2
      hgrid(icourant)=champ_courant%bmf_champ%hgrid
      vgrid(icourant)=champ_courant%bmf_champ%vgrid
      scat(icourant)=champ_courant%bmf_champ%scat
     if(indice.gt.length) then
         write(*,*) 'BMF_CATALOG: number of fields found .gt. expected'
         call abort
      endif

      if(.not.trouve) indice=indice+1
      champ_courant=>champ_courant%champ_suivant
   enddo
return
end subroutine bmf_catalog
