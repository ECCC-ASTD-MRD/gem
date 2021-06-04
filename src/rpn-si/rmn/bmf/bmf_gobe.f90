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
function bmf_gobe(filename) result(length)
!
!FUNCTION BMF_gobe(filename) result(length)  L. Corbeil
!
! ARGUMENTS
! filename     name of BMF file                   character*
! length       number of distinct variables read  integer
!
! DESCRIPTION
! Routine that reads the BMF file specified by "filename" 
! (previously written with bmf_write or bmf_splitwr*) and 
! keeps all the data in memory for future use by bmf_get.
! The "length" variable returned by the function corresponds
! to the number of distincts fields read. Two fields are 
! not distinct if one of "time1", "time2" and "nom" are
! the same for both (see bmf_write,  bmf_catalog and bmf_get).
! Furthermore, some checks are performed in order to assure
! that  ni, nj, nk and dtyp are the same for non-distincts 
! fields for data consistency. It is possible to perform
! many bmf_gobe between one bmf_init and bmf_clear. In this
! case, the return value "length" will be the number of new
! distinct fields found.
! 
! revision Aug.18,2005 V.Lee (removed extra allocate tab_temp)

  use bmf_mod

! now standard in f90...
  implicit none

  character*(*) filename
  integer :: length
  integer :: nbyt_head,nbyt_head2,nbyt_data,nbyt_data2

! Declarations

  type(bmf_liste),pointer :: champ_courant,champ_recherche
  integer :: fclos,fnom,longueur
  external fclos,fnom,longueur
  integer :: iun,ierr,status,i,j,k,size
  logical fichier,fin_liste,trouve
  character*4 nom_tempo
  integer ni,istart,iend,nj,jstart,jend,nk,kstart,kend
  integer time1,time2,hgrid,vgrid,dtyp,scat,ndata,indice
  integer sizeofint, toto(2)
  integer, dimension(:), pointer :: tab_temp

  fichier=.true.
  sizeofint = loc(toto(2))-loc(toto(1))
  nullify(champ_recherche)

! Is BMF  mode started?
  if(.not.bmf_started) STOP 'BMF_GOBE: BMF mode not started: use bmf_init'
  iun=0
  ierr = fnom(iun,filename(1:longueur(filename)),'SEQ/UNF',0)

  if(ierr.ne.0) then
     write(*,*)  'BMF_GOBE: error while opening ',filename
     length=-1
     return
  endif

  read(iun,err = 100, end = 9100 ) nbyt_head,nom_tempo
  ! Do we have a liste already?

  if(.not.bmf_liste_started) then
     nullify(liste)
     bmf_liste_started=.true.
     bmf_length=0      
  endif

! Read stuff in temp variables, we'll find later where we put them

  do while(fichier)

     read(iun,err = 100, end = 9100 ) ni,istart,iend
     read(iun,err = 100, end = 9100 ) nj,jstart,jend
     read(iun,err = 100, end = 9100 ) nk,kstart,kend
     read(iun,err = 100, end = 9100 ) time1,time2
     read(iun,err = 100, end = 9100 ) hgrid,vgrid
     read(iun,err = 100, end = 9100 ) dtyp, scat
     read(iun,err = 100, end = 9100 ) ndata
     read(iun,err = 100, end = 9100 ) nbyt_head2
     if(nbyt_head.ne.nbyt_head2) STOP 'probleme avec en-tete, nom=' 
     size=max(dtyp/10/4,1)

! check if we have a similar record already

     trouve=.false.	
     if(bmf_length.gt.0) then
        trouve=.false.
        fin_liste=.false.
        trouve=.false.
        champ_recherche=>liste
        do while((.not.fin_liste).and.(.not.trouve))
           ! Do we have a match?
           if((champ_recherche%bmf_champ%nom.eq.nom_tempo).and. &
                (champ_recherche%bmf_champ%time1.eq.time1).and. &
                (champ_recherche%bmf_champ%time2.eq.time2).and. &
                (champ_recherche%bmf_champ%vgrid.eq.vgrid).and. &
                (champ_recherche%bmf_champ%hgrid.eq.hgrid)) then
              champ_courant => champ_recherche
              trouve=.true.
           endif
           if(trouve) exit

! If not, let see the next record

           if(associated(champ_recherche%champ_suivant)) then
              champ_recherche=>champ_recherche%champ_suivant
           else
              fin_liste=.true.
           endif
        enddo
     endif

     if(.not.trouve) then
        bmf_length=bmf_length+1
        allocate(champ_courant,STAT=status)
        if (status > 0) STOP 'ERREUR ALLOCATION NOUVEAU CHAMP'
     endif

! Beginning of the storage

     champ_courant%bmf_champ%hgrid=hgrid
     champ_courant%bmf_champ%vgrid=vgrid
     champ_courant%bmf_champ%scat=scat

! Some checks for the dimensions

     if(trouve) then

        champ_courant%bmf_champ%istart=min(istart,champ_courant%bmf_champ%istart)
        champ_courant%bmf_champ%iend=max(iend,champ_courant%bmf_champ%iend)
        champ_courant%bmf_champ%jstart=min(jstart,champ_courant%bmf_champ%jstart)
        champ_courant%bmf_champ%jend=max(jend,champ_courant%bmf_champ%jend)
        champ_courant%bmf_champ%kstart=min(kstart,champ_courant%bmf_champ%kstart)
        champ_courant%bmf_champ%kend=max(kend,champ_courant%bmf_champ%kend)

        if(ni.ne.champ_courant%bmf_champ%ni) then
           write(*,*) 'BMF_gobe: inconsistent NI!'
           length=-1
           return
        endif

        if(nj.ne.champ_courant%bmf_champ%nj) then
           write(*,*) 'BMF_gobe: inconsistent NJ!'
           length=-1
           return
        endif

        if(nk.ne.champ_courant%bmf_champ%nk) then
           write(*,*) 'BMF_gobe: inconsistant NK!'
           length=-1
           return
        endif

        if(dtyp.ne.champ_courant%bmf_champ%dtyp) then
           write(*,*) 'BMF_gobe: inconsistent DTYP!'
           length=-1
           return
        endif

        champ_courant%bmf_champ%ndata=champ_courant%bmf_champ%ndata+ndata

     else

        champ_courant%bmf_champ%nom=nom_tempo
        champ_courant%bmf_champ%time1=time1
        champ_courant%bmf_champ%time2=time2
        champ_courant%bmf_champ%istart=istart
        champ_courant%bmf_champ%iend=iend
        champ_courant%bmf_champ%jstart=jstart
        champ_courant%bmf_champ%jend=jend
        champ_courant%bmf_champ%kstart=kstart
        champ_courant%bmf_champ%kend=kend
        champ_courant%bmf_champ%ni=ni
        champ_courant%bmf_champ%nj=nj
        champ_courant%bmf_champ%nk=nk
        champ_courant%bmf_champ%dtyp=dtyp
        champ_courant%bmf_champ%ndata=ndata

     endif

     read(iun,err = 100, end = 9100 ) nbyt_data
     if(dtyp.ne.bmf_character) then
        allocate(tab_temp(ndata*size),stat=ierr)
        read(iun,err = 100, end = 9100 ) tab_temp
     else
        allocate(tab_temp(2+(ndata-1)/sizeofint))
        read(iun,err = 100, end = 9100 ) tab_temp
     endif
     read(iun,err = 100, end = 9100 ) nbyt_data2

     if(nbyt_data.ne.nbyt_data2) STOP 'probleme avec les donnees, nom=' 

     if(scat.eq.0) then
        if(.not.trouve) then
           allocate(champ_courant%bmf_champ%tableau(size*ni,nj,nk),stat=ierr)
           champ_courant%bmf_champ%tableau=-1
        endif
        indice=0
        if(dtyp.ne.bmf_character) then
           do k=kstart,kend
           do j=jstart,jend
           do i=size*(istart-1)+1,size*iend
              indice=indice+1
              champ_courant%bmf_champ%tableau(i,j,k)=tab_temp(indice)
           enddo
           enddo
           enddo
        else
           do i=1,2+(ndata-1)/sizeofint
              champ_courant%bmf_champ%tableau(i,1,1)=tab_temp(i)
           enddo
        endif
     else if(scat.eq.1) then
        ! hole mode, stream
        if(.not.trouve) then
           allocate(champ_courant%bmf_champ%tableau(size*nbyt_data,1,1),stat=ierr)
           champ_courant%bmf_champ%tableau=-1
        endif
        do i=1,size*nbyt_data
           champ_courant%bmf_champ%tableau(i,1,1)=tab_temp(i)
        enddo
     else
        write(*,*) 'scat mode unimplemented yet'
        stop
     endif

   deallocate(tab_temp)

   if(.not.trouve) then
     champ_courant%champ_suivant=>liste
     liste=> champ_courant
   else
     champ_courant=>liste
     liste=> champ_courant
   endif

    read(iun,err = 100, end = 9200 ) nbyt_head,nom_tempo
  enddo

 9200  continue

!  write(*,*) 'BMF_GOBE: Fin fichier'

   length=bmf_length
   ierr = fclos(iun)

   return
 9100  write(*,*) 'BMF_GOBE: Fin prematuree fichier'
       length=-1
       stop
       return
 100   write(*,*) 'BMF_GOBE: erreur lecture'
       length=-1
       stop
       return

end function bmf_gobe
